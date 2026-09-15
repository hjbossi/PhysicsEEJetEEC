//=====================================================================-*-C++-*-
// ALEPH EEC unfolding with projected statistical covariance.
//
//   1. Nominal unfolding is done once with RooUnfoldBayes (kNoError) for the
//      central values, and a sparse reimplementation of the same D'Agostini
//      iteration is validated against it.
//   2. The exact adjoint of that iteration builds the composite operator
//      G = W * J, where W is the E1E2 weighted projection onto the angular axis
//      and J is the Jacobian of the unfolding. G has only kTotalBins rows, so
//      the flattened measured covariance is never constructed.
//   3. The statistical covariance of the 1D EEC is accumulated event by event,
//      with events as the independent statistical unit. Within event
//      correlations are preserved because each event contributes one short
//      vector g_e = G y_e.
//   4. Pre unfolding covariance is the same code path with G replaced by W.
//   5. A fully nonlinear Poisson(1) event bootstrap runs through the sparse
//      unfolder, parallelised over replicas.
//
// Authors: Tim Adye, Fergus Wilson (original RooUnfold example)
// Adapted by Hannah Bossi for the ALEPH EEC analysis.
//==============================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTree.h"

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include "EffCorrFactor.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//==============================================================================
// Configuration
//==============================================================================

#define MAXPAIR 10000

// The original macro computed a tracking efficiency weight 1/(eff1*eff2) but
// never applied it to any fill. That behaviour is preserved by default.
static const bool kApplyTrackingEfficiency = false;

// Memory strategy for the per event pair cache.
//
// If the fill weight of a pair is a pure function of its flattened cell, which
// is the case when the only weight is a binned fake correction, then the weight
// does not need to be stored per pair. It is folded once into the columns of the
// projection operators and the cache holds nothing but 4 byte cell indices.
// This is a factor 4 reduction, roughly 800 MB instead of 7 GB at 1e8 pairs per
// observable, and it removes a multiply from the inner loop.
//
// Set kForcePerPairWeights to true if the weight genuinely varies within a cell,
// which happens as soon as kApplyTrackingEfficiency is enabled. The code checks
// the assumption at run time and refuses to continue if it is violated.
static const bool kForcePerPairWeights = true;

// Relative tolerance for the within cell weight consistency check.
static const double kCellWeightTolerance = 1e-12;

// How many offending EffCorrFactor queries to print in full.
static const int kMaxWeightComplaints = 20;

//==============================================================================
// Binning helper
//==============================================================================

// 0 based bin index, -1 for underflow and nBins for overflow. The loop must
// inspect bins[nBins] as well, otherwise the last real bin is misclassified.
static int FindBin(double value, int nBins, const std::vector<double> &bins)
{
   for (int i = 0; i <= nBins; ++i)
      if (value < bins[i])
         return i - 1;
   return nBins;
}

//==============================================================================
// Sparse smearing matrix and the D'Agostini iteration
//
// Flattening, identical on the truth and measured sides:
//     flat = angularBin * nE + e1e2Bin
//==============================================================================

struct SparseResponse {
   int nT = 0;
   int nM = 0;
   std::vector<int>    tIdx;
   std::vector<int>    mIdx;
   std::vector<double> P;      // P(m|t) = R(m,t) / N_truth(t)
   std::vector<double> eff;    // eff[t] = sum_m P(m|t)
   std::vector<double> prior;  // MC truth counts, the initial Bayes prior

   void Fold(const std::vector<double> &x, std::vector<double> &out) const
   {
      out.assign(nM, 0.0);
      for (size_t k = 0; k < P.size(); ++k)
         out[mIdx[k]] += P[k] * x[tIdx[k]];
   }

   void FoldTranspose(const std::vector<double> &y, std::vector<double> &out) const
   {
      out.assign(nT, 0.0);
      for (size_t k = 0; k < P.size(); ++k)
         out[tIdx[k]] += P[k] * y[mIdx[k]];
   }
};

struct BayesTrace {
   std::vector<std::vector<double>> n;   // n[0] = prior, n[k] = result after iteration k
   std::vector<std::vector<double>> f;   // f[k-1] = P n[k-1]
};

// n^k_t = (n^{k-1}_t / eff_t) * sum_m P(m|t) d_m / f^{k-1}_m
static void UnfoldBayesSparse(const SparseResponse &R, const std::vector<double> &d, int nIter,
                              std::vector<double> &result, BayesTrace *trace)
{
   std::vector<double> nCur = R.prior;
   if (trace) {
      trace->n.clear();
      trace->f.clear();
      trace->n.push_back(nCur);
   }

   std::vector<double> f, ratio, back, nNew;
   for (int k = 1; k <= nIter; ++k) {
      R.Fold(nCur, f);
      if (trace)
         trace->f.push_back(f);

      ratio.assign(R.nM, 0.0);
      for (int m = 0; m < R.nM; ++m)
         if (f[m] > 0.0)
            ratio[m] = d[m] / f[m];

      R.FoldTranspose(ratio, back);

      nNew.assign(R.nT, 0.0);
      for (int t = 0; t < R.nT; ++t)
         if (R.eff[t] > 0.0)
            nNew[t] = nCur[t] / R.eff[t] * back[t];

      nCur.swap(nNew);
      if (trace)
         trace->n.push_back(nCur);
   }
   result = nCur;
}

// gRow = u^T J with J = dn^K/dd, including the prior derivative at every step.
// Forward recursion (Adye 2011): J^k = M^k + S^k J^{k-1}, J^0 = 0, with
//     M^k_{t,m} = P(m|t) n^{k-1}_t / (eff_t f^{k-1}_m)
//     S^k       = diag(n^k / n^{k-1}) - diag(n^{k-1}/eff) P^T diag(d/f^2) P
// Transposed: u^T J^K = sum_k (v^k)^T M^k, v^K = u, v^{k-1} = (S^k)^T v^k.
static void AdjointRowBayes(const SparseResponse &R, const std::vector<double> &d, const BayesTrace &tr,
                            int nIter, const std::vector<double> &u, std::vector<double> &gRow)
{
   gRow.assign(R.nM, 0.0);
   std::vector<double> v = u, a, y, cy, z;

   for (int k = nIter; k >= 1; --k) {
      const std::vector<double> &nPrev = tr.n[k - 1];
      const std::vector<double> &nCurr = tr.n[k];
      const std::vector<double> &f     = tr.f[k - 1];

      a.assign(R.nT, 0.0);
      for (int t = 0; t < R.nT; ++t)
         if (R.eff[t] > 0.0)
            a[t] = nPrev[t] / R.eff[t] * v[t];

      R.Fold(a, y);

      for (int m = 0; m < R.nM; ++m)
         if (f[m] > 0.0)
            gRow[m] += y[m] / f[m];

      if (k == 1)
         break; // J^0 vanishes

      cy.assign(R.nM, 0.0);
      for (int m = 0; m < R.nM; ++m)
         if (f[m] > 0.0)
            cy[m] = d[m] * y[m] / (f[m] * f[m]);

      R.FoldTranspose(cy, z);

      for (int t = 0; t < R.nT; ++t) {
         const double r = (nPrev[t] > 0.0) ? nCurr[t] / nPrev[t] : 0.0;
         v[t] = r * v[t] - z[t];
      }
   }
}

//==============================================================================
// Projection helpers
//==============================================================================

static std::vector<double> ProjectVector(const std::vector<double> &x, const std::vector<double> &centers,
                                         int nX, int nE)
{
   std::vector<double> out(nX, 0.0);
   for (int i = 0; i < nX; ++i) {
      double s = 0.0;
      for (int j = 0; j < nE; ++j)
         s += x[i * nE + j] * centers[j];
      out[i] = s;
   }
   return out;
}

static void WeightedProjection(const TH2D *h2D, TH1D *h1D)
{
   h1D->Reset();
   for (int i = 1; i <= h2D->GetNbinsX(); ++i) {
      double sum = 0.0;
      for (int j = 1; j <= h2D->GetNbinsY(); ++j)
         sum += h2D->GetBinContent(i, j) * h2D->GetYaxis()->GetBinCenter(j);
      h1D->SetBinContent(i, sum);
      h1D->SetBinError(i, 0.0);
   }
}

//==============================================================================
// One observable, either theta or z
//==============================================================================

struct Observable {
   std::string name;
   std::vector<double> xBins;
   TH2D *hRaw = nullptr;
   TH2D *hSmeared = nullptr;
   TH2D *hTrue = nullptr;
   TH1D *hClosure = nullptr;
   RooUnfoldResponse response;
   SparseResponse sparse;
   EffCorrFactor *fake = nullptr;

   // Per event pair cache in CSR layout. cells[evtOffset[e] .. evtOffset[e+1])
   // are the flattened measured cells filled by event e.
   std::vector<int>    cells;
   std::vector<size_t> evtOffset;
   std::vector<float>  pairWeights;  // parallel to cells, only if kForcePerPairWeights

   std::vector<double> cellWeight;   // per cell fill weight, 0 for never filled cells
   std::vector<double> dNominal;

   std::vector<double> Wt, Gt;       // raw operators, column major op[flat * nRows + row]
   std::vector<double> WtEff, GtEff; // same, with cellWeight folded into each column

   std::vector<std::unordered_map<int, double>> Racc;
   std::vector<double> truthCount;

   // EffCorrFactor diagnostics
   long long nBadWeight = 0;
   double    badWeightSum = 0.0;
   int       nComplaints = 0;
   long long nCellMismatch = 0;
   double    maxCellMismatch = 0.0;
   double    minWeight =  std::numeric_limits<double>::max();
   double    maxWeight = -std::numeric_limits<double>::max();
};

//==============================================================================
// Event level moment accumulation
//
// sumG  = sum_e g_e,  sumGG = sum_e g_e g_e^T,  with g_e = op * y_e.
//==============================================================================

static void AccumulateEventMoments(const Observable &o, const std::vector<double> &op, int nRows,
                                   bool usePairWeights, std::vector<double> &sumG,
                                   std::vector<double> &sumGG)
{
   sumG.assign(nRows, 0.0);
   sumGG.assign((size_t)nRows * nRows, 0.0);
   const int nEvents = (int)o.evtOffset.size() - 1;

   auto processEvent = [&](int e, std::vector<double> &g, std::vector<double> &sg, std::vector<double> &sgg) {
      std::fill(g.begin(), g.end(), 0.0);
      const size_t kBegin = o.evtOffset[e], kEnd = o.evtOffset[e + 1];
      for (size_t k = kBegin; k < kEnd; ++k) {
         const double *col = &op[(size_t)o.cells[k] * nRows];
         if (usePairWeights) {
            const double w = o.pairWeights[k];
            for (int i = 0; i < nRows; ++i)
               g[i] += w * col[i];
         } else {
            for (int i = 0; i < nRows; ++i)
               g[i] += col[i];
         }
      }
      for (int i = 0; i < nRows; ++i) {
         if (g[i] == 0.0)
            continue;
         sg[i] += g[i];
         double *row = &sgg[(size_t)i * nRows];
         for (int j = 0; j <= i; ++j)
            row[j] += g[i] * g[j];
      }
   };

#ifdef _OPENMP
#pragma omp parallel
   {
      std::vector<double> g(nRows, 0.0), locSumG(nRows, 0.0), locSumGG((size_t)nRows * nRows, 0.0);
#pragma omp for schedule(static)
      for (int e = 0; e < nEvents; ++e)
         processEvent(e, g, locSumG, locSumGG);
#pragma omp critical
      {
         for (int i = 0; i < nRows; ++i)
            sumG[i] += locSumG[i];
         for (size_t k = 0; k < sumGG.size(); ++k)
            sumGG[k] += locSumGG[k];
      }
   }
#else
   std::vector<double> g(nRows, 0.0);
   for (int e = 0; e < nEvents; ++e)
      processEvent(e, g, sumG, sumGG);
#endif

   for (int i = 0; i < nRows; ++i)
      for (int j = i + 1; j < nRows; ++j)
         sumGG[(size_t)i * nRows + j] = sumGG[(size_t)j * nRows + i];
}

//==============================================================================
// Histogram writers
//==============================================================================

static TH2D *MakeCovHist(const std::string &name, const std::vector<double> &cov, int nRows)
{
   TH2D *h = new TH2D(name.c_str(), (name + ";bin index;bin index").c_str(), nRows, 0, nRows, nRows, 0, nRows);
   for (int i = 0; i < nRows; ++i)
      for (int j = 0; j < nRows; ++j)
         h->SetBinContent(i + 1, j + 1, cov[(size_t)i * nRows + j]);
   return h;
}

static TH2D *MakeCorrHist(const std::string &name, const std::vector<double> &cov, int nRows)
{
   TH2D *h = new TH2D(name.c_str(), (name + ";bin index;bin index").c_str(), nRows, 0, nRows, nRows, 0, nRows);
   for (int i = 0; i < nRows; ++i) {
      const double cii = cov[(size_t)i * nRows + i];
      for (int j = 0; j < nRows; ++j) {
         const double cjj = cov[(size_t)j * nRows + j];
         if (cii > 0.0 && cjj > 0.0)
            h->SetBinContent(i + 1, j + 1, cov[(size_t)i * nRows + j] / std::sqrt(cii * cjj));
      }
   }
   return h;
}

// If `edges` is provided (size nRows+1, physical bin edges matching the
// observable's xBins), the error is divided by the physical bin width
// Delta_i = edges[i+1]-edges[i]. If `nEvents` is provided (!= 1.0), the error
// is additionally divided by it. Together these reproduce EXACTLY the
// sigma / (width * nEvents) convention used for stat_error_Z/stat_error_Theta
// in preUnfoldingCovariance.cxx: since sqrt(cov) and the central value must
// be rescaled by the identical factor to remain comparable, and the
// non-uniform double-log binning here varies in width by orders of magnitude
// across the range. Passing edges=nullptr and/or leaving nEvents at its
// default of 1.0 reproduces the previous (un-normalized) behaviour exactly.
static TH1D *MakeErrorHist(const std::string &name, const std::vector<double> &cov, int nRows,
                           const std::vector<double> *edges = nullptr, double nEvents = 1.0)
{
   if (edges && (int)edges->size() != nRows + 1) {
      std::cout << "MakeErrorHist: edges size " << edges->size() << " does not match nRows+1 = " << (nRows + 1)
                << " for " << name << " -- falling back to un-normalized error." << std::endl;
      edges = nullptr;
   }

   TH1D *h = new TH1D(name.c_str(), (name + ";bin index;stat. uncertainty").c_str(), nRows, 0, nRows);
   for (int i = 0; i < nRows; ++i) {
      const double c   = cov[(size_t)i * nRows + i];
      double       err = c > 0.0 ? std::sqrt(c) : 0.0;
      if (edges) {
         const double width = (*edges)[i + 1] - (*edges)[i];
         err                = (width > 0.0) ? err / width : 0.0;
      }
      if (nEvents != 1.0 && nEvents > 0.0)
         err /= nEvents;
      h->SetBinContent(i + 1, err);
   }
   return h;
}

static TH1D *MakeValueHist(const std::string &name, const std::vector<double> &v, int nRows)
{
   TH1D *h = new TH1D(name.c_str(), (name + ";bin index;EEC").c_str(), nRows, 0, nRows);
   for (int i = 0; i < nRows; ++i)
      h->SetBinContent(i + 1, v[i]);
   return h;
}

//==============================================================================
// Main routine
//==============================================================================

void RooEEC_ProjectedCovariance(std::string date = "09122026", int iter = 4,
                                int nNonlinearReplicas = 200, UInt_t seed = 42,
                                bool writeOperator = false)
{
   //---------------------------------------------------------------
   // Binning
   //---------------------------------------------------------------
   const int BinCount   = 100;
   const int kTotalBins = 2 * BinCount;

   const double BinMin  = 0.002;
   const double BinMax  = M_PI / 2;
   const double zBinMin = (1 - cos(0.002)) / 2;
   const double zBinMax = 0.5;

   std::vector<double> thetaBins(kTotalBins + 1, 0.0), zBins(kTotalBins + 1, 0.0);
   for (int i = 0; i <= BinCount; ++i) {
      thetaBins[i]                = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      thetaBins[2 * BinCount - i] = BinMax * 2 - thetaBins[i];
      zBins[i]                    = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2 * BinCount - i]     = zBinMax * 2 - zBins[i];
   }

   const std::vector<double> e1e2Bins = {0.0,     0.0001,  0.0002,  0.0005, 0.00075, 0.001, 0.00125, 0.0015,
                                         0.00175, 0.002,   0.00225, 0.0025, 0.00275, 0.003, 0.0035,  0.004,
                                         0.005,   0.007,   0.01,    0.02,   0.03,    0.04,  0.05,    0.07,
                                         0.10,    0.15,    0.20,    0.3};
   const int nE = (int)e1e2Bins.size() - 1;
   const int nM = kTotalBins * nE;
   const int nT = nM;

   std::cout << "Flattened space: " << kTotalBins << " angular x " << nE << " E1E2 = " << nM << " cells"
             << std::endl;
   std::cout << "Pair weight storage: " << (kForcePerPairWeights ? "per pair" : "folded per cell") << std::endl;
#ifdef _OPENMP
   std::cout << "OpenMP threads: " << omp_get_max_threads() << std::endl;
#else
   std::cout << "OpenMP: disabled at compile time" << std::endl;
#endif

   //---------------------------------------------------------------
   // Observables
   //---------------------------------------------------------------
   EffCorrFactor fakeZ, fakeTheta;
   fakeZ.init("/home/hbossi/PhysicsEEJetEEC/Unfolding/20250317_Unfolding/matchingScheme2/FakeCorr.root", "z");
   fakeTheta.init("/home/hbossi/PhysicsEEJetEEC/Unfolding/20250317_Unfolding/matchingScheme2/FakeCorr.root",
                  "theta");

   Observable obsTheta, obsZ;
   obsTheta.name  = "Theta";
   obsTheta.xBins = thetaBins;
   obsTheta.fake  = &fakeTheta;
   obsZ.name      = "Z";
   obsZ.xBins     = zBins;
   obsZ.fake      = &fakeZ;

   std::vector<Observable *> obs = {&obsTheta, &obsZ};

   for (Observable *o : obs) {
      const std::string &n = o->name;
      o->hRaw     = new TH2D(("r_" + n).c_str(), ("raw_" + n).c_str(), kTotalBins, 0, kTotalBins, nE,
                             e1e2Bins.data());
      o->hSmeared = new TH2D(("smeared_" + n).c_str(), ("smeared_" + n).c_str(), kTotalBins, 0, kTotalBins, nE,
                             e1e2Bins.data());
      o->hTrue    = new TH2D(("true_" + n).c_str(), ("true_" + n).c_str(), kTotalBins, 0, kTotalBins, nE,
                             e1e2Bins.data());
      o->hClosure = new TH1D(("h1MCGen_" + n).c_str(), ("h1MCGen_" + n).c_str(), kTotalBins, 0, kTotalBins);
      o->hRaw->Sumw2();
      o->hSmeared->Sumw2();
      o->hTrue->Sumw2();
      o->response.Setup(o->hSmeared, o->hTrue);
      o->dNominal.assign(nM, 0.0);
      o->cellWeight.assign(nM, 0.0);
      o->Racc.assign(nT, {});
      o->truthCount.assign(nT, 0.0);
   }

   // E1E2 centres taken from the axis, so the projection operator is bit for
   // bit the one used by WeightedProjection.
   std::vector<double> centers(nE, 0.0);
   for (int j = 0; j < nE; ++j)
      centers[j] = obsTheta.hTrue->GetYaxis()->GetBinCenter(j + 1);

   //---------------------------------------------------------------
   // Data pass
   //---------------------------------------------------------------
   TFile *inputData = TFile::Open("UnfoldingInputData_03192025.root");
   TTree *smeared   = (TTree *)inputData->Get("UnmatchedPairTree");
   const Int_t nEv  = smeared->GetEntries();
   std::cout << "Data events: " << nEv << std::endl;

   // Pre scan of the pair multiplicity only, so the CSR cache can be reserved
   // exactly and never reallocates a multi hundred MB buffer.
   {
      Int_t nPairsScan = 0;
      smeared->SetBranchStatus("*", 0);
      smeared->SetBranchStatus("NUnmatchedPair", 1);
      smeared->SetBranchAddress("NUnmatchedPair", &nPairsScan);
      long long totalPairs = 0;
      for (Int_t iEntry = 0; iEntry < nEv; ++iEntry) {
         smeared->GetEntry(iEntry);
         totalPairs += nPairsScan;
      }
      smeared->ResetBranchAddresses();
      smeared->SetBranchStatus("*", 1);
      std::cout << "Total data pairs: " << totalPairs << std::endl;
      for (Observable *o : obs) {
         o->cells.reserve(totalPairs);
         if (kForcePerPairWeights)
            o->pairWeights.reserve(totalPairs);
         o->evtOffset.reserve(nEv + 1);
         o->evtOffset.push_back(0);
      }
   }

   Double_t e1e2Data[MAXPAIR], thetaData[MAXPAIR];
   Double_t eff1[MAXPAIR], eff2[MAXPAIR], recoE1Data[MAXPAIR], recoE2Data[MAXPAIR];
   Int_t    nPairsData = 0;

   smeared->SetBranchAddress("NUnmatchedPair", &nPairsData);
   smeared->SetBranchAddress("E1E2RecoUnmatched", &e1e2Data);
   smeared->SetBranchAddress("DistanceUnmatchedReco", &thetaData);
   smeared->SetBranchAddress("RecoE1Unmatched", &recoE1Data);
   smeared->SetBranchAddress("RecoE2Unmatched", &recoE2Data);
   smeared->SetBranchAddress("RecoEfficiency1", &eff1);
   smeared->SetBranchAddress("RecoEfficiency2", &eff2);

   std::vector<char> cellSet[2] = {std::vector<char>(nM, 0), std::vector<char>(nM, 0)};
   long long droppedData = 0;

   for (Int_t iEntry = 0; iEntry < nEv; ++iEntry) {
      smeared->GetEntry(iEntry);

      for (int i = 0; i < nPairsData; ++i) {
         const int binE = obsTheta.hRaw->GetYaxis()->FindBin(e1e2Data[i]);
         if (binE < 1 || binE > nE) {
            droppedData += 2;
            continue;
         }

         double trackingWeight = 1.0;
         if (kApplyTrackingEfficiency) {
            const double e1 = std::min(eff1[i], 1.0), e2 = std::min(eff2[i], 1.0);
            if (e1 > 0.0 && e2 > 0.0)
               trackingWeight = 1.0 / (e1 * e2);
         }

         const double z    = (1 - cos(thetaData[i])) / 2;
         const int binX[2] = {FindBin(thetaData[i], kTotalBins, thetaBins),
                              FindBin(z, kTotalBins, zBins)};

         for (int k = 0; k < 2; ++k) {
            Observable *o = obs[k];
            if (binX[k] < 0 || binX[k] >= kTotalBins) {
               ++droppedData;
               continue;
            }

            const double fakeW = o->fake->efficiency(binX[k], e1e2Data[i]);
            const double w     = trackingWeight * fakeW;

            // EffCorrFactor diagnostics. A non finite or non positive
            // correction silently corrupts both the central value and the
            // covariance, so it is counted and reported rather than ignored.
            if (!std::isfinite(fakeW) || fakeW <= 0.0) {
               ++o->nBadWeight;
               o->badWeightSum += std::isfinite(w) ? std::fabs(w) : 0.0;
               if (o->nComplaints < kMaxWeightComplaints) {
                  std::cout << "  BAD WEIGHT " << o->name << ": binX = " << binX[k]
                            << ", E1E2 = " << e1e2Data[i] << ", E1E2 bin = " << binE
                            << ", efficiency() = " << fakeW << std::endl;
                  ++o->nComplaints;
               }
            }
            o->minWeight = std::min(o->minWeight, fakeW);
            o->maxWeight = std::max(o->maxWeight, fakeW);

            const int flat = binX[k] * nE + (binE - 1);

            // Within cell weight consistency, the assumption that lets the
            // weight be folded into the operator columns.
            if (!cellSet[k][flat]) {
               cellSet[k][flat]   = 1;
               o->cellWeight[flat] = w;
            } else if (std::fabs(w - o->cellWeight[flat]) >
                       kCellWeightTolerance * std::max(1.0, std::fabs(o->cellWeight[flat]))) {
               ++o->nCellMismatch;
               o->maxCellMismatch =
                   std::max(o->maxCellMismatch,
                            std::fabs(w - o->cellWeight[flat]) / std::max(1e-300, std::fabs(o->cellWeight[flat])));
            }

            o->hRaw->Fill(binX[k], e1e2Data[i], w);
            o->cells.push_back(flat);
            if (kForcePerPairWeights)
               o->pairWeights.push_back((float)w);
            o->dNominal[flat] += w;
         }
      }

      for (Observable *o : obs)
         o->evtOffset.push_back(o->cells.size());
   }

   std::cout << "Raw integrals, theta " << obsTheta.hRaw->Integral() << ", z " << obsZ.hRaw->Integral()
             << ", pair fills dropped as out of range " << droppedData << std::endl;

   //---------------------------------------------------------------
   // EffCorrFactor report and the folding guard
   //---------------------------------------------------------------
   bool foldingValid = true;
   for (Observable *o : obs) {
      double totalWeight = 0.0;
      for (int m = 0; m < nM; ++m)
         totalWeight += o->dNominal[m];

      std::cout << o->name << " fake correction: range [" << o->minWeight << ", " << o->maxWeight
                << "], bad queries " << o->nBadWeight;
      if (o->nBadWeight > 0)
         std::cout << " (affected weight fraction " << o->badWeightSum / totalWeight << ")";
      std::cout << ", cached pairs " << o->cells.size() << std::endl;

      if (o->nBadWeight > 0)
         std::cout << "  WARNING: EffCorrFactor returned non finite or non positive values. This affects the "
                      "central values as well as the covariance. Locate the print statement in "
                      "EffCorrFactor and fix the out of range behaviour."
                   << std::endl;

      if (o->nCellMismatch > 0) {
         foldingValid = false;
         std::cout << "  ERROR: " << o->nCellMismatch
                   << " pairs disagree with the cell weight, max relative deviation " << o->maxCellMismatch
                   << ". The fill weight is not a pure function of the cell, so it cannot be folded into the "
                      "operator. Set kForcePerPairWeights = true and rerun."
                   << std::endl;
      }
   }
   if (!kForcePerPairWeights && !foldingValid) {
      std::cout << "Aborting before the covariance is built, since folding would give a wrong answer."
                << std::endl;
      inputData->Close();
      return;
   }

   //---------------------------------------------------------------
   // MC pass, building the RooUnfold response and the sparse copy
   //---------------------------------------------------------------
   TFile *inputMC = TFile::Open(
       "/home/hbossi/PhysicsEEJetEEC/Unfolding/20250317_Unfolding/matchingScheme2/LEP1MCMerged_Matched.root");
   TTree *mc = (TTree *)inputMC->Get("PairTree");

   Double_t e1e2Reco[MAXPAIR], e1e2Gen[MAXPAIR], thetaReco[MAXPAIR], thetaGen[MAXPAIR];
   Double_t recoE1[MAXPAIR], recoE2[MAXPAIR];
   Int_t    nPairsMC = 0;

   mc->SetBranchAddress("NPair", &nPairsMC);
   mc->SetBranchAddress("E1E2Reco", &e1e2Reco);
   mc->SetBranchAddress("E1E2Gen", &e1e2Gen);
   mc->SetBranchAddress("DistanceReco", &thetaReco);
   mc->SetBranchAddress("DistanceGen", &thetaGen);
   mc->SetBranchAddress("RecoE1", &recoE1);
   mc->SetBranchAddress("RecoE2", &recoE2);

   const Int_t nEvMC = mc->GetEntries();
   std::cout << "MC events: " << nEvMC << std::endl;

   for (Int_t iEntry = 0; iEntry < nEvMC; ++iEntry) {
      mc->GetEntry(iEntry);
      for (int i = 0; i < nPairsMC; ++i) {
         if (recoE1[i] < 0 || recoE2[i] < 0)
            continue;
         if (thetaReco[i] < BinMin || thetaGen[i] < BinMin)
            continue;

         const int binEReco = obsTheta.hTrue->GetYaxis()->FindBin(e1e2Reco[i]);
         const int binEGen  = obsTheta.hTrue->GetYaxis()->FindBin(e1e2Gen[i]);

         const double zReco = (1 - cos(thetaReco[i])) / 2;
         const double zGen  = (1 - cos(thetaGen[i])) / 2;

         const int binXReco[2] = {FindBin(thetaReco[i], kTotalBins, thetaBins),
                                  FindBin(zReco, kTotalBins, zBins)};
         const int binXGen[2]  = {FindBin(thetaGen[i], kTotalBins, thetaBins),
                                  FindBin(zGen, kTotalBins, zBins)};

         for (int k = 0; k < 2; ++k) {
            Observable *o = obs[k];
            o->hTrue->Fill(binXGen[k], e1e2Gen[i]);
            o->hClosure->Fill(binXGen[k], e1e2Gen[i]);
            o->hSmeared->Fill(binXReco[k], e1e2Reco[i]);
            o->response.Fill(binXReco[k], e1e2Reco[i], binXGen[k], e1e2Gen[i]);

            // Truth normalisation counts every generated pair with a valid
            // truth cell, so pairs whose measured cell is out of range act as
            // misses and lower the efficiency, matching RooUnfold.
            if (binXGen[k] < 0 || binXGen[k] >= kTotalBins || binEGen < 1 || binEGen > nE)
               continue;
            const int tFlat = binXGen[k] * nE + (binEGen - 1);
            o->truthCount[tFlat] += 1.0;

            if (binXReco[k] < 0 || binXReco[k] >= kTotalBins || binEReco < 1 || binEReco > nE)
               continue;
            o->Racc[tFlat][binXReco[k] * nE + (binEReco - 1)] += 1.0;
         }
      }
   }

   // Finalise the sparse responses.
   for (Observable *o : obs) {
      SparseResponse &R = o->sparse;
      R.nT = nT;
      R.nM = nM;
      R.eff.assign(nT, 0.0);
      R.prior.assign(nT, 0.0);
      size_t nnz = 0;
      for (int t = 0; t < nT; ++t)
         nnz += o->Racc[t].size();
      R.tIdx.reserve(nnz);
      R.mIdx.reserve(nnz);
      R.P.reserve(nnz);

      for (int t = 0; t < nT; ++t) {
         const double N = o->truthCount[t];
         R.prior[t] = N;
         if (N <= 0.0) {
            o->Racc[t].clear();
            continue;
         }
         for (const auto &entry : o->Racc[t]) {
            const double p = entry.second / N;
            R.tIdx.push_back(t);
            R.mIdx.push_back(entry.first);
            R.P.push_back(p);
            R.eff[t] += p;
         }
         o->Racc[t].clear();
      }
      o->Racc.clear();
      o->Racc.shrink_to_fit();

      // Unsupported cells, quantified by data weight rather than by cell count.
      // Weight is what biases the result, cell count is not.
      std::vector<double> fPrior;
      R.Fold(R.prior, fPrior);
      int    unsupportedCells  = 0;
      double unsupportedWeight = 0.0, totalWeight = 0.0;
      for (int m = 0; m < nM; ++m) {
         totalWeight += o->dNominal[m];
         if (o->dNominal[m] > 0.0 && fPrior[m] <= 0.0) {
            ++unsupportedCells;
            unsupportedWeight += o->dNominal[m];
         }
      }

      std::cout << o->name << " sparse response: " << R.P.size() << " nonzeros, occupancy "
                << (double)R.P.size() / ((double)nT * nM) << std::endl;
      std::cout << "  unsupported cells " << unsupportedCells << ", unsupported data weight fraction "
                << unsupportedWeight / totalWeight << std::endl;
      if (unsupportedWeight / totalWeight > 1e-4)
         std::cout << "  WARNING: this is above 1e-4. Those cells are dropped by the iteration and bias the "
                      "central value. Consider merging MC bins or excluding them explicitly."
                   << std::endl;
   }

   //---------------------------------------------------------------
   // Output file and MC truth projections
   //---------------------------------------------------------------
   TFile *fout = new TFile(Form("unfoldingE2C_ProjectedCovariance_%s.root", date.c_str()), "RECREATE");
   fout->cd();

   for (Observable *o : obs) {
      o->hRaw->Write();
      o->hSmeared->Write();
      o->hTrue->Write();
      o->hClosure->Write();

      TH1D *hTrueProj = new TH1D(("h1True_" + o->name + "_ProjectionX").c_str(),
                                 ("h1True_" + o->name + "_ProjectionX").c_str(), kTotalBins, 0, kTotalBins);
      WeightedProjection(o->hTrue, hTrueProj);
      hTrueProj->Write();
   }

   //---------------------------------------------------------------
   // Nominal unfolding, central values, and the operator G
   //---------------------------------------------------------------
   for (Observable *o : obs) {
      std::cout << "Nominal RooUnfoldBayes for " << o->name << ", iter = " << iter << std::endl;
      TStopwatch t;
      t.Start();
      RooUnfoldBayes unfolder(&o->response, o->hRaw, iter);
      unfolder.SetVerbose(0);
      TH2D *hUnf  = dynamic_cast<TH2D *>(unfolder.Hreco(RooUnfold::kNoError));
      TH1  *hFold = o->response.ApplyToTruth(hUnf, "");
      t.Stop();
      std::cout << "  wall " << t.RealTime() << " s, CPU " << t.CpuTime() << " s" << std::endl;

      ((TH2D *)hUnf->Clone(Form("Bayesian_Unfoldediter%d_%s", iter, o->name.c_str())))->Write();
      ((TH2D *)hFold->Clone(Form("Bayesian_Foldediter%d_%s", iter, o->name.c_str())))->Write();

      TH1D *central = new TH1D(("central_Unfolded_" + o->name).c_str(),
                               ("central_Unfolded_" + o->name + ";bin index;EEC").c_str(), kTotalBins, 0,
                               kTotalBins);
      WeightedProjection(hUnf, central);
      central->Write();

      TH1D *measured = new TH1D(("central_Measured_" + o->name).c_str(),
                                ("central_Measured_" + o->name + ";bin index;EEC").c_str(), kTotalBins, 0,
                                kTotalBins);
      WeightedProjection(o->hRaw, measured);
      measured->Write();

      // Sparse forward pass. Provides the per iteration trace the adjoint
      // consumes, and its projection is compared against RooUnfold to confirm
      // that the prior, efficiency and overflow conventions agree. If they do
      // not, G is built from a different operator than the central values and
      // the covariance is wrong with no other symptom.
      std::vector<double> nOwn;
      BayesTrace trace;
      UnfoldBayesSparse(o->sparse, o->dNominal, iter, nOwn, &trace);

      const std::vector<double> projOwn = ProjectVector(nOwn, centers, kTotalBins, nE);

      double maxRel = 0.0;
      int    maxBin = -1;
      for (int i = 0; i < kTotalBins; ++i) {
         const double ref = central->GetBinContent(i + 1);
         if (ref == 0.0)
            continue;
         const double rel = std::fabs(projOwn[i] - ref) / std::fabs(ref);
         if (rel > maxRel) {
            maxRel = rel;
            maxBin = i;
         }
      }
      std::cout << "  sparse vs RooUnfold: max relative difference " << maxRel << " at bin " << maxBin
                << std::endl;
      if (maxRel > 1e-6)
         std::cout << "  WARNING: sparse iteration does not reproduce RooUnfold. Check the prior, the "
                      "efficiency definition and the overflow handling before using the covariance."
                   << std::endl;

      MakeValueHist("sparseCheck_Unfolded_" + o->name, projOwn, kTotalBins)->Write();

      // W, the weighted projection.
      o->Wt.assign((size_t)nM * kTotalBins, 0.0);
      for (int i = 0; i < kTotalBins; ++i)
         for (int j = 0; j < nE; ++j)
            o->Wt[(size_t)(i * nE + j) * kTotalBins + i] = centers[j];

      // G = W J, one adjoint pass per angular bin.
      o->Gt.assign((size_t)nM * kTotalBins, 0.0);
      TStopwatch tG;
      tG.Start();
      std::vector<double> u(nT, 0.0), gRow;
      for (int i = 0; i < kTotalBins; ++i) {
         std::fill(u.begin(), u.end(), 0.0);
         for (int j = 0; j < nE; ++j)
            u[i * nE + j] = centers[j];
         AdjointRowBayes(o->sparse, o->dNominal, trace, iter, u, gRow);
         for (int m = 0; m < nM; ++m)
            o->Gt[(size_t)m * kTotalBins + i] = gRow[m];
      }
      tG.Stop();
      std::cout << "  built G (" << kTotalBins << " x " << nM << ") in " << tG.RealTime() << " s" << std::endl;

      // Weight folded copies, used by the covariance loop when the cache holds
      // bare cell indices. Kept separate so the raw operators stay available.
      if (!kForcePerPairWeights) {
         o->WtEff = o->Wt;
         o->GtEff = o->Gt;
         for (int m = 0; m < nM; ++m) {
            const double cw = o->cellWeight[m];
            for (int i = 0; i < kTotalBins; ++i) {
               o->WtEff[(size_t)m * kTotalBins + i] *= cw;
               o->GtEff[(size_t)m * kTotalBins + i] *= cw;
            }
         }
      }

      if (writeOperator) {
         TH2D *hG = new TH2D(("operatorG_" + o->name).c_str(), ("operatorG_" + o->name).c_str(), kTotalBins, 0,
                             kTotalBins, nM, 0, nM);
         for (int i = 0; i < kTotalBins; ++i)
            for (int m = 0; m < nM; ++m)
               hG->SetBinContent(i + 1, m + 1, o->Gt[(size_t)m * kTotalBins + i]);
         hG->Write();
      }
   }

   //---------------------------------------------------------------
   // Event level covariance, before and after unfolding
   //---------------------------------------------------------------
   const double N = (double)nEv;
   for (Observable *o : obs) {
      struct Stage {
         std::string tag;
         const std::vector<double> *op;
      };
      const Stage stages[2] = {{"PreUnfolded", kForcePerPairWeights ? &o->Wt : &o->WtEff},
                               {"Unfolded",    kForcePerPairWeights ? &o->Gt : &o->GtEff}};

      for (const Stage &s : stages) {
         TStopwatch t;
         t.Start();
         std::vector<double> sumG, sumGG;
         AccumulateEventMoments(*o, *s.op, kTotalBins, kForcePerPairWeights, sumG, sumGG);
         t.Stop();
         std::cout << o->name << " " << s.tag << " event moments in " << t.RealTime() << " s" << std::endl;

         // Poisson(1) event resampling: Var(sum_e w_e g_e) = sum_e g_e g_e^T,
         // the exact analytic limit of the bootstrap, including fluctuation of
         // the total event count.
         const std::vector<double> &covPoisson = sumGG;

         // Fixed sample size alternative, Bessel corrected. Differs by the rank
         // one term (sum g)(sum g)^T / N, that is by whether the overall
         // normalisation is allowed to fluctuate.
         std::vector<double> covFixedN((size_t)kTotalBins * kTotalBins, 0.0);
         for (int i = 0; i < kTotalBins; ++i)
            for (int j = 0; j < kTotalBins; ++j)
               covFixedN[(size_t)i * kTotalBins + j] =
                   (N * sumGG[(size_t)i * kTotalBins + j] - sumG[i] * sumG[j]) / (N - 1.0);

         MakeCovHist("cov_" + s.tag + "_Poisson_" + o->name, covPoisson, kTotalBins)->Write();
         MakeCorrHist("corr_" + s.tag + "_Poisson_" + o->name, covPoisson, kTotalBins)->Write();
         MakeErrorHist("stat_error_" + s.tag + "_Poisson_" + o->name, covPoisson, kTotalBins, &o->xBins,
                       (double)nEv)
             ->Write();

         MakeCovHist("cov_" + s.tag + "_FixedN_" + o->name, covFixedN, kTotalBins)->Write();
         MakeCorrHist("corr_" + s.tag + "_FixedN_" + o->name, covFixedN, kTotalBins)->Write();
         MakeErrorHist("stat_error_" + s.tag + "_FixedN_" + o->name, covFixedN, kTotalBins, &o->xBins,
                       (double)nEv)
             ->Write();
      }
   }

   //---------------------------------------------------------------
   // Fully nonlinear Poisson(1) event bootstrap through the sparse unfolder,
   // parallelised over replicas. Each replica seeds its own generator from
   // seed + r, so results are reproducible independently of thread count.
   // Set nNonlinearReplicas to 0 or 1 to skip.
   //---------------------------------------------------------------
   if (nNonlinearReplicas > 1) {
      std::cout << "Running " << nNonlinearReplicas << " nonlinear Poisson replicas..." << std::endl;

      const int nObs = (int)obs.size();
      std::vector<std::vector<double>> sumRep(nObs, std::vector<double>(kTotalBins, 0.0));
      std::vector<std::vector<double>> sumRepSq(nObs,
                                                std::vector<double>((size_t)kTotalBins * kTotalBins, 0.0));

      TStopwatch tBoot;
      tBoot.Start();

#ifdef _OPENMP
#pragma omp parallel
#endif
      {
         std::vector<std::vector<double>> cnt(nObs, std::vector<double>(nM, 0.0));
         std::vector<std::vector<double>> locSum(nObs, std::vector<double>(kTotalBins, 0.0));
         std::vector<std::vector<double>> locSumSq(nObs,
                                                   std::vector<double>((size_t)kTotalBins * kTotalBins, 0.0));
         std::vector<double> d(nM, 0.0), nRep;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
         for (int r = 0; r < nNonlinearReplicas; ++r) {
            TRandom3 rng(seed + (UInt_t)r + 1);

            for (int k = 0; k < nObs; ++k)
               std::fill(cnt[k].begin(), cnt[k].end(), 0.0);

            // One Poisson(1) weight per event, applied to every pair from that
            // event and shared by both observables.
            for (Int_t e = 0; e < nEv; ++e) {
               const double w = rng.Poisson(1.0);
               if (w == 0.0)
                  continue;
               for (int k = 0; k < nObs; ++k) {
                  const Observable *o = obs[k];
                  const size_t kBegin = o->evtOffset[e], kEnd = o->evtOffset[e + 1];
                  if (kForcePerPairWeights) {
                     for (size_t q = kBegin; q < kEnd; ++q)
                        cnt[k][o->cells[q]] += w * o->pairWeights[q];
                  } else {
                     for (size_t q = kBegin; q < kEnd; ++q)
                        cnt[k][o->cells[q]] += w;
                  }
               }
            }

            for (int k = 0; k < nObs; ++k) {
               const Observable *o = obs[k];
               if (kForcePerPairWeights) {
                  d = cnt[k];
               } else {
                  for (int m = 0; m < nM; ++m)
                     d[m] = cnt[k][m] * o->cellWeight[m];
               }

               UnfoldBayesSparse(o->sparse, d, iter, nRep, nullptr);
               const std::vector<double> proj = ProjectVector(nRep, centers, kTotalBins, nE);

               for (int i = 0; i < kTotalBins; ++i) {
                  locSum[k][i] += proj[i];
                  for (int j = 0; j <= i; ++j)
                     locSumSq[k][(size_t)i * kTotalBins + j] += proj[i] * proj[j];
               }
            }
         }

#ifdef _OPENMP
#pragma omp critical
#endif
         {
            for (int k = 0; k < nObs; ++k) {
               for (int i = 0; i < kTotalBins; ++i)
                  sumRep[k][i] += locSum[k][i];
               for (size_t q = 0; q < sumRepSq[k].size(); ++q)
                  sumRepSq[k][q] += locSumSq[k][q];
            }
         }
      }

      tBoot.Stop();
      std::cout << "  bootstrap wall time " << tBoot.RealTime() << " s, CPU " << tBoot.CpuTime() << " s"
                << std::endl;

      const double R = (double)nNonlinearReplicas;
      for (int k = 0; k < nObs; ++k) {
         std::vector<double> cov((size_t)kTotalBins * kTotalBins, 0.0), mean(kTotalBins, 0.0);
         for (int i = 0; i < kTotalBins; ++i) {
            mean[i] = sumRep[k][i] / R;
            for (int j = 0; j <= i; ++j) {
               const double c =
                   (sumRepSq[k][(size_t)i * kTotalBins + j] - sumRep[k][i] * sumRep[k][j] / R) / (R - 1.0);
               cov[(size_t)i * kTotalBins + j] = c;
               cov[(size_t)j * kTotalBins + i] = c;
            }
         }
         const std::string n = obs[k]->name;
         MakeValueHist("bootstrapMean_" + n, mean, kTotalBins)->Write();
         MakeCovHist("cov_Unfolded_Bootstrap_" + n, cov, kTotalBins)->Write();
         MakeCorrHist("corr_Unfolded_Bootstrap_" + n, cov, kTotalBins)->Write();
         MakeErrorHist("stat_error_Unfolded_Bootstrap_" + n, cov, kTotalBins, &obs[k]->xBins, (double)nEv)
             ->Write();
      }
   }

   fout->Close();
   inputData->Close();
   inputMC->Close();
}

int main()
{
   RooEEC_ProjectedCovariance("09122026", 4, 200, 42, false);
   return 0;
}