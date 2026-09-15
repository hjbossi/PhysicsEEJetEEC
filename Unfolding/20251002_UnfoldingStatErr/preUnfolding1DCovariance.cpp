//==============================================================================
//  preUnfoldingCovariance.cxx
//
//  Builds the statistical covariance matrices of the reconstructed (smeared)
//  EEC distributions in theta and in z, before unfolding.
//
//  Build:  g++ -O2 -std=c++17 preUnfoldingCovariance.cxx -o preUnfoldingCovariance \
//              $(root-config --cflags --libs)
//==============================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"

#include "EffCorrFactor.h"
#include "alephTrkEfficiency.h"

namespace {

//==============================================================================
//  Configuration
//==============================================================================

constexpr int    kMaxPair  = 10000;   // maximum number of pairs per event
constexpr int    kBinCount = 100;     // half of the (symmetric) angular binning
constexpr int    kTotalBins = 2 * kBinCount;

const TString kInputFileName  = "UnfoldingInputData_03192025.root";
const TString kInputTreeName  = "UnmatchedPairTree";
const TString kOutputFileName = "cov_matrix_09122026.root";
const TString kFakeCorrFile =
    "/home/hbossi/PhysicsEEJetEEC/Unfolding/20250317_Unfolding/matchingScheme2/FakeCorr.root";

//==============================================================================
//  Helper functions
//==============================================================================

// Returns the ROOT-style bin index (1..NBins) of the bin containing Value,
// for the array Bins of NBins+1 monotonically increasing edges
// (Bins[0] .. Bins[NBins]). Returns 0 for underflow and NBins+1 for overflow,
// matching TH1's own SetBinContent/GetBinContent convention.
//
// NOTE: the loop must run over i = 0 .. NBins (inclusive) so that the final
// edge Bins[NBins] is actually inspected. The previous version looped only
// to i < NBins, so any value belonging to the last real bin
// (Bins[NBins-1] <= Value < Bins[NBins]) never satisfied "Value < Bins[i]"
// for any tested i, and fell through to the overflow branch. That silently
// zeroed out the last bin of h2RawTheta/h2RawZ (and therefore of the
// covariance matrices) in every single event.
int FindBin(double Value, int NBins, const double Bins[])
{
   for (int i = 0; i <= NBins; ++i)
      if (Value < Bins[i])
         return i - 1;
   return NBins;
}

// Collapses a (angle, E1E2) histogram onto the angular axis, weighting each
// entry by the E1E2 bin center, and propagates the bin errors in quadrature.
void Projection(const TH2D *h2D, TH1D *h1D)
{
   h1D->Reset();

   for (int i = 1; i <= h2D->GetNbinsX(); ++i) {
      double weight = 0.0;
      double error2 = 0.0;

      for (int j = 1; j <= h2D->GetNbinsY(); ++j) {
         const double binContent = h2D->GetBinContent(i, j);
         const double binError   = h2D->GetBinError(i, j);
         const double binCenter  = h2D->GetYaxis()->GetBinCenter(j);

         weight += binContent * binCenter;
         error2 += (binError * binCenter) * (binError * binCenter);
      }

      h1D->SetBinContent(i, weight);
      h1D->SetBinError(i, std::sqrt(error2));
   }
}

// Fills the symmetric double logarithmic binning used for theta and z.
void BuildBinning(double binMin, double binMax, int binCount, std::vector<double> &bins)
{
   bins.assign(2 * binCount + 1, 0.0);

   const double logMin = std::log(binMin);
   const double logMax = std::log(binMax);

   for (int i = 0; i <= binCount; ++i) {
      const double edge = std::exp(logMin + (logMax - logMin) / binCount * i);
      bins[i]                = edge;
      bins[2 * binCount - i] = binMax * 2.0 - edge;
   }
}

} // namespace

//==============================================================================
//  Main analysis
//==============================================================================

void preUnfoldingCovariance()
{
   // ---------------------------------------------------------------- binning
   const double thetaBinMin = 0.002;
   const double thetaBinMax = M_PI / 2.0;
   const double zBinMin     = (1.0 - std::cos(0.002)) / 2.0;
   const double zBinMax     = 0.5;

   std::vector<double> thetaBins;
   std::vector<double> zBins;
   BuildBinning(thetaBinMin, thetaBinMax, kBinCount, thetaBins);
   BuildBinning(zBinMin, zBinMax, kBinCount, zBins);

   // Binning option #1 for the pair energy weight.
   const std::vector<double> e1e2BinsUnfolded = {
       0.0,     0.0001,  0.0002,  0.0005, 0.00075, 0.001, 0.00125, 0.0015,
       0.00175, 0.002,   0.00225, 0.0025, 0.00275, 0.003, 0.0035,  0.004,
       0.005,   0.007,   0.01,    0.02,   0.03,    0.04,  0.05,    0.07,
       0.10,    0.15,    0.20,    0.3};

   // ------------------------------------------------------ fake corrections
   EffCorrFactor fakeCorrFactorZ;
   fakeCorrFactorZ.init(kFakeCorrFile.Data(), "z");

   EffCorrFactor fakeCorrFactorTheta;
   fakeCorrFactorTheta.init(kFakeCorrFile.Data(), "theta");

   // ---------------------------------------------------------------- input
   TFile *inputSmeared = TFile::Open(kInputFileName);
   if (inputSmeared == nullptr || inputSmeared->IsZombie()) {
      std::cerr << "Error: cannot open " << kInputFileName << std::endl;
      return;
   }

   TTree *smeared = dynamic_cast<TTree *>(inputSmeared->Get(kInputTreeName));
   if (smeared == nullptr) {
      std::cerr << "Error: cannot find tree " << kInputTreeName << std::endl;
      inputSmeared->Close();
      return;
   }

   Int_t    nPairsData = 0;
   Double_t e1e2Data[kMaxPair];
   Double_t thetaData[kMaxPair];
   Double_t eff1[kMaxPair];
   Double_t eff2[kMaxPair];

   smeared->SetBranchAddress("NUnmatchedPair", &nPairsData);
   smeared->SetBranchAddress("E1E2RecoUnmatched", e1e2Data);
   smeared->SetBranchAddress("DistanceUnmatchedReco", thetaData);
   smeared->SetBranchAddress("RecoEfficiency1", eff1);
   smeared->SetBranchAddress("RecoEfficiency2", eff2);

   const Long64_t nEv = smeared->GetEntries();
   std::cout << "Number of entries in the smeared tree: " << nEv << std::endl;

   // ------------------------------------------------------------ containers
   //
   // Bin-index convention used throughout the rest of this file:
   //   - All histograms (h1RawTheta, h1RawZ, covZ, covTheta, ...) use ROOT's
   //     native 1-indexed convention: bin 0 = underflow, bins 1..kTotalBins
   //     = real bins, bin kTotalBins+1 = overflow.
   //   - sumsTheta/sumsZ and tempTheta/tempZ are plain std::vector's, so
   //     they are declared with kTotalBins+1 or kTotalBins+2 entries and are
   //     always indexed with the SAME 1..kTotalBins range as the ROOT bins,
   //     leaving element 0 (and, for sums*, element kTotalBins+1) unused.
   //     This keeps every loop below consistent with GetBinContent/
   //     SetBinContent's own indexing, so there is no off-by-one translation
   //     needed anywhere.
   TH2D *covZ = new TH2D("cov_Z", "cov_Z", kTotalBins, 0, kTotalBins,
                         kTotalBins, 0, kTotalBins);
   TH2D *covTheta = new TH2D("cov_Theta", "cov_Theta", kTotalBins, 0, kTotalBins,
                             kTotalBins, 0, kTotalBins);

   // Sized kTotalBins + 1 (indices 0..kTotalBins) so that ROOT bin indices
   // 1..kTotalBins are all valid and index 0 is simply left unused.
   std::vector<std::vector<double>> tempTheta(kTotalBins + 1, std::vector<double>(kTotalBins + 1, 0.0));
   std::vector<std::vector<double>> tempZ(kTotalBins + 1, std::vector<double>(kTotalBins + 1, 0.0));

   // Sized kTotalBins + 2 so that ROOT style bin indices 1 .. kTotalBins are valid.
   std::vector<double> sumsTheta(kTotalBins + 2, 0.0);
   std::vector<double> sumsZ(kTotalBins + 2, 0.0);

   TH2D *h2RawTheta = new TH2D("r_Theta", "r_Theta", kTotalBins, 0, kTotalBins,
                               static_cast<int>(e1e2BinsUnfolded.size()) - 1,
                               e1e2BinsUnfolded.data());
   TH2D *h2RawZ = new TH2D("r_Z", "r_Z", kTotalBins, 0, kTotalBins,
                           static_cast<int>(e1e2BinsUnfolded.size()) - 1,
                           e1e2BinsUnfolded.data());

   // Created once and reused for every event to avoid leaking one histogram
   // per entry.
   TH1D *h1RawTheta = new TH1D("h1raw_Theta", "h1raw_Theta", kTotalBins, 0, kTotalBins);
   TH1D *h1RawZ     = new TH1D("h1raw_Z", "h1raw_Z", kTotalBins, 0, kTotalBins);
   h1RawTheta->Sumw2();
   h1RawZ->Sumw2();

   // ==================== loop over events ====================
   for (Long64_t iEntry = 0; iEntry < nEv; ++iEntry) {
      smeared->GetEntry(iEntry);

      for (int i = 0; i < nPairsData; ++i) {
         // FindBin now returns a proper 0..kTotalBins-1 "which real bin"
         // index that is then used as the x-coordinate fed into Fill()
         // below, which itself follows TH2::Fill's own (0,kTotalBins)-range
         // axis convention -- consistent with the fix in FindBin() above.
         const int    binTheta = FindBin(thetaData[i], kTotalBins, thetaBins.data());
         const double z        = (1.0 - std::cos(thetaData[i])) / 2.0;
         const int    binZ     = FindBin(z, kTotalBins, zBins.data());

         eff1[i] = std::min(eff1[i], 1.0);
         eff2[i] = std::min(eff2[i], 1.0);

         const double trackingEff = 1.0 / (eff1[i] * eff2[i]);
         if (trackingEff < 1.0)
            std::cout << "Unexpected value of the tracking efficiency " << trackingEff << std::endl;

         const double fakeCorrZ     = fakeCorrFactorZ.efficiency(binZ, e1e2Data[i]);
         const double fakeCorrTheta = fakeCorrFactorTheta.efficiency(binTheta, e1e2Data[i]);

         h2RawTheta->Fill(binTheta, e1e2Data[i], fakeCorrTheta);
         h2RawZ->Fill(binZ, e1e2Data[i], fakeCorrZ);
      }

      // Collapse this event onto the angular axes.
      Projection(h2RawTheta, h1RawTheta);
      Projection(h2RawZ, h1RawZ);

      // Running sums of the per event yields. ROOT bins 1..kTotalBins.
      for (int iBin = 1; iBin <= kTotalBins; ++iBin) {
         const double contentTheta = h1RawTheta->GetBinContent(iBin);
         const double contentZ     = h1RawZ->GetBinContent(iBin);

         if (contentTheta > 0.0)
            sumsTheta[iBin] += contentTheta;
         if (contentZ > 0.0)
            sumsZ[iBin] += contentZ;
      }

      // Running sums of the cross products. Loop bounds now match the
      // sums* loop above exactly (ROOT bins 1..kTotalBins), so every real
      // bin -- including the last one -- is accumulated, and neither the
      // underflow bin (0) nor an out-of-range index is touched.
      for (int iBin = 1; iBin <= kTotalBins; ++iBin) {
         for (int jBin = 1; jBin <= kTotalBins; ++jBin) {
            tempTheta[iBin][jBin] +=
                h1RawTheta->GetBinContent(iBin) * h1RawTheta->GetBinContent(jBin);
            tempZ[iBin][jBin] += h1RawZ->GetBinContent(iBin) * h1RawZ->GetBinContent(jBin);
         }
      }

      h2RawTheta->Reset();
      h2RawZ->Reset();
   }
   // ==================== end loop over events ====================

   // ------------------------------------------------------------ covariance
   //
   // Cov(H_i, H_j) = Sum_k x_i^k x_j^k - (Sum_k x_i^k)(Sum_k x_j^k) / N
   // is already the (unnormalized) sample covariance of the SUMMED yield
   // H_i = Sum_k x_i^k, i.e. of the raw histogram you actually plot as your
   // EEC distribution. The previous version divided this by an additional
   // factor of nEv*nEv, which instead gives the covariance of the per-event
   // MEAN H_i/N -- a mismatch in overall normalization relative to the raw
   // yield histogram (off by a factor of nEv on the variance, sqrt(nEv) on
   // the uncertainty). Only the /nEv term inside the parentheses belongs
   // here; there should be no further division by nEv outside of it.
   //
   // If your nominal EEC histogram is itself normalized by nEv (i.e. you
   // plot H_i/N rather than H_i), divide covZVal/covThetaVal by nEv*nEv
   // instead of leaving them as the raw-yield covariance -- just make sure
   // whichever convention you pick here matches the convention of the
   // central-value histogram you're attaching these uncertainties to.
   for (int iBin = 1; iBin <= kTotalBins; ++iBin) {
      for (int jBin = 1; jBin <= kTotalBins; ++jBin) {
         double covZVal     = tempZ[iBin][jBin]     - sumsZ[iBin]     * sumsZ[jBin]     / nEv;
         double covThetaVal = tempTheta[iBin][jBin] - sumsTheta[iBin] * sumsTheta[jBin] / nEv;

         if (std::abs(covZVal) > 0.0) {
            std::cout << "iBin: " << iBin << ", jBin: " << jBin
                      << ", sums_z[iBin] " << sumsZ[iBin]
                      << ", sums_z[jBin] " << sumsZ[jBin]
                      << ", temp_Z[iBin][jBin]: " << tempZ[iBin][jBin]
                      << ", covZVal " << covZVal << std::endl;
         }

         covZ->SetBinContent(iBin, jBin, covZVal);
         covTheta->SetBinContent(iBin, jBin, covThetaVal);
      }
   }

   TH2D *corr_Z = (TH2D *)covZ->Clone("corr_Z");
   corr_Z->SetTitle("Correlation matrix;bin i;bin j");

   int nBins = covZ->GetNbinsX();

   for (int i = 1; i <= nBins; ++i) {
      double Cii = covZ->GetBinContent(i, i);
      if (Cii <= 0)
         continue;

      for (int j = 1; j <= nBins; ++j) {
         double Cjj = covZ->GetBinContent(j, j);
         if (Cjj <= 0)
            continue;

         double Cij = covZ->GetBinContent(i, j);
         double rho = Cij / std::sqrt(Cii * Cjj);

         corr_Z->SetBinContent(i, j, rho);
      }
   }

   TH2D *corr_theta = (TH2D *)covTheta->Clone("corrTheta");
   corr_theta->SetTitle("Correlation matrix;bin i;bin j");

   int nBinsTheta = covTheta->GetNbinsX();

   for (int i = 1; i <= nBinsTheta; ++i) {
      double Cii = covTheta->GetBinContent(i, i);
      if (Cii <= 0)
         continue;

      for (int j = 1; j <= nBinsTheta; ++j) {
         double Cjj = covTheta->GetBinContent(j, j);
         if (Cjj <= 0)
            continue;

         double Cij = covTheta->GetBinContent(i, j);
         double rho = Cij / std::sqrt(Cii * Cjj);

         corr_theta->SetBinContent(i, j, rho);
      }
   }

   // ------------------------------------------------------- stat. uncertainty
   //
   // sqrt(diag(cov)) is an ABSOLUTE uncertainty on the raw summed bin content
   // sumsTheta[i]/sumsZ[i], i.e. it lives in the same (un-normalized) units
   // as sumsTheta[i]/sumsZ[i] themselves. To compare it against a nominal
   // differential distribution that has been turned into a per-event,
   // per-unit-z (or per-unit-theta) density -- via DivideByBin(...) followed
   // by Scale(1./nEvents), as in the "without considering correlations"
   // cross-check -- sqrt(diag(cov)) needs to go through the SAME two
   // normalization steps, not be divided by the central value. Dividing by
   // the central value would give a genuinely different quantity (a
   // relative/fractional error), which is not what's being compared here:
   // the non-uniform bin widths of the symmetric double-log binning (narrow
   // at the two extremes, wide in the middle) are themselves a major driver
   // of the difference in shape between the two curves.
   TH1D *statErrZ = new TH1D("stat_error_Z", "stat_error_Z;bin index;Stat. Error",
                              kTotalBins, 0, kTotalBins);
   TH1D *statErrTheta = new TH1D("stat_error_Theta", "stat_error_Theta;bin index;Stat. Error",
                                 kTotalBins, 0, kTotalBins);

   for (int iBin = 1; iBin <= kTotalBins; ++iBin) {
      // zBins/thetaBins hold kTotalBins+1 edges (indices 0..kTotalBins), so
      // ROOT bin iBin (1..kTotalBins) spans [Bins[iBin-1], Bins[iBin]] --
      // the same physical interval used by FindBin() to fill this bin.
      const double widthZ     = zBins[iBin]     - zBins[iBin - 1];
      const double widthTheta = thetaBins[iBin] - thetaBins[iBin - 1];

      const double sigmaZ     = std::sqrt(covZ->GetBinContent(iBin, iBin));
      const double sigmaTheta = std::sqrt(covTheta->GetBinContent(iBin, iBin));

      statErrZ->SetBinContent(iBin, sigmaZ / (widthZ * nEv));
      statErrTheta->SetBinContent(iBin, sigmaTheta / (widthTheta * nEv));
   }

   // ---------------------------------------------------------------- output
   TFile outFile(kOutputFileName, "RECREATE");
   covZ->Write();
   covTheta->Write();
   corr_theta->Write();
   corr_Z->Write();
   statErrZ->Write();
   statErrTheta->Write();
   outFile.Close();

   delete h1RawTheta;
   delete h1RawZ;
   delete h2RawTheta;
   delete h2RawZ;
   delete covZ;
   delete covTheta;
   delete corr_Z;
   delete corr_theta;
   delete statErrZ;
   delete statErrTheta;

   inputSmeared->Close();
   delete inputSmeared;
}

int main()
{
   preUnfoldingCovariance();
   return 0;
}