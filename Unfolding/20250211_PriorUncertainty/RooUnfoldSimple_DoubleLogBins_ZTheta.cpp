//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldExample.cxx 279 2011-02-11 18:23:44Z T.J.Adye $
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
// Adapted by Hannah Bossi for the ALEPH EEC analysis.
//==============================================================================


//#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
#include <fstream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "TFile.h"
#include "TVectorD.h"

#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TPostScript.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TMath.h"
#include "TF1.h"

#include "RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
//#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;
#define MAXPAIR 10000

int FindBin(double Value, int NBins, double Bins[]);
void RooSimplenSDPbPb_pp2_DoubleLogBins(std::string date);
//==============================================================================
// Helper Functions
//==============================================================================

int FindBin(double Value, int NBins, double Bins[])
{
   for(int i = 0; i < NBins; i++)
      if(Value < Bins[i])
         return i - 1;
   return NBins;
}


//==============================================================================
// Unfolding for Z and Theta
//==============================================================================

void RooUnfoldSimple_DoubleLogBins_ZTheta(std::string date = "02122025"){

    //------------------------------------
   // define the binning
   //------------------------------------

   // theta binning
   const int BinCount = 100;
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;

   // z binning
   double zBins[2*BinCount+1];
   double zBinMin = (1- cos(0.002))/2;
   double zBinMax = 0.5;

   // energy binning
   const int EnergyBinCount = 15;
   double EnergyBins[EnergyBinCount+1];
   double EnergyBinMin = 4e-6;
   double EnergyBinMax = 0.2;
   double logMin = std::log10(EnergyBinMin);
   double logMax = std::log10(EnergyBinMax);
   double logStep = (logMax - logMin) / (EnergyBinCount);

   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);

   }

  //  for(int e = 0; e <= EnergyBinCount; e++){
  //     double logValue = logMin + e * logStep;
  //     EnergyBins[e] =  std::pow(10, logValue);
  //  }

  // here is what I denote as binning option #1
  std::vector<double> e1e2BinsUnfolded = {0.0, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.00125, 0.0015, 0.00175, 0.002, 0.00225, 0.0025, 0.00275, 0.003, 0.0035, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.3};
  // binning in the original version used for the preliminary result
  // std::vector<double> e1e2BinsUnfolded = {0.0, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.00125, 0.0015, 0.00175, 0.002, 0.0025, 0.003, 0.004, 0.01, 0.04, 0.07, 0.15, 0.3};

  // binning option #3
  //std::vector<double> e1e2BinsUnfolded = {0.0, 0.0001, 0.00020, 0.0005, 0.00075, 0.001, 0.00125, 0.0015, 0.00175, 0.002, 0.00225, 0.0025,  0.00275, 0.003, 0.004, 0.01, 0.04, 0.07, 0.15, 0.3};


  // create the raw histograms
  TH2D *h2raw_Theta = new TH2D("r_Theta","raw_Theta",2 * BinCount, 0, 2 * BinCount , e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data() );
  TH2D *h2raw_Z     = new TH2D("r_Z","raw_Z",2 * BinCount, 0, 2 * BinCount ,e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data() );

  // create the smeared histograms
  TH2D *h2smeared_Theta = new TH2D("smeared_Theta","smeared_Theta", 2 * BinCount, 0, 2 * BinCount , e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data() );
  TH2D *h2smeared_Z     = new TH2D("smeared_Z","smeared_Z", 2 * BinCount, 0, 2 * BinCount ,e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data());

  // create the response histograms
  TH2D* h2resp_Theta = new TH2D("h2resp_Theta", "h2resp_Z", 2 * BinCount, 0, 2 * BinCount, 2 * BinCount, 0, 2 * BinCount);
  TH2D* h2resp_Z     = new TH2D("h2resp_Z", "h2resp_Z", 2 * BinCount, 0, 2 * BinCount, 2 * BinCount, 0, 2 * BinCount);

  // create the true histograms
  TH2D *h2true_Theta = new TH2D("true_Theta","true_Theta",2 * BinCount, 0, 2 * BinCount,e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data() );
  TH2D *h2true_Z = new TH2D("true_Z","true_Z",2 * BinCount, 0, 2 * BinCount,e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data() );

  TH1D* closureCheck_Theta = new TH1D("h1MCGen_Theta", "h1MCGen_Theta", 2*BinCount, 0, 2*BinCount);
  TH1D* closureCheck_Z = new TH1D("h1MCGen_Z", "h1MCGen_Z", 2*BinCount, 0, 2*BinCount);


  // do a sumw2 on the histograms
  h2raw_Theta->Sumw2();
  h2raw_Z->Sumw2();
  h2smeared_Theta->Sumw2();
  h2smeared_Z->Sumw2();
  h2resp_Theta->Sumw2();
  h2resp_Z->Sumw2();
  h2true_Theta->Sumw2();
  h2true_Z->Sumw2();

  TString fnamesmeared = "UnfoldingInputData_12282024.root";
  TFile *inputsmeared =TFile::Open(fnamesmeared);
  TTree *smeared=(TTree*)inputsmeared->Get("UnmatchedPairTree");
  Int_t nEv=smeared->GetEntries();
  Double_t e1e2data[MAXPAIR], thetaData[MAXPAIR];
  double eff1[MAXPAIR], eff2[MAXPAIR], recoE1Data[MAXPAIR], recoE2Data[MAXPAIR];
  int nPairsData;
  smeared->SetBranchAddress("NUnmatchedPair",&nPairsData);
  smeared->SetBranchAddress("E1E2RecoUnmatched", &e1e2data);
  smeared->SetBranchAddress("DistanceUnmatchedReco", &thetaData);
  smeared->SetBranchAddress("RecoE1Unmatched", &recoE1Data);
  smeared->SetBranchAddress("RecoE2Unmatched", &recoE2Data);

  std::cout << "Number of entries in the smeared tree: " << nEv << std::endl;

  for(int iEntry=0; iEntry< nEv; iEntry++){
    smeared->GetEntry(iEntry);
    for(int i=0; i<nPairsData; i++){
      if(recoE1Data[i] < 0 || recoE2Data[i] < 0) continue; // skip over the unmatched pairs
      if(thetaData[i] <  BinMin) continue;
      int BinTheta = FindBin(thetaData[i], 2 * BinCount, Bins);
      double z = (1-cos(thetaData[i]))/2;
      h2raw_Theta->Fill(BinTheta,e1e2data[i], 1.0);
      int BinZ = FindBin(z, 2*BinCount, zBins);
      h2raw_Z->Fill(BinZ,e1e2data[i], 1.0);
    }
  }
  std::cout << "Integral in theta " << h2raw_Theta->Integral() << " Integral in Z " << h2raw_Z->Integral() << std::endl;



  RooUnfoldResponse response_Theta;
  response_Theta.Setup(h2smeared_Theta,h2true_Theta);

  RooUnfoldResponse response_Z;
  response_Z.Setup(h2smeared_Z,h2true_Z);



  double e1e2recoMC[MAXPAIR], e1e2gen[MAXPAIR], thetaRecoMC[MAXPAIR], thetaGen[MAXPAIR];
  double recoE1[MAXPAIR], recoE2[MAXPAIR], genE1[MAXPAIR], genE2[MAXPAIR], genE[MAXPAIR];
  double recoEfficiency1[MAXPAIR], recoEfficiency2[MAXPAIR];
  int nPairsMC;


  TString fnamemc = "LEP1MC1994_recons_aftercut-001_Matched.root";
  TFile *inputmc =TFile::Open(fnamemc);
  TTree *mc=(TTree*)inputmc->Get("PairTree");

  TFile* reweightFile = TFile::Open("ReweightingUncertainty_01092025.root");
  TH2D* reweightFactors_Theta = (TH2D*)reweightFile->Get("reweightFactors_Theta");

  std::vector<TF1*> reweightingFits;

  int nBinsX = reweightFactors_Theta->GetNbinsX();
  for (int i = 1; i <= nBinsX; ++i) { // Bins start from 1 in ROOT
      // Get projection onto Y-axis for bin i in X
      TH1D *hProj = reweightFactors_Theta->ProjectionY(Form("py_%d", i), i, i);
      TF1 *fitFunc = new TF1(Form("fit_%d", i), "pol1", 1e-6, 0.3);
      hProj->Fit(fitFunc, "QR"); // "Q" option suppresses output
      reweightingFits.push_back(fitFunc);
      gStyle->SetOptStat(0);
      gStyle->SetOptTitle(0);
      TCanvas *canvas = new TCanvas(Form("canvas_%d", i), Form("Fit Projection X-bin %d", i), 600, 600);
      canvas->SetLogx();
      hProj->Draw(); // Draw histogram
      fitFunc->SetLineColor(kRed);
      fitFunc->Draw("same"); // Overlay fit
        // Add a legend
      TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.85);
      legend->AddEntry(hProj, "Data Projection", "lep");
      legend->AddEntry(fitFunc, Form("Fit: %.3f + %.3f*x", fitFunc->GetParameter(0), fitFunc->GetParameter(1)), "l");
      legend->Draw();

      // Save the canvas as an image and a ROOT file
      canvas->SaveAs(Form("./fits/fit_projection_bin_%d.png", i));  // Save as PNG
      delete canvas;
      delete hProj; // Clean up histogram after use
  }

  Int_t nEv2=mc->GetEntries();
  std::cout << "nEvents in the mc " << nEv2 << std::endl;
  //------------------------------------------------
  mc->SetBranchAddress("NPair",&nPairsMC);
  mc->SetBranchAddress("E1E2Reco", &e1e2recoMC);
  mc->SetBranchAddress("E1E2Gen", &e1e2gen);
  mc->SetBranchAddress("DistanceReco", &thetaRecoMC);
  mc->SetBranchAddress("DistanceGen", &thetaGen);
  mc->SetBranchAddress("RecoE1", &recoE1);
  mc->SetBranchAddress("RecoE2", &recoE2);

  Int_t countm=0;
  int countLargeReweight = 0;
  for(int iEntry=0; iEntry< nEv2; iEntry++){
    mc->GetEntry(iEntry);
    for(int i=0; i<nPairsMC; i++){
      if(recoE1[i] < 0 || recoE2[i] < 0) continue; // skip over the unmatched pairs
      if(thetaRecoMC[i] < BinMin || thetaGen[i] < BinMin )continue;
      int BinThetaMeasuredMC = FindBin(thetaRecoMC[i], 2 * BinCount, Bins);
      int BinThetaGenMC = FindBin(thetaGen[i], 2 * BinCount, Bins);
      double zMeasuredMC = (1-cos(thetaRecoMC[i]))/2;
      double zGenMC = (1-cos(thetaGen[i]))/2;
      int BinZMeasured = FindBin(zMeasuredMC, 2 * BinCount, zBins);
      int BinZGen = FindBin(zGenMC, 2*BinCount,zBins);
      int BinEnergyGen = FindBin(e1e2gen[i], EnergyBinCount, EnergyBins);
      int BinEnergyRecoMC = FindBin(e1e2recoMC[i], EnergyBinCount, EnergyBins);
      if(BinZMeasured < 0 || BinZGen < 0 || BinEnergyGen < 0 || BinEnergyRecoMC < 0 || BinThetaGenMC < 0 || BinThetaMeasuredMC < 0){
        std::cout << "Theta " << thetaRecoMC[i] << " ZData " << zMeasuredMC << " EnergyData " << e1e2recoMC[i] << std::endl;
      }
      if(BinZMeasured > 200 || BinZGen > 200 || BinEnergyGen > 200 || BinEnergyRecoMC > 200 || BinThetaGenMC > 200 || BinThetaMeasuredMC > 200){
        std::cout << "Theta " << thetaRecoMC[i] << " ZData " << zMeasuredMC << " EnergyData " << e1e2recoMC[i] << std::endl;
      }

      closureCheck_Theta->Fill(BinThetaGenMC,e1e2gen[i]);
      h2true_Theta->Fill(BinThetaGenMC, e1e2gen[i]);
      h2true_Z->Fill(BinZGen, e1e2gen[i]);
      closureCheck_Z->Fill(BinZGen, e1e2gen[i]);
      h2smeared_Theta->Fill(BinThetaMeasuredMC,e1e2recoMC[i]);
      h2smeared_Z->Fill(BinZMeasured, e1e2recoMC[i]);
      // reweight the response
      // std vectors are 0 indexed while binning is index 1
     //std::cout << "Checking reweighting factor for bin " << BinThetaMeasuredMC-1 << " and  energy " << e1e2recoMC[i] << std::endl;

      double reweightFactor;
      if(BinThetaMeasuredMC > 0 && BinThetaMeasuredMC < 199) reweightFactor = reweightingFits.at(BinThetaMeasuredMC-1)->Eval(e1e2recoMC[i]);//reweightFactors_Theta->GetBinContent(BinThetaMeasuredMC, BinEnergyRecoMC);
      else reweightFactor = 1;
      if(reweightFactor < 0)reweightFactor = abs(reweightFactor);
      if(abs(reweightFactor-1) > 0.5){
        // std::cout << "Large Reweighting factor of " << reweightFactor << " for bin " << BinThetaMeasuredMC-1 << " and  energy " << e1e2recoMC[i] <<  std::endl;
        countLargeReweight++;
      }
      response_Theta.Fill(BinThetaMeasuredMC, e1e2recoMC[i],BinThetaGenMC,e1e2gen[i]*reweightFactor);
      response_Z.Fill(BinZMeasured, e1e2recoMC[i], BinZGen, e1e2gen[i]*reweightFactor);
    }
  }

  std::cout << "Entries in h2true_Z " << h2true_Z->GetEntries() << " entries in h2true_Theta " << h2true_Theta->GetEntries() << std::endl;
  std::cout << "Large Reweight Count: " << countLargeReweight << std::endl;
  //------------------------------------------------




  TFile *fout = new TFile(Form("unfoldingE2C_DataUnfolding_DoubleLogBinning_BinningOption1_Reweighted_FittingMethod%s.root",date.c_str()),"RECREATE");
  fout->cd();
  h2raw_Theta->Write();
  h2raw_Z->Write();
  h2smeared_Theta->Write();
  h2smeared_Z->Write();
  h2true_Theta->Write();
  h2true_Z->Write();
  closureCheck_Theta->Write();
  closureCheck_Z->Write();

   TH1D *HTrue_Theta = (TH1D*)h2true_Theta->ProjectionX("h1True_Theta_ProjectionX");
   for (int i = 1; i <= h2true_Theta->GetNbinsX(); ++i) {
        double weight = 0;
        double error = 0;
        for (int j = 1; j <= h2true_Theta->GetNbinsY(); ++j) {
            double binContent = h2true_Theta->GetBinContent(i, j);
            double binError= h2true_Theta->GetBinError(i,j);
            double binCenter = h2true_Theta->GetYaxis()->GetBinCenter(j);
            weight += binContent*((binCenter));
            error += pow(binError*binCenter, 2);;
        }
        HTrue_Theta->SetBinContent(i, weight);
        HTrue_Theta->SetBinError(i, sqrt(error));
    }

    HTrue_Theta->Write();

       TH1D *HTrue_Z = (TH1D*)h2true_Z->ProjectionX("h1True_Z_ProjectionX");
   for (int i = 1; i <= h2true_Z->GetNbinsX(); ++i) {
        double weight = 0;
        double error = 0;
        for (int j = 1; j <= h2true_Z->GetNbinsY(); ++j) {
            double binContent = h2true_Z->GetBinContent(i, j);
            double binError= h2true_Z->GetBinError(i,j);
            double binCenter = h2true_Z->GetYaxis()->GetBinCenter(j);
            weight += binContent*((binCenter));
            error += pow(binError*binCenter, 2);;
        }
        HTrue_Z->SetBinContent(i, weight);
        HTrue_Z->SetBinError(i, sqrt(error));
    }
    HTrue_Z->Write();

  // int iter = 2;
  for(int iter = 4; iter < 5; iter++){

    std::cout << "Unfolding for theta iter " << iter << std::endl;
    // -------------------------------
    // unfolding for the theta
    RooUnfoldBayes  unfold_Theta(&response_Theta, h2raw_Theta, iter);
    // unfold_Theta.SetNToys(0);
    TH2D* hunf_Theta =  dynamic_cast<TH2D*>(unfold_Theta.Hreco());
    TH1* hfold_Theta = response_Theta.ApplyToTruth(hunf_Theta, "");
    TH2D *htempUnf_Theta=(TH2D*)hunf_Theta->Clone("htempUnf_Theta");
    htempUnf_Theta->SetName(Form("Bayesian_Unfoldediter%d_Theta",iter));

    TH2D *htempFold_Theta=(TH2D*)hfold_Theta->Clone("htempFold_Theta");
    htempFold_Theta->SetName(Form("Bayesian_Foldediter%d_Theta",iter));
    // -------------------------------

    std::cout << "Unfolding for z iter " << iter << std::endl;


    // -------------------------------
    // unfolding for the Z
    RooUnfoldBayes  unfold_Z(&response_Z, h2raw_Z, iter);
    TH2D* hunf_Z =  dynamic_cast<TH2D*>(unfold_Z.Hreco());

    TH1* hfold_Z = response_Z.ApplyToTruth(hunf_Z, "");
    TH2D *htempUnf_Z=(TH2D*)hunf_Theta->Clone("htempUnf_Z");
    htempUnf_Z->SetName(Form("Bayesian_Unfoldediter%d_Z",iter));

    TH2D *htempFold_Z=(TH2D*)hfold_Z->Clone("htempFold_Z");
    htempFold_Z->SetName(Form("Bayesian_Foldediter%d_Z",iter));

    // -------------------------------

    htempUnf_Theta->Write();
    htempFold_Theta->Write();
    htempUnf_Z->Write();
    htempFold_Z->Write();

    TH1D *H_Theta = (TH1D*)htempUnf_Theta->ProjectionX("");
    TH1D *H_Z = (TH1D*)htempUnf_Z->ProjectionX();
    H_Theta->Reset();
    H_Z->Reset();
    // handle the errors
    for (int i = 1; i <= htempUnf_Z->GetNbinsX(); ++i) {
        double weight = 0;
        double error = 0;
        for (int j = 1; j <= htempUnf_Z->GetNbinsY(); ++j) {
            double binContent = htempUnf_Z->GetBinContent(i, j);
            double binError= htempUnf_Z->GetBinError(i,j);
            double binCenter = htempUnf_Z->GetYaxis()->GetBinCenter(j);
            weight += binContent*((binCenter));
            error += pow(binError*binCenter, 2);;
        }
        H_Z->SetBinContent(i, weight);
        H_Z->SetBinError(i, sqrt(error));
    }
    // handle the errors
    for (int i = 1; i <= htempUnf_Theta->GetNbinsX(); ++i) {
        double weight = 0;
        double error = 0;
        for (int j = 1; j <= htempUnf_Theta->GetNbinsY(); ++j) {
            double binContent = htempUnf_Theta->GetBinContent(i, j);
            double binError= htempUnf_Theta->GetBinError(i,j);
            double binCenter = htempUnf_Theta->GetYaxis()->GetBinCenter(j);
            weight += binContent*((binCenter));
            error += pow(binError*binCenter, 2);;
        }
        H_Theta->SetBinContent(i, weight);
        H_Theta->SetBinError(i, sqrt(error));
    }
    H_Theta->SetName(Form("Bayesian_Unfoldediter%d_Theta_ProjectionX",iter));
    H_Z->SetName(Form("Bayesian_Unfoldediter%d_Z_ProjectionX",iter));

    H_Theta->Write();
    H_Z->Write();

  }
  fout->Close();
}




//#ifndef __CINT__
int main () {  RooUnfoldSimple_DoubleLogBins_ZTheta(); return 0; }  // Main program when run stand-alone
//#endif
