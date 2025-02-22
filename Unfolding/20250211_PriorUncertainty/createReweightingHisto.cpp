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
#include "TPad.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TPostScript.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TMath.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TF1.h"

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;
#define MAXPAIR 10000

int FindBin(double Value, int NBins, double Bins[]);
void createReweightingHisto(std::string date);
void SetPad(TPad &P);
void MakeCanvas2D(TH2D* Histogram, std::string Output, std::string X, std::string Y, double WorldMin, double WorldMax,  bool LogX);
void DivideByBin(TH1D &H, double Bins[]);
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

void createReweightingHisto(std::string date = "01092025"){

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

  // here is what I denote as binning option #1
  std::vector<double> e1e2BinsUnfolded = {0.000001, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.00125, 0.0015, 0.00175, 0.002, 0.00225, 0.0025, 0.00275, 0.003, 0.0035, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.3};

  // create the raw histograms
  TH2D *h2raw_Theta = new TH2D("r_Theta","raw_Theta",2 * BinCount, 0, 2 * BinCount , e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data() );
  TH2D *h2raw_Z     = new TH2D("r_Z","raw_Z",2 * BinCount, 0, 2 * BinCount ,e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data() );

  // create the smeared histograms
  TH2D *h2smeared_Theta = new TH2D("smeared_Theta","smeared_Theta", 2 * BinCount, 0, 2 * BinCount , e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data() );
  TH2D *h2smeared_Z     = new TH2D("smeared_Z","smeared_Z", 2 * BinCount, 0, 2 * BinCount ,e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data());

  // do a sumw2 on the histograms
  h2raw_Theta->Sumw2();
  h2raw_Z->Sumw2();
  h2smeared_Theta->Sumw2();
  h2smeared_Z->Sumw2();

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

  double e1e2recoMC[MAXPAIR], e1e2gen[MAXPAIR], thetaRecoMC[MAXPAIR], thetaGen[MAXPAIR];
  double recoE1[MAXPAIR], recoE2[MAXPAIR], genE1[MAXPAIR], genE2[MAXPAIR], genE[MAXPAIR];
  double recoEfficiency1[MAXPAIR], recoEfficiency2[MAXPAIR];
  int nPairsMC;


  TString fnamemc = "LEP1MC1994_recons_aftercut-001_Matched.root";
  TFile *inputmc =TFile::Open(fnamemc);
  TTree *mc=(TTree*)inputmc->Get("PairTree");

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
      h2smeared_Theta->Fill(BinThetaMeasuredMC,e1e2recoMC[i]);
      h2smeared_Z->Fill(BinZMeasured, e1e2recoMC[i]);
    }
  }

  //------------------------------------------------
  h2raw_Theta->Scale(1./nEv);
  h2smeared_Theta->Scale(1./nEv2);
  h2raw_Theta->Divide(h2smeared_Theta);
  h2raw_Theta->SetName("reweightFactors_Theta");

  TFile *fout = new TFile(Form("ReweightingUncertainty_%s.root",date.c_str()),"RECREATE");
  fout->cd();
  h2raw_Theta->Write();

  fout->Close();



  //------------------------------------------------
  // now do the plotting
  //------------------------------------------------
  TLatex* cms = new TLatex(0.12,0.92,"#bf{ALEPH} Work in Progress e^{+}e^{-}");
  cms->SetNDC();
  cms->SetTextSize(0.05);
  cms->SetTextFont(42);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  // --------------------------------------------
  // test draw a projection onto a canvas
  // --------------------------------------------
  TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
  c2->SetTickx(1);
  c2->SetTicky(1);
  //c2->SetLogy();
  c2->SetLogx();
  // c->SetLogz();
  c2->SetRightMargin(0.13);
  c2->SetLeftMargin(0.13);

  TH1D* projectionTest = (TH1D*)h2raw_Theta->ProjectionY("projY", 20, 20);
  projectionTest->Draw();

  TF1 *fpol = new TF1("fpol", "pol1", 1e-5, 1e-2);
  fpol->SetLineWidth(2);
  projectionTest->Fit(fpol, "R");

  c2->SaveAs("BinningTest.pdf");


  // --------------------------------------------
  // draw divisional reweighting factor
  // --------------------------------------------

  MakeCanvas2D(h2raw_Theta,  "ReweightingFactors_EEC2","#theta_{L}", "E_{i}E_{j}/E^{2}", 2e-4, 100, true);
}

void MakeCanvas2D(TH2D* Histogram, std::string Output,
   std::string X, std::string Y, double WorldMin, double WorldMax, bool LogX)
{
   int N = Histogram->GetNbinsX();

   double MarginL = 180;
   double MarginR = 60;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 0 : 0;
   double WorldXMax = LogX ? N : M_PI;

   double PadWidth = 1200;
   double PadHeight =  640 + 240;
   double PadRHeight = 0.001;

   double CanvasWidth = MarginL + PadWidth + MarginR;
   double CanvasHeight = MarginT + PadHeight + PadRHeight + MarginB;

   MarginL = MarginL / CanvasWidth;
   MarginR = MarginR / CanvasWidth;
   MarginT = MarginT / CanvasHeight;
   MarginB = MarginB / CanvasHeight;

   PadWidth   = PadWidth / CanvasWidth;
   PadHeight  = PadHeight / CanvasHeight;
   PadRHeight = PadRHeight / CanvasHeight;

   TCanvas Canvas("Canvas", "", CanvasWidth, CanvasHeight);

   TPad Pad("Pad", "", MarginL, MarginB + PadRHeight, MarginL + PadWidth, MarginB + PadHeight + PadRHeight);
   Pad.SetLogy();
   SetPad(Pad);

   Pad.cd();


   double WorldRMin = Histogram->GetYaxis()->GetXmin();
   double WorldRMax = Histogram->GetYaxis()->GetXmax();

   TH2D HWorld("HWorld", "", N, WorldXMin, WorldXMax, 100, WorldRMin, WorldRMax);
   HWorld.SetStats(0);
   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);
   HWorld.SetMaximum(Histogram->GetMaximum());
   HWorld.SetMinimum(Histogram->GetMinimum());

   HWorld.Draw("axis");
   Histogram->Draw("colz same");

   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : M_PI / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : M_PI / 2, 1000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kBlack);
   G.SetLineWidth(1);
   G.Draw("l");


   TH2D HWorldR("HWorldR", "", N, WorldXMin, WorldXMax, 100, WorldRMin, WorldRMax);
   TGraph G2;
    HWorldR.SetMaximum(Histogram->GetMaximum()/100);
    HWorldR.SetMinimum(Histogram->GetMinimum());


   double BinMin    = 0.002;
   double BinMiddle = M_PI / 2;
   double BinMax    = M_PI - 0.002;

   Canvas.cd();
   // this is the axis from 0 to pi/2
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth*0.9/ 2, MarginB, BinMin, BinMiddle, 510, "GS");
   // this is the axis from pi/2 to pi
   TGaxis X2(MarginL + 0.9*PadWidth, MarginB, MarginL +(0.9*PadWidth) / 2, MarginB, BinMin, BinMiddle, 510, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldRMin, WorldRMax, 510, "G");
   std::cout << "WorldR min" << WorldRMin << " WorldRMax " << WorldRMax  << std::endl;
   // used if log X == false
   TGaxis XL1(MarginL, MarginB, MarginL + PadWidth, MarginB, 0, M_PI, 510, "S");
   TGaxis XL2(MarginL, MarginB + PadRHeight, MarginL + PadWidth, MarginB + PadRHeight, 0, M_PI, 510, "+-S");

   Y2.SetLabelFont(42);
   XL1.SetLabelFont(42);
   XL2.SetLabelFont(42);

   X1.SetLabelSize(0);
   X2.SetLabelSize(0);
   X3.SetLabelSize(0);
   X4.SetLabelSize(0);
   // XL1.SetLabelSize(0);
   XL2.SetLabelSize(0);
  //Y2.SetLabelSize(0);

   X1.SetTickSize(0.06);
   X2.SetTickSize(0.06);
   X3.SetTickSize(0.06);
   X4.SetTickSize(0.06);
   XL1.SetTickSize(0.03);
   XL2.SetTickSize(0.03);

   if(LogX == true)
   {
      X1.Draw();
      X2.Draw();
   }
   if(LogX == false)
   {
      XL1.Draw();
   }

   Y2.Draw();

   TLatex Latex;
   Latex.SetNDC();
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.035);
   Latex.SetTextAlign(23);
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.115, MarginB - 0.01, "0.01");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.270, MarginB - 0.01, "0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.425, MarginB - 0.01, "1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.505, MarginB - 0.01, "#pi - 1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.650, MarginB - 0.01, "#pi - 0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.8, MarginB - 0.01, "#pi - 0.01");

   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.47, 1 - MarginT - 0.015, "#theta_{L} = #pi/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.5*0.9, MarginB * 0.3, X.c_str());

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(90);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV, Work-in-progress");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "Work-in-progress, 2025 Feb 12th, Hannah Bossi");


   Canvas.SaveAs((Output + ".pdf").c_str());
}

void SetPad(TPad &P){
   P.SetLeftMargin(0);
   P.SetTopMargin(0);
   P.SetRightMargin(0.1);
   P.SetBottomMargin(0);
   P.SetTickx();
   P.SetTicky();
   P.Draw();
}




void DivideByBin(TH1D &H, double Bins[])
{
   int N = H.GetNbinsX();
   for(int i = 1; i <= N; i++)
   {
      double L = Bins[i-1];
      double R = Bins[i];
      H.SetBinContent(i, H.GetBinContent(i) / (R - L));
      H.SetBinError(i, H.GetBinError(i) / (R - L));
   }
}




//#ifndef __CINT__
int main () {  createReweightingHisto(); return 0; }  // Main program when run stand-alone
//#endif
