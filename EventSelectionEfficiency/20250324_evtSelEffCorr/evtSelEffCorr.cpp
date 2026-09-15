// ./EvtSelEffCorr.exe --Input v0/LEP1MC1994_recons_aftercut-001_Matched.root --Matched PairTree --Unmatched UnmatchedPairTree
#include <iostream>
#include <vector>
#include <map>
#include <filesystem>
using namespace std;

// root includes
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"

#include "Messenger.h"
#include "CommandLine.h"
#include "ProgressBar.h"
#include "TauHelperFunctions3.h"
#include "SetStyle.h"
#include "EffCorrFactor.h"

int main(int argc, char *argv[]);
int FindBin(double Value, int NBins, double Bins[]);
void MakeCanvasZ(vector<TH1D>& Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX);
void MakeCanvas(vector<TH1D>& Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX);
void SetPad(TPad &P);
void DivideByBin(TH1D &H, double Bins[]);

// can we put this in the CommonCode/ ?
void FillChain(TChain &chain, const vector<string> &files) {
  for (auto file : files) {
    chain.Add(file.c_str());
  }
}

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   SetThesisStyle();
   static vector<int> Colors = GetCVDColors6();

   string InputFileName          = CL.Get("Input");
   string EvtSelEffFileName      = CL.Get("EvtSelEffName", "EvtSelEff.root");
   string EvtSelEffArgName       = CL.Get("EvtSelEffArgName", "z");
   bool MakeEvtSelEffCorrFactor  = CL.GetBool("MakeEvtSelEffCorrFactor", false);
   string GenTreeName            = CL.Get("Gen", "tgen");
   string GenBeforeTreeName      = CL.Get("GenBefore", "tgenBefore");
   TFile InputFile(InputFileName.c_str());

   double TotalE = 91.1876;
   ParticleTreeMessenger MReco(InputFile, "t");
   ParticleTreeMessenger MGen(InputFile, GenTreeName);
   ParticleTreeMessenger MGenBefore(InputFile, GenBeforeTreeName);

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

   // [Warning] A future todo improvement task after HP2024
   //           to make the bin boundary configurable in the future
   // energy binning
   // const int EnergyBinCount = 10;
   // double EnergyBins[EnergyBinCount+1];
   // double EnergyBinMin = 4e-6;
   // double EnergyBinMax = 0.2;
   // double logMin = std::log10(EnergyBinMin);
   // double logMax = std::log10(EnergyBinMax);
   // double logStep = (logMax - logMin) / (EnergyBinCount);

   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);

   }

   // [Warning] A future todo improvement task after HP2024
   //           to make the bin boundary configurable in the future
   // EnergyBins[0] = 0;
   // for(int e = 1; e <= EnergyBinCount; e++){
   //    double logValue = logMin + e * logStep;
   //    EnergyBins[e] =  std::pow(10, logValue);
   //    printf("%f ", EnergyBins[e]);
   // }
   vector<double> EnergyBins =  {0.0, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.00125, 0.0015, 0.00175, 0.002, 0.00225, 0.0025, 0.00275, 0.003, 0.0035, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.3};


   // -------------------------------------------
   // allocate the histograms
   // -------------------------------------------

   // 2D histograms
   TH2D h2_EvtSelBefore_Theta("h2_EvtSelBefore_Theta", "h2_EvtSelBefore_Theta", 2 * BinCount, 0, 2 * BinCount, EnergyBins.size()-1, EnergyBins.data());
   TH2D h2_EvtSelBefore_Z("h2_EvtSelBefore_Z", "h2_EvtSelBefore_Z", 2 * BinCount, 0, 2 * BinCount, EnergyBins.size()-1, EnergyBins.data());
   TH2D h2_EvtSel_Theta("h2_EvtSel_Theta", "h2_EvtSel_Theta", 2 * BinCount, 0, 2 * BinCount, EnergyBins.size()-1, EnergyBins.data());
   TH2D h2_EvtSel_Z("h2_EvtSel_Z", "h2_EvtSel_Z", 2 * BinCount, 0, 2 * BinCount, EnergyBins.size()-1, EnergyBins.data());

   // 1D histograms
   TH1D h1_EvtSel_Z("h1_EvtSel_Z", "h1_EvtSel_Z", 2 * BinCount, 0, 2 * BinCount);
   TH1D h1_EvtSelBefore_Z("h1_EvtSelBefore_Z", "h1_EvtSelBefore_Z", 2 * BinCount, 0, 2 * BinCount);
   TH1D h1_EvtSel_Theta("h1_EvtSel_Theta", "h1_EvtSel_Theta", 2 * BinCount, 0, 2 * BinCount);
   TH1D h1_EvtSelBefore_Theta("h1_EvtSelBefore_Theta", "h1_EvtSelBefore_Theta", 2 * BinCount, 0, 2 * BinCount);

   // corrected 1D histograms
   TH1D h1_EvtSelCorrected_Z("h1_EvtSelCorrected_Z", "h1_EvtSelCorrected_Z", 2 * BinCount, 0, 2 * BinCount);
   TH1D h1_EvtSelCorrected_Theta("h1_EvtSelCorrected_Theta", "h1_EvtSelCorrected_Theta", 2 * BinCount, 0, 2 * BinCount);

   //------------------------------------
   // define the event selection efficiency correction factor
   //------------------------------------
   EffCorrFactor EvtSelEffCorrFactor;
   if (!filesystem::exists(EvtSelEffFileName.c_str())) MakeEvtSelEffCorrFactor = true;

   if (MakeEvtSelEffCorrFactor)
   {
      printf("[INFO] produce event selection efficiency correction factor (%s), EvtSelEffArgName=%s\n", EvtSelEffFileName.c_str(), EvtSelEffArgName.c_str());
   } else {
      printf("[INFO] applying event selection efficiency correction factor (%s), EvtSelEffArgName=%s\n", EvtSelEffFileName.c_str(), EvtSelEffArgName.c_str());
      EvtSelEffCorrFactor.init(EvtSelEffFileName.c_str(), EvtSelEffArgName.c_str());
   }

   // -------------------------------------
   // loop over the tree after event selections
   // -------------------------------------
   int EntryCount = MGen.GetEntries();
   int nAcceptedEventMGen = 0; 
   for(int iE = 0; iE < EntryCount; iE++)
   {
      MGen.GetEntry(iE);
      MReco.GetEntry(iE); 
      
      // add reco event selection
      if(MReco.passesSTheta < 0.5) continue; // require that the distribution passSTheta
      if(MReco.passesNTrkMin < 0.5) continue; // require at least 5 tracks
      if(MReco.passesTotalChgEnergyMin < 0.5) continue; // require that the total energy is at least 15 GeV

      nAcceptedEventMGen++; 
      // fill the four vector
      vector<FourVector> PGen;
      vector<int> pwflagVec; 
      for(int i = 0; i < MGen.nParticle; i++){
         // charged particle selection
         if(MGen.charge[i] == 0) continue;
         if(MGen.highPurity[i] == false) continue;
         if(abs(cos(MGen.theta[i])) > 0.94) continue; 
         
         PGen.push_back(MGen.P[i]);
         pwflagVec.push_back(MGen.pwflag[i]);
      } // end loop over the particles


      // now calculate and fill the EECs
      for(int i = 0; i < PGen.size(); i++)
      {
        for(int j = i+1; j < PGen.size();j++)
        {
            FourVector Gen1 = PGen.at(i);
            FourVector Gen2 = PGen.at(j);

            // get the proper bins
            int BinThetaGen  = FindBin(GetAngle(Gen1,Gen2), 2 * BinCount, Bins);
            // int BinEnergyGen = FindBin(Gen1[0]*Gen2[0]/(TotalE*TotalE), EnergyBinCount, EnergyBins);
            double zGen = (1-cos(GetAngle(Gen1,Gen2)))/2;
            int BinZGen = FindBin(zGen, 2*BinCount, zBins);

            // calculate the EEC
            double EEC =  Gen1[0]*Gen2[0]/(TotalE*TotalE);

            // fill the histograms
            h2_EvtSel_Theta.Fill(BinThetaGen, EEC, EEC);
            h2_EvtSel_Z.Fill(BinZGen, EEC, EEC);
            h1_EvtSel_Z.Fill(BinZGen, EEC);
            h1_EvtSel_Theta.Fill(BinThetaGen, EEC);
            if (!MakeEvtSelEffCorrFactor)
            {
               double efficiency = EvtSelEffCorrFactor.efficiency((EvtSelEffArgName=="z")? BinZGen: BinThetaGen,
                                                                    EEC);
               // printf("argBin: %d, z: %.3f, eff: %.3f \n", BinZGen, zGen, efficiency);
               // printf("argBin: %d, normEEBin:%d, z: %.3f, normEE: %.3f, eff: %.3f \n", BinZGen, BinEnergyGen, zGen, EEC, efficiency);
               h1_EvtSelCorrected_Z.Fill(BinZGen, EEC/efficiency);
               h1_EvtSelCorrected_Theta.Fill(BinThetaGen, EEC/efficiency);
            }
         }
      }
   } // end loop over the number of events'

   // -------------------------------------
   // loop over the tree before event selection
   // -------------------------------------
   int EntryCountBefore = MGenBefore.GetEntries();
   for(int iE = 0; iE < EntryCountBefore; iE++)
   {
      MGenBefore.GetEntry(iE);
      // fill the four vector
      vector<FourVector> PGenBefore;
      vector<int> pwflagVec; 
      for(int i = 0; i < MGenBefore.nParticle; i++){
        // charged particle selection
       if(MGenBefore.charge[i] == 0) continue;
       // remove conversion electrons
       bool areBothElectrons = (i > 0 && MGenBefore.pwflag[i] == 2 && MGenBefore.pwflag[i-1]  == 2); 
       bool areOppositeCharge = (i > 0 && MGenBefore.charge[i] != MGenBefore.charge[i-1]); 
       float conversionDPhi =  0.05;
       float conversionDTheta = 0.05;
       double deltaPhi = (TMath::Abs(MGenBefore.theta[i] - MGenBefore.theta[i-1])); 
       double deltaTheta = (TMath::ACos(TMath::Cos(MGenBefore.phi[i] - MGenBefore.phi[i-1])) ); 
       bool meetsDeltaThetaReq = 0.0001 < conversionDTheta-deltaTheta; 
       bool meetsDeltaPhiReq = 0.00001 < conversionDPhi- deltaPhi; 
       bool isConversionElectron = areBothElectrons && areOppositeCharge && meetsDeltaThetaReq && meetsDeltaPhiReq; 
       if(isConversionElectron){
          if(PGenBefore.size() > 0 ){
             FourVector poppedParticle = PGenBefore.back(); 
             if(GetAngle(poppedParticle, MGenBefore.P[i]) < 0.1 && pwflagVec.back() == 2 )PGenBefore.pop_back(); 
          }
          continue; 
       }
       
       PGenBefore.push_back(MGenBefore.P[i]);
       pwflagVec.push_back(MGenBefore.pwflag[i]);
       
      } // end loop over the particles

      // now calculate and fill the EECs
      for(int i = 0; i < PGenBefore.size(); i++){
        for(int j = i+1; j < PGenBefore.size();j++){
            FourVector Gen1 = PGenBefore.at(i);
            FourVector Gen2 = PGenBefore.at(j);

            // get the proper bins
            int BinThetaGen  = FindBin(GetAngle(Gen1,Gen2), 2 * BinCount, Bins);
            // int BinEnergyGen = FindBin(Gen1[0]*Gen2[0]/(TotalE*TotalE), EnergyBinCount, EnergyBins);
            double zGen = (1-cos(GetAngle(Gen1,Gen2)))/2;
            int BinZGen = FindBin(zGen, 2*BinCount, zBins);

            // calculate the EEC
            double EEC =  Gen1[0]*Gen2[0]/(TotalE*TotalE);

            // fill the histograms
            h2_EvtSelBefore_Theta.Fill(BinThetaGen, EEC, EEC);
            h2_EvtSelBefore_Z.Fill(BinZGen, EEC, EEC);
            h1_EvtSelBefore_Z.Fill(BinZGen, EEC);
            h1_EvtSelBefore_Theta.Fill(BinThetaGen, EEC);
         }
      }
   } // end loop over the number of events
   // EEC is per-event so scale by the event number
   printf( "h1_EvtSel_Z.GetEntries(): %.3f, EntryCount: %d, h1_EvtSelBefore_Z.GetEntries(): %.3f, EntryCountBefore: %d\n",
            h1_EvtSel_Z.GetEntries(), EntryCount,
            h1_EvtSelBefore_Z.GetEntries(), EntryCountBefore );
   h1_EvtSel_Z.Scale(1.0/nAcceptedEventMGen);
   h1_EvtSelBefore_Z.Scale(1.0/EntryCountBefore);
   h1_EvtSel_Theta.Scale(1.0/nAcceptedEventMGen);
   h1_EvtSelBefore_Theta.Scale(1.0/EntryCountBefore);
   h2_EvtSel_Z.Scale(1.0/nAcceptedEventMGen);
   h2_EvtSelBefore_Z.Scale(1.0/EntryCountBefore);
   h2_EvtSel_Theta.Scale(1.0/nAcceptedEventMGen);
   h2_EvtSelBefore_Theta.Scale(1.0/EntryCountBefore);
   if (!MakeEvtSelEffCorrFactor)
   {
      h1_EvtSelCorrected_Z.Scale(1.0/nAcceptedEventMGen);
      h1_EvtSelCorrected_Theta.Scale(1.0/nAcceptedEventMGen);
   }

   // divide by the bin width
   DivideByBin(h1_EvtSel_Z, zBins);
   DivideByBin(h1_EvtSelBefore_Z, zBins);
   DivideByBin(h1_EvtSel_Theta, Bins);
   DivideByBin(h1_EvtSelBefore_Theta, Bins);
   if (!MakeEvtSelEffCorrFactor)
   {
      DivideByBin(h1_EvtSelCorrected_Z,zBins);
      DivideByBin(h1_EvtSelCorrected_Theta,Bins);
   }

   // set the style for the plots
   h1_EvtSel_Z.SetMarkerColor(Colors[2]);
   h1_EvtSelBefore_Z.SetMarkerColor(Colors[3]);
   h1_EvtSel_Theta.SetMarkerColor(Colors[2]);
   h1_EvtSelBefore_Theta.SetMarkerColor(Colors[3]);
   h1_EvtSel_Z.SetLineColor(Colors[2]);
   h1_EvtSelBefore_Z.SetLineColor(Colors[3]);
   h1_EvtSel_Theta.SetLineColor(Colors[2]);
   h1_EvtSelBefore_Theta.SetLineColor(Colors[3]);
   h1_EvtSel_Z.SetMarkerStyle(20);
   h1_EvtSelBefore_Z.SetMarkerStyle(20);
   h1_EvtSel_Theta.SetMarkerStyle(20);
   h1_EvtSelBefore_Theta.SetMarkerStyle(20);
   h1_EvtSel_Z.SetLineWidth(2);
   h1_EvtSelBefore_Z.SetLineWidth(2);
   h1_EvtSel_Theta.SetLineWidth(2);
   h1_EvtSelBefore_Theta.SetLineWidth(2);
   if (!MakeEvtSelEffCorrFactor)
   {
      h1_EvtSelCorrected_Z.SetMarkerColor(Colors[4]);
      h1_EvtSelCorrected_Z.SetLineColor(Colors[4]);
      h1_EvtSelCorrected_Z.SetMarkerStyle(20);
      h1_EvtSelCorrected_Z.SetLineWidth(2);
      h1_EvtSelCorrected_Theta.SetMarkerColor(Colors[4]);
      h1_EvtSelCorrected_Theta.SetLineColor(Colors[4]);
      h1_EvtSelCorrected_Theta.SetMarkerStyle(20);
      h1_EvtSelCorrected_Theta.SetLineWidth(2);

   }

   if (!MakeEvtSelEffCorrFactor)
   {
      std::vector<TH1D> hists = {h1_EvtSelBefore_Z, h1_EvtSel_Z, h1_EvtSelCorrected_Z};
      std::vector<TH1D> hists_theta = {h1_EvtSelBefore_Theta, h1_EvtSel_Theta, h1_EvtSelCorrected_Theta};

      system("mkdir -p plot/");
      string plotPrefix(EvtSelEffFileName);
      plotPrefix.replace(plotPrefix.find(".root"), 5, "");
      MakeCanvasZ(hists, {"Before Evt. Sel.", "After Evt. Sel.", "Evt. Sel. Corrected"},
                  Form("plot/%s_Closure", plotPrefix.c_str()), "#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",2e-3,1e3, true, true);
      MakeCanvas(hists_theta, {"Before Evt. Sel.", "After Evt. Sel.", "Evt. Sel. Corrected"},
                  Form("plot/%s_Theta_Closure", plotPrefix.c_str()),"#theta", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{#theta}}",2e-3, 5e0, true, true);
   }
   else
   {
      std::vector<TH1D> hists = {h1_EvtSelBefore_Z, h1_EvtSel_Z};
      std::vector<TH1D> hists_theta = {h1_EvtSelBefore_Theta, h1_EvtSel_Theta};

      system("mkdir -p plot/");
      string plotPrefix(EvtSelEffFileName);
      plotPrefix.replace(plotPrefix.find(".root"), 5, "");
      MakeCanvasZ(hists, {"Before Evt. Sel.", "After Evt. Sel."},
                  Form("plot/%s_EvtSel", plotPrefix.c_str()), "#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",2e-3,1e3, true, true);
      MakeCanvas(hists_theta, {"Before Evt. Sel.", "After Evt. Sel."},
                  Form("plot/%s_EvtSel_Theta", plotPrefix.c_str()),"#theta", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{#theta}}",2e-3, 2e0, true, true);
   }

   // now handle the 2D plots

   // write to the output file
   if (MakeEvtSelEffCorrFactor)
   {
      TFile OutputFile(EvtSelEffFileName.c_str(), "RECREATE");
      EvtSelEffCorrFactor.write(OutputFile, "theta", &h1_EvtSel_Theta, &h1_EvtSelBefore_Theta);
      EvtSelEffCorrFactor.write(OutputFile, "z", &h1_EvtSel_Z, &h1_EvtSelBefore_Z);
      EvtSelEffCorrFactor.write(OutputFile, "theta", &h2_EvtSel_Theta, &h2_EvtSelBefore_Theta);
      EvtSelEffCorrFactor.write(OutputFile, "z", &h2_EvtSel_Z, &h2_EvtSelBefore_Z);

      OutputFile.Close();
   }
   // cleanup
   InputFile.Close();

   return 0;
}

int FindBin(double Value, int NBins, double Bins[])
{
   for(int i = 0; i < NBins; i++)
      if(Value < Bins[i])
         return i - 1;
   return NBins;
}

void MakeCanvasZ(vector<TH1D>& Histograms, vector<string> Labels, string Output,
   string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX){


   int NLine = Histograms.size();
   int N = Histograms[0].GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 0 : 0;
   double WorldXMax = LogX ? N : 1;

   double PadWidth = 1200;
   double PadHeight = DoRatio ? 640 : 640 + 240;
   double PadRHeight = DoRatio ? 240 : 0.001;

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
   Canvas.cd();

   TPad Pad("Pad", "", MarginL, MarginB + PadRHeight, MarginL + PadWidth, MarginB + PadHeight + PadRHeight);
   Pad.SetLogy();
   SetPad(Pad);

   TPad PadR("PadR", "", MarginL, MarginB, MarginL + PadWidth, MarginB + PadRHeight);
   if(DoRatio)
      SetPad(PadR);

   Pad.cd();

   TH2D HWorld("HWorld", "", N, WorldXMin, WorldXMax, 100, WorldMin, WorldMax);
   HWorld.SetStats(0);
   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);

   HWorld.Draw("axis");
   for(TH1D H : Histograms){
      TH1D *HClone = (TH1D *)H.Clone();
      HClone->Draw("exp same");
   }

   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : 1 / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : 1/ 2, 10000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   if(DoRatio)
      PadR.cd();

   double WorldRMin = 0.75;
   double WorldRMax = 1.25;

   TH2D HWorldR("HWorldR", "", N, WorldXMin, WorldXMax, 100, WorldRMin, WorldRMax);
   TGraph G2;

   if(DoRatio)
   {
      HWorldR.SetStats(0);
      HWorldR.GetXaxis()->SetTickLength(0);
      HWorldR.GetXaxis()->SetLabelSize(0);
      HWorldR.GetYaxis()->SetNdivisions(505);

      HWorldR.Draw("axis");
      for(int i = 1; i < NLine; i++)
      {
         TH1D *H = (TH1D *)Histograms[i].Clone();
         H->Divide(&Histograms[0]);
         H->Draw("same");
      }

      G.Draw("l");

      G2.SetPoint(0, 0, 1);
      G2.SetPoint(1, 99999, 1);
      G2.Draw("l");
   }

   double BinMin    = (1- cos(0.002))/2;
   double BinMiddle = 0.5;
   double BinMax    = 1 - BinMin;

   Canvas.cd();
   std::cout << "Here 1" << std::endl;
   int nDiv = 505;
   std::cout << "Bin Min " << BinMin << " Bin Middle " << BinMiddle << std::endl;
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   // axis on the x axis on the left hand side for the top of the axis 
   TGaxis X5(MarginL, MarginB + PadHeight + PadRHeight, MarginL + PadWidth / 2, MarginB + PadHeight + PadRHeight, BinMin, BinMiddle, 510, "-GS"); // - in the draw options means we only draw axis on the "negative" side
   // axis on the x axis on the right hand side for the top of the plot
   TGaxis X6(MarginL + PadWidth, MarginB + PadHeight + PadRHeight, MarginL + PadWidth / 2, MarginB + PadHeight + PadRHeight, BinMin, BinMiddle, 510, "+GS");// - in the draw options means we only draw axis on the "negative" side

   std::cout << "Here 2" << std::endl;
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");

   std::cout << "Here 3" << std::endl;
   TGaxis XL1(MarginL, MarginB, MarginL + PadWidth, MarginB,  (1- cos(0.002))/2, 1, 210, "S");
   TGaxis XL2(MarginL, MarginB + PadRHeight, MarginL + PadWidth, MarginB + PadRHeight,  (1- cos(0.002))/2, 1, 210, "+-S");

   Y1.SetLabelFont(42);
   Y2.SetLabelFont(42);
   XL1.SetLabelFont(42);
   XL2.SetLabelFont(42);

   X1.SetLabelSize(0);
   X2.SetLabelSize(0);
   X3.SetLabelSize(0);
   X4.SetLabelSize(0);
   X5.SetLabelSize(0);
   X6.SetLabelSize(0);
   // XL1.SetLabelSize(0);
   XL2.SetLabelSize(0);

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
      if(DoRatio) X3.Draw();
      if(DoRatio) X4.Draw();
      X5.Draw(); 
      X6.Draw(); 
   }
   if(LogX == false)
   {
      XL1.Draw();
      if(DoRatio)
         XL2.Draw();
   }
   if(DoRatio)
      Y1.Draw();
   Y2.Draw();

   TLatex Latex;
   Latex.SetNDC();
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.035);
   Latex.SetTextAlign(23);
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.02, MarginB - 0.01, "10^{-6} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.180, MarginB - 0.01, "10^{-4} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.365, MarginB - 0.01, "10^{-2} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.500, MarginB - 0.01, "1/2");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.650, MarginB - 0.01, "1 - 10^{-2}");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.823, MarginB - 0.01, "1 - 10^{-4}");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.995, MarginB - 0.01, "1 - 10^{-6}");

   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kGray);
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.04, "#it{z} = 1/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.9, MarginB * 0.4, X.c_str());

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(90);
   Latex.SetTextColor(kBlack);
   if(DoRatio)
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "Event Selection Eff.");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV, Work-in-progress");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "Finalization of EEC Results - April 24th (HB)");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.04 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
      Legend.AddEntry(&Histograms[i], Labels[i].c_str(), "pl");
   Legend.Draw();

   TLegend Legend2(0.55, 0.90, 0.8, 0.90 - 0.04 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(&Histograms[i], Labels[i].c_str(), "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}




void SetPad(TPad &P){
   P.SetLeftMargin(0);
   P.SetTopMargin(0);
   P.SetRightMargin(0);
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


void MakeCanvas(vector<TH1D>& Histograms, vector<string> Labels, string Output,
   string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0].GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 0 : 0;
   double WorldXMax = LogX ? N : M_PI;

   double PadWidth = 1200;
   double PadHeight = DoRatio ? 640 : 640 + 240;
   double PadRHeight = DoRatio ? 240 : 0.001;

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

   TPad PadR("PadR", "", MarginL, MarginB, MarginL + PadWidth, MarginB + PadRHeight);
   if(DoRatio)
      SetPad(PadR);

   Pad.cd();

   TH2D HWorld("HWorld", "", N, WorldXMin, WorldXMax, 100, WorldMin, WorldMax);
   HWorld.SetStats(0);
   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);

   HWorld.Draw("axis");
   for(TH1D H : Histograms){
      TH1D *H1 = (TH1D *)H.Clone();
      H1->Draw("same");
   }

   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : M_PI / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : M_PI / 2, 1000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   if(DoRatio)
      PadR.cd();

   double WorldRMin = 0.75;
   double WorldRMax = 1.25;

   TH2D HWorldR("HWorldR", "", N, WorldXMin, WorldXMax, 100, WorldRMin, WorldRMax);
   TGraph G2;

   if(DoRatio)
   {
      HWorldR.SetStats(0);
      HWorldR.GetXaxis()->SetTickLength(0);
      HWorldR.GetXaxis()->SetLabelSize(0);
      HWorldR.GetYaxis()->SetNdivisions(505);

      HWorldR.Draw("axis");
      for(int i = 1; i < NLine; i++)
      {
         TH1D *H = (TH1D *)Histograms[i].Clone();
         H->Divide(&Histograms[0]);
         H->Draw("same");
      }

      G.Draw("l");

      G2.SetPoint(0, 0, 1);
      G2.SetPoint(1, 99999, 1);
      G2.Draw("l");
   }

   double BinMin    = 0.002;
   double BinMiddle = M_PI / 2;
   double BinMax    = M_PI - 0.002;

   Canvas.cd();
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, 510, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, 510, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   // axis on the x axis on the left hand side for the top of the axis 
   TGaxis X5(MarginL, MarginB + PadHeight + PadRHeight, MarginL + PadWidth / 2, MarginB + PadHeight + PadRHeight, BinMin, BinMiddle, 510, "-GS"); // - in the draw options means we only draw axis on the "negative" side
   // axis on the x axis on the right hand side for the top of the plot
   TGaxis X6(MarginL + PadWidth, MarginB + PadHeight + PadRHeight, MarginL + PadWidth / 2, MarginB + PadHeight + PadRHeight, BinMin, BinMiddle, 510, "+GS");// - in the draw options means we only draw axis on the "negative" side
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");

   TGaxis XL1(MarginL, MarginB, MarginL + PadWidth, MarginB, 0, M_PI, 510, "S");
   TGaxis XL2(MarginL, MarginB + PadRHeight, MarginL + PadWidth, MarginB + PadRHeight, 0, M_PI, 510, "+-S");

   Y1.SetLabelFont(42);
   Y2.SetLabelFont(42);
   XL1.SetLabelFont(42);
   XL2.SetLabelFont(42);

   X1.SetLabelSize(0);
   X2.SetLabelSize(0);
   X3.SetLabelSize(0);
   X4.SetLabelSize(0);
   X5.SetLabelSize(0);
   X6.SetLabelSize(0);
   // XL1.SetLabelSize(0);
   XL2.SetLabelSize(0);

   X1.SetTickSize(0.06);
   X2.SetTickSize(0.06);
   X3.SetTickSize(0.06);
   X4.SetTickSize(0.06);
   X5.SetTickSize(0.06);
   X6.SetTickSize(0.06);
   XL1.SetTickSize(0.03);
   XL2.SetTickSize(0.03);

   if(LogX == true)
   {
      X1.Draw();
      X2.Draw();
      if(DoRatio) X3.Draw();
      if(DoRatio) X4.Draw();
      X5.Draw(); 
      X6.Draw(); 
   }
   if(LogX == false)
   {
      XL1.Draw();
      if(DoRatio)
         XL2.Draw();
   }
   if(DoRatio)
      Y1.Draw();
   Y2.Draw();

   TLatex Latex;
   Latex.SetNDC();
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.035);
   Latex.SetTextAlign(23);
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.115, MarginB - 0.01, "0.01");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.290, MarginB - 0.01, "0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.465, MarginB - 0.01, "1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.535, MarginB - 0.01, "#pi - 1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.710, MarginB - 0.01, "#pi - 0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.885, MarginB - 0.01, "#pi - 0.01");

   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kGray);
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.04, "#theta_{L} = #pi/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.5, MarginB * 0.3, X.c_str());

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(90);
   Latex.SetTextColor(kBlack);
   if(DoRatio)
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "Ratio");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV, Work-in-progress");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "Finalization of Results - April 24th (HB)");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.06 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
      Legend.AddEntry(&Histograms[i], Labels[i].c_str(), "pl");
   Legend.Draw();

   TLegend Legend2(0.7, 0.90, 0.9, 0.90 - 0.06 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(&Histograms[i], Labels[i].c_str(), "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}