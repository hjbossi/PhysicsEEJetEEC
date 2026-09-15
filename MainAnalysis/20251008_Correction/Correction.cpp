#include <iostream>
#include <vector>
#include <map>
using namespace std;

// root includes
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"

#include "Messenger.h"
#include "CommandLine.h"
#include "ProgressBar.h"
#include "TauHelperFunctions3.h"
#include "SetStyle.h"
#include "EffCorrFactor.h"

int FindBin(double Value, int NBins, double Bins[]);
double FindBinFraction(double Value, int NBins, double Bins[]); 
void MakeCanvasZ(vector<TH1D > Histograms, TGraphErrors DataSyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX);
void SetPad(TPad &P);
void DivideByBin(TH1D &H, double Bins[]);
TGraphAsymmErrors getTheoryPlot();
TGraphAsymmErrors getTheoryPlotUpdated();
TGraphAsymmErrors getTheoryPlotTheta(); 
void MakeCanvasZTheory(vector<TH1D > Histograms, TGraphErrors DataSyst,  TGraphAsymmErrors TheorySyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX);
void MakeCanvas(vector<TH1D> Histograms,  TGraphErrors DataSyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX);
void MakeCanvasRatioOnly(vector<TH1D > Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool LogX);
void MakeCanvasTheory(vector<TH1D > Histograms, TGraphErrors DataSyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX); 
// JC: the projection is copied here just so a quick check of the projected 1D histogram is accessible
void projection(TH2D* h_2D, TH1D* h_1D)
{
   h_1D->Reset();
   for (int i = 1; i <= h_2D->GetNbinsX(); ++i) {
      double weight = 0;
      double error = 0;
      for (int j = 1; j <= h_2D->GetNbinsY(); ++j) {
         double binContent = h_2D->GetBinContent(i, j);
         double binError= h_2D->GetBinError(i,j);
         double binCenter = h_2D->GetYaxis()->GetBinCenter(j);
         weight += binContent*((binCenter));
         error += pow(binError*binCenter, 2);;

      }
      h_1D->SetBinContent(i, weight);
      h_1D->SetBinError(i, sqrt(error));
   }
}

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   SetThesisStyle();
   static vector<int> Colors = GetCVDColors10();

   string InputDataPath          = CL.Get("InputData");
   string HistoName              = CL.Get("HistoName", "Bayesian_Unfoldediter4_Z");
   string HistoNameTheta         = CL.Get("HistoNameTheta", "Bayesian_Unfoldediter4_Theta");
   string CovBaseName            = CL.Get("CovBaseName", "cov_Unfolded_FixedN_");
   string OutputFileName         = CL.Get("Output", "CorrectedData.root");
   // For tgenBefore comparison
   string InputMCPath            = CL.Get("InputMC", "/data/ALEPH/MC/LEP1MC/LEP1MCMerged.root");
   string InputSHERPAPath          = CL.Get("InputHERWIG", "/data/ALEPH/MC/SHERPA/Sherpa_RNG100_0_0.root");
   string InputHERWIGPath          = CL.Get("InputHERWIG", "/data/ALEPH/MC/HERWIG/LEP-Matchbox-S1000-1_0_0.root");
   string InputPYTHIA8Path           = CL.Get("InputPYTHIA8", "/data/ALEPH/MC/PYTHI8Gen/LEP1_PYTHIA8_MC_TGenBefore.root");
   string InputPYTHIA8VinciaPath     = CL.Get("InputPYTHIA8Vincia", "/data/ALEPH/MC/PYTHI8Gen/LEP1_PYTHIA8_MC_TGenBefore_VINCIA.root");
   string InputPYTHIA8DirePath       = CL.Get("InputPYTHIA8Dire", "/data/ALEPH/MC/PYTHI8Gen/LEP1_PYTHIA8_MC_TGenBefore_DIRE.root");

   string GenBeforeTreeName      = CL.Get("GenBefore", "tgenBefore");
   string GenTreeName            = CL.Get("Gen","tgen");
   string RecoTreeName           = CL.Get("Reco", "t"); // used for herwig and sherpa

   TFile InputData(InputDataPath.c_str(), "READ");
   TH2D HDataBfMatchingCorr( *((TH2D*) InputData.Get(HistoName.c_str())) );
   TH2D HDataAfMatchingCorr( *((TH2D*) HDataBfMatchingCorr.Clone(Form("%s_afCorrZ", HistoName.c_str()))) );
   TH2D HDataBfMatchingCorr_Theta( *((TH2D*) InputData.Get(HistoNameTheta.c_str())) );
   TH2D HDataAfMatchingCorr_Theta( *((TH2D*) HDataBfMatchingCorr_Theta.Clone(Form("%s_afCorr", HistoNameTheta.c_str()))) );

   // HDataBfCorr.Print("all");

   int applyEffCorrOnHistoErrorStatus = 0;
   int applyEffCorrOnHistoErrorStatus_theta = 0;
   // apply the matching efficiency correction
   std::cout << "------------- Efficiency Correction: [Applying the matching efficiency] ----------" << std::endl;
   EffCorrFactor matchingEffCorrFactor;
   matchingEffCorrFactor.init("/home/hbossi/PhysicsEEJetEEC/Unfolding/20250317_Unfolding/matchingScheme2/MatchingEff.root", "z");
   applyEffCorrOnHistoErrorStatus += matchingEffCorrFactor.applyEffCorrOnHisto(&HDataBfMatchingCorr, &HDataAfMatchingCorr);
   EffCorrFactor matchingEffCorrFactorTheta;
   matchingEffCorrFactorTheta.init("/home/hbossi/PhysicsEEJetEEC/Unfolding/20250317_Unfolding/matchingScheme2/MatchingEff.root", "theta");
   applyEffCorrOnHistoErrorStatus_theta += matchingEffCorrFactorTheta.applyEffCorrOnHisto(&HDataBfMatchingCorr_Theta, &HDataAfMatchingCorr_Theta);
   // done applying the matchig efficiency correction


   TH2D HDataBfCorr( *((TH2D*) HDataAfMatchingCorr.Clone(Form("%s_afMatchingCorrZ", HistoName.c_str()))) );
   TH2D HDataAfCorr( *((TH2D*) HDataAfMatchingCorr.Clone(Form("%s_afAllZ", HistoName.c_str()))) );
   TH2D HDataBfCorr_Theta( *((TH2D*) HDataAfMatchingCorr_Theta.Clone(Form("%s_afMatchingCorr", HistoNameTheta.c_str()))) );
   TH2D HDataAfCorr_Theta( *((TH2D*) HDataAfMatchingCorr_Theta.Clone(Form("%s_afAll", HistoNameTheta.c_str()))) );

   
   TH2D *HCovZ     = (TH2D*) InputData.Get(Form("%sZ",     CovBaseName.c_str()));
   TH2D *HCovTheta = (TH2D*) InputData.Get(Form("%sTheta", CovBaseName.c_str()));

   // print out error message if the covariance matrix is not there
   if (!HCovZ)     { std::cerr << "[Error] Could not load " << CovBaseName << "Z from "     << InputDataPath << std::endl; return 1; }
   if (!HCovTheta) { std::cerr << "[Error] Could not load " << CovBaseName << "Theta from " << InputDataPath << std::endl; return 1; }

   // Save pre-correction 1D projections for z and theta.
   // These are used later to compute the total per-bin correction
   // factor needed to propagate the covariance-derived stat uncertainty.
   TH1D *HDataRaw1D_z     = (TH1D*) HDataBfMatchingCorr.ProjectionX("HDataRaw1D_z");
   TH1D *HDataRaw1D_Theta = (TH1D*) HDataBfMatchingCorr_Theta.ProjectionX("HDataRaw1D_Theta");
   projection(&HDataBfMatchingCorr,       HDataRaw1D_z);
   projection(&HDataBfMatchingCorr_Theta, HDataRaw1D_Theta);

   
   
   
   std::cout << "------------- Efficiency Correction: [Applying the event selection efficiency] ----------" << std::endl;
   EffCorrFactor EvtSelEffCorrFactor;
   EvtSelEffCorrFactor.init("/home/hbossi/PhysicsEEJetEEC/EventSelectionEfficiency/20250324_evtSelEffCorr/EvtSelEff.root", "z");
   applyEffCorrOnHistoErrorStatus += EvtSelEffCorrFactor.applyEffCorrOnHisto(&HDataBfCorr, &HDataAfCorr);

   EffCorrFactor EvtSelEffCorrFactor_Theta;
   EvtSelEffCorrFactor_Theta.init("/home/hbossi/PhysicsEEJetEEC/EventSelectionEfficiency/20250324_evtSelEffCorr/EvtSelEff.root", "theta");
   applyEffCorrOnHistoErrorStatus_theta += EvtSelEffCorrFactor_Theta.applyEffCorrOnHisto(&HDataBfCorr_Theta, &HDataAfCorr_Theta);

   TFile Output(OutputFileName.c_str(), "RECREATE");
   Output.cd();
   HDataBfCorr.Write();
   HDataBfCorr_Theta.Write();
   HDataAfCorr.Write();
   HDataAfCorr_Theta.Write();

   // JC: the projection is copied here just so a quick check of the projected 1D histogram is accessible
   TH1D* HDataBfCorr1D = (TH1D*) HDataBfCorr.ProjectionX();
   TH1D* HDataAfCorr1D = (TH1D*) HDataAfCorr.ProjectionX();
   projection(&HDataBfCorr, HDataBfCorr1D);
   projection(&HDataAfCorr, HDataAfCorr1D);
   Output.cd();
   HDataBfCorr1D->Write();
   HDataAfCorr1D->Write();

   TH1D* HDataBfCorr1D_Theta = (TH1D*) HDataBfCorr_Theta.ProjectionX();
   TH1D* HDataAfCorr1D_Theta = (TH1D*) HDataAfCorr_Theta.ProjectionX();
   projection(&HDataBfCorr_Theta, HDataBfCorr1D_Theta);
   projection(&HDataAfCorr_Theta, HDataAfCorr1D_Theta);
   Output.cd();
   HDataBfCorr1D_Theta->Write();
   HDataAfCorr1D_Theta->Write();

   TH1D* HDataBfCorr1D_unfoldBinCorr = (TH1D*) HDataAfCorr1D->Clone(Form("%s_unfoldBinCorrx", HDataBfCorr1D->GetName()));
   TH1D* HDataAfCorr1D_unfoldBinCorr = (TH1D*) HDataAfCorr1D->Clone(Form("%s_unfoldBinCorr", HDataAfCorr1D->GetName()));
   EffCorrFactor UnfoldingBinCorrFactor;
   UnfoldingBinCorrFactor.init("/home/hbossi/PhysicsEEJetEEC/Unfolding/20250324_UnfoldingBinningCorrection/UnfoldingBinCorr_with_z.root", "z");
   //applyEffCorrOnHistoErrorStatus += UnfoldingBinCorrFactor.applyEffCorrOnHisto(HDataBfCorr1D, HDataBfCorr1D_unfoldBinCorr);
   applyEffCorrOnHistoErrorStatus += UnfoldingBinCorrFactor.applyEffCorrOnHisto(HDataBfCorr1D_unfoldBinCorr, HDataAfCorr1D_unfoldBinCorr);
   Output.cd();
   HDataBfCorr1D_unfoldBinCorr->Write();
   HDataAfCorr1D_unfoldBinCorr->Write();


   TH1D* HDataBfCorr1D_unfoldBinCorr_Theta = (TH1D*) HDataAfCorr1D_Theta->Clone(Form("%s_unfoldBinCorrTheta", HDataBfCorr1D_Theta->GetName()));
   TH1D* HDataAfCorr1D_unfoldBinCorr_Theta = (TH1D*) HDataAfCorr1D_Theta->Clone(Form("%s_unfoldBinCorrTheta", HDataAfCorr1D_Theta->GetName()));
   EffCorrFactor UnfoldingBinCorrFactor_Theta;
   UnfoldingBinCorrFactor_Theta.init("/home/hbossi/PhysicsEEJetEEC/Unfolding/20250324_UnfoldingBinningCorrection/UnfoldingBinCorr_with_theta.root", "theta");
   //applyEffCorrOnHistoErrorStatus_theta += UnfoldingBinCorrFactor_Theta.applyEffCorrOnHisto(HDataBfCorr1D_Theta, HDataBfCorr1D_unfoldBinCorr_Theta);
   applyEffCorrOnHistoErrorStatus_theta += UnfoldingBinCorrFactor_Theta.applyEffCorrOnHisto(HDataBfCorr1D_unfoldBinCorr_Theta, HDataAfCorr1D_unfoldBinCorr_Theta);
   Output.cd();
   HDataBfCorr1D_unfoldBinCorr_Theta->Write();
   HDataAfCorr1D_unfoldBinCorr_Theta->Write();

   //------------------------------------
   // plotting
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

   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);

      //std::cout << "i = " << i << " Bins[i] = " << Bins[i] << " zBins[i] " << zBins[i] << " " <<  (1- cos(Bins[i]))/2 << std::endl;

   }

   // [Warning] JC: this is a quick hack to get the normalization, we would need to parse that normalization value
  TString fnamesmeared = "/home/hbossi/PhysicsEEJetEEC/Unfolding/20250303_FakeCorrection/UnfoldingInputData_03192025.root";
  TFile *inputsmeared =TFile::Open(fnamesmeared);
  TTree *smeared=(TTree*)inputsmeared->Get("UnmatchedPairTree");
  Int_t nEv=smeared->GetEntries();
  
  std::cout << "Number of events " << nEv << std::endl;

   HDataBfCorr1D->Scale(1.0/nEv);
   HDataAfCorr1D->Scale(1.0/nEv);
   HDataBfCorr1D_unfoldBinCorr->Scale(1.0/nEv);
   HDataAfCorr1D_unfoldBinCorr->Scale(1.0/nEv);

   HDataBfCorr1D_Theta->Scale(1.0/nEv);
   HDataAfCorr1D_Theta->Scale(1.0/nEv);
   HDataBfCorr1D_unfoldBinCorr_Theta->Scale(1.0/nEv);
   HDataAfCorr1D_unfoldBinCorr_Theta->Scale(1.0/nEv);

   //------------------------------------
   // getting the EEC(z) at the gen level before the event selection
   //------------------------------------
   TFile InputMC(InputMCPath.c_str(), "READ");
   TH1D HzMCGenBeforeRef("HzMCGenBeforeRef", "HzMCGenBeforeRef", 2 * BinCount, 0, 2 * BinCount);
   TH1D HthetaMCGenBeforeRef("HthetaMCGenBeforeRef", "HthetaMCGenBeforeRef", 2 * BinCount, 0, 2 * BinCount);

   double TotalE = 91.1876;
   ParticleTreeMessenger MGenBefore(InputMC, GenBeforeTreeName);

   int EntryCountBefore = MGenBefore.GetEntries();
   for(int iE = 0; iE < EntryCountBefore; iE++){
      MGenBefore.GetEntry(iE);
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

         //if(MGenBefore.highPurity[i] == false) continue;
         // place cut on the reco energy, not included at gen level
         PGenBefore.push_back(MGenBefore.P[i]);
         pwflagVec.push_back(MGenBefore.pwflag[i]); 
      }


      for(int i = 0; i < PGenBefore.size(); i++){
         for(int j = i+1; j < PGenBefore.size();j++){
            if(i == j) continue; // don't fill the EEC with self correlations particles with themselves
            FourVector Gen1 = PGenBefore.at(i);
            FourVector Gen2 = PGenBefore.at(j);
            double genTheta = GetAngle(Gen1,Gen2);
            double genZ = (1-cos(genTheta))/2;
            int BinThetaGen = FindBin(genTheta, 2 * BinCount, Bins);
            int BinZGen = FindBin(genZ, 2*BinCount, zBins);
            double  genEEC  = Gen1[0]*Gen2[0]/(TotalE*TotalE);
            HthetaMCGenBeforeRef.Fill(BinThetaGen, genEEC);
            HzMCGenBeforeRef.Fill(BinZGen, genEEC);

         }
      }
    } // end loop over the reco tree

   // EEC is per-event so scale by the event number
   HzMCGenBeforeRef.Scale(1.0/EntryCountBefore);
   HthetaMCGenBeforeRef.Scale(1.0/EntryCountBefore);
   
   // write out the error for the sake of the systematic uncertainty
   Output.cd();
   HthetaMCGenBeforeRef.Write(); 
   HzMCGenBeforeRef.Write(); 
   
   TH1D HzMCGenPYTHIA8Dire("HzMCGenPYTHIA8Dire", "HzMCGenPYTHIA8Dire", 2 * BinCount, 0, 2 * BinCount);
   TH1D HthetaMCGenPYTHIA8Dire("HthetaMCGenPYTHIA8Dire", "HthetaMCGenPYTHIA8Dire", 2 * BinCount, 0, 2 * BinCount);
     TH1D HzMCGenSHERPA("HzMCGenSHERPA", "HzMCGenSHERPA", 2 * BinCount, 0, 2 * BinCount);
   TH1D HthetaMCGenSHERPA("HthetaMCGenSHERPA", "HthetaMCGenSHERPA", 2 * BinCount, 0, 2 * BinCount);
      TH1D HzMCGenHERWIG("HzMCGenHERWIG", "HzMCGenHERWIG", 2 * BinCount, 0, 2 * BinCount);
   TH1D HthetaMCGenHERWIG("HthetaMCGenHERWIG", "HthetaMCGenHERWIG", 2 * BinCount, 0, 2 * BinCount);
      TH1D HzMCGenPYTHIA8("HzMCGenPYTHIA8", "HzMCGenPYTHIA8", 2 * BinCount, 0, 2 * BinCount);
   TH1D HthetaMCGenPYTHIA8("HthetaMCGenPYTHIA8", "HthetaMCGenPYTHIA8", 2 * BinCount, 0, 2 * BinCount);
      TH1D HzMCGenPYTHIA8Vincia("HzMCGenPYTHIA8Vincia", "HzMCGenPYTHIA8Vincia", 2 * BinCount, 0, 2 * BinCount);
   TH1D HthetaMCGenPYTHIA8Vincia("HthetaMCGenPYTHIA8Vincia", "HthetaMCGenPYTHIA8Vincia", 2 * BinCount, 0, 2 * BinCount);
   
   //------------------ SHERPA  ---------------------------
   TFile InputSHERPA(InputSHERPAPath.c_str(), "READ");

   ParticleTreeMessenger MGenSHERPA(InputSHERPA, RecoTreeName);
   int EntryCountSHERPA = MGenSHERPA.GetEntries();
   for(int iE = 0; iE < EntryCountSHERPA; iE++)
   {
      MGenSHERPA.GetEntry(iE);
      // fill the four vector
      vector<FourVector> PGenBefore;
      //std::cout << MGenSHERPA.nParticle << std::endl;
      for(int i = 0; i < MGenSHERPA.nParticle; i++){
        // charged particle selection
        //std::cout << "MGenSHERPA.charge[i]: " << MGenSHERPA.charge[i] << std::endl;
       if(MGenSHERPA.charge[i] == 0) continue;
         PGenBefore.push_back(MGenSHERPA.P[i]);
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
            HzMCGenSHERPA.Fill(BinZGen, EEC);
            HthetaMCGenSHERPA.Fill(BinThetaGen, EEC);
         }
      }
   } // end loop over the number of events
   // EEC is per-event so scale by the event number

   HzMCGenSHERPA.Scale(1.0/EntryCountSHERPA);
   HthetaMCGenSHERPA.Scale(1.0/EntryCountSHERPA);

   
   //------------------ HERWIG  ---------------------------
   TFile InputHERWIG(InputHERWIGPath.c_str(), "READ");

   ParticleTreeMessenger MGenHERWIG(InputHERWIG, RecoTreeName);
   int EntryCountHERWIG = MGenHERWIG.GetEntries();
   for(int iE = 0; iE < EntryCountHERWIG; iE++)
   {
      MGenHERWIG.GetEntry(iE);
      // fill the four vector
      vector<FourVector> PGenBefore;
      //std::cout << MGenHERWIG.nParticle << std::endl;
      for(int i = 0; i < MGenHERWIG.nParticle; i++){
        // charged particle selection
        //std::cout << "MGenHERWIG.charge[i]: " << MGenHERWIG.charge[i] << std::endl;
       if(MGenHERWIG.charge[i] == 0) continue;
         PGenBefore.push_back(MGenHERWIG.P[i]);
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
            HzMCGenHERWIG.Fill(BinZGen, EEC);
            HthetaMCGenHERWIG.Fill(BinThetaGen, EEC);
         }
      }
   } // end loop over the number of events
   // EEC is per-event so scale by the event number

   HzMCGenHERWIG.Scale(1.0/EntryCountHERWIG);
   HthetaMCGenHERWIG.Scale(1.0/EntryCountHERWIG);

   //------------------ PYTHIA 8 ---------------------------
   TFile InputPYTHIA8(InputPYTHIA8Path.c_str(), "READ");

   ParticleTreeMessenger MGenPYTHIA8(InputPYTHIA8, GenBeforeTreeName);
   int EntryCountPYTHIA8 = MGenPYTHIA8.GetEntries();
   for(int iE = 0; iE < EntryCountPYTHIA8; iE++)
   {
      MGenPYTHIA8.GetEntry(iE);
      // fill the four vector
      vector<FourVector> PGenBefore;
      for(int i = 0; i < MGenPYTHIA8.nParticle; i++){
        // charged particle selection
       if(MGenPYTHIA8.isCharged[i] == 0) continue;
         PGenBefore.push_back(MGenPYTHIA8.P[i]);
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
            HzMCGenPYTHIA8.Fill(BinZGen, EEC);
            HthetaMCGenPYTHIA8.Fill(BinThetaGen, EEC);
         }
      }
   } // end loop over the number of events
   // EEC is per-event so scale by the event number

   HzMCGenPYTHIA8.Scale(1.0/EntryCountPYTHIA8);
   HthetaMCGenPYTHIA8.Scale(1.0/EntryCountPYTHIA8);


   //------------------ PYTHIA 8 Vincia ---------------------------
   TFile InputPYTHIA8Vincia(InputPYTHIA8VinciaPath.c_str(), "READ");

   ParticleTreeMessenger MGenPYTHIA8Vincia(InputPYTHIA8Vincia, GenBeforeTreeName);
   int EntryCountPYTHIA8Vincia = MGenPYTHIA8Vincia.GetEntries();
   for(int iE = 0; iE < EntryCountPYTHIA8Vincia; iE++)
   {
      MGenPYTHIA8Vincia.GetEntry(iE);
      // fill the four vector
      vector<FourVector> PGenBefore;
      for(int i = 0; i < MGenPYTHIA8Vincia.nParticle; i++){
        // charged particle selection
       if(MGenPYTHIA8Vincia.isCharged[i] == 0) continue;
         PGenBefore.push_back(MGenPYTHIA8Vincia.P[i]);
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
            HzMCGenPYTHIA8Vincia.Fill(BinZGen, EEC);
            HthetaMCGenPYTHIA8Vincia.Fill(BinThetaGen, EEC);
         }
      }
   } // end loop over the number of events
   // EEC is per-event so scale by the event number

   HzMCGenPYTHIA8Vincia.Scale(1.0/EntryCountPYTHIA8Vincia);
   HthetaMCGenPYTHIA8Vincia.Scale(1.0/EntryCountPYTHIA8Vincia);

   //-------------------------------------------------

    //------------------ PYTHIA 8 Dire ---------------------------
   TFile InputPYTHIA8Dire(InputPYTHIA8DirePath.c_str(), "READ");
   ParticleTreeMessenger MGenPYTHIA8Dire(InputPYTHIA8Dire, GenBeforeTreeName);
   int EntryCountPYTHIA8Dire = MGenPYTHIA8Dire.GetEntries();
   for(int iE = 0; iE < EntryCountPYTHIA8Dire; iE++)
   {
      MGenPYTHIA8Dire.GetEntry(iE);
      // fill the four vector
      vector<FourVector> PGenBefore;
      for(int i = 0; i < MGenPYTHIA8Dire.nParticle; i++){
        // charged particle selection
       if(MGenPYTHIA8Dire.isCharged[i] == 0) continue;
         PGenBefore.push_back(MGenPYTHIA8Dire.P[i]);
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
            HzMCGenPYTHIA8Dire.Fill(BinZGen, EEC);
            HthetaMCGenPYTHIA8Dire.Fill(BinThetaGen, EEC);
         }
      }
   } // end loop over the number of events
   // EEC is per-event so scale by the event number

   HzMCGenPYTHIA8Dire.Scale(1.0/EntryCountPYTHIA8Dire);
   HthetaMCGenPYTHIA8Dire.Scale(1.0/EntryCountPYTHIA8Dire);
   //-------------------------------------------------

   std::cout << "zBins 17: " << zBins[17] << " zbins 183: " << zBins[183] << std::endl;
   std::cout << "theta 17: " << Bins[17] << " theta Bins 183: " << Bins[183] << std::endl;
   std::cout << "Sanity check:  " << (1-cos(Bins[17]))/2 <<std::endl;

   // divide by the bin width
   DivideByBin(*HDataBfCorr1D, zBins);
   DivideByBin(*HDataAfCorr1D, zBins);
   DivideByBin(*HDataBfCorr1D_unfoldBinCorr, zBins);
   DivideByBin(*HDataAfCorr1D_unfoldBinCorr, zBins);
   DivideByBin(HzMCGenBeforeRef, zBins);
   DivideByBin(HzMCGenPYTHIA8, zBins);
   DivideByBin(HzMCGenPYTHIA8Vincia, zBins);
   DivideByBin(HzMCGenPYTHIA8Dire, zBins);
   DivideByBin(HzMCGenHERWIG, zBins);
   DivideByBin(HzMCGenSHERPA, zBins);


   DivideByBin(*HDataBfCorr1D_Theta, Bins);
   DivideByBin(*HDataAfCorr1D_Theta, Bins);
   DivideByBin(*HDataBfCorr1D_unfoldBinCorr_Theta, Bins);
   DivideByBin(*HDataAfCorr1D_unfoldBinCorr_Theta, Bins);
   DivideByBin(HthetaMCGenBeforeRef, Bins);
   DivideByBin(HthetaMCGenPYTHIA8, Bins);
   DivideByBin(HthetaMCGenPYTHIA8Vincia, Bins);
   DivideByBin(HthetaMCGenPYTHIA8Dire, Bins);
   DivideByBin(HthetaMCGenHERWIG, Bins);
   DivideByBin(HthetaMCGenSHERPA, Bins);

   // set the style for the plots
   HDataBfCorr1D->SetMarkerColor(Colors[5]);
   HDataAfCorr1D->SetMarkerColor(Colors[3]);
   HDataBfCorr1D_unfoldBinCorr->SetMarkerColor(Colors[4]);
   HDataAfCorr1D_unfoldBinCorr->SetMarkerColor(kBlack);
   HzMCGenBeforeRef.SetMarkerColor(kAzure);
   HzMCGenPYTHIA8.SetMarkerColor(Colors[1]);
   HzMCGenHERWIG.SetMarkerColor(Colors[6]);
   HzMCGenSHERPA.SetMarkerColor(kMagenta-3);
   HzMCGenPYTHIA8Vincia.SetMarkerColor(Colors[4]);
   HzMCGenPYTHIA8Dire.SetMarkerColor(kGreen-2);
   HDataBfCorr1D->SetLineColor(Colors[5]);
   HDataAfCorr1D->SetLineColor(Colors[3]);
   HDataBfCorr1D_unfoldBinCorr->SetLineColor(Colors[4]);
   HDataAfCorr1D_unfoldBinCorr->SetLineColor(kBlack);
   HzMCGenBeforeRef.SetLineColor(kAzure);
   HzMCGenPYTHIA8.SetLineColor(Colors[1]);
   HzMCGenPYTHIA8Vincia.SetLineColor(Colors[4]);
   HzMCGenPYTHIA8Dire.SetLineColor(kGreen-2);
   HzMCGenHERWIG.SetLineColor(Colors[6]);
   HzMCGenSHERPA.SetLineColor(kMagenta-3);

   HDataBfCorr1D->SetMarkerStyle(20);
   HDataAfCorr1D->SetMarkerStyle(20);
   HDataBfCorr1D_unfoldBinCorr->SetMarkerStyle(20);
   HDataAfCorr1D_unfoldBinCorr->SetMarkerStyle(20);
   HzMCGenBeforeRef.SetMarkerStyle(20);
   HzMCGenPYTHIA8.SetMarkerStyle(21);
   HzMCGenPYTHIA8Vincia.SetMarkerStyle(33);
   HzMCGenPYTHIA8Vincia.SetMarkerSize(2);
   HzMCGenPYTHIA8Dire.SetMarkerStyle(25);
   HzMCGenHERWIG.SetMarkerStyle(27);
   HzMCGenHERWIG.SetMarkerSize(2);

   HzMCGenSHERPA.SetMarkerStyle(24);


   HDataBfCorr1D->SetLineWidth(2);
   HDataAfCorr1D->SetLineWidth(2);
   HDataBfCorr1D_unfoldBinCorr->SetLineWidth(2);
   HDataAfCorr1D_unfoldBinCorr->SetLineWidth(2);
   HzMCGenBeforeRef.SetLineWidth(2);
   HzMCGenPYTHIA8.SetLineWidth(2);
   HzMCGenPYTHIA8Vincia.SetLineWidth(2);
   HzMCGenPYTHIA8Dire.SetLineWidth(2);
   HzMCGenHERWIG.SetLineWidth(2);
   HzMCGenSHERPA.SetLineWidth(2);

  
   HDataBfCorr1D_Theta->SetMarkerColor(Colors[5]);
   HDataAfCorr1D_Theta->SetMarkerColor(Colors[3]);
   HDataBfCorr1D_unfoldBinCorr_Theta->SetMarkerColor(Colors[4]);
   HDataAfCorr1D_unfoldBinCorr_Theta->SetMarkerColor(kBlack);
   HthetaMCGenBeforeRef.SetMarkerColor(kAzure);
   HthetaMCGenPYTHIA8.SetMarkerColor(Colors[1]);
   HthetaMCGenPYTHIA8Vincia.SetMarkerColor(Colors[4]);
   HthetaMCGenPYTHIA8Dire.SetMarkerColor(kGreen-2);
   HthetaMCGenHERWIG.SetMarkerColor(Colors[6]);
   HthetaMCGenSHERPA.SetMarkerColor(kMagenta-3);


   HDataBfCorr1D_Theta->SetLineColor(Colors[5]);
   HDataAfCorr1D_Theta->SetLineColor(Colors[3]);
   HDataBfCorr1D_unfoldBinCorr_Theta->SetLineColor(Colors[4]);
   HDataAfCorr1D_unfoldBinCorr_Theta->SetLineColor(kBlack);
   HthetaMCGenBeforeRef.SetLineColor(kAzure);
   HthetaMCGenPYTHIA8.SetLineColor(Colors[1]);
   HthetaMCGenPYTHIA8Vincia.SetLineColor(Colors[4]);
   HthetaMCGenPYTHIA8Dire.SetLineColor(kGreen-2);
   HthetaMCGenHERWIG.SetLineColor(Colors[6]);
   HthetaMCGenSHERPA.SetLineColor(kMagenta-3);

   HDataBfCorr1D_Theta->SetMarkerStyle(20);
   HDataAfCorr1D_Theta->SetMarkerStyle(20);
   HDataBfCorr1D_unfoldBinCorr_Theta->SetMarkerStyle(20);
   HDataAfCorr1D_unfoldBinCorr_Theta->SetMarkerStyle(20);
   HthetaMCGenBeforeRef.SetMarkerStyle(20);
   HthetaMCGenPYTHIA8.SetMarkerStyle(21);
   HthetaMCGenHERWIG.SetMarkerStyle(27);
   HthetaMCGenHERWIG.SetMarkerSize(2);
   HthetaMCGenSHERPA.SetMarkerStyle(24);
   HthetaMCGenPYTHIA8Vincia.SetMarkerStyle(33);
   HthetaMCGenPYTHIA8Vincia.SetMarkerSize(2);
   HthetaMCGenPYTHIA8Dire.SetMarkerStyle(25);
   HDataBfCorr1D_Theta->SetLineWidth(2);
   HDataAfCorr1D_Theta->SetLineWidth(2);
   HDataBfCorr1D_unfoldBinCorr_Theta->SetLineWidth(2);
   HDataAfCorr1D_unfoldBinCorr_Theta->SetLineWidth(2);
   HthetaMCGenBeforeRef.SetLineWidth(2);
   HthetaMCGenPYTHIA8.SetLineWidth(2);
   HthetaMCGenPYTHIA8Vincia.SetLineWidth(2);
   HthetaMCGenPYTHIA8Dire.SetLineWidth(2);
   HthetaMCGenHERWIG.SetLineWidth(2);
   HthetaMCGenSHERPA.SetLineWidth(2);
   
   


   // Override bin errors for z using covariance matrix diagonal.
   for (int i = 1; i <= HDataAfCorr1D_unfoldBinCorr->GetNbinsX(); i++) {
      double raw     = HDataRaw1D_z->GetBinContent(i);
      double covDiag = HCovZ->GetBinContent(i, i);
      if (raw <= 0.0 || covDiag < 0.0) {
         HDataAfCorr1D_unfoldBinCorr->SetBinError(i, 0.0);
         continue;
      }
      double corrFactor = HDataAfCorr1D_unfoldBinCorr->GetBinContent(i) / raw;
      HDataAfCorr1D_unfoldBinCorr->SetBinError(i, std::sqrt(covDiag) * corrFactor);
   }
   
   std::cout << "[Debug] Z bin errors AFTER covariance override (first 10 bins):" << std::endl;
   for (int i = 1; i <= std::min(10, HDataAfCorr1D_unfoldBinCorr->GetNbinsX()); i++) {
      double raw      = HDataRaw1D_z->GetBinContent(i);
      double covDiag  = HCovZ->GetBinContent(i, i);
      double content  = HDataAfCorr1D_unfoldBinCorr->GetBinContent(i);
      double errAfter = HDataAfCorr1D_unfoldBinCorr->GetBinError(i);
      std::cout << "  bin " << i
               << "  raw=" << raw
               << "  covDiag=" << covDiag
               << "  content=" << content
               << "  err=" << errAfter
               << "  rel_err=" << (content > 0 ? errAfter/content : -1)
               << std::endl;
   }
      
   
   // do the same thing for theta
   for (int i = 1; i <= HDataAfCorr1D_unfoldBinCorr_Theta->GetNbinsX(); i++) {
      double raw     = HDataRaw1D_Theta->GetBinContent(i);
      double covDiag = HCovTheta->GetBinContent(i, i);
      if (raw <= 0.0 || covDiag < 0.0) {
         HDataAfCorr1D_unfoldBinCorr_Theta->SetBinError(i, 0.0);
         continue;
      }
      double corrFactor = HDataAfCorr1D_unfoldBinCorr_Theta->GetBinContent(i) / raw;
      HDataAfCorr1D_unfoldBinCorr_Theta->SetBinError(i, std::sqrt(covDiag) * corrFactor);
   }


   // sys Z
   TFile* sysFile = TFile::Open("/home/hbossi/PhysicsEEJetEEC/Systematics/20251030_SystematicsUpdate/SystematicsE2C_Z_10302025.root");
   TH2D* h2SysFromFile = (TH2D*)sysFile->Get("Systematics_Z_Total");
   TH2D* h2Sys = (TH2D*) HDataBfCorr.Clone("h2Sys");
   TH2D* h2SysPlus =(TH2D*) HDataBfCorr.Clone("h2SysPlus");


   h2Sys->Add(h2SysFromFile, -1);
   h2SysPlus->Add(h2SysFromFile);

   TH2D HSysAfCorr( *((TH2D*) h2Sys->Clone("h2Sys_afCorr") ));
   TH2D HSysAfCorrPlus( *((TH2D*) h2SysPlus->Clone("h2SysPlus_afCorr") ));

   matchingEffCorrFactor.applyEffCorrOnHisto(h2Sys, &HSysAfCorr);
   matchingEffCorrFactor.applyEffCorrOnHisto(h2SysPlus, &HSysAfCorrPlus);

   EvtSelEffCorrFactor.applyEffCorrOnHisto(h2Sys, &HSysAfCorr);
   EvtSelEffCorrFactor.applyEffCorrOnHisto(h2SysPlus, &HSysAfCorrPlus);

   // sys theta
   TFile* sysFileTheta = TFile::Open("/home/hbossi/PhysicsEEJetEEC/Systematics/20251030_SystematicsUpdate/SystematicsE2C_Theta_10302025.root");
   TH2D* h2SysFromFileTheta = (TH2D*)sysFileTheta->Get("Systematics_Theta_Total");
   TH2D* h2SysTheta = (TH2D*) HDataBfCorr_Theta.Clone("h2Sys");
   TH2D* h2SysPlusTheta =(TH2D*) HDataBfCorr_Theta.Clone("h2SysPlus");




   h2SysTheta->Add(h2SysFromFileTheta, -1);
   h2SysPlusTheta->Add(h2SysFromFileTheta);

   TH2D HSysAfCorrTheta( *((TH2D*) h2SysTheta->Clone("h2SysTheta_afCorr") ));
   TH2D HSysAfCorrPlusTheta( *((TH2D*) h2SysPlusTheta->Clone("h2SysPlusTheta_afCorr") ));

   matchingEffCorrFactorTheta.applyEffCorrOnHisto(h2SysTheta, &HSysAfCorrTheta);
   matchingEffCorrFactorTheta.applyEffCorrOnHisto(h2SysPlusTheta, &HSysAfCorrPlusTheta);

   EvtSelEffCorrFactor_Theta.applyEffCorrOnHisto(h2SysTheta, &HSysAfCorrTheta);
   EvtSelEffCorrFactor_Theta.applyEffCorrOnHisto(h2SysPlusTheta, &HSysAfCorrPlusTheta);


   // JC: the projection is copied here just so a quick check of the projected 1D histogram is accessible
   TH1D* HSysBfCorr1D = (TH1D*) h2Sys->ProjectionX();
   TH1D* HSysBfCorr1DPlus = (TH1D*) h2SysPlus->ProjectionX();
   TH1D* HSysAfCorr1D = (TH1D*) HSysAfCorr.ProjectionX();
   TH1D* HSysAfCorr1DPlus = (TH1D*) HSysAfCorrPlus.ProjectionX();

   TH1D* HSysBfCorr1DTheta = (TH1D*) h2SysTheta->ProjectionX();
   TH1D* HSysBfCorr1DPlusTheta = (TH1D*) h2SysPlusTheta->ProjectionX();
   TH1D* HSysAfCorr1DTheta = (TH1D*) HSysAfCorrTheta.ProjectionX();
   TH1D* HSysAfCorr1DPlusTheta = (TH1D*) HSysAfCorrPlusTheta.ProjectionX();

   projection(h2Sys, HSysBfCorr1D);
   projection(h2SysPlus, HSysBfCorr1DPlus);

   projection(&HSysAfCorr, HSysAfCorr1D);
   projection(&HSysAfCorrPlus, HSysAfCorr1DPlus);

   projection(h2SysTheta, HSysBfCorr1DTheta);
   projection(h2SysPlusTheta, HSysBfCorr1DPlusTheta);

   projection(&HSysAfCorrTheta, HSysAfCorr1DTheta);
   projection(&HSysAfCorrPlusTheta, HSysAfCorr1DPlusTheta);



   TH1D* HSysBfCorr1D_unfoldBinCorr = (TH1D*) HSysBfCorr1D->Clone(Form("%s_unfoldBinCorr", HSysBfCorr1D->GetName()));
   TH1D* HSysAfCorr1D_unfoldBinCorr = (TH1D*) HSysAfCorr1D->Clone(Form("%s_unfoldBinCorr", HSysAfCorr1D->GetName()));
   TH1D* HSysBfCorr1DPlus_unfoldBinCorr = (TH1D*) HSysBfCorr1DPlus->Clone(Form("%s_unfoldBinCorr", HSysBfCorr1DPlus->GetName()));
   TH1D* HSysAfCorr1DPlus_unfoldBinCorr = (TH1D*) HSysAfCorr1DPlus->Clone(Form("%s_unfoldBinCorr", HSysAfCorr1DPlus->GetName()));


   TH1D* HSysBfCorr1D_unfoldBinCorrTheta = (TH1D*) HSysBfCorr1DTheta->Clone(Form("%s_unfoldBinCorrTheta", HSysBfCorr1DTheta->GetName()));
   TH1D* HSysAfCorr1D_unfoldBinCorrTheta = (TH1D*) HSysAfCorr1DTheta->Clone(Form("%s_unfoldBinCorrTheta", HSysAfCorr1DTheta->GetName()));
   TH1D* HSysBfCorr1DPlus_unfoldBinCorrTheta = (TH1D*) HSysBfCorr1DPlusTheta->Clone(Form("%s_unfoldBinCorrTheta", HSysBfCorr1DPlusTheta->GetName()));
   TH1D* HSysAfCorr1DPlus_unfoldBinCorrTheta = (TH1D*) HSysAfCorr1DPlusTheta->Clone(Form("%s_unfoldBinCorrTheta", HSysAfCorr1DPlusTheta->GetName()));


   UnfoldingBinCorrFactor.applyEffCorrOnHisto(HSysBfCorr1D, HSysBfCorr1D_unfoldBinCorr);
   UnfoldingBinCorrFactor.applyEffCorrOnHisto(HSysAfCorr1D, HSysAfCorr1D_unfoldBinCorr);
   UnfoldingBinCorrFactor.applyEffCorrOnHisto(HSysBfCorr1DPlus, HSysBfCorr1DPlus_unfoldBinCorr);
   UnfoldingBinCorrFactor.applyEffCorrOnHisto(HSysAfCorr1DPlus, HSysAfCorr1DPlus_unfoldBinCorr);

   UnfoldingBinCorrFactor_Theta.applyEffCorrOnHisto(HSysBfCorr1DTheta, HSysBfCorr1D_unfoldBinCorrTheta);
   UnfoldingBinCorrFactor_Theta.applyEffCorrOnHisto(HSysAfCorr1DTheta, HSysAfCorr1D_unfoldBinCorrTheta);
   UnfoldingBinCorrFactor_Theta.applyEffCorrOnHisto(HSysBfCorr1DPlusTheta, HSysBfCorr1DPlus_unfoldBinCorrTheta);
   UnfoldingBinCorrFactor_Theta.applyEffCorrOnHisto(HSysAfCorr1DPlusTheta, HSysAfCorr1DPlus_unfoldBinCorrTheta);


   HSysAfCorr1D_unfoldBinCorr->Scale(1.0/nEv);
   DivideByBin(*HSysAfCorr1D_unfoldBinCorr, zBins);
   HSysAfCorr1DPlus_unfoldBinCorr->Scale(1.0/nEv);
   DivideByBin(*HSysAfCorr1DPlus_unfoldBinCorr, zBins);

   HSysAfCorr1D_unfoldBinCorrTheta->Scale(1.0/nEv);
   DivideByBin(*HSysAfCorr1D_unfoldBinCorrTheta, Bins);
   HSysAfCorr1DPlus_unfoldBinCorrTheta->Scale(1.0/nEv);
   DivideByBin(*HSysAfCorr1DPlus_unfoldBinCorrTheta, Bins);


   TGraphErrors GzDataSyst( HDataAfCorr1D_unfoldBinCorr);
   GzDataSyst.SetName("HzDataSyst");




   for(int i = 1; i <= HDataAfCorr1D_unfoldBinCorr->GetNbinsX(); i++)
   {
      int iGraph = i-1;
      double err = HSysAfCorr1DPlus_unfoldBinCorr->GetBinContent(i) - HDataAfCorr1D_unfoldBinCorr->GetBinContent(i);
      GzDataSyst.SetPoint(iGraph, HDataAfCorr1D_unfoldBinCorr->GetBinCenter(i), HDataAfCorr1D_unfoldBinCorr->GetBinContent(i));
      GzDataSyst.SetPointError(iGraph, HDataAfCorr1D_unfoldBinCorr->GetBinWidth(i)/2, err);
      printf("index: %d { %.7e , %.7f , sys: %.7f stat: %.7f }\n", iGraph, zBins[iGraph], HDataAfCorr1D_unfoldBinCorr->GetBinContent(i), err,  HDataAfCorr1D_unfoldBinCorr->GetBinError(i));
   }
   GzDataSyst.SetLineWidth(0);
   GzDataSyst.SetFillStyle(1001);
   GzDataSyst.SetFillColorAlpha(kBlack, 0.3);

   std::cout << "HERWIG integral: " << HzMCGenHERWIG.Integral() << std::endl;

   std::vector<TH1D> hists = {*HDataAfCorr1D_unfoldBinCorr, HzMCGenBeforeRef, HzMCGenPYTHIA8, HzMCGenPYTHIA8Vincia, HzMCGenPYTHIA8Dire, HzMCGenHERWIG, HzMCGenSHERPA};
   std::vector<TH1D> histsPYTHIA8 = {*HDataAfCorr1D_unfoldBinCorr, HzMCGenPYTHIA8, HzMCGenPYTHIA8Vincia, HzMCGenPYTHIA8Dire};

   std::vector<TH1D> histsTheory = {*HDataAfCorr1D_unfoldBinCorr};
   


   system("mkdir -p plot/");
   MakeCanvasZ(histsPYTHIA8,GzDataSyst, {"Fully Corrected Data","PYTHIA8", "PYTHIA8 Vincia", "PYTHIA8 Dire"}, Form("plot/CorrectedData_PYTHIA8_z_september15"), "#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",1e-2,1e3, true, true);
   MakeCanvasZ(hists,GzDataSyst, {"Fully Corrected Data", "Archived MC","PYTHIA8", "PYTHIA8 Vincia", "PYTHIA8 Dire","HERWIG", "SHERPA"}, Form("plot/CorrectedData_z_september15"), "#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",1e-2,1e3, true, true);
   MakeCanvasZ(histsTheory,GzDataSyst, {"Fully Corrected Data"}, Form("plot/CorrectedDataONLY_z_september15"), "#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",1e-2,1e3, false, true);
   MakeCanvasZ({*HDataAfCorr1D_unfoldBinCorr,HzMCGenBeforeRef},GzDataSyst, {"Archived MC","Fully Corrected Data"}, Form("plot/CorrectedDataMC_z_september15"), "#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",1e-2,1e3, true, true);

   TGraphAsymmErrors theory = getTheoryPlotUpdated();
   //MakeCanvasZTheory(histsTheory,GzDataSyst, theory, {"Fully Corrected Data"}, Form("plot/CorrectedData_z_theory_May28th"), "#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",1e-2,1e3, false, true);

   TGraphErrors GthetaDataSyst( HDataAfCorr1D_unfoldBinCorr_Theta );
   GthetaDataSyst.SetName("GthetaDataSyst");

   for(int i = 1; i <= HDataAfCorr1D_unfoldBinCorr_Theta->GetNbinsX(); i++)
   {
      int iGraph = i-1;
      double err = HSysAfCorr1DPlus_unfoldBinCorrTheta->GetBinContent(i) - HDataAfCorr1D_unfoldBinCorr_Theta->GetBinContent(i);
      GthetaDataSyst.SetPoint(iGraph,HDataAfCorr1D_unfoldBinCorr_Theta->GetBinCenter(i), HDataAfCorr1D_unfoldBinCorr_Theta->GetBinContent(i));
      GthetaDataSyst.SetPointError(iGraph, HDataAfCorr1D_unfoldBinCorr_Theta->GetBinWidth(i)/2, err);
   }
   GthetaDataSyst.SetLineWidth(0);
   GthetaDataSyst.SetFillStyle(1001);
   GthetaDataSyst.SetFillColorAlpha(kBlack, 0.3);

   std::vector<TH1D> histsTheta = {*HDataAfCorr1D_unfoldBinCorr_Theta, HthetaMCGenBeforeRef, HthetaMCGenPYTHIA8, HthetaMCGenPYTHIA8Vincia, HthetaMCGenPYTHIA8Dire, HthetaMCGenHERWIG, HthetaMCGenSHERPA};
   std::vector<TH1D> histsThetaPYTHIA8 = {*HDataAfCorr1D_unfoldBinCorr_Theta, HthetaMCGenPYTHIA8, HthetaMCGenPYTHIA8Vincia, HthetaMCGenPYTHIA8Dire};
   MakeCanvas(histsThetaPYTHIA8,GthetaDataSyst, {"Fully Corrected Data", "PYTHIA8", "PYTHIA8 Vincia", "PYTHIA8 Dire"}, Form("plot/CorrectedData_PYTHIA8_theta_september15"), "#theta_{L}", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#theta_{L}}",1e-2,3, true, true);

   MakeCanvas(histsTheta,GthetaDataSyst, {"Fully Corrected Data", "Archived MC","PYTHIA8", "PYTHIA8 Vincia", "PYTHIA8 Dire", "HERWIG", "SHERPA"}, Form("plot/CorrectedData_theta_september15"), "#theta_{L}", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#theta_{L}}",1e-2,3, true, true);
   MakeCanvas({*HDataAfCorr1D_unfoldBinCorr_Theta},GthetaDataSyst, {"Fully Corrected Data"}, Form("plot/CorrectedDataONLY_theta_september15"), "#theta_{L}", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#theta_{L}}",1e-2,0.5, false, true);
   MakeCanvas({*HDataAfCorr1D_unfoldBinCorr_Theta, HthetaMCGenBeforeRef}, GthetaDataSyst, {"Fully Corrected Data","Archived MC"}, Form("plot/CorrectedDataMC_theta_September15th"), "#theta_{L}", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#theta_{L}}",1e-3,5, true, true);
   
   //MakeCanvasTheory({*HDataAfCorr1D_unfoldBinCorr_Theta},GthetaDataSyst, {"Fully Corrected Data"}, Form("plot/CorrectedData_theta_Theory_May28th"), "#theta_{L}", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#theta_{L}}",3e-3,0.5, false, true);

   Output.Close();

   if (applyEffCorrOnHistoErrorStatus>0)
   {
      printf("[Error] Something wrong with applyEffCorrOnHisto.\n");
      return 1;
   } else return 0;
}

int FindBin(double Value, int NBins, double Bins[])
{
   for(int i = 0; i < NBins; i++)
      if(Value < Bins[i])
         return i - 1;
   return NBins;
}

/*
* Make the canvas for the results as a function of z. Note that these results include a cutoff at 1e-6.
* Additionally, this is intended for MC comparisons.
*/
void MakeCanvasZ(vector<TH1D > Histograms, TGraphErrors DataSyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0].GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 17 : 0;
   double WorldXMax = LogX ? 183: 1;

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

   const int BinCount = 100;


   // z binning
   double zBins[2*BinCount+1];
   double zBinMin = (1- cos(0.002))/2;
   double zBinMax = 0.5;

   for(int i = 0; i <= BinCount; i++){
      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
   }

   TH2D HWorld("HWorld", "", N, WorldXMin, WorldXMax, 100, WorldMin, WorldMax);
   HWorld.SetStats(0);
   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);

   HWorld.Draw("axis");
   for(TH1D H : Histograms){
      TH1D *HClone = (TH1D *)H.Clone();
      HClone->Draw("exp same");
   }
   DataSyst.DrawClone("2 same");


   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : 1 / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : 1/ 2, 10000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   if(DoRatio)
      PadR.cd();

   double WorldRMin = 0.59999999;
   double WorldRMax = 1.39;//99999;

   TH2D HWorldR("HWorldR", "", N, WorldXMin, WorldXMax, 100, WorldRMin, WorldRMax);
   TGraph G2;

   if(DoRatio)
   {
      HWorldR.SetStats(0);
      HWorldR.GetXaxis()->SetTickLength(0);
      HWorldR.GetXaxis()->SetLabelSize(0);
      HWorldR.GetYaxis()->SetNdivisions(505);

      HWorldR.Draw("axis");
      for(int i = 1; i <= Histograms[0].GetNbinsX(); i++)
      {
         int iGraph = i-1;
         DataSyst.SetPoint(iGraph,
                                 DataSyst.GetPointX(iGraph),
                                 DataSyst.GetPointY(iGraph)/Histograms[0].GetBinContent(i));
         DataSyst.SetPointError( iGraph,
                                 DataSyst.GetErrorX(iGraph),
                                 DataSyst.GetErrorY(iGraph)/Histograms[0].GetBinContent(i));
      }
      DataSyst.DrawClone("2 same");
      
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

   double BinMin    = zBins[17];//(1- cos(0.002))/2;
   double BinMiddle = 0.5;
   double BinMax    = 1 - BinMin;

   Canvas.cd();
   int nDiv = 505;
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X5(MarginL, MarginB + PadHeight + PadRHeight, MarginL + PadWidth / 2, MarginB + PadHeight + PadRHeight, BinMin, BinMiddle, 510, "-GS"); // - in the draw options means we only draw axis on the "negative" side
   // axis on the x axis on the right hand side for the top of the plot
   TGaxis X6(MarginL + PadWidth, MarginB + PadHeight + PadRHeight, MarginL + PadWidth / 2, MarginB + PadHeight + PadRHeight, BinMin, BinMiddle, 510, "+GS");// - in the draw options means we only draw axis on the "negative" side
   
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");

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
   // if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.02, MarginB - 0.01, "10^{-6} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.1, MarginB - 0.01, "10^{-4} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.32, MarginB - 0.01, "10^{-2} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.500, MarginB - 0.01, "1/2");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.68, MarginB - 0.01, "1 - 10^{-2}");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.9, MarginB - 0.01, "1 - 10^{-4}");
   // if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.995, MarginB - 0.01, "1 - 10^{-6}");



   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kGray);
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.015, "#it{z} = 1/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.9, MarginB * 0.4, X.c_str());

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(90);
   Latex.SetTextColor(kBlack);
   if(DoRatio)
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "MC/Data");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "2025 HB - Finalization of Result");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.04 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
   {
      if (Labels[i]=="Data")
      {
         Histograms[i].SetFillStyle(DataSyst.GetFillStyle());
         Histograms[i].SetFillColor(DataSyst.GetFillColor());
      }
      Legend.AddEntry(&Histograms[i], Labels[i].c_str(),
                      (Labels[i]=="Data")? "plf": "pl");
   }
   Legend.Draw();

   TLegend Legend2(0.55, 0.90, 0.8, 0.90 - 0.04 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(&Histograms[i], Labels[i].c_str(),
                      (Labels[i]=="Data")? "plf": "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}

// make the theory comparison canvas
void MakeCanvasZTheory(vector<TH1D > Histograms, TGraphErrors DataSyst,  TGraphAsymmErrors TheorySyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0].GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 17 : 0;
   double WorldXMax = LogX ? 183: 1;

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

   const int BinCount = 100;


   // z binning
   double zBins[2*BinCount+1];
   double zBinMin = (1- cos(0.002))/2;
   double zBinMax = 0.5;

   for(int i = 0; i <= BinCount; i++){
      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
   }

   TH2D HWorld("HWorld", "", N, WorldXMin, WorldXMax, 100, WorldMin, WorldMax);
   HWorld.SetStats(0);
   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);

   HWorld.Draw("axis");
   TheorySyst.DrawClone("3 l same");

   for(TH1D H : Histograms){
      TH1D *HClone = (TH1D *)H.Clone();
      HClone->Draw("exp same");
   }
   DataSyst.DrawClone("2 same");



   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : 1 / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : 1/ 2, 10000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   if(DoRatio)
      PadR.cd();

   double WorldRMin = 0.0;
   double WorldRMax = 1.99999;

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

      for(int i = 1; i <= Histograms[0].GetNbinsX(); i++)
      {
         int iGraph = i-1;
         TheorySyst.SetPoint(iGraph,
                                 TheorySyst.GetPointX(iGraph),
                                 TheorySyst.GetPointY(iGraph)/Histograms[0].GetBinContent(i));
         TheorySyst.SetPointError( iGraph,
                                 TheorySyst.GetErrorXlow(iGraph), TheorySyst.GetErrorXhigh(iGraph),
                                 TheorySyst.GetErrorYlow(iGraph)/Histograms[0].GetBinContent(i), TheorySyst.GetErrorYhigh(iGraph)/Histograms[0].GetBinContent(i));
      }
      TheorySyst.DrawClone("2 same");

      G.Draw("l");

      G2.SetPoint(0, 0, 1);
      G2.SetPoint(1, 99999, 1);
      G2.Draw("l");
   }

   double BinMin    = zBins[17];//(1- cos(0.002))/2;
   double BinMiddle = 0.5;
   double BinMax    = 1 - BinMin;

   Canvas.cd();
   int nDiv = 505;
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X5(MarginL, MarginB + PadHeight + PadRHeight, MarginL + PadWidth / 2, MarginB + PadHeight + PadRHeight, BinMin, BinMiddle, 510, "-GS"); // - in the draw options means we only draw axis on the "negative" side
   // axis on the x axis on the right hand side for the top of the plot
   TGaxis X6(MarginL + PadWidth, MarginB + PadHeight + PadRHeight, MarginL + PadWidth / 2, MarginB + PadHeight + PadRHeight, BinMin, BinMiddle, 510, "+GS");// - in the draw options means we only draw axis on the "negative" side
   
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");

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
   // if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.02, MarginB - 0.01, "10^{-6} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.1, MarginB - 0.01, "10^{-4} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.32, MarginB - 0.01, "10^{-2} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.500, MarginB - 0.01, "1/2");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.68, MarginB - 0.01, "1 - 10^{-2}");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.9, MarginB - 0.01, "1 - 10^{-4}");
   // if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.995, MarginB - 0.01, "1 - 10^{-6}");

   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kGray);
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.015, "#it{z} = 1/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.9, MarginB * 0.4, X.c_str());

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(90);
   Latex.SetTextColor(kBlack);
   if(DoRatio)
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "MC/Data");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "2025 HB - Finalization of Result");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.04 * min(3, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.03);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
   {
      if (Labels[i]=="Data")
      {
         Histograms[i].SetFillStyle(DataSyst.GetFillStyle());
         Histograms[i].SetFillColor(DataSyst.GetFillColor());
      }
      Legend.AddEntry(&Histograms[i], Labels[i].c_str(),
                      (Labels[i]=="Data")? "plf": "pl");
   }
   Legend.AddEntry(&TheorySyst, "Track Function Theory Calculation", "pl");
   Legend.AddEntry((TObject*)0, "(NNLL Collinear + NNNLL Sudakov)", "");
   Legend.Draw();

   TLegend Legend2(0.55, 0.90, 0.8, 0.90 - 0.04 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.03);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(&Histograms[i], Labels[i].c_str(),
                      (Labels[i]=="Data")? "plf": "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}

// create the settings for the tpad
void SetPad(TPad &P){
   P.SetLeftMargin(0);
   P.SetTopMargin(0);
   P.SetRightMargin(0);
   P.SetBottomMargin(0);
   P.SetTickx();
   P.SetTicky();
   P.Draw();
}


// normalize by the bin width
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



// takes in theory points from Ian and returns TGraphAsymmErrors of the points
TGraphAsymmErrors getTheoryPlotUpdated(){
      // points for the central values of the theory curves
      double centralvals[100][2] = {{1.5708, 1.69714}, {1.69409, 1.75999}, {1.80373, 1.8824}, {1.90226, 
                                    1.99367}, {1.99148, 2.15729}, {2.07275, 2.31996}, {2.14712, 
                                    2.54475}, {2.21543, 2.78266}, {2.27836, 3.03587}, {2.33648, 
                                    3.35711}, {2.39026, 3.7039}, {2.44012, 4.13161}, {2.48639, 
                                    4.52754}, {2.5294, 5.07051}, {2.56942, 5.56135}, {2.60668, 
                                    6.11072}, {2.6414, 6.74856}, {2.67378, 7.38398}, {2.70399, 
                                    8.06813}, {2.73219, 8.86919}, {2.75852, 9.72367}, {2.78313, 
                                    10.529}, {2.80611, 11.3214}, {2.8276, 12.364}, {2.84769, 
                                    13.4936}, {2.86648, 14.5506}, {2.88405, 15.5715}, {2.90048, 
                                    16.5704}, {2.91587, 17.7279}, {2.93026, 18.9662}, {2.94373, 
                                    20.214}, {2.95633, 21.4576}, {2.96812, 22.7318}, {2.97917, 
                                    24.07}, {2.9895, 25.4153}, {2.99918, 26.7228}, {3.00823, 
                                    27.9559}, {3.01671, 29.1038}, {3.02465, 30.1469}, {3.03208, 
                                    31.0205}, {3.03904, 31.6504}, {3.04556, 32.1013}, {3.05167, 
                                    32.3991}, {3.05738, 32.5427}, {3.06273, 32.5418}, {3.06774, 
                                    32.4245}, {3.07244, 32.1746}, {3.07684, 31.7905}, {3.08095, 
                                    31.272}, {3.0848, 30.6208}, {3.08841, 29.8399}, {3.09179, 
                                    28.9327}, {3.09494, 27.9158}, {3.09791, 26.7947}, {3.1007, 
                                    25.5876}, {3.10328, 24.3392}, {3.1057, 23.0468}, {3.10801, 
                                    21.719}, {3.11016, 20.3944}, {3.11213, 19.118}, {3.11402, 
                                    17.84}, {3.11575, 16.6394}, {3.11743, 15.4388}, {3.11896, 
                                    14.317}, {3.12033, 13.3054}, {3.12169, 12.287}, {3.12294, 
                                    11.346}, {3.12416, 10.4166}, {3.12522, 9.59842}, {3.12636, 
                                    8.71581}, {3.12731, 7.97389}, {3.12818, 7.28946}, {3.1291, 
                                    6.5485}, {3.12976, 6.01506}, {3.13064, 5.29115}, {3.1312, 
                                    4.81831}, {3.132, 4.12892}, {3.13265, 3.55408}, {3.13311, 
                                    3.13465}, {3.13359, 2.67804}, {3.13411, 2.17383}, {3.13466, 
                                    1.60629}, {3.13496, 1.29143}, {3.13559, 0.576117}, {3.13594, 
                                    0.153484}, {3.1363, -0.00265485}, {3.13669, 
                                    1.40345e-7}, {3.13669, 1.40345e-7}, {3.13712, 
                                    1.28117e-7}, {3.13759, 1.14591e-7}, {3.13759, 
                                    1.14591e-7}, {3.13813, 9.92391e-8}, {3.13813, 
                                    9.92391e-8}, {3.13813, 9.92391e-8}, {3.13876, 
                                    8.10284e-8}, {3.13876, 8.10284e-8}, {3.13876, 
                                    8.10284e-8}, {3.13959, 5.72958e-8}, {3.13959, 
                                    5.72958e-8}, {3.13959, 5.72958e-8}};
                                    
         double  err_low[100][2] ={{1.5708, 1.62255}, {1.69409, 1.68275}, {1.80373, 1.79461}, {1.90226, 
                                    1.90247}, {1.99148, 2.05483}, {2.07275, 2.21023}, {2.14712, 
                                    2.42057}, {2.21543, 2.64203}, {2.27836, 2.8943}, {2.33648, 
                                    3.20007}, {2.39026, 3.53116}, {2.44012, 3.93114}, {2.48639, 
                                    4.31198}, {2.5294, 4.82583}, {2.56942, 5.31425}, {2.60668, 
                                    5.85832}, {2.6414, 6.47311}, {2.67378, 7.08382}, {2.70399, 
                                    7.7374}, {2.73219, 8.49264}, {2.75852, 9.29809}, {2.78313, 
                                    10.0696}, {2.80611, 10.8358}, {2.8276, 11.8155}, {2.84769, 
                                    12.8731}, {2.86648, 13.8767}, {2.88405, 14.8568}, {2.90048, 
                                    15.8248}, {2.91587, 16.9312}, {2.93026, 18.1111}, {2.94373, 
                                    19.3048}, {2.95633, 20.501}, {2.96812, 21.7292}, {2.97917, 
                                    23.0184}, {2.9895, 24.3211}, {2.99918, 25.6005}, {3.00823, 
                                    26.8215}, {3.01671, 27.5946}, {3.02465, 28.1413}, {3.03208, 
                                    28.4435}, {3.03904, 28.4616}, {3.04556, 28.2845}, {3.05167, 
                                    27.9559}, {3.05738, 27.5668}, {3.06273, 27.1859}, {3.06774, 
                                    26.886}, {3.07244, 26.6662}, {3.07684, 26.5113}, {3.08095, 
                                    26.3792}, {3.0848, 26.2028}, {3.08841, 25.8983}, {3.09179, 
                                    25.0487}, {3.09494, 24.0527}, {3.09791, 22.9775}, {3.1007, 
                                    21.8369}, {3.10328, 20.6725}, {3.1057, 19.4804}, {3.10801, 
                                    18.2678}, {3.11016, 17.0683}, {3.11213, 15.9203}, {3.11402, 
                                    14.7771}, {3.11575, 13.7077}, {3.11743, 12.6418}, {3.11896, 
                                    11.6457}, {3.12033, 10.7484}, {3.12169, 9.84665}, {3.12294, 
                                    9.01449}, {3.12416, 8.19338}, {3.12522, 7.46964}, {3.12636, 
                                    6.68553}, {3.12731, 6.02379}, {3.12818, 5.4112}, {3.1291, 
                                    4.74534}, {3.12976, 4.26342}, {3.13064, 3.60455}, {3.1312, 
                                    3.17013}, {3.132, 2.5291}, {3.13265, 1.98579}, {3.13311, 
                                    1.5832}, {3.13359, 1.13802}, {3.13411, 0.634698}, {3.13466, 
                                    0.113149}, {3.13496, -0.00825052}, {3.13559, 
                                    1.71887e-7}, {3.13594, 1.62056e-7}, {3.1363, 
                                    1.5159e-7}, {3.13669, 1.40345e-7}, {3.13669, 
                                    1.40345e-7}, {3.13712, 1.28117e-7}, {3.13759, 
                                    1.14591e-7}, {3.13759, 1.14591e-7}, {3.13813, 
                                    9.92391e-8}, {3.13813, 9.92391e-8}, {3.13813, 
                                    9.92391e-8}, {3.13876, 8.10284e-8}, {3.13876, 
                                    8.10284e-8}, {3.13876, 8.10284e-8}, {3.13959, 
                                    5.72958e-8}, {3.13959, 5.72958e-8}, {3.13959, 5.72958e-8}}; 
                                    
           double  err_high[100][2] ={{1.5708, 1.76402}, {1.69409, 1.82836}, {1.80373, 1.96297}, {1.90226, 
                                 2.07583}, {1.99148, 2.2531}, {2.07275, 2.42262}, {2.14712, 
                                 2.66296}, {2.21543, 2.90983}, {2.27836, 3.16607}, {2.33648, 
                                 3.49746}, {2.39026, 3.84577}, {2.44012, 4.28083}, {2.48639, 
                                 4.66392}, {2.5294, 5.21625}, {2.56942, 5.71031}, {2.60668, 
                                 6.27395}, {2.6414, 6.93337}, {2.67378, 7.58344}, {2.70399, 
                                 8.28404}, {2.73219, 9.11351}, {2.75852, 9.99839}, {2.78313, 
                                 10.819}, {2.80611, 11.6218}, {2.8276, 12.7042}, {2.84769, 
                                 13.8837}, {2.86648, 14.9731}, {2.88405, 16.0163}, {2.90048, 
                                 17.0355}, {2.91587, 18.2354}, {2.93026, 19.5251}, {2.94373, 
                                 20.8266}, {2.95633, 22.105}, {2.96812, 23.3952}, {2.97917, 
                                 24.738}, {2.9895, 26.0575}, {2.99918, 27.2981}, {3.00823, 
                                 28.4208}, {3.01671, 29.4543}, {3.02465, 30.4157}, {3.03208, 
                                 31.2937}, {3.03904, 32.035}, {3.04556, 32.6703}, {3.05167, 
                                 33.1974}, {3.05738, 33.595}, {3.06273, 33.8726}, {3.06774, 
                                 34.1495}, {3.07244, 34.3953}, {3.07684, 34.5571}, {3.08095, 
                                 34.6216}, {3.0848, 34.5826}, {3.08841, 34.4358}, {3.09179, 
                                 34.1806}, {3.09494, 33.8183}, {3.09791, 33.3714}, {3.1007, 
                                 32.8567}, {3.10328, 32.2495}, {3.1057, 31.5451}, {3.10801, 
                                 30.7404}, {3.11016, 29.8542}, {3.11213, 28.9182}, {3.11402, 
                                 27.897}, {3.11575, 26.8572}, {3.11743, 25.7359}, {3.11896, 
                                 24.6122}, {3.12033, 23.5326}, {3.12169, 22.3794}, {3.12294, 
                                 21.2537}, {3.12416, 20.0848}, {3.12522, 19.0087}, {3.12636, 
                                 17.7995}, {3.12731, 16.746}, {3.12818, 15.7462}, {3.1291, 
                                 14.6367}, {3.12976, 13.8226}, {3.13064, 12.701}, {3.1312, 
                                 11.9603}, {3.132, 10.8734}, {3.13265, 9.96475}, {3.13311, 
                                 9.30289}, {3.13359, 8.58561}, {3.13411, 7.80001}, {3.13466, 
                                 6.9274}, {3.13496, 6.45003}, {3.13559, 5.38735}, {3.13594, 
                                 4.78628}, {3.1363, 4.12287}, {3.13669, 3.37815}, {3.13669, 
                                 3.37815}, {3.13712, 2.52213}, {3.13759, 1.50609}, {3.13759, 
                                 1.50609}, {3.13813, 0.246218}, {3.13813, 0.246218}, {3.13813, 
                                 0.246218}, {3.13876, -1.42451}, {3.13876, -1.42451}, {3.13876, -1.42451}, {3.13959, -3.947}, {3.13959, -3.947}, {3.13959, -3.947}};
  

   const int BinCount = 100;
   // z binning
   double zBins[2*BinCount+1];
   double zBinMin = (1- cos(0.002))/2;
   double zBinMax = 0.5;

   const int thetaBinCount = 100;
   double thetaBins[2*thetaBinCount+1];
   double thetaBinMin = 0.002;
   double thetaBinMax = M_PI / 2;

   for(int i = 0; i <= BinCount; i++){
      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
       // theta double log binning
      thetaBins[i] = exp(log(thetaBinMin) + (log(thetaBinMax) - log(thetaBinMin)) / thetaBinCount * i);
      thetaBins[2*thetaBinCount-i] = thetaBinMax * 2 - exp(log(thetaBinMin) + (log(thetaBinMax) - log(thetaBinMin)) / thetaBinCount * i);
   }
                                    
   TH1D HTemp("HTemp", "HTemp", 100, 0, 100);

   TGraphAsymmErrors graph(&HTemp);



   for(int i = 1; i <= 100; i++)
   {
      int iGraph = i-1;
      double centralValueRadians = centralvals[iGraph][0];
      double centralValueZ = (1- cos(centralValueRadians))/2; 
      
      double centralValueRadiansUp = centralvals[iGraph+1][0];
      double centralValueZUp = (1- cos(centralValueRadiansUp))/2; 
      
      double centralValueRadiansDown = centralvals[iGraph-1][0];
      double centralValueZDown = (1- cos(centralValueRadiansDown))/2; 
      
      int bin     = FindBin(centralValueRadians,  2*thetaBinCount, thetaBins);
      int binUp   = FindBin(centralValueRadiansUp, 2*thetaBinCount, thetaBins);
      int binDown = FindBin(centralValueRadiansDown,  2*thetaBinCount, thetaBins); 

      double binFraction = FindBinFraction(centralValueRadians,  2*thetaBinCount, thetaBins); 
      double binFractionUp = FindBinFraction(centralValueRadiansUp, 2*thetaBinCount, thetaBins); 
      double binFractionDown = FindBinFraction(centralValueRadiansDown,  2*thetaBinCount, thetaBins); 

      std::cout << "Setting point: " << iGraph << " with a theta value of " << centralValueRadians << " with a z value of " << centralValueZ << " with a bin of " << bin + binFraction << " with a y value of " << centralvals[iGraph][1]*(TMath::Pi()/180)*2 << std::endl;
      // set the point, not worrying about the errors. 
      graph.SetPoint(iGraph, HTemp.GetBinLowEdge(bin+1)+binFraction, centralvals[iGraph][1]*(TMath::Pi()/180)*0.5); 
      graph.SetPointError(iGraph, HTemp.GetBinWidth(i)/2, HTemp.GetBinWidth(i)/2, centralvals[iGraph][1]*(TMath::Pi()/180)*0.5 - err_low[iGraph][1]*(TMath::Pi()/180)*0.5, err_high[iGraph][1]*(TMath::Pi()/180)*0.5 - centralvals[iGraph][1]*(TMath::Pi()/180)*0.5);

   }

   static vector<int> Colors = GetCVDColors6();

   graph.SetFillStyle(1001);
   graph.SetMarkerStyle(20);
   graph.SetMarkerColor(kRed);
   graph.SetLineColor(kRed);
   graph.SetLineWidth(2);
   graph.SetFillColorAlpha(kRed, 0.3);


   return graph;

}

// takes in theory points from Ian and returns TGraphAsymmErrors of the points
TGraphAsymmErrors getTheoryPlotTheta(){
      // points for the central values of the theory curves
      double centralvals[100][2] =  {{1.5708, 1.69714}, {1.69409, 1.75999}, {1.80373, 1.8824}, {1.90226, 
  1.99367}, {1.99148, 2.15729}, {2.07275, 2.31996}, {2.14712, 
  2.54475}, {2.21543, 2.78266}, {2.27836, 3.03587}, {2.33648, 
  3.35711}, {2.39026, 3.7039}, {2.44012, 4.13161}, {2.48639, 
  4.52754}, {2.5294, 5.07051}, {2.56942, 5.56135}, {2.60668, 
  6.11072}, {2.6414, 6.74856}, {2.67378, 7.38398}, {2.70399, 
  8.06813}, {2.73219, 8.86919}, {2.75852, 9.72367}, {2.78313, 
  10.529}, {2.80611, 11.3214}, {2.8276, 12.364}, {2.84769, 
  13.4936}, {2.86648, 14.5506}, {2.88405, 15.5715}, {2.90048, 
  16.5704}, {2.91587, 17.7279}, {2.93026, 18.9662}, {2.94373, 
  20.214}, {2.95633, 21.4576}, {2.96812, 22.7318}, {2.97917, 
  24.07}, {2.9895, 25.4153}, {2.99918, 26.7228}, {3.00823, 
  27.9559}, {3.01671, 29.1038}, {3.02465, 30.1469}, {3.03208, 
  31.0205}, {3.03904, 31.6504}, {3.04556, 32.1013}, {3.05167, 
  32.3991}, {3.05738, 32.5427}, {3.06273, 32.5418}, {3.06774, 
  32.4245}, {3.07244, 32.1746}, {3.07684, 31.7905}, {3.08095, 
  31.272}, {3.0848, 30.6208}, {3.08841, 29.8399}, {3.09179, 
  28.9327}, {3.09494, 27.9158}, {3.09791, 26.7947}, {3.1007, 
  25.5876}, {3.10328, 24.3392}, {3.1057, 23.0468}, {3.10801, 
  21.719}, {3.11016, 20.3944}, {3.11213, 19.118}, {3.11402, 
  17.84}, {3.11575, 16.6394}, {3.11743, 15.4388}, {3.11896, 
  14.317}, {3.12033, 13.3054}, {3.12169, 12.287}, {3.12294, 
  11.346}, {3.12416, 10.4166}, {3.12522, 9.59842}, {3.12636, 
  8.71581}, {3.12731, 7.97389}, {3.12818, 7.28946}, {3.1291, 
  6.5485}, {3.12976, 6.01506}, {3.13064, 5.29115}, {3.1312, 
  4.81831}, {3.132, 4.12892}, {3.13265, 3.55408}, {3.13311, 
  3.13465}, {3.13359, 2.67804}, {3.13411, 2.17383}, {3.13466, 
  1.60629}, {3.13496, 1.29143}, {3.13559, 0.576117}, {3.13594, 
  0.153484}, {3.1363, -0.00265485}, {3.13669, 
  1.40345e-7}, {3.13669, 1.40345e-7}, {3.13712, 
  1.28117e-7}, {3.13759, 1.14591e-7}, {3.13759, 
  1.14591e-7}, {3.13813, 9.92391e-8}, {3.13813, 
  9.92391e-8}, {3.13813, 9.92391e-8}, {3.13876, 
  8.10284e-8}, {3.13876, 8.10284e-8}, {3.13876, 
  8.10284e-8}, {3.13959, 5.72958e-8}, {3.13959, 
  5.72958e-8}, {3.13959, 5.72958e-8}}; 

      // points for the higher values of the error bands
      double err_high[100][2] ={{1.5708, 1.76402}, {1.69409, 1.82836}, {1.80373, 1.96297}, {1.90226, 
  2.07583}, {1.99148, 2.2531}, {2.07275, 2.42262}, {2.14712, 
  2.66296}, {2.21543, 2.90983}, {2.27836, 3.16607}, {2.33648, 
  3.49746}, {2.39026, 3.84577}, {2.44012, 4.28083}, {2.48639, 
  4.66392}, {2.5294, 5.21625}, {2.56942, 5.71031}, {2.60668, 
  6.27395}, {2.6414, 6.93337}, {2.67378, 7.58344}, {2.70399, 
  8.28404}, {2.73219, 9.11351}, {2.75852, 9.99839}, {2.78313, 
  10.819}, {2.80611, 11.6218}, {2.8276, 12.7042}, {2.84769, 
  13.8837}, {2.86648, 14.9731}, {2.88405, 16.0163}, {2.90048, 
  17.0355}, {2.91587, 18.2354}, {2.93026, 19.5251}, {2.94373, 
  20.8266}, {2.95633, 22.105}, {2.96812, 23.3952}, {2.97917, 
  24.738}, {2.9895, 26.0575}, {2.99918, 27.2981}, {3.00823, 
  28.4208}, {3.01671, 29.4543}, {3.02465, 30.4157}, {3.03208, 
  31.2937}, {3.03904, 32.035}, {3.04556, 32.6703}, {3.05167, 
  33.1974}, {3.05738, 33.595}, {3.06273, 33.8726}, {3.06774, 
  34.1495}, {3.07244, 34.3953}, {3.07684, 34.5571}, {3.08095, 
  34.6216}, {3.0848, 34.5826}, {3.08841, 34.4358}, {3.09179, 
  34.1806}, {3.09494, 33.8183}, {3.09791, 33.3714}, {3.1007, 
  32.8567}, {3.10328, 32.2495}, {3.1057, 31.5451}, {3.10801, 
  30.7404}, {3.11016, 29.8542}, {3.11213, 28.9182}, {3.11402, 
  27.897}, {3.11575, 26.8572}, {3.11743, 25.7359}, {3.11896, 
  24.6122}, {3.12033, 23.5326}, {3.12169, 22.3794}, {3.12294, 
  21.2537}, {3.12416, 20.0848}, {3.12522, 19.0087}, {3.12636, 
  17.7995}, {3.12731, 16.746}, {3.12818, 15.7462}, {3.1291, 
  14.6367}, {3.12976, 13.8226}, {3.13064, 12.701}, {3.1312, 
  11.9603}, {3.132, 10.8734}, {3.13265, 9.96475}, {3.13311, 
  9.30289}, {3.13359, 8.58561}, {3.13411, 7.80001}, {3.13466, 
  6.9274}, {3.13496, 6.45003}, {3.13559, 5.38735}, {3.13594, 
  4.78628}, {3.1363, 4.12287}, {3.13669, 3.37815}, {3.13669, 
  3.37815}, {3.13712, 2.52213}, {3.13759, 1.50609}, {3.13759, 
  1.50609}, {3.13813, 0.246218}, {3.13813, 0.246218}, {3.13813, 
  0.246218}, {3.13876, -1.42451}, {3.13876, -1.42451}, {3.13876, -1.42451}, {3.13959, -3.947}, {3.13959, -3.947}, {3.13959, -3.947}};

         double err_low[100][2] = {{1.5708, 1.62255}, {1.69409, 1.68275}, {1.80373, 1.79461}, {1.90226, 
  1.90247}, {1.99148, 2.05483}, {2.07275, 2.21023}, {2.14712, 
  2.42057}, {2.21543, 2.64203}, {2.27836, 2.8943}, {2.33648, 
  3.20007}, {2.39026, 3.53116}, {2.44012, 3.93114}, {2.48639, 
  4.31198}, {2.5294, 4.82583}, {2.56942, 5.31425}, {2.60668, 
  5.85832}, {2.6414, 6.47311}, {2.67378, 7.08382}, {2.70399, 
  7.7374}, {2.73219, 8.49264}, {2.75852, 9.29809}, {2.78313, 
  10.0696}, {2.80611, 10.8358}, {2.8276, 11.8155}, {2.84769, 
  12.8731}, {2.86648, 13.8767}, {2.88405, 14.8568}, {2.90048, 
  15.8248}, {2.91587, 16.9312}, {2.93026, 18.1111}, {2.94373, 
  19.3048}, {2.95633, 20.501}, {2.96812, 21.7292}, {2.97917, 
  23.0184}, {2.9895, 24.3211}, {2.99918, 25.6005}, {3.00823, 
  26.8215}, {3.01671, 27.5946}, {3.02465, 28.1413}, {3.03208, 
  28.4435}, {3.03904, 28.4616}, {3.04556, 28.2845}, {3.05167, 
  27.9559}, {3.05738, 27.5668}, {3.06273, 27.1859}, {3.06774, 
  26.886}, {3.07244, 26.6662}, {3.07684, 26.5113}, {3.08095, 
  26.3792}, {3.0848, 26.2028}, {3.08841, 25.8983}, {3.09179, 
  25.0487}, {3.09494, 24.0527}, {3.09791, 22.9775}, {3.1007, 
  21.8369}, {3.10328, 20.6725}, {3.1057, 19.4804}, {3.10801, 
  18.2678}, {3.11016, 17.0683}, {3.11213, 15.9203}, {3.11402, 
  14.7771}, {3.11575, 13.7077}, {3.11743, 12.6418}, {3.11896, 
  11.6457}, {3.12033, 10.7484}, {3.12169, 9.84665}, {3.12294, 
  9.01449}, {3.12416, 8.19338}, {3.12522, 7.46964}, {3.12636, 
  6.68553}, {3.12731, 6.02379}, {3.12818, 5.4112}, {3.1291, 
  4.74534}, {3.12976, 4.26342}, {3.13064, 3.60455}, {3.1312, 
  3.17013}, {3.132, 2.5291}, {3.13265, 1.98579}, {3.13311, 
  1.5832}, {3.13359, 1.13802}, {3.13411, 0.634698}, {3.13466, 
  0.113149}, {3.13496, -0.00825052}, {3.13559, 
  1.71887e-7}, {3.13594, 1.62056e-7}, {3.1363, 
  1.5159e-7}, {3.13669, 1.40345e-7}, {3.13669, 
  1.40345e-7}, {3.13712, 1.28117e-7}, {3.13759, 
  1.14591e-7}, {3.13759, 1.14591e-7}, {3.13813, 
  9.92391e-8}, {3.13813, 9.92391e-8}, {3.13813, 
  9.92391e-8}, {3.13876, 8.10284e-8}, {3.13876, 
  8.10284e-8}, {3.13876, 8.10284e-8}, {3.13959, 
  5.72958e-8}, {3.13959, 5.72958e-8}, {3.13959, 5.72958e-8}};



   TH1D HTemp("HTemp", "HTemp", 100, 100, 200);

   TGraphAsymmErrors graph(&HTemp);


   for(int i = 1; i <= HTemp.GetNbinsX(); i++)
   {
      int iGraph = i-1;
      graph.SetPoint(iGraph, HTemp.GetBinCenter(i), centralvals[iGraph][1]);
      graph.SetPointError(iGraph, HTemp.GetBinWidth(i)/2, HTemp.GetBinWidth(i)/2, centralvals[iGraph][1] - err_low[iGraph][1], err_high[iGraph][1] - centralvals[iGraph][1]);
   }

   static vector<int> Colors = GetCVDColors6();

   graph.SetFillStyle(1001);
   graph.SetMarkerStyle(20);
   graph.SetMarkerColor(kRed);
   graph.SetLineColor(kRed);
   graph.SetLineWidth(2);
   graph.SetFillColorAlpha(kRed, 0.3);


   return graph;

}


// takes in theory points from Ian and returns TGraphAsymmErrors of the points
TGraphAsymmErrors getTheoryPlot(){
      // points for the central values of the theory curves
      double centralvals[200][2] =  {{1.14022*10e-6, 10.6999}, {1.30011*10e-6, 10.6999}, {1.48241*10e-6,
      10.6999}, {1.69028*10e-6, 10.6998}, {1.9273*10e-6,
      10.6998}, {2.19755*10e-6, 10.6997}, {2.50569*10e-6,
      10.6997}, {2.85705*10e-6, 10.6996}, {3.25767*10e-6,
      10.6995}, {3.71447*10e-6, 10.6994}, {4.23532*10e-6,
      10.6993}, {4.82921*10e-6, 10.6992}, {5.50638*10e-6,
      10.699}, {6.2785*10e-6, 10.6988}, {7.15888*10e-6,
      10.6985}, {8.16272*10e-6, 10.6982}, {9.30732*10e-6,
      10.6978}, {0.0000106124, 10.6973}, {0.0000121005,
      10.6967}, {0.0000137973, 10.6959}, {0.000015732,
      10.6951}, {0.000017938, 10.694}, {0.0000204533,
      10.6927}, {0.0000233213, 10.6911}, {0.0000265915,
      10.6892}, {0.0000303202, 10.6868}, {0.0000345718,
      10.684}, {0.0000394195, 10.6805}, {0.0000449471,
      10.6762}, {0.0000512496, 10.6711}, {0.000058436,
      10.6648}, {0.0000666301, 10.6572}, {0.0000759731,
      10.6479}, {0.0000866263, 10.6366}, {0.0000987733,
      10.6229}, {0.000112624, 10.6063}, {0.000128416,
      10.5862}, {0.000146423, 10.5617}, {0.000166955,
      10.5321}, {0.000190365, 10.4963}, {0.000217059,
      10.453}, {0.000247496, 10.4007}, {0.0002822, 10.3378}, {0.000321771,
         10.2623}, {0.00036689, 10.1717}, {0.000418337,
      10.0637}, {0.000476997, 9.93512}, {0.000543883,
      9.78298}, {0.000620148, 9.60392}, {0.000707107,
      9.39457}, {0.000806259, 9.15168}, {0.000919315,
      8.87238}, {0.00104822, 8.55452}, {0.00119521, 8.19696}, {0.0013628,
      7.80003}, {0.0015539, 7.36573}, {0.00177179, 6.89812}, {0.00202024,
      6.40316}, {0.00230352, 5.88873}, {0.00262653, 5.36402}, {0.00299483,
         4.83905}, {0.00341477, 4.32383}, {0.0038936, 3.82763}, {0.00443957,
         3.35839}, {0.0050621, 2.92221}, {0.00577192, 2.52321}, {0.00658127,
         2.16354}, {0.00750412, 1.84357}, {0.00855636,
      1.56227}, {0.00975616, 1.3175}, {0.0111242, 1.10644}, {0.0126841,
      0.925849}, {0.0144627, 0.772364}, {0.0164906, 0.641358}, {0.018803,
      0.531175}, {0.0214396, 0.448236}, {0.0244459, 0.386252}, {0.0278738,
         0.331888}, {0.0317824, 0.28539}, {0.036239, 0.245715}, {0.0413205,
      0.211811}, {0.0471146, 0.182846}, {0.0537211, 0.158078}, {0.061254,
      0.136891}, {0.0698433, 0.118765}, {0.0796369, 0.103247}, {0.0908038,
         0.0899689}, {0.103537, 0.0783422}, {0.118055,
      0.0681593}, {0.134609, 0.0594393}, {0.153484, 0.0521564}, {0.175006,
         0.0462156}, {0.199546, 0.0413229}, {0.227526,
      0.0372613}, {0.259431, 0.033881}, {0.295809, 0.0311329}, {0.337288,
      0.0289907}, {0.384583, 0.0274659}, {0.438511, 0.0266299}, {0.5,
      0.0266588}, {0.561489, 0.0266287}, {0.615417, 0.0285988}, {0.662712,
         0.0314269}, {0.704191, 0.0350517}, {0.740569,
      0.0395662}, {0.772474, 0.0450902}, {0.800454, 0.0517807}, {0.824994,
         0.0598343}, {0.846516, 0.0694895}, {0.865391,
      0.0810323}, {0.881945, 0.0948032}, {0.896463, 0.111205}, {0.909196,
      0.130712}, {0.920363, 0.153883}, {0.930157, 0.181372}, {0.938746,
      0.213943}, {0.946279, 0.252494}, {0.952885, 0.298203}, {0.958679,
      0.351671}, {0.963761, 0.409324}, {0.968218, 0.477187}, {0.972126,
      0.554457}, {0.975554, 0.642385}, {0.97856, 0.742408}, {0.981197,
      0.855938}, {0.983509, 0.984361}, {0.985537, 1.12917}, {0.987316,
      1.29184}, {0.988876, 1.47371}, {0.990244, 1.67615}, {0.991444,
      1.90046}, {0.992496, 2.14754}, {0.993419, 2.41836}, {0.994228,
      2.71303}, {0.994938, 3.03188}, {0.99556, 3.3737}, {0.996106,
      3.73788}, {0.996585, 4.12231}, {0.997005, 4.5241}, {0.997373,
      4.93954}, {0.997696, 5.36533}, {0.99798, 5.7979}, {0.998228,
      6.22944}, {0.998446, 6.65777}, {0.998637, 7.07661}, {0.998805,
      7.48292}, {0.998952, 7.87025}, {0.999081, 8.23594}, {0.999194,
      8.57615}, {0.999293, 8.88846}, {0.99938, 9.17203}, {0.999456,
      9.42418}, {0.999523, 9.64688}, {0.999582, 9.83995}, {0.999633,
      10.0011}, {0.999678, 10.1357}, {0.999718, 10.2459}, {0.999753,
      10.3323}, {0.999783, 10.3962}, {0.99981, 10.4435}, {0.999833,
      10.4743}, {0.999854, 10.4931}, {0.999872, 10.501}, {0.999887,
      10.5007}, {0.999901, 10.4941}, {0.999913, 10.4829}, {0.999924,
      10.4675}, {0.999933, 10.451}, {0.999942, 10.4304}, {0.999949,
      10.4113}, {0.999955, 10.3926}, {0.999961, 10.3713}, {0.999965,
      10.3556}, {0.99997, 10.3341}, {0.999973, 10.32}, {0.999977,
      10.2996}, {0.99998, 10.2829}, {0.999982, 10.2711}, {0.999984,
      10.2585}, {0.999986, 10.2451}, {0.999988, 10.2308}, {0.999989,
      10.2232}, {0.999991, 10.2071}, {0.999992, 10.1985}, {0.999993,
      10.1895}, {0.999994, 10.1801}, {0.999994, 10.1801}, {0.999995,
      10.1703}, {0.999996, 10.1603}, {0.999996, 10.1603}, {0.999997,
      10.1508}, {0.999997, 10.1508}, {0.999997, 10.1508}, {0.999998,
      10.1443}, {0.999998, 10.1443}, {0.999998, 10.1443}, {0.999999,
      10.1568}, {0.999999, 10.1568}, {0.999999, 10.1568}};

      // points for the higher values of the error bands
      double err_high[200][2] = {{1.14022*10e-6, 15.6}, {1.30011*10e-6, 15.6}, {1.48241*10e-6,
      15.6}, {1.69028*10e-6, 15.6}, {1.9273*10e-6, 15.6}, {2.19755*10e-6,
      15.5999}, {2.50569*10e-6, 15.5999}, {2.85705*10e-6,
      15.5999}, {3.25767*10e-6, 15.5999}, {3.71447*10e-6,
      15.5999}, {4.23532*10e-6, 15.5998}, {4.82921*10e-6,
      15.5998}, {5.50638*10e-6, 15.5997}, {6.2785*10e-6,
      15.5997}, {7.15888*10e-6, 15.5996}, {8.16272*10e-6,
      15.5995}, {9.30732*10e-6, 15.5993}, {0.0000106124,
      15.5992}, {0.0000121005, 15.599}, {0.0000137973,
      15.5987}, {0.000015732, 15.5984}, {0.000017938,
      15.5979}, {0.0000204533, 15.5974}, {0.0000233213,
      15.5967}, {0.0000265915, 15.5959}, {0.0000303202,
      15.5948}, {0.0000345718, 15.5935}, {0.0000394195,
      15.5918}, {0.0000449471, 15.5897}, {0.0000512496,
      15.587}, {0.000058436, 15.5837}, {0.0000666301,
      15.5795}, {0.0000759731, 15.5742}, {0.0000866263,
      15.5675}, {0.0000987733, 15.5592}, {0.000112624,
      15.5487}, {0.000128416, 15.5355}, {0.000146423,
      15.5189}, {0.000166955, 15.4981}, {0.000190365,
      15.472}, {0.000217059, 15.4393}, {0.000247496, 15.3983}, {0.0002822,
         15.3471}, {0.000321771, 15.2831}, {0.00036689,
      15.2034}, {0.000418337, 15.1043}, {0.000476997,
      14.9814}, {0.000543883, 14.8296}, {0.000620148,
      14.643}, {0.000707107, 14.4147}, {0.000806259,
      14.1374}, {0.000919315, 13.8033}, {0.00104822,
      13.4047}, {0.00119521, 12.9348}, {0.0013628, 12.3885}, {0.0015539,
      11.7633}, {0.00177179, 11.061}, {0.00202024, 10.2883}, {0.00230352,
      9.45702}, {0.00262653, 8.58438}, {0.00299483, 7.69144}, {0.00341477,
         6.80135}, {0.0038936, 5.93693}, {0.00443957, 5.11847}, {0.0050621,
      4.3619}, {0.00577192, 3.67792}, {0.00658127, 3.07189}, {0.00750412,
      2.54439}, {0.00855636, 2.09235}, {0.00975616, 1.7101}, {0.0111242,
      1.3905}, {0.0126841, 1.12578}, {0.0144627, 0.908237}, {0.0164906,
      0.727761}, {0.018803, 0.581304}, {0.0214396, 0.482044}, {0.0244459,
      0.415775}, {0.0278738, 0.356327}, {0.0317824, 0.305523}, {0.036239,
      0.262287}, {0.0413205, 0.225431}, {0.0471146, 0.194018}, {0.0537211,
         0.167316}, {0.061254, 0.14454}, {0.0698433, 0.125099}, {0.0796369,
      0.108495}, {0.0908038, 0.094326}, {0.103537, 0.0819665}, {0.118055,
      0.0711902}, {0.134609, 0.061999}, {0.153484, 0.0543494}, {0.175006,
      0.0481302}, {0.199546, 0.0430273}, {0.227526, 0.0388008}, {0.259431,
         0.0352873}, {0.295809, 0.0324343}, {0.337288, 0.030217}, {0.384583,
         0.0286472}, {0.438511, 0.0278006}, {0.5, 0.0278633}, {0.561489,
      0.0281555}, {0.615417, 0.0302384}, {0.662712, 0.0332074}, {0.704191,
         0.0370211}, {0.740569, 0.0417724}, {0.772474,
      0.0475863}, {0.800454, 0.0546263}, {0.824994, 0.0630977}, {0.846516,
         0.0732486}, {0.865391, 0.0853764}, {0.881945,
      0.0998347}, {0.896463, 0.11704}, {0.909196, 0.137484}, {0.920363,
      0.161742}, {0.930157, 0.190486}, {0.938746, 0.2245}, {0.946279,
      0.264701}, {0.952885, 0.312421}, {0.958679, 0.367707}, {0.963761,
      0.425388}, {0.968218, 0.494427}, {0.972126, 0.573641}, {0.975554,
      0.664938}, {0.97856, 0.772384}, {0.981197, 0.892013}, {0.983509,
      1.02653}, {0.985537, 1.17852}, {0.987316, 1.3496}, {0.988876,
      1.54126}, {0.990244, 1.75476}, {0.991444, 1.99078}, {0.992496,
      2.24881}, {0.993419, 2.52751}, {0.994228, 2.82376}, {0.994938,
      3.13375}, {0.99556, 3.45088}, {0.996106, 3.78308}, {0.996585,
      4.19338}, {0.997005, 4.68023}, {0.997373, 5.23372}, {0.997696,
      5.83477}, {0.99798, 6.47909}, {0.998228, 7.15381}, {0.998446,
      7.85235}, {0.998637, 8.56014}, {0.998805, 9.26649}, {0.998952,
      9.95397}, {0.999081, 10.6115}, {0.999194, 11.2265}, {0.999293,
      11.7897}, {0.99938, 12.296}, {0.999456, 12.7385}, {0.999523,
      13.1198}, {0.999582, 13.4397}, {0.999633, 13.6958}, {0.999678,
      13.8984}, {0.999718, 14.0534}, {0.999753, 14.1636}, {0.999783,
      14.2344}, {0.99981, 14.2754}, {0.999833, 14.2905}, {0.999854,
      14.2863}, {0.999872, 14.2675}, {0.999887, 14.2406}, {0.999901,
      14.206}, {0.999913, 14.1691}, {0.999924, 14.1248}, {0.999933,
      14.1202}, {0.999942, 14.1818}, {0.999949, 14.2509}, {0.999955,
      14.3081}, {0.999961, 14.3652}, {0.999965, 14.4031}, {0.99997,
      14.4502}, {0.999973, 14.4781}, {0.999977, 14.5148}, {0.99998,
      14.5418}, {0.999982, 14.5594}, {0.999984, 14.5766}, {0.999986,
      14.5932}, {0.999988, 14.6091}, {0.999989, 14.6168}, {0.999991,
      14.6312}, {0.999992, 14.6379}, {0.999993, 14.6441}, {0.999994,
      14.65}, {0.999994, 14.65}, {0.999995, 14.6553}, {0.999996,
      14.6602}, {0.999996, 14.6602}, {0.999997, 14.6651}, {0.999997,
      14.6651}, {0.999997, 14.6651}, {0.999998, 14.6723}, {0.999998,
      14.6723}, {0.999998, 14.6723}, {0.999999, 14.696}, {0.999999,
      14.696}, {0.999999, 14.696}};

         double err_low[200][2] = {{1.14022*10e-6, 5.79982}, {1.30011*10e-6, 5.79979}, {1.48241*10e-6,
         5.79974}, {1.69028*10e-6, 5.79969}, {1.9273*10e-6,
         5.79962}, {2.19755*10e-6, 5.79954}, {2.50569*10e-6,
         5.79944}, {2.85705*10e-6, 5.79932}, {3.25767*10e-6,
         5.79917}, {3.71447*10e-6, 5.799}, {4.23532*10e-6,
         5.79879}, {4.82921*10e-6, 5.79853}, {5.50638*10e-6,
         5.79822}, {6.2785*10e-6, 5.79784}, {7.15888*10e-6,
         5.79739}, {8.16272*10e-6, 5.79683}, {9.30732*10e-6,
         5.79617}, {0.0000106124, 5.79536}, {0.0000121005,
         5.79438}, {0.0000137973, 5.7932}, {0.000015732,
         5.79178}, {0.000017938, 5.79006}, {0.0000204533,
         5.78798}, {0.0000233213, 5.78548}, {0.0000265915,
         5.78245}, {0.0000303202, 5.77881}, {0.0000345718,
         5.77442}, {0.0000394195, 5.76913}, {0.0000449471,
         5.76277}, {0.0000512496, 5.75511}, {0.000058436,
         5.74591}, {0.0000666301, 5.73486}, {0.0000759731,
         5.72161}, {0.0000866263, 5.70572}, {0.0000987733,
         5.68671}, {0.000112624, 5.66398}, {0.000128416,
         5.63685}, {0.000146423, 5.60453}, {0.000166955,
         5.56611}, {0.000190365, 5.52055}, {0.000217059,
         5.46666}, {0.000247496, 5.40313}, {0.0002822,
         5.32855}, {0.000321771, 5.24138}, {0.00036689,
         5.14005}, {0.000418337, 5.02299}, {0.000476997,
         4.88881}, {0.000543883, 4.73634}, {0.000620148,
         4.56489}, {0.000707107, 4.37447}, {0.000806259,
         4.16599}, {0.000919315, 3.94149}, {0.00104822,
         3.70432}, {0.00119521, 3.45909}, {0.0013628, 3.21159}, {0.0015539,
         2.96817}, {0.00177179, 2.73519}, {0.00202024, 2.51806}, {0.00230352,
            2.32044}, {0.00262653, 2.14367}, {0.00299483,
         1.98666}, {0.00341477, 1.84631}, {0.0038936, 1.71833}, {0.00443957,
         1.59831}, {0.0050621, 1.48251}, {0.00577192, 1.36849}, {0.00658127,
         1.25518}, {0.00750412, 1.14275}, {0.00855636, 1.03218}, {0.00975616,
            0.924896}, {0.0111242, 0.822379}, {0.0126841, 0.72592}, {0.0144627,
            0.636492}, {0.0164906, 0.554956}, {0.018803, 0.481047}, {0.0214396,
            0.414429}, {0.0244459, 0.356729}, {0.0278738,
         0.307449}, {0.0317824, 0.265257}, {0.036239, 0.229142}, {0.0413205,
         0.19819}, {0.0471146, 0.171674}, {0.0537211, 0.148839}, {0.061254,
         0.129242}, {0.0698433, 0.112431}, {0.0796369,
         0.0979989}, {0.0908038, 0.0856118}, {0.103537,
         0.0747178}, {0.118055, 0.0651285}, {0.134609, 0.0568796}, {0.153484,
            0.0499634}, {0.175006, 0.044301}, {0.199546, 0.0396184}, {0.227526,
            0.0357218}, {0.259431, 0.0324747}, {0.295809,
         0.0298314}, {0.337288, 0.0277645}, {0.384583, 0.0262847}, {0.438511,
            0.0254593}, {0.5, 0.0254542}, {0.561489, 0.0251019}, {0.615417,
         0.0269592}, {0.662712, 0.0296464}, {0.704191, 0.0330823}, {0.740569,
            0.03736}, {0.772474, 0.0425942}, {0.800454, 0.0489351}, {0.824994,
         0.056571}, {0.846516, 0.0657305}, {0.865391, 0.0766881}, {0.881945,
         0.0897716}, {0.896463, 0.105369}, {0.909196, 0.12394}, {0.920363,
         0.146024}, {0.930157, 0.172258}, {0.938746, 0.203386}, {0.946279,
         0.240286}, {0.952885, 0.283985}, {0.958679, 0.335635}, {0.963761,
         0.393259}, {0.968218, 0.459947}, {0.972126, 0.535273}, {0.975554,
         0.619832}, {0.97856, 0.712431}, {0.981197, 0.819863}, {0.983509,
         0.942188}, {0.985537, 1.07983}, {0.987316, 1.23407}, {0.988876,
         1.40617}, {0.990244, 1.59755}, {0.991444, 1.81015}, {0.992496,
         2.04627}, {0.993419, 2.30921}, {0.994228, 2.6023}, {0.994938,
         2.93002}, {0.99556, 3.29652}, {0.996106, 3.69268}, {0.996585,
         4.05124}, {0.997005, 4.36798}, {0.997373, 4.64536}, {0.997696,
         4.89588}, {0.99798, 5.11671}, {0.998228, 5.30506}, {0.998446,
         5.46319}, {0.998637, 5.59308}, {0.998805, 5.69934}, {0.998952,
         5.78653}, {0.999081, 5.86035}, {0.999194, 5.92579}, {0.999293,
         5.98723}, {0.99938, 6.04808}, {0.999456, 6.10989}, {0.999523,
         6.17399}, {0.999582, 6.24024}, {0.999633, 6.30652}, {0.999678,
         6.3729}, {0.999718, 6.43849}, {0.999753, 6.50094}, {0.999783,
         6.55795}, {0.99981, 6.61152}, {0.999833, 6.65799}, {0.999854,
         6.69998}, {0.999872, 6.73445}, {0.999887, 6.76087}, {0.999901,
         6.78225}, {0.999913, 6.79668}, {0.999924, 6.81033}, {0.999933,
         6.78183}, {0.999942, 6.67909}, {0.999949, 6.57172}, {0.999955,
         6.47705}, {0.999961, 6.37749}, {0.999965, 6.30819}, {0.99997,
         6.21802}, {0.999973, 6.16187}, {0.999977, 6.08438}, {0.99998,
         6.02408}, {0.999982, 5.98275}, {0.999984, 5.94044}, {0.999986,
         5.89704}, {0.999988, 5.85246}, {0.999989, 5.82968}, {0.999991,
         5.78307}, {0.999992, 5.7592}, {0.999993, 5.73494}, {0.999994,
         5.71032}, {0.999994, 5.71032}, {0.999995, 5.68541}, {0.999996,
         5.66051}, {0.999996, 5.66051}, {0.999997, 5.63646}, {0.999997,
         5.63646}, {0.999997, 5.63646}, {0.999998, 5.61636}, {0.999998,
         5.61636}, {0.999998, 5.61636}, {0.999999, 5.61749}, {0.999999,
         5.61749}, {0.999999, 5.61749}};


   TH1D HTemp("HTemp", "HTemp", 200, 0, 200);

   TGraphAsymmErrors graph(&HTemp);


   for(int i = 1; i <= HTemp.GetNbinsX(); i++)
   {
      int iGraph = i-1;
      graph.SetPoint(iGraph, HTemp.GetBinCenter(i), centralvals[iGraph][1]);
      graph.SetPointError(iGraph, HTemp.GetBinWidth(i)/2, HTemp.GetBinWidth(i)/2, centralvals[iGraph][1] - err_low[iGraph][1], err_high[iGraph][1] - centralvals[iGraph][1]);
   }

   static vector<int> Colors = GetCVDColors6();

   graph.SetFillStyle(1001);
   graph.SetMarkerStyle(20);
   graph.SetMarkerColor(Colors[4]);
   graph.SetLineColor(Colors[4]);
   graph.SetLineWidth(2);
   graph.SetFillColorAlpha(Colors[4], 0.3);


   return graph;

}

void MakeCanvasTheory(vector<TH1D > Histograms, TGraphErrors DataSyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0].GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 17 : 0.5;
   double WorldXMax = LogX ? 183: 1;

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
   DataSyst.DrawClone("2 same");
   TGraphAsymmErrors theory =  (TGraphAsymmErrors) getTheoryPlotUpdated(); 
   theory.DrawClone("3 l same");


   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : 1 / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : 1/ 2, 10000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   if(DoRatio)
      PadR.cd();

   double WorldRMin = 0.5999999;
   double WorldRMax = 1.39;//99999

   TH2D HWorldR("HWorldR", "", N, WorldXMin, WorldXMax, 100, WorldRMin, WorldRMax);
   TGraph G2;

   if(DoRatio)
   {
      HWorldR.SetStats(0);
      HWorldR.GetXaxis()->SetTickLength(0);
      HWorldR.GetXaxis()->SetLabelSize(0);
      HWorldR.GetYaxis()->SetNdivisions(505);

      HWorldR.Draw("axis");
      
      for(int i = 1; i <= Histograms[0].GetNbinsX(); i++)
      {
         int iGraph = i-1;
         DataSyst.SetPoint(iGraph,
                                 DataSyst.GetPointX(iGraph),
                                 DataSyst.GetPointY(iGraph)/Histograms[0].GetBinContent(i));
         DataSyst.SetPointError( iGraph,
                                 DataSyst.GetErrorX(iGraph),
                                 DataSyst.GetErrorY(iGraph)/Histograms[0].GetBinContent(i));
      }
      DataSyst.DrawClone("2 same");
      
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

   const int thetaBinCount = 100;
   double thetaBins[2*thetaBinCount+1];
   double thetaBinMin = 0.002;
   double thetaBinMax = M_PI / 2;

   for(int i = 0; i <= thetaBinCount; i++){
      // theta double log binning
      thetaBins[i] = exp(log(thetaBinMin) + (log(thetaBinMax) - log(thetaBinMin)) / thetaBinCount * i);
      thetaBins[2*thetaBinCount-i] = thetaBinMax * 2 - exp(log(thetaBinMin) + (log(thetaBinMax) - log(thetaBinMin)) / thetaBinCount * i);
   }

   std::cout << "Theta bins 17 is " << thetaBins[17] << std::endl;
   double BinMin    = thetaBins[17];
   double BinMiddle = M_PI/2;
   double BinMax    = M_PI-BinMin;

   Canvas.cd();
   int nDiv = 505;
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X5(MarginL, MarginB + PadHeight + PadRHeight, MarginL + PadWidth / 2, MarginB + PadHeight + PadRHeight, BinMin, BinMiddle, 510, "-GS"); // - in the draw options means we only draw axis on the "negative" side
   // axis on the x axis on the right hand side for the top of the plot
   TGaxis X6(MarginL + PadWidth, MarginB + PadHeight + PadRHeight, MarginL + PadWidth / 2, MarginB + PadHeight + PadRHeight, BinMin, BinMiddle, 510, "+GS");// - in the draw options means we only draw axis on the "negative" side
   
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL+ PadWidth / 2, MarginB + PadRHeight, MarginL+ PadWidth / 2, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");

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
   X5.SetTickSize(0.06);
   X6.SetTickSize(0.06);
   XL1.SetTickSize(0.03);
   XL2.SetTickSize(0.03);

   if(LogX == true)
   {
      //X1.Draw();
      X2.Draw();
      if(DoRatio) X3.Draw();
      if(DoRatio) X4.Draw();
      //X5.Draw();
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
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.05, MarginB - 0.01, "0.01");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.25, MarginB - 0.01, "0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.500, MarginB - 0.01, "#pi/2");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.75, MarginB - 0.01, "#pi - 0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.95, MarginB - 0.01, "#pi - 0.01");

   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kGray);
   //Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.015, "#theta_{L} = #pi/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.9, MarginB * 0.4, X.c_str());

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(90);
   Latex.SetTextColor(kBlack);
   if(DoRatio)
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "MC/Data");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL+PadWidth/2, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "2025 HB - Finalization of Result");

   TLegend Legend(0.60, 0.25, 0.70, 0.25 - 0.04 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
   {
      if (Labels[i]=="Data")
      {
         Histograms[i].SetFillStyle(DataSyst.GetFillStyle());
         Histograms[i].SetFillColor(DataSyst.GetFillColor());
      }
      Legend.AddEntry(&Histograms[i], Labels[i].c_str(),
                      (Labels[i]=="Data")? "plf": "pl");
   }
   Legend.Draw();

   TLegend Legend2(0.55, 0.90, 0.8, 0.90 - 0.04 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(&Histograms[i], Labels[i].c_str(),
                      (Labels[i]=="Data")? "plf": "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}

void MakeCanvas(vector<TH1D > Histograms, TGraphErrors DataSyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0].GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 17 : 0;
   double WorldXMax = LogX ? 183: 1;

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
   DataSyst.DrawClone("2 same");


   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : 1 / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : 1/ 2, 10000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   if(DoRatio)
      PadR.cd();

   double WorldRMin = 0.5999999;
   double WorldRMax = 1.39;//99999

   TH2D HWorldR("HWorldR", "", N, WorldXMin, WorldXMax, 100, WorldRMin, WorldRMax);
   TGraph G2;

   if(DoRatio)
   {
      HWorldR.SetStats(0);
      HWorldR.GetXaxis()->SetTickLength(0);
      HWorldR.GetXaxis()->SetLabelSize(0);
      HWorldR.GetYaxis()->SetNdivisions(505);

      HWorldR.Draw("axis");
      
      for(int i = 1; i <= Histograms[0].GetNbinsX(); i++)
      {
         int iGraph = i-1;
         DataSyst.SetPoint(iGraph,
                                 DataSyst.GetPointX(iGraph),
                                 DataSyst.GetPointY(iGraph)/Histograms[0].GetBinContent(i));
         DataSyst.SetPointError( iGraph,
                                 DataSyst.GetErrorX(iGraph),
                                 DataSyst.GetErrorY(iGraph)/Histograms[0].GetBinContent(i));
      }
      DataSyst.DrawClone("2 same");
      
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

   const int thetaBinCount = 100;
   double thetaBins[2*thetaBinCount+1];
   double thetaBinMin = 0.002;
   double thetaBinMax = M_PI / 2;

   for(int i = 0; i <= thetaBinCount; i++){
      // theta double log binning
      thetaBins[i] = exp(log(thetaBinMin) + (log(thetaBinMax) - log(thetaBinMin)) / thetaBinCount * i);
      thetaBins[2*thetaBinCount-i] = thetaBinMax * 2 - exp(log(thetaBinMin) + (log(thetaBinMax) - log(thetaBinMin)) / thetaBinCount * i);
   }

   std::cout << "Theta bins 17 is " << thetaBins[17] << std::endl;
   double BinMin    = thetaBins[17];
   double BinMiddle = M_PI/2;
   double BinMax    = M_PI-BinMin;

   Canvas.cd();
   int nDiv = 505;
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X5(MarginL, MarginB + PadHeight + PadRHeight, MarginL + PadWidth / 2, MarginB + PadHeight + PadRHeight, BinMin, BinMiddle, 510, "-GS"); // - in the draw options means we only draw axis on the "negative" side
   // axis on the x axis on the right hand side for the top of the plot
   TGaxis X6(MarginL + PadWidth, MarginB + PadHeight + PadRHeight, MarginL + PadWidth / 2, MarginB + PadHeight + PadRHeight, BinMin, BinMiddle, 510, "+GS");// - in the draw options means we only draw axis on the "negative" side
   
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");

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
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.05, MarginB - 0.01, "0.01");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.25, MarginB - 0.01, "0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.500, MarginB - 0.01, "#pi/2");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.75, MarginB - 0.01, "#pi - 0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.95, MarginB - 0.01, "#pi - 0.01");

   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kGray);
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.015, "#theta_{L} = #pi/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.9, MarginB * 0.4, X.c_str());

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(90);
   Latex.SetTextColor(kBlack);
   if(DoRatio)
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "MC/Data");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "2025 HB - Finalization of Result");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.04 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
   {
      if (Labels[i]=="Data")
      {
         Histograms[i].SetFillStyle(DataSyst.GetFillStyle());
         Histograms[i].SetFillColor(DataSyst.GetFillColor());
      }
      Legend.AddEntry(&Histograms[i], Labels[i].c_str(),
                      (Labels[i]=="Data")? "plf": "pl");
   }
   Legend.Draw();

   TLegend Legend2(0.55, 0.90, 0.8, 0.90 - 0.04 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(&Histograms[i], Labels[i].c_str(),
                      (Labels[i]=="Data")? "plf": "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}


double FindBinFraction(double Value, int NBins, double Bins[])
{
   int binIndex = FindBin(Value, NBins, Bins);

   // Handle out-of-range values
   if (binIndex < 0 || binIndex >= NBins)
      return -1.0; // Indicates out-of-range

   // Calculate the fraction within the bin
   double BinStart = Bins[binIndex];
   double BinEnd = Bins[binIndex + 1];
   return (Value - BinStart) / (BinEnd - BinStart);
}

