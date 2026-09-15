
// Full Analysis Crosscheck
// Hannah Bossi, <hannah.bossi@cern.ch>, March 6th, 2025

//#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
#include <fstream>
// using std::cout;
// using std::endl;
#include <vector>
#include <map>
using namespace std;

// ROOT include statements
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
#include "TGraph.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

// include statements for RooUnfold
#include "RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

// include statements for the tracking efficiency
#include "alephTrkEfficiency.h"
#include "SetStyle.h"
#include "ProgressBar.h"
#include "CommandLine.h"
#include "Messenger.h"
#include "JetCorrector.h"
#include "EffCorrFactor.h"
//#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;
// total center of mass energy
const double TotalE = 91.1876;
#define MAXPAIR 10000

static vector<int> Colors = GetCVDColors10();





// main function that controls everything
void FullAnalysisClosureCheck(std::string date);
// plots the distributions before matching
void preMatchingClosureCheck(std::string InputFileName, std::string GenTreeName, std::string RecoTreeName);
// plots the fake correction
void fakeCorrectionClosureCheck(std::string InputFileNameFull, std::string RecoTreeName, std::string InputFileNameMatched);
// plots the matching efficiency
void matchingEfficiencyClosureCheck(std::string InputFileNameFull, std::string GenTreeName, std::string InputFileNameMatched);
// function for the unfolding procedure (mode == 0 is trivial check, mode == 1 is the split closure test, mode == 2 is the unfolding of the data)
void unfoldingClosureCheck(std::string InputFileNameData, std::string InputFileNameFull, std::string InputFileNameMatched, int mode, int iterLow, int iterHigh, std::string date);
// function to apply the matching efficiency and see that this returns the pre-matching distribution
void applyMatchingEfficiencyClosure(std::string InputFileNameUnfolded, std::string InputFileNameFull, int iterLow, int iterHigh, std::string tag); 
// function to plot jingyu's results with ours
TH1D plotJingyuSameBinning(std::string jsonFileName); 
// function to plot the evis check
void eVisCheck(std::string InputFileName, std::string RecoTreeName);
// check of the binning issue
void binningCheck(std::string InputFileName); 
// check of the highPurity
void highPurityCheck(std::string InputFileName); 
// helper functions
int FindBin(double Value, int NBins, double Bins[]);
void DivideByBin(TH1D &H, double Bins[]);
void SetPad(TPad &P);
double FindBinFraction(double Value, int NBins, double Bins[]);
void energyCheck(std::string InputFileName, std::string RecoTreeName); 
// pythia 8 conversion electron test
void conversionElectronCheck(std::string InputFileNameWConversionElectrons, std::string InputFileNameWOConversionElectrons );
// plotting
void MakeCanvasPointerVec(vector<TH1D *> Histograms, vector<string> Labels, string Output,string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX);
void MakeCanvas(vector<TH1D> Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX);
void MakeCanvasZ(vector<TH1D> Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX); 

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

int FindBinLong(long double Value, int NBins, long double Bins[])
{
   for(int i = 0; i < NBins; i++)
      if(Value < Bins[i])
         return i - 1;
   return NBins;
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

void SetPad(TPad &P){
   P.SetLeftMargin(0);
   P.SetTopMargin(0);
   P.SetRightMargin(0);
   P.SetBottomMargin(0);
   P.SetTickx();
   P.SetTicky();
   P.Draw();
}


// ------------------------------------------------------
// function to do the weighted projection from 2D to 1D
// ------------------------------------------------------
TH1D* Project2Dto1D(TH2D* h2Temp){

    // create 1D histogram with the correct binning then reset it so that we can do things the proper way
    TH1D* h1Temp = (TH1D*)h2Temp->ProjectionX(Form("h1Temp_%s", h2Temp->GetName()));
    h1Temp->Reset();

    for (int i = 1; i <= h2Temp->GetNbinsX(); ++i) {
        double weight = 0;
        double error = 0;
        for (int j = 1; j <= h2Temp->GetNbinsY(); ++j) {
            double binContent = h2Temp->GetBinContent(i, j);
            double binError= h2Temp->GetBinError(i,j);
            double binCenter = h2Temp->GetYaxis()->GetBinCenter(j);
            weight += binContent*((binCenter));
            error += pow(binError*binCenter, 2);;
        }
        h1Temp->SetBinContent(i, weight);
        h1Temp->SetBinError(i, sqrt(error));
    }

    return h1Temp;

}


//==============================================================================
// Full Analysis Closure Check for Z and Theta
//==============================================================================

void FullAnalysisClosureCheck(std::string date = "04182025"){

   // define all the input files
   string InputFileNameFull = "/data/ALEPH/MC/LEP1MC/LEP1MCMerged.root";
   string InputFileNameMatched = "/home/hbossi/PhysicsEEJetEEC/Unfolding/20250317_Unfolding/matchingScheme2/LEP1MCMerged_Matched.root";
   string InputFileNameData = "UnfoldingInputData_03192025.root";
   string InputFileNameUnfolded = "unfoldingCheckE2C_Mode0_IterLow1_IterHigh2_04182025.root"; 
   string InputFileNameUnfoldedData = "unfoldingE2C_DataUnfolding_03242025.root";
   string JingyuResultsJSON = "results_WithGenCut.json";
   string InputFileRawData = "/data/ALEPH/Data/LEP1Data/20190602/LEP1Data1994All_recons_aftercut-MERGED.root";
   string InputFileNameWithConvEle = "/data/ALEPH/MC/PYTHI8Gen/LEP1_PYTHIA8_MC_TGenBefore.root";
   string InputFileNameWithoutConvEle = "/data/ALEPH/MC/PYTHI8Gen/LEP1_PYTHIA8_MC_TGenBefore_NoQEDConversions.root";

   
   // step by step abalysis checks
   //preMatchingClosureCheck(InputFileNameFull, "tgen", "t");
   //fakeCorrectionClosureCheck(InputFileNameFull, "t", InputFileNameMatched);
   //matchingEfficiencyClosureCheck(InputFileNameFull, "tgen", InputFileNameMatched);
   // run the unfolding
   //unfoldingClosureCheck(InputFileNameData,InputFileNameFull, InputFileNameMatched, 0, 1, 2, date);
   //applyMatchingEfficiencyClosure(InputFileNameUnfolded, InputFileNameFull, 1, 2, "ClosureCheck");
   applyMatchingEfficiencyClosure(InputFileNameUnfoldedData, InputFileNameFull, 4, 5, "DataUnfoldingJingyuComp");
   //debug(InputFileNameUnfoldedData, InputFileNameMatched, "Debugg"); 
   //eVisCheck(InputFileNameFull, "tgen");
   // TH1D test = (TH1D) plotJingyuSameBinning(JingyuResultsJSON); 
   //binningCheck(InputFileRawData); 
   //energyCheck(InputFileRawData, "t"); 
   //highPurityCheck(InputFileNameFull);

   //conversionElectronCheck(InputFileNameWithConvEle, InputFileNameWithoutConvEle ); 
}




//#ifndef __CINT__
int main () {  FullAnalysisClosureCheck(); return 0; }  // Main program when run stand-alone
//#endif


void energyCheck(std::string InputFileName, std::string RecoTreeName){

   //------------------------------------
   // define the binning
   //------------------------------------

   // theta binning
   const int BinCount = 100;
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;
   
   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

   }

   std::cout << "Reading from " << InputFileName.c_str() << std::endl;
   TFile InputFile(InputFileName.c_str());

   // create messengers for the ttrees
   ParticleTreeMessenger MReco(InputFile, RecoTreeName);

   TH1D HEEC2_reco("HEEC2_reco", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // reco level EEC after event selections
   TH1D HEEC2_recoWithECut("HEEC2_recoWithEcut", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // reco level EEC after event selections with Evis cut
   TH1D HNChargedHadrons("HNChargedHadrons", "", 60, 0, 60); // number of charged hadrons
   TH1D HNChargedHadrons_withECut("HNChargedHadrons_withECut", "", 60, 0, 60); // number of charged hadrons with energy cut 

   //------------------------------
   // loop over the reco tree
   //------------------------------
   int EntryCount = MReco.GetEntries();
   int numAcceptedRuns = 0; 
   int numHighETrackRuns = 0; 
   for(int iE = 0; iE < EntryCount; iE++){
      MReco.GetEntry(iE);
      double Evis = 0; 
      
      if(MReco.passesSTheta < 0.5) continue; // require that the distribution passSTheta
      if(MReco.passesNTrkMin < 0.5) continue; // require at least 5 tracks
      if(MReco.passesTotalChgEnergyMin < 0.5) continue; // require that the total energy is at least 15 GeV
      
      vector<FourVector> PReco;
      vector<FourVector> PReco45; 
      bool foundHighETrack = false; 
      for(int i = 0; i < MReco.nParticle; i++){
         if(MReco.charge[i] == 0) continue;
         if(MReco.highPurity[i] == false) continue;
         if(abs(cos(MReco.theta[i])) > 0.94) continue; 
         // add the energy to the total visible energy
         Evis = Evis + MReco.P[i][0]; 
         PReco.push_back(MReco.P[i]);
         if(MReco.P[i][0] < 45) PReco45.push_back(MReco.P[i]);
         if(MReco.P[i][0] > 45 && MReco.P[i][0] < 46){
            foundHighETrack = true; 
            //std::cout << "Found a track of high energy " << MReco.P[i][0] << " in an event with mult " << MReco.nChargedHadronsHP << std::endl;
         }
      }

      if(Evis > 200) continue;
      numAcceptedRuns++; 
      
      HNChargedHadrons.Fill(MReco.nChargedHadronsHP); 
      if(foundHighETrack){
         HNChargedHadrons_withECut.Fill(MReco.nChargedHadronsHP); 
         numHighETrackRuns++; 
      }
   

      for(int i = 0; i < PReco.size(); i++){
         for(int j = i+1; j < PReco.size();j++){
            if(i == j) continue; // don't fill the EEC with self correlations particles with themselves
            FourVector Reco1 = PReco.at(i);
            FourVector Reco2 = PReco.at(j);
            double recoTheta = GetAngle(Reco1,Reco2);
            int BinThetaReco = FindBin(recoTheta, 2 * BinCount, Bins);
            double  recoEEC  = Reco1[0]*Reco2[0]/(TotalE*TotalE);
            HEEC2_reco.Fill(BinThetaReco, recoEEC);
         }
      }
      
      for(int i = 0; i < PReco45.size(); i++){
         for(int j = i+1; j < PReco45.size();j++){
            if(i == j) continue; // don't fill the EEC with self correlations particles with themselves
            FourVector Reco1 = PReco45.at(i);
            FourVector Reco2 = PReco45.at(j);
            double recoTheta = GetAngle(Reco1,Reco2);
            int BinThetaReco = FindBin(recoTheta, 2 * BinCount, Bins);
            double  recoEEC  = Reco1[0]*Reco2[0]/(TotalE*TotalE);
            HEEC2_recoWithECut.Fill(BinThetaReco, recoEEC);
         }
      }
      
      
      
    } // end loop over the gen tree
    
    std::cout << "In total there were " << numAcceptedRuns << " accepted runs for a fraction of " << float(numAcceptedRuns)/EntryCount << std::endl; 

   // now do all the plotting
   // set the color and drawing style
   HEEC2_reco.SetMarkerColor(Colors[0]);
   HEEC2_reco.SetMarkerStyle(20);
   HEEC2_reco.SetLineColor(Colors[0]);
   HEEC2_reco.SetLineWidth(2);

   HEEC2_recoWithECut.SetMarkerColor(Colors[1]);
   HEEC2_recoWithECut.SetMarkerStyle(20);
   HEEC2_recoWithECut.SetLineColor(Colors[1]);
   HEEC2_recoWithECut.SetLineWidth(2);

   // divide by the bin width
   DivideByBin(HEEC2_reco, Bins);
   DivideByBin(HEEC2_recoWithECut, Bins);
   
   // scale by the number of events
   HEEC2_reco.Scale(1.0/numAcceptedRuns);
   HEEC2_recoWithECut.Scale(1.0/numAcceptedRuns);


   std::vector<TH1D> hists_baseline = {HEEC2_reco, HEEC2_recoWithECut};
   MakeCanvas(hists_baseline, {"Nominal", "With E < 45 GeV"}, "EnergyCheck",  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 1, true, true);

   TCanvas* c = new TCanvas("c", "c", 600, 600);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetLogy(); 
   c->SetRightMargin(0.05);
   c->SetLeftMargin(0.13);
   gStyle->SetOptStat(0); 
   gStyle->SetOptTitle(0); 
   
   TLegend* leg = new TLegend(0.16, 0.16, 0.5, 0.30);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.035);
   
   HNChargedHadrons.SetMarkerColor(Colors[0]); 
   HNChargedHadrons.SetMarkerStyle(20); 
   HNChargedHadrons.SetLineColor(Colors[0]); 
   
   HNChargedHadrons_withECut.SetMarkerColor(Colors[1]); 
   HNChargedHadrons_withECut.SetMarkerStyle(20); 
   HNChargedHadrons_withECut.SetLineColor(Colors[1]); 
   
   HNChargedHadrons.Scale(1.0/numAcceptedRuns);
   HNChargedHadrons_withECut.Scale(1.0/numHighETrackRuns);
   
   std::cout << "Total events " << numAcceptedRuns << " num events with a high E track " << numHighETrackRuns << std::endl;
      
   leg->AddEntry(&HNChargedHadrons, "All Events"); 
   leg->AddEntry(&HNChargedHadrons_withECut, "Events with High Energy Track");
   
   HNChargedHadrons.GetXaxis()->SetTitle("Charged Hadron Multiplicity"); 
   HNChargedHadrons.GetYaxis()->SetTitle("Per-Event Norm."); 
   HNChargedHadrons.Draw(); 
   HNChargedHadrons_withECut.Draw("same"); 

   leg->Draw();
   c->SaveAs("ChargedHadronMultTest.pdf");

   // close the input files
   InputFile.Close();

}

// ----------------------------------------------
// Check to see the impact of the binning
// ----------------------------------------------
void binningCheck(std::string InputFileName){

   const int BinCount = 100;
   long double jingyuBins[2*BinCount+1];  
     std::cout << std::setprecision (17);
   // define the jingyu bins in a hard coded way
   std::vector< long double> jingyuBinsVec = {0.0020000000000000005, 0.002137868038148267, 0.0022852398742679594, 0.0024427706433497163, 0.0026111606414721184, 0.0027911584389369833, 0.0029835642080055934, 0.003189233280029149, 0.0034090799477865387, 0.003644081529932501, 0.0038952827156245647, 0.00416380020864257, 0.004450827691646017, 0.004757641132637624, 0.005085604457222746, 0.005436175611880437, 0.0058109130452001405, 0.006211482635896096, 0.0066396650983975986, 0.007097363898936396, 0.007586613717321741, 0.008109589492019679, 0.008668616088745954, 0.009266178635553906, 0.009904933570361504, 0.010587720450028827, 0.011317574573482704, 0.012097740475004087, 0.012931686347661932, 0.013823119461012367, 0.014776002641601817, 0.01579457188953743, 0.016883355209438573, 0.01804719273948138, 0.019291258268019335, 0.020621082233431023, 0.022042576309439614, 0.02356205968519756, 0.025186287156962837, 0.02692247915624751, 0.02877835384792723, 0.0307621614410024, 0.03288272086453801, 0.03514945897182347, 0.03757245244703261, 0.040162472600678324, 0.04293103325299784, 0.045890441918132234, 0.049053854516637144, 0.05243533385954677, 0.0560499121639793, 0.05991365787819455, 0.0640437471131711, 0.0684585399982494, 0.07317766230027604, 0.07822209266908371, 0.08361425589715295, 0.08937812260803674, 0.09553931581670937, 0.10212522483554805, 0.10916512703231186, 0.11669031798138736, 0.12473425058688299, 0.1333326837960369, 0.14252384156403827, 0.15234858277693245, 0.16285058288799467, 0.1740765280750294, 0.18607632278171227, 0.1989033115655914, 0.21261451623896224, 0.22727088935681625, 0.24293758517873423, 0.259684249309269, 0.2775853283044061, 0.2967204006204415, 0.3171745303764955, 0.3390386455032981, 0.36240994195929077, 0.387392315810968, 0.41409682509825396, 0.442642183538115, 0.4731552882611474, 0.5057717839271684, 0.5406366657275616, 0.5779049239550011, 0.6177422330059507, 0.6603256878788806, 0.7058445914422634, 0.7545012959721182, 0.8065121027001184, 0.8621082233711677, 0.9215368080850028, 0.9850620439910504, 1.0529663297207341, 1.1255515307781232, 1.203140321469701, 1.2860776193387522, 1.3747321184810655, 1.4694979285582628, 1.5707963267948974, 1.6720947250315303, 1.7668605351087276, 1.855515034251041, 1.9384523321200922, 2.0160411228116697, 2.088626323869059, 2.1565306095987427, 2.22005584550479, 2.2794844302186252, 2.335080550889675, 2.387091357617675, 2.4357480621475296, 2.4812669657109128, 2.5238504205838423, 2.563687729634792, 2.6009559878622315, 2.635820869662625, 2.6684373653286455, 2.698950470051678, 2.727495828491539, 2.7542003377788253, 2.7791827116305026, 2.802554008086495, 2.8244181232132974, 2.8448722529693518, 2.864007325285387, 2.8819084042805243, 2.898655068411059, 2.914321764232977, 2.928978137350831, 2.9426893420242015, 2.9555163308080807, 2.9675161255147637, 2.9787420707017986, 2.9892440708128607, 2.999068812025755, 3.008259969793756, 3.01685840300291, 3.024902335608406, 3.0324275265574814, 3.039467428754245, 3.046053337773084, 3.0522145309817565, 3.05797839769264, 3.0633705609207094, 3.068414991289517, 3.0731341135915438, 3.077548906476622, 3.0816789957115986, 3.085542741425814, 3.0891573197302464, 3.092538799073156, 3.095702211671661, 3.098661620336795, 3.101430180989115, 3.1040202011427604, 3.10644319461797, 3.108709932725255, 3.1108304921487906, 3.1128142997418657, 3.1146701744335457, 3.1164063664328303, 3.1180305939045954, 3.1195500772803535, 3.120971571356362, 3.1223013953217738, 3.1235454608503117, 3.1247092983803544, 3.1257980817002555, 3.126816650948191, 3.1277695341287806, 3.1286609672421313, 3.129494913114789, 3.1302750790163105, 3.1310049331397645, 3.1316877200194315, 3.132326474954239, 3.132924037501047, 3.1334830640977733, 3.1340060398724714, 3.134495289690857, 3.1349529884913956, 3.135381170953897, 3.135781740544593, 3.1361564779779125, 3.1365070491325704, 3.1368350124571553, 3.137141825898147, 3.1374288533811505, 3.1376973708741684, 3.1379485720598606, 3.1381835736420065, 3.138403420309764, 3.1386090893817875, 3.1388014951508563, 3.138981492948321, 3.1391498829464433, 3.139307413715525, 3.139454785551645, 3.1395926535897933};
   for(int s = 0; s < jingyuBinsVec.size(); s++){
    jingyuBins[s] = jingyuBinsVec.at(s);
    std::cout << "jingyuBins[" << s << "]: " << jingyuBins[s] << std::endl;
   }
   
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;
   
   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
   }
   
    //std::cout << "hannahBins[" << i << "]: " << Bins[i] << std::endl;
   std::cout << "Reading from " << InputFileName.c_str() << std::endl;
   TFile InputFile(InputFileName.c_str());

   ParticleTreeMessenger MReco(InputFile, "t");
   TH1D HEEC2_jingyuBins("HEEC2_jingyuBins", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // jingyu binning
   TH1D HEEC2_hannahBins("HEEC2_hannahBins", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // hannah binning

   //------------------------------
   // loop over the reco tree
   //------------------------------
   int EntryCount = MReco.GetEntries();
   int numAcceptedRuns = 0; 
   for(int iE = 0; iE < EntryCount; iE++){
      MReco.GetEntry(iE);
      double Evis = 0; 
      
      if(MReco.passesSTheta < 0.5) continue; // require that the distribution passSTheta
      if(MReco.passesNTrkMin < 0.5) continue; // require at least 5 tracks
      if(MReco.passesTotalChgEnergyMin < 0.5) continue; // require that the total energy is at least 15 GeV
      
      vector<FourVector> PReco;
      for(int i = 0; i < MReco.nParticle; i++){
         if(MReco.charge[i] == 0) continue;
         if(MReco.highPurity[i] == false) continue;
         if(abs(cos(MReco.theta[i])) > 0.94) continue; 
         // add the energy to the total visible energy
         Evis = Evis + MReco.P[i][0]; 
         PReco.push_back(MReco.P[i]);
      }

      if(Evis > 200) continue; 

      numAcceptedRuns++;

      for(int i = 0; i < PReco.size(); i++){
         for(int j = i+1; j < PReco.size();j++){
            if(i == j) continue; // don't fill the EEC with self correlations particles with themselves
            FourVector Reco1 = PReco.at(i);
            FourVector Reco2 = PReco.at(j);
            double recoTheta = GetAngle(Reco1,Reco2);
            int BinThetaReco = FindBin(recoTheta, 2 * BinCount, Bins);
            int BinJingyu = FindBinLong((long double)recoTheta, 2 * BinCount, jingyuBins);
            double  recoEEC  = Reco1[0]*Reco2[0]/(TotalE*TotalE);
            HEEC2_hannahBins.Fill(BinThetaReco, recoEEC);
            HEEC2_jingyuBins.Fill(BinJingyu, recoEEC);
            if(BinThetaReco != BinJingyu) std::cout << "Warning:  Found a case where the bins were different" << std::endl;

         }
      }
    } // end loop over the gen tree
    

   // now do all the plotting
   // set the color and drawing style
   HEEC2_hannahBins.SetMarkerColor(Colors[0]);
   HEEC2_hannahBins.SetMarkerStyle(20);
   HEEC2_hannahBins.SetLineColor(Colors[0]);
   HEEC2_hannahBins.SetLineWidth(2);

   HEEC2_jingyuBins.SetMarkerColor(Colors[1]);
   HEEC2_jingyuBins.SetMarkerStyle(20);
   HEEC2_jingyuBins.SetLineColor(Colors[1]);
   HEEC2_jingyuBins.SetLineWidth(2);

   // divide by the bin width
   DivideByBin(HEEC2_hannahBins, Bins);
   DivideByBin(HEEC2_jingyuBins, Bins);
   
   // scale by the number of events
   HEEC2_hannahBins.Scale(1.0/numAcceptedRuns);
   HEEC2_jingyuBins.Scale(1.0/numAcceptedRuns);


   std::vector<TH1D> hists_baseline = {HEEC2_hannahBins, HEEC2_jingyuBins};
   MakeCanvas(hists_baseline, {"Nominal Binning", "Crosscheck Binning"}, "BinningCheck",  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 1, true, true);

   // close the input files
   InputFile.Close();
   
   
}

// ----------------------------------------------
// conversionElectronCheck
// ----------------------------------------------
void conversionElectronCheck(std::string InputFileNameWConversionElectrons, std::string InputFileNameWOConversionElectrons ){

   //------------------------------------
   // define the binning
   //------------------------------------

   // theta binning
   const int BinCount = 100;
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;
   
   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

   }

   std::cout << "Opening file " << InputFileNameWConversionElectrons.c_str() << " with conversion electrons" << std::endl;
   TFile InputFileWithConversionElectrons(InputFileNameWConversionElectrons.c_str());

   std::cout << "Opening file " << InputFileNameWOConversionElectrons.c_str() << " without conversion electrons" << std::endl;
   TFile InputFileWithoutConversionElectrons(InputFileNameWOConversionElectrons.c_str());
   
   // create messengers for the ttrees
   ParticleTreeMessenger MGenWithConvEle(InputFileWithConversionElectrons, "tgenBefore");
   ParticleTreeMessenger MGenWithoutConvEle(InputFileWithoutConversionElectrons, "tgenBefore");

   TH1D HEEC2_genWithConvEle("HEEC2_genWithConvEle", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // gen level without conversion electrons
   TH1D HEEC2_genWithoutConvEle("HEEC2_genWithoutHighPurityReq", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // gen level without conversion electrons

   //------------------------------
   // loop over the gen tree
   //------------------------------
   int EntryCount = MGenWithConvEle.GetEntries();
   for(int iE = 0; iE < EntryCount; iE++){
      MGenWithConvEle.GetEntry(iE);
      if(iE % 10000 == 0) std::cout << "On with Conversion Electron Case, Event: " << iE << std::endl;
      vector<FourVector> PGen;
      for(int i = 0; i < MGenWithConvEle.nParticle; i++){
         if(MGenWithConvEle.isCharged[i] == 0) continue;
         PGen.push_back(MGenWithConvEle.P[i]);
      }
      
      // loop over the collection of particles with the high purity cut
      for(int i = 0; i < PGen.size(); i++){
         for(int j = i+1; j < PGen.size();j++){
            if(i == j) continue; // don't fill the EEC with self correlations particles with themselves
            FourVector Gen1 = PGen.at(i);
            FourVector Gen2 = PGen.at(j);
            double genTheta = GetAngle(Gen1,Gen2);
            int BinThetaGen = FindBin(genTheta, 2 * BinCount, Bins);
            double  genEEC  = Gen1[0]*Gen2[0]/(TotalE*TotalE);
            HEEC2_genWithConvEle.Fill(BinThetaGen, genEEC);
         }
      }
    } // end loop over the gen tree
    
 

   // now do all the plotting
   // set the color and drawing style
   HEEC2_genWithConvEle.SetMarkerColor(Colors[0]);
   HEEC2_genWithConvEle.SetMarkerStyle(20);
   HEEC2_genWithConvEle.SetLineColor(Colors[0]);
   HEEC2_genWithConvEle.SetLineWidth(2);
   
   std::cout << "~~~~~~ Moving on to the without conversion Electron case ~~~~~~" << std::endl;
   
      //------------------------------
   // loop over the gen tree
   //------------------------------
   int EntryCount2 = MGenWithoutConvEle.GetEntries();
   for(int iE = 0; iE < EntryCount2; iE++){
      MGenWithoutConvEle.GetEntry(iE);
      
      if(iE % 10000 == 0) std::cout << "On without Conversion Electron Case, Event: " << iE << std::endl;

      vector<FourVector> PGenWithoutEle;
      for(int i = 0; i < MGenWithoutConvEle.nParticle; i++){
         if(MGenWithoutConvEle.isCharged[i] == 0) continue;
         PGenWithoutEle.push_back(MGenWithoutConvEle.P[i]);
      }
      
      // loop over the collection of particles with the high purity cut
      for(int i = 0; i < PGenWithoutEle.size(); i++){
         for(int j = i+1; j < PGenWithoutEle.size();j++){
            if(i == j) continue; // don't fill the EEC with self correlations particles with themselves
            FourVector Gen1 = PGenWithoutEle.at(i);
            FourVector Gen2 = PGenWithoutEle.at(j);
            double genTheta = GetAngle(Gen1,Gen2);
            int BinThetaGen = FindBin(genTheta, 2 * BinCount, Bins);
            double  genEEC  = Gen1[0]*Gen2[0]/(TotalE*TotalE);
            HEEC2_genWithoutConvEle.Fill(BinThetaGen, genEEC);
         }
      }
    } // end loop over the gen tree

   HEEC2_genWithoutConvEle.SetMarkerColor(Colors[1]);
   HEEC2_genWithoutConvEle.SetMarkerStyle(20);
   HEEC2_genWithoutConvEle.SetLineColor(Colors[1]);
   HEEC2_genWithoutConvEle.SetLineWidth(2);

   // divide by the bin width
   DivideByBin(HEEC2_genWithConvEle, Bins);
   DivideByBin(HEEC2_genWithoutConvEle, Bins);
   
   // scale by the number of events
   HEEC2_genWithConvEle.Scale(1.0/EntryCount);
   HEEC2_genWithoutConvEle.Scale(1.0/EntryCount2);


   std::vector<TH1D> hists_baseline = {HEEC2_genWithConvEle, HEEC2_genWithoutConvEle};
   MakeCanvas(hists_baseline, {"PYTHIA 8 With QED Conv. Ele.", "PYTHIA 8 Without QED Conv. Ele."}, "ConvElectron",  "#theta_{L, gen}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 3, true, true);

   // close the input files
   InputFileWithConversionElectrons.Close();
   InputFileWithoutConversionElectrons.Close(); 
}

// ----------------------------------------------
// HighPurity Check - Check to see the impact of a high purity selection at the gen level
// this implicitly also includes 
// ----------------------------------------------
void highPurityCheck(std::string InputFileName){

   //------------------------------------
   // define the binning
   //------------------------------------

   // theta binning
   const int BinCount = 100;
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;
   
   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

   }

   std::cout << "Reading from " << InputFileName.c_str() << std::endl;
   TFile InputFile(InputFileName.c_str());

   // create messengers for the ttrees
   ParticleTreeMessenger MGen(InputFile, "tgen");

   TH1D HEEC2_gen("HEEC2_gen", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // gen level EEC as we have it in the nominal framework (with high purity requirement)
   TH1D HEEC2_genWithoutHighPurityReq("HEEC2_genWithoutHighPurityReq", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // gen level without the high purity requirement

   //------------------------------
   // loop over the reco tree
   //------------------------------
   int EntryCount = MGen.GetEntries();
   int numAcceptedEvents = 0; 
   int nEventsWithSelectionError = 0; 
   for(int iE = 0; iE < EntryCount; iE++){
      MGen.GetEntry(iE);
      double Evis = 0; 
      
      //if(MGen.passesSTheta < 0.5) continue; // require that the distribution passSTheta

      vector<FourVector> PGen;
      vector<FourVector> PGenNoHighPurity; 
      vector<int> pwflagVec; 
      vector<int> indices; 
      int numConvElectrons = 0; 
      for(int i = 0; i < MGen.nParticle; i++){
         Evis = Evis + MGen.P[i][0]; 
         if(MGen.charge[i] == 0) continue;
         if(abs(cos(MGen.theta[i])) > 0.94) continue; 
         // see if the particle may be a conversion electron and try to also remove that
         bool areBothElectrons = (i > 0 && MGen.pwflag[i] == 2 && MGen.pwflag[i-1]  == 2); 
         bool areOppositeCharge = (i > 0 && MGen.charge[i] != MGen.charge[i-1]); 
         float conversionDPhi =  0.05;
         float conversionDTheta = 0.05;
         double deltaPhi = (TMath::Abs(MGen.theta[i] - MGen.theta[i-1])); 
         double deltaTheta = (TMath::ACos(TMath::Cos(MGen.phi[i] - MGen.phi[i-1])) ); 
         bool meetsDeltaThetaReq = 0.0001 < conversionDTheta-deltaTheta; 
         bool meetsDeltaPhiReq = 0.00001 < conversionDPhi- deltaPhi; 
         bool isConversionElectron = areBothElectrons && areOppositeCharge && meetsDeltaThetaReq && meetsDeltaPhiReq; 
         if(isConversionElectron && MGen.highPurity[i] == true){
            std::cout << "Something is wrong with the conversion electron definition " << std::endl;
            std::cout << "Delta Phi " << deltaPhi << " delta theta " << deltaTheta << std::endl; 
            std::cout << "meets delta phi " << meetsDeltaPhiReq << "meets delta theta " << meetsDeltaThetaReq << std::endl;
           // std::cout << "PWFlag: " << MGen.pwflag[i] << " and charge " << MGen.charge[i-1] << " delta phi " << (TMath::Abs(float(MGen.theta[i]) - float(MGen.theta[i-1])))  << " delta theta " << TMath::ACos(TMath::Cos(float(MGen.phi[i]) - float(MGen.phi[i-1])))  << std::endl;
         }
     
         // keep the theta cut
         float thetaCutLow = 20.*TMath::Pi()/180.;       //currently not used in highPurity
         float thetaCutHigh = 160.*TMath::Pi()/180.;     //currently not used in highPurity
         if(MGen.theta[i] > thetaCutHigh || MGen.theta[i] < thetaCutLow  ) continue;
         if (!(MGen.pwflag[i] == 0 || MGen.pwflag[i] == 1 || MGen.pwflag[i] == 2)) continue;
         if(isConversionElectron){
            numConvElectrons = numConvElectrons + 1; 
            // remove previous element from conversion electron list
            if(PGenNoHighPurity.size() > 0){
               //std::cout << "Event " << iE <<  " Found conversion electron. Size before removal is " <<  PGenNoHighPurity.size(); 
               //std::cout << "Found a conversion electron and the previous particle in the list has high purity " << MGen.highPurity[i-1] << std::endl;
               FourVector poppedParticle = PGenNoHighPurity.back(); 
               if(GetAngle(poppedParticle, MGen.P[i]) < 0.1 && pwflagVec.back() == 2 ){
                  PGenNoHighPurity.pop_back(); 
                  pwflagVec.pop_back();
                  indices.pop_back(); 
                  //std::cout << " Removing previous particle which has a pT " << poppedParticle[0] << " does it match previous? " << MGen.pt[i-1] << std::endl;
               }
               else{
                  //std::cout << " Not Removing previous particle which has a pT " << poppedParticle[0] << " does it match previous? " << MGen.pt[i-1] << " angle is " << GetAngle(poppedParticle, MGen.P[i]) << std::endl;
                }
            }
            //std::cout << " after removal size is " << PGenNoHighPurity.size() << std::endl;
            continue; 
         }
     
         // if(MGen.highPurity[i] == false){
         //    //std::cout << " ------- Event " << iE << " ------- " << std::endl;
         //    //std::cout << "Found a particle [" << i << "] " << " that was not high purity, but was not removed by the cuts. "  << " with pt: " << MGen.pt[i] << " pwflag: " << MGen.pwflag[i] << std::endl;
         //    continue; 
            
         // }
         pwflagVec.push_back(MGen.pwflag[i]);
         PGenNoHighPurity.push_back(MGen.P[i]);
         indices.push_back(i);
         if(MGen.pt[i] < 0.2) continue;
         // add the energy to the total visible energy
         PGen.push_back(MGen.P[i]);
      }
      
   
      // if(PGenNoHighPurity.size() != PGen.size()){
      //    nEventsWithSelectionError++; 
      //    std::cout << " ------- Event " << iE << " ------- " << std::endl;
      //    std::cout << "PGenNoHighPurity.size() " << PGenNoHighPurity.size() << " PGen.size() " << PGen.size() << std::endl;
      //    int size = PGenNoHighPurity.size(); 
      //    int BiggerSize = PGen.size(); 
      //    if (PGen.size() < size){
      //       size = PGen.size(); 
      //       BiggerSize = PGenNoHighPurity.size(); 
      //    }
      //    for(int p = 0; p < size; p++){
      //       FourVector GenHand = PGenNoHighPurity.at(p);
      //       FourVector GenHP = PGen.at(p);
      //       std::cout << "On particle " << p << " pT by hand : " << GenHand[0] << " pT w/ HP " << GenHP[0] << " pwflag: " << pwflagVec.at(p) << " original index: " << indices.at(p) << std::endl;
      //    }
      //    for(int ind = size; ind < BiggerSize; ind++){
      //       if(PGenNoHighPurity.size() >  PGen.size()){
      //          std::cout << "Extra element " << ind << " pT by hand : " << PGenNoHighPurity.at(ind)[0] <<  " pwflag: " << pwflagVec.at(ind) << " original index: " << indices.at(ind) << " theta: " << abs(cos(MGen.theta[indices.at(ind)])) << " Phi " << MGen.phi[indices.at(ind)] << std::endl;
      //       }
      //    }
      //    std::cout << "-------------------------------" << std::endl;
      // }
      // remove events for the laser calibration set (won't do anything for gen)
      if (Evis > 200) continue; 
      numAcceptedEvents++; 
      
      
      
      // loop first over the collection of particles without the high purity cut
      for(int i = 0; i < PGenNoHighPurity.size(); i++){
         for(int j = i+1; j < PGenNoHighPurity.size();j++){
            if(i == j) continue; // don't fill the EEC with self correlations particles with themselves
            FourVector Gen1 = PGenNoHighPurity.at(i);
            FourVector Gen2 = PGenNoHighPurity.at(j);
            double genTheta = GetAngle(Gen1,Gen2);
            int BinThetaGen = FindBin(genTheta, 2 * BinCount, Bins);
            double  genEEC  = Gen1[0]*Gen2[0]/(Evis*Evis);
            HEEC2_genWithoutHighPurityReq.Fill(BinThetaGen, genEEC);
         }
      }

      // loop over the collection of particles with the high purity cut
      for(int i = 0; i < PGen.size(); i++){
         for(int j = i+1; j < PGen.size();j++){
            if(i == j) continue; // don't fill the EEC with self correlations particles with themselves
            FourVector Gen1 = PGen.at(i);
            FourVector Gen2 = PGen.at(j);
            double genTheta = GetAngle(Gen1,Gen2);
            int BinThetaGen = FindBin(genTheta, 2 * BinCount, Bins);
            double  genEEC  = Gen1[0]*Gen2[0]/(TotalE*TotalE);
            HEEC2_gen.Fill(BinThetaGen, genEEC);
         }
      }
    } // end loop over the gen tree
    
    std::cout << "In total there were " << numAcceptedEvents << " accepted runs for a fraction of " << float(numAcceptedEvents)/EntryCount << " where " << nEventsWithSelectionError << " had the error " << std::endl; 

   // now do all the plotting
   // set the color and drawing style
   HEEC2_gen.SetMarkerColor(Colors[0]);
   HEEC2_gen.SetMarkerStyle(20);
   HEEC2_gen.SetLineColor(Colors[0]);
   HEEC2_gen.SetLineWidth(2);

   HEEC2_genWithoutHighPurityReq.SetMarkerColor(Colors[1]);
   HEEC2_genWithoutHighPurityReq.SetMarkerStyle(20);
   HEEC2_genWithoutHighPurityReq.SetLineColor(Colors[1]);
   HEEC2_genWithoutHighPurityReq.SetLineWidth(2);

   // divide by the bin width
   DivideByBin(HEEC2_gen, Bins);
   DivideByBin(HEEC2_genWithoutHighPurityReq, Bins);
   
   // scale by the number of events
   HEEC2_gen.Scale(1.0/numAcceptedEvents);
   HEEC2_genWithoutHighPurityReq.Scale(1.0/numAcceptedEvents);


   std::vector<TH1D> hists_baseline = {HEEC2_gen, HEEC2_genWithoutHighPurityReq};
   MakeCanvas(hists_baseline, {"Gen MC w/ pT cut", "Gen MC w/ pT cut"}, "PtCut",  "#theta_{L, gen}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 3, true, true);

   // close the input files
   InputFile.Close();

}

// ----------------------------------------------
// Evis Check - Check to see the impact of an Evis cut
// ----------------------------------------------
void eVisCheck(std::string InputFileName, std::string RecoTreeName){

   //------------------------------------
   // define the binning
   //------------------------------------

   // theta binning
   const int BinCount = 100;
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;
   
   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

   }

   std::cout << "Reading from " << InputFileName.c_str() << std::endl;
   TFile InputFile(InputFileName.c_str());

   // create messengers for the ttrees
   ParticleTreeMessenger MReco(InputFile, RecoTreeName);

   TH1D HEEC2_reco("HEEC2_reco", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // reco level EEC after event selections
   TH1D HEEC2_recoWithEvisCut("HEEC2_recoWithEvis", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // reco level EEC after event selections with Evis cut

   //------------------------------
   // loop over the reco tree
   //------------------------------
   int EntryCount = MReco.GetEntries();
   int numAcceptedRuns = 0; 
   for(int iE = 0; iE < EntryCount; iE++){
      MReco.GetEntry(iE);
      double Evis = 0; 
      
      if(MReco.passesSTheta < 0.5) continue; // require that the distribution passSTheta
      if(MReco.passesNTrkMin < 0.5) continue; // require at least 5 tracks
      if(MReco.passesTotalChgEnergyMin < 0.5) continue; // require that the total energy is at least 15 GeV
      
      vector<FourVector> PReco;
      for(int i = 0; i < MReco.nParticle; i++){
         if(MReco.charge[i] == 0) continue;
         if(MReco.highPurity[i] == false) continue;
         if(abs(cos(MReco.theta[i])) > 0.94) continue; 
         // add the energy to the total visible energy
         Evis = Evis + MReco.P[i][0]; 
         PReco.push_back(MReco.P[i]);
      }

      if(MReco.passesSTheta > 0.5 && Evis < 200){
               numAcceptedRuns++; 
      }

      for(int i = 0; i < PReco.size(); i++){
         for(int j = i+1; j < PReco.size();j++){
            if(i == j) continue; // don't fill the EEC with self correlations particles with themselves
            FourVector Reco1 = PReco.at(i);
            FourVector Reco2 = PReco.at(j);
            double recoTheta = GetAngle(Reco1,Reco2);
            int BinThetaReco = FindBin(recoTheta, 2 * BinCount, Bins);
            double  recoEEC  = Reco1[0]*Reco2[0]/(TotalE*TotalE);
            HEEC2_reco.Fill(BinThetaReco, recoEEC);
            if(MReco.passesSTheta > 0.5 && Evis < 200){
               HEEC2_recoWithEvisCut.Fill(BinThetaReco, recoEEC);
            }

         }
      }
    } // end loop over the gen tree
    
    std::cout << "In total there were " << numAcceptedRuns << " accepted runs for a fraction of " << float(numAcceptedRuns)/EntryCount << std::endl; 

   // now do all the plotting
   // set the color and drawing style
   HEEC2_reco.SetMarkerColor(Colors[0]);
   HEEC2_reco.SetMarkerStyle(20);
   HEEC2_reco.SetLineColor(Colors[0]);
   HEEC2_reco.SetLineWidth(2);

   HEEC2_recoWithEvisCut.SetMarkerColor(Colors[1]);
   HEEC2_recoWithEvisCut.SetMarkerStyle(20);
   HEEC2_recoWithEvisCut.SetLineColor(Colors[1]);
   HEEC2_recoWithEvisCut.SetLineWidth(2);

   // divide by the bin width
   DivideByBin(HEEC2_reco, Bins);
   DivideByBin(HEEC2_recoWithEvisCut, Bins);
   
   // scale by the number of events
   HEEC2_reco.Scale(1.0/EntryCount);
   HEEC2_recoWithEvisCut.Scale(1.0/numAcceptedRuns);


   std::vector<TH1D> hists_baseline = {HEEC2_reco, HEEC2_recoWithEvisCut};
   MakeCanvas(hists_baseline, {"Gen", "Gen w/ Passes STheta"}, "GenSthetaCheck",  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 1, true, true);

   // close the input files
   InputFile.Close();

}
// -------------------------------------------------------------------


// ----------------------------------------------
// plot jingyu's results
// ----------------------------------------------
TH1D plotJingyuSameBinning(std::string jsonFileName){
   
   const int BinCount = 100;
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;

   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

   }
   
   
   TH1D HJingyu("HEEC2_Jingyu", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount);
   
   // Read the JSON file
    std::ifstream file(jsonFileName.c_str());
    if (!file.is_open()) {
        std::cerr << "Failed to open file!" << std::endl;
        return HJingyu;
    }

    json data;
    file >> data;
    file.close();

    // Extract points
    std::vector<double> x_values, y_values;
    int iGraph = 1; 
    for (const auto& point : data) {
      // get the central value for the data-point
      double thetaRadians = point[0];
      int bin     = FindBin(thetaRadians, 2*BinCount, Bins); 
      double hannahThetaRadians = (Bins[bin]+Bins[bin+1])/2; 
      std::cout << std::setprecision (17);
      std::cout << "Bin Center Jingyu: " << thetaRadians << " Bin Center Hannah " << hannahThetaRadians << std::endl;
      HJingyu.SetBinContent(bin+1, point[1]); 
      HJingyu.SetBinError(bin+1, point[2]);
      iGraph++;
    }

    HJingyu.SetMarkerStyle(20); 
    HJingyu.SetLineColor(kBlack); 
    HJingyu.SetLineWidth(2); 
 

    return HJingyu; 
}
// ----------------------------------------------



// ----------------------------------------------
// applyMatchingEfficiencyClosure
// ----------------------------------------------
void applyMatchingEfficiencyClosure(std::string InputFileNameUnfolded, std::string InputFileNameFull, int iterLow, int iterHigh, std::string tag){
   
   // define the binning
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

   }

   // here is what I denote as binning option #1
   std::vector<double> e1e2BinsUnfolded = {0.0, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.00125, 0.0015, 0.00175, 0.002, 0.00225, 0.0025, 0.00275, 0.003, 0.0035, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.3};

   
   TFile InputFile(InputFileNameFull.c_str());

   // create messengers for the ttrees
   ParticleTreeMessenger MGen(InputFile, "tgen");
   ParticleTreeMessenger MReco(InputFile, "t"); 
   TH1D HEEC2_gen("HEEC2_gen", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // gen level EEC
   TH2D HEEC2_gen2D("HEEC2_gen2D", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount, e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data()); // gen level EEC
   TH1D HEEC2_gen_Z("HEEC2_gen_Z", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // gen level EEC as a function of Z
   TH2D HEEC2_gen2D_Z("HEEC2_gen2D_Z", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount, e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data()); // gen level EEC as a function of Z
   //------------------------------
   // loop over the gen tree
   //------------------------------
   int EntryCount = MGen.GetEntries();
   int numAcceptedEvents = 0; 
   for(int iE = 0; iE < EntryCount; iE++){
      MGen.GetEntry(iE);
      MReco.GetEntry(iE);
      
      if(MReco.passesSTheta < 0.5) continue; // require that the distribution passSTheta
      if(MReco.passesNTrkMin < 0.5) continue; // require at least 5 tracks
      if(MReco.passesTotalChgEnergyMin < 0.5) continue; // require that the total energy is at least 15 GeV
      numAcceptedEvents++; 
      
      vector<FourVector> PGen;
      for(int i = 0; i < MGen.nParticle; i++){
         if(MGen.charge[i] == 0) continue;
         if(MGen.highPurity[i] == false) continue;
         PGen.push_back(MGen.P[i]);
      }


      for(int i = 0; i < PGen.size(); i++){
         for(int j = i+1; j < PGen.size();j++){
            if(i == j) continue; // don't fill the EEC with self correlations particles with themselves
            FourVector Gen1 = PGen.at(i);
            FourVector Gen2 = PGen.at(j);
            double genTheta = GetAngle(Gen1,Gen2);
            double genZ = (1-cos(genTheta))/2;
            int BinThetaGen = FindBin(genTheta, 2 * BinCount, Bins);
            int BinZGen = FindBin(genZ, 2*BinCount, zBins);
            double  genEEC  = Gen1[0]*Gen2[0]/(TotalE*TotalE);
            HEEC2_gen.Fill(BinThetaGen, genEEC);
            HEEC2_gen2D.Fill(BinThetaGen, genEEC);
            HEEC2_gen_Z.Fill(BinZGen, genEEC);
            HEEC2_gen2D_Z.Fill(BinZGen, genEEC);
         }
      }
    } // end loop over the reco tree

   ParticleTreeMessenger MGenBefore(InputFile, "tgenBefore");
   TH1D HEEC2_genBefore("HEEC2_genBefore", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // gen level EEC
   TH1D HEEC2_genBefore_Z("HEEC2_genBeforeZ", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // gen level EEC
   TH1D HEEC2_genBeforeAccCorr("HEEC2_genBeforeAccCorr", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // gen level EEC

   //------------------------------
   // loop over the gen tree
   //------------------------------
   int EntryCountBefore = MGenBefore.GetEntries();
   for(int iE = 0; iE < EntryCount; iE++){
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
            HEEC2_genBefore.Fill(BinThetaGen, genEEC);
            HEEC2_genBefore_Z.Fill(BinZGen, genEEC);

         }
      }
    } // end loop over the reco tree


   TH1D* MCPreMatch = (TH1D*)Project2Dto1D(&HEEC2_gen2D);
   // now do all the plotting
   // set the color and drawing style
   MCPreMatch->SetMarkerColor(Colors[1]);
   MCPreMatch->SetMarkerStyle(20);
   MCPreMatch->SetLineColor(Colors[1]);
   MCPreMatch->SetLineWidth(2);
   // divide by the bin width
   DivideByBin(*MCPreMatch, Bins);

   // scale by the number of events
   MCPreMatch->Scale(1.0/numAcceptedEvents);
   
   
   TH1D* MCPreMatch_Z = (TH1D*)Project2Dto1D(&HEEC2_gen2D_Z);
   // now do all the plotting
   // set the color and drawing style
   MCPreMatch_Z->SetMarkerColor(Colors[1]);
   MCPreMatch_Z->SetMarkerStyle(20);
   MCPreMatch_Z->SetLineColor(Colors[1]);
   MCPreMatch_Z->SetLineWidth(2);
   // divide by the bin width
   DivideByBin(*MCPreMatch_Z, zBins);

   // scale by the number of events
   MCPreMatch_Z->Scale(1.0/EntryCount);
   
   
   // create the vectors that will later be used for plotting
   std::vector<TH1D> hists;
   std::vector<std::string> labels;
   std::vector<TH1D> hists_z;
   std::vector<std::string> labels_z;
   hists.push_back(*MCPreMatch); 
   labels.push_back("Archived MC - Gen");
   // hists_z.push_back(*MCPreMatch_Z); 
   // labels_z.push_back("Archived MC - Gen");
   
   // divide by the bin width
   HEEC2_gen.SetMarkerColor(Colors[0]);
   HEEC2_gen.SetMarkerStyle(20);
   HEEC2_gen.SetLineColor(Colors[0]);
   HEEC2_gen.SetLineWidth(2);
   DivideByBin(HEEC2_gen, Bins);
   // scale by the number of events
   HEEC2_gen.Scale(1.0/numAcceptedEvents);
   
   HEEC2_gen_Z.SetMarkerColor(Colors[0]);
   HEEC2_gen_Z.SetMarkerStyle(20);
   HEEC2_gen_Z.SetLineColor(Colors[0]);
   HEEC2_gen_Z.SetLineWidth(2);
   DivideByBin(HEEC2_gen_Z, zBins);
   // scale by the number of events
   HEEC2_gen_Z.Scale(1.0/EntryCount);
   
   HEEC2_gen.SetMarkerColor(Colors[0]);
   HEEC2_gen.SetMarkerStyle(20);
   HEEC2_gen.SetLineColor(Colors[0]);
   HEEC2_gen.SetLineWidth(2);
   DivideByBin(HEEC2_gen, Bins);
   // scale by the number of events
   HEEC2_gen.Scale(1.0/numAcceptedEvents);
   
   std::vector<TH1D> histsAfterBinningCorrection; 
   std::vector<std::string> labelsAfterBinningCorrection; 

   histsAfterBinningCorrection.push_back(HEEC2_gen); 
   labelsAfterBinningCorrection.push_back("Archived MC - Gen 1D");
   
   std::vector<TH1D> histsAfterBinningCorrection_Z; 
   std::vector<std::string> labelsAfterBinningCorrection_Z; 

   histsAfterBinningCorrection_Z.push_back(HEEC2_gen_Z); 
   labelsAfterBinningCorrection_Z.push_back("Archived MC - Gen 1D");
   
      
   HEEC2_genBefore_Z.SetMarkerColor(Colors[0]);
   HEEC2_genBefore_Z.SetMarkerStyle(20);
   HEEC2_genBefore_Z.SetLineColor(Colors[0]);
   HEEC2_genBefore_Z.SetLineWidth(2);
   DivideByBin(HEEC2_genBefore_Z, zBins);
   // scale by the number of events
   HEEC2_genBefore_Z.Scale(1.0/EntryCount);
   
   HEEC2_genBefore.SetMarkerColor(Colors[0]);
   HEEC2_genBefore.SetMarkerStyle(20);
   HEEC2_genBefore.SetLineColor(Colors[0]);
   HEEC2_genBefore.SetLineWidth(2);
   DivideByBin(HEEC2_genBefore, Bins);
   // scale by the number of events
   HEEC2_genBefore.Scale(1.0/EntryCount);
   
   std::vector<TH1D> histsAfterEventCorrection; 
   std::vector<std::string> labelsAfterEventCorrection; 

  
   
   std::vector<TH1D> histsAfterEventCorrection_Z; 
   std::vector<std::string> labelsAfterEventCorrection_Z; 

   histsAfterEventCorrection_Z.push_back(HEEC2_genBefore_Z); 
   labelsAfterEventCorrection_Z.push_back("Archived MC Before Event Selection");
   
   TFile* unfoldedFile = TFile::Open(InputFileNameUnfolded.c_str()); 
   
   

   
   // get the MC true histograms
   TH2D* h2true_Theta = (TH2D*)unfoldedFile->Get("true_Theta"); 
   TH1D* MCtrue = (TH1D*)Project2Dto1D(h2true_Theta);
   DivideByBin(*MCtrue, Bins);
   MCtrue->Scale(1.0/753725);
   MCtrue->SetMarkerStyle(20);
   MCtrue->SetMarkerColor(Colors[0]);
   MCtrue->SetLineColor(Colors[0]);
   MCtrue->SetLineWidth(2);
   
   // TH2D* h2true_Z = (TH2D*)unfoldedFile->Get("true_Z"); 
   // TH1D* MCtrue_Z = (TH1D*)Project2Dto1D(h2true_Z);
   // DivideByBin(*MCtrue_Z, zBins);
   // MCtrue_Z->Scale(1.0/771597);
   // MCtrue_Z->SetMarkerStyle(20);
   // MCtrue_Z->SetMarkerColor(Colors[0]);
   // MCtrue_Z->SetLineColor(Colors[0]);
   // MCtrue_Z->SetLineWidth(2);

   hists.push_back(*MCtrue);
   labels.push_back("Matched MC 2D to 1D");
   
   
   
   // hists_z.push_back(*MCtrue_Z);
   // labels_z.push_back("Matched MC 2D to 1D");
   
   for(int iter = iterLow; iter < iterHigh; iter++){
      TH2D* hunf_Theta =  (TH2D*)unfoldedFile->Get(Form("Bayesian_Unfoldediter%d_Theta", iter));
      TH2D* hUnfoldedPostCorrection = (TH2D*)hunf_Theta->Clone(Form("hUnfoldedPostCorrection_iter%d", iter));
      hUnfoldedPostCorrection->Reset(); 
      int applyEffCorrOnHistoErrorStatus_theta = 0;
      std::cout << "------------- Efficiency Correction: [Applying the matching efficiency] ----------" << std::endl;
      EffCorrFactor matchingEffCorrFactorTheta;
      matchingEffCorrFactorTheta.init("/home/hbossi/PhysicsEEJetEEC/Unfolding/20250317_Unfolding/matchingScheme2/MatchingEff.root", "theta");
      applyEffCorrOnHistoErrorStatus_theta += matchingEffCorrFactorTheta.applyEffCorrOnHisto(hunf_Theta , hUnfoldedPostCorrection);
      
      
      // TH2D* hunf_Z =  (TH2D*)unfoldedFile->Get(Form("Bayesian_Unfoldediter%d_Z", iter));
      // TH2D* hUnfoldedPostCorrection_Z = (TH2D*)hunf_Z->Clone(Form("hUnfoldedPostCorrectionZ_iter%d", iter));
      // hUnfoldedPostCorrection_Z->Reset(); 
      // int applyEffCorrOnHistoErrorStatus_Z = 0;
      // EffCorrFactor matchingEffCorrFactorZ;
      // matchingEffCorrFactorZ.init("/home/hbossi/PhysicsEEJetEEC/Unfolding/20250317_Unfolding/matchingScheme2/MatchingEff.root", "z");
      // applyEffCorrOnHistoErrorStatus_Z += matchingEffCorrFactorZ.applyEffCorrOnHisto(hunf_Z , hUnfoldedPostCorrection_Z);
      
      TH1D* hUnfoldedPreMatching =  Project2Dto1D(hunf_Theta);
      hUnfoldedPreMatching->SetName(Form("Bayesian_Unfolded1Diter%d_Theta_PreMatching",iter));
      DivideByBin(*hUnfoldedPreMatching, Bins);
      hUnfoldedPreMatching->Scale(1.0/1333529);//1333529//753725
      // probably need to also add correction for the matching efficiency here??
      hUnfoldedPreMatching->SetMarkerStyle(20);
      hUnfoldedPreMatching->SetMarkerColor(Colors[iter+1]);
      hUnfoldedPreMatching->SetLineColor(Colors[iter+1]);
      hUnfoldedPreMatching->SetLineWidth(2);
      hists.push_back(*hUnfoldedPreMatching);
      labels.push_back(Form("Unfolded Iteration %d w/o Matching Correction", iter));
      
      // TH1D* hUnfoldedPreMatching_Z =  Project2Dto1D(hunf_Z);
      // hUnfoldedPreMatching_Z->SetName(Form("Bayesian_Unfolded1Diter%d_Zsadfsdf",iter));
      // DivideByBin(*hUnfoldedPreMatching_Z, zBins);
      // hUnfoldedPreMatching_Z->Scale(1.0/1333529);
      // // probably need to also add correction for the matching efficiency here??
      // hUnfoldedPreMatching_Z->SetMarkerStyle(20);
      // hUnfoldedPreMatching_Z->SetMarkerColor(Colors[iter+1]);
      // hUnfoldedPreMatching_Z->SetLineColor(Colors[iter+1]);
      // hUnfoldedPreMatching_Z->SetLineWidth(2);
      // hists_z.push_back(*hUnfoldedPreMatching_Z);
      // labels_z.push_back(Form("Unfolded Iteration %d w/o Matching Correction", iter));
   
      TH1D* hUnfolded =  Project2Dto1D(hUnfoldedPostCorrection);
      hUnfolded->SetName(Form("Bayesian_Unfolded1Diter%d_Theta",iter));
      DivideByBin(*hUnfolded, Bins);
      hUnfolded->Scale(1.0/1333529);
      // probably need to also add correction for the matching efficiency here??
      hUnfolded->SetMarkerStyle(20);
      hUnfolded->SetMarkerColor(Colors[iter+2]);
      hUnfolded->SetLineColor(Colors[iter+2]);
      hUnfolded->SetLineWidth(2);
      hists.push_back(*hUnfolded);
      labels.push_back(Form("Unfolding Iteration %d w/ matching correction", iter));

      //labels.push_back(Form("Unfolded Iteration %d w/ Matching Correction", iter));
      
      // TH1D* hUnfolded_Z =  Project2Dto1D(hUnfoldedPostCorrection_Z);
      // hUnfolded_Z->SetName(Form("Bayesian_Unfolded1Diter%d_Z",iter));
      // DivideByBin(*hUnfolded_Z,zBins);
      // hUnfolded_Z->Scale(1.0/1333529);
      // // probably need to also add correction for the matching efficiency here??
      // hUnfolded_Z->SetMarkerStyle(20);
      // hUnfolded_Z->SetMarkerColor(Colors[iter+2]);
      // hUnfolded_Z->SetLineColor(Colors[iter+2]);
      // hUnfolded_Z->SetLineWidth(2);
      // hists_z.push_back(*hUnfolded_Z);
      // labels_z.push_back(Form("Unfolded Iteration %d w/ Matching Correction", iter));
      
      // ---------------------------------------------
      // binning correction
      // -----------------------------------------------
      
      TH1D* hUnfoldedPreBinningCorrection = (TH1D*)hUnfolded->Clone(Form("hUnfoldedPreBinningCorrection%d", iter));
      TH1D* hUnfoldedPostBinningCorrection = (TH1D*)hUnfoldedPreBinningCorrection->Clone(Form("hUnfoldedPostBinningCorrection%d", iter));
      hUnfoldedPostBinningCorrection->Reset(); 
      EffCorrFactor UnfoldingBinCorrFactor_Theta;
      UnfoldingBinCorrFactor_Theta.init("/home/hbossi/PhysicsEEJetEEC/Unfolding/20250324_UnfoldingBinningCorrection/UnfoldingBinCorr_with_theta.root", "theta");
      applyEffCorrOnHistoErrorStatus_theta += UnfoldingBinCorrFactor_Theta.applyEffCorrOnHisto(hUnfoldedPreBinningCorrection, hUnfoldedPostBinningCorrection);

      histsAfterBinningCorrection.push_back(*hUnfoldedPostBinningCorrection);
      labelsAfterBinningCorrection.push_back(Form("Unfolded Iteration %d w/ Matching and Binning Correction", iter));
      
            
      // TH1D* hUnfoldedPreBinningCorrection_Z = (TH1D*)hUnfolded_Z->Clone(Form("hUnfoldedPreBinningCorrectionZ%d", iter));
      // TH1D* hUnfoldedPostBinningCorrection_Z = (TH1D*)hUnfoldedPreBinningCorrection_Z->Clone(Form("hUnfoldedPostBinningCorrectionZ%d", iter));
      // hUnfoldedPostBinningCorrection_Z->Reset(); 
      // EffCorrFactor UnfoldingBinCorrFactor_Z;
      // UnfoldingBinCorrFactor_Z.init("/home/hbossi/PhysicsEEJetEEC/Unfolding/20250222_UnfoldingBinningCorrection/UnfoldingBinCorr_with_z.root", "z");
      // applyEffCorrOnHistoErrorStatus_theta += UnfoldingBinCorrFactor_Z.applyEffCorrOnHisto(hUnfoldedPreBinningCorrection_Z, hUnfoldedPostBinningCorrection_Z);

      // histsAfterBinningCorrection_Z.push_back(*hUnfoldedPostBinningCorrection_Z);
      // labelsAfterBinningCorrection_Z.push_back(Form("Unfolded Iteration %d w/ Matching and Binning Correction", iter));

      // ---------------------------------------------
      // event selection efficiency correction
      // -----------------------------------------------
      
      TH1D* hUnfoldedPreEventCorrection = (TH1D*)hUnfoldedPostBinningCorrection->Clone(Form("hUnfoldedPreEventCorrection%d", iter));
      TH1D* hUnfoldedPostEventCorrection = (TH1D*)hUnfoldedPostBinningCorrection->Clone(Form("hUnfoldedPostEventCorrection%d", iter));
      hUnfoldedPostEventCorrection->Reset(); 
      
      EffCorrFactor EvtSelEffCorrFactor_Theta;
      EvtSelEffCorrFactor_Theta.init("/home/hbossi/PhysicsEEJetEEC/EventSelectionEfficiency/20250324_evtSelEffCorr/EvtSelEff.root", "theta");
      applyEffCorrOnHistoErrorStatus_theta += EvtSelEffCorrFactor_Theta.applyEffCorrOnHisto(hUnfoldedPreEventCorrection, hUnfoldedPostEventCorrection); 
      // now take the post correction version and project that
  

      hUnfoldedPostEventCorrection->SetMarkerStyle(20);
      hUnfoldedPostEventCorrection->SetMarkerColor(Colors[iter+1]);
      hUnfoldedPostEventCorrection->SetLineColor(Colors[iter+1]);
      hUnfoldedPostEventCorrection->SetLineWidth(2);
      

      
      histsAfterEventCorrection.push_back(*hUnfoldedPostEventCorrection);
      labelsAfterEventCorrection.push_back(Form("Unfolded Iteration %d w/ All Corrections", iter));

      
      // TH1D* hUnfoldedPreEventCorrection_Z = (TH1D*)hUnfoldedPostBinningCorrection_Z->Clone(Form("hUnfoldedPreEventCorrectionZ%d", iter));
      // TH1D* hUnfoldedPostEventCorrection_Z = (TH1D*)hUnfoldedPostBinningCorrection_Z->Clone(Form("hUnfoldedPostEventCorrectionZ%d", iter));
      // hUnfoldedPostEventCorrection_Z->Reset(); 
      
      // EffCorrFactor EvtSelEffCorrFactor_Z;
      // EvtSelEffCorrFactor_Z.init("/home/hbossi/PhysicsEEJetEEC/EventSelectionEfficiency/20250222_evtSelEffCorr/EvtSelEff.root", "z");
      // applyEffCorrOnHistoErrorStatus_Z += EvtSelEffCorrFactor_Z.applyEffCorrOnHisto(hUnfoldedPreEventCorrection_Z, hUnfoldedPostEventCorrection_Z); 
      // // now take the post correction version and project that
  

      // hUnfoldedPostEventCorrection_Z->SetMarkerStyle(20);
      // hUnfoldedPostEventCorrection_Z->SetMarkerColor(Colors[iter+1]);
      // hUnfoldedPostEventCorrection_Z->SetLineColor(Colors[iter+1]);
      // hUnfoldedPostEventCorrection_Z->SetLineWidth(2);
      

      
      // histsAfterEventCorrection_Z.push_back(*hUnfoldedPostEventCorrection_Z);
      // labelsAfterEventCorrection_Z.push_back(Form("Unfolded Iteration %d w/ All Corrections", iter));
   }
   

   TH1D jingyu = plotJingyuSameBinning("results_050125_WithRecoStheta.json"); 
   histsAfterEventCorrection.push_back(jingyu); 
   labelsAfterEventCorrection.push_back("Jingyu Comparison - 05/01/25"); 
   
   // make the archived MC last
   // histsAfterEventCorrection.push_back(HEEC2_genBefore); 
   // labelsAfterEventCorrection.push_back("Archived MC Before Event Selection");
   
   MakeCanvas(hists, labels, Form("appyMatchingEffCheck_Theta_%s", tag.c_str()),  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 10, true, true);
   MakeCanvas(histsAfterBinningCorrection, labelsAfterBinningCorrection, Form("applyBinningCorrCheck_%s", tag.c_str()),  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 10, true, true);
   MakeCanvas(histsAfterEventCorrection, labelsAfterEventCorrection, Form("applyEventCorrCheck_%s", tag.c_str()),  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 10, true, true);

   //MakeCanvasZ(hists_z, labels_z, Form("appyMatchingEffCheck_Z_%s", tag.c_str()),  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 10, true, true);
   // MakeCanvasZ(histsAfterBinningCorrection_Z, labelsAfterBinningCorrection_Z, Form("applyBinningCorrCheck_Z_%s", tag.c_str()),  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 10, true, true);
   // MakeCanvasZ(histsAfterEventCorrection_Z, labelsAfterEventCorrection_Z, Form("applyEventCorrCheck_Z_%s", tag.c_str()),  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 10, true, true);



   
}
// 



// ----------------------------------------------
// unfoldingClosureCheck
// ----------------------------------------------
void unfoldingClosureCheck(std::string InputFileNameData, std::string InputFileNameFull, std::string InputFileNameMatched, int mode, int iterLow, int iterHigh, std::string date){

   // first print the mode to make sure that this is the intended mode by the user
   if(mode == 0){
      std::cout << "........ Running Unfolding Full Analysis Closure Check (mode 0)........ " << std::endl;
   }
   else if(mode == 1){
      std::cout << "........ Running Unfolding Trivial Closure Check (mode 1)........ " << std::endl;
   }
   else if(mode == 2){
      std::cout << "........ Running Unfolding Split Closure Check (mode 1)........ " << std::endl;
   }
   else if(mode == 3){
      std::cout << "........ Running Unfolding of the data (mode 2) ................" << std::endl;
   }
   else{
      std::cout << "[Error]: Mode not recognized" << std::endl;
   }

   // define the binning
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

   }

   // here is what I denote as binning option #1
   std::vector<double> e1e2BinsUnfolded = {0.0, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.00125, 0.0015, 0.00175, 0.002, 0.00225, 0.0025, 0.00275, 0.003, 0.0035, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.3};

   // define the histograms
   TH2D *h2raw_Theta = new TH2D("r_Theta","raw_Theta",2 * BinCount, 0, 2 * BinCount , e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data() );
   TH2D *h2smeared_Theta = new TH2D("smeared_Theta","smeared_Theta", 2 * BinCount, 0, 2 * BinCount , e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data() );
   TH2D *h2true_Theta = new TH2D("true_Theta","true_Theta",2 * BinCount, 0, 2 * BinCount,e1e2BinsUnfolded.size()-1, e1e2BinsUnfolded.data() );

   h2raw_Theta->Sumw2();
   h2smeared_Theta->Sumw2();
   h2true_Theta->Sumw2();

   // set up the fake correction factor
   EffCorrFactor fakeCorrFactorTheta;
   fakeCorrFactorTheta.init("/home/hbossi/PhysicsEEJetEEC/Unfolding/20250317_Unfolding/matchingScheme2/FakeCorr.root", "theta");

   // set up the RooUnfold Response
   RooUnfoldResponse response_Theta;
   response_Theta.Setup(h2smeared_Theta,h2true_Theta);

   // if the mode is 0, we want to do the full analysis crosscheck, so we will unfold reco MC as if it were data
   int nEventsRaw = 0;
   if (mode == 0){
      TFile InputFile(InputFileNameFull.c_str());
      // create messengers for the ttrees
      ParticleTreeMessenger MReco(InputFile, "t");
      alephTrkEfficiency efficiencyCorrector;

      //------------------------------
      // loop over the reco tree
      //------------------------------
      int EntryCount = MReco.GetEntries();
      // for the case of this mode, nEvents reflects the number of MC pseudodata events
      //nEventsRaw = EntryCount;
      for(int iE = 0; iE < EntryCount; iE++){
         MReco.GetEntry(iE);
         
         if(MReco.passesSTheta < 0.5) continue; // require that the distribution passSTheta
         if(MReco.passesNTrkMin < 0.5) continue; // require at least 5 tracks
         if(MReco.passesTotalChgEnergyMin < 0.5) continue; // require that the total energy is at least 15 GeV
      
         vector<FourVector> PReco;
         vector<double> efficiencyVec;
         double Evis = 0; 
         for(int i = 0; i < MReco.nParticle; i++){
            if(MReco.charge[i] == 0) continue;
            if(MReco.highPurity[i] == false) continue;
            // place cut on the reco energy, not included at gen level
            if(MReco.P[i][0] < 0.2) continue;
            double Efficiency = efficiencyCorrector.efficiency(MReco.P[i].GetTheta(), MReco.P[i].GetPhi(), MReco.P[i].GetPT(), MReco.nChargedHadronsHP);
            PReco.push_back(MReco.P[i]);
            Evis = Evis + MReco.P[i][0];
            efficiencyVec.push_back(Efficiency);
         }
         if(Evis > 200) continue; 
         nEventsRaw++;

         // now fill the unmatched tree
         //std::cout << "Preco.size " << PReco.size() << std::endl;
         for(int i = 0; i < PReco.size(); i++){
            for(int j = i+1; j < PReco.size();j++){
               if(i == j) continue; // don't include self correlations in the EEC
               // fill the tree
               FourVector Reco1 = PReco.at(i);
               FourVector Reco2 = PReco.at(j);
               double recoTheta = GetAngle(Reco1,Reco2);
               int BinThetaReco = FindBin(recoTheta, 2 * BinCount, Bins);
               double recoEEC   = Reco1[0]*Reco2[0]/(TotalE*TotalE);
               double eff1 = efficiencyVec.at(i);
               double eff2 = efficiencyVec.at(j);
               double fakeFactor = fakeCorrFactorTheta.efficiency(BinThetaReco, recoEEC);
               h2raw_Theta->Fill(BinThetaReco, recoEEC, fakeFactor);
            }
         }
      } // end loop over the reco tree
   }// end the mode 0 for the full analysis crosscheck


   // if the mode is 3 this is the only time we need to to loop over the data to fill the raw histogram
   if (mode == 3){
      TFile *inputsmeared =TFile::Open(InputFileNameData.c_str());
      TTree *smeared=(TTree*)inputsmeared->Get("UnmatchedPairTree");
      Int_t nEvData=smeared->GetEntries();
      nEventsRaw = nEvData;
      Double_t e1e2data[MAXPAIR], thetaData[MAXPAIR];
      double eff1[MAXPAIR], eff2[MAXPAIR], recoE1Data[MAXPAIR], recoE2Data[MAXPAIR];
      int nPairsData;
      smeared->SetBranchAddress("NUnmatchedPair",&nPairsData);
      smeared->SetBranchAddress("E1E2RecoUnmatched", &e1e2data);
      smeared->SetBranchAddress("DistanceUnmatchedReco", &thetaData);
      smeared->SetBranchAddress("RecoE1Unmatched", &recoE1Data);
      smeared->SetBranchAddress("RecoE2Unmatched", &recoE2Data);
      smeared->SetBranchAddress("RecoEfficiency1", &eff1);
      smeared->SetBranchAddress("RecoEfficiency2", &eff2);
      std::cout << "Number of entries in the data tree: " << nEvData << std::endl;
      for(int iEntry=0; iEntry< nEvData; iEntry++){
         smeared->GetEntry(iEntry);
         for(int i=0; i<nPairsData; i++){
            int BinTheta = FindBin(thetaData[i], 2 * BinCount, Bins);
            double z = (1-cos(thetaData[i]))/2;
            int BinZ = FindBin(z, 2*BinCount, zBins);

            if(eff1[i] > 1.0)eff1[i] = 1.0;
            if(eff2[i] > 1.0)eff2[i] = 1.0;
            double trackingEff = 1.0/(eff1[i]*eff2[i]);
            if(trackingEff < 1.0) std::cout << "Unexpected value of the tracking efficiency " << trackingEff << std::endl;


            //double fakeCorrZ = fakeCorrFactorZ.efficiency(BinZ, e1e2data[i]);
            double fakeCorrTheta = fakeCorrFactorTheta.efficiency(BinTheta, e1e2data[i]);
            h2raw_Theta->Fill(BinTheta,e1e2data[i], fakeCorrTheta*trackingEff);
            // h2raw_Z->Fill(BinZ,e1e2data[i], fakeCorrZ*trackingEff);
         }
      }
   } // end mode 3 loop over the data file

   // no matter the mode, we need to loop over the MC
   double e1e2recoMC[MAXPAIR], e1e2gen[MAXPAIR], thetaRecoMC[MAXPAIR], thetaGen[MAXPAIR];
   double recoE1[MAXPAIR], recoE2[MAXPAIR], genE1[MAXPAIR], genE2[MAXPAIR], genE[MAXPAIR];
   double recoEfficiency1[MAXPAIR], recoEfficiency2[MAXPAIR];
   int nPairsMC;


   TFile *inputmc =TFile::Open(InputFileNameMatched.c_str());
   TTree *mc=(TTree*)inputmc->Get("PairTree");
   Int_t nEv2=mc->GetEntries();
   if(mode ==1 || mode == 2) nEventsRaw = nEv2;
   std::cout << "nEvents in the mc " << nEv2 << std::endl;
   //------------------------------------------------
   mc->SetBranchAddress("NPair",&nPairsMC);
   mc->SetBranchAddress("E1E2Reco", &e1e2recoMC);
   mc->SetBranchAddress("E1E2Gen", &e1e2gen);
   mc->SetBranchAddress("DistanceReco", &thetaRecoMC);
   mc->SetBranchAddress("DistanceGen", &thetaGen);
   mc->SetBranchAddress("RecoE1", &recoE1);
   mc->SetBranchAddress("RecoE2", &recoE2);
   mc->SetBranchAddress("RecoEfficiency1", &recoEfficiency1);
   mc->SetBranchAddress("RecoEfficiency2", &recoEfficiency2);

   Int_t countm=0;
   int countLargeReweight = 0;
   for(int iEntry=0; iEntry< nEv2; iEntry++){
      mc->GetEntry(iEntry);
      for(int i=0; i<nPairsMC; i++){

         if(recoE1[i] < 0 || recoE2[i] < 0) continue; // skip over the unmatched pairs

         int BinThetaMeasured = FindBin(thetaRecoMC[i], 2 * BinCount, Bins);
         int BinThetaGenMC = FindBin(thetaGen[i], 2 * BinCount, Bins);

         double zMeasuredMC = (1-cos(thetaRecoMC[i]))/2;
         double zGenMC = (1-cos(thetaGen[i]))/2;
         int BinZMeasured = FindBin(zMeasuredMC, 2 * BinCount, zBins);
         int BinZGen = FindBin(zGenMC, 2*BinCount,zBins);

         // for trivial or split closure test also fill the raw histogram from MC
         if(mode == 1 || mode == 2)h2raw_Theta->Fill(BinThetaMeasured,e1e2recoMC[i]);

         // fill the true distributions
         h2true_Theta->Fill(BinThetaGenMC, e1e2gen[i]);

         // fill the smeared distributions
         h2smeared_Theta->Fill(BinThetaMeasured,e1e2recoMC[i]);

         double totalTrackingEff = recoEfficiency1[i]*recoEfficiency2[i];

         // now fill the response with the tracking efficiency weight
         response_Theta.Fill(BinThetaMeasured, e1e2recoMC[i],BinThetaGenMC,e1e2gen[i]);

      }
   } // end loop over the mc events


   TFile *fout = new TFile(Form("unfoldingCheckE2C_Mode%d_IterLow%d_IterHigh%d_%s.root",mode, iterLow, iterHigh, date.c_str()),"RECREATE");
   fout->cd();
   h2raw_Theta->Write();
   h2smeared_Theta->Write();
   h2true_Theta->Write();


   // create a vector for the unfolded histograms
   std::vector<TH1D> hists;
   std::vector<std::string> labels;
   TH1D* MCtrue = (TH1D*)Project2Dto1D(h2true_Theta);

   DivideByBin(*MCtrue, Bins);
   MCtrue->Scale(1.0/nEv2);
   MCtrue->SetMarkerStyle(20);
   MCtrue->SetMarkerColor(Colors[0]);
   MCtrue->SetLineColor(Colors[0]);
   MCtrue->SetLineWidth(2);

   hists.push_back(*MCtrue);
   labels.push_back("MC True");


   // now do the actual unfolding
   for(int iter = iterLow; iter < iterHigh; iter++){

      std::cout << "Unfolding for theta iter " << iter << std::endl;
      RooUnfoldBayes  unfold_Theta(&response_Theta, h2raw_Theta, iter);
      TH2D* hunf_Theta =  dynamic_cast<TH2D*>(unfold_Theta.Hreco());
      TH1* hfold_Theta = response_Theta.ApplyToTruth(hunf_Theta, "");

      TH2D *htempUnf_Theta=(TH2D*)hunf_Theta->Clone("htempUnf_Theta");
      htempUnf_Theta->SetName(Form("Bayesian_Unfoldediter%d_Theta",iter));

      TH2D *htempFold_Theta=(TH2D*)hfold_Theta->Clone("htempFold_Theta");
      htempFold_Theta->SetName(Form("Bayesian_Foldediter%d_Theta",iter));

      htempUnf_Theta->Write();
      htempFold_Theta->Write();

      TH1D* hUnfolded =  Project2Dto1D(htempUnf_Theta);
      hUnfolded->SetName(Form("Bayesian_Unfolded1Diter%d_Theta",iter));
      DivideByBin(*hUnfolded, Bins);
      hUnfolded->Scale(1.0/nEventsRaw);
      // probably need to also add correction for the matching efficiency here??
      hUnfolded->SetMarkerStyle(20);
      hUnfolded->SetMarkerColor(Colors[iter]);
      hUnfolded->SetLineColor(Colors[iter]);
      hUnfolded->SetLineWidth(2);
      hists.push_back(*hUnfolded);
      labels.push_back(Form("Iteration %d", iter));
  }
  fout->Close();
  
  

  MakeCanvas(hists, labels, "UnfoldedCheck",  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 1, true, true);



} // end the unfolding check





// ----------------------------------------------
// matchingEfficiencyClosureCheck
// ----------------------------------------------
void matchingEfficiencyClosureCheck(std::string InputFileNameFull, std::string GenTreeName, std::string InputFileNameMatched){

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

   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);

   }

  // here is what I denote as binning option #1
  std::vector<double> e1e2BinsUnfolded = {0.0, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.00125, 0.0015, 0.00175, 0.002, 0.00225, 0.0025, 0.00275, 0.003, 0.0035, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.3};

   //------------------------------------

   // ---------------------------------------------------------------------
   // Handle the case of the unmatched distributions
   // ---------------------------------------------------------------------
   std::cout << "Reading from " << InputFileNameFull.c_str() << std::endl;
   TFile InputFile(InputFileNameFull.c_str());

   // create messengers for the ttrees
   ParticleTreeMessenger MGen(InputFile, GenTreeName);
   ParticleTreeMessenger MReco(InputFile, "t");
   TH1D HEEC2_gen("HEEC2_gen", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // gen level EEC


   //------------------------------
   // loop over the gen tree
   //------------------------------
   int EntryCount = MGen.GetEntries();
   int numAcceptedEvents = 0; 
   for(int iE = 0; iE < EntryCount; iE++){
      MGen.GetEntry(iE);
      MReco.GetEntry(iE); 
      
      if(MReco.passesSTheta < 0.5) continue; // require that the distribution passSTheta
      if(MReco.passesNTrkMin < 0.5) continue; // require at least 5 tracks
      if(MReco.passesTotalChgEnergyMin < 0.5) continue; // require that the total energy is at least 15 GeV
      
      numAcceptedEvents++; 
      
      vector<FourVector> PGen;
      for(int i = 0; i < MGen.nParticle; i++){
         if(MGen.charge[i] == 0) continue;
         if(MGen.highPurity[i] == false) continue;
         // place cut on the reco energy, not included at gen level
         PGen.push_back(MGen.P[i]);
      }


      for(int i = 0; i < PGen.size(); i++){
         for(int j = i+1; j < PGen.size();j++){
            if(i == j) continue; // don't fill the EEC with self correlations particles with themselves
            FourVector Gen1 = PGen.at(i);
            FourVector Gen2 = PGen.at(j);
            double genTheta = GetAngle(Gen1,Gen2);
            int BinThetaGen = FindBin(genTheta, 2 * BinCount, Bins);
            double  genEEC  = Gen1[0]*Gen2[0]/(TotalE*TotalE);
            HEEC2_gen.Fill(BinThetaGen, genEEC);

         }
      }
    } // end loop over the reco tree



   // now do all the plotting
   // set the color and drawing style
   HEEC2_gen.SetMarkerColor(Colors[1]);
   HEEC2_gen.SetMarkerStyle(20);
   HEEC2_gen.SetLineColor(Colors[1]);
   HEEC2_gen.SetLineWidth(2);




   // divide by the bin width
   DivideByBin(HEEC2_gen, Bins);

   // scale by the number of events
   HEEC2_gen.Scale(1.0/numAcceptedEvents);

   // ---------------------------------------------------------------------
   // Handle the case of the matched distributions
   // ---------------------------------------------------------------------
   TH1D HEEC2_gen_matched("HEEC2_gen_matched", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount);


   double e1e2recoMC[MAXPAIR], e1e2gen[MAXPAIR], thetaRecoMC[MAXPAIR], thetaGen[MAXPAIR];
   double recoE1[MAXPAIR], recoE2[MAXPAIR], genE1[MAXPAIR], genE2[MAXPAIR], genE[MAXPAIR];
   double recoEfficiency1[MAXPAIR], recoEfficiency2[MAXPAIR];
   int nPairsMC;


   TFile *inputmc =TFile::Open(InputFileNameMatched.c_str());
   TTree *mc=(TTree*)inputmc->Get("PairTree");


   Int_t nEv2=mc->GetEntries();
   std::cout << "nEvents in the mc " << nEv2 << std::endl;
   //------------------------------------------------
   mc->SetBranchAddress("NPair",&nPairsMC);
   mc->SetBranchAddress("E1E2Gen", &e1e2gen);
   mc->SetBranchAddress("DistanceGen", &thetaGen);
   mc->SetBranchAddress("RecoE1", &recoE1);
   mc->SetBranchAddress("RecoE2", &recoE2);

   for(int iEntry=0; iEntry< nEv2; iEntry++){
      mc->GetEntry(iEntry);
      for(int i=0; i<nPairsMC; i++){

      if(recoE1[i] < 0 || recoE2[i] < 0) continue; // skip over the unmatched pairs still in MC

      int BinThetaGen = FindBin(thetaGen[i], 2 * BinCount, Bins);
      //std::cout << "Bin Theta Measured: " << BinThetaMeasured << " EEC: " << e1e2recoMC[i] << std::endl;
      HEEC2_gen_matched.Fill(BinThetaGen,e1e2gen[i]);


      }
   } // end loop over the mc events

   HEEC2_gen_matched.SetMarkerColor(Colors[5]);
   HEEC2_gen_matched.SetMarkerStyle(20);
   HEEC2_gen_matched.SetLineColor(Colors[5]);
   HEEC2_gen_matched.SetLineWidth(2);
   std::cout << HEEC2_gen_matched.Integral() << std::endl;
   DivideByBin(HEEC2_gen_matched, Bins);
   HEEC2_gen_matched.Scale(1.0/nEv2);
   std::cout << HEEC2_gen_matched.Integral() << std::endl;

   std::vector<TH1D> hists = {HEEC2_gen, HEEC2_gen_matched};
   MakeCanvas(hists, {"Archived MC Gen", "Archived MC Gen Matched"}, "MatchingEfficiencyCheck",  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 1, true, true);

   // close the input files
   InputFile.Close();
   inputmc->Close();


}
// ----------------------------------------------
// fakeCorrectionClosureCheck
// ----------------------------------------------
void fakeCorrectionClosureCheck(std::string InputFileNameFull, std::string RecoTreeName, std::string InputFileNameMatched){

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

   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);

   }

  // here is what I denote as binning option #1
  std::vector<double> e1e2BinsUnfolded = {0.0, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.00125, 0.0015, 0.00175, 0.002, 0.00225, 0.0025, 0.00275, 0.003, 0.0035, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.3};

   //------------------------------------

   // ---------------------------------------------------------------------
   // Handle the case of the unmatched distributions
   // ---------------------------------------------------------------------
   std::cout << "Reading from " << InputFileNameFull.c_str() << std::endl;
   TFile InputFile(InputFileNameFull.c_str());

   // create messengers for the ttrees
   ParticleTreeMessenger MGen(InputFile, "tgen");
   ParticleTreeMessenger MReco(InputFile, RecoTreeName);
   TH1D HEEC2_reco("HEEC2_reco", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // reco level EEC no corrections with pT selection

   //------------------------------
   // loop over the reco tree
   //------------------------------
   int EntryCount = MReco.GetEntries();
   int EntryCountMC = MGen.GetEntries();
   std::cout << "Initial Entry count in the MC reco is " << EntryCount << " as compared to the MC gen is " << EntryCountMC << std::endl;
   int numAcceptedEvents = 0; 
   for(int iE = 0; iE < EntryCount; iE++){
      MReco.GetEntry(iE);
      MGen.GetEntry(iE);
      
      // event selection 
      //if(MGen.passesSTheta < 0.5) continue; // require that the distribution passSTheta
      if(MReco.passesSTheta < 0.5) continue; // require that the distribution passSTheta
      if(MReco.passesNTrkMin < 0.5) continue; // require at least 5 tracks
      if(MReco.passesTotalChgEnergyMin < 0.5) continue; // require that the total energy is at least 15 GeV
     // if(MGen.passesSTheta != MReco.passesSTheta) std::cout << "MGen.passesStheta " << MGen.passesSTheta << " and  MReco.passesStheta " << MReco.passesSTheta << " do not agree. " << std::endl;

   
      vector<FourVector> PReco;
      vector<double> efficiencyVec;
      double EvisReco = 0; 
      for(int i = 0; i < MReco.nParticle; i++){
         if(MReco.charge[i] == 0) continue;
         if(MReco.highPurity[i] == false) continue;
         if(abs(cos(MReco.theta[i])) > 0.94) continue; 
         // place cut on the reco energy, not included at gen level
         if(MReco.P[i][0] < 0.2) continue;
         EvisReco = EvisReco + MReco.P[i][0]; 
         PReco.push_back(MReco.P[i]);
      }
      
      if(EvisReco > 200) continue; 
      numAcceptedEvents++; 

      // now fill the unmatched tree
      //std::cout << "Preco.size " << PReco.size() << std::endl;
      for(int i = 0; i < PReco.size(); i++){
         for(int j = i+1; j < PReco.size();j++){
            if(i == j) continue; // don't include self correlations in the EEC
            // fill the tree
            FourVector Reco1 = PReco.at(i);
            FourVector Reco2 = PReco.at(j);
            double recoTheta = GetAngle(Reco1,Reco2);
            int BinThetaReco = FindBin(recoTheta, 2 * BinCount, Bins);
            double recoEEC   = Reco1[0]*Reco2[0]/(TotalE*TotalE);
            HEEC2_reco.Fill(BinThetaReco, recoEEC);
         }
      }
    } // end loop over the reco tree



   // now do all the plotting
   // set the color and drawing style
   HEEC2_reco.SetMarkerColor(Colors[0]);
   HEEC2_reco.SetMarkerStyle(20);
   HEEC2_reco.SetLineColor(Colors[0]);
   HEEC2_reco.SetLineWidth(2);



   // divide by the bin width
   DivideByBin(HEEC2_reco, Bins);

   // scale by the number of events
   HEEC2_reco.Scale(1.0/numAcceptedEvents);

   // ---------------------------------------------------------------------
   // Handle the case of the matched distributions
   // ---------------------------------------------------------------------
   TH1D HEEC2_reco_matched("HEEC2_reco_matched", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount);


   double e1e2recoMC[MAXPAIR], e1e2gen[MAXPAIR], thetaRecoMC[MAXPAIR], thetaGen[MAXPAIR];
   double recoE1[MAXPAIR], recoE2[MAXPAIR], genE1[MAXPAIR], genE2[MAXPAIR], genE[MAXPAIR];
   double recoEfficiency1[MAXPAIR], recoEfficiency2[MAXPAIR];
   int nPairsMC;


   TFile *inputmc =TFile::Open(InputFileNameMatched.c_str());
   TTree *mc=(TTree*)inputmc->Get("PairTree");


   Int_t nEv2=mc->GetEntries();
   std::cout << "nEvents in the mc " << nEv2 << std::endl;
   //------------------------------------------------
   mc->SetBranchAddress("NPair",&nPairsMC);
   mc->SetBranchAddress("E1E2Reco", &e1e2recoMC);
   mc->SetBranchAddress("DistanceReco", &thetaRecoMC);
   mc->SetBranchAddress("RecoE1", &recoE1);
   mc->SetBranchAddress("RecoE2", &recoE2);

   for(int iEntry=0; iEntry< nEv2; iEntry++){
      mc->GetEntry(iEntry);
      for(int i=0; i<nPairsMC; i++){

      if(recoE1[i] < 0 || recoE2[i] < 0) continue; // skip over the unmatched pairs

      int BinThetaMeasured = FindBin(thetaRecoMC[i], 2 * BinCount, Bins);
      //std::cout << "Bin Theta Measured: " << BinThetaMeasured << " EEC: " << e1e2recoMC[i] << std::endl;
      HEEC2_reco_matched.Fill(BinThetaMeasured,e1e2recoMC[i]);


      }
   } // end loop over the mc events

   HEEC2_reco_matched.SetMarkerColor(Colors[5]);
   HEEC2_reco_matched.SetMarkerStyle(20);
   HEEC2_reco_matched.SetLineColor(Colors[5]);
   HEEC2_reco_matched.SetLineWidth(2);
   DivideByBin(HEEC2_reco_matched, Bins);
   HEEC2_reco_matched.Scale(1.0/nEv2);

   std::vector<TH1D> hists = {HEEC2_reco, HEEC2_reco_matched};
   MakeCanvas(hists, {"Archived MC Reco", "Archived MC Reco Matched"}, "FakeFractionCheck",  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 1, true, true);

   // close the input files
   InputFile.Close();
   inputmc->Close();


}



// ----------------------------------------------
// checks before the matching occurs
// ----------------------------------------------
void preMatchingClosureCheck(std::string InputFileName, std::string GenTreeName, std::string RecoTreeName){

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

   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);

   }

  // here is what I denote as binning option #1
  std::vector<double> e1e2BinsUnfolded = {0.0, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.00125, 0.0015, 0.00175, 0.002, 0.00225, 0.0025, 0.00275, 0.003, 0.0035, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.3};



   std::cout << "Reading from " << InputFileName.c_str() << std::endl;
   TFile InputFile(InputFileName.c_str());

   // create messengers for the ttrees
   ParticleTreeMessenger MGen(InputFile, GenTreeName);
   ParticleTreeMessenger MReco(InputFile, RecoTreeName);

   TH1D HEEC2_gen("HEEC2_gen", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // gen level EEC after event selections
   TH1D HEEC2_reco("HEEC2_reco", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // reco level EEC no corrections with pT selection
   TH1D HEEC2_recoWithTrackEff("HEEC2_recoWithTrackEff", ";EEC_{2};", 2 * BinCount, 0, 2 * BinCount); // reco level EEC with pT selection with tracking efficiency

   alephTrkEfficiency efficiencyCorrector;
   //------------------------------
   // loop over the reco tree
   //------------------------------
   int EntryCount = MReco.GetEntries();
   int nAcceptedRecoEvents = 0; 
   for(int iE = 0; iE < EntryCount; iE++){
      MReco.GetEntry(iE);
      
      if(MReco.passesSTheta < 0.5) continue; // require that the distribution passSTheta
      if(MReco.passesNTrkMin < 0.5) continue; // require at least 5 tracks
      if(MReco.passesTotalChgEnergyMin < 0.5) continue; // require that the total energy is at least 15 GeV
   
      
      vector<FourVector> PReco;
      vector<double> efficiencyVec;
      double EvisReco = 0; 
      for(int i = 0; i < MReco.nParticle; i++){
         if(MReco.charge[i] == 0) continue;
         if(MReco.highPurity[i] == false) continue;
         if(abs(cos(MReco.theta[i])) > 0.94) continue; 
         // place cut on the reco energy, not included at gen level
         if(MReco.P[i][0] < 0.2) continue;
         EvisReco = EvisReco + MReco.P[i][0]; 
         double Efficiency = efficiencyCorrector.efficiency(MReco.P[i].GetTheta(), MReco.P[i].GetPhi(), MReco.P[i].GetPT(), MReco.nChargedHadronsHP);
         PReco.push_back(MReco.P[i]);
         efficiencyVec.push_back(Efficiency);
      }
      
      // reject events with reco Evis > 200 to remove laser calibration runs
      if(EvisReco > 200) continue; 
      
      nAcceptedRecoEvents++; 


      // now fill the unmatched tree
      //std::cout << "Preco.size " << PReco.size() << std::endl;
      for(int i = 0; i < PReco.size(); i++){
         for(int j = i+1; j < PReco.size();j++){
            if(i == j) continue; // don't include self correlations in the EEC
            // fill the tree
            FourVector Reco1 = PReco.at(i);
            FourVector Reco2 = PReco.at(j);
            double recoTheta = GetAngle(Reco1,Reco2);
            int BinThetaReco = FindBin(recoTheta, 2 * BinCount, Bins);
            double recoEEC   = Reco1[0]*Reco2[0]/(TotalE*TotalE);
            double eff1 = efficiencyVec.at(i);
            double eff2 = efficiencyVec.at(j);

            HEEC2_reco.Fill(BinThetaReco, recoEEC);
            HEEC2_recoWithTrackEff.Fill(BinThetaReco, recoEEC/(eff1*eff2));
         }
      }
    } // end loop over the reco tree




   //------------------------------
   // loop over the gen tree
   //------------------------------
   int EntryCountGen = MGen.GetEntries();
   int nAcceptedGenEvents = 0; 
   for(int iE = 0; iE < EntryCountGen; iE++){
      MGen.GetEntry(iE);
      
      if(MGen.passesSTheta < 0.5) continue; // require that the distribution passSTheta
      
      double EvisGen = 0; 
      vector<FourVector> PGen;
      for(int i = 0; i < MGen.nParticle; i++){
         if(MGen.charge[i] == 0) continue;
         if(MGen.highPurity[i] == false) continue;
         if(abs(cos(MGen.theta[i])) > 0.94) continue; 
         EvisGen = EvisGen + MGen.P[i][0];
         PGen.push_back(MGen.P[i]);
      }
      
      // trivial cut to remove laser calibrations at the Gen Level, should effectively do nothing.
      if(EvisGen > 200) continue; 
      nAcceptedGenEvents++; 


      for(int i = 0; i < PGen.size(); i++){
         for(int j = i+1; j < PGen.size();j++){
            if(i == j) continue; // don't fill the EEC with self correlations particles with themselves
            FourVector Gen1 = PGen.at(i);
            FourVector Gen2 = PGen.at(j);
            double genTheta = GetAngle(Gen1,Gen2);
            int BinThetaGen = FindBin(genTheta, 2 * BinCount, Bins);
            double  genEEC  = Gen1[0]*Gen2[0]/(TotalE*TotalE);
            HEEC2_gen.Fill(BinThetaGen, genEEC);

         }
      }
    } // end loop over the gen tree

   // now do all the plotting
   // set the color and drawing style
   HEEC2_reco.SetMarkerColor(Colors[0]);
   HEEC2_reco.SetMarkerStyle(20);
   HEEC2_reco.SetLineColor(Colors[0]);
   HEEC2_reco.SetLineWidth(2);

   HEEC2_recoWithTrackEff.SetMarkerColor(Colors[3]);
   HEEC2_recoWithTrackEff.SetMarkerStyle(20);
   HEEC2_recoWithTrackEff.SetLineColor(Colors[3]);
   HEEC2_recoWithTrackEff.SetLineWidth(2);

   HEEC2_gen.SetMarkerColor(Colors[1]);
   HEEC2_gen.SetMarkerStyle(20);
   HEEC2_gen.SetLineColor(Colors[1]);
   HEEC2_gen.SetLineWidth(2);

   // divide by the bin width
   DivideByBin(HEEC2_reco, Bins);
   DivideByBin(HEEC2_recoWithTrackEff, Bins);
   DivideByBin(HEEC2_gen, Bins);

   // scale by the number of events
   HEEC2_reco.Scale(1.0/nAcceptedRecoEvents);
   HEEC2_gen.Scale(1.0/nAcceptedGenEvents);
   HEEC2_recoWithTrackEff.Scale(1.0/nAcceptedRecoEvents);

   std::vector<TH1D> hists_baseline = {HEEC2_gen, HEEC2_reco};
   MakeCanvas(hists_baseline, {"Archived MC - Gen", "Archived MC Reco"}, "PreMatchingCheck_NewEventSel",  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 1, true, true);


   std::vector<TH1D> hists = {HEEC2_gen, HEEC2_recoWithTrackEff};
   MakeCanvas(hists, {"Archived MC - Gen", "Archived MC Reco w/ Tracking Eff"}, "PreMatchingCheck_withTrackEff_NewEventSel",  "#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-3, 1, true, true);


   // -------------------------------------------------------------------

   // close the input files
   InputFile.Close();

}


//==============================================================================
// Plotting for Z and Theta
//==============================================================================

void MakeCanvasPointerVec(vector<TH1D *> Histograms, vector<string> Labels, string Output,
   string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0]->GetNbinsX();

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
   // Canvas.SetLogy();
   // Canvas.SetRightMargin(MarginR);
   // Canvas.SetLeftMargin(MarginL);
   // Canvas.SetTopMargin(MarginT);
   // Canvas.SetBottomMargin(MarginB);

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
   for(TH1D *H : Histograms){
        H->Draw("E same");
        std::cout << "On histogram " << H->GetName() << " which has an x axis range of " << H->GetXaxis()->GetXmin() << " to " << H->GetXaxis()->GetXmax() << std::endl;
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

   double WorldRMin = 0.6;
   double WorldRMax = 1.4;

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
         TH1D *H = (TH1D *)Histograms[i]->Clone();
         H->Divide(Histograms[0]);
         H->Draw("E same");
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
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.015, "#theta_{L} = #pi/2");

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
   Latex.DrawLatex(0.01, 0.01, "Work-in-progress, 2024 December 31st, HB");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.035 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++){
      Legend.AddEntry(Histograms[i], Labels[i].c_str(), "pl");
   }
   Legend.Draw();

   TLegend Legend2(0.7, 0.90, 0.9, 0.90 - 0.035 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);

   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(Histograms[i], Labels[i].c_str(), "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}


void MakeCanvas(vector<TH1D> Histograms, vector<string> Labels, string Output,
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

   double WorldRMin = 0.9;
   double WorldRMax = 1.1;

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
         H->Draw(" ep same");
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
   // axis on the x axis on the left hand side 
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, 510, "GS");
   // axis on the x axis on the right hand side 
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, 510, "-GS");
   // axis on the x axis on the left hand side for the ratio
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   // axis on the x axis on the right hand side for the ratio
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
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "Full Analysis Closure Studies - April 24th");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.06 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   //Legend.AddEntry(&jingyuGraph, "Crosscheck Results"); 
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


void MakeCanvasZ(vector<TH1D> Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX){


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

   double WorldRMin = 0.6;
   double WorldRMax = 1.4;

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
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.015, "#it{z} = 1/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.9, MarginB * 0.4, X.c_str());

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
   Latex.DrawLatex(0.01, 0.01, "Work-in-progress, 2024 August 7th HB");

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




