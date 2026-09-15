#include <iostream>
#include <vector>
#include <string>
#include <map>
using namespace std;

#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "Messenger.h"
#include "CommandLine.h"
#include "Matching.h"
#include "ProgressBar.h"
#include "TauHelperFunctions3.h"
#include "alephTrkEfficiency.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLatex.h"

#define MAX 1000
#define MAXPAIR 10000

int main(int argc, char *argv[]);

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName  = CL.Get("Input");
   string GenTreeName    = CL.Get("Gen", "tgen");
   string RecoTreeName   = CL.Get("Reco", "t");
   string OutputFileName = CL.Get("Output");
   double Fraction       = CL.GetDouble("Fraction", 1.00);

    std::cout << "Reading from " << InputFileName.c_str() << std::endl;
   TFile InputFile(InputFileName.c_str());
   TFile OutputFile(OutputFileName.c_str(), "RECREATE");


   //------------------------------------
   // define the trees
   //------------------------------------
   // tree that just takes pairs, not matching taken into account
   TTree OutputUnmatchedTree("UnmatchedPairTree", "");



     // variables for the unmatched pair tree
   int NUnmatchedPair;
   double GenE1Unmatched[MAXPAIR], GenX1Unmatched[MAXPAIR], GenY1Unmatched[MAXPAIR], GenZ1Unmatched[MAXPAIR];
   double GenE2Unmatched[MAXPAIR], GenX2Unmatched[MAXPAIR], GenY2Unmatched[MAXPAIR], GenZ2Unmatched[MAXPAIR];
   double RecoE1Unmatched[MAXPAIR], RecoX1Unmatched[MAXPAIR], RecoY1Unmatched[MAXPAIR], RecoZ1Unmatched[MAXPAIR];
   double RecoE2Unmatched[MAXPAIR], RecoX2Unmatched[MAXPAIR], RecoY2Unmatched[MAXPAIR], RecoZ2Unmatched[MAXPAIR];
   double DistanceUnmatchedGen[MAXPAIR], DistanceUnmatchedReco[MAXPAIR], DeltaPhiUnmatchedGen[MAXPAIR], DeltaPhiUnmatchedReco[MAXPAIR], DeltaEUnmatchedReco[MAXPAIR], DeltaEUnmatchedGen[MAXPAIR], DeltaThetaUnmatchedGen[MAXPAIR], DeltaThetaUnmatchedReco[MAXPAIR];
   double E1E2GenUnmatched[MAXPAIR], E1E2RecoUnmatched[MAXPAIR];
   double RecoEfficiency1[MAXPAIR], RecoEfficiency2[MAXPAIR];
   OutputUnmatchedTree.Branch("NUnmatchedPair", &NUnmatchedPair);
   OutputUnmatchedTree.Branch("RecoE1Unmatched", &RecoE1Unmatched, "RecoE1Unmatched[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("RecoE2Unmatched", &RecoE2Unmatched, "RecoE2Unmatched[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("DistanceUnmatchedReco", &DistanceUnmatchedReco, "DistanceUnmatchedReco[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("DeltaPhiUnmatchedReco", &DeltaPhiUnmatchedReco, "DeltaPhiUnmatchedReco[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("DeltaEUnmatchedReco", &DeltaEUnmatchedReco, "DeltaEUnmatchedReco[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("DeltaThetaUnmatchedReco", &DeltaThetaUnmatchedReco, "DeltaThetaUnmatchedReco[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("E1E2RecoUnmatched", &E1E2RecoUnmatched, "E1E2RecoUnmatched[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("RecoEfficiency1", &RecoEfficiency1, "RecoEfficiency1[NUnmatchedPair]/D");
   OutputUnmatchedTree.Branch("RecoEfficiency2", &RecoEfficiency2, "RecoEfficiency2[NUnmatchedPair]/D");


   ParticleTreeMessenger MReco(InputFile, RecoTreeName);
   alephTrkEfficiency efficiencyCorrector;

   int EntryCount = MReco.GetEntries() * Fraction;
   std::cout << "Entry count " << EntryCount << std::endl;
   int nAcceptedEvents = 0;
   ProgressBar Bar(cout, EntryCount);
   Bar.SetStyle(-1);
   for(int iE = 0; iE < EntryCount; iE++)
   {

      MReco.GetEntry(iE);

      // add some explicit event selection for the hadronic events
      if(MReco.passesSTheta < 0.5) continue; // require that the distribution passSTheta
      if(MReco.passesNTrkMin < 0.5) continue; // require at least 5 tracks
      if(MReco.passesTotalChgEnergyMin < 0.5) continue; // require that the total energy is at least 15 GeV


      vector<FourVector> PReco;
      vector<double> efficiencyVec;
      double Evis = 0; 
      for(int i = 0; i < MReco.nParticle; i++){
         // apply the charged particle selection
         if(MReco.charge[i] == 0) continue;
         if(MReco.highPurity[i] == false) continue;
         if(abs(cos(MReco.theta[i])) > 0.94) continue;   
         
         // place cut on the reco energy, not included at gen level
         if(MReco.P[i][0] < 0.2) continue;
         double Efficiency = efficiencyCorrector.efficiency(MReco.P[i].GetTheta(), MReco.P[i].GetPhi(), MReco.P[i].GetPT(), MReco.nChargedHadronsHP);
         PReco.push_back(MReco.P[i]);
          Evis = Evis + MReco.P[i][0]; 
         efficiencyVec.push_back(Efficiency);
      }

      if (Evis > 200.0) continue; // skip laser calibration runs
      nAcceptedEvents++; 
      // now fill the unmatched tree
      int index_counter = 0;
      double TotalE = 91.1876;
      //std::cout << "Preco.size " << PReco.size() << std::endl;
      for(int i = 0; i < PReco.size(); i++){
         for(int j = i+1; j < PReco.size();j++){
            if(i == j) continue; // don't match particles with themselves
            // fill the tree
            FourVector Reco1 = PReco.at(i);
            FourVector Reco2 = PReco.at(j);
            RecoE1Unmatched[index_counter] = Reco1[0];
            RecoE2Unmatched[index_counter] = Reco2[0];
            DistanceUnmatchedReco[index_counter] = GetAngle(Reco1,Reco2);
            DeltaPhiUnmatchedReco[index_counter] = GetDPhi(Reco1, Reco2);
            DeltaEUnmatchedReco[index_counter] = Reco1[0]-Reco2[0];
            DeltaThetaUnmatchedReco[index_counter] = Reco1.GetTheta() - Reco2.GetTheta();
            E1E2RecoUnmatched[index_counter] =  Reco1[0]*Reco2[0]/(TotalE*TotalE);
            RecoEfficiency1[index_counter] = efficiencyVec.at(i);
            RecoEfficiency2[index_counter] = efficiencyVec.at(j);
            index_counter++;
         }
      }
      NUnmatchedPair = index_counter;

      OutputUnmatchedTree.Fill();
      

    }
    
    std::cout << "Number of accepted events " << nAcceptedEvents << std::endl;
   Bar.Update(EntryCount);
   Bar.Print();
   Bar.PrintLine();

   OutputFile.cd();

   // write the trees to the output file
   OutputUnmatchedTree.Write();
   // -------------------------------------------------------------------

   // write the output files
   OutputFile.Close();
   InputFile.Close();


   return 0;
}