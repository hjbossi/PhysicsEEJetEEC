/* plotFromUnfolding.cpp: Macro to plot the final result from the unfolded result.
* Hannah Bossi, <hannah.bossi@cern.ch>
* 05/28/2024
 */

#include <iostream>
using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TPad.h"
#include "TGraph.h"

#include "SetStyle.h"
#include "ProgressBar.h"
#include "CommandLine.h"
#include "Messenger.h"
#include "JetCorrector.h"
#include "alephTrkEfficiency.h"


void MakeCanvas(vector<TH1D *> Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX);
void SetPad(TPad &P);
void DivideByBin(TH1D &H, double Bins[]);
int main(int argc, char *argv[]);

int main(int argc, char *argv[]){

    // take the parameters from the command line
    CommandLine CL(argc, argv);
    vector<string> FileNames = CL.GetStringVector("Input");
    string output = CL.Get("Output", "plots");
    vector<string> Labels = CL.GetStringVector("Label");
    string Prefix = CL.Get("Prefix", "");
    bool DoRatio = CL.GetBool("DoRatio", true);
    int iterNom = CL.GetInt("Iter", 3);

    // SetThesisStyle();
    static vector<int> Colors = GetCVDColors6();


   const int BinCount = 100;
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;
   for(int i = 0; i <= BinCount; i++)
   {
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
   }


    int nEvents = 592388;
    double TotalE = 91.1876; // GeV


    int N = FileNames.size();
    vector<TFile *> Files(N);
    vector<TH1D *> Histograms;
    for(int i = 0; i < N; i++){
      Files[i] = new TFile(FileNames[i].c_str());
    }

    int nIters = 4;

    for(int i = 0; i < N; i++){
        std::cout << "On file " << i << std::endl;
        TH2D *hUnfolded = (TH2D *)Files[i]->Get(Form("Bayesian_Unfoldediter%d_Theta", iterNom));
        if(hUnfolded == NULL){
          cerr << "Histogram not found" << endl;
          return 0;
        }

        TH1D *H = (TH1D*)hUnfolded->ProjectionX();
        H->Reset();
        for (int i = 1; i <= hUnfolded->GetNbinsX(); ++i) {
            double weight = 0;
            double error = 0;
            for (int j = 1; j <= hUnfolded->GetNbinsY(); ++j) {
               double binContent = hUnfolded->GetBinContent(i, j);
               double binError= hUnfolded->GetBinError(i,j); // manully enforce poisson errors //hUnfolded->GetBinError(i,j);
               double binCenter = hUnfolded->GetYaxis()->GetBinCenter(j);
               weight += binContent*((binCenter));
                error += pow(binError*binCenter, 2);

            }
            H->SetBinContent(i, weight);
            H->SetBinError(i, sqrt(error));
        }
        DivideByBin(*H, Bins);
        H->Scale(1.0 / H->Integral());
        H->SetDirectory(0);
        H->SetStats(0);
        H->SetTitle("");
        H->GetXaxis()->SetTitle("#theta_{L}");
        H->GetYaxis()->SetTitle("#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}");
        H->SetMarkerColor(Colors[i]);
        H->SetMarkerStyle(20);
        H->SetLineColor(Colors[i]);
        H->SetLineWidth(2);

      //   TH1D* H2 = (TH1D*)Files[i]->Get("h1MCGen_Theta");
      //   DivideByBin(*H2, Bins);
      //   H2->Scale(1.0 / H2->Integral());
      //   H2->SetDirectory(0);
      //   H2->SetStats(0);
      //   H2->SetTitle("");
      //   H2->GetXaxis()->SetTitle("#theta_{L}");
      //   H2->GetYaxis()->SetTitle("#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}");
      //   H2->SetMarkerColor(Colors[1]);
      //   H2->SetMarkerStyle(20);
      //   H2->SetLineColor(Colors[1]);
      //   H2->SetLineWidth(2);

      //   Histograms.push_back(H2);
        Histograms.push_back(H);

    }


    double YMax = DoRatio ? 1 : 1;
   std::cout << "Plotting the distribution" << std::endl;
    MakeCanvas(Histograms, Labels, output + Prefix + "EEC2","#theta_{L}", "#frac{1}{N_{event}} #frac{d(Sum E_{i}E_{j}/E^{2})}{d #theta_{L}}", 2e-4, 0.05,
         DoRatio, true);


    for(int i = 0; i < N; i++){
        Files[i]->Close();
        delete Files[i];
    }

    return 0;

}


void MakeCanvas(vector<TH1D *> Histograms, vector<string> Labels, string Output,
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

