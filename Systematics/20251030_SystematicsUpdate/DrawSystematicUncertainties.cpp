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

void MakeCanvasRatioOnlyZ(vector<TH1D*> Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool LogX);
void MakeCanvasRatioOnly(vector<TH1D*> Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool LogX);
void SetPad(TPad &P);
void DivideByBin(TH1D &H, double Bins[]);
int main(int argc, char *argv[]);

int main(int argc, char *argv[]){

   // take the parameters from the command line
   CommandLine CL(argc, argv);
   vector<string> FileNames = CL.GetStringVector("Input");
   string output = CL.Get("Output", "plots");
   string Variable = CL.Get("Variable", "Theta");
   vector<string> Labels = CL.GetStringVector("Label");
   string Prefix = CL.Get("Prefix", "");
   bool DoRatio = CL.GetBool("DoRatio", true);
   int iterNom = CL.GetInt("Iter", 3);

   // SetThesisStyle();
   static vector<int> Colors = GetCVDColors10();
   const int BinCount = 100;
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;
   for(int i = 0; i <= BinCount; i++)
   {
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
   }

    int N = FileNames.size();
    vector<TFile *> Files(N);
    vector<TH1D*> Histograms;
    for(int i = 0; i < N; i++){
      Files[i] = new TFile(FileNames[i].c_str());
    }

    int nIters = 1;

    for(int i = 0; i < N; i++){
         TH2D *hNom = (TH2D *)Files[i]->Get("NominalHistogram");
        if(hNom == NULL){
            cerr << "Histogram not found" << endl;
            return 0;
        }
        TH1D *H0 = (TH1D*)hNom->ProjectionX();
        H0->Reset();
        for (int h = 1; h <= hNom->GetNbinsX(); ++h) {
            double weight = 0;
            double error = 0;
            for (int k = 1; k <= hNom->GetNbinsY(); ++k) {
                double binContent = hNom->GetBinContent(h, k);
                double binError= hNom->GetBinError(h,k); // manully enforce poisson errors //hUnfolded->GetBinError(i,j);
                double binCenter = hNom->GetYaxis()->GetBinCenter(k);
                weight += binContent*((binCenter));
                error += pow(binError*binCenter, 2);
            }
            H0->SetBinContent(h, weight);
            H0->SetBinError(h, sqrt(error));
        }
        DivideByBin(*H0, Bins);
        //H0->Scale(1.0 / H->Integral());
        H0->SetDirectory(0);
        H0->SetStats(0);
        H0->SetTitle("");
        H0->SetMarkerColor(Colors[0]);
        H0->SetMarkerStyle(20);
        H0->SetLineColor(Colors[0]);
        H0->SetLineWidth(2);
        Histograms.push_back(H0);
        for(int sysIndex = 0; sysIndex < Labels.size(); sysIndex++){
            std::cout << "Processing Systematic " << Labels.at(sysIndex) << std::endl;
            TH2D *hSysNom = (TH2D *)hNom->Clone("hSysNom");
            TH2D *hSys;
            if(Labels.at(sysIndex) != "Total"){
                hSys= (TH2D *)Files[i]->Get(Form("Systematics_%s_%s_Variant", Variable.c_str(), Labels.at(sysIndex).c_str()));
               for (int x = 1; x <= hSysNom->GetNbinsX(); x++) {
                  for (int y = 1; y <= hSysNom->GetNbinsY(); y++) {
                        double value1 = hSysNom->GetBinContent(x, y);
                        double value2 = hSys->GetBinContent(x, y);
                        hSys->SetBinContent(x, y, abs(value1 - value2)); // Store the difference
                  }
               }
            }
            else{
               hSys= (TH2D *)Files[i]->Get(Form("Systematics_%s_%s", Variable.c_str(), Labels.at(sysIndex).c_str()));
            }


            std::cout << "hSys integral " << hSys->Integral() << std::endl;

            if(hSys == NULL){
                cerr << "Histogram not found" << endl;
                return 0;
            }
            TH1D *H = (TH1D*)hSys->ProjectionX(Form("hSysProj_%s", Labels.at(sysIndex)));
            H->Reset();
            for (int h = 1; h <= hSys->GetNbinsX(); ++h) {
                double weight = 0;
                double error = 0;
                for (int k = 1; k <= hSys->GetNbinsY(); ++k) {
                double binContent = hSys->GetBinContent(h, k);
                double binError= hSys->GetBinError(h,k); // manully enforce poisson errors //hUnfolded->GetBinError(i,j);
                double binCenter = hSys->GetYaxis()->GetBinCenter(k);
                weight += binContent*((binCenter));
                    error += pow(binError*binCenter, 2);

                }
                H->SetBinContent(h, weight);
                H->SetBinError(h, sqrt(error));
            }
            DivideByBin(*H, Bins);
            std::cout << "H integral " << hSys->Integral() << std::endl;
            H->SetDirectory(0);
            H->SetStats(0);
            H->SetTitle("");
            H->SetMarkerColor(Colors[sysIndex+1]);
            H->SetMarkerStyle(20);
            H->SetLineColor(Colors[sysIndex+1]);
            H->SetLineWidth(2);

            Histograms.push_back(H);
            delete hSys;
            delete hSysNom;
            //delete H;
        }// end loop over systematics
    }

   //if(Variable.c_str() == "Theta"){
   //MakeCanvasRatioOnly(Histograms, Labels, output + Prefix + "SystematicsTheta","#theta_{L}", "Fractional Relative Error", 0, 0.25, true);
   // }
   // else{
   MakeCanvasRatioOnlyZ(Histograms, Labels, output + Prefix + "Systematics_" + Variable,"#it{z} = (1- cos(#theta))/2", "Fractional Relative Error", 0, 0.25, true);
   // }

    return 0;

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


void MakeCanvasRatioOnly(vector<TH1D* > Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0]->GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 17 : 0;
   double WorldXMax = LogX ? 183: 1;

   double PadWidth = 1200;
   double PadHeight = 880;
   double PadRHeight = 0;

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

   TPad Pad("Pad", "", MarginL, MarginB , MarginL + PadWidth, MarginB + PadHeight);
   SetPad(Pad);

   Pad.cd();

   TH2D HWorld("HWorld", "", N, WorldXMin, WorldXMax, 100, WorldMin, WorldMax);
   HWorld.SetStats(0);
   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);
   HWorld.Draw("axis");
   for(int i = 1; i < NLine; i++)
   {
      std::cout << "Histograms[i]: " << Histograms[i]->Integral() << std::endl;
      TH1D *H = (TH1D *)Histograms[i]->Clone();
      H->Divide(Histograms[0]);
      H->Draw("hist p l same");
   }


   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : 1 / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : 1/ 2, 10000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");


   TGraph G2;

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
   
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadHeight, WorldMin, WorldMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadHeight, WorldMin, WorldMax, 510, "");

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
      X5.Draw();
      X6.Draw();
   }
   if(LogX == false)
   {
      XL2.Draw();
   }
   // Y1.Draw();
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
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.04, "#theta_{L} = #pi/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.9, MarginB * 0.4, X.c_str());

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(90);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "Finalization of Result April 24 (HB)");

   TLegend Legend(0.15, 0.88, 0.35, 0.88 - 0.04 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 1; i < NLine && i < 4; i++)
   {
      Legend.AddEntry(Histograms[i], Labels[i-1].c_str(),"pl");
   }
   Legend.Draw();

   TLegend Legend2(0.55, 0.88, 0.8, 0.88 - 0.04 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(Histograms[i], Labels[i-1].c_str(),"pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}


void MakeCanvasRatioOnlyZ(vector<TH1D* > Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0]->GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 17 : 0;
   double WorldXMax = LogX ? 183: 1;

   double PadWidth = 1200;
   double PadHeight = 880;
   double PadRHeight = 0;

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

   TPad Pad("Pad", "", MarginL, MarginB , MarginL + PadWidth, MarginB + PadHeight);
   SetPad(Pad);

   Pad.cd();

   TH2D HWorld("HWorld", "", N, WorldXMin, WorldXMax, 100, WorldMin, WorldMax);
   HWorld.SetStats(0);
   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);
   HWorld.Draw("axis");
   for(int i = 1; i < NLine; i++)
   {
      std::cout << "Histograms[i]: " << Histograms[i]->Integral() << std::endl;
      TH1D *H = (TH1D *)Histograms[i]->Clone();
      H->Divide(Histograms[0]);
      H->Draw("hist p l same");
   }


   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : 1 / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : 1/ 2, 10000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");


   TGraph G2;

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
   
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadHeight, WorldMin, WorldMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadHeight, WorldMin, WorldMax, 510, "");

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
      X5.Draw();
      X6.Draw();
   }
   if(LogX == false)
   {
      XL2.Draw();
   }
   // Y1.Draw();
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
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.04, "#it{z} = 1/2");


   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.9, MarginB * 0.4, X.c_str());

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(90);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "Finalization of Result April 24 (HB)");

   TLegend Legend(0.15, 0.88, 0.35, 0.88 - 0.04 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 1; i < NLine && i < 4; i++)
   {
      Legend.AddEntry(Histograms[i], Labels[i-1].c_str(),"pl");
   }
   Legend.Draw();

   TLegend Legend2(0.55, 0.88, 0.8, 0.88 - 0.04 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(Histograms[i], Labels[i-1].c_str(),"pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}

