// drawCov.C: Code to draw the covariance matrices pre-unfolding.
// Hannah Bossi, <hannah.bossi@cern.ch>

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

void drawCov1D_PostUnfolding(){
    
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

   
    
    TFile* inFile = TFile::Open("unfoldingE2C_ProjectedCovariance_09122026.root");
    
    TH2D* cov_Z = (TH2D*)inFile->Get("cov_Unfolded_FixedN_Z"); 
    
    // now draw and save the covariances
    gStyle->SetOptStat(0); 
    gStyle->SetOptTitle(0);
    TCanvas *c1 = new TCanvas("", "", 600, 600);
    c1->SetRightMargin(0.17);
    c1->SetLeftMargin(0.15);
    gStyle->SetPalette(109);
    //c1->SetLogz(); 
    // make the visualization symmetric
    
    std::cout << "Cov(10,10) = "
          << cov_Z->GetBinContent(11,11)
          << std::endl;
    double min = cov_Z->GetMinimum(); 
    double max = cov_Z->GetMaximum();
    double mult =  0.000001;
    if(abs(min) > abs(max)){
    cov_Z->GetZaxis()->SetRangeUser(mult*min, mult*abs(min));
    cov_Z->SetMaximum(mult*abs(min)); 
    cov_Z->SetMinimum(mult*min);
    }
    else{
    cov_Z->GetZaxis()->SetRangeUser(-1*mult*max, mult*abs(max));
    cov_Z->SetMaximum(mult*abs(max)); 
    cov_Z->SetMinimum(-1*mult*max); 
    }
    // cov_Z->SetMinimum(0);
    cov_Z->GetXaxis()->SetTitle("#theta_{L} bin index");
    cov_Z->GetYaxis()->SetTitle("#theta_{L} bin index");
    cov_Z->Draw("colz0"); 
    // for some reason the last bin has 0 counts so we will cut it off here. 
    cov_Z->GetXaxis()->SetRangeUser(0, 199); 
    cov_Z->GetYaxis()->SetRangeUser(0,199); 
    c1->SaveAs("cov_Z_PostUnfolding_09132026.pdf"); 
    
    
    TH2D* corr_Z = (TH2D*)inFile->Get("corr_Unfolded_FixedN_Z");
    corr_Z->SetTitle("Correlation matrix;bin i;bin j");
    corr_Z->SetMaximum(1.0); 
    corr_Z->SetMinimum(-1.0); 
    // corr_Z->Reset(); 
    
    // int nBins = cov_Z->GetNbinsX();
    // std::cout << "Number of bins is " << nBins << std::endl; 

    // for (int i = 1; i <= nBins; ++i) {
    //     double Cii = cov_Z->GetBinContent(i,i);

    //     for (int j = 1; j <= nBins; ++j) {
    //         double Cjj = cov_Z->GetBinContent(j,j);
    //         double Cij = cov_Z->GetBinContent(i,j);
    //         double rho = Cij / std::sqrt(Cii * Cjj);
    //         if(i == j){
    //             std::cout << "i = " << i << " j = " << j << " Cii: " << Cii << " Cij: " << Cij << " Cjj: " << Cjj  << " cov  = " << rho << std::endl; 
    //         }

    //         corr_Z->SetBinContent(i,j, rho);
    //     }
    // }
    
    TCanvas *c2 = new TCanvas("", "", 600, 600);
    c2->SetRightMargin(0.17);
    c2->SetLeftMargin(0.15);
   
    corr_Z->GetXaxis()->SetTitle("#theta_{L} bin index");
    corr_Z->GetYaxis()->SetTitle("#theta_{L} bin index");
    corr_Z->Draw("colz0"); 
    c2->SaveAs("corr_Z_PostUnfolding_09122026.pdf"); 
    
    // now make a plot to compare the statistical uncertainties 
    TH1D* statErr = (TH1D*)inFile->Get("stat_error_Unfolded_FixedN_Z");
    int nEvents = 1333529; 

    
    // now for the purposes of comparison, we want the statistical error on the raw data
    TFile* unfoldingFile = TFile::Open("unfoldingE2C_DataUnfolding_StatErrCheck_10022025.root"); 
    TH2D* raw = (TH2D*)unfoldingFile->Get("Bayesian_Unfoldediter4_Z"); 
    TH1D *raw1D= (TH1D*)raw->ProjectionX("statErrNorm");
    raw1D->Reset(); 
    for (int i = 1; i <= raw->GetNbinsX(); ++i) {
        double weight = 0;
        double error2 = 0;
        for (int j = 1; j <= raw->GetNbinsY(); ++j) {
            double binContent = raw->GetBinContent(i, j);
            double binError= raw->GetBinError(i,j);
            double binCenter = raw->GetYaxis()->GetBinCenter(j);
            weight += binContent*((binCenter));
            error2 += pow(binError*binCenter, 2);;
        }
        raw1D->SetBinContent(i, weight);
        raw1D->SetBinError(i,sqrt(error2));
        std::cout << "Setting Bin error i = " << i << " to be " << error2 << std::endl;
    }
    
    // now do the proper normalization 
    DivideByBin(*raw1D, zBins); 
    raw1D->Scale(1.0/nEvents); 
    
    TH1D* statErrNorm = (TH1D*)raw1D->Clone("statErrNorm"); 
    statErrNorm->Reset(); 
    for (int i = 1; i <= raw1D->GetNbinsX(); ++i) {
        statErrNorm->SetBinContent(i, raw1D->GetBinError(i)); 
        statErrNorm->SetBinError(i, 0.0); 
    }
    
    
    TCanvas *c3 = new TCanvas("", "", 800, 600);
    c3->SetRightMargin(0.05);
    c3->SetLeftMargin(0.12);
    c3->SetTopMargin(0.05); 
    c3->SetLogy();
    c3->SetTicks(1, 1); 

 
  
    statErr->SetMarkerColor(kBlack);
    statErr->SetMarkerStyle(20); 
    statErr->GetXaxis()->SetTitle("#theta_{L} bin index");
    statErr->GetYaxis()->SetTitle("Stat. Error");
    statErrNorm->SetMarkerColor(kBlack);
    statErrNorm->SetMarkerStyle(24); 
    
    TLegend *leg = new TLegend(0.40, 0.70, 0.68, 0.88);
    leg->AddEntry(statErr, "considering correlations", "p");
    leg->AddEntry(statErrNorm, "w/o considering correlations", "p");
    leg->SetTextSize(0.035);
    // no border + transparent background
    leg->SetBorderSize(0);
    leg->SetLineColor(0);
    leg->SetLineWidth(0);
    leg->SetFillStyle(0);
    
    statErr->Draw("p"); 
    statErrNorm->Draw(" p same"); 
    leg->Draw();
    c3->SaveAs("EEC_statErr_PostUnfolding_09132026.pdf");
    
    
}