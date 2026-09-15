void MCStatSyst(const char* infile = "CorrectedData.root",
                          const char* histname = "HzMCGenBeforeRef",
                          const char* outfile = "MCStatSyst_Z.root") {
  // Open input file
  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    Error("MCStatSyst", "Cannot open input file %s", infile);
    return;
  }

  // Retrieve histogram
  TH1* h = dynamic_cast<TH1*>(fin->Get(histname));

  // Clone histogram for output
  TH1* h_shifted = (TH1*)h->Clone(Form("%s_plus1sigma", histname));
  h_shifted->SetTitle(Form("%s (central + 1σ stat)", h->GetTitle()));

  // Shift bin contents by +1 × stat. error
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double val = h->GetBinContent(i);
    double err = h->GetBinError(i);
    h_shifted->SetBinContent(i, val + err);
    // keep the same error bar
    h_shifted->SetBinError(i, err);
  }

  // Open output file and write
  TFile* fout = TFile::Open(outfile, "RECREATE");
  h->Write(); 
  h_shifted->Write();
  fout->Close();
  fin->Close();

}