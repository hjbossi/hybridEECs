void plotPtDependence(bool doEEC = true) {
  TFile *wakeR040 = TFile::Open("Hybrid_Projhadron_passedjets_yeswake_3101_n1.root"); 
  wakeR040->ls();
  TH2D* eec_pt_hist_wake; 
  if(doEEC){
    eec_pt_hist_wake = (TH2D*)wakeR040->Get("eec_pt_hist");
  } 
  else{
    eec_pt_hist_wake = (TH2D*)wakeR040->Get("e3c_pt_hist");
  }
  eec_pt_hist_wake->SetName("eec_pt_hist_wake");
  // print out the y axis range
  cout << "y axis range: " << eec_pt_hist_wake->GetYaxis()->GetXmin() << " to " << eec_pt_hist_wake->GetYaxis()->GetXmax() << endl;
  TH1D* eec_projection_40_100 = eec_pt_hist_wake->ProjectionX("eec_projection_40_100", eec_pt_hist_wake->GetYaxis()->FindBin(40), eec_pt_hist_wake->GetYaxis()->FindBin(100));
  eec_projection_40_100->SetName("eec_projection_40_100");
  TH1D* eec_projection_100_200 = eec_pt_hist_wake->ProjectionX("eec_projection_100_200", eec_pt_hist_wake->GetYaxis()->FindBin(100), eec_pt_hist_wake->GetYaxis()->FindBin(200));
  eec_projection_100_200->SetName("eec_projection_100_200");
  TH1D* eec_projection_200_300 = eec_pt_hist_wake->ProjectionX("eec_projection_200_300", eec_pt_hist_wake->GetYaxis()->FindBin(200), eec_pt_hist_wake->GetYaxis()->FindBin(300));
  eec_projection_200_300->SetName("eec_projection_200_300");
  TH1D* eec_projection_300_400 = eec_pt_hist_wake->ProjectionX("eec_projection_300_400", eec_pt_hist_wake->GetYaxis()->FindBin(300), eec_pt_hist_wake->GetYaxis()->FindBin(400));
  eec_projection_300_400->SetName("eec_projection_300_400");
  TH1D* eec_projection_400_500 = eec_pt_hist_wake->ProjectionX("eec_projection_400_500", eec_pt_hist_wake->GetYaxis()->FindBin(400), eec_pt_hist_wake->GetYaxis()->FindBin(500));
  eec_projection_400_500->SetName("eec_projection_400_500");
  TH1D* eec_projection_500_1000 = eec_pt_hist_wake->ProjectionX("eec_projection_500_1000", eec_pt_hist_wake->GetYaxis()->FindBin(500), eec_pt_hist_wake->GetYaxis()->FindBin(1000));
  eec_projection_500_1000->SetName("eec_projection_500_600");

  // scale by the area and width
  eec_projection_40_100->Scale(1.0/eec_projection_40_100->Integral(), "width");
  eec_projection_100_200->Scale(1.0/eec_projection_100_200->Integral(), "width");
  eec_projection_200_300->Scale(1.0/eec_projection_200_300->Integral(), "width");
  eec_projection_300_400->Scale(1.0/eec_projection_300_400->Integral(), "width");
  eec_projection_400_500->Scale(1.0/eec_projection_400_500->Integral(), "width");
  eec_projection_500_1000->Scale(1.0/eec_projection_500_1000->Integral(), "width");


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->SetLogy();
  c1->SetLogx();
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.05);
  c1->SetBottomMargin(0.15);
  c1->SetTopMargin(0.05);

  // create a legend in the bottom left corner 
  TLegend *legend = new TLegend(0.20, 0.25, 0.5, 0.6);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.045);
  legend->AddEntry(eec_projection_40_100, "40 < #it{p}_{T, jet} < 100 GeV/#it{c}", "lp");
  legend->AddEntry(eec_projection_100_200, "100 < #it{p}_{T, jet} < 200 GeV/#it{c}", "lp");
  legend->AddEntry(eec_projection_200_300, "200 < #it{p}_{T, jet} < 300 GeV/#it{c}", "lp");
  legend->AddEntry(eec_projection_300_400, "300 < #it{p}_{T, jet} < 400 GeV/#it{c}", "pl");
  legend->AddEntry(eec_projection_400_500, "400 < #it{p}_{T, jet} < 500 GeV/#it{c}", "pl");
  legend->AddEntry(eec_projection_500_1000, "500 < #it{p}_{T, jet} < 1000 GeV/#it{c}", "pl");

  //include tlatex above legend that says "Hybrid model w/ wake"
  TLatex *latex = new TLatex();
  latex->SetTextSize(0.045);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);



  // restreec_projection_40_100
  eec_projection_40_100->GetXaxis()->SetTitle("#it{R}_{L}");
  if(doEEC){
    eec_projection_40_100->GetYaxis()->SetTitle("Normalized EEC");
  }
  else{
    eec_projection_40_100->GetYaxis()->SetTitle("Normalized EEEC");
  }
  eec_projection_40_100->GetXaxis()->SetTitleSize(0.05);
  eec_projection_40_100->GetYaxis()->SetTitleSize(0.05);
  eec_projection_40_100->GetXaxis()->SetRangeUser(0.005, 1); 
  eec_projection_100_200->GetXaxis()->SetRangeUser(0.005, 1);
  eec_projection_200_300->GetXaxis()->SetRangeUser(0.005, 1);
  eec_projection_300_400->GetXaxis()->SetRangeUser(0.005, 1);
  eec_projection_400_500->GetXaxis()->SetRangeUser(0.005, 1);
  eec_projection_500_1000->GetXaxis()->SetRangeUser(0.005, 1);
  eec_projection_40_100->GetYaxis()->SetRangeUser(1e-6, 500);

  eec_projection_40_100->Draw();
  eec_projection_40_100->SetLineColor(kRed);
  eec_projection_40_100->SetMarkerColor(kRed);
  eec_projection_40_100->SetMarkerStyle(20);
  eec_projection_100_200->Draw("same");
  eec_projection_100_200->SetLineColor(kBlue);
  eec_projection_100_200->SetMarkerColor(kBlue);
  eec_projection_100_200->SetMarkerStyle(20);
  eec_projection_200_300->Draw("same");
  eec_projection_200_300->SetLineColor(kGreen);
  eec_projection_200_300->SetMarkerColor(kGreen);
  eec_projection_200_300->SetMarkerStyle(20);
  eec_projection_300_400->Draw("same");
  eec_projection_300_400->SetLineColor(kYellow);
  eec_projection_300_400->SetMarkerColor(kYellow);
  eec_projection_300_400->SetMarkerStyle(20);
  eec_projection_400_500->Draw("same");
  eec_projection_400_500->SetLineColor(kMagenta);
  eec_projection_400_500->SetMarkerColor(kMagenta);
  eec_projection_400_500->SetMarkerStyle(20);
  eec_projection_500_1000->Draw("same");  
  eec_projection_500_1000->SetLineColor(kCyan);
  eec_projection_500_1000->SetMarkerColor(kCyan);
  eec_projection_500_1000->SetMarkerStyle(20);
  legend->Draw();
  latex->DrawLatexNDC(0.6, 0.88, "Hybrid model w/ wake");
  latex->DrawLatexNDC(0.6, 0.83, "Hadrons, n = 1.0");
  latex->DrawLatexNDC(0.6, 0.78, "Full anti-#it{k}_{T} jets, #it{R} = 0.8");
  if(doEEC){
    c1->SaveAs("eec_projection_pt_dependence_R080.png");
  }
  else{
    c1->SaveAs("eeec_projection_pt_dependence_R080.png");
  }


}