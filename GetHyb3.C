#include <iostream>
#include <fstream>
#include <math.h>
#include "TMath.h"

void GetHyb3(int pt1_bin_dat, int pt2_bin_dat,bool ifgamma, bool ifvac, bool ifwake, bool ifnowake)
{
    
    
    //    //xxxxxxxxxxxxxxxxxxxxx------------------DATA HISTOGRAMS-------------------------------xxxxxxxxxxxxxxxxxxxxxxxxx
    TFile* fvac;
    TFile* fnw;
    TFile* fwake;
    
    TH1D* jet_pt_dat;
    TH1D* jet_pt_dat_wake;
    TH1D* jet_pt_dat_nw;
    
    TH3D* eeec_pt_hist;
    TH3D* eeec_pt_hist_wake;
    TH3D* eeec_pt_hist_nw;
    
    TH2D *eeec_r_slice;
    TH2D *eeec_r_slice_wake;
    TH2D *eeec_r_slice_nw;
    TH2D *eeec_diff;
    
    TH3D* eeec_diff_3;
    
    
    Double_t fromR = 0;
    Double_t toR = 1;
    Int_t binsR = 20;
    Double_t widthR = (toR-fromR)/binsR;
    Double_t new_binsR[21] = {};
    for (Int_t i = 0; i <= binsR; i++)
    {
        new_binsR[i] = (fromR + i * widthR);
    }
    
    Double_t fromRs = 0;
    Double_t toRs = M_PI/2;
    Int_t binsRs = 20;
    Double_t widthRs = (toRs-fromRs)/binsRs;
    Double_t new_binsRs[21] = {};
    for (Int_t i = 0; i <= binsRs; i++)
    {
        new_binsRs[i] = (fromRs + i * widthRs);
    }


    Double_t from_const = 40;
    Double_t to_const = 1000;
    Int_t bins_const = 48;
    Double_t width_const = (to_const-from_const)/bins_const;
    Double_t new_bins_const[49] = {};
    for (Int_t i = 0; i <= bins_const; i++)
    {
    new_bins_const[i] = (from_const + i * width_const);
    }
    
//    eeec_diff_3 = (TH3D*)eeec_pt_hist->Clone("diff");
//    eeec_diff = (TH2D*)eeec_pt_hist->Project3D("yx");
    eeec_diff_3 = new TH3D("eeec_diff","Vac-Wake",20, new_binsR, 20, new_binsRs, 48, new_bins_const);

    
    if (ifgamma == false)
    {
        cout<<"Inclusive"<<endl;
        fvac = TFile::Open("./Hybrid_Vac_hadron_passedjets_Vac_1203_67.root");
        fnw = TFile::Open("./Hybrid_hadron_passedjets_nowake_1203_67.root");
        fwake = TFile::Open("./Hybrid_hadron_passedjets_yeswake_1211_67.root");
    }
    else{
        cout<<"Gamma"<<endl;
        fvac = TFile::Open("./Hybrid_Vac_hadron_gamma_Vac_1129_23.root");
        fnw = TFile::Open("./Hybrid_hadron_gamma_nowake_1129_23.root");
        fwake = TFile::Open("./Hybrid_hadron_gamma_yeswake_1129_23.root");
    }
    
     cout<<"low  "<<endl;
    
    if (fvac)
    {
        
        if(ifgamma==true)
        {
            if(ifvac==true)
            {
             cout<<"low  "<<endl;
                jet_pt_dat = (TH1D*)fvac->Get("gam_jet_pt_hist");
                eeec_pt_hist = (TH3D*)fvac->Get("eeec_pt_hist");
                eeec_pt_hist->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_r_slice= (TH2D*)eeec_pt_hist->Project3D("yx");
                 cout<<"low  "<<endl;
                 
            }
            else if(ifwake==true)
            {
             cout<<"low  "<<endl;
                jet_pt_dat_wake = (TH1D*)fwake->Get("gam_jet_pt_hist");
                eeec_pt_hist_wake = (TH3D*)fwake->Get("eeec_pt_hist");
                eeec_pt_hist_wake->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_r_slice_wake= (TH2D*)eeec_pt_hist_wake->Project3D("yx");
                 cout<<"low  "<<endl;
            }
            else
            {
             cout<<"low  "<<endl;
                jet_pt_dat_nw = (TH1D*)fnw->Get("gam_jet_pt_hist");
                eeec_pt_hist_nw = (TH3D*)fnw->Get("eeec_pt_hist");
                eeec_pt_hist_nw->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_r_slice_nw = (TH2D*)eeec_pt_hist_nw->Project3D("yx");
                 cout<<"low  "<<endl;
            }
        }
        else
        {
            if(ifvac==true)
            {
            cout<<"low  "<<endl;
                jet_pt_dat = (TH1D*)fvac->Get("jet_pt_hist");
                jet_pt_dat->SetMarkerStyle(kOpenSquare);
                
                eeec_pt_hist = (TH3D*)fvac->Get("eeec_pt_hist");
                eeec_pt_hist->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_r_slice= (TH2D*)eeec_pt_hist->Project3D("yx");
                cout<<"low  "<<endl;
            }
            if(ifwake==true)
            {
            cout<<"low  "<<endl;
            cout<<"here"<<endl;
                jet_pt_dat_wake = (TH1D*)fwake->Get("jet_pt_hist");
                jet_pt_dat_wake->SetLineColor(kRed);
                jet_pt_dat_wake->SetMarkerColor(kRed);
                jet_pt_dat_wake->SetMarkerStyle(kFullCircle);
                
                eeec_pt_hist_wake = (TH3D*)fwake->Get("eeec_pt_hist");
                eeec_pt_hist_wake->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_r_slice_wake= (TH2D*)eeec_pt_hist_wake->Project3D("yx");
                cout<<"low  "<<endl;
            }
            if(ifnowake==true)
            {
            cout<<"low  "<<endl;
                jet_pt_dat_nw = (TH1D*)fnw->Get("jet_pt_hist");
                jet_pt_dat_nw->SetLineColor(kGreen+2);
                jet_pt_dat_nw->SetMarkerColor(kGreen+2);
                jet_pt_dat_nw->SetMarkerStyle(22);
                
                eeec_pt_hist_nw = (TH3D*)fnw->Get("eeec_pt_hist");
                eeec_pt_hist_nw->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_r_slice_nw = (TH2D*)eeec_pt_hist_nw->Project3D("yx");
                cout<<"low  "<<endl;
            }
        }
        cout<<"heregain"<<endl;
//        jet_pt_dat->RebinX(20);
        TCanvas *cpt = new TCanvas();
//        jet_pt_dat->SetMaximum(1);
//        jet_pt_dat->SetMinimum(1e-8);
        jet_pt_dat->SetStats(0);
        jet_pt_dat->Draw();
        cout<<"heregain"<<endl;
    
         ////Getting bin edges
        float low_pt = eeec_pt_hist->GetZaxis()->GetBinLowEdge(pt1_bin_dat);
        float high_pt = eeec_pt_hist->GetZaxis()->GetBinUpEdge(pt2_bin_dat);
        
        float low_pt_wake = eeec_pt_hist_wake->GetZaxis()->GetBinLowEdge(pt1_bin_dat);
        float high_pt_wake = eeec_pt_hist_wake->GetZaxis()->GetBinUpEdge(pt2_bin_dat);
        
        float low_pt_nw = eeec_pt_hist_nw->GetZaxis()->GetBinLowEdge(pt1_bin_dat);
        float high_pt_nw = eeec_pt_hist_nw->GetZaxis()->GetBinUpEdge(pt2_bin_dat);
        
        float low_pt_jet = jet_pt_dat->GetXaxis()->GetBinLowEdge(pt1_bin_dat);
        float high_pt_jet = jet_pt_dat->GetXaxis()->GetBinUpEdge(pt2_bin_dat);
        
        //checking if axes are same on both plots
        cout<<"low  "<<low_pt<<endl;
        cout<<"high  "<<high_pt<<endl;
        cout<<"low  "<<low_pt_jet<<endl;
        cout<<"high  "<<high_pt_jet<<endl;
        
        
        double num_jets_dat = 0;
        double num_entries_dat = 0;
       
        for (int i=pt1_bin_dat; i<=pt2_bin_dat; i++)
        {
            num_entries_dat = jet_pt_dat->GetBinContent(i);
            num_jets_dat += num_entries_dat;
            
        }
        cout <<"num dat no overflow "<<num_jets_dat<<endl;
        //        TCanvas *cfnew = new TCanvas();
        //        TH1D *eeec_r_slice_phi= eeec_r_slice->ProjectionX("eeec_phi",1,2);
        //        TH1D *eeec_r_slice_phi1= eeec_r_slice->ProjectionX("eeec_phi1",3,4);
        //        eeec_r_slice_phi->Draw();
        //        eeec_r_slice_phi1->SetLineColor(kRed);
        //        eeec_r_slice_phi1->Draw("SAME");
        //        eeec_r_slice_phi1->Draw("SAME");
        
        
        
        //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
        //DEFINING ALL PLOTS WITH CORRECT AXES AND SCALINGS
        //2D plots if needed
        //                TCanvas *c = new TCanvas();
        //                eec_pt_2d_dat->Draw("LEGO1");
        //                //Jet spectrum
        //                        TCanvas *cpt = new TCanvas();
        //                        jet_pt_dat->Draw();
        jet_pt_dat->GetXaxis()->SetTitle("pT (GeV/c)");
        jet_pt_dat->GetYaxis()->SetTitle("Entries");
        double entries = jet_pt_dat->GetEntries();
        cout<<"all entries data "<<entries<<endl;
        
        
       
        //-------------Defining global attributes of plots------------------------------
 
        //Writing pT range on plots
         TLatex latex, latex1, latex2, latex3, latex4, latex5;
        latex.SetTextColor(kWhite);
        latex1.SetTextColor(kWhite);
        latex2.SetTextColor(kWhite);
        latex3.SetTextColor(kWhite);
        latex4.SetTextColor(kWhite);
        latex5.SetTextColor(kWhite);
         latex.SetNDC ();
        //               const char *str = "p^{ch}_{T,jet}";
        const char *str = "p_{T,jet}";
        latex.SetTextSize(0.04);
        latex.SetTextFont(42);
        latex.SetTextAlign(50);
        //Writing particle track cut range on plots
        
        latex1.SetTextFont(42);
        latex1.SetTextAlign(50);
        latex1.SetNDC ();
        const char *str1 = "p^{ch}_{T,min}";
        latex1.SetTextSize(0.04);
        //Writing jet algorithm
       
        latex2.SetTextFont(42);
        latex2.SetTextAlign(50);
        latex2.SetNDC ();
        const char *str2 = "anti-k_{T}";
        latex2.SetTextSize(0.04);
        //Writing collision energy
        latex3.SetTextFont(42);
        latex3.SetTextAlign(50);
        latex3.SetNDC ();
        const char *str3 = "#sqrt{s} = 5.02 TeV";
        latex3.SetTextSize(0.04);
        //Writing R_L for the 3 point
        const char *str4 = "R_{L}";
        //Writing num jets and if data is corrected or uncorrected
        const char *str5 = "Uncorrected";
        const char *str6 = "#gamma";
        const char *str7 = "p_{T}^{#gamma}";
        latex4.SetTextFont(42);
        latex4.SetTextAlign(50);
        latex4.SetNDC ();
        latex4.SetTextSize(0.04);
        
        latex5.SetNDC ();
        latex5.SetTextFont(42);
        latex5.SetTextAlign(50);
        latex5.SetTextSize(0.04);
        
        
        //Format of application:               latex.DrawLatex(0.8,0.7 ,Form("%.1d GeV < %s < %.1d GeV", low, str, high));
        // Adding a line at 1 to guide the eye
        TLine *line = new TLine(0.0 ,1 ,1,1);
        line->SetLineColorAlpha(kBlue,0.35);
        line->SetLineWidth(2);
        line->SetLineStyle(7);
        
        //xxxxxxxxxxxxxxxxxxxxxxxxxxxx------------------DATA/MC PLOTS------------------------------xxxxxxxxxxxxxxxxxxxxxxxxxx//
        //        eeec_pt_hist->RebinY(2);
        //        eeec_pt_hist->RebinX(2);
        const char *strphi = "#phi";
        const char *strx = "#xi";
        eeec_pt_hist->GetYaxis()->SetTitle(strphi);
        eeec_pt_hist->GetYaxis()->LabelsOption("h");
        eeec_pt_hist->GetXaxis()->SetTitle(strx);
        eeec_pt_hist->GetXaxis()->LabelsOption("h");
        
        
        TCanvas *cf = new TCanvas();
        eeec_r_slice->SetStats(0);
        eeec_r_slice->Scale(1./jet_pt_dat->Integral(pt1_bin_dat,pt2_bin_dat),"width");
        eeec_r_slice_wake->Scale(1./jet_pt_dat_wake->Integral(pt1_bin_dat,pt2_bin_dat),"width");
        eeec_r_slice_nw->Scale(1./jet_pt_dat_nw->Integral(pt1_bin_dat,pt2_bin_dat),"width");
        
//        for(int i = 1; i<=eeec_r_slice->GetXaxis()->GetNbins(); i++)
//        {
//            for(int j = 1; j<=eeec_r_slice->GetYaxis()->GetNbins(); j++){
//                cout<<i<<j<<endl;
//                double value_vac = eeec_r_slice->GetBinContent(i,j);
//                double value_wake = eeec_r_slice_wake->GetBinContent(i,j);
//                double value = ((value_vac - value_wake)/value_vac);
////                eeec_diff->SetBinContent(i,j,value);
//                eeec_diff->SetBinContent(i,j,value);
//
//                cout<<eeec_r_slice->GetBinContent(i,j)<<endl;
//                cout<<value_vac<<endl;
//                cout<<"picl"<<endl;
//                cout<<eeec_r_slice_wake->GetBinContent(i,j)<<endl;
//                cout<<value_wake<<endl;
////                cout<<"value in set"<<eeec_diff->GetBinContent(i,j)<<endl;
////                 cout<<"value in set"<<eeec_r_slice_wake->GetBinContent(i,j)<<endl;
//
//            }
//        }
//         cout<<eeec_r_slice->GetMaximum()<<endl;
//        cout<<eeec_r_slice->GetMinimum()<<endl;
//        eeec_r_slice_wake->Add(eeec_r_slice, -1.0);
//        eeec_r_slice->SetMaximum(0.0103458);
//        eeec_r_slice->SetMinimum(0.00141669);
        
//        eeec_r_slice->SetMaximum(eeec_r_slice->GetMaximum());
//        eeec_r_slice->SetMinimum(eeec_r_slice->GetMinimum());
//
        cout<<eeec_r_slice->GetMaximum()<<endl;
        cout<<eeec_r_slice->GetMinimum()<<endl;
        eeec_r_slice_wake->Draw("colz");
//        eeec_diff->SetMaximum(eeec_diff->GetMaximum());
//        eeec_diff->SetMinimum(eeec_diff->GetMinimum());
//        eeec_diff->Draw("box");
       
       
       //xxxxxxxxxxxxxxxxxxxxxxxxxxxx------------------Legends&Latex------------------------------xxxxxxxxxxxxxxxxxxxxxxxxxx//
        latex.SetTextSize(0.03);latex1.SetTextSize(0.03);latex2.SetTextSize(0.03);latex3.SetTextSize(0.03);latex4.SetTextSize(0.03);latex5.SetTextSize(0.03);
        
        //    latex1.DrawLatex(0.60,0.70 ,Form("%s > 1 GeV/c", str1));
        latex2.DrawLatex(0.60,0.75 ,Form("%s , R = 0.8", str2));
        latex4.DrawLatex(0.60,0.65 ,Form("0.6 < %s < 0.7 ", str4));
        
        if(ifgamma==true)
        {
            latex4.DrawLatex(0.60,0.85 ,Form("%s - tagged", str6));
            latex.DrawLatex(0.60,0.70 ,Form("%.1f GeV/c < %s < %.1f GeV/c", low_pt, str7, high_pt));
        }
        else{
            latex4.DrawLatex(0.60,0.85 ,Form("Inclusive"));
            latex.DrawLatex(0.60,0.70 ,Form("%.1f GeV/c < %s < %.1f GeV/c", low_pt, str, high_pt));
        }
        
        if(ifvac==true)
        {
            latex3.DrawLatex(0.60,0.80 ,Form("Hybrid Model, hadrons, vacuum"));
        }
        else if(ifwake==true)
        {
            latex3.DrawLatex(0.60,0.80 ,Form("Hybrid Model, hadrons, wake"));
        }
        else
        {
            latex3.DrawLatex(0.60,0.80 ,Form("Hybrid Model, hadrons, no wake"));
        }
        
        TLegend *legend_dat_rat = new TLegend(0.15,0.84,0.40,0.90);
        
        //        legend_dat_rat->SetBorderSize(0);
        //        legend_dat_rat->SetTextFont(42);
        //        legend_dat_rat->SetTextSize(0.035);
        //legend_dat_rat->Draw();
        
        
        if (!eeec_pt_hist)
            printf("No such histogram found!\n");
        
    }
    else
        printf("No such file found!\n");
}
