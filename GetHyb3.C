#include <iostream>
#include <fstream>
#include <math.h>
#include "TMath.h"

void GetHyb3(int pt1_bin_dat, int pt2_bin_dat,bool ifgamma, bool ifvac, bool ifwake, bool ifnowake, bool ifdivide)
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
    

    if (ifgamma == false)
    {
        cout<<"Inclusive"<<endl;
        fvac = TFile::Open("/Users/ar2545/Desktop/hybdat/Hybrid_Vac_hadron_passedjets_Vac_1211_67.root");
        fwake = TFile::Open("/Users/ar2545/Desktop/hybdat/Hybrid_hadron_passedjets_yeswake_1211_67.root");
        fnw = TFile::Open("/Users/ar2545/Desktop/hybdat/Hybrid_hadron_passedjets_nowake_1203_67.root");
    }
    else{
        cout<<"Gamma"<<endl;
        fvac = TFile::Open("/Users/ar2545/Desktop/hybdat/Hybrid_Vac_hadron_gamma_Vac_1129_23.root");
        fnw = TFile::Open("/Users/ar2545/Desktop/hybdat/Hybrid_hadron_gamma_nowake_1129_23.root");
        fwake = TFile::Open("/Users/ar2545/Desktop/hybdat/Hybrid_hadron_gamma_yeswake_1129_23.root");
    }
    
     cout<<"low  "<<endl;
     TCanvas *cpt = new TCanvas();//for drawing jet pt
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
                eeec_pt_hist_wake->SetName("eeec_pt_hist_wake");
                eeec_pt_hist_wake->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_r_slice_wake= (TH2D*)eeec_pt_hist_wake->Project3D("yx");
                cout<<"low  "<<endl;
            }
            else
            {
                cout<<"low  "<<endl;
                jet_pt_dat_nw = (TH1D*)fnw->Get("gam_jet_pt_hist");
                eeec_pt_hist_nw = (TH3D*)fnw->Get("eeec_pt_hist");
                eeec_pt_hist_nw->SetName("eeec_pt_hist_nw");
                eeec_pt_hist_nw->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_r_slice_nw = (TH2D*)eeec_pt_hist_nw->Project3D("yx");
                cout<<"low  "<<endl;
            }
        }
        else
        {
            
                jet_pt_dat = (TH1D*)fvac->Get("jet_pt_hist");
                jet_pt_dat->SetMarkerStyle(kOpenSquare);
                jet_pt_dat->SetStats(0);
                
                eeec_pt_hist = (TH3D*)fvac->Get("eeec_pt_hist");
                eeec_pt_hist->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_r_slice= (TH2D*)eeec_pt_hist->Project3D("yx");
                eeec_r_slice->Scale(1./jet_pt_dat->Integral(pt1_bin_dat,pt2_bin_dat),"width");
                
            
                cout<<"wake "<<endl;
                
                jet_pt_dat_wake = (TH1D*)fwake->Get("jet_pt_hist");
                jet_pt_dat_wake->SetName("jet_pt_hist_wake");
                jet_pt_dat_wake->SetLineColor(kRed);
                jet_pt_dat_wake->SetMarkerColor(kRed);
                jet_pt_dat_wake->SetMarkerStyle(kFullCircle);
                
                jet_pt_dat_wake->SetStats(0);
//                jet_pt_dat_wake->Draw();
                
                eeec_pt_hist_wake = (TH3D*)fwake->Get("eeec_pt_hist");
                eeec_pt_hist_wake->SetName("eeec_pt_hist_wake");
                eeec_pt_hist_wake->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_r_slice_wake= (TH2D*)eeec_pt_hist_wake->Project3D("yx");
                eeec_r_slice_wake->Scale(1./jet_pt_dat_wake->Integral(pt1_bin_dat,pt2_bin_dat),"width");
            
                cout<<"no wake "<<endl;
                jet_pt_dat_nw = (TH1D*)fnw->Get("jet_pt_hist");
                jet_pt_dat_nw->SetName("jet_pt_hist_nw");
                jet_pt_dat_nw->SetLineColor(kGreen+2);
                jet_pt_dat_nw->SetMarkerColor(kGreen+2);
                jet_pt_dat_nw->SetMarkerStyle(22);
                
                jet_pt_dat_nw->SetStats(0);
//                jet_pt_dat_nw->Draw();
                
                eeec_pt_hist_nw = (TH3D*)fnw->Get("eeec_pt_hist");
                eeec_pt_hist_nw->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_pt_hist_nw->SetName("eeec_pt_hist_nw");
                eeec_r_slice_nw = (TH2D*)eeec_pt_hist_nw->Project3D("yx");
                eeec_r_slice_nw->Scale(1./jet_pt_dat_nw->Integral(pt1_bin_dat,pt2_bin_dat),"width");

           
        }
        
        
        float low_pt, high_pt, low_pt_wake,high_pt_wake,low_pt_nw,high_pt_nw,low_pt_jet,high_pt_jet,low_pt_jet_wake,high_pt_jet_wake,low_pt_jet_nw,high_pt_jet_nw;
        double num_entries_dat, num_entries_dat_wake, num_entries_dat_nw;
        double num_jets_dat, num_jets_dat_wake, num_jets_dat_nw;
        ////Getting bin edges
        if(ifvac==true)
        {
        jet_pt_dat->Draw();
            low_pt = eeec_pt_hist->GetZaxis()->GetBinLowEdge(pt1_bin_dat);
            high_pt = eeec_pt_hist->GetZaxis()->GetBinUpEdge(pt2_bin_dat);
            low_pt_jet = jet_pt_dat->GetXaxis()->GetBinLowEdge(pt1_bin_dat);
            high_pt_jet = jet_pt_dat->GetXaxis()->GetBinUpEdge(pt2_bin_dat);
            
            cout<<"low 2d "<<low_pt<<endl;
            cout<<"high 2d "<<high_pt<<endl;
            
            cout<<"low jet "<<low_pt_jet<<endl;
            cout<<"high jet  "<<high_pt_jet<<endl;
            
            for (int i=pt1_bin_dat; i<=pt2_bin_dat; i++)
            {
                num_entries_dat = jet_pt_dat->GetBinContent(i);
                num_jets_dat += num_entries_dat;
                
            }
            cout <<"num dat no overflow "<<num_jets_dat<<endl;
        }
        if(ifwake==true)
        {
        jet_pt_dat_wake->Draw();
        jet_pt_dat->Draw("SAME");
            low_pt_wake = eeec_pt_hist_wake->GetZaxis()->GetBinLowEdge(pt1_bin_dat);
            high_pt_wake = eeec_pt_hist_wake->GetZaxis()->GetBinUpEdge(pt2_bin_dat);
            low_pt_jet_wake = jet_pt_dat_wake->GetXaxis()->GetBinLowEdge(pt1_bin_dat);
            high_pt_jet_wake = jet_pt_dat_wake->GetXaxis()->GetBinUpEdge(pt2_bin_dat);
            
            cout<<"low 2d "<<low_pt_wake<<endl;
            cout<<"high 2d "<<high_pt_wake<<endl;
            
            cout<<"low jet "<<low_pt_jet_wake<<endl;
            cout<<"high jet  "<<high_pt_jet_wake<<endl;
            
            for (int i=pt1_bin_dat; i<=pt2_bin_dat; i++)
            {
                num_entries_dat_wake = jet_pt_dat_wake->GetBinContent(i);
                num_jets_dat_wake += num_entries_dat_wake;
                
            }
            cout <<"num dat no overflow "<<num_jets_dat_wake<<endl;
        }
        if(ifnowake==true)
        {
        jet_pt_dat_nw->Draw();
            low_pt_nw = eeec_pt_hist_nw->GetZaxis()->GetBinLowEdge(pt1_bin_dat);
            high_pt_nw = eeec_pt_hist_nw->GetZaxis()->GetBinUpEdge(pt2_bin_dat);
            low_pt_jet_nw = jet_pt_dat_nw->GetXaxis()->GetBinLowEdge(pt1_bin_dat);
            high_pt_jet_nw = jet_pt_dat_nw->GetXaxis()->GetBinUpEdge(pt2_bin_dat);
            
            cout<<"low 2d "<<low_pt_nw<<endl;
            cout<<"high 2d "<<high_pt_nw<<endl;
            
            cout<<"low jet "<<low_pt_jet_nw<<endl;
            cout<<"high jet  "<<high_pt_jet_nw<<endl;
            
            for (int i=pt1_bin_dat; i<=pt2_bin_dat; i++)
            {
                num_entries_dat_nw = jet_pt_dat_nw->GetBinContent(i);
                num_jets_dat_nw += num_entries_dat_nw;
                
            }
            cout <<"num dat no overflow "<<num_jets_dat_nw<<endl;
        }
    
        
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
        
        cout<<"made it here"<<endl;
         //xxxxxxxxxxxxxxxxxxxxxxxxxxxx------------------Legends&Latex------------------------------xxxxxxxxxxxxxxxxxxxxxxxxxx//
        latex.SetTextSize(0.03);latex1.SetTextSize(0.03);latex2.SetTextSize(0.03);latex3.SetTextSize(0.03);latex4.SetTextSize(0.03);latex5.SetTextSize(0.03);
        //xxxxxxxxxxxxxxxxxxxxxxxxxxxx------------------DATA/MC PLOTS------------------------------xxxxxxxxxxxxxxxxxxxxxxxxxx//
        
        const char *strphi = "#phi";
        const char *strx = "#xi";
        
        
        TCanvas *cvac = new TCanvas();
        TCanvas *cwake = new TCanvas();
        TCanvas *cnw = new TCanvas();
        
        if(ifvac==true)
        {
            cvac->cd();
            eeec_r_slice->GetYaxis()->SetTitle(strphi);
            eeec_r_slice->GetYaxis()->LabelsOption("h");
            eeec_r_slice->GetXaxis()->SetTitle(strx);
            eeec_r_slice->GetXaxis()->LabelsOption("h");
            
            eeec_r_slice->SetStats(0);
            eeec_r_slice->Draw("colz");
            latex3.DrawLatex(0.60,0.80 ,Form("Hybrid Model, hadrons, vacuum"));
            latex.DrawLatex(0.60,0.70 ,Form("%.1f GeV/c < %s < %.1f GeV/c", low_pt, str, high_pt));
            latex2.DrawLatex(0.60,0.75 ,Form("%s , R = 0.8", str2));
            latex4.DrawLatex(0.60,0.65 ,Form("0.6 < %s < 0.7 ", str4));
             if(ifgamma==true)
            {
                latex4.DrawLatex(0.60,0.85 ,Form("%s - tagged", str6));
               }
            else{
                latex4.DrawLatex(0.60,0.85 ,Form("Inclusive"));
            }
        }
        if(ifwake==true)
        {
            cwake->cd();
            eeec_r_slice_wake->GetYaxis()->SetTitle(strphi);
            eeec_r_slice_wake->GetYaxis()->LabelsOption("h");
            eeec_r_slice_wake->GetXaxis()->SetTitle(strx);
            eeec_r_slice_wake->GetXaxis()->LabelsOption("h");
            cout<<"here"<<endl;
            
            eeec_r_slice_wake->SetStats(0);
            eeec_r_slice_wake->Draw("colz");
            latex3.DrawLatex(0.60,0.80 ,Form("Hybrid Model, hadrons, jet + wake"));
            latex.DrawLatex(0.60,0.70 ,Form("%.1f GeV/c < %s < %.1f GeV/c", low_pt_wake, str, high_pt_wake));
            latex2.DrawLatex(0.60,0.75 ,Form("%s , R = 0.8", str2));
            latex4.DrawLatex(0.60,0.65 ,Form("0.6 < %s < 0.7 ", str4));
             if(ifgamma==true)
            {
                latex4.DrawLatex(0.60,0.85 ,Form("%s - tagged", str6));
               }
            else{
                latex4.DrawLatex(0.60,0.85 ,Form("Inclusive"));
            }
        }
        if(ifnowake==true)
        {
            eeec_pt_hist_nw->GetYaxis()->SetTitle(strphi);
            eeec_pt_hist_nw->GetYaxis()->LabelsOption("h");
            eeec_pt_hist_nw->GetXaxis()->SetTitle(strx);
            eeec_pt_hist_nw->GetXaxis()->LabelsOption("h");
            
            cnw->cd();
            eeec_r_slice_nw->SetStats(0);
            eeec_r_slice_nw->Draw("colz");
            latex3.DrawLatex(0.60,0.80 ,Form("Hybrid Model, hadrons, no wake"));
            latex.DrawLatex(0.60,0.70 ,Form("%.1f GeV/c < %s < %.1f GeV/c", low_pt_nw, str, high_pt_nw));
            latex2.DrawLatex(0.60,0.75 ,Form("%s , R = 0.8", str2));
            latex4.DrawLatex(0.60,0.65 ,Form("0.6 < %s < 0.7 ", str4));
            
            if(ifgamma==true)
            {
                latex4.DrawLatex(0.60,0.85 ,Form("%s - tagged", str6));
               }
            else{
                latex4.DrawLatex(0.60,0.85 ,Form("Inclusive"));
            }
        }
        
        if(ifgamma==true)
        {
            latex4.DrawLatex(0.60,0.85 ,Form("%s - tagged", str6));
            //            latex.DrawLatex(0.60,0.70 ,Form("%.1f GeV/c < %s < %.1f GeV/c", low_pt, str7, high_pt));
        }
        else{
            latex4.DrawLatex(0.60,0.85 ,Form("Inclusive"));
            latex2.DrawLatex(0.60,0.75 ,Form("%s , R = 0.8", str2));
            latex4.DrawLatex(0.60,0.65 ,Form("0.6 < %s < 0.7 ", str4));
            //            latex.DrawLatex(0.60,0.70 ,Form("%.1f GeV/c < %s < %.1f GeV/c", low_pt, str, high_pt));
        }
        
        TCanvas *cdiv = new TCanvas();
        if(ifdivide == true)
        {
            cdiv->cd();
            TH2D* eeec_r_slice_wake_cop =(TH2D*)eeec_r_slice_wake->Clone();
            eeec_r_slice_wake_cop->Divide(eeec_r_slice);
            cout<<eeec_r_slice_wake->GetEntries()<<endl;
            eeec_r_slice_wake_cop->Draw("colz");
            latex3.DrawLatex(0.60,0.80 ,Form("Hybrid Model, hadrons, jet+wake/jet"));
            latex.DrawLatex(0.60,0.70 ,Form("%.1f GeV/c < %s < %.1f GeV/c", low_pt_wake, str, high_pt_wake));
            latex2.DrawLatex(0.60,0.75 ,Form("%s , R = 0.8", str2));
            latex4.DrawLatex(0.60,0.65 ,Form("0.6 < %s < 0.7 ", str4));
            if(ifgamma==true)
            {
                latex4.DrawLatex(0.60,0.85 ,Form("%s - tagged", str6));
            }
            else{
                latex4.DrawLatex(0.60,0.85 ,Form("Inclusive"));
            }
        }
        
    
        
        
        if (!eeec_pt_hist)
            printf("No such histogram found!\n");
        
    }
    else
        printf("No such file found!\n");
}
