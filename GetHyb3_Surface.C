#include <iostream>
#include <fstream>
#include <math.h>
#include "TMath.h"

void GetHyb3_Surface(int pt1_bin_dat, int pt2_bin_dat,bool ifgamma, bool ifvac, bool ifwake, bool ifnowake, bool ifdivide,bool ifdividenw, bool ifsave)
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
    TH2D* eeec_r_slice_wake_cop;
    TH2D *eeec_r_slice_nw;
    TH2D* eeec_r_slice_nw_cop;
    TH2D *eeec_diff;

    
    
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
    

    Double_t rL_low = 0.6;
    Double_t rL_high = 0.7;

    
    if (ifgamma == false)
    {
        cout<<"Inclusive"<<endl;
	// Hybrid_Vac_hadron_passedjets_Vac_0301_%d%d_n1.root
        fvac = TFile::Open(Form("Hybrid_Vac_hadron_passedjets_Vac_1812_%d%d.root", (int)(rL_low*10), (int)(rL_high*10)));
        fwake = TFile::Open(Form("Hybrid_hadron_passedjets_yeswake_1812_%d%d.root", (int)(rL_low*10), (int)(rL_high*10)));
        fnw = TFile::Open(Form("Hybrid_hadron_passedjets_nowake_1812_%d%d.root", (int)(rL_low*10), (int)(rL_high*10)));
    }
    else{
        cout<<"Gamma"<<endl;
        fvac = TFile::Open("./Hybrid_Vac_hadron_gamma_Vac_1129_23.root");
        fnw = TFile::Open("./Hybrid_hadron_gamma_nowake_1129_23.root");
        fwake = TFile::Open("./Hybrid_hadron_gamma_yeswake_1129_23.root");
    }
    
     cout<<"low  "<<endl;
     TCanvas *cpt = new TCanvas();//for drawing jet pt
     TCanvas *surfcanv = new TCanvas("surfcanv"); 
    if (fvac)
    {
        
        if(ifgamma==true)
        {
           
                cout<<"gamma  "<<endl;
                jet_pt_dat = (TH1D*)fvac->Get("gam_jet_pt_hist");
                eeec_pt_hist = (TH3D*)fvac->Get("eeec_pt_hist");
                eeec_pt_hist->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_r_slice= (TH2D*)eeec_pt_hist->Project3D("yx");
            
                
                jet_pt_dat_wake = (TH1D*)fwake->Get("gam_jet_pt_hist");
                jet_pt_dat_wake->SetName("gam_jet_pt_hist_wake");
                eeec_pt_hist_wake = (TH3D*)fwake->Get("eeec_pt_hist");
                eeec_pt_hist_wake->SetName("eeec_pt_hist_wake");
                eeec_pt_hist_wake->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                
             
                jet_pt_dat_nw = (TH1D*)fnw->Get("gam_jet_pt_hist");
                jet_pt_dat_nw->SetName("gam_jet_pt_hist_nw");
                eeec_pt_hist_nw = (TH3D*)fnw->Get("eeec_pt_hist");
                eeec_pt_hist_nw->SetName("eeec_pt_hist_nw");
                eeec_pt_hist_nw->GetZaxis()->SetRange(pt1_bin_dat,pt2_bin_dat);
                eeec_r_slice_nw = (TH2D*)eeec_pt_hist_nw->Project3D("yx");
      
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
                eeec_r_slice_wake_cop =(TH2D*)eeec_r_slice_wake->Clone("wake_cop");
                eeec_r_slice_wake_cop->Divide(eeec_r_slice);
                
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
                eeec_r_slice_nw_cop =(TH2D*)eeec_r_slice_nw->Clone("nowake_cop");
                eeec_r_slice_nw_cop->Divide(eeec_r_slice);

           
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
        
        TLatex latex, latex1, latex2, latex3, latex4, latex5;
        latex.SetTextColor(kWhite);latex1.SetTextColor(kWhite);latex2.SetTextColor(kWhite);latex3.SetTextColor(kWhite);latex4.SetTextColor(kWhite);
        latex5.SetTextColor(kWhite);
        
        latex.SetTextAlign(33);latex1.SetTextAlign(33);latex2.SetTextAlign(33); latex3.SetTextAlign(33);latex4.SetTextAlign(33);latex5.SetTextAlign(33);
        latex.SetTextSize(0.03);latex1.SetTextSize(0.03);latex2.SetTextSize(0.03);latex3.SetTextSize(0.03);latex4.SetTextSize(0.03);latex5.SetTextSize(0.03);
        
        latex1.SetTextFont(42);
        latex1.SetTextAlign(50);
        latex1.SetNDC ();
        
        latex.SetNDC ();
        latex.SetTextFont(42);
        latex.SetTextAlign(50);
        latex2.SetTextFont(42);
        latex2.SetTextAlign(50);
        latex2.SetNDC ();
 
        latex3.SetTextFont(42);
        latex3.SetTextAlign(50);
        latex3.SetNDC ();
 
        latex4.SetTextFont(42);
        latex4.SetTextAlign(50);
        latex4.SetNDC ();
        
        latex5.SetNDC ();
        latex5.SetTextFont(42);
        latex5.SetTextAlign(50);
        
        const char *str = "#it{p}_{T,jet}";
        const char *str1 = "#it{p}^{ch}_{T,min}";
        const char *str2 = "anti-#it{k}_{T}";   //Writing jet algorithm
        const char *str3 = "#sqrt{s} = 5.02 TeV";        //Writing collision energy
        const char *str4 = "#it{R}_{L}";    //Writing R_L for the 3 point
        const char *str5 = "Uncorrected";
        const char *str6 = "#gamma";
        const char *str7 = "#it{p_{T}^{#gamma}}";
        const char *str8 = "Jet_{#it{Wake=On}}^{#it{Med}}";
        const char *str9 = "Jet_{#it{Wake=Off}}^{#it{Med}}";
        const char *str10 = "Jet^{#it{Vac}}";
        const char *strphi = "#phi";
        const char *strx = "#xi";
        
    
        // Adding a line at 1 to guide the eye
        TLine *line = new TLine(0.0 ,1 ,1,1);
        line->SetLineColorAlpha(kBlue,0.75);
        line->SetLineWidth(2);
        line->SetLineStyle(7);
        
        cout<<"made it here"<<endl;
        
       
    
        
        
        //xxxxxxxxxxxxxxxxxxxxxxxx------------------Handle min & max values of plots------------------------------xxxxxxxxxxxxxxxxxxxxxxxxxx//
        std::vector <double> min_val, min_val_rat;
        std::vector <double> max_val, max_val_rat;
        
        max_val.push_back(eeec_r_slice->GetMaximum()); max_val.push_back(eeec_r_slice_wake->GetMaximum()); max_val.push_back(eeec_r_slice_nw->GetMaximum());
        min_val.push_back(eeec_r_slice->GetMinimum());  min_val.push_back(eeec_r_slice_wake->GetMinimum()); min_val.push_back(eeec_r_slice_nw->GetMinimum());
        
        max_val_rat.push_back(eeec_r_slice_wake_cop->GetMaximum()); max_val_rat.push_back(eeec_r_slice_nw_cop->GetMaximum());
        min_val_rat.push_back(eeec_r_slice_wake_cop->GetMinimum()); min_val_rat.push_back(eeec_r_slice_nw_cop->GetMinimum());
        
        int max = distance(max_val.begin(), max_element(max_val.begin(), max_val.end()));
        int min = distance(min_val.begin(), min_element(min_val.begin(), min_val.end()));
        
        int max_rat = distance(max_val_rat.begin(), max_element(max_val_rat.begin(), max_val_rat.end()));
        int min_rat = distance(min_val_rat.begin(), min_element(min_val_rat.begin(), min_val_rat.end()));
        
        cout<<max_val[0]<<" "<<max_val[1]<<" "<<max_val[2]<< endl;
        cout<<max_val[max]<<endl;
        cout<<min_val[min]<<endl;
        eeec_r_slice->SetMaximum(max_val[max]); eeec_r_slice_wake->SetMaximum(max_val[max]); eeec_r_slice_nw->SetMaximum(max_val[max]);
        eeec_r_slice->SetMinimum(min_val[min]); eeec_r_slice_wake->SetMinimum(min_val[min]); eeec_r_slice_nw->SetMinimum(min_val[min]);
        
        eeec_r_slice_wake_cop->SetMaximum(max_val_rat[max_rat]); eeec_r_slice_nw_cop->SetMaximum(max_val_rat[max_rat]);
        eeec_r_slice_wake_cop->SetMinimum(min_val_rat[min_rat]); eeec_r_slice_nw_cop->SetMinimum(min_val_rat[min_rat]);
        
        
        //xxxxxxxxxxxxxxxxxxxxxxxx-------------------Create Plots-------------------------xxxxxxxxxxxxxxxxxxxxxxxxxx//
        std::stringstream svac;std::stringstream swake;std::stringstream snw;std::stringstream sdiv;std::stringstream sdivnw;
        TCanvas *cvac = new TCanvas("cvac","Vac",900,650);
        cvac->SetRightMargin(0.12);
        TCanvas *cwake = new TCanvas("cwake","Wake",900,650);
        cwake->SetRightMargin(0.12);
        TCanvas *cnw = new TCanvas("cnwake","No Wake",900,650);
        cnw->SetRightMargin(0.12);
        TCanvas *cdiv = new TCanvas("cdiv","WakeVac",900,650);
        cdiv->SetRightMargin(0.12);
        TCanvas *cdivnw = new TCanvas("cdivnw","NoWakeVac",900,650);
        cdivnw->SetRightMargin(0.12);
        
        if(ifvac==true)
        {
            cvac->cd();
            
            eeec_r_slice->SetTitle("");
            eeec_r_slice->GetYaxis()->SetTitle(strphi);
            eeec_r_slice->GetYaxis()->LabelsOption("h");
            eeec_r_slice->GetXaxis()->SetTitle(strx);
            eeec_r_slice->GetXaxis()->LabelsOption("h");
	    eeec_r_slice->GetXaxis()->SetTitleOffset(1.4);
            eeec_r_slice->GetYaxis()->SetTitleOffset(1.4);
            eeec_r_slice->GetXaxis()->SetTitleSize(0.05);
            eeec_r_slice->GetYaxis()->SetTitleSize(0.05);
	    eeec_r_slice->SetMinimum(0.0);
	    eeec_r_slice->SetMaximum(0.011); 
	    
	    TCanvas* surfCanv0 = new TCanvas("surfCanv0", "", 800, 800);
            surfCanv0->cd();
            surfCanv0->SetRightMargin(0.12);
            surfCanv0->SetLeftMargin(0.12);
            surfCanv0->SetBottomMargin(0.1);
            surfCanv0->SetTopMargin(0.1);
           
            eeec_r_slice->SetStats(0);
            eeec_r_slice->Draw("LEGO2");
	    latex3.SetTextColor(kBlack);latex.SetTextColor(kBlack);latex2.SetTextColor(kBlack);latex4.SetTextColor(kBlack);
            latex3.DrawLatex(0.16,0.93 ,Form("Hadrons, Vacuum"));
            latex.DrawLatex(0.58,0.89 ,Form("%.0f GeV/#it{c} < %s < %.0f GeV/#it{c}", low_pt, str, high_pt));
            latex2.DrawLatex(0.58,0.93 ,Form("Full %s jets, #it{R} = 0.8", str2));
            latex4.DrawLatex(0.16,0.85 ,Form("n = 1, %0.1f < %s < %0.1f ", rL_low, str4, rL_high));
             if(ifgamma==true)
            {
                latex4.DrawLatex(0.58,0.85 ,Form("%s - tagged", str6));
               }
            else{
                latex4.DrawLatex(0.16,0.89 ,Form("Hybrid Model (Inclusive sample)"));
            }
	     surfCanv0->SaveAs(Form("vacuumSurface_%d%d.pdf", (int)(rL_low*10), (int)(rL_high*10))); 
        }
        if(ifwake==true)
        {
            cwake->cd();
            eeec_r_slice_wake->SetTitle("");
            eeec_r_slice_wake->GetYaxis()->SetTitle(strphi);
            eeec_r_slice_wake->GetYaxis()->LabelsOption("h");
            eeec_r_slice_wake->GetXaxis()->SetTitle(strx);
            eeec_r_slice_wake->GetXaxis()->LabelsOption("h");
            cout<<"here"<<endl;
            TCanvas* surfCanv2 = new TCanvas("surfCanv2", "", 800, 800);
            surfCanv2->cd();
            surfCanv2->SetRightMargin(0.12);
            surfCanv2->SetLeftMargin(0.12);
            surfCanv2->SetBottomMargin(0.1);
            surfCanv2->SetTopMargin(0.1);

            eeec_r_slice_wake->SetStats(0);
            eeec_r_slice_wake->GetXaxis()->SetTitleOffset(1.4);
            eeec_r_slice_wake->GetYaxis()->SetTitleOffset(1.4);
            eeec_r_slice_wake->GetXaxis()->SetTitleSize(0.05);
            eeec_r_slice_wake->GetYaxis()->SetTitleSize(0.05);
	    eeec_r_slice_wake->SetMinimum(0.0);
	    eeec_r_slice_wake->SetMaximum(0.011); 
            eeec_r_slice_wake->Draw("LEGO2");
          
            latex3.SetTextColor(kBlack);latex.SetTextColor(kBlack);latex2.SetTextColor(kBlack);latex4.SetTextColor(kBlack);
            latex3.DrawLatex(0.16,0.93 ,Form("Hadrons, Wake = ON"));
            latex.DrawLatex(0.58,0.89 ,Form("%.0f GeV/#it{c} < %s < %.0f GeV/#it{c}", low_pt_wake, str, high_pt_wake));
            latex2.DrawLatex(0.58,0.93 ,Form("Full %s jets, #it{R} = 0.8", str2));
            latex4.DrawLatex(0.16,0.85 ,Form("n = 1, %0.1f < %s < %0.1f ", rL_low, str4, rL_high));
             if(ifgamma==true)
            {
                latex4.DrawLatex(0.58,0.85 ,Form("%s - tagged", str6));
               }
            else{
                latex4.DrawLatex(0.16,0.89 ,Form("Hybrid Model (Inclusive Sample)"));
            }
	     surfCanv2->SaveAs(Form("surface_wakeOn_%d%d.pdf", (int)(rL_low*10), (int)(rL_high*10) ));
        }
        if(ifnowake==true)
        {
            cnw->cd();
            eeec_r_slice_nw->SetTitle("");
            eeec_r_slice_nw->GetYaxis()->SetTitle(strphi);
            eeec_r_slice_nw->GetYaxis()->LabelsOption("h");
            eeec_r_slice_nw->GetXaxis()->SetTitle(strx);
            eeec_r_slice_nw->GetXaxis()->LabelsOption("h");
            
            eeec_r_slice_nw->SetStats(0);
            eeec_r_slice_nw->Draw("colz");
	    TCanvas* surfCanv = new TCanvas("surfCanv", "", 800, 800);
	    surfCanv->cd();
	    surfCanv->SetRightMargin(0.12);
	    surfCanv->SetLeftMargin(0.12);
	    surfCanv->SetBottomMargin(0.1);
	    surfCanv->SetTopMargin(0.1);
	    std::cout << eeec_pt_hist_wake->Integral() << std::endl;
	    eeec_r_slice_nw->GetXaxis()->SetTitleOffset(1.4);
	    eeec_r_slice_nw->GetYaxis()->SetTitleOffset(1.4);
	    eeec_r_slice_nw->GetXaxis()->SetTitleSize(0.05);
	    eeec_r_slice_nw->GetYaxis()->SetTitleSize(0.05);
	    eeec_r_slice_nw->SetMinimum(0.0);
	    eeec_r_slice_nw->SetMaximum(0.011); 
	    eeec_r_slice_nw->Draw("LEGO2");
	   
	  
            latex3.SetTextColor(kBlack);latex.SetTextColor(kBlack);latex2.SetTextColor(kBlack);latex4.SetTextColor(kBlack);
            latex3.DrawLatex(0.16,0.93 ,Form("Hadrons, Wake = OFF"));
            latex.DrawLatex(0.58,0.89 ,Form("%.0f GeV/#it{c} < %s < %.0f GeV/#it{c}", low_pt_nw, str, high_pt_nw));
            latex2.DrawLatex(0.58,0.93 ,Form("Full %s jets, #it{R} = 0.8", str2));
            latex4.DrawLatex(0.16,0.85 ,Form("n = 1, %0.1f < %s  < %0.1f ", rL_low, str4, rL_high));
            
            if(ifgamma==true)
            {
                latex4.DrawLatex(0.58,0.85 ,Form("%s - tagged", str6));
               }
            else{
                latex4.DrawLatex(0.16,0.89 ,Form("Hybrid Model (Inclusive Sample)"));
            }
	    surfCanv->SaveAs(Form("wakeoff_surface_%d%d.pdf", (int)(rL_low*10), (int)(rL_high*10)));

        }

         cout<<"max val"<<eeec_r_slice_wake_cop->GetMaximum()<<endl;
         cout<<"min val"<<eeec_r_slice_wake_cop->GetMinimum()<<endl;
        if(ifdivide == true)
        {
            cdiv->cd();
	    TCanvas* surfCanv3 = new TCanvas("surfCanv3", "", 800, 800);
            surfCanv3->cd();
	    surfCanv3->SetRightMargin(0.12);
	    surfCanv3->SetLeftMargin(0.12);
	    surfCanv3->SetBottomMargin(0.1);
	    surfCanv3->SetTopMargin(0.1);
            eeec_r_slice_wake_cop->SetStats(0);
            eeec_r_slice_wake_cop->SetTitle("");
            eeec_r_slice_wake_cop->GetYaxis()->SetTitle(strphi);
            eeec_r_slice_wake_cop->GetYaxis()->LabelsOption("h");
            eeec_r_slice_wake_cop->GetXaxis()->SetTitle(strx);
            eeec_r_slice_wake_cop->GetXaxis()->LabelsOption("h");
	    eeec_r_slice_wake_cop->GetXaxis()->SetTitleOffset(1.4);
            eeec_r_slice_wake_cop->GetYaxis()->SetTitleOffset(1.4);
            eeec_r_slice_wake_cop->GetXaxis()->SetTitleSize(0.05);
            eeec_r_slice_wake_cop->GetYaxis()->SetTitleSize(0.05);
	    
//            cout<<eeec_r_slice_wake->GetEntries()<<endl;
//            cout<<eeec_r_slice_wake_cop->GetMaximum()<<endl;
//            cout<<eeec_r_slice_wake_cop->GetMinimum()<<endl;

            eeec_r_slice_wake_cop->SetMaximum(2.0);
            eeec_r_slice_wake_cop->SetMinimum(0);
           

            eeec_r_slice_wake_cop->Draw("LEGO2");
    
            latex3.SetTextColor(kBlack);latex.SetTextColor(kBlack);latex2.SetTextColor(kBlack);latex4.SetTextColor(kBlack);
            latex3.DrawLatex(0.16,0.93 ,Form("Hadrons, %s / %s",str8,str10));
            latex.DrawLatex(0.58,0.89 ,Form("%.0f GeV/#it{c} < %s < %.0f GeV/#it{c}", low_pt_wake, str, high_pt_wake));
            latex2.DrawLatex(0.58,0.93 ,Form("Full %s jets, #it{R} = 0.8", str2));
            latex4.DrawLatex(0.16,0.85 ,Form("n = 1, %0.1f < %s < %0.1f ", rL_low,  str4, rL_high));
            if(ifgamma==true)
            {
                latex4.DrawLatex(0.58,0.85 ,Form("%s - tagged", str6));
            }
            else{
                latex4.DrawLatex(0.16,0.89 ,Form("Hybrid Model (Inclusive Sample)"));
            }
	    surfCanv3->SaveAs(Form("wakeovernowake_%d%d.pdf", (int)(rL_low*10), (int)(rL_high*10))); 
        }
        
      
        if(ifdividenw == true)
        {
            cdivnw->cd();
	    TCanvas* surfCanv4 = new TCanvas("surfCanv4", "", 800, 800);
            surfCanv4->cd();
            surfCanv4->SetRightMargin(0.12);
            surfCanv4->SetLeftMargin(0.12);
            surfCanv4->SetBottomMargin(0.1);
            surfCanv4->SetTopMargin(0.1);
            eeec_r_slice_nw_cop->SetTitle(0);
            eeec_r_slice_nw_cop->SetStats(0);
            eeec_r_slice_nw_cop->SetTitle("");
            eeec_r_slice_nw_cop->GetYaxis()->SetTitle(strphi);
            eeec_r_slice_nw_cop->GetYaxis()->LabelsOption("h");
            eeec_r_slice_nw_cop->GetXaxis()->SetTitle(strx);
            eeec_r_slice_nw_cop->GetXaxis()->LabelsOption("h");
	    eeec_r_slice_nw_cop->GetXaxis()->SetTitleOffset(1.4);
            eeec_r_slice_nw_cop->GetYaxis()->SetTitleOffset(1.4);
            eeec_r_slice_nw_cop->GetXaxis()->SetTitleSize(0.05);
            eeec_r_slice_nw_cop->GetYaxis()->SetTitleSize(0.05);
    
            cout<<eeec_r_slice_nw->GetEntries()<<endl;
    
           eeec_r_slice_nw_cop->SetMaximum(2.0);
	   eeec_r_slice_nw_cop->SetMinimum(0.0); 
            
            eeec_r_slice_nw_cop->Draw("LEGO2");
	    latex3.SetTextColor(kBlack);latex.SetTextColor(kBlack);latex2.SetTextColor(kBlack);latex4.SetTextColor(kBlack);
            latex3.DrawLatex(0.16,0.93 ,Form("Hadrons, %s / %s",str9,str10));
            latex.DrawLatex(0.58,0.89 ,Form("%.0f GeV/#it{c} < %s < %.0f GeV/#it{c}", low_pt_nw, str, high_pt_nw));
            latex2.DrawLatex(0.58,0.93 ,Form("%s , #it{R} = 0.8", str2));
            latex4.DrawLatex(0.16,0.85 ,Form("n = 1, %0.1f < %s < %0.1f ", rL_low, str4, rL_high));
            if(ifgamma==true)
            {
                latex4.DrawLatex(0.58,0.85 ,Form("%s - tagged", str6));
            }
            else{
                latex4.DrawLatex(0.16,0.89 ,Form("Hybrid Model (Inclusive Sample)"));
            }
	    surfCanv4->SaveAs(Form("nowake_over_vacuum_%d%d.pdf", (int)(rL_low*10), (int)(rL_high*10))); 
        }
    
    if(ifsave)
    {
        svac << "EEEC_V_34" << "_" << low_pt <<  "to"<< high_pt << "_R8_n1_Dec18.pdf";
        cvac->SaveAs(svac.str().c_str());
        swake << "EEEC_W_34" << "_" <<  low_pt_wake  <<  "to"<< high_pt_wake << "_R8_n1_Dec18.pdf";
        cwake->SaveAs(swake.str().c_str());
        snw << "EEEC_NW_34" << "_" <<  low_pt_nw  <<  "to"<< high_pt_nw << "_R8_n1_Dec18.pdf";
        cnw->SaveAs(snw.str().c_str());
        
        if(ifdivide == true)
        {
        sdiv <<"EEEC_WV_34" << "_" << low_pt <<  "to"<< high_pt << "_R8_n1_Dec18.pdf";
        cdiv->SaveAs(sdiv.str().c_str());
        }
        if(ifdividenw == true)
        {
        sdivnw <<"EEEC_NWV_34" << "_" << low_pt <<  "to"<< high_pt << "_R8_n1_Dec18.pdf";
        cdivnw->SaveAs(sdivnw.str().c_str());
        }
    }
    
     if (!eeec_pt_hist)
            printf("No such histogram found!\n");
    }
    else
        printf("No such file found!\n");
}
