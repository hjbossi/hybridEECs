#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <iterator>
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include "TColor.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

using namespace Pythia8;
using namespace fastjet;
using namespace std;

Double_t delR(fastjet::PseudoJet ps1, fastjet::PseudoJet ps2)
{
    double dphi = abs(ps1.phi()-ps2.phi());
    if(dphi>TMath::Pi()){dphi = (2.*TMath::Pi() - dphi);}
    
    double deta = std::fabs(ps1.eta() - ps2.eta());
    Double_t delR = std::sqrt(dphi*dphi + deta*deta);
    
    return delR;
}

int main(int argc, char* argv[]) {
    
    //Its useful to have options for what I am trying to run
    std::string argu1 = argv[1]; //parton or hadron level
    std::string argu2 = argv[2]; //gamma jets or no
    std::string argu3 = argv[3]; //3 point corr or no
    std::string argu4 = argv[4]; //name of the input file
    std::string argu5 = argv[5]; //wake or no wake
    char* argu6 = argv[6]; //adds date to file
    std::string argu7 = argv[7]; //R_low
    std::string argu8 = argv[8]; //R_high
    std::string argu9 = argv[9]; //power n of weighting
    //   std::string argu10 = argv[10]; //jet radius
    
    
    std::string parton = "parton";
    std::string hadron = "hadron";
    std::string gamma = "gamma";
    std::string passedjets = "passedjets";
    std::string gammaNoJet = "nojet";
    std::string three = "3"; //3 point corr
    std::string vac = "Vac"; //vacuum or no
    std::string nowake = "nowake"; //wake  no
    std::string yeswake = "yeswake"; //wake yes or no
    
    TH3D *R_sizes;
    TH3D *R_sizes_www;
    TH3D *R_sizes_jjw;
    TH3D *Rdistribution;
    TH3D *Rdistribution_wt;
    TH3D *R_sizes_gam;
    TH3D *R_sizes_gam1;
    TH3D *R_sizes_gam2;
    TH3D *R_sizes_gam3;
    TH2D *eec_pt_hist;
    TH2D *gam_pt_jet_pt;
    TH2D *eec_pt_hist_gam;
    TH1D *eec_hist;
    TH1D *e3c_hist;
    TH2D *e3c_pt_hist;
    TH2D *e3c_pt_hist1;
    TH2D *e3c_pt_hist2;
    TH2D *e3c_pt_hist_gam;
    TH1 *jet_pt_hist;
    TH1 *gam_jet_pt_hist;
    TFile *fout_nom;
    TH2D *e3c_pt_hist_gam_jet;
    TH3D *neg_map;
    TH3D *pos_map;
    TH3D *jet_map;
    
    Double_t from = -4;
    Double_t to = 0;
    Int_t bins = 100;
    Double_t width = (to-from)/bins;
    Double_t new_bins[101] = {};
    for (Int_t i = 0; i <= bins; i++)
    {
        new_bins[i] = TMath::Power(10.0, from + i * width);
    }
    
    //xi bins
    Double_t fromR = 0;
    Double_t toR = 1;
    Int_t binsR = 20;
    Double_t widthR = (toR-fromR)/binsR;
    Double_t new_binsR[21] = {};
    for (Int_t i = 0; i <= binsR; i++)
    {
        new_binsR[i] = (fromR + i * widthR);
    }
    
    //bins for ratio of R distribution of particles in jjw case
    Double_t fromRrat = 0;
    Double_t toRrat = 5;
    Int_t binsRrat = 100;
    Double_t widthRrat = (toRrat-fromRrat)/binsRrat;
    Double_t new_binsRrat[101] = {};
    for (Int_t i = 0; i <= binsRrat; i++)
    {
        new_binsRrat[i] = (fromRrat + i * widthRrat);
    }
    
    
    //phi bins
    Double_t fromRs = 0;
    Double_t toRs = M_PI/2;
    Int_t binsRs = 20;
    Double_t widthRs = (toRs-fromRs)/binsRs;
    Double_t new_binsRs[21] = {};
    for (Int_t i = 0; i <= binsRs; i++)
    {
        new_binsRs[i] = (fromRs + i * widthRs);
    }
    
    //jet pT bins
    Double_t from_const = 40;
    Double_t to_const = 1000;
    Int_t bins_const = 48;
    Double_t width_const = (to_const-from_const)/bins_const;
    Double_t new_bins_const[49] = {};
    for (Int_t i = 0; i <= bins_const; i++)
    {
        new_bins_const[i] = (from_const + i * width_const);
    }
    
    
    
    jet_pt_hist = new TH1D("jet_pt_hist", "Jet Pt", 48, 40, 1000);
    
    gam_pt_jet_pt = new TH2D("gamma_pt_jet_pt", "Gamma pt vs jet pt", 50, 100, 1000,48,40,1000);
    
    gam_jet_pt_hist = new TH1D("gam_jet_pt_hist", "Gam-Jet Pt", 48, 40, 1000);
    
    eec_hist = new TH1D("eec_hist","EEC", 100, new_bins);
    
    eec_pt_hist = new TH2D("eec_pt_hist", "EEC and jet_pt 2D", 100, new_bins, 48, 40, 1000);
    
    eec_pt_hist_gam = new TH2D("eec_pt_hist_gam", "EEC and jet_pt 2D gam", 100, new_bins, 48, 40, 1000);
    
    R_sizes = new TH3D("eeec_pt_hist", "EEEC and jet_pt 3D", 20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    
    R_sizes_www = new TH3D("eeec_pt_hist_www", "EEEC and jet_pt 3D", 20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    
    R_sizes_jjw = new TH3D("eeec_pt_hist_jjw", "EEEC and jet_pt 3D", 20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    
    Rdistribution = new TH3D("Rdistribution", "R distribution of jjw", 20, new_binsR, 100, new_binsRrat, 100, new_binsRrat);
    
    // Rdistribution_wt = new TH3D("Rdistribution_wt", "EEEC and jet_pt 3D", 20, new_binsR, 20, new_binsR, 20, new_binsR);
    
    Rdistribution_wt = new TH3D("Rdistribution_wt", "EEEC and jet_pt 3D", 20, new_binsR, 100, new_binsRrat, 100, new_binsRrat);
    
    
    R_sizes_gam = new TH3D("eeec_pt_hist_gam", "EEEC and jet_pt 3D gam", 20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    
    R_sizes_gam1 = new TH3D("eeec_pt_hist_gam1", "EEEC and jet_pt 3D gam 20", 20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    
    R_sizes_gam2 = new TH3D("eeec_pt_hist_gam2", "EEEC and jet_pt 3D gam 50", 20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    
    R_sizes_gam3 = new TH3D("eeec_pt_hist_gam3", "EEEC and jet_pt 3D gam 70", 20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    
    neg_map = new TH3D("neg_map", "phi,eta,pt map of negatives", 100, -3.14, 3.14, 200, -3, 3, 200, 0, 200);
    
    pos_map = new TH3D("pos_map", "phi,eta,pt map of positives", 100, -3.14, 3.14, 200, -3, 3, 200, 0, 200);
    
    jet_map = new TH3D("jeT_map", "phi,eta,pt map of jet", 100, -3.14, 3.14, 200, -3, 3, 200, 0, 200);
    
    
    if(argu3.find(three) != std::string::npos)
    {
        e3c_hist = new TH1D("e3c_hist","E3C", 100, new_bins);
        
        e3c_pt_hist = new TH2D("e3c_pt_hist", "E3C and jet_pt 2D", 100, new_bins, 48, 40, 1000);
        
        e3c_pt_hist1 = new TH2D("e3c_pt_hist1", "E3C and jet_pt 2D", 100, new_bins, 48, 40, 1000);
        
        e3c_pt_hist2 = new TH2D("e3c_pt_hist2", "E3C and jet_pt 2D", 100, new_bins, 48, 40, 1000);
        
        e3c_pt_hist_gam = new TH2D("e3c_pt_hist_gam", "E3C and jet_pt 2D", 100, new_bins, 48, 40, 1000);
        
        // e3c_pt_hist_gam_jet = new TH2D("e3c_pt_hist_gam_jet", "E3C and gam_pt and jet_pt 3D", 100, new_bins, 48, new_bins_const, 48, new_bins_const));
    }
    
    const char* charPtrIn = argu4.c_str();
    TFile *f = new TFile(charPtrIn);//import the root file with the trees
    if(!f){exit(0);}
    
    
    std::stringstream fname;
    if(argu4.find(vac) != std::string::npos)
    {
        fname << "Hybrid_Vac_" <<  argu1  <<  "_"<< argu2 <<"_"<<argu5<<"_"<<argu6<<".root";
    }
    else
    {
        fname << "Hybrid_" <<  argu1  <<  "_"<< argu2 << "_"<<argu5<<"_"<<argu6<<".root";
    }
    std::string fname_str =  fname.str();
    const char* charPtr = fname_str.c_str();
    fout_nom = TFile::Open(charPtr,"recreate");
    
    TTree *t1 = (TTree*)f->Get("T");
    Long64_t event_id;
    Double_t wt;
    Double_t crossx;
    Double_t px, py, pz, m;
    Int_t id, label;
    t1->SetBranchAddress("px",&px);
    t1->SetBranchAddress("py",&py);
    t1->SetBranchAddress("pz",&pz);
    t1->SetBranchAddress("m",&m);
    t1->SetBranchAddress("Event",&event_id);
    t1->SetBranchAddress("id",&id);
    t1->SetBranchAddress("label",&label);
    t1->SetBranchAddress("Weight",&wt);
    
    
    //// Select common parameters for FastJet analyses.
    int    power   = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT.
    double R       = 0.8;    // Jet size.
    double pTMin   = 40;    // Min jet pT.
    double pTMax   = 1000;    // Max jet pT.
    float track_pt = 1.0 ; //Min track pT
    double R_low = std::stod(argu7);
    double R_high = std::stod(argu8);
    double power_wt = std::stod(argu9);
    // double R = std::stod(argu10); //jet size
    //  double etaMax  = 0.9;    // Pseudorapidity range of detector.
    
    // cout<<R_low<<"and "<<R_high<<endl;
    
    // Defining a class to help store input particle PIDs to access later to tag heavy quark hadrons
    class MyInfo: public fastjet::PseudoJet::UserInfoBase //class name my info, inherits from Pseudojet::UserInfoBase
    {
    public:
        MyInfo(int id, int label, float pT): _pdg_id(id),_label_id(label),_pT(pT){} //MyInfo(int id) is a parametrized class constructor. Has to have the same name as the class but no return type
        int pdg_id() const {return _pdg_id;}//here we are defining a member function pdg_id() of the class that returns an integer value
        int label_id() const {return _label_id;}
        float pT() const {return _pT;}
    protected:
        int _pdg_id;
        int _label_id;
        float _pT;
    };
    
    //Define a class to store gamma pT so that you can change jet pT to gamma pT
    class MyInfoPt: public fastjet::PseudoJet::UserInfoBase //class name my info, inherits from Pseudojet::UserInfoBase
    {
    public:
        MyInfoPt(double gam_pt): _gam_pt(gam_pt){} //MyInfo(int id) is a parametrized class constructor. Has to have the same name as the class but no return type
        double gam_pt() const {return _gam_pt;}//here we are defining a member function pdg_id() of the class that returns an integer value
    protected:
        double _gam_pt;
    };
    
    //Define a class to store gamma pT so that you can change jet pT to gamma pT
    class MyInfoPtJet: public fastjet::PseudoJet::UserInfoBase //class name my info, inherits from Pseudojet::UserInfoBase
    {
    public:
        MyInfoPtJet(double jet_pt_included): _jet_pt_included(jet_pt_included){} //MyInfo(int id) is a parametrized class constructor. Has to have the same name as the class but no return type
        double jet_pt_included() const {return _jet_pt_included;}//here we are defining a member function pdg_id() of the class that returns an integer value
    protected:
        double _jet_pt_included;
    };
    
    // // Defining a class to help store input particle labels to access wake or jet particles
    // class MyInfoLabel: public fastjet::PseudoJet::UserInfoBase //class name my info, inherits from Pseudojet::UserInfoBase
    // {
    // public:
    //     MyInfoLabel(int label): _label_id(label){}
    //     int label_id() const {return _label_id;}
    // protected:
    //     int _label_id;
    // };
    ////Variables
    // double weight,crossx, x, y, mass, label, index, mom_x, mom_y,mom_z;
    // int eventID, pdg_id;
    Vec4 particleTemp;
    // vector<double> event_id, x_pos, y_pos, px, py, pz,m, unique_event_id, pdg_ID ;
    vector<double> unique_event_id;
    
    cout <<unique_event_id.size()<<" after set"<<endl;
    
    std::vector <fastjet::PseudoJet> fjInputs;
    std::vector <double> weights;
    int n = 0;
    
    // Set the condition
    Int_t conditionValue = 0;
    
    // // Index variable to be incremented so you can check if all entries are read
    Int_t a = 0;
    
    // Loop over entries in the tree
    for (Int_t i = 0; i < t1->GetEntries(); ++i) {
        t1->GetEntry(i);
        
        // Check the condition (e.g., event_id equals the condition value)
        if (event_id == conditionValue) {
            // cout<< event_id<<" here "<<wt<<endl;
            //add particles to fjinput list
            if(label != -2) //This excludes the hard parton that scatters from the jet
            {
                if(argu4.find(vac) != std::string::npos) //in case of vac include everything
                {
                    //  cout<<"a "<<a<<" label "<<label<<endl;
                    double E = sqrt((m*m)+(px*px+py*py+pz*pz));
                    PseudoJet pTemp(px,py,pz,E);
                    float pT = pTemp.pt();
                    pTemp.set_user_info(new MyInfo(id,label,pT));
                    fjInputs.push_back(pTemp);
                    
                }
                else //in case of hybrid data
                {
                    if(argu5.find(nowake) != std::string::npos && label != 1) //do not include the wake particles
                    {
                        //  cout<<"no wake"<<endl;
                        double E = sqrt((m*m)+(px*px+py*py+pz*pz));
                        PseudoJet pTemp(px,py,pz,E);
                        float pT = pTemp.pt();
                        pTemp.set_user_info(new MyInfo(id,label,pT));
                        fjInputs.push_back(pTemp);
                        
                    }
                    else //include all particles from jets + wake
                    {
                        // if(argu5.find(yeswake) != std::string::npos && label !=2) //excluding negative wake particles
                        if(argu5.find(yeswake) != std::string::npos) //for pT subtracted files
                        {
                            //  cout<<"yes wake"<< " " <<a<<endl;
                            double E = sqrt((m*m)+(px*px+py*py+pz*pz));
                            PseudoJet pTemp(px,py,pz,E);
                            
                            float pT = pTemp.pt();
                            pTemp.set_user_info(new MyInfo(id,label,pT));
                            
                            
                            fjInputs.push_back(pTemp);
                        }
                    }
                }
            }
        }
        else {
            //When the condition is not met, process the jets because now list of particles is complete
            t1->GetEntry(i-1);//define this so you can get the event weight of the right event (which is the previous event)
            //  cout<<i<<" and "<< i-1 << " and "<<event_id<<" and here "<<wt<<endl;
            ////Jet definition
            vector <fastjet::PseudoJet> inclusiveJets, sortedJets,passedJets,constit_i,constit_j;
            vector <fastjet::PseudoJet> gammaJets, EtPhoton, PhotonParticle;
            fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, R, power);
            fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
            inclusiveJets = clustSeq.inclusive_jets(pTMin);
            sortedJets    = sorted_by_pt(inclusiveJets);
            
            ////Jet selections
            double rap_cut = 2.0;//Jet axis rapidity cut from Dani
            //  double rap_cut = 0.9 - R;//Jet axis rapidity cut
            //  double rap_cut = 2.4 - R;//Jet axis rapidity cut
            fastjet::Selector select_pt_min = fastjet::SelectorPtMin(pTMin);
            fastjet::Selector select_pt_max = fastjet::SelectorPtMax(pTMax);
            fastjet::Selector select_rap = fastjet::SelectorAbsRapMax(rap_cut);//For setting min jet pt and rapidity cut (rapidity cut should be (0.9-jetR)
            fastjet::Selector select_all = select_pt_min && select_rap && select_pt_max;//order of operation doesn't matter in how its defined here
            passedJets = select_all(sortedJets);
            
            vector <fastjet::PseudoJet> constit;
            vector<double> R_dist, constit_pt, R_dist_jjw, R_dist_www,R_dist_new;
            double eee_jsm,max_pt_part_id;
            double jet_pt_included;
            double p_x, p_y, p_z, p_t;
            double Rj1j2;
           
            if(argu2.find(passedjets) != std::string::npos)
            {
                for(int i=0; i< int (passedJets.size()); ++i)
                {
                    constit = sorted_by_pt(passedJets[i].constituents());
                    // constit = passedJets[i].constituents();
                    double jet_pt = passedJets[i].pt();
                    jet_pt_hist->Fill(jet_pt,wt);
                    
                    for(int j=0; j < int(constit.size());j++)
                    {
                        // if(j==0){
                        //     for(int j=0; j < int(constit.size());j++){
                        //         constit_pt.push_back(constit[j].pt());
                        //     }
                        
                        //      max_pt_part_id =  distance(constit_pt.begin(), max_element(constit_pt.begin(), constit_pt.end()));
                        
                        //     }
                        
                        if(argu3.find(three) != std::string::npos)
                        {
                            
                            for (int s=0; s<int(constit.size());s++)
                            {
                                if(s==j) continue;
                                
                                
                                for( int m=0; m<j && m<s; m++)
                                {
                                    if(s>j) continue;
                                    
                                    if(argu5.find(yeswake) != std::string::npos)
                                    {
                                        //-------------------------------jet-jet-wake correlation-----------------------------------------
                                        if((constit[j].user_info<MyInfo>().label_id() == 0 && constit[s].user_info<MyInfo>().label_id() == 0 && constit[m].user_info<MyInfo>().label_id() == 1) ||
                                           (constit[j].user_info<MyInfo>().label_id() == 1 && constit[s].user_info<MyInfo>().label_id() == 0  && constit[m].user_info<MyInfo>().label_id() == 0) ||
                                           (constit[j].user_info<MyInfo>().label_id() == 0 && constit[s].user_info<MyInfo>().label_id() == 1 && constit[m].user_info<MyInfo>().label_id() == 0)
                                           )
                                        {
                                            
                                            if(constit[j].user_info<MyInfo>().label_id() == 0 && constit[s].user_info<MyInfo>().label_id() == 0 && constit[m].user_info<MyInfo>().label_id() == 1)
                                            {
                                                Rj1j2=constit[j].delta_R(constit[s]);
                                            }
                                            else if(constit[j].user_info<MyInfo>().label_id() == 1 && constit[s].user_info<MyInfo>().label_id() == 0  && constit[m].user_info<MyInfo>().label_id() == 0)
                                            {
                                                Rj1j2=constit[s].delta_R(constit[m]);
                                            }
                                            else
                                            {
                                                Rj1j2=constit[j].delta_R(constit[m]);
                                            }
                                            
                                            //                         // if(constit[j].pt()<1 || constit[s].pt()<1 || constit[m].pt()<1) continue;
                                            if(Rj1j2<R_low || Rj1j2>R_high) continue; //make sure the jet-jet distance is between R_low and R_high, don't enforce that it is the longest distance
                                            
                                            //----------------------------wake-wake-wake correlation-----------------------------------------
                                            // if(constit[j].user_info<MyInfo>().label_id() == 1 && constit[s].user_info<MyInfo>().label_id() == 1 && constit[m].user_info<MyInfo>().label_id() == 1)
                                            // {
                                            //     double eee_jsmnew = ((6*constit[j].pt()*constit[s].pt()*constit[m].pt())/(pow(jet_pt,3)));
                                            //     double Rjs = constit[j].delta_R(constit[s]);
                                            //     double Rjm = constit[j].delta_R(constit[m]);
                                            //     double Rsm = constit[s].delta_R(constit[m]);
                                            
                                            //     R_dist_www.push_back(Rjs);
                                            //     R_dist_www.push_back(Rjm);
                                            //     R_dist_www.push_back(Rsm);
                                            
                                            //     int max_R_www = distance(R_dist_www.begin(), max_element(R_dist_www.begin(), R_dist_www.end()));//pick the longest side to compute the correlators with
                                            //     int min_R_www = distance(R_dist_www.begin(), min_element(R_dist_www.begin(), R_dist_www.end()));
                                            //     int mid_R_www;
                                            //     for(int q=0; q<3; q++)
                                            //         {
                                            //             if(q!=max_R_www && q!=min_R_www){mid_R_www==q;}
                                            //         }
                                            
                                            //     double x_www = R_dist_www[min_R_www]/R_dist_www[mid_R_www];
                                            //     double diff_www = (R_dist_www[max_R_www]-R_dist_www[mid_R_www]);
                                            //     double rat_www = (diff_www*diff_www)/ (R_dist_www[min_R_www]*R_dist_www[min_R_www]);
                                            //     double phi_www = asin(sqrt(1-rat_www));
                                            
                                            //     if(R_low<R_dist_www[max_R_www] && R_dist_www[max_R_www]<R_high)
                                            //     {
                                            //         R_sizes_www->Fill(x_www,phi_www,jet_pt,eee_jsmnew*wt);
                                            //     }
                                            // }
                                            
                                            
                                            //-------------------------------basic computation-----------------------------
                                            if(power_wt == 1.5)
                                            {
                                                eee_jsm = ((6*(pow(constit[j].pt(),1.5))*(pow(constit[s].pt(),1.5))*(pow(constit[m].pt(),1.5)))/(pow(jet_pt,4.5)));
                                            }
                                            else if(power_wt == 0.5)
                                            {
                                                eee_jsm = ((6*(pow(constit[j].pt(),0.5))*(pow(constit[s].pt(),0.5))*(pow(constit[m].pt(),0.5)))/(pow(jet_pt,1.5)));
                                            }
                                            else
                                            {
                                                eee_jsm = ((6*constit[j].pt()*constit[s].pt()*constit[m].pt())/(pow(jet_pt,3)));
                                            }
                                            
                                           
                                            
                                            double deltaR_js = constit[j].delta_R(constit[s]);
                                            double deltaR_jm = constit[j].delta_R(constit[m]);
                                            double deltaR_sm = constit[s].delta_R(constit[m]);
                                            
                                            
                                            R_dist.push_back(deltaR_js);
                                            R_dist.push_back(deltaR_jm);
                                            R_dist.push_back(deltaR_sm);
                                            
                                            int max_R = distance(R_dist.begin(), max_element(R_dist.begin(), R_dist.end()));//pick the longest side to compute the correlators with
                                            int min_R = distance(R_dist.begin(), min_element(R_dist.begin(), R_dist.end()));//pick the shortest side to compute the correlators with
                                            
                                            int mid_R;
                                            for(int q=0; q<3; q++)
                                            {
                                                if(q!=max_R && q!=min_R){mid_R==q;}
                                            }
                                            
                                            
                                            double x = R_dist[min_R]/R_dist[mid_R];
                                            double diff = (R_dist[max_R]-R_dist[mid_R]);
                                            double rat = (diff*diff)/ (R_dist[min_R]*R_dist[min_R]);
                                            double phi = asin(sqrt(1-rat));
                                            
                                            
                                            // if(R_low<R_dist[max_R] && R_dist[max_R]<R_high)
                                            // {
                                            R_sizes->Fill(x,phi,jet_pt,eee_jsm*wt);
                                            // }
                                           
                                            
                                            // e3c_pt_hist->Fill(R_dist[max_R],jet_pt,eee_jsm*wt);
                                            
                                            R_dist.clear();
                                            //                         R_dist_jjw.clear();
                                            //                         R_dist_www.clear();
                                            
                                        }
                                    }
                                    else
                                    {
                                        Rj1j2=constit[s].delta_R(constit[m]);
                                        if(Rj1j2<R_low || Rj1j2>R_high) continue; //make sure the jet-jet distance is between R_low and R_high, s & m are the higher pT particles because my code only allows for jsm where j is the highest index
                                        
                                        if(power_wt == 1.5)
                                        {
                                            eee_jsm = ((6*(pow(constit[j].pt(),1.5))*(pow(constit[s].pt(),1.5))*(pow(constit[m].pt(),1.5)))/(pow(jet_pt,4.5)));
                                        }
                                        else if(power_wt == 0.5)
                                        {
                                            eee_jsm = ((6*(pow(constit[j].pt(),0.5))*(pow(constit[s].pt(),0.5))*(pow(constit[m].pt(),0.5)))/(pow(jet_pt,1.5)));
                                        }
                                        else
                                        {
                                            eee_jsm = ((6*constit[j].pt()*constit[s].pt()*constit[m].pt())/(pow(jet_pt,3)));
                                        }
                                        
                                        double deltaR_js = constit[j].delta_R(constit[s]);
                                        double deltaR_jm = constit[j].delta_R(constit[m]);
                                        double deltaR_sm = constit[s].delta_R(constit[m]);
                                        
                                        
                                        R_dist.push_back(deltaR_js);
                                        R_dist.push_back(deltaR_jm);
                                        R_dist.push_back(deltaR_sm);
                                        
                                        int max_R = distance(R_dist.begin(), max_element(R_dist.begin(), R_dist.end()));//pick the longest side to compute the correlators with
                                        int min_R = distance(R_dist.begin(), min_element(R_dist.begin(), R_dist.end()));//pick the shortest side to compute the correlators with
                                        
                                        int mid_R;
                                        for(int q=0; q<3; q++)
                                        {
                                            if(q!=max_R && q!=min_R){mid_R==q;}
                                        }
                                        
                                        
                                        double x = R_dist[min_R]/R_dist[mid_R];
                                        double diff = (R_dist[max_R]-R_dist[mid_R]);
                                        double rat = (diff*diff)/ (R_dist[min_R]*R_dist[min_R]);
                                        double phi = asin(sqrt(1-rat));
                                        
                                        
                                        
                                        // if(R_low<R_dist[max_R] && R_dist[max_R]<R_high)
                                        // {
                                        R_sizes->Fill(x,phi,jet_pt,eee_jsm*wt);
                                        // }
                                        
                                        R_dist.clear();
                                    }
                                }//close m loop
                            }//close s loop
                        } //close if loop for 3 point corr
                    } //close j loop
                } //close i loop
            }//close if loop for passed jets
            
            if(argu2.find(gamma) != std::string::npos)
            {
                //This condition should hold true for any gamma-jet pair because fastjet should see the recoil photon as a jet
                if(passedJets.size()>=2)
                {
                    
                    //         ///Finding the direct photon
                    //          ///The cone of the direct photon is 0.4. It can only have a certain amount of transverse energy in this cone to be considered isolated.
                    //          //It also has certain phase space cuts to be considered isolated (delta_phi)
                    double transE=0;
                    double delPhi_iso = 0;
                    
                    //Loop over all jets in an event
                    for(int i=0; i< int(passedJets.size()); i++)
                    {
                        
                        //Loop over all the particles in an event
                        for(int j=0; j< int (fjInputs.size()); ++j)
                        {
                            
                            //Find a photon that is "opposite" to a jet.
                            //Opposite here means: distance between photon and jet > 0.4. Angle between photon and jet > (2/3)pi
                            if(fjInputs[j].user_info<MyInfo>().pdg_id() == 22 && fjInputs[j].pt()>100 && fjInputs[j].pt()<1000  && fjInputs[j].rap()<2.)
                            {
                                double delR_iso = fjInputs[j].delta_R(passedJets[i]);
                                double delPhi = abs(fjInputs[j].phi()-passedJets[i].phi());
                                
                                if (delPhi>pi){delPhi_iso = 2*pi - delPhi;}
                                else{delPhi_iso = delPhi;}
                                if (delR_iso > 0.4 && delPhi_iso > ((2*pi)/3) )
                                {
                                    //Here I check if the transverse energy within the photon cone is less than 5GeV. This is an isolation cut.
                                    //First I find all the particles that lie within the photon cone of radius 0.4.
                                    //Then I sum the transverse energy of all particles that are within the photon cone.
                                    for(int k=0; k< int (fjInputs.size()); ++k)
                                    {
                                        if (j!=k && fjInputs[j].delta_R(fjInputs[k])<0.4)
                                        {
                                            // cout<<"photon cone "<<fjInputs[j].delta_R(fjInputs[k])<<endl;
                                            // cout<<"particle energy "<<fjInputs[k].Et()<<endl;
                                            transE = transE + fjInputs[k].Et();
                                        }
                                    }
                                    
                                    //If the Et_iso condition is satisfied, I have found a recoil jet to the direct photon. Add this to the list of jets
                                    if(transE <= 5.)
                                    {
                                        gam_pt_jet_pt->Fill(fjInputs[j].pt(),passedJets[i].pt());
                                        double gam_pt=fjInputs[j].pt();
                                        //  cout<<"value "<<gam_pt<<endl;
                                        passedJets[i].set_user_info(new MyInfoPt(gam_pt));
                                        gammaJets.push_back(passedJets[i]);
                                        
                                    }
                                    
                                    transE = 0;
                                }
                            }
                        }
                    }
                }
                // if(gammaJets.size()==0) break;
                for(int i=0; i< int (gammaJets.size()); ++i)
                {
                    
                    constit = gammaJets[i].constituents();
                    double quench_jet_pt = gammaJets[i].pt();
                    double jet_pt = gammaJets[i].user_info<MyInfoPt>().gam_pt();
                    double aj = quench_jet_pt/jet_pt;
                    gam_jet_pt_hist->Fill(jet_pt,wt);
                    // cout<<aj<<endl;
                    
                    for(int j=0; j < int(constit.size());j++)
                    {
                        // if(j==0){
                        //     for(int j=0; j < int(constit.size());j++){
                        //         constit_pt.push_back(constit[j].pt());
                        //     }
                        
                        //      max_pt_part_id =  distance(constit_pt.begin(), max_element(constit_pt.begin(), constit_pt.end()));
                        
                        //     }
                        if(argu3.find(three) != std::string::npos)
                        {
                            
                            for (int s=0; s<int(constit.size());s++)
                            {
                                if(s==j) continue;
                                
                                //   for( int m=0; m!=j && m!=s; m++)
                                for( int m=0; m<j && m<s; m++)
                                {
                                    if(s>j) continue;
                                    // if(constit[j].pt()<1 || constit[s].pt()<1 || constit[m].pt()<1) continue;
                                    //To compute correlation over wake particles
                                    // if(constit[j].user_info<MyInfo>().label_id() == 0 || constit[s].user_info<MyInfo>().label_id() == 0 || constit[m].user_info<MyInfoLabel>().label_id() == 0){continue;}
                                    
                                    //jet-wake=wake correlation
                                    
                                    // if(j == max_pt_part_id && constit[s].user_info<MyInfoLabel>().label_id() == 1 && constit[m].user_info<MyInfoLabel>().label_id() == 1)
                                    // {
                                    if(power_wt == 1.5)
                                    {
                                        eee_jsm = ((6*(pow(constit[j].pt(),1.5))*(pow(constit[s].pt(),1.5))*(pow(constit[m].pt(),1.5)))/(pow(jet_pt,4.5)));
                                    }
                                    else if(power_wt == 0.5)
                                    {
                                        eee_jsm = ((6*(pow(constit[j].pt(),0.5))*(pow(constit[s].pt(),0.5))*(pow(constit[m].pt(),0.5)))/(pow(jet_pt,1.5)));
                                    }
                                    else
                                    {
                                        eee_jsm = ((6*constit[j].pt()*constit[s].pt()*constit[m].pt())/(pow(jet_pt,3)));
                                    }
                                    
                                    
                                    double deltaR_js = constit[j].delta_R(constit[s]);
                                    double deltaR_jm = constit[j].delta_R(constit[m]);
                                    double deltaR_sm = constit[s].delta_R(constit[m]);
                                    
                                    
                                    R_dist.push_back(deltaR_js);
                                    R_dist.push_back(deltaR_jm);
                                    R_dist.push_back(deltaR_sm);
                                    
                                    int max_R = distance(R_dist.begin(), max_element(R_dist.begin(), R_dist.end()));//pick the longest side to compute the correlators with
                                    int min_R = distance(R_dist.begin(), min_element(R_dist.begin(), R_dist.end()));//pick the shortest side to compute the correlators with
                                    
                                    int mid_R;
                                    for(int q=0; q<3; q++)
                                    {
                                        if(q!=max_R && q!=min_R){mid_R==q;}
                                    }
                                    
                                    
                                    double x = R_dist[min_R]/R_dist[mid_R];
                                    double diff = (R_dist[max_R]-R_dist[mid_R]);
                                    double rat = (diff*diff)/ (R_dist[min_R]*R_dist[min_R]);
                                    double phi = asin(sqrt(1-rat));
                                    
                                    if(R_low<R_dist[max_R] && R_dist[max_R]<R_high)
                                        // if(R_low<R_dist[max_R] && R_dist[max_R]<R_high && 1.4<phi && phi<1.5 && 0.8<x && x<1.0)
                                    {
                                        R_sizes_gam->Fill(x,phi,jet_pt,eee_jsm*wt);
                                        
                                        // if(0.2 < aj && aj < 0.3){R_sizes_gam1->Fill(x,phi,jet_pt,eee_jsm*wt);}
                                        // if(0.5 < aj && aj < 0.6){R_sizes_gam2->Fill(x,phi,jet_pt,eee_jsm*wt);}
                                        // if(0.7 < aj && aj < 0.8){R_sizes_gam3->Fill(x,phi,jet_pt,eee_jsm*wt);}
                                        
                                    }
                                    
                                    // e3c_pt_hist_gam->Fill(R_dist[max_R],jet_pt,eee_jsm*wt);
                                    
                                    
                                    R_dist.clear();
                                    
                                    // }// jet-wake-wake corr
                                }//close m loop
                            }//close s loop
                        } //close if loop for 3 point corr
                    } //close j loop
                } //close i loop
            }//close if loop for gamma jets
            
            
            //Once jet processing is over, increase the condition value to match the next event id
            fjInputs.resize(0);//reshape the fjInputs
            sortedJets.resize(0);//reshape the sortedJets
            passedJets.resize(0);//reshape the passedJets
            t1->GetEntry(i);//getting entry again so we can increment the event id
            conditionValue = event_id;
        }
        a++;
        if(a%1000 == 0){cout<<a<<endl;}
    }
    
    cout<<"a = "<<a<<" and tree entries = " << t1->GetEntries()<<endl; //this makes sure I am going through all entries in the tree
    cout<<conditionValue<<endl; // this makes sure i hit all the "events"
    
    if(argu2.find(passedjets) != std::string::npos)
    {
        jet_pt_hist->Write();
        R_sizes->Write();
        if(argu3.find(three) != std::string::npos)
        {
            e3c_pt_hist->Write();
            // e3c_pt_hist1->Write();
            // e3c_pt_hist2->Write();
            
            // R_sizes_www->Write();
            R_sizes_jjw->Write();
            Rdistribution->Write();
            Rdistribution_wt->Write();
            // neg_map->Write();
            // pos_map->Write();
            // jet_map->Write();
        }
    }
    
    
    if(argu2.find(gamma) != std::string::npos)
    {
        gam_jet_pt_hist->Write();
        gam_pt_jet_pt->Write();
        if(argu3.find(three) != std::string::npos)
        {
            // e3c_pt_hist_gam->Write();
            // e3c_pt_hist_gam_jet->Write();
            R_sizes_gam->Write();
            // R_sizes_gam1->Write();
            // R_sizes_gam2->Write();
            // R_sizes_gam3->Write();
        }
    }
    
    
    delete fout_nom;
    f->Close();
    
}

