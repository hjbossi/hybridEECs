
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

int main(int argc, char* argv[]) {
    
//Its useful to have options for what I am trying to run
  std::string argu1 = argv[1]; //parton or hadron level
  std::string argu2 = argv[2]; //gamma jets or no
  std::string argu3 = argv[3]; //3 point corr or no
  std::string argu4 = argv[4]; //name of the input file 
  std::string argu5 = argv[5]; //wake or no wake
  char* argu6 = argv[6]; //adds date to file
 
 
  
//   int start = std::stoi(argu5);
//   int end = std::stoi(argu6);
  
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

    Double_t from = -4;
    Double_t to = 0;
    Int_t bins = 100;
    Double_t width = (to-from)/bins;
    Double_t new_bins[101] = {};
    for (Int_t i = 0; i <= bins; i++)
    {
        new_bins[i] = TMath::Power(10.0, from + i * width);
    }

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

    // const int numLabels = 4;
    // double step = (maxAngleRad - minAngleRad) / nDivisions;
    // for (int i = 0; i <= nDivisions; i++) {
    //     double angleRad = minAngleRad + i * step;
    //     double anglePi = angleRad / TMath::Pi();
    //     std::stringstream label;
    //     label << anglePi << "π";
    //     angleAxis->SetBinLabel(i + 1, label.str().c_str());
    // }

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
    
    R_sizes = new TH3D("eeec_pt_hist", "EEEC and jet_pt 3D gam", 20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    
    
    if(argu3.find(three) != std::string::npos)
    {
        e3c_hist = new TH1D("e3c_hist","E3C", 100, new_bins);

        e3c_pt_hist = new TH2D("e3c_pt_hist", "E3C and jet_pt 2D", 100, new_bins, 48, 40, 1000);
        
        e3c_pt_hist1 = new TH2D("e3c_pt_hist1", "E3C and jet_pt 2D", 100, new_bins, 48, 40, 1000);
        
        e3c_pt_hist2 = new TH2D("e3c_pt_hist2", "E3C and jet_pt 2D", 100, new_bins, 48, 40, 1000);
    
        e3c_pt_hist_gam = new TH2D("e3c_pt_hist_gam", "E3C and jet_pt 2D", 100, new_bins, 48, 40, 1000);
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
   Float_t px, py, pz, m, wt;
   Int_t event_id, id, label;
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
    float track_pt = 1 ; //Min track pT
    double R_low = 0.6;
    double R_high = 0.7;
    //  double etaMax  = 0.9;    // Pseudorapidity range of detector.

    
    // Defining a class to help store input particle PIDs to access later to tag heavy quark hadrons
    class MyInfo: public fastjet::PseudoJet::UserInfoBase //class name my info, inherits from Pseudojet::UserInfoBase
    {
    public:
        MyInfo(int id): _pdg_id(id){} //MyInfo(int id) is a parametrized class constructor. Has to have the same name as the class but no return type
        int pdg_id() const {return _pdg_id;}//here we are defining a member function pdg_id() of the class that returns an integer value
    protected:
        int _pdg_id;
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
////Variables
    // double weight,crossx, x, y, mass, label, index, mom_x, mom_y,mom_z;
    // int eventID, pdg_id;
    Vec4 particleTemp;
    // vector<double> event_id, x_pos, y_pos, px, py, pz,m, unique_event_id, pdg_ID ;
     vector<double> unique_event_id, ev_wt;
    
 
////Remove repeated event ids so can loop over events
    Int_t nentries = (Int_t)t1->GetEntries();
    set<int> s;
    for(unsigned i = 0; i < nentries; ++i){t1->GetEntry(i); s.insert(event_id);}
    unique_event_id.assign(s.begin(), s.end());
    
    cout <<unique_event_id.size()<<" after set"<<endl;

    std::vector <fastjet::PseudoJet> fjInputs;
    std::vector <double> weights;
    int a = 0;
    int n = 0;

//This for loop loops over events
     for(int k = 0; k<unique_event_id.size(); k++)
             {
                 for (Int_t i=n; i<nentries; i++) 
                 {
                     t1->GetEntry(i);
                    //  cout<<"n "<<n<<" i "<<i<<endl; 
                     if(a == event_id)
                     {   
                        if(label != -2) //This excludes the hard parton that scatters from the jet 
                         {
                            if(argu4.find(vac) != std::string::npos) //in case of vac include everything 
                                {
                                //  cout<<"a "<<a<<" label "<<label<<endl;
                                 double E = sqrt((m*m)+(px*px+py*py+pz*pz));
                                 PseudoJet pTemp(px,py,pz,E);
                                 pTemp.set_user_info(new MyInfo(id));
                                 fjInputs.push_back(pTemp);
                                 weights.push_back(wt);
                                }
                            else //in case of hybrid data 
                            {
                            if(argu5.find(nowake) != std::string::npos && label != 1 && label != 2) //do not include the wake particles 
                                {
                                //  cout<<"no wake"<<endl;
                                 double E = sqrt((m*m)+(px*px+py*py+pz*pz));
                                 PseudoJet pTemp(px,py,pz,E);
                                 pTemp.set_user_info(new MyInfo(id));
                                 fjInputs.push_back(pTemp);
                                 weights.push_back(wt);
                                }
                            else //include all particles from jets + wake 
                                 {
                                    if(argu5.find(yeswake) != std::string::npos)
                                    {
                                    //  cout<<"yes wake"<< " " <<a<<endl;
                                     double E = sqrt((m*m)+(px*px+py*py+pz*pz));
                                     PseudoJet pTemp(px,py,pz,E);
                                     pTemp.set_user_info(new MyInfo(id));
                                     fjInputs.push_back(pTemp);
                                     weights.push_back(wt);
                                    }
                                 }
                            }
                         }
                       n = n+1;
                     }
                     else break;
                   }
              
            weights.erase(unique(weights.begin(), weights.end()), weights.end()); //removing repeated weights from the vector 
            // cout<<"weight size "<<weights.size()<<endl;
            // cout<<"here"<<endl;
            double event_wt = weights[k];
            // cout<<k<<" associated wt "<<event_wt<<endl;
            
             ////Jet definition
             vector <fastjet::PseudoJet> inclusiveJets, sortedJets,passedJets,constit_i,constit_j;
             vector <fastjet::PseudoJet> gammaJets, EtPhoton, PhotonParticle;
             fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, R, power);
             fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
             inclusiveJets = clustSeq.inclusive_jets(pTMin);
             sortedJets    = sorted_by_pt(inclusiveJets);
             
              //cout statements to check jets
            //  if(sortedJets.size()==0){cout<<"no jets"<<endl;}
            //  else{cout<<a<<" "<<sortedJets.size()<<endl;}
            //  for(int j = 0; j<sortedJets.size(); j++)
            //         {
            //             cout<< a <<" "<<sortedJets[j].pt()<<endl;
            //         }
            //  //
            
             ////Jet selections
             double rap_cut = 2;//Jet axis rapidity cut
             fastjet::Selector select_pt_min = fastjet::SelectorPtMin(pTMin);
             fastjet::Selector select_pt_max = fastjet::SelectorPtMax(pTMax);
             fastjet::Selector select_rap = fastjet::SelectorAbsRapMax(rap_cut);//For setting min jet pt and rapidity cut (rapidity cut should be (0.9-jetR)
             fastjet::Selector select_all = select_pt_min && select_rap && select_pt_max;//order of operation doesn't matter in how its defined here
             
             passedJets = select_all(sortedJets);
             
             //This condition should hold true for any gamma-jet pair because fastjet should see the recoil photon as a jet
             if(argu2.find(gamma) != std::string::npos)
             {
             if(passedJets.size()>=2)
             {
            
    //         ///Finding the direct photon
    //          ///The cone of the direct photon is 0.4. It can only have a certain amount of transverse energy in this cone to be considered isolated. 
    //          //It also has certain phase space cuts to be considered isolated (delta_phi)
             double transE=0;
             double delPhi_iso = 0;
             
                 //Loop over all jets in an event
                 for(int i=0; i< passedJets.size(); i++)
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
                                //  cout<<delR_iso<<" "<<delPhi_iso<<endl;
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
                                //  cout<<"et "<<transE<<setprecision(7)<<endl;
                                 //If the Et_iso condition is satisfied, I have found a recoil jet to the direct photon. Add this to the list of jets
                                 if(transE <= 5.)
                                 {
                                     gam_pt_jet_pt->Fill(fjInputs[j].pt(),passedJets[i].pt());
                                     double gam_pt=fjInputs[j].pt();
                                     cout<<"value "<<gam_pt<<endl;
                                     passedJets[i].set_user_info(new MyInfoPt(gam_pt));
                                     gammaJets.push_back(passedJets[i]);
                                     
                                 }
                                 
                                 transE = 0;
                    
                            }
                         }
                         
                     }
                 }
            }
             if(gammaJets.size() != 0){cout<<gammaJets.size()<<" and passjets "<< passedJets.size()<<endl;}
             }
            
             //
             //////ENC CODE BEGINS----------------------------------------------------------------------------------------------------------------------------------------------------------
             /////This has to be multiplying pairwise WITHIN a single jet so we first need to loop over the vector of Jets (done above)
             //////Initializing some variables to use
             vector<double> energy_pairs_e2c, jetE_e2c, Rdistvec_e2c, logR_distvec_e2c;
             vector<double> energy_pairs, Rdistvec, jetE, logR_distvec, jet_spectra, max_logR_distvec, max_R_distvec, min_R_distvec, mid_R_distvec, R_dist, logR_dist, const_size, num_pairs;
             vector <fastjet::PseudoJet> constit, constit_2;
    
              ////For all passed jets
    if(argu2.find(passedjets) != std::string::npos)
    {
        ////Loop for jets with more than 2 constituents
        for(int i=0; i< int (passedJets.size()); ++i)
        {
            constit = passedJets[i].constituents();
            if(int(constit.size()) == 2) continue;//get out of the loop if the jet has only 2 constituents
            double jet_pt = passedJets[i].pt();
            jet_pt_hist->Fill(jet_pt,event_wt);
            
            for(int j=0; j < int(constit.size());j++)
            {
                if(argu3.find(three) != std::string::npos)
                {
                    for (int s=0; s<int(constit.size());s++)
                    {
                        for( int m=0; m<int(constit.size()); m++)
                        {

                            double eee_jsm = ((constit[j].pt()*constit[s].pt()*constit[m].pt())/(pow(jet_pt,3)));
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
                                 if(q!=max_R && q!=min_R){mid_R=q;}
                                 break;
                                }
                                     
                                     
                             double x = R_dist[min_R]/R_dist[mid_R];
                             double diff = (R_dist[max_R]-R_dist[mid_R]);
                             double rat = (diff*diff)/ (R_dist[min_R]*R_dist[min_R]);
                             double phi = asin(sqrt(1-rat));
                             
                             if(0.2<R_dist[max_R] && R_dist[max_R]<0.3)
                             {
                                R_sizes->Fill(x,phi,jet_pt,eee_jsm*event_wt);
                             }
                             
                             
                            if((1.4<phi && phi<1.5) && (0.9<x && x<1.0))
                            {
                             e3c_pt_hist->Fill(R_dist[max_R],jet_pt,eee_jsm*event_wt);
                            }
                            
                            if((1.4<phi && phi<1.5) && (0.05<x && x<0.1))
                            {
                            e3c_pt_hist1->Fill(R_dist[max_R],jet_pt,eee_jsm*event_wt);
                            }
                            
                            if((0.7<phi && phi<0.8) && (0.9<x && x<1.0) )
                            {
                            e3c_pt_hist1->Fill(R_dist[max_R],jet_pt,eee_jsm*event_wt);
                            }
                             R_dist.clear();
                            
                        }//close m loop
                    }//close s loop
                } //close if loop for 3 point corr
            } //close j loop
        } //close i loop
    }//close if loop for passed jets
             
             
    // // // //         ////FOR GAMMA JETS -------------------------------------------------
             if(argu2.find(gamma) != std::string::npos)
             {

    //              ////Loop for jets with more than 2 constituents
                 for(int i=0; i< int (gammaJets.size()); ++i)
                 {
                     constit = gammaJets[i].constituents();
                     
                    //  double jet_pt_gam = gammaJets[i].pt();
                     double jet_pt_gam =  gammaJets[i].user_info<MyInfoPt>().gam_pt();
                     gam_jet_pt_hist->Fill(jet_pt_gam,event_wt);
                     
    
                    if(argu3.find(three) != std::string::npos)
                    {
                        for(int j=0; j < int(constit.size());j++)
                        {
                             for (int s=0; s<int(constit.size());s++)
                             {
                                 for( int m=0; m<int(constit.size()); m++)
    
                                 {
                                    
                                     double eee_jsm = ((constit[j].pt()*constit[s].pt()*constit[m].pt())/(pow(jet_pt_gam,3)));
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
                                        if(q!=max_R && q!=min_R){mid_R=q;}
                                        break;
                                    }
                                     
                                     
                                     double x = R_dist[min_R]/R_dist[mid_R];
                                     double diff = (R_dist[max_R]-R_dist[mid_R]);
                                     double rat = (diff*diff)/ (R_dist[min_R]*R_dist[min_R]);
                                     double phi = asin(sqrt(1-rat));
                                     
                                    
                                     e3c_pt_hist_gam->Fill(R_dist[max_R],jet_pt_gam,eee_jsm*event_wt);
                                     
                                     if(R_low<R_dist[max_R] && R_dist[max_R]<R_high)
                                     {
                                     R_sizes->Fill(x,phi,jet_pt_gam,eee_jsm*event_wt);
                                     }
                                     R_dist.clear();
    
                                 }//close m loop
                             } //close s loop for the 2 point
                         }//close if loop for 3 point
                    //  }//close min track pT loop
                     }//close j loop
                 }//close jet loop
             }//close if loop for gamma jets
   
    // // ///PhotonNoJet----------------------------------------------------------------
    // //          vector<double>  particle_pt;
    
    // //          if(argu2.find(gammaNoJet) != std::string::npos)
    // //          {
                 
    // //              //Loop over all the particles in an event
    // //              for(int j=0; j< int (fjInputs.size()); ++j)
    // //              {
    // //                  //Find a photon
    // //                  if(fjInputs[j].user_info<MyInfo>().pdg_id() == 22 && fjInputs[j].pt()>100  && fjInputs[j].rap()<2.)
    // //                  {
    // //                      cout<<"event "<<k << " found"<<endl;
    // //                      //Find particles that are greater than pi/2 away from the photon.
    // //                      for(int i=0; i< int (fjInputs.size()); ++i)
    // //                      {
    // //                          double delPhi = fjInputs[j].delta_phi_to(fjInputs[i]);
    // //                          if (delPhi > (1/2)*pi )
    // //                          {
    // //                              PhotonParticle.push_back(fjInputs[i]);
    // //                              particle_pt.push_back(fjInputs[i].pt());
    // //                              cout<<fjInputs[i].pt()<<endl;;
                                 
    // //                          }
    // //                      }
    // //                      double pt = std::accumulate(particle_pt.begin(), particle_pt.end(), 0.0f);
                         
    // //                      for(int j=0; j< int (PhotonParticle.size()); ++j)
    // //                      {
                             
    // //                          for (int s=0; s< j;s++)
    // //                          {
    // //                              //                          cout<<"S "<<s<< "J "<<j<<endl;
    // //                              double ee_js = ((2*PhotonParticle[j].pt()*PhotonParticle[s].pt())/(pow(pt,2)));
    // //                              double deltaR_e2c_js = PhotonParticle[j].delta_R(PhotonParticle[s]);
                                 
    // //                              //Appending the vectors with this information
    // //                              energy_pairs_e2c.push_back(ee_js);
    // //                              Rdistvec_e2c.push_back(deltaR_e2c_js);
                                 
    // //                          } //close s loop for the 2 point
    // //                      }//close j loop for the 2 point
    // //                      PhotonParticle.clear();
    // //                      particle_pt.clear();
    // //                  }//close if loop
    // //              }//close particle loop
    // //          }//close no jet loop
    
             

           
             fjInputs.resize(0);//reshape the fjInputs
             a++; //increment the value of a to move on to next event
             if(a%1000 == 0){cout<<a<<endl;}
             if(a>unique_event_id.size()) break;
            }
        
   
    if(argu2.find(passedjets) != std::string::npos)
         { 
            jet_pt_hist->Write();
            // eec_pt_hist->Write();
            if(argu3.find(three) != std::string::npos)
             {
                e3c_pt_hist->Write();
                e3c_pt_hist1->Write();
                e3c_pt_hist2->Write();
                R_sizes->Write();
             }
         }
   
    
    if(argu2.find(gamma) != std::string::npos)
         {
             gam_jet_pt_hist->Write();
            //  eec_pt_hist_gam->Write();
             gam_pt_jet_pt->Write();
             if(argu3.find(three) != std::string::npos)
             {
                // e3c_pt_hist_gam->Write();
                R_sizes->Write();
             }
         }
    cout<<"success"<<endl;
    delete fout_nom;
    
    f->Close();

}

