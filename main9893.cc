//in this script i compute hyb model correlators in new coordinates
//but in this case i am not fixing to the highest pT particle in the jet
#include <chrono>
#include <ctime>  
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
               
struct MyCoordinates {
double x;
double y;
};

struct XiPhiCoordinates {
    double xi;
    double phi;
};

struct ParticleTriplet {
PseudoJet baseParticle1;
PseudoJet baseParticle2;
PseudoJet thirdParticle;
int tripletType; // jjj = 0, jjw = 1, jww = 2, www = 3
MyCoordinates myCoordinates;
XiPhiCoordinates xiPhiCoordinates;
};  

// Gets the angle of rotation for constructing MyCoordinates for a particle-triplet.
double getTheta(double deltaY, double deltaPhi, double deltaR) {
    double theta = 0.0;
    if (deltaY > 0 && deltaPhi > 0) {
        theta = asin(deltaPhi/deltaR);
    } else if (deltaY < 0 && deltaPhi > 0) {
        theta = pi - asin(deltaPhi/deltaR);
    } else if (deltaY < 0 && deltaPhi < 0) {
        theta = pi - asin(deltaPhi/deltaR);
    } else if (deltaY > 0 && deltaPhi < 0) {
        theta = 2*pi + asin(deltaPhi/deltaR);
    }
    return theta;
};

MyCoordinates getMyCoordinates(ParticleTriplet particleTriplet) {
    // get distance to shift and angle to rotate between the highest pt and second particle
    // double deltaR = particleTriplet.baseParticle1.delta_R(particleTriplet.baseParticle2);
    double deltaR = delR(particleTriplet.baseParticle1, particleTriplet.baseParticle2);
    // double deltaY = particleTriplet.baseParticle2.rapidity() - particleTriplet.baseParticle1.rapidity();
    double deltaEta = particleTriplet.baseParticle2.eta() - particleTriplet.baseParticle1.eta();
    double deltaPhi = particleTriplet.baseParticle1.delta_phi_to(particleTriplet.baseParticle2);
    double theta = getTheta(deltaEta, deltaPhi, deltaR);

    // shift and rotate the third particle to to where it should lie in the new x-y plane
    // double y = particleTriplet.thirdParticle.rapidity() - particleTriplet.baseParticle1.rapidity();
    double eta_diff = particleTriplet.thirdParticle.eta() - particleTriplet.baseParticle1.eta();
    double phi = particleTriplet.baseParticle1.delta_phi_to(particleTriplet.thirdParticle);
    // double yRotated = y * cos(theta) + phi * sin(theta);
    double etaRotated = eta_diff * cos(theta) + phi * sin(theta);
    // double phiRotated = -1.0 * y * sin(theta) + phi * cos(theta);
    double phiRotated = -1.0 * eta_diff * sin(theta) + phi * cos(theta);

    // normalize coordinate choice so that 0 < y < sqrt(3)/2 and 0 < x < 1/2
    MyCoordinates coordinates;
    // coordinates.x = yRotated * 1.0/deltaR;
    coordinates.x = etaRotated * 1.0/deltaR;
    // coordinates.y = fabs(phiRotated) * 1.0/deltaR;
    coordinates.y = fabs(phiRotated) * 1.0/deltaR;
    return coordinates;
};

XiPhiCoordinates getXiPhiCoordinates(ParticleTriplet particleTriplet) {
    // get distance to shift and angle to rotate between the highest pt and second particle
    double RL = delR(particleTriplet.baseParticle1, particleTriplet.baseParticle2); 
    double RM = delR(particleTriplet.baseParticle1, particleTriplet.thirdParticle); 
    double RS = delR(particleTriplet.baseParticle2, particleTriplet.thirdParticle); 
    if (RS > RM) swap(RM, RS);
    if (RL < RM || RL < RS) cout << "error" << endl;

    XiPhiCoordinates coordinates;
    coordinates.xi = RS/RM;
    double difference = RL - RM;
    double ratio = difference*difference / (RS*RS);
    coordinates.phi = asin(sqrt(1 - ratio));
    return coordinates;
}

// Gets tripletType based on number of (positive) wake particles (user_index of 1) -- 0 = jjj, 1 = jjw, 2 = jww, 3 = www //modified from Arjun's where index of wake is 2
int getTripletType(PseudoJet p1, PseudoJet p2, PseudoJet p3) {
    int tripletType = 0;
    if (p1.user_info<MyInfo>().label_id() == 1) tripletType++;
    if (p2.user_info<MyInfo>().label_id() == 1) tripletType++;
    if (p3.user_info<MyInfo>().label_id() == 1) tripletType++;
    
    // if (p1.user_index() == 1) tripletType++;
    // if (p2.user_index() == 1) tripletType++;
    // if (p3.user_index() == 1) tripletType++;
   
    return tripletType;
};

// Orders the particles so that two base particles of a triangle are the two particles that define RL
// ParticleTriplet getTriplet(PseudoJet p1, PseudoJet p2, PseudoJet p3) {
//     double maxDistance = 0.0;
//     PseudoJet firstParticle; PseudoJet secondParticle; PseudoJet thirdParticle;

//     if (delR(p1,p2) > maxDistance) {
//         maxDistance = delR(p1,p2);
//         firstParticle = p1; secondParticle = p2; thirdParticle = p3;
//     }
//     if (delR(p1,p3) > maxDistance) {
//         maxDistance = delR(p1,p3);
//         firstParticle = p1; secondParticle = p3; thirdParticle = p2;
//     }
//     if (delR(p2,p3) > maxDistance) {
//         maxDistance = delR(p2,p3);
//         firstParticle = p2; secondParticle = p3; thirdParticle = p1;
//     }

//     double randomNumber = RANDOM_NUMBER_GENERATOR->Uniform(0, 1);
//     if (randomNumber >= 0.5) swap(firstParticle, secondParticle);

//     ParticleTriplet triplet;
//     triplet.baseParticle1 = firstParticle;
//     triplet.baseParticle2 = secondParticle;
//     triplet.thirdParticle = thirdParticle;
//     triplet.tripletType = getTripletType(p1, p2, p3);
//     triplet.myCoordinates = getMyCoordinates(triplet);

//     return triplet;
// }

struct TripletResult {
    ParticleTriplet triplet;
    double maxDistance;
    // int tripletType;
    bool isValid;
};

TripletResult getTriplet(PseudoJet p1, PseudoJet p2, PseudoJet p3, double Rlow, double Rhigh) {
    double maxDistance = 0.0;
    PseudoJet firstParticle, secondParticle, thirdParticle;

    // Determine the pair of particles with the maximum distance
    if (delR(p1, p2) > maxDistance) {
        maxDistance = delR(p1, p2);
        firstParticle = p1; secondParticle = p2; thirdParticle = p3;
    }
    if (delR(p1, p3) > maxDistance) {
        maxDistance = delR(p1, p3);
        firstParticle = p1; secondParticle = p3; thirdParticle = p2;
    }
    if (delR(p2, p3) > maxDistance) {
        maxDistance = delR(p2, p3);
        firstParticle = p2; secondParticle = p3; thirdParticle = p1;
    }
    
    if(maxDistance<Rlow || Rhigh<maxDistance)
    {
           return { ParticleTriplet(), maxDistance, false };
    }
    

    // Randomly swap firstParticle and secondParticle
    double randomNumber = ((double) rand() / (RAND_MAX));
    if (randomNumber >= 0.5) std::swap(firstParticle, secondParticle);

    // Fill the ParticleTriplet
    ParticleTriplet triplet;
    triplet.baseParticle1 = firstParticle;
    triplet.baseParticle2 = secondParticle;
    triplet.thirdParticle = thirdParticle;
    triplet.tripletType = getTripletType(p1, p2, p3);
    triplet.myCoordinates = getMyCoordinates(triplet);
    triplet.xiPhiCoordinates = getXiPhiCoordinates(triplet);

    // Return both the triplet and maxDistance
    // return { triplet, maxDistance, true };
     return { triplet, maxDistance, true };
};

double ComputeEEEC(ParticleTriplet particleTriplet, double jetpT, double n)
{
    
    double eee_jsm;
        if(n == 1.5)
        {
            eee_jsm = ((6*(pow(particleTriplet.baseParticle1.pt(),1.5))*(pow(particleTriplet.baseParticle2.pt(),1.5))*(pow(particleTriplet.thirdParticle.pt(),1.5)))/(pow(jetpT,4.5)));
        }
        else if(n == 0.5)
        {
            eee_jsm = ((6*(pow(particleTriplet.baseParticle1.pt(),0.5))*(pow(particleTriplet.baseParticle2.pt(),0.5))*(pow(particleTriplet.thirdParticle.pt(),0.5)))/(pow(jetpT,1.5)));
        }
        else
        {
            eee_jsm = ((6*particleTriplet.baseParticle1.pt()*particleTriplet.baseParticle2.pt()*particleTriplet.thirdParticle.pt())/(pow(jetpT,3)));
        }
 return eee_jsm;
};



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
  std::string argu10 = argv[10]; //tree_start_value
  std::string argu11 = argv[11]; //tree_stop_value
  std::string argu12 = argv[12]; //jet radius
 
 auto start = std::chrono::system_clock::now();
 
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
  TH3D *R_sizes_gam;
  TH3D *R_sizes_gam1;
  TH3D *R_sizes_gam2;
  TH3D *R_sizes_gam3;
  TH2D *eec_pt_hist;
  TH2D *gam_pt_jet_pt;
  TH2D *eec_pt_hist_gam;
  TH2D *e3c_pt_hist_gam;
  TH1 *jet_pt_hist;
  TH1 *gam_jet_pt_hist;
  TFile *fout_nom;
  TH2D *e3c_pt_hist_gam_jet;
  TH3D *neg_map;
  TH3D *pos_map;
  TH3D *jet_map;
  TH3D *eeec_xy;
  TH3D *eeec_jjj_xy;
  TH3D *eeec_jjw_xy;
  TH3D *eeec_jww_xy;
  TH3D *eeec_www_xy;
  TH3D *eeec;
  TH3D *eeec_jjj;
  TH3D *eeec_jjw;
  TH3D *eeec_jww;
  TH3D *eeec_www;
  

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
    
    // R_sizes = new TH3D("eeec_pt_hist", "EEEC and jet_pt 3D", 20, new_binsR, 20, new_binsRs, 48, new_bins_const);
   
    // R_sizes_gam = new TH3D("eeec_pt_hist_gam", "EEEC and jet_pt 3D gam", 20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    
    double y_max = sqrt(3.0)/2.0;
    eeec_xy = new TH3D("eeec_pt_hist","x,y coordinates of all",40,0,1,40,0,y_max,48,40,1000);
    eeec_jjj_xy = new TH3D("eeec_pt_hist_jjj","x,y coordinates of jjj",40,0,1,40,0,y_max,48,40,1000);
    eeec_jjw_xy = new TH3D("eeec_pt_hist_jjw","x,y coordinates of jjw",40,0,1,40,0,y_max,48,40,1000);
    eeec_jww_xy = new TH3D("eeec_pt_hist_jww","x,y coordinates of jww",40,0,1,40,0,y_max,48,40,1000);
    eeec_www_xy = new TH3D("eeec_pt_hist_www","x,y coordinates of jww",40,0,1,40,0,y_max,48,40,1000);
    
    eeec = new TH3D("eeec_pt_his_nomt","xi,phi coordinates of all",20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    eeec_jjj = new TH3D("eeec_pt_hist_jjj_nom","xi,phi coordinates of jjj",20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    eeec_jjw = new TH3D("eeec_pt_hist_jjw_nom","xi,phi coordinates of jjw",20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    eeec_jww = new TH3D("eeec_pt_hist_jww_nom","xi,phi coordinates of jww",20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    eeec_www = new TH3D("eeec_pt_hist_www_nom","xi,phi coordinates of jww",20, new_binsR, 20, new_binsRs, 48, new_bins_const);
    
   
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
    
    if(argu4.find(vac) != std::string::npos)
         {
          cout<< "Hybrid_Vac_" <<  argu1  <<  "_"<< argu2 << "_"<<argu5<<"_"<<argu6<<" "<<argu9<<endl;
         }
    else{
        cout<< "Hybrid_" <<  argu1  <<  "_"<< argu2 << "_"<<argu5<<"_"<<argu6<<" "<<argu9<<endl;
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
    // double R       = 0.8;    // Jet size.
    double pTMin   = 40;    // Min jet pT.
    double pTMax   = 1000;    // Max jet pT.
    float track_pt = 1.0 ; //Min track pT
    double R_low = std::stod(argu7);
    double R_high = std::stod(argu8);
    double power_wt = std::stod(argu9);
    int startvalue = std::stoi(argu10);
    int stopvalue = std::stod(argu11); //jet size
    double R = std::stod(argu12); //jet size
    //  double etaMax  = 0.9;    // Pseudorapidity range of detector.

// cout<<R_low<<"and "<<R_high<<endl;
    
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
    
  
    Vec4 particleTemp;
  
    cout<<"About to loop over trees"<<endl;

    std::vector <fastjet::PseudoJet> fjInputs;
    std::vector <double> weights;
    int n = 0;

// Set the condition so that you can run over a subset of events 
    Int_t conditionValue = startvalue;

    // // Index variable to be incremented so you can check if all entries are read
    Int_t a = 0;
    cout<<"Start value is "<<startvalue<<endl;
    cout<<"##############"<<endl;
    
    t1->GetEntry(startvalue);
    int event_start = event_id;
    int event_stop = event_start + 999999;
    // int event_stop = event_start + 100000;
    
    //For parsing the tree
    //   for (int i = 517064518; i < t1->GetEntries(); ++i) {
    //      t1->GetEntry(i);
    //     if(event_id==1800000){cout<<event_id<<" at "<<i<<endl;break; 
    //   }
    // cout<<"Stop when you get to event number: "<<event_stop<<endl;
    
    // cout<<"condition value is "<<conditionValue<<endl;
    // Loop over entries in the tree
    for (int i = startvalue; i < stopvalue; ++i) {
         t1->GetEntry(i);
        //  if(event_stop == event_id) break;
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
                    if(argu5.find(nowake) != std::string::npos && label != 1 && label != 2) //do not include the wake (positive or negative) particles
                        {
                        //  cout<<"no wake"<<endl;
                        double E = sqrt((m*m)+(px*px+py*py+pz*pz));
                        PseudoJet pTemp(px,py,pz,E);
                        float pT = pTemp.pt();
                        pTemp.set_user_info(new MyInfo(id,label,pT));
                        fjInputs.push_back(pTemp);
                        
                        }
                        if(argu5.find(yeswake) != std::string::npos && label !=2) //for pT subtracted files
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
            
             vector <fastjet::PseudoJet> constit, sorted_constit;
             vector<double> R_dist, constit_pt, R_dist_jjw, R_dist_www,R_dist_new;
             double eee_jsm,max_pt_part_id, deltaR_js;
             ParticleTriplet triplet;
             TripletResult result;
             int type = -1;
            
             if(argu2.find(passedjets) != std::string::npos)
             {
                for(int i=0; i< int (passedJets.size()); ++i)
                {
                    constit = passedJets[i].constituents();
                    double jet_pt = passedJets[i].pt();
                    jet_pt_hist->Fill(jet_pt,wt);
                    
                      for (int j=0; j<int(constit.size());j++)
                          {
                  
                            for (int s=j+1; s<int(constit.size()) ;s++)
                            {
                               
                                 for( int m=s+1; m<int(constit.size()); m++){
                                                
                                        
                                        result = getTriplet(constit[j],constit[s],constit[m], R_low, R_high);
                                        // cout<<"here"<<endl;
                                        
                                        if (result.isValid) {
                                            ParticleTriplet partriplet = result.triplet;
                                            double maxDistance = result.maxDistance;
                                            type = partriplet.tripletType;
                                            
                                            // cout<<constit[j].user_info<MyInfo>().label_id()<<" "<<constit[s].user_info<MyInfo>().label_id()<<" "<<constit[m].user_info<MyInfo>().label_id()<<endl;
                                            // cout<<type<<endl;
                                            
                                            double energy_weight = ComputeEEEC(partriplet, jet_pt, power_wt);
                                          
                                            if(argu5.find(yeswake) != std::string::npos)
                                            { 
                                                eeec_xy->Fill(partriplet.myCoordinates.x, partriplet.myCoordinates.y, jet_pt,energy_weight*wt);
                                                eeec->Fill(partriplet.xiPhiCoordinates.xi, partriplet.xiPhiCoordinates.phi,jet_pt,energy_weight*wt);
                                                if(type == 0)
                                                {
                                                eeec_jjj_xy->Fill(partriplet.myCoordinates.x, partriplet.myCoordinates.y, jet_pt,energy_weight*wt);
                                                 eeec_jjj->Fill(partriplet.xiPhiCoordinates.xi, partriplet.xiPhiCoordinates.phi,jet_pt,energy_weight*wt);
                                                }
                                                else if(type == 1)
                                                {
                                                 eeec_jjw_xy->Fill(partriplet.myCoordinates.x, partriplet.myCoordinates.y, jet_pt,energy_weight*wt); 
                                                 eeec_jjw->Fill(partriplet.xiPhiCoordinates.xi, partriplet.xiPhiCoordinates.phi,jet_pt,energy_weight*wt);
                                                }
                                                else if(type == 2)
                                                {
                                                eeec_jww_xy->Fill(partriplet.myCoordinates.x, partriplet.myCoordinates.y, jet_pt,energy_weight*wt);
                                                eeec_jww->Fill(partriplet.xiPhiCoordinates.xi, partriplet.xiPhiCoordinates.phi,jet_pt,energy_weight*wt);
                                                }
                                                else{
                                                 eeec_www_xy->Fill(partriplet.myCoordinates.x, partriplet.myCoordinates.y,jet_pt,energy_weight*wt); 
                                                 eeec_www->Fill(partriplet.xiPhiCoordinates.xi, partriplet.xiPhiCoordinates.phi,jet_pt,energy_weight*wt);
                                                }
                                            }
                                            else{
                                                 eeec_xy->Fill(partriplet.myCoordinates.x, partriplet.myCoordinates.y, jet_pt,energy_weight*wt); 
                                                 eeec->Fill(partriplet.xiPhiCoordinates.xi, partriplet.xiPhiCoordinates.phi,jet_pt,energy_weight*wt);
                                                
                                            }
                                        }
                                        else continue;
                                     }  
                                    }
                                 }
                            }
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
                    
                      for (int j=0; j<int(constit.size());j++)
                          {
                  
                            for (int s=j+1; s<int(constit.size()) ;s++)
                            {
                               
                                 for( int m=s+1; m<int(constit.size()); m++){
                                                
                                        
                                        result = getTriplet(constit[j],constit[s],constit[m], R_low, R_high);
                                        // cout<<"here"<<endl;
                                        
                                        if (result.isValid) {
                                            ParticleTriplet partriplet = result.triplet;
                                            double maxDistance = result.maxDistance;
                                            type = partriplet.tripletType;
                                            
                                            // cout<<constit[j].user_info<MyInfo>().label_id()<<" "<<constit[s].user_info<MyInfo>().label_id()<<" "<<constit[m].user_info<MyInfo>().label_id()<<endl;
                                            // cout<<type<<endl;
                                            
                                            double energy_weight = ComputeEEEC(partriplet, jet_pt, power_wt);
                                          
                                            if(argu5.find(yeswake) != std::string::npos)
                                            { 
                                                eeec_xy->Fill(partriplet.myCoordinates.x, partriplet.myCoordinates.y, jet_pt,energy_weight*wt);
                                                eeec->Fill(partriplet.xiPhiCoordinates.xi, partriplet.xiPhiCoordinates.phi,jet_pt,energy_weight*wt);
                                                if(type == 0)
                                                {
                                                eeec_jjj_xy->Fill(partriplet.myCoordinates.x, partriplet.myCoordinates.y, jet_pt,energy_weight*wt);
                                                 eeec_jjj->Fill(partriplet.xiPhiCoordinates.xi, partriplet.xiPhiCoordinates.phi,jet_pt,energy_weight*wt);
                                                }
                                                else if(type == 1)
                                                {
                                                 eeec_jjw_xy->Fill(partriplet.myCoordinates.x, partriplet.myCoordinates.y, jet_pt,energy_weight*wt); 
                                                 eeec_jjw->Fill(partriplet.xiPhiCoordinates.xi, partriplet.xiPhiCoordinates.phi,jet_pt,energy_weight*wt);
                                                }
                                                else if(type == 2)
                                                {
                                                eeec_jww_xy->Fill(partriplet.myCoordinates.x, partriplet.myCoordinates.y, jet_pt,energy_weight*wt);
                                                eeec_jww->Fill(partriplet.xiPhiCoordinates.xi, partriplet.xiPhiCoordinates.phi,jet_pt,energy_weight*wt);
                                                }
                                                else{
                                                 eeec_www_xy->Fill(partriplet.myCoordinates.x, partriplet.myCoordinates.y,jet_pt,energy_weight*wt); 
                                                 eeec_www->Fill(partriplet.xiPhiCoordinates.xi, partriplet.xiPhiCoordinates.phi,jet_pt,energy_weight*wt);
                                                }
                                            }
                                            else{
                                                eeec_xy->Fill(partriplet.myCoordinates.x, partriplet.myCoordinates.y, jet_pt,energy_weight*wt); 
                                                 eeec->Fill(partriplet.xiPhiCoordinates.xi, partriplet.xiPhiCoordinates.phi,jet_pt,energy_weight*wt);
                                                
                                            }
                                        }
                                        else continue;
                                     }  
                                    }
                                 }
                            }
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
    cout<<"Computation was from "<<event_start<<" to "<<event_stop<< " which matches the 999999"<<endl;
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "finished computation at " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s"<< std::endl;
    
    
    
    if(argu2.find(passedjets) != std::string::npos)
         { 
            jet_pt_hist->Write();
            eeec_xy->Write();
            eeec_jjj_xy->Write();
            eeec_jjw_xy->Write();
            eeec_jww_xy->Write();
            eeec_www_xy->Write();
            eeec->Write();
            eeec_jjj->Write();
            eeec_jjw->Write();
            eeec_jww->Write();
            eeec_www->Write();
         }
   
    
    if(argu2.find(gamma) != std::string::npos)
         {
             gam_jet_pt_hist->Write();
             gam_pt_jet_pt->Write();
             if(argu3.find(three) != std::string::npos)
             {
                eeec_xy->Write();
                eeec_jjj_xy->Write();
                eeec_jjw_xy->Write();
                eeec_jww_xy->Write();
                eeec_www_xy->Write();
                eeec->Write();
                eeec_jjj->Write();
                eeec_jjw->Write();
                eeec_jww->Write();
                eeec_www->Write();
             }
         }
        
    
    delete fout_nom;
    f->Close();

}

