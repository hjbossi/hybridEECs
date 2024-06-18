#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Recluster.hh"
#include "../../utils/NegaRecombiner.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace fastjet;
using namespace std;

// ROOT Classes
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TRandom3.h"

char* PARTICLES_FILE_NAME = "../../common-data/HYBRID_Hadrons_NoElastic_HighPt.out";
char* OUT_FILE_NAME = "../data-unbiased-140-160-gev-jets/E3C_PBPB_NO_WAKE.root";
char* OUT_INFO_FILE_NAME = "../data-unbiased-140-160-gev-jets/E3C_PbPb_No_Wake_Info.out";
int MAX_EVENTS = 10000000;
bool NO_ANGLES = true;

bool DO_NOT_INCLUDE_NEGATIVES = true;
bool DO_NOT_INCLUDE_ANY_WAKE_PARTICLES = true;

double JET_RADIUS = 0.8;
double MIN_PT = 140.0;
double MAX_PT = 160.0;
double MAX_ETA = 2.0;

TRandom3* RANDOM_NUMBER_GENERATOR = new TRandom3();

vector<double> RLS{
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
    0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8
};

int getRLIndex(double RL) {
    bool indexFound = false;
    for (int i = 0; i < RLS.size() - 1; i++) {
        if (RLS[i] <= RL && RL < RLS[i+1]) return i;
    }
    return RLS.size() - 2;
}

char* getHistName(int RLIndex, int tripletType, char* coordinates) {
    double RL_low = RLS[RLIndex];
    double RL_high = RLS[RLIndex + 1];

    char* histType;
    if (tripletType == 0) histType = "jjj";
    if (tripletType == 1) histType = "jjw";
    if (tripletType == 2) histType = "jww";
    if (tripletType == 3) histType = "www";
    if (tripletType == 4) histType = "all";

    char* histName;
    asprintf(&histName, "%s%s%s%s%.1f%s%.1f", histType, "_", coordinates, "_", RL_low, "_", RL_high);
    return histName;
}

double N_EVT = 0.0;
double AVG_CROSS = 0.0;
double N_JETS = 0.0;

const Int_t NUM_X_BINS = 20;
Double_t X_MIN = 0.0;
Double_t X_MAX = 1.0;

const Int_t NUM_Y_BINS = 20;
Double_t Y_MIN = 0.0;
Double_t Y_MAX = sqrt(3.0)/2.0;

const Int_t NUM_XI_BINS = 20;
Double_t XI_MIN = 0.0;
Double_t XI_MAX = 1.0;

const Int_t NUM_PHI_BINS = 20;
Double_t PHI_MIN = 0.0;
Double_t PHI_MAX = TMath::Pi()/2.0;

struct MyParticle {
    PseudoJet particle;
    int pdg_id;
    int label;
    double weight;
    double cross;
    double X;
    double Y;
    bool shouldDisregard;
};

struct JetData {
    PseudoJet jet;
    double weight; // this is the event weight
    vector<PseudoJet> constituents;
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

// Collects particles in the event
vector< MyParticle > collectParticles(ifstream& file) {
    // Clear all particles from past event
    vector<MyParticle> particles; particles.clear();
    string temp;

    // Get the weight and cross section of the event
    while (temp != "weight") {
        file >> temp;
    }
    file >> temp; double weight = stod(temp);
    file >> temp; file >> temp; double cross = stod(temp);
    file >> temp; file >> temp; double X = stod(temp);
    file >> temp; file >> temp; double Y = stod(temp);
    file >> temp; // Go to next line

    while (temp != "end") {
        // Get particle in the line
        double px = stod(temp);
        file >> temp; double py = stod(temp);
        file >> temp; double pz = stod(temp);
        file >> temp; double m = stod(temp);
        double E = sqrt(px*px + py*py + pz*pz + m*m);
        PseudoJet particle = PseudoJet(px, py, pz, E);

        // Get pdg_id and label of particle in the line
        file >> temp; int pdg_id = stoi(temp);
        file >> temp; int label = stoi(temp);

        // Add data to myParticle, and add myParticle to particles vector
        MyParticle myParticle;
        myParticle.particle = PseudoJet(px, py, pz, E);
        myParticle.pdg_id = pdg_id;
        myParticle.label = label;
        myParticle.weight = weight;
        myParticle.cross = cross;
        myParticle.X = X;
        myParticle.Y = Y;
        myParticle.shouldDisregard = false;

        if (!(myParticle.label == 2 && DO_NOT_INCLUDE_NEGATIVES) && !((label == 1 || label == 2) && DO_NOT_INCLUDE_ANY_WAKE_PARTICLES)) particles.push_back(myParticle);

        file >> temp; // Go to next line
    }

    // Check if file has ended
    file >> temp;
    if (temp != "#") file.close();

    return particles;
}

// Reconstruct jets from particles
vector<JetData> reconstructJets(vector<MyParticle> myParticles, double jetRadius) {
    double weight; double crossSection;
    vector<PseudoJet> particles; particles.clear();
    for(MyParticle myParticle : myParticles) {
        if (myParticle.label == -2) continue; // dont recluster particles in hard scattering (jet originators)
        PseudoJet particle = myParticle.particle;
        
        // set user_index to each particle so that we can distinguish between wake and nonwake particles even after jet reconstruction
        // label = 0: positive particle from jet fragmentation -> positive user_index
        // label = 1: positive particle in wake -> positive user_index
        // label = 2: negative particle in wake (depletion of plasma in direction opposite jet) -> negative user_index
        // label = 3: negative particle (hole) from elastic scattering -> negative user_index
        if (myParticle.label == 0) {
            particle.set_user_index(1);
        } else if (myParticle.label == 1) {
            particle.set_user_index(2);
        } else if (myParticle.label == 2 ) {
            particle.set_user_index(-1);
        } else if (myParticle.label == 3) {
            particle.set_user_index(-2);
        }

        if (myParticle.label != 0 && myParticle.label != 1 && myParticle.label != 2 && myParticle.label != 3) cout << "Error: Found particle with label " << myParticle.label << endl;
        particles.push_back(particle);

        weight = myParticle.weight;
        crossSection = myParticle.cross;
    }

    // choose a jet definition with antikt_algorithm and negative recombiner (negative recombiner only has an effect on negative particles)
    int ui = -123456;
    NegativeEnergyRecombiner uir(ui);
    JetDefinition jet_def(antikt_algorithm, jetRadius);
    jet_def.set_recombiner(&uir);

    // run the clustering, extract the jets
    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jetsInThisEvent = sorted_by_pt(cs.inclusive_jets());

    vector<JetData> jetsData; jetsData.clear();
    for (PseudoJet jet : jetsInThisEvent) {
        if (jet.pt() > MIN_PT && jet.pt() < MAX_PT && fabs(jet.eta()) < MAX_ETA) {
            JetData jetData;
            jetData.jet = jet;
            jetData.weight = weight;
            jetData.constituents = jet.constituents();
            jetsData.push_back(jetData);
            N_JETS += weight;
        }
    }

    // Add weight to N_evt
    N_EVT += weight;
    AVG_CROSS += crossSection;
    
    // cout << "Moving on to next event" << endl;
    return jetsData;
}

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
}

MyCoordinates getMyCoordinates(ParticleTriplet particleTriplet) {
    // get distance to shift and angle to rotate between the highest pt and second particle
    double deltaR = particleTriplet.baseParticle1.delta_R(particleTriplet.baseParticle2);
    double deltaY = particleTriplet.baseParticle2.rapidity() - particleTriplet.baseParticle1.rapidity();
    double deltaPhi = particleTriplet.baseParticle1.delta_phi_to(particleTriplet.baseParticle2);
    double theta = getTheta(deltaY, deltaPhi, deltaR);

    // shift and rotate the third particle to to where it should lie in the new x-y plane
    double y = particleTriplet.thirdParticle.rapidity() - particleTriplet.baseParticle1.rapidity();
    double phi = particleTriplet.baseParticle1.delta_phi_to(particleTriplet.thirdParticle);
    double yRotated = y * cos(theta) + phi * sin(theta);
    double phiRotated = -1.0 * y * sin(theta) + phi * cos(theta);

    // normalize coordinate choice so that 0 < y < sqrt(3)/2 and 0 < x < 1/2
    MyCoordinates coordinates;
    coordinates.x = yRotated * 1.0/deltaR;
    coordinates.y = fabs(phiRotated) * 1.0/deltaR;
    return coordinates;
}

XiPhiCoordinates getXiPhiCoordinates(ParticleTriplet particleTriplet) {
    // get distance to shift and angle to rotate between the highest pt and second particle
    double RL = particleTriplet.baseParticle1.delta_R(particleTriplet.baseParticle2);
    double RM = particleTriplet.baseParticle1.delta_R(particleTriplet.thirdParticle);
    double RS = particleTriplet.baseParticle2.delta_R(particleTriplet.thirdParticle);
    if (RS > RM) swap(RM, RS);
    if (RL < RM || RL < RS) cout << "error" << endl;

    XiPhiCoordinates coordinates;
    coordinates.xi = RS/RM;
    double difference = RL - RM;
    double ratio = difference*difference / (RS*RS);
    coordinates.phi = asin(sqrt(1 - ratio));
    return coordinates;
}

// Gets tripletType based on number of (positive) wake particles (user_index of 2) -- 0 = jjj, 1 = jjw, 2 = jww, 3 = www
int getTripletType(PseudoJet p1, PseudoJet p2, PseudoJet p3) {
    int tripletType = 0;
    if (p1.user_index() == 2) tripletType++;
    if (p2.user_index() == 2) tripletType++;
    if (p3.user_index() == 2) tripletType++;
    return tripletType;
}

// Orders the particles so that two base particles of a triangle are the two particles that define RL
ParticleTriplet getTriplet(PseudoJet p1, PseudoJet p2, PseudoJet p3) {
    double maxDistance = 0.0;
    PseudoJet firstParticle; PseudoJet secondParticle; PseudoJet thirdParticle;

    if (p1.delta_R(p2) > maxDistance) {
        maxDistance = p1.delta_R(p2);
        firstParticle = p1; secondParticle = p2; thirdParticle = p3;
    }
    if (p1.delta_R(p3) > maxDistance) {
        maxDistance = p1.delta_R(p3);
        firstParticle = p1; secondParticle = p3; thirdParticle = p2;
    }
    if (p2.delta_R(p3) > maxDistance) {
        maxDistance = p2.delta_R(p3);
        firstParticle = p2; secondParticle = p3; thirdParticle = p1;
    }

    double randomNumber = RANDOM_NUMBER_GENERATOR->Uniform(0, 1);
    if (randomNumber >= 0.5) swap(firstParticle, secondParticle);

    ParticleTriplet triplet;
    triplet.baseParticle1 = firstParticle;
    triplet.baseParticle2 = secondParticle;
    triplet.thirdParticle = thirdParticle;
    triplet.tripletType = getTripletType(p1, p2, p3);
    triplet.myCoordinates = getMyCoordinates(triplet);
    triplet.xiPhiCoordinates = getXiPhiCoordinates(triplet);

    return triplet;
}

// Constructs triplets of particles within the specified RL range
vector<ParticleTriplet> createParticleTriplets(vector<PseudoJet> particles) {
    vector<ParticleTriplet> triplets; triplets.clear();
    for (int i = 0; i < particles.size(); i++) {
        for (int j = i+1; j < particles.size(); j++) {
            for (int k = j+1; k < particles.size(); k++) {
                ParticleTriplet triplet = getTriplet(particles[i], particles[j], particles[k]);
                triplets.push_back(triplet);
            }
        }
    }
    return triplets;
}

int main(int argc, char **argv) {
    // Open the files to read from
    ifstream particlesFile; particlesFile.open(PARTICLES_FILE_NAME);

    // Create histograms
    vector< vector<TH2D*> > histsXY; histsXY.clear();
    for (int rlIndex = 0; rlIndex < RLS.size() - 1; rlIndex++) {
        vector<TH2D*> histsAtRL; histsAtRL.clear();
        for (int i = 0; i <= 4; i++) {
            char* histID; asprintf(&histID, "%s%d%d", "histXY", rlIndex, i);
            TH2D* hist = new TH2D(histID, "E3C in x, y coordinates", NUM_X_BINS, X_MIN, X_MAX, NUM_Y_BINS, Y_MIN, Y_MAX);
            histsAtRL.push_back(hist);
        }
        histsXY.push_back(histsAtRL);
    }
    vector< vector<TH2D*> > histsXiPhi; histsXiPhi.clear();
    for (int rlIndex = 0; rlIndex < RLS.size() - 1; rlIndex++) {
        vector<TH2D*> histsAtRL; histsAtRL.clear();
        for (int i = 0; i <= 4; i++) {
            char* histID; asprintf(&histID, "%s%d%d", "histXiPhi", rlIndex, i);
            TH2D* hist = new TH2D(histID, "E3C in xi, phi coordinates", NUM_XI_BINS, XI_MIN, XI_MAX, NUM_PHI_BINS, PHI_MIN, PHI_MAX);
            histsAtRL.push_back(hist);
        }
        histsXiPhi.push_back(histsAtRL);
    }

    // Construct the jets using anti-kt per event and add them to the file
    int eventNumber = 0;
    while (particlesFile.is_open() && eventNumber < MAX_EVENTS) {
        vector<MyParticle> particles = collectParticles(particlesFile);
        vector<JetData> jets = reconstructJets(particles, JET_RADIUS);

        for (JetData jetData : jets) {
            vector<ParticleTriplet> triplets = createParticleTriplets(jetData.constituents);
            for (ParticleTriplet triplet : triplets) {
                double histWeight;
                double pt_product = triplet.baseParticle1.pt() * triplet.baseParticle2.pt() * triplet.thirdParticle.pt();
                double distance_product = triplet.baseParticle1.delta_R(triplet.baseParticle2) * triplet.baseParticle1.delta_R(triplet.thirdParticle) * triplet.baseParticle2.delta_R(triplet.thirdParticle);
                if (NO_ANGLES) histWeight = jetData.weight * pt_product / (pow(jetData.jet.pt(), 3));
                else histWeight = jetData.weight * pt_product * distance_product / (pow(jetData.jet.pt(), 3));

                int RLIndex = getRLIndex(triplet.baseParticle1.delta_R(triplet.baseParticle2));
                histsXY[RLIndex][triplet.tripletType]->Fill(triplet.myCoordinates.x, triplet.myCoordinates.y, histWeight);
                histsXY[RLIndex][4]->Fill(triplet.myCoordinates.x, triplet.myCoordinates.y, histWeight);
                histsXiPhi[RLIndex][triplet.tripletType]->Fill(triplet.xiPhiCoordinates.xi, triplet.xiPhiCoordinates.phi, histWeight);
                histsXiPhi[RLIndex][4]->Fill(triplet.xiPhiCoordinates.xi, triplet.xiPhiCoordinates.phi, histWeight);
            }
        }

        cout << "Got observables in event " << eventNumber << endl;
        eventNumber++;
    } 
    int numEvents = eventNumber;

    TFile* ROOTFile = TFile::Open(OUT_FILE_NAME, "RECREATE");
    for (int RLIndex = 0; RLIndex < histsXY.size(); RLIndex++) {
        for (int tripletType = 0; tripletType <= 4; tripletType++) {
            char* histName = getHistName(RLIndex, tripletType, "xy");
            histsXY[RLIndex][tripletType]->SetName(histName);
            histsXY[RLIndex][tripletType]->Scale(1.0/N_JETS);
            histsXY[RLIndex][tripletType]->Write();
        }
    }
    for (int RLIndex = 0; RLIndex < histsXY.size(); RLIndex++) {
        for (int tripletType = 0; tripletType <= 4; tripletType++) {
            char* histName = getHistName(RLIndex, tripletType, "xiphi");
            histsXiPhi[RLIndex][tripletType]->SetName(histName);
            histsXiPhi[RLIndex][tripletType]->Scale(1.0/N_JETS);
            histsXiPhi[RLIndex][tripletType]->Write();
        }
    }
    ROOTFile->Write();
    ROOTFile->ls();
    ROOTFile->Close();

    ofstream outInfoFile; outInfoFile.open(OUT_INFO_FILE_NAME);
    outInfoFile << "particlesFile" << " " << PARTICLES_FILE_NAME << endl;
    outInfoFile << "numEvents" << " " << numEvents << endl;
    outInfoFile << "jetRadius" << " " << JET_RADIUS << endl;
    outInfoFile << "minJetpt" << " " << MIN_PT << endl;
    outInfoFile << "maxJetpt" << " " << MAX_PT << endl;
    outInfoFile << "maxJetEta" << " " << MAX_ETA << endl;
    outInfoFile << "N_evt" << " " << N_EVT << endl;
    outInfoFile << "avg_cross" << " " << AVG_CROSS/numEvents << endl;
    outInfoFile << "N_jets" << " " << N_JETS << endl;

    // Close the files
    if (particlesFile.is_open()) particlesFile.close();
    if (outInfoFile.is_open()) outInfoFile.close();
}
