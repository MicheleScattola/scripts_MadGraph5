#include "TH2D.h"
#include "ROOT/RDataFrame.hxx"
#include <iostream>

using namespace std;
 
void Zmumu(){

    string file = "/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/myExercises/Zmumu/Events/ATLAS/tag_1_delphes_events.root"
    // retrieve data and filter by RP_costheta
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(treeName, infile_data);

    if (!df.HasColumn(dataColName)) {
    std::cerr << "[Fitter] Data column missing: " << dataColName << std::endl;
    return result;
    }

}

