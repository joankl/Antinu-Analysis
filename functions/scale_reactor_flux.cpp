/***
Reactor scaling tool
Author: Charlie Mills <Charlie.Mills@sussex.ac.uk>
Revision History: 08/01/20: Tool created 
                  21/02/20: Oscillation added (using Stefan Alexandru Nae's rat code)
                  06/10/20: Added flag for rat versions between 6.18.4 and 6.18.7
                  05/10/21: Added matter effects option and other small tweaks - James Page <James.Page@sussex.ac.uk>
Usage:  example command to compile: g++ -g -std=c++1y -o new_oscillation_scale_reactor_flux.exe new_oscillation_scale_reactor_flux.cpp `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux g++ -g -std=c++1y -o scale_reactor_flux.exe scale_reactor_flux.cpp `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux -I${RATROOT}/include/external
        run with ./scale_reactor_flux.exe FILE OFFLINE OSCILLATION REACTOR_INPUT Ne
            where FILE is an ntuple, OFFLINE is a bool to use the tables in the rat/data table - this is faster but less accurate,
            OSCILLATION is a bool to apply oscillation to the antineutrinos, REACTOR_INPUT is a bool on whether to alter the reactor names,
            and Ne is the electron density which can be: left blank (defaulting to zero), given a float, or given the string "crust" for the
            average crust electron density.
***/

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TRandom3.h>
#include <RAT/DB.hh>
#include <TMath.h>
#include <TVector3.h>

//define functions

TVector3 LLAtoECEF(double longitude, double latitude, double altitude) {
  // reference http://www.mathworks.co.uk/help/aeroblks/llatoecefposition.html
  static double toRad = TMath::Pi()/180.;
  static double Earthradius = 6378137.0; //Radius of the Earth (in meters)
  static double f = 1./298.257223563; //Flattening factor WGS84 Model
  static double L, rs, x, y, z;
  L = atan( pow((1. - f),2)*tan(latitude*toRad))*180./TMath::Pi();
  rs = sqrt( pow(Earthradius,2)/(1. + (1./pow((1. - f),2) - 1.)*pow(sin(L*toRad),2)));
  x = (rs*cos(L*toRad)*cos(longitude*toRad) + altitude*cos(latitude*toRad)*cos(longitude*toRad))/1000; // in km
  y = (rs*cos(L*toRad)*sin(longitude*toRad) + altitude*cos(latitude*toRad)*sin(longitude*toRad))/1000; // in km
  z = (rs*sin(L*toRad) + altitude*sin(latitude*toRad))/1000; // in km

  TVector3 ECEF = TVector3(x,y,z);

  return ECEF;
}

double GetReactorDistanceLLA(const double &longitude, const double&latitude, const double &altitude) {
    const TVector3 SNO_ECEF_coord_ = TVector3(672.87,-4347.18,4600.51);
    double dist = (LLAtoECEF(longitude, latitude,altitude) - SNO_ECEF_coord_).Mag();
  return dist;
}

std::vector<std::string> SplitString(std::string str, bool Osc){
    std::istringstream buf(str);
    std::istream_iterator<std::string> beg(buf), end;
    std::vector<std::string> tokens(beg, end); //each word of string now in vector
    std::vector<std::string> info;
    std::string dummy = "";

    for(int i=0;i<tokens.size();i++){ //combine back to reactor name, core number
        if(i==0){
            dummy = tokens.at(i);
            if(i!=tokens.size()-2){
                dummy += " ";
            }
        }
        else if(i!=tokens.size()-1){
            dummy += tokens.at(i);
            if(i!=tokens.size()-2){
                dummy += " ";
            }
        }
        else{
            info.push_back(dummy);
            info.push_back(tokens.at(i));
        }
    }

    return info;
}

int ScaleReactorFlux(std::string inputFile, bool offline_mode, bool oscillations, bool rat_version_corr, Double_t Ne, std::string ratDBTag){

    if(offline_mode){
        std::cout << "WARNING: Working in offline mode, this is less accurate as default tables will be used." << std::endl;
    }

    if(ratDBTag == "liveDB"){
        std::cout << "WARNING: Using the live version of ratdb, these tables may change! Use a ratdb tag for a consistent set of tables." << std::endl;
    }

    TFile *reactorFile;

    //load in ntuple
    try{
        reactorFile = TFile::Open(inputFile.c_str());
    }
    catch(...){
        std::cout << "Could not open input file." << std::endl;
        return 1;
    }

    std::string output_string = "scaled_" + inputFile;
    std::string output_hist_string = "scaled_histogram_" + inputFile;
    TFile *scaled_hist = new TFile(output_hist_string.c_str(),"RECREATE");
    TFile *scaled_ntuple = new TFile(output_string.c_str(),"RECREATE");

    //clone the ntuple, we want to keep the original ntuple untouched and produce a new file after scaling
    TTree *reactorEventInfo = (TTree *) reactorFile->Get("output"); 
    TTree *scaledReactorEventInfo = reactorEventInfo->CloneTree(0);

    //Set branch addresses for relevant parameters
    Int_t evIndex;
    Int_t mcIndex;
    Double_t parentKE1; // antineutrino energy
    Double_t reconEnergy; //recon energy
    Int_t runNum; //run number associated with MC
    Int_t lastRunNum = -999; //last run number loaded in DB
    TString *originReactor = NULL; // string for reactor where antineutrino was produced
    Bool_t *valid;

    reactorEventInfo->SetBranchAddress("mcIndex",&mcIndex);
    reactorEventInfo->SetBranchAddress("evIndex",&evIndex);
    reactorEventInfo->SetBranchAddress("parentMeta1",&originReactor);
    reactorEventInfo->SetBranchAddress("parentKE1",&parentKE1);
    reactorEventInfo->SetBranchAddress("energy",&reconEnergy);
    reactorEventInfo->SetBranchAddress("runID",&runNum);
    reactorEventInfo->SetBranchAddress("fitValid",&valid);

    //loop over events and decide whether to keep the event or not, initialise variables and histograms first

    TH1 *spectraRaw = new TH1D("raw", "Antinu Unscaled Spectra", 100, 0, 10);
    TH1 *spectraScaled = new TH1D("scaled", "Antinu Scaled Spectra", 100, 0, 10);
    TH1 *spectraRawRecon = new TH1D("raw_recon", "Antinu Unscaled Spectra", 100, 0, 10);
    TH1 *spectraScaledRecon = new TH1D("scaled_recon", "Antinu Scaled Spectra", 100, 0, 10);

    Int_t nentries = reactorEventInfo->GetEntries();
    TRandom3 *rndm_dummy = new TRandom3();
    Double_t rndm;
    Double_t efficiency;
    std::vector<std::string> reacts_eff_names;
    std::vector<Double_t> reacts_eff_vals;
    int num_excess_eff_evts = 0;
    int num_no_table_evts = 0;
    RAT::DBLinkPtr linkdb;
    RAT::DB *db = RAT::DB::Get();
    RAT::DS::Run run;

    if(offline_mode){
        db->SetAirplaneModeStatus(true);
        db->LoadDefaults();
    }
    else{
        try{
            db->LoadDefaults();
            db->SetServer("postgres://snoplus@pgsql.snopl.us:5400/ratdb");
            run.SetRunID(runNum);
            db->BeginOfRun(run);
            if(ratDBTag != "liveDB") db->SetDbTag(ratDBTag);
            std::cout << "RAT DB tag: " << db->GetDBTag() << std::endl;
        }
        catch(...){
            std::cout << "Invalid ratDBTag!" << std::endl;
            return 1;            
        }
    }

    Int_t previousmcIndex = -999;
    //go through entries of ttree
    for(Int_t a=0;a<nentries;a++){
        reactorEventInfo->GetEntry(a);

	// Save all the retriggers and delayed event of the saved prompt event
	if( mcIndex == previousmcIndex ){ scaledReactorEventInfo->Fill(); }

	// Only apply scaling to the first trigger of the mcIndex, which is the most likely to be the prompt
	// Skip all events that are not propmt
	if( evIndex != 0 ){ continue; }
	
        const char *originReactorString = originReactor->Data();
        std::vector<std::string> originReactorVect = SplitString(originReactorString, 0);
        spectraRaw->Fill(parentKE1);

        if(valid){
            spectraRawRecon->Fill(reconEnergy);
        }
        rndm = rndm_dummy->Rndm();

        if(!offline_mode && lastRunNum != runNum){
            run.SetRunID(runNum);
            db->BeginOfRun(run);
            lastRunNum = runNum;
        }
        if(originReactorVect[0] == "LENINGRAD-2" and rat_version_corr){ //hacks for rat versions 6.18.4 < 6.18.7
            originReactorVect[0] = "LENINGRAD 2";
        }
        if(originReactorVect[0] == "SHIN-KORI" and rat_version_corr){
            originReactorVect[0] = "SHIN";
        }
        if (originReactorVect[0] == "SHIN-WOLSONG" and rat_version_corr){
            originReactorVect[0] = "SHIN";
            if(originReactorVect[1] == "0"){
                originReactorVect[1] = "3";
            }
            else if(originReactorVect[1] == "1"){
                originReactorVect[1] = "4";
            }
        }
        if(originReactorVect[0] == "FUKUSHIMA-DAINI" and rat_version_corr){
            originReactorVect[0] = "FUKUSHIMA";
        }
        if(originReactorVect[0] == "HIGASHI DORI-1 (TOHOKU)" and rat_version_corr){
            originReactorVect[0] = "HIGASHI DORI";
        }
        if((originReactorVect[0] == "GUNDREMMINGEN-C" or originReactorVect[0] == "GUNDREMMINGEN-B") and rat_version_corr){
            originReactorVect[0] = "GUNDREMMINGEN";
        }

        linkdb = db->GetLink("REACTOR_STATUS",originReactorVect[0]);
        try{
            efficiency = linkdb->GetDArray("core_power_scale_factor")[std::stoi(originReactorVect[1])];
        }
        catch(...){
            num_no_table_evts++;
            std::cout << "Reactors status table does not exist for reactor " << originReactorString << ". Removing this event." << std::endl;
            continue;
        }
        if(efficiency > 1){ //we cant add events, so for now make a note of which reactors this occurs, SOLUTION: simulate higher than design power
            num_excess_eff_evts++;
            if(std::find(reacts_eff_names.begin(), reacts_eff_names.end(),originReactorString)==reacts_eff_names.end()){
                reacts_eff_names.push_back(originReactorString);
            }
        }
        if(rndm<efficiency){ //if keeping the event, write to new ttree
	    previousmcIndex = mcIndex;
            scaledReactorEventInfo->Fill();
            spectraScaled->Fill(parentKE1);
            if(valid){
                spectraScaledRecon->Fill(reconEnergy);
            }
        }
    }

    //write and close everything, give user info

    Int_t nentries_scaled = scaledReactorEventInfo->GetEntries();
    std::cout << "Raw entries: " << nentries << ", scaled entries: " << nentries_scaled << ", removed entries: " << nentries-nentries_scaled << std::endl;
    if(num_excess_eff_evts != 0){
        std::cout << "The following reactors had efficiencies > 1, events were not removed:" << std::endl;
        for(int i=0;i<reacts_eff_names.size();i++){
            std::cout << reacts_eff_names.at(i) << std::endl;
        }
        std::cout << "This was the case for " << num_excess_eff_evts << " events, equal to " << (float(num_excess_eff_evts) / nentries) * 100 << " percent of events" << std::endl;
    }
    if(num_no_table_evts != 0){
        std::cout << "Could not find a table for " << num_no_table_evts << " events, equal to " << (float(num_no_table_evts) / nentries) * 100 << " percent of events" << std::endl;
    }

    //If oscillating, do it here
    if(oscillations){
        std::cout << "Applying oscillations" << std::endl;
        Double_t nuE_parent;
        Double_t nuE;
        TString *originReactorOsc = NULL;
        Bool_t *valid_osc;
        //create ntuple
	std::string output_string = "scaled_oscillated_" + inputFile;
        TFile *scaled_oscillated_ntuple = new TFile(output_string.c_str(),"RECREATE");
        TH1 *spectraScaledOscillated = new TH1D("oscillated", "Antinu Scaled & Oscillated Spectra", 100, 0, 10);
        TH1 *spectraScaledOscillatedRecon = new TH1D("oscillated_recon", "Antinu Scaled & Oscillated Spectra", 100, 0, 10);
        TTree *scaledOscillatedReactorEventInfo = scaledReactorEventInfo->CloneTree(0);
        scaledReactorEventInfo->SetBranchAddress("parentKE1",&nuE_parent);
        scaledReactorEventInfo->SetBranchAddress("energy",&nuE);
        scaledReactorEventInfo->SetBranchAddress("parentMeta1",&originReactorOsc);
        scaledReactorEventInfo->SetBranchAddress("fitValid",&valid_osc);
        //get parameters
        linkdb = db->GetLink("OSCILLATIONS");
        Double_t fDmSqr21 = linkdb->GetD("deltamsqr21");
        Double_t fDmSqr32 = linkdb->GetD("deltamsqr32");
        Double_t fSSqrTheta12 = linkdb->GetD("sinsqrtheta12");
        Double_t fSSqrTheta13 = linkdb->GetD("sinsqrtheta13");

        // Declare quantities to use in loops
        Double_t scale;
        Double_t fOscProb = 0.0;

        // If electron density is zero, use vacuum calculation
        if(Ne == 0.0){
            Double_t fSSqr2Theta12 = 4.0*fSSqrTheta12*(1.0-fSSqrTheta12);
            Double_t fSSqr2Theta13 = 4.0*fSSqrTheta13*(1.0-fSSqrTheta13);
            Double_t fC4 = pow(1.0-fSSqrTheta13, 2.0);

            // Declare paramters outside loop
            Double_t sSqrDm21BE;
            Double_t sSqrDm31BE;
            Double_t sSqrDm32BE;

	    previousmcIndex = -999;
            //loop and oscillate
            for(int a=0; a<nentries_scaled;a++){
                scaledReactorEventInfo->GetEntry(a);

		// Save all the retriggers and delayed event of the saved prompt event
		if( mcIndex == previousmcIndex ){ scaledOscillatedReactorEventInfo->Fill(); }

		// Only apply scaling to the first trigger of the mcIndex, which is the most likely to be the prompt
		// Skip all events that are not propmt
		if( evIndex != 0 ){ continue; }
		
                //get baseline
                const char *originReactorString = originReactorOsc->Data();
                std::vector<std::string> originReactorVect = SplitString(originReactorString, 1);

                if(originReactorVect[0] == "LENINGRAD-2" and rat_version_corr){ //hacks for rat versions 6.18.4 < 6.18.7
                    originReactorVect[0] = "LENINGRAD 2";
                }
                if(originReactorVect[0] == "SHIN-KORI" and rat_version_corr){
                    originReactorVect[0] = "SHIN";
                }
                if (originReactorVect[0] == "SHIN-WOLSONG" and rat_version_corr){
                    originReactorVect[0] = "SHIN";
                    if(originReactorVect[1] == "0"){
                        originReactorVect[1] = "3";
                    }
                    else if(originReactorVect[1] == "1"){
                        originReactorVect[1] = "4";
                    }
                }
                if(originReactorVect[0] == "FUKUSHIMA-DAINI" and rat_version_corr){
                    originReactorVect[0] = "FUKUSHIMA";
                }
                if(originReactorVect[0] == "HIGASHI DORI-1 (TOHOKU)" and rat_version_corr){
                    originReactorVect[0] = "HIGASHI DORI";
                }
                if((originReactorVect[0] == "GUNDREMMINGEN-C" or originReactorVect[0] == "GUNDREMMINGEN-B") and rat_version_corr){
                    originReactorVect[0] = "GUNDREMMINGEN";
                }

                linkdb = db->GetLink("REACTOR",originReactorVect[0]);
                std::vector<Double_t> fLatitude  = linkdb->GetDArray("latitude");
                std::vector<Double_t> fLongitute = linkdb->GetDArray("longitude");
                std::vector<Double_t> fAltitude = linkdb->GetDArray("altitude");
                const double baseline = GetReactorDistanceLLA(fLongitute[std::stoi(originReactorVect[1])], fLatitude[std::stoi(originReactorVect[1])], fAltitude[std::stoi(originReactorVect[1])]);

                //do calculation
                scale = 1.267e3 * baseline / nuE_parent; // for nuE in [MeV] and baseline in [km]
                sSqrDm21BE = pow(sin(scale * fDmSqr21), 2.0);
                sSqrDm31BE = pow(sin(scale * (fDmSqr32-fDmSqr21)), 2.0);
                sSqrDm32BE = pow(sin(scale * fDmSqr32), 2.0);
                fOscProb = fC4 * fSSqr2Theta12 * sSqrDm21BE + fSSqr2Theta13 * ( (1.0-fSSqrTheta12)*sSqrDm31BE + fSSqrTheta12*sSqrDm32BE );
                //decide whether to keep or not
                rndm = rndm_dummy->Rndm();
                if(rndm > fOscProb){
		    previousmcIndex = mcIndex;
                    scaledOscillatedReactorEventInfo->Fill();
                    spectraScaledOscillated->Fill(nuE_parent);
                    if(valid_osc){
                        spectraScaledOscillatedRecon->Fill(nuE);
                    }
                }
            }
        // Otherswsie take into account matter effects
        } else {
            // Compute fixed parameters for oscillation calculation
            Double_t fDmSqr31 = fDmSqr32 + fDmSqr21;

            Double_t H_ee_vac = fDmSqr21 * (fSSqrTheta12 * (1 - fSSqrTheta13) - (1.0/3.0)) + fDmSqr31 * (fSSqrTheta13 - (1.0/3.0));
            Double_t H_neq2 = (1 - fSSqrTheta13) * (fDmSqr21*fDmSqr21 * fSSqrTheta12 * (1 + fSSqrTheta12 * (fSSqrTheta13 - 1))
                            + fDmSqr31*fDmSqr31 * fSSqrTheta13 - 2.0 * fDmSqr21 * fDmSqr31 * fSSqrTheta12 * fSSqrTheta13);

            Double_t a0_vac = - (2.0/27.0) * (fDmSqr21*fDmSqr21*fDmSqr21 + fDmSqr31*fDmSqr31*fDmSqr31)
                                + (1.0/9.0) * (fDmSqr21*fDmSqr21 * fDmSqr31 + fDmSqr21 * fDmSqr31*fDmSqr31);
            Double_t a1_vac = (1.0/3.0) * (fDmSqr21 * fDmSqr31 - fDmSqr21*fDmSqr21 - fDmSqr31*fDmSqr31);
            Double_t Y_ee_vac = (2.0/3.0) * a1_vac + H_ee_vac*H_ee_vac + H_neq2;

            // Declare quantities to use in loop
            Double_t alpha = - 2.535e-31 * Ne;  // conversion factor in eV2/MeV
            Double_t A_CC;
            Double_t alpha_1;
            Double_t a0;
            Double_t a1;
            Double_t Y_ee;
            Double_t H_ee;
            Double_t eigen[3];
            Double_t X[3];
            Double_t arcCos;
            Double_t preFact;
            Double_t s_10;
            Double_t s_20;
            Double_t s_21;

            //loop and oscillate
	    previousmcIndex = -999;
            for(int a=0; a<nentries_scaled;a++){
                scaledReactorEventInfo->GetEntry(a);

		// Save all the retriggers and delayed event of the saved prompt event
		if( mcIndex == previousmcIndex ){ scaledOscillatedReactorEventInfo->Fill(); }

		// Only apply scaling to the first trigger of the mcIndex, which is the most likely to be the prompt
		// Skip all events that are not propmt
		if( evIndex != 0 ){ continue; }
		
                //get baseline
                const char *originReactorString = originReactorOsc->Data();
                std::vector<std::string> originReactorVect = SplitString(originReactorString, 1);

                if(originReactorVect[0] == "LENINGRAD-2" and rat_version_corr){ //hacks for rat versions 6.18.4 < 6.18.7
                    originReactorVect[0] = "LENINGRAD 2";
                }
                if(originReactorVect[0] == "SHIN-KORI" and rat_version_corr){
                    originReactorVect[0] = "SHIN";
                }
                if (originReactorVect[0] == "SHIN-WOLSONG" and rat_version_corr){
                    originReactorVect[0] = "SHIN";
                    if(originReactorVect[1] == "0"){
                        originReactorVect[1] = "3";
                    }
                    else if(originReactorVect[1] == "1"){
                        originReactorVect[1] = "4";
                    }
                }
                if(originReactorVect[0] == "FUKUSHIMA-DAINI" and rat_version_corr){
                    originReactorVect[0] = "FUKUSHIMA";
                }
                if(originReactorVect[0] == "HIGASHI DORI-1 (TOHOKU)" and rat_version_corr){
                    originReactorVect[0] = "HIGASHI DORI";
                }
                if((originReactorVect[0] == "GUNDREMMINGEN-C" or originReactorVect[0] == "GUNDREMMINGEN-B") and rat_version_corr){
                    originReactorVect[0] = "GUNDREMMINGEN";
                }

                linkdb = db->GetLink("REACTOR",originReactorVect[0]);
                std::vector<Double_t> fLatitude  = linkdb->GetDArray("latitude");
                std::vector<Double_t> fLongitute = linkdb->GetDArray("longitude");
                std::vector<Double_t> fAltitude = linkdb->GetDArray("altitude");
                const double baseline = GetReactorDistanceLLA(fLongitute[std::stoi(originReactorVect[1])], fLatitude[std::stoi(originReactorVect[1])], fAltitude[std::stoi(originReactorVect[1])]);

                //do calculation
                scale = 1.267e3 * baseline / nuE_parent; // for nuE in [MeV] and baseline in [km]
                A_CC = alpha * nuE_parent; // for A_CC in [eV^2] and nuE in [MeV]
                
                // Compute new values for H_ee, Y, a0 and a1 (make sure and Y are updated after their use by others)
                alpha_1 = H_ee_vac * A_CC + (1.0/3.0) * A_CC*A_CC;

                a0 = a0_vac - Y_ee_vac * A_CC - (1.0/3.0) * H_ee_vac * A_CC*A_CC - (2.0/27.0) * A_CC*A_CC*A_CC;
                a1 = a1_vac - alpha_1;
                Y_ee = Y_ee_vac + (2.0/3.0) * alpha_1;
                H_ee = H_ee_vac + (2.0/3.0) * A_CC;

                // Get eigenvalues of H, and constants X and theta
                arcCos = (1.0/3.0) * acos(1.5 * (a0/a1) * sqrt(- 3.0 / a1));
                preFact = 2.0 * sqrt(- a1 / 3.0);

                for(int i=0; i<3; ++i){
                    eigen[i] = preFact * cos(arcCos - (2.0/3.0) * M_PI * i);
                    X[i] = (1.0/3.0) + (eigen[i] * H_ee + Y_ee) / (3.0 * eigen[i]*eigen[i] + a1);
                }

                s_10 = sin(scale * (eigen[1] - eigen[0]));
                s_20 = sin(scale * (eigen[2] - eigen[0]));
                s_21 = sin(scale * (eigen[2] - eigen[1]));

                // Compute probability
                fOscProb = 4.0 * (X[1]*X[0]*s_10*s_10 + X[2]*X[0]*s_20*s_20 + X[2]*X[1]*s_21*s_21);

                //decide whether to keep or not
                rndm = rndm_dummy->Rndm();
                if(rndm > fOscProb){
		    previousmcIndex = mcIndex;
                    scaledOscillatedReactorEventInfo->Fill();
                    spectraScaledOscillated->Fill(nuE_parent);
                    if(valid_osc){
                        spectraScaledOscillatedRecon->Fill(nuE);
                    }
                }
            }
        }

        
        Int_t nentries_oscillated = scaledOscillatedReactorEventInfo->GetEntries();
        std::cout << "Scaled entries: " << nentries_scaled << ", oscillated entries: " << nentries_oscillated << ", removed entries: " << nentries_scaled-nentries_oscillated << std::endl;
        scaled_oscillated_ntuple->cd();
        scaledOscillatedReactorEventInfo->Write();
        scaled_hist->cd();
        spectraScaledOscillated->Write();
        spectraScaledOscillatedRecon->Write();
        delete scaledOscillatedReactorEventInfo;
        scaled_oscillated_ntuple->Close();
    }

    scaled_ntuple->cd();
    scaledReactorEventInfo->Write();
    scaled_hist->cd();
    spectraRaw->Write();
    spectraScaled->Write();
    spectraRawRecon->Write();
    spectraScaledRecon->Write();

    delete reactorEventInfo;
    delete scaledReactorEventInfo;
    reactorFile->Close();
    scaled_ntuple->Close();
    scaled_hist->Close();
    return 0;
}

//main

int main(int argv, char** argc){
    if(argv<6){
        std::cout << "Not enough input arguments. Input the file name, a bool for working in offline mode, a bool for oscillation, a bool to specify if using rat versions between 6.18.4 and 6.18.7, a value for electron density (0 for vacuum oscillation) and a ratdb tag." << std::endl;
        return 1;
    }
    bool offline_input;
    if(std::string(argc[2])=="true" || std::string(argc[2]) == "True" || std::string(argc[2]) == "TRUE" || std::string(argc[2]) == "1"){ 
        offline_input = true;
    }
    else{ //don't work in offline mode by default
        offline_input = false;
    }
    bool oscillation_input;
    Double_t Ne;
    if(std::string(argc[3])=="true" || std::string(argc[3]) == "True" || std::string(argc[3]) == "TRUE" || std::string(argc[3]) == "1"){ 
        oscillation_input = true;
        if(std::string(argc[5])=="crust") {
            // Average crust electron density (cm-3), given by density 2.7g/cm^3, and <proton/nucleon>=0.5.
            // (If you wish modify these, use the fact that A_CC is proportional to both quantities)
            Ne = 8.13e23;
        }
        else{
            Ne = atof(argc[5]);  // Otherwise, electron density is inputted
        }
    }
    else{ //don't oscillate by default
        oscillation_input = false;
    }
    bool reactor_name_input;
    if(std::string(argc[4])=="true" || std::string(argc[4]) == "True" || std::string(argc[4]) == "TRUE" || std::string(argc[4]) == "1"){ 
        reactor_name_input = true;
    }
    else{ //don't alter reactor names by default
        reactor_name_input = false;
    }
    std::string ratDBTag_input;
    if(argv >= 7){
        ratDBTag_input = argc[6];
    }
    else{
        std::cout << "Please specify a ratdb tag from the list at snopl.us/ratdb or pass 'liveDB' to use the current state ratdb (untagged)." << std::endl;
        std::cout << "Exiting..." << std::endl;
        return 1;
    }
    std::cout << "Scaling reactor fluxes." << std::endl;
    int return_code = ScaleReactorFlux(argc[1], offline_input, oscillation_input, reactor_name_input, Ne, ratDBTag_input);
    if(return_code != 0){
        std::cout << "Exiting..." << std::endl;
        return 1;
    }
    else{
        std::cout << "Reactor fluxes scaled sucessfully." << std::endl;
        return 0;
    }
} 
