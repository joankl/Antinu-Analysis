#ifndef __CINT__ 

#include "TROOT.h"
#include "TObjArray.h"
#include "TString.h"
#include "TRegexp.h"
#include "TObjString.h"
#include "TH1D.h"
#include "TLine.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TSystem.h"
// include used RAT classes below
#include <RAT/DS/EV.hh> 
//#include <RAT/DS/Root.hh> 
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCParticle.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DB.hh>
#include <RAT/DBTable.hh>
#include <RAT/DBLink.hh>
// include STL classes below (C++ headers)
#include <iostream>
#include <cmath>
#include <stdint.h>
#include <string>
#include <RAT/GeoUtils.hh>
#include <sstream>
#include <fstream>
#include <string>

#endif

#include <stdint.h>

void Analysis_antinu(TString filename, string prefix, string folder, Double_t FV, TString outputname)
{
  //define common pattern in name of files to be merged

  string dirs = "/share/neutrino/snoplus/Data/FullFill_2p2/"+prefix+"/"+folder+"/";
  const char *indir = dirs.c_str();

  cout<<indir<<endl;
  ifstream input;
  char name[200];
  input.open(filename);	
  int counter = 0;
  TObjString *fn = NULL;
  TObjArray *flist = new TObjArray();
  const char *file;
  TString fname;
  // int size = pass;
  int j = 0;

  while(1){
    input>>name;
    j++;
    string namen = name;
    cout<<namen<<endl;
    string pattern = prefix+"_r0000"+namen;
    
    if (!indir) {
      cerr << "Failed to get directory name." << endl;
      return;
    }
    gSystem->ChangeDirectory(indir);
    cout << "Working in directory : " << gSystem->WorkingDirectory() << endl;
    
    void *dir = gSystem->OpenDirectory(".");
    if (!dir) {
      cerr << "Couldn't open directory : " << indir << endl;
      return;
    }
    
    TRegexp regexp(pattern.c_str());
    
    while ((file = gSystem->GetDirEntry(dir))) {
      fname = file;
      //ignore subfolders and existing merged files
      
      if (!fname.CompareTo(".") || !fname.CompareTo("..") ||  (fname.Contains("cc"))) {
      //      std::cout << "Ignoring file " << fname << std::endl;
	continue;
      }
      if (fname.Index(regexp) != kNPOS) {
	std::cout << "Adding file " << fname << std::endl;
	flist->Add(new TObjString(fname));
      }
    }
    if (!input.good()) break;
  }

  gSystem->FreeDirectory(dir);
  flist->Print("*");
  
  //loop over array and add file to the TChain
  tc = new TChain("output");
  TIter next(flist);
  while ((fn = (TObjString*)next())){
    cout << "__ Adding file " << (fn->GetString()).Data() << endl;
    tc->Add((fn->GetString()).Data());
    counter++;
  } 

  input.close();

  Long_t iEntries = tc->GetEntries();
  cout << "Created chain of " << counter<< " with: "<< iEntries << " entries" << endl;
  
 // TString outputname;
  TFile *PF = new TFile(outputname+=".root", "recreate");
  
  RAT::DU::Utility::Get()->LoadDBAndBeginRun();

  RAT::DB* db = RAT::DB::Get();
  
  RAT::DS::Run Run;
  Run.SetRunID(300999);

  db->BeginOfRun(Run);

  RAT::DBLinkPtr dblink = db->GetLink("AV_OFFSET_RUNTIME");

  std::vector<Double_t> AVpos = dblink->GetDArray("position");

  Double_t Zoffset = AVpos[2];
  cout << "AV Offset is " << Zoffset << " mm" << endl;
  
 //// --------  Plots

  TH1D *Total_nocuts = new TH1D("Total_nocuts","Energy", 600,0.,6.0);
  Total_nocuts->SetXTitle("Energy [MeV]");
  Total_nocuts->SetYTitle("events/10keV");

  TH1D *Total_dc = new TH1D("Total_dc","Energy", 600,0.,6.0);
  Total_dc->SetXTitle("Energy [MeV]");
  Total_dc->SetYTitle("events/10keV");

  TH1D *Total_hscut = new TH1D("Total_hscut","Energy", 600,0.,6.0);
  Total_hscut->SetXTitle("Energy [MeV]");
  Total_hscut->SetYTitle("events/10keV");

  TH1D *Total_nhighcut = new TH1D("Total_nhighcut","Energy", 600,0.,6.0);
  Total_nhighcut->SetXTitle("Energy [MeV]");
  Total_nhighcut->SetYTitle("events/10keV");

  TH1D *Total_dtcut = new TH1D("Total_dtcut","Energy", 600,0.,6.0);
  Total_dtcut->SetXTitle("Energy [MeV]");
  Total_dtcut->SetYTitle("events/10keV");

  TH1D *Total_prompt = new TH1D("Total_prompt","Energy", 600,0.,6.0);
  Total_prompt->SetXTitle("Energy [MeV]");
  Total_prompt->SetYTitle("events/10keV");

  TH1D *Total_delay = new TH1D("Total_delay","Energy", 600,0.,6.0);
  Total_delay->SetXTitle("Energy [MeV]");
  Total_delay->SetYTitle("events/10keV");

  TH1D *Total_dt = new TH1D("Total_dt","Time", 2000,0.,2000.);
  Total_dt->SetXTitle("Time diff [us]");
  Total_dt->SetYTitle("events/1 us");

  TH1D *Total_dr = new TH1D("Total_dr","Position diff", 1500,0.,1500.0);
  Total_dr->SetXTitle("Pos diff [mm]");
  Total_dr->SetYTitle("events/1mm");




  //////////////////////
  // Declare Variables
  //////////////////////

  Double_t Zoffset;        // AV zOffset for the run
  Double_t RunDuration;    // Duration of the run in seconds
  Double_t TotalRunTime;   // Total duration of the runs in seconds (different from vRunDuration when the number of input runs is more than one)
  Double_t nhits1;
  Double_t nhits0;
  Double_t itr0;
  Double_t posr0;
  Double_t posx0;
  Double_t posy0;
  Double_t posz0;
  Double_t Dr;
  Double_t rho;
  Double_t rho0;
  Double_t r3;
  Double_t r30;
  Double_t energy0, energy1;
  Double_t energyF;
  Double_t scaledLLH;

  Bool_t nhexists;     
  Bool_t vetocut, vetocutR;    
  Bool_t pair;        
  Bool_t hs;      
      
  ULong64_t timepev;   
  ULong64_t time0; 
  ULong64_t time;      
  ULong64_t timeF;    
  ULong64_t timeFj;  
  ULong64_t timeDiffj;    
  ULong64_t time1;
  ULong64_t Dt;
  ULong64_t timeDiff, timeDiffR, timeDiffRF;

  Int_t IDc;
  Int_t ID;

  // Counters
  Int_t k_e, nh, rt, tag;  // Counter for the number of tagged bipos

  // From the ttree
  Int_t runID;
  Int_t  evIndex;
  Int_t  mcIndex;
  Int_t eventID;
  Bool_t mc;
  Int_t pdg1;
  Int_t pdg2;
  Int_t parentpdg1;
  Int_t parentpdg2;
  Double_t mcPosr;
  Double_t mcPosx;
  Double_t mcPosy;
  Double_t mcPosz;
  Int_t     nhits;
  Int_t     nhitsCleaned;
  Double_t correctedNhits;
  Int_t necknhits;
  Int_t necknhits0;
  Int_t runID;
  ULong64_t clockCount50;
  ULong64_t dcApplied;
  ULong64_t dcFlagged;
  Int_t triggerWord;
  Bool_t   fitValid;
  Bool_t   scintFit;
  Double_t itr;
  Double_t beta14;
  Double_t beta140;
  Double_t posx;
  Double_t posy;
  Double_t posz;
  Double_t posr;
  Double_t energy;
  Double_t alphaBeta212;
  Double_t alphaBeta214;
  Double_t posFOM;
  Int_t ROIID, ID0;

  tc->SetBranchAddress("evIndex",&evIndex);
  tc->SetBranchAddress("nhits",&nhits);
  tc->SetBranchAddress("energy",&energy);
  tc->SetBranchAddress("scintFit",&scintFit);
  tc->SetBranchAddress("beta14",&beta14);
  tc->SetBranchAddress("itr",&itr);
  tc->SetBranchAddress("posr",&posr);
  tc->SetBranchAddress("posx",&posx);
  tc->SetBranchAddress("posy",&posy);
  tc->SetBranchAddress("posz",&posz);
  tc->SetBranchAddress("fitValid",&fitValid); 
  tc->SetBranchAddress("dcFlagged",&dcFlagged);
  tc->SetBranchAddress("triggerWord",&triggerWord);
  tc->SetBranchAddress("dcApplied",&dcApplied);
  tc->SetBranchAddress("nhitsCleaned",&nhitsCleaned);
  tc->SetBranchAddress("clockCount50",&clockCount50);
  tc->SetBranchAddress("eventID",&eventID);
  tc->SetBranchAddress("necknhits",&necknhits);
  tc->SetBranchAddress("correctedNhits",&correctedNhits);
  tc->SetBranchAddress("runID",&runID);
  tc->SetBranchAddress("posFOM",&posFOM);

  nhexists = false;
  hs = false;
  vetocut = false;
  vetocutR = false;
  nh = rt = tag = 0;
  pair = false;

  // now the code

  cout<<" now start the code "<<endl;

  // Loop through all of them
  for(Long_t i=0; i<iEntries;++i)
    {
      //previous event
      if ( i > 0 )
	{
	  tc->GetEntry(i-1);
	  timepev = clockCount50 * 20.0;
	}

      //the actual event

      tc->GetEntry(i);

      //define quantities

      Zoffset = 183.75;
      ROIID = 7141281;

      posz = posz - Zoffset;
      posr = sqrt(posx*posx+posy*posy+posz*posz);
      rho = sqrt(posx*posx + posy*posy)/1000;
      r3 = posr*posr*posr/(6005*6005*6005.);
      ID = eventID;
      time = clockCount50 * 20.0;
      pair = false;

      //  if ( ID == IDc ) {IDc = 0; continue;}

      if ( nhitsCleaned < 20 ) continue;

      if (  posr < FV ) Total_nocuts->Fill( energy ); //to study the dbd

      if (((0x2100000042C2) & dcFlagged ) != (0x2100000042C2)) continue; 

      if ( posr < FV ) Total_dc->Fill( energy );

      //---- veto ----
      /*     if ( nhits > 3000 & !nhexists ) //was 3000
	{
	  nh++;
	  if ( ( posx > 500 && posx < 900 && posy > -900 && posy < -500 && posz > 5300 ) || (posx > 0 && posx < 500 && posy > -1250 && posy < -800 && posz > 5300) ) 
	    {
	      hs = true;
	    }
      
	  nhexists = true;
	  timeF = clockCount50 * 20.0; 
	}
  
      if ( nhexists )
	{
	  timeDiff = ( time - timeF ) & 0x7FFFFFFFFFF;
	  if ( hs == true )
	    {
	      if ( timeDiff < 1.0e9 ) {vetocut = true;}
	      else { vetocut = false; nhexists = false; hs = false; continue; }
	    }
	  if ( hs == false )
	    {
	      if ( timeDiff < 20.0e9 ) vetocut = true; //mine was 20s, Sam is 5ms
	      else { vetocut = false; nhexists = false; }
	    }
	}

      else if ( !nhexists )
      	{
	  
      	  timeDiffR = ( time - timepev ) & 0x7FFFFFFFFFF;
      	  if ( timeDiffR < 1000 ) {vetocutR = true; rt++;} //mine was 500
      	  else {vetocutR = false;}
      	}

      
      if ( vetocut ) continue;

      if ( ID == ROIID ) cout<<" still here "<<endl;

      if ( posr < FV ) Total_hscut->Fill( energy );

      if ( vetocutR ) continue;
      */
      //check if a retrigger

      if ( fitValid != 1 ) continue;
      if ( posr >= FV ) continue;

      if ( energy < 0.9 || energy > 8.0 ) continue;

      Total_dtcut->Fill( energy );

      posx0 = posx;
      posy0 = posy;
      posz0 = posz;
      rho0 = rho;
      energy0 = energy;
      itr0 = itr;
      beta140 = beta14;
      necknhits0 = necknhits;
      nhits0 = nhits;
      time0 = time;
      ID0 = ID;


      //cout<<"here ene cut"<<endl;

      //select the delay

      for( Long_t n = i+1; n<iEntries; n++)
	{
	  if (pair) break;
	  tc->GetEntry(n);
	  if ( energy < 1.8 || energy > 3.0 ) continue;// || energy > 10.0 ) continue; //mine only <0.1, Sam 0.5
          time1 = clockCount50 * 20.0;
	  posz = posz - Zoffset;
	  posr = sqrt(posx*posx+posy*posy+posz*posz);
	  Dt = ( ( time1 - time0 ) & 0x7FFFFFFFFFF );
	  Dr = sqrt( (posx-posx0)*(posx-posx0) + (posy-posy0)*(posy-posy0) + (posz-posz0)*(posz-posz0) );
	 // if ( nhits > 40 && Dt <=460 ) { pair = true; IDc = eventID; break; } //I didn't have this line
	  if ( Dt < 500 ) continue;
          if ( Dt > 2.0e6 ) break;
          if ( Dr > 1000 ) continue; //mine 2.5m, Sam 1.5m
	      tag++;
	      pair = true;
	      IDc = eventID;
	      energy1 = energy;
	      cout<<energy1<<endl;
	      Total_prompt->Fill( energy0 );
	      Total_delay->Fill( energy1 );
	      Total_dt->Fill( Dt/1e3 );
	      Total_dr->Fill( Dr );
	      //   if ( ID0 == ROIID ) cout<<" energy: "<<energy<<" Dr "<<Dr<<" Dt "<<Dt<<" ID: "<<IDc<<endl;
	}
      

	  // scaledLLH = posFOM/nhits;
	  //    cout<<"here"<<endl;
    }
  
  cout<<" tagged: "<<tag<<endl; 
  PF->Write();
  PF->Close();

}


