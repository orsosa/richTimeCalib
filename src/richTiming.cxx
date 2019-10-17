#include "Riostream.h"
#include "TApplication.h"
#include "TBenchmark.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
using namespace std;

//#include<TFile.h>
//#include<TMath.h>
#include<TTree.h>
#include<TChain.h>
#include<TH2F.h>
#include<TF1.h>
#include<TVector3.h>
#include<TLorentzVector.h>
#include<TCanvas.h>
#include<TApplication.h>
#include<TRint.h>
#include<TStyle.h>
//#include <TROOT.h>

/****************************************/
/* CLAS12 Bank definition */
#include <Clas12Banks4.h>


/* ================================================================= */
void PrintUsage(char *processName);

/* ================================================================ */
/* MAX number of hipo files */
#define MAXFILES 10000

/* For time correction */
int iTimeCorr = 0;

/* For RICH event solution */
int iRayTracing = 0;

/* For the event start time correction */
int iSTtimeCorr = 0;

/* Average uncorrect DTime for histograms */
double HistoAveDT = -90;

/* Flag for Zero Filed run analysis */
int IsZeroField = 0;

/* ================================================================ */
#include <RichHW.h>
#include "Clas12Detectors.h"
#include "Clas12Pid.h"
#include <RichTimeCorr.h>

/* ================================================================ */

/* Histograms */
void makeHistos();
TObjArray Hlist(0);

/******************************************/
int main(int argc, char** argv) {

  /* number of input files */
  int nFiles = 0;
  /* max number of entries to read */
  int nEntries = 0;

  int RunNumber = 0;
  
  /* list of input file names */
  char inputFile[256];
  char *inputFiles[MAXFILES];

  /* output root file */
  char rootFile[256];

  /* ======================================== */
  /* Scanning the command line */
  char *argptr = NULL;
  
  if(argc == 1) {
    PrintUsage(argv[0]);
    exit(0);
  }
  
  
  for (int i=1; i<argc; i++) {
    argptr = argv[i];
    //cout << argv[i] << endl;
    
    if (*argptr == '-') {
      argptr++;

      switch (*argptr) {
	
	
      case 'n':
	nEntries = atoi(++argptr);
	printf("Reading %d entries per file\n", nEntries);
	break;
      case 'R':
	RunNumber = atoi(++argptr);
	printf("Run number: %d\n", RunNumber);
	break;
      case 'r':
	iRayTracing = 1;
	printf("Using Ray Tracing solution\n");
	break;
      case 's':
	iSTtimeCorr = 1;
	printf("Correcting the event start time by hand\n");
	break;
      case 'T':
	iTimeCorr = atoi(++argptr);
	printf("Applying time correction with flag %d\n", iTimeCorr);
	break;
      case 't':
	HistoAveDT = atof(++argptr);
	printf("DTime histogram average: %lf\n", HistoAveDT);
	break;
      case 'P':
	RichParticleID = atoi(++argptr);
	printf("Using particles with PID %d\n", RichParticleID);
	break;
      case 'Z':
	IsZeroField = 1;
	printf("Analysis of Zero Field data\n");
	break;

      default:
	fprintf(stderr, "Unrecognized argument: [-%s]\n\n", argptr);
	PrintUsage(argv[0]);
	break;
      }

    }
    else {
      
      sprintf(inputFile,"%s",argptr);

      inputFiles[nFiles] = (char*)malloc(256*sizeof(char));
      sprintf(inputFiles[nFiles], "%s", inputFile);
      nFiles++;
    }
    
//printf("argptr: %s   runNumber=%d\n", argptr, RunNumber);
  }

  if (RunNumber == 0) {
    printf("Please give the run number\n");
    return -1;
  }
  
  sprintf(rootFile, "RichTimeCalib_%d_%d.root", RunNumber, iTimeCorr);
  printf("Reading %d files\n", nFiles);
  printf("Root file name: %s\n", rootFile);


  /* ===================================== */
  /* Histograms */
  makeHistos();
  
  /* ===================================== */
  /* Time corrections */
  InitTimeCorrections();
  if (iTimeCorr) {
    LoadTimeOffsets();
    if (iTimeCorr == 2) LoadTimeWalkPars();
  }

  /* ===================================== */
  /* Some counters */
  int entry = 0;
  int nTriggers = 0;
  int nRec = 0;
  int nRichBadClusters = 0;
  int nRichTracks = 0;
  int nRich = 0;
  int nphotons = 0;


  // int tile;
  int pmt, anode, absChannel;
  //TH1F *h1;
  TH2F *h2;
  char name[200];

  
  /* ===================================== */
  /* hipo4 object inizializations */
  fEvent = new hipo::event();
  fFactory = new hipo::dictionary();
  fReader = new hipo::reader();
 

  
  /* ====================================== */
  /* LOOP OVER THE HIPO FILES */
  for (int l=0; l<nFiles; l++) {
    printf("==>> READING HIPO FILE: %s\n", inputFiles[l]);
 
    /* opening the hipo4 file */
    //fReader->open(flist[f]);  TString input
    fReader->open(inputFiles[l]);

    /* Bank definition */
    fReader->readDictionary(*fFactory);
    InitBanks();


    /* looping over the current file */
    while( (fReader->next() )  && ( (entry < nEntries)||(nEntries == 0) ) ) {
      fReader->read(*fEvent);

       /* looking for CLAS12 banks */
      if (FillBanks() ) {


	/* looking for good triggers */
	if (RUN__config->getRows()) {
	  get_RUN__config(0);
	  nTriggers++;

	  if (REC__Event->getRows()) {
	    get_REC__Event(0);
	    nRec++;

	    /* time stamp correction */
	    double ts = 0;
	    if ( (RUN__config_timestamp % 2) == 1) ts = 4.;


	    /* Loading the pindex Map */
	    LoadPindexMap(REC__Particle->getRows());


	    
	    /* Selecting events with 1 track in the RICH */
	    if (RICH__hadrons->getRows() == 1) {
	      int richhadron_index = 0;

	      /* Selecting the good particle in the RICH */
	      int RichParticleQ = 0;
	      if (IsZeroField == 0) {
		RichParticleQ = GoodRichParticle(richhadron_index);
	      }
	      else if (IsZeroField == 1) {
		RichParticleQ = GoodRichParticle_ZeroField(richhadron_index);
	      }
	      if (RichParticleQ > 0) {
		nRichTracks++;
	      
		/* ------------------------------ */
		/* loop over the RICH photons */
		nphotons = 0;
		for (int f=0; f<RICH__photons->getRows(); f++) {
		  get_RICH__photons(f);
		  int richhit_pindex = RICH__photons_hit_index;
		  //int richhadron_pindex = RICH__photons_hadron_index;

		  /* good ray tracing photon, if selected */
		  int GoodPhoton = 1;
		  if (iRayTracing) {
		    if ( (RICH__photons_traced_the == 0) || (RICH__photons_traced_phi == 0) || (RICH__photons_traced_EtaC == 0) ) GoodPhoton = 0; 
		  }

		  if ( (RICH__photons_type == 0) && GoodPhoton) {
		    nphotons++;

		    /* channel info */
		    get_RICH__hits(richhit_pindex);
		    //tile = RICH__hits_tile;
		    pmt = RICH__hits_pmt;
		    anode = RICH__hits_anode;
		    absChannel = anode + (pmt-1)*nANODES;

		    /* signal duration */
		    double duration = RICH__hits_duration;

		    /* Measured photon time */
		    double MeasPhotonTime = RICH__hits_rawtime + ts;

		    /* Photon path time from production to the MAPMT */
		    double PhotonPathTime = RICH__photons_analytic_time;
		    if (iRayTracing) PhotonPathTime = RICH__photons_traced_time;

		    /* Calculated photon time with respect to the event start time */
		    double PhotonStartTime = RICH__photons_start_time;
		    if (iSTtimeCorr) PhotonStartTime = RICH__photons_start_time + REC__Event_startTime;
		    double CalcPhotonTime = PhotonStartTime + PhotonPathTime;
	 	
		    /* Delta T (measured-calculated) */
		    double DTime = MeasPhotonTime - CalcPhotonTime;
		
		
		    //printf("i=%d  rawT=%f  rawTc=%f  ts=%f\n", RICH__hits_rawtime, MeasPhotonTime, ts);
		    //printf("  stT=%f   anT=%f  DT=%f\n", RICH__photons_start_time, RICH__photons_analytic_time, DTime);
	


		
		    sprintf(name, "hDTime");
		    gDirectory->GetObject(name, h2);
		    h2->Fill(absChannel, DTime);
		  
		    sprintf(name, "hDTime2");
		    gDirectory->GetObject(name, h2);
		    h2->Fill(absChannel, DTime);

		    sprintf(name, "hDTimeVsDuration_Pmt%d", pmt);
		    gDirectory->GetObject(name, h2);
		    h2->Fill(duration, DTime);

		    /* Corrected delta T */
		    double DTimeCorr = DTime;
		    if (iTimeCorr) {
		      /* Corrected measured time */
		      double MeasPhotonTimeCorr = GetCorrectedTime(pmt, anode, MeasPhotonTime, duration);
		      DTimeCorr = MeasPhotonTimeCorr - CalcPhotonTime;
		    }

		    sprintf(name, "hDTimeCorr");
		    gDirectory->GetObject(name, h2);
		    h2->Fill(absChannel, DTimeCorr);

		    sprintf(name, "hDLenCorr");
		    gDirectory->GetObject(name, h2);
		    h2->Fill(absChannel, DTimeCorr*vLight);
		    
		    sprintf(name, "hDTimeCorrVsDuration_Pmt%d", pmt);
		    gDirectory->GetObject(name, h2);
		    h2->Fill(duration, DTimeCorr);
		
		    
		  }/* END of good photon type */


		}/* END loop over RICH photons */
		  
		if (nphotons) nRich++;
	      }/* Good RICH particle */
	      else {
		if (RichParticleQ == 0) {
		  nRichBadClusters++;
		}
	      }
	    

	    }/* END events with one track in the RICH */
	    

	  }


	}
      }/* END event with CLAS12 banks */

      entry++;
      if ( (entry%1000) == 0) {
	fprintf(stdout, "%d  %d  %d  %d  %d\r", entry, nTriggers, nRec, nRichTracks, nRich);
	fflush(stdout);	
      }


    }/* END loop over the current file */
	  
  }
  printf("Total number of events:  %d \n", entry);
  printf("Number of triggers found: %d\n", nTriggers);
  printf("Number of rec. events: %d\n", nRec);
  printf("Number of RICH clusters: %d\n", (nRichTracks+nRichBadClusters));
  printf("Number of RICH good tracks: %d\n", nRichTracks);
  printf("Number of RICH events: %d\n", nRich);

 
  /*******************************************************/
  /* Out histograms */
  TFile f(&rootFile[0], "recreate");
  //Hlist.ls();
  Hlist.Write();
  f.Close();
  
  return 1;
}
/* ------------------------------------------------------ */
void makeHistos()
{
  /* Histograms for the time calibration */

  char name[200];
  char title[200];


  /* channel binning */
  int nbins = nPMTS*nANODES;
  double bmin = 0.5;
  double bmax = nbins + 0.5;

  /* DeltaT binning */
  int na = 200;
  double amin = -120;
  double amax = -20;
  amin = HistoAveDT - 50;
  amax = HistoAveDT + 50;

  /* ====================================== */
  /* larger time window to accomodate all the possible timing variations */
  na = 1000;
  amin = -700.;
  amax = 300;
  /* ======== */
  

  /* DeltaTcorr binning */
  int ndt = 400;
  double dt1 = -50;
  double dt2 = 50;

  /* duration binning */
  int nc = 100;
  double cmin = -0.5;
  double cmax = cmin + nc;


  sprintf(title, "#Delta T vs channel");
  sprintf(name, "hDTime");
  TH2F *hDTime = new TH2F(name, title, nbins, bmin, bmax, na, amin, amax);
  Hlist.Add(hDTime);
  
  sprintf(title, "#Delta T vs channel");
  sprintf(name, "hDTime2");
  TH2F *hDTime2 = new TH2F(name, title, nbins, bmin, bmax, na, -1500, 1500);
  Hlist.Add(hDTime2);

  sprintf(title, "#Delta T_{corr} vs channel");
  sprintf(name, "hDTimeCorr");
  TH2F *hDTimeCorr = new TH2F(name, title, nbins, bmin, bmax, ndt, dt1, dt2);
  Hlist.Add(hDTimeCorr);

  sprintf(title, "#Delta L vs channel");
  sprintf(name, "hDLenCorr");
  TH2F *hDLenCorr = new TH2F(name, title, nbins, bmin, bmax, ndt, 30*dt1, 30*dt2);
  Hlist.Add(hDLenCorr);


  /* Delta T vs Duration per PMT */
  for (Int_t p=0; p<nPMTS; p++) {

      
    sprintf(name, "hDTimeVsDuration_Pmt%d", p+1);
    sprintf(title, "#Delta T vs Duration, PMT %d", p+1);
    TH2F *hDTimeVsDuration_Pmt = new TH2F(name, title, nc, cmin, cmax, na, amin, amax);
    Hlist.Add(hDTimeVsDuration_Pmt);

    sprintf(name, "hDTimeCorrVsDuration_Pmt%d", p+1);
    sprintf(title, "#Delta T_{corr} vs Duration, PMT %d", p+1);
    TH2F *hDTimeCorrVsDuration_Pmt = new TH2F(name, title, nc, cmin, cmax, ndt, dt1, dt2);
    Hlist.Add(hDTimeCorrVsDuration_Pmt);

  }
  
}
/* ===================================================== */
void PrintUsage(char *processName)
{
  fprintf(stderr,"Usage: %s [-Options] \n\n",processName);
  fprintf(stderr,"  Options:\n");

  fprintf(stderr,"\t-n[#]\t\tMax number of entries\n");
  fprintf(stderr,"\t-R[#]\t\tRun number\n");
  fprintf(stderr,"\t-T[#]\t\tApplying the time corrections: 1=Only offset; 2=Also time-walk\n");
  fprintf(stderr,"\t-s\t\tCorrecting by hand for the event start time (def: already done in the RICH rec.)\n");
  fprintf(stderr,"\t-r\t\tUsing the ray tracing solution (def.: using the analytic solution)\n");
  fprintf(stderr,"\t-P[#]\t\tUsing particle with PID=# (def.: 11->electrons, etc; -1->AllNeg; +1->AllPos; 0->All; 99->Straight tracks)\n");
  fprintf(stderr,"\t-t[#]\t\tSetting the average of the DTime window in the histograms (default: -90 ns)\n");
  fprintf(stderr,"\t-Z\t\tAnalysis of Zero Field data\n");

  exit(0);
}
