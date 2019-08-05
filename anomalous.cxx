/// Standard C++ Header Files
#include <iostream>
#include <cstring>

/// ROOT Header Files
#include <root.h>
#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <TVector3.h>

/// AMS Offline Software Header Files
#include "TrdKCluster.h"
#include "amschain.h"
#include "bcorr.h"
#include "TrCharge.h"
#include "TrSim.h"
#include "TrExtAlignDB.h"

bool debug = false;

AMSEventR*    pAMSEvent;
HeaderR*      pHeader;
DaqEventR*    pDaqEvent;
Level1R*      pLevel1;
ParticleR*    pParticle;
BetaHR*       pBetaH;
BetaR*        pBeta;
TrTrackR*     pTrTrack;
TrRecHitR*    pTrRecHit;
TrdTrackR*    pTrdTrack;
TrdSegmentR*  pTrdSegment;
EcalShowerR*  pEcalShower;
RichRingR*    pRichRing;
TrdKCluster*  pTrdKCluster;
TrdKHit*      pTrdKHit;
AntiClusterR* pAntiCluster;

TH1F* hEvtCounter;

char* outputFileName;
char* inputFileName;
const char* skirmishRunPath = "root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6/1319745123.00000001.root";
unsigned int nTreesInChain = 0;
unsigned int nCut;
unsigned int nEvents;
bool goodhw;
bool hwerror;
bool unbiased;
unsigned int nLevel1;
unsigned int nExaminedParticles;
unsigned int nSelected;
unsigned int nSelectedParticles;
unsigned int nRun;
unsigned long int nEvent;
unsigned int nLevel;
unsigned int nParticle;
unsigned int nCharge;
unsigned int nTrTrack;
unsigned int nTrdTrack;                  // Number of the TRD tracks which are successfully reconstructed
unsigned int nAntiCluster;               // Number of clusters on the ACC
unsigned int nTofClustersInTime;         // Number of in-time clusters on the TOF
unsigned int nRichRing;                  // Number of successfully reconstructed the RICH rings
unsigned int nRichRingB;                 // Number of successfully reconstructed the RICH rings with algorithm B
unsigned int nBeta;                      // Number of successfully estimated beta(v/c) values
unsigned int nBetaB;                     // Number of successfully estimated beta(v/c) values with algorithm B
unsigned int nBetaH;                     // Number of successfully estimated beta(v/c) values with algorithm H
unsigned int nEcalShower;                    // Number of the ECAL shower objects
unsigned int nVertex;                    // Number of vertices in this event

unsigned int nAntiMCCluster;             // Number of MC ACC clusters
unsigned int nTrMCCluster;               // Number of MC Tracker clusters
unsigned int nTofMCCluster;              // Number of MC TOF clusters
unsigned int nTofMCPmtHit;               // Number of MC TOF PMT hits
unsigned int nEcalMCHit;                 // Number of MC ECAL hits
unsigned int nTrdMCCluster;              // Number of MC TRD clusters
unsigned int nRichMCCluster;             // Number of MC RICH clusters
unsigned int nMCTrack;                   // Number of MC Tracks
unsigned int nMCEventg;                  // Number of MCEventg objects

float         liveTime;                   // Livetime fraction
float         utcTime;                    // UTC time
float         orbitAltitude;              // (cm) in GTOD coordinates system.
float         orbitLatitude;              // (rad) in GTOD coordinates system.
float         orbitLongitude;             // (rad) in GTOD coordinates system.
float         orbitLatitudeM;             // (rad) in eccentric dipole coordinate system.
float         orbitLongitudeM;            // (rad) in eccentric dipole coordinate system.
float         velR;                       // Speed of the ISS in radial direction
float         velTheta;                   // Angular speed of the ISS in polar angle direction
float         velPhi;                     // Angular speed of the ISS in azimuthal angle direction
float         yaw;                        // A parameter describing ISS attitude (tilted angle with respect to the x-axis)
float         pitch;                      // A parameter describing ISS attitude (tilted angle with respect to the y-axis)
float         roll;                       // A parameter describing ISS attitude (tilted angle with respect to the z-axis)
float         gLongitude;                 // Galactic longitude of the incoming particle.
float         gLatitude;                  // Galactic latitude of the incoming particle.
int           gCoordCalcResult;           // Return value for galactic coordinate calculation.
int           gCoordCalcResultBit;        // Result of GetGalCoo written in bit form
float         sunPosAzimuth;              // Azimuthal angle of the position of the Sun.
float         sunPosElevation;            // Elevation angle of the position of the Sun.
int           sunPosCalcResult;           // Return value for the Sun's position calculation.
unsigned int  unixTime;                   // UNIX time
float         acEventTime;
float         solarArrayCoord[3];
float         ptlCutOffStoermer;
float         ptlCutOffDipole;
float         ptlCutOffIGRF;
int           isInShadow;                 // Value for check whether the AMS is in ISS solar panel shadow or not.
float         zenithAngle;
int           isInSAA;
// TOF variables
int           tofNCluster;
int           tofNClusterH;
int           tofNUsedHits;
int           tofNUnusedHits;
int           tofNUsedLayersForQ;
float         tofBeta;
float         tofBetaH;
float         tofInvBetaErr;
float         tofInvBetaErrC;
float         tofInvBetaErrH;
float         tofNormEBetaV;
float         tofBetaSH;
float         tofBetaC;
float         tofBetaCH;
float         tofEBetaCV;
float         tofMass;
float         tofMassError;
int           isGoodBeta;
int           isTkTofMatch;
float         tofChisqT;
float         tofChisqC;
float         tofChisqTH;
float         tofChisqCH;
float         tofReducedChisqT;
float         tofReducedChisqC;
int           tofSumHit;
int           tofNUsedClusterH;
int           tofPatternOnLayer[4];
float         tofDepositedEnergyOnLayer[4];
float         tofEstimatedChargeOnLayer[4];
float         tofCharge;
float         tofChargeRMS;
float         tofUpperCharge;
float         tofLowerCharge;
float         tofChargeOnLayer[4];
float         tofTimeOnLayer[4];
float         tofETimeOnLayer[4];
float         tofETCooOnLayer[4];
float         tofTResidualOnLayer[4];
float         tofT0;
float         tofTkTFLenOnLayer[4];
float         tofCResidualXOnLayer[4];
float         tofCResidualYOnLayer[4];
int           nTofLUsedForZ;
float         probTOFZ;
int           tofZ;

/// Shower related variables
int           isEcalAvailable;
float         showerEnergyD;
float         showerEnergyDL[18];
float         showerEnergyE;
float         showerEnergyCorrected;
float         showerBDT;
float         showerCofG[3];
float         showerCofGDist;
float         showerCofGdX;
float         showerCofGdY;

// Track related variables
float         q_inn;
int           z_inn;
float         mass;
int           algo;
int           refit;
bool          xside[9];
bool          yside[9];
int           id_inner;
int           trkFitCodeL1Inner;
float         trkRigidityL1Inner;
float         trkRigidityInverseErrorL1Inner;
float         trkReducedChisquareL1InnerX;
float         trkReducedChisquareL1InnerY;
int           trkFitCodeL9Inner;
float         trkRigidityL9Inner;
float         trkRigidityInverseErrorL9Inner;
float         trkReducedChisquareL9InnerX;
float         trkReducedChisquareL9InnerY;
int           trkFitCodeUpperInner;
float         trkRigidityUpperInner;
float         trkRigidityInverseErrorUpperInner;
float         trkReducedChisquareUpperInnerX;
float         trkReducedChisquareUpperInnerY;
int           trkFitCodeLowerInner;
float         trkRigidityLowerInner;
float         trkRigidityInverseErrorLowerInner;
float         trkReducedChisquareLowerInnerX;
float         trkReducedChisquareLowerInnerY;
int           trkFitCodeFS;
float         trkRigidityFS;
float         trkRigidityInverseErrorFS;
float         trkReducedChisquareFSX;
float         trkReducedChisquareFSY;
int           trkFitCodeMS;
float         trkRigidityMS;
float         trkRigidityInverseErrorMS;
float         trkReducedChisquareMSX;
float         trkReducedChisquareMSY;
int           trkFitCodeInner;
float         trkRigidityInner;
float         trkRigidityInverseErrorInner;
float         trkReducedChisquareInnerX;
float         trkReducedChisquareInnerY;
float         trkLayerJQ[9];
int           trkEdepLayerJXSideOK[9];
int           trkEdepLayerJYSideOK[9];
float         trkEdepLayerJX[9];
float         trkEdepLayerJY[9];
float         trkPosXLJ[9];
float         trkPosYLJ[9];
float         trkPosZLJ[9];
float         trkHitCooXLJ[9];
float         trkHitCooYLJ[9];
float         trkDirThetaLJ[9];
float         trkDirPhiLJ[9];
float         trkCharge;
float         trkChargeRMS;
float         trkInnerCharge;
float         trkInnerChargeRMS;
int           trkZ;
int           trkInnerZ;
int           trkLayerJZ[9];
int           trkHasExtLayers;
int           isRichAvailable;
int           richRebuild;
int           richIsGood;
int           richIsClean;
int           richIsNaF;
int           richTileIndex;
float         richDistanceTileBorder;
int           richChargeCorrections;
int           richPMTCorrectionsFailed;
int           richUsedHits;
int           richUsedTileIndex;
float         richRingWidth;
float         richWidth;
int           richNHits;
int           richNPMTsOnRing;
float         richBeta;
float         richBetaError;
float         richChargeSquared;
float         richKolmogorovProbability;
float         richPhotoelectrons;                         // Return value from RichRingR::getPhotoelectrons(bool corr=true)
float         richPhotoElectrons;                         // Return value from RichRingR::getPhotoElectrons(int pmt, bool corr)
float         richExpectedPhotoelectrons;                 // Return value from RichRingR::getExpectedPhotoElectrons(bool corr=true)
float         richExpectedPhotoElectrons;                 // Return value from RichRingR::getExpectedPhotoElectrons(int pmt, bool corr)
int           richNUsedHits;
float         richTrackEmissionPoint[5];
float         richTheta;
float         richPhi;
float         richBetaConsistency;
int           richReflectedHits;                          // Return value from RichRingR::getReflectedHits()
float         richPMTChargeConsistency;                   // Return value from RichRingR::getPMTChargeConsistency()
float         richTrackerTrackOnRadiatorX;
float         richTrackerTrackOnRadiatorY;
int           trdNClusters;
int           trdNUnusedHits;
int           trdNTracks;
float         trdTrackPhi;
float         trdTrackTheta;
float         trdTrackChi2;
int           trdTrackPattern;
float         trdTrackCharge;
float         trdTrackEdepL[20];
float         trdTrackMeanDepositedEnergy;
float         trdTrackTotalDepositedEnergy;
float         trdTrackDeviationXWithInnerTrk;
float         trdTrackDeviationYWithInnerTrk;
float         trdTrackDistanceBetweenInnerTrk;
float         ResidualXBetweenInnerTrackAndSplineTrack;
float         ResidualYBetweenInnerTrackAndSplineTrack;
int           trdNVertex;
int           trdKNRawHits;
int           trdKIsReadAlignmentOK;
int           trdKIsReadCalibOK;
int           trdKNHits;
int           trdKIsValid;
float         trdKElectronToProtonLogLikelihoodRatio;
float         trdKHeliumToElectronLogLikelihoodRatio;
float         trdKHeliumToProtonLogLikelihoodRatio;
float         trdKCharge;
float         trdKChargeError;
int           trdKNUsedHitsForCharge;
float         trdKAmpLayer[20];
float         trdKTotalPathLength;
float         trdKTotalAmp;
float         trdKElectronLikelihood;
float         trdKProtonLikelihood;
float         trdKHeliumLikelihood;
float         trdPElectronToProtonLogLikelihoodRatio;
float         trdPHeliumToElectronLogLikelihoodRatio;
float         trdPHeliumToProtonLogLikelihoodRatio;
float         trdPElectronLikelihood;
float         trdPProtonLikelihood;
float         trdPHeliumLikelihood;
int           trdPNHybridHits;
int           trdPNActiveLayers;
int           trdPNLowAdcHits;
int           trdPNLowDxHits;
int           accNHits;
int           accNRecoClPG;

//char skirmishRunPath[] = "root://eosams.cern.ch//eos/ams/Data/AMS02/2014/ISS.B950/pass6/1366261348.00000001.root";  // Bad run example


//==============================================================================//
// ---------------------- Physical Variables Declaration -----------------------//
//==============================================================================//

const std::string CurrentDateTime()
{
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

  return buf;
}

int main(int argc, char* argv[])
{
  std::cout << "Initiating He3/He4 Data Reduction Software ................" << std::endl;
  std::cout << "RUN MODE: ISS" << std::endl;

  TStopwatch clock;
  clock.Start();
  std::cout << "TStopwatch is started." << std::endl;
  std::cout << "Current time: " << CurrentDateTime() << std::endl;

  std::cout << "Initializing AMSChain instance ..................." << std::endl;
  AMSChain chain;

  switch( argc )
  {
    case 1:
      std::cout << "RUN MODE: Single Test Run (Cat. 1)" << std::endl;
      outputFileName = new char[ strlen( "testrun.root" ) + 1 ];
      strncpy( outputFileName, "testrun.root", strlen( "testrun.root" ) + 1 );

      if( chain.Add( skirmishRunPath ) == 1 )
        std::cout << "File [ " << skirmishRunPath << " ] opened successfully." << std::endl;
      else
      {
        std::cout << " ERROR: File open error, [ " << skirmishRunPath << " ] can not be found!" << std::endl;
        return 1;
      }
      break;
    case 2:
      std::cout << "RUN MODE: Single Test Run (Cat. 2)" << std::endl;
      outputFileName = new char[ strlen( "testrun.root" ) + 1 ];
      strncpy( outputFileName, "testrun.root", strlen( "testrun.root" ) + 1 );

      if( chain.Add( argv[2] ) == 1 )
        std::cout << "File [ " << skirmishRunPath << " ] opened successfully." << std::endl;
      else
      {
        std::cout << " ERROR: File open error, [ " << skirmishRunPath << " ] can not be found!" << std::endl;
        return 1;
      }
      break;
    case 3:
      std::cout << "RUN MODE: Batch-job Mode" << std::endl;
      outputFileName = new char[ strlen(argv[2]) + 1 ];
      strncpy( outputFileName, argv[2], strlen( argv[2] ) + 1 );

      char listFileName[256];
      char inputFileName[256];

      FILE* fp_listfile;
      strcpy( listFileName, argv[1] );

      if( ( fp_listfile = fopen( listFileName, "r" ) ) == NULL )
      {
        std::cerr << " ERROR: Failed to open file [ " << listFileName << " ]!" << std::endl;
        return 2;
      }

      char* line_p;
      while( fgets( inputFileName, 256, fp_listfile ) != NULL )
      {
        if( ( line_p = strchr( inputFileName, '\n') ) != NULL )
          *line_p = 0;

        if( chain.Add( inputFileName ) == 1 )
        {
          std::cout << "The file [ " << inputFileName << " ] is added to the chain." << std::endl;
          std::cout << "Currently loaded events: " << chain.GetEntries() << std::endl;
          nTreesInChain++;
        }
        else
        {
          std::cout << " ERROR: Failed to open file [ " << inputFileName << " ]!" << std::endl;
          return 1;
        }
      }
      fclose(fp_listfile);
      break;
    case 4:
      std::cout << "RUN MODE: Single Test Run (Cat. 3)" << std::endl;
      outputFileName = new char[ strlen( argv[2] ) + 1 ];
      strncpy( outputFileName, argv[2], strlen( argv[2] ) + 1 );
      strncpy( inputFileName, argv[1], strlen( argv[1] ) + 1 );
      if( chain.Add(inputFileName) == 1 )
      {
        std::cout << "The file [ " << inputFileName << " ] is added to the chain." << std::endl;
        std::cout << argv[3] << " events will be processed for test." << std::endl;
      }
      else
      {
        std::cerr << " ERROR: File open error, [ " << inputFileName << " ] can not be found!" << std::endl;
        return 1;
      }
      break;
    default:
      break;
  }

  TFile* pOutputTFile = new TFile( outputFileName, "RECREATE" );
  TTree* outTree = new TTree("KNUTree", "KNUTree");

  //===========================================================================
  //
  //
  // Initialize
  //
  //
  //===========================================================================

  std::cout << "Initializing variables ... " << std::endl;

  nEvents                           = chain.GetEntries();
  nCut                              = 0;
  nExaminedParticles                = 0;
  nSelected                         = 0;
  nSelectedParticles                = 0;
  nRun                              = 0;
  nEvent                            = 0;
  nLevel1                           = 0;
  nParticle                         = 0;
  nCharge                           = 0;
  nTrTrack                          = 0;
  nTrdTrack                         = 0;
  nAntiCluster                      = 0;
  nTofClustersInTime                = 0;
  nRichRing                         = 0;
  nRichRingB                        = 0;
  nBeta                             = 0;
  nBetaB                            = 0;
  nBetaH                            = 0;
  nEcalShower                       = 0;
  nVertex                           = 0;
  nAntiMCCluster                    = 0;
  nTrMCCluster                      = 0;
  nTofMCCluster                     = 0;
  nTofMCPmtHit                      = 0;
  nEcalMCHit                        = 0;
  nTrdMCCluster                     = 0;
  nRichMCCluster                    = 0;
  nMCTrack                          = 0;
  nMCEventg                         = 0;
  liveTime                          = 0;
  utcTime                           = 0;
  orbitAltitude                     = 0;
  orbitLatitude                     = 0;
  orbitLongitude                    = 0;
  orbitLatitudeM                    = 0;
  orbitLongitudeM                   = 0;
  velR                              = 0;
  velTheta                          = 0;
  velPhi                            = 0;
  yaw                               = 0;
  pitch                             = 0;
  roll                              = 0;
  gLongitude                        = 0;
  gLatitude                         = 0;
  gCoordCalcResult                  = 0;
  gCoordCalcResultBit               = 0;
  sunPosAzimuth                     = 0;
  sunPosElevation                   = 0;
  sunPosCalcResult                  = 0;
  unixTime                          = 0;
  acEventTime                       = 0;
  solarArrayCoord[0]                = 0;
  solarArrayCoord[1]                = 0;
  solarArrayCoord[2]                = 0;
  isInShadow                        = 0;
  zenithAngle                       = 0;
  isInSAA                           = 0;
  tofNCluster                       = 0;
  tofNClusterH                      = 0;
  tofNUsedHits                      = 0;
  tofNUnusedHits                    = 0;
  tofNUsedLayersForQ                = 0;
  tofBeta                           = 0;
  tofBetaC                          = 0;
  tofBetaH                          = 0;
  tofBetaCH                         = 0;
  tofBetaSH                         = 0;
  tofInvBetaErr                     = 0;
  tofInvBetaErrC                    = 0;
  tofInvBetaErrH                    = 0;
  tofMass                           = 0;
  tofMassError                      = 0;
  isGoodBeta                        = 0;
  isTkTofMatch                      = 0;
  tofReducedChisqC                  = 0;
  tofReducedChisqT                  = 0;
  tofDepositedEnergyOnLayer[0]      = 0;
  tofDepositedEnergyOnLayer[1]      = 0;
  tofDepositedEnergyOnLayer[2]      = 0;
  tofDepositedEnergyOnLayer[3]      = 0;
  tofEstimatedChargeOnLayer[0]      = 0;
  tofEstimatedChargeOnLayer[1]      = 0;
  tofEstimatedChargeOnLayer[2]      = 0;
  tofEstimatedChargeOnLayer[3]      = 0;
  tofCharge                         = -1;
  tofChargeRMS                      = 0;
  tofUpperCharge                    = -1;
  tofLowerCharge                    = -1;
  tofChargeOnLayer[0]               = -1;
  tofChargeOnLayer[1]               = -1;
  tofChargeOnLayer[2]               = -1;
  tofChargeOnLayer[3]               = -1;
  tofZ                              = -1;
  isEcalAvailable                   = 0;
  showerEnergyD                     = 0;
  showerEnergyDL[0]                 = 0;
  showerEnergyDL[1]                 = 0;
  showerEnergyDL[2]                 = 0;
  showerEnergyDL[3]                 = 0;
  showerEnergyDL[4]                 = 0;
  showerEnergyDL[5]                 = 0;
  showerEnergyDL[6]                 = 0;
  showerEnergyDL[7]                 = 0;
  showerEnergyDL[8]                 = 0;
  showerEnergyDL[9]                 = 0;
  showerEnergyDL[10]                = 0;
  showerEnergyDL[11]                = 0;
  showerEnergyDL[12]                = 0;
  showerEnergyDL[13]                = 0;
  showerEnergyDL[14]                = 0;
  showerEnergyDL[15]                = 0;
  showerEnergyDL[16]                = 0;
  showerEnergyDL[17]                = 0;
  showerEnergyE                     = 0;
  showerEnergyCorrected             = 0;
  showerBDT                         = 0;
  showerCofG[0]                     = 0;
  showerCofG[1]                     = 0;
  showerCofG[2]                     = 0;
  showerCofGDist                    = 0;
  showerCofGdX                      = 0;
  showerCofGdY                      = 0;
  trkFitCodeL1Inner                 = -1;
  trkRigidityL1Inner                = -9999;
  trkRigidityInverseErrorL1Inner    = -9999;
  trkReducedChisquareL1InnerX       = -9999;
  trkReducedChisquareL1InnerY       = -9999;
  trkFitCodeFS                      = -1;
  trkRigidityFS                     = -9999;
  trkRigidityInverseErrorFS         = -9999;
  trkReducedChisquareFSX            = -9999;
  trkReducedChisquareFSY            = -9999;
  trkFitCodeMS                      = -1;
  trkRigidityMS                     = -9999;
  trkRigidityInverseErrorMS         = -9999;
  trkReducedChisquareMSY            = -9999;
  trkReducedChisquareMSX            = -9999;
  trkFitCodeInner                   = -1;
  trkRigidityInverseErrorInner      = -9999;
  trkReducedChisquareInnerX         = -9999;
  trkReducedChisquareInnerY         = -9999;
  trkFitCodeUpperInner              = -1;
  trkRigidityInverseErrorUpperInner = -9999;
  trkReducedChisquareUpperInnerX    = -9999;
  trkReducedChisquareUpperInnerY    = -9999;
  trkFitCodeLowerInner              = -1;
  trkRigidityInverseErrorLowerInner = -9999;
  trkReducedChisquareLowerInnerX    = -9999;
  trkReducedChisquareLowerInnerY    = -9999;
  trkLayerJQ[0]                     = -1;
  trkLayerJQ[1]                     = -1;
  trkLayerJQ[2]                     = -1;
  trkLayerJQ[3]                     = -1;
  trkLayerJQ[4]                     = -1;
  trkLayerJQ[5]                     = -1;
  trkLayerJQ[6]                     = -1;
  trkLayerJQ[7]                     = -1;
  trkLayerJQ[8]                     = -1;
  trkEdepLayerJXSideOK[0]           = 0;
  trkEdepLayerJXSideOK[1]           = 0;
  trkEdepLayerJXSideOK[2]           = 0;
  trkEdepLayerJXSideOK[3]           = 0;
  trkEdepLayerJXSideOK[4]           = 0;
  trkEdepLayerJXSideOK[5]           = 0;
  trkEdepLayerJXSideOK[6]           = 0;
  trkEdepLayerJXSideOK[7]           = 0;
  trkEdepLayerJXSideOK[8]           = 0;
  trkEdepLayerJYSideOK[0]           = 0;
  trkEdepLayerJYSideOK[1]           = 0;
  trkEdepLayerJYSideOK[2]           = 0;
  trkEdepLayerJYSideOK[3]           = 0;
  trkEdepLayerJYSideOK[4]           = 0;
  trkEdepLayerJYSideOK[5]           = 0;
  trkEdepLayerJYSideOK[6]           = 0;
  trkEdepLayerJYSideOK[7]           = 0;
  trkEdepLayerJYSideOK[8]           = 0;
  trkEdepLayerJX[0]                  = 0;
  trkEdepLayerJX[1]                  = 0;
  trkEdepLayerJX[2]                  = 0;
  trkEdepLayerJX[3]                  = 0;
  trkEdepLayerJX[4]                  = 0;
  trkEdepLayerJX[5]                  = 0;
  trkEdepLayerJX[6]                  = 0;
  trkEdepLayerJX[7]                  = 0;
  trkEdepLayerJX[8]                  = 0;
  trkEdepLayerJY[0]                  = 0;
  trkEdepLayerJY[1]                  = 0;
  trkEdepLayerJY[2]                  = 0;
  trkEdepLayerJY[3]                  = 0;
  trkEdepLayerJY[4]                  = 0;
  trkEdepLayerJY[5]                  = 0;
  trkEdepLayerJY[6]                  = 0;
  trkEdepLayerJY[7]                  = 0;
  trkEdepLayerJY[8]                  = 0;
  for( int i = 0; i < 9; i++ )
  {
    trkPosXLJ[i]                    = 0;
    trkPosYLJ[i]                    = 0;
    trkPosZLJ[i]                    = 0;
    trkHitCooXLJ[i]                 = 0;
    trkHitCooYLJ[i]                 = 0;
    trkDirThetaLJ[i]                = 0;
    trkDirPhiLJ[i]                  = 0;
  }
  trkCharge                         = -1;
  trkInnerCharge                    = -1;
  trkZ                              = -1;
  trkInnerZ                         = -1;
  for( int i = 0; i < 9; i++ ) trkLayerJZ[i] = -1;
  trkHasExtLayers                   = 0;
  isRichAvailable                   = 0;
  richRebuild                       = 0;
  richIsGood                        = 0;
  richIsClean                       = 0;
  richIsNaF                         = 0;
  richTileIndex                     = -1;
  richDistanceTileBorder            = -1;
  richChargeCorrections             = 0;
  richPMTCorrectionsFailed          = -10;
  richWidth                         = -1;
  richBetaConsistency               = -999;
  richReflectedHits                 = -1;
  richPMTChargeConsistency          = -999;
  richUsedHits                      = 0;
  richRingWidth                     = 0;
  richNHits                         = 0;
  richNPMTsOnRing                   = 0;
  richBeta                          = 0;
  richBetaError                     = 0;
  richChargeSquared                 = 0;
  richKolmogorovProbability         = 0;
  richPhotoelectrons                = 0;
  richPhotoElectrons                = 0;
  richExpectedPhotoelectrons        = 0;
  richExpectedPhotoElectrons        = 0;
  richNUsedHits                     = 0;
  richTheta                         = 0;
  richPhi                           = 0;
  richTrackerTrackOnRadiatorX       = -999;
  richTrackerTrackOnRadiatorY       = -999;
  trdNClusters                      = 0;
  trdNUnusedHits                    = 0;
  trdNTracks                        = 0;
  trdNVertex                        = 0;
  trdTrackPhi                       = 0;
  trdTrackTheta                     = 0;
  trdTrackChi2                      = 0;
  trdTrackPattern                   = 0;
  trdTrackCharge                    = 0;
  trdTrackEdepL[0]                  = 0;
  trdTrackEdepL[1]                  = 0;
  trdTrackEdepL[2]                  = 0;
  trdTrackEdepL[3]                  = 0;
  trdTrackEdepL[4]                  = 0;
  trdTrackEdepL[5]                  = 0;
  trdTrackEdepL[6]                  = 0;
  trdTrackEdepL[7]                  = 0;
  trdTrackEdepL[8]                  = 0;
  trdTrackEdepL[9]                  = 0;
  trdTrackEdepL[10]                 = 0;
  trdTrackEdepL[11]                 = 0;
  trdTrackEdepL[12]                 = 0;
  trdTrackEdepL[13]                 = 0;
  trdTrackEdepL[14]                 = 0;
  trdTrackEdepL[15]                 = 0;
  trdTrackEdepL[16]                 = 0;
  trdTrackEdepL[17]                 = 0;
  trdTrackEdepL[18]                 = 0;
  trdTrackEdepL[19]                 = 0;
  trdTrackMeanDepositedEnergy       = 0;
  trdTrackTotalDepositedEnergy      = 0;
  trdTrackDeviationXWithInnerTrk    = 0;
  trdTrackDeviationYWithInnerTrk    = 0;
  trdTrackDistanceBetweenInnerTrk   = 0;
  ResidualXBetweenInnerTrackAndSplineTrack = -9999;
  ResidualYBetweenInnerTrackAndSplineTrack = -9999;
  trdKNRawHits                             = 0;
  trdKIsReadAlignmentOK                    = 0;
  trdKIsReadCalibOK                        = 0;
  trdKNHits                                = 0;
  trdKIsValid                              = 0;
  trdKElectronToProtonLogLikelihoodRatio   = 0;
  trdKHeliumToElectronLogLikelihoodRatio   = 0;
  trdKHeliumToProtonLogLikelihoodRatio     = 0;
  trdKCharge                               = -1;
  trdKChargeError                          = 0;
  trdKNUsedHitsForCharge                   = 0;
  trdKAmpLayer[0]                          = 0;
  trdKAmpLayer[1]                          = 0;
  trdKAmpLayer[2]                          = 0;
  trdKAmpLayer[3]                          = 0;
  trdKAmpLayer[4]                          = 0;
  trdKAmpLayer[5]                          = 0;
  trdKAmpLayer[6]                          = 0;
  trdKAmpLayer[7]                          = 0;
  trdKAmpLayer[8]                          = 0;
  trdKAmpLayer[9]                          = 0;
  trdKAmpLayer[10]                         = 0;
  trdKAmpLayer[11]                         = 0;
  trdKAmpLayer[12]                         = 0;
  trdKAmpLayer[13]                         = 0;
  trdKAmpLayer[14]                         = 0;
  trdKAmpLayer[15]                         = 0;
  trdKAmpLayer[16]                         = 0;
  trdKAmpLayer[17]                         = 0;
  trdKAmpLayer[18]                         = 0;
  trdKAmpLayer[19]                         = 0;
  trdKTotalPathLength                      = 0;
  trdKTotalAmp                             = 0;
  trdKElectronLikelihood                   = 0;
  trdKHeliumLikelihood                     = 0;
  accNHits                                     = 0;
  accNRecoClPG                                 = 0;

  outTree->Branch("nRun",                                    &nRun,                                    "nRun/i");
  outTree->Branch("nEvent",                                  &nEvent,                                  "nEvent/i");
  outTree->Branch("nLevel1",                                 &nLevel1,                                 "nLevel1/i");
  outTree->Branch("nParticle",                               &nParticle,                               "nParticle/i");
  outTree->Branch("nCharge",                                 &nCharge,                                 "nCharge/i");
  outTree->Branch("nTrTrack",                                &nTrTrack,                                "nTrTrack/i");
  outTree->Branch("nTrdTrack",                               &nTrdTrack,                               "nTrdTrack/i");
  outTree->Branch("nAntiCluster",                            &nAntiCluster,                            "nAntiCluster/i");
  outTree->Branch("nTofClustersInTime",                      &nTofClustersInTime,                      "nTofClustersInTime/I");
  outTree->Branch("nRichRing",                               &nRichRing,                               "nRichRing/i");
  outTree->Branch("nRichRingB",                              &nRichRingB,                              "nRichRingB/i");
  outTree->Branch("nBeta",                                   &nBeta,                                   "nBeta/i");
  outTree->Branch("nBetaB",                                  &nBetaB,                                  "nBetaB/i");
  outTree->Branch("nBetaH",                                  &nBetaH,                                  "nBetaH/i");
  outTree->Branch("nShower",                                 &nEcalShower,                             "nShower/i");
  outTree->Branch("nVertex",                                 &nVertex,                                 "nVertex/i");
  outTree->Branch("liveTime",                                &liveTime,                                "liveTime/F");
  outTree->Branch("utcTime",                                 &utcTime,                                 "utcTime/F");
  outTree->Branch("orbitAltitude",                           &orbitAltitude,                           "orbitAltitude/F");
  outTree->Branch("orbitLatitude",                           &orbitLatitude,                           "orbitLatitude/F");
  outTree->Branch("orbitLongitude",                          &orbitLongitude,                          "orbitLongitude/F");
  outTree->Branch("orbitLatitudeM",                          &orbitLatitudeM,                          "orbitLatitudeM/F");
  outTree->Branch("orbitLongitudeM",                         &orbitLongitudeM,                         "orbitLongitudeM/F");
  outTree->Branch("velR",                                    &velR,                                    "velR/F");
  outTree->Branch("velTheta",                                &velTheta,                                "velTheta/F");
  outTree->Branch("velPhi",                                  &velPhi,                                  "velPhi/F");
  outTree->Branch("yaw",                                     &yaw,                                     "yaw/F");
  outTree->Branch("pitch",                                   &pitch,                                   "pitch/F");
  outTree->Branch("roll",                                    &roll,                                    "roll/F");
  outTree->Branch("gLongitude",                              &gLongitude,                              "gLongitude/F");
  outTree->Branch("gLatitude",                               &gLatitude,                               "gLatitude/F");
  outTree->Branch("gCoordCalcResult",                        &gCoordCalcResult,                        "gCoordCalcResult/I");
  outTree->Branch("gCoordCalcResultBit",                     &gCoordCalcResultBit,                     "gCoordCalcResultBit/I");
  outTree->Branch("sunPosAzimuth",                           &sunPosAzimuth,                           "sunPosAzimuth/F");
  outTree->Branch("sunPosElevation",                         &sunPosElevation,                         "sunPosElevation/F");
  outTree->Branch("sunPosCalcResult",                        &sunPosCalcResult,                        "sunPosCalcResult/I");
  outTree->Branch("unixTime",                                &unixTime,                                "unixTime/i");
  outTree->Branch("acEventTime",                             &acEventTime,                             "acEventTime/F");
  outTree->Branch("solarArrayCoord",                         &solarArrayCoord,                         "solarArrayCoord[3]/F");
  outTree->Branch("ptlCutOffStoermer",                       &ptlCutOffStoermer,                       "ptlCutOffStoermer/F");
  outTree->Branch("ptlCutOffDipole",                         &ptlCutOffDipole,                         "ptlCutOffDipole/F");
  outTree->Branch("ptlCutOffIGRF",                           &ptlCutOffIGRF,                           "ptlCutOffIGRF/F");
  outTree->Branch("isInShadow",                              &isInShadow,                              "isInShadow/I");
  outTree->Branch("zenithAngle",                             &zenithAngle,                             "zenithAngle/F");
  outTree->Branch("isInSAA",                                 &isInSAA,                                 "isInSAA/I");
  outTree->Branch("tofNCluster",                             &tofNCluster,                             "tofNCluster/I");
  outTree->Branch("tofNClusterH",                            &tofNClusterH,                            "tofNClusterH/I");
  outTree->Branch("tofNUsedHits",                            &tofNUsedHits,                            "tofNUsedHits/I");
  outTree->Branch("tofNUnusedHits",                          &tofNUnusedHits,                          "tofNUnusedHits/I");
  outTree->Branch("tofNUsedLayersForQ",                      &tofNUsedLayersForQ,                      "tofNUsedLayersForQ/I");
  outTree->Branch("tofBeta",                                 &tofBeta,                                 "tofBeta/F");
  outTree->Branch("tofBetaH",                                &tofBetaH,                                "tofBetaH/F");
  outTree->Branch("tofInvBetaErr",                           &tofInvBetaErr,                           "tofInvBetaErr/F");
  outTree->Branch("tofInvBetaErrC",                          &tofInvBetaErrC,                          "tofInvBetaErrC/F");
  outTree->Branch("tofInvBetaErrH",                          &tofInvBetaErrH,                          "tofInvBetaErrH/F");
  outTree->Branch("tofNormEBetaV",                           &tofNormEBetaV,                           "tofNormEBetaV/F");
  outTree->Branch("tofBetaSH",                               &tofBetaSH,                               "tofBetaSH/F");
  outTree->Branch("tofBetaC",                                &tofBetaC,                                "tofBetaC/F");
  outTree->Branch("tofBetaCH",                               &tofBetaCH,                               "tofBetaCH/F");
  outTree->Branch("tofEBetaCV",                              &tofEBetaCV,                              "tofEBetaCV/F");
  outTree->Branch("tofMass",                                 &tofMass,                                 "tofMass/F");
  outTree->Branch("tofMassError",                            &tofMassError,                            "tofMassError/F");
  outTree->Branch("isGoodBeta",                              &isGoodBeta,                              "isGoodBeta/I");
  outTree->Branch("isTkTofMatch",                            &isTkTofMatch,                            "isTkTofMatch/I");
  outTree->Branch("tofChisqT",                               &tofChisqT,                               "tofChisqT/F");
  outTree->Branch("tofChisqC",                               &tofChisqC,                               "tofChisqC/F");
  outTree->Branch("tofReducedChisqT",                        &tofReducedChisqT,                        "tofReducedChisqT/F");
  outTree->Branch("tofReducedChisqC",                        &tofReducedChisqC,                        "tofReducedChisqC/F");
  outTree->Branch("tofSumHit",                               &tofSumHit,                               "tofSumHit/I");
  outTree->Branch("tofNUsedClusterH",                        &tofNUsedClusterH,                        "tofNUsedClusterH/I");
  outTree->Branch("tofPatternOnLayer",                       &tofPatternOnLayer,                       "tofPatternOnLayer[4]/I");
  outTree->Branch("tofDepositedEnergyOnLayer",               &tofDepositedEnergyOnLayer,               "tofDepositedEnergyOnLayer[4]/F");
  outTree->Branch("tofEstimatedChargeOnLayer",               &tofEstimatedChargeOnLayer,               "tofEstimatedChargeOnLayer[4]/F");
  outTree->Branch("tofCharge",                               &tofCharge,                               "tofCharge/F");
  outTree->Branch("tofChargeRMS",                            &tofChargeRMS,                            "tofChargeRMS/F");
  outTree->Branch("tofChargeOnLayer",                        &tofChargeOnLayer,                        "tofChargeOnLayer[4]/F");
  outTree->Branch("tofTimeOnLayer",                          &tofTimeOnLayer,                          "tofTimeOnLayer[4]/F");
  outTree->Branch("tofETimeOnLayer",                         &tofETimeOnLayer,                         "tofETimeOnLayer[4]/F");
  outTree->Branch("tofETCooOnLayer",                         &tofETCooOnLayer,                         "tofETCooOnLayer[4]/F");
  outTree->Branch("tofTResidualOnLayer",                     &tofTResidualOnLayer,                     "tofTResidualOnLayer[4]/F");
  outTree->Branch("tofT0",                                   &tofT0,                                   "tofT0/F");
  outTree->Branch("tofTkTFLenOnLayer",                       &tofTkTFLenOnLayer,                       "tofTkTFLenOnLayer[4]/F");
  outTree->Branch("tofCResidualXOnLayer",                    &tofCResidualXOnLayer,                    "tofCResidualXOnLayer[4]/F");
  outTree->Branch("tofCResidualYOnLayer",                    &tofCResidualYOnLayer,                    "tofCResidualYOnLayer[4]/F");
  outTree->Branch("nTofLUsedForZ",                           &nTofLUsedForZ,                           "nTofLUsedForZ/I");
  outTree->Branch("probTOFZ",                                &probTOFZ,                                "probTOFZ/F");
  outTree->Branch("tofZ",                                    &tofZ,                                    "tofZ/I");
  outTree->Branch("isEcalAvailable",                         &isEcalAvailable,                         "isEcalAvailable/I");
  outTree->Branch("showerEnergyD",                           &showerEnergyD,                           "showerEnergyD/F");
  outTree->Branch("showerEnergyDL",                          &showerEnergyDL,                          "showerEnergyDL[18]/F");
  outTree->Branch("showerEnergyE",                           &showerEnergyE,                           "showerEnergyE/F");
  outTree->Branch("showerEnergyCorrected",                   &showerEnergyCorrected,                   "showerEnergyCorrected/F");
  outTree->Branch("showerBDT",                               &showerBDT,                               "showerBDT/F");
  outTree->Branch("showerCofG",                              &showerCofG,                              "showerCofG[3]/F");
  outTree->Branch("showerCofGDist",                          &showerCofGDist,                          "showerCofGDist/F");
  outTree->Branch("showerCofGdX",                            &showerCofGdX,                            "showerCofGdX/F");
  outTree->Branch("showerCofGdY",                            &showerCofGdY,                            "showerCofGdY/F");
  outTree->Branch("trkFitCodeFS",                            &trkFitCodeFS,                            "trkFitCodeFS/I");
  outTree->Branch("trkRigidityFS",                           &trkRigidityFS,                           "trkRigidityFS/F");
  outTree->Branch("trkRigidityInverseErrorFS",               &trkRigidityInverseErrorFS,               "trkRigidityInverseErrorFS/F");
  outTree->Branch("trkReducedChisquareFSX",                  &trkReducedChisquareFSX,                  "trkReducedChisquareFSX/F");
  outTree->Branch("trkReducedChisquareFSY",                  &trkReducedChisquareFSY,                  "trkReducedChisquareFSY/F");

  outTree->Branch("trkFitCodeL1Inner",                       &trkFitCodeL1Inner,                       "trkFitCodeL1Inner/I");
  outTree->Branch("trkRigidityL1Inner",                      &trkRigidityL1Inner,                      "trkRigidityL1Inner/F");
  outTree->Branch("trkRigidityInverseErrorL1Inner",          &trkRigidityInverseErrorL1Inner,          "trkRigidityInverseErrorL1Inner/F");
  outTree->Branch("trkReducedChisquareL1InnerX",             &trkReducedChisquareL1InnerX,             "trkReducedChisquareL1InnerX/F");
  outTree->Branch("trkReducedChisquareL1InnerY",             &trkReducedChisquareL1InnerY,             "trkReducedChisquareL1InnerY/F");

  outTree->Branch("trkFitCodeL9Inner",                       &trkFitCodeL9Inner,                       "trkFitCodeL9Inner/I");
  outTree->Branch("trkRigidityL9Inner",                      &trkRigidityL9Inner,                      "trkRigidityL9Inner/F");
  outTree->Branch("trkRigidityInverseErrorL9Inner",          &trkRigidityInverseErrorL9Inner,          "trkRigidityInverseErrorL9Inner/F");
  outTree->Branch("trkReducedChisquareL9InnerX",             &trkReducedChisquareL9InnerX,             "trkReducedChisquareL9InnerX/F");
  outTree->Branch("trkReducedChisquareL9InnerY",             &trkReducedChisquareL9InnerY,             "trkReducedChisquareL9InnerY/F");

  outTree->Branch("trkFitCodeUpperInner",                    &trkFitCodeUpperInner,                    "trkFitCodeUpperInner/I");
  outTree->Branch("trkRigidityUpperInner",                   &trkRigidityUpperInner,                   "trkRigidityUpperInner/F");
  outTree->Branch("trkRigidityInverseErrorUpperInner",       &trkRigidityInverseErrorUpperInner,       "trkRigidityInverseErrorUpperInner/F");
  outTree->Branch("trkReducedChisquareUpperInnerX",          &trkReducedChisquareUpperInnerX,          "trkReducedChisquareUpperInnerX/F");
  outTree->Branch("trkReducedChisquareUpperInnerY",          &trkReducedChisquareUpperInnerY,          "trkReducedChisquareUpperInnerY/F");

  outTree->Branch("trkFitCodeLowerInner",                    &trkFitCodeLowerInner,                    "trkFitCodeLowerInner/I");
  outTree->Branch("trkRigidityLowerInner",                   &trkRigidityLowerInner,                   "trkRigidityLowerInner/F");
  outTree->Branch("trkRigidityInverseErrorLowerInner",       &trkRigidityInverseErrorLowerInner,       "trkRigidityInverseErrorLowerInner/F");
  outTree->Branch("trkReducedChisquareLowerInnerX",          &trkReducedChisquareLowerInnerX,          "trkReducedChisquareLowerInnerX/F");
  outTree->Branch("trkReducedChisquareLowerInnerY",          &trkReducedChisquareLowerInnerY,          "trkReducedChisquareLowerInnerY/F");

  outTree->Branch("trkFitCodeMS",                            &trkFitCodeMS,                            "trkFitCodeMS/I");
  outTree->Branch("trkRigidityMS",                           &trkRigidityMS,                           "trkRigidityMS/F");
  outTree->Branch("trkRigidityInverseErrorMS",               &trkRigidityInverseErrorMS,               "trkRigidityInverseErrorMS/F");
  outTree->Branch("trkReducedChisquareMSX",                  &trkReducedChisquareMSX,                  "trkReducedChisquareMSX/F");
  outTree->Branch("trkReducedChisquareMSY",                  &trkReducedChisquareMSY,                  "trkReducedChisquareMSY/F");

  outTree->Branch("trkFitCodeInner",                         &trkFitCodeInner,                         "trkFitCodeInner/I");
  outTree->Branch("trkRigidityInner",                        &trkRigidityInner,                        "trkRigidityInner/F");
  outTree->Branch("trkRigidityInverseErrorInner",            &trkRigidityInverseErrorInner,            "trkRigidityInverseErrorInner/F");
  outTree->Branch("trkReducedChisquareInnerX",               &trkReducedChisquareInnerX,               "trkReducedChisquareInnerX/F");
  outTree->Branch("trkReducedChisquareInnerY",               &trkReducedChisquareInnerY,               "trkReducedChisquareInnerY/F");

  outTree->Branch("trkLayerJQ",                              &trkLayerJQ,                              "trkLayerJQ[9]/F");
  outTree->Branch("trkEdepLayerJXSideOK",                    &trkEdepLayerJXSideOK,                    "trkEdepLayerJXSideOK[9]/I");
  outTree->Branch("trkEdepLayerJYSideOK",                    &trkEdepLayerJYSideOK,                    "trkEdepLayerJYSideOK[9]/I");
  outTree->Branch("trkEdepLayerJX",                          &trkEdepLayerJX,                          "trkEdepLayerJX[9]/F");
  outTree->Branch("trkEdepLayerJY",                          &trkEdepLayerJY,                          "trkEdepLayerJY[9]/F");
  outTree->Branch("trkPosXLJ",                               &trkPosXLJ,                               "trkPosXLJ[9]/F");
  outTree->Branch("trkPosYLJ",                               &trkPosYLJ,                               "trkPosYLJ[9]/F");
  outTree->Branch("trkPosZLJ",                               &trkPosZLJ,                               "trkPosZLJ[9]/F");
  outTree->Branch("trkHitCooXLJ",                            &trkHitCooXLJ,                            "trkHitCooXLJ[9]/F");
  outTree->Branch("trkHitCooYLJ",                            &trkHitCooYLJ,                            "trkHitCooYLJ[9]/F");
  outTree->Branch("trkDirThetaLJ",                           &trkDirThetaLJ,                           "trkDirThetaLJ[9]/F");
  outTree->Branch("trkDirPhiLJ",                             &trkDirPhiLJ,                             "trkDirPhiLJ[9]/F");
  outTree->Branch("trkCharge",                               &trkCharge,                               "trkCharge/F");
  outTree->Branch("trkChargeRMS",                            &trkChargeRMS,                            "trkChargeRMS/F");
  outTree->Branch("trkInnerCharge",                          &trkInnerCharge,                          "trkInnerCharge/F");
  outTree->Branch("trkInnerChargeRMS",                       &trkInnerChargeRMS,                       "trkInnerChargeRMS/F");
  outTree->Branch("trkZ",                                    &trkZ,                                    "trkZ/I");
  outTree->Branch("trkInnerZ",                               &trkInnerZ,                               "trkInnerZ/I");
  outTree->Branch("trkLayerJZ",                              &trkLayerJZ,                              "trkLayerJZ[9]/I");
  outTree->Branch("trkHasExtLayers",                         &trkHasExtLayers,                         "trkHasExtLayers/F");

  outTree->Branch("isRichAvailable",                         &isRichAvailable,                         "isRichAvailable/I");
  outTree->Branch("richRebuild",                             &richRebuild,                             "richRebuild/I");
  outTree->Branch("richIsGood",                              &richIsGood,                              "richIsGood/I");
  outTree->Branch("richIsClean",                             &richIsClean,                             "richIsClean/I");
  outTree->Branch("richIsNaF",                               &richIsNaF,                               "richIsNaF/I");
  outTree->Branch("richTileIndex",                           &richTileIndex,                           "richTileIndex/I");
  outTree->Branch("richDistanceTileBorder",                  &richDistanceTileBorder,                  "richDistanceTileBorder/F");
  outTree->Branch("richChargeCorrections",                   &richChargeCorrections,                   "richChargeCorrections/I");
  outTree->Branch("richPMTCorrectionsFailed",                &richPMTCorrectionsFailed,                "richPMTCorrectionsFailed/I");
  outTree->Branch("richTrackEmissionPoint",                  &richTrackEmissionPoint,                  "richTrackEmissionPoint[5]/F");
  outTree->Branch("richUsedHits",                            &richUsedHits,                            "richUsedHits/I");
  outTree->Branch("richRingWidth",                           &richRingWidth,                           "richRingWidth/F");
  outTree->Branch("richWidth",                               &richWidth,                               "richRingWidth/F");
  outTree->Branch("richNHits",                               &richNHits,                               "richNHits/I");
  outTree->Branch("richNPMTsOnRing",                         &richNPMTsOnRing,                         "richNPMTsOnRing/I");
  outTree->Branch("richBeta",                                &richBeta,                                "richBeta/F");
  outTree->Branch("richBetaError",                           &richBetaError,                           "richBetaError/F");
  outTree->Branch("richChargeSquared",                       &richChargeSquared,                       "richChargeSquared/F");
  outTree->Branch("richKolmogorovProbability",               &richKolmogorovProbability,               "richKolmogorovProbability/F");
  outTree->Branch("richPhotoelectrons",                      &richPhotoelectrons,                      "richPhotoelectrons/F");
  outTree->Branch("richExpectedPhotoelectrons",              &richExpectedPhotoelectrons,              "richExpectedPhotoelectrons/F");
  outTree->Branch("richTheta",                               &richTheta,                               "richTheta/F");
  outTree->Branch("richPhi",                                 &richPhi,                                 "richPhi/F");
  outTree->Branch("richBetaConsistency",                     &richBetaConsistency,                     "richBetaConsistency/F");
  outTree->Branch("richReflectedHits",                       &richReflectedHits,                       "richReflectedHits/I");
  outTree->Branch("richPMTChargeConsistency",                &richPMTChargeConsistency,                "richPMTChargeConsistency/F");
  outTree->Branch("richTrackerTrackOnRadiatorX",             &richTrackerTrackOnRadiatorX,             "richTrackerTrackOnRadiatorX/F");
  outTree->Branch("richTrackerTrackOnRadiatorY",             &richTrackerTrackOnRadiatorY,             "richTrackerTrackOnRadiatorY/F");
  outTree->Branch("trdNClusters",                            &trdNClusters,                            "trdNClusters/I");
  outTree->Branch("trdNUnusedHits",                          &trdNUnusedHits,                          "trdNUnusedHits/I");
  outTree->Branch("trdNTracks",                              &trdNTracks,                              "trdNTracks/I");
  outTree->Branch("trdNVertex",                              &trdNVertex,                              "trdNVertex/I");
  outTree->Branch("trdTrackTheta",                           &trdTrackTheta,                           "trdTrackTheta/F");
  outTree->Branch("trdTrackPhi",                             &trdTrackPhi,                             "trdTrackPhi/F");
  outTree->Branch("trdTrackChi2",                            &trdTrackChi2,                            "trdTrackChi2/F");
  outTree->Branch("trdTrackPattern",                         &trdTrackPattern,                         "trdTrackPattern/I");
  outTree->Branch("trdTrackCharge",                          &trdTrackCharge,                          "trdTrackCharge/F");
  outTree->Branch("trdTrackEdepL",                           &trdTrackEdepL,                           "trdTrackEdepL[20]/F");
  outTree->Branch("trdTrackMeanDepositedEnergy",             &trdTrackMeanDepositedEnergy,             "trdTrackMeanDepositedEnergy/F");
  outTree->Branch("trdTrackTotalDepositedEnergy",            &trdTrackTotalDepositedEnergy,            "trdTrackTotalDepositedEnergy/F");
  outTree->Branch("trdTrackDeviationXWithInnerTrk",          &trdTrackDeviationXWithInnerTrk,          "trdTrackDeviationXWithInnerTrk/F");
  outTree->Branch("trdTrackDeviationYWithInnerTrk",          &trdTrackDeviationYWithInnerTrk,          "trdTrackDeviationYWithInnerTrk/F");
  outTree->Branch("trdTrackDistanceBetweenInnerTrk",         &trdTrackDistanceBetweenInnerTrk,         "trdTrackDistanceBetweenInnerTrk/F");
  outTree->Branch("ResidualXBetweenInnerTrackAndSplineTrack",&ResidualXBetweenInnerTrackAndSplineTrack,"ResidualXBetweenInnerTrackAndSplineTrack/F");
  outTree->Branch("ResidualYBetweenInnerTrackAndSplineTrack",&ResidualYBetweenInnerTrackAndSplineTrack,"ResidualYBetweenInnerTrackAndSplineTrack/F");
  outTree->Branch("trdKNRawHits",                            &trdKNRawHits,                            "trdKNRawHits/I");
  outTree->Branch("trdKIsReadAlignmentOK",                   &trdKIsReadAlignmentOK,                   "trdKIsReadAlignmentOK/I");
  outTree->Branch("trdKIsReadCalibOK",                       &trdKIsReadCalibOK,                       "trdKIsReadCalibOK/I");
  outTree->Branch("trdKNHits",                               &trdKNHits,                               "trdKNHits/I");
  outTree->Branch("trdKIsValid",                             &trdKIsValid,                             "trdKIsValid/I");
  outTree->Branch("trdKElectronToProtonLogLikelihoodRatio",  &trdKElectronToProtonLogLikelihoodRatio,  "trdKElectronToProtonLogLikelihoodRatio/F");
  outTree->Branch("trdKHeliumToElectronLogLikelihoodRatio",  &trdKHeliumToElectronLogLikelihoodRatio,  "trdKHeliumToElectronLogLikelihoodRatio/F");
  outTree->Branch("trdKHeliumToProtonLogLikelihoodRatio",    &trdKHeliumToProtonLogLikelihoodRatio,    "trdKHeliumToProtonLogLikelihoodRatio/F");
  outTree->Branch("trdKCharge",                              &trdKCharge,                              "trdCharge/F");
  outTree->Branch("trdKChargeError",                         &trdKChargeError,                         "trdKChargeError/F");
  outTree->Branch("trdKNUsedHitsForCharge",                  &trdKNUsedHitsForCharge,                  "trdKNUsedHitsForCharge/I");
  outTree->Branch("trdKAmpLayer",                            &trdKAmpLayer,                            "trdKAmpLayer[20]/F");
  outTree->Branch("trdKTotalPathLength",                     &trdKTotalPathLength,                     "trdKTotalPathLength/F");
  outTree->Branch("trdKTotalAmp",                            &trdKTotalAmp,                            "trdKTotalAmp/F");
  outTree->Branch("trdKElectronLikelihood",                  &trdKElectronLikelihood,                  "trdKElectronLikelihood/F");
  outTree->Branch("trdKProtonLikelihood",                    &trdKProtonLikelihood,                    "trdKProtonLikelihood/F");
  outTree->Branch("trdKHeliumLikelihood",                    &trdKHeliumLikelihood,                    "trdKHeliumLikelihood/F");
  outTree->Branch("trdPElectronToProtonLogLikelihoodRatio",  &trdPElectronToProtonLogLikelihoodRatio,  "trdPElectronToProtonLogLikelihoodRatio/F");
  outTree->Branch("trdPHeliumToElectronLogLikelihoodRatio",  &trdPHeliumToElectronLogLikelihoodRatio,  "trdPHeliumToElectronLogLikelihoodRatio/F");
  outTree->Branch("trdPHeliumToProtonLogLikelihoodRatio",    &trdPHeliumToProtonLogLikelihoodRatio,    "trdPHeliumToProtonLogLikelihoodRatio/F");
  outTree->Branch("trdPElectronLikelihood",                  &trdPElectronLikelihood,                  "trdPElectronLikelihood/F");
  outTree->Branch("trdPProtonLikelihood",                    &trdPProtonLikelihood,                    "trdPProtonLikelihood/F");
  outTree->Branch("trdPHeliumLikelihood",                    &trdPHeliumLikelihood,                    "trdPHeliumLikelihood/F");
  outTree->Branch("trdPNHybridHits",                         &trdPNHybridHits,                         "trdPNHybridHits/I");
  outTree->Branch("trdPNActiveLayers",                       &trdPNActiveLayers,                       "trdPNActiveLayers/I");
  outTree->Branch("trdPNLowAdcHits",                         &trdPNLowAdcHits,                         "trdPNLowAdcHits/I");
  outTree->Branch("trdPNLowDxHits",                          &trdPNLowDxHits,                          "trdPNLowDxHits/I");
  outTree->Branch("accNHits",                                     &accNHits,                                     "accNHits/I");
  outTree->Branch("accNRecoClPG",                                 &accNRecoClPG,                                 "accNRecoClPG/I");
  std::cout << "Initialization Complete." << std::endl;


  hEvtCounter = new TH1F("hEvtCounter", "Number of survived events after cut", 19, 0., 19.);
  hEvtCounter->GetXaxis()->SetBinLabel(1, "No AMSEventR");
  hEvtCounter->GetXaxis()->SetBinLabel(2, "Bad hardware");
  hEvtCounter->GetXaxis()->SetBinLabel(3, "Unbiased physics trigger");
  hEvtCounter->GetXaxis()->SetBinLabel(4, "SAA test");
  hEvtCounter->GetXaxis()->SetBinLabel(5, "No RTI object");
  hEvtCounter->GetXaxis()->SetBinLabel(6, "Livetime > 0.7");
  hEvtCounter->GetXaxis()->SetBinLabel(7, "Tracker alignment");
  hEvtCounter->GetXaxis()->SetBinLabel(8, "Single particle event");
  hEvtCounter->GetXaxis()->SetBinLabel(9, "Particle pointer exist");
  hEvtCounter->GetXaxis()->SetBinLabel(10, "BetaH object test");
  hEvtCounter->GetXaxis()->SetBinLabel(11, "TrTrack object test");
  hEvtCounter->GetXaxis()->SetBinLabel(12, "TrdTrack object test");
  hEvtCounter->GetXaxis()->SetBinLabel(13, "Inner track fit test");
  hEvtCounter->GetXaxis()->SetBinLabel(14, "Hits on L2, L3/4, L5/6, L7/8");
  hEvtCounter->GetXaxis()->SetBinLabel(15, "Single track event");
  hEvtCounter->GetXaxis()->SetBinLabel(16, "Track Chi2XInner");
  hEvtCounter->GetXaxis()->SetBinLabel(17, "Beta > 0");
  hEvtCounter->GetXaxis()->SetBinLabel(18, "IsGoodBeta?");
  hEvtCounter->GetXaxis()->SetBinLabel(19, "IsTkTofMatch?");
  //===========================================================================
  //
  //
  //
  //   Initialization Finished
  //
  //
  //
  //===========================================================================

  AMSSetupR::RTI::UseLatest(6);
  TkDBc::UseFinal();
  TRMCFFKEY.ReadFromFile = 0;
  TRFITFFKEY.ReadFromFile = 0;
  TRFITFFKEY.magtemp = 0;

  for( unsigned long int e = 0; e < nEvents; e++ )
  {
    ///// Initialize
    nRun = 0;
    nEvent = 0;
    nLevel = 0;
    nParticle = 0;
    nCharge = 0;
    nTrTrack = 0;
    nTrdTrack = 0;
    nAntiCluster = 0;
    nTofClustersInTime = 0;
    nRichRing = 0;
    nRichRingB = 0;
    nBeta = 0;
    nBetaB = 0;
    nBetaH = 0;
    nEcalShower = 0;
    nVertex = 0;

    liveTime = 0;
    utcTime = 0;
    orbitAltitude = 0;
    orbitLatitude = 0;
    orbitLongitude = 0;
    orbitLatitudeM = 0;
    orbitLongitudeM = 0;
    velR = -30;
    velTheta = -30;
    velPhi = -30;
    yaw = -30;
    pitch = -30;
    roll = -30;
    gLongitude = -30;
    gLatitude = -30;
    sunPosAzimuth = -30;
    sunPosElevation = -30;
    sunPosCalcResult = -30;
    unixTime = -1;
    solarArrayCoord[0] = -999;
    solarArrayCoord[1] = -999;
    solarArrayCoord[2] = -999;
    ptlCutOffStoermer = -3333333;
    ptlCutOffDipole = -3333333;
    ptlCutOffIGRF = -3333333;
    isInShadow = -1;
    zenithAngle = -30;
    isInSAA = -1;

    tofNCluster = -1;
    tofNClusterH = -1;
    tofNUsedHits = -1;
    tofNUnusedHits = -1;
    tofNUsedLayersForQ = -1;
    tofBeta = -99;
    tofBetaH = -99;
    tofInvBetaErr = -99;
    tofInvBetaErrH = -99;
    tofInvBetaErrC = -99;
    tofNormEBetaV = -99;
    tofBetaSH = -99;
    tofBetaC = -99;
    tofBetaCH = -99;
    tofEBetaCV = -99;
    tofMass = -99;
    tofMassError = -99;
    isGoodBeta = -1;
    isTkTofMatch = -1;
    tofChisqT = -1;
    tofChisqC = -1;
    tofReducedChisqT = -1;
    tofReducedChisqC = -1;
    tofSumHit = -1;
    tofNUsedClusterH = -1;
    tofPatternOnLayer[0] = -1;
    tofPatternOnLayer[1] = -1;
    tofPatternOnLayer[2] = -1;
    tofPatternOnLayer[3] = -1;
    tofDepositedEnergyOnLayer[0] = -1;
    tofDepositedEnergyOnLayer[1] = -1;
    tofDepositedEnergyOnLayer[2] = -1;
    tofDepositedEnergyOnLayer[3] = -1;
    tofEstimatedChargeOnLayer[0] = -1;
    tofEstimatedChargeOnLayer[1] = -1;
    tofEstimatedChargeOnLayer[2] = -1;
    tofEstimatedChargeOnLayer[3] = -1;
    tofCharge = -1;
    tofChargeRMS = -1;
    tofUpperCharge = -1;
    tofLowerCharge = -1;
    tofChargeOnLayer[0] = -1;
    tofChargeOnLayer[1] = -1;
    tofChargeOnLayer[2] = -1;
    tofChargeOnLayer[3] = -1;
    tofTimeOnLayer[0] = -9999;
    tofTimeOnLayer[1] = -9999;
    tofTimeOnLayer[2] = -9999;
    tofTimeOnLayer[3] = -9999;
    tofETimeOnLayer[0] = -9999;
    tofETimeOnLayer[1] = -9999;
    tofETimeOnLayer[2] = -9999;
    tofETimeOnLayer[3] = -9999;
    tofETCooOnLayer[0] = -999;
    tofETCooOnLayer[1] = -999;
    tofETCooOnLayer[2] = -999;
    tofETCooOnLayer[3] = -999;
    tofTResidualOnLayer[0] = -999;
    tofTResidualOnLayer[1] = -999;
    tofTResidualOnLayer[2] = -999;
    tofTResidualOnLayer[3] = -999;
    tofT0 = -99;
    tofTkTFLenOnLayer[0] = -999;
    tofTkTFLenOnLayer[1] = -999;
    tofTkTFLenOnLayer[2] = -999;
    tofTkTFLenOnLayer[3] = -999;
    isEcalAvailable = 0;
    showerEnergyD = -1.;
    showerEnergyDL[0] = 0;
    showerEnergyDL[1] = 0;
    showerEnergyDL[2] = 0;
    showerEnergyDL[3] = 0;
    showerEnergyDL[4] = 0;
    showerEnergyDL[5] = 0;
    showerEnergyDL[6] = 0;
    showerEnergyDL[7] = 0;
    showerEnergyDL[8] = 0;
    showerEnergyDL[9] = 0;
    showerEnergyDL[10] = 0;
    showerEnergyDL[11] = 0;
    showerEnergyDL[12] = 0;
    showerEnergyDL[13] = 0;
    showerEnergyDL[14] = 0;
    showerEnergyDL[15] = 0;
    showerEnergyDL[16] = 0;
    showerEnergyDL[17] = 0;
    showerEnergyE = -1.;
    showerEnergyCorrected = -1.;
    showerBDT = -99.;
    showerCofG[0] = -99;
    showerCofG[1] = -99;
    showerCofG[2] = -99;
    showerCofGDist = -99;
    showerCofGdX = -99;
    showerCofGdY = -99;
    trkFitCodeFS = -1;
    trkRigidityFS = -9999;
    trkRigidityInverseErrorFS = -9999;
    trkReducedChisquareFSX = -9999;
    trkReducedChisquareFSY = -9999;
    trkFitCodeMS = -1;
    trkRigidityMS = -9999;
    trkRigidityInverseErrorMS = -9999;
    trkReducedChisquareMSY = -9999;
    trkReducedChisquareMSX = -9999;
    trkFitCodeInner = -1;
    trkRigidityInverseErrorInner = -9999;
    trkReducedChisquareInnerX = -9999;
    trkReducedChisquareInnerY = -9999;
    trkLayerJQ[0] = -1;
    trkLayerJQ[1] = -1;
    trkLayerJQ[2] = -1;
    trkLayerJQ[3] = -1;
    trkLayerJQ[4] = -1;
    trkLayerJQ[5] = -1;
    trkLayerJQ[6] = -1;
    trkLayerJQ[7] = -1;
    trkLayerJQ[8] = -1;
    trkEdepLayerJXSideOK[0] = -1;
    trkEdepLayerJXSideOK[1] = -1;
    trkEdepLayerJXSideOK[2] = -1;
    trkEdepLayerJXSideOK[3] = -1;
    trkEdepLayerJXSideOK[4] = -1;
    trkEdepLayerJXSideOK[5] = -1;
    trkEdepLayerJXSideOK[6] = -1;
    trkEdepLayerJXSideOK[7] = -1;
    trkEdepLayerJXSideOK[8] = -1;
    trkEdepLayerJYSideOK[0] = -1;
    trkEdepLayerJYSideOK[1] = -1;
    trkEdepLayerJYSideOK[2] = -1;
    trkEdepLayerJYSideOK[3] = -1;
    trkEdepLayerJYSideOK[4] = -1;
    trkEdepLayerJYSideOK[5] = -1;
    trkEdepLayerJYSideOK[6] = -1;
    trkEdepLayerJYSideOK[7] = -1;
    trkEdepLayerJYSideOK[8] = -1;
    trkEdepLayerJX[0] = -12345;
    trkEdepLayerJX[1] = -12345;
    trkEdepLayerJX[2] = -12345;
    trkEdepLayerJX[3] = -12345;
    trkEdepLayerJX[4] = -12345;
    trkEdepLayerJX[5] = -12345;
    trkEdepLayerJX[6] = -12345;
    trkEdepLayerJX[7] = -12345;
    trkEdepLayerJX[8] = -12345;
    trkEdepLayerJY[0] = -12345;
    trkEdepLayerJY[1] = -12345;
    trkEdepLayerJY[2] = -12345;
    trkEdepLayerJY[3] = -12345;
    trkEdepLayerJY[4] = -12345;
    trkEdepLayerJY[5] = -12345;
    trkEdepLayerJY[6] = -12345;
    trkEdepLayerJY[7] = -12345;
    trkEdepLayerJY[8] = -12345;
    trkCharge = -1;
    trkInnerCharge = -1;
    trkHasExtLayers = -1;
    isRichAvailable = 0;
    richRebuild = -1;
    richIsGood = -1;
    richIsClean = -1;
    richIsNaF = -1;
    richUsedHits = -1;
    richRingWidth = -1;
    richNHits = -1;
    richNPMTsOnRing = -1;
    richBeta = -1;
    richBetaError = -1;
    richChargeSquared = -1;
    richKolmogorovProbability = -1;
    richPhotoelectrons = -1;
    richPhotoElectrons = -1;
    richExpectedPhotoelectrons = -1;
    richExpectedPhotoElectrons = -1;
    richNUsedHits = -1;
    richTheta = -9;
    richPhi = -9;
    trdNClusters = 0;
    trdNUnusedHits = 0;
    trdNTracks = 0;
    trdNVertex = 0;
    trdTrackPhi = -9.;
    trdTrackTheta = -9.;
    trdTrackChi2 = -9.;
    trdTrackPattern = -9.;
    trdTrackCharge = -9.;
    trdTrackEdepL[0] = -1;
    trdTrackEdepL[1] = -1;
    trdTrackEdepL[2] = -1;
    trdTrackEdepL[3] = -1;
    trdTrackEdepL[4] = -1;
    trdTrackEdepL[5] = -1;
    trdTrackEdepL[6] = -1;
    trdTrackEdepL[7] = -1;
    trdTrackEdepL[8] = -1;
    trdTrackEdepL[9] = -1;
    trdTrackEdepL[10] = -1;
    trdTrackEdepL[11] = -1;
    trdTrackEdepL[12] = -1;
    trdTrackEdepL[13] = -1;
    trdTrackEdepL[14] = -1;
    trdTrackEdepL[15] = -1;
    trdTrackEdepL[16] = -1;
    trdTrackEdepL[17] = -1;
    trdTrackEdepL[18] = -1;
    trdTrackEdepL[19] = -1;
    trdTrackMeanDepositedEnergy  = -1;
    trdTrackTotalDepositedEnergy = -1;
    trdKNRawHits                                = -1;
    trdKIsReadAlignmentOK                       = -1;
    trdKIsReadCalibOK                           = -1;
    trdKNHits                                   = -1;
    trdKIsValid                                 = -1;
    trdKElectronToProtonLogLikelihoodRatio      = -1;
    trdKHeliumToElectronLogLikelihoodRatio      = -1;
    trdKHeliumToProtonLogLikelihoodRatio        = -1;
    trdKCharge                                  = -1;
    trdKChargeError                             = -1;
    trdKNUsedHitsForCharge                      = -1;
    for(int k = 0; k < 20; k++) trdKAmpLayer[k] = -1;
    trdKTotalPathLength                         = -1;
    trdKElectronLikelihood                      = -1;
    trdKProtonLikelihood                        = -1;
    trdKHeliumLikelihood                        = -1;
    trdPElectronToProtonLogLikelihoodRatio = 0;
    trdPHeliumToElectronLogLikelihoodRatio = 0;
    trdPHeliumToProtonLogLikelihoodRatio = 0;
    trdPElectronLikelihood = 0;
    trdPProtonLikelihood = 0;
    trdPHeliumLikelihood = 0;
    trdPNHybridHits = 0;
    if( (int)e % 10000 == 0 || e == nEvents - 1 )
    {
      std::cout << "Processing " << e << " out of " << nEvents << " (" << (float)e/nEvents*100. << "%), "
        << "nSelected = " << nSelected << ", currently, selection efficiency is (" << (float)nSelected/e*100. << "%) " << std::endl;
    }

    // Get event
    pAMSEvent = chain.GetEvent(e);
    pHeader = &( pAMSEvent->fHeader );

    /// =======================================================================
    //
    //
    //            Cuts
    //
    //
    //
    /// =======================================================================
    if( !pAMSEvent ) continue;
    hEvtCounter->Fill(0);

    // Hardware status check
    goodhw  = true;
    hwerror = false;

    for( int i = 0; i < pAMSEvent->nDaqEvent(); ++i )
    {
      pDaqEvent = pAMSEvent->pDaqEvent(i);
      for( int iJINJ = 0; iJINJ < 4;  iJINJ++ ) hwerror |= (bool)( ( pDaqEvent->JINJStatus[iJINJ]>>8 ) & 0x7F );
      for( int iJErr = 0; iJErr < 24; iJErr++ ) hwerror |= (bool)(   pDaqEvent->JError[iJErr] & 0x7F );
      if (hwerror) goodhw &= false;
    }
    if( goodhw == false ) continue;
    hEvtCounter->Fill(1);


    // Unbiased Physics Trigger Check
    nLevel1 = pAMSEvent->nLevel1();
    unbiased = false;

    for( int i = 0; i < nLevel1; ++i )
    {
      pLevel1 = pAMSEvent->pLevel1(i);
      bitset<8> physicsBitPattern(pLevel1->PhysBPatt);

      if( physicsBitPattern.test(1) || physicsBitPattern.test(2) || physicsBitPattern.test(3) || physicsBitPattern.test(4) || physicsBitPattern.test(5))
        unbiased |= false;
      else
        unbiased |= true;
    }
    if( unbiased == true ) continue;
    hEvtCounter->Fill(2);

    /// SAA test
    if( pAMSEvent->IsInSAA() ) continue;
    hEvtCounter->Fill(3);

    /// Tracker alignment test
    AMSSetupR::RTI rti;
    if( pAMSEvent->GetRTI(rti) != 0 ) continue;
    hEvtCounter->Fill(4);

    if( rti.lf < 0.7 ) continue;
    hEvtCounter->Fill(5);
    AMSPoint pn1, pn9, pd1, pd9;
    pAMSEvent->GetRTIdL1L9( 0, pn1, pd1, pAMSEvent->UTime(), 60);
    pAMSEvent->GetRTIdL1L9( 1, pn9, pd9, pAMSEvent->UTime(), 60);
    if( pd1.y() > 35 || pd9.y() > 45) continue;
    hEvtCounter->Fill(6);

    /// Single particle test
    if( pAMSEvent->nParticle() != 1 ) continue;
    hEvtCounter->Fill(7);
    pParticle = pAMSEvent->pParticle(0);
    if( pParticle == NULL ) continue;
    hEvtCounter->Fill(8);

    /// Object test
    pBetaH = pParticle->pBetaH();
    if( pBetaH == NULL ) continue;
    hEvtCounter->Fill(9);

    pBeta = pParticle->pBeta();
    if( pBeta == NULL ) continue;
    hEvtCounter->Fill(10);

    pTrTrack = pParticle->pTrTrack();
    if( pTrTrack == NULL ) continue;
    hEvtCounter->Fill(11);

    /// Good track test
    q_inn = pTrTrack->GetInnerQ();
    z_inn = floor( q_inn + 0.5);

    if( z_inn < 1 )
    {
      z_inn = 1;
      mass = TrFit::Mproton;
    }
    if( z_inn >= 2 ) mass = TrFit::Mhelium;
    if( z_inn >= 3 ) mass = TrFit::Mproton*2*z_inn;
    refit = 3;

    // For ISS, enable ion leinearity (charge Z >= 3 ) and dZ correction
    TRCLFFKEY.UseNonLinearity = 0; // Off default linearity correction option in datacard (it speed up refit)
    TRFITFFKEY.Zshift = 2; // Ion dZ correction.
    if( z_inn >= 3 )
      TrClusterR::SetLinearityCorrection(); // Linearity correction (important for ions with 2<Z<13)
    else
      TrClusterR::UnsetLinearityCorrection(); // Off linearity correction (Z< 3)
    refit = 3;

    id_inner = 0;
    id_inner = pTrTrack->iTrTrackPar(1, 3, 20+refit, mass, z_inn);
    if(  id_inner < 0 )
    {
      if( debug )
      {
        switch( id_inner )
        {
          case -1: cerr << "[IN]The requested fit cannot be performed on this track." << endl; break;
          case -2: cerr << "[IN]The requested fit it is not available without refitting." << endl; break;
          case -3: cerr << "[IN]The refit failed." << endl; break;
          case -4: cerr << "[IN]Should not happen!! Contact the developers!." << endl; break;
          case -5: cerr << "[IN]The refit failed because there was a problem retrieving dynamic alignment for Ext Layers" << endl; break;
          default: break;
        }
      }
      continue;       //// Skip this event when there is no reconstructed track for inner span
    }
    hEvtCounter->Fill(12);

    for( int ilayer = 0; ilayer < 9; ilayer++ )
    {
      xside[ilayer] = false;
      yside[ilayer] = false;
    }
    for( int ilayer = 0; ilayer < 9; ilayer++ )
    {
      pTrRecHit = pTrTrack->GetHitLJ(ilayer+1);       // argument should be in the range from [1,9]
      if( !pTrRecHit ) continue;
      else
      {
        if( pTrRecHit->GetEdep(0) != 0 ) { xside[ilayer] = true; }
        else xside[ilayer] = false;
        if( pTrRecHit->GetEdep(1) != 0 ) { yside[ilayer] = true; }
        else yside[ilayer] = false;
      }
    }

    if( !(yside[1])             ) { if( debug ) { std::cerr << "Event: " << e << ", No hit on Tracker layer 2" << std::endl; }   continue;}
    if( !(yside[2] || yside[3]) ) { if( debug ) { std::cerr << "Event: " << e << ", No hit on Tracker layer 3/4" << std::endl; } continue;}
    if( !(yside[4] || yside[5]) ) { if( debug ) { std::cerr << "Event: " << e << ", No hit on Tracker layer 5/6" << std::endl; } continue;}
    if( !(yside[7] || yside[8]) ) { if( debug ) { std::cerr << "Event: " << e << ", No hit on Tracker layer 7/8" << std::endl; } continue;}
    //if( !(xside[0] && yside[0]) ) { std::cerr << "Event: " << e << ", No hit on L1X and L1Y!" << std::endl;      continue;}
    hEvtCounter->Fill(13);
    trkReducedChisquareInnerX = pTrTrack->GetNormChisqX( id_inner );
    trkReducedChisquareInnerY = pTrTrack->GetNormChisqY( id_inner );


    // Single track event only
    if( pAMSEvent->nTrTrack() != 1 ) continue;
    hEvtCounter->Fill(14);

    // Track quality
    if( trkReducedChisquareInnerX > 20 )    continue;
    hEvtCounter->Fill(15);
    // TOF quality
    if( pBetaH->GetBeta()         < 0 )     continue;
    hEvtCounter->Fill(16);
    if( pBetaH->IsGoodBeta()     == false ) continue;
    hEvtCounter->Fill(17);
    if( pBetaH->IsTkTofMatch()   == false ) continue;
    hEvtCounter->Fill(18);

    /// =======================================================================
    ///
    //
    //    End of cut
    //
    //
    //
    /// =======================================================================

    nSelected++;
    // End of Cut

    // Save variables
    nRun            = pAMSEvent->Run();
    nEvent          = pAMSEvent->Event();
    nLevel1         = pAMSEvent->nLevel1();
    nParticle       = pAMSEvent->nParticle();
    nCharge         = pAMSEvent->nCharge();
    nTrTrack        = pAMSEvent->nTrTrack();
    nTrdTrack       = pAMSEvent->nTrdTrack();
    nAntiCluster    = pAMSEvent->nAntiCluster();
    nRichRing       = pAMSEvent->nRichRing();
    nRichRingB      = pAMSEvent->nRichRingB();
    nBeta           = pAMSEvent->nBeta();
    nBetaB          = pAMSEvent->nBetaB();
    nBetaH          = pAMSEvent->nBetaH();
    nEcalShower     = pAMSEvent->nEcalShower();
    nVertex         = pAMSEvent->nVertex();
    liveTime        = pAMSEvent->LiveTime();
    utcTime         = (float)pHeader->UTCTime(0);
    orbitAltitude   = pHeader->RadS;
    orbitLatitude   = pHeader->ThetaS;
    orbitLongitude  = pHeader->PhiS;
    orbitLatitudeM  = pHeader->ThetaM;
    orbitLongitudeM = pHeader->PhiM;
    velR            = pHeader->VelocityS;
    velTheta        = pHeader->VelTheta;
    velPhi          = pHeader->VelPhi;
    yaw             = pHeader->Yaw;
    pitch           = pHeader->Pitch;
    roll            = pHeader->Roll;
    zenithAngle     = pHeader->Zenith();
    isInSAA         = pAMSEvent->IsInSAA();

    double tmpglong;
    double tmpglat;
    gCoordCalcResult  = pAMSEvent->GetGalCoo(
        gCoordCalcResultBit,                           // Galactic coordinate calculation result bit
        tmpglong,                                      // Galactic longitude
        tmpglat);                                      // Galactic latitude
    gLongitude  = tmpglong;
    gLatitude   = tmpglat;

    double tmpSunPosAzimuth;
    double tmpSunPosElevation;
    sunPosCalcResult   = pHeader->getSunAMS(tmpSunPosAzimuth, tmpSunPosElevation);
    sunPosAzimuth      = tmpSunPosAzimuth;
    sunPosElevation    = tmpSunPosElevation;

    AMSPoint solarArray;
    isInShadow         = pAMSEvent->isInShadow(solarArray, 0);
    solarArrayCoord[0] = solarArray.x();
    solarArrayCoord[1] = solarArray.y();
    solarArrayCoord[2] = solarArray.z();

    ptlCutOffStoermer = pParticle->CutoffS;
    ptlCutOffDipole = pParticle->Cutoff;

    // Object preparation
    //pTrTrack = pParticle->pTrTrack();
    //int id_inner = pTrTrack->iTrTrackPar(1, 3, 20);
    //pBetaH = pParticle->pBetaH();
    pEcalShower = pParticle->pEcalShower();


    // Tracker
    if( id_inner < 0 ) cout << "ERROR!" << endl;
    // Inner tracker
    trkFitCodeInner              = id_inner;
    trkRigidityInner             = pTrTrack->GetRigidity(id_inner);
    trkRigidityInverseErrorInner = pTrTrack->GetErrRinv(id_inner);
    trkReducedChisquareInnerX    = pTrTrack->GetNormChisqX(id_inner);
    trkReducedChisquareInnerY    = pTrTrack->GetNormChisqY(id_inner);

    // Upper half of inner tracker
    trkFitCodeUpperInner = pTrTrack->iTrTrackPar(1, 1, 23);
    if( trkFitCodeUpperInner > 0 )
    {
      trkRigidityUpperInner             = pTrTrack->GetRigidity(trkFitCodeUpperInner);
      trkRigidityInverseErrorUpperInner = pTrTrack->GetErrRinv(trkFitCodeUpperInner);
      trkReducedChisquareUpperInnerX    = pTrTrack->GetNormChisqX(trkFitCodeUpperInner);
      trkReducedChisquareUpperInnerY    = pTrTrack->GetNormChisqY(trkFitCodeUpperInner);
    }

    // Lower half of inner tracker
    trkFitCodeLowerInner = pTrTrack->iTrTrackPar(1, 2, 23);
    if( trkFitCodeLowerInner > 0 )
    {
      trkRigidityLowerInner             = pTrTrack->GetRigidity(trkFitCodeLowerInner);
      trkRigidityInverseErrorLowerInner = pTrTrack->GetErrRinv(trkFitCodeLowerInner);
      trkReducedChisquareLowerInnerX    = pTrTrack->GetNormChisqX(trkFitCodeLowerInner);
      trkReducedChisquareLowerInnerY    = pTrTrack->GetNormChisqY(trkFitCodeLowerInner);
    }

    // L1 + inner tracker
    trkFitCodeL1Inner = pTrTrack->iTrTrackPar(1, 5, 20);
    if( trkFitCodeL1Inner > 0 )
    {
      trkRigidityL1Inner             = pTrTrack->GetRigidity(trkFitCodeL1Inner);
      trkRigidityInverseErrorL1Inner = pTrTrack->GetErrRinv(trkFitCodeL1Inner);
      trkReducedChisquareL1InnerX    = pTrTrack->GetNormChisqX(trkFitCodeL1Inner);
      trkReducedChisquareL1InnerY    = pTrTrack->GetNormChisqY(trkFitCodeL1Inner);
    }

    // L9 + inner tracker
    trkFitCodeL9Inner = pTrTrack->iTrTrackPar(1, 6, 20);
    if( trkFitCodeL9Inner > 0 )
    {
      trkRigidityL9Inner             = pTrTrack->GetRigidity(trkFitCodeL9Inner);
      trkRigidityInverseErrorL9Inner = pTrTrack->GetErrRinv(trkFitCodeL9Inner);
      trkReducedChisquareL9InnerX    = pTrTrack->GetNormChisqX(trkFitCodeL9Inner);
      trkReducedChisquareL9InnerY    = pTrTrack->GetNormChisqY(trkFitCodeL9Inner);
    }

    // Full span tracker
    trkFitCodeFS = pTrTrack->iTrTrackPar(1, 7, 20);
    if( trkFitCodeFS > 0 )
    {
      trkRigidityFS             = pTrTrack->GetRigidity(trkFitCodeFS);
      trkRigidityInverseErrorFS = pTrTrack->GetErrRinv(trkFitCodeFS);
      trkReducedChisquareFSX    = pTrTrack->GetNormChisqX(trkFitCodeFS);
      trkReducedChisquareFSY    = pTrTrack->GetNormChisqY(trkFitCodeFS);
    }

    // Maximum span tracker
    trkFitCodeMS = pTrTrack->iTrTrackPar(1, 0, 20);
    if( trkFitCodeMS > 0 )
    {
      trkRigidityMS             = pTrTrack->GetRigidity(trkFitCodeMS);
      trkRigidityInverseErrorMS = pTrTrack->GetErrRinv(trkFitCodeMS);
      trkReducedChisquareMSX    = pTrTrack->GetNormChisqX(trkFitCodeMS);
      trkReducedChisquareMSY    = pTrTrack->GetNormChisqY(trkFitCodeMS);
    }

    // Get IGRF Cutoff
    AMSDir inc_dir;
    inc_dir.SetTheta(pParticle->Theta);
    inc_dir.SetPhi(pParticle->Phi);
    double cf;
    int rsign=(trkRigidityInner>=0)?1:-1;
    int flag_igrf = pAMSEvent->GetIGRFCutoff(cf, rsign, inc_dir);

    if( debug )
    {
      std::cout << "Event: " << e << ", R(Inner): " << trkRigidityInner <<
        ", rsign: " << rsign <<
        ", inc_dir.Theta: " << inc_dir.gettheta() <<
        ", inc_dir.Phi: " << inc_dir.getphi() <<
        ", flag_igrf: " << flag_igrf <<
        ", cf: " << cf << std::endl;
    }
    ptlCutOffIGRF = (float)cf;

    // Now we calculate Tracker dE/dx related information

    // Initialize
    for(int ii = 0; ii < 9; ii++)
    {
      trkLayerJQ[ii] = 0.;
      trkEdepLayerJX[ii] = 0.;
      trkEdepLayerJY[ii] = 0.;
      trkEdepLayerJXSideOK[ii] = -1.;
      trkEdepLayerJYSideOK[ii] = -1.;
    }

    AMSPoint amspnt;
    AMSDir amsdir;
    for( int ilayer = 0; ilayer < 9; ilayer++)     // loops over from ilayer = 0 to ilayer = 8, total 9 loops.
    {
      trkLayerJQ[ilayer] = pTrTrack->GetLayerJQ(ilayer+1);
      pTrRecHit = pTrTrack->GetHitLJ(ilayer+1);
      if( !pTrRecHit )
      {
        trkEdepLayerJXSideOK[ilayer] = -2;
        trkEdepLayerJYSideOK[ilayer] = -2;
      }
      else
      {
        if( pTrRecHit->GetEdep(0) > 0 ) trkEdepLayerJXSideOK[ilayer] = 1;
        else trkEdepLayerJXSideOK[ilayer] = 0;
        if( pTrRecHit->GetEdep(1) > 0 ) trkEdepLayerJYSideOK[ilayer] = 1;
        else trkEdepLayerJYSideOK[ilayer] = 0;
        trkEdepLayerJX[ilayer] += pTrRecHit->GetEdep(0);
        trkEdepLayerJY[ilayer] += pTrRecHit->GetEdep(1);
      }

      // Tracker coordinate information
      pTrTrack->InterpolateLayerJ(ilayer+1, amspnt, amsdir);
      trkPosXLJ[ilayer]     = amspnt.x();
      trkPosYLJ[ilayer]     = amspnt.y();
      trkPosZLJ[ilayer]     = amspnt.z();
      trkHitCooXLJ[ilayer]  = pTrTrack->GetHitCooLJ(ilayer).x();
      trkHitCooYLJ[ilayer]  = pTrTrack->GetHitCooLJ(ilayer).y();
      trkDirThetaLJ[ilayer] = amsdir.gettheta();
      trkDirPhiLJ[ilayer]   = amsdir.getphi();
    }

    // 
    //
    // TOF Section
    //
    //

    tofNCluster         = pAMSEvent->nTofCluster();
    tofNClusterH        = pAMSEvent->nTofClusterH();
    if(!pBetaH) { std::cerr << "Event: " << e << " NULL pBetaH pointer! This should not be happen! Exit the program!" << std::endl; exit(1); }
    if(!pBeta)  { std::cerr << "Event: " << e << " NULL pBeta pointer!  This should not be happen! Exit the program!" << std::endl; exit(1); }
    tofBeta             = pBeta->Beta;
    tofInvBetaErr       = pBeta->Error;
    tofBetaC            = pBeta->BetaC;
    tofInvBetaErrC      = pBeta->ErrorC;
    tofBetaH            = pBetaH->GetBeta();
    tofInvBetaErrH      = pBetaH->GetEBetaV();
    tofNormEBetaV       = pBetaH->GetNormEBetaV();
    tofBetaSH           = pBetaH->GetBetaS();
    tofBetaCH           = pBetaH->GetBetaC();
    tofEBetaCV          = pBetaH->GetEBetaCV();
    tofMass             = pBetaH->GetMass();
    tofMassError        = pBetaH->GetEMass();
    tofNUsedHits        = pBetaH->NTofClusterH();
    tofNUnusedHits      = tofNCluster - tofNUsedHits;
    if( pBetaH->IsGoodBeta()   == true ) isGoodBeta = 1;
    else isGoodBeta     = 0;
    if( pBetaH->IsTkTofMatch() == true ) isTkTofMatch = 1;
    else isTkTofMatch   = 0;
    tofChisqT           = pBeta->Chi2;
    tofChisqC           = pBeta->Chi2S;
    tofChisqTH          = pBetaH->GetChi2T();
    tofChisqCH          = pBetaH->GetChi2C();
    tofReducedChisqT    = pBetaH->GetNormChi2T();
    tofReducedChisqC    = pBetaH->GetNormChi2C();
    tofSumHit           = pBetaH->GetSumHit();
    tofNUsedClusterH    = pBetaH->GetUseHit();
    for( int k = 0; k < 4; k++ )
    {
      tofPatternOnLayer[k] = pBetaH->GetPattern(k);
      tofTimeOnLayer[k] = pBetaH->GetTime(k);
      tofETimeOnLayer[k] = pBetaH->GetETime(k);
      tofETCooOnLayer[k] = pBetaH->GetETCoo(k);
      tofTResidualOnLayer[k] = pBetaH->GetTResidual(k);
      tofTkTFLenOnLayer[k] = pBetaH->GetTkTFLen(k);
      tofCResidualXOnLayer[k] = pBetaH->GetCResidual(k, 0);
      tofCResidualYOnLayer[k] = pBetaH->GetCResidual(k, 1);
    }

    tofCharge = pBetaH->GetQ(tofNUsedLayersForQ, tofChargeRMS);
    tofChargeOnLayer[0] = pBetaH->GetQL(0);
    tofChargeOnLayer[1] = pBetaH->GetQL(1);
    tofChargeOnLayer[2] = pBetaH->GetQL(2);
    tofChargeOnLayer[3] = pBetaH->GetQL(3);
    tofZ = pBetaH->GetZ(nTofLUsedForZ, probTOFZ);

    for(int i = 0; i < 4; ++i)
    {
      TofClusterHR* pTofClusterH = pBetaH->GetClusterHL(i);
      if( pTofClusterH != 0 )
      {
        tofDepositedEnergyOnLayer[i] = pTofClusterH->GetEdep();
      }
      else
        tofDepositedEnergyOnLayer[i] = -1;
    }

    int ncls[4] = {0, 0, 0, 0};
    nTofClustersInTime  = pAMSEvent->GetNTofClustersInTime(pBetaH, ncls);

    // Tracker charge information
    // !! Tracker charge information can be acquired correctly after beta information given!
    trkCharge = (pTrTrack->GetQ_all(tofBeta, trkFitCodeInner)).Mean;
    trkChargeRMS = (pTrTrack->GetQ_all(tofBeta, trkFitCodeInner)).RMS;
    std::vector<like_t> like;
    if( fabs(trkRigidityInner) < 3 && tofBeta < 1)
    {
      trkZ = pTrTrack->GetZ(like, tofBeta, trkFitCodeL1Inner);
      trkInnerZ = pTrTrack->GetInnerZ(like, tofBeta, trkFitCodeInner);
      for(int k = 0; k < 9; k++) trkLayerJZ[k] = pTrTrack->GetLayerJZ(like, k, tofBeta, trkFitCodeInner);
    }
    else
    {
      trkZ = pTrTrack->GetZ(like);
      trkInnerZ = pTrTrack->GetInnerZ(like);
      for(int k = 0; k < 9; k++) trkLayerJZ[k] = pTrTrack->GetLayerJZ(like, k, tofBeta, trkFitCodeInner);
    }
    trkInnerCharge = (pTrTrack->GetInnerQ_all(tofBeta, trkFitCodeInner)).Mean;
    trkInnerChargeRMS = (pTrTrack->GetInnerQ_all(tofBeta, trkFitCodeInner)).RMS;

    trkHasExtLayers = pTrTrack->HasExtLayers();

    //
    //
    // ECAL shower Section
    //
    //

    // Save ECAL shower related variables in case of EcalShowerR object exists.
    isEcalAvailable = 0;
    if( pEcalShower )
    {
      isEcalAvailable = 1;
      showerEnergyCorrected = pEcalShower->GetCorrectedEnergy(2, 2);
      if( showerEnergyCorrected < 0.5 ) continue;
      showerEnergyD = (pEcalShower->EnergyD)/1000.;
      for(int i = 0; i < pEcalShower->NEcal2DCluster(); ++i)
      {
        for(int j = 0; j < pEcalShower->pEcal2DCluster(i)->NEcalCluster(); ++j)
        {
          if (pEcalShower->pEcal2DCluster(i)->pEcalCluster(j)->Edep > 0)
            showerEnergyDL[pEcalShower->pEcal2DCluster(i)->pEcalCluster(j)->Plane] += pEcalShower->pEcal2DCluster(i)->pEcalCluster(j)->Edep;
          else
            showerEnergyDL[pEcalShower->pEcal2DCluster(i)->pEcalCluster(j)->Plane] = -1;
        }
      }
      showerEnergyE = pEcalShower->EnergyE;
      showerBDT     = pEcalShower->GetEcalBDT();
      showerCofG[0] = pEcalShower->CofG[0];
      showerCofG[1] = pEcalShower->CofG[1];
      showerCofG[2] = pEcalShower->CofG[2];
      AMSPoint trackPoint;
      AMSDir   trackDirection;
      pTrTrack->Interpolate(pEcalShower->CofG[2], trackPoint, trackDirection, id_inner);
      pEcalShower->NormaliseVariableLAPP();
      showerCofGDist = std::sqrt(std::pow(trackPoint[0] - showerCofG[0], 2) + std::pow( trackPoint[1] - showerCofG[1], 2));
      showerCofGdX   = std::fabs(trackPoint[0] - showerCofG[0]);
      showerCofGdY   = std::fabs(trackPoint[1] - showerCofG[1]);
    }

    // RICH Section
    
    pRichRing = pParticle->pRichRing();
    if( pRichRing )
    {
      isRichAvailable            = 1;
      richRebuild                =   (int)pRichRing->Rebuild();
      richIsGood                 =   (int)pRichRing->IsGood();
      richIsClean                =   (int)pRichRing->IsClean();
      richIsNaF                  =   (int)pRichRing->IsNaF();
      richTileIndex              =   (int)pRichRing->getTileIndex();
      richUsedHits               =   (int)pRichRing->getUsedHits();
      richRingWidth              = (float)pRichRing->RingWidth();
      richDistanceTileBorder     = (float)pRichRing->DistanceTileBorder();
      richChargeCorrections      =   (int)pRichRing->buildChargeCorrections();
      richPMTCorrectionsFailed   =   (int)pRichRing->PmtCorrectionsFailed();
      richWidth                  =        pRichRing->getWidth();
      richNHits                  =        pRichRing->getHits();
      richNPMTsOnRing            =        pRichRing->getPMTs();
      richBeta                   =        pRichRing->getBeta();
      richBetaError              =        pRichRing->getBetaError();
      richChargeSquared          =        pRichRing->getCharge2Estimate();
      richKolmogorovProbability  =        pRichRing->getProb();
      richPhotoelectrons         =        pRichRing->getPhotoElectrons();
      richExpectedPhotoelectrons =        pRichRing->getExpectedPhotoelectrons();
      const float* tmpRichTrackEmissionPoint = pRichRing->getTrackEmissionPoint();
      richTrackEmissionPoint[0]  = tmpRichTrackEmissionPoint[0];  // X
      richTrackEmissionPoint[1]  = tmpRichTrackEmissionPoint[1];  // Y
      richTrackEmissionPoint[2]  = tmpRichTrackEmissionPoint[2];  // Z
      richTrackEmissionPoint[3]  = tmpRichTrackEmissionPoint[3];  // Theta
      richTrackEmissionPoint[4]  = tmpRichTrackEmissionPoint[4];  // Phi

      AMSPoint trackPointRICH;
      AMSDir   trackDirAtRICH;
      richTheta                  =        pRichRing->getTrackTheta();
      richPhi                    =        pRichRing->getTrackPhi();
      richBetaConsistency        =        pRichRing->getBetaConsistency();
      richReflectedHits          =        pRichRing->getReflectedHits();
      richPMTChargeConsistency   =        pRichRing->getPMTChargeConsistency();
    }

    // TRD vertex cut
    pTrdTrack = pParticle->pTrdTrack();
    int trdNUsedHits = 0;
    int trdNUsedSegment = 0;
    if(!pTrdTrack)
    {
      trdNClusters      = -1;
      trdNUsedHits      = -1;
      trdNUsedSegment   = -1;
      trdNUsedHits      = -1;
      trdNTracks        = -1;
      trdTrackTheta     = -9.;
      trdTrackPhi       = -9.;
      trdTrackPattern   = -9;
      trdTrackCharge    = -9;
      for(int i = 0; i < 20; ++i) trdTrackEdepL[i] = -9.;
    }
    else
    {
      trdNClusters = pAMSEvent->nTrdCluster();
      trdNUsedSegment = pTrdTrack->NTrdSegment();
      for(int i = 0; i < trdNUsedSegment; i++)
      {
        pTrdSegment = pTrdTrack->pTrdSegment(i);
        trdNUsedHits += pTrdSegment->NTrdCluster();
      }
      trdNUnusedHits = trdNClusters - trdNUsedHits;

      trdNTracks      = pAMSEvent->nTrdTrack();
      trdTrackTheta   = pTrdTrack->Theta;
      trdTrackPhi     = pTrdTrack->Phi;
      trdTrackPattern = pTrdTrack->Pattern;
      trdTrackCharge  = pTrdTrack->Q;
      trdTrackChi2    = pTrdTrack->Chi2;
      trdTrackMeanDepositedEnergy = 0.;
      int ntrdlayers = 0;
      for(int i = 0; i < pTrdTrack->NTrdSegment(); i++)
      {
        for(int j = 0; j < pTrdTrack->pTrdSegment(i)->NTrdCluster(); j++)
        {
          int trdLayer  = pTrdTrack->pTrdSegment(i)->pTrdCluster(j)->Layer;
          float trdEdep = pTrdTrack->pTrdSegment(i)->pTrdCluster(j)->EDep;
          trdTrackEdepL[trdLayer] = trdEdep;
          trdTrackTotalDepositedEnergy += trdEdep;
          ntrdlayers++;
        }
      }
      trdTrackMeanDepositedEnergy = trdTrackTotalDepositedEnergy/(float)ntrdlayers;

    }

    // Below of this for TrdK

    trdKNRawHits = pAMSEvent->NTrdRawHit();
    if(trdKNRawHits > 0)
    {
      double trdKLikelihoodRatio[3] = {-1., -1., -1.};
      double trdKLikelihood[3]      = {-1., -1., -1.};

      pTrdKCluster = new TrdKCluster(pAMSEvent, pTrTrack, id_inner);
      if(!pTrdKCluster) continue;

      trdKIsReadAlignmentOK = pTrdKCluster->IsReadAlignmentOK;
      trdKIsReadCalibOK     = pTrdKCluster->IsReadCalibOK;
      trdKIsValid           = pTrdKCluster->GetLikelihoodRatio_TrTrack(15, trdKLikelihoodRatio, trdKLikelihood, trdKNHits, trdKTotalPathLength, trdKTotalAmp, -1, 0);
      if(trdKIsValid != 0 && trdKNHits != 0)
      {
        pTrdKCluster->CalculateTRDCharge();
        trdKCharge                              = pTrdKCluster->GetTRDCharge();
        trdKChargeError                         = pTrdKCluster->GetTRDChargeError();
        trdKNUsedHitsForCharge                  = pTrdKCluster->GetQTRDHitCollectionNuclei().size();
        trdKElectronToProtonLogLikelihoodRatio  = trdKLikelihoodRatio[0];
        trdKHeliumToElectronLogLikelihoodRatio  = trdKLikelihoodRatio[1];
        trdKHeliumToProtonLogLikelihoodRatio    = trdKLikelihoodRatio[2];
        trdKElectronLikelihood                  = trdKLikelihood[0];
        trdKProtonLikelihood                    = trdKLikelihood[1];
        trdKHeliumLikelihood                    = trdKLikelihood[2];

        AMSPoint trExtraP0;
        AMSDir   trExtraDir;
        pTrdKCluster->GetTrTrackExtrapolation(trExtraP0, trExtraDir);

        for(int l = 0; l < pTrdKCluster->NHits(); l++)
        {
          TrdKHit* pTrdKHit = pTrdKCluster->GetHit(l);
          int tmpL = 0;
          tmpL = pTrdKHit->TRDHit_Layer;
          float tmpAmp = 0;
          tmpAmp = pTrdKHit->TRDHit_Amp;
          float tmpPathLength = 0;
          tmpPathLength = pTrdKHit->Tube_Track_3DLength(&trExtraP0, &trExtraDir);

          if( tmpAmp < 15 || tmpPathLength <= 0 ) continue;
          if( pTrdKHit->IsAligned == 0 || pTrdKHit->IsCalibrated == 0 ) continue;
          trdKAmpLayer[tmpL] += tmpAmp;
        }
      }
      delete pTrdKCluster;
    }
    outTree->Fill();
  } // End of for loop

  std::cout << "Finalizing.." << std::endl << std::endl;
  pOutputTFile->cd("/");
  hEvtCounter->Write();
  outTree->Write();
  delete outTree;
  pOutputTFile->Write();
  pOutputTFile->Close();
  delete pOutputTFile;
  std::cout << "The program exited successfully" << std::endl << std::endl;
  std::cout << nEvents << " events were analyzed and " << nSelected << " events are selected" << std::endl;
  std::cout << "Output is saved as: " << outputFileName << std::endl;
  delete [] outputFileName;

  clock.Stop();
  std::cout << "TStopwatch report: " << std::endl;
  std::cout << "Current time: " << CurrentDateTime() << std::endl;
  clock.Print();

  return 0;
}
