#include "CrossSectionManager.hh"
#include <typeinfo>

using namespace std;

#include "CrossCalc.hxx"
#include "CrossCalc_exp.hxx"
#ifdef __MAKECINT__
#pragma link C++ class ResonanceCalculator+;
#endif

class ResonanceCalculator;

// Constants
namespace {
  const Double_t amuSI = 1.6726e-27;//kg
  //  const Double_t amu = 931.494;
  const Double_t clight = 2.99792458E+10; // cm/sec
  const Double_t avo = 6.022e+23; // molecules per mole
  const Double_t Coul = 1.602176E-19; // Coulombs
  const Double_t perm = 8.854188e-12; //SI
  // Central temperature (Centigrade)
  //  const Double_t temp = 21.0;
  // Gas constant (L.atm/mol.K)
  //  const Double_t R = 0.08206;
  // Loschmidt constant
  const Double_t Losch = 2.6867805e19; //cm^-3
  //polar angle of beam in zx plane (rad)
  const Double_t beam_thetax = 0. * 1.e-3;
  //polar angle of beam in zy plane (rad)
  const Double_t beam_thetay = 0. * 1.e-3;
}

CrossSectionManager::CrossSectionManager(InputManager* aInMgr) {

  InMgr = aInMgr;
  CSfile = new TFile("CS.root","RECREATE");

  Event = false;
  GenerateCrossSection();

}

CrossSectionManager::~CrossSectionManager() {

  CSfile->Write();
  CSfile->Close();

}

//---------------------------------------------------------------------------------
//Calculates scattering cross section at specific energy & angle across whole range
//---------------------------------------------------------------------------------

void CrossSectionManager::GenerateCrossSection() {

  EventCount=0;

  double mcFac;

  InMgr->GetVariable("mcFactor",mcFac);//pnA of beam on target
  InMgr->GetVariable("beamA",a[BEAM]);
  InMgr->GetVariable("beamZ",z[BEAM]);
  InMgr->GetVariable("targetA",a[TARGET]);
  InMgr->GetVariable("targetZ",z[TARGET]);

  KinMgr = new KinManager(a[BEAM],a[TARGET]);

  double elab, eloss_tot;//incident lab beam energy & total loss through target
  InMgr->GetVariable("beamE",elab);
  InMgr->GetVariable("eloss",eloss_tot);

  Double_t elabmin = elab-eloss_tot;

  string targetType;
  double P/*torr*/, temp/*°C*/,Nvol/*atom/cm3*/;
  int nu;
  InMgr->GetVariable("targetTemp",temp);
  InMgr->GetVariable("targetP",P);

  if      (z[TARGET] == 2 && a[TARGET] == 4) targetType = "4He";
  else if (z[TARGET] == 2 && a[TARGET] == 3) targetType = "3He";
  else if (z[TARGET] == 1 && a[TARGET] == 1) targetType = "H";
  else {
    cerr << "Error in <CrossSectionManager::CrossSectionManager>: Invalid target type, only 4He, 3He & H allowed" << endl;
    exit(1);
  }

  if (targetType == "H") nu = 2;
  else nu = 1;

  Nvol = nu*Losch*(P / 760)*(273.15 / (temp+273.15) );

  string GeomType;
  InMgr->GetVariable("GeomType",GeomType);

  if (GeomType == "Regular") TargetLength = 12.3;//cm
  if (GeomType == "SONIK" ){
    TargetLength = 37;//cm
    // if (targetType == "H") TargetLength = 21.2776;//cm
    // else TargetLength = 22.0043;//cm
  }

  if (GeomType == "Regular") TargetOffset = -TargetLength/2.;//cm
  if (GeomType == "SONIK" ) TargetOffset = 7.125;//cm

  Double_t Narea = Nvol * TargetLength;//Total atoms per cm**2

  int nLayers;
  InMgr->GetVariable("nLayers",nLayers);

  //Trace beam through from entrance collimator & define energy for
  //each target layer.
  vector <Double_t> elay(nLayers);
  Double_t ecurr = elab;
  Double_t eloss = (elab-elabmin)/nLayers;
  for (Int_t i=0; i<nLayers; i++) {
    ecurr = ecurr - eloss;
    elay.at(i) = ecurr;
  }

  double sigma_e;//beam energy sd
  InMgr->GetVariable("sigma",sigma_e);

  sigma_e = sigma_e*elab;// Beam energy standard deviation (MeV)
  double e_min = elay.at(nLayers-1)-5.*sigma_e;//min and max beam energy limits in target
  double e_max = elab+5.*sigma_e;

  InMgr->GetVariable("z_min", z_min);//beam axis limits (cm)
  InMgr->GetVariable("z_max", z_max);

  //converts z position limits from detector geomertry to that relative to start of target
  z_min -= TargetOffset;//cm
  z_max -= TargetOffset;

  int zbin_max = int(z_max/TargetLength*nLayers+1.);//min and max z bins in z limits
  int zbin_min = int(z_min/TargetLength*nLayers-1.);

  //recalculates z limits to be consistent with bin limits
  z_min = double(zbin_min)/double(nLayers)*TargetLength;
  z_max = double(zbin_max)/double(nLayers)*TargetLength;

  double e_min_lim = elay.at(zbin_max)-5.*sigma_e;//min and max beam energy limits in z limits
  double e_max_lim = elay.at(zbin_min)+5.*sigma_e;

  int n_ebin;//number of energy bins
  InMgr->GetVariable("n_ebin",n_ebin);

  ebin_min = int((e_min_lim-e_min)/(e_max-e_min)*n_ebin-1);//min and max energy bins in limits
  int ebin_max = int((e_max_lim-e_min)/(e_max-e_min)*n_ebin+1);
  int n_ebin_lim = ebin_max-ebin_min;// number of energy bins in limits

  //recalculates e limits to be consistent with bin limits
  e_min_lim = (e_max-e_min)/double(n_ebin)*double(ebin_min)+e_min;
  e_max_lim = (e_max-e_min)/double(n_ebin)*double(ebin_max)+e_min;

  //Creates histo detailing beam energy pdf through target length
  hE = TH2F("E vs Z","E vs Z",nLayers,0.,TargetLength,n_ebin,e_min,e_max);
  hE.SetOption("COLZ");

  for (Int_t j=0; j<nLayers; j++) { // Loops through target layers

    pdone = j*100/nLayers;
    cout << '\r' << "Generating beam energy histogram  " << pdone << "%";

    double z = hE.ProjectionX()->GetBinCenter(j);
    double strag = 0;//(2.039*pow(z+TargetLength/2.,0.4997))/1000.;// Straggling in MeV (taken from program bohr)
    double sigma_tot = pow(sigma_e*sigma_e+strag*strag,0.5);

    for (int ii=0;ii<1e5;ii++) {
      hE.Fill(z,gRandom->Gaus(elay.at(j),sigma_tot));
    }

  }

  cout << '\r' << "Generating beam energy histogram  " << "100%" << endl;

  //Creates histo of beam energy pdf
  TH1F hEpdf = TH1F("E prob","E prob",n_ebin,e_min,e_max);

  for (int j=0; j<n_ebin; j++) {

    pdone = j*100/nLayers;
    cout << '\r' << "Generating beam energy pdf  " << pdone << "%";

    int n_ecount = hE.ProjectionX("hEprojX",j,j)->Integral();
    hEpdf.SetBinContent(j,n_ecount);

  }

  cout << '\r' << "Generating beam energy pdf  " << "100%" << endl;

  //Normalized pdf to 1
  double norm = hEpdf.Integral();
  for (int j=0; j<n_ebin; j++) {
    hEpdf.SetBinContent(j, hEpdf.GetBinContent(j)/norm);
  }

  double theta_min, theta_max;//limits on lab angle
  int n_theta, n_phi;//number of bins in angle range

  InMgr->GetVariable("theta_max", theta_max);//deg
  InMgr->GetVariable("theta_min", theta_min);

  InMgr->GetVariable("phi_max", phi_max);//deg
  InMgr->GetVariable("phi_min", phi_min);

  InMgr->GetVariable("n_theta", n_theta);//number of angle bins
  InMgr->GetVariable("n_phi", n_phi);

  Double_t theta_mincm = TMath::Pi() - (2.*theta_max*TMath::Pi()/180.);//convert to cm
  Double_t theta_maxcm = TMath::Pi() - (2.*theta_min*TMath::Pi()/180.);
  phi_min              = phi_min*TMath::Pi()/180.;//rad
  phi_max              = phi_max*TMath::Pi()/180.;
  Double_t dtheta      = (theta_maxcm-theta_mincm)/n_theta;
  Double_t dphi        = (phi_max-phi_min)/n_phi;

  if (theta_min < -90 || theta_max > 90) {
    cerr << "Illegal theta range chosen. Probabilty of scattering < 0" << endl;
    exit(1);
  }

  //Creates histo of cs for every E & theta within z limits
  hEtheta = TH2F("E vs theta","E vs theta",n_ebin_lim,e_min_lim,e_max_lim,n_theta,theta_mincm,theta_maxcm);
  hEtheta.SetOption("COLZ");

  //Creates histo of scattering probability for every beam energy within z limits
  hEscat = TH1F("E vs scat prob","E vs scat prob",n_ebin_lim,e_min_lim,e_max_lim);

  string cx_type;
  InMgr->GetVariable("cx_type", cx_type);

  if (cx_type != "rutherford") InMgr->GetVariable("cx_file",cxFile);//cs input file

  pdone=0;

  //Calculates absolute cs for given energy and theta bin, used to fill hEtheta
  for (int j=0; j<=n_ebin_lim; j++) {//Energy bin loop

    Double_t xadd;//cs for given energy and theta bin

    pdone_prev = pdone;
    pdone = j*100/n_ebin_lim;

    ecurr = hEpdf.GetBinCenter(j+ebin_min);
    Double_t ecm = ecurr*a[TARGET]/(a[BEAM]+a[TARGET]); //  MeV

    for (int i=0; i<n_theta; i++) {//theta bin loop

      if (pdone%1 == 0 && pdone>pdone_prev) cout << '\r' << "Generating cross sections  " << pdone << "%";

      Double_t thetacm = i*dtheta + theta_mincm;

      //Caluclates absolute cs in this energy & angle bin using user defined method
      if      (cx_type=="resonance")     xadd = res(ecm,thetacm)* dtheta* TMath::Sin(thetacm)*dphi;//fm2, 10mb
      else if (cx_type=="resonance_exp") xadd = res_exp(ecm,thetacm)* dtheta* TMath::Sin(thetacm)*dphi;//fm2, 10mb
      else if (cx_type=="rutherford") {
        xadd = ruth(z[BEAM],z[TARGET],ecm,thetacm)* dtheta* TMath::Sin(thetacm)*dphi;//fm2, 10mb
        //        cout << "\t" << ruth(z[BEAM],z[TARGET],ecm,thetacm) << endl;
      }
      else {
        cerr << "Error in <CrossSectionManager::CrossSectionManager>: Invalid cross section type" << endl;
        exit(1);
      }

      xadd=xadd*n_phi;//phi integral
      hEtheta.SetBinContent(j,i,xadd);//fills histogram with cs data

    }

    //Calculates scattering probability for given lab energy bin based on cs
    double dcsdE = hEtheta.ProjectionY("hEthetaprojY",j,j)->Integral();//cs at this energy (10mb)
    //number of expeted scattering events from one incident beam ion
    double pscat = dcsdE * Narea * hEpdf.GetBinContent(j+ebin_min) * 1.e-26;
    hEscat.SetBinContent(j,pscat);

  }

  time = clock(); //for calculating ETA
  prev_count = 0;

  n_scat = hEscat.Integral()*6.25e9*mcFac;//total number of scattering events in limits

  gRandom->SetSeed(0);

  //cout << "\n" << elab-(((z_max+z_min)/2)*eloss_tot/TargetLength) << endl; //beam energy at observation point, used in kinematics calculations
  KinMgr->SetEbeam(elab-(((z_max+z_min)/2)*eloss_tot/TargetLength));

}

//---------------------------------------------------------------------------------
//Generates position, angle and energy for one scattering event
//---------------------------------------------------------------------------------

bool CrossSectionManager::GenerateEvent() {

  EventCount+= 1;

  pdone_prev = pdone;
  pdone = int (100 * (EventCount/n_scat));


  if (EventCount%10000 == 0)
    {
      if ((EventCount+3)%100000 == 0)
        {
          time = clock();
          prev_count = EventCount;
        }
      cout << "\n" << "\033[J" << "Estimated time to completion: "<< (int)((clock()-time)*(n_scat-EventCount)/((EventCount-prev_count)*CLOCKS_PER_SEC*60)) << " minutes" << "\033[F";
    }

  if (pdone%10 == 0 && pdone>pdone_prev) cout << pdone << "%" << "\n" << "\033[J";
  if (EventCount == int(n_scat)) cout << "100%" << "\n" << "\033[J";

  //generates random lab energy weighted to scattering pdf
  double E_lab = hEscat.GetRandom(); //  MeV

  Beam_E_cm = E_lab * a[TARGET]/(a[BEAM]+a[TARGET]);//converts to cm

  int Ebin = hEscat.FindBin(E_lab);//energy bin number

  //finds z (beam axis) position
  pos[2] = hE.ProjectionX("Z dist",Ebin+ebin_min,Ebin+ebin_min)->GetRandom();

  if (pos[2]<z_min || pos[2]>z_max) {//discards event if outside z limits
    Event = false;
    return false;
  }

  //generates random cm angle weighted to angular disribution at this energy
  theta_cm = hEtheta.ProjectionY("Ang dist",Ebin,Ebin)->GetRandom();

  double sigx, sigy;//s.d if beam in x & y axis (mm)
  InMgr->GetVariable("sig_x",sigx);
  InMgr->GetVariable("sig_y",sigy);

  double x_offset, y_offset;
  InMgr->GetVariable("off_x",x_offset);
  InMgr->GetVariable("off_y",y_offset);

  //Generates random x & y position weighted to gaussian
  pos[0] = gRandom->Gaus(x_offset/10.,sigx/10.) + (pos[2]+TargetLength/2.)*TMath::Tan(beam_thetax);
  pos[1] = gRandom->Gaus(y_offset/10.,sigy/10.) + (pos[2]+TargetLength/2.)*TMath::Tan(beam_thetay);

  //Generate random phi direction between phimin and phimax
  phi = gRandom->Uniform(phi_min,phi_max);

  //finds coresponding lab angle and energy for recoil
  Rcl_theta_lab = KinMgr->GetPhi(theta_cm);//recoil
  Rcl_E_lab = KinMgr->GetE_B(theta_cm);

  //finds corresponding lab angle and energy for beam particle
  Ejc_theta_lab = KinMgr->GetTheta(theta_cm);//ejectile
  Ejc_E_lab = KinMgr->GetE_b(theta_cm);

  Event = true;

  return true;

}

//---------------------------------------------------------------------------------
// Definition of cross section calculation functions.
//---------------------------------------------------------------------------------

// Calculates a resonant cross section using R matrix fit input
double CrossSectionManager::res(Double_t ecm, Double_t theta_cm) {

  Double_t thetaDeg = theta_cm * 180. / TMath::Pi();

  static ResonanceCalculator rc(cxFile.c_str());

  if(rc.good())
    return rc.Calculate(ecm, thetaDeg) / 10.;//fm^2/sr
  else
    exit(1);

}

// Calculates a resonant cross section using experimental data
double CrossSectionManager::res_exp(Double_t ecm, Double_t theta_cm) {

  Double_t thetaDeg = theta_cm * 180. / TMath::Pi();
  static ResonanceCalculator_exp rc_exp(cxFile.c_str());

  if(rc_exp.good())
    return rc_exp.Calculate(ecm, thetaDeg) / 10.; // n.b. divide by 10 to convert mb to fm^2
  else
    exit(1);

}

//Rutherford cross-section
double CrossSectionManager::ruth(Double_t z_a,Double_t z_b, Double_t ecm, Double_t theta_cm) {
  //  cout << z_a << "\t" << z_b << "\t" << ecm << "\t" << theta_cm << endl;
  double dsig = 1.296*pow(z_a*z_b/ecm,2)*TMath::Power(TMath::Sin(theta_cm/2.),-4)/10.;//fm^2/sr
  return dsig;
}

double CrossSectionManager::GetTargetLength() {
  return TargetLength;
}

double CrossSectionManager::GetTargetOffset() {
  return TargetOffset;
}

bool CrossSectionManager::GetEvent() {
  return Event;
}

double CrossSectionManager::GetN_scat() {
  return(n_scat);
}

double CrossSectionManager::GetX() {
  return(pos[0]);
}

double CrossSectionManager::GetY() {
  return(pos[1]);
}

double CrossSectionManager::GetZ() {
  return(pos[2]);
}

double CrossSectionManager::GetRcl_Theta_lab() {
  return(Rcl_theta_lab);
}

double CrossSectionManager::GetPhi() {
  return(phi);
}

double CrossSectionManager::GetRcl_E_lab() {
  return(Rcl_E_lab);
}

double CrossSectionManager::GetBeam_E_cm() {
  return(Beam_E_cm);
}

double CrossSectionManager::GetIonA() {
  return(a[BEAM]+a[TARGET]);
}

double CrossSectionManager::GetIonZ() {
  return(z[BEAM]+z[TARGET]);
}

double CrossSectionManager::GetTheta_cm() {
  return(theta_cm);
}

double CrossSectionManager::GetEjc_E_lab() {
  return(Ejc_E_lab);
}

double CrossSectionManager::GetEjc_Theta_lab() {
  return(Ejc_theta_lab);
}

//////////////////////////////////////////////////////////////////////////////
