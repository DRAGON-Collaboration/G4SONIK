#include <memory>
#include "KinManager.hh"

#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TSpline.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TF2.h>

TGraph* g_theta;
TGraph* g_phi;
TGraph* g_E_b;
TGraph* g_E_b1;
TGraph* g_E_b2;
TGraph* g_E_B;

std::auto_ptr<TSpline3> s_theta(0) ;
std::auto_ptr<TSpline3> s_phi(0) ;
std::auto_ptr<TSpline3> s_E_b(0) ;
std::auto_ptr<TSpline3> s_E_B(0) ;

double m_a, m_A;

KinManager::KinManager(double n_a, double n_A) {

  //graphs relate COM angle to lab energy and lab angle
  g_theta = new TGraph();//Ejectile
  g_phi = new TGraph();//Recoil
  g_E_b = new TGraph();//Ejectile
  g_E_b1 = new TGraph();//Ejectile+
  g_E_b2 = new TGraph();//Ejectile-
  g_E_B = new TGraph();//Recoil

  g_theta->SetName("theta_ejectile");
  g_phi->SetName("phi_recoil");
  g_E_b->SetName("E_ejectile");
  g_E_b1->SetName("E_ejectile+");
  g_E_b2->SetName("E_ejectile-");
  g_E_B->SetName("E_recoil");

  m_a = n_a;//Beam
  m_A = n_A;//Target

}

//----------------------------------------------------------------------------------
//fills graphs relating COM angle to lab energy and lab angle
//----------------------------------------------------------------------------------

//Ejc_Theta_lab = Theta
//Rcl_Theta_lab = Phi
//Ejc_E_lab = E_b
//Rcl_E_lab = E_B

void KinManager::SetEbeam(double E_a) {

  double const pi = 3.14159;
  double E_b1, E_b2, min;

  double diff[1800];

  gamma = m_a/m_A;

  double s = E_a*(m_A-m_a)/(m_a+m_A);

  double theta_cmm = 0;

  for (int i = 0; i<1800; i+=1) {

    theta_cmm+=0.1;
    theta_cm = theta_cmm/180.*pi;

    theta = TMath::ACos( (gamma+cos(theta_cm)) / pow(1+gamma*gamma+2*gamma*cos(theta_cm),0.5) );

    double r = pow(m_a*m_a*E_a,0.5)*cos(theta)/(m_a+m_A);

    E_b1 = pow(r+pow(r*r+s,0.5),2);
    E_b2 = pow(r-pow(r*r+s,0.5),2);

    diff[i] = E_b1-E_b2;

  }

  min = diff[0];
  theta_cmm = 0;

  for (int i = 0; i<1800; i+=1) {
    theta_cmm += 0.1;
    if (diff[i] < min) {
      min = diff[i];
      change = theta_cmm;
    }
  }

  theta_cmm =0;

  for (int i = 0; i<1800; i+=1) {

    theta_cmm+=0.1;
    theta_cm = theta_cmm/180.*pi;
    theta =  TMath::ACos( (gamma+cos(theta_cm)) / pow( 1+gamma*gamma+2*gamma*cos(theta_cm),0.5 ) );

    double r = pow(m_a*m_a*E_a,0.5)*cos(theta)/(m_a+m_A);

    if (theta_cmm<=change) E_b = pow(r+pow(r*r+s,0.5),2);
    if (theta_cmm>change) E_b = pow(r-pow(r*r+s,0.5),2);

    E_B = E_a-E_b;

    double phi = TMath::ASin(pow(2*m_a*E_b,0.5)*sin(theta)/pow(2*m_A*E_B,0.5));

    //    double phi = (3.14159-theta_cm)/2.;
    /*
      theta = pow(2*m_a*E_a,0.5) - pow(2*m_A*E_B,0.5)*cos(phi);
      theta = theta/pow(2*m_a*E_b,0.5);
      theta = TMath::ACos(theta);
    */

    g_theta->SetPoint(i,theta_cm,theta);
    g_phi->SetPoint(i,theta_cm,phi);
    g_E_b->SetPoint(i,theta_cm,E_b);
    g_E_B->SetPoint(i,theta_cm,E_B);

    E_b = pow(r+pow(r*r+s,0.5),2);
    g_E_b1->SetPoint(i,theta_cm,E_b);

    E_b = pow(r-pow(r*r+s,0.5),2);
    g_E_b2->SetPoint(i,theta_cm,E_b);

  }


  s_phi.reset(new TSpline3("g_phi",g_phi));
  s_E_b.reset(new TSpline3("g_E_b",g_E_b));
  s_E_B.reset(new TSpline3("s_E_B",g_E_B));
  s_theta.reset(new TSpline3("g_theta",g_theta));

  TFile File("kin.root","RECREATE");

  g_theta->SetTitle("theta;theta_{cm} (rad);theta_{lab} (rad)");
  g_phi->SetTitle("phi;theta_{cm} (rad);phi_{lab} (rad)");
  g_E_b->SetTitle("E_{b};theta_{cm} (rad);E_{b}_{lab} (MeV)");
  g_E_B->SetTitle("E_{B};theta_{cm} (rad);E_{B}_{lab} (MeV)");

  g_theta->Write();
  g_phi->Write();
  g_E_b->Write();
  g_E_B->Write();

  TMultiGraph *mg = new TMultiGraph();

  mg->Add(g_E_b1);
  mg->Add(g_E_b2);
  mg->Add(g_E_b);
  mg->Write();

}


double KinManager::GetTheta (double theta_cm) {

  return (s_theta->Eval(theta_cm));

}

double KinManager::GetPhi (double theta_cm) {

  return (s_phi->Eval(theta_cm));

}

double KinManager::GetE_b (double theta_cm) {

  return (s_E_b->Eval(theta_cm));

}

double KinManager::GetE_B (double theta_cm) {

  return (s_E_B->Eval(theta_cm));

}
