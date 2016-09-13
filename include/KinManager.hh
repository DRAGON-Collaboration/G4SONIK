#ifndef KinManager_h
#define KinManager_h 1

#include <string>

using namespace std;

class KinManager {
  public:
  
  KinManager(double m_a, double m_A);
  ~KinManager();

  void SetEbeam(double E_a);
  double GetTheta(double theta_cm); //Ejc_Theta_lab = Theta
  double GetPhi(double theta_cm);   //Rcl_Theta_lab = Phi
  double GetE_b(double theta_cm);   //Ejc_E_lab = E_b
  double GetE_B(double theta_cm);   //Rcl_E_lab = E_B

  private:

  double E_b, change, gamma, s, theta, theta_cm, E_B;

};

#endif
