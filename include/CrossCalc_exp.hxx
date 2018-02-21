#include <iostream>
#include <fstream>
#include <TGraph.h>
#include <TMath.h>
#include <vector>

class ResonanceCalculator_exp {

private:

  TGraph2D* gcs;
  bool isGood; //< tells whether the constructor succeded in reading the cross-section data file

public:
  bool good() { return isGood; }

  /// The filename should be 3-collumn:
  /// cm energy (MeV)    cm angle (degrees)    cm cross section (mb/sr)
  ResonanceCalculator_exp(const char* filename) {

    gcs = new TGraph2D();

    std::ifstream inFile(filename);
    int nlines;

    double buffer[3];

    vector<double> E_cm(0);
    vector<double> theta_cm(0);
    vector<double> cs_cm(0);

    for (int i=0; i>=0; i++) {

      inFile >> buffer[0] >> buffer[1] >> buffer[2];

      if (!inFile.good() || buffer[0]==-1){
        nlines = i;
        break;
      }

      if (i==0) {
        E_cm.push_back(buffer[0]);
        theta_cm.push_back(buffer[1]);
        cs_cm.push_back(buffer[2]);
        continue;
      }

      for (int j=0; j<(int)E_cm.size(); j++) {

        if (E_cm.at(j) == buffer[0] && theta_cm.at(j) == buffer[1]) {
          cs_cm.at(j) = (buffer[2]+cs_cm.at(j))/2;
          break;
        }

        if (j==(int)E_cm.size()-1) {
          E_cm.push_back(buffer[0]);
          theta_cm.push_back(buffer[1]);
          cs_cm.push_back(buffer[2]);
        }

      }

    }

    for (int i=0; i<(int)E_cm.size(); i++) {
      theta_cm.at(i) = theta_cm.at(i)*3.14159/180.;
      gcs->SetPoint(i,E_cm.at(i),theta_cm.at(i),cs_cm.at(i));
    }

  }

  ~ResonanceCalculator_exp() {

  }

  Double_t Calculate(Double_t energy, Double_t angle) {
    return (gcs->Interpolate(energy,angle));
  }

};
