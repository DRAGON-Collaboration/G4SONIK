#include <map>
#include <iostream>
#include <fstream>
#include <TGraph.h>
#include <TMath.h>

/// A class to facilitate calculation of resonant cross sections using interpolation.
/// We are provided with cross section data at specific cm energies and angles, and want
/// to calculate cross sections in between those points. The method used here is to use two
/// spline interpolations, with the calculation handled by ROOT's TGraph::Eval() method.
class ResonanceCalculator {

private:
  /// These maps basically allow look-up of plots of either energy vs. cx (at fixed angles)
  /// or angle vs. cx (at fixed energies).
  std::map<Double_t, TGraph*> energyGraphs;
  std::map<Double_t, TGraph*> angleGraphs;
  bool isGood; //< tells whether the constructor succeded in reading the cross-section data file

public:
  bool good() {
    return isGood;
  }

  /// Constructor, read in a file and populate two sets of graphs:
  ///    1) energy vs. cross section at fixed angles
  ///    2) angle vs. cross section at fixed energies
  /// The filename should be 3-collumn:
  /// energy (MeV)    angle (degrees)    cross section (mb/sr)
  ResonanceCalculator(const char* filename) {
    std::ifstream ifs(filename);
    if(!ifs.good()) {
      std::cerr << "Error in <ResonanceCalculator::ResonanceCalculator>: Invalid input file " << filename <<"\n\n";
      isGood = false;
    }
    else {
      isGood = true;
      double energy, angle, cross;
      std::map<Double_t, TGraph*>::iterator it;
      while(1) {
	ifs >> energy >> angle >> cross;
//        cout << energy << "\t" << angle << "\t" << cross << endl;
	if(!ifs.good()) break;

	it = energyGraphs.find(energy);
	if(it == energyGraphs.end()) {
	  energyGraphs.insert(std::make_pair(energy, new TGraph()));
	  energyGraphs[energy]->SetPoint(0, angle, cross);
	} else {
	  int nnn = energyGraphs[energy]->GetN();
	  energyGraphs[energy]->SetPoint(nnn, angle, cross);
	}

	it = angleGraphs.find(angle);
	if(it == angleGraphs.end()) {
	  angleGraphs.insert(std::make_pair(angle, new TGraph()));
	  angleGraphs[angle]->SetPoint(0, energy, cross);
	} else {
	  int nnn = angleGraphs[angle]->GetN();
	  angleGraphs[angle]->SetPoint(nnn, energy, cross);
	}
      }
      ifs.close();
    }
  }

  /// Destructor, free dynamically allocated memory in maps
  ~ResonanceCalculator() {
    map<Double_t, TGraph*>::iterator it;
    it = energyGraphs.begin();
    while(it!=energyGraphs.end()) {
      delete it->second;
      ++it;
    }
    it = angleGraphs.begin();
    while(it != angleGraphs.end()) {
      delete it->second;
      ++it;
    }
  }

  /// The actual calculation function.
  /// Method:
  ///   1) For each fixed-angle graph (angleGraphs), use spline interpolation to get the
  ///      cross section at the input energy.  Take these values and populate a new graph
  ///      of angle vs. cross section.
  ///   2) Using spline interpolation on the graph created in (1), calculate the cross section
  ///      at the input angle.
  Double_t Calculate(Double_t energy, Double_t angle) {
    // populate a TGraph with (angle, cx) pairs @ given energy
    TGraph graph;
    std::map<Double_t, TGraph*>::iterator it;
    for(it = angleGraphs.begin(); it != angleGraphs.end(); ++it) {
      graph.SetPoint(graph.GetN(), it->first, it->second->Eval(energy, 0, "S"));  /// spline interp
    }

    return graph.Eval(angle, 0, "S");
  }



};
