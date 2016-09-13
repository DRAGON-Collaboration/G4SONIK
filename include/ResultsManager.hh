#ifndef ResultsManager_h
#define ResultsManager_h 1

#include "InputManager.hh"

#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

class ResultsManager {

  public:
  ResultsManager(InputManager* aInMgr);
 ~ResultsManager();
  void CreateTable();

  private:
  TFile *f;  
  InputManager *InMgr;

};

#endif
