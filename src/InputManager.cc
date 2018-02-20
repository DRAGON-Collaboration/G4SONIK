#include "InputManager.hh"

InputManager::InputManager() {

}

InputManager::~InputManager() {

}

void InputManager::ReadFile(const char* filename) {

  ifstream ifs(filename);
  if(!ifs.good()) {
    cerr << "Invalid config file: \"" << filename << "\"\n\n";
    exit(1);
  }
  string line;
  int nlines=0;

  cout << endl;

  while(1) {

    nlines +=1;
    getline(ifs, line);
    if(!ifs.good()) break;
    line = line.substr(0, line.find("#")); // # = comment

    for (int i=0; i< line.size(); ++i) { // tab -> space
      if (line[i] == '\t') {
        line[i]  = ' ';
      }
    }

    for(int i=line.size()-1; i >0; --i) { // get rid of duplicate spaces
      if(line[i] == ' ' && line[i-1] == ' ')
        line.erase(i,1);
    }
    if(line.size() <= 1) continue; // empty line
    if(line.find(" ") > line.size()) { // no break for option to be specified
      cerr << "Invalid input file at line " << nlines << "n\n";
      exit(1);
    }

    string s0 = line.substr(0, line.find(" "));//key
    string s1 = line.substr(line.find(" ")+1);//value

    config_var[s0] = s1;//keys value to map
    cout << s0 << "\t\t" << s1 << endl;

  }

  cout << endl;

  ifs.close();

}
