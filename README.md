# SONIK Scattering Simulation

- To compile simulation use commands(in /local/astro/scat) { source env.sh } and the { make }
- To run use command { scat config.dat SONIK.dat }
- Run settings are input from the config.dat and SONIK.dat text files

- Output files are Events0.root, Table0.tsv (or whatever specified in config.dat), CS.root (generated cross-section and beam profiles), Kin.root (generated kinematics solutions), and currentEvent.rndm, currentRun.rndm (pseudoRNG seeds)

- To disable multiple scattering, comment out line 243, 248 in src/PhysicsList.cc { pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1); }

___
## Compile on your machine: 

To compile on your machine
1. source G4SONIKenv.sh (making sure the directory paths to ROOT and Geant4 simulation match that of your machine).
2. delete GNUmakefile, this is no longer needed.
3. enter the following command: cmake -DGeant4_DIR=/path/To/geant4.9.6 /path/To/G4SONIK
4. compile with make.
___

### Explained in the Appendix section of Charlie Akers' Thesis
- Commences by generating a beam histogram through the target, and using the input cross section data it calulcates the number of scattering events
- Creates a scattering distribution to get a weighted energy and angle(cm) for each generated scattering event
- Kinematics are calculated by evaluating a histogram for a given &theta;<sub>cm</sub>, this histogram is calculated at the beginning for a single energy(set to the beam energy at 20cm)
- The Cross Section and Kinematics histograms are saved in the files CS.root and Kin.root

- For each scattering event, Geant4 will generate a recoil particle and then an ejectile particle
- If a particle is detected it will fill the TotalDetector branch and the Recoil(or Ejectile) Detector Branch

- Estimated time remaining is for current run only
- Run time depends on number of scattering events, which depends on probability of scattering, which depends on theta range and energy
- Target isotope may only be H, <sup>3</sup>He, or <sup>4</sup>He (although this could be changed in CrossSectionManager.cc and DetectorConstruction.cc)
___

Source files are in directory src, corresponding headers in directory include. Some important source files are:

- scat.cc (in main directory)
  - Entry point of the application

- CrossCalc.hxx (in src directory)
	- Calculates the cross section for a given angle and energy from the input cross section data data
	- Different CrossCalc(1 or 2).hxx are included in the src directory for different format of input cross section data (format specified in source comments)
	- Input cross section is specified in the config.dat file
	- Default format is: energy(MeV)	angle(degrees)	cross section(mb/sr)

- CrossSectionManager.cc
	- Generates beam histogram and calculates number of scattering events 
	- Generates CS.root and calls the kinematics manager to generate Kin.root as well
	- GenerateEvent() called by PrimaryActionGenerator to generate events details (angles, energies, positions) using kinematics manager

- DetectorConstruction.cc
	- Constructs SONIK geometry for use by Geant4
	- This can be viewed by enabling the visualizer in config.dat (use indirect rendering { export LIBGL_ALWAYS_INDIRECT = 1 } and define +iglx when starting xTerm if problems are experienced)

- InputManager.cc
	- Reads config.dat and SONIK.dat variables into a string map
	- Can be called to copy individual variables to a declared variable

- KinManager.cc
	- Generates kinematics histogram for a given energy
	- Is called with a given theta(cm) to evaluate histograms and return desired value(energy, lab angle)

- PrimartyActionGenerator.cc
	- Calls CrossSectionManager to generate event position, energy, and direction and use G4ParticleGun to generate particles
	- Assigns values to Tree Manager to be filled if particle is detected

- TrackerSD.cc
	- Processes detector hits, and calls TreeManager to fill corresponding branches

- TreeManager.cc
	- ROOT file is created, filled, and written here

- ResultsManager.cc
	- Opens ROOT file at end of run, and creates a table with summarized results


<!-- Local Variables:  -->
<!-- mode: gfm  -->
<!-- End:  -->
