#!/bin/bash

# BASH Shell Script for running the SONIK scattering simulations
# Place config.dat (with desired parameters), SONIK.dat, and 3He4He(2).dat in the energy folder
# Script will correct SONIK.dat (col size, det angle, and theta range only) for each detector
# HOWEVER, theta range (SNtm, and SNtM) must be modified here for different energies (recommend using CatKIN)
# Ignore recoils in theta range for detector angles larger than 75deg
# Run the script from the energy folder (e.g. /local/astro/scat/SONIK_He3/0.5MeV)


declare -a col=("1" "2" "3" "5" "7")
declare -a deg=("22.5" "25" "35" "40" "45" "55" "60" "65" "75" "90" "120")

declare -a SNtm=("17" "20" "30" "35" "40" "39" "36" "32" "26" "18" "8") #Ejectile + Recoil for He3->He4 at 0.5MeV
declare -a SNtM=("74" "73" "64" "60" "55" "60" "65" "70" "80" "24" "11") 

#declare -a SNtm=("64" "63" "55" "51" "47" "39" "35" "32" "26" "18" "8") #Ejectile only for He3->He4 at 0.5MeV
#declare -a SNtM=("74" "73" "64" "60" "56" "47" "44" "40" "33" "24" "11")

#declare -a SNtm=("17.5" "20" "30" "35" "40" "50" "55" "60" "70" "85" "115") #Recoil only for He3->He4 at 0.5MeV
#declare -a SNtM=("27.5" "30" "40" "45" "50" "60" "65" "70" "80" "95" "125") 

for i in "${col[@]}";do

	mkdir "$i"mm_10k_pnA
	cd "$i"mm_10k_pnA

	l=0

	for j in "${deg[@]}";do
 
		mkdir "$j"deg
  		cp "/local/astro/scat/SONIK_He3/2.0MeV/3He(a,a)_abinitio.dat" "$j""deg/3He(a,a)_abinitio.dat"
		cp "/local/astro/scat/SONIK_He3/2.0MeV/config.dat" "$j""deg/config.dat"

        	k=1
		while read line;do
			if((k==3));then
				echo 'colsize		'"$i"'.'
			elif((k==4));then
				echo 'detAngle	'"$j"'.'
			elif((k==6));then
				echo 'theta_min	'${SNtm[l]}'.'
			elif((k==7));then
				echo 'theta_max	'${SNtM[l]}'.'
			else
				echo "$line"
			fi
			((k++))
		done < "/local/astro/scat/SONIK_He3/2.0MeV/SONIK.dat" > "$j""deg/SONIK.dat"

		cd "$j"deg
	
		scat config.dat SONIK.dat

		cd ..

		((l++))		

	done

	cd ..
done



