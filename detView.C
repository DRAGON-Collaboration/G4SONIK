#include <iostream>

#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TFile.h>
#include <TH1F.h>

#include <TMatrixD.h>
#include <TArrayD.h>

//Compile in ROOT with .L detView.C+, run with DetView()
//Returns effective length Z of the segment of target that each detector sees 
//Format: Collimator Diameter(mm) , Detector Angle , effective Z (cm)
//Written by Arthur Firmino, firmino@ualberta.ca when trying to calculate back the cross-section from the number of counts in each detector in the simulation using the formula 
//Double_t cx = (Double_t) conv * nevents/(SolidAngle*pnA*nVol*z);

/*
-ROOT C++ script for calculating based on geometry an effective z for each collimator, and angle setup
-Compile in ROOT with the command { .L detView.C+ } and run using { detView() }

-Integrates a weighing factor sigma(z) over z (between 19cm and 21cm) 
-sigma(z) is the double integral of the visibility ratio of the collimator multiplied by the relative intensity of the beam over the x,y plane...
-divided by the double intergal of the relative intensity of the beam over the x,y plane

-Visibility ratio is calculated by the following proccess
--Setup slit (4 vertices) and collimator (20 vertices and a center) with actual sizes and locations
--A camera position x,y,z is specified
--Entire setup is transformed to coordinates centered at the camera, with the z-axis now aligned between the camera and collimator center
--A camera projection matrix is used to transform the scene to a 2D image from the camera's perspective
---This can be viewed by commenting out the first return in the IonViewRatio function and running it individually after compiling, IonViewRatio(x,y,z,colRadius(cm),angle(deg))
--The 2D projection of the collimator is split into 400 pieces (20 divides radially, 20 divides angularly)
---Each piece is checked if within slit (within slit if any ray originating from piece intersects slit lines an odd number of times)
--Each visible piece area is summed up and divided by total area

-Also takes into account how much larger(or smaller) the solid angle is at (x,y,z) than at (0,0,20) 

-Relative intensity of the beam is a 2D Gaussian centered at x=0,y=0 with std deviations 1mm in both x and y (this is what is used in the simulation)

-Program will print Collimator Diameter, Detector Angle, and effective z for all setups (Runtime ~5-10mins, depending on zbins and xybins)
*/

double IonViewRatio(double cameraX, double cameraY, double cameraZ, double colRadius, double thetadeg)
{
	/////////////////////////////////////////////////////////////////////////// SET SLIT AND COLLIMATOR
	const double pi = 3.14159265;
	const int nvertices = 20;
	const double slitWidth =  0.2; //cm
	const double slitLength = 1.; //cm
	double theta = thetadeg * pi / 180.;
	double collimator[nvertices+1][3]; //x,y,z
	double collimatorTemp[nvertices+1][3]; //for operations		
	double phi;

	for (int j = 0; j < nvertices; j++)
	{
		phi = j * 2 * pi / nvertices;
		collimator[j][0] = colRadius * cos(phi);
		collimator[j][1] = colRadius * sin(phi);
		collimator[j][2] = 0.;
	}
	collimator[nvertices][0] = 0.;
	collimator[nvertices][1] = 0.;
	collimator[nvertices][2] = 0.;		

	double slit[4][3];
	double slitTemp[4][3];
	slit[0][0] = slitWidth / 2.;
	slit[0][1] = slitLength / 2.;
	slit[0][2] = 0.;
	slit[1][0] = -slitWidth / 2.;
	slit[1][1] = slitLength / 2.;
	slit[1][2] = 0.;
	slit[2][0] = -slitWidth / 2.;
	slit[2][1] = -slitLength / 2.;
	slit[2][2] = 0.;
	slit[3][0] = slitWidth / 2.;
	slit[3][1] = -slitLength / 2.;
	slit[3][2] = 0.;
	
	double slitTranslate[3] = {0, 0, 6.1588};
	double colTranslate[3] = {0, 0, 10.84 + 6.1588};

	double RotateToPos[3][3] =  { {cos(theta), 0., -sin(theta)}, {0., 1., 0.}, {sin(theta), 0., cos(theta)} }; //Rotation matrix
	
	for (int j = 0; j < nvertices+1; j++) { collimator[j][2] = collimator[j][2] + colTranslate[2]; } 
	for (int j = 0; j < 4; j++) { slit[j][2] = slit[j][2] + slitTranslate[2]; } 

	for (int j = 0; j < nvertices + 1; j++) //collimator vectors * Rotation Matrix ...Sorry I didn't use a linear algebra library
	{
		for (int k = 0; k < 3; k++)
		{
			collimatorTemp[j][k] = collimator[j][0]*RotateToPos[0][k] + collimator[j][1]*RotateToPos[1][k] + collimator[j][2]*RotateToPos[2][k];
		}
		for (int k = 0; k < 3; k++)
		{
			collimator[j][k] = collimatorTemp[j][k];
		}
	}

	for (int j = 0; j < 4; j++) //slit vectors * Rotation Matrix
	{
		for (int k = 0; k < 3; k++)
		{
			slitTemp[j][k] = slit[j][0]*RotateToPos[0][k] + slit[j][1]*RotateToPos[1][k] + slit[j][2]*RotateToPos[2][k];
		}
		for (int k = 0; k < 3; k++)
		{
			slit[j][k] = slitTemp[j][k];
		}
	}
	
	for (int j = 0; j < nvertices + 1; j++) { collimator[j][2] = collimator[j][2] + 20.; }
	for (int j = 0; j < 4; j++) { slit[j][2] = slit[j][2] + 20.; }

	/////////////////////////////////////////////////////////////////////////// SET CAMERA POSITION

	double cameraPos[3] = {cameraX, cameraY, cameraZ};

	/////////////////////////////////////////////////////////////////////////// TRANSFORM TO CAMERA CENTERED COORDINATES

	for (int j = 0; j < nvertices + 1; j++)
	{
		for (int k = 0; k < 3; k++)
		{
			collimator[j][k] = collimator[j][k] - cameraPos[k];
		}
	}

	for (int j = 0; j < 4; j++)
	{
		for (int k = 0; k < 3; k++)
		{
			slit[j][k] = slit[j][k] - cameraPos[k];
		}
	}

	double cameraPhi = -atan( collimator[nvertices][1] / collimator[nvertices][0] );
	//double cameraPhi = (pi / 2) - atan( collimator[nvertices][1] / collimator[nvertices][0] ); 
	double RotatePhi[3][3] =  { {cos(cameraPhi), sin(cameraPhi), 0.}, {-sin(cameraPhi), cos(cameraPhi), 0.}, {0., 0., 1.} }; //Rotation Matrix

	for (int j = 0; j < nvertices + 1; j++) //collimator vertices * Rotation Matrix
	{
		for (int k = 0; k < 3; k++)
		{
			collimatorTemp[j][k] = collimator[j][0]*RotatePhi[0][k] + collimator[j][1]*RotatePhi[1][k] + collimator[j][2]*RotatePhi[2][k];
		}
		for (int k = 0; k < 3; k++)
		{
			collimator[j][k] = collimatorTemp[j][k];
		}
	}

	for (int j = 0; j < 4; j++) //slit vectors * Rotation Matrix
	{
		for (int k = 0; k < 3; k++)
		{
			slitTemp[j][k] = slit[j][0]*RotatePhi[0][k] + slit[j][1]*RotatePhi[1][k] + slit[j][2]*RotatePhi[2][k];
		}
		for (int k = 0; k < 3; k++)
		{
			slit[j][k] = slitTemp[j][k];
		}
	}

	double cameraTheta = -atan( (collimator[nvertices][0]) / (collimator[nvertices][2]) ); 
	double RotateTheta[3][3] =  { {cos(cameraTheta), 0., -sin(cameraTheta)}, {0., 1., 0.}, {sin(cameraTheta), 0., cos(cameraTheta)} }; //Rotation Matrix

	for (int j = 0; j < nvertices + 1; j++) //collimator vectors * Rotation Matrix
	{
		for (int k = 0; k < 3; k++)
		{
			collimatorTemp[j][k] = collimator[j][0]*RotateTheta[0][k] + collimator[j][1]*RotateTheta[1][k] + collimator[j][2]*RotateTheta[2][k];
		}
		for (int k = 0; k < 3; k++)
		{
			collimator[j][k] = collimatorTemp[j][k];
		}
	}

	for (int j = 0; j < 4; j++) //slit vectors * Rotation Matrix
	{
		for (int k = 0; k < 3; k++)
		{
			slitTemp[j][k] = slit[j][0]*RotateTheta[0][k] + slit[j][1]*RotateTheta[1][k] + slit[j][2]*RotateTheta[2][k];
		}
		for (int k = 0; k < 3; k++)
		{
			slit[j][k] = slitTemp[j][k];
		}
	}

	/////////////////////////////////////////////////////////////////////////// PROJECT TO 2 DIMENSIONAL IMAGE

	double collimatorImageX[nvertices + 1];
	double collimatorImageY[nvertices + 1];
	double slitImageX[4];
	double slitImageY[4];
	
	double cameraFocus = 1.;

	for (int j = 0; j < nvertices+1; j++)
	{
		collimatorImageX[j] = (cameraFocus / collimator[j][2]) * collimator[j][0];
		collimatorImageY[j] = (cameraFocus / collimator[j][2]) * collimator[j][1];
	}

	for (int j = 0; j < 4; j++)
	{
		slitImageX[j] = (cameraFocus / slit[j][2]) * slit[j][0];
		slitImageY[j] = (cameraFocus / slit[j][2]) * slit[j][1];
	}

	/////////////////////////////////////////////////////////////////////////// FIND VISIBLE AREA RATIO

	int n = nvertices;
	double areaSeg, R, x, y;
	double colArea = 0.;
	double vis_colArea = 0.;

	double distance = sqrt(pow(collimator[nvertices][0],2)+pow(collimator[nvertices][1],2)+pow(collimator[nvertices][2],2));

	double p[4][2]; //p is first point on line, p+r is second
	double r[4][2];

	for (int j = 0; j < 4; j++)
	{
		p[j][0] = slitImageX[j];
		p[j][1] = slitImageY[j];
		r[j][0] = slitImageX[(j+1)%4]-p[j][0];
		r[j][1] = slitImageY[(j+1)%4]-p[j][1];
	}

	for (int g = 0; g < n; g++)
	{
		for (int h = 0; h < n; h++)
		{
			R = sqrt(pow(collimatorImageX[nvertices]-collimatorImageX[h],2) + pow(collimatorImageY[nvertices]-collimatorImageY[h],2));
			areaSeg = pi*(pow((g+1)*(R/n),2)-pow(g*R/n,2))/n;
			x = collimatorImageX[nvertices] + ((g+0.5)*(R/n)*cos(h*2*pi/n));
			y = collimatorImageY[nvertices] + ((g+0.5)*(R/n)*sin(h*2*pi/n));
			int m = 0;
			for (int j = 0; j < 4; j++)
			{		
				double time0 = (-p[j][0]*y + p[j][1]*x) / (r[j][0]*y - r[j][1]*x);
				double time1 = (-p[j][0]*r[j][1] + p[j][1]*r[j][0]) / (r[j][0]*y - r[j][1]*x);	
				if (time0 >= 0. && time0 <= 1. && time1 >= 1.) {m++;}
			}
			if (m % 2 != 0) {vis_colArea += areaSeg;}
			colArea += areaSeg;
		}
	}

	return vis_colArea * pow(16.9988/distance,2) / colArea;

	cout << "The ratio of visible to total area is: " << vis_colArea/colArea << endl;

	/////////////////////////////////////////////////////////////////////////// DRAW 2 DIMENSIONAL IMAGE

	double ImageX[nvertices+5];
	double ImageY[nvertices+5];
	
	for (int j = 0; j < nvertices+1; j++)
	{
		ImageX[j] = collimatorImageX[j];
		ImageY[j] = collimatorImageY[j];
	}

	for (int j = 0; j < 4; j++)
	{
		ImageX[j+nvertices+1] = slitImageX[j];
		ImageY[j+nvertices+1] = slitImageY[j];
	}

	TGraph *g = new TGraph(nvertices+5, ImageY, ImageX);
	g->SetMarkerColor(6);
	g->SetMarkerStyle(21);
	g->SetTitle("Ion's Eye View");
	g->Draw("ALP");

	return 0.;
}

double BeamGauss(double x, double y)
{
	const double sigmaxy = 	0.1; //cm
	double f = exp(-( (pow(x/sigmaxy,2)/2.) + (pow(y/sigmaxy,2)/2.) ) );
	return f;
}

void detView()
{
	const Double_t colr[] = {0.5, 1.0, 1.5, 2.5, 3.5};//radius of colimator sizes
	const Double_t deg[] = {22.5, 25., 35., 40., 45., 55., 60., 65., 75., 90., 120.};

	for (int p = 0; p < 5; p++)
	{
		for (int q = 0; q < 11; q++)
		{
			int zbins = 80;
			int xybins = 60;
			double Z = 0.;
			double Zi = (21.-19.) / zbins;
			double XYi = (1.) / xybins;
			for (int i = 0; i < zbins; i++)
			{
				double z = 19. + (i * Zi);
				double sigma = 0;
				double integral1 = 0.;
				double integral2 = 0.;
				for (int j = 0; j < xybins; j++)
				{
					double x = -0.5 + (j * XYi);
					for (int k = 0; k < xybins; k++)
					{
						double y = -0.5 + (k * XYi);
						integral1 += (IonViewRatio(x,y,z,colr[p]/10.,deg[q]) * BeamGauss(x,y) * XYi);
						integral2 += (BeamGauss(x,y) * XYi);
					}
				}
				sigma = integral1 / integral2;
				Z = Z + (sigma * Zi );
			}
			cout << 2.*colr[p] << "  " << deg[q] << "  " << Z << endl;
		}
		cout << endl;
	}
}




















