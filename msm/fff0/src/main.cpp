#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>

#include "header.h"
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

int main(int argc, char* argv[])
{		

	int nstxout = 1; //prints every nstxout steps; if 0 does not print the trj

	std::string folder(argv[1]); 			//Directory where to save all data
	const char * folderChar = folder.c_str(); 	//Convert from "std::string" to "const char *" --> necessary for the function mkdir
	int status = mkdir(folderChar, S_IRWXU); 	//Create a new directory, if it does not exist

	//parameters
	int nsteps = 100000000; //number of timesteps
	double dt = 0.001, sdt=sqrt(dt);
	double temperature = 3;
	double kboltz = 1;
	double mass = 1;
	double sigma = 15;
//cout << "ok" << endl;
	//M S M
	double delta=1;	//move the window every delta timesteps
	const int nTau = 11; //number of timelags to test
	double tau[nTau] = {5, 10, 20, 30, 50, 70, 100, 200, 300, 400, 500};	//lagtime

	int mx = 50; //number of states (bins)
	int my = 50;

    double aX[mx]; //axes
	double aY[my];

	//Integration algorithm
	double Fx, Fy; 			
	double xt, yt; //position 
	double etaX, etaY; //noise	
	
	double *X = new double[nsteps];
	double *Y = new double[nsteps];
    

    
	//others
	int i,j,k,t,u;
	par parameters;
	string outputTrj  = folder + "trajectory.txt";
	ofstream fileTrj, fileLog; // to print the trajectory
	string outputLog  = folder + "log.txt"; //save here all the parameters

/////////////////////////
cout << "//////////////////////////////////////////// S T A R T ////////////////////////////////////////////"<< endl;

	cout << endl;
	cout << endl;
	cout << endl;
	cout << "# # S i m u l a t i o n   p a r a m e t e r s: # #" << endl;
	cout << "Timesteps= " << nsteps << endl;
	cout << "dt= " << dt << endl;
 	cout << "sigma= " << sigma << endl;
	cout << endl;

// C R E A T E   L O G   F I L E
	fileLog.open(outputLog.c_str());
	fileLog << "Saving directory: " << folder << endl;
	fileLog << "Simulation parameters:" << endl;
	fileLog << "K boltz= " << kboltz << endl;
	fileLog << "Timesteps= " << nsteps << endl;
	fileLog << "dt= " << dt << endl;
	fileLog << "Output trj= " << nstxout << endl;  
	fileLog << "MSM parameters:" << endl;
	fileLog << "mx= " << mx << endl; //number of states (bins)
	fileLog << "my= " << my << endl; //number of states (bins)

	fileLog.close();

////////initialize
	//random seed: 
	srand(time(NULL));
	
	xt = normRand(1);
	yt = normRand(1);

	X[0] = xt;
	Y[0] = yt;

////////Simulation
	cout << "# # S i m u l a t i o n   s t a r t # #" << endl;
	fileTrj.open(outputTrj.c_str());

	for (t=1; t<nsteps; t++)	//At time t....
	{	
		if( fmod((double)(t+1) , (double)((nsteps)/100.))==0)
			cout << "\r" << "Integrating... "<< 100*(double)(t+1)/(double)(nsteps) << "% " <<flush;


		Fx = - DVx(xt, yt , 0);
		Fy = - DVy(xt, yt , 0);

		etaX = normRand(1);			//Random force 
		etaY = normRand(1);		
		xt = xt + Fx * dt + sigma * etaX * sdt;
		yt = yt + Fy * dt + sigma * etaY * sdt;

		X[t] = xt;
		Y[t] = yt;
		
////////Save trajectory
		if (fmod((double)(t+1),(double)(nstxout)) == 0)
		{
			fileTrj << xt << " ";
			fileTrj << yt << " " << endl;
		}
	}//close time loop


	cout << endl;
	cout << "Simulation, done." << endl;



/////////////// M S M 
	double xmin = *min_element(X, X+nsteps) - 0.001;
	double xmax = *max_element(X, X+nsteps) + 0.001;

	double xdelta = (xmax - xmin) / mx;
    
    double ymin = *min_element(Y, Y+nsteps) - 0.001;
	double ymax = *max_element(Y, Y+nsteps) + 0.001;
    
	double ydelta = (ymax - ymin) / my;

    for (i=0; i<mx; i++)
        aX[i] = xmin + i*xdelta;

    for (i=0; i<my; i++)
        aY[i] = ymin + i*ydelta;
    
	parameters.nsteps = nsteps;
	parameters.mx = mx;	
	parameters.xmin = xmin;	
	parameters.xmax = xmax;
	parameters.xdelta = xdelta;
	parameters.my = my;	
	parameters.ymin = ymin;
	parameters.ymax = ymax;	
	parameters.ydelta = ydelta;
	parameters.delta = delta;


	cout << endl;
	cout << "# # M S M   s t a r t # #" << endl;
	for (i=0; i <nTau; i++)
	{
		parameters.tau = tau[i];
		cout << "Lagtime= "<< tau[i] <<endl; 

		transitionMatrix1d(0, X, parameters, folder); 

		transitionMatrix1d(1, Y, parameters, folder); 

		transitionMatrix2d(X, Y, parameters, folder); 
	}

	saveVector(tau, nTau, folder + "tau.txt");
	saveVector(aX, mx, folder + "Ax.txt");
	saveVector(aY, my, folder + "Ay.txt");

cout << endl;
cout << endl;

cout << "////////////////////////////////////////////// E N D //////////////////////////////////////////////"<< endl;

	delete [] X;
	delete [] Y;

	return 0;
}//end main   

