#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include "header.h"
using namespace std;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void transitionMatrix1d(int f, double *x, par parameters, const std::string& folder) 
{

////////PARAMETERS
	int T = parameters.nsteps;		
	int tau = parameters.tau;
	int delta = parameters.delta;

	int mx; 
    double xmin, xmax;
	double xdelta;

	if (f==0)
	{	
		mx = parameters.mx;	
		xmin = parameters.xmin;	
		xmax = parameters.xmax;
		xdelta = parameters.xdelta;
	}
	else if (f==1)
	{	
		mx = parameters.my;	
		xmin = parameters.ymin;	
		xmax = parameters.ymax;
		xdelta = parameters.ydelta;
	}

/////////  COUNTING MATRIX
	int i,j;
	double C[mx][mx]; 
    int start, end;

	int k;

	for (i=0; i<mx; i++)
		for (j=0; j<mx; j++)
			C[i][j]=0;


	for (k=0; k<T-tau; k=k+delta)
	{	
		if( fmod((double)(k+1), (double)((T-tau)/100.))==0)
			cout << "\r" << "Counting matrix..." <<  100*(double)(k+1)/(double)(T-tau) << "% " <<std::flush;

        //Get index of start
        start = floor((x[k] - xmin) / xdelta);

        //Get index of end
        end = floor((x[k+tau] - xmin) / xdelta);

        //cout << x[k] << "  " <<  xmin << endl;
        // increase count in the corresponding transition
        C[start][end]++;
        // enforce detailed balance
        C[end][start]++ ;  

	}
	cout << endl;

//////// SAVE MATRIX
	string name;
	if (f==0)
		name = ".//" + folder + "Cx" + doubleToString(tau) + ".txt";  
	else if (f==1)
		name = ".//" + folder + "Cy" + doubleToString(tau) + ".txt";  

	ofstream file;
	file.open(name.c_str());
	for (int i=0 ; i<mx ; i++)
	{
		for (int j=0 ; j<mx; j++)
			file<< C[i][j] << "  ";
		file << endl;
	}

	file.close();
	file.clear();	
	
	cout <<"Transition matrix " << doubleToString(tau) << ", done." <<endl;
	cout << endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void transitionMatrix2d(double *x, double *y, par parameters, const std::string& folder) 
{
////////PARAMETERS
	int T = parameters.nsteps;		 
	int tau = parameters.tau;
	int delta = parameters.delta;

	int mx, my;
    double xmin, xmax, ymin, ymax;
	double xdelta, ydelta;


	mx = parameters.mx;	
	xmin = parameters.xmin;	
	xmax = parameters.xmax;
	xdelta = parameters.xdelta;	
	my = parameters.my;	
	ymin = parameters.ymin;	
	ymax = parameters.ymax;
	ydelta = parameters.ydelta;
	
///////// COUNTING MATRIX

	int i,j;	
    int m=mx*my;
    
       
	double** C = new double*[m];
    for (i = 0; i < m; ++i)
	    C[i] = new double[m];	
         
         
	int x_start, y_start, start;
	int x_end, y_end, end;
	int k;

                

	for (i=0; i<m; i++)
		for (j=0; j<m; j++)
			C[i][j]=0;
        


	for (k=0; k<T-tau; k=k+delta)
	{	
	//	if( fmod((double)(k+1), (double)((T-tau)/100.))==0)
                //cout << "\r" << "Counting matrix..." <<  100*(double)(k+1)/(double)(T-tau) << "% " <<std::flush;

        /*
        //Get index of start
        x_start   = x[k] < xmin ? 0 
                : x[k] >= xmax ? mx-1 
		        : floor((x[k] - xmin) / xdelta);
        
        y_start   = y[k] < ymin ? 0 
                : y[k] >= ymax ? my-1 
		        : floor((y[k] - ymin) / ydelta);
		         
         //Get index of end
        x_end     = x[k+tau] < xmin ? 0
                : x[k+tau] >= xmax ? mx-1 
		        : floor((x[k+tau] - xmin) / xdelta);
        
        y_end     = y[k+tau] < ymin ? 0 
                : y[k+tau] >= ymax ? my-1 
		        : floor((y[k+tau] - ymin) / ydelta);
		         
		start = x_start + my*y_start;
		end   = x_end + my*y_end; 
        */
            x_start = floor((x[k] - xmin) / xdelta); //cout << xmin << "  " << x[k] << "  " <<  x_start << endl;
            x_end = floor((x[k+tau] - xmin) / xdelta); //cout << xmin << "  " << x[k+tau] << "  " <<  x_end << endl;

            y_start = floor((y[k] - ymin) / ydelta); //cout << ymin << "  " << y[k] << "  " <<  y_start << endl;
            y_end = floor((y[k+tau] - ymin) / ydelta); //cout << ymin << "  " << y[k+tau] << "  " <<  y_end << endl;
            
            start = x_start + mx * y_start; 
            end = x_end + my * y_end;
        

        

        
        //cout << start << "  " << end << endl;
        
                    //cout << k << "   " << T-tau << endl;    
       
        // increase count in the corresponding transition
        C[start][end]++;
        // enforce detailed balance
        C[end][start]++ ;   
 //cout << " dsaas " <<endl;
	}
	cout << endl;


//////// SAVE MATRIX
	string name;
	name = ".//" + folder + "Cxy" + doubleToString(tau) + ".txt";  

	ofstream file;
	file.open(name.c_str());
	for (int i=0 ; i<m ; i++)
	{
		for (int j=0 ; j<m ; j++)
			file<< C[i][j] << "  ";
		file << endl;
	}

	file.close();
	file.clear();	
 

	for (i = 0; i < m; ++i)
		delete [] C[i];
	delete [] C;
 
	cout <<"Transition matrix2d " << doubleToString(tau) << ", done." <<endl;
 	cout << endl;
 	cout << endl;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
