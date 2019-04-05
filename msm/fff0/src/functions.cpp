#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
using namespace std;

// Functions: //

// doubleToString
// RMSD
// norm of a vector (used also for euclidean distance)
// dot -> Scalar product
// cross -> Vector product
// dihedral angle (only if N>=5)
// histogram 2d
// KISS algorithm 
// Generate random number from uniform distribution
// Generate random number from gaussian distribution
// Save number
// Save vector
// Save matrix

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string doubleToString(double t) // Convert value into a string.
{
	std::string ch;
	ostringstream outs;
	outs << t; 
	ch = outs.str();
	return ch;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double RMSD(double *x, double *y, double *z, int N, int Np)
{
	int i,j;
	double d=0;
	
	for (i=0; i<N; i++)
		for (j=0; j<N; j++)
		{
			d = d + (x[i] - x[j])*(x[i] - x[j]) + (y[i] - y[j])*(y[i] - y[j]) + (z[i] - z[j])*(z[i] - z[j]);			
		}
	d = sqrt( d / (2*Np) );
	return d;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double norm(double *v, int n)
{
	double s=0;
	for (int i=0; i<n; i++)
		s = s + v[i]*v[i];
	
	return sqrt(s);
} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double distance(double *xt, double *yt, double *zt, int p1, int p2) //Indicate the position of the atom in the molecule
{
	double distVect[3] = {xt[p1] - xt[p2], yt[p1] - yt[p2], zt[p1] - zt[p2]};
	return norm(distVect, 3);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double dot(double *a, double *b, int n) //scalar product
{
	double s=0;
	int i;

	for (i=0; i<n; i++)
		s = s + a[i] * b[i];
	
	return s;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cross(double *a, double *b, double *c) //cross product (vettoriale)
{
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];	
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double dihedralAngle(double *x, double *y, double *z)  
{
//x=[x1 x2 x3 x4]
	double Y, X, angle;
	double c1[3], c2[3], c3[3];
	
	double b1[3] = {x[1]-x[0], y[1]-y[0], z[1]-z[0]};
	double b2[3] = {x[2]-x[1], y[2]-y[1], z[1]-z[0]};
	double b3[3] = {x[3]-x[2], y[3]-y[2], z[3]-z[2]};

	double b = norm(b2, 3);
	double b2n[3] = {b2[0]/b, b2[1]/b, b2[2]/b};

	cross(b1, b2, c1);
	cross(b2, b3, c2);
	cross(c1, c2, c3);

	Y = dot(c3, b2n, 3);
	X = dot(c1, c2, 3);

	angle = atan2(Y,X);
	//angle = fmod(angle,360.);
			  
	//angle = angle <0 ?  angle + 360. : angle;
	
	return angle;
	//atan2(dot(cross(cross(b1,b2),cross(b2,b3)),b2/norm(b2))  ,    dot(cross(b1,b2),cross(b2,b3))  );
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void hist2(double x, double y, double xmin, double delta_x, double ymin, double delta_y, int *id)
{
	int x_id, y_id;
	x_id = (int)((x-1e-6 - xmin) / delta_x)+1;
	y_id = (int)((y-1e-6 - ymin) / delta_y)+1;

	id[0] = (int)(x_id-1);	
	id[1] = (int)(y_id-1);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static unsigned int x=123456789,y=987654321,z=43219876,c=6543217; 
void init_KISS()
{
	x = rand();
	while (!(y = rand())); 
		z = rand();
	c = rand() %698769068 + 1; /*Should be less than 698769069 */
}
unsigned int JKISS() //KISS period 2^127
{
        init_KISS();
	unsigned long long t, a=4294584393ULL;
	x = 314527869*x+1234567;
	y ^= (y<<5);
	y ^= (y>>7);
	y ^= (y<<22);
	t = a*z+c;
	c = (t>>32);
	z = t;      
	return x+y+z; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double unifRand()
{
	double x;
	unsigned long long a;
	a = ((unsigned long long)JKISS()<<32) + JKISS();
	a = (a >> 12) | 0x3FF0000000000000ULL; 
	*((unsigned long long *)&x) = a; 
	return x-1.0;
	
	//return rand() / double(RAND_MAX);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double normRand(const double &variance) //box muller
{
	static bool hasSpare = false;
	static double rand1, rand2;
 
	if(hasSpare)
	{
		hasSpare = false;
		return sqrt(variance * rand1) * sin(rand2);
	}
 
	hasSpare = true;
 
	rand1 = unifRand();
	if(rand1 < 1e-100) 
		rand1 = 1e-100;
	rand1 = -2 * log(rand1);
	rand2 = unifRand() * 2 * M_PI;
 
	return sqrt(variance * rand1) * cos(rand2);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void saveValue(double x, const std::string& name)
{
	ofstream file;
	file.open(name.c_str());

	file << x << endl;

	file.close();
	file.clear();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void saveVector(double *x, int N, const std::string& name)
{
	ofstream file;
	file.open(name.c_str());

	for (int i=0 ; i<N ; i++)
		file << x[i] << endl;

	file.close();
	file.clear();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void saveMatrix(double **x, int N, int M, const std::string& name)
{
	ofstream file;
	file.open(name.c_str());
	for (int i=0 ; i<N ; i++)
	{
		for (int j=0 ; j<M; j++)
			file << x[i][j] << "  ";
		file << endl;
	}

	file.close();
	file.clear();	

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
