#include <iostream>
#include <stdlib.h>
#include <math.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*### Potential 1
k1=0.003
a1=3.3
### Potential 2
k2=1
a2=3.5
###
k12 =-50
c=1
*/
static double k1=0.003, a1=3.3, k2=1.0, a2=3.5, k12=-50, c=1;  						

double V(double x, double y, double z) //actually never used
{
	return k1 * pow( (x*x - a1*a1), 4) + k2 * pow(y*y - a2*a2, 2) + k12 / sqrt((x - y)*(x - y) + c*c) ; 
}
double DVx(double x, double y, double z)
{	
	return 8 * k1 * x * pow( (x*x - a1*a1), 3) - k12 * (x - y) / pow((x - y)*(x - y) + c*c, 1.5); 
}
double DVy(double x, double y, double z)
{
	return 4 * k2 * y * (y*y - a2*a2) + k12 * (x - y) / pow((x - y)*(x - y) + c*c, 1.5); 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* DOUBLE DOUBLE WELL POTENTIAL (WORKS FINE)
double V(double x, double y, double z) //actually never used
{
	return (x*x-1)*(x*x-1) + (y*y-1)*(y*y-1) + Kx*x + Ky*y; 
}
double DVx(double x, double y, double z)
{	
	return 4*x*(x*x-1) + Kx; 
}
double DVy(double x, double y, double z)
{
	return 4*y*(y*y-1) + Ky;
}
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
