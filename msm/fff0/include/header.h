std::string doubleToString(double);
//S T R U C T U R E S 
struct par
{
	double nsteps, dt, sigma, beta; //simulation parameters	
	int delta, tau; //lagtime for MSM
	int mx, my; //number of states for MSM
	double xmin, ymin, xmax, ymax, xdelta, ydelta; //intervals for MSM
};
double RMSD(double *, double *, double *, int , int );
double distance(double *, double *, double *, int , int );
double norm(double *, int );
double dot(double *, double *, int);
void cross(double *, double *, double *);
double dihedralAngle(double *, double *, double *);  
void hist2(double , double , double , double , double , double , int *);

void init_KISS();
unsigned int JKISS(); 
double unifRand();
double normRand(const double &);

void saveValue(double , const std::string&);
void saveVector(double *, int, const std::string&);
void saveMatrix(double **, int, int, const std::string&);


void transitionMatrix1d(int , double *, par, const std::string&);
void transitionMatrix2d(double *, double *, par, const std::string&);

double V(double, double, double);
double DVx(double, double, double);
double DVy(double, double, double);


 
