#include <cmath>

const double q = 1.0;
const double m = 1.0;
const double c = 1.0;
const double mc = m * c;
const double q_c = q / c;
const double s = 2.0;

const double Bx = 0.0;
const double By = 0.0;
const double Bz = 0.5;

const double Ex = 0.0;
const double Ey = 0.0;
const double Ez = 0.0;

const double Rx = -10.0;
const double Ry = 0.0;
const double Rz = 0.0;

const double Px = 1.0;
const double Py = 0.0;
const double Pz = 1.0;

const double dmax = 0.1;
const double tmax = 500.0;
const double CFL = 0.1;

// Square function
inline double Sqr(double x)
{
   return (x * x);
};

// Norm function
inline double Norm(double* v)
{
   return sqrt(Sqr(v[0])+Sqr(v[1])+Sqr(v[2]));
};

// Cross product a x b = c
inline void Cross(double* a, double* b, double* c)
{
   c[0] = a[1]*b[2]-a[2]*b[1];
   c[1] = a[2]*b[0]-a[0]*b[2];
   c[2] = a[0]*b[1]-a[1]*b[0];
};

// Function to convert momentum to velocity in magnitude
inline double Vel(double mom)
{
   return mom / (m * sqrt(1.0 + Sqr(mom / mc)));
};

// Function to convert momentum to velocity in vector form
inline void Vel(double* mom, double* vel)
{
   double mmag = Norm(mom);
   double ratio = Vel(mmag) / mmag;
 
   vel[0] = mom[0] * ratio;
   vel[1] = mom[1] * ratio;
   vel[2] = mom[2] * ratio;
};

// Function to compute gyroradius
inline double LarmorRadius(double mmag, double Bmag)
{
   return mmag * c / (q * Bmag);
};
