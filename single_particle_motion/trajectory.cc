// Include packages/dependencies here
#include <iostream>
#include <iomanip>

#include "trajectory.hh"

// Function for magnetic field
inline void Bfield(double t, double* pos, double* B)
{
// Constant field
   B[0] = Bx;
   B[1] = By;
   B[2] = Bz;
};

// Function for electric field
inline void Efield(double t, double* pos, double* E)
{
// Constant field
	E[0] = Ex;
	E[1] = Ey;
	E[2] = Ez;
};

// Lorentz force function
inline void Lorentz(double t, double* pos, double* vel, double* force)
{
   double B[3], E[3];

// Magnetic (v x B) force
   Bfield(t, pos, B);
   Cross(vel, B, force); 
   force[0] *= q_c;
   force[1] *= q_c;
   force[2] *= q_c;

// Electric force
   Efield(t, pos, E);
   force[0] += q*E[0];
   force[1] += q*E[1];
   force[2] += q*E[2];
};

// Initial condions function
inline void SetInitialConditions(double* pos, double* mom)
{
// Position
   pos[0] = Rx;
   pos[1] = Ry;
   pos[2] = Rz;
// Momentum
   mom[0] = Px;
   mom[1] = Py;
   mom[2] = Pz;
};

// Timestep function
inline double TimeStep(double t, double* pos, double* mom)
{
   double B[3];
   Bfield(t, pos, B);
   double rg = LarmorRadius(Norm(mom),Norm(B));
   double dr = 2.0 * rg * M_PI / 200.0;
// Two constraints: physical gyro-radius and absolute maximum distance traveled
   return CFL * fmin(dr, dmax) / Vel(Norm(mom));
};

// Advance position and momentum using Euler's method
void AdvanceEuler(double t, double* pos, double* mom, double dt)
{
   double vel[3], F[3];

// Find velocity (position slope) and force (momentum slope)
   Vel(mom, vel);
   Lorentz(t, pos, vel, F);
   
// Advance position and momentum to next step
   pos[0] += vel[0] * dt;
   pos[1] += vel[1] * dt;
   pos[2] += vel[2] * dt;
   mom[0] += F[0] * dt;
   mom[1] += F[1] * dt;
   mom[2] += F[2] * dt;
};

// Advance position and momentum using Heun's method
void AdvanceHeun(double t, double* pos, double* mom, double dt)
{
// TODO
};

// Advance position and momentum using RK4 method
void AdvanceRK4(double t, double* pos, double* mom, double dt)
{
// TODO
};

int main(int argc, char** argv)
{
   double t, dt, pos[3], mom[3];

// Set your initial position and momentum
   t = 0.0;
   SetInitialConditions(pos, mom);

// Iterate over time
   std::setprecision(8);
   while(t < tmax) {
// Output trajectory
      std::cerr << "\r" << std::setw(16) << t;
      std::cout << std::setw(16) << t
                << std::setw(16) << pos[0]
                << std::setw(16) << pos[1]
                << std::setw(16) << pos[2]
                << std::setw(16) << mom[0]
                << std::setw(16) << mom[1]
                << std::setw(16) << mom[2]
                << std::endl;

// Figure out a suitable timestep
      dt = TimeStep(t, pos, mom);

// Advance one step using the uncommented method
      AdvanceEuler(t, pos, mom, dt);
      t += dt;
   };
   std::cerr << "\r" << std::setw(16) << t;

   return 0;
};
