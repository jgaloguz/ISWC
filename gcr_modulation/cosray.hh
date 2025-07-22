#include <cmath>
#include <cstring>
#include <string>

// Physical scales
#define unit_length_fluid 1.496e+13
#define unit_velocity_fluid 1.000e+07
#define unit_frequency_fluid (unit_velocity_fluid / unit_length_fluid)
#define unit_diffusion_fluid (unit_velocity_fluid * unit_length_fluid)
#define unit_density_fluid 1.673e-24
#define unit_magnetic_fluid (unit_velocity_fluid * sqrt(unit_density_fluid))
#define unit_energy_particle 1.602e-12
#define unit_mass_particle (unit_energy_particle / unit_velocity_fluid / unit_velocity_fluid)
#define unit_charge_particle 4.803e-10

// Properties of the simulation
const double c_code = 2.998e+10 / unit_velocity_fluid;
const double c2_code = c_code * c_code;
const double m_p = 1.673e-24 / unit_mass_particle;
const double q_p = 4.803e-10 / unit_charge_particle;
const double B_0 = 5.0E-5 / unit_magnetic_fluid;
const double u_0 = 4.0E7 / unit_velocity_fluid;
const double omega = 2.0 * M_PI / (25.0 * 24.0 * 3600.0) / unit_frequency_fluid;
const double kappa_0 = 1.5E22 / unit_diffusion_fluid;
const double eta = 0.05;
const double f_0 = 1.0;
const double J_0 = 1.0;
const double AU_cgs = 1.496e+13;
const double r_in = 1.0 * AU_cgs / unit_length_fluid;
const double lambda = 0.1 * AU_cgs / unit_length_fluid;
const double r_0 = 1.0 * AU_cgs / unit_length_fluid;
const double r_min = 0.01 * AU_cgs / unit_length_fluid;
const double r_max = 80.0 * AU_cgs / unit_length_fluid;
const double MeV_cgs = 1.602e-06;
const double T_b = 150.0 * MeV_cgs / unit_energy_particle;
const double T_0 = 1000.0 * MeV_cgs / unit_energy_particle;
const double T_min = 10.0 * MeV_cgs / unit_energy_particle;
const double T_max = 5000.0 * MeV_cgs / unit_energy_particle;
const double cfl = 0.5;
const double delta = 1.0E-5 * AU_cgs / unit_length_fluid;
const int nbins = 100;
const int ntraj = 20000;

const std::string trajectory_fname = "trajectory.dat";
const std::string intensity_fname = "intensity_1au.dat";

// Square function
inline double Sqr(double x) {return x * x;};

// Kinetic energy from momentum
inline double EnrKin(double mom) {return c_code * sqrt(Sqr(mom) + Sqr(m_p * c_code)) - m_p * c2_code;};

// Momentum from kinetic energy
inline double Mom(double T) {return sqrt(T * (T + 2.0 * m_p * c2_code)) / c_code;};

const double p_0 = Mom(T_0);

// Velocity from momentum
inline double Vel(double mom) {return mom / (m_p * sqrt(1.0 + Sqr(mom / (m_p * c_code))));};

const double lnp_min = log(Mom(T_min));
const double lnp_max = log(Mom(T_max));
const double dlnp = (lnp_max - lnp_min) / nbins;

// Norm function
inline double Norm(double* v) {return sqrt(Sqr(v[0]) + Sqr(v[1]) + Sqr(v[2]));};

// Cross product a x b = c
inline void Cross(double* a, double* b, double* c)
{
   c[0] = a[1]*b[2]-a[2]*b[1];
   c[1] = a[2]*b[0]-a[0]*b[2];
   c[2] = a[0]*b[1]-a[1]*b[0];
};

// Unit vector function
inline void UnitVec(double* v_in, double* v_out)
{
   double mag = Norm(v_in);
   v_out[0] = v_in[0] / mag;
   v_out[1] = v_in[1] / mag;
   v_out[2] = v_in[2] / mag;
};

// Get perpendicular vector function
inline void GetPerpVec(double* v_in, double* v_out)
{
   int max_ang;
   double unit_vec[3] = {0.0};

// Find the unit vector making the largest angle with "first" (smallest component of "first")
   if(fabs(v_in[0]) < fabs(v_in[1])) {
      if(fabs(v_in[0]) < fabs(v_in[2])) max_ang = 0;
      else max_ang = 2;
   }
   else {
      if(fabs(v_in[1]) < fabs(v_in[2])) max_ang = 1;
      else max_ang = 2;
   };

   unit_vec[max_ang] = 1.0;
   Cross(v_in, unit_vec, v_out);
};

// Change vector from basis to cartesian
void ChangeFromBasis(double basis[][3], double* v)
{
   double v_new[3];
   v_new[0] = v[0] * basis[0][0] + v[1] * basis[1][0] + v[2] * basis[2][0];
   v_new[1] = v[0] * basis[0][1] + v[1] * basis[1][1] + v[2] * basis[2][1];
   v_new[2] = v[0] * basis[0][2] + v[1] * basis[1][2] + v[2] * basis[2][2];
   std::memcpy(v, v_new, 3 * sizeof(double));
};
