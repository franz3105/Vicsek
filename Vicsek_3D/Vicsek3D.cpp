//compile with: g++ Vicsek3D.cpp nrutil.cpp random_mars.cpp -o MD -O2 -w

#include <fstream>
#include <iomanip>
#include <iostream>


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string>
#include <string.h>
#include <sstream>



#include "random_mars.h"
#include "Vicsek3D.h"
#include "nrutil.h"

#define PI 3.1415926


int main (void) {

  // system parameters
  //den = 2.0; //number density
  //eta = 0.3;
  delta = 0.5;
  eta = 0.5;
  Tinit = 1.0; //temperature in units of epsilon/kB.
  v = 0.5;
  N = 10000;
  //simulation parameters
  dt = 0.1; //time-step length
  //Nsim = 2000; //8001; //8001; //16001; //number of timesteps in Eq = 80001
  Neq = 10000; //8000; //number of timesteps to equilibrate = 2000
  Nsim = 100000;
  freq = 500;
  den = 0.5;
  std::fstream parameters("Vparam", std::ios_base::in);
  parameters >> N >> den >> eta;
  delta=eta;
  N = 10000;


  N_MSD=3200;

  // calculate constants for cell list calculation
  rc = 1;
  L = pow(N/den, 1/3.);//pow(N/den,1./3.); //box length
  Ncell = (int) L/rc;
  Nmax = 500;

  printf("Ncell = %d \n",Ncell);

  // create arrays
  x = dvector(0,N-1);
  y = dvector(0,N-1);
  z = dvector(0,N-1);

  x_accum = dvector(0,N-1);
  y_accum = dvector(0,N-1);
  z_accum = dvector(0,N-1);

  x_init = dvector(0,N-1);
  y_init = dvector(0,N-1);
  z_init = dvector(0,N-1);


  vx = dvector(0,N-1);
  vy = dvector(0,N-1);
  vz = dvector(0,N-1);

  vx_init = dvector(0,N-1);
  vy_init = dvector(0,N-1);
  vz_init = dvector(0,N-1);


  fx = dvector(0,N-1);
  fy = dvector(0,N-1);
  fz = dvector(0,N-1);

  theta = dvector(0,N-1);
  alpha = dvector(0,N-1);
  SF = dvector(0,NSF-1);
  SFimg = dvector(0,NSF-1);
  realsum = dvector(0,3*NSF-1);
  imgsum = dvector(0,3*NSF-1);
  cell_list_index = ivector(0,N-1);
  cell_list = imatrix(0,(Ncell*Ncell*Ncell)-1,0,Nmax-1);


  x_tracked =  dvector(0,N);

  cvvx_accum_array = dvector(0,N_MSD);
  cvvx_init_array = dvector(0,N);


  ratio = (double) Nsim/ (double) N_MSD;



  // initialize random number generator
  random_mars = new RanMars(time(NULL));

  // initialize particles
  lattice_packing();

  // initialize velocities
  random_velocities();



  //print_positions(0);
  //printf("Ciao!\n");
  create_cell_lists();
  //printf("Ciao!\n");


  // perform simulation
  struct timeval tp;
  gettimeofday(&tp, NULL);
  long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
  printf("Start Simulation! \n");
  printf("step, m\n");

  eq_temp_sum = 0.0;
  eq_temp_count = 0;

  //equilibration:
  for (int i=0; i<Neq; i++) {

    md_step();

  }

  std::stringstream sstr;
  sstr << "output.dat";
  const std::string tmp = sstr.str();
  const char* cstr = tmp.c_str();

  FILE * fout = fopen(cstr,"w");


  for (int i=0; i<Nsim; i++) {

    md_step();
    //printf("%f \n", i*dt);
    if (i % freq == 0){
        double m = calculate_order();
        printf("%f %f \n", i*dt, m);
        fprintf(fout, "%f %f \n", i*dt, m);
        //print_positions(i);
    }
  }
  fclose(fout);
}
/*********************************************************************************************************************************/




void output_MSD_VCC(){

  std::stringstream sstr;
  sstr << "MSD_VCC.dat";
  const std::string tmp = sstr.str();
  const char* cstr = tmp.c_str();

  FILE * fout = fopen(cstr,"w");

  for (int k=0; k<N_MSD; k++)
    {
      msdx_accum_array[k] = msdx_accum_array[k]/(ratio * (float)  N);
      cvvx_accum_array[k] = cvvx_accum_array[k]/(ratio * (float)  N);
      fprintf(fout,"%f %f %f\n",k*dt, msdx_accum_array[k], cvvx_accum_array[k]);
    }

  fclose(fout);
}

/*********************************************************************************************************************************/

void lattice_packing() {
  int K = (int) (pow(N,1./3.)); //number of particles per row
  if (pow(K,3)<N){ K++;} //correct K if not all the particles will fit
  double a = L/((double) K); //spacing between particles

  for (int i=0; i<K; i++) {
    for (int j=0; j<K; j++) {
      for (int n=0; n<K; n++) {
        int index = i*K*K+j*K+n;
        if(index < N){
          x[index] = i*a-0.5*L;
          y[index] = j*a-0.5*L;
          z[index] = n*a-0.5*L;
        }
      }
    }
  }

  printf("Initialized...  Density %f, Number of particles %d, Box length %f\n",\
	 den,N,L);
}

/*********************************************************************************************************************************/

void random_velocities() {

  //assign each angle randomly in a uniform distribution
  for (int i=0; i<N; i++) {
    theta[i] = (random_mars->uniform() - 0.5)*2*PI;
    alpha[i] = (random_mars->uniform() - 0.5)*PI;
    vx[i] = v*cos(theta[i])*sin(alpha[i]);
    vy[i] = v*sin(theta[i])*sin(alpha[i]);
    vz[i] = v*cos(alpha[i]);
  }

}


/*********************************************************************************************************************************/

//calculate the potential energy of the system
// (as the sum of all pair-wise interactions)
double calculate_potential() {

  double potential_sum = 0.0;
  double r_cut = pow(2,1./6.);
  double r_cut2 = r_cut*r_cut;

  //loop through all particles
  for (int i=0; i<N; i++) {
    //determine x,y,z cell index
    int cell_list_indexi = cell_list_index[i];
    int cell_list_z  = cell_list_indexi % Ncell;
    int cell_list_y  = ((cell_list_indexi - cell_list_z)/Ncell) % Ncell;
    int cell_list_x  = (cell_list_indexi - cell_list_z - cell_list_y*Ncell)\
      /(Ncell*Ncell);

    //loop through neighboring cells
    for (int a=cell_list_x - 1; a<= cell_list_x + 1; a++) {
      for (int b=cell_list_y - 1; b<= cell_list_y + 1; b++) {
        for (int c=cell_list_z - 1; c<= cell_list_z + 1; c++) {
          //account for periodic boundary conditions
          int aloc = a;
          int bloc = b;
          int cloc = c;
          if (aloc < 0) aloc += Ncell;
          if (aloc >= Ncell) aloc -= Ncell;
          if (bloc < 0) bloc += Ncell;
          if (bloc >= Ncell) bloc -= Ncell;
          if (cloc < 0) cloc += Ncell;
          if (cloc >= Ncell) cloc -= Ncell;

          //loop through all particles in the cell
          for (int j=0; j< Nmax; j++) {
	    int jloc = cell_list[aloc*Ncell*Ncell + bloc*Ncell + cloc][j];
	    if (jloc != -1 && jloc != i) {
              //calculate interaction between particles
	      double dxloc = x[i] - x[jloc];
              double dyloc = y[i] - y[jloc];
              double dzloc = z[i] - z[jloc];
	      pbc(dxloc, dyloc, dzloc);
	      double r2 = dxloc*dxloc + dyloc*dyloc + dzloc*dzloc;
              if(r2<r_cut2){
                double r2i = 1.0/r2;
                double r6i = r2i * r2i * r2i;
                double r12i = r6i * r6i;
                potential_sum += 4.*(r12i - r6i) + 1.;
              }
            }
          }

        }
      }
    }

  }

  //account for double counting every pair of particles
  potential_sum /= 2.0;

  return potential_sum;

}


void update_angles() {
  int nint;
  double r_cut = 1;
  double r_cut2 = r_cut*r_cut;
  double * cos_theta_new;
  double * sin_theta_new;
  double * cos_alpha_new;
  double * sin_alpha_new;
  cos_theta_new = dvector(0, N-1);
  sin_theta_new = dvector(0, N-1);
  cos_alpha_new = dvector(0, N-1);
  sin_alpha_new = dvector(0, N-1);

  //loop through all particles
  for (int i=0; i<N; i++) {
    //determine x,y,z cell index
    int cell_list_indexi = cell_list_index[i];
    int cell_list_z  = cell_list_indexi % Ncell;
    int cell_list_y  = ((cell_list_indexi - cell_list_z)/Ncell) % Ncell;
    int cell_list_x  = (cell_list_indexi - cell_list_z - cell_list_y*Ncell)\
      /(Ncell*Ncell);

    //loop through neighboring cells
    for (int a=cell_list_x - 1; a<= cell_list_x + 1; a++) {
      for (int b=cell_list_y - 1; b<= cell_list_y + 1; b++) {
        for (int c=cell_list_z - 1; c<= cell_list_z + 1; c++) {
          //account for periodic boundary conditions
          int aloc = a;
          int bloc = b;
          int cloc = c;
          if (aloc < 0) aloc += Ncell;
          if (aloc >= Ncell) aloc -= Ncell;
          if (bloc < 0) bloc += Ncell;
          if (bloc >= Ncell) bloc -= Ncell;
          if (cloc < 0) cloc += Ncell;
          if (cloc >= Ncell) cloc -= Ncell;

          //loop through all particles in the cell
          for (int j=0; j< Nmax; j++) {
	    int jloc = cell_list[aloc*Ncell*Ncell + bloc*Ncell + cloc][j];
	    if (jloc != -1 && jloc != i) {
              //calculate interaction between particles
	      double dxloc = x[i] - x[jloc];
          double dyloc = y[i] - y[jloc];
          double dzloc = z[i] - z[jloc];
	      pbc(dxloc, dyloc, dzloc);
	      double r2 = dxloc*dxloc + dyloc*dyloc + dzloc*dzloc;
              if(r2<r_cut2){
              cos_theta_new[i] += cos(theta[jloc]);
              sin_theta_new[i] += sin(theta[jloc]);
              cos_alpha_new[i] += cos(alpha[jloc]);
              sin_alpha_new[i] += sin(alpha[jloc]);
              }
            }
          }

        }
      }
    }
  }
  for (int i=0; i<N; i++) {
    theta[i]=atan2(sin_theta_new[i], cos_theta_new[i])+eta*(random_mars->uniform() -0.5)*sqrt(dt)*2*PI;
    alpha[i]=atan2(sin_alpha_new[i], cos_alpha_new[i])+delta*(random_mars->uniform() -0.5)*sqrt(dt)*PI;
    vx[i] = v*cos(theta[i])*sin(alpha[i]);
    vy[i] = v*sin(theta[i])*sin(alpha[i]);
    vz[i] = v*cos(alpha[i]);
  }


  free_dvector(cos_theta_new,0,N-1);
  free_dvector(sin_theta_new,0,N-1);
  free_dvector(cos_alpha_new,0,N-1);
  free_dvector(sin_alpha_new,0,N-1);

}
/*********************************************************************************************************************************/

//update the positions of all particles
//(with one iteration of the verlet integrator.)
void md_step(){

  for (int i=0; i<N; i++) {
    //second half step: update positions
    x[i] = x[i] + dt*vx[i];
    y[i] = y[i] + dt*vy[i];
    z[i] = z[i] + dt*vz[i];

    x_accum[i] = x_accum[i] + dt*vx[i];
    y_accum[i] = y_accum[i] + dt*vy[i];
    z_accum[i] = z_accum[i] + dt*vz[i];

    x_tracked[i] = x_tracked[i]+dt*vx[i];

    //apply periodic boundary conditions
    pbc(x[i],y[i],z[i]);
  }

  //recalculate cell list
  create_cell_lists();

  //reculculate forces, now that particles have moved
  update_angles();

}

/*********************************************************************************************************************************/

// calculate and store the force on every particle
//(from the derivative of the interaction potential,
//summed over neighbouring particles)
void calculate_forces() {

  double r_cut = pow(2,1./6.);
  double r_cut2 = r_cut*r_cut;
  virial = 0.0;

  //loop through all particles
  for (int i=0; i<N; i++) {
    //initilize arrays to hold forces
    fx[i] = fy[i] = fz[i] = 0.0;

    //determine x,y,z cell index
    int cell_list_indexi = cell_list_index[i];
    int cell_list_z  = cell_list_indexi % Ncell;
    int cell_list_y  = ((cell_list_indexi - cell_list_z)/Ncell) % Ncell;
    int cell_list_x  = (cell_list_indexi - cell_list_z - cell_list_y*Ncell)\
      /(Ncell*Ncell);

    //loop through neighboring cells
    for (int a=cell_list_x - 1; a<= cell_list_x + 1; a++) {
      for (int b=cell_list_y - 1; b<= cell_list_y + 1; b++) {
        for (int c=cell_list_z - 1; c<= cell_list_z + 1; c++) {
          //account for periodic boundary conditions
          int aloc = a;
          int bloc = b;
          int cloc = c;
          if (aloc < 0) aloc += Ncell;
          if (aloc >= Ncell) aloc -= Ncell;
          if (bloc < 0) bloc += Ncell;
          if (bloc >= Ncell) bloc -= Ncell;
          if (cloc < 0) cloc += Ncell;
          if (cloc >= Ncell) cloc -= Ncell;

          //loop through all particles in the cell
          for (int j=0; j< Nmax; j++) {
	    int jloc = cell_list[aloc*Ncell*Ncell + bloc*Ncell + cloc][j];
	    if (jloc != -1 && jloc != i) {
              //calculate force between particles
              //(and add it to the sum for particle i)
	      double dxloc = x[i] - x[jloc];
              double dyloc = y[i] - y[jloc];
              double dzloc = z[i] - z[jloc];
	      pbc(dxloc, dyloc, dzloc);
	      double r2 = dxloc*dxloc + dyloc*dyloc + dzloc*dzloc;
              if(r2<r_cut2){
                double r2i = 1.0/r2;
                double r6i = r2i * r2i * r2i;
                double r12i = r6i * r6i;
                double fc = 48. * (r12i - 0.5*r6i) * r2i;
                fx[i] += fc * dxloc;
                fy[i] += fc * dyloc;
                fz[i] += fc * dzloc;
                virial += fc*dxloc*dxloc + fc*dyloc*dyloc + fc*dzloc*dzloc;
              }
            }
          }

        }
      }
    }

  }

  virial /= 2.0;

}

/*********************************************************************************************************************************/

//calculate the kinetic energy of the system
double calculate_kinetic(){
  double ke_sum = 0.0;
  for (int i=0; i<N; i++) {
    ke_sum += vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  double ke = ke_sum/2.0;
  return ke;
}

/*********************************************************************************************************************************/

// calculate and print g(r)
void print_gr(int step){

  // constants for histogram
  int Nhisto = 100;
  double drhisto = 0.05;

  // declare and initialize histogram
  double * histo = dvector(0,Nhisto-1);
  for (int i=0; i<Nhisto; i++) {
    histo[i] = 0.0;
  }

  for (int i=0; i<N; i++) {
    double xi = x[i];
    double yi = y[i];
    double zi = z[i];

    for (int j=0; j<N; j++) {
      if (j!=i) {
	double dxloc = xi - x[j];
	double dyloc = yi - y[j];
	double dzloc = zi - z[j];
	pbc(dxloc, dyloc, dzloc);

	double dr = sqrt(dxloc*dxloc + dyloc*dyloc + dzloc*dzloc);
	int drint = (int) ( dr/drhisto + 0.5);

	//printf("dr %f drint %d\n",dr,drint);
	// bin the particle distance into a histogram
	if (drint < Nhisto) histo[drint]++;
      }
    }

  }

  // write histogram data
  std::stringstream sstr;
  sstr << "gr_step" << step << ".dat";
  const std::string tmp = sstr.str();
  const char* cstr = tmp.c_str();

  FILE * out = fopen(cstr,"w");

  for (int i=0; i<Nhisto; i++) {
    double V = 4.* PI * ( pow(i*drhisto + drhisto,3)- pow(i*drhisto,3) )/3.;
    fprintf(out,"%f %f\n",i*drhisto,histo[i]*L*L*L/((double) N * N *V));

  }

  fclose(out);
}

/*********************************************************************************************************************************/

// print particle positions
void print_positions(int i){

  std::stringstream sstr;
  sstr << "positions_step_" << i << ".dat";
  const std::string tmp = sstr.str();
  const char* cstr = tmp.c_str();

  FILE * out = fopen(cstr,"w");

  for (int i=0; i<N; i++) {
    double xloc = x[i];
    double yloc = y[i];
    double zloc = z[i];
    pbc(xloc, yloc, zloc);
    fprintf(out,"%f %f %f %f %f %f\n",xloc,yloc,zloc,vx[i],vy[i],vz[i]);
  }

  fclose(out);

}

/*********************************************************************************************************************************/

void create_cell_lists(){
  // initialize cell lists
  for (int a=0; a<Ncell; a++) {
    for (int b=0; b<Ncell; b++) {
      for (int c=0; c<Ncell; c++) {
        for (int j=0; j< Nmax; j++) {
	  cell_list[a*Ncell*Ncell + b*Ncell + c][j] = -1;
        }
      }
    }
  }

  // find cell lists for particles
  for (int i=0; i<N; i++) {
    int xint = (int) ( Ncell* (x[i]/L + 0.5) );
    int yint = (int) ( Ncell* (y[i]/L + 0.5) );
    int zint = (int) ( Ncell* (z[i]/L + 0.5) );

    if(xint<0) xint=0;
    if(yint<0) yint=0;
    if(zint<0) zint=0;
    if(xint>(Ncell-1)) xint=Ncell-1;
    if(yint>(Ncell-1)) yint=Ncell-1;
    if(zint>(Ncell-1)) zint=Ncell-1;

    cell_list_index[i] = xint * Ncell*Ncell + yint*Ncell + zint;


    bool particle_stored = false;
    for (int j=0; j< Nmax; j++) {
      if (cell_list[xint*Ncell*Ncell + yint*Ncell + zint][j] == -1) {
	cell_list[xint*Ncell*Ncell + yint*Ncell + zint][j] = i;
        particle_stored = true;
	break;
      }
    }
    if (!particle_stored){
      printf("ERROR: cell list storage exceeded!!! Nmax must be increased! \n");
    }
  }

}


/*********************************************************************************************************************************/
// consider periodic boundary conditions
void pbc(double &x, double &y, double &z){
  if (x >= L/2.0) x -= L;
  if (x < -L/2.0) x += L;
  if (y >= L/2.0) y -= L;
  if (y < -L/2.0) y += L;
  if (z >= L/2.0) z -= L;
  if (z < -L/2.0) z += L;
}

double calculate_msd()
{
  double msd_now = 0;

  for (int i=0; i<N; i++)
    {
      msd_now +=
	(x_accum[i]-x_init[i])*(x_accum[i]-x_init[i])+
	(y_accum[i]-y_init[i])*(y_accum[i]-y_init[i])+
	(z_accum[i]-z_init[i])*(z_accum[i]-z_init[i]);
    }

  msd_now = msd_now / (float)N;

  return msd_now;
}

double calculate_fourth()
{
  double fourth = 0;

  double dx,dy,dz;

  for (int i=0; i<N; i++)
    {
      dx = (x_accum[i]-x_init[i]);
      dy = (y_accum[i]-y_init[i]);
      dz = (z_accum[i]-z_init[i]);
      fourth += dx*dx*dx*dx+dy*dy*dy*dy+dz*dz*dz*dz;

    }

  fourth = fourth / (float)N;

  return fourth;
}


void initialize_msd_variables()
{
  for (int i=0; i<N; i++)
    {
      x_init[i] = x_accum[i] = x[i];
      y_init[i] = y_accum[i] = y[i];
      z_init[i] = z_accum[i] = z[i];

    }
}

void initialize_vel_auto_variables()
{
  for (int i=0; i<N; i++)
    {
      vx_init[i] = vx[i];
      vy_init[i] = vy[i];
      vz_init[i] = vz[i];

    }
}


double calculate_vel_auto()
{

  double corr = 0;

  for (int i = 0; i<N; i++ )
    {
      corr += vx[i]*vx_init[i] + vy[i]*vy_init[i]+ vz[i]*vz_init[i];
    }

  corr = corr / (float) N;

  return corr;

}

void calculate_msd_2(double &msd_x, double &msd_y, double &msd_z)
{

  for (int i=0; i<N; i++)
    {
      msd_x += (x_accum[i] - x_init[i])*(x_accum[i] - x_init[i]);
      msd_y += (y_accum[i] - y_init[i])*(y_accum[i] - y_init[i]);
      msd_z += (z_accum[i] - z_init[i])*(z_accum[i] - z_init[i]);
    }

  msd_x = msd_x/(float)N;
  msd_y = msd_y/(float)N;
  msd_z = msd_z/(float)N;

}

void calculate_vel_auto_2(double &corr_vx,double &corr_vy, double &corr_vz)
{
  for (int i=0; i<N; i++)
    {


      corr_vx += vx[i] * vx_init[i];
      corr_vy += vy[i] * vy_init[i];
      corr_vz += vz[i] * vz_init[i];
    }

  corr_vx = corr_vx/N;
  corr_vy = corr_vy/N;
  corr_vz = corr_vz/N;
}

double calculate_order()
{
  double m_x=0.0;
  double m_y=0.0;
  double m_z = 0.0;
  double m=0.0;

  for (int i=0; i<N; i++)
  {
    m_x+=vx[i];
    m_y+=vy[i];
    m_z+=vz[i];
  }

  m = sqrt(m_x*m_x+m_y*m_y + m_z*m_z)/v;

  return m/((double)N);
}

