//compile with: g++ Vicsek.cpp nrutil.cpp random_mars.cpp -o MD -O2 -w //-static-libgcc -static-libstdc++

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

#include "Vicsek.h"
#include "nrutil.h"
#include "random_mars.h"

#define PI 3.1415926

 
int main () {
  
  // system parameters
  //den = 0.8; number density
  //N = 1000; // 1000; //number of particles
  //eta = 0.1;
  v = 0.5;
  std::fstream parameters("Vparam", std::ios_base::in);
  parameters >> N >> den >> eta;
  N=1000;
  
  //simulation parameters
  
  N_blocks = 10; // 10; // number of blocks
  N_levels = 6; //5; // 8;  4;  // number of levels 

  Nsim = (int) pow(N_blocks,N_levels); //8001; //8001; //16001; //number of timesteps in Eq = 80001
  
  // declare an array for blockings
  blocking_sumx = f3tensor(0,N_blocks,0,N_levels,0,N);
  blocking_sumy = f3tensor(0,N_blocks,0,N_levels,0,N);
  
  msd_blocking = dmatrix(0,N_blocks,0,N_levels);
  count_block = dmatrix(0,N_blocks,0,N_levels);
  
  // initialize array with zeros
  for (int i=0; i<N_blocks; i++) 
  {
    for (int j=0; j<N_levels; j++)
      {
	      msd_blocking[i][j] = 0;

	      count_block[i][j]  = 0;
	      for (int k=0; k<N; k++) {
	        blocking_sumx[i][j][k]=0;
	        blocking_sumy[i][j][k]=0;
	      }
      }
  }
  
  dt = 0.1; //time-step length
  Neq = 50000; //8000; //number of timesteps to equilibrate = 2000
  int freqeq = 1; //frequency of sampling (in steps) = 400
  
  N_MSD=3200;
  int MSDyn=0;

  // initialize random number generator
  random_mars = new RanMars(time(NULL));

  
  // calculate constants for cell list calculation
  rc = 2; 
  L = pow(N/den,1./2.); //box length
  Ncell = (int)round(L/rc);
  Nmax = 400;
  
  printf("Ncell = %d \n",Ncell);
  
  // create arrays
  x = dvector(0,N-1);
  y = dvector(0,N-1);

  x_accum = dvector(0,N-1);
  y_accum = dvector(0,N-1);

  x_init = dvector(0,N-1);
  y_init = dvector(0,N-1);

  theta = dvector(0,N-1);

  vx = dvector(0,N-1);
  vy = dvector(0,N-1);

  vx_init = dvector(0,N-1);
  vy_init = dvector(0,N-1);


  cell_list_index = ivector(0,N-1);
  cell_list = imatrix(0,(Ncell*Ncell)-1,0,Nmax-1);


  msdx_accum_array = dvector(0,N_MSD);
  msdx_init_array = dvector(0,N);

  x_tracked =  dvector(0,N);

  cvvx_accum_array = dvector(0,N_MSD);
  cvvx_init_array = dvector(0,N);
  

  ratio = (double) Nsim/ (double) N_MSD;
  
  
  std::stringstream sstre;
  sstre << "outputEq.dat";
  const std::string tmpe = sstre.str();
  const char* cstre = tmpe.c_str();
  
  FILE * foute = fopen(cstre,"w");
  
  // perform simulation
  struct timeval tp;
  gettimeofday(&tp, NULL);
  long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
  printf("Start Simulation! \n");
  printf("step KE m\n");
  fprintf(foute, "#time KE m\n");
  
  char buf[80];
  sprintf(buf, "%s_%d.dat", "positions_step", Nsim-1);
  std::ifstream ifile(buf);
  //std::cout << "positions_step_"+s+".dat" << "\n";
  if (ifile){
    std::fstream myfile(buf, std::ios_base::in);
    float a;
    int i=0;
    while (myfile >> a)
    {
      if ((i%4)==0) x[(int)floor(i/4)]=a;
      else if ((i%4)==1) y[(int)floor(i/4)]=a;
      else if ((i%4)==2) vx[(int)floor(i/4)]=a;
      else vy[(int)floor(i/4)]=a;
      i++;
    }
    print_positions(0);
    
  }
  else {
    // initialize particles
    lattice_packing();
    
    // initialize velocities
    random_velocities();
    
    print_positions(0);

  



    eq_temp_sum = 0.0;
    eq_temp_count = 0;
  
    //equilibration:
    for (int i=0; i<Neq; i++) {
  
      md_step();
  
      if (i%100==0) {
        double m = calculate_order();
        //      double msd = calculate_msd();
  
        //printf("%f %f\n",i*dt, m);
        fprintf(foute,"%f %f\n",i*dt, m);
      }
  
    }
  
    fclose(foute);
  }
  
  
  create_cell_lists();

  initialize_msd_variables();
  

  std::stringstream sstr;
  sstr << "output.dat";
  const std::string tmp = sstr.str();
  const char* cstr = tmp.c_str();
  
  FILE * fout = fopen(cstr,"w");

  std::stringstream sstr_drs;
  sstr_drs << "displacements.dat";
  const std::string tmp_drs = sstr_drs.str();
  const char* cstr_drs = tmp_drs.c_str();
 
  FILE * fout_displ = fopen(cstr_drs,"w");
 
  std::stringstream sstr_corr_vv;
  sstr_corr_vv << "velocity_corrs.dat";
  const std::string tmp_corr_vv = sstr_corr_vv.str();
  const char* cstr_corr_vv = tmp_corr_vv.c_str();
 
  FILE * fout_corvv = fopen(cstr_corr_vv,"w");
 
  


  initialize_msd_variables();
  initialize_vel_auto_variables();
  double m_sum = 0.0;
  double m2_sum = 0.0;
  int m_count = 0;


  fprintf(foute, "#time KE m\n");

  //main run:
	 
  int pow2p0 = 0;
  int tick_tock = 100*pow2p0; //  pow(2.0,pow2p0);
  double msd_x, msd_y;
  double corr_vx, corr_vy;
	 
  for (int i=0; i<Nsim; i++) {

    
    if (i==tick_tock && MSDyn!=0)
    {

	    pow2p0++;
	
      tick_tock = 100*pow2p0; //  pow(2.0,pow2);
        
      printf("%d\n", tick_tock );
	
      msd_x = msd_y = 0;

      calculate_msd_2(msd_x, msd_y);

      printf("%f %f %f\n",i*dt, msd_x, msd_y);
      fprintf(fout_displ,"%f %f %f\n",i*dt, msd_x, msd_y);

      corr_vx = corr_vy = 0;
	
      calculate_vel_auto_2(corr_vx, corr_vy);


      fprintf(fout_corvv,"%f %f %f\n",i*dt, corr_vx, corr_vy);

      printf("%f %f %f\n",i*dt,  corr_vx, corr_vy);



    }



    md_step();
    
    //perform_blocking_operation(i, );
    if (MSDyn!=0) {
      for (int n=0; n<N; n++) {
        double delx, dely;
        int j0; // = i % N_blocks;
        int max_levels = (int) round((log(i)/log(N_blocks)));
      
        for (int k=0; k<(max_levels+1); k++) {
      
	        if (k==0) {
	          delx = vx[n];
	          dely = vy[n];
	        } else {
	          delx = blocking_sumx[N_blocks-1][k-1][n];
	          dely = blocking_sumy[N_blocks-1][k-1][n];
	        }
      
          int nblocks_to_k = 1 ; // (int) round(pow(N_blocks,k));
      
          for (int kk=0; kk<k; kk++)
            nblocks_to_k *= N_blocks;
      
            if (i % nblocks_to_k == 0 ) {
      
              j0 =  ((i) / nblocks_to_k-1 ) % N_blocks;  // was -1
      
              if (j0==0) {
                blocking_sumx[j0][k][n] = delx; //blocking_sum[N_blocks-1][k-1][n] ;
                blocking_sumy[j0][k][n] = dely;
              } else {
                blocking_sumx[j0][k][n] = blocking_sumx[j0-1][k][n] + delx; 
                blocking_sumy[j0][k][n] = blocking_sumy[j0-1][k][n] + dely; 
              }
	    
              double dr2x = blocking_sumx[j0][k][n]*blocking_sumx[j0][k][n]*dt*dt;
              double dr2y = blocking_sumy[j0][k][n]*blocking_sumy[j0][k][n]*dt*dt;
              msd_blocking[j0][k] +=  dr2x+dr2y;
      
              count_block[j0][k] += 1;
            }
        }	
      }
      
      int delta_t = i % N_MSD;
      
      if (delta_t == 0)
      {
        for (int j=0; j<N; j++)
        {
          msdx_init_array[j]=x[j];
          x_tracked[j] = x[j];
          cvvx_init_array[j]=vx[j];
        }
      }
      
      
      for (int j=0; j<N; j++)
      {
        msdx_accum_array[delta_t] += (x_tracked[j]-msdx_init_array[j])*(x_tracked[j]-msdx_init_array[j]);
        cvvx_accum_array[delta_t] += vx[j]*cvvx_init_array[j];
      }
    }

    
    if (i%freqeq==0) {
      double m = calculate_order();
      //double msd = calculate_msd();
      //      double fourth=calculate_fourth();

      //double vel_corr = calculate_vel_auto();


      //      double nonGaussian = 3.0/5.0*fourth/(msd*msd)-1.0;
      

      //printf("%f %f %f\n",i*dt, ke, m);//, msd, fourth, vel_corr, nonGaussian);

      fprintf(fout,"%f %f\n",i*dt, m);

      m_sum += m;
      m2_sum += m*m;
      m_count ++;
      //printf("%f %f\n",i*dt, m);
      //print_gr(i);
      //print_positions(i);
      
    }
    
    if ((i+1)%((int)Nsim/4)==0) {
      print_positions(i);
    }

  }
  
  if (MSDyn!=0) {
    std::stringstream sstr_blocks2;
    sstr_blocks2 << "output_msd_blocks2.dat";
    const std::string tmp_blocks2 = sstr_blocks2.str();
    const char* cstr_b2 = tmp_blocks2.c_str();
    
    FILE * fout_msds_blocks2 = fopen(cstr_b2,"w");
  
    for (int i=0; i<N_levels; i++)
      for (int j=0; j<N_blocks; j++)
        {
    double time_here = (j+1) * (pow(N_blocks,i))*dt;
    double msd_here = msd_blocking[j][i]/count_block[j][i];
  
    fprintf(fout_msds_blocks2,"%d %d %f %f %.10f\n",i, j, count_block[j][i], time_here, msd_here);
  
        }
  
    fclose(fout_msds_blocks2);
  }

  double m_mean = m_sum/m_count;
  double m_var = m2_sum/m_count - pow(m_sum/m_count,2);
  printf("m = %f +/- %f , count = %d \n", m_mean, sqrt(m_var), m_count);
  

  struct timeval tpend;
  gettimeofday(&tpend, NULL);
  long int msend = tpend.tv_sec * 1000 + tpend.tv_usec / 1000;
  long int difftime = msend - ms;
  printf("Finished Simulation! Time: %d ms \n",difftime);

  fclose(fout);

  fclose(fout_displ);
  fclose(fout_corvv);

  output_MSD_VCC();


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
  int K = (int) (pow(N,1./2.)); //number of particles per row
  if (pow(K,2)<N){ K++;} //correct K if not all the particles will fit
  double a = L/((double) K); //spacing between particles  
  
  for (int i=0; i<K; i++) {
    for (int j=0; j<K; j++) {
      int index = i*K+j;
      if(index < N){
        x[index] = i*a-0.5*L;
        y[index] = j*a-0.5*L;
      }
    }
  }
  
  printf("Initialized...  Density %f, Noise %f, Number of particles %d, Box length %f\n",\
	 den,eta,N,L);
}

/*********************************************************************************************************************************/

void random_velocities() {


  //assign each angle randomly in a uniform distribution 
  for (int i=0; i<N; i++) {
    theta[i] = random_mars->uniform()*2.*PI;
    vx[i] = v*cos(theta[i]);
    vy[i] = v*sin(theta[i]);
  }

}


/*********************************************************************************************************************************/

//update the positions of all particles 
//(with one iteration of the verlet integrator.) 
void md_step(){

  for (int i=0; i<N; i++) {
    x[i] = x[i] + dt*vx[i];
    y[i] = y[i] + dt*vy[i];

    x_accum[i] = x_accum[i] + dt*vx[i];
    y_accum[i] = y_accum[i] + dt*vy[i];
    
    x_tracked[i] = x_tracked[i]+dt*vx[i];

    //apply periodic boundary conditions
    pbc(x[i],y[i]);
  }

  //recalculate cell list
  create_cell_lists();

  //reculculate angles, now that particles have moved
  update_angles();

}

/*********************************************************************************************************************************/

// calculate and store the force on every particle
//(from the derivative of the interaction potential,
//summed over neighbouring particles)
void update_angles() {
  
  int nint;
  double r_cut = 1;
  double r_cut2 = r_cut*r_cut;
  double * sin_theta_new;
  double * cos_theta_new;
  
  sin_theta_new = dvector(0,N-1);
  cos_theta_new = dvector(0,N-1);

  //loop through all particles
  for (int i=0; i<N; i++) {
    //initilize arrays to hold forces
    sin_theta_new[i] = 0.0;
    cos_theta_new[i] = 0.0;
    nint = 0;

    //determine x,y,z cell index
    int cell_list_indexi = cell_list_index[i];
    int cell_list_y  = cell_list_indexi % Ncell;
    int cell_list_x  = ((cell_list_indexi - cell_list_y)/Ncell) % Ncell;

    //loop through neighboring cells
    for (int b=cell_list_x - 1; b<= cell_list_x + 1; b++) {
      for (int c=cell_list_y - 1; c<= cell_list_y + 1; c++) {
        //account for periodic boundary conditions
        int bloc = b;
        int cloc = c;
        if (bloc < 0) bloc += Ncell;
        if (bloc >= Ncell) bloc -= Ncell;
        if (cloc < 0) cloc += Ncell;
        if (cloc >= Ncell) cloc -= Ncell;

        //loop through all particles in the cell
        for (int j=0; j< Nmax; j++) {
	  int jloc = cell_list[bloc*Ncell + cloc][j];
	  if (jloc != -1) {
            //calculate force between particles
            //(and add it to the sum for particle i)
	    double dxloc = x[i] - x[jloc];
            double dyloc = y[i] - y[jloc];
	    pbc(dxloc, dyloc);
	    double r2 = dxloc*dxloc + dyloc*dyloc;
            if(r2<r_cut2){
              sin_theta_new[i] += vy[jloc];
              cos_theta_new[i] += vx[jloc];
              nint++;
            }
          }
        }

      }
    }
    
    //sin_theta_new[i]/=(nint*v);
    //cos_theta_new[i]/=(nint*v);

  }

  //update theta
  for (int i=0; i<N; i++) {
    theta[i]=atan2(sin_theta_new[i],cos_theta_new[i])+eta*sqrt(dt)*2*PI*(random_mars->uniform()-0.5);//PI*random_mars->gaussian();
    vx[i] = v*cos(theta[i]);
    vy[i] = v*sin(theta[i]);
  }
  
  free_dvector(sin_theta_new,0,N-1);
  free_dvector(cos_theta_new,0,N-1);
}

/*********************************************************************************************************************************/

//calculate the kinetic energy of the system
double calculate_kinetic(){
  double ke_sum = 0.0;
  for (int i=0; i<N; i++) {
    ke_sum += vx[i]*vx[i]+vy[i]*vy[i];
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
    
    for (int j=0; j<N; j++) {
      if (j!=i) {
	double dxloc = xi - x[j];
	double dyloc = yi - y[j];
	pbc(dxloc, dyloc);
	
	double dr = sqrt(dxloc*dxloc + dyloc*dyloc);
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
    double V = 4.* PI * ( pow(i*drhisto + drhisto,2)- pow(i*drhisto,2) )/2.;
    fprintf(out,"%f %f\n",i*drhisto,histo[i]*L*L/((double) N * N *V));
    
  }
  
  fclose(out);
}

/*********************************************************************************************************************************/

// print particle positions
void print_positions(int step){
  
  std::stringstream sstr;
  sstr << "positions_step_" << step << ".dat";
  const std::string tmp = sstr.str();
  const char* cstr = tmp.c_str();
  
  FILE * out = fopen(cstr,"w");
  
  for (int i=0; i<N; i++) {
    double xloc = x[i];
    double yloc = y[i];
    pbc(xloc, yloc);
    fprintf(out,"%f %f %f %f\n",xloc,yloc,vx[i],vy[i]);
  }
  
  fclose(out);
    
}

/*********************************************************************************************************************************/

void create_cell_lists(){



  // initialize cell lists
  for (int b=0; b<Ncell; b++) {
    for (int c=0; c<Ncell; c++) {
      for (int j=0; j< Nmax; j++) {
	cell_list[b*Ncell + c][j] = -1;
      }
    }
  }


  

  // find cell lists for particles
  for (int i=0; i<N; i++) {
    int xint = (int) ( Ncell* (x[i]/L + 0.5) );
    int yint = (int) ( Ncell* (y[i]/L + 0.5) );

    if(xint<0) xint=0;
    if(yint<0) yint=0;
    if(xint>(Ncell-1)) xint=Ncell-1;
    if(yint>(Ncell-1)) yint=Ncell-1;

    cell_list_index[i] = xint * Ncell + yint;


    bool particle_stored = false;
    for (int j=0; j< Nmax; j++) {
      if (cell_list[xint*Ncell + yint][j] == -1) {
	cell_list[xint*Ncell + yint][j] = i;
        particle_stored = true;
	break;
      }
    }
    if (!particle_stored){ 
      //printf("ERROR: cell list storage exceeded!!! Nmax must be increased! \n");
    }
  }

}


/*********************************************************************************************************************************/
// consider periodic boundary conditions
void pbc(double &x, double &y){
  if (x >= L/2.0) x -= L;
  if (x < -L/2.0) x += L;
  if (y >= L/2.0) y -= L;
  if (y < -L/2.0) y += L;
}
    
double calculate_msd()
{
  double msd_now = 0;

  for (int i=0; i<N; i++)
    {
      msd_now += 
	(x_accum[i]-x_init[i])*(x_accum[i]-x_init[i])+
	(y_accum[i]-y_init[i])*(y_accum[i]-y_init[i]);
    }

  msd_now = msd_now / (float)N;

  return msd_now;
}

double calculate_fourth()
{
  double fourth = 0;
  
  double dx,dy;

  for (int i=0; i<N; i++)
    {
      dx = (x_accum[i]-x_init[i]);
      dy = (y_accum[i]-y_init[i]);
      fourth += dx*dx*dx*dx+dy*dy*dy*dy;

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
      
    }
}

void initialize_vel_auto_variables()
{
  for (int i=0; i<N; i++)
    {
      vx_init[i] = vx[i];
      vy_init[i] = vy[i];
      
    }
}

double calculate_vel_auto()
{

  double corr = 0;

  for (int i = 0; i<N; i++ )
    {      
      corr += vx[i]*vx_init[i] + vy[i]*vy_init[i];
    }

  corr = corr / (float) N;

  return corr;

}

void calculate_msd_2(double &msd_x, double &msd_y)
{

  for (int i=0; i<N; i++)
    {
      msd_x += (x_accum[i] - x_init[i])*(x_accum[i] - x_init[i]);
      msd_y += (y_accum[i] - y_init[i])*(y_accum[i] - y_init[i]);
    }

  msd_x = msd_x/(float)N;
  msd_y = msd_y/(float)N;

}

void calculate_vel_auto_2(double &corr_vx,double &corr_vy)
{
  for (int i=0; i<N; i++)
    {
      corr_vx += vx[i] * vx_init[i];
      corr_vy += vy[i] * vy_init[i];
    }

  corr_vx = corr_vx/(float)N;
  corr_vy = corr_vy/(float)N;
}

double calculate_order()
{
  double m_x=0;
  double m_y=0;
  double m=0;
  
  for (int i=0; i<N; i++) 
  {
    m_x+=vx[i];
    m_y+=vy[i];
  }
  
  m=sqrt(m_x*m_x+m_y*m_y)/v;
  
  return m/((double)N);
}

