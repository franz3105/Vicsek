

double L;
double den;
double mu;
double eta;
double delta;
double v;
int freq;
int order;

int N;				// number of particles
double Tinit;

double * x;			// x, y and z positions of the particles
double * y;
double * z;

double * x_accum;               // x,y,z positions but not counting periodic boundaries
double * y_accum;
double * z_accum;

double x_save;
double y_save;
double z_save;

double * x_init;              //x,y,z store initial postitions
double * y_init;
double * z_init;

double * vx;			// x, y and z velocities of the particles
double * vy;
double * vz;

double * vx_init;
double * vy_init;
double * vz_init;

double * theta;
double *alpha;
double * fx;			// x, y and z components of force on the particles
double * fy;
double * fz;
double * SF;
double * SFimg;
double * realsum;
double * imgsum;
int * cell_list_index;
int **cell_list;


double * msdx_accum_array;
double * msdx_init_array;

double * x_tracked;

double * cvvx_accum_array;
double * cvvx_init_array;

double ratio;

double rc;			// constants for cell list calculation
int Ncell;
int Nmax;

double dt;		// step_size for MC moves


int Neq;			// Number of equillibrations steps, number of simulation steps
int Nsim;

int NSF;
int countSF;

int N_MSD;
int count_MSD;

double press_sum;
double virial;
int press_count;

double eq_temp_sum;
int eq_temp_count, freqEqTemp;

void initialize_msd_variables();
void initialize_SF();
void sample_SF();
void sample_SF_fast();
void output_SF();

void adjust_temp(int i, double ke);

void lattice_packing();
void random_velocities();
void initialize_oldpos();

void pbc(double &x, double &y, double &z);

double calculate_potential();
void md_step();
void calculate_forces();
double calculate_kinetic();
double calculate_msd();
void calculate_msd_2(double &msd_x, double &msd_y, double &msd_z);
double calculate_fourth();
double calculate_order();

double calculate_vel_auto();
void calculate_vel_auto_2(double &corr_vx,double &corr_vy, double &corr_vz);
void initialize_vel_auto_variables();
double gaussnormal(double m, double s);
void sum_pressure();

void print_gr(int step);
void print_positions(int i);

void create_cell_lists();

void output_MSD_VCC();
void update_angles();

class RanMars *random_mars;

