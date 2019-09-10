

double L;
double den;

int N;				// number of particles
double v;
double eta;

int N_blocks; // number of blocks
int N_levels; // nimber of levels

double * x;			// x, y and z positions of the particles
double * y;

double * x_accum;               // x,y,z positions but not counting periodic boundaries
double * y_accum;

double * x_init;              //x,y,z store initial postitions
double * y_init; 

double * theta;

double * vx;			// x, y and z velocities of the particles
double * vy;

double * vx_init;
double * vy_init;

double * fx;			// x, y and z components of force on the particles
double * fy;
int * cell_list_index;
int **cell_list;

float *** blocking_sumx;
float *** blocking_sumy;

double ** msd_blocking;
double ** count_block;

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

int N_MSD;
int count_MSD;


double eq_temp_sum;
int eq_temp_count, freqEqTemp;

void initialize_msd_variables();

void adjust_temp(int i, double ke);

void lattice_packing();
void random_velocities();
void initialize_oldpos();

void pbc(double &x, double &y);

void md_step();
void update_angles();
double calculate_kinetic();
double calculate_msd();
void calculate_msd_2(double &msd_x, double &msd_y);
double calculate_fourth();

double calculate_vel_auto();
void calculate_vel_auto_2(double &corr_vx,double &corr_vy);
void initialize_vel_auto_variables();

double calculate_order();

void print_gr(int step);
void print_positions(int step);

void create_cell_lists();

void output_MSD_VCC();

void perform_blocking_operation(int time_step, double );

class RanMars *random_mars;
