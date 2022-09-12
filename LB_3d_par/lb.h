///////////////////////////////////////////////
//libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>
#include <math.h>

typedef enum {false, true} bool;
typedef char *string;

// The dimension of the lattice space
#define 	 SDIM 		 	3

// The number of different links: D3Q19 Model
#define 	 NDIM 			19

//The number of ples that have positive or neg. component to a given direction
#define   	 COL			5

// The sound velocity = 1/sqrt(3)
//#define 	 CS 		 	0.5773502692

// The sound velocity squared
#define 	 CS2 		 	(1.0/3.0)
#define 	 CS4             	(CS2*CS2)

// The eq. coeff. has been copied from Y.H. Qian et al. Europhys. Lett. 17 (6) 479
//#define  	 t0		 	(1.0/3.0)
//#define  	 t1 		 	(1.0/18.0)
//#define   	 t2		 	(1.0/36.0)


///////////////////////////////////////////////
//structs

//struct contends macroscopic information
typedef struct{
	int t_max; //Maximum number of iterations
	double density; //Fluid density per link
	double accel; //Accelleration
	double omega; //Relaxation parameter
	double r_rey; //Linear dimension for Reynolds number
} s_properties;

//lattice structure
typedef struct {
	long int lx; //nodes number in axis x
	long int ly; //nodes number in axis y
	long int lz; //nodes number in axis z
	long int lx_total;
	long int ly_total;
	long int lz_total;
	long int lx_first; //initial position of the sub-lattice
	long int lx_last; //final position of the sub-lattice
	long int ly_first; //initial position of the sub-lattice
	long int ly_last; //final position of the sub-lattice
	long int lz_first; //initial position of the sub-lattice
	long int lz_last; //final position of the sub-lattice

	int d; //spatial dimensions
	int n; //lattice dimension elements
	bool ***obst; //obstacle array lx * ly * lz
	double ****node; //n-speed lattice  lx * ly * lz * n
	double ****temp; //temporarely storage of fluid densities
} s_lattice;

//assistent structure for preprocessing
typedef struct {
	int **e; //edges
	//the coefficients for equilibrium distrution. A is for the constant term, B is for e*u-term, C is for (e*u)^2-term and D is for u^2-term
	double *A, *B, *C, *D;
	//The six tables that tell the numbers of ples that have positive or negative x, y or z component. Each of these tables have 5 (CLNBR) components
	int *pos_x, *pos_y, *pos_z;
	int *neg_x, *neg_y, *neg_z;
	// tables that give the population numbers when the direction is known: 0 = -component, 1 = no component 2 = +component
	int ***points;
} s_aux;

//structure for partitioning and communications process
typedef struct {
	MPI_Comm grid_comm;
	int myid;
	int myid_x;
	int myid_y;
	int myid_z;
	int proc_x; // number of processors at x-axis
	int proc_y; // number of processors at y-axis
	int proc_z; // number of processors at z-axis
} s_comm;


///////////////////////////////////////////////
//Functions

s_comm *initialize_comm(string dim1, string dim2, string dim3);

void define_sublattice(s_comm *comm, s_lattice *l);

void send_pos_x(s_aux *aux, s_lattice *l, s_comm *comm);

void send_neg_x(s_aux *aux, s_lattice *l, s_comm *comm);

void send_pos_y(s_aux *aux, s_lattice *l, s_comm *comm);

void send_neg_y(s_aux *aux, s_lattice *l, s_comm *comm);

void send_pos_z(s_aux *aux, s_lattice *l, s_comm *comm);

void send_neg_z(s_aux *aux, s_lattice *l, s_comm *comm);

void copy_x(s_lattice *l, s_aux *aux); 

void copy_y(s_lattice *l, s_aux *aux);

void copy_z(s_lattice *l, s_aux *aux);

void corner_send_x(s_lattice *l, s_aux *aux, s_comm *comm);

void corner_send_y(s_lattice *l, s_aux *aux, s_comm *comm);

void corner_send_z(s_lattice *l, s_aux *aux, s_comm *comm);

void corner_receive_x(s_lattice *l, s_aux *aux, s_comm *comm);

void corner_receive_y(s_lattice *l, s_aux *aux, s_comm *comm);

void corner_receive_z(s_lattice *l, s_aux *aux, s_comm *comm);

void corner_z(s_lattice *l, s_aux *aux, s_comm *comm);

void corner_y(s_lattice *l, s_aux *aux, s_comm *comm);

void corner_x(s_lattice *l, s_aux *aux, s_comm *comm);

void sincronization(s_aux *aux, s_lattice *l, s_comm *comm);

double check_density_par(s_lattice *l, s_comm *comm, int time);

double calc_velocity_par(s_lattice *l, s_comm *comm, int time);

s_properties *read_parametrs(string file);

void alloc_lattice(s_lattice *l);

s_lattice *create_lattice(string file, s_comm *comm);

s_aux *init_constants(s_lattice *l, s_properties *p);

void init_density(s_lattice * l, double *A);

double check_density(s_lattice *l);

void redistribute(s_lattice *l, double accel, double density);

void propagate(s_lattice *l);

void bounceback(s_lattice *l);

void relaxation(s_lattice *l, double density, double omega, s_aux *aux);

double *calc_velocity(s_lattice *l, int time);

void write_results(string file, s_lattice *l, double density);

void comp_rey(double vel, s_properties *p, int time, double execution, double communication, string dir, int proc_x, int proc_y, int proc_z);

double crono();

void parallel_corner_y(s_lattice *l, s_aux *aux, s_comm *comm);

void parallel_corner_x(s_lattice *l, s_aux *aux, s_comm *comm);

void parallel_corner_z(s_lattice *l, s_aux *aux, s_comm *comm);

void corner_send_y_neg_z_neg(int y_neg_z_neg, double *buffer_x, s_lattice *l, s_comm *comm);

void corner_send_y_pos_z_neg(int y_pos_z_neg, double *buffer_x, s_lattice *l, s_comm *comm);

void corner_send_y_neg_z_pos(int y_neg_z_pos, double *buffer_x, s_lattice *l, s_comm *comm);

void corner_send_y_pos_z_pos(int y_pos_z_pos, double *buffer_x, s_lattice *l, s_comm *comm);


void corner_send_x_neg_z_neg(int x_neg_z_neg, double *buffer_y, s_lattice *l, s_comm *comm);

void corner_send_x_pos_z_neg(int x_pos_z_neg, double *buffer_y, s_lattice *l, s_comm *comm);

void corner_send_x_neg_z_pos(int x_neg_z_pos, double *buffer_y, s_lattice *l, s_comm *comm);

void corner_send_x_pos_z_pos(int x_pos_z_pos, double *buffer_y, s_lattice *l, s_comm *comm);


void corner_send_x_neg_y_neg(int x_neg_y_neg, double *buffer_z, s_lattice *l, s_comm *comm);

void corner_send_x_pos_y_neg(int x_pos_y_neg, double *buffer_z, s_lattice *l, s_comm *comm);

void corner_send_x_neg_y_pos(int x_neg_y_pos, double *buffer_z, s_lattice *l, s_comm *comm);

void corner_send_x_pos_y_pos(int x_pos_y_pos, double *buffer_z, s_lattice *l, s_comm *comm);


void corner_recv_y_neg_z_neg(double *buffer_x, s_lattice *l, s_comm *comm);

void corner_recv_y_pos_z_neg(double *buffer_x, s_lattice *l, s_comm *comm);

void corner_recv_y_neg_z_pos(double *buffer_x, s_lattice *l, s_comm *comm);

void corner_recv_y_pos_z_pos(double *buffer_x, s_lattice *l, s_comm *comm);

void corner_recv_x_neg_z_neg(double *buffer_y, s_lattice *l, s_comm *comm);

void corner_recv_x_pos_z_neg(double *buffer_y, s_lattice *l, s_comm *comm);

void corner_recv_x_neg_z_pos(double *buffer_y, s_lattice *l, s_comm *comm);

void corner_recv_x_pos_z_pos(double *buffer_y, s_lattice *l, s_comm *comm);

void corner_recv_x_neg_y_neg(double *buffer_z, s_lattice *l, s_comm *comm);

void corner_recv_x_pos_y_neg(double *buffer_z, s_lattice *l, s_comm *comm);

void corner_recv_x_neg_y_pos(double *buffer_z, s_lattice *l, s_comm *comm);

void corner_recv_x_pos_y_pos(double *buffer_z, s_lattice *l, s_comm *comm);


void corner_y_neg_z_neg(double *buffer_x, s_lattice *l, s_comm *comm);

void corner_y_pos_z_neg(double *buffer_x, s_lattice *l, s_comm *comm);

void corner_y_neg_z_pos(double *buffer_x, s_lattice *l, s_comm *comm);

void corner_y_pos_z_pos(double *buffer_x, s_lattice *l, s_comm *comm);


void corner_x_neg_z_neg(double *buffer_y, s_lattice *l, s_comm *comm);

void corner_x_pos_z_neg(double *buffer_y, s_lattice *l, s_comm *comm);

void corner_x_neg_z_pos(double *buffer_y, s_lattice *l, s_comm *comm);

void corner_x_pos_z_pos(double *buffer_y, s_lattice *l, s_comm *comm);


void corner_x_neg_y_neg(double *buffer_z, s_lattice *l, s_comm *comm);

void corner_x_pos_y_neg(double *buffer_z, s_lattice *l, s_comm *comm);

void corner_x_neg_y_pos(double *buffer_z, s_lattice *l, s_comm *comm);

void corner_x_pos_y_pos(double *buffer_z, s_lattice *l, s_comm *comm);


