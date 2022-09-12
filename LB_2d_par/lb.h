/* The densities are numbered as follows:
 *               6   2   5
 *                 \ | /
 *               3 - 0 - 1
 *                 / | \
 *               7   4   8
 *
 * The lattice nodes are numbered as follows:
 *      ^
 *      |
 *      y
 *           :    :    :
 *
 *      3    *    *    *  ..
 *
 *           *    *    *  ..
 *
 *      1    *    *    *  ..
 *                            x ->
 *           1    2    3 
 *      
*/


///////////////////////////////////////////////
//libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <sys/time.h>
#include<math.h>

typedef enum {false, true} bool;
typedef char *string;

#define COL 3
//int numprocs; //number of processors
//int myid; // processor identification


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
	long int lx_first; //initial position of the sub-lattice
	long int lx_last; //final position of the sub-lattice
	long int ly_first; //initial position of the sub-lattice
	long int ly_last; //final position of the sub-lattice
	long int lx; //vertex number in axis x
	long int ly; //vertex number in axis y
	long int lx_total;
	long int ly_total;
	int n; //lattice dimension - number of directions
	int d; //spatial dimensions
	bool **obst; //obstacle Array lx * ly
	double ***node; //n-speed lattice  n * lx * ly
	double ***temp; //temporarely storage of fluid densities
} s_lattice;


//structure for partitioning and communications process
typedef struct {
	MPI_Comm grid_comm;
	int myid;
	int myid_x;
	int myid_y;
	int proc_x; // number of processors at x-axis
	int proc_y; // number of processors at y-axis
} s_comm;

	
///////////////////////////////////////////////
//Functions

s_properties* read_parametrs(string file);

void alloc_lattice(s_lattice *l);

s_lattice* read_obstacles(string file, s_comm *comm);

void init_density(s_lattice * l, double density);

double check_density(s_lattice *l);

void redistribute(s_lattice *l, double accel, double density);

void propagate(s_lattice *l);

void bounceback(s_lattice *l);

void relaxation(s_lattice *l, double density, double omega);

double *calc_velocity(s_lattice *l);

void write_results(string dir, string file, s_lattice *l, double density);

void comp_rey(double vel, s_properties *p, int time, double execution, double communication, string dir, int proc_x, int proc_y);

double crono();



s_lattice *create_lattice(string file, s_comm *comm);

void define_sublattice(s_comm *comm, s_lattice *l);

s_comm *initialize_comm(string dim1, string dim2);

double check_density_par(s_lattice *l, s_comm *comm);

double calc_velocity_par(s_lattice *l, s_comm *comm);



double *alloc_border(int buffer_size);

void send_border(double *buffer, int buffer_size, int *coord, int label, MPI_Comm grid_comm);

void send_corner(double *buffer, int buffer_size, int *coord, int label, MPI_Comm grid_comm);


void send_pos_x(s_lattice *l, s_comm *comm);
void send_neg_x(s_lattice *l, s_comm *comm);
void send_pos_y(s_lattice *l, s_comm *comm);
void send_neg_y(s_lattice *l, s_comm *comm);

void recv_pos_x(s_lattice *l, s_comm *comm);
void recv_neg_x(s_lattice *l, s_comm *comm);
void recv_pos_y(s_lattice *l, s_comm *comm);
void recv_neg_y(s_lattice *l, s_comm *comm);

void copy_x(s_lattice *l);
void copy_y(s_lattice *l);
void sincronization(s_lattice *l, s_comm *comm);

void copy_corners(s_lattice *l, s_comm *comm);

void corner_max_max(s_lattice *l, s_comm *comm);
void corner_max_min(s_lattice *l, s_comm *comm);
void corner_min_max(s_lattice *l, s_comm *comm);
void corner_min_min(s_lattice *l, s_comm *comm);

