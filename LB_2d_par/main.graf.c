#include "lb.h"

///////////////////////////////////////////////
int main(int argc, char **argv) {
	
	//Begin MPI initialization
	MPI_Init(&argc,&argv);
	
	//MPI initialization structures
	s_comm *comm = initialize_comm(argv[4], argv[5]);

	//Iteration counter
	int time;
	
	//Tolerance
	//double tolerance = 0.00001;
	
	//Density
	//double den;
	
	//Time counters
	double start_time;
	double start_comm;
	double time_exec;
	double time_comm;

	//Average velocity
	double vel = 0;
	//double vel_b = 0;

	//Input structure
	s_properties *properties;

	//Lattice structure
	s_lattice *lattice;

	//startup information message	
	/* printf("Lattice Boltzmann Method\n");
	printf("Claudio Schepke\n");
	printf("Instituto de Informática - UFRGS\n");
	printf("Date: 2006, January 09\n\n"); */
	
	//Checking arguments
	/* if (comm->myid == 0 && argc != 7) {
		fprintf(stderr, "Usage: %s [file_configuration] [file_colision] [file_plot] [proc_dimension_x] [proc_dimension_y] [Directory]\n\n", argv[0]);	
		exit(1);
	}*/
	
	//Read parameter file
	//properties->t_max
	//properties->density
	//properties->accel
	//properties->omega
	//properties->r_rey
	properties = (s_properties*) read_parametrs(argv[1]);

	//read obstacle file
	//<x> <y> <n directions> <number of obstacles> 
	//x- and y-coordinates of any obstacles
	//wall boundaries are also defined here by adding single obstacles
	//Define to each node the first and the last position of the sublattice
	lattice = (s_lattice*) create_lattice(argv[2], comm);
	
	//Initialize Density
	init_density(lattice, properties->density);

	//Begin timer
	if (comm->myid == (comm->proc_x * comm->proc_y)/2)
		start_time = MPI_Wtime();

	//Begin of the main loop
	//do {
	for (time = 0; time < properties->t_max; time++) {

		//den = check_density_par(lattice, comm);
		//if(comm->myid == 0 && !(time%(properties->t_max/1))) {
			fprintf("%d %f\n", time, den);
		//}
		
		if (comm->myid_x == 0)
			redistribute(lattice, properties->accel, properties->density);
		
		propagate(lattice);

		//Communications
		start_comm = MPI_Wtime();
		sincronization(lattice, comm); 
		time_comm = time_comm + MPI_Wtime() - start_comm;
		
		bounceback(lattice);
		
		relaxation(lattice, properties->density, properties->omega);
		
		//Calcule velocity
		//vel_b = vel;
		//vel = calc_velocity_par(lattice, comm);
		//if (comm->myid == 0)
		//	printf("%d %f\n", time, vel);
		
		//time ++;
	//} while(sqrt((vel - vel_b) * (vel - vel_b)) > tolerance);
	
		if (!(time%10000)) {
			string tt = (string) malloc(100*sizeof(string));
			string it = (string) malloc(100*sizeof(string));
			strcpy(tt, argv[3]);
			sprintf(it, "%d", time);
			strcat(tt, it);
			write_results(argv[6], tt, lattice, properties->density);
			free(it);
			free(tt);
		}
	}
	
	//End counter
	if (comm->myid == (comm->proc_x * comm->proc_y)/2)
		time_exec = MPI_Wtime();

	//Calcule velocity
	vel = calc_velocity_par(lattice, comm);
	
	//Macroscopic properties and time resume
	if (comm->myid == (comm->proc_x * comm->proc_y)/2) {
		comp_rey(vel, properties, time, time_exec - start_time, time_comm, argv[6], comm->proc_x, comm->proc_y);
		write_results(argv[6], argv[3], lattice, properties->density);
	} 
		

	MPI_Finalize();
	//printf("End of the execution\n\n");
	return 1;
}
