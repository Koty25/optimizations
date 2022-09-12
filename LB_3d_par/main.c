#include "lb.h"

///////////////////////////////////////////////
int main(int argc, char **argv) {
	
	//MPI cartesian communications	
	MPI_Init(&argc, &argv);

	//MPI initialization structures
	s_comm *comm = initialize_comm(argv[4], argv[5], argv[6]); 
	//Parameters

	//Iteration counter
	int time = 0;

	//Tolerance
	double tolerance = 0.0001;

	//Density
	//double den;
	
	//Execution Time
	double execution_time;
	double communication_time = 0;
	double iter_time;

	//Average velocity
	double vel = 0, vel_b = 0;

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
	if (comm->myid == 0 && argc != 8) {
		fprintf(stderr, "Usage: %s [file_configuration] [file_colision] [file_plot] [proc_dimension_x] [proc_dimension_y] [proc_dimension_z] [Directory]\n\n", argv[0]);
		exit(1);
	}
	
	//Begin initialization
	
	//Read parameter file
	//properties->t_max
	//properties->density
	//properties->accel
	//properties->omega
	//properties->r_rey
	properties = (s_properties*) read_parametrs(argv[1]);

	//read obstacle file
	//<x> <y> <z> <d dimensions> <n directions> <number of obstacles> 
	//x-, y-, and z-coordinates of any obstacles
	//wall boundaries are also defined here by adding single obstacles

	lattice = (s_lattice*) create_lattice(argv[2], comm);
	
	//alloc memory
	s_aux *aux = init_constants(lattice, properties);
	
	init_density(lattice, aux->A);
	
	//Execution Time
	execution_time = crono();
	//printf("Passou\n");
		

	//Begin of the main loop
	//for (time = 0; time < properties->t_max; time++) {
	do {
	//while (time < 5) {
		//den = check_density_par(lattice, comm, time);
		//if(comm->myid == 0) {
		//	printf("%f\n", den);
		//}
	
		redistribute(lattice, properties->accel, properties->density);

		propagate(lattice);
	
		iter_time = crono();
		sincronization(aux, lattice, comm);
		iter_time = crono() - iter_time;
		communication_time += iter_time;

		bounceback(lattice);

		relaxation(lattice, properties->density, properties->omega, aux);

		vel_b = vel;
		vel = calc_velocity_par(lattice, comm, time);
		//if (comm->myid == 0)
		//	printf("%d %f\n", time, vel);

		time++;

		//if (time%10000)
		//	write_results(strcat(argv[3], time), lattice, properties->density);
		
	} while(sqrt((vel - vel_b) * (vel - vel_b)) > tolerance);
	//}
		
	execution_time = crono() - execution_time;
	//if (comm->myid == 0)
	//	printf("Execution %d %f\n", time, execution_time);
		
	if (comm->myid == 0)
		comp_rey(vel, properties, time, execution_time, communication_time, argv[7], comm->proc_x, comm->proc_y, comm->proc_z);

	if (comm->proc_x == 1 && comm->proc_y == 1 && comm->proc_z == 1)
		write_results(argv[3], lattice, properties->density);

	//printf("End of the execution\n\n");
	MPI_Finalize();
	return true;
}



/*
if (comm->myid == 0) {
	int x = 0, y = 0, z, n;
	for (x = 0; x < lattice->lx; x++) {
		for (y = 0; y < lattice->ly; y++) {
			for (z = 0; z < lattice->lz; z++) {
				//for (n = 0; n < lattice->n; n++) {
					//printf("%d %d %d %d\n", lattice->obst[x][y][z], x, y, z);
					printf("%d ", lattice->obst[x][y][z]);
				//}
			} printf("\n");
		}printf("\n");
	}
}*/

/*		int x, y, z;
		for (x = 0; x < lattice->lx; x++)
			for (y = 0; y < lattice->ly; y++)
				for (z = 0; z < lattice->lz; z++)
					printf("%d %d %d %d %f\n", x, y, z, lattice->obst[x][y][z], lattice->node[x][y][z][2]); 
*/

