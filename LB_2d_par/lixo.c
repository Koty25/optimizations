#include "lb.h"

//////////////////////////////////////////////////////////////////////
// comm initialization and sublattice definition
//////////////////////////////////////////////////////////////////////


//Initialize MPI communications structure
s_comm *initialize_comm(string dim1, string dim2) {
	//Initialize comm structure
	s_comm *comm = (s_comm *) malloc(sizeof(s_comm));
	//Check for memory utilization
	if (comm == NULL) {
		printf("Error in memory allocation in function initialize_comm\n");
		exit(-10);
	}
	
	//Local variables
	int dimensions[2];
	int wrap_around[2];
	int coordinates[2];
	int reorder = false;
	
	//Get dimensions
	dimensions[0] = atoi(dim1);
	dimensions[1] = atoi(dim2);

	//Define wrap_around type
	wrap_around[0] = true;
	wrap_around[1] = true;
	
	//Create cartesian communication structure
	MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, reorder, &comm->grid_comm);
	
	//Define rank an coordinates for each processor
	MPI_Comm_rank(comm->grid_comm, &comm->myid);
	MPI_Cart_coords(comm->grid_comm, comm->myid, 2, coordinates);
	
	//Save number of processors for each dimension in comm structure
	comm->proc_x = dimensions[0];
	comm->proc_y = dimensions[1];

	//Save rank identification for each dimension in comm structure
	comm->myid_x = coordinates[0];
	comm->myid_y = coordinates[1];
	
	//Return the structure
	return comm;
}


//Define subdimensions of the lattice
void define_sublattice(s_comm *comm, s_lattice *l) {
	//Initialize variables
	int size_x = l->lx_total/comm->proc_x;
	int size_y = l->ly_total/comm->proc_y;
	int remain_x = l->lx_total%comm->proc_x;
	int remain_y = l->ly_total%comm->proc_y;
	
	//x-axis
	if(comm->myid_x > comm->proc_x - 1 - remain_x) {
		size_x += 1;
		l->lx_last = l->lx_total - (comm->proc_x - 1 - comm->myid_x)*size_x;
		l->lx_first = l->lx_last - size_x;
	} 
	else {
		l->lx_first = size_x * comm->myid_x;
		l->lx_last = size_x * (comm->myid_x + 1);
	}
	
	//y-axis
	if(comm->myid_y > comm->proc_y - 1 - remain_y) {
		size_y += 1;
		l->ly_last = l->ly_total - (comm->proc_y - 1 - comm->myid_y)*size_y;
		l->ly_first = l->ly_last - size_y;
	}
	else {
		l->ly_first = size_y * comm->myid_y;
		l->ly_last = size_y * (comm->myid_y + 1);
	}
	
	//Define the real large size of each sublattice
	//Overlapping borders
	l->lx = l->lx_last - l->lx_first + 2;
	l->ly = l->ly_last - l->ly_first + 2; 
}


//////////////////////////////////////////////////////////////////////
// Parallel Velocity and Density
//////////////////////////////////////////////////////////////////////

//Check the density for the parallel implementation
double check_density_par(s_lattice *l, s_comm *comm) {
	//local variables
	double dens = check_density(l);
	double buf;
	int i;
	int procs = comm->proc_x * comm->proc_y - 1;
	MPI_Status s;

	if (comm->myid == 0) {
		for(i = 0; i < procs; i++) {
			MPI_Recv(&buf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 30, MPI_COMM_WORLD, &s);
			dens += buf;
		}
	}
	else {
		MPI_Send(&dens, 1, MPI_DOUBLE, 0, 30, MPI_COMM_WORLD);
	}
	
	//printf("Time: %d\t Density: %f\n",time, dens);
	return dens;
}


//Calcule the velocity for the parallel code
double calc_velocity_par(s_lattice *l, s_comm *comm) {
	//local variables
	int i;
	int procs = comm->proc_y;
	
	double vel[2]; 
	vel[0] = 0;
	vel[1] = 0;
		
	double *buffer= malloc(2 * sizeof(double));
	//check memory allocation
	if (buffer == NULL) {
		printf("Erro calc_velocity_par \n");
		exit(11);
	}
	
	//If sublattice has l->lx_total/2 than calculate velocity dimension for this point
	if (l->lx_total/2 >= l->lx_first  && l->lx_total/2 < l->lx_last) {
		buffer = calc_velocity(l);
		vel[0] = buffer[0];
		vel[1] = buffer[1];
	}

	//Master processor receive velocities
	if (comm->myid == 0) {
		MPI_Status s;
		//If master processor has l->lx_total/2 is not necessary communicate it
		if (l->lx_total/2 >= l->lx_first && l->lx_total/2 < l->lx_last)
			procs -= 1;
		for(i = 0; i < procs; i++) {
			MPI_Recv(buffer, 2, MPI_DOUBLE, MPI_ANY_SOURCE, 25, MPI_COMM_WORLD, &s);
			vel[0] += buffer[0];
			vel[1] += buffer[1];
		}
	}
	//Other processors send his velocities if it has l->lx_total/2
	else {
		if (l->lx_total/2 >= l->lx_first && l->lx_total/2 < l->lx_last) 
			MPI_Send(vel, 2, MPI_DOUBLE, 0, 25, MPI_COMM_WORLD);
	}

	//All processors receive the result
	MPI_Bcast(vel, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	//And finally, the velocity at axis l->lx_total/2
	return vel[0]/vel[1];
}


//////////////////////////////////////////////////////////////////////
// Send
//////////////////////////////////////////////////////////////////////

//Memory allocate for communications data and check if the allocation is correct
double *alloc_border(int buffer_size) {
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
	if (buffer == NULL) {
		printf("error: alloc_border\n");
		exit(11);
	}
	return buffer;
}


//Sends a data structure through MPI_Isend
void send_border(double *buffer, int buffer_size, int *coord, int label) {
	//local variables
	int neighbour;
	MPI_Request request;
	MPI_Status stat;
	MPI_Cart_rank(comm->grid_comm, coord, &neighbour);
	MPI_Isend(buffer, buffer_size, MPI_DOUBLE, neighbour, label, comm->grid_comm, &request);
	//printf("Eu pos_x %d Enviando... %d\n", comm->myid, neighbour);
	//fflush (stdout);
	MPI_Wait(&request, &stat);
	free(buffer);
}


//Sends a point through MPI_Isend
void send_border(double *buffer, int buffer_size, int *coord, int label) {
	//local variables
	int neighbour;
	MPI_Request request;
	MPI_Status stat;
	MPI_Cart_rank(comm->grid_comm, coord, &neighbour);
	MPI_Isend(buffer, buffer_size, MPI_DOUBLE, neighbour, label, comm->grid_comm, &request);
	//printf("Eu pos_x %d Enviando... %d\n", comm->myid, neighbour);
	//fflush (stdout);
	MPI_Wait(&request, &stat);
}


//Communication function to send positive x neighbours
void send_pos_x(s_lattice *l, s_comm *comm) {
	//local variables
	long int buffer_size = COL * (l->ly - 2);
	double *buffer = alloc_border(buffer_size);
	
	int y;
	long int k = 0;

	for(y = 1; y < l->ly - 1; y++) {
		buffer[k++] = l->temp[l->lx - 1][y][1];
		buffer[k++] = l->temp[l->lx - 1][y][5];
		buffer[k++] = l->temp[l->lx - 1][y][8];
	}
	
	int neighbour;
	int coord[2];
	coord[0] = (comm->myid_x + 1)%comm->proc_x;
	coord[1] = comm->myid_y;

	send_border(buffer, buffer_size, coord, 1);
}


//Communication function to send negative x neighbours
void send_neg_x(s_lattice *l, s_comm *comm) {
	//local variables
	long int buffer_size = COL * (l->ly - 2);
	double *buffer = alloc_border(buffer_size);
	
	int y;
	long int k = 0;

	for(y = 1; y < l->ly - 1; y++) {
		buffer[k++] = l->temp[0][y][3];
		buffer[k++] = l->temp[0][y][6];
		buffer[k++] = l->temp[0][y][7];
	}
	
	int neighbour;
	int coord[2];
	coord[0] = (comm->myid_x - 1 + comm->proc_x)%comm->proc_x;
	coord[1] = comm->myid_y;

	send_border(buffer, buffer_size, coord, 2);
}


//Communication function to send positive y neighbours
void send_pos_y(s_lattice *l, s_comm *comm) {
	//local variables
	long int buffer_size = COL * (l->lx - 2);
	double *buffer = alloc_border(buffer_size);
	
	int x;
	long int k = 0;

	for(x = 1; x < l->lx - 1; x++) {
		buffer[k++] = l->temp[x][l->ly - 1][2];
		buffer[k++] = l->temp[x][l->ly - 1][5];
		buffer[k++] = l->temp[x][l->ly - 1][6];
	}
	
	int neighbour;
	int coord[2];
	coord[1] = (comm->myid_y + 1)%comm->proc_y;
	coord[0] = comm->myid_x;

	send_border(buffer, buffer_size, coord, 3);
}


//Communication function to send negative y neighbours
void send_neg_y(s_lattice *l, s_comm *comm) {
	//local variables
	long int buffer_size = COL * (l->lx - 2);
	double *buffer = alloc_border(buffer_size);
	
	int x;
	long int k = 0;

	for(x = 1; x < l->lx - 1; x++) {
		buffer[k++] = l->temp[x][0][4];
		buffer[k++] = l->temp[x][0][7];
		buffer[k++] = l->temp[x][0][8];
	}
	
	int neighbour;
	int coord[2];
	coord[1] = (comm->myid_y - 1 + comm->proc_y)%comm->proc_y;
	coord[0] = comm->myid_x;

	send_border(buffer, buffer_size, coord, 4);
}


void send_corners(s_lattice *l, s_comm *comm) {
	buffer_size = 1;
	int coord[2];

	coord[0] = (comm->myid_x - 1 + comm->proc_x)%comm->proc_x;
	coord[1] = (comm->myid_y - 1 + comm->proc_y)%comm->proc_y;
	send_border(&l->temp[0][0][7], buffer_size, coord, 5);
	
	coord[0] = (comm->myid_x - 1 + comm->proc_x)%comm->proc_x;
	coord[1] = (comm->myid_y + 1)%comm->proc_y;
	send_border(&l->temp[0][l->ly - 1][6], buffer_size, coord, 6);

	coord[0] = (comm->myid_x + 1)%comm->proc_x;
	coord[1] = (comm->myid_y - 1 + comm->proc_y)%comm->proc_y;
	send_border(&l->temp[l->lx - 1][0][8], buffer_size, coord, 7);

	coord[0] = (comm->myid_x + 1)%comm->proc_x;
	coord[1] = (comm->myid_y + 1)%comm->proc_y;
	send_border(&l->temp[l->lx - 1][l->ly - 1][5], buffer_size, coord, 8);

}
