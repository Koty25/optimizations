#include "lb.h"


//////////////////////////////////////////////////////////////////////
// comm initialization and sublattice definition
//////////////////////////////////////////////////////////////////////


//Initialize MPI communications structure
s_comm *initialize_comm(string dim1, string dim2, string dim3) {
	
	//Initialize comm structure
	s_comm *comm = (s_comm *) malloc(sizeof(s_comm));
	
	//Check for memory utilization
	if (comm == NULL) {
		printf("Error in memory allocation in function initialize_comm\n");
		exit(-10);
	}
	
	//Local variables
	int dimensions[3];
	int wrap_around[3];
	int coordinates[3];
	int reorder = false;
	
	//Get dimensions
	dimensions[0] = atoi(dim1);
	dimensions[1] = atoi(dim2);
	dimensions[2] = atoi(dim3);

	//Define wrap_around type
	wrap_around[0] = true;
	wrap_around[1] = true;
	wrap_around[2] = true;

	//Create cartesian communication structure
	MPI_Cart_create(MPI_COMM_WORLD, 3, dimensions, wrap_around, reorder, &comm->grid_comm);
	
	//Define rank an coordinates for each processor
	MPI_Comm_rank(comm->grid_comm, &comm->myid);
	MPI_Cart_coords(comm->grid_comm, comm->myid, 3, coordinates);
	
	//Save number of processors for each dimension in comm structure
	comm->proc_x = dimensions[0];
	comm->proc_y = dimensions[1];
	comm->proc_z = dimensions[2];

	//Save rank identification for each dimension in comm structure
	comm->myid_x = coordinates[0];
	comm->myid_y = coordinates[1];
	comm->myid_z = coordinates[2];
	
	//Return the structure
	return comm;
}


//Define subdimensions of the lattice
void define_sublattice(s_comm *comm, s_lattice *l) {
	
	//Initialize variables
	int size_x = l->lx_total/comm->proc_x;
	int size_y = l->ly_total/comm->proc_y;
	int size_z = l->lz_total/comm->proc_z;
	int remain_x = l->lx_total%comm->proc_x;
	int remain_y = l->ly_total%comm->proc_y;
	int remain_z = l->lz_total%comm->proc_z;
	
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
	
	//z-axis
	if(comm->myid_z > comm->proc_z - 1 - remain_z) {
		size_z += 1;
		l->lz_last = l->lz_total - (comm->proc_z - 1 - comm->myid_z)*size_z;
		l->lz_first = l->lz_last - size_z;
	}
	else {
		l->lz_first = size_z * comm->myid_z;
		l->lz_last = size_z * (comm->myid_z + 1);
	}
	
	//Define the real large size of each sublattice - Overlapping borders
	l->lx = l->lx_last - l->lx_first + 2;
	l->ly = l->ly_last - l->ly_first + 2; 
	l->lz = l->lz_last - l->lz_first + 2;
}


//////////////////////////////////////////////////////////////////////
// Parallel Velocity and Density
//////////////////////////////////////////////////////////////////////

//Check the density for the parallel implementation
double check_density_par(s_lattice *l, s_comm *comm, int time) {
	
	//local variables
	double dens = check_density(l);
	double buf;
	int i;
	int procs = comm->proc_x * comm->proc_y * comm->proc_z - 1;
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
double calc_velocity_par(s_lattice *l, s_comm *comm, int time) {
	
	//local variables
	int i;
	int procs = comm->proc_y * comm->proc_z;
	
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
		buffer = calc_velocity(l, time);
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
	
	//Memory release
	free(buffer);
	
	//And finally, the velocity at axis l->lx_total/2
	return vel[0]/vel[1];
}


//////////////////////////////////////////////////////////////////////
// Send
//////////////////////////////////////////////////////////////////////


//Communication function to send positive x neighbours
void send_pos_x(s_aux *aux, s_lattice *l, s_comm *comm) {
	
	//local variables
	long int buffer_size = COL * (l->ly - 1) * (l->lz - 1);
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
	if (buffer == NULL) {
		printf("Erro send_pos_x\n");
		exit(11);
	}

	MPI_Request request;
	MPI_Status stat;

	int neighbour;
	int coord[3];
	int y, z, j;
	long int k = 0;
	
	for(y = 1; y < l->ly - 1; y++) {
		for(z = 1; z < l->lz - 1; z++) {
			for(j = 0; j < COL; j++) {
				buffer[k] = l->temp[l->lx - 1][y][z][aux->pos_x[j]];
				k++;
			} 
		}
	}
	
	coord[0] = (comm->myid_x + 1)%comm->proc_x;
	coord[1] = comm->myid_y;
	coord[2] = comm->myid_z;
	
	MPI_Cart_rank(comm->grid_comm, coord, &neighbour);
	MPI_Isend(buffer, buffer_size, MPI_DOUBLE, neighbour, 1, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
	free(buffer);
}


//Communication function to send negative x neighbours
void send_neg_x(s_aux *aux, s_lattice *l, s_comm *comm) {
	
	//local variables
	long int buffer_size = COL * (l->ly - 1) * (l->lz - 1);
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
		if (buffer == NULL) {
		printf("Erro send_neg_x\n");
		exit(11);
	}

	MPI_Request request;
	MPI_Status stat;

	int neighbour;
	int coord[3];
	int y, z, j;
	long int k = 0;
	
	for(y = 1; y < l->ly - 1; y++) {
		for(z = 1; z < l->lz - 1; z++) {
			for(j = 0; j < COL; j++) {
				buffer[k] = l->temp[0][y][z][aux->neg_x[j]];
				k++;
			}
		}
	}

	coord[0] = (comm->myid_x - 1 + comm->proc_x)%comm->proc_x;
	coord[1] = comm->myid_y;
	coord[2] = comm->myid_z;
	
	MPI_Cart_rank(comm->grid_comm, coord, &neighbour);
	MPI_Isend(buffer, buffer_size, MPI_DOUBLE, neighbour, 2, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
	free(buffer);
}


//Communication function to send positive y neighbours
void send_pos_y(s_aux *aux, s_lattice *l, s_comm *comm) {
	
	//local variables
	long int buffer_size = COL * (l->lx - 1) * (l->lz - 1);
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
	if (buffer == NULL) {
		printf("Erro send_pos_y\n");
		exit(11);
	}

	MPI_Request request;
	MPI_Status stat;

	int neighbour;
	int coord[3];
	int x, z, j;
	long int k = 0;
	
	for(x = 1; x < l->lx - 1; x++) {
		for(z = 1; z < l->lz - 1; z++) {
			for(j = 0; j < COL; j++) {
				buffer[k] = l->temp[x][l->ly - 1][z][aux->pos_y[j]];
				k++;
			}
		}
	}
	
	coord[0] = comm->myid_x;
	coord[1] = (comm->myid_y + 1)%comm->proc_y;
	coord[2] = comm->myid_z;
	
	MPI_Cart_rank(comm->grid_comm, coord, &neighbour);
	MPI_Isend(buffer, buffer_size, MPI_DOUBLE, neighbour, 3, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
	free(buffer);
}


//Communication function to send negative y neighbours
void send_neg_y(s_aux *aux, s_lattice *l, s_comm *comm) {
	
	//local variables
	long int buffer_size = COL * (l->lx - 1) * (l->lz - 1);
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
	if (buffer == NULL) {
		printf("Erro send_neg_y\n");
		exit(11);
	}

	MPI_Request request;
	MPI_Status stat;

	int neighbour;
	int coord[3];
	int x, z, j;
	long int k = 0;
	
	for(x = 1; x < l->lx - 1; x++) {
		for(z = 1; z < l->lz - 1; z++) {
			for(j = 0; j < COL; j++) {
				buffer[k] = l->temp[x][0][z][aux->neg_y[j]];
				k++;
			}
		}
	}

	coord[0] = comm->myid_x;
	coord[1] = (comm->myid_y - 1 + comm->proc_y)%comm->proc_y;
	coord[2] = comm->myid_z;
	
	MPI_Cart_rank(comm->grid_comm, coord, &neighbour);
	MPI_Isend(buffer, buffer_size, MPI_DOUBLE, neighbour, 4, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
	free(buffer);
}


//Communication function to send positive z neighbours
void send_pos_z(s_aux *aux, s_lattice *l, s_comm *comm) {
	
	//local variables		
	long int buffer_size = COL * (l->lx - 1) * (l->ly - 1);
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
	if (buffer == NULL) {
		printf("Erro send_pos_z\n");
		exit(11);
	}

	MPI_Request request;
	MPI_Status stat;

	int neighbour;
	int coord[3];
	int x, y, j;
	long int k = 0;
	
	for(x = 1; x < l->lx - 1; x++) {
		for(y = 1; y < l->ly - 1; y++) {
			for(j = 0; j < COL; j++) {
				buffer[k] = l->temp[x][y][l->lz - 1][aux->pos_z[j]];
				k++;
			}
		}
	}
	
	coord[0] = comm->myid_x;
	coord[1] = comm->myid_y;
	coord[2] = (comm->myid_z + 1)%comm->proc_z;

	MPI_Cart_rank(comm->grid_comm, coord, &neighbour);
	MPI_Isend(buffer, buffer_size, MPI_DOUBLE, neighbour, 5, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
	free(buffer);
}


//Communication function to send negative z neighbours
void send_neg_z(s_aux *aux, s_lattice *l, s_comm *comm) {
	
	//local variables
	long int buffer_size = COL * (l->lx - 1) * (l->ly - 1);
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
	if (buffer == NULL) {
		printf("Erro send_neg_z\n");
		exit(11);
	}

	MPI_Request request;
	MPI_Status stat;
	
	int neighbour;
	int coord[3];
	int x, y, j;
	long int k = 0;
	
	for(x = 1; x < l->lx - 1; x++) {
		for(y = 1; y < l->ly - 1; y++) {
			for(j = 0; j < COL; j++) {
				buffer[k] = l->temp[x][y][0][aux->neg_z[j]];
				k++;
			}
		}
	}

	coord[0] = comm->myid_x;
	coord[1] = comm->myid_y;
	coord[2] = (comm->myid_z - 1 + comm->proc_z)%comm->proc_z;

	MPI_Cart_rank(comm->grid_comm, coord, &neighbour);
	MPI_Isend(buffer, buffer_size, MPI_DOUBLE, neighbour, 6, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
	free(buffer);
}


//////////////////////////////////////////////////////////////////////
// Receive
//////////////////////////////////////////////////////////////////////


//Communication function to recv negative z neighbours
void recv_pos_x(s_aux *aux, s_lattice *l, s_comm *comm) {
	
	//local variables
	long int buffer_size = COL * (l->ly - 1) * (l->lz - 1);
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
	if (buffer == NULL) {
		printf("Erro recv \n");
		exit(11);
	}
	
	long int k = 0;
	int y, z, j;
	
	MPI_Status status;
	
	MPI_Recv(buffer, buffer_size, MPI_DOUBLE, MPI_ANY_SOURCE, 1, comm->grid_comm, &status);

	for(y = 1; y < l->ly - 1; y++) {
		for(z = 1; z < l->lz - 1; z++) {
			for(j = 0; j < COL; j++) {
				l->temp[1][y][z][aux->pos_x[j]] = buffer[k];
				k++;
			}
		}
	}
	free(buffer);
}


//Communication function to recv negative z neighbours
void recv_neg_x(s_aux *aux, s_lattice *l, s_comm *comm) {
	
	//local variables
	long int buffer_size = COL * (l->ly - 1) * (l->lz - 1);
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
	if (buffer == NULL) {
		printf("Erro recv \n");
		exit(11);
	}
	
	long int k = 0;
	int y, z, j;

	MPI_Status status;

	MPI_Recv(buffer, buffer_size, MPI_DOUBLE, MPI_ANY_SOURCE, 2, comm->grid_comm, &status);

	for(y = 1; y < l->ly - 1; y++) {
		for(z = 1; z < l->lz - 1; z++) {
			for(j = 0; j < COL; j++) {
				l->temp[l->lx - 2][y][z][aux->neg_x[j]] = buffer[k];
				k++;
			}
		}
	}
	free(buffer);
}


//Communication function to recv negative z neighbours
void recv_pos_y(s_aux *aux, s_lattice *l, s_comm *comm) {
	
	//local variables
	long int buffer_size = COL * (l->lx - 1) * (l->lz - 1);
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
	if (buffer == NULL) {
		printf("Erro recv \n");
		exit(11);
	}

	long int k = 0;
	int x, z, j;
	
	MPI_Status status;
	
	MPI_Recv(buffer, buffer_size, MPI_DOUBLE, MPI_ANY_SOURCE, 3, comm->grid_comm, &status);
	
	for(x = 1; x < l->lx - 1; x++) {
		for(z = 1; z < l->lz - 1; z++) {
			for(j = 0; j < COL; j++) {
				l->temp[x][1][z][aux->pos_y[j]] = buffer[k];
				k++;
			}
		}
	}
	free(buffer);
}


//Communication function to recv negative z neighbours
void recv_neg_y(s_aux *aux, s_lattice *l, s_comm *comm) {
	
	//local variables
	long int buffer_size = COL * (l->lx - 1) * (l->lz - 1);
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
	if (buffer == NULL) {
		printf("Erro recv \n");
		exit(11);
	}

	long int k = 0;
	int x, z, j;

	MPI_Status status;
	
	MPI_Recv(buffer, buffer_size, MPI_DOUBLE, MPI_ANY_SOURCE, 4, comm->grid_comm, &status);

	for(x = 1; x < l->lx - 1; x++) {
		for(z = 1; z < l->lz - 1; z++) {
			for(j = 0; j < COL; j++) {
				l->temp[x][l->ly - 2][z][aux->neg_y[j]] = buffer[k];
				k++;
			}
		}
	}
	free(buffer);
}



//Communication function to recv negative z neighbours
void recv_pos_z(s_aux *aux, s_lattice *l, s_comm *comm) {
	
	//local variables
	long int buffer_size = COL * (l->lx - 1) * (l->ly - 1);
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
	if (buffer == NULL) {
		printf("Erro recv \n");
		exit(11);
	}

	long int k = 0;
	int x, y, j;
	
	MPI_Status status;

	MPI_Recv(buffer, buffer_size, MPI_DOUBLE, MPI_ANY_SOURCE, 5, comm->grid_comm, &status);
	
	for(x = 1; x < l->lx - 1; x++) {
		for(y = 1; y < l->ly - 1; y++) {
			for(j = 0; j < COL; j++) {
				l->temp[x][y][1][aux->pos_z[j]] = buffer[k];
				k++;
			}
		}
	}
	free(buffer);
}


//Communication function to recv negative z neighbours
void recv_neg_z(s_aux *aux, s_lattice *l, s_comm *comm) {
	
	//local variables
	long int buffer_size = COL * (l->lx - 1) * (l->ly - 1);
	double *buffer = (double *) malloc(buffer_size * sizeof(double));
	if (buffer == NULL) {
		printf("Erro recv \n");
		exit(11);
	}

	long int k = 0;
	int x, y, j;

	MPI_Status status;
	
	MPI_Recv(buffer, buffer_size, MPI_DOUBLE, MPI_ANY_SOURCE, 6, comm->grid_comm, &status);

	for(x = 1; x < l->lx - 1; x++) {
		for(y = 1; y < l->ly - 1; y++) {
			for(j = 0; j < COL; j++) {
				l->temp[x][y][l->lz - 2][aux->neg_z[j]] = buffer[k];
				k++;
			}
		}
	}
	free(buffer);
}


///////////////////////////////////////////////////////////////////////////////
// Copies
///////////////////////////////////////////////////////////////////////////////

//Instead communicate, this function make a local data copy of x
void copy_x(s_lattice *l, s_aux *aux) {
	
	//local variables
	int i, j, k;

	for(i = 1; i < l->ly - 1; i++) {
		for(j = 1; j < l->lz - 1; j++) {
			for(k = 0; k < COL; k++) {
				l->temp[1][i][j][aux->pos_x[k]] = l->temp[l->lx - 1][i][j][aux->pos_x[k]];
				l->temp[l->lx - 2][i][j][aux->neg_x[k]] = l->temp[0][i][j][aux->neg_x[k]];
			}
		}
	}
}


//Instead communicate, this function make a local data copy of y
void copy_y(s_lattice *l, s_aux *aux) {
	
	//local variables
	int i, j, k;
	
	for(i = 1; i < l->lx - 1; i++) {
		for(j = 1; j < l->lz - 1; j++) {
			for(k = 0; k < COL; k++) {
				l->temp[i][1][j][aux->pos_y[k]] = l->temp[i][l->ly - 1][j][aux->pos_y[k]];
				l->temp[i][l->ly - 2][j][aux->neg_y[k]] = l->temp[i][0][j][aux->neg_y[k]];
			}
		}
	}
}


//Instead communicate, this function make a local data copy of z
void copy_z(s_lattice *l, s_aux *aux) {
	
	//local variables
	int i, j, k;
	
	for(i = 1; i < l->lx - 1; i++) {
		for(j = 1; j < l->ly - 1; j++) {
			for(k = 0; k < COL; k++) {
				l->temp[i][j][1][aux->pos_z[k]] = l->temp[i][j][l->lz - 1][aux->pos_z[k]];
				l->temp[i][j][l->lz - 2][aux->neg_z[k]] = l->temp[i][j][0][aux->neg_z[k]];
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////
// Corners 
// Receive diagonal points when its are at the border
//////////////////////////////////////////////////////////////////////

void corner_recv_y_neg_z_neg(double *buffer_x, s_lattice *l, s_comm *comm) {
	MPI_Status status;
	long int k;
	MPI_Recv(buffer_x, l->lx-2, MPI_DOUBLE, MPI_ANY_SOURCE, 19, comm->grid_comm, &status);
	for(k = 1; k < l->lx - 1; k++) 
		l->temp[k][l->ly-2][l->lz-2][17] = buffer_x[k-1];
}
	

void corner_recv_y_pos_z_neg(double *buffer_x, s_lattice *l, s_comm *comm) {
	MPI_Status status;
	long int k;
	MPI_Recv(buffer_x, l->lx-2, MPI_DOUBLE, MPI_ANY_SOURCE, 20, comm->grid_comm, &status);
	for(k = 1; k < l->lx - 1; k++) 
		l->temp[k][1][l->lz-2][18] = buffer_x[k-1];
}


void corner_recv_y_neg_z_pos(double *buffer_x, s_lattice *l, s_comm *comm) {
	MPI_Status status;
	long int k;
	MPI_Recv(buffer_x, l->lx-2, MPI_DOUBLE, MPI_ANY_SOURCE, 21, comm->grid_comm, &status);
	for(k = 1; k < l->lx - 1; k++) 
		l->temp[k][l->ly-2][1][16] = buffer_x[k-1];
}


void corner_recv_y_pos_z_pos(double *buffer_x, s_lattice *l, s_comm *comm) {
	MPI_Status status;
	long int k;
	MPI_Recv(buffer_x, l->lx-2, MPI_DOUBLE, MPI_ANY_SOURCE, 22, comm->grid_comm, &status);
	for(k = 1; k < l->lx - 1; k++) 
		l->temp[k][1][1][15] = buffer_x[k-1];
}


///////////////////////////////////////////////

void corner_recv_x_neg_z_neg(double *buffer_y, s_lattice *l, s_comm *comm) {
	MPI_Status status;
	long int k;
	MPI_Recv(buffer_y, l->ly-2, MPI_DOUBLE, MPI_ANY_SOURCE, 15, comm->grid_comm, &status);
	for(k = 1; k < l->ly - 1; k++) 
		l->temp[l->lx-2][k][l->lz-2][12] = buffer_y[k-1];
}
	

void corner_recv_x_pos_z_neg(double *buffer_y, s_lattice *l, s_comm *comm) {
	MPI_Status status;
	long int k;
	MPI_Recv(buffer_y, l->ly-2, MPI_DOUBLE, MPI_ANY_SOURCE, 16, comm->grid_comm, &status);
	for(k = 1; k < l->ly - 1; k++) 
		l->temp[1][k][l->lz-2][14] = buffer_y[k-1];
}


void corner_recv_x_neg_z_pos(double *buffer_y, s_lattice *l, s_comm *comm) {
	MPI_Status status;
	long int k;
	MPI_Recv(buffer_y, l->ly-2, MPI_DOUBLE, MPI_ANY_SOURCE, 17, comm->grid_comm, &status);
	for(k = 1; k < l->ly - 1; k++) 
		l->temp[l->lx-2][k][1][11] = buffer_y[k-1];
}


void corner_recv_x_pos_z_pos(double *buffer_y, s_lattice *l, s_comm *comm) {
	MPI_Status status;
	long int k;
	MPI_Recv(buffer_y, l->ly-2, MPI_DOUBLE, MPI_ANY_SOURCE, 18, comm->grid_comm, &status);
	for(k = 1; k < l->ly - 1; k++) 
		l->temp[1][k][1][9] = buffer_y[k-1];
}


///////////////////////////////////////////////

void corner_recv_x_neg_y_neg(double *buffer_z, s_lattice *l, s_comm *comm) {
	MPI_Status status;
	long int k;
	MPI_Recv(buffer_z, l->lz-2, MPI_DOUBLE, MPI_ANY_SOURCE, 11, comm->grid_comm, &status);
	for(k = 1; k < l->lz - 1; k++) 
		l->temp[l->lx-2][l->ly-2][k][6] = buffer_z[k-1];
}


void corner_recv_x_pos_y_neg(double *buffer_z, s_lattice *l, s_comm *comm) {
	MPI_Status status;
	long int k;
	MPI_Recv(buffer_z, l->lz-2, MPI_DOUBLE, MPI_ANY_SOURCE, 12, comm->grid_comm, &status);
	for(k = 1; k < l->lz - 1; k++) 
		l->temp[1][l->ly-2][k][8] = buffer_z[k-1];
}
	

void corner_recv_x_neg_y_pos(double *buffer_z, s_lattice *l, s_comm *comm) {
	MPI_Status status;
	long int k;
	MPI_Recv(buffer_z, l->lz-2, MPI_DOUBLE, MPI_ANY_SOURCE, 13, comm->grid_comm, &status);
	for(k = 1; k < l->lz - 1; k++)
		l->temp[l->lx-2][1][k][4] = buffer_z[k-1];
}
	

void corner_recv_x_pos_y_pos(double *buffer_z, s_lattice *l, s_comm *comm) {
	MPI_Status status;
	long int k;
	MPI_Recv(buffer_z, l->lz-2, MPI_DOUBLE, MPI_ANY_SOURCE, 14, comm->grid_comm, &status);
	for(k = 1; k < l->lz - 1; k++)
		l->temp[1][1][k][2] = buffer_z[k-1];
}
	

//////////////////////////////////////////////////////////////////////
// Corners
//////////////////////////////////////////////////////////////////////

//Send corners from x face
void parallel_corner_x(s_lattice *l, s_aux *aux, s_comm *comm) {
	
	//local variables
	double *buffer_x = (double *) malloc((l->lx - 2)*sizeof(double));
	if (buffer_x == NULL) {
		printf("Erro recv \n");
		exit(11);
	}

	corner_y_neg_z_neg(buffer_x, l, comm);
	corner_y_pos_z_neg(buffer_x, l, comm);
	corner_y_neg_z_pos(buffer_x, l, comm);
	corner_y_pos_z_pos(buffer_x, l, comm);
	
	free(buffer_x);
}


//k 0 0 17
void corner_y_neg_z_neg(double *buffer_x, s_lattice *l, s_comm *comm) {
	
	//local variables
	int coord[3];
	coord[0] = comm->myid_x;
	coord[1] = (comm->myid_y - 1 + comm->proc_y)%comm->proc_y;
	coord[2] = (comm->myid_z - 1 + comm->proc_z)%comm->proc_z;
	
	int y_neg_z_neg;

	MPI_Cart_rank(comm->grid_comm, coord, &y_neg_z_neg);
	if(comm->myid == y_neg_z_neg) {
		long int k;
		for(k = 1; k < l->lx - 1; k++) {
			l->temp[k][l->ly - 2][l->lz - 2][17] = l->temp[k][0][0][17];
		}
	} else {
		if(comm->myid_y%2 == 0) {
			corner_send_y_neg_z_neg(y_neg_z_neg, buffer_x, l, comm);
			corner_recv_y_neg_z_neg(buffer_x, l, comm);
		} else {			
			corner_recv_y_neg_z_neg(buffer_x, l, comm);
			corner_send_y_neg_z_neg(y_neg_z_neg, buffer_x, l, comm);
		}
	}
}


//k y 0 18
void corner_y_pos_z_neg(double *buffer_x, s_lattice *l, s_comm *comm) {
	
	//local variables
	int coord[3];
	coord[0] = comm->myid_x;
	coord[1] = (comm->myid_y + 1)%comm->proc_y;
	coord[2] = (comm->myid_z - 1 + comm->proc_z)%comm->proc_z;

	int y_pos_z_neg;

	MPI_Cart_rank(comm->grid_comm, coord, &y_pos_z_neg);
	if(comm->myid == y_pos_z_neg) {
		long int k;
		for(k = 1; k < l->lx - 1; k++) {
			l->temp[k][1][l->lz - 2][18] = l->temp[k][l->ly - 1][0][18];
		}
	} else {
		if(comm->myid_y%2 == 0) {
			corner_send_y_pos_z_neg(y_pos_z_neg, buffer_x, l, comm);
			corner_recv_y_pos_z_neg(buffer_x, l, comm);
		} else {
			corner_recv_y_pos_z_neg(buffer_x, l, comm);
			corner_send_y_pos_z_neg(y_pos_z_neg, buffer_x, l, comm);
		}
	}
}


//k 0 z 16
void corner_y_neg_z_pos(double *buffer_x, s_lattice *l, s_comm *comm) {
	
	//local variables
	int coord[3];
	coord[0] = comm->myid_x;
	coord[1] = (comm->myid_y - 1 + comm->proc_y)%comm->proc_y;
	coord[2] = (comm->myid_z + 1)%comm->proc_z;

	int y_neg_z_pos;

	MPI_Cart_rank(comm->grid_comm, coord, &y_neg_z_pos);
	if(comm->myid == y_neg_z_pos) {
		long int k;
		for(k = 1; k < l->lx - 1; k++) {
			l->temp[k][l->ly - 2][1][16] = l->temp[k][0][l->lz - 1][16];
		}
	} else {
		if(comm->myid_y%2 == 0) {
			corner_send_y_neg_z_pos(y_neg_z_pos, buffer_x, l, comm);
			corner_recv_y_neg_z_pos(buffer_x, l, comm);
		} else {
			corner_recv_y_neg_z_pos(buffer_x, l, comm);
			corner_send_y_neg_z_pos(y_neg_z_pos, buffer_x, l, comm);
		}
	}
}


//k y z 15
void corner_y_pos_z_pos(double *buffer_x, s_lattice *l, s_comm *comm) {
	
	//local variables
	int coord[3];
	coord[0] = comm->myid_x;
	coord[1] = (comm->myid_y + 1)%comm->proc_y;
	coord[2] = (comm->myid_z + 1)%comm->proc_z;

	int y_pos_z_pos;

	MPI_Cart_rank(comm->grid_comm, coord, &y_pos_z_pos);
	if(comm->myid == y_pos_z_pos) {
		long int k;
		for(k = 1; k < l->lx - 1; k++) {
			l->temp[k][1][1][15] = l->temp[k][l->ly - 1][l->lz - 1][15];
		}
	} else {
		if(comm->myid_y%2 == 0) {
			corner_send_y_pos_z_pos(y_pos_z_pos, buffer_x, l, comm);
			corner_recv_y_pos_z_pos(buffer_x, l, comm);
		} else {
			corner_recv_y_pos_z_pos(buffer_x, l, comm);
			corner_send_y_pos_z_pos(y_pos_z_pos, buffer_x, l, comm);
		}
	}
}


//
void corner_send_y_neg_z_neg(int y_neg_z_neg, double *buffer_x, s_lattice *l, s_comm *comm) {
	
	//local variables
	MPI_Request request;
	MPI_Status stat;
	long int k;
	for(k = 1; k < l->lx - 1; k++) {
		buffer_x[k-1] = l->temp[k][0][0][17];
	}
	MPI_Isend(buffer_x, l->lx-2, MPI_DOUBLE, y_neg_z_neg, 19, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
}


//
void corner_send_y_pos_z_neg(int y_pos_z_neg, double *buffer_x, s_lattice *l, s_comm *comm) {
	
	//local variables
	MPI_Request request;
	MPI_Status stat;
	long int k;
	for(k = 1; k < l->lx - 1; k++) {
		buffer_x[k-1] = l->temp[k][l->ly-1][0][18];
	}
	MPI_Isend(buffer_x, l->lx-2, MPI_DOUBLE, y_pos_z_neg, 20, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
}


//
void corner_send_y_neg_z_pos(int y_neg_z_pos, double *buffer_x, s_lattice *l, s_comm *comm) {
	
	//local variables
	MPI_Request request;
	MPI_Status stat;
	long int k;
	for(k = 1; k < l->lx - 1; k++) {
		buffer_x[k-1] = l->temp[k][0][l->lz-1][16];
	}
	MPI_Isend(buffer_x, l->lx-2, MPI_DOUBLE, y_neg_z_pos, 21, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
}


//
void corner_send_y_pos_z_pos(int y_pos_z_pos, double *buffer_x, s_lattice *l, s_comm *comm) {
	
	//local variables
	MPI_Request request;
	MPI_Status stat;
	long int k;
	for(k = 1; k < l->lx - 1; k++) {
		buffer_x[k-1] = l->temp[k][l->ly-1][l->lz-1][15];
	}
	MPI_Isend(buffer_x, l->lx-2, MPI_DOUBLE, y_pos_z_pos, 22, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
}


//Send corners from y face
void parallel_corner_y(s_lattice *l, s_aux *aux, s_comm *comm) {
	
	//local variables
	double *buffer_y = (double *) malloc((l->ly - 2)*sizeof(double));
	if (buffer_y == NULL) {
		printf("Erro recv \n");
		exit(11);
	}
	
	corner_x_neg_z_neg(buffer_y, l, comm);
	corner_x_pos_z_neg(buffer_y, l, comm);
	corner_x_neg_z_pos(buffer_y, l, comm);
	corner_x_pos_z_pos(buffer_y, l, comm);

	free(buffer_y);
}


//0 k z 12
void corner_x_neg_z_neg(double *buffer_y, s_lattice *l, s_comm *comm) {	
	
	//local variables
	int coord[3];
	coord[1] = comm->myid_y;
	coord[0] = (comm->myid_x - 1 + comm->proc_x)%comm->proc_x;
	coord[2] = (comm->myid_z - 1 + comm->proc_z)%comm->proc_z;
	
	int x_neg_z_neg;
	
	MPI_Cart_rank(comm->grid_comm, coord, &x_neg_z_neg);
	if (comm->myid == x_neg_z_neg) {
		long int k;
		for(k = 1; k < l->ly - 1; k++) {
			l->temp[l->lx - 2][k][l->lz - 2][12] = l->temp[0][k][0][12];
		}
	} else {
		if(comm->myid_z%2 == 0) {
			corner_send_x_neg_z_neg(x_neg_z_neg, buffer_y, l, comm);
			corner_recv_x_neg_z_neg(buffer_y, l, comm);
		} else {
			corner_recv_x_neg_z_neg(buffer_y, l, comm);
			corner_send_x_neg_z_neg(x_neg_z_neg, buffer_y, l, comm);
		}
	}
}


//x k 0 14
void corner_x_pos_z_neg(double *buffer_y, s_lattice *l, s_comm *comm) {	
	
	//local variables
	int coord[3];
	coord[1] = comm->myid_y;
	coord[0] = (comm->myid_x + 1)%comm->proc_x;
	coord[2] = (comm->myid_z - 1 + comm->proc_z)%comm->proc_z;
	
	int x_pos_z_neg;
	
	MPI_Cart_rank(comm->grid_comm, coord, &x_pos_z_neg);
	if (comm->myid == x_pos_z_neg) {
		long int k;
		for(k = 1; k < l->ly - 1; k++) {
			l->temp[1][k][l->lz - 2][14] = l->temp[l->lx - 1][k][0][14];
		}
	} else {
		if(comm->myid_z%2 == 0) {
			corner_send_x_pos_z_neg(x_pos_z_neg, buffer_y, l, comm);
			corner_recv_x_pos_z_neg(buffer_y, l, comm);
		} else {
			corner_recv_x_pos_z_neg(buffer_y, l, comm);
			corner_send_x_pos_z_neg(x_pos_z_neg, buffer_y, l, comm);
		}
	}
}


//0 k z 11
void corner_x_neg_z_pos(double *buffer_y, s_lattice *l, s_comm *comm) {	
	
	//local variables
	int coord[3];
	coord[1] = comm->myid_y;
	coord[0] = (comm->myid_x - 1 + comm->proc_x)%comm->proc_x;
	coord[2] = (comm->myid_z + 1)%comm->proc_z;
	
	int x_neg_z_pos;
	
	MPI_Cart_rank(comm->grid_comm, coord, &x_neg_z_pos);
	if (comm->myid == x_neg_z_pos) {
		long int k;
		for(k = 1; k < l->ly - 1; k++) {
			l->temp[l->lx - 2][k][1][11] = l->temp[0][k][l->lz - 1][11];
		}
	} else {
		if(comm->myid_z%2 == 0) {
			corner_send_x_neg_z_pos(x_neg_z_pos, buffer_y, l, comm);
			corner_recv_x_neg_z_pos(buffer_y, l, comm);
		} else {
			corner_recv_x_neg_z_pos(buffer_y, l, comm);
			corner_send_x_neg_z_pos(x_neg_z_pos, buffer_y, l, comm);
		}
	}
}


//x k z 9
void corner_x_pos_z_pos(double *buffer_y, s_lattice *l, s_comm *comm) {	
	
	//local variables
	int coord[3];
	coord[1] = comm->myid_y;
	coord[0] = (comm->myid_x + 1)%comm->proc_x;
	coord[2] = (comm->myid_z + 1)%comm->proc_z;
	
	int x_pos_z_pos;
	
	MPI_Cart_rank(comm->grid_comm, coord, &x_pos_z_pos);
	if (comm->myid == x_pos_z_pos) {
		long int k;
		for(k = 1; k < l->ly - 1; k++) {
			l->temp[1][k][1][9] = l->temp[l->lx - 1][k][l->lz - 1][9];
		}
	} else {
		if(comm->myid_z%2 == 0) {
			corner_send_x_pos_z_pos(x_pos_z_pos, buffer_y, l, comm);
			corner_recv_x_pos_z_pos(buffer_y, l, comm);
		} else {
			corner_recv_x_pos_z_pos(buffer_y, l, comm);
			corner_send_x_pos_z_pos(x_pos_z_pos, buffer_y, l, comm);
		}
	}
}


//
void corner_send_x_neg_z_neg(int x_neg_z_neg, double *buffer_y, s_lattice *l, s_comm *comm) {
	
	//local variables
	MPI_Request request;
	MPI_Status stat;
	long int k;
	for(k = 1; k < l->ly - 1; k++) {
		buffer_y[k-1] = l->temp[0][k][0][12];
	}
	MPI_Isend(buffer_y, l->ly-2, MPI_DOUBLE, x_neg_z_neg, 15, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
}


//
void corner_send_x_pos_z_neg(int x_pos_z_neg, double *buffer_y, s_lattice *l, s_comm *comm) {
	
	//local variables
	MPI_Request request;
	MPI_Status stat;
	long int k;
	for(k = 1; k < l->ly - 1; k++) {
		buffer_y[k-1] = l->temp[l->lx-1][k][0][14];
	}
	MPI_Isend(buffer_y, l->ly-2, MPI_DOUBLE, x_pos_z_neg, 16, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
}


//
void corner_send_x_neg_z_pos(int x_neg_z_pos, double *buffer_y, s_lattice *l, s_comm *comm) {
	
	//local variables
	MPI_Request request;
	MPI_Status stat;
	long int k;
	for(k = 1; k < l->ly - 1; k++) {
		buffer_y[k-1] = l->temp[0][k][l->lz-1][11];
	}
	MPI_Isend(buffer_y, l->ly-2, MPI_DOUBLE, x_neg_z_pos, 17, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
}
	

//
void corner_send_x_pos_z_pos(int x_pos_z_pos, double *buffer_y, s_lattice *l, s_comm *comm) {
	
	//local variables
	MPI_Request request;
	MPI_Status stat;
	long int k;
	for(k = 1; k < l->ly - 1; k++) {
		buffer_y[k-1] = l->temp[l->lx-1][k][l->lz-1][9];
	}
	MPI_Isend(buffer_y, l->ly-2, MPI_DOUBLE, x_pos_z_pos, 18, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
}
	


//Send corners from z face
void parallel_corner_z(s_lattice *l, s_aux *aux, s_comm *comm) {
	
	//local variables
	double *buffer_z = (double *) malloc((l->lz - 2) * sizeof(double));
	if (buffer_z == NULL) {
		printf("Erro recv \n");
		exit(11);
	}
	
	corner_x_neg_y_neg(buffer_z, l, comm);
	corner_x_pos_y_neg(buffer_z, l, comm);
	corner_x_neg_y_pos(buffer_z, l, comm);
	corner_x_pos_y_pos(buffer_z, l, comm);
	
	free(buffer_z);
}


//0 0 k 6
void corner_x_neg_y_neg(double *buffer_z, s_lattice *l, s_comm *comm) {
	
	//local variables
	int coord[3];
	coord[2] = comm->myid_z;
	coord[0] = (comm->myid_x - 1 + comm->proc_x)%comm->proc_x;
	coord[1] = (comm->myid_y - 1 + comm->proc_y)%comm->proc_y;

	int x_neg_y_neg;

	MPI_Cart_rank(comm->grid_comm, coord, &x_neg_y_neg);
	if(comm->myid == x_neg_y_neg) {
		long int k;
		for(k = 1; k < l->lz - 1; k++) {
			l->temp[l->lx - 2][l->ly - 2][k][6] = l->temp[0][0][k][6];
		}
	} else {
		if(comm->myid_x%2 == 0) {
			corner_send_x_neg_y_neg(x_neg_y_neg, buffer_z, l, comm);
			corner_recv_x_neg_y_neg(buffer_z, l, comm);
		} else {
			corner_recv_x_neg_y_neg(buffer_z, l, comm);
			corner_send_x_neg_y_neg(x_neg_y_neg, buffer_z, l, comm);
		}
	}
}


//x 0 k 8
void corner_x_pos_y_neg(double *buffer_z, s_lattice *l, s_comm *comm) {
	
	//local variables
	int coord[3];
	coord[2] = comm->myid_z;
	coord[0] = (comm->myid_x + 1)%comm->proc_x;
	coord[1] = (comm->myid_y - 1 + comm->proc_y)%comm->proc_y;

	int x_pos_y_neg;

	MPI_Cart_rank(comm->grid_comm, coord, &x_pos_y_neg);
	if(comm->myid == x_pos_y_neg) {
		long int k;
		for(k = 1; k < l->lz - 1; k++) {
			l->temp[1][l->ly - 2][k][8] = l->temp[l->lx - 1][0][k][8];
		}
	} else {
		if(comm->myid_x%2 == 0) {
			corner_send_x_pos_y_neg(x_pos_y_neg, buffer_z, l, comm);
			corner_recv_x_pos_y_neg(buffer_z, l, comm);
		} else {
			corner_recv_x_pos_y_neg(buffer_z, l, comm);
			corner_send_x_pos_y_neg(x_pos_y_neg, buffer_z, l, comm);
		}
	}
}


//0 y k 4
void corner_x_neg_y_pos(double *buffer_z, s_lattice *l, s_comm *comm) {
	
	//local variables
	int coord[3];
	coord[2] = comm->myid_z;
	coord[0] = (comm->myid_x - 1 + comm->proc_x)%comm->proc_x;
	coord[1] = (comm->myid_y + 1)%comm->proc_y;

	int x_neg_y_pos;

	MPI_Cart_rank(comm->grid_comm, coord, &x_neg_y_pos);
	if(comm->myid == x_neg_y_pos) {
		long int k;
		for(k = 1; k < l->lz - 1; k++) {
			l->temp[l->lx - 2][1][k][4] = l->temp[0][l->ly - 1][k][4];
		}
	} else {
		if(comm->myid_x%2 == 0) {
			corner_send_x_neg_y_pos(x_neg_y_pos, buffer_z, l, comm);
			corner_recv_x_neg_y_pos(buffer_z, l, comm);
		} else {
			corner_recv_x_neg_y_pos(buffer_z, l, comm);
			corner_send_x_neg_y_pos(x_neg_y_pos, buffer_z, l, comm);
		}
	}
}


//x y k 2
void corner_x_pos_y_pos(double *buffer_z, s_lattice *l, s_comm *comm) {
	
	//local variables
	int coord[3];
	coord[2] = comm->myid_z;
	coord[0] = (comm->myid_x + 1)%comm->proc_x;
	coord[1] = (comm->myid_y + 1)%comm->proc_y;

	int x_pos_y_pos;

	MPI_Cart_rank(comm->grid_comm, coord, &x_pos_y_pos);
	if(comm->myid == x_pos_y_pos) {
		long int k;
		for(k = 1; k < l->lz - 1; k++) {
			l->temp[1][1][k][2] = l->temp[l->lx - 1][l->ly - 1][k][2];
		}
	} else {
		if(comm->myid_x%2 == 0) {
			corner_send_x_pos_y_pos(x_pos_y_pos, buffer_z, l, comm);
			corner_recv_x_pos_y_pos(buffer_z, l, comm);
		} else {
			corner_recv_x_pos_y_pos(buffer_z, l, comm);
			corner_send_x_pos_y_pos(x_pos_y_pos, buffer_z, l, comm);
		}
	}
}


//
void corner_send_x_neg_y_neg(int x_neg_y_neg, double *buffer_z, s_lattice *l, s_comm *comm) {
	
	//local variables
	MPI_Request request;
	MPI_Status stat;
	long int k;
	for(k = 1; k < l->lz - 1; k++) {
		buffer_z[k-1] = l->temp[0][0][k][6];
	}
	MPI_Isend(buffer_z, l->lz-2, MPI_DOUBLE, x_neg_y_neg, 11, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
}


//
void corner_send_x_pos_y_neg(int x_pos_y_neg, double *buffer_z, s_lattice *l, s_comm *comm) {
	
	//local variables
	MPI_Request request;
	MPI_Status stat;
	long int k;

	for(k = 1; k < l->lz - 1; k++) {
		buffer_z[k-1] = l->temp[l->lx-1][0][k][8];
	}
	MPI_Isend(buffer_z, l->lz-2, MPI_DOUBLE, x_pos_y_neg, 12, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
}


//
void corner_send_x_neg_y_pos(int x_neg_y_pos, double *buffer_z, s_lattice *l, s_comm *comm) {
	
	//local variables
	MPI_Request request;
	MPI_Status stat;
	long int k;

	for(k = 1; k < l->lz - 1; k++) {
		buffer_z[k-1] = l->temp[0][l->ly-1][k][4];
	}
	MPI_Isend(buffer_z, l->lz-2, MPI_DOUBLE, x_neg_y_pos, 13, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
}


//
void corner_send_x_pos_y_pos(int x_pos_y_pos, double *buffer_z, s_lattice *l, s_comm *comm) {
	
	//local variables
	MPI_Request request;
	MPI_Status stat;
	long int k;

	for(k = 1; k < l->lz - 1; k++) {
		buffer_z[k-1] = l->temp[l->lx-1][l->ly-1][k][2];
	}
	MPI_Isend(buffer_z, l->lz-2, MPI_DOUBLE, x_pos_y_pos, 14, comm->grid_comm, &request);
	MPI_Wait(&request, &stat);
}


//////////////////////////////////////////////////////////////////////
// Communication/Sincronization
//////////////////////////////////////////////////////////////////////


//makes the sends and receive operations 
void sincronization(s_aux *aux, s_lattice *l, s_comm *comm) {

	//Send boundary points
	if (comm->proc_x == 1) 
		copy_x(l, aux);
	else {
		if(comm->myid_x%2 == 0){ 
			send_neg_x(aux, l, comm);
			recv_neg_x(aux, l, comm);
			send_pos_x(aux, l, comm);
			recv_pos_x(aux, l, comm);
	} else {
			recv_neg_x(aux, l, comm);
			send_neg_x(aux, l, comm);
			recv_pos_x(aux, l, comm);
			send_pos_x(aux, l, comm);
		}
	}

	if (comm->proc_y == 1) 
		copy_y(l, aux);
	else {
		if(comm->myid_y%2 == 0) {
			send_neg_y(aux, l, comm);
			recv_neg_y(aux, l, comm);
			send_pos_y(aux, l, comm);
			recv_pos_y(aux, l, comm);
		} else {
			recv_neg_y(aux, l, comm);
			send_neg_y(aux, l, comm);
			recv_pos_y(aux, l, comm);
			send_pos_y(aux, l, comm);
		}
	}

	if (comm->proc_z == 1)
		copy_z(l, aux);
	else {
		if(comm->myid_z%2 == 0) {
			send_neg_z(aux, l, comm);
			recv_neg_z(aux, l, comm);
			send_pos_z(aux, l, comm);
			recv_pos_z(aux, l, comm);
		} else {
			recv_neg_z(aux, l, comm);
			send_neg_z(aux, l, comm);
			recv_pos_z(aux, l, comm);
			send_pos_z(aux, l, comm);
		}
	}

	//Parallel corner
	parallel_corner_x(l, aux, comm);
	parallel_corner_y(l, aux, comm);
	parallel_corner_z(l, aux, comm);

}

