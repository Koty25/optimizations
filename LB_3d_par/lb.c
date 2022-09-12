#include "lb.h"

//////////////////////////////////////////
// Functions
//////////////////////////////////////////

//////////////////////////////////////////
// Read_parametrs
//////////////////////////////////////////
s_properties *read_parametrs(string file) {
	//struct with macroscopic properties of a fluid
	s_properties *p = (s_properties*) malloc(sizeof(s_properties));
	if (p == NULL) {
		printf("In: read_parametrs \n");
		exit(11);
	}
	
	//Open input file
	FILE *archive = fopen(file, "r");
	//Checking if the file is correct
	if (archive == NULL) {
		printf("Could not open input file\n\n");
		exit(-2);
	}
	//number of iterations
	fscanf(archive, "%d", &p->t_max);
	//fluid density per link
	fscanf(archive, "%lf", &p->density);
	//density redistribution
	fscanf(archive, "%lf", &p->accel);
	//relaxation parameter
	fscanf(archive, "%lf", &p->omega);
	//linear dimension (for reynolds number)
	fscanf(archive, "%lf", &p->r_rey);
	
	//Close file
	fclose(archive);

	//printf("Parameters read sucessful\n"); //All right
	return p;
}


//////////////////////////////////////////
// Alloc memory to the lattices structures
//////////////////////////////////////////
void alloc_lattice(s_lattice *l) {
	//local variables
	int x, y, z;
	
	//Alloc memory space to the grid 
	//Obstacles matrix
	l->obst = (bool ***) calloc(l->lx, sizeof(bool **));
	if (l->obst == NULL) {
		printf("alloc_lattice \n");
		exit(11);
	}
	for(x = 0; x < l->lx; x++) {
		l->obst[x] = (bool **) calloc(l->ly, sizeof(bool *));
		for(y = 0; y < l->ly; y++) 
			l->obst[x][y] = (bool *) calloc(l->lz, sizeof(bool));
	}
	
	//Lattice and temporary
	l->node = (double ****) calloc(l->lx, sizeof(double***));
	l->temp = (double ****) calloc(l->lx, sizeof(double***));
	for(x = 0; x < l->lx; x++) {
		l->node[x] = (double ***) calloc(l->ly, sizeof(double**));
		l->temp[x] = (double ***) calloc(l->ly, sizeof(double**));
		for(y = 0; y < l->ly; y++) {
			l->node[x][y] = (double **) calloc(l->lz, sizeof(double*));
			l->temp[x][y] = (double **) calloc(l->lz, sizeof(double*));
			for(z = 0; z < l->lz; z++) {
				l->node[x][y][z] = (double *) calloc(l->n, sizeof(double));
				l->temp[x][y][z] = (double *) calloc(l->n, sizeof(double));
			}
		}
	}
}


//////////////////////////////////////////
// Read_obstacles
//////////////////////////////////////////
s_lattice *create_lattice(string file, s_comm *comm) { 
	//local variables
	int max;
	int c = 0;
	int i, j, k;
	//int ii, jj, kk;
	s_lattice *l = (s_lattice *) malloc(sizeof(s_lattice));
	if (l == NULL)
		printf("Erro in create_lattice\n");
		
	//Open input file
	FILE *archive = fopen(file, "r");
	if (archive == NULL) {
		printf("Could not read colision input file\n\n");
                exit(-2);	
	}

	//Reading headers
	fscanf(archive, "%ld", &l->lx_total);
	fscanf(archive, "%ld", &l->ly_total);
	fscanf(archive, "%ld", &l->lz_total);
	fscanf(archive, "%d", &l->d);
	fscanf(archive, "%d", &l->n);
	fscanf(archive, "%d", &max);

	define_sublattice(comm, l);
	//printf("%d %d %d %d %d\n", l->lx, l->ly, l->lz, l->n, max);
	
	//alloc memory
	alloc_lattice(l);

	//Reading obstacle points
	//obstacle[l->lx_total][l->ly_total][l->lz_total];
	while(c < max) {
		fscanf(archive, "%d %d %d", &i, &j, &k);
		
		//Check if i and j are less then x_max and y max
		//if(i > l->lx_total || j > l->ly_total || k > l->lz_total)
			//printf("Obstacle input file is not valid\n\n");
		
		//Very important - All extreme point of the obstacle
		if(l->lx_first <= i && l->lx_last >= i && l->ly_first <= j && l->ly_last >= j && l->lz_first <= k && l->lz_last >= k)
			l->obst[i - l->lx_first + 1][j - l->ly_first + 1][k - l->lz_first + 1] = true;

		c++;
	}
	
	//close archive
	fclose(archive);
	
	return l;
}

//////////////////////////////////////////
// Init_constants
//////////////////////////////////////////
s_aux *init_constants(s_lattice * l, s_properties *p) {
	
	s_aux *aux = (s_aux *) malloc(sizeof(s_aux));
	if (aux == NULL) {
		printf("init_constants\n");
		exit(22);
	}
	int i, j;
	aux->e = (int **) calloc(l->n, sizeof(int *));
	if (aux->e == NULL) {
		printf("init_constants\n");
		exit(22);
	}
	for(i = 0; i < l->n; i++) {
		aux->e[i] = (int *) calloc(l->d, sizeof(int));
		if (aux->e[i] == NULL) {
			printf("init_constants\n");
			exit(22);
		}
	}
	//Lattice vectors for the particles in XY-plane
	aux->e[0][0] =  0; aux->e[0][1] =  0;  aux->e[0][2] = 0;  

	aux->e[1][0] =  1; aux->e[1][1] =  0;  aux->e[1][2] = 0;	  
	aux->e[2][0] =  1; aux->e[2][1] =  1;  aux->e[2][2] = 0;
	aux->e[3][0] =  0; aux->e[3][1] =  1;  aux->e[3][2] = 0;
	aux->e[4][0] = -1; aux->e[4][1] =  1;  aux->e[4][2] = 0;
	aux->e[5][0] = -1; aux->e[5][1] =  0;  aux->e[5][2] = 0;
	aux->e[6][0] = -1; aux->e[6][1] = -1;  aux->e[6][2] = 0;
	aux->e[7][0] =  0; aux->e[7][1] = -1;  aux->e[7][2] = 0;
	aux->e[8][0] =  1; aux->e[8][1] = -1;  aux->e[8][2] = 0;

	// Lattice vectors for the particles in XZ-plane
	aux->e[9][0]  =  1; aux->e[9][1]  =  0;  aux->e[9][2]  =  1;
	aux->e[10][0] =  0; aux->e[10][1] =  0;  aux->e[10][2] =  1;
	aux->e[11][0] = -1; aux->e[11][1] =  0;  aux->e[11][2] =  1;
	aux->e[12][0] = -1; aux->e[12][1] =  0;  aux->e[12][2] = -1;
	aux->e[13][0] =  0; aux->e[13][1] =  0;  aux->e[13][2] = -1;
	aux->e[14][0] =  1; aux->e[14][1] =  0;  aux->e[14][2] = -1;
 
	//Lattice vectors for the particles in YZ-plane

	aux->e[15][0] =  0; aux->e[15][1] =  1;  aux->e[15][2] =  1;
	aux->e[16][0] =  0; aux->e[16][1] = -1;  aux->e[16][2] =  1;
	aux->e[17][0] =  0; aux->e[17][1] = -1;  aux->e[17][2] = -1;
	aux->e[18][0] =  0; aux->e[18][1] =  1;  aux->e[18][2] = -1;

	double *t = (double *) calloc(l->d, sizeof(double));
	t[0] = p->density/3;
	t[1] = p->density/18;
	t[2] = p->density/36;
	
	//printf("%4.3f %4.3f %4.3f\n", t[0], t[1], t[2]);
	
	//The coefficients for eq. distr. A is for the constant term, 
	//B is for e*u-term, C is for (e*u)^2-term and D is for u^2-term 
	aux->A = (double *) calloc(l->n, sizeof(double));
	aux->B = (double *) calloc(l->n, sizeof(double));
	aux->C = (double *) calloc(l->n, sizeof(double));
	aux->D = (double *) calloc(l->n, sizeof(double));


	//Coefficients for the particles in XY-plane
	aux->A[0] = t[0]; aux->B[0] = 0.0;    aux->C[0] = 0.0;        aux->D[0] = -(t[0]/2.0)/CS2;
	
	aux->A[1] = t[1]; aux->B[1] = t[1]/CS2; aux->C[1] = (t[1]/2.0)/CS4; aux->D[1] = (-t[1]/2.0)/CS2;
	aux->A[2] = t[2]; aux->B[2] = t[2]/CS2; aux->C[2] = (t[2]/2.0)/CS4; aux->D[2] = (-t[2]/2.0)/CS2;
	aux->A[3] = t[1]; aux->B[3] = t[1]/CS2; aux->C[3] = (t[1]/2.0)/CS4; aux->D[3] = (-t[1]/2.0)/CS2;
	aux->A[4] = t[2]; aux->B[4] = t[2]/CS2; aux->C[4] = (t[2]/2.0)/CS4; aux->D[4] = (-t[2]/2.0)/CS2;
	aux->A[5] = t[1]; aux->B[5] = t[1]/CS2; aux->C[5] = (t[1]/2.0)/CS4; aux->D[5] = (-t[1]/2.0)/CS2;
	aux->A[6] = t[2]; aux->B[6] = t[2]/CS2; aux->C[6] = (t[2]/2.0)/CS4; aux->D[6] = (-t[2]/2.0)/CS2;
	aux->A[7] = t[1]; aux->B[7] = t[1]/CS2; aux->C[7] = (t[1]/2.0)/CS4; aux->D[7] = (-t[1]/2.0)/CS2;
	aux->A[8] = t[2]; aux->B[8] = t[2]/CS2; aux->C[8] = (t[2]/2.0)/CS4; aux->D[8] = (-t[2]/2.0)/CS2;

	//Coefficients for the particles in XZ-plane
	aux->A[9]  = t[2]; aux->B[9]  = t[2]/CS2; aux->C[9]  = (t[2]/2.0)/CS4; aux->D[9]  = (-t[2]/2.0)/CS2;
	aux->A[10] = t[1]; aux->B[10] = t[1]/CS2; aux->C[10] = (t[1]/2.0)/CS4; aux->D[10] = (-t[1]/2.0)/CS2;
	aux->A[11] = t[2]; aux->B[11] = t[2]/CS2; aux->C[11] = (t[2]/2.0)/CS4; aux->D[11] = (-t[2]/2.0)/CS2;
	aux->A[12] = t[2]; aux->B[12] = t[2]/CS2; aux->C[12] = (t[2]/2.0)/CS4; aux->D[12] = (-t[2]/2.0)/CS2;
	aux->A[13] = t[1]; aux->B[13] = t[1]/CS2; aux->C[13] = (t[1]/2.0)/CS4; aux->D[13] = (-t[1]/2.0)/CS2;
	aux->A[14] = t[2]; aux->B[14] = t[2]/CS2; aux->C[14] = (t[2]/2.0)/CS4; aux->D[14] = (-t[2]/2.0)/CS2;

	//Coefficients for the particles in YZ-plane
	aux->A[15] = t[2]; aux->B[15] = t[2]/CS2; aux->C[15] = (t[2]/2.0)/CS4; aux->D[15] = (-t[2]/2.0)/CS2;
	aux->A[16] = t[2]; aux->B[16] = t[2]/CS2; aux->C[16] = (t[2]/2.0)/CS4; aux->D[16] = (-t[2]/2.0)/CS2;
	aux->A[17] = t[2]; aux->B[17] = t[2]/CS2; aux->C[17] = (t[2]/2.0)/CS4; aux->D[17] = (-t[2]/2.0)/CS2;
	aux->A[18] = t[2]; aux->B[18] = t[2]/CS2; aux->C[18] = (t[2]/2.0)/CS4; aux->D[18] = (-t[2]/2.0)/CS2;
	
	//The six tables that tell the numbers of ples that have positive or negative x, y or z component. Each of these tables have 5 (CLNBR) components
	aux->pos_x = (int *) calloc(COL, sizeof(int));
	aux->neg_x = (int *) calloc(COL, sizeof(int));
	aux->pos_y = (int *) calloc(COL, sizeof(int));
	aux->neg_y = (int *) calloc(COL, sizeof(int));
	aux->pos_z = (int *) calloc(COL, sizeof(int));
	aux->neg_z = (int *) calloc(COL, sizeof(int));
	
	i = 0;
	for(j=1; j<NDIM; j++) {
		if(aux->e[j][0] > 0) {
			aux->pos_x[i] = j;
			i++;
		}
	}
	i = 0;
	for(j=1; j<NDIM; j++) {
		if(aux->e[j][0] < 0) {
			aux->neg_x[i] = j;
			i++;
		}
	}
	i = 0;
	for(j=1; j<NDIM; j++) {
		if(aux->e[j][1] > 0) {
			aux->pos_y[i] = j;
			i++;
		}
	}
	i = 0;
	for(j=1; j<NDIM; j++) {
		if(aux->e[j][1] < 0) {
			aux->neg_y[i] = j;
			i++;
		}
	}
	i = 0;
	for(j=1; j<NDIM; j++) {
		if(aux->e[j][2] > 0) {
			aux->pos_z[i] = j;
			i++;
		}
	}
	i = 0;
	for(j=1; j<NDIM; j++) {
		if(aux->e[j][2] < 0) {
			aux->neg_z[i] = j;
			i++;
		}
	}
	aux->points = (int ***) calloc(l->d, sizeof(int **));
	for(i = 0; i < l->d; i++) {
		aux->points[i] = (int **) calloc(l->d, sizeof(int *));
		for(j = 0; j < l->d; j++) 
			aux->points[i][j] = (int *) calloc(l->d, sizeof(int));
	}

	//Tables that give the population numbers when the direction is known
	//0 = -component, 1 = no component 2 = +component
	aux->points[0+1][0+1][0+1]   = 0;  aux->points[1+1][0+1][0+1]   = 1;
	aux->points[1+1][1+1][0+1]   = 2;  aux->points[0+1][1+1][0+1]   = 3; 
	aux->points[-1+1][1+1][0+1]  = 4;  aux->points[-1+1][0+1][0+1]  = 5; 
	aux->points[-1+1][-1+1][0+1] = 6;  aux->points[0+1][-1+1][0+1]  = 7; 
	aux->points[1+1][-1+1][0+1]  = 8;  aux->points[1+1][0+1][1+1]   = 9;
	aux->points[0+1][0+1][1+1]   = 10; aux->points[-1+1][0+1][1+1]  = 11; 
	aux->points[-1+1][0+1][-1+1] = 12; aux->points[0+1][0+1][-1+1]  = 13; 
	aux->points[1+1][0+1][-1+1]  = 14; aux->points[0+1][1+1][1+1]   = 15; 
	aux->points[0+1][-1+1][1+1]  = 16; aux->points[0+1][-1+1][-1+1] = 17; 
	aux->points[0+1][1+1][-1+1]  = 18; 

	return aux;
}


//////////////////////////////////////////
// Init_density
//////////////////////////////////////////
void init_density(s_lattice * l, double *A) {
	//local variables
	int x, y, z, n;
	
	//loop over computational domain
	for (x = 1; x < l->lx - 1; x++) {
		for (y = 1; y < l->ly - 1; y++) {
			for (z = 1; z < l->lz - 1; z++) {
				for (n = 0; n < l->n; n++) {
					l->node[x][y][z][n] = A[n];
				}
			}
		}
	}
}



//////////////////////////////////////////
// Check_density
//////////////////////////////////////////
double check_density(s_lattice *l) { //allreduce: sum of all densities, but not consider de boundary densities of the sublattices
	//local variables
	int x, y, z, n;
	double n_sum = 0;
	
	for (x = 1; x < l->lx - 1; x++) {
		for (y = 1; y < l->ly - 1; y++) {
			for (z = 1; z < l->lz - 1; z++) {
				for (n = 0; n < l->n; n++) {
					n_sum += l->node[x][y][z][n];
				} 
			} 
		} 
	} 
	
	return n_sum;
}


//////////////////////////////////////////
// Redistribute
//////////////////////////////////////////
// It is interesting to redistribute de forces to all points
void redistribute(s_lattice *l, double accel, double density) {
	//local variables
	int x, y, z, n;
	double t_1 = density * accel / 18.0;
	double t_2 = density * accel / 36.0;
	int v_1[NDIM], v_2[NDIM];
	int length = 1;
	
	v_1[0] = 1;
	v_2[0] = 5;
	x=0;
	for (x = 1; x < l->lx - 1; x++) {
		for (y = 1; y < l->ly - 1; y++) {
			for (z = 1; z < l->lz - 1; z++) {
				//check to avoid negative densities
				//check false | true
				if (l->obst[x][y][z] == false) {
					for (n = 0; n < length; n++) {
						l->node[x][y][z][v_1[n]] += t_1;
						l->node[x][y][z][v_2[n]] -= t_1;
					}
				}
			}
		}
	}

	v_1[0] = 2;
	v_1[1] = 8;
	v_1[2] = 9;
	v_1[3] = 14;

	v_2[0] = 4;
	v_2[1] = 6;
	v_2[2] = 11;
	v_2[3] = 12;

	length = 4;
	for (x = 0; x < l->lx; x++) {
		for (y = 0; y < l->ly; y++) {
			for (z = 0; z < l->lz; z++) {
				//check to avoid negative densities
				//check false | true
				if (l->obst[x][y][z] == false) {
					for (n = 0; n < length; n++) {
						l->node[x][y][z][v_1[n]] += t_2;
						l->node[x][y][z][v_2[n]] -= t_2;
					}
				}				
			}
		}
	}
}


//////////////////////////////////////////
// Propagate
//////////////////////////////////////////
void propagate(s_lattice *l) {
        //local variables
	int x, y, z;
		
	for (x = 1; x < l->lx - 1; x++) {
		for(y = 1; y < l->ly - 1; y++) {
			for(z = 1; z < l->lz - 1; z++) {
							
				l->temp[x][y][z][0] = l->node[x][y][z][0];

				l->temp[x+1][y][z][1] = l->node[x][y][z][1];
				l->temp[x+1][y+1][z][2] = l->node[x][y][z][2];
				l->temp[x][y+1][z][3] = l->node[x][y][z][3];
				l->temp[x-1][y+1][z][4] = l->node[x][y][z][4];
				l->temp[x-1][y][z][5] = l->node[x][y][z][5];
				l->temp[x-1][y-1][z][6] = l->node[x][y][z][6];
				l->temp[x][y-1][z][7] = l->node[x][y][z][7];
				l->temp[x+1][y-1][z][8] = l->node[x][y][z][8];

				l->temp[x+1][y][z+1][9] = l->node[x][y][z][9];
				l->temp[x][y][z+1][10] = l->node[x][y][z][10];
				l->temp[x-1][y][z+1][11] = l->node[x][y][z][11];
				l->temp[x-1][y][z-1][12] = l->node[x][y][z][12];
				l->temp[x][y][z-1][13] = l->node[x][y][z][13];
				l->temp[x+1][y][z-1][14] = l->node[x][y][z][14];
				
				l->temp[x][y+1][z+1][15] = l->node[x][y][z][15];
				l->temp[x][y-1][z+1][16] = l->node[x][y][z][16];
				l->temp[x][y-1][z-1][17] = l->node[x][y][z][17];
				l->temp[x][y+1][z-1][18] = l->node[x][y][z][18];
			}
		}
	}
}


//////////////////////////////////////////
// Bounceback
//////////////////////////////////////////
void bounceback(s_lattice *l) {
	//local variables
	int x, y, z;

	for (x = 1; x < l->lx - 1; x++) {
		for(y = 1; y < l->ly - 1; y++) {
			for(z = 1; z < l->lz - 1; z++) {
				if (l->obst[x][y][z] == true) {
					
					l->node[x][y][z][1] = l->temp[x][y][z][5];
					l->node[x][y][z][2] = l->temp[x][y][z][6];
					l->node[x][y][z][3] = l->temp[x][y][z][7];
					l->node[x][y][z][4] = l->temp[x][y][z][8];
					l->node[x][y][z][5] = l->temp[x][y][z][1];
					l->node[x][y][z][6] = l->temp[x][y][z][2];
					l->node[x][y][z][7] = l->temp[x][y][z][3];
					l->node[x][y][z][8] = l->temp[x][y][z][4];

					l->node[x][y][z][9] = l->temp[x][y][z][12];
					l->node[x][y][z][10] = l->temp[x][y][z][13];
					l->node[x][y][z][11] = l->temp[x][y][z][14];
					l->node[x][y][z][12] = l->temp[x][y][z][9];
					l->node[x][y][z][13] = l->temp[x][y][z][10];
					l->node[x][y][z][14] = l->temp[x][y][z][11];
					
					l->node[x][y][z][15] = l->temp[x][y][z][17];
					l->node[x][y][z][16] = l->temp[x][y][z][18];
					l->node[x][y][z][17] = l->temp[x][y][z][15];
					l->node[x][y][z][18] = l->temp[x][y][z][16];
				}
			}
		}
	}
} 


//////////////////////////////////////////
// Relaxation
//////////////////////////////////////////
void relaxation(s_lattice *l, double density, double omega, s_aux *aux) {
	//local variables
	int x, y, z, i;
	double u_x, u_y, u_z;
	double u_n[l->n], n_equ[l->n], u_squ, d_loc;

	for (x = 1; x < l->lx - 1; x++) {
		for(y = 1; y < l->ly - 1; y++) {
			for(z = 1; z < l->lz - 1; z++) {
				if (l->obst[x][y][z] == false) {
					//d_loc = 0.0;	
					//for (i = 0; i < l->n; i++) {
					//	d_loc += l->node[x][y][z][i];
					//}
					d_loc = l->temp[x][y][z][0] + l->temp[x][y][z][1] + l->temp[x][y][z][2] + l->temp[x][y][z][3] + l->temp[x][y][z][4] + l->temp[x][y][z][5] + l->temp[x][y][z][6] + l->temp[x][y][z][7] + l->temp[x][y][z][8] + l->temp[x][y][z][9] + l->temp[x][y][z][10] + l->temp[x][y][z][11] + l->temp[x][y][z][12] + l->temp[x][y][z][13] + l->temp[x][y][z][14] + l->temp[x][y][z][15] + l->temp[x][y][z][16] + l->temp[x][y][z][17] + l->temp[x][y][z][18];
					
					//x-, y- and z- velocity components
					u_x = l->temp[x][y][z][1] + l->temp[x][y][z][2] + l->temp[x][y][z][8] + l->temp[x][y][z][9] + l->temp[x][y][z][14]  - (l->temp[x][y][z][4] + l->temp[x][y][z][5] + l->temp[x][y][z][6] + l->temp[x][y][z][11] + l->temp[x][y][z][12]);
					u_y = l->temp[x][y][z][2] + l->temp[x][y][z][3] + l->temp[x][y][z][4] + l->temp[x][y][z][15] + l->temp[x][y][z][18]  - (l->temp[x][y][z][6] + l->temp[x][y][z][7] + l->temp[x][y][z][8] + l->temp[x][y][z][16] + l->temp[x][y][z][17]); 
					u_z = l->temp[x][y][z][9] + l->temp[x][y][z][10] + l->temp[x][y][z][11] + l->temp[x][y][z][15] + l->temp[x][y][z][16]  - (l->temp[x][y][z][12] + l->temp[x][y][z][13] + l->temp[x][y][z][14] + l->temp[x][y][z][17] + l->temp[x][y][z][18]);
					u_x = u_x / d_loc;
					u_y = u_y / d_loc;
					u_z = u_z / d_loc;
					
					//square velocity
					u_squ = u_x * u_x + u_y * u_y + u_z * u_z;
					//n- velocity compnents
					//only 3 speeds would be necessary
					u_n[1] = u_x;
					u_n[2] = u_x + u_y;
					u_n[3] = u_y;
					u_n[4] = -u_x + u_y;
					u_n[5] = -u_x;
					u_n[6] = -u_x - u_y;
					u_n[7] = -u_y;
					u_n[8] = u_x - u_y;

					u_n[9] = u_x + u_z;
					u_n[10] = u_z;
					u_n[11] = -u_x + u_z;
					u_n[12] = -u_x - u_z;
					u_n[13] = -u_z;
					u_n[14] = u_x - u_z;
					u_n[15] = u_y + u_z;
					u_n[16] = -u_y + u_z;
					u_n[17] = -u_y - u_z;
					u_n[18] = u_y - u_z;

					//Equilibrium distribution function
					for (i = 0; i < l->n; i++)
						n_equ[i] = d_loc * (aux->A[i] + u_n[i] * aux->B[i] + u_n[i] * u_n[i] * aux->C[i] + u_squ * aux->D[i]) / density;

					//relaxation step
					for (i = 0; i < l->n; i++) {
						l->node[x][y][z][i] = l->temp[x][y][z][i] + omega * (n_equ[i] - l->temp[x][y][z][i]);
					}
				}
			}
		}
	}
}


//////////////////////////////////////////
// Calc_velocity
//////////////////////////////////////////
double *calc_velocity(s_lattice *l, int time) { 
	//local variables
	int x, y, z, i, n_free;
	double u_x, d_loc;

	x = l->lx_total/2 - l->lx_first + 1;
	
	n_free = 0;
	u_x = 0;

	for(y = 1; y < l->ly - 1; y++) {
		for(z = 1; z < l->lz - 1; z++) {
			if (l->obst[x][y][z] == false) {
				d_loc = 0;	
				for (i = 0; i < l->n; i++)
					d_loc = d_loc + l->node[x][y][z][i];
				u_x = u_x + (l->node[x][y][z][1] + l->node[x][y][z][2] + l->node[x][y][z][8] + l->node[x][y][z][9] + l->node[x][y][z][14]  - (l->node[x][y][z][4] + l->node[x][y][z][5] + l->node[x][y][z][6] + l->node[x][y][z][11] + l->node[x][y][z][12])) / d_loc;
				n_free++;
			}
		}
	}
	//Optional
	//printf("%d %lf\n", time, u_x / n_free);
	//return u_x / n_free;
	double *ret = malloc(2 * sizeof(double));
	ret[0] = u_x;
	ret[1] = n_free;
	return ret;
}


////////////////////////////////////////////
//// Write_results
////////////////////////////////////////////
void write_results(string file, s_lattice *l, double density) {
	//local variables
	int x, y, z, i;
	bool obsval;
	double u_x, u_y, u_z, d_loc, press;

	//Square speed of sound
	double c_squ = 1.0 / 3.0;

	//Open results output file
	FILE *archive = fopen(file, "w");

	//write results
	fprintf(archive,"VARIABLES = X, Y, Z, VX, VY, VZ, PRESS, OBST\n");
	fprintf(archive,"ZONE I= %ld, J= %ld, K= %ld F=POINT\n", l->lx - 2, l->ly - 2, l->lz - 2);

	for(z = 1; z < l->lz - 1; z++) {
		for(y = 1; y < l->ly - 1; y++) {
			for(x = 1; x < l->lx - 1; x++) {
				//if obstacle node, nothing is to do
				if (l->obst[x][y][z] == true) {
					//obstacle indicator
					obsval = true;
					//velocity components = 0
					u_x = 0.0;
					u_y = 0.0;
					u_z = 0.0;
					//pressure = average pressure
					press = density * c_squ;
				} else {
					//integral local density
					//initialize variable d_loc
					d_loc = 0.0;
					for (i = 0; i < l->n; i++) {
						d_loc += l->node[x][y][z][i];
					}
					
					//x-, y- and z- velocity components
					u_x = l->node[x][y][z][1] + l->node[x][y][z][2] + l->node[x][y][z][8] + l->node[x][y][z][9] + l->node[x][y][z][14]  - (l->node[x][y][z][4] + l->node[x][y][z][5] + l->node[x][y][z][6] + l->node[x][y][z][11] + l->node[x][y][z][13]);
					u_x /= d_loc;

					u_y = l->node[x][y][z][2] + l->node[x][y][z][3] + l->node[x][y][z][4] + l->node[x][y][z][15] + l->node[x][y][z][18]  - (l->node[x][y][z][6] + l->node[x][y][z][7] + l->node[x][y][z][8] + l->node[x][y][z][16] + l->node[x][y][z][17]);
					u_y /= d_loc;
					
					u_z = l->node[x][y][z][9] + l->node[x][y][z][10] + l->node[x][y][z][11] + l->node[x][y][z][15] + l->node[x][y][z][16]  - (l->node[x][y][z][12] + l->node[x][y][z][13] + l->node[x][y][z][14] + l->node[x][y][z][17] + l->node[x][y][z][18]);
					u_z /= d_loc;
	
					//pressure
					press = d_loc * c_squ;
					obsval = false;
				}
				fprintf(archive,"%d %d %d %f %f %f %f %d\n", x - 1, y - 1, z - 1, u_x, u_y, u_z, press, obsval);
			}
		}
	}
	
	fclose(archive);
}


//////////////////////////////////////////
// Compute Reynolds number
//////////////////////////////////////////
void comp_rey(double vel, s_properties *p, int time, double execution, double communication, string dir, int proc_x, int proc_y, int proc_z) {

	//Local variables
	double visc, rey;

	//Compute viscosity
	visc = 1.0 / 6.0 * (2.0 / p->omega - 1.0);
	
	//Compute Reynolds number
	rey = vel * p->r_rey / visc;

	//Messages
	//printf("Calculations finished, results:\n");
	
	//Create directory if it don't exist
	string command = (string) malloc(20*sizeof(string));
	strcpy(command, "mkdir ");
	strcat(command, dir);
	system(command);
	
	//Create name of output files directory + /name + .out
	//directory is a string argument from main function
	//Generical informations
	string tt = (string) malloc(100*sizeof(string));
	strcpy(tt, dir);
	strcat(tt, "/time.out");
	
	//Communication time
	string cc = (string) malloc(100*sizeof(string));
	strcpy(cc, dir);
	strcat(cc, "/comm.out");
	
	//Execution time
	string ee = (string) malloc(100*sizeof(string));
	strcpy(ee, dir);
	strcat(ee, "/exec.out");
	
	//printf("Creating files...\n%s \n%s \n%s\n", tt, cc, ee);
	
	FILE *t = fopen(tt, "a");
	FILE *e = fopen(ee, "a");
	FILE *c = fopen(cc, "a");
	if(t == NULL)
		printf("Erro\n");
	//printf("%s\n", dir);

	fprintf(t, "%d %d %d %d %lg %lg %lg\n", proc_x, proc_y, proc_z, time, visc, vel, rey);
	fprintf(e, "%lg\n", execution);
	fprintf(c, "%lg\n", communication);
	
	//Close files
	fclose(t);
	fclose(e);
	fclose(c);

	free(tt);
	free(ee);
	free(cc);
}


//////////////////////////////////////////
// Return de time in a especific moment
//////////////////////////////////////////
double crono() {
	struct timeval tv;
	gettimeofday(&tv, 0);
	return tv.tv_sec + tv.tv_usec / 1e6;
}
