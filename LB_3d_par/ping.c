#include <mpi.h>
#include <stdio.h>
//#include <MPI-rastro.h>

//#define SIZE 0
//#define SIZE 512
#define SIZE 16384
//#define SIZE 8192
//#define SIZE 1048576


int main (int argc, char **argv)
{
	int rank, size, i;
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Request *req = (MPI_Request *) malloc ((size)*sizeof (MPI_Request));
	MPI_Status sta;

	double *buf = (double *) malloc (SIZE * sizeof (double));
	size = 10000;
	
	double tempo = 0;
	
	double menor = 1.0;

	int iter = -1;
	
	if (rank == 0) {
		for (i = 1; i < size; i++){
			double t1 = MPI_Wtime ();
			MPI_Isend (buf,SIZE,MPI_DOUBLE,1,0,MPI_COMM_WORLD, req);
			MPI_Wait (req, &sta);
			//MPI_Send(buf, SIZE,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
			MPI_Recv (buf,SIZE,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&sta);
			double t2 = MPI_Wtime ();
			tempo += t2-t1;
			printf ("MESTRE %d dif: %f\n", i, t2-t1);
			fflush (stdout);
			if(menor > t2-t1) {
				menor = t2-t1;
				iter = i;
			}
			//sleep(1);
		}
		printf("M %d %lf\n", iter, menor/2);
		printf("Tempo %lf\n", (tempo/size)/2);// - menor/2);
	} else {
		for (i = 1; i < size; i++){	
			double t1 = MPI_Wtime ();
			MPI_Recv (buf,SIZE,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&sta);
			MPI_Isend (buf,SIZE,MPI_DOUBLE,0,0,MPI_COMM_WORLD,req);
			MPI_Wait (req, &sta);
			//MPI_Send(buf, SIZE,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
			
			//double t2 = MPI_Wtime ();
			//printf ("ESCRAVO %d dif: %f\n", i, t2-t1);
			//fflush (stdout);
		
			//if(menor > t2-t1) {
			//	menor = t2-t1;
			//	iter = i;
			//}
			
			//sleep(1);
		}
		//printf("E %d %lf\n", iter, menor/2);
	}
	MPI_Finalize ();
	return 0;
}
