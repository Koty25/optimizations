#/usr/local/mpich-gm/bin/
mpicc -o lb comm.c lb.c main.c -lm -Wall -Lefence
