NEC C/C++ Compiler (3.0.4) for Vector Engine     Thu Jul 23 22:24:20 2020
FILE NAME: c_print_results.c

FUNCTION NAME: c_print_results
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    38: par(1801): Parallel routine generated.: c_print_results$1
    40: par(1805): Master-construct generated.


NEC C/C++ Compiler (3.0.4) for Vector Engine     Thu Jul 23 22:24:20 2020
FILE NAME: c_print_results.c

FUNCTION NAME: c_print_results
FORMAT LIST

 LINE   LOOP      STATEMENT

    10:           void c_print_results( char   *name,
    11:                                 char   class,
    12:                                 int    n1,
    13:                                 int    n2,
    14:                                 int    n3,
    15:                                 int    niter,
    16:                                 double t,
    17:                                 double mops,
    18:           		      char   *optype,
    19:                                 int    passed_verification,
    20:                                 char   *npbversion,
    21:                                 char   *compiletime,
    22:                                 char   *cc,
    23:                                 char   *clink,
    24:                                 char   *c_lib,
    25:                                 char   *c_inc,
    26:                                 char   *cflags,
    27:                                 char   *clinkflags )
    28:           {
    29:               int num_threads, max_threads;
    30:           
    31:           
    32:               max_threads = 0;
    33:               num_threads = 0;
    34:           
    35:           /*   figure out number of threads used */
    36:           #ifdef _OPENMP
    37:               max_threads = omp_get_max_threads();
    38:           #pragma omp parallel shared(num_threads)
    39:           {
    40:               #pragma omp master
    41:               num_threads = omp_get_num_threads();
    42:           }
    43:           #endif
    44:           
    45:           
    46:               printf( "\n\n %s Benchmark Completed\n", name );
    47:           
    48:               printf( " Class           =                        %c\n", class );
    49:           
    50:               if( n3 == 0 ) {
    51:                   long nn = n1;
    52:                   if ( n2 != 0 ) nn *= n2;
    53:                   printf( " Size            =             %12ld\n", nn );   /* as in IS */
    54:               }
    55:               else
    56:                   printf( " Size            =             %4dx%4dx%4d\n", n1,n2,n3 );
    57:           
    58:               printf( " Iterations      =             %12d\n", niter );
    59:            
    60:               printf( " Time in seconds =             %12.2f\n", t );
    61:           
    62:               if (num_threads > 0)
    63:                   printf( " Total threads   =             %12d\n", num_threads);
    64:           
    65:               if (max_threads > 0)
    66:                   printf( " Avail threads   =             %12d\n", max_threads);
    67:           
    68:               if (num_threads != max_threads)
    69:                   printf( " Warning: Threads used differ from threads available\n");
    70:           
    71:               printf( " Mop/s total     =             %12.2f\n", mops );
    72:           
    73:               if (num_threads > 0)
    74:                   printf( " Mop/s/thread    =             %12.2f\n",
    75:                          mops/(double)num_threads );
    76:           
    77:               printf( " Operation type  = %24s\n", optype);
    78:           
    79:               if( passed_verification < 0 )
    80:                   printf( " Verification    =            NOT PERFORMED\n" );
    81:               else if( passed_verification )
    82:                   printf( " Verification    =               SUCCESSFUL\n" );
    83:               else
    84:                   printf( " Verification    =             UNSUCCESSFUL\n" );
    85:           
    86:               printf( " Version         =             %12s\n", npbversion );
    87:           
    88:               printf( " Compile date    =             %12s\n", compiletime );
    89:           
    90:               printf( "\n Compile options:\n" );
    91:           
    92:               printf( "    CC           = %s\n", cc );
    93:           
    94:               printf( "    CLINK        = %s\n", clink );
    95:           
    96:               printf( "    C_LIB        = %s\n", c_lib );
    97:           
    98:               printf( "    C_INC        = %s\n", c_inc );
    99:           
   100:               printf( "    CFLAGS       = %s\n", cflags );
   101:           
   102:               printf( "    CLINKFLAGS   = %s\n", clinkflags );
   103:           
   104:               printf( "\n\n" );
   105:               printf( " Please send all errors/feedbacks to:\n\n" );
   106:               printf( " NPB Development Team\n" );
   107:               printf( " npb@nas.nasa.gov\n\n\n" );
   108:           /*    printf( " Please send the results of this run to:\n\n" );
   109:               printf( " NPB Development Team\n" );
   110:               printf( " Internet: npb@nas.nasa.gov\n \n" );
   111:               printf( " If email is not available, send this to:\n\n" );
   112:               printf( " MS T27A-1\n" );
   113:               printf( " NASA Ames Research Center\n" );
   114:               printf( " Moffett Field, CA  94035-1000\n\n" );
   115:               printf( " Fax: 650-604-3957\n\n" ); */
   116:           }


