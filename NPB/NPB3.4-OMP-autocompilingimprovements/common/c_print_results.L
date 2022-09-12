NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME : c_print_results.c

  COMPILER OPTIONS : -c -O4 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2 -finline-functions -fouterloop-unroll

  PARAMETER : 

      Optimization Options : 
        -On                                      : 4
        -fargument-alias                         : disable
        -fargument-noalias                       : enable
        -fassociative-math                       : enable
        -faggressive-associative-math            : disable
        -fcse-after-vectorization                : disable
        -ffast-math                              : enable
        -fignore-induction-variable-overflow     : enable
        -fignore-volatile                        : enable
        -fivdep                                  : disable
        -floop-collapse                          : enable
        -floop-count                             : 5000
        -floop-fusion                            : enable
        -floop-interchange                       : enable
        -floop-normalize                         : enable
        -floop-split                             : disable
        -floop-strip-mine                        : enable
        -floop-unroll                            : enable
        -floop-unroll-complete                   : 4
        -floop-unroll-max-times                  : auto
        -fmatrix-multiply                        : enable
        -fmove-loop-invariants                   : enable
        -fmove-loop-invariants-if                : enable
        -fmove-loop-invariants-unsafe            : enable
        -fmove-nested-loop-invariants-outer      : enable
        -fnaked-ivdep                            : disable
        -fnamed-alias                            : disable
        -fnamed-noalias                          : enable
        -fouterloop-unroll                       : enable
        -fouterloop-unroll-max-size              : 4
        -fouterloop-unroll-max-times             : 4
        -fprecise-math                           : disable
        -freciprocal-math                        : enable
        -freplace-loop-equation                  : enable
        -fstrict-aliasing                        : enable
        -fthis-pointer-alias                     : disable
        -fthis-pointer-noalias                   : enable
        -mlist-vector                            : disable
        -mretain-all                             : enable
        -mretain-list-vector                     : disable
        -mretain-none                            : disable
        -msched-none                             : disable
        -msched-insns                            : enable
        -msched-block                            : disable
        -mvector                                 : enable
        -mvector-advance-gather                  : enable
        -mvector-advance-gather-limit            : 56
        -mvector-dependency-test                 : enable
        -mvector-floating-divide-instruction     : disable
        -mvector-fma                             : enable
        -mvector-intrinsic-check                 : disable
        -mvector-iteration                       : enable
        -mvector-iteration-unsafe                : disable
        -mvector-loop-count-test                 : disable
        -mvector-low-precise-divide-function     : disable
        -mvector-merge-conditional               : enable
        -mvector-packed                          : disable
        -mvector-power-to-explog                 : disable
        -mvector-power-to-sqrt                   : enable
        -mvector-reduction                       : enable
        -mvector-shortloop-reduction             : disable
        -mvector-sqrt-instruction                : disable
        -mvector-threshold                       : 5
        -mwork-vector-kind=none                  : disable

      Parallelization Options : 
        -fopenmp                                 : enable
        -mparallel                               : disable
        -mparallel-innerloop                     : disable
        -mparallel-omp-routine                   : enable
        -mparallel-outerloop-strip-mine          : disable
        -mparallel-sections                      : disable
        -mparallel-threshold                     : 2000
        -mschedule-dynamic                       : disable
        -mschedule-runtime                       : disable
        -mschedule-static                        : disable
        -mschedule-chunk-size                    : 0
        -pthread                                 : disable

      Inlining Options : 
        -finline                                 : enable
        -finline-copy-arguments                  : disable
        -finline-functions                       : enable
        -finline-max-depth                       : 2
        -finline-max-function-size               : 50
        -finline-max-times                       : 6

      Code Generation Options : 
        -finstrument-functions                   : disable
        -fpic | -fPIC                            : disable
        -ftrace                                  : disable
        -p                                       : disable

      Debugging Options : 
        -g                                       : disable
        -minit-stack                             : not specified
        -traceback                               : disable

      Language Options : 
        -fallow-keyword-macros                   : enable
        -fgnu89-inline                           : disable
        -frestrict                               : enable
        -fsigned-char                            : enable
        -std                                     : gnu11

      Message Options : 
        -Wall                                    : disable
        -Wcomment                                : disable
        -Werror                                  : disable
        -Wno-div-by-zero                         : disable
        -Wunknown-pragma                         : disable
        -Wunused                                 : disable
        -Wunused-but-set-parameter               : disable
        -Wunused-but-set-variable                : disable
        -Wunused-parameter                       : disable
        -Wunused-value                           : disable
        -Wunused-variable                        : disable
        -fdiag-inline                            : 2
        -fdiag-parallel                          : 2
        -fdiag-vector                            : 2
        -fdiag-system-header                     : disable
        -pedantic                                : disable
        -pedantic-errors                         : disable
        -w                                       : disable

      List Output Options : 
        -report-file                             : not specified
        -report-append-mode                      : disable
        -report-all                              : enable
        -report-cg                               : enable
        -report-diagnostics                      : enable
        -report-format                           : enable
        -report-inline                           : enable
        -report-option                           : enable
        -report-vector                           : enable

      Miscellaneous Options : 
        -v                                       : disable


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_print_results.c

FUNCTION NAME: c_print_results
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

    38: par(1801): Parallel routine generated.: c_print_results$1
    40: par(1805): Master-construct generated.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_print_results.c

FUNCTION NAME: c_print_results
FORMAT LIST

 Loop Mark                              Line Mark
 ---------                              ---------
 C - Conditionally Vectorized           C - Vector Scatter
 P - Parallelized                       F - Fused-multiply-add
 S - Partially Vectorized               G - Vector Gather
 U - Unrolled                           I - Inlined
 V - Vectorized                         M - Vector Matrix Multiply
 W - Collapsed and Vectorized           R - Retain
 Y - Parallelized and Vectorized        V - Vreg
 X - Interchanged and Vectorized        
 + - Not Vectorized                     
 * - Expanded                           

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


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_print_results.c

FUNCTION NAME: c_print_results
INLINE LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_print_results.c

FUNCTION NAME: c_print_results
VECTORIZATION LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Fri May 21 01:36:50 2021
FILE NAME: c_print_results.c

FUNCTION NAME: c_print_results
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 13 [s0 s2 s8-s12 s15-s16 s23 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 208 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :  16 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


