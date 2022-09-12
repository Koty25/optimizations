NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME : is.c

  COMPILER OPTIONS : -c -O3 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2

  PARAMETER : 

      Optimization Options : 
        -On                                      : 3
        -fargument-alias                         : enable
        -fargument-noalias                       : disable
        -fassociative-math                       : enable
        -faggressive-associative-math            : disable
        -fcse-after-vectorization                : disable
        -ffast-math                              : enable
        -fignore-induction-variable-overflow     : disable
        -fignore-volatile                        : disable
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
        -fmove-loop-invariants-unsafe            : disable
        -fmove-nested-loop-invariants-outer      : enable
        -fnaked-ivdep                            : disable
        -fnamed-alias                            : disable
        -fnamed-noalias                          : enable
        -fouterloop-unroll                       : enable
        -fouterloop-unroll-max-size              : 4
        -fouterloop-unroll-max-times             : 4
        -fprecise-math                           : disable
        -freciprocal-math                        : enable
        -freplace-loop-equation                  : disable
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
        -finline-copy-arguments                  : enable
        -finline-functions                       : disable
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


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: randlc
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   333: vec( 101): Vectorized loop.
   335: vec( 126): Idiom detected.: Product
   336: vec( 126): Idiom detected.: Product
   338: vec( 101): Vectorized loop.
   340: vec( 126): Idiom detected.: Product
   341: vec( 126): Idiom detected.: Product


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: randlc
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

   316:           double	randlc( double *X, double *A )
   317:           {
   318:                 double		T1, T2, T3, T4;
   319:                 double		A1;
   320:                 double		A2;
   321:                 double		X1;
   322:                 double		X2;
   323:                 double		Z;
   324:                 int     		i, j;
   325:           
   326:                 if (KS == 0)
   327:                 {
   328:                   R23 = 1.0;
   329:                   R46 = 1.0;
   330:                   T23 = 1.0;
   331:                   T46 = 1.0;
   332:            
   333: V------>          for (i=1; i<=23; i++)
   334: |                 {
   335: |                   R23 = 0.50 * R23;
   336: |                   T23 = 2.0 * T23;
   337: V------           }
   338: V------>          for (i=1; i<=46; i++)
   339: |                 {
   340: |                   R46 = 0.50 * R46;
   341: |                   T46 = 2.0 * T46;
   342: V------           }
   343:                   KS = 1;
   344:                 }
   345:           
   346:           /*  Break A into two parts such that A = 2^23 * A1 + A2 and set X = N.  */
   347:           
   348:                 T1 = R23 * *A;
   349:                 j  = T1;
   350:                 A1 = j;
   351:                 A2 = *A - T23 * A1;
   352:           
   353:           /*  Break X into two parts such that X = 2^23 * X1 + X2, compute
   354:               Z = A1 * X2 + A2 * X1  (mod 2^23), and then
   355:               X = 2^23 * Z + A2 * X2  (mod 2^46).                            */
   356:           
   357:                 T1 = R23 * *X;
   358:                 j  = T1;
   359:                 X1 = j;
   360:                 X2 = *X - T23 * X1;
   361:                 T1 = A1 * X2 + A2 * X1;
   362:            
   363:                 j  = R23 * T1;
   364:                 T2 = j;
   365:                 Z = T1 - T23 * T2;
   366:                 T3 = T23 * Z + A2 * X2;
   367:                 j  = R46 * T3;
   368:                 T4 = j;
   369:                 *X = T3 - T46 * T4;
   370:                 return(R46 * *X);
   371:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: randlc
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (is.c:333)
    <Vectorized loop.>
    **  Runtime-tested short loop. (is.c:333)
    *** The number of VGT,   VSC.    :  0,  0. (is.c:333)
    *** The number of VLOAD, VSTORE. :  0,  0. (is.c:333)
    *** Idiom detected. : Product (is.c:335)
    *** Idiom detected. : Product (is.c:336)
  LOOP END

  LOOP BEGIN: (is.c:338)
    <Vectorized loop.>
    **  Runtime-tested short loop. (is.c:338)
    *** The number of VGT,   VSC.    :  0,  0. (is.c:338)
    *** The number of VLOAD, VSTORE. :  0,  0. (is.c:338)
    *** Idiom detected. : Product (is.c:340)
    *** Idiom detected. : Product (is.c:341)
  LOOP END


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: randlc
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 31 [s0-s1 s8-s11 s23-s28 s45-s63]
      Vector registers         :  8 [v56-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 16640 bytes
      Register spill area      :  8224 bytes
      Parameter area           :     0 bytes
      Register save area       :   176 bytes
      User data area           :  8224 bytes
      Others                   :     8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: find_my_seed
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   436: par(1801): Parallel routine generated.: create_seq$1
   408: vec( 103): Unvectorized loop.
   408: vec( 110): Vectorization obstructive function reference.: randlc


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: find_my_seed
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

   390:           double   find_my_seed( int kn,        /* my processor rank, 0<=kn<=num procs */
   391:                                  int np,        /* np = num procs                      */
   392:                                  long nn,       /* total num of ran numbers, all procs */
   393:                                  double s,      /* Ran num seed, for ex.: 314159265.00 */
   394:                                  double a )     /* Ran num gen mult, try 1220703125.00 */
   395:           {
   396:           
   397:                 double t1,t2;
   398:                 long   mq,nq,kk,ik;
   399:           
   400:                 if ( kn == 0 ) return s;
   401:           
   402:                 mq = (nn/4 + np - 1) / np;
   403:                 nq = mq * 4 * kn;               /* number of rans to be skipped */
   404:           
   405:                 t1 = s;
   406:                 t2 = a;
   407:                 kk = nq;
   408: +------>        while ( kk > 1 ) {
   409: |               	 ik = kk / 2;
   410: |                  if( 2 * ik ==  kk ) {
   411: |                     (void)randlc( &t2, &t2 );
   412: |         	    kk = ik;
   413: |         	 }
   414: |         	 else {
   415: |                     (void)randlc( &t1, &t2 );
   416: |         	    kk = kk - 1;
   417: |         	 }
   418: +------         }
   419:                 (void)randlc( &t1, &t2 );
   420:           
   421:                 return( t1 );
   422:           
   423:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: find_my_seed
VECTORIZATION LIST

  LOOP BEGIN: (is.c:408)
    <Unvectorized loop.>
    *** Vectorization obstructive function reference. : randlc (is.c:408)
  LOOP END


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: find_my_seed
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 20 [s0-s4 s8-s12 s15-s16 s18-s19 s23-s25 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 208 bytes
      Register spill area      :   0 bytes
      Parameter area           :  16 bytes
      Register save area       : 176 bytes
      User data area           :  16 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (is.c:408)
    *** Estimated execution cycle                       : 22
    *** The number of SCALAR REGISTER TRANSFER          : 6
  LOOP END


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: create_seq
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   459: vec( 103): Unvectorized loop.
   459: vec( 110): Vectorization obstructive function reference.: randlc
   461: opt(1025): Reference to this function inhibits optimization.: randlc


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: create_seq
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

   431:           void	create_seq( double seed, double a )
   432:           {
   433:           	double x, s;
   434:           	INT_TYPE i, k;
   435:           
   436:           #pragma omp parallel private(x,s,i,k)
   437:               {
   438:           	INT_TYPE k1, k2;
   439:           	double an = a;
   440:           	int myid = 0, num_threads = 1;
   441:                   INT_TYPE mq;
   442:           
   443:           #ifdef _OPENMP
   444:           	myid = omp_get_thread_num();
   445:           	num_threads = omp_get_num_threads();
   446:           #endif
   447:           
   448:           	mq = (NUM_KEYS + num_threads - 1) / num_threads;
   449:           	k1 = mq * myid;
   450:           	k2 = k1 + mq;
   451:           	if ( k2 > NUM_KEYS ) k2 = NUM_KEYS;
   452:           
   453:           	KS = 0;
   454:           	s = find_my_seed( myid, num_threads,
   455:           			  (long)4*NUM_KEYS, seed, an );
   456:           
   457:                   k = MAX_KEY/4;
   458:           
   459: +------>  	for (i=k1; i<k2; i++)
   460: |         	{
   461: |         	    x = randlc(&s, &an);
   462: |         	    x += randlc(&s, &an);
   463: |             	    x += randlc(&s, &an);
   464: |         	    x += randlc(&s, &an);
   465: |         
   466: |                     key_array[i] = k*x;
   467: +------   	}
   468:               } /*omp parallel*/
   469:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: create_seq
VECTORIZATION LIST

  LOOP BEGIN: (is.c:459)
    <Unvectorized loop.>
    *** Vectorization obstructive function reference. : randlc (is.c:459)
  LOOP END


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: create_seq
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 22 [s0-s4 s8-s12 s15-s16 s18 s23-s28 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 256 bytes
      Register spill area      :   0 bytes
      Parameter area           :  40 bytes
      Register save area       : 176 bytes
      User data area           :  32 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (is.c:459)
    *** Estimated execution cycle                       : 37
    *** The number of SCALAR REGISTER TRANSFER          : 10
  LOOP END


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: alloc_mem
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: alloc_mem
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

   476:           void *alloc_mem( size_t size )
   477:           {
   478:               void *p;
   479:           
   480:               p = (void *)malloc(size);
   481:               if (!p) {
   482:                   perror("Memory allocation error");
   483:                   exit(1);
   484:               }
   485:               return p;
   486:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: alloc_mem
VECTORIZATION LIST

  No message.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: alloc_mem
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 12 [s0 s8-s12 s15-s16 s23 s61-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 192 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :   0 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: alloc_key_buff
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   501: vec( 103): Unvectorized loop.
   501: vec( 110): Vectorization obstructive function reference.: alloc_mem
   505: par(1801): Parallel routine generated.: alloc_key_buff$1
   502: opt(1025): Reference to this function inhibits optimization.: alloc_mem
   544: par(1801): Parallel routine generated.: full_verify$1
   506: par(1803): Parallelized by "for".
   506: vec( 101): Vectorized loop.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: alloc_key_buff
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

   488:           void alloc_key_buff( void )
   489:           {
   490:               INT_TYPE i;
   491:               int      num_threads = 1;
   492:           
   493:           
   494:           #ifdef _OPENMP
   495:               num_threads = omp_get_max_threads();
   496:           #endif
   497:           
   498:           #ifdef USE_BUCKETS
   499:               bucket_size = (INT_TYPE **)alloc_mem(sizeof(INT_TYPE *) * num_threads);
   500:           
   501: +------>      for (i = 0; i < num_threads; i++) {
   502: |                 bucket_size[i] = (INT_TYPE *)alloc_mem(sizeof(INT_TYPE) * NUM_BUCKETS);
   503: +------       }
   504:           
   505:               #pragma omp parallel for
   506: Y------>      for( i=0; i<NUM_KEYS; i++ )
   507: Y------           key_buff2[i] = 0;
   508:           
   509:           #else /*USE_BUCKETS*/
   510:           
   511:               key_buff1_aptr = (INT_TYPE **)alloc_mem(sizeof(INT_TYPE *) * num_threads);
   512:           
   513:               key_buff1_aptr[0] = key_buff1;
   514:               for (i = 1; i < num_threads; i++) {
   515:                   key_buff1_aptr[i] = (INT_TYPE *)alloc_mem(sizeof(INT_TYPE) * MAX_KEY);
   516:               }
   517:           
   518:           #endif /*USE_BUCKETS*/
   519:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: alloc_key_buff
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (is.c:506)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (is.c:506)
    *** The number of VLOAD, VSTORE. :  0,  1. (is.c:506)
  LOOP END


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: alloc_key_buff
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 18 [s0-s1 s8-s12 s15-s16 s18 s23-s25 s59-s63]
      Vector registers         :  1 [v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 192 bytes
      Register spill area      :   0 bytes
      Parameter area           :   8 bytes
      Register save area       : 176 bytes
      User data area           :   0 bytes
      Others                   :   8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (is.c:506)
    *** Estimated execution cycle                       : 38
  LOOP END


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: full_verify
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   546: par(1803): Parallelized by "for".
   546: opt(3008): Reference within a conditional branch moved outside loop.
   546: vec( 103): Unvectorized loop.
   546: vec( 109): Vectorization obstructive statement.
   546: vec( 113): Overhead of loop division is too large.
   549: vec( 102): Partially vectorized loop.
   549: vec( 113): Overhead of loop division is too large.
   590: par(1801): Parallel routine generated.: full_verify$2
   550: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: *(key_buff_ptr_global)
   591: par(1803): Parallelized by "for".
   591: vec( 101): Vectorized loop.
   593: par(1809): Barrier synchronization.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: full_verify
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

   528:           void full_verify( void )
   529:           {
   530:               INT_TYPE   i, j;
   531:               INT_TYPE   k, k1, k2;
   532:           
   533:           
   534:           /*  Now, finally, sort the keys:  */
   535:           
   536:           /*  Copy keys into work array; keys in key_array will be reassigned. */
   537:           
   538:           #ifdef USE_BUCKETS
   539:           
   540:               /* Buckets are already sorted.  Sorting keys within each bucket */
   541:           #ifdef SCHED_CYCLIC
   542:               #pragma omp parallel for private(i,j,k,k1) schedule(static,1)
   543:           #else
   544:               #pragma omp parallel for private(i,j,k,k1) schedule(dynamic)
   545:           #endif
   546: P------>      for( j=0; j< NUM_BUCKETS; j++ ) {
   547: |         
   548: |                 k1 = (j > 0)? bucket_ptrs[j-1] : 0;
   549: |S----->          for ( i = k1; i < bucket_ptrs[j]; i++ ) {
   550: ||                    k = --key_buff_ptr_global[key_buff2[i]];
   551: ||                    key_array[k] = key_buff2[i];
   552: |S-----           }
   553: P------       }
   554:           
   555:           #else
   556:           
   557:           #pragma omp parallel private(i,j,k,k1,k2)
   558:             {
   559:               #pragma omp for
   560:               for( i=0; i<NUM_KEYS; i++ )
   561:                   key_buff2[i] = key_array[i];
   562:           
   563:               /* This is actual sorting. Each thread is responsible for
   564:                  a subset of key values */
   565:           #ifdef _OPENMP
   566:               j = omp_get_num_threads();
   567:               j = (MAX_KEY + j - 1) / j;
   568:               k1 = j * omp_get_thread_num();
   569:           #else
   570:               j = MAX_KEY;
   571:               k1 = 0;
   572:           #endif
   573:               k2 = k1 + j;
   574:               if (k2 > MAX_KEY) k2 = MAX_KEY;
   575:           
   576:               for( i=0; i<NUM_KEYS; i++ ) {
   577:                   if (key_buff2[i] >= k1 && key_buff2[i] < k2) {
   578:                       k = --key_buff_ptr_global[key_buff2[i]];
   579:                       key_array[k] = key_buff2[i];
   580:                   }
   581:               }
   582:             } /*omp parallel*/
   583:           
   584:           #endif
   585:           
   586:           
   587:           /*  Confirm keys correctly sorted: count incorrectly sorted keys, if any */
   588:           
   589:               j = 0;
   590:               #pragma omp parallel for reduction(+:j)
   591: Y------>      for( i=1; i<NUM_KEYS; i++ )
   592: |                 if( key_array[i-1] > key_array[i] )
   593: Y------               j++;
   594:           
   595:               if( j != 0 )
   596:                   printf( "Full_verify: number of keys out of sort: %ld\n", (long)j );
   597:               else
   598:                   passed_verification++;
   599:           
   600:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: full_verify
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (is.c:591)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (is.c:591)
    *** The number of VLOAD, VSTORE. :  2,  0. (is.c:591)
  LOOP END


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: full_verify
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 21 [s0-s2 s8-s12 s15-s16 s23-s26 s57-s63]
      Vector registers         :  3 [v61-v63]
      Vector mask registers    :  1 [vm15]
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

  LOOP BEGIN: (is.c:591)
    *** Estimated execution cycle                       : 105
    *** The number of SCALAR REGISTER TRANSFER          : 1
  LOOP END


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: rank
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   641: par(1801): Parallel routine generated.: rank$1
   628: opt(1418): Constant-length loop is expanded.
   660: vec( 101): Vectorized loop.
   665: par(1803): Parallelized by "for".
   665: vec( 102): Partially vectorized loop.
   666: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: *(work_buff)
   671: vec( 101): Vectorized loop.
   672: vec( 126): Idiom detected.: Sum
   675: opt(1037): Feedback of array elements.: bucket_ptrs
   676: vec( 101): Vectorized loop.
   677: vec( 126): Idiom detected.: Sum
   678: vec( 101): Vectorized loop.
   679: vec( 126): Idiom detected.: Sum
   685: par(1803): Parallelized by "for".
   685: vec( 102): Partially vectorized loop.
   689: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: bucket_ptrs
   693: opt(1408): Loop interchanged.
   693: opt(1589): Outer loop moved inside inner loop(s).: i
   693: vec( 101): Vectorized loop.
   694: opt(1592): Outer loop unrolled inside inner loop.: k
   709: par(1803): Parallelized by "for".
   709: opt(3008): Reference within a conditional branch moved outside loop.
   709: vec( 109): Vectorization obstructive statement.
   714: vec( 101): Vectorized loop.
   723: vec( 102): Partially vectorized loop.
   724: vec( 122): Dependency unknown. Unvectorizable dependency is assumed.: *(key_buff_ptr)
   731: vec( 101): Vectorized loop.
   732: opt(1037): Feedback of array elements.: key_buff_ptr
   732: vec( 126): Idiom detected.: Iteration
   780: opt(3008): Reference within a conditional branch moved outside loop.
   780: vec( 103): Unvectorized loop.
   780: vec( 110): Vectorization obstructive function reference.: printf
   791: opt(1394): Moved invariant if outside of an inner loop.
   847: opt(1118): This I/O statement inhibits optimization of loop.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: rank
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

   610:           void rank( int iteration )
   611:           {
   612:           
   613:               INT_TYPE    i, k;
   614:               INT_TYPE    *key_buff_ptr, *key_buff_ptr2;
   615:           
   616:           #ifdef USE_BUCKETS
   617:               int shift = MAX_KEY_LOG_2 - NUM_BUCKETS_LOG_2;
   618:               INT_TYPE num_bucket_keys = (1L << shift);
   619:           #endif
   620:           
   621:           
   622:               key_array[iteration] = iteration;
   623:               key_array[iteration+MAX_ITERATIONS] = MAX_KEY - iteration;
   624:           
   625:           
   626:           /*  Determine where the partial verify test keys are, load into  */
   627:           /*  top of array bucket_size                                     */
   628: *------>      for( i=0; i<TEST_ARRAY_SIZE; i++ )
   629: *------           partial_verify_vals[i] = key_array[test_index_array[i]];
   630:           
   631:           
   632:           /*  Setup pointers to key buffers  */
   633:           #ifdef USE_BUCKETS
   634:               key_buff_ptr2 = key_buff2;
   635:           #else
   636:               key_buff_ptr2 = key_array;
   637:           #endif
   638:               key_buff_ptr = key_buff1;
   639:           
   640:           
   641:           #pragma omp parallel private(i, k)
   642:             {
   643:               INT_TYPE *work_buff, m, k1, k2;
   644:               int myid = 0, num_threads = 1;
   645:           
   646:           #ifdef _OPENMP
   647:               myid = omp_get_thread_num();
   648:               num_threads = omp_get_num_threads();
   649:           #endif
   650:           
   651:           
   652:           /*  Bucket sort is known to improve cache performance on some   */
   653:           /*  cache based systems.  But the actual performance may depend */
   654:           /*  on cache size, problem size. */
   655:           #ifdef USE_BUCKETS
   656:           
   657:               work_buff = bucket_size[myid];
   658:           
   659:           /*  Initialize */
   660: V------>      for( i=0; i<NUM_BUCKETS; i++ )
   661: V------           work_buff[i] = 0;
   662:           
   663:           /*  Determine the number of keys in each bucket */
   664:               #pragma omp for schedule(static)
   665: Y------>      for( i=0; i<NUM_KEYS; i++ )
   666: Y------           work_buff[key_array[i] >> shift]++;
   667:           
   668:           /*  Accumulative bucket sizes are the bucket pointers.
   669:               These are global sizes accumulated upon to each bucket */
   670:               bucket_ptrs[0] = 0;
   671: V------>      for( k=0; k< myid; k++ )
   672: V------ G         bucket_ptrs[0] += bucket_size[k][0];
   673:           
   674: +------>      for( i=1; i< NUM_BUCKETS; i++ ) {
   675: |                 bucket_ptrs[i] = bucket_ptrs[i-1];
   676: |V----->          for( k=0; k< myid; k++ )
   677: |V----- G             bucket_ptrs[i] += bucket_size[k][i];
   678: |V----->          for( k=myid; k< num_threads; k++ )
   679: |V----- G             bucket_ptrs[i] += bucket_size[k][i-1];
   680: +------       }
   681:           
   682:           
   683:           /*  Sort into appropriate bucket */
   684:               #pragma omp for schedule(static)
   685: Y------>      for( i=0; i<NUM_KEYS; i++ )
   686: |             {
   687: |                 k = key_array[i];
   688: |                 key_buff2[bucket_ptrs[k >> shift]++] = k;
   689: Y------       }
   690:           
   691:           /*  The bucket pointers now point to the final accumulated sizes */
   692:               if (myid < num_threads-1) {
   693: X------>          for( i=0; i< NUM_BUCKETS; i++ )
   694: |U----->              for( k=myid+1; k< num_threads; k++ )
   695: XU-----                   bucket_ptrs[i] += bucket_size[k][i];
   696:               }
   697:           
   698:           
   699:           /*  Now, buckets are sorted.  We only need to sort keys inside
   700:               each bucket, which can be done in parallel.  Because the distribution
   701:               of the number of keys in the buckets is Gaussian, the use of
   702:               a dynamic schedule should improve load balance, thus, performance     */
   703:           
   704:           #ifdef SCHED_CYCLIC
   705:               #pragma omp for schedule(static,1)
   706:           #else
   707:               #pragma omp for schedule(dynamic)
   708:           #endif
   709: P------>      for( i=0; i< NUM_BUCKETS; i++ ) {
   710: |         
   711: |         /*  Clear the work array section associated with each bucket */
   712: |                 k1 = i * num_bucket_keys;
   713: |                 k2 = k1 + num_bucket_keys;
   714: |V----->          for ( k = k1; k < k2; k++ )
   715: |V-----               key_buff_ptr[k] = 0;
   716: |         
   717: |         /*  Ranking of all keys occurs in this section:                 */
   718: |         
   719: |         /*  In this section, the keys themselves are used as their
   720: |             own indexes to determine how many of each there are: their
   721: |             individual population                                       */
   722: |                 m = (i > 0)? bucket_ptrs[i-1] : 0;
   723: |S----->          for ( k = m; k < bucket_ptrs[i]; k++ )
   724: |S-----               key_buff_ptr[key_buff_ptr2[k]]++;  /* Now they have individual key   */
   725: |                                                /* population                     */
   726: |         
   727: |         /*  To obtain ranks of each key, successively add the individual key
   728: |             population, not forgetting to add m, the total of lesser keys,
   729: |             to the first key population                                          */
   730: |                 key_buff_ptr[k1] += m;
   731: |V----->          for ( k = k1+1; k < k2; k++ )
   732: |V-----               key_buff_ptr[k] += key_buff_ptr[k-1];
   733: |         
   734: P------       }
   735:           
   736:           #else /*USE_BUCKETS*/
   737:           
   738:           
   739:               work_buff = key_buff1_aptr[myid];
   740:           
   741:           
   742:           /*  Clear the work array */
   743:               for( i=0; i<MAX_KEY; i++ )
   744:                   work_buff[i] = 0;
   745:           
   746:           
   747:           /*  Ranking of all keys occurs in this section:                 */
   748:           
   749:           /*  In this section, the keys themselves are used as their
   750:               own indexes to determine how many of each there are: their
   751:               individual population                                       */
   752:           
   753:               #pragma omp for nowait schedule(static)
   754:               for( i=0; i<NUM_KEYS; i++ )
   755:                   work_buff[key_buff_ptr2[i]]++;  /* Now they have individual key   */
   756:                                                  /* population                     */
   757:           
   758:           /*  To obtain ranks of each key, successively add the individual key
   759:               population                                          */
   760:           
   761:               for( i=0; i<MAX_KEY-1; i++ )
   762:                   work_buff[i+1] += work_buff[i];
   763:           
   764:               #pragma omp barrier
   765:           
   766:           /*  Accumulate the global key population */
   767:               for( k=1; k<num_threads; k++ ) {
   768:                   #pragma omp for nowait schedule(static)
   769:                   for( i=0; i<MAX_KEY; i++ )
   770:                       key_buff_ptr[i] += key_buff1_aptr[k][i];
   771:               }
   772:           
   773:           #endif /*USE_BUCKETS*/
   774:           
   775:             } /*omp parallel*/
   776:           
   777:           /* This is the partial verify test section */
   778:           /* Observe that test_rank_array vals are   */
   779:           /* shifted differently for different cases */
   780: +------>      for( i=0; i<TEST_ARRAY_SIZE; i++ )
   781: |             {
   782: |                 k = partial_verify_vals[i];          /* test vals were put here */
   783: |                 if( 0 < k  &&  k <= NUM_KEYS-1 )
   784: |                 {
   785: |                     INT_TYPE key_rank = key_buff_ptr[k-1];
   786: |                     INT_TYPE test_rank = test_rank_array[i];
   787: |                     int failed = 0;
   788: |         
   789: |                     switch( CLASS )
   790: |                     {
   791: |                         case 'S':
   792: |                             if( i <= 2 )
   793: |                                 test_rank += iteration;
   794: |                             else
   795: |                                 test_rank -= iteration;
   796: |                             break;
   797: |                         case 'W':
   798: |                             if( i < 2 )
   799: |                                 test_rank += iteration - 2;
   800: |                             else
   801: |                                 test_rank -= iteration;
   802: |                             break;
   803: |                         case 'A':
   804: |                             if( i <= 2 )
   805: |                                 test_rank += iteration - 1;
   806: |                             else
   807: |                                 test_rank -= iteration - 1;
   808: |                             break;
   809: |                         case 'B':
   810: |                             if( i == 1 || i == 2 || i == 4 )
   811: |                                 test_rank += iteration;
   812: |                             else
   813: |                                 test_rank -= iteration;
   814: |                             break;
   815: |                         case 'C':
   816: |                             if( i <= 2 )
   817: |                                 test_rank += iteration;
   818: |                             else
   819: |                                 test_rank -= iteration;
   820: |                             break;
   821: |                         case 'D':
   822: |                             if( i < 2 )
   823: |                                 test_rank += iteration;
   824: |                             else
   825: |                                 test_rank -= iteration;
   826: |                             break;
   827: |                         case 'E':
   828: |                             if( i < 2 )
   829: |                                 test_rank += iteration - 2;
   830: |                             else if( i == 2 )
   831: |                             {
   832: |                                 test_rank += iteration - 2;
   833: |                                 if (iteration > 4)
   834: |                                     test_rank -= 2;
   835: |                                 else if (iteration > 2)
   836: |                                     test_rank -= 1;
   837: |                             }
   838: |                             else
   839: |                                 test_rank -= iteration - 2;
   840: |                             break;
   841: |                     }
   842: |                     if( key_rank != test_rank )
   843: |                         failed = 1;
   844: |                     else
   845: |                         passed_verification++;
   846: |                     if( failed == 1 )
   847: |                         printf( "Failed partial verification: "
   848: |                                 "iteration %d, test key %d\n",
   849: |                                  iteration, (int)i );
   850: |                 }
   851: +------       }
   852:           
   853:           
   854:           
   855:           
   856:           /*  Make copies of rank info for use by full_verify: these variables
   857:               in rank are local; making them global slows down the code, probably
   858:               since they cannot be made register by compiler                        */
   859:           
   860:               if( iteration == MAX_ITERATIONS )
   861:                   key_buff_ptr_global = key_buff_ptr;
   862:           
   863:           }


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: rank
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (is.c:660)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (is.c:660)
    *** The number of VLOAD, VSTORE. :  0,  1. (is.c:660)
  LOOP END

  LOOP BEGIN: (is.c:665)
    <Partially vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (is.c:665)
    *** The number of VLOAD, VSTORE. :  1,  1. (is.c:665)
    *** Dependency unknown. Unvectorizable dependency is assumed. : *(work_buff) (is.c:666)
  LOOP END

  LOOP BEGIN: (is.c:671)
    <Vectorized loop.>
    *** The number of VGT,   VSC.    :  1,  0. (is.c:671)
    *** The number of VLOAD, VSTORE. :  1,  0. (is.c:671)
    *** VGT generated (is.c:672)
    *** Idiom detected. : Sum (is.c:672)
  LOOP END

  LOOP BEGIN: (is.c:674)
    <Unvectorized loop.>

    LOOP BEGIN: (is.c:676)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  1,  0. (is.c:676)
      *** The number of VLOAD, VSTORE. :  1,  0. (is.c:676)
      *** VGT generated (is.c:677)
      *** Idiom detected. : Sum (is.c:677)
    LOOP END

    LOOP BEGIN: (is.c:678)
      <Vectorized loop.>
      *** The number of VGT,   VSC.    :  1,  0. (is.c:678)
      *** The number of VLOAD, VSTORE. :  1,  0. (is.c:678)
      *** VGT generated (is.c:679)
      *** Idiom detected. : Sum (is.c:679)
    LOOP END
  LOOP END

  LOOP BEGIN: (is.c:685)
    <Partially vectorized loop.>
    *** The number of VGT,   VSC.    :  0,  0. (is.c:685)
    *** The number of VLOAD, VSTORE. :  1,  2. (is.c:685)
    *** Dependency unknown. Unvectorizable dependency is assumed. : bucket_ptrs (is.c:689)
  LOOP END

  LOOP BEGIN: (is.c:694)
    <Unvectorized loop.>
    **  Copied loop. (is.c:694)

    LOOP BEGIN: (is.c:693)
      <Vectorized loop.>
      **  Copied loop. (is.c:693)
      **  Interchanged loop. (is.c:693)
      *** The number of VGT,   VSC.    :  0,  0. (is.c:693)
      *** The number of VLOAD, VSTORE. :  2,  1. (is.c:693)
    LOOP END
  LOOP END

  LOOP BEGIN: (is.c:694)
    <Unvectorized loop.>
    **  Outer loop unrolled loop. : The unroll time is 4. (is.c:694)

    LOOP BEGIN: (is.c:693)
      <Vectorized loop.>
      **  Interchanged loop. (is.c:693)
      *** The number of VGT,   VSC.    :  0,  0. (is.c:693)
      *** The number of VLOAD, VSTORE. :  5,  1. (is.c:693)
    LOOP END
  LOOP END

  LOOP BEGIN: (is.c:709)
    <Unvectorized loop.>
    *** Vectorization obstructive statement. (is.c:709)

    LOOP BEGIN: (is.c:709)
      <Unvectorized loop.>

      LOOP BEGIN: (is.c:714)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (is.c:714)
        *** The number of VLOAD, VSTORE. :  0,  1. (is.c:714)
      LOOP END

      LOOP BEGIN: (is.c:723)
        <Partially vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (is.c:723)
        *** The number of VLOAD, VSTORE. :  1,  1. (is.c:723)
        *** Dependency unknown. Unvectorizable dependency is assumed. : *(key_buff_ptr) (is.c:724)
      LOOP END

      LOOP BEGIN: (is.c:731)
        <Vectorized loop.>
        *** The number of VGT,   VSC.    :  0,  0. (is.c:731)
        *** The number of VLOAD, VSTORE. :  1,  1. (is.c:731)
        *** Idiom detected. : Iteration (is.c:732)
      LOOP END
    LOOP END
  LOOP END


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: rank
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 49 [s0-s4 s8-s12 s15-s16 s18-s33 s43-s63]
      Vector registers         : 50 [v14-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 8608 bytes
      Register spill area      : 4168 bytes
      Parameter area           :   40 bytes
      Register save area       :  176 bytes
      User data area           : 4216 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (is.c:660)
    *** Estimated execution cycle                       : 38
  LOOP END

  LOOP BEGIN: (is.c:665)
    *** Estimated execution cycle                       : 136
  LOOP END

  LOOP BEGIN: (is.c:665)
    *** Estimated execution cycle                       : 9
  LOOP END

  LOOP BEGIN: (is.c:671)
    *** Estimated execution cycle                       : 158
  LOOP END

  LOOP BEGIN: (is.c:674)
    *** Estimated execution cycle                       : 130
    *** The number of SCALAR REGISTER TRANSFER          : 18

    LOOP BEGIN: (is.c:676)
      *** Estimated execution cycle                     : 158
    LOOP END

    LOOP BEGIN: (is.c:678)
      *** Estimated execution cycle                     : 158
    LOOP END
  LOOP END

  LOOP BEGIN: (is.c:685)
    *** Estimated execution cycle                       : 136
  LOOP END

  LOOP BEGIN: (is.c:685)
    *** Estimated execution cycle                       : 13
  LOOP END

  LOOP BEGIN: (is.c:694)
    *** Estimated execution cycle                       : 7

    LOOP BEGIN: (is.c:693)
      *** Estimated execution cycle                     : 106
    LOOP END
  LOOP END

  LOOP BEGIN: (is.c:694)
    *** Estimated execution cycle                       : 13

    LOOP BEGIN: (is.c:693)
      *** Estimated execution cycle                     : 299
    LOOP END
  LOOP END

  LOOP BEGIN: (is.c:709)
    *** Estimated execution cycle                       : 27
    *** The number of SCALAR REGISTER SPILL
          Total                                         : 4
            Over basic blocks                           : 4
    *** The number of SCALAR REGISTER RESTORE
          Total                                         : 2
            Over basic blocks                           : 2
    *** The number of SCALAR REGISTER TRANSFER          : 9

    LOOP BEGIN: (is.c:709)
      *** Estimated execution cycle                     : 29
      *** The number of SCALAR REGISTER SPILL
            Total                                       : 1
              Across calls                              : 1
      *** The number of SCALAR REGISTER RESTORE
            Total                                       : 5
              Across calls                              : 1
              Over basic blocks                         : 4
      *** The number of SCALAR REGISTER TRANSFER        : 61

      LOOP BEGIN: (is.c:714)
        *** Estimated execution cycle                   : 38
      LOOP END

      LOOP BEGIN: (is.c:723)
        *** Estimated execution cycle                   : 104
      LOOP END

      LOOP BEGIN: (is.c:723)
        *** Estimated execution cycle                   : 9
      LOOP END

      LOOP BEGIN: (is.c:731)
        *** Estimated execution cycle                   : 379
        *** The number of SCALAR REGISTER TRANSFER      : 6
      LOOP END
    LOOP END
  LOOP END


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: main
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   892: vec( 101): Vectorized loop.
   895: opt(1394): Moved invariant if outside of an inner loop.
   961: vec( 103): Unvectorized loop.
   961: vec( 110): Vectorization obstructive function reference.: printf
   961: vec( 110): Vectorization obstructive function reference.: rank
   963: opt(1118): This I/O statement inhibits optimization of loop.
   963: opt(1317): Logical clause simplified.


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: main
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

   870:           int main( int argc, char **argv )
   871:           {
   872:           
   873:               int             i, iteration, timer_on;
   874:           
   875:               double          timecounter;
   876:           
   877:           
   878:           /*  Initialize timers  */
   879:               timer_on = check_timer_flag();
   880:           
   881:               timer_clear( 0 );
   882:               if (timer_on) {
   883:                   timer_clear( 1 );
   884:                   timer_clear( 2 );
   885:                   timer_clear( 3 );
   886:               }
   887:           
   888:               if (timer_on) timer_start( 3 );
   889:           
   890:           
   891:           /*  Initialize the verification arrays if a valid class */
   892: V------>      for( i=0; i<TEST_ARRAY_SIZE; i++ )
   893: |                 switch( CLASS )
   894: |                 {
   895: |                     case 'S':
   896: |                         test_index_array[i] = S_test_index_array[i];
   897: |                         test_rank_array[i]  = S_test_rank_array[i];
   898: |                         break;
   899: |                     case 'A':
   900: |                         test_index_array[i] = A_test_index_array[i];
   901: |                         test_rank_array[i]  = A_test_rank_array[i];
   902: |                         break;
   903: |                     case 'W':
   904: |                         test_index_array[i] = W_test_index_array[i];
   905: |                         test_rank_array[i]  = W_test_rank_array[i];
   906: |                         break;
   907: |                     case 'B':
   908: |                         test_index_array[i] = B_test_index_array[i];
   909: |                         test_rank_array[i]  = B_test_rank_array[i];
   910: |                         break;
   911: |                     case 'C':
   912: |                         test_index_array[i] = C_test_index_array[i];
   913: |                         test_rank_array[i]  = C_test_rank_array[i];
   914: |                         break;
   915: |                     case 'D':
   916: |                         test_index_array[i] = D_test_index_array[i];
   917: |                         test_rank_array[i]  = D_test_rank_array[i];
   918: |                         break;
   919: |                     case 'E':
   920: |                         test_index_array[i] = E_test_index_array[i];
   921: |                         test_rank_array[i]  = E_test_rank_array[i];
   922: |                         break;
   923: V------           };
   924:           
   925:            
   926:           
   927:           /*  Printout initial NPB info */
   928:               printf
   929:                 ( "\n\n NAS Parallel Benchmarks (NPB3.4-OMP) - IS Benchmark\n\n" );
   930:               printf( " Size:  %ld  (class %c)\n", (long)TOTAL_KEYS, CLASS );
   931:               printf( " Iterations:  %d\n", MAX_ITERATIONS );
   932:           #ifdef _OPENMP
   933:               printf( " Number of available threads:  %d\n", omp_get_max_threads() );
   934:           #endif
   935:               printf( "\n" );
   936:           
   937:               if (timer_on) timer_start( 1 );
   938:           
   939:           /*  Generate random number sequence and subsequent keys on all procs */
   940:               create_seq( 314159265.00,                    /* Random number gen seed */
   941:                           1220703125.00 );                 /* Random number gen mult */
   942:           
   943:               alloc_key_buff();
   944:               if (timer_on) timer_stop( 1 );
   945:           
   946:           
   947:           /*  Do one interation for free (i.e., untimed) to guarantee initialization of
   948:               all data and code pages and respective tables */
   949:               rank( 1 );
   950:           
   951:           /*  Start verification counter */
   952:               passed_verification = 0;
   953:           
   954:               if( CLASS != 'S' ) printf( "\n   iteration\n" );
   955:           
   956:           /*  Start timer  */
   957:               timer_start( 0 );
   958:           
   959:           
   960:           /*  This is the main iteration */
   961: +------>      for( iteration=1; iteration<=MAX_ITERATIONS; iteration++ )
   962: |             {
   963: |                 if( CLASS != 'S' ) printf( "        %d\n", iteration );
   964: |                 rank( iteration );
   965: +------       }
   966:           
   967:           
   968:           /*  End of timing, obtain maximum time of all processors */
   969:               timer_stop( 0 );
   970:               timecounter = timer_read( 0 );
   971:           
   972:           
   973:           /*  This tests that keys are in sequence: sorting of last ranked key seq
   974:               occurs here, but is an untimed operation                             */
   975:               if (timer_on) timer_start( 2 );
   976:               full_verify();
   977:               if (timer_on) timer_stop( 2 );
   978:           
   979:               if (timer_on) timer_stop( 3 );
   980:           
   981:           
   982:           /*  The final printout  */
   983:               if( passed_verification != 5*MAX_ITERATIONS + 1 )
   984:                   passed_verification = 0;
   985:               c_print_results( "IS",
   986:                                CLASS,
   987:                                TOTAL_KS1,
   988:                                TOTAL_KS2,
   989:                                0,
   990:                                MAX_ITERATIONS,
   991:                                timecounter,
   992:                                1.0e-6*(double)(TOTAL_KEYS)*MAX_ITERATIONS
   993:                                                             /timecounter,
   994:                                "keys ranked",
   995:                                passed_verification,
   996:                                NPBVERSION,
   997:                                COMPILETIME,
   998:                                CC,
   999:                                CLINK,
  1000:                                C_LIB,
  1001:                                C_INC,
  1002:                                CFLAGS,
  1003:                                CLINKFLAGS );
  1004:           
  1005:           
  1006:           /*  Print additional timers  */
  1007:               if (timer_on) {
  1008:                  double t_total, t_percent;
  1009:           
  1010:                  t_total = timer_read( 3 );
  1011:                  printf("\nSECTION   Time (secs)\n");
  1012:                  printf(" Total : %8.3f\n", t_total);
  1013:                  if (t_total == 0.0) t_total = 1.0;
  1014:                  timecounter = timer_read(1);
  1015:                  t_percent = timecounter/t_total * 100.;
  1016:                  printf(" Initialization : %8.3f ( %5.2f %%)\n", timecounter, t_percent);
  1017:                  timecounter = timer_read(0);
  1018:                  t_percent = timecounter/t_total * 100.;
  1019:                  printf(" Benchmarking   : %8.3f ( %5.2f %%)\n", timecounter, t_percent);
  1020:                  timecounter = timer_read(2);
  1021:                  t_percent = timecounter/t_total * 100.;
  1022:                  printf(" Sorting        : %8.3f ( %5.2f %%)\n", timecounter, t_percent);
  1023:               }
  1024:           
  1025:               return 0;
  1026:                    /**************************/
  1027:           }        /*  E N D  P R O G R A M  */


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: main
VECTORIZATION LIST

  NOTE: The number of operation does not include the instruction
        for register spill and restore.

  LOOP BEGIN: (is.c:892)
    <Vectorized loop.>
    **  Runtime-tested short loop. (is.c:892)
    *** The number of VGT,   VSC.    :  0,  0. (is.c:892)
    *** The number of VLOAD, VSTORE. :  2,  2. (is.c:892)
  LOOP END

  LOOP BEGIN: (is.c:961)
    <Unvectorized loop.>
    *** Vectorization obstructive function reference. : printf (is.c:961)
    *** Vectorization obstructive function reference. : rank (is.c:961)
  LOOP END


NEC C/C++ Compiler (3.2.0) for Vector Engine     Thu Apr 15 12:28:59 2021
FILE NAME: is.c

FUNCTION NAME: main
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 27 [s0-s12 s15-s16 s18 s23-s30 s61-s63]
      Vector registers         :  2 [v62-v63]
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 320 bytes
      Register spill area      :   0 bytes
      Parameter area           : 144 bytes
      Register save area       : 176 bytes
      User data area           :   0 bytes
      Others                   :   0 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.

  LOOP BEGIN: (is.c:961)
    *** Estimated execution cycle                       : 16
    *** The number of SCALAR REGISTER TRANSFER          : 3
  LOOP END


