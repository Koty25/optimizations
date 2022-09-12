NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:58 2021
FILE NAME : ua_data.f

  COMPILER OPTIONS : -c -O4 -fopenmp -report-all -fdiag-inline=2 -fdiag-parallel=2 -fdiag-vector=2 -finline-functions -fouterloop-unroll

  PARAMETER : 

      Optimization Options : 
        -On                                      : 4
        -fargument-alias                         : disable
        -fargument-noalias                       : enable
        -fassociative-math                       : enable
        -faggressive-associative-math            : disable
        -fassume-contiguous                      : disable
        -fcopyin-intent-out                      : enable
        -fcse-after-vectorization                : disable
        -ffast-formatted-io                      : enable
        -ffast-math                              : enable
        -fignore-asynchronous                    : disable
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
        -floop-unroll-complete-nest              : 3
        -floop-unroll-max-times                  : auto
        -fmatrix-multiply                        : enable
        -fmove-loop-invariants                   : enable
        -fmove-loop-invariants-if                : enable
        -fmove-loop-invariants-unsafe            : enable
        -fmove-nested-loop-invariants-outer      : enable
        -fnamed-alias                            : disable
        -fnamed-noalias                          : enable
        -fnamed-noalias-aggressive               : disable
        -fouterloop-unroll                       : enable
        -fouterloop-unroll-max-size              : 4
        -fouterloop-unroll-max-times             : 4
        -fprecise-math                           : disable
        -freciprocal-math                        : enable
        -freorder-logical-expression             : enable
        -freplace-loop-equation                  : enable
        -marray-io                               : enable
        -mlist-vector                            : disable
        -mretain-all                             : enable
        -mretain-list-vector                     : disable
        -mretain-none                            : disable
        -msched-none                             : disable
        -msched-insns                            : enable
        -msched-block                            : disable
        -mstack-arrays                           : enable
        -mno-stack-arrays                        : disable
        -mvector                                 : enable
        -mvector-advance-gather                  : enable
        -mvector-advance-gather-limit            : 56
        -mvector-dependency-test                 : enable
        -mvector-floating-divide-instruction     : disable
        -mvector-fma                             : enable
        -mvector-intrinsic-check                 : disable
        -mvector-iteration                       : enable
        -mvector-iteration-unsafe                : enable
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
        -fbounds-check                           : disable
        -fcheck=all                              : disable
        -fcheck=alias                            : disable
        -fcheck=bits                             : disable
        -fcheck=bounds                           : disable
        -fcheck=dangling                         : disable
        -fcheck=do                               : disable
        -fcheck=iovf                             : disable
        -fcheck=pointer                          : disable
        -fcheck=present                          : disable
        -fcheck=recursion                        : disable
        -g                                       : disable
        -minit-stack                             : not specified
        -mmemory-trace                           : disable
        -mmemory-trace-full                      : disable
        -traceback                               : disable

      Language Options : 
        -bss                                     : disable
        -fdefault-integer                        : 4
        -fdefault-double                         : 8
        -fdefault-real                           : 4
        -fextend-source                          : disable
        -ffree-form                              : disable
        -ffixed-form                             : enable
        -ff90-sign                               : disable
        -fmax-continuation-lines                 : 1023
        -frealloc-lhs                            : enable
        -frealloc-lhs-array                      : enable
        -frealloc-lhs-scalar                     : enable
        -save                                    : disable
        -std                                     : Fortran 2008

      Message Options : 
        -Wall                                    : disable
        -Werror                                  : disable
        -Wextension                              : enable
        -Wobsolescent                            : enable
        -Woverflow                               : disable
        -Woverflow-errors                        : disable
        -fdiag-inline                            : 2
        -fdiag-parallel                          : 2
        -fdiag-vector                            : 2
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

      Directory Options : 
        -fintrinsic-modules-path                 : /opt/nec/ve/nfort/3.2.0/include
        -module                                  : not specified

      Miscellaneous Options : 
        -v                                       : disable


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:58 2021
FILE NAME: ua_data.f

PROCEDURE NAME: ALLOC_SPACE
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:58 2021
FILE NAME: ua_data.f

PROCEDURE NAME: ALLOC_SPACE
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

   295:                 subroutine alloc_space
   296:           
   297:           c---------------------------------------------------------------------
   298:           c---------------------------------------------------------------------
   299:           
   300:           c---------------------------------------------------------------------
   301:           c allocate space dynamically for data arrays
   302:           c---------------------------------------------------------------------
   303:           
   304:                 use ua_data
   305:                 implicit none
   306:           
   307:                 integer ios
   308:           
   309:           
   310:                 allocate (
   311:                &        ta1  (lx1,lx1,lx1,lelt), ta2   (lx1,lx1,lx1,lelt),
   312:                &        trhs (lx1,lx1,lx1,lelt), t     (lx1,lx1,lx1,lelt),
   313:                &        tmult(lx1,lx1,lx1,lelt), dpcelm(lx1,lx1,lx1,lelt),
   314:                &        pdiff(lx1,lx1,lx1,lelt), pdiffp(lx1,lx1,lx1,lelt),
   315:                >        stat = ios)
   316:           
   317:                 if (ios .eq. 0) allocate (
   318:                &        umor(lmor), tmmor(lmor),
   319:                &        rmor(lmor), dpcmor (lmor), pmorx(lmor), ppmor(lmor),
   320:                &        mormult(lmor), tmort(lmor),
   321:                >        stat = ios)
   322:           
   323:                 if (ios .eq. 0) allocate (
   324:                &        idmo    (lx1,lx1,lnje,lnje,nsides,lelt),
   325:                &        idel    (lx1,lx1,          nsides,lelt),
   326:                &        sje     (2,2,              nsides,lelt),
   327:                &        sje_new (2,2,              nsides,lelt),
   328:                &        ijel    (2,                nsides,lelt),
   329:                &        ijel_new(2,                nsides,lelt),
   330:                &        cbc     (                  nsides,lelt),
   331:                &        cbc_new (                  nsides,lelt),
   332:                &        vassign (8,lelt),       emo(2,8,8*lelt),
   333:                &        nemo    (8*lelt),
   334:                &        diagn   (2,12,lelt),
   335:                >        stat = ios)
   336:           
   337:                 if (ios .eq. 0) allocate (
   338:                &        tree   (lelt), mt_to_id    (lelt),
   339:                &        newc   (lelt), mt_to_id_old(lelt),
   340:                &        newi   (lelt), id_to_mt    (lelt),
   341:                &        newe   (lelt), ref_front_id(lelt),
   342:                &        front  (lelt), action      (lelt),
   343:                &        ich    (lelt), size_e      (lelt),
   344:                &        treenew(lelt),
   345:                >        stat = ios)
   346:           
   347:                 if (ios .eq. 0) allocate (
   348:                &        ifpcmor  (8* lelt),
   349:                &        eassign  (12,lelt),  if_1_edge(12,lelt),
   350:                &        ncon_edge(12,lelt),
   351:                &        skip (lelt), ifcoa (lelt), ifcoa_id(lelt),
   352:                &        fassign(nsides,lelt), edgevis(4,nsides,lelt),
   353:                >        stat = ios)
   354:           
   355:           c.....coordinate of element vertices
   356:                 if (ios .eq. 0) allocate (
   357:                &        xc    (8,lelt),yc    (8,lelt),zc    (8,lelt),
   358:                &        xc_new(8,lelt),yc_new(8,lelt),zc_new(8,lelt),
   359:                >        stat = ios)
   360:           
   361:                 if (ios .ne. 0) then
   362:                    write(*,*) 'Error encountered in allocating space'
   363:                    stop
   364:                 endif
   365:           
   366:                 return
   367:                 end


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:58 2021
FILE NAME: ua_data.f

PROCEDURE NAME: ALLOC_SPACE
INLINE LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:58 2021
FILE NAME: ua_data.f

PROCEDURE NAME: ALLOC_SPACE
VECTORIZATION LIST

  No message.


NEC Fortran Compiler (3.2.0) for Vector Engine   Fri May 21 01:36:58 2021
FILE NAME: ua_data.f

PROCEDURE NAME: ALLOC_SPACE
CODE GENERATION LIST

  Hardware registers
    Reserved                   : 10 [sl fp lr sp s12 s13 tp got plt s17]
    Callee-saved               : 16 [s18-s33]
    Assigned
      Scalar registers         : 22 [s0-s12 s15-s16 s18 s23-s24 s60-s63]
      Vector registers         :  0
      Vector mask registers    :  0
      VREG directive           :  0

  Routine stack
    Total size                 : 1312 bytes
      Register spill area      :    0 bytes
      Parameter area           :   72 bytes
      Register save area       :  176 bytes
      User data area           : 1056 bytes
      Others                   :    8 bytes

    Note: Total size of Routine stack does not include
          the size extended by alloca() and so on.


