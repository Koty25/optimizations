NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:38 2020
FILE NAME: lu.f

PROCEDURE NAME: APPLU
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

   171: vec( 103): Unvectorized loop.
   171: vec( 110): Vectorization obstructive procedure reference.: TIMER_READ
   179: vec( 103): Unvectorized loop.
   179: vec( 180): I/O statement obstructs vectorization.
   181: opt(1118): This I/O statement inhibits optimization of loop.


NEC Fortran Compiler (3.0.4) for Vector Engine   Sat Jul 18 20:59:38 2020
FILE NAME: lu.f

PROCEDURE NAME: APPLU
FORMAT LIST

 LINE   LOOP      STATEMENT

    49:                 program applu
    50:           c---------------------------------------------------------------------
    51:           
    52:           c---------------------------------------------------------------------
    53:           c
    54:           c   driver for the performance evaluation of the solver for
    55:           c   five coupled parabolic/elliptic partial differential equations.
    56:           c
    57:           c---------------------------------------------------------------------
    58:           
    59:                 use lu_data
    60:                 implicit none
    61:           
    62:                 character class
    63:                 logical verified
    64:                 double precision mflops
    65:           
    66:                 double precision t, tmax, timer_read, trecs(t_last)
    67:                 external timer_read
    68:                 integer i
    69:                 character t_names(t_last)*8
    70:           
    71:           c---------------------------------------------------------------------
    72:           c     Setup info for timers
    73:           c---------------------------------------------------------------------
    74:           
    75:                 call check_timer_flag( timeron )
    76:                 if (timeron) then
    77:                    t_names(t_total) = 'total'
    78:                    t_names(t_rhsx) = 'rhsx'
    79:                    t_names(t_rhsy) = 'rhsy'
    80:                    t_names(t_rhsz) = 'rhsz'
    81:                    t_names(t_rhs) = 'rhs'
    82:                    t_names(t_jacld) = 'jacld'
    83:                    t_names(t_blts) = 'blts'
    84:                    t_names(t_jacu) = 'jacu'
    85:                    t_names(t_buts) = 'buts'
    86:                    t_names(t_add) = 'add'
    87:                    t_names(t_l2norm) = 'l2norm'
    88:                 endif
    89:           
    90:           c---------------------------------------------------------------------
    91:           c   read input data
    92:           c---------------------------------------------------------------------
    93:                 call read_input()
    94:           
    95:           c---------------------------------------------------------------------
    96:           c   set up domain sizes
    97:           c---------------------------------------------------------------------
    98:                 call domain()
    99:           
   100:                 call alloc_space
   101:           
   102:           c---------------------------------------------------------------------
   103:           c   set up coefficients
   104:           c---------------------------------------------------------------------
   105:                 call setcoeff()
   106:           
   107:           c---------------------------------------------------------------------
   108:           c   set the boundary values for dependent variables
   109:           c---------------------------------------------------------------------
   110:                 call setbv()
   111:           
   112:           c---------------------------------------------------------------------
   113:           c   set the initial values for dependent variables
   114:           c---------------------------------------------------------------------
   115:                 call setiv()
   116:           
   117:           c---------------------------------------------------------------------
   118:           c   compute the forcing term based on prescribed exact solution
   119:           c---------------------------------------------------------------------
   120:                 call erhs()
   121:           
   122:           c---------------------------------------------------------------------
   123:           c   perform one SSOR iteration to touch all data pages
   124:           c---------------------------------------------------------------------
   125:                 call ssor(1)
   126:           
   127:           c---------------------------------------------------------------------
   128:           c   reset the boundary and initial values
   129:           c---------------------------------------------------------------------
   130:                 call setbv()
   131:                 call setiv()
   132:           
   133:           c---------------------------------------------------------------------
   134:           c   perform the SSOR iterations
   135:           c---------------------------------------------------------------------
   136:                 call ssor(itmax)
   137:           
   138:           c---------------------------------------------------------------------
   139:           c   compute the solution error
   140:           c---------------------------------------------------------------------
   141:                 call error()
   142:           
   143:           c---------------------------------------------------------------------
   144:           c   compute the surface integral
   145:           c---------------------------------------------------------------------
   146:                 call pintgr()
   147:           
   148:           c---------------------------------------------------------------------
   149:           c   verification test
   150:           c---------------------------------------------------------------------
   151:                 call verify ( rsdnm, errnm, frc, class, verified )
   152:                 mflops = 1.0d-6*dble(itmax)*(1984.77*dble( nx0 )
   153:                >     *dble( ny0 )
   154:                >     *dble( nz0 )
   155:                >     -10923.3*(dble( nx0+ny0+nz0 )/3.)**2
   156:                >     +27770.9* dble( nx0+ny0+nz0 )/3.
   157:                >     -144010.)
   158:                >     / maxtime
   159:           
   160:                 call print_results('LU', class, nx0,
   161:                >  ny0, nz0, itmax,
   162:                >  maxtime, mflops, '          floating point', verified,
   163:                >  npbversion, compiletime, cs1, cs2, cs3, cs4, cs5, cs6,
   164:                >  '(none)')
   165:           
   166:           c---------------------------------------------------------------------
   167:           c      More timers
   168:           c---------------------------------------------------------------------
   169:                 if (.not.timeron) goto 999
   170:           
   171: +------>        do i=1, t_last
   172: |                  trecs(i) = timer_read(i)
   173: +------         end do
   174:                 tmax = maxtime
   175:                 if ( tmax .eq. 0. ) tmax = 1.0
   176:           
   177:                 write(*,800)
   178:            800  format('  SECTION     Time (secs)')
   179: +------>        do i=1, t_last
   180: |                  if (i.ne.t_jacld .and. i.ne.t_jacu) then
   181: |                     write(*,810) t_names(i), trecs(i), trecs(i)*100./tmax
   182: |                  endif
   183: |                  if (i.eq.t_rhs) then
   184: |                     t = trecs(t_rhsx) + trecs(t_rhsy) + trecs(t_rhsz)
   185: |                     write(*,820) 'sub-rhs', t, t*100./tmax
   186: |                     t = trecs(i) - t
   187: |                     write(*,820) 'rest-rhs', t, t*100./tmax
   188: |                  endif
   189: |          810     format(2x,a8,':',f9.3,'  (',f6.2,'%)')
   190: |          820     format(5x,'--> ',a8,':',f9.3,'  (',f6.2,'%)')
   191: +------         end do
   192:           
   193:            999  continue
   194:                 end


