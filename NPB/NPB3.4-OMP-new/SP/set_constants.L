NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:56 2021
FILE NAME: set_constants.f

PROCEDURE NAME: SET_CONSTANTS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Feb 11 10:50:56 2021
FILE NAME: set_constants.f

PROCEDURE NAME: SET_CONSTANTS
FORMAT LIST

 LINE   LOOP      STATEMENT

     5:                  subroutine  set_constants
     6:           
     7:           c---------------------------------------------------------------------
     8:           c---------------------------------------------------------------------
     9:           
    10:                  use sp_data
    11:                  implicit none
    12:            
    13:                  ce(1,1)  = 2.0d0
    14:                  ce(1,2)  = 0.0d0
    15:                  ce(1,3)  = 0.0d0
    16:                  ce(1,4)  = 4.0d0
    17:                  ce(1,5)  = 5.0d0
    18:                  ce(1,6)  = 3.0d0
    19:                  ce(1,7)  = 0.5d0
    20:                  ce(1,8)  = 0.02d0
    21:                  ce(1,9)  = 0.01d0
    22:                  ce(1,10) = 0.03d0
    23:                  ce(1,11) = 0.5d0
    24:                  ce(1,12) = 0.4d0
    25:                  ce(1,13) = 0.3d0
    26:            
    27:                  ce(2,1)  = 1.0d0
    28:                  ce(2,2)  = 0.0d0
    29:                  ce(2,3)  = 0.0d0
    30:                  ce(2,4)  = 0.0d0
    31:                  ce(2,5)  = 1.0d0
    32:                  ce(2,6)  = 2.0d0
    33:                  ce(2,7)  = 3.0d0
    34:                  ce(2,8)  = 0.01d0
    35:                  ce(2,9)  = 0.03d0
    36:                  ce(2,10) = 0.02d0
    37:                  ce(2,11) = 0.4d0
    38:                  ce(2,12) = 0.3d0
    39:                  ce(2,13) = 0.5d0
    40:           
    41:                  ce(3,1)  = 2.0d0
    42:                  ce(3,2)  = 2.0d0
    43:                  ce(3,3)  = 0.0d0
    44:                  ce(3,4)  = 0.0d0
    45:                  ce(3,5)  = 0.0d0
    46:                  ce(3,6)  = 2.0d0
    47:                  ce(3,7)  = 3.0d0
    48:                  ce(3,8)  = 0.04d0
    49:                  ce(3,9)  = 0.03d0
    50:                  ce(3,10) = 0.05d0
    51:                  ce(3,11) = 0.3d0
    52:                  ce(3,12) = 0.5d0
    53:                  ce(3,13) = 0.4d0
    54:           
    55:                  ce(4,1)  = 2.0d0
    56:                  ce(4,2)  = 2.0d0
    57:                  ce(4,3)  = 0.0d0
    58:                  ce(4,4)  = 0.0d0
    59:                  ce(4,5)  = 0.0d0
    60:                  ce(4,6)  = 2.0d0
    61:                  ce(4,7)  = 3.0d0
    62:                  ce(4,8)  = 0.03d0
    63:                  ce(4,9)  = 0.05d0
    64:                  ce(4,10) = 0.04d0
    65:                  ce(4,11) = 0.2d0
    66:                  ce(4,12) = 0.1d0
    67:                  ce(4,13) = 0.3d0
    68:           
    69:                  ce(5,1)  = 5.0d0
    70:                  ce(5,2)  = 4.0d0
    71:                  ce(5,3)  = 3.0d0
    72:                  ce(5,4)  = 2.0d0
    73:                  ce(5,5)  = 0.1d0
    74:                  ce(5,6)  = 0.4d0
    75:                  ce(5,7)  = 0.3d0
    76:                  ce(5,8)  = 0.05d0
    77:                  ce(5,9)  = 0.04d0
    78:                  ce(5,10) = 0.03d0
    79:                  ce(5,11) = 0.1d0
    80:                  ce(5,12) = 0.3d0
    81:                  ce(5,13) = 0.2d0
    82:           
    83:                  c1 = 1.4d0
    84:                  c2 = 0.4d0
    85:                  c3 = 0.1d0
    86:                  c4 = 1.0d0
    87:                  c5 = 1.4d0
    88:           
    89:                  bt = dsqrt(0.5d0)
    90:           
    91:                  dnxm1 = 1.0d0 / dble(grid_points(1)-1)
    92:                  dnym1 = 1.0d0 / dble(grid_points(2)-1)
    93:                  dnzm1 = 1.0d0 / dble(grid_points(3)-1)
    94:           
    95:                  c1c2 = c1 * c2
    96:                  c1c5 = c1 * c5
    97:                  c3c4 = c3 * c4
    98:                  c1345 = c1c5 * c3c4
    99:           
   100:                  conz1 = (1.0d0-c1c5)
   101:           
   102:                  tx1 = 1.0d0 / (dnxm1 * dnxm1)
   103:                  tx2 = 1.0d0 / (2.0d0 * dnxm1)
   104:                  tx3 = 1.0d0 / dnxm1
   105:           
   106:                  ty1 = 1.0d0 / (dnym1 * dnym1)
   107:                  ty2 = 1.0d0 / (2.0d0 * dnym1)
   108:                  ty3 = 1.0d0 / dnym1
   109:            
   110:                  tz1 = 1.0d0 / (dnzm1 * dnzm1)
   111:                  tz2 = 1.0d0 / (2.0d0 * dnzm1)
   112:                  tz3 = 1.0d0 / dnzm1
   113:           
   114:                  dx1 = 0.75d0
   115:                  dx2 = 0.75d0
   116:                  dx3 = 0.75d0
   117:                  dx4 = 0.75d0
   118:                  dx5 = 0.75d0
   119:           
   120:                  dy1 = 0.75d0
   121:                  dy2 = 0.75d0
   122:                  dy3 = 0.75d0
   123:                  dy4 = 0.75d0
   124:                  dy5 = 0.75d0
   125:           
   126:                  dz1 = 1.0d0
   127:                  dz2 = 1.0d0
   128:                  dz3 = 1.0d0
   129:                  dz4 = 1.0d0
   130:                  dz5 = 1.0d0
   131:           
   132:                  dxmax = dmax1(dx3, dx4)
   133:                  dymax = dmax1(dy2, dy4)
   134:                  dzmax = dmax1(dz2, dz3)
   135:           
   136:                  dssp = 0.25d0 * dmax1(dx1, dmax1(dy1, dz1) )
   137:           
   138:                  c4dssp = 4.0d0 * dssp
   139:                  c5dssp = 5.0d0 * dssp
   140:           
   141:                  dttx1 = dt*tx1
   142:                  dttx2 = dt*tx2
   143:                  dtty1 = dt*ty1
   144:                  dtty2 = dt*ty2
   145:                  dttz1 = dt*tz1
   146:                  dttz2 = dt*tz2
   147:           
   148:                  c2dttx1 = 2.0d0*dttx1
   149:                  c2dtty1 = 2.0d0*dtty1
   150:                  c2dttz1 = 2.0d0*dttz1
   151:           
   152:                  dtdssp = dt*dssp
   153:           
   154:                  comz1  = dtdssp
   155:                  comz4  = 4.0d0*dtdssp
   156:                  comz5  = 5.0d0*dtdssp
   157:                  comz6  = 6.0d0*dtdssp
   158:           
   159:                  c3c4tx3 = c3c4*tx3
   160:                  c3c4ty3 = c3c4*ty3
   161:                  c3c4tz3 = c3c4*tz3
   162:           
   163:                  dx1tx1 = dx1*tx1
   164:                  dx2tx1 = dx2*tx1
   165:                  dx3tx1 = dx3*tx1
   166:                  dx4tx1 = dx4*tx1
   167:                  dx5tx1 = dx5*tx1
   168:            
   169:                  dy1ty1 = dy1*ty1
   170:                  dy2ty1 = dy2*ty1
   171:                  dy3ty1 = dy3*ty1
   172:                  dy4ty1 = dy4*ty1
   173:                  dy5ty1 = dy5*ty1
   174:            
   175:                  dz1tz1 = dz1*tz1
   176:                  dz2tz1 = dz2*tz1
   177:                  dz3tz1 = dz3*tz1
   178:                  dz4tz1 = dz4*tz1
   179:                  dz5tz1 = dz5*tz1
   180:           
   181:                  c2iv  = 2.5d0
   182:                  con43 = 4.0d0/3.0d0
   183:                  con16 = 1.0d0/6.0d0
   184:            
   185:                  xxcon1 = c3c4tx3*con43*tx3
   186:                  xxcon2 = c3c4tx3*tx3
   187:                  xxcon3 = c3c4tx3*conz1*tx3
   188:                  xxcon4 = c3c4tx3*con16*tx3
   189:                  xxcon5 = c3c4tx3*c1c5*tx3
   190:           
   191:                  yycon1 = c3c4ty3*con43*ty3
   192:                  yycon2 = c3c4ty3*ty3
   193:                  yycon3 = c3c4ty3*conz1*ty3
   194:                  yycon4 = c3c4ty3*con16*ty3
   195:                  yycon5 = c3c4ty3*c1c5*ty3
   196:           
   197:                  zzcon1 = c3c4tz3*con43*tz3
   198:                  zzcon2 = c3c4tz3*tz3
   199:                  zzcon3 = c3c4tz3*conz1*tz3
   200:                  zzcon4 = c3c4tz3*con16*tz3
   201:                  zzcon5 = c3c4tz3*c1c5*tz3
   202:           
   203:                  return
   204:                  end


