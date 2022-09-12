NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:16 2020
FILE NAME: set_constants.f

PROCEDURE NAME: SET_CONSTANTS
DIAGNOSTIC LIST

 LINE              DIAGNOSTIC MESSAGE

  No diagnostic message.


NEC Fortran Compiler (3.0.4) for Vector Engine   Thu Jul 23 22:24:16 2020
FILE NAME: set_constants.f

PROCEDURE NAME: SET_CONSTANTS
FORMAT LIST

 LINE   LOOP      STATEMENT

     4:                 subroutine  set_constants
     5:           
     6:           c---------------------------------------------------------------------
     7:           c---------------------------------------------------------------------
     8:           
     9:                 use bt_data
    10:                 implicit none
    11:            
    12:                 ce(1,1)  = 2.0d0
    13:                 ce(1,2)  = 0.0d0
    14:                 ce(1,3)  = 0.0d0
    15:                 ce(1,4)  = 4.0d0
    16:                 ce(1,5)  = 5.0d0
    17:                 ce(1,6)  = 3.0d0
    18:                 ce(1,7)  = 0.5d0
    19:                 ce(1,8)  = 0.02d0
    20:                 ce(1,9)  = 0.01d0
    21:                 ce(1,10) = 0.03d0
    22:                 ce(1,11) = 0.5d0
    23:                 ce(1,12) = 0.4d0
    24:                 ce(1,13) = 0.3d0
    25:            
    26:                 ce(2,1)  = 1.0d0
    27:                 ce(2,2)  = 0.0d0
    28:                 ce(2,3)  = 0.0d0
    29:                 ce(2,4)  = 0.0d0
    30:                 ce(2,5)  = 1.0d0
    31:                 ce(2,6)  = 2.0d0
    32:                 ce(2,7)  = 3.0d0
    33:                 ce(2,8)  = 0.01d0
    34:                 ce(2,9)  = 0.03d0
    35:                 ce(2,10) = 0.02d0
    36:                 ce(2,11) = 0.4d0
    37:                 ce(2,12) = 0.3d0
    38:                 ce(2,13) = 0.5d0
    39:           
    40:                 ce(3,1)  = 2.0d0
    41:                 ce(3,2)  = 2.0d0
    42:                 ce(3,3)  = 0.0d0
    43:                 ce(3,4)  = 0.0d0
    44:                 ce(3,5)  = 0.0d0
    45:                 ce(3,6)  = 2.0d0
    46:                 ce(3,7)  = 3.0d0
    47:                 ce(3,8)  = 0.04d0
    48:                 ce(3,9)  = 0.03d0
    49:                 ce(3,10) = 0.05d0
    50:                 ce(3,11) = 0.3d0
    51:                 ce(3,12) = 0.5d0
    52:                 ce(3,13) = 0.4d0
    53:           
    54:                 ce(4,1)  = 2.0d0
    55:                 ce(4,2)  = 2.0d0
    56:                 ce(4,3)  = 0.0d0
    57:                 ce(4,4)  = 0.0d0
    58:                 ce(4,5)  = 0.0d0
    59:                 ce(4,6)  = 2.0d0
    60:                 ce(4,7)  = 3.0d0
    61:                 ce(4,8)  = 0.03d0
    62:                 ce(4,9)  = 0.05d0
    63:                 ce(4,10) = 0.04d0
    64:                 ce(4,11) = 0.2d0
    65:                 ce(4,12) = 0.1d0
    66:                 ce(4,13) = 0.3d0
    67:           
    68:                 ce(5,1)  = 5.0d0
    69:                 ce(5,2)  = 4.0d0
    70:                 ce(5,3)  = 3.0d0
    71:                 ce(5,4)  = 2.0d0
    72:                 ce(5,5)  = 0.1d0
    73:                 ce(5,6)  = 0.4d0
    74:                 ce(5,7)  = 0.3d0
    75:                 ce(5,8)  = 0.05d0
    76:                 ce(5,9)  = 0.04d0
    77:                 ce(5,10) = 0.03d0
    78:                 ce(5,11) = 0.1d0
    79:                 ce(5,12) = 0.3d0
    80:                 ce(5,13) = 0.2d0
    81:           
    82:                 c1 = 1.4d0
    83:                 c2 = 0.4d0
    84:                 c3 = 0.1d0
    85:                 c4 = 1.0d0
    86:                 c5 = 1.4d0
    87:           
    88:                 dnxm1 = 1.0d0 / dble(grid_points(1)-1)
    89:                 dnym1 = 1.0d0 / dble(grid_points(2)-1)
    90:                 dnzm1 = 1.0d0 / dble(grid_points(3)-1)
    91:           
    92:                 c1c2 = c1 * c2
    93:                 c1c5 = c1 * c5
    94:                 c3c4 = c3 * c4
    95:                 c1345 = c1c5 * c3c4
    96:           
    97:                 conz1 = (1.0d0-c1c5)
    98:           
    99:                 tx1 = 1.0d0 / (dnxm1 * dnxm1)
   100:                 tx2 = 1.0d0 / (2.0d0 * dnxm1)
   101:                 tx3 = 1.0d0 / dnxm1
   102:           
   103:                 ty1 = 1.0d0 / (dnym1 * dnym1)
   104:                 ty2 = 1.0d0 / (2.0d0 * dnym1)
   105:                 ty3 = 1.0d0 / dnym1
   106:            
   107:                 tz1 = 1.0d0 / (dnzm1 * dnzm1)
   108:                 tz2 = 1.0d0 / (2.0d0 * dnzm1)
   109:                 tz3 = 1.0d0 / dnzm1
   110:           
   111:                 dx1 = 0.75d0
   112:                 dx2 = 0.75d0
   113:                 dx3 = 0.75d0
   114:                 dx4 = 0.75d0
   115:                 dx5 = 0.75d0
   116:           
   117:                 dy1 = 0.75d0
   118:                 dy2 = 0.75d0
   119:                 dy3 = 0.75d0
   120:                 dy4 = 0.75d0
   121:                 dy5 = 0.75d0
   122:           
   123:                 dz1 = 1.0d0
   124:                 dz2 = 1.0d0
   125:                 dz3 = 1.0d0
   126:                 dz4 = 1.0d0
   127:                 dz5 = 1.0d0
   128:           
   129:                 dxmax = dmax1(dx3, dx4)
   130:                 dymax = dmax1(dy2, dy4)
   131:                 dzmax = dmax1(dz2, dz3)
   132:           
   133:                 dssp = 0.25d0 * dmax1(dx1, dmax1(dy1, dz1) )
   134:           
   135:                 c4dssp = 4.0d0 * dssp
   136:                 c5dssp = 5.0d0 * dssp
   137:           
   138:                 dttx1 = dt*tx1
   139:                 dttx2 = dt*tx2
   140:                 dtty1 = dt*ty1
   141:                 dtty2 = dt*ty2
   142:                 dttz1 = dt*tz1
   143:                 dttz2 = dt*tz2
   144:           
   145:                 c2dttx1 = 2.0d0*dttx1
   146:                 c2dtty1 = 2.0d0*dtty1
   147:                 c2dttz1 = 2.0d0*dttz1
   148:           
   149:                 dtdssp = dt*dssp
   150:           
   151:                 comz1  = dtdssp
   152:                 comz4  = 4.0d0*dtdssp
   153:                 comz5  = 5.0d0*dtdssp
   154:                 comz6  = 6.0d0*dtdssp
   155:           
   156:                 c3c4tx3 = c3c4*tx3
   157:                 c3c4ty3 = c3c4*ty3
   158:                 c3c4tz3 = c3c4*tz3
   159:           
   160:                 dx1tx1 = dx1*tx1
   161:                 dx2tx1 = dx2*tx1
   162:                 dx3tx1 = dx3*tx1
   163:                 dx4tx1 = dx4*tx1
   164:                 dx5tx1 = dx5*tx1
   165:            
   166:                 dy1ty1 = dy1*ty1
   167:                 dy2ty1 = dy2*ty1
   168:                 dy3ty1 = dy3*ty1
   169:                 dy4ty1 = dy4*ty1
   170:                 dy5ty1 = dy5*ty1
   171:            
   172:                 dz1tz1 = dz1*tz1
   173:                 dz2tz1 = dz2*tz1
   174:                 dz3tz1 = dz3*tz1
   175:                 dz4tz1 = dz4*tz1
   176:                 dz5tz1 = dz5*tz1
   177:           
   178:                 c2iv  = 2.5d0
   179:                 con43 = 4.0d0/3.0d0
   180:                 con16 = 1.0d0/6.0d0
   181:            
   182:                 xxcon1 = c3c4tx3*con43*tx3
   183:                 xxcon2 = c3c4tx3*tx3
   184:                 xxcon3 = c3c4tx3*conz1*tx3
   185:                 xxcon4 = c3c4tx3*con16*tx3
   186:                 xxcon5 = c3c4tx3*c1c5*tx3
   187:           
   188:                 yycon1 = c3c4ty3*con43*ty3
   189:                 yycon2 = c3c4ty3*ty3
   190:                 yycon3 = c3c4ty3*conz1*ty3
   191:                 yycon4 = c3c4ty3*con16*ty3
   192:                 yycon5 = c3c4ty3*c1c5*ty3
   193:           
   194:                 zzcon1 = c3c4tz3*con43*tz3
   195:                 zzcon2 = c3c4tz3*tz3
   196:                 zzcon3 = c3c4tz3*conz1*tz3
   197:                 zzcon4 = c3c4tz3*con16*tz3
   198:                 zzcon5 = c3c4tz3*c1c5*tz3
   199:           
   200:                 return
   201:                 end


