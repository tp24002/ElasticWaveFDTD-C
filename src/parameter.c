#include "../header/struct.h"
#include "../header/init.h"
#include "../header/parameter.h"

#include <stdio.h>

void parameter(Range *ran,Diff *dif, Pml *pml,Medium *med,Object *con,Object *clack,Coord *out,Inpaluse *ip){
   initMedium(med);
   initDiff(dif, med);
   initPml(pml, med, *dif);//pml 32*2
   Coord region;
   //計算領域
   // initCoord(&region, 202, 202, 202);
   initCoord(&region, 15, 15, 15);
   initRange(ran, region, *pml);
   // //空気
   // Coord air_st, air_size;
   // initCoord(&air_st, 0, 0, 0);
   // initCoord(&air_size, ran->sr.Txx.x, ran->sr.Txx.y, ran->sr.Txx.z);
   //コンクリート
   Coord con_st, con_size;
   // initCoord(&con_st, 2, 2, 2);
   // initCoord(&con_size, 200, 200, 200);
   initCoord(&con_st, 3, 3, 3);
   initCoord(&con_size, 11, 11, 11);
   // 1 < clack_st < con - 1 (表面欠陥は今は考えない)
   //欠陥
   Coord clack_st, clack_size;
   // initCoord(&clack_st, 77, 77, 102);
   // initCoord(&clack_size, 50, 50, 100);
   initCoord(&clack_st, 1, 1, 1);
   initCoord(&clack_size, 1, 1, 1);
   ip->freq = 1.0e7;
   ip->mode = E_SINE;//E_SINE,E_RCOS

   Coord halfcon;
   Coord center;
   initCoord(&halfcon,(con_size.x - 1) / 2,(con_size.y - 1) / 2,(con_size.z - 1) / 2);
   initCoord(&center,(region.x + 2 * pml->pl1.x - 1) / 2,(region.y + 2 * pml->pl1.y - 1) / 2,(region.z + 2 * pml->pl1.z - 1) / 2);
   //上面中心
   initCoord(&(ip->in),center.x,center.y,center.z + halfcon.z);
   
   // 一旦コンクリートのxyの小さい方の半分に合わせて計測地点変更
   // 一つのデバイスとしてまとめたい->できるだけ小型化(ある一定値でできるように)
   int halfconmin;
   if(halfcon.x < halfcon.y) {
      halfconmin = (halfcon.x / 2); 
   } else {
      halfconmin = (halfcon.y / 2); 
   }
   halfconmin = 3;
   //  int halfconmin2 = 25;
   //  int halfconmin2 = 25;
   initCoord(&out[0],center.x - halfconmin, center.y, ip->in.z);
   initCoord(&out[1],center.x + halfconmin, center.y, ip->in.z);
   initCoord(&out[2],center.x, center.y - halfconmin, ip->in.z);
   initCoord(&out[3],center.x, center.y + halfconmin, ip->in.z);
   
   initObject(con, med[E_CON], *pml, con_st, con_size);
   initObject(clack, med[E_AIR], *pml, clack_st, clack_size);
   // for(int i = 0; i < clacknum; i++) {
   //     initObject(&clack[i], med[E_AIR], *pml, );
   // }
}
