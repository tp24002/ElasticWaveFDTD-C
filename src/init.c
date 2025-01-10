#define _USE_MATH_DEFINES
#include "../header/init.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../header/struct.h"

#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define MAX(i, j) (((i) > (j)) ? (i) : (j))

void initMedium(Medium *med) {
  for (int mednum = 0; mednum < E_M_END; mednum++) {
    switch (mednum) {
      case E_AIR:
        med[mednum].rho = 1.205;  // 密度
        med[mednum].K = 1.422e5;  // 体積弾性率
        med[mednum].E = 0.;       // ヤング率
        med[mednum].G = 0.;       // 剛性率
        med[mednum].nu = 0.;     // ポアソン比
        med[mednum].ramda = med[mednum].K - 2. / 3. * med[mednum].G;  // 第1ラメ定数(1.422000e+05)
        med[mednum].zeta = 5.;    //セル間摩擦係数
        med[mednum].gamma = 1.8e-5;   //第1粘性
        med[mednum].khi = 0.;         //第2粘性
        med[mednum].eta = 0.;
        med[mednum].omega = 0.;
        break;
      case E_CON:
        med[mednum].rho = 2400.;  // 密度
        med[mednum].E = 2.4e10;   // ヤング率
        med[mednum].nu = 0.2;     // ポアソン比
        // 第2ラメ定数，剛性率（1.000000e+10）mu
        med[mednum].G = med[mednum].E / 2. / (1. + med[mednum].nu);

        // 体積弾性率（1.333333e+10）
        med[mednum].K = med[mednum].E / 3. / (1. - 2. * med[mednum].nu);

        // 第1ラメ定数（6.666667e+09）
        med[mednum].ramda = med[mednum].E * med[mednum].nu / (1. + med[mednum].nu) / (1. - 2. * med[mednum].nu);

        med[mednum].zeta = 2.5e4;//セル間摩擦係数
        med[mednum].eta = 0.005;//粘性定数算出係数(損失係数)
        med[mednum].omega = 2. * M_PI * 32.;//粘性定数算出係数(角周波数)

        // 第1粘性定数（2.486796e+05）
        med[mednum].gamma = med[mednum].eta * med[mednum].G / med[mednum].omega;
        
        // 第2粘性定数（1.657864e+05）
        med[mednum].khi = med[mednum].eta * med[mednum].ramda / med[mednum].omega;
        break;
      default:
        break;
    }
  }
}

void initCoord(Coord *co, int x, int y, int z) {
  co->x = x;
  co->y = y;
  co->z = z;
}

int getId(Coord co, int x, int y, int z) {
  int id = co.x * co.y * z + co.x * y + x;
  return id;
}

//差分間隔
void initDiff(Diff *dif, Medium *med) {
  dif->dx = 0.001;
  dif->dy = 0.001;
  dif->dz = 0.001;
  double tmp;
  double ramda;
  
  for(int i = E_AIR; i < E_M_END - 1; i++){
    tmp = MAX(sqrt((med[i].K + 4. / 3. * med[i].G) / med[i].rho),tmp);
    ramda = MAX(ramda, med[i].ramda);
  }
  printf("v = %lf\n", tmp);
  double n = 0.5;
  dif->dt = MIN(0.9 * dif->dx / (tmp * sqrt(3)), 1 / ramda);
  // dif->dt = n * 3e-9/20;
  // dif->dt = 0.9 * dif->dx / (tmp * sqrt(3));
  // dif->dt = dif->dx / tmp / 100.;
}

void initPml(Pml *pml, Medium *med, Diff dif) {
  pml->ta = 3.;
  pml->fm = 3.574e4;
  double R = 1.e-10;
  double tmp,tmp_v;//max
  initCoord(&pml->pl1, 32, 32, 32);
  initCoord(&pml->pl2, 32, 32, 32);
  //計算領域内最高速度
  for(int i = E_AIR; i < E_M_END - 1; i++){
    tmp_v = MAX(sqrt((med[i].K + 4. / 3. * med[i].G) / med[i].rho),tmp);
  }
  //減衰係数最大値(PML層)
  for (int i = E_AIR + 1; i < E_M_END; i++) {
    tmp = tmp_v * (pml->ta + 1) / (2. * (double)pml->pl1.x * dif.dx) * log(1/R);
    pml->fm = MAX(tmp, pml->fm);
  }
}

void initObject(Object *ob, Medium med, Pml pml, Coord sp, Coord ran) {
  ob->med = med;
  sp.x = sp.x + pml.pl1.x - 1, sp.y = sp.y + pml.pl1.y - 1, sp.z = sp.z + pml.pl1.z - 1;
  initCoord(&ob->sp, sp.x, sp.y, sp.z);
  initCoord(&ob->range, ran.x, ran.y, ran.z);
}

void initRange(Range *ran, Coord region, Pml pml) {
  int x, y, z;
  x = region.x, y = region.y, z = region.z; 
  initCoord(&ran->sr.Txx, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->sr.Tyy, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->sr.Tzz, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->tr.Txy, x + pml.pl1.x + pml.pl2.x + 1, y + pml.pl1.y + pml.pl2.y + 1, z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->tr.Tyz, x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y + 1, z + pml.pl1.z + pml.pl2.z + 1);
  initCoord(&ran->tr.Tzx, x + pml.pl1.x + pml.pl2.x + 1, y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z + 1);
  initCoord(&ran->vr.Vx , x + pml.pl1.x + pml.pl2.x + 1, y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->vr.Vy , x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y + 1, z + pml.pl1.z + pml.pl2.z    );
  initCoord(&ran->vr.Vz , x + pml.pl1.x + pml.pl2.x    , y + pml.pl1.y + pml.pl2.y    , z + pml.pl1.z + pml.pl2.z + 1);
} 

void initSigArr(SigArr *sa, SigRan sr) {
  int i, j, k;
  sa->Txx = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  sa->Txxx = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  sa->Txxy = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  sa->Txxz = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);

  sa->Tyy = malloc(sizeof(double **) * sr.Tyy.x * sr.Tyy.y * sr.Tyy.z);
  sa->Tyyx = malloc(sizeof(double **) * sr.Tyy.x * sr.Tyy.y * sr.Tyy.z);
  sa->Tyyy = malloc(sizeof(double **) * sr.Tyy.x * sr.Tyy.y * sr.Tyy.z);
  sa->Tyyz = malloc(sizeof(double **) * sr.Tyy.x * sr.Tyy.y * sr.Tyy.z);

  sa->Tzz = malloc(sizeof(double **) * sr.Tzz.x * sr.Tzz.y * sr.Tzz.z);
  sa->Tzzx = malloc(sizeof(double **) * sr.Tzz.x * sr.Tzz.y * sr.Tzz.z);
  sa->Tzzy = malloc(sizeof(double **) * sr.Tzz.x * sr.Tzz.y * sr.Tzz.z);
  sa->Tzzz = malloc(sizeof(double **) * sr.Tzz.x * sr.Tzz.y * sr.Tzz.z);
}

void initTauArr(TauArr *ta, TauRan tr) {
  int i, j, k;
  ta->Txy = malloc(sizeof(double **) * tr.Txy.x * tr.Txy.y * tr.Txy.z);
  ta->Txyx = malloc(sizeof(double **) * tr.Txy.x * tr.Txy.y * tr.Txy.z);
  ta->Txyy = malloc(sizeof(double **) * tr.Txy.x * tr.Txy.y * tr.Txy.z);

  ta->Tyz = malloc(sizeof(double **) * tr.Tyz.x * tr.Tyz.y * tr.Tyz.z);
  ta->Tyzy = malloc(sizeof(double **) * tr.Tyz.x * tr.Tyz.y * tr.Tyz.z);
  ta->Tyzz = malloc(sizeof(double **) * tr.Tyz.x * tr.Tyz.y * tr.Tyz.z);

  ta->Tzx = malloc(sizeof(double **) * tr.Tzx.x * tr.Tzx.y * tr.Tzx.z);
  ta->Tzxz = malloc(sizeof(double **) * tr.Tzx.x * tr.Tzx.y * tr.Tzx.z);
  ta->Tzxx = malloc(sizeof(double **) * tr.Tzx.x * tr.Tzx.y * tr.Tzx.z);
}

void initVelArr(VelArr *va, VelRan vr) {
  int i, j, k;
  va->Vx = malloc(sizeof(double **) * vr.Vx.x * vr.Vx.y * vr.Vx.z);
  va->Vxx = malloc(sizeof(double **) * vr.Vx.x * vr.Vx.y * vr.Vx.z);
  va->Vxy = malloc(sizeof(double **) * vr.Vx.x * vr.Vx.y * vr.Vx.z);
  va->Vxz = malloc(sizeof(double **) * vr.Vx.x * vr.Vx.y * vr.Vx.z);

  va->Vy = malloc(sizeof(double **) * vr.Vy.x * vr.Vy.y * vr.Vy.z);
  va->Vyx = malloc(sizeof(double **) * vr.Vy.x * vr.Vy.y * vr.Vy.z);
  va->Vyy = malloc(sizeof(double **) * vr.Vy.x * vr.Vy.y * vr.Vy.z);
  va->Vyz = malloc(sizeof(double **) * vr.Vy.x * vr.Vy.y * vr.Vy.z);

  va->Vz = malloc(sizeof(double **) * vr.Vz.x * vr.Vz.y * vr.Vz.z);
  va->Vzx = malloc(sizeof(double **) * vr.Vz.x * vr.Vz.y * vr.Vz.z);
  va->Vzy = malloc(sizeof(double **) * vr.Vz.x * vr.Vz.y * vr.Vz.z);
  va->Vzz = malloc(sizeof(double **) * vr.Vz.x * vr.Vz.y * vr.Vz.z);

}

void initBefAft(BefAft *ba, Range ran) {
  initSigArr(&ba->sa, ran.sr);
  initTauArr(&ba->ta, ran.tr);
  initVelArr(&ba->va, ran.vr);
}

void initInpalse(Inpaluse *ip, SigRan sr, Pml pml) {
  int i, j, k;
  // initCoord(&ip->in, x + pml.pl1.x - 1, y + pml.pl1.y - 1, z + pml.pl1.z - 1);//ok
  ip->Txx = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ip->Tyy = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ip->Tzz = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
}

void initMedArr(MedArr *ma, SigRan sr) {
  int i, j, k;
  ma->ramda = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->mu = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->c11 = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->rho = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->zetaxx = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->zetaxy = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->zetaxz = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->zetayx = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->zetayy = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->zetayz = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->zetazx = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->zetazy = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->zetazz = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->gamma = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->khi = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->xi11 = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->zetadx = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->zetady = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
  ma->zetadz = malloc(sizeof(double **) * sr.Txx.x * sr.Txx.y * sr.Txx.z);
}

void initClack(Object *clack, Medium med, Pml *pml, int spx, int spy, int spz, int ranx, int rany, int ranz) {
  clack->med = med;
  spx = spx + pml->pl1.x - 1, spy = spy + pml->pl1.y - 1, spz = spz + pml->pl1.z - 1;//ok
  initCoord(&clack->sp, spx, spy, spz);
  initCoord(&clack->range, ranx, rany, ranz);
}
