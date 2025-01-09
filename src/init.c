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
  x = region.x - 1, y = region.y - 1, z = region.z - 1;//ok
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
  sa->Txx = malloc(sizeof(double **) * (sr.Txx.x + 1));
  sa->Txxx = malloc(sizeof(double **) * (sr.Txx.x + 1));
  sa->Txxy = malloc(sizeof(double **) * (sr.Txx.x + 1));
  sa->Txxz = malloc(sizeof(double **) * (sr.Txx.x + 1));
  for (i = 0; i <= sr.Txx.x; i++) {
    sa->Txx[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    sa->Txxx[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    sa->Txxy[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    sa->Txxz[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
  }
  for (i = 0; i <= sr.Txx.x; i++) {
    for (j = 0; j <= sr.Txx.y; j++) {
      sa->Txx[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      sa->Txxx[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      sa->Txxy[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      sa->Txxz[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
    }
  }
  sa->Tyy = malloc(sizeof(double **) * (sr.Tyy.x + 1));
  sa->Tyyx = malloc(sizeof(double **) * (sr.Tyy.x + 1));
  sa->Tyyy = malloc(sizeof(double **) * (sr.Tyy.x + 1));
  sa->Tyyz = malloc(sizeof(double **) * (sr.Tyy.x + 1));
  for (i = 0; i <= sr.Tyy.x; i++) {
    sa->Tyy[i] = malloc(sizeof(double *) * (sr.Tyy.y + 1));
    sa->Tyyx[i] = malloc(sizeof(double *) * (sr.Tyy.y + 1));
    sa->Tyyy[i] = malloc(sizeof(double *) * (sr.Tyy.y + 1));
    sa->Tyyz[i] = malloc(sizeof(double *) * (sr.Tyy.y + 1));
  }
  for (i = 0; i <= sr.Tyy.x; i++) {
    for (j = 0; j <= sr.Tyy.y; j++) {
      sa->Tyy[i][j] = malloc(sizeof(double) * (sr.Tyy.z + 1));
      sa->Tyyx[i][j] = malloc(sizeof(double) * (sr.Tyy.z + 1));
      sa->Tyyy[i][j] = malloc(sizeof(double) * (sr.Tyy.z + 1));
      sa->Tyyz[i][j] = malloc(sizeof(double) * (sr.Tyy.z + 1));
    }
  }
  sa->Tzz = malloc(sizeof(double **) * (sr.Tzz.x + 1));
  sa->Tzzx = malloc(sizeof(double **) * (sr.Tzz.x + 1));
  sa->Tzzy = malloc(sizeof(double **) * (sr.Tzz.x + 1));
  sa->Tzzz = malloc(sizeof(double **) * (sr.Tzz.x + 1));
  for (i = 0; i <= sr.Tzz.x; i++) {
    sa->Tzz[i] = malloc(sizeof(double *) * (sr.Tzz.y + 1));
    sa->Tzzx[i] = malloc(sizeof(double *) * (sr.Tzz.y + 1));
    sa->Tzzy[i] = malloc(sizeof(double *) * (sr.Tzz.y + 1));
    sa->Tzzz[i] = malloc(sizeof(double *) * (sr.Tzz.y + 1));
  }
  for (i = 0; i <= sr.Tzz.x; i++) {
    for (j = 0; j <= sr.Tzz.y; j++) {
      sa->Tzz[i][j] = malloc(sizeof(double) * (sr.Tzz.z + 1));
      sa->Tzzx[i][j] = malloc(sizeof(double) * (sr.Tzz.z + 1));
      sa->Tzzy[i][j] = malloc(sizeof(double) * (sr.Tzz.z + 1));
      sa->Tzzz[i][j] = malloc(sizeof(double) * (sr.Tzz.z + 1));
    }
  }
}

void initTauArr(TauArr *ta, TauRan tr) {
  int i, j, k;
  ta->Txy = malloc(sizeof(double **) * (tr.Txy.x + 1));
  ta->Txyx = malloc(sizeof(double **) * (tr.Txy.x + 1));
  ta->Txyy = malloc(sizeof(double **) * (tr.Txy.x + 1));
  for (i = 0; i <= tr.Txy.x; i++) {
    ta->Txy[i] = malloc(sizeof(double *) * (tr.Txy.y + 1));
    ta->Txyx[i] = malloc(sizeof(double *) * (tr.Txy.y + 1));
    ta->Txyy[i] = malloc(sizeof(double *) * (tr.Txy.y + 1));
  }
  for (i = 0; i <= tr.Txy.x; i++) {
    for (j = 0; j <= tr.Txy.y; j++) {
      ta->Txy[i][j] = malloc(sizeof(double) * (tr.Txy.z + 1));
      ta->Txyx[i][j] = malloc(sizeof(double) * (tr.Txy.z + 1));
      ta->Txyy[i][j] = malloc(sizeof(double) * (tr.Txy.z + 1));
    }
  }
  ta->Tyz = malloc(sizeof(double **) * (tr.Tyz.x + 1));
  ta->Tyzy = malloc(sizeof(double **) * (tr.Tyz.x + 1));
  ta->Tyzz = malloc(sizeof(double **) * (tr.Tyz.x + 1));
  for (i = 0; i <= tr.Tyz.x; i++) {
    ta->Tyz[i] = malloc(sizeof(double *) * (tr.Tyz.y + 1));
    ta->Tyzy[i] = malloc(sizeof(double *) * (tr.Tyz.y + 1));
    ta->Tyzz[i] = malloc(sizeof(double *) * (tr.Tyz.y + 1));
  }
  for (i = 0; i <= tr.Tyz.x; i++) {
    for (j = 0; j <= tr.Tyz.y; j++) {
      ta->Tyz[i][j] = malloc(sizeof(double) * (tr.Tyz.z + 1));
      ta->Tyzy[i][j] = malloc(sizeof(double) * (tr.Tyz.z + 1));
      ta->Tyzz[i][j] = malloc(sizeof(double) * (tr.Tyz.z + 1));
    }
  }
  ta->Tzx = malloc(sizeof(double **) * (tr.Tzx.x + 1));
  ta->Tzxz = malloc(sizeof(double **) * (tr.Tzx.x + 1));
  ta->Tzxx = malloc(sizeof(double **) * (tr.Tzx.x + 1));
  for (i = 0; i <= tr.Tzx.x; i++) {
    ta->Tzx[i] = malloc(sizeof(double *) * (tr.Tzx.y + 1));
    ta->Tzxz[i] = malloc(sizeof(double *) * (tr.Tzx.y + 1));
    ta->Tzxx[i] = malloc(sizeof(double *) * (tr.Tzx.y + 1));
  }
  for (i = 0; i <= tr.Tzx.x; i++) {
    for (j = 0; j <= tr.Tzx.y; j++) {
      ta->Tzx[i][j] = malloc(sizeof(double) * (tr.Tzx.z + 1));
      ta->Tzxz[i][j] = malloc(sizeof(double) * (tr.Tzx.z + 1));
      ta->Tzxx[i][j] = malloc(sizeof(double) * (tr.Tzx.z + 1));
    }
  }
}

void initVelArr(VelArr *va, VelRan vr) {
  int i, j, k;
  va->Vx = malloc(sizeof(double **) * (vr.Vx.x + 1));
  va->Vxx = malloc(sizeof(double **) * (vr.Vx.x + 1));
  va->Vxy = malloc(sizeof(double **) * (vr.Vx.x + 1));
  va->Vxz = malloc(sizeof(double **) * (vr.Vx.x + 1));
  for (i = 0; i <= vr.Vx.x; i++) {
    va->Vx[i] = malloc(sizeof(double *) * (vr.Vx.y + 1));
    va->Vxx[i] = malloc(sizeof(double *) * (vr.Vx.y + 1));
    va->Vxy[i] = malloc(sizeof(double *) * (vr.Vx.y + 1));
    va->Vxz[i] = malloc(sizeof(double *) * (vr.Vx.y + 1));
  }
  for (i = 0; i <= vr.Vx.x; i++) {
    for (j = 0; j <= vr.Vx.y; j++) {
      va->Vx[i][j] = malloc(sizeof(double) * (vr.Vx.z + 1));
      va->Vxx[i][j] = malloc(sizeof(double) * (vr.Vx.z + 1));
      va->Vxy[i][j] = malloc(sizeof(double) * (vr.Vx.z + 1));
      va->Vxz[i][j] = malloc(sizeof(double) * (vr.Vx.z + 1));
    }
  }
  va->Vy = malloc(sizeof(double **) * (vr.Vy.x + 1));
  va->Vyx = malloc(sizeof(double **) * (vr.Vy.x + 1));
  va->Vyy = malloc(sizeof(double **) * (vr.Vy.x + 1));
  va->Vyz = malloc(sizeof(double **) * (vr.Vy.x + 1));
  for (i = 0; i <= vr.Vy.x; i++) {
    va->Vy[i] = malloc(sizeof(double *) * (vr.Vy.y + 1));
    va->Vyx[i] = malloc(sizeof(double *) * (vr.Vy.y + 1));
    va->Vyy[i] = malloc(sizeof(double *) * (vr.Vy.y + 1));
    va->Vyz[i] = malloc(sizeof(double *) * (vr.Vy.y + 1));
  }
  for (i = 0; i <= vr.Vy.x; i++) {
    for (j = 0; j <= vr.Vy.y; j++) {
      va->Vy[i][j] = malloc(sizeof(double) * (vr.Vy.z + 1));
      va->Vyx[i][j] = malloc(sizeof(double) * (vr.Vy.z + 1));
      va->Vyy[i][j] = malloc(sizeof(double) * (vr.Vy.z + 1));
      va->Vyz[i][j] = malloc(sizeof(double) * (vr.Vy.z + 1));
    }
  }
  va->Vz = malloc(sizeof(double **) * (vr.Vz.x + 1));
  va->Vzx = malloc(sizeof(double **) * (vr.Vz.x + 1));
  va->Vzy = malloc(sizeof(double **) * (vr.Vz.x + 1));
  va->Vzz = malloc(sizeof(double **) * (vr.Vz.x + 1));
  for (i = 0; i <= vr.Vz.x; i++) {
    va->Vz[i] = malloc(sizeof(double *) * (vr.Vz.y + 1));
    va->Vzx[i] = malloc(sizeof(double *) * (vr.Vz.y + 1));
    va->Vzy[i] = malloc(sizeof(double *) * (vr.Vz.y + 1));
    va->Vzz[i] = malloc(sizeof(double *) * (vr.Vz.y + 1));
  }
  for (i = 0; i <= vr.Vz.x; i++) {
    for (j = 0; j <= vr.Vz.y; j++) {
      va->Vz[i][j] = malloc(sizeof(double) * (vr.Vz.z + 1));
      va->Vzx[i][j] = malloc(sizeof(double) * (vr.Vz.z + 1));
      va->Vzy[i][j] = malloc(sizeof(double) * (vr.Vz.z + 1));
      va->Vzz[i][j] = malloc(sizeof(double) * (vr.Vz.z + 1));
    }
  }
}

void initBefAft(BefAft *ba, Range ran) {
  initSigArr(&ba->sa, ran.sr);
  initTauArr(&ba->ta, ran.tr);
  initVelArr(&ba->va, ran.vr);
}

void initInpalse(Inpaluse *ip, SigRan sr, Pml pml) {
  int i, j, k;
  // initCoord(&ip->in, x + pml.pl1.x - 1, y + pml.pl1.y - 1, z + pml.pl1.z - 1);//ok
  ip->Txx = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ip->Tyy = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ip->Tzz = malloc(sizeof(double **) * (sr.Txx.x + 1));
  for (i = 0; i <= sr.Txx.x; i++) {
    ip->Txx[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ip->Tyy[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ip->Tzz[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
  }
  for (i = 0; i <= sr.Txx.x; i++) {
    for (j = 0; j <= sr.Txx.y; j++) {
      ip->Txx[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ip->Tyy[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ip->Tzz[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
    }
  }
  for (k = 0; k <= sr.Txx.z; k++) {
    for (j = 0; j <= sr.Txx.y; j++) {
      for (i = 0; i <= sr.Txx.x; i++) {
        ip->Txx[i][j][k] = 0.;
        ip->Tyy[i][j][k] = 0.;
        ip->Tzz[i][j][k] = 0.;
      }
    }
  }
}

void initMedArr(MedArr *ma, SigRan sr) {
  int i, j, k;
  ma->ramda = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->mu = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->c11 = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->rho = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetaxx = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetaxy = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetaxz = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetayx = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetayy = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetayz = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetazx = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetazy = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetazz = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->gamma = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->khi = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->xi11 = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetadx = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetady = malloc(sizeof(double **) * (sr.Txx.x + 1));
  ma->zetadz = malloc(sizeof(double **) * (sr.Txx.x + 1));
  for (i = 0; i <= sr.Txx.x; i++) {
    ma->ramda[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->mu[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->c11[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->rho[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetaxx[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetaxy[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetaxz[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetayx[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetayy[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetayz[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetazx[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetazy[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetazz[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->gamma[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->khi[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->xi11[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetadx[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetady[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
    ma->zetadz[i] = malloc(sizeof(double *) * (sr.Txx.y + 1));
  }
  for (i = 0; i <= sr.Txx.x; i++) {
    for (j = 0; j <= sr.Txx.y; j++) {
      ma->ramda[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->mu[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->c11[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->rho[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetaxx[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetaxy[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetaxz[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetayx[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetayy[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetayz[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetazx[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetazy[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetazz[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->gamma[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->khi[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->xi11[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetadx[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetady[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
      ma->zetadz[i][j] = malloc(sizeof(double) * (sr.Txx.z + 1));
    }
  }
}

void initClack(Object *clack, Medium med, Pml *pml, int spx, int spy, int spz, int ranx, int rany, int ranz) {
  clack->med = med;
  spx = spx + pml->pl1.x - 1, spy = spy + pml->pl1.y - 1, spz = spz + pml->pl1.z - 1;//ok
  initCoord(&clack->sp, spx, spy, spz);
  initCoord(&clack->range, ranx, rany, ranz);
}
