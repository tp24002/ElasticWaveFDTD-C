#define _USE_MATH_DEFINES

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "../header/struct.h"
#include "../header/init.h"
#include "../header/update.h"

void Txx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran, Inpaluse ip, int t) {
  int i, j, k;
  int Txximax = ran.sr.Txx.x, Txxjmax = ran.sr.Txx.y, Txxkmax = ran.sr.Txx.z;
  int id = getId(ran.sr.Txx, ip.in.x, ip.in.y, ip.in.z);
  if (ip.mode == E_SINE) {
    ip.Txx[id] = 0;
  } else if (ip.mode == E_RCOS) {
    if (t < 1. / ip.freq / dif.dt) {
      ip.Txx[id] = 0;///* 8.e3 * 0.5 * */(1. - cos(2. * M_PI * ip.freq * (double)t * dif.dt)) / 2.;
    } else {
      ip.Txx[id] = 0.;
    }
  } else {
    ip.Txx[id] = 0.;
  }
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        int idtxx = getId(ran.sr.Txx, i, j, k);
        int idvx  = getId(ran.vr.Vx, i, j, k);
        int idvxi = getId(ran.vr.Vx, i + 1, j, k);
        int idvy  = getId(ran.vr.Vy, i, j, k);
        int idvyj = getId(ran.vr.Vy, i, j + 1, k);
        int idvz  = getId(ran.vr.Vz, i, j, k);
        int idvzk = getId(ran.vr.Vz, i, j, k + 1);
        aft->sa.Txxx[idtxx] = (2. - ma.zetadx[idtxx] * dif.dt) / (2. + ma.zetadx[idtxx] * dif.dt) * bef->sa.Txxx[idtxx]
         + 2. * (ma.c11[idtxx] * dif.dt + ma.xi11[idtxx]) / (2. + ma.zetadx[idtxx] * dif.dt) * (aft->va.Vx[idvxi] - aft->va.Vx[idvx]) / dif.dx 
          - 2. * ma.xi11[idtxx] / (2. + ma.zetadx[idtxx] * dif.dt) * (bef->va.Vx[idvxi] - bef->va.Vx[idvx]) / dif.dx;

        aft->sa.Txxy[idtxx] = (2. - ma.zetady[idtxx] * dif.dt) / (2. + ma.zetady[idtxx] * dif.dt) * bef->sa.Txxy[idtxx]
         + 2. * (ma.ramda[idtxx] * dif.dt + ma.khi[idtxx]) / (2. + ma.zetady[idtxx] * dif.dt) * (aft->va.Vy[idvyj] - aft->va.Vy[idvy]) / dif.dy
          - 2. * ma.khi[idtxx] / (2. + ma.zetady[idtxx] * dif.dt) * (bef->va.Vy[idvyj] - bef->va.Vy[idvy]) / dif.dy;

        aft->sa.Txxz[idtxx] = (2. - ma.zetadz[idtxx] * dif.dt) / (2. + ma.zetadz[idtxx] * dif.dt) * bef->sa.Txxz[idtxx]
         + 2. * (ma.ramda[idtxx] * dif.dt + ma.khi[idtxx]) / (2. + ma.zetadz[idtxx] * dif.dt) * (aft->va.Vz[idvzk] - aft->va.Vz[idvz]) / dif.dz
          - 2. * ma.khi[idtxx] / (2. + ma.zetadz[idtxx] * dif.dt) * (bef->va.Vz[idvzk] - bef->va.Vz[idvz]) / dif.dz;
      }
    }
  }

//全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        int idtxx = getId(ran.sr.Txx, i, j, k);
        aft->sa.Txx[idtxx] = aft->sa.Txxx[idtxx] + aft->sa.Txxy[idtxx] + aft->sa.Txxz[idtxx] + ip.Txx[idtxx];
      }
    }
  }
}

void Tyy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran, Inpaluse ip, int t) {
  int i, j, k;
  int id;
  int Tyyimax = ran.sr.Tyy.x, Tyyjmax = ran.sr.Tyy.y, Tyykmax = ran.sr.Tyy.z;
  id = getId(ran.sr.Tyy, ip.in.x, ip.in.y, ip.in.z);
  if (ip.mode == E_SINE) {
    ip.Tyy[id] = 0;
  } else if (ip.mode == E_RCOS) {
    if (t < 1. / ip.freq / dif.dt) {
      ip.Tyy[id] = 0;///* 8.e3 * 0.5 * */(1. - cos(2. * M_PI * ip.freq * (double)t * dif.dt)) / 2.;
    } else {
      ip.Tyy[id] = 0.;
    }
  } else {
    ip.Tyy[id] = 0.;
  }
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Tyykmax; k++) {
    for (j = 0; j < Tyyjmax; j++) {
      for (i = 0; i < Tyyimax; i++) {
        int idtyy = getId(ran.sr.Tyy, i, j, k);
        int idvx  = getId(ran.vr.Vx, i, j, k);
        int idvxi = getId(ran.vr.Vx, i + 1, j, k);
        int idvy  = getId(ran.vr.Vy, i, j, k);
        int idvyj = getId(ran.vr.Vy, i, j + 1, k);
        int idvz  = getId(ran.vr.Vz, i, j, k);
        int idvzk = getId(ran.vr.Vz, i, j, k + 1);

        aft->sa.Tyyx[idtyy] = (2. - ma.zetadx[idtyy] * dif.dt) / (2. + ma.zetadx[idtyy] * dif.dt) * bef->sa.Tyyx[idtyy]
        + 2. * (ma.ramda[idtyy] * dif.dt + ma.khi[idtyy]) / (2. + ma.zetadx[idtyy] * dif.dt) * (aft->va.Vx[idvxi] - aft->va.Vx[idvx]) / dif.dx
          - 2. * ma.khi[idtyy] / (2. + ma.zetadx[idtyy] * dif.dt) * (bef->va.Vx[idvxi] - bef->va.Vx[idvx]) / dif.dx;

        aft->sa.Tyyy[idtyy] = (2. - ma.zetady[idtyy] * dif.dt) / (2. + ma.zetady[idtyy] * dif.dt) * bef->sa.Tyyy[idtyy]
        + 2. * (ma.c11[idtyy] * dif.dt + ma.xi11[idtyy]) / (2. + ma.zetady[idtyy] * dif.dt) * (aft->va.Vy[idvyj] - aft->va.Vy[idvy]) / dif.dy
          - 2. * ma.xi11[idtyy] / (2. + ma.zetady[idtyy] * dif.dt) * (bef->va.Vy[idvyj] - bef->va.Vy[idvy]) / dif.dy;

        aft->sa.Tyyz[idtyy] = (2. - ma.zetadz[idtyy] * dif.dt) / (2. + ma.zetadz[idtyy] * dif.dt) * bef->sa.Tyyz[idtyy]
        + 2. * (ma.ramda[idtyy] * dif.dt + ma.khi[idtyy]) / (2. + ma.zetadz[idtyy] * dif.dt) * (aft->va.Vz[idvzk] - aft->va.Vz[idvz]) / dif.dz
          - 2. * ma.khi[idtyy] / (2. + ma.zetadz[idtyy] * dif.dt) * (bef->va.Vz[idvzk] - bef->va.Vz[idvz]) / dif.dz;
      }
    }
  }

// 全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Tyykmax; k++) {
    for (j = 0; j < Tyyjmax; j++) {
      for (i = 0; i < Tyyimax; i++) {
        int idtyy = getId(ran.sr.Tyy, i, j, k);
        aft->sa.Tyy[idtyy] = aft->sa.Tyyx[idtyy] + aft->sa.Tyyy[idtyy] + aft->sa.Tyyz[idtyy] + ip.Tyy[idtyy];
      }
    }
  }
}

void Tzz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran, Inpaluse ip, int t) {
  int i, j, k;
  int id;
  int Tzzimax = ran.sr.Tzz.x, Tzzjmax = ran.sr.Tzz.y, Tzzkmax = ran.sr.Tzz.z;
  id = getId(ran.sr.Tzz, ip.in.x, ip.in.y, ip.in.z);
  if (ip.mode == E_SINE) {
    ip.Tzz[id] = (-1) * sin(2. * M_PI * ip.freq * (double)t * dif.dt) / 2.;
  } else if (ip.mode == E_RCOS) {
    if (t < 1. / ip.freq / dif.dt) {
      ip.Tzz[id] = 8.e3 * 0.5 * (-1) * (1. - cos(2. * M_PI * ip.freq * (double)t * dif.dt)) / 2.;
    } else {
      ip.Tzz[id] = 0.;
    }
  } else {
    ip.Tzz[id] = 0.;
  }
  // Tzzの更新式
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Tzzkmax; k++) {
    for (j = 0; j < Tzzjmax; j++) {
      for (i = 0; i < Tzzimax; i++) {
        int idtzz = getId(ran.sr.Tzz, i, j, k);
        int idvx  = getId(ran.vr.Vx, i, j, k);
        int idvxi = getId(ran.vr.Vx, i + 1, j, k);
        int idvy  = getId(ran.vr.Vy, i, j, k);
        int idvyj = getId(ran.vr.Vy, i, j + 1, k);
        int idvz  = getId(ran.vr.Vz, i, j, k);
        int idvzk = getId(ran.vr.Vz, i, j, k + 1);

        aft->sa.Tzzx[idtzz] = (2. - ma.zetadx[idtzz] * dif.dt) / (2. + ma.zetadx[idtzz] * dif.dt) * bef->sa.Tzzx[idtzz]
        + 2. * (ma.ramda[idtzz] * dif.dt + ma.khi[idtzz]) / (2. + ma.zetadx[idtzz] * dif.dt) * (aft->va.Vx[idvxi] - aft->va.Vx[idvx]) / dif.dx
          - 2. * ma.khi[idtzz] / (2. + ma.zetadx[idtzz] * dif.dt) * (bef->va.Vx[idvxi] - bef->va.Vx[idvx]) / dif.dx;

        aft->sa.Tzzy[idtzz] = (2. - ma.zetady[idtzz] * dif.dt) / (2. + ma.zetady[idtzz] * dif.dt) * bef->sa.Tzzy[idtzz]
        + 2. * (ma.ramda[idtzz] * dif.dt + ma.khi[idtzz]) / (2. + ma.zetady[idtzz] * dif.dt) * (aft->va.Vy[idvyj] - aft->va.Vy[idvy]) / dif.dy
          - 2. * ma.khi[idtzz] / (2. + ma.zetady[idtzz] * dif.dt) * (bef->va.Vy[idvyj] - bef->va.Vy[idvy]) / dif.dy;

        aft->sa.Tzzz[idtzz] = (2. - ma.zetadz[idtzz] * dif.dt) / (2. + ma.zetadz[idtzz] * dif.dt) * bef->sa.Tzzz[idtzz]
        + 2. * (ma.c11[idtzz] * dif.dt + ma.xi11[idtzz]) / (2. + ma.zetadz[idtzz] * dif.dt) * (aft->va.Vz[idvzk] - aft->va.Vz[idvz]) / dif.dz
          - 2. * ma.xi11[idtzz] / (2. + ma.zetadz[idtzz] * dif.dt) * (bef->va.Vz[idvzk] - bef->va.Vz[idvz]) / dif.dz;
      }
    }
  }

// 全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Tzzkmax; k++) {
    for (j = 0; j < Tzzjmax; j++) {
      for (i = 0; i < Tzzimax; i++) {
        int idtzz = getId(ran.sr.Tzz, i, j, k);
        aft->sa.Tzz[idtzz] = aft->sa.Tzzx[idtzz] + aft->sa.Tzzy[idtzz] + aft->sa.Tzzz[idtzz] + ip.Tzz[idtzz];
      }
    }
  }
  i = 39, j = 39, k = 43;
  // printf("tzza:%le\n", 2. * (ma.c11[i][j][k] * dif.dt + ma.xi11[i][j][k]) / (2. + ma.zetadz[i][j][k] * dif.dt));
  // printf("tzzb:%le\n", 2. * ma.xi11[i][j][k] / (2. + ma.zetadz[i][j][k] * dif.dt));

}
//垂直応力計算
void Sig(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran, Inpaluse ip, int t) {
  Txx(aft, bef, ma, dif, ran, ip, t);
  Tyy(aft, bef, ma, dif, ran, ip, t);
  Tzz(aft, bef, ma, dif, ran, ip, t);
}

void Txy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran) {
  int i, j, k;
  int Txyimax = ran.tr.Txy.x, Txyjmax = ran.tr.Txy.y, Txykmax = ran.tr.Txy.z;
  double Hzetadx, Hzetady, Hzetadz, Hmu, Hgamma;
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Txykmax; k++) {
    for (j = 1; j < Txyjmax - 1; j++) {
      for (i = 1; i < Txyimax - 1; i++) {
        int idtxy = getId(ran.tr.Txy, i, j, k);
        int idt   = getId(ran.sr.Txx, i, j, k);
        int idti  = getId(ran.sr.Txx, i - 1, j, k);
        int idtj  = getId(ran.sr.Txx, i, j - 1, k);
        int idtij = getId(ran.sr.Txx, i - 1, j - 1, k);
        int idvx  = getId(ran.vr.Vx , i, j, k);
        int idvxj = getId(ran.vr.Vx , i, j - 1, k);
        int idvy  = getId(ran.vr.Vy , i, j, k);
        int idvyi = getId(ran.vr.Vy , i - 1, j, k);
        //PML:減衰係数
        //計算領域:摩擦定数
        Hzetadx = 4. * pow((1. / ma.zetadx[idtij]) + (1. / ma.zetadx[idtj]) + (1. / ma.zetadx[idti]) + (1. / ma.zetadx[idt]), -1.);
        Hzetady = 4. * pow((1. / ma.zetady[idtij]) + (1. / ma.zetady[idtj]) + (1. / ma.zetady[idti]) + (1. / ma.zetady[idt]), -1.);
        //第2ラメ，横弾性係数(剛性率)
        Hmu     = 4. * pow((1. /     ma.mu[idtij]) + (1. /     ma.mu[idtj]) + (1. /     ma.mu[idti]) + (1. /     ma.mu[idt]), -1.);
        //第１粘性定数
        Hgamma  = 4. * pow((1. /  ma.gamma[idtij]) + (1. /  ma.gamma[idtj]) + (1. /  ma.gamma[idti]) + (1. /  ma.gamma[idt]), -1.);

        aft->ta.Txyx[idtxy] = (2. - Hzetadx * dif.dt) / (2. + Hzetadx * dif.dt) * bef->ta.Txyx[idtxy]
         + 2. * (Hmu * dif.dt + Hgamma) / (2. + Hzetadx * dif.dt) * (aft->va.Vy[idvy] - aft->va.Vy[idvyi]) / dif.dx
          - 2. * Hgamma / (2. + Hzetadx * dif.dt) * (bef->va.Vy[idvy] - bef->va.Vy[idvyi]) / dif.dx;

        aft->ta.Txyy[idtxy] = (2. - Hzetady * dif.dt) / (2. + Hzetady * dif.dt) * bef->ta.Txyy[idtxy]
         + 2. * (Hmu * dif.dt + Hgamma) / (2. + Hzetady * dif.dt) * (aft->va.Vx[idvx] - aft->va.Vx[idvxj]) / dif.dy
          - 2. * Hgamma / (2. + Hzetady * dif.dt) * (bef->va.Vx[idvx] - bef->va.Vx[idvxj]) / dif.dy;
      }
    }
  }
//全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Txykmax; k++) {
    for (j = 0; j < Txyjmax; j++) {
      for (i = 0; i < Txyimax; i++) {
        int idtxy = getId(ran.tr.Txy, i, j, k);
        aft->ta.Txy[idtxy] = aft->ta.Txyx[idtxy] + aft->ta.Txyy[idtxy];
      }
    }
  }
}

void Tyz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran) {
  int i, j, k;
  int Tyzimax = ran.tr.Tyz.x, Tyzjmax = ran.tr.Tyz.y, Tyzkmax = ran.tr.Tyz.z;
  double Hzetadx, Hzetady, Hzetadz, Hmu, Hgamma;
#pragma omp parallel for private(i, j, k)
  for (k = 1; k < Tyzkmax - 1; k++) {
    for (j = 1; j < Tyzjmax - 1; j++) {
      for (i = 0; i < Tyzimax; i++) {
        int idtyz = getId(ran.tr.Tyz, i, j, k);
        int idt   = getId(ran.sr.Tyy, i, j, k);
        int idtj  = getId(ran.sr.Tyy, i, j - 1, k);
        int idtk  = getId(ran.sr.Tyy, i, j, k - 1);
        int idtjk = getId(ran.sr.Tyy, i, j - 1, k - 1);
        int idvz  = getId(ran.vr.Vz, i, j, k);
        int idvzj = getId(ran.vr.Vz, i, j - 1, k);
        int idvy  = getId(ran.vr.Vy, i, j, k);
        int idvyk = getId(ran.vr.Vy, i, j, k - 1);

        // PML: 減衰係数
        Hzetady = 4. * pow((1. / ma.zetady[idtjk]) + (1. / ma.zetady[idtj]) + (1. / ma.zetady[idtk]) + (1. / ma.zetady[idt]), -1.);
        Hzetadz = 4. * pow((1. / ma.zetadz[idtjk]) + (1. / ma.zetadz[idtj]) + (1. / ma.zetadz[idtk]) + (1. / ma.zetadz[idt]), -1.);
        // 第2ラメ，横弾性係数(剛性率)
        Hmu     = 4. * pow((1. / ma.mu[idtjk]) + (1. / ma.mu[idtj]) + (1. / ma.mu[idtk]) + (1. / ma.mu[idt]), -1.);
        // 第１粘性定数
        Hgamma  = 4. * pow((1. / ma.gamma[idtjk]) + (1. / ma.gamma[idtj]) + (1. / ma.gamma[idtk]) + (1. / ma.gamma[idt]), -1.);

        aft->ta.Tyzy[idtyz] = (2. - Hzetady * dif.dt) / (2. + Hzetady * dif.dt) * bef->ta.Tyzy[idtyz]
        + 2. * (Hmu * dif.dt + Hgamma) / (2. + Hzetady * dif.dt) * (aft->va.Vz[idvz] - aft->va.Vz[idvzj]) / dif.dy
          - 2. * Hgamma / (2. + Hzetady * dif.dt) * (bef->va.Vz[idvz] - bef->va.Vz[idvzj]) / dif.dy;

        aft->ta.Tyzz[idtyz] = (2. - Hzetadz * dif.dt) / (2. + Hzetadz * dif.dt) * bef->ta.Tyzz[idtyz]
        + 2. * (Hmu * dif.dt + Hgamma) / (2. + Hzetadz * dif.dt) * (aft->va.Vy[idvy] - aft->va.Vy[idvyk]) / dif.dz 
          - 2. * Hgamma / (2. + Hzetadz * dif.dt) * (bef->va.Vy[idvy] - bef->va.Vy[idvyk]) / dif.dz;
      }
    }
  }

// 全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Tyzkmax; k++) {
    for (j = 0; j < Tyzjmax; j++) {
      for (i = 0; i < Tyzimax; i++) {
        int idtyz = getId(ran.tr.Tyz, i, j, k);
        aft->ta.Tyz[idtyz] = aft->ta.Tyzy[idtyz] + aft->ta.Tyzz[idtyz];
      }
    }
  }
}

void Tzx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran) {
  int i, j, k;
  int Tzximax = ran.tr.Tzx.x, Tzxjmax = ran.tr.Tzx.y, Tzxkmax = ran.tr.Tzx.z;
  double Hzetadx, Hzetady, Hzetadz, Hmu, Hgamma;
#pragma omp parallel for private(i, j, k)
  for (k = 1; k < Tzxkmax - 1; k++) {
    for (j = 0; j < Tzxjmax; j++) {
      for (i = 1; i < Tzximax - 1; i++) {
        int idtzx = getId(ran.tr.Tzx, i, j, k);
        int idt   = getId(ran.sr.Txx, i, j, k);
        int idti  = getId(ran.sr.Txx, i - 1, j, k);
        int idtk  = getId(ran.sr.Txx, i, j, k - 1);
        int idtik = getId(ran.sr.Txx, i - 1, j, k - 1);
        int idvx  = getId(ran.vr.Vx, i, j, k);
        int idvxk = getId(ran.vr.Vx, i, j, k - 1);
        int idvz  = getId(ran.vr.Vz, i, j, k);
        int idvzi = getId(ran.vr.Vz, i - 1, j, k);

        // PML: 減衰係数
        Hzetadx = 4. * pow((1. / ma.zetadx[idtik]) + (1. / ma.zetadx[idti]) + (1. / ma.zetadx[idtk]) + (1. / ma.zetadx[idt]), -1.);
        Hzetadz = 4. * pow((1. / ma.zetadz[idtik]) + (1. / ma.zetadz[idti]) + (1. / ma.zetadz[idtk]) + (1. / ma.zetadz[idt]), -1.);
        // 第2ラメ，横弾性係数(剛性率)
        Hmu     = 4. * pow((1. / ma.mu[idtik]) + (1. / ma.mu[idti]) + (1. / ma.mu[idtk]) + (1. / ma.mu[idt]), -1.);
        // 第１粘性定数
        Hgamma  = 4. * pow((1. / ma.gamma[idtik]) + (1. / ma.gamma[idti]) + (1. / ma.gamma[idtk]) + (1. / ma.gamma[idt]), -1.);

        aft->ta.Tzxz[idtzx] = (2. - Hzetadz * dif.dt) / (2. + Hzetadz * dif.dt) * bef->ta.Tzxz[idtzx]
        + 2. * (Hmu * dif.dt + Hgamma) / (2. + Hzetadz * dif.dt) * (aft->va.Vx[idvx] - aft->va.Vx[idvxk]) / dif.dz
          - 2. * Hgamma / (2. + Hzetadz * dif.dt) * (bef->va.Vx[idvx] - bef->va.Vx[idvxk]) / dif.dz;

        aft->ta.Tzxx[idtzx] = (2. - Hzetadx * dif.dt) / (2. + Hzetadx * dif.dt) * bef->ta.Tzxx[idtzx]
        + 2. * (Hmu * dif.dt + Hgamma) / (2. + Hzetadx * dif.dt) * (aft->va.Vz[idvz] - aft->va.Vz[idvzi]) / dif.dx
          - 2. * Hgamma / (2. + Hzetadx * dif.dt) * (bef->va.Vz[idvz] - bef->va.Vz[idvzi]) / dif.dx;
      }
    }
  }

// 全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Tzxkmax; k++) {
    for (j = 0; j < Tzxjmax; j++) {
      for (i = 0; i < Tzximax; i++) {
        int idtzx = getId(ran.tr.Tzx, i, j, k);
        aft->ta.Tzx[idtzx] = aft->ta.Tzxx[idtzx] + aft->ta.Tzxz[idtzx];
      }
    }
  }

}
//せん断応力計算
void Tau(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran) {
  Txy(aft, bef, ma, dif, ran);
  Tyz(aft, bef, ma, dif, ran);
  Tzx(aft, bef, ma, dif, ran);
}

void Vx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran) {
  int i, j, k;
  int Vximax = ran.vr.Vx.x, Vxjmax = ran.vr.Vx.y, Vxkmax = ran.vr.Vx.z;
  double Azetaxx, Azetaxy, Azetaxz, Arho;
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Vxkmax; k++) {
    for (j = 0; j < Vxjmax; j++) {
      for (i = 1; i < Vximax - 1; i++) {
        int idvx   = getId(ran.vr.Vx, i, j, k);
        int idtxx  = getId(ran.sr.Txx, i, j, k);
        int idtxxi = getId(ran.sr.Txx, i - 1, j, k);
        int idtxy  = getId(ran.tr.Txy, i, j, k);
        int idtxyj = getId(ran.tr.Txy, i, j + 1, k);
        int idtzx  = getId(ran.tr.Tzx, i, j, k);
        int idtzxk = getId(ran.tr.Tzx, i, j, k + 1);
        Azetaxx = (ma.zetaxx[idtxxi] + ma.zetaxx[idtxx]) / 2.;
        Azetaxy = (ma.zetaxy[idtxxi] + ma.zetaxy[idtxx]) / 2.;
        Azetaxz = (ma.zetaxz[idtxxi] + ma.zetaxz[idtxx]) / 2.;
        Arho    = (   ma.rho[idtxxi] +    ma.rho[idtxx]) / 2.;
        aft->va.Vxx[idvx] = (2. * Arho - Azetaxx * dif.dt) / (2. * Arho + Azetaxx * dif.dt) * bef->va.Vxx[idvx]
         + 2. * dif.dt / (2. * Arho + Azetaxx * dif.dt) * (bef->sa.Txx[idtxx] - bef->sa.Txx[idtxxi]) / dif.dx;

        aft->va.Vxy[idvx] = (2. * Arho - Azetaxy * dif.dt) / (2. * Arho + Azetaxy * dif.dt) * bef->va.Vxy[idvx]
         + 2. * dif.dt / (2. * Arho + Azetaxy * dif.dt) * (bef->ta.Txy[idtxyj] - bef->ta.Txy[idtxy]) / dif.dy;

        aft->va.Vxz[idvx] = (2. * Arho - Azetaxz * dif.dt) / (2. * Arho + Azetaxz * dif.dt) * bef->va.Vxz[idvx]
         + 2. * dif.dt / (2. * Arho + Azetaxz * dif.dt) * (bef->ta.Tzx[idtzxk] - bef->ta.Tzx[idtzx]) / dif.dz;
      }
    }
  }
//全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Vxkmax; k++) {
    for (j = 0; j < Vxjmax; j++) {
      for (i = 0; i < Vximax; i++) {
        int idvx = getId(ran.vr.Vx, i, j, k);
        aft->va.Vx[idvx] = aft->va.Vxx[idvx] + aft->va.Vxy[idvx] + aft->va.Vxz[idvx];
      }
    }
  }
}

void Vy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran) {
  int i, j, k;
  int Vyimax = ran.vr.Vy.x, Vyjmax = ran.vr.Vy.y, Vykmax = ran.vr.Vy.z;
  double Azetayx, Azetayy, Azetayz, Arho;
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Vykmax; k++) {
    for (j = 1; j < Vyjmax - 1; j++) {
      for (i = 0; i < Vyimax; i++) {
        int idvy   = getId(ran.vr.Vy, i, j, k);
        int idtyy  = getId(ran.sr.Tyy, i, j, k);
        int idtyyj = getId(ran.sr.Tyy, i, j - 1, k);
        int idtxy  = getId(ran.tr.Txy, i, j, k);
        int idtxyi = getId(ran.tr.Txy, i + 1, j, k);
        int idtyz  = getId(ran.tr.Tyz, i, j, k);
        int idtyzk = getId(ran.tr.Tyz, i, j, k + 1);

        Azetayx = (ma.zetayx[idtyyj] + ma.zetayx[idtyy]) / 2.;
        Azetayy = (ma.zetayy[idtyyj] + ma.zetayy[idtyy]) / 2.;
        Azetayz = (ma.zetayz[idtyyj] + ma.zetayz[idtyy]) / 2.;
        Arho    = (   ma.rho[idtyyj] +    ma.rho[idtyy]) / 2.;

        aft->va.Vyx[idvy] = (2. * Arho - Azetayx * dif.dt) / (2. * Arho + Azetayx * dif.dt) * bef->va.Vyx[idvy]
          + 2. * dif.dt / (2. * Arho + Azetayx * dif.dt) * (bef->ta.Txy[idtxyi] - bef->ta.Txy[idtxy]) / dif.dx;

        aft->va.Vyy[idvy] = (2. * Arho - Azetayy * dif.dt) / (2. * Arho + Azetayy * dif.dt) * bef->va.Vyy[idvy]
          + 2. * dif.dt / (2. * Arho + Azetayy * dif.dt) * (bef->sa.Tyy[idtyy] - bef->sa.Tyy[idtyyj]) / dif.dy;

        aft->va.Vyz[idvy] = (2. * Arho - Azetayz * dif.dt) / (2. * Arho + Azetayz * dif.dt) * bef->va.Vyz[idvy]
          + 2. * dif.dt / (2. * Arho + Azetayz * dif.dt) * (bef->ta.Tyz[idtyzk] - bef->ta.Tyz[idtyz]) / dif.dz;
      }
    }
  }

// 全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Vykmax; k++) {
    for (j = 0; j < Vyjmax; j++) {
      for (i = 0; i < Vyimax; i++) {
        int idvy = getId(ran.vr.Vy, i, j, k);
        aft->va.Vy[idvy] = aft->va.Vyx[idvy] + aft->va.Vyy[idvy] + aft->va.Vyz[idvy];
      }
    }
  }
}

void Vz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran) {
  int i, j, k;
  int Vzimax = ran.vr.Vz.x, Vzjmax = ran.vr.Vz.y, Vzkmax = ran.vr.Vz.z;
  double Azetazx, Azetazy, Azetazz, Arho;
#pragma omp parallel for private(i, j, k)
  for (k = 1; k < Vzkmax - 1; k++) {
    for (j = 0; j < Vzjmax; j++) {
      for (i = 0; i < Vzimax; i++) {
        int idvz   = getId(ran.vr.Vz, i, j, k);
        int idtzz  = getId(ran.sr.Tzz, i, j, k);
        int idtzzk = getId(ran.sr.Tzz, i, j, k - 1);
        int idtzx  = getId(ran.tr.Tzx, i, j, k);
        int idtzxi = getId(ran.tr.Tzx, i + 1, j, k);
        int idtyz  = getId(ran.tr.Tyz, i, j, k);
        int idtyzj = getId(ran.tr.Tyz, i, j + 1, k);

        Azetazx = (ma.zetazx[idtzzk] + ma.zetazx[idtzz]) / 2.;
        Azetazy = (ma.zetazy[idtzzk] + ma.zetazy[idtzz]) / 2.;
        Azetazz = (ma.zetazz[idtzzk] + ma.zetazz[idtzz]) / 2.;
        Arho    = (   ma.rho[idtzzk] +    ma.rho[idtzz]) / 2.;

        aft->va.Vzx[idvz] = (2. * Arho - Azetazx * dif.dt) / (2. * Arho + Azetazx * dif.dt) * bef->va.Vzx[idvz]
          + 2. * dif.dt / (2. * Arho + Azetazx * dif.dt) * (bef->ta.Tzx[idtzxi] - bef->ta.Tzx[idtzx]) / dif.dx;

        aft->va.Vzy[idvz] = (2. * Arho - Azetazy * dif.dt) / (2. * Arho + Azetazy * dif.dt) * bef->va.Vzy[idvz]
          + 2. * dif.dt / (2. * Arho + Azetazy * dif.dt) * (bef->ta.Tyz[idtyzj] - bef->ta.Tyz[idtyz]) / dif.dy;

        aft->va.Vzz[idvz] = (2. * Arho - Azetazz * dif.dt) / (2. * Arho + Azetazz * dif.dt) * bef->va.Vzz[idvz]
          + 2. * dif.dt / (2. * Arho + Azetazz * dif.dt) * (bef->sa.Tzz[idtzz] - bef->sa.Tzz[idtzzk]) / dif.dz;
      }
    }
  }

// 全方向加算
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Vzkmax; k++) {
    for (j = 0; j < Vzjmax; j++) {
      for (i = 0; i < Vzimax; i++) {
        int idvz = getId(ran.vr.Vz, i, j, k);
        aft->va.Vz[idvz] = aft->va.Vzx[idvz] + aft->va.Vzy[idvz] + aft->va.Vzz[idvz];
      }
    }
  }

}

//粒子速度計算
void Vel(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran) {
  Vx(aft, bef, ma, dif, ran);
  Vy(aft, bef, ma, dif, ran);
  Vz(aft, bef, ma, dif, ran);
}

void Acceleration(CoordD *Acc,BefAft *aft, BefAft *bef, VelRan vr, Diff dif, Coord out){
  int idvx  = getId(vr.Vx, out.x, out.y, out.z);
  int idvxi = getId(vr.Vx, out.x + 1, out.y, out.z);
  int idvy  = getId(vr.Vy, out.x, out.y, out.z);
  int idvyj = getId(vr.Vy, out.x, out.y + 1, out.z);
  int idvz  = getId(vr.Vz, out.x, out.y, out.z);
  int idvzk = getId(vr.Vz, out.x, out.y, out.z + 1);
  Acc->x = ((aft->va.Vx[idvxi] - bef->va.Vx[idvxi]) / dif.dt  + (aft->va.Vx[idvx] - bef->va.Vx[idvx]) / dif.dt) / 2;
  Acc->y = ((aft->va.Vy[idvyj] - bef->va.Vy[idvyj]) / dif.dt  + (aft->va.Vy[idvy] - bef->va.Vy[idvy]) / dif.dt) / 2;
  Acc->z = ((aft->va.Vz[idvzk] - bef->va.Vz[idvzk]) / dif.dt  + (aft->va.Vz[idvz] - bef->va.Vz[idvz]) / dif.dt) / 2;
}

//更新
void swapBefAft(BefAft *aft, BefAft *bef, Range ran) {
  int i, j, k;
  int Txximax = ran.sr.Txx.x, Txxjmax = ran.sr.Txx.y, Txxkmax = ran.sr.Txx.z;
  int Tyyimax = ran.sr.Tyy.x, Tyyjmax = ran.sr.Tyy.y, Tyykmax = ran.sr.Tyy.z;
  int Tzzimax = ran.sr.Tzz.x, Tzzjmax = ran.sr.Tzz.y, Tzzkmax = ran.sr.Tzz.z;
  int Txyimax = ran.tr.Txy.x, Txyjmax = ran.tr.Txy.y, Txykmax = ran.tr.Txy.z;
  int Tyzimax = ran.tr.Tyz.x, Tyzjmax = ran.tr.Tyz.y, Tyzkmax = ran.tr.Tyz.z;
  int Tzximax = ran.tr.Tzx.x, Tzxjmax = ran.tr.Tzx.y, Tzxkmax = ran.tr.Tzx.z;
  int Vximax = ran.vr.Vx.x, Vxjmax = ran.vr.Vx.y, Vxkmax = ran.vr.Vx.z;
  int Vyimax = ran.vr.Vy.x, Vyjmax = ran.vr.Vy.y, Vykmax = ran.vr.Vy.z;
  int Vzimax = ran.vr.Vz.x, Vzjmax = ran.vr.Vz.y, Vzkmax = ran.vr.Vz.z;
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        int idtxx = getId(ran.sr.Txx, i, j, k);
        bef->sa.Txx[idtxx] = aft->sa.Txx[idtxx];
        bef->sa.Txxx[idtxx] = aft->sa.Txxx[idtxx];
        bef->sa.Txxy[idtxx] = aft->sa.Txxy[idtxx];
        bef->sa.Txxz[idtxx] = aft->sa.Txxz[idtxx];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Tyykmax; k++) {
    for (j = 0; j < Tyyjmax; j++) {
      for (i = 0; i < Tyyimax; i++) {
        int idtyy = getId(ran.sr.Tyy, i, j, k);
        bef->sa.Tyy[idtyy] = aft->sa.Tyy[idtyy];
        bef->sa.Tyyx[idtyy] = aft->sa.Tyyx[idtyy];
        bef->sa.Tyyy[idtyy] = aft->sa.Tyyy[idtyy];
        bef->sa.Tyyz[idtyy] = aft->sa.Tyyz[idtyy];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Tzzkmax; k++) {
    for (j = 0; j < Tzzjmax; j++) {
      for (i = 0; i < Tzzimax; i++) {
        int idtzz = getId(ran.sr.Tzz, i, j, k);
        bef->sa.Tzz[idtzz] = aft->sa.Tzz[idtzz];
        bef->sa.Tzzx[idtzz] = aft->sa.Tzzx[idtzz];
        bef->sa.Tzzy[idtzz] = aft->sa.Tzzy[idtzz];
        bef->sa.Tzzz[idtzz] = aft->sa.Tzzz[idtzz];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (i = 0; i < Txyimax; i++) {
    for (j = 0; j < Txyjmax; j++) {
      for (k = 0; k < Txykmax; k++) {
        int idtxy = getId(ran.tr.Txy, i, j, k);
        bef->ta.Txy[idtxy] = aft->ta.Txy[idtxy];
        bef->ta.Txyx[idtxy] = aft->ta.Txyx[idtxy];
        bef->ta.Txyy[idtxy] = aft->ta.Txyy[idtxy];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (i = 0; i < Tyzimax; i++) {
    for (j = 0; j < Tyzjmax; j++) {
      for (k = 0; k < Tyzkmax; k++) {
        int idtyz = getId(ran.tr.Tyz, i, j, k);
        bef->ta.Tyz[idtyz] = aft->ta.Tyz[idtyz];
        bef->ta.Tyzy[idtyz] = aft->ta.Tyzy[idtyz];
        bef->ta.Tyzz[idtyz] = aft->ta.Tyzz[idtyz];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (i = 0; i < Tzximax; i++) {
    for (j = 0; j < Tzxjmax; j++) {
      for (k = 0; k < Tzxkmax; k++) {
        int idtzx = getId(ran.tr.Tzx, i, j, k);
        bef->ta.Tzx[idtzx] = aft->ta.Tzx[idtzx];
        bef->ta.Tzxz[idtzx] = aft->ta.Tzxz[idtzx];
        bef->ta.Tzxx[idtzx] = aft->ta.Tzxx[idtzx];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Vxkmax; k++) {
    for (j = 0; j < Vxjmax; j++) {
      for (i = 0; i < Vximax; i++) {
        int idvx = getId(ran.vr.Vx, i, j, k);
        bef->va.Vx[idvx] = aft->va.Vx[idvx];
        bef->va.Vxx[idvx] = aft->va.Vxx[idvx];
        bef->va.Vxy[idvx] = aft->va.Vxy[idvx];
        bef->va.Vxz[idvx] = aft->va.Vxz[idvx];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Vykmax; k++) {
    for (j = 0; j < Vyjmax; j++) {
      for (i = 0; i < Vyimax; i++) {
        int idvy = getId(ran.vr.Vy, i, j, k);
        bef->va.Vy[idvy] = aft->va.Vy[idvy];
        bef->va.Vyx[idvy] = aft->va.Vyx[idvy];
        bef->va.Vyy[idvy] = aft->va.Vyy[idvy];
        bef->va.Vyz[idvy] = aft->va.Vyz[idvy];
      }
    }
  }
#pragma omp parallel for private(i, j, k)
  for (k = 0; k < Vzkmax; k++) {
    for (j = 0; j < Vzjmax; j++) {
      for (i = 0; i < Vzimax; i++) {
        int idvz = getId(ran.vr.Vz, i, j, k);
        bef->va.Vz[idvz] = aft->va.Vz[idvz];
        bef->va.Vzx[idvz] = aft->va.Vzx[idvz];
        bef->va.Vzy[idvz] = aft->va.Vzy[idvz];
        bef->va.Vzz[idvz] = aft->va.Vzz[idvz];
      }
    }
  }
}