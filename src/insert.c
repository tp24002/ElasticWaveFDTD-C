#include "../header/insert.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../header/struct.h"

void insertAir(MedArr *ma, SigRan sr, Medium air) {
  int i, j, k;
  for (k = 0; k < sr.Txx.z; k++) {
    for (j = 0; j < sr.Txx.y; j++) {
      for (i = 0; i < sr.Txx.x; i++) {
        ma->ramda[i][j][k] = air.ramda;
        ma->mu[i][j][k] = air.G;
        ma->c11[i][j][k] = air.ramda + 2. * air.G;
        ma->rho[i][j][k] = air.rho;
        ma->zetaxx[i][j][k] = air.zeta;
        ma->zetaxy[i][j][k] = air.zeta;
        ma->zetaxz[i][j][k] = air.zeta;
        ma->zetayx[i][j][k] = air.zeta;
        ma->zetayy[i][j][k] = air.zeta;
        ma->zetayz[i][j][k] = air.zeta;
        ma->zetazx[i][j][k] = air.zeta;
        ma->zetazy[i][j][k] = air.zeta;
        ma->zetazz[i][j][k] = air.zeta;
        ma->gamma[i][j][k] = air.gamma;
        ma->khi[i][j][k] = air.khi;
        ma->xi11[i][j][k] = air.khi + 2. * air.gamma;
        ma->zetadx[i][j][k] = 0.;
        ma->zetady[i][j][k] = 0.;
        ma->zetadz[i][j][k] = 0.;
      }
    }
  }
}

void insertObject(MedArr *ma, Object ob) {
  int i, j, k;
  Medium objmed = ob.med;
  for (k = ob.sp.z; k < ob.sp.z + ob.range.z; k++) {
    for (j = ob.sp.y; j < ob.sp.y + ob.range.y; j++) {
      for (i = ob.sp.x; i < ob.sp.x + ob.range.x; i++) {
        ma->ramda[i][j][k] = objmed.ramda;
        ma->mu[i][j][k] = objmed.G;
        ma->c11[i][j][k] = objmed.ramda + 2. * objmed.G;
        ma->rho[i][j][k] = objmed.rho;
        ma->zetaxx[i][j][k] = objmed.zeta;
        ma->zetaxy[i][j][k] = objmed.zeta;
        ma->zetaxz[i][j][k] = objmed.zeta;
        ma->zetayx[i][j][k] = objmed.zeta;
        ma->zetayy[i][j][k] = objmed.zeta;
        ma->zetayz[i][j][k] = objmed.zeta;
        ma->zetazx[i][j][k] = objmed.zeta;
        ma->zetazy[i][j][k] = objmed.zeta;
        ma->zetazz[i][j][k] = objmed.zeta;
        ma->gamma[i][j][k] = objmed.gamma;
        ma->khi[i][j][k] = objmed.khi;
        ma->xi11[i][j][k] = objmed.khi + 2. * objmed.gamma;
      }
    }
  }
}

void insertPml(MedArr *ma, SigRan sr, Pml pml) {
  int plx1 = pml.pl1.x, plx2 = pml.pl2.x;
  int ply1 = pml.pl1.y, ply2 = pml.pl2.y;
  int plz1 = pml.pl1.z, plz2 = pml.pl2.z;
  int Txximax = sr.Txx.x, Txxjmax = sr.Txx.y, Txxkmax = sr.Txx.z;
  double zeta_max = pml.fm, ta = pml.ta;
  int i, j, k;
  //x方向
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < plx1; i++) {
        ma->zetaxx[i][j][k] += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        ma->zetayx[i][j][k] += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        ma->zetazx[i][j][k] += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        ma->zetadx[i][j][k] = ma->zetaxx[i][j][k] / ma->rho[i][j][k];
      }
    }
  }
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = Txximax - 1; i > Txximax - 1 - plx2; i--) {
        ma->zetaxx[i][j][k] += zeta_max * pow(((double)plx2 - (double)(Txximax - i)) / (double)plx2, ta);
        ma->zetayx[i][j][k] += zeta_max * pow(((double)plx2 - (double)(Txximax - i)) / (double)plx2, ta);
        ma->zetazx[i][j][k] += zeta_max * pow(((double)plx2 - (double)(Txximax - i)) / (double)plx2, ta);
        ma->zetadx[i][j][k] = ma->zetaxx[i][j][k] / ma->rho[i][j][k];
      }
    }
  }
  //y方向
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < ply1; j++) {
      for (i = 0; i < Txximax; i++) {
        ma->zetaxy[i][j][k] += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        ma->zetayy[i][j][k] += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        ma->zetazy[i][j][k] += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        ma->zetady[i][j][k] = ma->zetaxy[i][j][k] / ma->rho[i][j][k];
      }
    }
  }
  for (k = 0; k < Txxkmax; k++) {
    for (j = Txxjmax - 1; j > Txxjmax - 1 - ply2; j--) {
      for (i = 0; i < Txximax; i++) {
        ma->zetaxy[i][j][k] += zeta_max * pow(((double)ply2 - (double)(Txxjmax - j)) / (double)ply2, ta);
        ma->zetayy[i][j][k] += zeta_max * pow(((double)ply2 - (double)(Txxjmax - j)) / (double)ply2, ta);
        ma->zetazy[i][j][k] += zeta_max * pow(((double)ply2 - (double)(Txxjmax - j)) / (double)ply2, ta);
        ma->zetady[i][j][k] = ma->zetaxy[i][j][k] / ma->rho[i][j][k];
      }
    }
  }
  //z方向
  for (k = 0; k < plz1; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        ma->zetaxz[i][j][k] += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        ma->zetayz[i][j][k] += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        ma->zetazz[i][j][k] += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        ma->zetadz[i][j][k] = ma->zetaxz[i][j][k] / ma->rho[i][j][k];
      }
    }
  }
  for (k = Txxkmax - 1; k > Txxkmax - 1 - plz2; k--) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        ma->zetaxz[i][j][k] += zeta_max * pow(((double)plz2 - (double)(Txxkmax - k)) / (double)plz2, ta);
        ma->zetayz[i][j][k] += zeta_max * pow(((double)plz2 - (double)(Txxkmax - k)) / (double)plz2, ta);
        ma->zetazz[i][j][k] += zeta_max * pow(((double)plz2 - (double)(Txxkmax - k)) / (double)plz2, ta);
        ma->zetadz[i][j][k] = ma->zetaxz[i][j][k] / ma->rho[i][j][k];
      }
    }
  }
}

void zeroPadSig(SigArr *sa, SigRan sr) {
  int i, j, k;
  for (k = 0; k < sr.Txx.z; k++) {
    for (j = 0; j < sr.Txx.y; j++) {
      for (i = 0; i < sr.Txx.x; i++) {
        sa->Txx[i][j][k] = 0.;
        sa->Txxx[i][j][k] = 0.;
        sa->Txxy[i][j][k] = 0.;
        sa->Txxz[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < sr.Tyy.z; k++) {
    for (j = 0; j < sr.Tyy.y; j++) {
      for (i = 0; i < sr.Tyy.x; i++) {
        sa->Tyy[i][j][k] = 0.;
        sa->Tyyx[i][j][k] = 0.;
        sa->Tyyy[i][j][k] = 0.;
        sa->Tyyz[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < sr.Tzz.z; k++) {
    for (j = 0; j < sr.Tzz.y; j++) {
      for (i = 0; i < sr.Tzz.x; i++) {
        sa->Tzz[i][j][k] = 0.;
        sa->Tzzx[i][j][k] = 0.;
        sa->Tzzy[i][j][k] = 0.;
        sa->Tzzz[i][j][k] = 0.;
      }
    }
  }
}

void zeroPadTau(TauArr *ta, TauRan tr) {
  int i, j, k;
  for (k = 0; k < tr.Txy.z; k++) {
    for (j = 0; j < tr.Txy.y; j++) {
      for (i = 0; i < tr.Txy.x; i++) {
        ta->Txy[i][j][k] = 0.;
        ta->Txyx[i][j][k] = 0.;
        ta->Txyy[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < tr.Tyz.z; k++) {
    for (j = 0; j < tr.Tyz.y; j++) {
      for (i = 0; i < tr.Tyz.x; i++) {
        ta->Tyz[i][j][k] = 0.;
        ta->Tyzy[i][j][k] = 0.;
        ta->Tyzz[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < tr.Tzx.z; k++) {
    for (j = 0; j < tr.Tzx.y; j++) {
      for (i = 0; i < tr.Tzx.x; i++) {
        ta->Tzx[i][j][k] = 0.;
        ta->Tzxz[i][j][k] = 0.;
        ta->Tzxx[i][j][k] = 0.;
      }
    }
  }
}

void zeroPadVel(VelArr *va, VelRan vr) {
  int i, j, k;
  for (k = 0; k < vr.Vx.z; k++) {
    for (j = 0; j < vr.Vx.y; j++) {
      for (i = 0; i < vr.Vx.x; i++) {
        va->Vx[i][j][k] = 0.;
        va->Vxx[i][j][k] = 0.;
        va->Vxy[i][j][k] = 0.;
        va->Vxz[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k <= vr.Vy.z; k++) {
    for (j = 0; j < vr.Vy.y; j++) {
      for (i = 0; i < vr.Vy.x; i++) {
        va->Vy[i][j][k] = 0.;
        va->Vyx[i][j][k] = 0.;
        va->Vyy[i][j][k] = 0.;
        va->Vyz[i][j][k] = 0.;
      }
    }
  }
  for (k = 0; k < vr.Vz.z; k++) {
    for (j = 0; j < vr.Vz.y; j++) {
      for (i = 0; i < vr.Vz.x; i++) {
        va->Vz[i][j][k] = 0.;
        va->Vzx[i][j][k] = 0.;
        va->Vzy[i][j][k] = 0.;
        va->Vzz[i][j][k] = 0.;
      }
    }
  }
}

void zeroPadding(BefAft *ba, Range ran) {
  zeroPadSig(&ba->sa, ran.sr);
  zeroPadTau(&ba->ta, ran.tr);
  zeroPadVel(&ba->va, ran.vr);
}