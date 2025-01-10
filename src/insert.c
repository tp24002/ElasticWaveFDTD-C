#include "../header/insert.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../header/struct.h"
#include "../header/init.h"

void insertAir(MedArr *ma, SigRan sr, Medium air) {
  int i, j, k;
  int id;
  for (k = 0; k < sr.Txx.z; k++) {
    for (j = 0; j < sr.Txx.y; j++) {
      for (i = 0; i < sr.Txx.x; i++) {
        id = getId(sr.Txx, i, j, k);
        ma->ramda[id] = air.ramda;
        ma->mu[id] = air.G;
        ma->c11[id] = air.ramda + 2. * air.G;
        ma->rho[id] = air.rho;
        ma->zetaxx[id] = air.zeta;
        ma->zetaxy[id] = air.zeta;
        ma->zetaxz[id] = air.zeta;
        ma->zetayx[id] = air.zeta;
        ma->zetayy[id] = air.zeta;
        ma->zetayz[id] = air.zeta;
        ma->zetazx[id] = air.zeta;
        ma->zetazy[id] = air.zeta;
        ma->zetazz[id] = air.zeta;
        ma->gamma[id] = air.gamma;
        ma->khi[id] = air.khi;
        ma->xi11[id] = air.khi + 2. * air.gamma;
        ma->zetadx[id] = 0.;
        ma->zetady[id] = 0.;
        ma->zetadz[id] = 0.;
      }
    }
  }
}

void insertObject(MedArr *ma, SigRan sr, Object ob) {
  int i, j, k;
  int id;
  Medium objmed = ob.med;
  for (k = ob.sp.z; k < ob.sp.z + ob.range.z; k++) {
    for (j = ob.sp.y; j < ob.sp.y + ob.range.y; j++) {
      for (i = ob.sp.x; i < ob.sp.x + ob.range.x; i++) {
        id = getId(sr.Txx, i, j, k);
        ma->ramda[id] = objmed.ramda;
        ma->mu[id] = objmed.G;
        ma->c11[id] = objmed.ramda + 2. * objmed.G;
        ma->rho[id] = objmed.rho;
        ma->zetaxx[id] = objmed.zeta;
        ma->zetaxy[id] = objmed.zeta;
        ma->zetaxz[id] = objmed.zeta;
        ma->zetayx[id] = objmed.zeta;
        ma->zetayy[id] = objmed.zeta;
        ma->zetayz[id] = objmed.zeta;
        ma->zetazx[id] = objmed.zeta;
        ma->zetazy[id] = objmed.zeta;
        ma->zetazz[id] = objmed.zeta;
        ma->gamma[id] = objmed.gamma;
        ma->khi[id] = objmed.khi;
        ma->xi11[id] = objmed.khi + 2. * objmed.gamma;
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
  int id;
  //x方向
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < plx1; i++) {
        id = getId(sr.Txx, i, j, k);
        ma->zetaxx[id] += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        ma->zetayx[id] += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        ma->zetazx[id] += zeta_max * pow(((double)plx1 - (double)i) / (double)plx1, ta);
        ma->zetadx[id] = ma->zetaxx[id] / ma->rho[id];
      }
    }
  }
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = Txximax - 1; i > Txximax - 1 - plx2; i--) {
        id = getId(sr.Txx, i, j, k);
        ma->zetaxx[id] += zeta_max * pow(((double)plx2 - (double)(Txximax - i)) / (double)plx2, ta);
        ma->zetayx[id] += zeta_max * pow(((double)plx2 - (double)(Txximax - i)) / (double)plx2, ta);
        ma->zetazx[id] += zeta_max * pow(((double)plx2 - (double)(Txximax - i)) / (double)plx2, ta);
        ma->zetadx[id] = ma->zetaxx[id] / ma->rho[id];
      }
    }
  }
  //y方向
  for (k = 0; k < Txxkmax; k++) {
    for (j = 0; j < ply1; j++) {
      for (i = 0; i < Txximax; i++) {
        id = getId(sr.Txx, i, j, k);
        ma->zetaxy[id] += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        ma->zetayy[id] += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        ma->zetazy[id] += zeta_max * pow(((double)ply1 - (double)j) / (double)ply1, ta);
        ma->zetady[id] = ma->zetaxy[id] / ma->rho[id];
      }
    }
  }
  for (k = 0; k < Txxkmax; k++) {
    for (j = Txxjmax - 1; j > Txxjmax - 1 - ply2; j--) {
      for (i = 0; i < Txximax; i++) {
        id = getId(sr.Txx, i, j, k);
        ma->zetaxy[id] += zeta_max * pow(((double)ply2 - (double)(Txxjmax - j)) / (double)ply2, ta);
        ma->zetayy[id] += zeta_max * pow(((double)ply2 - (double)(Txxjmax - j)) / (double)ply2, ta);
        ma->zetazy[id] += zeta_max * pow(((double)ply2 - (double)(Txxjmax - j)) / (double)ply2, ta);
        ma->zetady[id] = ma->zetaxy[id] / ma->rho[id];
      }
    }
  }
  //z方向
  for (k = 0; k < plz1; k++) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        id = getId(sr.Txx, i, j, k);
        ma->zetaxz[id] += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        ma->zetayz[id] += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        ma->zetazz[id] += zeta_max * pow(((double)plz1 - (double)k) / (double)plz1, ta);
        ma->zetadz[id] = ma->zetaxz[id] / ma->rho[id];
      }
    }
  }
  for (k = Txxkmax - 1; k > Txxkmax - 1 - plz2; k--) {
    for (j = 0; j < Txxjmax; j++) {
      for (i = 0; i < Txximax; i++) {
        id = getId(sr.Txx, i, j, k);
        ma->zetaxz[id] += zeta_max * pow(((double)plz2 - (double)(Txxkmax - k)) / (double)plz2, ta);
        ma->zetayz[id] += zeta_max * pow(((double)plz2 - (double)(Txxkmax - k)) / (double)plz2, ta);
        ma->zetazz[id] += zeta_max * pow(((double)plz2 - (double)(Txxkmax - k)) / (double)plz2, ta);
        ma->zetadz[id] = ma->zetaxz[id] / ma->rho[id];
      }
    }
  }
}

void zeroPadSig(SigArr *sa, SigRan sr) {
  int i, j, k;
  int id;
  for (k = 0; k < sr.Txx.z; k++) {
    for (j = 0; j < sr.Txx.y; j++) {
      for (i = 0; i < sr.Txx.x; i++) {
        id = getId(sr.Txx, i, j, k);
        sa->Txx[id] = 0.;
        sa->Txxx[id] = 0.;
        sa->Txxy[id] = 0.;
        sa->Txxz[id] = 0.;
      }
    }
  }
  for (k = 0; k < sr.Tyy.z; k++) {
    for (j = 0; j < sr.Tyy.y; j++) {
      for (i = 0; i < sr.Tyy.x; i++) {
        id = getId(sr.Tyy, i, j, k);
        sa->Tyy[id] = 0.;
        sa->Tyyx[id] = 0.;
        sa->Tyyy[id] = 0.;
        sa->Tyyz[id] = 0.;
      }
    }
  }
  for (k = 0; k < sr.Tzz.z; k++) {
    for (j = 0; j < sr.Tzz.y; j++) {
      for (i = 0; i < sr.Tzz.x; i++) {
        id = getId(sr.Tzz, i, j, k);
        sa->Tzz[id] = 0.;
        sa->Tzzx[id] = 0.;
        sa->Tzzy[id] = 0.;
        sa->Tzzz[id] = 0.;
      }
    }
  }
}

void zeroPadTau(TauArr *ta, TauRan tr) {
  int i, j, k;
  int id;
  for (k = 0; k < tr.Txy.z; k++) {
    for (j = 0; j < tr.Txy.y; j++) {
      for (i = 0; i < tr.Txy.x; i++) {
        id = getId(tr.Txy, i, j, k);
        ta->Txy[id] = 0.;
        ta->Txyx[id] = 0.;
        ta->Txyy[id] = 0.;
      }
    }
  }
  for (k = 0; k < tr.Tyz.z; k++) {
    for (j = 0; j < tr.Tyz.y; j++) {
      for (i = 0; i < tr.Tyz.x; i++) {
        id = getId(tr.Tyz, i, j, k);
        ta->Tyz[id] = 0.;
        ta->Tyzy[id] = 0.;
        ta->Tyzz[id] = 0.;
      }
    }
  }
  for (k = 0; k < tr.Tzx.z; k++) {
    for (j = 0; j < tr.Tzx.y; j++) {
      for (i = 0; i < tr.Tzx.x; i++) {
        id = getId(tr.Tzx, i, j, k);
        ta->Tzx[id] = 0.;
        ta->Tzxz[id] = 0.;
        ta->Tzxx[id] = 0.;
      }
    }
  }
}

void zeroPadVel(VelArr *va, VelRan vr) {
  int i, j, k;
  int id;
  for (k = 0; k < vr.Vx.z; k++) {
    for (j = 0; j < vr.Vx.y; j++) {
      for (i = 0; i < vr.Vx.x; i++) {
        id = getId(vr.Vx, i, j, k);
        va->Vx[id] = 0.;
        va->Vxx[id] = 0.;
        va->Vxy[id] = 0.;
        va->Vxz[id] = 0.;
      }
    }
  }
  for (k = 0; k <= vr.Vy.z; k++) {
    for (j = 0; j < vr.Vy.y; j++) {
      for (i = 0; i < vr.Vy.x; i++) {
        id = getId(vr.Vy, i, j, k);
        va->Vy[id] = 0.;
        va->Vyx[id] = 0.;
        va->Vyy[id] = 0.;
        va->Vyz[id] = 0.;
      }
    }
  }
  for (k = 0; k < vr.Vz.z; k++) {
    for (j = 0; j < vr.Vz.y; j++) {
      for (i = 0; i < vr.Vz.x; i++) {
        id = getId(vr.Vz, i, j, k);
        va->Vz[id] = 0.;
        va->Vzx[id] = 0.;
        va->Vzy[id] = 0.;
        va->Vzz[id] = 0.;
      }
    }
  }
}

void zeroPadding(BefAft *ba, Range ran) {
  zeroPadSig(&ba->sa, ran.sr);
  zeroPadTau(&ba->ta, ran.tr);
  zeroPadVel(&ba->va, ran.vr);
}