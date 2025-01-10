#pragma once

typedef struct {
  int x;
  int y;
  int z;
} Coord;

typedef struct {
  double x;
  double y;
  double z;
} CoordD;

typedef struct {
  double rho;
  double K;
  double E;
  double G;
  double nu;
  double ramda;
  double khi;
  double gamma;
  double zeta;
  double eta;   // loss factor
  double omega; // first-order natural angular frequency
} Medium;

typedef struct {
  Medium med;
  Coord range;
  Coord sp;
} Object;

typedef struct {
  Coord Txx;
  Coord Tyy;
  Coord Tzz;
} SigRan;

typedef struct {
  Coord Txy;
  Coord Tyz;
  Coord Tzx;
} TauRan;

typedef struct {
  Coord Vx;
  Coord Vy;
  Coord Vz;
} VelRan;

typedef struct {
  SigRan sr;
  TauRan tr;
  VelRan vr;
} Range;

typedef struct {
  Coord pl1;
  Coord pl2;
  double ta;
  double fm;
} Pml;

typedef struct {
  double dx;
  double dy;
  double dz;
  double dt;
} Diff;

typedef struct {
  double *ramda;
  double *mu;
  double *c11;
  double *rho;
  double *zetaxx;
  double *zetaxy;
  double *zetaxz;
  double *zetayx;
  double *zetayy;
  double *zetayz;
  double *zetazx;
  double *zetazy;
  double *zetazz;
  double *gamma;
  double *khi;
  double *xi11;
  double *zetadx;
  double *zetady;
  double *zetadz;
} MedArr;

typedef struct {
  double *Txx;
  double *Txxx;
  double *Txxy;
  double *Txxz;
  double *Tyy;
  double *Tyyx;
  double *Tyyy;
  double *Tyyz;
  double *Tzz;
  double *Tzzx;
  double *Tzzy;
  double *Tzzz;
} SigArr;

typedef struct {
  double *Txy;
  double *Txyx;
  double *Txyy;
  double *Tyz;
  double *Tyzy;
  double *Tyzz;
  double *Tzx;
  double *Tzxz;
  double *Tzxx;
} TauArr;

typedef struct {
  double *Vx;
  double *Vxx;
  double *Vxy;
  double *Vxz;
  double *Vy;
  double *Vyx;
  double *Vyy;
  double *Vyz;
  double *Vz;
  double *Vzx;
  double *Vzy;
  double *Vzz;
} VelArr;

typedef struct {
  SigArr sa;
  TauArr ta;
  VelArr va;
  SigArr noin;
} BefAft;

typedef struct {
  double *Txx;
  double *Tyy;
  double *Tzz;
  double freq;
  int mode;
  Coord in;
} Inpaluse;

typedef enum { E_AIR = 0, E_CON, E_STEEL, E_M_END } E_Mednum;
typedef enum { E_SINE = 100, E_RCOS } E_KIND_OF_IP;
typedef enum { E_HOLE = 200, E_CLACK } E_KIND_OF_DEFECT;