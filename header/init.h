#pragma once
#include "./struct.h"

void initMedium(Medium *med);
void initCoord(Coord *co, int x, int y, int z);
void initDiff(Diff *dif, Medium *med);
void initPml(Pml *pml, Medium *med, Diff dif);
void initObject(Object *con, Medium med, Pml pml, Coord sp, Coord ran);
void initRange(Range *ran, Coord region, Pml pml);
void initSigArr(SigArr *sa, SigRan sr);
void initTauArr(TauArr *ta, TauRan tr);
void initVelArr(VelArr *va, VelRan vr);
void initBefAft(BefAft *ba, Range ran);
void initInpalse(Inpaluse *ip, SigRan sr, Pml pml);
void initMedArr(MedArr *ma, SigRan sr);
void initClack(Object *clack, Medium med, Pml *pml, int spx, int spy, int spz, int ranx, int rany, int ranz);