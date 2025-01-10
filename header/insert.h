#pragma once
#include "./struct.h"

void insertAir(MedArr *ma, SigRan sr, Medium air);
void insertObject(MedArr *ma, SigRan sr, Object ob);
void insertPml(MedArr *ma, SigRan sr, Pml pml);
void zeroPadSig(SigArr *sa, SigRan sr);
void zeroPadTau(TauArr *ta, TauRan tr);
void zeroPadVel(VelArr *va, VelRan vr);
void zeroPadding(BefAft *ba, Range ran);
