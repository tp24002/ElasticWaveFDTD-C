#pragma once
#include "./struct.h"

void Vx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
void Vy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
void Vz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
void Vel(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, VelRan vr);
void Txx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t);
void Tyy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t);
void Tzz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t);
void Sig(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, SigRan sr, Inpaluse ip, int t);
void Txy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
void Tyz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
void Tzx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
void Tau(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, TauRan tr);
void Acceleration(CoordD *Acc,BefAft *aft, BefAft *bef, Diff dif, Coord out);
void swapBefAft(BefAft *aft, BefAft *bef, Range ran);