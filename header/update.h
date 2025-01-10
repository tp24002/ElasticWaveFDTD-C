#pragma once
#include "./struct.h"

void Vx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran);
void Vy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran);
void Vz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran);
void Vel(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran);
void Txx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran, Inpaluse ip, int t);
void Tyy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran, Inpaluse ip, int t);
void Tzz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran, Inpaluse ip, int t);
void Sig(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran, Inpaluse ip, int t);
void Txy(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran);
void Tyz(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran);
void Tzx(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran);
void Tau(BefAft *aft, BefAft *bef, MedArr ma, Diff dif, Range ran);
void Acceleration(CoordD *Acc,BefAft *aft, BefAft *bef, VelRan vr, Diff dif, Coord out);
void swapBefAft(BefAft *aft, BefAft *bef, Range ran);