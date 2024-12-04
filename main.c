#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "./header/init.h"
#include "./header/insert.h"
#include "./header/struct.h"
#include "./header/update.h"
#include "./header/parameter.h"

#define clacknum 1
#define outnum 4
#define tmax 8192

void progressBar(int now, int max);

int main(void) {
  Medium med[E_M_END];
  // Object air; // 領域全体に空気を設定するため今は必要なし
  Object con;
  Object *clack;
  Range ran;
  Pml pml;
  Diff dif;
  MedArr ma;
  BefAft bef, aft;
  Inpaluse ip;
  FILE *fp1;
  char fn1[256];
  Coord *out;
  CoordD *acc;

  clack = (Object*)malloc(clacknum * sizeof(Object));
  out   = (Coord*)malloc(outnum * sizeof(Coord));
  acc   = (CoordD*)malloc(outnum * sizeof(CoordD));

  initMedium(med);

  parameter(&ran,&dif,&pml,med,&con,clack,out,&ip);
  // //ここをfor文で
  // initCoord(&clack_st, clack_st.x, clack_st.y, clack_st.z);
  initBefAft(&bef, ran);
  initBefAft(&aft, ran);
  initMedArr(&ma, ran.sr);
  initInpalse(&ip, ran.sr, pml);

  //出力
  printf("range:%d,%d,%d\n", ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
  printf("in:%d,%d,%d\n", ip.in.x, ip.in.y, ip.in.z);
  printf("med[E_CON].gamma = %le\n", med[E_CON].gamma);
  printf("med[E_CON].khi = %le\n", med[E_CON].khi);
  printf("dif.dt = %le\n", dif.dt);
  printf("pml.fm = %le\n", pml.fm);
  printf("start:%d,%d,%d\n",con.sp.x,con.sp.y,con.sp.z);
  printf("size:%d,%d,%d\n",con.range.x,con.range.y,con.range.z);
  printf("in :%d,%d,%d\n", ip.in.x, ip.in.y, ip.in.z);
  for(int i = 0; i < outnum; i++) {
    printf("out[%d]:%d,%d,%d\n", i, out[i].x,out[i].y,out[i].z);
  }
  if(ip.mode == E_SINE){
    printf("sin:%f\n", ip.freq);
  } else if(ip.mode == E_RCOS){
    printf("cos:%f\n", ip.freq);
  }

  // 媒質配置
  insertAir(&ma, ran.sr, med[E_AIR]);
  insertObject(&ma, con);
  for(int i = 0; i < clacknum; i++) {
    insertObject(&ma, clack[i]);
    printf("clack[%d]:%d,%d,%d\n", i, clack->sp.x, clack->sp.y, clack->sp.z);
  }

  insertPml(&ma, ran.sr, pml);
  zeroPadding(&bef, ran);
  zeroPadding(&aft, ran);

  //ファイル名出力
  sprintf(fn1, "./concrete_input.csv");
  printf("%.*s\n", (int)sizeof fn1, fn1);
  
  fp1 = fopen(fn1, "wb");

  for (int t = 0; t < tmax; t++) {
    Vel(&aft, &bef, ma, dif, ran.vr);
    Sig(&aft, &bef, ma, dif, ran.sr, ip, t);
    Tau(&aft, &bef, ma, dif, ran.tr);
    // 加速度算出＆書き込み
    for(int i = 0; i < outnum; i++){
      Acceleration(&acc[i],&aft, &bef, dif, out[i]);
      fprintf(fp1, "%le,%le,%le," , acc[i].x,acc[i].y,acc[i].z);
    }
    fprintf(fp1,"\n");

    swapBefAft(&aft, &bef, ran);
    progressBar(t, tmax);
  }
  fclose(fp1);
  printf("loop end.\n");
  return 0;
}

void progressBar(int now, int max) {
  int bar_width = 50;
  double progress = (double)(now + 1) / (double)max;
  int bar_length = (int)(progress * bar_width);
  printf("Progress: [");
  for (int j = 0; j < bar_length; j++) {
    printf("=");
  }
  for (int j = bar_length; j < bar_width; j++) {
    printf(" ");
  }
  printf("] %.2f%%\r", progress * 100);
  fflush(stdout);
}
