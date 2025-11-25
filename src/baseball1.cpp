#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

struct Params {
  double g;   // acceleration [m/s^2]
  double m;   // mass [kg]
  double d;   // diameter of ball [m]
  double b;   // linear drag coefficient
  double c;   // quadratic drag coefficient
};

// --- ODE system ---
// y[0] = x, y[1] = vx
// y[2] = z, y[3] = vz

double f_rx(double t, const vector<double> &y, void *params=0){
  return y[1];
}

double f_rz(double t, const vector<double> &y, void *params=0){
  return y[3];
}

double f_vx(double t, const vector<double> &y, void *params=0){
  Params *p = (Params*)params;

  double vx = y[1];
  double vz = y[3];
  double v  = sqrt(vx*vx + vz*vz);

  double Fdrag = p->b * v + p->c * v * v;

  return -(Fdrag * vx / v) / p->m;
}

double f_vz(double t, const vector<double> &y, void *params=0){
  Params *p = (Params*)params;

  double vx = y[1];
  double vz = y[3];
  double v  = sqrt(vx*vx + vz*vz);

  double Fdrag = p->b * v + p->c * v * v;

  return -(Fdrag * vz / v) / p->m - p->g;
}


double f_stop(double t, const vector<double> &y, void *params=0){
  double x = y[0];
  if (x >= 18.5) return 1;
  return 0;
}


double simulatePitch(double v0, double theta, Params *p){

  vector<pfunc_t> funcs(4);
  funcs[0] = f_rx;
  funcs[1] = f_vx;
  funcs[2] = f_rz;
  funcs[3] = f_vz;

  vector<double> y(4);
  y[0] = 0.0;
  y[1] = v0 * cos(theta);
  y[2] = 1.4;
  y[3] = v0 * sin(theta);

  double t0 = 0;
  double tmax = 3.0;
  int nsteps = 1000;

  auto sol = RK4SolveN(funcs, y, nsteps, t0, tmax, (void*)p, f_stop);

  double x,z;
  sol[2].GetPoint(sol[2].GetN()-1, x, z);

  return z;
}

int main(int argc, char **argv){

  Params pars;
  pars.g = 9.81;
  pars.m = 0.145;
  pars.d = 0.075;   

  //air resistance coefs
  pars.b = 1.6e-4 * pars.d;
  pars.c = 0.25    * pars.d * pars.d;

  double xend = 18.5;
  double z0 = 1.4;
  double zend = 0.9;
  double theta0 = 1.0;

  bool showPlot = false; //i don't see any popups in vscode so no, dont show me the plot here lol

  int c;
  while ((c = getopt(argc, argv, "x:z:t:p")) != -1)
    switch(c){
      case 'x': xend = atof(optarg); break;
      case 'z': z0 = atof(optarg); break;
      case 't': theta0 = atof(optarg); break;
      case 'p': showPlot = true; break;
    }

  theta0 *= M_PI/180.0;

  TApplication theApp("App", &argc, argv);

  double vLow = 10.0;
  double vHigh = 60.0;
  double vMid = 0;
  double zFinal;

  for(int i=0; i<40; i++){
    vMid = 0.5 * (vLow + vHigh);
    zFinal = simulatePitch(vMid, theta0, &pars);

    if(zFinal > zend)
      vHigh = vMid;
    else
      vLow = vMid;
  }

  double vPitch = vMid;

  printf("********************************\n");
  printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n", xend, z0, theta0*180/M_PI);
  printf("v_pitch = %lf m/s\n", vPitch);
  printf("********************************\n");

  if(showPlot){
    theApp.SetIdleTimer(30, ".q");
    theApp.Run();
  }

  return 0;
}
