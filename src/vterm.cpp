#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <TString.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>

using namespace std;

struct Params {
  double g;
  double m;
  double air_k;
};

// ODE 
double f_ri(double x, const vector<double> &y, void *params=0){ 
  (void) x;
  return y[1];
}

double f_vi(double x, const vector<double> &y, void *params=0){ 
  (void) x;
  (void) y;
  return 0.0;  // no air resistance
}

double f_rj(double x, const vector<double> &y, void *params=0){  
  (void) x;
  return y[3];
}

double f_vj(double x, const vector<double> &y, void *params=0){  
  (void) x;
  Params *p = (Params*)params;
  return -p->g;
}

double f_stop(double x, const vector<double> &y, void *params=0){
  (void) x;
  if (y[2] < 0) return 1;
  return 0;
}

// Energy
double total_energy(const vector<double> &y, Params *p){
  double vx = y[1];
  double vy = y[3];
  double h  = y[2];

  double T = 0.5 * p->m * (vx*vx + vy*vy);
  double U = p->m * p->g * h;

  return T + U;
}

int main(int argc, char **argv){

  Params pars;
  pars.g = 9.81;
  pars.m = 10.0;
  pars.air_k = 0.0;

  void *p_par = (void*) &pars;

  double theta = 45;
  double v0 = 100;

  int c;
  while ((c = getopt(argc, argv, "v:t:m:")) != -1){
    if (c == 'v') v0 = atof(optarg);
    if (c == 't') theta = atof(optarg);
    if (c == 'm') pars.m = atof(optarg);
  }


  TApplication theApp("App", &argc, argv);

  UInt_t dh = gClient->GetDisplayHeight()/2;
  UInt_t dw = 1.1*dh;

  vector<pfunc_t> v_fun(4);
  v_fun[0] = f_ri;
  v_fun[1] = f_vi;
  v_fun[2] = f_rj;
  v_fun[3] = f_vj;

  vector<double> y(4);
  y[0] = 0;
  y[1] = v0*cos(theta*M_PI/180);
  y[2] = 0;
  y[3] = v0*sin(theta*M_PI/180);

  double x = 0;
  double xmax = 20;
  int nsteps = 200;

  auto tgN = RK4SolveN(v_fun, y, nsteps, x, xmax, p_par, f_stop);

  //baseline energy is before we mess with number of steps
  TGraph *gEnergy = new TGraph();
  gEnergy->SetName("energy_baseline");
  gEnergy->SetTitle("Energy vs Time (200 steps); t [s]; Energy [J]");

  int N = tgN[0].GetN();

  for (int i = 0; i < N; i++){

    double t, xpos, vx, ypos, vy;

    tgN[0].GetPoint(i, t, xpos);
    tgN[1].GetPoint(i, t, vx);
    tgN[2].GetPoint(i, t, ypos);
    tgN[3].GetPoint(i, t, vy);

    vector<double> ytemp(4);
    ytemp[1] = vx;
    ytemp[2] = ypos;
    ytemp[3] = vy;

    double E = total_energy(ytemp, &pars);
    gEnergy->SetPoint(i, t, E);
  }

  //i use JSROOT extension in VSCode to open ROOT files but these energy plots show up blank
  //so I tested to see that my energy plot wasn't empty (it wasn't). my README says a bit more
  cout << "Baseline energy points: " << gEnergy->GetN() << endl;

  //more step sizes
  vector<int> stepList = {50, 100, 200, 400};
  vector<TGraph*> energyGraphs;

  for (int steps : stepList){

    double x_local = 0;

    vector<double> y_local(4);
    y_local[0] = 0;
    y_local[1] = v0*cos(theta*M_PI/180);
    y_local[2] = 0;
    y_local[3] = v0*sin(theta*M_PI/180);

    auto tgLocal = RK4SolveN(v_fun, y_local, steps, x_local, xmax, p_par, f_stop);

    TGraph *gE = new TGraph();
    gE->SetName(Form("energy_%dsteps", steps));
    gE->SetTitle(Form("Energy vs Time (%d steps); t [s]; Energy [J]", steps));

    int NL = tgLocal[0].GetN();

    for (int i = 0; i < NL; i++){

      double t, xpos, vx, ypos, vy;

      tgLocal[0].GetPoint(i, t, xpos);
      tgLocal[1].GetPoint(i, t, vx);
      tgLocal[2].GetPoint(i, t, ypos);
      tgLocal[3].GetPoint(i, t, vy);

      vector<double> ytemp(4);
      ytemp[1] = vx;
      ytemp[2] = ypos;
      ytemp[3] = vy;

      double E = total_energy(ytemp, &pars);
      gE->SetPoint(i, t, E);
    }

    energyGraphs.push_back(gE);
  }

  //RKnDemo puts everything in a ROOT file, so this script does too. Plotting the different energies isn't done here. 
  //look at test.cpp for that plot
  TFile *tf = new TFile("vterm.root", "RECREATE");

  for (unsigned i = 0; i < v_fun.size(); i++)
    tgN[i].Write();

  gEnergy->Write();

  for (auto g : energyGraphs)
    g->Write();

  tf->Close();

  cout << "Saved everything to vterm.root" << endl;
  cout << "Press Ctrl+C to exit\n";

  theApp.SetIdleTimer(30, ".q");
  theApp.Run();
}
