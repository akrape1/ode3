#include "RKn.hpp"
#include "TROOT.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>   

using namespace std;

struct Params {
  double g;      // gravity
  double B;      // Magnus coefficient (dimensionless)
  double omega;  // spin rate [rad/s]
  double phi;    // spin axis angle [rad]
};

// some unit conversions because I do everything in m and s but problem wants ft
double mph_to_mps(double v){ return v * 0.44704; }
double rpm_to_rad(double w){ return w * 2.0 * M_PI / 60.0; }
double ft_to_m(double x){ return x * 0.3048; }
double m_to_ft(double x){ return x / 0.3048; }

// drag law 
double F_drag(double v){
  const double v_d = 35.0;     // m/s
  const double Delta = 5.0;    // m/s
  return 0.0039 + 0.0058 / (1.0 + exp((v - v_d)/Delta));
}

// ODEs: 
// y[0] = x, y[1] = v_x
// y[2] = y, y[3] = v_y
// y[4] = z, y[5] = v_z

double f_x(double t, const vector<double>& y, void *p){ 
  (void)t; (void)p;
  return y[1]; 
}

double f_y(double t, const vector<double>& y, void *p){ 
  (void)t; (void)p;
  return y[3]; 
}

double f_z(double t, const vector<double>& y, void *p){ 
  (void)t; (void)p;
  return y[5]; 
}

double f_vx(double t, const vector<double>& y, void *pp){
  (void)t;
  Params *p = (Params*)pp;
  double vx=y[1], vy=y[3], vz=y[5];
  double v = sqrt(vx*vx + vy*vy + vz*vz);
  double F = F_drag(v);
  return -F * v * vx + p->B * p->omega * ( vz * sin(p->phi) - vy * cos(p->phi) );
}

double f_vy(double t, const vector<double>& y, void *pp){
  (void)t;
  Params *p = (Params*)pp;
  double vx=y[1], vy=y[3], vz=y[5];
  double v = sqrt(vx*vx + vy*vy + vz*vz);
  double F = F_drag(v);
  return -F * v * vy + p->B * p->omega * vx * cos(p->phi);
}

double f_vz(double t, const vector<double>& y, void *pp){
  (void)t;
  Params *p = (Params*)pp;
  double vx=y[1], vy=y[3], vz=y[5];
  double v = sqrt(vx*vx + vy*vy + vz*vz);
  double F = F_drag(v);
  return -p->g - F * v * vz - p->B * p->omega * vx * sin(p->phi);
}

// Stop when ball reaches plate at x = 18.44 m (60 ft)
double f_stop(double t, const vector<double>& y, void *pp){
  (void)t; (void)pp;
  const double plate_x = 18.44;  // m
  if (y[0] >= plate_x) return 1.0;
  return 0.0;
}

// Setup initial conditions given v0, theta (deg)
void SetupInitialConditions(vector<double> &y0, double v0_mph, double theta_deg){
  double v0 = mph_to_mps(v0_mph);
  double theta = theta_deg * M_PI/180.0;

  y0.resize(6);
  // initial position (x=y=z=0)
  y0[0] = 0.0;  // x
  y0[2] = 0.0;  // y
  y0[4] = 0.0;  // z

  // initial velocity components
  y0[1] = v0 * cos(theta);   // v_x
  y0[3] = 0.0;               // v_y
  y0[5] = v0 * sin(theta);   // v_z
}

int main(int argc, char **argv){

  gROOT->SetBatch(kTRUE); //i run in vscode so it's like batch mode

  // default: curveball (ip=1)
  int ip = 1;
  int opt;
  while ((opt = getopt(argc, argv, "p:")) != -1){
    switch(opt){
      case 'p':
        ip = atoi(optarg);
        break;
      default:
        break;
    }
  }

  // Set physical + pitch parameters
  Params pars;
  pars.g = 9.81;
  pars.B = 4.1e-4;   // Magnus coefficient (dimensionless)

  double v0_mph;
  double theta0_deg = 1.0;
  double omega_rpm = 1800.0;
  double phi_deg;
  std::string pitchName;

  if (ip == 0){
    pitchName = "slider";
    v0_mph = 85.0;
    phi_deg = 0.0;
  }
  else if (ip == 1){
    pitchName = "curveball";
    v0_mph = 85.0;
    phi_deg = 45.0;
  }
  else if (ip == 2){
    pitchName = "screwball";
    v0_mph = 85.0;
    phi_deg = 135.0;
  }
  else {
    pitchName = "fastball";
    v0_mph = 95.0;
    phi_deg = 225.0;
  }

  pars.omega = rpm_to_rad(omega_rpm);      // rad/s
  pars.phi   = phi_deg * M_PI/180.0;      // rad

  // Initial conditions
  vector<double> y0;
  SetupInitialConditions(y0, v0_mph, theta0_deg);

  vector<pfunc_t> funcs(6);
  funcs[0] = f_x;
  funcs[1] = f_vx;
  funcs[2] = f_y;
  funcs[3] = f_vy;
  funcs[4] = f_z;
  funcs[5] = f_vz;

  double l = 18.44;                        // m
  double v0_mps = mph_to_mps(v0_mph);
  double T_est = l / v0_mps;              // estimated time-of-flight
  double h = 1e-4 * T_est;                // time step [s]
  int nmax = 20000;                       // max steps (backup)

  double t = 0.0;                         // start time

  auto sol = RK4SolveN(funcs, y0, h, t, &pars, f_stop, nmax);

  int N = sol[0].GetN();

 
  double tt, x_m, y_m, z_m, vx_m, vy_m, vz_m;
  sol[0].GetPoint(N-1, tt, x_m);
  sol[2].GetPoint(N-1, tt, y_m);
  sol[4].GetPoint(N-1, tt, z_m);
  sol[1].GetPoint(N-1, tt, vx_m);
  sol[3].GetPoint(N-1, tt, vy_m);
  sol[5].GetPoint(N-1, tt, vz_m);


  //this is changed from the one in the original script but it's the same info
  printf("********************************\n");
  printf("Pitch type: %s\n", pitchName.c_str());
  printf("Coordinates when x = 60 ft:\n");
  printf("(x,y,z) = (%lf, %lf, %lf) ft\n",
         m_to_ft(x_m), m_to_ft(y_m), m_to_ft(z_m));
  printf("(vx,vy,vz) = (%lf, %lf, %lf) ft/s\n",
         m_to_ft(vx_m), m_to_ft(vy_m), m_to_ft(vz_m));
  printf("********************************\n");

  TGraph *gY = new TGraph();
  TGraph *gZ = new TGraph();

  for (int i = 0; i < N; i++){
    double t_i, x_i, y_i, z_i;

    sol[0].GetPoint(i, t_i, x_i);  // x(t)
    sol[2].GetPoint(i, t_i, y_i);  // y(t)
    sol[4].GetPoint(i, t_i, z_i);  // z(t)

    double x_ft = m_to_ft(x_i);
    double y_ft = m_to_ft(y_i);
    double z_ft = m_to_ft(z_i);

    gY->SetPoint(i, x_ft, y_ft);
    gZ->SetPoint(i, x_ft, z_ft);
  }

  TCanvas *c1 = new TCanvas("c1","Trajectory",800,600);

  gZ->SetLineWidth(2);
  gZ->SetLineColor(kBlue);
  gZ->SetTitle((pitchName + " pitch; x (ft); displacement (ft)").c_str());

  gY->SetLineWidth(2);
  gY->SetLineColor(kRed);

  gZ->Draw("AL");
  gY->Draw("L SAME");

  TLegend *leg = new TLegend(0.65,0.75,0.9,0.9);
  leg->AddEntry(gZ,"Vertical displacement z(x)","l");
  leg->AddEntry(gY,"Horizontal displacement y(x)","l");
  leg->Draw();

  std::string outfile = pitchName + ".png";
  c1->SaveAs(outfile.c_str());

  return 0;
}
