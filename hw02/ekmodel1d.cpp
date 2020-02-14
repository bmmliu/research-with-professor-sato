#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

int main(void) {
  const int num=400;
  
  //array
  double v[num];//membrane voltage
  double h[num];//Na channel inactivation gate
  double f[num];//Ca channel inactivation gate
  double stim[num];
  double tmp[num];//temporary variable to solve diffusion eqn

  //initial values
  for (int i=0;i<num;i++) {
    v[i]=0;
    h[i]=1.0;
    f[i]=0.9;
    stim[i]=0;
  }

  //constants
  const double dt=0.1;// time step (0.1 ms)
  const double pcl=200;//pacing cycle length (300 ms)
  const int itr=10;//# of beats
  double tmax=pcl*itr;//time to simulate

  int tnmax=tmax/dt;//convert to int
  int pcln=pcl/dt;//convert to int
  int durn=1.0/dt;//duration of stimulation (1.0 ms)

  ofstream os("result.txt");

  //main loop
  for (int tn=0;tn<tnmax;tn++) {
    if (tn%100==0) {//write results every 10 ms
      double t=tn*dt;
      for (int i=0;i<num;i++) {
        os<<v[i]<<"\t";
      }
      os<<endl;
    }

    //stimlate the end of the cable
    if (tn%pcln < durn) {
      for (int i=0;i<5;i++) {
        stim[i]=0.3;
      }
    }
    else {
      for (int i=0;i<5;i++) {
        stim[i]=0;
      }
    }

    //Euler method
    for (int i=0;i<num;i++) {
      const double tauso=15;
      const double taufi=0.8;
      const double tauh1=4.8;
      const double tauh2=10.0;
      const double tausi=4.0;
      const double tauf1=100;
      const double tauf2=30;
      double minf=pow((v[i]/0.2),6)/(1+pow((v[i]/0.2),6));
      double hinf=1/(1+pow((v[i]/0.1),6));
      double dinf=pow((v[i]/0.4),4)/(1+pow((v[i]/0.4),4));
      double finf=1/(1+pow((v[i]/0.1),4));
      double tauh=tauh1+tauh2*exp(-20*pow((v[i]-0.1),2));
      double tauf=tauf2+(tauf1-tauf2)*v[i]*v[i]*v[i];

      double jfi=h[i]*minf*(v[i]-1.3)/taufi;//Fast inward current (Na current)
      double jsi=f[i]*dinf*(v[i]-1.4)/tausi;//Slow inward current (Ca current)
      double jso=(1-exp(-4*v[i]))/tauso;//outward current (K current)
      double ion=-(jfi+jsi+jso-stim[i]);//total transmembrane current

      double dh=(hinf-h[i])/tauh;
      double df=(finf-f[i])/tauf;

      //update variables
      v[i]+=ion*dt;
      h[i]+=dh*dt;
      f[i]+=df*dt;
    }
    //solve Diffusion
    const double dfu=0.0005;//diffusion coefficient
    const double dx=0.015;//0.015 cm
    //non-flux boundary
    v[0]=v[2];
    v[num-1]=v[num-3];
    for (int i=1;i<num-1;i++) {
      tmp[i]=v[i]+(v[i-1]+v[i+1]-2*v[i])*dfu*dt/(dx*dx);
    }
    for (int i=1;i<num-1;i++) {
      v[i]=tmp[i];
    }
  }
  return 0;
}
