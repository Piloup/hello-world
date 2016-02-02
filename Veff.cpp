#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "vecteur.h"
#define m 1
#define mu 1
#define a0 -10
#define b0 2
#define sigmax 0.5
#define sigmay 0.5
#define PI 3.14159265359
#define epsilon 1
#define alpha 2
using namespace std;

double gaus(double x,double x0,double sigma){
return exp(-pow(x-x0,2)/(2*sigma*sigma));
}


double Vpot(double x,double y,double xeq,double yeq){
return a0*(gaus(x,xeq,sigmax)*gaus(y,yeq,sigmay)+epsilon*gaus(x,-xeq,sigmax)*gaus(y,-yeq,sigmay));
}

double Vbord(double x,double y){
return b0*pow(pow(pow(x,2)+pow(y,2),0.5),alpha);
}

double Vbord_gradx(double x,double y){
return b0*alpha*2*x*pow(pow(pow(x,2)+pow(y,2),0.5),alpha-1);
}

double Vbord_grady(double x,double y){
return  b0*alpha*2*y*pow(pow(pow(x,2)+pow(y,2),0.5),alpha-1);
}

double V_gradx(double x,double y,double xeq,double yeq){
return -(x-xeq)*Vpot(x,y,xeq,yeq)/pow(sigmax,2);
}

double V_grady(double x,double y,double xeq,double yeq){
return -(y-yeq)*Vpot(x,y,xeq,yeq)/pow(sigmay,2);
}

double Verlet_x(Vector &pos,Vector vit,Vector acc,double dt,int ite){
vector <double> x=pos.V_GetX();
vector <double> vx=vit.V_GetX();
vector <double> ax=acc.V_GetX();
//cout<<x[ite]<<" "<<vx[ite]<<" "<<ax[ite]/2;
double xpdt=x[ite]+dt*vx[ite]+pow(dt,2)*ax[ite]/2;
return xpdt;
}

double Verlet_y(Vector pos,Vector vit,Vector acc,double dt,int ite){
vector <double> y=pos.V_GetY();
vector <double> vy=vit.V_GetY();
vector <double> ay=acc.V_GetY();
double ypdt=y[ite]+dt*vy[ite]+pow(dt,2)*ay[ite]/2;
return ypdt;
}


double Verlet_vx(Vector eq,Vector pos,Vector vit,Vector acc,double dt,int ite,double w){
vector <double> x=pos.V_GetX();
vector <double> y=pos.V_GetY();
vector <double> xeq=eq.V_GetX();
vector <double> yeq=eq.V_GetY();
vector <double> vx=vit.V_GetX();
vector <double> ax=acc.V_GetX();
double vxpdt=(vx[ite]+mu*w*(y[ite]+y[ite+1])+dt*ax[ite]/2-dt*Vbord_gradx(x[ite+1],y[ite+1])-dt*V_gradx(x[ite+1],y[ite+1],xeq[ite+1],yeq[ite+1]))/(1+dt*mu/2);
return vxpdt;
}


double Verlet_vy(Vector eq,Vector pos, Vector vit, Vector acc,double dt,int ite,double w){
vector <double> x=pos.V_GetX();
vector <double> y=pos.V_GetY();
vector <double> xeq=eq.V_GetX();
vector <double> yeq=eq.V_GetY();
vector <double> vy=vit.V_GetY();
vector <double> ay=acc.V_GetY();
double vypdt=(vy[ite]-mu*w*(x[ite]+x[ite+1])+dt*ay[ite]/2-dt*Vbord_grady(x[ite+1],y[ite+1])-dt*V_grady(x[ite+1],y[ite+1],xeq[ite+1],yeq[ite+1]))/(1+dt*mu/2);
return vypdt;
}

double Verlet_vit(Vector eq,Vector pos,Vector vit, Vector acc,double dt,double ite,double w){



}


double Verlet_ax(Vector eq,Vector pos,Vector vit,Vector acc,double dt,int ite,double w){
double x=pos.V_GetX(ite+1);
double y=pos.V_GetY(ite+1);
double xeq=eq.V_GetX(ite+1);
double yeq=eq.V_GetY(ite+1);
double vx=vit.V_GetX(ite+1);
double vy=vit.V_GetY(ite+1);
double axpdt=-mu*vx/m+mu*w*y/m-Vbord_gradx(x,y)/m-V_gradx(x,y,xeq,yeq)/m+x*w*w/m+2*w*vy/m;
//cout<<axpdt<<endl;
return axpdt;
}

double Verlet_ay(Vector eq,Vector pos,Vector vit,Vector acc,double dt,int ite,double w){
double x=pos.V_GetX(ite+1);
double y=pos.V_GetY(ite+1);
double xeq=eq.V_GetX(ite+1);
double yeq=eq.V_GetY(ite+1);
double vy=vit.V_GetY(ite+1);
double vx=vit.V_GetX(ite+1);
double aypdt=-mu*vy/m-mu*w*x/m-V_grady(x,y,xeq,yeq)/m-Vbord_grady(x,y)/m-w*w*y+2*w*vx;
//cout<<-Vpot(x,y,xeq,yeq)<<endl;
return aypdt;
}

void Verlet(Vector eq,Vector pos,Vector vit,Vector acc,double dt,double t,double w){
ofstream data("Rtournant.out");
ofstream pot("potentiel.out");
vector <double> x=pos.V_GetX(); vector <double> y=pos.V_GetY();
vector <double> xeq=eq.V_GetX(); vector <double> yeq=eq.V_GetY();
vector <double> vx=vit.V_GetX(); vector <double> vy=vit.V_GetY();
vector <double> ay=acc.V_GetY(); vector <double> ax=acc.V_GetX();
data<<eq.V_GetX(0)<<" "<<eq.V_GetY(0)<<endl;
for(double i=0;i<t;i+=dt){
int ite=i/dt;

data<<pos.V_GetX(ite)<<" "<<pos.V_GetY(ite)<<" "<<vit.V_GetX(ite)<<" "<<vit.V_GetY(ite)<<" "<<acc.V_GetX(ite)<<" "<<acc.V_GetY(ite)<<endl;
pos.V_SetX(Verlet_x(pos,vit,acc,dt,ite),ite+1);
pos.V_SetY(Verlet_y(pos,vit,acc,dt,ite),ite+1);
vit.V_SetX(Verlet_vx(eq,pos,vit,acc,dt,ite,w),ite+1);
vit.V_SetY(Verlet_vy(eq,pos,vit,acc,dt,ite,w),ite+1);
acc.V_SetX(Verlet_ax(eq,pos,vit,acc,dt,ite,w),ite+1);
acc.V_SetY(Verlet_ay(eq,pos,vit,acc,dt,ite,w),ite+1);
//cout<<"("<<x[ite+1]<<";"<<y[ite+1]<<"): "<<Vpot(x[ite],y[ite],xeq[ite],yeq[ite])<<" at ("<<xeq[ite]<<";"<<yeq[ite]<<")"<<endl;
pot<<Vpot(x[ite],y[ite],xeq[ite],yeq[ite])<<" "<<V_gradx(x[ite],y[ite],xeq[ite],yeq[ite])<<" "<<V_grady(x[ite],y[ite],xeq[ite],yeq[ite])<<endl;
	}

}

void rk_1(Vector eq,Vector &pos,Vector &vit,Vector &acc, double dt,double t){

vector <double> xeq=eq.V_GetX();vector <double> yeq=eq.V_GetY();
vector <double> x=pos.V_GetX();vector <double> y=pos.V_GetY();
vector <double> vx=vit.V_GetX();vector <double> vy=vit.V_GetY();
vector <double> ax=acc.V_GetX();vector <double> ay=acc.V_GetY();
for(double i=dt;i<t;i+=dt){
	x[i/dt]=x[i/dt-1]+dt*vx[i/dt-1]+dt*dt*ax[i/dt-1]/2;
	y[i/dt]=y[i/dt-1]+dt*vy[i/dt-1]+dt*dt*ay[i/dt-1]/2;
	vx[i/dt]=vx[i/dt-1]+dt*ax[i/dt-1];
	vy[i/dt]=vy[i/dt-1]+dt*ay[i/dt-1];
	ax[i/dt]=-mu*vx[i/dt]-V_gradx(x[i/dt],y[i/dt],xeq[i/dt],yeq[i/dt])-Vbord_grady(x[i/dt],y[i/dt]);
	ay[i/dt]=-mu*vy[i/dt]/m-V_grady(x[i/dt],y[i/dt],xeq[i/dt],yeq[i/dt])/m-Vbord_grady(x[i/dt],y[i/dt])/m;
	}
pos.V_SetX(x);pos.V_SetY(y);vit.V_SetX(vx);vit.V_SetY(vy);acc.V_SetX(ax);acc.V_SetY(ay);
}



int main(){
Vector eq; Vector pos; Vector vit;Vector acc;
double t=10, dt=0.01;
double w=0.1;
double xeq=0,yeq=1,vxin=0,vyin=0;
eq.V_resize(t/dt);

for(double i=0;i<t;i+=dt){
eq.V_SetX(xeq,int(i/dt)); eq.V_SetY(yeq,int(i/dt));
}

pos.V_resize(t/dt);pos.V_SetX(xeq+0.2,0);pos.V_SetY(yeq,0);
vit.V_resize(t/dt);vit.V_SetX(vxin,0);vit.V_SetY(vyin,0);
acc.V_resize(t/dt);acc.V_SetX(-V_gradx(pos.V_GetX(0),pos.V_GetY(0),xeq,yeq),0);acc.V_SetY(-V_grady(pos.V_GetX(0),pos.V_GetY(0),xeq,yeq),0);

ofstream data("Rtournant.out");
Verlet(eq,pos,vit,acc,dt,t,w);
//rk_1(eq,pos,vit,acc,dt,t);
for(double i=0;i<t;i+=dt){
//cout<<i/dt<<" ";
data<<i/dt<<" "<<pos.V_GetX(int(i/dt))<<" "<<pos.V_GetY(int(i/dt))<<" : "<<vit.V_norme(i/dt)<<" ;  "<<Vpot(pos.V_GetX(i/dt),pos.V_GetY(i/dt),xeq,yeq)<<endl;
}
system("gnuplot script.gnu");


ofstream potentiel("3Dpotentiel.out");
ofstream fx("fx.out");
ofstream fy("fy.out");
for(double x=-3;x<3;x+=0.1){
	for(double y=-3;y<3;y+=0.1){
		potentiel<<x<<" "<<y<<" "<<Vpot(x,y,xeq,yeq)+Vbord(x,y)<<endl;
		fx<<x<<" "<<y<<" "<<(x-xeq)*gaus(x,xeq,sigmax)*gaus(y,yeq,sigmay)*a0/pow(sigmax,2)+epsilon*a0*(x+xeq)*gaus(x,-xeq,sigmax)*gaus(y,-yeq,sigmay)/pow(sigmax,2)<<endl;
		fy<<x<<" "<<y<<" "<<(y-yeq)*gaus(x,xeq,sigmax)*gaus(y,yeq,sigmay)*a0/pow(sigmay,2)+epsilon*a0*(y+yeq)*gaus(x,-xeq,sigmax)*gaus(y,-yeq,sigmay)/pow(sigmay,2)<<endl;
	}
}
system("gnuplot potentiel_script.gnu");

return 0;
}	
