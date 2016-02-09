#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstdlib>
#define m 1
#define sigmax 0.2
#define sigmay 0.2
#define epsilon 1
#define a0 -10
#define b0 0
#define alpha 4
#define mu 0.1
using namespace std;


double gaus(double x,double x0,double sigma){
return exp(-pow(x-x0,2)/(2*sigma*sigma));
}


double Vpot(double x,double y,double xeq,double yeq){
return a0*(gaus(x,xeq,sigmax)*gaus(y,yeq,sigmay));
}

double Vbord(double x,double y){
return b0*pow(pow(pow(x,2)+pow(y,2),0.5),alpha);
}

double Vbord_gradx(double x,double y){
return b0*alpha*2.0*x*pow(pow(pow(x,2)+pow(y,2),0.5),alpha-1);
}

double Vbord_grady(double x,double y){

return  b0*alpha*2.0*y*pow(pow(pow(x,2)+pow(y,2),0.5),alpha-1);
}

double V_gradx(double x,double y,double xeq,double yeq){
return -(x-xeq)*Vpot(x,y,xeq,yeq)/pow(sigmax,2)-(x+xeq)*epsilon*Vpot(x,y,-xeq,-yeq)/pow(sigmax,2);
}

double V_grady(double x,double y,double xeq,double yeq){
return -(y-yeq)*Vpot(x,y,xeq,yeq)/pow(sigmay,2)-(y+yeq)*epsilon*Vpot(x,y,-xeq,-yeq)/pow(sigmay,2);
}

double Fx(double x,double y,double xeq,double yeq,double w){
double f;
	if((pow(x,2)+pow(y,2))>4*(pow(xeq,2)+pow(yeq,2))){
f=m*w*w*x+w*mu*y-V_gradx(x,y,xeq,yeq)-Vbord_gradx(x-xeq,y-yeq);
	}
	else{f=m*w*w*x+w*mu*y-V_gradx(x,y,xeq,yeq);}
	return f;
}

double Fy(double x,double y,double xeq,double yeq,double w){
double f;
	if((pow(x,2)+pow(y,2))>4*(pow(xeq,2)+pow(yeq,2))){
f=-w*mu*x+m*w*w*y-V_grady(x,y,xeq,yeq)-Vbord_grady(x-xeq,y-yeq);
	
	}
	else{	
		f=-w*mu*x+m*w*w*y-V_grady(x,y,xeq,yeq);
		}
	return f;
}


int main(){
vector <vector <vector <double> > >gridx;
vector <vector <vector <double> > >gridy;
double r=1;double w=0;double dx=0.1,dy=0.1,dw=0.1;
double xeq=0,yeq=1; double xmin=-2,xmax=2,ymin=-2,ymax=2,wmin=0,wmax=10;
vector <double> temp((ymax-ymin)/dy,0);
gridx.resize((xmax-xmin)/dx);
gridy.resize((xmax-xmin)/dx);
for(int i=0;i<(xmax-xmin)/dx;i++){
	gridx[i].resize((ymax-ymin)/dy);
	gridy[i].resize((ymax-ymin)/dy);
	for(int j=0;j<(ymax-ymin)/dy;j++){
		gridx[i][j].resize((wmax-wmin)/dw);
		gridy[i][j].resize((wmax-wmin)/dw);
		}
	}


ofstream data("./force.out");

for(w=wmin;w<wmax;w+=dw){
//cin>>w;
for(double x=xmin;x<xmax;x+=dx){
	for(double y=ymin;y<ymax;y+=dy){
		//gridx[(x-xmin)/dx][(y-ymin)/dy][w/dw]=Fx(x,y,xeq,yeq,w);	
		//gridy[(x-xmin)/dx][(y-ymin)/dy][w/dw]=Fy(x,y,xeq,yeq,w);
		data<<x<<" "<<y<<" "<<pow(pow(Fx(x,y,xeq,yeq,w),2)+pow(Fy(x,y,xeq,yeq,w),2),0.5)<<endl;//<<pow(pow(Fx(x,y,xeq,yeq,w),2)+pow(Fy(x,y,xeq,yeq,w),2),1/2)<<endl;			
		}
		//grid[x/dx].push_back(temp);
	//grid.push_back(0);
	}
}
system("gnuplot force_w_0_10.gnu");

double fmin=100;
double xp,yp;
for(double x=0;x<xmax;x+=dx){
	for(double y=0;y<ymax;y+=dy){
	if(pow(pow(Fx(x,y,xeq,yeq,w),2)+pow(Fy(x,y,xeq,yeq,w),2),0.5)<fmin){
		fmin=pow(pow(Fx(x,y,xeq,yeq,w),2)+pow(Fy(x,y,xeq,yeq,w),2),1/2);xp=x;yp=y;
		}
	}
}
cout<<w<<" :("<<xp<<";"<<yp<<") "<<fmin<<endl;
return 0;
}



