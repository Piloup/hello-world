#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
using namespace std;

class Vector{
private:
public:
int size;
vector <vector <double> > coord;

void V_resize(int s){
coord.resize(2);
coord[0].resize(s);coord[1].resize(s);
}

void V_SetX(double x,int t){
	coord[0][t]=x;
}

void V_SetY(double y,int t){
	coord[1][t]=y;
}

void V_SetX(vector <double> x){
	coord[0]=x;
}

void V_SetY(vector <double> y){
	coord[1]=y;
}

double V_norme(int t){
	return pow(pow(coord[0][t],2)+pow(coord[1][t],2),1/2);	
}

double V_GetX(int t){
	return coord[0][t];
}

double V_GetY(int t){
	return  coord[1][t];
}

vector <double> V_GetX(){
	return coord[0];
}

vector <double> V_GetY(){
	return coord[1];
}
	
	};
