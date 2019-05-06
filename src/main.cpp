//  main.cpp
//  Integrative Experience
//
//  Created by Willie Beeson on 5/3/19.
//

#include <iostream>
#include <cmath>
#include <complex>

using namespace std;
typedef complex<double> dcomp;
typedef basic_string<char> string;

//Initialize variables
double wl[100] = {0};
dcomp t[100];
dcomp T[100];
double r;
double R;
double d[20];
double n[20];
double g[20];
dcomp m[20][2][2] = {{{0}}};
dcomp M[20][2][2] = {{{0}}};
double ps[20];
int nLayers;
string HL;

//Initialize constants
double n0 = 1.0;
double nSi = 1.45704;
double nTa = 2.1306;
double ns = 1.5;
double pi = acos(-1);
dcomp im (0,1);
double e0 = 8.85*pow(10,-12);
double mu0 = 4*pi*pow(10,-7);
double g0 = n0*pow(e0*mu0,0.5);
double gs = ns*pow(e0*mu0,0.5);
double peak;

int main(int argc, const char * argv[]) {
    //User input parameters (number of layers, interfered wavelength, high-low structure)
    cout << "Enter the number of layers (max 20): ";
    cin >> nLayers;
    cout << endl << "Enter target wavelength to be interfered in nm: ";
    cin >> peak;
    
    for (int i=1; i<=nLayers; i++){
        cout << endl << "Layer " << i << endl;
        
        cout << "High or low index of refraction? (enter H or L):" << endl;
        cin >> HL;
        if (HL == "H"){
            n[i] = nTa;
        }
        if (HL == "L"){
            n[i] = nSi;
        }
        d[i] = peak/(n[i]*4);
        cout << d[i] << endl;
    }
    
    //Loop over wavelengths
    for (int w=0; w<100; w++){
        wl[w] = 400 + w*4;
    
    //Gamma factors and phase shifts
    for (int l=1; l<=nLayers; l++){
        g[l] = n[l]*pow(e0*mu0,0.5);
        ps[l] = (2*pi/wl[w])*n[l]*d[l];
        //cout << l << " " << nl[l] << endl;
    }
    
    //Transfer matrix of each layer
    for (int l=1; l<=nLayers; l++){
        m[l][0][0] = cos(ps[l]);
        m[l][0][1] = im*sin(ps[l])/g[l];
        m[l][1][0] = im*g[l]*sin(ps[l]);
        m[l][1][1] = cos(ps[l]);
        //cout << l << " " << ps[l] << " " << m[l][0][1] << endl;
    }
    
    //Multiply transfer matrices
    int r1 = 2, r2 = 2, c2 = 2;
        for(int l = 2; l<=nLayers; l++){
            if (l == 2){
                for (int i = 0; i < r1; i++){
                    for (int j = 0; j < c2; j++){
                        M[l][i][j] = 0;
                        for (int k = 0; k < r2; k++){
                            M[l][i][j] += m[l-1][i][k] * m[l][k][j];
                            //cout << l << " " << M[l][0][0] << endl;
                        }
                    }
                }
            }
            if (l>2){
                for (int i = 0; i < r1; i++){
                    for (int j = 0; j < c2; j++){
                        M[l][i][j] = 0;
                        for (int k = 0; k < r2; k++){
                            M[l][i][j] += M[l-1][i][k] * m[l][k][j];
                            //cout << l << " " << k << " " << j << " " << m[l][k][j] << endl;
                        }
                    }
                }
            }
        }
    
    //Transmission coefficient
    t[w] = 2*g0/(g0*M[nLayers][0][0] + g0*gs*M[nLayers][0][1] + M[nLayers][1][0] + gs*M[nLayers][1][1]);
    
    //Transmission percentage
    T[w] = norm(t[w])*100;
    
    cout << real(T[w]) << endl;
    }
    return 0;
}
