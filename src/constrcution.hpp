#include <cmath>
#include <iostream>
#include "define.hpp"

class CSSWM {
public:
    class patch {
    public: 
        // constructer
        patch();

        double hp[NX][NY], h[NX][NY], hm[NX][NY];
        double up[NX][NY], u[NX][NY], um[NX][NY];
        double vp[NX][NY], v[NX][NY], vm[NX][NY];

        double lon[NX][NY], lat[NX][NY];
        double lon_L[NX][NY], lat_L[NX][NY];
        double lon_R[NX][NY], lat_R[NX][NY];
        double lon_D[NX][NY], lat_D[NX][NY];
        double lon_U[NX][NY], lat_U[NX][NY];

        double x[NX][NY], y[NX][NY];
        double x_L[NX][NY], y_L[NX][NY];
        double x_R[NX][NY], y_R[NX][NY];
        double x_U[NX][NY], y_U[NX][NY];
        double x_D[NX][NY], y_D[NX][NY];

        double A[NX][NY][4], IA[NX][NY][4];
        double A_L[NX][NY][4], IA_L[NX][NY][4];
        double A_R[NX][NY][4], IA_R[NX][NY][4];
        double A_D[NX][NY][4], IA_D[NX][NY][4];
        double A_U[NX][NY][4], IA_U[NX][NY][4];
    };

    CSSWM();
    patch csswm[6];
    double sqrtG[NX][NY], gamma[NX][NY], gLower[NX][NY][4], gUpper[NX][NY][4];
    double sqrtG_L[NX][NY], gamma_L[NX][NY], gLower_L[NX][NY][4], gUpper_L[NX][NY][4];
    double sqrtG_R[NX][NY], gamma_R[NX][NY], gLower_R[NX][NY][4], gUpper_R[NX][NY][4];
    double sqrtG_D[NX][NY], gamma_D[NX][NY], gLower_D[NX][NY][4], gUpper_D[NX][NY][4];
    double sqrtG_U[NX][NY], gamma_U[NX][NY], gLower_U[NX][NY][4], gUpper_U[NX][NY][4];

    // ***********************************************************************************
    // In transform.cpp
    double Cube2Sphere_U(CSSWM &, int, int, int);
    double Cube2Sphere_V(CSSWM &, int, int, int);
    double Sphere2Cube_U(CSSWM &, int, int, int);
    double Sphere2Cube_V(CSSWM &, int, int, int);
    void matrixMul(double firstMatrix[4], double secondMatrix[4], double mult[2][2]);
    // ***********************************************************************************


private:
    void Construct_gamma_sqrtG_GUpper(double **alpha2D, double **beta2D, double gamma[NX][NY], double sqrtG[NX][NY], double gUpper[NX][NY][4], double gLower[NX][NY][4]);
    void Construct_p0123_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]);
    void Construct_p4_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]);
    void Construct_p5_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]);
};