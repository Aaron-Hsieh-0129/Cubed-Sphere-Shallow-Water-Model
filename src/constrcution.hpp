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

    };

    CSSWM();
    
    patch csswm[6];
    

private:



};