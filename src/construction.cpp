#include "constrcution.hpp"

CSSWM::patch::patch() {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            hp[i][j] = h[i][j] = hm[i][j] = FILLVALUE;
            up[i][j] = u[i][j] = um[i][j] = FILLVALUE;
            vp[i][j] = v[i][j] = vm[i][j] = FILLVALUE;

            lon[i][j] = lat[i][j] = FILLVALUE;
            lon_L[i][j] = lat_L[i][j] = FILLVALUE;
            lon_R[i][j] = lat_R[i][j] = FILLVALUE;
            lon_U[i][j] = lat_U[i][j] = FILLVALUE;
            lon_D[i][j] = lat_D[i][j] = FILLVALUE;

            x[i][j] = y[i][j] = FILLVALUE;
            x_L[i][j] = y_L[i][j] = FILLVALUE;
            x_R[i][j] = y_R[i][j] = FILLVALUE;
            x_U[i][j] = y_U[i][j] = FILLVALUE;
            x_D[i][j] = y_D[i][j] = FILLVALUE;
        }
    }
}

CSSWM::CSSWM() {
    double *alpha = new double[NX], *beta = new double[NY];
    double *alpha_L = new double[NX], *alpha_R = new double[NX], *beta_U = new double[NY], *beta_D = new double[NY];

    for (int i = 0; i < NX; i++) {
        alpha[i] = -M_PI/4. + (M_PI/2.) / (NX-4) * (i-1.5);
        alpha_L[i] = -M_PI/4. + (M_PI/2.) / (NX-4) * (i-2);
        alpha_R[i] = -M_PI/4. + (M_PI/2.) / (NX-4) * (i-1);

    }

    delete[] alpha, delete[] beta, delete[] alpha_L, delete[] beta_U, delete[] alpha_R, delete[] beta_D;
}