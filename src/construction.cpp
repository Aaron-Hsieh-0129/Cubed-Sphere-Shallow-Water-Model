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

void CSSWM::Construct_gamma_sqrtG_GUpper(double **alpha2D, double **beta2D, double gamma[NX][NY], double sqrtG[NX][NY], double gUpper[NX][NY][4], double gLower[NX][NY][4]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            gamma[i][j] = sqrt(1 + pow(tan(alpha2D[i][j]), 2) + pow(tan(beta2D[i][j]), 2));
            sqrtG[i][j] = 1. / (pow(gamma[i][j], 3) * pow(cos(alpha2D[i][j]), 2) * pow(cos(beta2D[i][j]), 2));

            gUpper[i][j][0] = pow((gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j])), 2) * (1 + pow(tan(beta2D[i][j]), 2));
            gUpper[i][j][1] = pow((gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j])), 2) * (tan(alpha2D[i][j]) * tan(beta2D[i][j]));
            gUpper[i][j][2] = gUpper[i][j][1];
            gUpper[i][j][3] = pow((gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j])), 2) * (1 + pow(tan(alpha2D[i][j]), 2));

            gLower[i][j][0] = 1. / (pow(gamma[i][j], 4) * pow(cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2)) * (1 + pow(tan(alpha2D[i][j]), 2));
            gLower[i][j][1] = 1. / (pow(gamma[i][j], 4) * pow(cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2)) * (-tan(alpha2D[i][j]) * tan(beta2D[i][j]));
            gLower[i][j][2] = gLower[i][j][1];
            gLower[i][j][3] = 1. / (pow(gamma[i][j], 4) * pow(cos(alpha2D[i][j]) * cos(beta2D[i][j]), 2)) * (1 + pow(tan(beta2D[i][j]), 2));
        }
    }
    return;
}

void CSSWM::Construct_p0123_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            // lon/lat
            lon[i][j] = alpha2D[i][j] + p * M_PI/2.;
            lat[i][j] = atan(tan(beta2D[i][j]) * cos(alpha2D[i][j]));

            // x/y
            x[i][j] = RADIUS * (lon[i][j] - p * M_PI/2.);
            y[i][j] = RADIUS * atan(tan(lat[i][j]) / cos(lon[i][j] - p * M_PI/2.));

            // A/IA
            A[i][j][0] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * gamma[i][j] * cos(beta2D[i][j]);
            A[i][j][1] = 0;
            A[i][j][2] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-tan(alpha2D[i][j]) * sin(beta2D[i][j]));
            A[i][j][3] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) / cos(beta2D[i][j]);

            IA[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) / cos(beta2D[i][j]);
            IA[i][j][1] = 0;
            IA[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (tan(alpha2D[i][j]) * sin(beta2D[i][j]));
            IA[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]));

            if (lon[i][j] < 0) lon[i][j] += 2 * M_PI;
        }
    }
}

void CSSWM::Construct_p4_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            // lon/lat
            lon[i][j] = atan2(tan(alpha2D[i][j]), -tan(beta2D[i][j]));
            lat[i][j] = atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2)+pow(tan(beta2D[i][j]), 2)));

            // x/y
            x[i][j] = RADIUS * atan(sin(lon[i][j]) / tan(lat[i][j]));
            y[i][j] = RADIUS * atan(-cos(lon[i][j]) / tan(lat[i][j]));

            // A/AInverse
            A[i][j][0] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(lon[i][j]));
            A[i][j][1] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(lon[i][j]));
            A[i][j][2] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(lon[i][j]));
            A[i][j][3] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(lon[i][j]));

            IA[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(lon[i][j]));
            IA[i][j][1] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (-gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(lon[i][j]));
            IA[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(lon[i][j]));
            IA[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(lon[i][j]));

            if (lon[i][j] < 0) lon[i][j] += 2 * M_PI;
        }
    }
}

void CSSWM::Construct_p5_lonlat_xy_AIA(int p, double **alpha2D, double **beta2D, double gamma[NX][NY], double lon[NX][NY], double lat[NX][NY], double x[NX][NY], double y[NX][NY], double A[NX][NY][4], double IA[NX][NY][4]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            // lon/lat
            lon[i][j] = atan2(tan(alpha2D[i][j]), tan(beta2D[i][j]));
            lat[i][j] = -atan(1 / sqrt(pow(tan(alpha2D[i][j]), 2)+pow(tan(beta2D[i][j]), 2)));

            // x/y
            x[i][j] = RADIUS * atan(-sin(lon[i][j]) / tan(lat[i][j]));
            y[i][j] = RADIUS * atan(-cos(lon[i][j]) / tan(lat[i][j]));

            // A/AInverse
            A[i][j][0] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(lon[i][j]));
            A[i][j][1] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (-gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(lon[i][j]));
            A[i][j][2] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(lon[i][j]));
            A[i][j][3] = 1 / (pow(gamma[i][j], 2) * cos(alpha2D[i][j]) * cos(beta2D[i][j])) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(lon[i][j]));

            IA[i][j][0] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (cos(alpha2D[i][j]) / cos(beta2D[i][j]) * cos(lon[i][j]));
            IA[i][j][1] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(alpha2D[i][j]) / cos(beta2D[i][j]) * sin(lon[i][j]));
            IA[i][j][2] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (-cos(beta2D[i][j]) / cos(alpha2D[i][j]) * sin(lon[i][j]));
            IA[i][j][3] = gamma[i][j] * cos(alpha2D[i][j]) * cos(beta2D[i][j]) * (gamma[i][j] * cos(beta2D[i][j]) / cos(alpha2D[i][j]) * cos(lon[i][j]));

            if (lon[i][j] < 0) lon[i][j] += 2 * M_PI;
        }
    }
}

CSSWM::CSSWM() {
    // Init new 1D array
    double *alpha = new double[NX], *beta = new double[NY];
    double *alpha_L = new double[NX], *alpha_R = new double[NX], *beta_U = new double[NY], *beta_D = new double[NY];

    for (int i = 0; i < NX; i++) {
        alpha[i] = -M_PI/4. + (M_PI/2.) / (NX-4) * (i-1.5);
        alpha_L[i] = -M_PI/4. + (M_PI/2.) / (NX-4) * (i-2);
        alpha_R[i] = -M_PI/4. + (M_PI/2.) / (NX-4) * (i-1);
    }
    for (int j = 0; j < NY; j++) {
        beta[j] = -M_PI/4. + (M_PI/2.) / (NX-4) * (j-1.5);
        beta_D[j] = -M_PI/4. + (M_PI/2.) / (NX-4) * (j-2);
        beta_U[j] = -M_PI/4. + (M_PI/2.) / (NX-4) * (j-1);
    }

    // Init new 2D array
    double **alpha2D = new double *[NX], **beta2D = new double *[NY];
    double **alpha2D_L = new double *[NX], **alpha2D_R = new double *[NX];
    double **beta2D_U = new double *[NX], **beta2D_D = new double *[NX];
    for (int i = 0; i < NX; i++) {
        alpha2D[i] = new double[NY];
        alpha2D_L[i] = new double [NY];
        alpha2D_R[i] = new double [NY];
        beta2D[i] = new double [NY];
        beta2D_U[i] = new double [NY];
        beta2D_D[i] = new double [NY];
    }
    
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            alpha2D[i][j] = alpha[i];
            beta2D[i][j] = beta[j];
            alpha2D_L[i][j] = alpha_L[i];
            alpha2D_R[i][j] = alpha_R[i];
            beta2D_U[i][j] = beta_U[j];
            beta2D_D[i][j] = beta_D[j];
        }
    }

    Construct_gamma_sqrtG_GUpper(alpha2D, beta2D, gamma, sqrtG, gUpper, gLower);
    Construct_gamma_sqrtG_GUpper(alpha2D_L, beta2D, gamma_L, sqrtG_L, gUpper_L, gLower_L);
    Construct_gamma_sqrtG_GUpper(alpha2D_R, beta2D, gamma_R, sqrtG_R, gUpper_R, gLower_R);
    Construct_gamma_sqrtG_GUpper(alpha2D, beta2D_U, gamma_U, sqrtG_U, gUpper_U, gLower_U);
    Construct_gamma_sqrtG_GUpper(alpha2D, beta2D_D, gamma_D, sqrtG_D, gUpper_D, gLower_D);

    for (int p = 0; p < 6; p++) {
        if (p == 4) {
            Construct_p4_lonlat_xy_AIA(p, alpha2D, beta2D, gamma, csswm[p].lon, csswm[p].lat, csswm[p].x, csswm[p].y, csswm[p].A, csswm[p].IA);
            Construct_p4_lonlat_xy_AIA(p, alpha2D_L, beta2D, gamma_L, csswm[p].lon_L, csswm[p].lat_L, csswm[p].x_L, csswm[p].y_L, csswm[p].A_L, csswm[p].IA_L);
            Construct_p4_lonlat_xy_AIA(p, alpha2D_R, beta2D, gamma_R, csswm[p].lon_R, csswm[p].lat_R, csswm[p].x_R, csswm[p].y_R, csswm[p].A_R, csswm[p].IA_R);
            Construct_p4_lonlat_xy_AIA(p, alpha2D, beta2D_U, gamma_U, csswm[p].lon_U, csswm[p].lat_U, csswm[p].x_U, csswm[p].y_U, csswm[p].A_U, csswm[p].IA_U);
            Construct_p4_lonlat_xy_AIA(p, alpha2D, beta2D_D, gamma_D, csswm[p].lon_D, csswm[p].lat_D, csswm[p].x_D, csswm[p].y_D, csswm[p].A_D, csswm[p].IA_D);
            continue;
        }
        if (p == 5) {
            Construct_p5_lonlat_xy_AIA(p, alpha2D, beta2D, gamma, csswm[p].lon, csswm[p].lat, csswm[p].x, csswm[p].y, csswm[p].A, csswm[p].IA);
            Construct_p5_lonlat_xy_AIA(p, alpha2D_L, beta2D, gamma_L, csswm[p].lon_L, csswm[p].lat_L, csswm[p].x_L, csswm[p].y_L, csswm[p].A_L, csswm[p].IA_L);
            Construct_p5_lonlat_xy_AIA(p, alpha2D_R, beta2D, gamma_R, csswm[p].lon_R, csswm[p].lat_R, csswm[p].x_R, csswm[p].y_R, csswm[p].A_R, csswm[p].IA_R);
            Construct_p5_lonlat_xy_AIA(p, alpha2D, beta2D_U, gamma_U, csswm[p].lon_U, csswm[p].lat_U, csswm[p].x_U, csswm[p].y_U, csswm[p].A_U, csswm[p].IA_U);
            Construct_p5_lonlat_xy_AIA(p, alpha2D, beta2D_D, gamma_D, csswm[p].lon_D, csswm[p].lat_D, csswm[p].x_D, csswm[p].y_D, csswm[p].A_D, csswm[p].IA_D);
            continue;
        }
        Construct_p0123_lonlat_xy_AIA(p, alpha2D, beta2D, gamma, csswm[p].lon, csswm[p].lat, csswm[p].x, csswm[p].y, csswm[p].A, csswm[p].IA);
        Construct_p0123_lonlat_xy_AIA(p, alpha2D_L, beta2D, gamma_L, csswm[p].lon_L, csswm[p].lat_L, csswm[p].x_L, csswm[p].y_L, csswm[p].A_L, csswm[p].IA_L);
        Construct_p0123_lonlat_xy_AIA(p, alpha2D_R, beta2D, gamma_R, csswm[p].lon_R, csswm[p].lat_R, csswm[p].x_R, csswm[p].y_R, csswm[p].A_R, csswm[p].IA_R);
        Construct_p0123_lonlat_xy_AIA(p, alpha2D, beta2D_U, gamma_U, csswm[p].lon_U, csswm[p].lat_U, csswm[p].x_U, csswm[p].y_U, csswm[p].A_U, csswm[p].IA_U);
        Construct_p0123_lonlat_xy_AIA(p, alpha2D, beta2D_D, gamma_D, csswm[p].lon_D, csswm[p].lat_D, csswm[p].x_D, csswm[p].y_D, csswm[p].A_D, csswm[p].IA_D);
    }

    // delete dynamic array
    delete[] alpha, delete[] beta, delete[] alpha_L, delete[] beta_U, delete[] alpha_R, delete[] beta_D;
    for (int i = 0; i < NX; i++) {
        delete[] alpha2D[i], delete[] beta2D[i];
        delete[] alpha2D_L[i], delete[] alpha2D_R[i];
        delete[] beta2D_U[i], delete[] beta2D_D[i];
    }
    delete[] alpha2D, delete[] beta2D, delete[] alpha2D_L, delete[] alpha2D_R, delete[] beta2D_U, delete[] beta2D_D;
}