#include "initialField.hpp"

void Init::Init2d(CSSWM & model) {
    for (int p = 0; p < 6; p++) {
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                #ifdef Jung
                    model.csswm[p].hp[i][j] = JungH(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    double mult[2][2];
                    if (p == 0 || p == 1 || p == 5) {
                        model.matrixMul(model.gLower_L[i][j], model.csswm[p].IA_L[i][j], mult);
                        model.csswm[p].up[i][j] = mult[0][0] * JungU(model.csswm[p].lon_L[i][j], model.csswm[p].lat_L[i][j]) + 
                                                  mult[0][1] * JungV(model.csswm[p].lon_L[i][j]);

                        model.matrixMul(model.gLower_D[i][j], model.csswm[p].IA_D[i][j], mult);
                        model.csswm[p].vp[i][j] = mult[1][0] * JungU(model.csswm[p].lon_D[i][j], model.csswm[p].lat_D[i][j]) + 
                                                  mult[1][1] * JungV(model.csswm[p].lon_D[i][j]);
                    }
                    else if (p == 2 || p == 3) {
                        model.matrixMul(model.gLower_L[i][j], model.csswm[p].IA_L[i][j], mult);
                        model.csswm[p].up[i][j] = mult[0][0] * JungU(model.csswm[p].lon_L[i][j], model.csswm[p].lat_L[i][j]) + 
                                                  mult[0][1] * JungV(model.csswm[p].lon_L[i][j]);
                        model.matrixMul(model.gLower_U[i][j], model.csswm[p].IA_U[i][j], mult);                
                        model.csswm[p].vp[i][j] = mult[1][0] * JungU(model.csswm[p].lon_U[i][j], model.csswm[p].lat_U[i][j]) + 
                                                  mult[1][1] * JungV(model.csswm[p].lon_U[i][j]);
                    }
                    else {
                        model.matrixMul(model.gLower_R[i][j], model.csswm[p].IA_R[i][j], mult);
                        model.csswm[p].up[i][j] = mult[0][0] * JungU(model.csswm[p].lon_R[i][j], model.csswm[p].lat_R[i][j]) + 
                                                  mult[0][1] * JungV(model.csswm[p].lon_R[i][j]);
                        model.matrixMul(model.gLower_D[i][j], model.csswm[p].IA_D[i][j], mult);
                        model.csswm[p].vp[i][j] = mult[1][0] * JungU(model.csswm[p].lon_D[i][j], model.csswm[p].lat_D[i][j]) + 
                                                  mult[1][1] * JungV(model.csswm[p].lon_D[i][j]);
                    }
                #endif

                #ifdef GravityWave
                    model.csswm[p].hp[i][j] = JungH(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]) / 100.;
                    model.csswm[p].up[i][j] = 0.;
                    model.csswm[p].vp[i][j] = 0.;
                #endif

                #ifdef SteadyGeostrophy
                    model.csswm[p].hp[i][j] = JungH(model.csswm[p].lon[i][j], model.csswm[p].lat[i][j]);
                    if (p == 0 || p == 1 || p == 5) {
                        model.csswm[p].up[i][j] = (model.gLower_L[i][j][0] * model.csswm[p].IA_L[i][j][0] + model.gLower_L[i][j][1] * model.csswm[p].IA_L[i][j][2]) * JungU(model.csswm[p].lon_L[i][j], model.csswm[p].lat_L[i][j]) + 
                                                  (model.gLower_L[i][j][0] * model.csswm[p].IA_L[i][j][1] + model.gLower_L[i][j][1] * model.csswm[p].IA_L[i][j][3]) * JungV(model.csswm[p].lon_L[i][j]);
                        model.csswm[p].vp[i][j] = (model.gLower_D[i][j][2] * model.csswm[p].IA_D[i][j][0] + model.gLower_D[i][j][3] * model.csswm[p].IA_D[i][j][2]) * JungU(model.csswm[p].lon_D[i][j], model.csswm[p].lat_D[i][j]) + 
                                                  (model.gLower_D[i][j][2] * model.csswm[p].IA_D[i][j][1] + model.gLower_D[i][j][3] * model.csswm[p].IA_D[i][j][3]) * JungV(model.csswm[p].lon_D[i][j]);
                    }
                    else if (p == 2 || p == 3) {
                        model.csswm[p].up[i][j] = (model.gLower_L[i][j][0] * model.csswm[p].IA_L[i][j][0] + model.gLower_L[i][j][1] * model.csswm[p].IA_L[i][j][2]) * JungU(model.csswm[p].lon_L[i][j], model.csswm[p].lat_L[i][j]) + 
                                                  (model.gLower_L[i][j][0] * model.csswm[p].IA_L[i][j][1] + model.gLower_L[i][j][1] * model.csswm[p].IA_L[i][j][3]) * JungV(model.csswm[p].lon_L[i][j]);
                        model.csswm[p].vp[i][j] = (model.gLower_U[i][j][2] * model.csswm[p].IA_U[i][j][0] + model.gLower_U[i][j][3] * model.csswm[p].IA_U[i][j][2]) * JungU(model.csswm[p].lon_U[i][j], model.csswm[p].lat_U[i][j]) + 
                                                  (model.gLower_U[i][j][2] * model.csswm[p].IA_U[i][j][1] + model.gLower_U[i][j][3] * model.csswm[p].IA_U[i][j][3]) * JungV(model.csswm[p].lon_U[i][j]);
                    }
                    else {
                        model.csswm[p].up[i][j] = (model.gLower_R[i][j][0] * model.csswm[p].IA_R[i][j][0] + model.gLower_R[i][j][1] * model.csswm[p].IA_R[i][j][2]) * JungU(model.csswm[p].lon_R[i][j], model.csswm[p].lat_R[i][j]) + 
                                                  (model.gLower_R[i][j][0] * model.csswm[p].IA_R[i][j][1] + model.gLower_R[i][j][1] * model.csswm[p].IA_R[i][j][3]) * JungV(model.csswm[p].lon_R[i][j]);
                        model.csswm[p].vp[i][j] = (model.gLower_D[i][j][2] * model.csswm[p].IA_D[i][j][0] + model.gLower_D[i][j][3] * model.csswm[p].IA_D[i][j][2]) * JungU(model.csswm[p].lon_D[i][j], model.csswm[p].lat_D[i][j]) + 
                                                  (model.gLower_D[i][j][2] * model.csswm[p].IA_D[i][j][1] + model.gLower_D[i][j][3] * model.csswm[p].IA_D[i][j][3]) * JungV(model.csswm[p].lon_D[i][j]);
                    }
                #endif
            }
        }
    }

    for (int p = 0; p < 6; p++) {
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                model.csswm[p].hm[i][j] = model.csswm[p].hp[i][j]; 
                model.csswm[p].um[i][j] = model.csswm[p].up[i][j]; 
                model.csswm[p].vm[i][j] = model.csswm[p].vp[i][j];   

                model.csswm[p].h[i][j] = model.csswm[p].hp[i][j]; 
                model.csswm[p].u[i][j] = model.csswm[p].up[i][j]; 
                model.csswm[p].v[i][j] = model.csswm[p].vp[i][j];  
            }
        }
    }
}

double Init::JungH(double lon, double lat) {
    double h0 = 1000;
    double lonC = 0., latC = 0.;
    double rd = RADIUS * acos(sin(latC) * sin(lat) + cos(latC) * cos(lat) * cos(lon-lonC));
    double r0 = RADIUS / 3.;
    if (rd < r0) return h0 / 2. * (1 + cos(M_PI * rd / r0));
    else return 0.;
}

double Init::JungU(double lon, double lat) {
    double u0 = 2 * M_PI * RADIUS / (12. * 86400);
    double u = u0 * (cos(ALPHA0) * cos(lat) + sin(ALPHA0) * sin(lon) * sin(lat));
    return u;
}

double Init::JungV(double lon) {
    double u0 = 2 * M_PI * RADIUS / (12. * 86400);
    double v = u0 * sin(ALPHA0) * cos(lon);
    return v;
}

double Init::SteadyGeostrophyH(double lon, double lat) {
    double h0 = 2.94E4 / GRAVITY;
    double u0 = 2 * M_PI * RADIUS / (12. * 86400);
    return h0 - (RADIUS * OMEGA * u0 + u0 * u0 / 2.) * pow(-cos(lon) * cos(lat) * sin(ALPHA0) + sin(lat) * cos(ALPHA0), 2) / GRAVITY;
}