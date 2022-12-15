#include "iteration.hpp"

void Iteration::ph_pt(CSSWM &model) {
    for (int p = 0; p < 6; p++) {
        double psqrtGHU_px = 0, psqrtGHU_py = 0, dx_for_h = 0, dy_for_h = 0;
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                if (p == 0 || p == 1 || p == 5) {
                    dx_for_h = model.csswm[p].x_L[i+1][j] - model.csswm[p].x_L[i][j];
                    dy_for_h = model.csswm[p].y_D[i][j+1] - model.csswm[p].y_D[i][j];

                    psqrtGHU_px = (1. / (model.sqrtG[i][j] * dx_for_h)) * 
                                  ((model.sqrtG_L[i+1][j] * 0.5*(model.csswm[p].h[i+1][j]+model.csswm[p].h[i][j]) * 
                                    (model.gUpper_L[i+1][j][0] * model.csswm[p].u[i+1][j] + 
                                     model.gUpper_L[i+1][j][1] * 0.25*(model.csswm[p].v[i+1][j+1]+model.csswm[p].v[i][j+1]+model.csswm[p].v[i+1][j]+model.csswm[p].v[i][j]))) - 
                                   (model.sqrtG_L[i][j] * 0.5*(model.csswm[p].h[i][j]+model.csswm[p].h[i-1][j]) * 
                                    (model.gUpper_L[i][j][0] * model.csswm[p].u[i][j] + 
                                     model.gUpper_L[i][j][1] * 0.25*(model.csswm[p].v[i][j+1]+model.csswm[p].v[i][j]+model.csswm[p].v[i-1][j+1]+model.csswm[p].v[i-1][j]))));

                    psqrtGHU_py = (1. / (model.sqrtG[i][j] * dy_for_h)) * 
                                  ((model.sqrtG_D[i][j+1] * 0.5*(model.csswm[p].h[i][j+1]+model.csswm[p].h[i][j]) * 
                                    (model.gUpper_D[i][j+1][3] * model.csswm[p].v[i][j+1] + 
                                     model.gUpper_D[i][j+1][2] * 0.25*(model.csswm[p].u[i+1][j+1]+model.csswm[p].u[i][j+1]+model.csswm[p].u[i+1][j]+model.csswm[p].u[i][j]))) - 
                                   (model.sqrtG_D[i][j] * 0.5*(model.csswm[p].h[i][j]+model.csswm[p].h[i][j-1]) * 
                                    (model.gUpper_D[i][j][3] * model.csswm[p].v[i][j] + 
                                     model.gUpper_D[i][j][2] * 0.25*(model.csswm[p].u[i+1][j]+model.csswm[p].u[i][j]+model.csswm[p].u[i+1][j-1]+model.csswm[p].u[i][j-1]))));
                
                    model.csswm[p].hp[i][j] = model.csswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py);
                    
                    #ifdef DIFFUSION
                        model.csswm[p].hp[i][j] += D2T * Kx * (model.csswm[p].hm[i+1][j] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i-1][j]) / pow(dx_for_h, 2) + 
                                                   D2T * Ky  * (model.csswm[p].hm[i][j+1] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i][j-1]) / pow(dy_for_h, 2);
                    #endif
                }
                else if (p == 2 || p == 3) {
                    dx_for_h = model.csswm[p].x_L[i+1][j] - model.csswm[p].x_L[i][j];
                    dy_for_h = model.csswm[p].y_U[i][j] - model.csswm[p].y_U[i][j-1];

                    psqrtGHU_px = (1. / (model.sqrtG[i][j] * dx_for_h)) * 
                                  ((model.sqrtG_L[i+1][j] * 0.5*(model.csswm[p].h[i+1][j]+model.csswm[p].h[i][j]) * 
                                    (model.gUpper_L[i+1][j][0] * model.csswm[p].u[i+1][j] + 
                                     model.gUpper_L[i+1][j][1] * 0.25*(model.csswm[p].v[i+1][j]+model.csswm[p].v[i+1][j-1]+model.csswm[p].v[i][j]+model.csswm[p].v[i][j-1]))) - 
                                   (model.sqrtG_L[i][j] * 0.5*(model.csswm[p].h[i][j]+model.csswm[p].h[i-1][j]) * 
                                    (model.gUpper_L[i][j][0] * model.csswm[p].u[i][j] + 
                                     model.gUpper_L[i][j][1] * 0.25*(model.csswm[p].v[i][j]+model.csswm[p].v[i][j-1]+model.csswm[p].v[i-1][j]+model.csswm[p].v[i-1][j-1]))));

                    psqrtGHU_py = (1. / (model.sqrtG[i][j] * dy_for_h)) * 
                                  ((model.sqrtG_U[i][j] * 0.5*(model.csswm[p].h[i][j+1]+model.csswm[p].h[i][j]) * 
                                    (model.gUpper_U[i][j][3] * model.csswm[p].v[i][j] + 
                                    model.gUpper_U[i][j][2] * 0.25*(model.csswm[p].u[i+1][j+1]+model.csswm[p].u[i+1][j]+model.csswm[p].u[i][j+1]+model.csswm[p].u[i][j]))) - 
                                   (model.sqrtG_U[i][j-1] * 0.5*(model.csswm[p].h[i][j]+model.csswm[p].h[i][j-1]) * 
                                    (model.gUpper_U[i][j-1][3] * model.csswm[p].v[i][j-1] + 
                                     model.gUpper_U[i][j-1][2] * 0.25*(model.csswm[p].u[i+1][j]+model.csswm[p].u[i+1][j-1]+model.csswm[p].u[i][j]+model.csswm[p].u[i][j-1]))));

                    model.csswm[p].hp[i][j] = model.csswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py);

                    #ifdef DIFFUSION
                        model.csswm[p].hp[i][j] += D2T * Kx * (model.csswm[p].hm[i+1][j] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i-1][j]) / pow(dx_for_h, 2) + 
                                                   D2T * Ky  * (model.csswm[p].hm[i][j+1] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i][j-1]) / pow(dy_for_h, 2);
                    #endif
                }
                else {
                    dx_for_h = model.csswm[p].x_R[i][j] - model.csswm[p].x_R[i-1][j];
                    dy_for_h = model.csswm[p].y_D[i][j+1] - model.csswm[p].y_D[i][j];

                    psqrtGHU_px = (1. / (model.sqrtG[i][j] * dx_for_h)) * 
                                  ((model.sqrtG_R[i][j] * 0.5*(model.csswm[p].h[i+1][j]+model.csswm[p].h[i][j]) * 
                                    (model.gUpper_R[i][j][0] * model.csswm[p].u[i][j] + 
                                     model.gUpper_R[i][j][1] * 0.25*(model.csswm[p].v[i+1][j+1]+model.csswm[p].v[i+1][j]+model.csswm[p].v[i][j+1]+model.csswm[p].v[i][j]))) - 
                                   (model.sqrtG_R[i-1][j] * 0.5*(model.csswm[p].h[i][j]+model.csswm[p].h[i-1][j]) * 
                                    (model.gUpper_R[i-1][j][0] * model.csswm[p].u[i-1][j] + 
                                     model.gUpper_R[i-1][j][1] * 0.25*(model.csswm[p].v[i][j+1]+model.csswm[p].v[i][j]+model.csswm[p].v[i-1][j+1]+model.csswm[p].v[i-1][j]))));

                    psqrtGHU_py = (1. / (model.sqrtG[i][j] * dy_for_h)) * 
                                  ((model.sqrtG_D[i][j+1] * 0.5*(model.csswm[p].h[i][j+1]+model.csswm[p].h[i][j]) * 
                                    (model.gUpper_D[i][j+1][3] * model.csswm[p].v[i][j+1] + 
                                    model.gUpper_D[i][j+1][2] * 0.25*(model.csswm[p].u[i][j+1]+model.csswm[p].u[i][j]+model.csswm[p].u[i-1][j+1]+model.csswm[p].u[i-1][j]))) - 
                                   (model.sqrtG_D[i][j] * 0.5*(model.csswm[p].h[i][j]+model.csswm[p].h[i][j-1]) * 
                                    (model.gUpper_D[i][j][3] * model.csswm[p].v[i][j] + 
                                     model.gUpper_D[i][j][2] * 0.25*(model.csswm[p].u[i][j]+model.csswm[p].u[i][j-1]+model.csswm[p].u[i-1][j]+model.csswm[p].u[i-1][j-1]))));

                    model.csswm[p].hp[i][j] = model.csswm[p].hm[i][j] + D2T * (-psqrtGHU_px - psqrtGHU_py);

                    // diffusion
                    #ifdef DIFFUSION
                        model.csswm[p].hp[i][j] += D2T * Kx * (model.csswm[p].hm[i+1][j] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i-1][j]) / pow(dx_for_h, 2) + 
                                                   D2T * Ky * (model.csswm[p].hm[i][j+1] - 2. * model.csswm[p].hm[i][j] + model.csswm[p].hm[i][j-1]) / pow(dy_for_h, 2);
                    #endif
                }
            }
        }
    }
    return;
}

void Iteration::pu_pt(CSSWM &model) {
    double dx_for_u = 0, dy_for_u = 0;
    double pgH_px = 0, pU2_px = 0, pUV_px = 0, pV2_px = 0, rotationU = 0;

    for (int p = 0; p < 6; p++) {
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                if (p == 0 || p == 1 || p == 5) {
                    dx_for_u = model.csswm[p].x[i][j] - model.csswm[p].x[i-1][j];
                    dy_for_u = 0.5 * (model.csswm[p].y_L[i][j+1] - model.csswm[p].y_L[i][j-1]);

                    pgH_px = GRAVITY / dx_for_u * (model.csswm[p].h[i][j] - model.csswm[p].h[i-1][j]);

                    pU2_px = 0.5 / dx_for_u * 
                            (model.gUpper[i][j][0] * pow(0.5*(model.csswm[p].u[i+1][j] + model.csswm[p].u[i][j]), 2) -
                             model.gUpper[i-1][j][0] * pow(0.5*(model.csswm[p].u[i][j] + model.csswm[p].u[i-1][j]), 2));

                    pUV_px = 1 / dx_for_u * 
                         ((model.gUpper[i][j][1] * 0.5*(model.csswm[p].u[i+1][j]+model.csswm[p].u[i][j]) * 0.5*(model.csswm[p].v[i][j+1]+model.csswm[p].v[i][j])) - 
                          (model.gUpper[i-1][j][1] * 0.5*(model.csswm[p].u[i][j]+model.csswm[p].u[i-1][j]) * 0.5*(model.csswm[p].v[i-1][j+1]+model.csswm[p].v[i-1][j])));

                    pV2_px = 0.5 / dx_for_u * 
                            (model.gUpper[i][j][3] * pow(0.5*(model.csswm[p].v[i][j+1] + model.csswm[p].v[i][j]), 2) -
                             model.gUpper[i-1][j][3] * pow(0.5*(model.csswm[p].v[i-1][j+1] + model.csswm[p].v[i-1][j]), 2));

                    rotationU = (((0.5*(model.csswm[p].v[i][j+1] + model.csswm[p].v[i][j]) - 0.5*(model.csswm[p].v[i-1][j+1] + model.csswm[p].v[i-1][j])) / dx_for_u) - 
                                 (0.5*(model.csswm[p].u[i][j+1] - model.csswm[p].u[i][j-1]) / dy_for_u) + model.sqrtG_L[i][j] * CF) * 
                                (model.gUpper_L[i][j][2] * model.csswm[p].u[i][j] + 
                                 model.gUpper_L[i][j][3] * 0.25*(model.csswm[p].v[i][j+1] + model.csswm[p].v[i][j] + model.csswm[p].v[i-1][j+1] + model.csswm[p].v[i-1][j]));
                
                    model.csswm[p].up[i][j] = model.csswm[p].um[i][j] + D2T * (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);

                    #ifdef DIFFUSION
                        model.csswm[p].up[i][j] += D2T * Kx * (model.csswm[p].um[i+1][j] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i-1][j]) / pow(dx_for_u, 2) + 
                                                   D2T * Ky  * (model.csswm[p].um[i][j+1] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i][j-1]) / pow(dy_for_u, 2);
                    #endif
                }
                else if (p == 2 || p == 3) {
                    dx_for_u = model.csswm[p].x[i][j] - model.csswm[p].x[i-1][j];
                    dy_for_u = 0.5 * (model.csswm[p].y_L[i][j+1] - model.csswm[p].y_L[i][j-1]);

                    pgH_px = GRAVITY / dx_for_u * (model.csswm[p].h[i][j] - model.csswm[p].h[i-1][j]);

                    pU2_px = 0.5 / dx_for_u * 
                            (model.gUpper[i][j][0] * pow(0.5*(model.csswm[p].u[i+1][j] + model.csswm[p].u[i][j]), 2) -
                             model.gUpper[i-1][j][0] * pow(0.5*(model.csswm[p].u[i][j] + model.csswm[p].u[i-1][j]), 2));

                    pUV_px = 1 / dx_for_u * 
                         ((model.gUpper[i][j][1] * 0.5*(model.csswm[p].u[i+1][j]+model.csswm[p].u[i][j]) * 0.5*(model.csswm[p].v[i][j]+model.csswm[p].v[i][j-1])) - 
                          (model.gUpper[i-1][j][1] * 0.5*(model.csswm[p].u[i][j]+model.csswm[p].u[i-1][j]) * 0.5*(model.csswm[p].v[i-1][j]+model.csswm[p].v[i-1][j-1])));

                    pV2_px = 0.5 / dx_for_u * 
                            (model.gUpper[i][j][3] * pow(0.5*(model.csswm[p].v[i][j] + model.csswm[p].v[i][j-1]), 2) -
                             model.gUpper[i-1][j][3] * pow(0.5*(model.csswm[p].v[i-1][j] + model.csswm[p].v[i-1][j-1]), 2));

                    rotationU = (((0.5*(model.csswm[p].v[i][j] + model.csswm[p].v[i][j-1]) - 0.5*(model.csswm[p].v[i-1][j] + model.csswm[p].v[i-1][j-1])) / dx_for_u) - 
                                 (0.5*(model.csswm[p].u[i][j+1] - model.csswm[p].u[i][j-1]) / dy_for_u) + model.sqrtG_L[i][j] * CF) * 
                                (model.gUpper_L[i][j][2] * model.csswm[p].u[i][j] + 
                                 model.gUpper_L[i][j][3] * 0.25*(model.csswm[p].v[i][j] + model.csswm[p].v[i][j-1] + model.csswm[p].v[i-1][j] + model.csswm[p].v[i-1][j-1]));

                    model.csswm[p].up[i][j] = model.csswm[p].um[i][j] + D2T * (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);

                    #ifdef DIFFUSION
                        model.csswm[p].up[i][j] += D2T * Kx * (model.csswm[p].um[i+1][j] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i-1][j]) / pow(dx_for_u, 2) + 
                                                  D2T * Ky  * (model.csswm[p].um[i][j+1] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i][j-1]) / pow(dy_for_u, 2);
                    #endif
                }
                else {
                    dx_for_u = model.csswm[p].x[i+1][j] - model.csswm[p].x[i][j];
                    dy_for_u = 0.5 * (model.csswm[p].y_R[i][j+1] - model.csswm[p].y_R[i][j-1]);

                    pgH_px = GRAVITY / dx_for_u * (model.csswm[p].h[i+1][j] - model.csswm[p].h[i][j]);

                    pU2_px = 0.5 / dx_for_u * 
                            (model.gUpper[i+1][j][0] * pow(0.5*(model.csswm[p].u[i+1][j] + model.csswm[p].u[i][j]), 2) -
                             model.gUpper[i][j][0] * pow(0.5*(model.csswm[p].u[i][j] + model.csswm[p].u[i-1][j]), 2));

                    pUV_px = 1 / dx_for_u * 
                         ((model.gUpper[i+1][j][1] * 0.5*(model.csswm[p].u[i+1][j]+model.csswm[p].u[i][j]) * 0.5*(model.csswm[p].v[i+1][j+1]+model.csswm[p].v[i+1][j])) - 
                          (model.gUpper[i][j][1] * 0.5*(model.csswm[p].u[i][j]+model.csswm[p].u[i-1][j]) * 0.5*(model.csswm[p].v[i][j+1]+model.csswm[p].v[i][j])));

                    pV2_px = 0.5 / dx_for_u * 
                            (model.gUpper[i+1][j][3] * pow(0.5*(model.csswm[p].v[i+1][j+1] + model.csswm[p].v[i+1][j]), 2) -
                             model.gUpper[i][j][3] * pow(0.5*(model.csswm[p].v[i][j+1] + model.csswm[p].v[i][j]), 2));

                    rotationU = (((0.5*(model.csswm[p].v[i+1][j+1] + model.csswm[p].v[i+1][j]) - 0.5*(model.csswm[p].v[i][j+1] + model.csswm[p].v[i][j])) / dx_for_u) - 
                                 (0.5*(model.csswm[p].u[i][j+1] - model.csswm[p].u[i][j-1]) / dy_for_u) + model.sqrtG_R[i][j] * CF) * 
                                (model.gUpper_R[i][j][2] * model.csswm[p].u[i][j] + 
                                 model.gUpper_R[i][j][3] * 0.25*(model.csswm[p].v[i+1][j+1] + model.csswm[p].v[i+1][j] + model.csswm[p].v[i][j+1] + model.csswm[p].v[i][j]));

                    model.csswm[p].up[i][j] = model.csswm[p].um[i][j] + D2T * (-pgH_px - pU2_px - pUV_px - pV2_px + rotationU);

                    #ifdef DIFFUSION
                        model.csswm[p].up[i][j] += D2T * Kx * (model.csswm[p].um[i+1][j] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i-1][j]) / pow(dx_for_u, 2) + 
                                                   D2T * Ky  * (model.csswm[p].um[i][j+1] - 2. * model.csswm[p].um[i][j] + model.csswm[p].um[i][j-1]) / pow(dy_for_u, 2);
                    #endif
                }
            }
        }
    }
    return;
}

void Iteration::pv_pt(CSSWM &model) {
    double dx_for_v = 0, dy_for_v = 0;
    double pgH_py = 0, pU2_py = 0, pUV_py = 0, pV2_py = 0, rotationV = 0;

    for (int p = 0; p < 6; p++) {
        for (int i = 2; i < NX-2; i++) {
            for (int j = 2; j < NY-2; j++) {
                if (p == 0 || p == 1 || p == 5) {
                    dx_for_v = 0.5*(model.csswm[p].x_D[i+1][j] - model.csswm[p].x_D[i-1][j]);
                    dy_for_v = model.csswm[p].y[i][j] - model.csswm[p].y[i][j-1];

                    pgH_py = GRAVITY / dy_for_v * (model.csswm[p].h[i][j] - model.csswm[p].h[i][j-1]);

                    pU2_py = 0.5 / dy_for_v * (model.gUpper[i][j][0] * pow(0.5*(model.csswm[p].u[i+1][j] + model.csswm[p].u[i][j]), 2) -
                                               model.gUpper[i][j-1][0] * pow(0.5*(model.csswm[p].u[i+1][j-1] + model.csswm[p].u[i][j-1]), 2));

                    pUV_py = 1 / dy_for_v * 
                            ((model.gUpper[i][j][1] * 0.5*(model.csswm[p].u[i+1][j]+model.csswm[p].u[i][j]) * 0.5*(model.csswm[p].v[i][j+1]+model.csswm[p].v[i][j])) - 
                             (model.gUpper[i][j-1][1] * 0.5*(model.csswm[p].u[i+1][j-1]+model.csswm[p].u[i][j-1]) * 0.5*(model.csswm[p].v[i][j]+model.csswm[p].v[i][j-1])));

                    pV2_py = 0.5 / dy_for_v * (model.gUpper[i][j][3] * pow(0.5*(model.csswm[p].v[i][j+1] + model.csswm[p].v[i][j]), 2) -
                                               model.gUpper[i][j-1][3] * pow(0.5*(model.csswm[p].v[i][j] + model.csswm[p].v[i][j-1]), 2));

                    rotationV = ((0.5*(model.csswm[p].v[i+1][j] - model.csswm[p].v[i-1][j]) / dx_for_v) - 
                                 ((0.5*(model.csswm[p].u[i+1][j]+model.csswm[p].u[i][j]) - 0.5*(model.csswm[p].u[i+1][j-1]+model.csswm[p].u[i][j-1])) / dy_for_v) + model.sqrtG_D[i][j] * CF) * 
                                (model.gUpper_D[i][j][1] * model.csswm[p].v[i][j] + 
                                 model.gUpper_D[i][j][0] * 0.25*(model.csswm[p].u[i+1][j] + model.csswm[p].u[i][j] + model.csswm[p].u[i+1][j-1] + model.csswm[p].u[i][j-1]));

                    model.csswm[p].vp[i][j] = model.csswm[p].vm[i][j] + D2T * (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);
                    
                    #ifdef DIFFUSION
                        model.csswm[p].vp[i][j] += D2T * Kx * (model.csswm[p].vm[i+1][j] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i-1][j]) / pow(dx_for_v, 2) + 
                                                   D2T * Ky * (model.csswm[p].vm[i][j+1] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i][j-1]) / pow(dy_for_v, 2);
                    #endif
                }
                else if (p == 2 || p == 3) {
                    dx_for_v = 0.5*(model.csswm[p].x_U[i+1][j] - model.csswm[p].x_U[i-1][j]);
                    dy_for_v = model.csswm[p].y[i][j+1] - model.csswm[p].y[i][j];

                    pgH_py = GRAVITY / dy_for_v * (model.csswm[p].h[i][j+1] - model.csswm[p].h[i][j]);

                    pU2_py = 0.5 / dy_for_v * (model.gUpper[i][j+1][0] * pow(0.5*(model.csswm[p].u[i+1][j+1] + model.csswm[p].u[i][j+1]), 2) -
                                               model.gUpper[i][j][0] * pow(0.5*(model.csswm[p].u[i+1][j] + model.csswm[p].u[i][j]), 2));

                    pUV_py = 1 / dy_for_v * 
                            ((model.gUpper[i][j+1][1] * 0.5*(model.csswm[p].u[i+1][j+1]+model.csswm[p].u[i][j+1]) * 0.5*(model.csswm[p].v[i][j+1]+model.csswm[p].v[i][j])) - 
                             (model.gUpper[i][j][1] * 0.5*(model.csswm[p].u[i+1][j]+model.csswm[p].u[i][j]) * 0.5*(model.csswm[p].v[i][j]+model.csswm[p].v[i][j-1])));

                    pV2_py = 0.5 / dy_for_v * (model.gUpper[i][j+1][3] * pow(0.5*(model.csswm[p].v[i][j+1] + model.csswm[p].v[i][j]), 2) -
                                               model.gUpper[i][j][3] * pow(0.5*(model.csswm[p].v[i][j] + model.csswm[p].v[i][j-1]), 2));

                    rotationV = ((0.5*(model.csswm[p].v[i+1][j] - model.csswm[p].v[i-1][j]) / dx_for_v) - 
                                 ((0.5*(model.csswm[p].u[i+1][j+1]+model.csswm[p].u[i][j+1]) - 0.5*(model.csswm[p].u[i+1][j]+model.csswm[p].u[i][j])) / dy_for_v) + model.sqrtG_U[i][j] * CF) * 
                                (model.gUpper_U[i][j][1] * model.csswm[p].v[i][j] + 
                                 model.gUpper_U[i][j][0] * 0.25*(model.csswm[p].u[i+1][j+1] + model.csswm[p].u[i][j+1] + model.csswm[p].u[i+1][j] + model.csswm[p].u[i][j]));

                    model.csswm[p].vp[i][j] = model.csswm[p].vm[i][j] + D2T * (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);

                    #ifdef DIFFUSION
                        model.csswm[p].vp[i][j] += D2T * Kx * (model.csswm[p].vm[i+1][j] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i-1][j]) / pow(dx_for_v, 2) + 
                                                   D2T * Ky * (model.csswm[p].vm[i][j+1] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i][j-1]) / pow(dy_for_v, 2);
                    #endif
                }
                else {
                    dx_for_v = 0.5*(model.csswm[p].x_D[i+1][j] - model.csswm[p].x_D[i-1][j]);
                    dy_for_v = model.csswm[p].y[i][j] - model.csswm[p].y[i][j-1];

                    pgH_py = GRAVITY / dy_for_v * (model.csswm[p].h[i][j] - model.csswm[p].h[i][j-1]);

                    pU2_py = 0.5 / dy_for_v * (model.gUpper[i][j][0] * pow(0.5*(model.csswm[p].u[i][j] + model.csswm[p].u[i-1][j]), 2) -
                                               model.gUpper[i][j-1][0] * pow(0.5*(model.csswm[p].u[i][j-1] + model.csswm[p].u[i-1][j-1]), 2));

                    pUV_py = 1 / dy_for_v * 
                            ((model.gUpper[i][j][1] * 0.5*(model.csswm[p].u[i][j]+model.csswm[p].u[i-1][j]) * 0.5*(model.csswm[p].v[i][j+1]+model.csswm[p].v[i][j])) - 
                             (model.gUpper[i][j-1][1] * 0.5*(model.csswm[p].u[i][j-1]+model.csswm[p].u[i-1][j-1]) * 0.5*(model.csswm[p].v[i][j]+model.csswm[p].v[i][j-1])));

                    pV2_py = 0.5 / dy_for_v * (model.gUpper[i][j][3] * pow(0.5*(model.csswm[p].v[i][j+1] + model.csswm[p].v[i][j]), 2) -
                                               model.gUpper[i][j-1][3] * pow(0.5*(model.csswm[p].v[i][j] + model.csswm[p].v[i][j-1]), 2));

                    rotationV = ((0.5*(model.csswm[p].v[i+1][j] - model.csswm[p].v[i-1][j]) / dx_for_v) - 
                                 ((0.5*(model.csswm[p].u[i][j]+model.csswm[p].u[i-1][j]) - 0.5*(model.csswm[p].u[i][j-1]+model.csswm[p].u[i-1][j-1])) / dy_for_v) + model.sqrtG_D[i][j] * CF) * 
                                (model.gUpper_D[i][j][1] * model.csswm[p].v[i][j] + 
                                 model.gUpper_D[i][j][0] * 0.25*(model.csswm[p].u[i][j] + model.csswm[p].u[i-1][j] + model.csswm[p].u[i][j-1] + model.csswm[p].u[i-1][j-1]));
                
                    model.csswm[p].vp[i][j] = model.csswm[p].vm[i][j] + D2T * (-pgH_py - pU2_py - pUV_py - pV2_py - rotationV);

                    #ifdef DIFFUSION
                        model.csswm[p].vp[i][j] += D2T * Kx * (model.csswm[p].vm[i+1][j] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i-1][j]) / pow(dx_for_v, 2) + 
                                                   D2T * Ky  * (model.csswm[p].vm[i][j+1] - 2. * model.csswm[p].vm[i][j] + model.csswm[p].vm[i][j-1]) / pow(dy_for_v, 2);
                    #endif
                }
            }
        }
    }
    return;
}

void Iteration::leap_frog(CSSWM &model) {
    Outputs::output_parameter(model);
    int n = 0;
    double timenow = 0.;
    double temp = TIMEEND / DT;
    int nmax = (int) temp;

    while (n < nmax) {
        std::cout << n << std::endl;

        if (n % OUTPUTINTERVAL == 0) {
            Outputs::output_h(n, model);
            Outputs::output_u(n, model);
            Outputs::output_v(n, model);
        }

        n++;
        timenow = n * DT;

        // calculate
        ph_pt(model);
        // pu_pt(model);
        // pv_pt(model);

        // TODO: boundary process
        model.BP_h(model);
        model.BP_h2(model);

        // Time filter
        #ifdef TIMEFILTER
            for (int p = 0; p < 6; p++) {
                for (int i = 0; i < NX; i++) {
                    for (int j = 0; j < NY; j++) {
                        model.csswm[p].h[i][j] += TIMETS * (model.csswm[p].hp[i][j] - 2 * model.csswm[p].h[i][j] + model.csswm[p].hm[i][j]);
                        model.csswm[p].u[i][j] += TIMETS * (model.csswm[p].up[i][j] - 2 * model.csswm[p].u[i][j] + model.csswm[p].um[i][j]);
                        model.csswm[p].v[i][j] += TIMETS * (model.csswm[p].vp[i][j] - 2 * model.csswm[p].v[i][j] + model.csswm[p].vm[i][j]);
                    }
                }
            }
        #endif

        // next step
        for (int p = 0; p < 6; p++) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    model.csswm[p].hm[i][j] = model.csswm[p].h[i][j];
                    model.csswm[p].h[i][j] = model.csswm[p].hp[i][j];

                    model.csswm[p].um[i][j] = model.csswm[p].u[i][j];
                    model.csswm[p].u[i][j] = model.csswm[p].up[i][j];

                    model.csswm[p].vm[i][j] = model.csswm[p].v[i][j];
                    model.csswm[p].v[i][j] = model.csswm[p].vp[i][j];
                }
            }
        }
    }
    
    return;



}