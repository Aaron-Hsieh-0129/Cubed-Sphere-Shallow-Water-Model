#include "constrcution.hpp"


void CSSWM::matrixMul(double firstMatrix[4], double secondMatrix[4], double mult[2][2]) {
	double A[2][2], B[2][2];

    // Init
    int count = 0;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			mult[i][j] = 0;
            A[i][j] = firstMatrix[count];
            B[i][j] = secondMatrix[count];
            count++;
		}
	}
    // multiplication
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

double CSSWM::Cube2Sphere_U(CSSWM &model, int p, int i, int j) {
    double mult[2][2];
    if (p == 0 || p == 1 || p == 5) {
        model.matrixMul(model.gUpper_L[i][j], model.csswm[p].A_L[i][j], mult);
        return mult[0][0] * model.csswm[p].u[i][j] + 
               mult[0][1] * 0.25 * (model.csswm[p].v[i][j+1]+model.csswm[p].v[i][j]+model.csswm[p].v[i-1][j+1]+model.csswm[p].v[i-1][j]);
    }
    else if (p == 2 || p == 3) {
        model.matrixMul(model.gUpper_L[i][j], model.csswm[p].A_L[i][j], mult);
        return mult[0][0] * model.csswm[p].u[i][j] +
               mult[0][1] * 0.25 * (model.csswm[p].v[i][j]+model.csswm[p].v[i][j-1]+model.csswm[p].v[i-1][j]+model.csswm[p].v[i-1][j-1]);
    }
    else {
        model.matrixMul(model.gUpper_R[i][j], model.csswm[p].A_R[i][j], mult);
        return mult[0][0] * model.csswm[p].u[i][j] +
               mult[0][1] * 0.25 * (model.csswm[p].v[i+1][j+1]+model.csswm[p].v[i+1][j]+model.csswm[p].v[i][j+1]+model.csswm[p].v[i][j]);
    }
}

double CSSWM::Cube2Sphere_V(CSSWM &model, int p, int i, int j) {
    double mult[2][2];
    if (p == 0 || p == 1 || p == 5) {
        matrixMul(model.csswm[p].A_D[i][j], model.gUpper_D[i][j], mult);
        return mult[1][0] * 0.25 * (model.csswm[p].u[i+1][j]+model.csswm[p].u[i+1][j-1]+model.csswm[p].u[i][j]+model.csswm[p].u[i][j-1]) +
               mult[1][1] * model.csswm[p].v[i][j];
    }
    else if (p == 2 || p == 3) {
        matrixMul(model.csswm[p].A_U[i][j], model.gUpper_U[i][j], mult);
        return mult[1][0] * 0.25 * (model.csswm[p].u[i+1][j+1]+model.csswm[p].u[i+1][j]+model.csswm[p].u[i][j+1]+model.csswm[p].u[i][j]) +
               mult[1][1] * model.csswm[p].v[i][j];
    }
    else {
        matrixMul(model.csswm[p].A_D[i][j], model.gUpper_D[i][j], mult);
        return mult[1][0] * 0.25 * (model.csswm[p].u[i][j]+model.csswm[p].u[i][j-1]+model.csswm[p].u[i-1][j]+model.csswm[p].u[i-1][j-1]) +
               mult[1][1] * model.csswm[p].v[i][j];
    }
}

double CSSWM::Sphere2Cube_U(CSSWM &model, int p, int i, int j) {
    double mult[2][2];
    if (p == 0 || p == 1 || p == 5) {
        matrixMul(model.gLower_L[i][j], model.csswm[p].IA_L[i][j], mult);
        return mult[0][0] * model.csswm[p].u[i][j] +
               mult[0][1] * 0.25 * (model.csswm[p].v[i][j+1]+model.csswm[p].v[i][j]+model.csswm[p].v[i-1][j+1]+model.csswm[p].v[i-1][j]);
    }
    else if (p == 2 || p == 3) {
        matrixMul(model.gLower_L[i][j], model.csswm[p].IA_L[i][j], mult);
        return mult[0][0] * model.csswm[p].u[i][j] +
               mult[0][1] * 0.25 * (model.csswm[p].v[i][j]+model.csswm[p].v[i][j-1]+model.csswm[p].v[i-1][j]+model.csswm[p].v[i-1][j-1]);
    }
    else {
        matrixMul(model.gLower_R[i][j], model.csswm[p].IA_R[i][j], mult);
        return mult[0][0] * model.csswm[p].u[i][j] +
               mult[0][1] * 0.25 * (model.csswm[p].v[i+1][j+1]+model.csswm[p].v[i+1][j]+model.csswm[p].v[i][j+1]+model.csswm[p].v[i][j]);
    }
}

double CSSWM::Sphere2Cube_V(CSSWM &model, int p, int i, int j) {
    double mult[2][2];
    if (p == 0 || p == 1 || p == 5) {
        matrixMul(model.gLower_D[i][j], model.csswm[p].IA_D[i][j], mult);
        return mult[1][0] * 0.25 * (model.csswm[p].u[i+1][j]+model.csswm[p].u[i+1][j-1]+model.csswm[p].u[i][j]+model.csswm[p].u[i][j-1]) +
               mult[1][1] * model.csswm[p].v[i][j];
    }
    else if (p == 2 || p == 3) {
        matrixMul(model.gLower_U[i][j], model.csswm[p].IA_U[i][j], mult);
        return mult[1][0] * 0.25 * (model.csswm[p].u[i+1][j+1]+model.csswm[p].u[i+1][j]+model.csswm[p].u[i][j+1]+model.csswm[p].u[i][j]) +
               mult[1][1] * model.csswm[p].v[i][j];
    }
    else {
        matrixMul(model.gLower_D[i][j], model.csswm[p].IA_D[i][j], mult);
        return mult[1][0] * 0.25 * (model.csswm[p].u[i][j]+model.csswm[p].u[i][j-1]+model.csswm[p].u[i-1][j]+model.csswm[p].u[i-1][j-1]) +
               mult[1][1] * model.csswm[p].v[i][j];
    }
}
/*
double CSSWM::Cube2Cube_U(CSSWM &model, int p1, int p2, int i1, int i2, int j1, int j2) {
    double mult[2][2], A[2][2], B[2][2];
    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mult[i][j] = A[i][j] = B[i][j] = 0.;
        }
    }

    if (p1 == 4) {
        matrixMul(model.gLower_R[i1][j1], model.csswm[p1].IA_R[i1][j1], A);
        matrixMul(model.csswm[p2].A_R[i2][j2], model.gUpper_R[i2][j2], B);
    }
    else {
        matrixMul(model.gLower_L[i1][j1], model.csswm[p1].IA_L[i1][j1], A);
        matrixMul(model.csswm[p2].A_L[i2][j2], model.gUpper_L[i2][j2], B);
    }

    // multiply A & B
    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += A[i][k] * B[k][j];
			}
		}
	}

    // TODO: complete this function
}
*/