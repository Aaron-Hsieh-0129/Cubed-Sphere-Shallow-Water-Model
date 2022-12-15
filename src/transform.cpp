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

double CSSWM::Cube2Cube_U(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2) {
    // p1 is other patch and p2 is the patch who needs other's information
    double mult[2][2], A[2][2], B[2][2];
    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mult[i][j] = A[i][j] = B[i][j] = 0.;
        }
    }

    if (p1 == 4) matrixMul(model.gLower_R[i1][j1], model.csswm[p1].IA_R[i1][j1], A);
    else  matrixMul(model.gLower_L[i1][j1], model.csswm[p1].IA_L[i1][j1], A);

    if (p2 == 4)  matrixMul(model.csswm[p2].A_R[i2][j2], model.gUpper_R[i2][j2], B);
    else matrixMul(model.csswm[p2].A_L[i2][j2], model.gUpper_L[i2][j2], B);
    

    // multiply A & B
    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += A[i][k] * B[k][j];
			}
		}
	}

    if (p2 == 0 || p2 == 1 || p2 == 5) {
        return mult[0][0] * model.csswm[p2].up[i2][j2] + 
               mult[0][1] * 0.25 * (model.csswm[p2].vp[i2][j2+1] + model.csswm[p2].vp[i2][j2] + model.csswm[p2].vp[i2-1][j2+1] + model.csswm[p2].vp[i2-1][j2]);
    }
    else if (p2 == 2 || p2 == 3) {
        return mult[0][0] * model.csswm[p2].up[i2][j2] + 
               mult[0][1] * 0.25 * (model.csswm[p2].vp[i2][j2] + model.csswm[p2].vp[i2][j2-1] + model.csswm[p2].vp[i2-1][j2] + model.csswm[p2].vp[i2-1][j2-1]);
    }
    else {
        return mult[0][0] * model.csswm[p2].up[i2][j2] + 
               mult[0][1] * 0.25 * (model.csswm[p2].vp[i2+1][j2+1] + model.csswm[p2].vp[i2+1][j2] + model.csswm[p2].vp[i2][j2+1] + model.csswm[p2].vp[i2][j2]);
    }
}

double CSSWM::Cube2Cube_V(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2) {
    // p1 is other patch and p2 is the patch who needs other's information
    double mult[2][2], A[2][2], B[2][2];
    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mult[i][j] = A[i][j] = B[i][j] = 0.;
        }
    }

    if (p1 == 2 || p1 == 3) matrixMul(model.gLower_U[i1][j1], model.csswm[p1].IA_U[i1][j1], A);
    else  matrixMul(model.gLower_D[i1][j1], model.csswm[p1].IA_D[i1][j1], A);

    if (p2 == 2 || p2 == 3)  matrixMul(model.csswm[p2].A_U[i2][j2], model.gUpper_U[i2][j2], B);
    else matrixMul(model.csswm[p2].A_D[i2][j2], model.gUpper_D[i2][j2], B);
    

    // multiply A & B
    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += A[i][k] * B[k][j];
			}
		}
	}

    if (p2 == 0 || p2 == 1 || p2 == 5) {
        return mult[1][0] * 0.25 * (model.csswm[p2].up[i2+1][j2] + model.csswm[p2].up[i2][j2] + model.csswm[p2].up[i2+1][j2-1] + model.csswm[p2].up[i2][j2-1]) + 
               mult[1][1] * model.csswm[p2].vp[i2][j2];
    }
    else if (p2 == 2 || p2 == 3) {
        return mult[1][0] * 0.25 * (model.csswm[p2].up[i2+1][j2+1] + model.csswm[p2].up[i2][j2+1] + model.csswm[p2].up[i2+1][j2] + model.csswm[p2].up[i2][j2]) + 
               mult[1][1] * model.csswm[p2].vp[i2][j2];
    }
    else {
        return mult[1][0] * 0.25 * (model.csswm[p2].up[i2][j2] + model.csswm[p2].up[i2][j2-1] + model.csswm[p2].up[i2-1][j2] + model.csswm[p2].up[i2-1][j2-1]) + 
               mult[1][1] * model.csswm[p2].vp[i2][j2];
    }
}

double CSSWM::Cube2Cube_BV2AU(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2) {
    // p1 is other patch and p2 is the patch who needs other's information
    double mult[2][2], A[2][2], B[2][2];
    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mult[i][j] = A[i][j] = B[i][j] = 0.;
        }
    }

    if (p1 == 4) matrixMul(model.gLower_R[i1][j1], model.csswm[p1].IA_R[i1][j1], A);
    else  matrixMul(model.gLower_L[i1][j1], model.csswm[p1].IA_L[i1][j1], A);

    if (p2 == 2 || p2 == 3)  matrixMul(model.csswm[p2].A_U[i2][j2], model.gUpper_U[i2][j2], B);
    else matrixMul(model.csswm[p2].A_D[i2][j2], model.gUpper_D[i2][j2], B);
    

    // multiply A & B
    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += A[i][k] * B[k][j];
			}
		}
	}

    if (p2 == 0 || p2 == 1 || p2 == 5) {
        return mult[0][0] * 0.25 * (model.csswm[p2].up[i2+1][j2]+model.csswm[p2].up[i2][j2]+model.csswm[p2].up[i2+1][j2-1]+model.csswm[p2].up[i2][j2-1]) + 
               mult[0][1] * model.csswm[p2].vp[i2][j2];
    }
    else if (p2 == 2 || p2 == 3) {
        return mult[0][0] * 0.25 * (model.csswm[p2].up[i2+1][j2+1]+model.csswm[p2].up[i2+1][j2]+model.csswm[p2].up[i2][j2+1]+model.csswm[p2].up[i2][j2]) + 
               mult[0][1] * model.csswm[p2].vp[i2][j2];
    }
    else {
        return mult[0][0] * 0.25 * (model.csswm[p2].up[i2][j2]+model.csswm[p2].up[i2][j2-1]+model.csswm[p2].up[i2-1][j2]+model.csswm[p2].up[i2-1][j2-1]) + 
               mult[0][1] * model.csswm[p2].vp[i2][j2];
    }
}

double CSSWM::Cube2Cube_BU2AV(CSSWM &model, int p1, int p2, int i1, int j1, int i2, int j2) {
    // p1 is other patch and p2 is the patch who needs other's information
    double mult[2][2], A[2][2], B[2][2];
    // init
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            mult[i][j] = A[i][j] = B[i][j] = 0.;
        }
    }

    if (p1 == 2 || p1 == 3) matrixMul(model.gLower_U[i1][j1], model.csswm[p1].IA_U[i1][j1], A);
    else  matrixMul(model.gLower_D[i1][j1], model.csswm[p1].IA_D[i1][j1], A);

    if (p2 == 4)  matrixMul(model.csswm[p2].A_R[i2][j2], model.gUpper_R[i2][j2], B);
    else matrixMul(model.csswm[p2].A_L[i2][j2], model.gUpper_L[i2][j2], B);
    

    // multiply A & B
    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				mult[i][j] += A[i][k] * B[k][j];
			}
		}
	}

    if (p2 == 0 || p2 == 1 || p2 == 5) {
        return mult[1][0] * model.csswm[p2].up[i2][j2] + 
               mult[1][1] * 0.25 * (model.csswm[p2].vp[i2][j2+1]+model.csswm[p2].vp[i2][j2]+model.csswm[p2].vp[i2-1][j2+1]+model.csswm[p2].vp[i2-1][j2]);
    }
    else if (p2 == 2 || p2 == 3) {
        return mult[1][0] * model.csswm[p2].up[i2][j2] + 
               mult[1][1] * 0.25 * (model.csswm[p2].vp[i2][j2]+model.csswm[p2].vp[i2][j2-1]+model.csswm[p2].vp[i2-1][j2]+model.csswm[p2].vp[i2-1][j2-1]);
    }
    else {
        return mult[1][0] * model.csswm[p2].up[i2][j2] + 
               mult[1][1] * 0.25 * (model.csswm[p2].vp[i2+1][j2+1]+model.csswm[p2].vp[i2+1][j2]+model.csswm[p2].vp[i2][j2+1]+model.csswm[p2].vp[i2][j2]);
    }
}