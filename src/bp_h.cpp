#include "constrcution.hpp"

void CSSWM::BP_h(CSSWM &model) {
    for (int idx = 0; idx < NX; idx++) {
        if (idx < (int) NX / 2) {

        }
        else if (idx == (int) NX / 2 && NX % 2 != 0) {

        }
        else {

        }
    }
}

double CSSWM::interpolate(double A1, double A2, double V1, double V2, double B) {
    if (!((A1 < B && B < A2) || (A1 > B && B > A2))) {
        std::cout << "Error at interpolation4lat: " << A1 << " " << B << " " << A2 << std::endl;
    }
    return V1 + (V2-V1) * (B-A1) / (A2-A1);
}