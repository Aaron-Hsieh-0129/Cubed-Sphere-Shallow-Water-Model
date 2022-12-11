#include "outputs.hpp"

CSSWM model;
int main(void) {
    clock_t start, stop;
    start = clock();

    Init::Init2d(model);
    Outputs::output_parameter(model);
    Outputs::output_h(0, model);
    Outputs::output_u(0, model);
    Outputs::output_v(0, model);

    stop = clock();
    std::cout << double(stop - start) / CLOCKS_PER_SEC << " s" << std::endl;
    return 0;
}