#include <fstream>
#include <cstdlib>
#include <cstring>
#include "constrcution.hpp"

class Outputs {
public:
    static void output_parameter(CSSWM &);

private:
    static void create_directory(std::string);
};