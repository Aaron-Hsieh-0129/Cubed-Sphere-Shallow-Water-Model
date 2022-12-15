#include "outputs.hpp"

using std::fstream;
using std::ios;
using std::string;

void Outputs::create_directory(string directory_name) {
    string str = "mkdir -p " + directory_name;
    const char *command = str.c_str();
    const int dir_err = system(command);
    if (-1 == dir_err) {
        std::cout << "Error on creating directory!\n" << std::endl;
        return;
    }
    return;
}

void Outputs::output_parameter(CSSWM &model) {
    create_directory("../outputs/grids");

    fstream fout[20];
    string dir = "../outputs/grids/";
    string grid[20] = {"lon.txt", "lat.txt", "lon_L.txt", "lat_L.txt", "lon_R.txt", "lat_R.txt", "lon_U.txt", "lat_U.txt", "lon_D.txt", "lat_D.txt", 
                       "x.txt", "y.txt", "x_L.txt", "y_L.txt", "x_R.txt", "y_R.txt", "x_U.txt", "y_U.txt", "x_D.txt", "y_D.txt"};

    for (int i = 0; i < 20; i++) {
        fout[i].open(dir + grid[i], ios::out);
    }

    for (int p = 0; p < 6; p++) {
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                fout[0] << model.csswm[p].lon[i][j] << " ";
                fout[1] << model.csswm[p].lat[i][j] << " ";
                fout[2] << model.csswm[p].lon_L[i][j] << " ";
                fout[3] << model.csswm[p].lat_L[i][j] << " ";
                fout[4] << model.csswm[p].lon_R[i][j] << " ";
                fout[5] << model.csswm[p].lat_R[i][j] << " ";
                fout[6] << model.csswm[p].lon_U[i][j] << " ";
                fout[7] << model.csswm[p].lat_U[i][j] << " ";
                fout[8] << model.csswm[p].lon_D[i][j] << " ";
                fout[9] << model.csswm[p].lat_D[i][j] << " ";
                
                fout[10] << model.csswm[p].x[i][j] << " ";
                fout[11] << model.csswm[p].y[i][j] << " ";
                fout[12] << model.csswm[p].x_L[i][j] << " ";
                fout[13] << model.csswm[p].y_L[i][j] << " ";
                fout[14] << model.csswm[p].x_R[i][j] << " ";
                fout[15] << model.csswm[p].y_R[i][j] << " ";
                fout[16] << model.csswm[p].x_U[i][j] << " ";
                fout[17] << model.csswm[p].y_U[i][j] << " ";
                fout[18] << model.csswm[p].x_D[i][j] << " ";
                fout[19] << model.csswm[p].y_D[i][j] << " ";
            }
        }
    }
}

void Outputs::output_h(int n, CSSWM &model) {
    create_directory("../outputs/h");

    fstream fouth;
    string hname = "../outputs/h/h_" + std::to_string(n) + ".txt";
    fouth.open(hname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                fouth << model.csswm[p].h[i][j] << " ";
            }
        }
    }
    return;
}

void Outputs::output_u(int n, CSSWM &model) {
    create_directory("../outputs/u");

    fstream foutu;
    string uname = "../outputs/u/u_" + std::to_string(n) + ".txt";
    foutu.open(uname, std::ios::out);

    fstream foutu_lon_lat;
    string u_lon_latname = "../outputs/u_lon_lat/u_lon_lat_" + std::to_string(n) + ".txt";
    foutu_lon_lat.open(u_lon_latname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                foutu << model.csswm[p].u[i][j] << " ";
                // foutu_lon_lat << model.Cube2Sphere_U(model, p, i, j) << " ";
            }
        }
    }
    return;
}

void Outputs::output_v(int n, CSSWM &model) {
    create_directory("../outputs/v");

    fstream foutv;
    string vname = "../outputs/v/v_" + std::to_string(n) + ".txt";
    foutv.open(vname, std::ios::out);

    fstream foutv_lon_lat;
    string v_lon_latname = "../outputs/v_lon_lat/v_lon_lat_" + std::to_string(n) + ".txt";
    foutv_lon_lat.open(v_lon_latname, std::ios::out);
    for (int p = 0; p < 6; p++) {
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                foutv << model.csswm[p].v[i][j] << " ";
                // foutv_lon_lat << model.Cube2Sphere_V(model, p, i, j) << " ";
            }
        }
    }
    return;
}