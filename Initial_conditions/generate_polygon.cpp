#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream> 
#include <ctime>
#include <vector>
#include <functional>
#include <string>
#include <unordered_map>
#include <map>
#include <cassert>
#include <filesystem>
#include <algorithm>

int main() {
    int N = 25;
    double R = 1.0;
    double a = 4.0;
    std::ofstream outfile("Polygon_" + std::to_string(N) + ".txt");

    for (int j = 0; j < N; j++){
        outfile << 1.0 << ' ';
        outfile << R*cos(2*M_PI*j/N) << ' ' << R*sin(2*M_PI*j/N) << ' ' << 0.0 << ' ';
        outfile << -a*sin(2*M_PI*j/N) << ' ' << a*cos(2*M_PI*j/N) << ' ' << 0.0 << ' ' << std::endl;
    }

    std::cout << "Done" << std::endl;
    return 0;
}
