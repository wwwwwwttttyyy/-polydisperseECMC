#include "InitialConfig/2d2packing/initconfig.h"
#include <iostream>
#include <iomanip>
#include <sstream>

int main() {
    std::cout << "Testing ECMC configuration generation with volume fraction control..." << std::endl;
    for (int i=0;i<21;i++){
        double target_vf = 0.9-i*0.01;
        auto [data2, sigma2, box2] = ini_config(1, 13, 16, target_vf);
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3) << target_vf;
        std::string vf_str = oss.str();
        
        write_ecmc_config(data2, box2, "config1_vf_" + vf_str + ".dat");
    }
    
    return 0;
}
