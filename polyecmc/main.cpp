/**
 * ECMC - Event-Chain Monte Carlo for Hard Disks
 * Main program for running ECMC simulations
 */

#include "system.h"
#include "ecmc.h"
#include "para.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <map>

// Define M_PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Parameter structure
struct Parameters {
    std::string inputFile = "";
    std::string outputFile = "config_final.dat";
    std::string paramFile = "";
    long long n_chains = 100000;
    double chain_length = 1.0;
    long long sample_interval_pressure = 10000;
    long long sample_interval_snapshot = 0;  // 0 means no snapshot output
    long long n_equilibration_pressure = 0;
    unsigned int rng_seed = 0;  // 0 means random
    bool binary_output = false;
};

// Read parameters from parameter file
bool readParamFile(const std::string& filename, Parameters& params) {
    std::ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "Error: Cannot open parameter file: " << filename << "\n";
        return false;
    }
    
    std::string line;
    int lineNum = 0;
    while (std::getline(ifs, line)) {
        ++lineNum;
        // Remove comments
        size_t commentPos = line.find('#');
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }
        
        // Trim whitespace
        size_t start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) continue;  // Empty line
        size_t end = line.find_last_not_of(" \t\r\n");
        line = line.substr(start, end - start + 1);
        
        // Parse key = value
        size_t equalPos = line.find('=');
        if (equalPos == std::string::npos) continue;
        
        std::string key = line.substr(0, equalPos);
        std::string value = line.substr(equalPos + 1);
        
        // Trim key and value
        key = key.substr(0, key.find_last_not_of(" \t") + 1);
        key = key.substr(key.find_first_not_of(" \t"));
        value = value.substr(0, value.find_last_not_of(" \t") + 1);
        value = value.substr(value.find_first_not_of(" \t"));
        
        // Set parameters
        if (key == "input_file") {
            params.inputFile = value;
        } else if (key == "output_file") {
            params.outputFile = value;
        } else if (key == "n_chains") {
            params.n_chains = std::atoll(value.c_str());
        } else if (key == "chain_length") {
            params.chain_length = std::atof(value.c_str());
        } else if (key == "sample_interval_pressure") {
            params.sample_interval_pressure = std::atoll(value.c_str());
        } else if (key == "sample_interval_snapshot") {
            params.sample_interval_snapshot = std::atoll(value.c_str());
        } else if (key == "n_equilibration_pressure") {
            params.n_equilibration_pressure = std::atoll(value.c_str());
        } else if (key == "rng_seed") {
            params.rng_seed = static_cast<unsigned int>(std::atol(value.c_str()));
        } else if (key == "binary_output") {
            params.binary_output = (value == "true" || value == "1" || value == "yes");
        } else {
            std::cerr << "Warning: Unknown parameter '" << key << "' at line " << lineNum << "\n";
        }
    }
    
    return true;
}

void printUsage(const char* progName) {
    std::cout << "Usage: " << progName << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -p <file>      Parameter file (optional)\n";
    std::cout << "  -i <file>      Input configuration file (required)\n";
    std::cout << "  -o <file>      Output configuration file (default: config_final.dat)\n";
    std::cout << "  -n <chains>    Number of event chains (default: 100000)\n";
    std::cout << "  -l <length>    Chain length (default: 1.0)\n";
    std::cout << "  -s <interval>  Pressure sampling interval (default: 10000)\n";
    std::cout << "  -S <interval>  Snapshot interval (default: 0, disabled)\n";
    std::cout << "  -e <chains>    Equilibration chains for pressure (default: 0)\n";
    std::cout << "  -r <seed>      Random seed (default: random)\n";
    std::cout << "  -b             Use binary output format\n";
    std::cout << "  -h             Show this help\n\n";
    std::cout << "Parameter File Format:\n";
    std::cout << "  input_file = initial.dat\n";
    std::cout << "  output_file = final.dat\n";
    std::cout << "  n_chains = 1000000\n";
    std::cout << "  chain_length = 1.0\n";
    std::cout << "  sample_interval_pressure = 50000\n";
    std::cout << "  sample_interval_snapshot = 100000\n";
    std::cout << "  n_equilibration_pressure = 100000\n";
    std::cout << "  rng_seed = 42\n";
    std::cout << "  binary_output = true\n\n";
    std::cout << "Note: Command line options override parameter file settings.\n\n";
    std::cout << "Example:\n";
    std::cout << "  " << progName << " -p params.txt -i initial.dat\n";
    std::cout << "  " << progName << " -i initial.dat -n 1000000 -s 50000 -e 100000\n";
}

int main(int argc, char* argv[]) {
    // Default parameters
    Parameters params;
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-p" && i + 1 < argc) {
            params.paramFile = argv[++i];
            // Immediately read parameter file
            if (!readParamFile(params.paramFile, params)) {
                return 1;
            }
        } else if (arg == "-i" && i + 1 < argc) {
            params.inputFile = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            params.outputFile = argv[++i];
        } else if (arg == "-n" && i + 1 < argc) {
            params.n_chains = std::atoll(argv[++i]);
        } else if (arg == "-l" && i + 1 < argc) {
            params.chain_length = std::atof(argv[++i]);
        } else if (arg == "-s" && i + 1 < argc) {
            params.sample_interval_pressure = std::atoll(argv[++i]);
        } else if (arg == "-S" && i + 1 < argc) {
            params.sample_interval_snapshot = std::atoll(argv[++i]);
        } else if (arg == "-e" && i + 1 < argc) {
            params.n_equilibration_pressure = std::atoll(argv[++i]);
        } else if (arg == "-r" && i + 1 < argc) {
            params.rng_seed = static_cast<unsigned int>(std::atol(argv[++i]));
        } else if (arg == "-b") {
            params.binary_output = true;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }
    
    // If no random seed specified, use random value
    if (params.rng_seed == 0) {
        params.rng_seed = std::random_device{}();
    }
    
    // Check required parameters
    if (params.inputFile.empty()) {
        std::cerr << "Error: Input file required (-i or input_file in parameter file)\n\n";
        printUsage(argv[0]);
        return 1;
    }

    // Print configuration
    std::cout << "========================================\n";
    std::cout << "ECMC - Event-Chain Monte Carlo\n";
    std::cout << "========================================\n\n";
    std::cout << "Configuration:\n";
    if (!params.paramFile.empty()) {
        std::cout << "  Parameter file: " << params.paramFile << "\n";
    }
    std::cout << "  Input:  " << params.inputFile << "\n";
    std::cout << "  Output: " << params.outputFile << (params.binary_output ? " (binary)" : " (text)") << "\n";
    std::cout << "  Chains: " << params.n_chains << "\n";
    std::cout << "  Chain length: " << params.chain_length << "\n";
    std::cout << "  Pressure sampling: every " << params.sample_interval_pressure << " chains\n";
    if (params.sample_interval_snapshot > 0) {
        std::cout << "  Snapshot interval: " << params.sample_interval_snapshot << " chains\n";
    } else {
        std::cout << "  Snapshot interval: disabled\n";
    }
    std::cout << "  Equilibration: " << params.n_equilibration_pressure << " chains\n";
    std::cout << "  Random seed: " << params.rng_seed << "\n\n";

    // Load initial configuration
    System sys;
    try {
        std::cout << "Loading configuration from " << params.inputFile << "...\n";
        sys.loadConfig(params.inputFile);
        std::cout << "  Loaded " << sys.size() << " particles\n";
        std::cout << "  Box: " << sys.boxsize[0] << " x " << sys.boxsize[1] << "\n";
        std::cout << "  Radius: mean = " << sys.radiusMean 
                  << ", max = " << sys.radiusMax 
                  << ", min = " << sys.radiusMin << "\n";
        std::cout << "  Packing fraction: φ = " << sys.getPackingFraction() << "\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error loading configuration: " << e.what() << "\n";
        return 1;
    }

    // Check initial overlap
    if (sys.checkOverlap()) {
        std::cerr << "Warning: Initial configuration has overlaps!\n\n";
    }

    // Set algorithm parameters
    AlgorithmPara algo_para;
    algo_para.set_n_steps(params.n_chains);
    algo_para.set_chain_length(params.chain_length);
    algo_para.set_sample_interval_pressure(params.sample_interval_pressure);
    algo_para.set_sample_interval_snapshot(params.sample_interval_snapshot);
    algo_para.set_rng_seed(params.rng_seed);
    algo_para.set_n_equilibration_pressure(params.n_equilibration_pressure);

    PhysicalPara phys_para;
    phys_para.set_n_atoms(sys.size());
    phys_para.set_box({sys.boxsize[0], sys.boxsize[1]});

    // Run ECMC
    std::cout << "========================================\n";
    std::cout << "Running ECMC Simulation\n";
    std::cout << "========================================\n\n";

    try {
        EcmcEngine engine(sys, algo_para, phys_para);
        engine.run();
        
        std::cout << "\n========================================\n";
        std::cout << "Simulation Completed\n";
        std::cout << "========================================\n\n";
        std::cout << "Total collisions: " << engine.getCollisionCount() << "\n";
        std::cout << "Collisions per chain: " 
                  << static_cast<double>(engine.getCollisionCount()) / params.n_chains << "\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error during simulation: " << e.what() << "\n";
        return 1;
    }

    // Check final overlap
    if (sys.checkOverlap()) {
        std::cerr << "Error: Final configuration has overlaps!\n";
        return 1;
    }

    // Save final configuration
    try {
        std::cout << "Saving final configuration to " << params.outputFile << "...\n";
        if (params.binary_output) {
            sys.saveConfigBinary(params.outputFile);
        } else {
            sys.saveConfig(params.outputFile);
        }
        std::cout << "Done!\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error saving configuration: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
