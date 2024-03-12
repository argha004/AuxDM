//
// Created by Arghadwip Paul.
//

/*
 * Reads QuadFile and stores x,y,z in a flattened vector
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include "../include/readQuadWeightsFile.h"

std::pair<std::vector<double>, std::vector<double>>
readQuadWeightsFile(const std::string& quadweightsFile) {

    std::vector<double> quadData; // Stores the first three values of each line
    std::vector<double> quadWt;   // Stores the fifth value of each line
    std::ifstream file(quadweightsFile);
    std::string line;

    // Check if the file is open
    if (file.is_open()) {
        double number;
        // Read each line from the file
        while (getline(file, line)) {
            std::istringstream iss(line);
            // Read the first three numbers from the line and store them in quadData
            for (int i = 0; i < 3; ++i) {
                if (iss >> number) {
                    quadData.push_back(number);
                }
            }
            // Skip the fourth number
            iss >> number; // Read and discard the fourth number

            // Read the fifth number and store it in quadWt
            if (iss >> number) {
                quadWt.push_back(number);
            }
        }
        file.close(); // Close the file
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }

    // Return both vectors
    return std::make_pair(quadData, quadWt);
}