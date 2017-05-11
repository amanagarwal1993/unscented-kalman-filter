#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    // Always check inputs!
    if ((estimations.size() != ground_truth.size()) | (estimations.size() == 0)) {
        std::cout << "Invalid size of estimations data. Error!"<<"/n";
        return rmse;
    }
    
    
    for (unsigned int i = 0; i < estimations.size(); i++) {
        VectorXd error = estimations[i] - ground_truth[i];
        error = error.array() * error.array();
        
        rmse += error;
        //std::cout<<rmse<<"\n";
    };
    
    rmse = rmse / (estimations.size());
    
    return rmse;
}
