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
  VectorXd rmse(5);
	rmse << 0,0,0,0,0;

	if (estimations.size() == 0) {
	    cout << "Empty vector." << endl;
	    return rmse;
	}
	if (estimations.size() != ground_truth.size()) {
	    cout << "Size is different."  << endl;
	    return rmse;
	}

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i) {
    VectorXd res = estimations[i] - ground_truth[i];
    rmse(0) += pow(res(0),2);
    rmse(1) += pow(res(1),2);
    rmse(2) += pow(res(2),2);
    rmse(3) += pow(res(3),2);
    rmse(4) += pow(res(4),2);
	}

	//calculate the mean
	// ... your code here
	rmse = rmse / estimations.size();

	//calculate the squared root
	// ... your code here
	rmse(0) = sqrt(rmse(0));
	rmse(1) = sqrt(rmse(1));
	rmse(2) = sqrt(rmse(2));
  rmse(3) = sqrt(rmse(3));
  rmse(4) = sqrt(rmse(4));

  cout << "RMSE: " << rmse << endl;
	//return the result
	return rmse;
}