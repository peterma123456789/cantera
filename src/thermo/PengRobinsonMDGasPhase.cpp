/**
 *  @file PengRobinsonGasPhase.cpp
 *   ThermoPhase object for the Peng-Robinson gas equation of
 * state - workhorse for %Cantera (see \ref thermoprops
 * and class \link Cantera::PengRobinsonGasPhase PengRobinsonGasPhase\endlink).
 */

#include "cantera/thermo/PengRobinsonMDGasPhase.h"
#include "cantera/base/vec_functions.h"

#include "cantera/numerics/spline.h"

using namespace std;

namespace Cantera
{

void PengRobinsonMDGasPhase::setSpline(std::string fileName)
{
    std::vector<double> T_SPLINE, rho_SPLINE, cp_SPLINE;

    std::ifstream infile(fileName);
    if (!infile) {
        std::cout << "couldn't open " << fileName << std::endl;
        throw CanteraError("Error in PengRobinsonMDGasPhase.h",
                           "read file for spline error...\n\n");
    }
    // infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    double T_tmp, rho_tmp, cp_tmp;
    while (infile >> T_tmp >> rho_tmp >> cp_tmp) {
        T_SPLINE.push_back(T_tmp);
        rho_SPLINE.push_back(rho_tmp);
        cp_SPLINE.push_back(cp_tmp);
        std::cout << T_tmp << rho_tmp << cp_tmp << std::endl;
    }
    infile.close();

    mySplineRho.set_points(T_SPLINE, rho_SPLINE);
    mySplineCp.set_points(T_SPLINE, cp_SPLINE);
}

}
