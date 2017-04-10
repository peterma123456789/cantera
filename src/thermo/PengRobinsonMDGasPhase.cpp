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
    int i = 1;
    while (infile >> T_tmp >> rho_tmp >> cp_tmp) {
        if (i == 1) {
            T_SPLINE.push_back(10.0);
            rho_SPLINE.push_back(rho_tmp);
            cp_SPLINE.push_back(cp_tmp * 4184e3);
            T_SPLINE.push_back(T_tmp - 10.0);
            rho_SPLINE.push_back(rho_tmp);
            cp_SPLINE.push_back(cp_tmp * 4184e3);
        }
        T_SPLINE.push_back(T_tmp);
        rho_SPLINE.push_back(rho_tmp);
        cp_SPLINE.push_back(cp_tmp * 4184e3);

        std::cout << i << ": " <<  T_tmp << " " << rho_tmp << " " << cp_tmp << std::endl;
        i++;
    }
    T_SPLINE.push_back(T_tmp + 10.0);
    rho_SPLINE.push_back(rho_tmp);
    cp_SPLINE.push_back(cp_tmp * 4184e3);
    T_SPLINE.push_back(5000.0);
    rho_SPLINE.push_back(rho_tmp);
    cp_SPLINE.push_back(cp_tmp * 4184e3);

    infile.close();

    mySplineRho.set_points(T_SPLINE, rho_SPLINE);
    mySplineCp.set_points(T_SPLINE, cp_SPLINE);
}

void PengRobinsonMDGasPhase::setLinear(std::string fileName)
{
    myLinearT.clear();
    myLinearRho.clear();
    myLinearCp.clear();

    std::ifstream infile(fileName);
    if (!infile) {
        std::cout << "couldn't open " << fileName << std::endl;
        throw CanteraError("Error in PengRobinsonMDGasPhase.h",
                           "read file for spline error...\n\n");
    }
    // infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    double T_tmp, rho_tmp, cp_tmp;
    int i = 1;
    while (infile >> T_tmp >> rho_tmp >> cp_tmp) {

        myLinearT.push_back(T_tmp);
        myLinearRho.push_back(rho_tmp);
        myLinearCp.push_back(cp_tmp * 4184e3);

        std::cout << i << ": " <<  T_tmp << " " << rho_tmp << " " << cp_tmp << std::endl;
        i++;
    }

    infile.close();
}

double PengRobinsonMDGasPhase::linearInterp(const std::vector<double>& T,
                                            const std::vector<double>& var,
                                            double Tnow) const
{
    int N = T.size();

    // no extrapolation
    if (Tnow <= T[0]) return var[0];
    if (Tnow >= T[N - 1]) return var[N - 1];

    int i = 1;
    while ((Tnow > T[i]) && (i < N - 1)) i++;
    return var[i - 1] +
           (var[i] - var[i - 1]) / (T[i] - T[i - 1]) * (Tnow - T[i - 1]);
}

doublereal PengRobinsonMDGasPhase::pressure() const { return 65.0 * OneAtm; }

void PengRobinsonMDGasPhase::setPressure(doublereal p)
{
    // temperature is set already
    PengRobinsonGasPhase::setPressure(p); // need this for partial h
    //printf("%f: %f, %f\n", temperature(), density(), mySplineRho(temperature()));
    //if (mySplineRho(temperature()) < 0) printf("%f\n", temperature());

    //if (temperature() < 80) {
    //    setDensity(1551.0);
    //} else if (temperature() > 228) {
    //    setDensity(154.0);
    //} else {
    //    setDensity(mySplineRho(temperature()));
    //}

    //setDensity(10.0);
    //setDensity(mySplineRho(temperature()));
    setDensity(linearInterp(myLinearT, myLinearRho, temperature()));
}

doublereal PengRobinsonMDGasPhase::cp_mole() const
{
    // temperature is set already
    //return mySplineCp(temperature());
    //printf("%g: %g, %g\n", temperature(), PengRobinsonGasPhase::cp_mole(), mySplineCp(temperature()));

    //return PengRobinsonGasPhase::cp_mole();
    //return 4.0e4;
    return linearInterp(myLinearT, myLinearCp, temperature());
}

}
