/**
 *  @file PengRobinsonMDGasPhase.h
 *   ThermoPhase object for the ideal gas equation of
 * state - workhorse for %Cantera (see \ref thermoprops
 * and class \link Cantera::PengRobinsonMDGasPhase PengRobinsonMDGasPhase\endlink).
 */

//  Copyright 2001 California Institute of Technology
#ifndef CT_PENGROBINSONMDGASPHASE_H
#define CT_PENGROBINSONMDGASPHASE_H

#include "mix_defs.h"
#include "PengRobinsonGasPhase.h"

#include "cantera/numerics/spline.h"

#include <fstream>
#include <iostream>

namespace Cantera
{

//!  Derived from PengRobinsonGasPhase with density and cp functions overwritten
//!  by MD interpolations
class PengRobinsonMDGasPhase: public PengRobinsonGasPhase
{
public:
    //! Default empty Constructor
    PengRobinsonMDGasPhase() : PengRobinsonGasPhase()
    {
        std::vector<double> T_tmp, rho_tmp, cp_tmp;
        T_tmp.push_back(10.0);
        T_tmp.push_back(5000.0);
        rho_tmp.push_back(1000.0);
        rho_tmp.push_back(10.0);
        cp_tmp.push_back(4.0e4);
        cp_tmp.push_back(4.0e4);
        mySplineRho.set_points(T_tmp, rho_tmp);
        mySplineCp.set_points(T_tmp, cp_tmp);

        myLinearT.push_back(10.0);
        myLinearT.push_back(5000.0);
        myLinearRho.push_back(1000.0);
        myLinearRho.push_back(10.0);
        myLinearCp.push_back(4.0e4);
        myLinearCp.push_back(4.0e4);
    };

    //! Copy Constructor
    PengRobinsonMDGasPhase(const PengRobinsonMDGasPhase& right)
        : PengRobinsonGasPhase(right),
          mySplineRho(right.mySplineRho),
          mySplineCp(right.mySplineCp),
          myLinearT(right.myLinearT),
          myLinearRho(right.myLinearRho),
          myLinearCp(right.myLinearCp){};

    //! Assignment operator
    PengRobinsonMDGasPhase& operator=(const PengRobinsonMDGasPhase& right)
    {
        PengRobinsonGasPhase::operator=(right);
        mySplineRho = right.mySplineRho;
        mySplineCp = right.mySplineCp;
        myLinearT = right.myLinearT;
        myLinearRho = right.myLinearRho;
        myLinearCp = right.myLinearCp;
        return *this;
    };

    //! Equation of state flag.
    /*!
     *  Returns the value cPengRobinsonGas, defined in mix_defs.h.
     */
    virtual int eosType() const
    {
        return cPengRobinsonMDGas;
    }

    //! read in file and set up spline interpolation
    //! HARD-CODED for single species (AR) right now
    void setSpline(std::string fileName);

    //! read in file and set up linear interpolation
    //! HARD-CODED for single species (AR) right now
    void setLinear(std::string fileName);

    //! pressure()
    virtual doublereal pressure() const;

    //! Set the pressure at constant temperature and composition.
    virtual void setPressure(doublereal p);

    //! Get cp_mole
    virtual doublereal cp_mole() const;

private:
    // linear interpolation, used with setLinear()
    double linearInterp(const std::vector<double>& T,
                        const std::vector<double>& var, double Tnow) const;

    //! spline object
    tk::spline mySplineRho;
    tk::spline mySplineCp; // cp_mole

    std::vector<double> myLinearT;
    std::vector<double> myLinearRho;
    std::vector<double> myLinearCp;
};

}

#endif
