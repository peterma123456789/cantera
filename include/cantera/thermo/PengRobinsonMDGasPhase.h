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
    PengRobinsonMDGasPhase() : PengRobinsonGasPhase() {};

    //! Copy Constructor
    PengRobinsonMDGasPhase(const PengRobinsonMDGasPhase& right)
        : PengRobinsonGasPhase(right),
          mySplineRho(right.mySplineRho),
          mySplineCp(right.mySplineCp){};

    //! Assignment operator
    PengRobinsonMDGasPhase& operator=(const PengRobinsonMDGasPhase& right)
    {
        PengRobinsonGasPhase::operator=(right);
        mySplineRho = right.mySplineRho;
        mySplineCp = right.mySplineCp;
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

    //! read in file and set up spline
    //! HARD-CODED for single species (AR) right now
    void setSpline(std::string fileName);

private:
    //! spline object
    tk::spline mySplineRho;
    tk::spline mySplineCp;
};
}

#endif
