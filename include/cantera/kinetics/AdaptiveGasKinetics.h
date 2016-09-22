/**
 * @file AdaptiveGasKinetics.h
 * @ingroup chemkinetics
 */

// Copyright 2016  Hao Wu (wuhao@stanford.edu)

#ifndef CT_ADAPTGASKINETICS_H
#define CT_ADAPTGASKINETICS_H

#include "GasKinetics.h"
#include "ThirdBodyCalc.h"
#include "FalloffMgr.h"
#include "Reaction.h"
#include "RxnActiveMgr.h"

namespace Cantera
{
class AdaptiveGasKinetics : public GasKinetics
{
public:
  AdaptiveGasKinetics(thermo_t* thermo = 0);

  virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const;

  virtual int type() const {
      warn_deprecated("AdaptiveGas::type",
                      "To be removed after Cantera 2.3.");
      return 101;
  }

  virtual std::string kineticsType() const {
      return "AdaptiveGas";
  }

  //! @}
  //! @name Reaction Mechanism Setup Routines
  //! @{
  virtual bool addReaction(shared_ptr<Reaction> r);
  //@}

  /**
   * Given an array of species properties 'g', return in array 'dg' the
   * change in this quantity in the reversible reactions. Array 'g' must
   * have a length at least as great as the number of species, and array
   * 'dg' must have a length as great as the total number of reactions.
   * This method only computes 'dg' for the reversible reactions, and the
   * entries of 'dg' for the irreversible reactions are unaltered. This is
   * primarily designed for use in calculating reverse rate coefficients
   * from thermochemistry for reversible reactions.
   */
  virtual void getRevReactionDelta_adp(const std::vector<std::uint8_t>& iactive,
                                       const doublereal* g, doublereal* dg);

  /**
   * Species net production rates [kmol/m^3/s or kmol/m^2/s]. Return the
   * species net production rates (creation - destruction) in array wdot,
   * which must be dimensioned at least as large as the total number of
   * species. @see nTotalSpecies. Adpative kinetics is activated.
   *
   * @param wdot   Output vector of net production rates. Length: m_kk.
   */
  virtual void getNetProductionRates(doublereal* wdot);

  /**
   * Species net production rates [kmol/m^3/s or kmol/m^2/s]. Return the
   * species net production rates (creation - destruction) in array wdot,
   * which must be dimensioned at least as large as the total number of
   * species. @see nTotalSpecies. Adpative kinetics is deactivated.
   *
   * @param wdot   Output vector of net production rates. Length: m_kk.
   */
  virtual void getNetProductionRatesFull(doublereal* wdot);

  void updateROP_adp();

  //! Update temperature-dependent portions of reaction rates and falloff
  //! functions.
  void update_rates_T_adp();

  //! Update properties that depend on concentrations.
  //! Currently the enhanced collision partner concentrations are updated
  //! here, as well as the pressure-dependent portion of P-log and Chebyshev
  //! reactions.
  void update_rates_C_adp();

  //! @}
  //! @name Adaptive Reactions Routines
  //! @{
  //! Prepare adaptation manager
  void prepareAdaptation();

  //! Update adaptation manager
  void updateAdaptation(const double relTol = 1e-6, const double absTol = 1e-8);

  virtual void setActiveReaction(const size_t i, const bool flag = true);

  virtual const std::vector<std::uint8_t>& iActiveRractions()
    { return m_active_reactions; };
  //@}

protected:
  //! Note: std::vector<std::uint8_t> is prefered over std::vector<bool> based
  //! on test performed with clang++ on Mac OSX with Intel i7 cpu.
  //! One may achieve better results with std::vector<bool> with other
  //! computational environment.
  //! 
  //! Vector of flag for reaction muting
  std::vector<std::uint8_t> m_active_reactions;
  //! Vector of flag for falloff reaction muting
  std::vector<std::uint8_t> m_active_fall;

  // Reaction activation manager
  RxnActiveMgr m_rxnactmgr;
  bool m_ready_rxnactmgr;

  void processFalloffReactions_adp();

  //! Update the equilibrium constants in molar units.
  void updateKc_adp();

};
}

#endif
