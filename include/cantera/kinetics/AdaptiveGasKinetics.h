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

  virtual void updateROP();

  //! Update temperature-dependent portions of reaction rates and falloff
  //! functions.
  virtual void update_rates_T();

  //! Update properties that depend on concentrations.
  //! Currently the enhanced collision partner concentrations are updated
  //! here, as well as the pressure-dependent portion of P-log and Chebyshev
  //! reactions.
  virtual void update_rates_C();

  //! @}
  //! @name Adaptive Reactions Routines
  //! @{
  virtual void setMuteReaction(const size_t i, const bool flag = true);

  // NOTE: Shall not need the following two.
  // Shall be able to perform full-kinetics calculations by casting to
  // the base class (GasKinetics).
  //
  // virtual void updateROP_full();
  //
  // virtual void getNetProductionRates_full(doublereal* wdot);

  virtual const std::vector<bool>& iMuteRractions()
    { return m_imuted_reactions; };
  //@}

protected:
  //! Vector of flag for reaction muting
  std::vector<bool> m_imuted_reactions;
  //! Vector of flag for falloff reaction muting
  std::vector<bool> m_imuted_fall;

  virtual void processFalloffReactions();

  //! Update the equilibrium constants in molar units.
  virtual void updateKc();

};
}

#endif
