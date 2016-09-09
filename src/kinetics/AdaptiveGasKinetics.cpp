/**
 *  @file AdaptiveGasKinetics.cpp Adaptive homogeneous kinetics in ideal gases
 */

// Copyright 2016 Hao Wu (wuhao@stanford.edu)

#include "cantera/kinetics/AdaptiveGasKinetics.h"

using namespace std;

namespace Cantera {

AdaptiveGasKinetics::AdaptiveGasKinetics(thermo_t *thermo)
    : GasKinetics(thermo) {}

Kinetics *AdaptiveGasKinetics::duplMyselfAsKinetics(
    const std::vector<thermo_t *> &tpVector) const {
  GasKinetics *gK = new AdaptiveGasKinetics(*this);
  gK->assignShallowPointers(tpVector);
  return gK;
}

void AdaptiveGasKinetics::update_rates_T()
{
    doublereal T = thermo().temperature();
    doublereal P = thermo().pressure();
    m_logStandConc = log(thermo().standardConcentration());
    doublereal logT = log(T);

    if (T != m_temp) {
        if (!m_rfn.empty()) {
            m_rates.update(T, logT, m_imuted_reactions, m_rfn.data());
        }

        if (!m_rfn_low.empty()) {
            m_falloff_low_rates.update(T, logT, m_imuted_fall,
                                       m_rfn_low.data());
            m_falloff_high_rates.update(T, logT, m_imuted_fall,
                                        m_rfn_high.data());
        }
        if (!falloff_work.empty()) {
            m_falloffn.updateTemp(T, m_imuted_fall, falloff_work.data());
        }
        updateKc();
        m_ROP_ok = false;
    }

    if (T != m_temp || P != m_pres) {
        if (m_plog_rates.nReactions()) {
            m_plog_rates.update(T, logT, m_imuted_reactions, m_rfn.data());
            m_ROP_ok = false;
        }

        if (m_cheb_rates.nReactions()) {
            m_cheb_rates.update(T, logT, m_imuted_reactions, m_rfn.data());
            m_ROP_ok = false;
        }
    }
    m_pres = P;
    m_temp = T;
}

void AdaptiveGasKinetics::update_rates_C()
{
    thermo().getActivityConcentrations(m_conc.data());
    doublereal ctot = thermo().molarDensity();

    // 3-body reactions
    if (!concm_3b_values.empty()) {
        m_3b_concm.update(m_conc, ctot, concm_3b_values.data());
    }

    // Falloff reactions
    if (!concm_falloff_values.empty()) {
        m_falloff_concm.update(m_conc, ctot, m_imuted_reactions,
                               concm_falloff_values.data());
    }

    // P-log reactions
    if (m_plog_rates.nReactions()) {
        double logP = log(thermo().pressure());
        m_plog_rates.update_C(m_imuted_reactions, &logP);
    }

    // Chebyshev reactions
    if (m_cheb_rates.nReactions()) {
        double log10P = log10(thermo().pressure());
        m_cheb_rates.update_C(m_imuted_reactions, &log10P);
    }

    m_ROP_ok = false;
}

void AdaptiveGasKinetics::updateKc()
{
    thermo().getStandardChemPotentials(m_grt.data());
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);

    // compute Delta G^0 for all reversible reactions
    getRevReactionDelta(m_grt.data(), m_rkcn.data());

    doublereal rrt = 1.0 / thermo().RT();
    for (size_t i = 0; i < m_revindex.size(); i++) {
        size_t irxn = m_revindex[i];
        if (m_imuted_reactions[irxn]) {
          m_rkcn[irxn] = 0.0;
        } else {
          m_rkcn[irxn] = std::min(exp(m_rkcn[irxn]*rrt - m_dn[irxn]*m_logStandConc),
                                  BigNumber);
        }
    }

    for (size_t i = 0; i != m_irrev.size(); ++i) {
        m_rkcn[ m_irrev[i] ] = 0.0;
    }
}

void AdaptiveGasKinetics::updateROP()
{
    update_rates_C();
    update_rates_T();
    if (m_ROP_ok) {
        return;
    }

    // copy rate coefficients into ropf
    m_ropf = m_rfn;

    // multiply ropf by enhanced 3b conc for all 3b rxns
    if (!concm_3b_values.empty()) {
        m_3b_concm.multiply(m_ropf.data(), concm_3b_values.data());
    }

    if (m_falloff_high_rates.nReactions()) {
        processFalloffReactions();
    }

    // multiply by perturbation factor
    multiply_each(m_ropf.begin(), m_ropf.end(), m_perturb.begin());

    // copy the forward rates to the reverse rates
    m_ropr = m_ropf;

    // for reverse rates computed from thermochemistry, multiply the forward
    // rates copied into m_ropr by the reciprocals of the equilibrium constants
    multiply_each(m_ropr.begin(), m_ropr.end(), m_rkcn.begin());

    // multiply ropf by concentration products
    m_reactantStoich.multiply(m_conc.data(), m_ropf.data());

    // for reversible reactions, multiply ropr by concentration products
    m_revProductStoich.multiply(m_conc.data(), m_ropr.data());

    for (size_t j = 0; j != nReactions(); ++j) {
        m_ropnet[j] = (m_imuted_reactions[j]) ? 0.0 : m_ropf[j] - m_ropr[j];
    }

    for (size_t i = 0; i < m_rfn.size(); i++) {
        AssertFinite(m_rfn[i], "GasKinetics::updateROP",
                     "m_rfn[{}] is not finite.", i);
        AssertFinite(m_ropf[i], "GasKinetics::updateROP",
                     "m_ropf[{}] is not finite.", i);
        AssertFinite(m_ropr[i], "GasKinetics::updateROP",
                     "m_ropr[{}] is not finite.", i);
    }
    m_ROP_ok = true;
}

void AdaptiveGasKinetics::processFalloffReactions()
{
    // use m_ropr for temporary storage of reduced pressure
    vector_fp& pr = m_ropr;

    for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
        if (m_imuted_fall[i]) {
          pr[i] = 0.;
        } else {
          pr[i] = concm_falloff_values[i] * m_rfn_low[i] / (m_rfn_high[i] + SmallNumber);
        }
        AssertFinite(pr[i], "GasKinetics::processFalloffReactions",
                     "pr[{}] is not finite.", i);
    }

    m_falloffn.pr_to_falloff(pr.data(), m_imuted_fall, falloff_work.data());

    for (size_t i = 0; i < m_falloff_low_rates.nReactions(); i++) {
        if (m_imuted_fall[i]) {
          pr[i] = 0.;
        } else {
          if (reactionType(m_fallindx[i]) == FALLOFF_RXN) {
              pr[i] *= m_rfn_high[i];
          } else { // CHEMACT_RXN
              pr[i] *= m_rfn_low[i];
          }
        }
    }

    scatter_copy(pr.begin(), pr.begin() + m_falloff_low_rates.nReactions(),
                 m_ropf.begin(), m_fallindx.begin());
}

bool AdaptiveGasKinetics::addReaction(shared_ptr<Reaction> r)
{
  // operations common to all reaction types
  bool added = GasKinetics::addReaction(r);
  if (!added) {
      return false;
  }
  //  grow m_imuted_reactions
  m_imuted_reactions.push_back(false);
  if (m_imuted_reactions.size() != nReactions())
    throw CanteraError("AdaptiveGasKinetics::addReaction",
        "Size of the reaction muting flag does not match the number of reactions.");
  // special instruction for fallof reactions
  if (r->reaction_type == FALLOFF_RXN || r->reaction_type == CHEMACT_RXN) {
    m_imuted_fall.push_back(false);
    if (m_imuted_fall.size() != m_falloff_high_rates.nReactions())
      throw CanteraError("AdaptiveGasKinetics::addReaction",
          "Size of the falloff reaction muting flag does not match the number of falloff reactions.");
  }
  return true;
}

void AdaptiveGasKinetics::setMuteReaction(const size_t i, const bool flag)
{
  // Check if i is in range of reaction indices
  GasKinetics::checkReactionIndex(i);
  // Return if muting flag of reaction i does not change
  if (m_imuted_reactions[i] == flag) return;
  // Get reaction i
  shared_ptr<Reaction>& r = m_reactions[i];
  // Mark reaction i as muted
  m_imuted_reactions[i] = flag;
  if (r->reaction_type == FALLOFF_RXN || r->reaction_type == CHEMACT_RXN) {
    m_imuted_fall[m_rfallindx[i]] = flag;
  }
}

}
