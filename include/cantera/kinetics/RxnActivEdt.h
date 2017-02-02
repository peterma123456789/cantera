/**
 * @file RxnActivEdt.h
 * @ingroup chemkinetics
 */

// Copyright 2016  Hao Wu (wuhao@stanford.edu)
// #define CT_FAST_CHEM

#ifndef CT_RXNACTIVDT_H
#define CT_RXNACTIVDT_H

#include "cantera/kinetics/FalloffMgr.h"
#include "cantera/kinetics/RateCoeffMgr.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/ThirdBodyCalc.h"
#include "cantera/thermo/ThermoPhase.h"

#include "cantera/numerics/eigen_dense.h"
#include "cantera/numerics/eigen_sparse.h"

// #include "GasKinetics.h"

namespace Cantera {
class RxnActivEdt {
public:
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpCMat;
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpRMat;
  typedef Eigen::VectorXd VecType;
  typedef std::vector<shared_ptr<Reaction>> vector_rxnptr;

  template <typename T>
  static void editVecs(std::vector<T> &_leftVec,
                       const std::vector<T> &_rightVec,
                       const std::vector<size_t> &_idList) {
    _leftVec.resize(_idList.size());
    for (size_t i = 0; i < _leftVec.size(); ++i) {
      _leftVec[i] = _rightVec[_idList[i]];
    }
  }

  template <typename T>
  static void editStoichMng(StoichManagerN &_leftStMng,
                            const StoichManagerN &_rightStMng,
                            const std::vector<T> &_iActiv,
                            const std::vector<size_t> &_idMap) {
    _leftStMng.m_c1_list.clear();
    for (const auto &c : _rightStMng.m_c1_list) {
      if (!_iActiv[c.m_rxn]) continue;
      _leftStMng.m_c1_list.emplace_back(_idMap[c.m_rxn], c.m_ic0);
    }
    _leftStMng.m_c2_list.clear();
    for (const auto &c : _rightStMng.m_c2_list) {
      if (!_iActiv[c.m_rxn]) continue;
      _leftStMng.m_c2_list.emplace_back(_idMap[c.m_rxn], c.m_ic0, c.m_ic1);
    }
    _leftStMng.m_c3_list.clear();
    for (const auto &c : _rightStMng.m_c3_list) {
      if (!_iActiv[c.m_rxn]) continue;
      _leftStMng.m_c3_list.emplace_back(_idMap[c.m_rxn], c.m_ic0, c.m_ic1,
                                        c.m_ic2);
    }
    _leftStMng.m_cn_list.clear();
    for (const auto &c : _rightStMng.m_cn_list) {
      if (!_iActiv[c.m_rxn]) continue;
      _leftStMng.m_cn_list.emplace_back(_idMap[c.m_rxn], c.m_ic, c.m_order,
                                        c.m_stoich);
    }
  }

  template <typename T>
  static void
  editRevs(std::vector<size_t> &_leftRev, std::vector<size_t> &_leftIrrev,
           const std::vector<size_t> &_rightRev,
           const std::vector<size_t> &_rightIrrev,
           const std::vector<T> &_iActiv, const std::vector<size_t> &_idMap) {
    _leftRev.clear();
    _leftIrrev.clear();
    for (size_t i = 0; i < _rightRev.size(); ++i) {
      if (!_iActiv[_rightRev[i]]) continue;
      _leftRev.push_back(_idMap[_rightRev[i]]);
    }
    for (size_t i = 0; i < _rightIrrev.size(); ++i) {
      if (!_iActiv[_rightIrrev[i]]) continue;
      _leftIrrev.push_back(_idMap[_rightIrrev[i]]);
    }
  }

  template <typename R, typename T>
  static void editRates(Rate1<R> &_leftRates, const Rate1<R> &_rightRates,
                        const std::vector<T> &_iActiv,
                        const std::vector<size_t> &_idMap) {
    _leftRates.m_rxn.clear();
    _leftRates.m_rates.clear();
    for (size_t i = 0; i < _rightRates.m_rates.size(); ++i) {
      if (!_iActiv[_rightRates.m_rxn[i]]) continue;
      _leftRates.m_rxn.push_back(_idMap[_rightRates.m_rxn[i]]);
      _leftRates.m_rates.push_back(_rightRates.m_rates[i]);
    }
  }

  template <typename T>
  static void editThirdBody(ThirdBodyCalc &_leftTBC,
                            const ThirdBodyCalc &_rightTBC,
                            const std::vector<T> &_iActiv,
                            const std::vector<size_t> &_idMap) {
    _leftTBC.m_reaction_index.clear();
    _leftTBC.m_default.clear();
    _leftTBC.m_species.clear();
    _leftTBC.m_eff.clear();
    for (size_t i = 0; i < _rightTBC.m_reaction_index.size(); ++i) {
      if (!_iActiv[_rightTBC.m_reaction_index[i]]) continue;
      _leftTBC.m_reaction_index.push_back(
          _idMap[_rightTBC.m_reaction_index[i]]);
      _leftTBC.m_default.push_back(_rightTBC.m_default[i]);
      _leftTBC.m_species.push_back((_rightTBC.m_species[i]));
      _leftTBC.m_eff.push_back((_rightTBC.m_eff[i]));
    }
  }

  static void editFalloff(FalloffMgr &_leftMgr, const FalloffMgr &_rightMgr,
                          const std::vector<size_t> &_idList) {
    _leftMgr.m_rxn.clear();
    _leftMgr.m_falloff.clear();
    _leftMgr.m_offset.clear();
    _leftMgr.m_worksize = 0;
#ifndef CT_FAST_CHEM
    _leftMgr.m_indices.clear();
#endif
    for (size_t i = 0; i < _idList.size(); ++i) {
      _leftMgr.m_rxn.push_back(i);
      _leftMgr.m_offset.push_back(_leftMgr.m_worksize);
      _leftMgr.m_worksize += _rightMgr.m_falloff[_idList[i]]->workSize();
      _leftMgr.m_falloff.push_back(_rightMgr.m_falloff[_idList[i]]);
      _leftMgr.m_reactionType.push_back(_rightMgr.m_reactionType[_idList[i]]);
#ifndef CT_FAST_CHEM
      _leftMgr.m_indices.insert(
          _leftMgr.m_indices.end(),
          std::pair<int, int>(i, _leftMgr.m_falloff.size() - 1));
#endif
    }
  }

protected:
};
}

#endif
