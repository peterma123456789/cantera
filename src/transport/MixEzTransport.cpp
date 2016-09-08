/**
 *  @file MixEzTransport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */
// copyright 2016 Hao Wu (wuhao@stanford.edu)

#include "cantera/transport/MixEzTransport.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{
MixEzTransport::MixEzTransport() :
    m_lambda(0.0),
    m_spcond_ok(false),
    m_condmix_ok(false),
    m_lewis_ok(false),
    m_debug(false)
{
}

MixEzTransport::MixEzTransport(const MixEzTransport& right) :
    GasTransport(right),
    m_lambda(0.0),
    m_spcond_ok(false),
    m_condmix_ok(false),
    m_debug(false)
{
    *this = right;
}

MixEzTransport& MixEzTransport::operator=(const MixEzTransport& right)
{
    if (&right == this) {
        return *this;
    }
    GasTransport::operator=(right);

    m_cond = right.m_cond;
    m_lambda = right.m_lambda;
    m_spcond_ok = right.m_spcond_ok;
    m_condmix_ok = right.m_condmix_ok;
    m_debug = right.m_debug;
    m_sqrt_mw = right.m_sqrt_mw;
    m_sqrt_rmw = right.m_sqrt_rmw;
    return *this;
}

Transport* MixEzTransport::duplMyselfAsTransport() const
{
    return new MixEzTransport(*this);
}

doublereal MixEzTransport::viscosity()
{
    update_T();
    update_C();

    if (m_visc_ok) {
        return m_viscmix;
    }

    doublereal vismix = 0.0;
    // update m_visc and m_phi if necessary
    if (!m_viscwt_ok) {
        updateViscosity_T();
    }

    // \Phi = \sqrt(MW)^{-1} \sqrt(MW)^{T}
    // m_spwork = \sqrt(MW)^{-1} \sqrt(MW)^{T} X
    //          = \sqrt(MW)^{-1} (\sqrt(MW)^{T} X)
    // Thus, cost is O(m_nsp)
    Eigen::Map<Eigen::VectorXd>(m_spwork.data(), m_nsp) = m_sqrt_rmw *
      (m_sqrt_mw.dot(Eigen::Map<Eigen::VectorXd>(m_molefracs.data(), m_nsp)));

    for (size_t k = 0; k < m_nsp; k++) {
        vismix += m_molefracs[k] * m_visc[k]/m_spwork[k]; //denom;
    }
    m_viscmix = vismix;
    return vismix;
}

void MixEzTransport::updateViscosity_T()
{
    if (!m_spvisc_ok) {
        updateSpeciesViscosities();
    }
    m_viscwt_ok = true;
}

void MixEzTransport::init(ThermoPhase* thermo, int mode, int log_level)
{
    GasTransport::init(thermo, mode, log_level);
    m_cond.resize(m_nsp);

    // set flags all false
    m_spcond_ok = false;
    m_condmix_ok = false;

    // calculate sqrt(MW) and sqrt(1/MW)
    m_sqrt_mw = Eigen::Map<Eigen::VectorXd>(m_mw.data(), m_nsp);
    m_sqrt_mw = m_sqrt_mw.array().sqrt();
    m_sqrt_rmw = m_sqrt_mw.cwiseInverse();

    // initialze array of Lewis numers
    m_rLeDiffCoeffs = Eigen::VectorXd::Ones(m_nsp);
    m_rLeDiffCoeffsMole = Eigen::VectorXd::Ones(m_nsp);
    m_rLeDiffCoeffsMass = Eigen::VectorXd::Ones(m_nsp);

    // calculate Lewis numbers
    updateLewisNumber();

    // release memory init by GasTransport::init
    m_wratjk.resize(0, 0); m_wratjk.data().shrink_to_fit();
    m_wratkj1.resize(0, 0); m_wratjk.data().shrink_to_fit();
}

void MixEzTransport::getMobilities(doublereal* const mobil)
{
    getMixDiffCoeffs(m_spwork.data());
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil[k] = c1 * m_spwork[k];
    }
}

doublereal MixEzTransport::thermalConductivity()
{
    update_T();
    update_C();
    if (!m_spcond_ok) {
        updateCond_T();
    }
    if (!m_condmix_ok) {
        doublereal sum1 = 0.0, sum2 = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            sum1 += m_molefracs[k] * m_cond[k];
            sum2 += m_molefracs[k] / m_cond[k];
        }
        m_lambda = 0.5*(sum1 + 1.0/sum2);
        m_condmix_ok = true;
    }
    return m_lambda;
}

void MixEzTransport::getThermalDiffCoeffs(doublereal* const dt)
{
    for (size_t k = 0; k < m_nsp; k++) {
        dt[k] = 0.0;
    }
}

void MixEzTransport::getMixDiffCoeffs(doublereal* const d)
{
    if (!m_lewis_ok) updateLewisNumber();
    // Update thermal conductivity
    thermalConductivity();
    // Compute thermal diffusivity
    double thermalDiffusivity = m_lambda /
        (m_thermo->density() * m_thermo->cp_mass());
    // Compute MixDiff from Lewis number
    Eigen::Map<Eigen::VectorXd>(d, m_nsp) =
        m_rLeDiffCoeffs * thermalDiffusivity;
}

void MixEzTransport::getMixDiffCoeffsMole(doublereal* const d)
{
    if (!m_lewis_ok) updateLewisNumber();
    // Update thermal conductivity
    thermalConductivity();
    // Compute thermal diffusivity
    double thermalDiffusivity = m_lambda /
        (m_thermo->density() * m_thermo->cp_mass());
    // Compute MixDiff from Lewis number
    Eigen::Map<Eigen::VectorXd>(d, m_nsp) =
        m_rLeDiffCoeffsMole * thermalDiffusivity;
}

void MixEzTransport::getMixDiffCoeffsMass(doublereal* const d)
{
    if (!m_lewis_ok) updateLewisNumber();
    // Update thermal conductivity
    thermalConductivity();
    // Compute thermal diffusivity
    double thermalDiffusivity = m_lambda /
        (m_thermo->density() * m_thermo->cp_mass());
    // Compute MixDiff from Lewis number
    Eigen::Map<Eigen::VectorXd>(d, m_nsp) =
        m_rLeDiffCoeffsMass * thermalDiffusivity;
}

void MixEzTransport::getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                    size_t ldx, const doublereal* const grad_X,
                                    size_t ldf, doublereal* const fluxes)
{
    update_T();
    update_C();
    getMixDiffCoeffs(m_spwork.data());
    GasTransport::getMixDiffCoeffs(m_spwork.data());
    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* y = m_thermo->massFractions();
    doublereal rhon = m_thermo->molarDensity();
    vector_fp sum(ndim,0.0);
    for (size_t n = 0; n < ndim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            fluxes[n*ldf + k] = -rhon * mw[k] * m_spwork[k] * grad_X[n*ldx + k];
            sum[n] += fluxes[n*ldf + k];
        }
    }
    // add correction flux to enforce sum to zero
    for (size_t n = 0; n < ndim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            fluxes[n*ldf + k] -= y[k]*sum[n];
        }
    }
}

void MixEzTransport::update_T()
{
    doublereal t = m_thermo->temperature();
    if (t == m_temp) {
        return;
    }
    if (t < 0.0) {
        throw CanteraError("MixEzTransport::update_T",
                           "negative temperature {}", t);
    }
    GasTransport::update_T();
    // temperature has changed, so polynomial fits will need to be redone.
    m_spcond_ok = false;
    m_bindiff_ok = false;
    m_condmix_ok = false;
}

void MixEzTransport::update_C()
{
    // signal that concentration-dependent quantities will need to be recomputed
    // before use, and update the local mole fractions.
    m_visc_ok = false;
    m_condmix_ok = false;
    m_thermo->getMoleFractions(m_molefracs.data());

    // add an offset to avoid a pure species condition
    for (size_t k = 0; k < m_nsp; k++) {
        m_molefracs[k] = std::max(Tiny, m_molefracs[k]);
    }
}

void MixEzTransport::updateCond_T()
{
    if (m_mode == CK_Mode) {
        for (size_t k = 0; k < m_nsp; k++) {
            m_cond[k] = exp(dot4(m_polytempvec, m_condcoeffs[k]));
        }
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            m_cond[k] = m_sqrt_t * dot5(m_polytempvec, m_condcoeffs[k]);
        }
    }
    m_spcond_ok = true;
    m_condmix_ok = false;
}

void MixEzTransport::updateLewisNumber()
{
    // Update thermal conductivity
    thermalConductivity();
    // Compute thermal diffusivity
    // \alpha = \frac{\lambda}{\rho \C_p}
    double thermalDiffusivity = m_lambda /
        (m_thermo->density() * m_thermo->cp_mass());
    // Compute mixutre-averaged diffusion coefficients
    GasTransport::getMixDiffCoeffs(m_rLeDiffCoeffs.data());
    GasTransport::getMixDiffCoeffsMole(m_rLeDiffCoeffsMole.data());
    GasTransport::getMixDiffCoeffsMass(m_rLeDiffCoeffsMass.data());
    // Compute reciprocals of Lewis numbers
    m_rLeDiffCoeffs /= thermalDiffusivity;
    m_rLeDiffCoeffsMole /= thermalDiffusivity;;
    m_rLeDiffCoeffsMass /= thermalDiffusivity;
    // Update ok flag
    m_lewis_ok = true;
}

}
