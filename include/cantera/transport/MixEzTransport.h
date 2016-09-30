/**
 *  @file MixEzTransport.h
 *    Headers for the MixEzTransport object, which models transport properties
 *    in ideal gas solutions using a mixture averaged approximation
 *    (see \ref tranprops and \link Cantera::MixEzTransport MixEzTransport \endlink) .
 */

// Copyright 2016  Hao Wu (wuhao@stanford.edu)

#ifndef CT_MIXEZTRAN_H
#define CT_MIXEZTRAN_H

#include "GasTransport.h"
#include "cantera/numerics/DenseMatrix.h"
#include <Eigen/Core>

namespace Cantera
{
//! Class MixEzTransport implements mixture-averaged transport properties for
//! ideal gas mixtures.
/*!
 * The model is based on that described in: R. J. Kee, M. J. Coltrin, and P.
 * Glarborg, "Chemically Reacting Flow: Theory & Practice", John Wiley & Sons,
 * 2003 and Reid, Robert C., John M. Prausnitz, and Bruce E. Poling. "The
 * properties of gases and liquids." (1987).
 *
 * The viscosity is computed using the Herning and Zipperer mixture rule
 * (kg /m /s)
 *
 * \f[
 *     \mu = \sum_k \frac{\mu_k X_k}{\sum_j \Phi_{k,j} X_j}.
 * \f]
 *
 * Here \f$ \mu_k \f$ is the viscosity of pure species \e k, and
 *
 * \f[
 *     \Phi_{k,j} = \sqrt{\frac{M_j}{M_k}} = \sqrt{M}^{1/2} \sqrt{M}^T
 * \f]
 *
 * The thermal conductivity is computed from the following mixture rule:
 * \f[
 *     \lambda = 0.5 \left( \sum_k X_k \lambda_k  + \frac{1}{\sum_k X_k/\lambda_k} \right)
 * \f]
 *
 * It's used to compute the flux of energy due to a thermal gradient
 *
 * \f[
 *     j_T =  - \lambda  \nabla T
 * \f]
 *
 * The flux of energy has units of energy (kg m2 /s2) per second per area.
 *
 * The units of lambda are W / m K which is equivalent to kg m / s^3 K.
 * @ingroup tranprops
 */
class MixEzTransport : public GasTransport
{
public:
    //! Default constructor.
    MixEzTransport();

    MixEzTransport(const MixEzTransport& right);
    MixEzTransport& operator=(const MixEzTransport& right);
    virtual Transport* duplMyselfAsTransport() const;

    //! Return the model id for transport
    /*!
     * @return cMixtureAverage
     */
    virtual int model() const {
        warn_deprecated("MixEzTransport::model",
                        "To be removed after Cantera 2.3.");
        return cMixtureAveragedEz;
    }

    virtual std::string transportType() const {
        return "MixEz";
    }

    //! Viscosity of the mixture  (kg /m /s)
    /*!
     * The viscosity is computed using the Herning and Zipperer mixture rule
     * (kg /m /s)
     *
     * \f[
     *     \mu = \sum_k \frac{\mu_k X_k}{\sum_j \Phi_{k,j} X_j}.
     * \f]
     *
     * Here \f$ \mu_k \f$ is the viscosity of pure species \e k, and
     *
     * \f[
     *     \Phi_{k,j} = \sqrt{\frac{M_j}{M_k}} = \sqrt{M}^{1/2} \sqrt{M}^T
     * \f]
     *
     * @returns the viscosity of the mixture (units =  Pa s = kg /m /s)
     *
     * @see updateViscosity_T();
     */
    virtual doublereal viscosity();

    //! Return the thermal diffusion coefficients
    /*!
     * For this approximation, these are all zero.
     *
     * @param dt  Vector of thermal diffusion coefficients. Units = kg/m/s
     */
    virtual void getThermalDiffCoeffs(doublereal* const dt);

    //! Returns the mixture thermal conductivity (W/m /K)
    /*!
     * The thermal conductivity is computed from the following mixture rule:
     * \f[
     *     \lambda = 0.5 \left( \sum_k X_k \lambda_k  + \frac{1}{\sum_k X_k/\lambda_k} \right)
     * \f]
     *
     * It's used to compute the flux of energy due to a thermal gradient
     *
     * \f[
     *     j_T =  - \lambda  \nabla T
     * \f]
     *
     * The flux of energy has units of energy (kg m2 /s2) per second per area.
     *
     * The units of lambda are W / m K which is equivalent to kg m / s^3 K.
     *
     * @returns the mixture thermal conductivity, with units of W/m/K
     */
    virtual doublereal thermalConductivity();

    //! Get the Electrical mobilities (m^2/V/s).
    /*!
     * This function returns the mobilities. In some formulations this is equal
     * to the normal mobility multiplied by Faraday's constant.
     *
     * Here, the mobility is calculated from the diffusion coefficient using the
     * Einstein relation
     *
     * \f[
     *     \mu^e_k = \frac{F D_k}{R T}
     * \f]
     *
     * @param mobil  Returns the mobilities of the species in array \c mobil.
     *               The array must be dimensioned at least as large as the
     *               number of species.
     */
    virtual void getMobilities(doublereal* const mobil);

    //! Update the internal parameters whenever the temperature has changed
    /*!
     * This is called whenever a transport property is requested if the
     * temperature has changed since the last call to update_T().
     */
    virtual void update_T();

    //! Update the internal parameters whenever the concentrations have changed
    /*!
     * This is called whenever a transport property is requested if the
     * concentrations have changed since the last call to update_C().
     */
    virtual void update_C();

    //! Returns the Mixture-averaged diffusion coefficients [m^2/s].
    /*!
     * Returns the mixture averaged diffusion coefficients for a gas,
     * appropriate for calculating the mass averaged diffusive flux with respect
     * to the mass averaged velocity using gradients of the mole fraction.
     * See GasTransport::getMixDiffCoeffs for more info.
     *
     * \f[  D_{km}' = \frac{1}{Le'} \frac{\lambda}{\rho \C_p} \f]
     *
     * @param[out] d  Vector of mixture diffusion coefficients, \f$ D_{km}' \f$ ,
     *     for each species (m^2/s). length m_nsp
     */
    virtual void getMixDiffCoeffs(doublereal* const d);

    //! Returns the mixture-averaged diffusion coefficients [m^2/s].
    //! These are the coefficients for calculating the molar diffusive fluxes
    //! from the species mole fraction gradients, computed using the
    //! corresponding Lewis number. The definition of the diffusion coefficients
    //! in the Lewis number is according to Eq. 12.176 in
    //! "Chemically Reacting Flow":
    //!
    //! \f[  D_{km}^* = \frac{1}{Le^*} \frac{\lambda}{\rho \C_p} \f]
    //!
    //! @param[out] d vector of mixture-averaged diffusion coefficients for
    //!     each species, length m_nsp.
    virtual void getMixDiffCoeffsMole(doublereal* const d);

    //! Returns the mixture-averaged diffusion coefficients [m^2/s].
    /*!
     * These are the coefficients for calculating the diffusive mass fluxes
     * from the species mass fraction gradients, computed using the
     * corresponding Lewis number. The definition of the diffusion coefficients
     * in the Lewis number is according to Eq. 12.178 in
     * "Chemically Reacting Flow":
     *
     * \f[
     *     D_{km} = \frac{1}{Le} \frac{\lambda}{\rho \C_p}
     * \f]
     *
     * @param[out] d vector of mixture-averaged diffusion coefficients for
     *     each species, length m_nsp.
     */
    virtual void getMixDiffCoeffsMass(doublereal* const d);

    //! Get the species diffusive mass fluxes wrt to the mass averaged velocity,
    //! given the gradients in mole fraction and temperature
    /*!
     * Units for the returned fluxes are kg m-2 s-1.
     *
     * The diffusive mass flux of species \e k is computed from
     * \f[
     *     \vec{j}_k = -n M_k D_k \nabla X_k.
     * \f]
     *
     * @param ndim      Number of dimensions in the flux expressions
     * @param grad_T    Gradient of the temperature (length = ndim)
     * @param ldx       Leading dimension of the grad_X array
     *                  (usually equal to m_nsp but not always)
     * @param grad_X    Gradients of the mole fraction. Flat vector with the
     *                  m_nsp in the inner loop. length = ldx * ndim
     * @param ldf       Leading dimension of the fluxes array
     *                  (usually equal to m_nsp but not always)
     * @param fluxes    Output of the diffusive mass fluxes. Flat vector with
     *                  the m_nsp in the inner loop. length = ldx * ndim
     */
    virtual void getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                  size_t ldx, const doublereal* const grad_X,
                                  size_t ldf, doublereal* const fluxes);

    virtual void init(thermo_t* thermo, int mode=0, int log_level=0);

    //! Update the Lewis number
    /*!
     * These are evaluated from the mixture-averaged diffusion coefficients in
     * GasTransport and thermal conductivity;
     */
    virtual void updateLewisNumber();

    //! Reference to Diffusion Coefficient
    Eigen::VectorXd& rLeDiffCoeffs() { return m_rLeDiffCoeffs; };
    Eigen::VectorXd& rLeDiffCoeffsMole() { return m_rLeDiffCoeffsMole; };
    Eigen::VectorXd& rLeDiffCoeffsMass() { return m_rLeDiffCoeffsMass; };

private:
    //! Calculate the pressure from the ideal gas law
    doublereal pressure_ig() const {
        return (m_thermo->molarDensity() * GasConstant *
                m_thermo->temperature());
    }

    //! Update the temperature-dependent viscosity terms.
    /**
     * Updates the array of pure species viscosities, and the weighting
     * functions in the viscosity mixture rule. The flag m_visc_ok is set to true.
     *
     * The formula for the weighting function is from Poling and Prausnitz,
     * Eq. (9-5.17):
     * \f[
     *     \Phi_{i,j} = \sqrt{\frac{M_j}{M_i}} = \sqrt{M}^{1/2} \sqrt{M}^T
     * \f]
     */
    virtual void updateViscosity_T();

    //! Update the temperature dependent parts of the species thermal
    //! conductivities
    /*!
     * These are evaluated from the polynomial fits of the temperature and are
     * assumed to be independent of pressure
     */
    void updateCond_T();

    //! vector of species thermal conductivities (W/m /K)
    /*!
     * These are used in wilke's rule to calculate the viscosity of the
     * solution. units = W /m /K = kg m /s^3 /K. length = m_kk.
     */
    vector_fp m_cond;

    //! Internal storage for the calculated mixture thermal conductivity
    /*!
     *  Units = W /m /K
     */
    doublereal m_lambda;

    //! Update boolean for the species thermal conductivities
    bool m_spcond_ok;

    //! Update boolean for the mixture rule for the mixture thermal conductivity
    bool m_condmix_ok;

    //! Update boolean for Lewis number precomputation
    bool m_lewis_ok;

    //! Debug flag - turns on more printing
    bool m_debug;

    //! Local copy of square roots of the species molecular weights.
    Eigen::VectorXd m_sqrt_mw;

    //! Local copy of square roots of the reciprocal of the species molecular
    //! weights.
    Eigen::VectorXd m_sqrt_rmw;

    //! Reciprocal of Le for MixDiffCoeffs
    Eigen::VectorXd m_rLeDiffCoeffs;

    //! Reciprocal of Le for MixDiffCoeffsMole
    Eigen::VectorXd m_rLeDiffCoeffsMole;

    //! Reciprocal of Le for MixDiffCoeffsMass
    Eigen::VectorXd m_rLeDiffCoeffsMass;
};
}
#endif
