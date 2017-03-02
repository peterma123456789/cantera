/**
 *  @file TranscriticalTransport.cpp
 *  Implementation file for class TranscriticalTransport
 *
 *  Transport parameters are calculated using corresponding states models:
 *      Binary diffusion coefficients use the generalized chart described by
 *      Takahashi, et al. and viscosity calcualtions use the Lucas method.
 *      All methods are described in Reid, Prausnitz, and Polling, "The Properties
 *      of Gases and Liquids, 4th ed., 1987 (viscosity in Ch. 9, Thermal
 *      conductivity in Ch. 10, and Diffusion coefficients in Ch. 11).
 *
 **/
#include "cantera/transport/TranscriticalTransport.h"
#include "cantera/numerics/ctlapack.h"
#include "cantera/base/utilities.h"
#include "cantera/transport/TransportParams.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/transport/MultiTransport.h"

using namespace std;

namespace Cantera
{

TranscriticalTransport::TranscriticalTransport(thermo_t* thermo)
    : MultiTransport(thermo)
{
}

void TranscriticalTransport::init(ThermoPhase* thermo, int mode, int log_level)
{
    MultiTransport::init(thermo, mode, log_level);

    // get critical properties from thermo object
    // this is quick hack
    // ideally should use GasTranportData object to set omega and dipole
    IsCrit.resize(m_nsp);
    Tcrit.resize(m_nsp);
    Pcrit.resize(m_nsp);
    Vcrit.resize(m_nsp);
    Zcrit.resize(m_nsp);
    omega.resize(m_nsp);
    dipole.resize(m_nsp);
    kappa.resize(m_nsp);

//    m_thermo->getCritTemperature(&Tcrit[0]);
//    m_thermo->getCritPressure(&Pcrit[0]);
//    m_thermo->getCritVolume(&Vcrit[0]);
//    m_thermo->getCritCompressibility(&Zcrit[0]);
//    m_thermo->getAcentricFactor(&omega[0]);
//    m_thermo->getDipoleMoment(&dipole[0]);
    ReadCriticalProperties();

    // hard-coded, otherwise have issues with Chung_HP viscosity
    //for (size_t k = 0; k < m_nsp; k++)
        //if (m_thermo->speciesName(k) == "H2")
            //omega[k] = 0.0;

    sigma_IJ.resize(m_nsp * m_nsp);
    epsOverk_IJ.resize(m_nsp * m_nsp);
    omega_IJ.resize(m_nsp * m_nsp);
    MW_IJ.resize(m_nsp * m_nsp);
    kappa_IJ.resize(m_nsp * m_nsp);

    for (size_t k = 0; k < m_nsp; k++) {
        for (size_t l = 0; l < m_nsp; l++) {
            size_t apos = k * m_nsp + l;
            sigma_IJ[apos]    = sqrt(0.809 * pow(1000.0 * Vcrit[l], 1.0 / 3.0) * 0.809 * pow(1000.0 * Vcrit[k], 1.0 / 3.0));
            epsOverk_IJ[apos] = sqrt((Tcrit[l] / 1.2593) * (Tcrit[k] / 1.2593));
            omega_IJ[apos]    = 0.5 * (omega[l] + omega[k]);
            MW_IJ[apos]       = 2.0 * m_mw[l] * m_mw[k] / (m_mw[l] + m_mw[k]);
            kappa_IJ[apos]    = sqrt(kappa[k] * kappa[l]);
        }
    }
}

void TranscriticalTransport::ReadCriticalProperties()
{
    //! \brief Reads critical properties based on JP's papaer
    for (size_t k = 0; k < m_nsp; k++) {
        if (m_thermo->speciesName(k) == "H2") {
            IsCrit[k] = 1;
            Tcrit[k] = 33.145; // K
            Pcrit[k] = 1.2964e+6; // Pa
            Vcrit[k] = 64.4834e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = -0.216;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "O2") {
            IsCrit[k] = 1;
            Tcrit[k] = 154.5800; // K
            Pcrit[k] = 5.0430e+6; // Pa
            Vcrit[k] = 73.37e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0222;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "H2O") {
            IsCrit[k] = 1;
            Tcrit[k] = 647.10; // K
            Pcrit[k] = 22.064e+6; // Pa
            Vcrit[k] = 55.95e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.3443;
            dipole[k] = 1.855;
            kappa[k] = 0.076; // modified from Chung's value
        } else if (m_thermo->speciesName(k) == "O") {
            IsCrit[k] = 0;
            Tcrit[k] = 105.28; // K
            Pcrit[k] = 7.088e+6; // Pa
            Vcrit[k] = 41.21e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "H") {
            IsCrit[k] = 0;
            Tcrit[k] = 190.82; // K
            Pcrit[k] = 31.013e+6; // Pa
            Vcrit[k] = 17.07e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "OH") {
            IsCrit[k] = 0;
            Tcrit[k] = 105.28; // K
            Pcrit[k] = 7.088e+6; // Pa
            Vcrit[k] = 41.21e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "H2O2") {
            IsCrit[k] = 0;
            Tcrit[k] = 141.34; // K
            Pcrit[k] = 4.786e+6; // Pa
            Vcrit[k] = 81.93e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "HO2") {
            IsCrit[k] = 0;
            Tcrit[k] = 141.34; // K
            Pcrit[k] = 4.786e+6; // Pa
            Vcrit[k] = 81.93e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "N2") {
            IsCrit[k] = 1;
            Tcrit[k] = 126.19; // K
            Pcrit[k] = 3.3958e+6; // Pa
            Vcrit[k] = 89.41e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0372;
            dipole[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "AR") {
            IsCrit[k] = 1;
            Tcrit[k] = 150.687; // K
            Pcrit[k] = 4.8630e+6; // Pa
            Vcrit[k] = 74.59e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "KR") {
            IsCrit[k] = 1;
            Tcrit[k] = 209.48; // K
            Pcrit[k] = 5.5250e+6; // Pa
            Vcrit[k] = 92.17e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "NE") {
            IsCrit[k] = 1;
            Tcrit[k] = 44.4918; // K
            Pcrit[k] = 2.6786e+6; // Pa
            Vcrit[k] = 41.87e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
        } else {
            cout << " WARNING -> Unknown species : " << m_thermo->speciesName(k) << ". All critical properties set to zero." << endl;
            IsCrit[k] = 0;
            Tcrit[k] = 0.0; // K
            Pcrit[k] = 0.0; // Pa
            Vcrit[k] = 0.0; // m3/kmol
            Zcrit[k] = 0.0;
            omega[k] = 0.0;
            dipole[k] = 0.0;
        }
    }
}

// Mass diffusivity models
//void TranscriticalTransport::getBinaryDiffCoeffs(const size_t ld, doublereal* const d)
//{
//    doublereal P_corr_ij, Tr_ij, Pr_ij;
//
//    update_T();
//
//    updateDiff_T();
//
//    if (ld < m_nsp) {
//        throw CanteraError("TranscriticalTransport::getBinaryDiffCoeffs()", "ld is too small");
//    }
//
//    doublereal rp = 1.0/m_thermo->pressure();
//    for (size_t i = 0; i < m_nsp; i++) {
//        for (size_t j = 0; j < m_nsp; j++) {
//            // Add an offset to avoid a condition where x_i and x_j both equal
//            //   zero (this would lead to Pr_ij = Inf):
//            doublereal x_i_mix = std::max(Tiny, m_molefracs[i]);
//            doublereal x_j_mix = std::max(Tiny, m_molefracs[j]);
//
//            // Weight mole fractions of i and j so that X_i + X_j = 1.0:
//            doublereal x_i = x_i_mix/(x_i_mix + x_j_mix);
//            doublereal x_j = x_j_mix/(x_i_mix + x_j_mix);
//
//            //Calculate Tr and Pr based on mole-fraction-weighted crit constants:
//            Tr_ij = m_temp/(x_i*Tcrit[i] + x_j*Tcrit[j]);
//            Pr_ij = m_thermo->pressure()/(x_i*Pcrit[i] + x_j*Pcrit[j]);
//
//            if (Pr_ij < 0.1) {
//                // If pressure is low enough, no correction is needed:
//                P_corr_ij = 1;
//            }else {
//                // Otherwise, calculate the parameters for Takahashi correlation
//                //   by interpolating on Pr_ij:
//                P_corr_ij = setPcorr(Pr_ij, Tr_ij);
//
//                // If the reduced temperature is too low, the correction factor
//                //  P_corr_ij will be < 0:
//                if (P_corr_ij < 0) {
//                    P_corr_ij = Tiny;
//                }
//            }
//
//            // Multiply the standard low-pressure binary diffusion coefficient
//            //   (m_bdiff) by the Takahashi correction factor P_corr_ij:
//            d[ld * j + i] = P_corr_ij * rp * m_bdiff(i, j);
//        }
//    }
//}

//void TranscriticalTransport::updateDiff_T()
//{
//    update_T();
//
//    // evaluate binary diffusion coefficients at unit pressure
//    size_t ic = 0;
//    if (m_mode == CK_Mode) {
//        for (size_t i = 0; i < m_nsp; i++) {
//            for (size_t j = i; j < m_nsp; j++) {
//                m_bdiff(i, j) = exp(dot4(m_polytempvec, m_diffcoeffs[ic]));
//                m_bdiff(j, i) = m_bdiff(i, j);
//                ic++;
//            }
//        }
//    } else {
//        for (size_t i = 0; i < m_nsp; i++) {
//            for (size_t j = i; j < m_nsp; j++) {
//                m_bdiff(i, j) = m_temp * m_sqrt_t * dot5(m_polytempvec, m_diffcoeffs[ic]);
//                m_bdiff(j, i) = m_bdiff(i, j);
//                ic++;
//            }
//        }
//    }
//
//    // Correct the binary diffusion coefficients for high-pressure effects
//    doublereal P_corr_ij, Tr_ij, Pr_ij;
//    for (size_t i = 0; i < m_nsp; i++) {
//        for (size_t j = 0; j < m_nsp; j++) {
//            if (j != i) {
//                // Add an offset to avoid a condition where x_i and x_j both equal
//                //   zero (this would lead to Pr_ij = Inf):
//                doublereal x_i_mix = std::max(Tiny, m_molefracs[i]);
//                doublereal x_j_mix = std::max(Tiny, m_molefracs[j]);
//
//                doublereal x_i = x_i_mix / (x_i_mix + x_j_mix);
//                doublereal x_j = x_j_mix / (x_i_mix + x_j_mix);
//
//                Tr_ij = m_temp / (x_i * Tcrit[i] + x_j * Tcrit[j]);
//                Pr_ij = m_thermo->pressure() / (x_i * Pcrit[i] + x_j * Pcrit[j]);
//
//                if (Pr_ij < 0.1) {
//                    P_corr_ij = 1;
//                } else {
//                    P_corr_ij = setPcorr(Pr_ij, Tr_ij);
//                    //if (P_corr_ij < 0) {
//                    //    P_corr_ij = Tiny;
//                    //}
//                }
//
//                if (P_corr_ij > 1.07) printf("%s, %s, Pr = %g, Tr = %g, P_corr = %g\n", m_thermo->speciesName(i).c_str(), m_thermo->speciesName(j).c_str(), Pr_ij, Tr_ij, P_corr_ij);
//
//                m_bdiff(i, j) *= max(P_corr_ij, 0.2);
//            }
//        }
//    }
//
//    m_bindiff_ok = true;
//}

void TranscriticalTransport::updateThermal_T()
{
    MultiTransport::updateThermal_T();

    // Correct the binary diffusion coefficients for high-pressure effects
    doublereal P_corr_ij, Tr_ij, Pr_ij;
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            if (j != i) {
                // Add an offset to avoid a condition where x_i and x_j both equal
                //   zero (this would lead to Pr_ij = Inf):
                doublereal x_i_mix = std::max(Tiny, m_molefracs[i]);
                doublereal x_j_mix = std::max(Tiny, m_molefracs[j]);

                doublereal x_i = x_i_mix / (x_i_mix + x_j_mix);
                doublereal x_j = x_j_mix / (x_i_mix + x_j_mix);

                Tr_ij = m_temp / (x_i * Tcrit[i] + x_j * Tcrit[j]);
                Pr_ij = m_thermo->pressure() / (x_i * Pcrit[i] + x_j * Pcrit[j]);

                if (Pr_ij < 0.1) {
                    P_corr_ij = 1;
                } else {
                    P_corr_ij = setPcorr(Pr_ij, Tr_ij);
                    //if (P_corr_ij < 0) {
                    //    P_corr_ij = Tiny;
                    //}
                }

                if (P_corr_ij > 1.07) printf("%s, %s, Pr = %g, Tr = %g, P_corr = %g\n", m_thermo->speciesName(i).c_str(), m_thermo->speciesName(j).c_str(), Pr_ij, Tr_ij, P_corr_ij);

                // we limit the correction, b/c otherwise have numerical issues
                // for flamelet calculations
                m_bdiff(i, j) *= max(P_corr_ij, 0.4);
            }
        }
    }
}

//void TranscriticalTransport::getMultiDiffCoeffs_Takahashi(const size_t ld, doublereal* const d)
//{
//    // update the mole fractions
//    update_C();
//
//    // update the binary diffusion coefficients
//    update_T();
//    updateThermal_T();
//
//    // Correct the binary diffusion coefficients for high-pressure effects
//    doublereal P_corr_ij, Tr_ij, Pr_ij;
//    for (size_t i = 0; i < m_nsp; i++) {
//        for (size_t j = 0; j < m_nsp; j++) {
//            if (j != i) {
//                // Add an offset to avoid a condition where x_i and x_j both equal
//                //   zero (this would lead to Pr_ij = Inf):
//                doublereal x_i_mix = std::max(Tiny, m_molefracs[i]);
//                doublereal x_j_mix = std::max(Tiny, m_molefracs[j]);
//
//                doublereal x_i = x_i_mix / (x_i_mix + x_j_mix);
//                doublereal x_j = x_j_mix / (x_i_mix + x_j_mix);
//
//                Tr_ij = m_temp / (x_i * Tcrit[i] + x_j * Tcrit[j]);
//                Pr_ij = m_thermo->pressure() / (x_i * Pcrit[i] + x_j * Pcrit[j]);
//
//                if (Pr_ij < 0.1) {
//                    P_corr_ij = 1;
//                } else {
//                    P_corr_ij = setPcorr(Pr_ij, Tr_ij);
//                    //if (P_corr_ij < 0) {
//                    //    P_corr_ij = Tiny;
//                    //}
//                }
//
//                if (P_corr_ij > 1.07) printf("%s, %s, Pr = %g, Tr = %g, P_corr = %g\n", m_thermo->speciesName(i).c_str(), m_thermo->speciesName(j).c_str(), Pr_ij, Tr_ij, P_corr_ij);
//
//                m_bdiff(i, j) *= max(P_corr_ij, 0.2);
//            }
//        }
//    }
//    m_bindiff_ok = false;  // m_bdiff is overwritten by the above routine.
//
//    // Having corrected m_bdiff for pressure and concentration effects, the
//    // routine now procedes the same as in the low-pressure case:
//
//    // evaluate L0000 if the temperature or concentrations have
//    // changed since it was last evaluated.
//    if (!m_l0000_ok) {
//        eval_L0000(DATA_PTR(m_molefracs));
//    }
//
//    // invert L00,00
//    int ierr = invert(m_Lmatrix, m_nsp);
//    if (ierr != 0) {
//        throw CanteraError("TranscriticalTransport::getMultiDiffCoeffs",
//                           string(" invert returned ierr = ") + int2str(ierr));
//    }
//    m_l0000_ok = false;           // matrix is overwritten by inverse
//    m_lmatrix_soln_ok = false;
//
//    doublereal prefactor = 16.0 * m_temp * m_thermo->meanMolecularWeight() / (25.0 * m_thermo->pressure());
//    doublereal c;
//
//    for (size_t i = 0; i < m_nsp; i++) {
//        for (size_t j = 0; j < m_nsp; j++) {
//            c = prefactor / m_mw[j];
//            d[ld * j + i] = c * m_molefracs[i] *
//                            (m_Lmatrix(i, j) - m_Lmatrix(i, i));
//        }
//    }
//}

//void TranscriticalTransport::getMixDiffCoeffs_Takahashi(doublereal* const d)
//{
//    update_T();
//    update_C();
//
//    // update the binary diffusion coefficients if necessary
//    if (!m_bindiff_ok) {
//        updateDiff_T();
//    }
//
//    // Correct the binary diffusion coefficients for high-pressure effects
//    doublereal P_corr_ij, Tr_ij, Pr_ij;
//    for (size_t i = 0; i < m_nsp; i++) {
//        for (size_t j = 0; j < m_nsp; j++) {
//            if (j != i) {
//                // Add an offset to avoid a condition where x_i and x_j both equal
//                //   zero (this would lead to Pr_ij = Inf):
//                doublereal x_i_mix = std::max(Tiny, m_molefracs[i]);
//                doublereal x_j_mix = std::max(Tiny, m_molefracs[j]);
//
//                doublereal x_i = x_i_mix / (x_i_mix + x_j_mix);
//                doublereal x_j = x_j_mix / (x_i_mix + x_j_mix);
//
//                Tr_ij = m_temp / (x_i * Tcrit[i] + x_j * Tcrit[j]);
//                Pr_ij = m_thermo->pressure() / (x_i * Pcrit[i] + x_j * Pcrit[j]);
//
//                if (Pr_ij < 0.1) {
//                    P_corr_ij = 1;
//                } else {
//                    P_corr_ij = setPcorr(Pr_ij, Tr_ij);
//                    //if (P_corr_ij < 0) {
//                    //    P_corr_ij = Tiny;
//                    //}
//                }
//
//                if (P_corr_ij > 1.07) printf("%s, %s, Pr = %g, Tr = %g, P_corr = %g\n", m_thermo->speciesName(i).c_str(), m_thermo->speciesName(j).c_str(), Pr_ij, Tr_ij, P_corr_ij);
//
//                m_bdiff(i, j) *= max(P_corr_ij, 0.4);
//            }
//        }
//    }
//    m_bindiff_ok = false;  // m_bdiff is overwritten by the above routine.
//
//    // Having corrected m_bdiff for pressure and concentration effects, the
//    // routine now procedes the same as in the low-pressure case:
//
//    doublereal mmw = m_thermo->meanMolecularWeight();
//    doublereal sumxw = 0.0;
//    doublereal p = m_thermo->pressure();
//    if (m_nsp == 1) {
//        d[0] = m_bdiff(0,0) / p;
//    } else {
//        for (size_t k = 0; k < m_nsp; k++) {
//            sumxw += m_molefracs[k] * m_mw[k];
//        }
//        for (size_t k = 0; k < m_nsp; k++) {
//            double sum2 = 0.0;
//            for (size_t j = 0; j < m_nsp; j++) {
//                if (j != k) {
//                    sum2 += m_molefracs[j] / m_bdiff(j,k);
//                }
//            }
//            if (sum2 <= 0.0) {
//                d[k] = m_bdiff(k,k) / p;
//            } else {
//                d[k] = (sumxw - m_molefracs[k] * m_mw[k])/(p * mmw * sum2);
//            }
//        }
//    }
//}

// Viscosity models
doublereal TranscriticalTransport::viscosity_Lucas_HP()
{
    // Calculate the high-pressure mixture viscosity, based on the Lucas method.
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal Tc_mix = 0.;
    doublereal Pc_mix_n = 0.;
    doublereal Pc_mix_d = 0.;
    doublereal MW_mix = m_thermo->meanMolecularWeight();
    doublereal Tc, Zc, Tr, Afac, Z1m, Z2m;
    doublereal MW_H = m_mw[0];
    doublereal MW_L = m_mw[0];
    doublereal x_H = molefracs[0];
    doublereal FP_mix_o = 0.;
    doublereal FQ_mix_o = 0.;
    doublereal tKelvin = m_thermo->temperature();
    // for gas phase, we don't have this method implemented
    //doublereal Pvp_mix = m_thermo->satPressure(tKelvin);

    for (size_t i = 0; i < m_nsp; i++) {
        // Calculate pure-species critical constants and add their contribution
        // to the mole-fraction-weighted mixture averages:
        Tc = Tcrit[i];
        Tr = tKelvin / Tc;
        Zc = Zcrit[i];
        Tc_mix += Tc * molefracs[i];
        Pc_mix_n += molefracs[i] * Zc;         //numerator
        Pc_mix_d += molefracs[i] * Vcrit[i]; //denominator

        // Need to calculate ratio of heaviest to lightest species:
        if (m_mw[i] > MW_H) {
            MW_H = m_mw[i];
            x_H = molefracs[i];
        } else if (m_mw[i] < MW_L) {
            MW_L = m_mw[i];
        }

        // Calculate reduced dipole moment for polar correction term:
        doublereal mu_ri = 52.46 * 100000 * m_dipole(i, i) * m_dipole(i, i) * Pcrit[i] / (Tc * Tc);
        if (mu_ri < 0.022) {
            FP_mix_o += molefracs[i];
        } else if (mu_ri < 0.075) {
            FP_mix_o += molefracs[i]*(1. + 30.55*pow(0.292 - Zc, 1.72));
        } else {
            FP_mix_o += molefracs[i] * (1. + 30.55 * pow(0.292 - Zc, 1.72) * fabs(0.96 + 0.1 * (Tr - 0.7)));
        }

        // Calculate contribution to quantum correction term.
        // SCD Note:  This assumes the species of interest (He, H2, and D2) have
        //   been named in this specific way.  They are perhaps the most obvious
        //   names, but it would of course be preferred to have a more general
        //   approach, here.
        std::vector<std::string> spnames = m_thermo->speciesNames();
        if (spnames[i] == "He") {
            FQ_mix_o += molefracs[i] * FQ_i(1.38, Tr, m_mw[i]);
        } else if (spnames[i] == "H2") {
            FQ_mix_o += molefracs[i] * (FQ_i(0.76, Tr, m_mw[i]));
        } else if (spnames[i] == "D2") {
            FQ_mix_o += molefracs[i] * (FQ_i(0.52, Tr, m_mw[i]));
        } else {
            FQ_mix_o += molefracs[i];
        }
    }

    double Tr_mix = tKelvin/Tc_mix;
    double Pc_mix = GasConstant*Tc_mix*Pc_mix_n/Pc_mix_d;
    double Pr_mix = m_thermo->pressure()/Pc_mix;
    double ratio = MW_H/MW_L;

    double ksi = pow(GasConstant * Tc_mix * 3.6277 * pow(10.0, 53.0) / (pow(MW_mix, 3) * pow(Pc_mix, 4)), 1.0 / 6.0);

    if (ratio > 9 && x_H > 0.05 && x_H < 0.7) {
        Afac = 1 - 0.01*pow(ratio,0.87);
    } else {
        Afac = 1;
    }
    FQ_mix_o *= Afac;

    // Calculate Z1m
    Z1m = (0.807 * pow(Tr_mix, 0.618) - 0.357 * exp(-0.449 * Tr_mix) + 0.340 * exp(-4.058 * Tr_mix) + 0.018) * FP_mix_o * FQ_mix_o;

    // Calculate Z2m:
//    if (Tr_mix <= 1.0){
//        if (Pr_mix < Pvp_mix/Pc_mix) {
//            doublereal alpha = 3.262 + 14.98*pow(Pr_mix,5.508);
//            doublereal beta = 1.390 + 5.746*Pr_mix;
//            Z2m = 0.600 + 0.760 * pow(Pr_mix, alpha) + (0.6990 * pow(Pr_mix, beta) - 0.60) * (1 - Tr_mix);
//        } else {
//            throw CanteraError("TranscriticalTransport::viscosity",
//                               "State is outside the limits of the Lucas model, Tr <= 1");
//        }
//    } else if ((Tr_mix > 1.0) && (Tr_mix < 40.0)) {
//        if ((Pr_mix > 0.0) && (Pr_mix <= 100.0)) {
            doublereal a_fac = 0.001245*exp(5.1726*pow(Tr_mix,-0.3286))/Tr_mix;
            doublereal b_fac = a_fac*(1.6553*Tr_mix - 1.2723);
            doublereal c_fac = 0.4489*exp(3.0578*pow(Tr_mix,-37.7332))/Tr_mix;
            doublereal d_fac = 1.7368*exp(2.2310*pow(Tr_mix,-7.6351))/Tr_mix;
            doublereal f_fac = 0.9425*exp(-0.1853*pow(Tr_mix,0.4489));

            Z2m = Z1m * (1 + a_fac * pow(Pr_mix, 1.3088) / (b_fac * pow(Pr_mix, f_fac) + pow(1 + c_fac * pow(Pr_mix, d_fac), -1)));
//        } else {
//            throw CanteraError("TranscriticalTransport::viscosity",
//                               "State is outside the limits of the Lucas model, 1.0 < Tr < 40");
//        }
//    } else {
//        throw CanteraError("TranscriticalTransport::viscosity",
//                           "State is outside the limits of the Lucas model, Tr > 40");
//    }

    // Calculate Y:
    doublereal Y = Z2m/Z1m;

    // Return the viscosity:
    return Z2m * (1 + (FP_mix_o - 1) * pow(Y, -3)) * (1 + (FQ_mix_o - 1) * (1 / Y - 0.007 * pow(log(Y), 4))) / (ksi * FP_mix_o * FQ_mix_o);
}

doublereal TranscriticalTransport::viscosity_Chung_LP()
{
    m_thermo->getMoleFractions(DATA_PTR(m_molefracs));

    doublereal T_in = m_thermo->temperature();

    doublereal sigma_M    = 0.0;
    doublereal epsOverk_M = 0.0;
    doublereal omega_M    = 0.0;
    doublereal mu_M       = 0.0;
    doublereal kappa_M    = 0.0; // polar correction, neglected here for all species
    doublereal MW_M       = 0.0;

    for (size_t k = 0; k < m_nsp; k++) {
        for (size_t l = 0; l < m_nsp; l++) {
            size_t apos = k * m_nsp + l;
            doublereal X_X = m_molefracs[l] * m_molefracs[k];
            sigma_M += X_X * pow(sigma_IJ[apos], 3.0);
            MW_M += X_X * epsOverk_IJ[apos] * pow(sigma_IJ[apos], 2) * sqrt(MW_IJ[apos]);
            epsOverk_M += X_X * epsOverk_IJ[apos] * pow(sigma_IJ[apos], 3);
            omega_M += X_X * omega_IJ[apos] * pow(sigma_IJ[apos], 3.0);
            kappa_M += X_X * kappa_IJ[apos];
            if (fabs(sigma_IJ[apos]) > 1.0e-20) {
                mu_M += X_X * pow(dipole[k] * dipole[l], 2.0) / pow(sigma_IJ[apos], 3.0);
            }
        }
    }

    sigma_M = pow(sigma_M, 1.0 / 3.0);
    epsOverk_M = epsOverk_M / pow(sigma_M, 3);
    omega_M = omega_M / pow(sigma_M, 3.0);
    mu_M = pow(mu_M * pow(sigma_M, 3.0), 0.25);
    MW_M = pow(MW_M / (epsOverk_M * pow(sigma_M, 2)), 2);

    doublereal VCrit_M = pow((sigma_M / 0.809), 3.0);
    doublereal TCrit_M = 1.2593 * epsOverk_M;
    doublereal TStar_M = T_in / epsOverk_M;
    doublereal mu_RM = 131.3 * mu_M / sqrt(VCrit_M * TCrit_M);
    doublereal Neufeld = 1.16145 * pow(TStar_M, -0.14874) + 0.52487 * exp(-0.7732 * TStar_M) + 2.16178 * exp(-2.43787 * TStar_M); // Eq (9-4-3) in Poling
    doublereal Fc = 1.0 - 0.2756 * omega_M + 0.059035 * pow(mu_RM, 4) + kappa_M;

    doublereal mixViscosity = 40.785 * Fc * sqrt(MW_M * T_in) / (Neufeld * pow(VCrit_M, 2.0 / 3.0)) / 10.0e6;

    // Return the viscosity:
    return mixViscosity;
}

doublereal TranscriticalTransport::viscosity_Chung_HP()
{
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal T_in = m_thermo->temperature();
    doublereal rho_in = m_thermo->density();

    doublereal sigma_M    = 0.0;
    doublereal epsOverk_M = 0.0;
    doublereal omega_M    = 0.0;
    doublereal mu_M       = 0.0;
    doublereal kappa_M    = 0.0; // polar correction, neglected here for all species
    doublereal MW_M       = 0.0;

    for (size_t k = 0; k < nsp; k++) {
        for (size_t l = 0; l < nsp; l++) {
            size_t apos = k * nsp + l;
            doublereal X_X = molefracs[l] * molefracs[k];
            sigma_M += X_X * pow(sigma_IJ[apos], 3.0);
            MW_M += X_X * epsOverk_IJ[apos] * pow(sigma_IJ[apos], 2) * sqrt(MW_IJ[apos]);
            epsOverk_M += X_X * epsOverk_IJ[apos] * pow(sigma_IJ[apos], 3);
            omega_M += X_X * 0.0 * pow(sigma_IJ[apos], 3.0);
            kappa_M += X_X * 0.0;
            if (fabs(sigma_IJ[apos]) > 1.0e-20) {
                mu_M += X_X * pow(0.0 * 0.0, 2.0) / pow(sigma_IJ[apos], 3.0);
            }
            //omega_M += X_X * omega_IJ[apos] * pow(sigma_IJ[apos], 3.0);
            //kappa_M += X_X * kappa_IJ[apos];
            //if (fabs(sigma_IJ[apos]) > 1.0e-20) {
            //    mu_M += X_X * pow(dipole[k] * dipole[l], 2.0) / pow(sigma_IJ[apos], 3.0);
            //}
        }
    }

    sigma_M = pow(sigma_M, 1.0 / 3.0);
    epsOverk_M = epsOverk_M / pow(sigma_M, 3);
    omega_M = omega_M / pow(sigma_M, 3.0);
    mu_M = pow(mu_M * pow(sigma_M, 3.0), 0.25);
    MW_M = pow(MW_M / (epsOverk_M * pow(sigma_M, 2)), 2);

    doublereal VCrit_M = pow((sigma_M / 0.809), 3.0);
    doublereal TCrit_M = 1.2593 * epsOverk_M;
    doublereal TStar_M = T_in / epsOverk_M;
    doublereal mu_RM = 131.3 * mu_M / sqrt(VCrit_M * TCrit_M);
    doublereal Neufeld = 1.16145 * pow(TStar_M, -0.14874) + 0.52487 * exp(-0.7732 * TStar_M) + 2.16178 * exp(-2.43787 * TStar_M); // Eq (9-4-3) in Poling
    doublereal Fc = 1.0 - 0.2756 * omega_M + 0.059035 * pow(mu_RM, 4) + kappa_M;

    // Chung_HP part starts...
    const size_t nCoeff = 10;
    static const double a1[nCoeff] = {6.324, 1.210e-3, 5.283, 6.623, 19.745, -1.900, 24.275, 0.7972, -0.2382, 0.06863};
    static const double b1[nCoeff] = {50.412, -1.154e-3, 254.209, 38.096, 7.630, -12.537, 3.450, 1.1170, 0.0677, 0.3479};
    static const double c1[nCoeff] = {-51.680, -6.257e-3, -168.480, -8.464, -14.354, 4.985, -11.291, 0.01235, -0.8163, 0.5926};
    static const double d1[nCoeff] = {1189.000, 0.03728, 3898.0, 31.42, 31.53, -18.15, 69.35, -4.117, 4.025, -0.727};
    vector_fp E(nCoeff);
    for (size_t i = 0; i < nCoeff; i++)
        E[i] = a1[i] + b1[i] * omega_M + c1[i] * pow(mu_RM, 4.0) + d1[i] * kappa_M;

    doublereal y = rho_in / (1000 * MW_M) * VCrit_M / 6.0;
    //jph the units of density are strange density is mol/cm^3. In poling it is cm^3/mol but I think from Chung 1988 they are typical kg/m^3

    doublereal G1 = (1.0 - 0.5 * y) / pow((1.0 - y), 3);
    doublereal G2 = (E[0] / y * (1.0 - exp(-E[3] * y)) + E[1] * G1 * exp(E[4] * y) + E[2] * G1) / (E[0] * E[3] + E[1] + E[2]);
    //G2 = max(G2, 0.8);

    doublereal mustarstar = E[6] * pow(y, 2) * G2 * exp(E[7] + E[8] / TStar_M + E[9] / pow(TStar_M, 2));
    doublereal mustar = pow(TStar_M, 0.5) / Neufeld * (Fc * (pow(G2, -1) + E[5] * y)) + mustarstar;

//printf("----------->mustar = %f, TStar_M = %f, Neufeld = %f, Fc = %f, G2 = %f, E[5] = %f, y = %f, mustarstar = %f\n", mustar, TStar_M, Neufeld, Fc, G2, E[5], y, mustarstar);

//printf("mustar = %f, MW_M = %f, TCrit_M = %f, VCrit_M = %f\n", mustar, MW_M, TCrit_M, VCrit_M);
    doublereal mixViscosity = (36.644 * mustar * sqrt(MW_M * TCrit_M) / pow(VCrit_M, 2.0 / 3.0)) / 10.0e6;

    // Return the viscosity:
    return mixViscosity;
}

// Thermal conductivity models
doublereal TranscriticalTransport::thermalConductivity_ElyHanley_HP()
{
    //  Method of Ely and Hanley:
    update_T();

    doublereal Lprime_m = 0.0;

    const doublereal c1 = 1. / 16.04;

    vector_fp cp_0_R(m_nsp);
    m_thermo->getCp_R_ref(&cp_0_R[0]);

    std::vector<doublereal> L_i(m_nsp);
    std::vector<doublereal> f_i(m_nsp);
    std::vector<doublereal> h_i(m_nsp);
    std::vector<doublereal> V_k(m_nsp);

    m_thermo->getPartialMolarVolumes(&V_k[0]);

    doublereal L_i_min = BigNumber;

    for (size_t i = 0; i < m_nsp; i++) {
        doublereal Tc_i = Tcrit[i];
        doublereal Vc_i = Vcrit[i];
        doublereal T_r = m_thermo->temperature() / Tc_i;
        doublereal V_r = V_k[i] / Vc_i;
        doublereal T_p = std::min(T_r, 2.0);
        doublereal V_p = std::max(0.5, std::min(V_r, 2.0));

        // Calculate variables for density-independent component:
        doublereal theta_p = 1.0 + (m_w_ac[i] - 0.011)*(0.56553
            - 0.86276*log(T_p) - 0.69852/T_p);
        doublereal phi_p = (1.0 + (m_w_ac[i] - 0.011)*(0.38560
            - 1.1617*log(T_p)))*0.288/Zcrit[i];
        doublereal f_fac = Tc_i*theta_p/190.4;
        doublereal h_fac = 1000*Vc_i*phi_p/99.2;
        doublereal T_0 = m_temp/f_fac;
        doublereal mu_0 = 1e-7*(2.90774e6/T_0 - 3.31287e6*pow(T_0,-2./3.)
            + 1.60810e6*pow(T_0,-1./3.) - 4.33190e5 + 7.06248e4*pow(T_0,1./3.)
            - 7.11662e3*pow(T_0,2./3.) + 4.32517e2*T_0 - 1.44591e1*pow(T_0,4./3.)
            + 2.03712e-1*pow(T_0,5./3.));
        doublereal H = sqrt(f_fac*16.04/m_mw[i])*pow(h_fac,-2./3.);
        doublereal mu_i = mu_0*H*m_mw[i]*c1;
        L_i[i] = mu_i*1.32*GasConstant*(cp_0_R[i] - 2.5)/m_mw[i];
        L_i_min = min(L_i_min,L_i[i]);
        // Calculate variables for density-dependent component:
        doublereal theta_s = 1 + (m_w_ac[i] - 0.011)*(0.09057 - 0.86276*log(T_p)
            + (0.31664 - 0.46568/T_p)*(V_p - 0.5));
        doublereal phi_s = (1 + (m_w_ac[i] - 0.011)*(0.39490*(V_p - 1.02355)
            - 0.93281*(V_p - 0.75464)*log(T_p)))*0.288/Zcrit[i];
        f_i[i] = Tc_i*theta_s/190.4;
        h_i[i] = 1000*Vc_i*phi_s/99.2;
    }

    doublereal h_m = 0;
    doublereal f_m = 0;
    doublereal mw_m = 0;
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            // Density-independent component:
            doublereal L_ij = 2*L_i[i]*L_i[j]/(L_i[i] + L_i[j] + Tiny);
            Lprime_m += m_molefracs[i]*m_molefracs[j]*L_ij;
            // Additional variables for density-dependent component:
            doublereal f_ij = sqrt(f_i[i]*f_i[j]);
            doublereal h_ij = 0.125*pow(pow(h_i[i],1./3.) + pow(h_i[j],1./3.),3.);
            doublereal mw_ij_inv = (m_mw[i] + m_mw[j])/(2*m_mw[i]*m_mw[j]);
            f_m += m_molefracs[i]*m_molefracs[j]*f_ij*h_ij;
            h_m += m_molefracs[i]*m_molefracs[j]*h_ij;
            mw_m += m_molefracs[i]*m_molefracs[j]*sqrt(mw_ij_inv*f_ij)*pow(h_ij,-4./3.);
        }
    }

    f_m = f_m/h_m;
    mw_m = pow(mw_m,-2.)*f_m*pow(h_m,-8./3.);

    doublereal rho_0 = 16.04*h_m/(1000*m_thermo->molarVolume());
    doublereal T_0 = m_temp/f_m;
    doublereal mu_0 = 1e-7*(2.90774e6/T_0 - 3.31287e6*pow(T_0,-2./3.)
                + 1.60810e6*pow(T_0,-1./3.) - 4.33190e5 + 7.06248e4
                *pow(T_0,1./3.) - 7.11662e3*pow(T_0,2./3.) + 4.32517e2*T_0
                - 1.44591e1*pow(T_0,4./3.) + 2.03712e-1*pow(T_0,5./3.));
    doublereal L_1m = 1944*mu_0;
    doublereal L_2m = (-2.5276e-4 + 3.3433e-4*pow(1.12 - log(T_0/1.680e2),2))*rho_0;
    doublereal L_3m = exp(-7.19771 + 85.67822/T_0)*(exp((12.47183
                - 984.6252*pow(T_0,-1.5))*pow(rho_0,0.1) + (rho_0/0.1617 - 1)
                *sqrt(rho_0)*(0.3594685 + 69.79841/T_0 - 872.8833*pow(T_0,-2))) - 1.)*1e-3;
    doublereal H_m = sqrt(f_m*16.04/mw_m)*pow(h_m,-2./3.);
    doublereal Lstar_m = H_m*(L_1m + L_2m + L_3m);

    return Lprime_m + Lstar_m;
}

doublereal TranscriticalTransport::thermalConductivity_Chung_LP()
{
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal T_in = m_thermo->temperature();
    //doublereal rho_in = m_thermo->density();

    doublereal sigma_M    = 0.0;
    doublereal epsOverk_M = 0.0;
    doublereal omega_M    = 0.0;
    doublereal mu_M       = 0.0;
    doublereal kappa_M    = 0.0; // polar correction, neglected here for all species
    doublereal MW_M       = 0.0;

    for (size_t k = 0; k < nsp; k++) {
        for (size_t l = 0; l < nsp; l++) {
            size_t apos = k * nsp + l;
            doublereal X_X = molefracs[l] * molefracs[k];
            sigma_M += X_X * pow(sigma_IJ[apos], 3.0);
            MW_M += X_X * epsOverk_IJ[apos] * pow(sigma_IJ[apos], 2) * sqrt(MW_IJ[apos]);
            epsOverk_M += X_X * epsOverk_IJ[apos] * pow(sigma_IJ[apos], 3);
            omega_M += X_X * omega_IJ[apos] * pow(sigma_IJ[apos], 3.0);
            kappa_M += X_X * kappa_IJ[apos];
            if (fabs(sigma_IJ[apos]) > 1.0e-20) {
                mu_M += X_X * pow(dipole[k] * dipole[l], 2.0) / pow(sigma_IJ[apos], 3.0);
            }
        }
    }

    sigma_M = pow(sigma_M, 1.0 / 3.0);
    epsOverk_M = epsOverk_M / pow(sigma_M, 3);
    omega_M = omega_M / pow(sigma_M, 3.0);
    mu_M = pow(mu_M * pow(sigma_M, 3.0), 0.25);
    MW_M = pow(MW_M / (epsOverk_M * pow(sigma_M, 2)), 2);

    doublereal VCrit_M = pow((sigma_M / 0.809), 3.0);
    doublereal TCrit_M = 1.2593 * epsOverk_M;
    //doublereal TStar_M = T_in / epsOverk_M;
    doublereal mu_RM = 131.3 * mu_M / sqrt(VCrit_M * TCrit_M);
    //doublereal Neufeld = 1.16145 * pow(TStar_M, -0.14874) + 0.52487 * exp(-0.7732 * TStar_M) + 2.16178 * exp(-2.43787 * TStar_M); // Eq (9-4-3) in Poling
    //doublereal Fc = 1.0 - 0.2756 * omega_M + 0.059035 * pow(mu_RM, 4) + kappa_M;

    // above is repeated as in viscosity_Chung_HP

    doublereal mu = viscosity_Chung_LP();
    //doublereal mu = viscosity_Chung_HP();

    // Chung's conductivity model is used
    // Currently, only the low pressure correction is implemented and verified. High pressure contains the same inconsistency as the chung viscosity. This needs to be addressed.
    const size_t nCoeff = 7;
    static const double a1[nCoeff] = {2.4166, -0.50924, 6.6107, 14.543, 0.79274, -5.8634, 91.089};
    static const double b1[nCoeff] = {0.74824, -1.5094, 5.6207, -8.9139, 0.82019, 12.801, 128.11};
    static const double c1[nCoeff] = {-0.91858, -49.991, 64.760, -5.6379, -0.69369, 9.5893, -54.217};
    static const double d1[nCoeff] = {121.72, 69.983, 27.039, 74.344, 6.3173, 65.529, 523.81};
    // compute B_i coefficients
    vector_fp B_coef(nCoeff);
    for (size_t i = 0; i < nCoeff; i++)
        B_coef[i] = a1[i] + b1[i] * omega_M + c1[i] * pow(mu_RM, 4.0) + d1[i] * kappa_M;

    //doublereal y = rho_in / (1000 * MW_M) * VCrit_M / 6.0;
    //jph the units of density are strange density is mol/cm^3. In poling it is cm^3/mol but I think from Chung 1988 they are typical kg/m^3

    //doublereal G1 = (1.0 - 0.5 * y) / pow((1.0 - y), 3);
    //doublereal G2 = (B_coef[0] / y) * (1.0 - exp(-B_coef[3] * y) + B_coef[1] * G1 * exp(B_coef[4] * y) + B_coef[2] * G1) / (B_coef[0] * B_coef[3] + B_coef[1] + B_coef[2]);

    vector_fp cp_r_ref(m_nsp);
    m_thermo->getCp_R_ref(DATA_PTR(cp_r_ref));
    doublereal cp_ig = std::inner_product(molefracs.begin(), molefracs.end(), cp_r_ref.begin(), 0.0) * GasConstant;
    doublereal cv_ig = cp_ig - GasConstant; // Cv of the ideal gas

//    printf("cv_ig = %f\n", cv_ig);

    doublereal alpha = cv_ig / GasConstant - 3.0 / 2.0;
    doublereal beta = 0.7862 - 0.7109 * omega_M + 1.3168 * pow(omega_M, 2);
    doublereal Zr = 2.0 + 10.5 * pow(T_in / TCrit_M, 2);
    //doublereal q = 0.003586 * sqrt(TCrit_M / (MW_M * 1000.0)) / pow(VCrit_M, 2.0 / 3.0);

//    printf("alpha = %f, beta = %f, Zr = %f\n", alpha, beta, Zr);
    doublereal psi = 1.0 + alpha * ((0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * Zr) / (0.6366 + beta * Zr + 1.061 * alpha * beta));

//    printf("mu = %f, psi = %f, R = %f, MW_M = %f\n", mu, psi, GasConstant, MW_M);

    doublereal mixConductivity_ig = 3.75 * mu * psi * GasConstant / MW_M;
    doublereal mixConductivity = mixConductivity_ig;
    //doublereal mixConductivity = 31.2 * mu * psi / (MW_M) * GasConstant * (1.0 / G2 + B_coef[5] * y) + q * B_coef[6] * pow(y, 2) * sqrt(T_in / TCrit_M) * G2;

    return mixConductivity;
}

doublereal TranscriticalTransport::thermalConductivity_Chung_HP()
{
    //printf("we are inside TranscriticalTransport::thermalConductivity_Chung_HP()\n");
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    doublereal T_in = m_thermo->temperature();
    doublereal rho_in = m_thermo->density();

    doublereal sigma_M    = 0.0;
    doublereal epsOverk_M = 0.0;
    doublereal omega_M    = 0.0;
    doublereal mu_M       = 0.0;
    doublereal kappa_M    = 0.0; // polar correction, only set value to water
    doublereal MW_M       = 0.0;

    for (size_t k = 0; k < nsp; k++) {
        for (size_t l = 0; l < nsp; l++) {
            size_t apos = k * nsp + l;
            doublereal X_X = molefracs[l] * molefracs[k];
            sigma_M += X_X * pow(sigma_IJ[apos], 3.0);
            MW_M += X_X * epsOverk_IJ[apos] * pow(sigma_IJ[apos], 2) * sqrt(MW_IJ[apos]);
            epsOverk_M += X_X * epsOverk_IJ[apos] * pow(sigma_IJ[apos], 3);
            omega_M += X_X * omega_IJ[apos] * pow(sigma_IJ[apos], 3.0);
            kappa_M += X_X * kappa_IJ[apos];
            if (fabs(sigma_IJ[apos]) > 1.0e-20) {
                mu_M += X_X * pow(dipole[k] * dipole[l], 2.0) / pow(sigma_IJ[apos], 3.0);
            }
        }
    }

    sigma_M = pow(sigma_M, 1.0 / 3.0);
    epsOverk_M = epsOverk_M / pow(sigma_M, 3);
    omega_M = omega_M / pow(sigma_M, 3.0);
    mu_M = pow(mu_M * pow(sigma_M, 3.0), 0.25);
    MW_M = pow(MW_M / (epsOverk_M * pow(sigma_M, 2)), 2);

    doublereal VCrit_M = pow((sigma_M / 0.809), 3.0);
    doublereal TCrit_M = 1.2593 * epsOverk_M;
    doublereal TStar_M = T_in / epsOverk_M;
    doublereal mu_RM = 131.3 * mu_M / sqrt(VCrit_M * TCrit_M);
    doublereal Neufeld = 1.16145 * pow(TStar_M, -0.14874) + 0.52487 * exp(-0.7732 * TStar_M) + 2.16178 * exp(-2.43787 * TStar_M); // Eq (9-4-3) in Poling
    doublereal Fc = 1.0 - 0.2756 * omega_M + 0.059035 * pow(mu_RM, 4) + kappa_M;

    // above is repeated as in viscosity_Chung_LP

    doublereal mu = 40.785 * Fc * sqrt(MW_M * T_in) / (Neufeld * pow(VCrit_M, 2.0 / 3.0)) / 10.0e6;
    //doublereal mu = viscosity_Chung_LP();
    //doublereal mu = viscosity_Chung_HP();

    // Chung's conductivity model is used
    // Currently, only the low pressure correction is implemented and verified. High pressure contains the same inconsistency as the chung viscosity. This needs to be addressed.
    const size_t nCoeff = 7;
    static const double a1[nCoeff] = {2.4166, -0.50924, 6.6107, 14.543, 0.79274, -5.8634, 91.089};
    static const double b1[nCoeff] = {0.74824, -1.5094, 5.6207, -8.9139, 0.82019, 12.801, 128.11};
    static const double c1[nCoeff] = {-0.91858, -49.991, 64.760, -5.6379, -0.69369, 9.5893, -54.217};
    static const double d1[nCoeff] = {121.72, 69.983, 27.039, 74.344, 6.3173, 65.529, 523.81};
    // compute B_i coefficients
    vector_fp B_coef(nCoeff);
    for (size_t i = 0; i < nCoeff; i++) {
        B_coef[i] = a1[i] + b1[i] * omega_M + c1[i] * pow(mu_RM, 4.0) + d1[i] * kappa_M;
    }

    doublereal y = rho_in / (MW_M * 1000) * VCrit_M / 6.0;
    //jph the units of density are strange density is mol/cm^3. In poling it is cm^3/mol but I think from Chung 1988 they are typical kg/m^3

    doublereal G1 = (1.0 - 0.5 * y) / pow((1.0 - y), 3);
    doublereal G2 = ((B_coef[0] / y) * (1.0 - exp(-B_coef[3] * y)) + B_coef[1] * G1 * exp(B_coef[4] * y) + B_coef[2] * G1) / (B_coef[0] * B_coef[3] + B_coef[1] + B_coef[2]);

    vector_fp cp_r_ref(m_nsp);
    m_thermo->getCp_R_ref(DATA_PTR(cp_r_ref));
    doublereal cp_ig = std::inner_product(molefracs.begin(), molefracs.end(), cp_r_ref.begin(), 0.0) * GasConstant;
    doublereal cv_ig = cp_ig - GasConstant; // Cv of the ideal gas

    doublereal alpha = cv_ig / GasConstant - 3.0 / 2.0;
    doublereal beta = 0.7862 - 0.7109 * omega_M + 1.3168 * pow(omega_M, 2);
    doublereal Zr = 2.0 + 10.5 * pow(T_in / TCrit_M, 2);
    doublereal q = 0.003586 * sqrt(TCrit_M / (MW_M / 1000.0)) / pow(VCrit_M, 2.0 / 3.0);
    doublereal psi = 1.0 + alpha * ((0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * Zr) / (0.6366 + beta * Zr + 1.061 * alpha * beta));

    //doublereal mixConductivity_ig = 3.75 * mu * psi * GasConstant / MW_M;
    doublereal mixConductivity = 31.2 * mu * psi / (MW_M / 1000) * (1.0 / G2 + B_coef[5] * y) + q * B_coef[6] * pow(y, 2) * sqrt(T_in / TCrit_M) * G2;

    return mixConductivity;
}

//void TranscriticalTransport::getThermalDiffCoeffs(doublereal* const dt)
//{
//    // Method for MultiTransport class:
//    // solveLMatrixEquation();
//    // const doublereal c = 1.6/GasConstant;
//    // for (size_t k = 0; k < m_nsp; k++) {
//    // dt[k] = c * m_mw[k] * m_molefracs[k] * m_a[k];
//    // }
//    throw CanteraError("TranscriticalTransport::getThermalDiffCoeffs",
//                       "Not yet implemented.");
//}

// Pure species critical properties - Tc, Pc, Vc, Zc:
//doublereal TranscriticalTransport::Tcrit_i(size_t i)
//{
//    size_t nsp = m_thermo->nSpecies();
//
//    // Store current molefracs and set temp molefrac of species i to 1.0:
//    vector_fp molefracs = store(i, nsp);
//
//    double tc = m_thermo->critTemperature();
//
//    // Restore actual molefracs:
//    m_thermo->setMoleFractions(&molefracs[0]);
//
//    return tc;
//}

//doublereal TranscriticalTransport::Pcrit_i(size_t i)
//{
//    size_t nsp = m_thermo->nSpecies();
//
//    // Store current molefracs and set temp molefrac of species i to 1.0:
//    vector_fp molefracs = store(i, nsp);
//
//    double pc = m_thermo->critPressure();
//
//    // Restore actual molefracs:
//    m_thermo->setMoleFractions(&molefracs[0]);
//
//    return pc;
//}

//doublereal TranscriticalTransport::Vcrit_i(size_t i)
//{
//    size_t nsp = m_thermo->nSpecies();
//
//    // Store current molefracs and set temp molefrac of species i to 1.0:
//    vector_fp molefracs = store(i, nsp);
//
//    double vc = m_thermo->critVolume();
//
//    // Restore actual molefracs:
//    m_thermo->setMoleFractions(&molefracs[0]);
//
//    return vc;
//}

//doublereal TranscriticalTransport::Zcrit_i(size_t i)
//{
//    size_t nsp = m_thermo->nSpecies();
//
//    // Store current molefracs and set temp molefrac of species i to 1.0:
//    vector_fp molefracs = store(i, nsp);
//
//    double zc = m_thermo->critCompressibility();
//
//    // Restore actual molefracs:
//    m_thermo->setMoleFractions(&molefracs[0]);
//
//    return zc;
//}

//vector_fp TranscriticalTransport::store(size_t i, size_t nsp)
//{
//    vector_fp molefracs(nsp);
//    m_thermo->getMoleFractions(&molefracs[0]);
//
//    vector_fp mf_temp(nsp);
//    for (size_t j = 0; j < nsp; j++) {
//        if (j == i) {
//            mf_temp[j] = 1;
//        } else {
//            mf_temp[j] = 0;
//        }
//    }
//    m_thermo->setMoleFractions(&mf_temp[0]);
//
//    return molefracs;
//}

// Calculates quantum correction term for a species based on Tr and MW, used in
//   viscosity calculation:
doublereal TranscriticalTransport::FQ_i(doublereal Q, doublereal Tr, doublereal MW)
{
    return 1.22 * pow(Q, 0.15) * (1 + 0.00385 * pow(pow(Tr - 12., 2.), 1. / MW) * fabs(Tr - 12) / (Tr - 12));
}

// Set value of parameter values for Takahashi correlation, by interpolating
// table of constants vs. Pr:
doublereal TranscriticalTransport::setPcorr(doublereal Pr, doublereal Tr)
{
    const static double Pr_lookup[17] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0,
        1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0};
    const static double DP_Rt_lookup[17] = {1.01, 1.01, 1.01, 1.01, 1.01, 1.01,
        1.01, 1.02, 1.02, 1.02, 1.02, 1.03, 1.03, 1.04, 1.05, 1.06, 1.07};
    const static double A_ij_lookup[17] = {0.038042, 0.067433, 0.098317,
        0.137610, 0.175081, 0.216376, 0.314051, 0.385736, 0.514553, 0.599184,
        0.557725, 0.593007, 0.696001, 0.790770, 0.502100, 0.837452, 0.890390};
    const static double B_ij_lookup[17] = {1.52267, 2.16794, 2.42910, 2.77605,
        2.98256, 3.11384, 3.50264, 3.07773, 3.54744, 3.61216, 3.41882, 3.18415,
        3.37660, 3.27984, 3.39031, 3.23513, 3.13001};
    const static double C_ij_lookup[17] = {0., 0., 0., 0., 0., 0., 0., 0.141211,
        0.278407, 0.372683, 0.504894, 0.678469, 0.665702, 0., 0.602907, 0., 0.};
    const static double E_ij_lookup[17] = {1., 1., 1., 1., 1., 1., 1., 13.45454,
        14., 10.00900, 8.57519, 10.37483, 11.21674, 1., 6.19043, 1., 1.};

    // Interpolate Pr vs. those used in Takahashi table:
    int Pr_i = 0;
    double frac = 0.;

    if (Pr < 0.1) {
        //frac = (Pr - Pr_lookup[0])/(Pr_lookup[1] - Pr_lookup[0]);
        frac = 0.0;
    } else {
        for (int j = 1; j < 17; j++) {
            if (Pr_lookup[j] > Pr) {
                frac = (Pr - Pr_lookup[j-1])/(Pr_lookup[j] - Pr_lookup[j-1]);
                break;
            }
            Pr_i++;
        }
    }

    // If Pr is greater than the greatest value used by Takahashi (5.0), use the
    //   final table value.  Should eventually add in an extrapolation:
    if (Pr_i == 17) {
        frac = 1.0;
    }

    doublereal f_1 = C_ij_lookup[Pr_i] * pow(Tr, -E_ij_lookup[Pr_i]);
    if (f_1 > 1.0) f_1 = 1.0;
    doublereal P_corr_1 = DP_Rt_lookup[Pr_i] *
                          (1.0 - A_ij_lookup[Pr_i] * pow(Tr, -B_ij_lookup[Pr_i])) *
                          (1 - f_1);
    if (P_corr_1 < Tiny) P_corr_1 = Tiny;

    doublereal f_2 = C_ij_lookup[Pr_i + 1] * pow(Tr, -E_ij_lookup[Pr_i + 1]);
    if (f_2 > 1.0) f_2 = 1.0;
    doublereal P_corr_2 = DP_Rt_lookup[Pr_i + 1] *
                          (1.0 - A_ij_lookup[Pr_i + 1] * pow(Tr, -B_ij_lookup[Pr_i + 1])) *
                          (1 - f_2);
    if (P_corr_2 < Tiny) P_corr_2 = Tiny;

    return P_corr_1 * (1.0 - frac) + P_corr_2 * frac;
}

}
