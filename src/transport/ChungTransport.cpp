/**
 *  @file ChungTransport.cpp
 *  Implementation file for class ChungTransport
 *
 *  Transport parameters are calculated using corresponding states models:
 *      Binary diffusion coefficients use the generalized chart described by
 *      Takahashi, et al. and viscosity calcualtions use the Lucas method.
 *      All methods are described in Reid, Prausnitz, and Polling, "The Properties
 *      of Gases and Liquids, 4th ed., 1987 (viscosity in Ch. 9, Thermal
 *      conductivity in Ch. 10, and Diffusion coefficients in Ch. 11).
 *
 **/
#include "cantera/transport/ChungTransport.h"
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

ChungTransport::ChungTransport(thermo_t* thermo) : MultiTransport(thermo) {}

void ChungTransport::init(ThermoPhase* thermo, int mode, int log_level)
{
    MultiTransport::init(thermo, mode, log_level);

    // get critical properties from thermo object, this is quick hack
    // ideally should use GasTranportData object to set omega and dipole
    IsCrit.resize(m_nsp);
    Tcrit.resize(m_nsp);
    Pcrit.resize(m_nsp);
    Vcrit.resize(m_nsp);
    Zcrit.resize(m_nsp);
    omega.resize(m_nsp);
    dipole.resize(m_nsp);
    kappa.resize(m_nsp);
    // TODO: need to push_back zero

    ReadCriticalProperties();

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

void ChungTransport::ReadCriticalProperties()
{
    //! \brief Reads critical properties based on JP's papaer
    for (size_t k = 0; k < m_nsp; k++) {
        if (m_thermo->speciesName(k) == "H2") {
            IsCrit[k] = 1;
            Tcrit[k] = 33.145;     // K
            Pcrit[k] = 1.2964e+6;  // Pa
            Vcrit[k] = 64.4834e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = -0.219;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "O2") {
            IsCrit[k] = 1;
            Tcrit[k] = 154.5800;   // K
            Pcrit[k] = 5.0430e+6;  // Pa
            Vcrit[k] = 73.3682e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0222;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "H2O") {
            IsCrit[k] = 1;
            Tcrit[k] = 647.10;    // K
            Pcrit[k] = 22.064e+6; // Pa
            Vcrit[k] = 55.95e-3;  // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.3443;
            dipole[k] = 1.855;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "O") {
            IsCrit[k] = 0;
            Tcrit[k] = 105.28;   // K
            Pcrit[k] = 7.088e+6; // Pa
            Vcrit[k] = 41.21e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "H") {
            IsCrit[k] = 0;
            Tcrit[k] = 190.82;    // K
            Pcrit[k] = 31.013e+6; // Pa
            Vcrit[k] = 17.07e-3;  // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "OH") {
            IsCrit[k] = 0;
            Tcrit[k] = 105.28;   // K
            Pcrit[k] = 7.088e+6; // Pa
            Vcrit[k] = 41.21e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "H2O2") {
            IsCrit[k] = 0;
            Tcrit[k] = 141.34;   // K
            Pcrit[k] = 4.786e+6; // Pa
            Vcrit[k] = 81.93e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "HO2") {
            IsCrit[k] = 0;
            Tcrit[k] = 141.34;   // K
            Pcrit[k] = 4.786e+6; // Pa
            Vcrit[k] = 81.93e-3; // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else if (m_thermo->speciesName(k) == "N2") {
            IsCrit[k] = 0;
            Tcrit[k] = 126.19;    // K
            Pcrit[k] = 3.3958e+6; // Pa
            Vcrit[k] = 89.41e-3;  // m3/kmol
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0372;
            dipole[k] = 0.0;
            kappa[k] = 0.0;
        } else {
            cout << " Unknown or non-major species : " << m_thermo->speciesName(k)
                 << ". All critical properties were set to zero." << endl;
        }
    }
}

void ChungTransport::updateThermal_T()
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
                m_bdiff(i, j) *= max(P_corr_ij, 0.4);
            }
        }
    }
}

doublereal ChungTransport::viscosity()
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

doublereal ChungTransport::thermalConductivity()
{
    //printf("we are inside ChungTransport::thermalConductivity_Chung_HP()\n");
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
    m_thermo->getCp_R_ref(&cp_r_ref[0]);
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

doublereal ChungTransport::setPcorr(doublereal Pr, doublereal Tr)
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
