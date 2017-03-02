/**
 *  @file TranscriticalTransport.h
 *  Interface for class TranscriticalTransport
 */

#ifndef CT_TRANSCRITICALTRAN_H
#define CT_TRANSCRITICALTRAN_H

// Cantera includes
#include "GasTransport.h"
#include "cantera/numerics/DenseMatrix.h"
#include "cantera/transport/MultiTransport.h"

namespace Cantera
{

//! Class MultiTransport implements transport properties for
//! high pressure gas mixtures.
/*!
 * The implementation employs a method of corresponding states, using
 *  the Takahashi approach for binary diffusion coefficients, (using
 *  multicomponent averaging rules for the mixture properties, and the
 *  Lucas method for the viscosity of a high-pressure gas mixture.
 *
 * @ingroup tranprops
 */
class TranscriticalTransport : public MultiTransport
{
  protected:

    //! default constructor
    /*!
     *   @param thermo  Optional parameter for the pointer to the ThermoPhase object
     */
    TranscriticalTransport(thermo_t* thermo=0);

  public:

    virtual int model() const
    {
        if (m_mode == LP_Mode) {
            return cTrans_LP;
        } else if (m_mode == HP_Mode) {
            return cTrans_HP;
        } else {
            throw CanteraError("TranscriticalTransport::model",
                               "Unknown mode");
        }
    }

    virtual doublereal viscosity()
    {
        //return viscosity_Lucas_HP();
        //return MultiTransport::viscosity();
        if (m_mode == LP_Mode) {
            return viscosity_Chung_LP();
        } else if (m_mode == HP_Mode) {
            return viscosity_Chung_HP();
         } else {
            throw CanteraError("TranscriticalTransport::viscosity",
                               "Unknown mode");
        }
    }

    //! Return the thermal diffusion coefficients (kg/m/s)
    /*!
     *  Currently not implemented for this model
     */
    //virtual void getThermalDiffCoeffs(doublereal* const dt);

    virtual double thermalConductivity()
    {
        //return thermalConductivity_ElyHanley_HP();
        //return MultiTransport::thermalConductivity();
        if (m_mode == LP_Mode) {
            return thermalConductivity_Chung_LP();
        } else if (m_mode == HP_Mode) {
            return thermalConductivity_Chung_HP();
        } else {
            throw CanteraError("TranscriticalTransport::thermalConductivity",
                               "Unknown mode");
        }
    }

    /*! Returns the matrix of binary diffusion coefficients
     *
     *      d[ld*j +  i] = rp*m_bdiff(i,j)*(DP)_R;
     *
     * @param ld    offset of rows in the storage
     * @param d     output vector of diffusion coefficients.  Units of m**2 / s
     */
    //virtual void getBinaryDiffCoeffs(const size_t ld, doublereal* const d);

    virtual void getMultiDiffCoeffs(const size_t ld, doublereal* const d)
    {
        //MultiTransport::getMultiDiffCoeffs(ld, d);
        if (m_mode == LP_Mode) {
            MultiTransport::getMultiDiffCoeffs(ld, d);
        } else if (m_mode == HP_Mode) {
            MultiTransport::getMultiDiffCoeffs(ld, d);
            //getMultiDiffCoeffs_Takahashi(ld, d);
        } else {
            throw CanteraError("TranscriticalTransport::getMultiDiffCoeffs",
                               "Unknown mode");
        }
    }

    virtual void getMixDiffCoeffs(doublereal* const d)
    {
        if (m_mode == LP_Mode) {
            MultiTransport::getMixDiffCoeffs(d);
        } else if (m_mode == HP_Mode) {
            getMixDiffCoeffs_Takahashi(d);
        } else {
            throw CanteraError("TranscriticalTransport::getMultiDiffCoeffs",
                               "Unknown mode");
        }
    }

//    virtual void updateDiff_T();

    virtual void updateThermal_T();

    friend class TransportFactory;

    virtual void init(ThermoPhase* thermo, int mode = 0, int log_level = 0);

    void ReadCriticalProperties() const;

  protected:

    //virtual doublereal Tcrit_i(size_t i);

    //virtual doublereal Pcrit_i(size_t i);

    //virtual doublereal Vcrit_i(size_t i);

    //virtual doublereal Zcrit_i(size_t i);

    //vector_fp store(size_t i, size_t nsp);

    //virtual doublereal CT_i(doublereal T_0);

    virtual doublereal FQ_i(doublereal Q, doublereal Tr, doublereal MW);

    virtual doublereal setPcorr(doublereal Pr, doublereal Tr);

    // critical properties
    mutable vector_fp Tcrit;
    mutable vector_fp Pcrit;
    mutable vector_fp Vcrit;
    mutable vector_fp Zcrit;
    mutable vector_fp omega;
    mutable vector_fp dipole;
    mutable vector_fp kappa;
    mutable vector_fp sigma_IJ;
    mutable vector_fp epsOverk_IJ;
    mutable vector_fp omega_IJ;
    mutable vector_fp MW_IJ;
    mutable vector_fp kappa_IJ;

  private:

    doublereal viscosity_Lucas_HP();
    doublereal viscosity_Chung_LP();
    doublereal viscosity_Chung_HP();

    doublereal thermalConductivity_ElyHanley_HP();
    doublereal thermalConductivity_Chung_LP();
    doublereal thermalConductivity_Chung_HP();

    void getMultiDiffCoeffs_Takahashi(const size_t ld, doublereal* const d);
    void getMixDiffCoeffs_Takahashi(doublereal* const d);
};

}

#endif
