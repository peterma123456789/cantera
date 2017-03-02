/**
 *  @file BlendTransport.h
 *  Interface for class BlendTransport
 */

#ifndef CT_BLENDTRAN_H
#define CT_BLENDTRAN_H

// Cantera includes
#include "GasTransport.h"
#include "cantera/numerics/DenseMatrix.h"
#include "cantera/transport/TranscriticalTransport.h"

namespace Cantera
{

//! Class BlendTransport implements transport properties for
//! high pressure gas mixtures.
/*!
 * The implementation employs a method of corresponding states, using
 *  the Takahashi approach for binary diffusion coefficients, (using
 *  multicomponent averaging rules for the mixture properties, and the
 *  Lucas method for the viscosity of a high-pressure gas mixture.
 *
 * @ingroup tranprops
 */
class BlendTransport : public TranscriticalTransport
{
  protected:

    BlendTransport(thermo_t* thermo = 0) {}

  public:

    virtual int model() const
    {
        return cBlend;
    }

    virtual doublereal viscosity()
    {
        doublereal TcMix = getTcMix();
        doublereal frac = 1.0/(1.0 + exp(-4.0*log(m_thermo->temperature()/TcMix/2.0)));
        //return MultiTransport::viscosity();
        //return TranscriticalTransport::viscosity();
        return frac * MultiTransport::viscosity() + (1.0 - frac) * TranscriticalTransport::viscosity();
    }

    virtual double thermalConductivity()
    {
        doublereal TcMix = getTcMix();
        doublereal frac = 1.0/(1.0 + exp(-4.0*log(m_thermo->temperature()/TcMix/2.0)));
        //printf("TcMix = %f, T = %f, frac = %f\n", TcMix, m_thermo->temperature(), frac);
        //return MultiTransport::thermalConductivity();
        //return TranscriticalTransport::thermalConductivity();
        return frac * MultiTransport::thermalConductivity() + (1.0 - frac) * TranscriticalTransport::thermalConductivity();
    }

    virtual void getMultiDiffCoeffs(const size_t ld, doublereal* const d)
    {
        //MultiTransport::getMultiDiffCoeffs(ld, d);
        TranscriticalTransport::getMultiDiffCoeffs(ld, d);
    }

    friend class TransportFactory;

    virtual void init(ThermoPhase* thermo, int mode = 0, int log_level = 0)
    {
        TranscriticalTransport::init(thermo, HP_Mode, log_level);
    }

  private:

    doublereal getTcMix()
    {
        // update m_molefracs
        m_thermo->getMoleFractions(&m_molefracs[0]);

        size_t nsp = m_thermo->nSpecies();
        doublereal TcMix = 0.0;
        for (size_t k = 0; k < nsp; k++)
            TcMix += m_molefracs[k] * Tcrit[k];
        return 2*TcMix;
    }
};

}

#endif
