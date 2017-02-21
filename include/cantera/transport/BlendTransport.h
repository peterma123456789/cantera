/**
 *  @file BlendTransport.h
 *  Interface for class BlendTransport
 */

#ifndef CT_BLENDTRAN_H
#define CT_BLENDTRAN_H

// Cantera includes
#include "GasTransport.h"
#include "cantera/transport/ChungTransport.h"

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
class BlendTransport : public ChungTransport
{
public:
    //! default constructor
    /*!
     *   @param thermo  Optional parameter for the pointer to the ThermoPhase object
     */
    BlendTransport(thermo_t* thermo = 0) {}

    virtual int model() const
    {
        warn_deprecated("BlendTransport::model",
                        "To be removed after Cantera 2.3.");
        return cBlend;
    }

    virtual std::string transportType() const { return "Blend"; }

    virtual doublereal viscosity()
    {
        doublereal TcMix = getTcMix();
        doublereal frac = 1.0/(1.0 + exp(-4.0*log(m_thermo->temperature()/TcMix/2.0)));
        //return MultiTransport::viscosity();
        //return ChungTransport::viscosity();
        return frac * MultiTransport::viscosity() + (1.0 - frac) * ChungTransport::viscosity();
    }

    virtual double thermalConductivity()
    {
        doublereal TcMix = getTcMix();
        doublereal frac = 1.0/(1.0 + exp(-4.0*log(m_thermo->temperature()/TcMix/2.0)));
        //return MultiTransport::thermalConductivity();
        //return ChungTransport::thermalConductivity();
        return frac * MultiTransport::thermalConductivity() + (1.0 - frac) * ChungTransport::thermalConductivity();
    }

    virtual void init(ThermoPhase* thermo, int mode = 0, int log_level = 0)
    {
        ChungTransport::init(thermo, mode, log_level);
    }

private:
    doublereal getTcMix()
    {
        // update m_molefracs
        m_thermo->getMoleFractions(&m_molefracs[0]);

        doublereal TcMix = 0.0;
        for (size_t k = 0; k < m_nsp; k++)
            TcMix += m_molefracs[k] * Tcrit[k];
        return 2 * TcMix;
    }
};

}

#endif
