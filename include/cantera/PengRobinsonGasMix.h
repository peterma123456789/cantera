//! @file PengRobinsonGasMix.h

#ifndef CXX_PENGROBINSONGASMIX
#define CXX_PENGROBINSONGASMIX

#include "thermo/PengRobinsonGasPhase.h"
#include "kinetics/GasKinetics.h"
#include "kinetics/importKinetics.h"
#include "base/stringUtils.h"

namespace Cantera
{

class PengRobinsonGasMix : public PengRobinsonGasPhase, public GasKinetics
{
public:
    PengRobinsonGasMix() : m_ok(false), m_r(0) {}

    PengRobinsonGasMix(const std::string& infile, std::string id_ = "")
        : m_ok(false), m_r(0)
    {
        m_r = get_XML_File(infile);
        m_id = id_;
        if (id_ == "-") {
            id_ = "";
        }
        m_ok = buildSolutionFromXML(*m_r, m_id, "phase", this, this);
        if (!m_ok)
            throw CanteraError("PengRobinsonGasMix",
                               "Cantera::buildSolutionFromXML returned false");
    }

    PengRobinsonGasMix(XML_Node& root, std::string id_)
        : m_ok(false), m_r(&root), m_id(id_)
    {
        m_ok = buildSolutionFromXML(root, id_, "phase", this, this);
    }

    PengRobinsonGasMix(const PengRobinsonGasMix& other)
        : m_ok(false), m_r(other.m_r), m_id(other.m_id)
    {
        m_ok = buildSolutionFromXML(*m_r, m_id, "phase", this, this);
    }

    bool operator!()
    {
        return !m_ok;
    }

    bool ready() const
    {
        return m_ok;
    }

    friend std::ostream& operator<<(std::ostream& s, PengRobinsonGasMix& mix)
    {
        std::string r = mix.report(true);
        s << r;
        return s;
    }

protected:
    bool m_ok;
    XML_Node* m_r;
    std::string m_id;
};

}

#endif
