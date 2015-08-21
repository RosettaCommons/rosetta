#ifndef INCLUDED_core_scoring_vall_lookback_VallLookbackEnergyCreator_hh
#define INCLUDED_core_scoring_vall_lookback_VallLookbackEnergyCreator_hh

#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

class VallLookbackEnergyCreator : public EnergyMethodCreator
{
public:
	/// @brief Instantiate a new vallLookbackScore
	virtual
	methods::EnergyMethodOP
	create_energy_method(
		methods::EnergyMethodOptions const &
	) const;

	/// @brief Return the set of score types claimed by the EnergyMethod
	/// this EnergyMethodCreator creates in its create_energy_method() function
	virtual
	ScoreTypes
	score_types_for_method() const;

};

}
}
}

#endif

