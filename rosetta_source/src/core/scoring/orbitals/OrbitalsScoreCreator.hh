/*
 * OrbitalsScoreCreator.hh
 *
 *  Created on: Jun 3, 2010
 *      Author: combss
 */

#ifndef ORBITALSSCORECREATOR_HH_
#define ORBITALSSCORECREATOR_HH_

#include <core/scoring/methods/EnergyMethodCreator.hh>

#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace orbitals {

class OrbitalsScoreCreator : public methods::EnergyMethodCreator
{
public:
	/// @brief Instantiate a new OrbitalsScore
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


#endif /* ORBITALSSCORECREATOR_HH_ */
