// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/scoring/MEnvElectroAtomParams.hh
/// @brief A container for storing implicit membrane environment electrostatic parameters and derivatives
/// @author Rituparna Samanta (rsamant2@jhu.edu)

#ifndef INCLUDED_protocols_membrane_scoring_MEnvElectroAtomParams_hh
#define INCLUDED_protocols_membrane_scoring_MEnvElectroAtomParams_hh

#include <protocols/membrane/scoring/MEnvElectroAtomParams.fwd.hh>

// Basic headers
#include <core/types.hh>
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/VirtualBase.hh>

namespace protocols {
namespace membrane {
namespace scoring {

/// @brief A container for storing implicit membrane environment electrostatic parameters and derivatives
class MEnvElectroAtomParams : public utility::VirtualBase {

public:

	MEnvElectroAtomParams();

	MEnvElectroAtomParams(
		std::string const & atom_type_name,
		core::Real const charge,
		core::Real const elec_field,
		core::Real const elec_field_gradient,
		core::Vector const f1,
		core::Vector const f2
	);


	~MEnvElectroAtomParams() override;

	MEnvElectroAtomParamsOP
	clone() const;

	std::string atom_type_name() const { return atom_type_name_; }
	core::Real charge() const { return charge_; }
	core::Real elec_field() const { return elec_field_; }
	core::Real elec_field_gradient() const { return elec_field_gradient_; }
	core::Vector f1() const { return f1_; }
	core::Vector f2() const { return f2_; }

private:

	std::string atom_type_name_="";
	core::Real charge_=0.0;
	core::Real elec_field_=0.0;
	core::Real elec_field_gradient_=0.0;
	numeric::xyzVector< core::Real > f1_{0.0};
	numeric::xyzVector< core::Real > f2_{0.0};

};

} //protocols
} //membrane
} //scoring

#endif //INCLUDED_protocols_membrane_scoring_MEnvElectroAtomParams_hh

