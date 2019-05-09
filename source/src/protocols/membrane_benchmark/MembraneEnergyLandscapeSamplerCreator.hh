// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane_benchmark/MembraneEnergyLandscapeSamplerCreator.hh
/// @brief Sample all points on a 2D membrane energy landscape given implicit model and protein dimensions
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_benchmark_MembraneEnergyLandscapeSamplerCreator_HH
#define INCLUDED_protocols_membrane_benchmark_MembraneEnergyLandscapeSamplerCreator_HH

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane_benchmark {

class MembraneEnergyLandscapeSamplerCreator : public protocols::moves::MoverCreator {
public:
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //protocols
} //membrane_benchmark

#endif //INCLUDED_protocols_membrane_benchmark_MembraneEnergyLandscapeSamplerCreator_HH
