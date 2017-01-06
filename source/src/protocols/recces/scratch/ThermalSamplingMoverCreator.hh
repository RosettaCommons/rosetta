// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/scratch/ThermalSamplingMoverCreator.hh
/// @brief Use a simulated tempering simulation to refine a pose
/// @author Andy Watkins (amw579@nyu.edu)

#ifndef INCLUDED_protocols_recces_ThermalSamplingMoverCreator_hh
#define INCLUDED_protocols_recces_ThermalSamplingMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace recces {
namespace scratch {

class ThermalSamplingMoverCreator : public protocols::moves::MoverCreator {

public:

	// XRW TEMP  virtual protocols::moves::MoverOP
	// XRW TEMP  create_mover() const;

	// XRW TEMP  virtual std::string
	// XRW TEMP  keyname() const;
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

} //scratch
} //recces
} //protocols

#endif //INCLUDED_protocols/recces_ThermalSamplingMover_fwd_hh
