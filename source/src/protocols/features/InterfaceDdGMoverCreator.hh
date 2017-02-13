// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/InterfaceDdGMoverCreator.hh
/// @brief Calculates ddG of binding
/// @author Kyle Barlow (kb@kylebarlow.com)

#ifndef INCLUDED_protocols_features_InterfaceDdGMoverCreator_hh
#define INCLUDED_protocols_features_InterfaceDdGMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace features {

class InterfaceDdGMoverCreator : public protocols::moves::MoverCreator {
public:
	virtual protocols::moves::MoverOP create_mover() const override;
	virtual std::string keyname() const override;
	static std::string mover_name();
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //protocols
} //features

#endif //INCLUDED_protocols_features_InterfaceDdGMoverCreator_fwd_hh
