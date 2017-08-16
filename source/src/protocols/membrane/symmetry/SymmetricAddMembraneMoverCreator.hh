// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/symmetry/SymmetricAddMembraneMoverCreator.hh
///
/// @brief      Add Membrane Representation to a Symmetric starting pose
/// @details    Given a symmetrized pose, add the membrane components,
///             ensuring all data descriptions are for the asymmetric subunit
///             (in terms of scoring) and the membrane is the root of the
///             whole system. After applying SymmetricAddMembraneMover
///             pose.conformaiton().is_membrane() AND is_symmetric( pose )
///             should both return true
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (2/9/15)


#ifndef INCLUDED_protocols_membrane_SymmetricAddMembraneMoverCreator_hh
#define INCLUDED_protocols_membrane_SymmetricAddMembraneMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {
namespace symmetry {

/// @brief Mover Creator
class SymmetricAddMembraneMoverCreator : public protocols::moves::MoverCreator {

public:

	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

} // symmetry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_SymmetricAddMembraneMoverCreator_hh
