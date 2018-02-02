// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/RelaseConstraintFromResidueMoverCreator.hh
/// @brief Remove constraints from the residues specified by a ResidueSelector
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#ifndef INCLUDED_protocols_fold_from_loops_movers_RelaseConstraintFromResidueMoverCreator_hh
#define INCLUDED_protocols_fold_from_loops_movers_RelaseConstraintFromResidueMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace fold_from_loops {
namespace movers {

class ReleaseConstraintFromResidueMoverCreator : public protocols::moves::MoverCreator {
public:
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

}
} //protocols
} //fold_from_loops

#endif //INCLUDED_protocols_fold_from_loops_RelaseConstraintFromResidueMoverCreator_hh
