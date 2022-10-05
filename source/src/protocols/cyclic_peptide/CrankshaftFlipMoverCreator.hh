// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/CrankshaftFlipMoverCreator.hh
/// @brief This mover performs a cyclic-aware crankshaft flip at a residue position.
/// @author P. Douglas Renfrew (doug.renfrew@gmail.com)
/// @author James Eastwood (jre318@nyu.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_CrankshaftFlipMoverCreator_HH
#define INCLUDED_protocols_cyclic_peptide_CrankshaftFlipMoverCreator_HH

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace cyclic_peptide {

class CrankshaftFlipMoverCreator : public protocols::moves::MoverCreator {
public:
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} //cyclic_peptide
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_CrankshaftFlipMoverCreator_HH
