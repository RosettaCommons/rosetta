// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/LinkageConformerMoverCreator.hh
/// @brief This code changes all of the dihedrals of a particular glycosidic linkage based on database info, esentially sampling carbohydrate dihedral conformers of two residues.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_carbohydrates_LinkageConformerMoverCreator_hh
#define INCLUDED_protocols_carbohydrates_LinkageConformerMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>


namespace protocols {
namespace carbohydrates {

class LinkageConformerMoverCreator : public protocols::moves::MoverCreator {

public:

	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;



};

} //protocols
} //carbohydrates

#endif //INCLUDED_protocols/carbohydrates_LinkageConformerMover_fwd_hh



