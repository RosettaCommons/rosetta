// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/IdealizeAnomericHydrogens.hh
/// @brief This code changes all of the dihedrals of a particular glycosidic linkage based on database info,
///   esentially sampling carbohydrate dihedral conformers of two residues.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_protocols_carbohydrates_IdealizeAnomericHydrogens_hh
#define INCLUDED_protocols_carbohydrates_IdealizeAnomericHydrogens_hh

// Unit headers
#include <protocols/carbohydrates/IdealizeAnomericHydrogens.fwd.hh>
#include <protocols/simple_moves/BBDihedralSamplerMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/chemical/carbohydrates/LinkageConformers.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/id/types.hh>

#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>


namespace protocols {
namespace carbohydrates {

///@brief This code sets all the anomeric hydrogen positions based on the input structure
class IdealizeAnomericHydrogens : public protocols::moves::Mover {

public:

	///@brief Default constructor
	IdealizeAnomericHydrogens() {};

	// copy constructor
	//IdealizeAnomericHydrogens( IdealizeAnomericHydrogens const & src ){};

	~IdealizeAnomericHydrogens() override;

	void
	apply( core::pose::Pose & pose ) override;


public:
	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	//void parse_my_tag(
	// utility::tag::TagCOP tag,
	// basic::datacache::DataMap & data,
	// protocols::filters::Filters_map const & filters,
	// protocols::moves::Movers_map const & movers,
	// core::pose::Pose const & pose ) override;

	/// @brief required in the context of the parser/scripting scheme
	moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;


private:

};

} //carbohydrates
} //protocols

#endif  // protocols_carbohydrates_IdealizeAnomericHydrogens_hh
