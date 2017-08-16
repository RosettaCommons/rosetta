// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SwitchResidueTypeSetMover.hh
/// @brief switch between residue type sets (e.g. centroid and all atom)

#ifndef INCLUDED_protocols_simple_moves_SwitchResidueTypeSetMover_hh
#define INCLUDED_protocols_simple_moves_SwitchResidueTypeSetMover_hh

#include <protocols/simple_moves/SwitchResidueTypeSetMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

/// @brief A mover that switches a pose between residue type sets (e.g. centroid and fullatom)
///
/// examples:
///     switch = protocols::simple_moves::SwitchResidueTypeSetMover("centroid")
/// See also:
///     Pose
///     Residue
///     ResidueType
///     ResidueTypeSet
class SwitchResidueTypeSetMover : public moves::Mover {
public:
	SwitchResidueTypeSetMover();
	SwitchResidueTypeSetMover( std::string );

	std::string get_residue_type_set() const;

	/// @brief Applies ResidueTypeSet converion on the pose
	/// @note: a single protocols::moves::Mover only converts in ONE direction e.g. to centroid
	void apply( core::pose::Pose & pose ) override;

	void show(std::ostream & output=std::cout) const override;

	void type_set_tag( std::string const & type_set_tag_in ) { type_set_tag_ = type_set_tag_in; }

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string type_set_tag_;
};

std::ostream &operator<< (std::ostream &os, SwitchResidueTypeSetMover const &mover);

} // simple_moves
} // protocols

#endif
