// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ncbb/ResidueReplacementRebuildMover.hh
/// @brief A simple method to mutate a pose fairly destructively, and barely recover
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_protocols_ncbb_ResidueReplacementRebuildMover_HH
#define INCLUDED_protocols_ncbb_ResidueReplacementRebuildMover_HH

// Unit headers
#include <protocols/ncbb/ResidueReplacementRebuildMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace ncbb {

///@brief A simple method to mutate a pose fairly destructively, and barely recover
class ResidueReplacementRebuildMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	ResidueReplacementRebuildMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	ResidueReplacementRebuildMover( ResidueReplacementRebuildMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~ResidueReplacementRebuildMover() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	core::pose::Pose
	make_new_pose(
		core::pose::Pose const & pose
	) const;

	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;


public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	//ResidueReplacementRebuildMover & operator=( ResidueReplacementRebuildMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void set_resi_chosen( core::Size const setting ) { resi_chosen_ = setting; }
	void set_resi_len( core::Size const setting ) { resi_len_ = setting; }
	void set_resname( std::string const & setting ) { resname_ = setting; }
	void score_function( core::scoring::ScoreFunctionCOP sfxn ) { scorefxn_ = sfxn; }

	core::Size set_resi_chosen() const { return resi_chosen_; }
	core::Size set_resi_len() const { return resi_len_; }
	std::string set_resname() const { return resname_; }
	core::scoring::ScoreFunctionCOP score_function() const { return scorefxn_; }

private: // methods

private: // data

	Size resi_chosen_ = 0;
	Size resi_len_ = 0;
	std::string resname_ = "ALA";
	core::scoring::ScoreFunctionCOP scorefxn_;

};

std::ostream &
operator<<( std::ostream & os, ResidueReplacementRebuildMover const & mover );

} //protocols
} //ncbb

#endif //protocols_ncbb_ResidueReplacementRebuildMover_HH
