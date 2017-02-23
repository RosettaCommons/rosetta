// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/LocalBackboneMover.hh
/// @brief LocalBackboneMover moves a stretch of backbone locally.
/// @author xingjiepan (xingjiepan@gmail.com)

#ifndef INCLUDED_protocols_backbone_moves_local_backbone_mover_LocalBackboneMover_hh
#define INCLUDED_protocols_backbone_moves_local_backbone_mover_LocalBackboneMover_hh

// Unit headers
#include <protocols/backbone_moves/local_backbone_mover/types.hh>
#include <protocols/backbone_moves/local_backbone_mover/FreePeptide.fwd.hh>
#include <protocols/backbone_moves/local_backbone_mover/GapCloser.fwd.hh>
#include <protocols/backbone_moves/local_backbone_mover/LocalBackboneMover.fwd.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/FreePeptideMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {

///@brief LocalBackboneMover moves a stretch of backbone locally.
class LocalBackboneMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	LocalBackboneMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	LocalBackboneMover( LocalBackboneMover const &) = default;

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~LocalBackboneMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

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

	//LocalBackboneMover & operator=( LocalBackboneMover const & src );

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

	////////////////////////
	/// Specific Methods ///
	////////////////////////

public:

	/// @brief Clear all free peptide movers
	void clear_free_peptide_movers(){
		free_peptide_movers_.clear();
	}

	/// @brief Add a free peptide mover
	void add_free_peptide_mover(free_peptide_movers::FreePeptideMoverOP free_peptide_mover){
		free_peptide_movers_.push_back(free_peptide_mover);
	}

	/// @brief Set the free peptide
	void set_free_peptide(core::pose::Pose &pose, Size pivot1, Size pivot2);

	/// @brief Set the maximum number of move trials
	void set_max_trial_num(Size max_trial_num){
		max_trial_num_ = max_trial_num;
	}

private: // methods

private: // data

	GapCloserOP gap_closer_ = nullptr;
	vector1 <free_peptide_movers::FreePeptideMoverOP> free_peptide_movers_;

	FreePeptideOP free_peptide_;

	Size max_trial_num_ = 0;

	Size pivot1_ = 0;
	Size pivot2_ = 0;
};

std::ostream &
operator<<( std::ostream & os, LocalBackboneMover const & mover );

} //protocols
} //backbone_moves
} //local_backbone_mover

#endif //protocols/backbone_moves/local_backbone_mover_LocalBackboneMover_hh
