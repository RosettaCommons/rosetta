// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/movers/AddFlankingVirtualResiduesMover.hh
/// @brief adds virtual residues to either side of a given Pose
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_pose_sewing_movers_AddFlankingVirtualResiduesMover_HH
#define INCLUDED_protocols_pose_sewing_movers_AddFlankingVirtualResiduesMover_HH

// Unit headers
#include <protocols/pose_sewing/movers/AddFlankingVirtualResiduesMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace pose_sewing {
namespace movers {

///@brief adds virtual residues to either side of a given Pose
class AddFlankingVirtualResiduesMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	AddFlankingVirtualResiduesMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	AddFlankingVirtualResiduesMover( AddFlankingVirtualResiduesMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~AddFlankingVirtualResiduesMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;
	void
	add_flanking_virtual_residues( core::pose::Pose & pose );
	void
	show( std::ostream & output = std::cout ) const override;

	core::Size
	get_N_term_length() const;

	core::Size
	get_C_term_length() const;

	core::select::residue_selector::ResidueSelectorCOP
	get_vital_selector() const;

	core::Size
	get_chain_to_modify() const;

	bool
	get_remove_pre_pose() const;


	void
	set_N_term_length(core::Size);

	void
	set_C_term_length(core::Size);

	void
	set_vital_selector(core::select::residue_selector::ResidueSelectorCOP);

	void
	set_chain_to_modify(core::Size);

	void
	set_remove_pre_pose(bool);

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	//AddFlankingVirtualResiduesMover & operator=( AddFlankingVirtualResiduesMover const & src );

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

private: // methods

private: // data
	core::Size N_term_length_ = 0;
	core::Size C_term_length_ = 0;
	core::Size chain_to_modify_ = 1;
	core::select::residue_selector::ResidueSelectorCOP vital_selector_;
	bool remove_pre_pose_ = false;
};

std::ostream &
operator<<( std::ostream & os, AddFlankingVirtualResiduesMover const & mover );

} //protocols
} //pose_sewing
} //movers

#endif //protocols_pose_sewing_movers_AddFlankingVirtualResiduesMover_HH
