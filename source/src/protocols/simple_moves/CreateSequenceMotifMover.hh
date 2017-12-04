// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/CreateSequenceMotifMover.hh
/// @brief Create a sequence motif in a region of protein using the SequenceMotifTaskOperation.  Uses psueo-regular expressions to define the motif.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_CreateSequenceMotifMover_HH
#define INCLUDED_protocols_simple_moves_CreateSequenceMotifMover_HH

// Unit headers
#include <protocols/simple_moves/CreateSequenceMotifMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/toolbox/task_operations/SequenceMotifTaskOperation.fwd.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <core/pack/task/operation/TaskOperations.fwd.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>





//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace simple_moves {

///@brief Simple mover to Create a sequence motif in a region of protein using the SequenceMotifTaskOperation.  Uses psueo-regular expressions to define the motif.
///
///@details
/// Simply calls the packer using the operation, optionally packing neighbor residues as it does so.
/// If you need something more complex, use the SequenceMotifTaskOperation directly.
///
///MOTIF:
///
///  This is slightly similar to a regex, but not quite. We are not matching a sequence, we are designing in a motif regardless of the current sequence, anywhere in a protein.
///
///   - Each letter corresponds to a position. Using [ ] indicates a more complicated expression for that position.
///   - An X indicates it can be anything, and that we are designing here.
///   - An AA Letter, like V, indicates that that position will be designed to a V.
///   - A - charactor indicates that that position stays with whatever it is currently.  We essentially skip this position.
///   - An expression like: [^PAV] indicates that we will design anything except Proline, Alanine, and Valine
///   - An expression like: [NTS] indicates that that position can be Asparigine, Threonine, or Serine and only of these will be enabled during the design.
///   - RESFILE commands are accepted as well. FOr example [POLAR] is totally cool.  Separate these by commas.
///
/// EXAMPLE:
///  Glycosylation N-Linked motif design: N[^P][ST]
///
class CreateSequenceMotifMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	CreateSequenceMotifMover();

	CreateSequenceMotifMover( core::select::residue_selector::ResidueSelectorCOP selector );

	/// @brief Copy constructor (not needed unless you need deep copies)
	CreateSequenceMotifMover( CreateSequenceMotifMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~CreateSequenceMotifMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

public:

	///@brief Set a string that tells this operation how to design.
	///@details
	///  This is slightly similar to a regex, but not quite. We are not matching a sequence, we are designing in a motif regardless of the current sequence, anywhere in a protein.
	///
	///   - Each letter corresponds to a position. Using [ ] indicates a more complicated expression for that position.
	///   - An X indicates it can be anything, and that we are designing here.
	///   - An AA Letter, like V, indicates that that position will be designed to a V.
	///   - A - charactor indicates that that position stays with whatever it is currently.  We essentially skip this position.
	///   - An expression like: [^PAV] indicates that we will design anything except Proline, Alanine, and Valine
	///   - An expression like: [NTS] indicates that that position can be Asparigine, Threonine, or Serine and only of these will be enabled during the design.
	///   - RESFILE commands are accepted as well. A % is required in front of them.
	///     For example [%POLAR] is totally cool. These are the same as in a resfile line
	///     So, Non-cannonicals work the same:
	///       [%EMPTY NC R2 NC T6 NC OP5]
	///
	/// EXAMPLE:
	///  Glycosylation N-Linked motif design: N[^P][ST]
	///
	void
	set_motif( std::string motif );

	///@brief Set a residue selector where each position returned from the selector is a place in which we create the motif.
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	///@brief Should we pack neighbors?  Default True.
	void
	set_pack_neighbors( bool pack_neighbors );

	///@brief Should we design neighbors?  Default False.
	void
	set_design_neighbors( bool design_neighbors );

	///@brief Set the distance for any neighbor packing set.
	/// Default is 6A
	void
	set_neighbor_distance( core::Real neighbor_distance );

	void
	set_score_function( core::scoring::ScoreFunctionCOP scorefxn );

	///@brief Set the number of packing/design rounds.  Default is 5!
	void
	set_pack_rounds( core::Size pack_rounds );

public:
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

	//CreateSequenceMotifMover & operator=( CreateSequenceMotifMover const & src );

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

	void
	initialize_objects(core::pose::Pose & pose);

private: // data

	core::select::residue_selector::ResidueSelectorCOP selector_ = nullptr;
	std::string motif_;
	bool pack_neighbors_ = true;
	bool design_neighbors_ = false;
	core::Real pack_distance_ = 6.0;

	core::pack::task::TaskFactoryOP tf_ = nullptr;
	utility::vector1< core::pack::task::operation::OperateOnResidueSubsetOP > neighbor_operations_;
	core::pack::task::operation::InitializeFromCommandlineOP cmd_line_operation_ = nullptr;

	protocols::toolbox::task_operations::SequenceMotifTaskOperationOP motif_operation_ = nullptr;
	core::scoring::ScoreFunctionCOP scorefxn_ = nullptr;

	core::Size pack_rounds_ = 5;

};

std::ostream &
operator<<( std::ostream & os, CreateSequenceMotifMover const & mover );

} //protocols
} //simple_moves

#endif //protocols_simple_moves_CreateSequenceMotifMover_HH
