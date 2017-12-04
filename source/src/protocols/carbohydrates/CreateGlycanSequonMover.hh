// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/CreateGlycanSequonMover.hh
/// @brief Mutates residues to create a potential glycosylation site using known sequence motifs of N- or C- linked glycans.  Includes options for Enhanced Sequons for N-linked glycans that have been shown to have higher rates of glycosylation as well as other positions that have been shown to influence the glycosylation chemistry.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_carbohydrates_CreateGlycanSequonMover_HH
#define INCLUDED_protocols_carbohydrates_CreateGlycanSequonMover_HH

// Unit headers
#include <protocols/carbohydrates/CreateGlycanSequonMover.fwd.hh>
#include <protocols/simple_moves/CreateSequenceMotifMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace carbohydrates {


enum GlycanSequon{
	n_linked_typical = 1,
	n_linked_basic_enhanced,
	n_linked_best_enhanced,
	c_linked_NxC,
	c_linked_WxxW,
	c_linked_WSTxC

};

/// @brief Mutates residues to create a potential glycosylation site using known sequence motifs of N- or C- linked glycans.
///  Includes options for Enhanced Sequons for N-linked glycans that have been shown to have higher rates of glycosylation
///   as well as other positions that have been shown to influence the glycosylation chemistry.
///
/// @details
///  Creates the glycan sequence motif around (and including) the potential glycosylation site.
///  If the site could not be created due to the position being too close to the beginning or end of the sequence,
///  Will set the mover status to fail, do not retry.
///
///  Creates an N-Linked Sequence Motif by default using the non-enhanced motif.  This can be changed in options.
class CreateGlycanSequonMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	CreateGlycanSequonMover();

	CreateGlycanSequonMover( core::select::residue_selector::ResidueSelectorCOP selector );

	/// @brief Copy constructor (not needed unless you need deep copies)
	CreateGlycanSequonMover( CreateGlycanSequonMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~CreateGlycanSequonMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	///
	/// @details
	///  If the sequon could not be created (for example the length of the sequon is longer than the protein/etc.
	///  It will set mover_failed_do_not retry.  It is your responsiblity to catch this and do whatever you wish with it.
	///
	void
	apply( core::pose::Pose & pose ) override;

public:


	///@brief Set the target glycosylation position for which the sequence motif will be created around.
	/// If this is in beginning or end of protein and could not be created, will set fail do not retry mover status.
	///
	void
	set_glycosylation_position( core::Size position, core::pose::Pose const & pose );

	///@brief Set a number of positions using a residue selector.
	///  The positions are where glycosylation is intended to occur.
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	///@brief Instead of having each position be the glycosylation position, have each position be the start of the sequon.
	///  This just makes using ResidueSelectors with this class a bit easier if needed.
	///
	/// This only effects the basic_enhanced and best_enhanced sequons!
	void
	set_positions_as_start_of_sequon( bool positions_as_start);

	///@brief Set to use the enhanced sequon which has been shown to result in higher rates of glycosylation and slightly lower complex types as the glycan site. (Aromatics at n+2)
	///
	/// DEFAULT: False
	///
	///REFS:
	///  "Enhanced Aromatic Sequons Increase Oligosaccharyltransferase Glycosylation Efficiency and Glycan Homogeneity"
	///   Murray et al., 2015, Chemistry & Biology 22, 1052–1062 http://dx.doi.org/10.1016/j.chembiol.2015.06.017
	///
	///  "Residues Comprising the Enhanced Aromatic Sequon Influence Protein N‐Glycosylation Efficiency"
	///   Yen-Wen Huang, J. Am. Chem. Soc. 2017, 139, 12947−12955 DOI: 10.1021/jacs.7b03868
	///
	void
	set_use_basic_enhanced_n_linked_sequon( bool enhanced);

	///@brief Set any sequon type
	///
	void
	set_glycan_sequon_type( GlycanSequon sequon );


	void
	set_score_function( core::scoring::ScoreFunctionCOP scorefxn );

	///@brief If there is an X in the motif, should we design it or leave it alone.
	/// Default is to leave it alone.
	void
	set_design_x_positions( bool const design_x_positions );

public:

	///@brief Number of rounds to run packing/design.  Default is 5
	void
	set_pack_rounds( core::Size pack_rounds );

	///@brief Should we pack neighbors of the motif?  Default is True
	void
	set_pack_neighbors( bool const pack_neighbors );

	///@brief Should we design neighbors of the motif? Default is False.
	void
	set_design_neighbors( bool const design_neighbors );

	///@brief Set the neighbor detection distance for any packing/design of neighbors.
	void
	set_pack_distance( core::Real const pack_distance );

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

	//CreateGlycanSequonMover & operator=( CreateGlycanSequonMover const & src );

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
	std::map< GlycanSequon, std::string > sequons_;
	core::select::residue_selector::ResidueSelectorCOP selector_ = nullptr;
	GlycanSequon sequon_type_ = n_linked_typical;
	bool positions_as_start_of_sequon_ = false;
	bool design_x_positions_;

	simple_moves::CreateSequenceMotifMoverOP motif_mover_ ;
	bool pack_neighbors_ = true;
	bool design_neighbors_ = false;
	core::Real pack_distance_ = 6.0;
	core::Size pack_rounds_ = 5;

	core::scoring::ScoreFunctionCOP scorefxn_ = nullptr;

};


///@brief Create a map of the name and full sequon
std::map< GlycanSequon, std::string >
create_sequons();

std::ostream &
operator<<( std::ostream & os, CreateGlycanSequonMover const & mover );

} //protocols
} //carbohydrates

#endif //protocols_carbohydrates_CreateGlycanSequonMover_HH
