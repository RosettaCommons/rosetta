// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanTreeRelax.hh
/// @brief A protocol for optimizing glycan trees using GlycanRelax from the base of the tree out to the leaves.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Sebastian Raemisch (raemisch@scripps.edu)

#ifndef INCLUDED_protocols_carbohydrates_GlycanTreeRelax_HH
#define INCLUDED_protocols_carbohydrates_GlycanTreeRelax_HH

// Unit headers
#include <protocols/carbohydrates/GlycanTreeRelax.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace carbohydrates {

enum QuenchDirection {
	forward = 1,
	backward
};

///@brief A protocol for optimizing glycan trees using GlycanRelax from the base of the tree out to the leaves.
///@details Works by making all other residues virtual except the ones it is working on (current Layer).
/// A virtual residue is not scored.
/// It will start at the first glycan residues, and then move out to the edges.
///
/// GENERAL ALGORITHM
///
/// We start at the roots, and make all other glycan residues virtual.
/// We first model towards the leaves and this is considered the forward direction.
/// GlycanRelax is used for the actual modeling, we only model a layer at a time, until we reach the tips.
/// If more than one round is set, the protocol will move backwards on the next round, from the leafs to the roots.
/// A third round will involve relaxation again in the forward direction.
/// So we go forward, back, forward, etc. for how ever many rounds you set.
///
/// QUECHING
///
/// By default, we model all glycans simultaneously. First, all glycan roots (the start of the tree), and slowly unvirtualize
/// all glycan residues, while only modeling each layer.
/// Alternatively, we can choose a particular glycan tree, run the algorithm, and then choose another glycan tree randomly until all
/// glycan trees have been optimized.
/// Here, we call this quenching.
///
/// GLYCAN LAYERS
///
/// Draw a tree on a paper.  We start with the beginning N residues, and work our way out towards the leaves.
/// Layers are defined by the glycan residue distance to the rooot.  This enables branching residues to be considered the same
/// layer conceptually and computationally, and allows them to be modeled together.
///
/// --LAYER SIZE--
///
///  The distance that make up a layer.  If we have a distance of 2,
///  we first model all glycans that are equal to or less than 2 residue distance to the root.
///  We then slide this layer up.  So we take all residues that have a distance between 3 and 1, and so on.
///
///  --WINDOW SIZE--
///
///  The layers are slid down throught the tree of the glycan.  The window size represents the overlap in the layers.
///  A window size of 1, means that the last residue (or residues of layer 1) from the last modeling effort, will be used again as
///  part of the next layer.  A window size of 0, means that no residues will be re-modeled.
///  Typically, we would want at least a window size of 1.
///
class GlycanTreeRelax : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	GlycanTreeRelax();

	/// @brief Copy constructor (not needed unless you need deep copies)
	GlycanTreeRelax( GlycanTreeRelax const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~GlycanTreeRelax() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

public:

	void
	set_rounds( core::Size const rounds );

	///@brief Set the layer size we will be using.  A layer is a set of glycan residues that we will be optimizing.
	///  We work our way through the layers, while the rest of the residues are virtual (not scored).
	///
	///@details
	///
	///  The distance that make up a layer.  If we have a distance of 2,
	///  we first model all glycans that are equal to or less than 2 residue distance to the root.
	///  We then slide this layer up.  So we take all residues that have a distance between 3 and 1, and so on.
	void
	set_layer_size( core::Size const layer_size );

	///@brief Set the window size.  This is the overlap of the layers during modeling.
	///
	///@details
	///
	///  The layers are slid down throught the tree of the glycan.  The window size represents the overlap in the layers.
	///  A window size of 1, means that the last residue (or residues of layer 1) from the last modeling effort, will be used again as
	///  part of the next layer.  A window size of 0, means that no residues will be re-modeled.
	///  Typically, we would want at least a window size of 1.
	void
	set_window_size( core::Size const window_size );

	///@brief Set the protocol to use a quench-like algorithm, where we work on a single glycan tree at a time until
	/// all are modeled.
	void
	set_quench_mode( bool quench_mode );

public:
	//GlcyanRelax options

	///@brief Set a boolean that we will be refining instead of de-novo modeling.
	void
	set_refine( bool const refine );

	///@brief Change the setting to do a final min/pack/min of all glycan residues or glycan residues set by the selector
	/// at the end of the protocol.
	void
	set_final_min_pack_min( bool const minpackmin );

public:

	///@brief Set the scorefunction used for modeling.
	void
	set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );

	///@brief Set a residue selector to limit the residues we will be modeling.
	/// If you are using quench mode, the selector will limit the trees to model.
	/// So the residues that are true should correspond to the start of the trees you wish to model.
	///
	/// If you wish to limit further, this is not currently supported.  So email me and I can add it.
	void
	set_selector( core::select::residue_selector::ResidueSelectorCOP selector );

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

	//GlycanTreeRelax & operator=( GlycanTreeRelax const & src );

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

	///@brief Setup classes, selectors, etc.
	void
	setup( core::pose::Pose & pose );

	bool
	is_quenched() const;

private: // data

	core::Size layer_size_ = 2;
	core::Size window_size_ = 1;
	core::Size rounds_ = 1;
	core::Size completed_quenches_ = 0;
	core::Size trees_to_model_ = 0;

	bool refine_ = false;
	bool quench_mode_ = false;
	bool final_min_pack_min_ = true;

	core::scoring::ScoreFunctionOP scorefxn_ = nullptr;
	core::select::residue_selector::ResidueSelectorCOP selector_ = nullptr;

};

std::ostream &
operator<<( std::ostream & os, GlycanTreeRelax const & mover );

} //protocols
} //carbohydrates

#endif //protocols_carbohydrates_GlycanTreeRelax_HH
