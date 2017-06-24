// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/devel/denovo_design/filters/FoldabilityFilter.hh
/// @brief Tom's Denovo Protocol. This is freely mutable and used for playing around with stuff
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)


#ifndef INCLUDED_devel_denovo_design_filters_FoldabilityFilter_hh
#define INCLUDED_devel_denovo_design_filters_FoldabilityFilter_hh

// Unit headers
#include <devel/denovo_design/filters/FoldabilityFilter.fwd.hh>

// Project headers
#include <devel/denovo_design/calculators/CavityCalculator.hh>
#include <protocols/denovo_design/components/Picker.fwd.hh>

// Protocol headers
#include <protocols/filters/Filter.hh>
#include <protocols/fldsgn/topology/HelixPairing.fwd.hh>
#include <protocols/fldsgn/topology/HSSTriplet.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <protocols/forge/components/VarLengthBuild.fwd.hh>
#include <protocols/parser/BluePrint.fwd.hh>

// Core headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/sc/MolecularSurfaceCalculator.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.fwd.hh>

// C++ headers
#include <string>

#include <core/io/silent/silent.fwd.hh>
#include <utility/vector1.hh>


namespace devel {
namespace denovo_design {
namespace filters {

class FoldabilityFilter : public protocols::filters::Filter {
	typedef protocols::denovo_design::components::Picker Picker;
	typedef protocols::denovo_design::components::PickerOP PickerOP;
public:
	/// @brief Initialize FoldabilityFilter
	FoldabilityFilter();

	/// @brief virtual constructor to allow derivation
	virtual ~FoldabilityFilter();

	/// @brief Parses the FoldabilityFilter tags
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;

	/// @brief Return the name of this mover.
	virtual std::string get_name() const;

	/// @brief return a fresh instance of this class in an owning pointer
	protocols::filters::FilterOP clone() const override;

	/// @brief Apply the FoldabilityFilter. Overloaded apply function from filter base class.
	protocols::filters::FilterOP fresh_instance() const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	bool apply( core::pose::Pose const & pose ) const override;

	core::Real compute( core::pose::Pose const & pose ) const;
	core::Real compute_segment(
		core::pose::Pose const & pose,
		core::select::residue_selector::ResidueRange const & segment ) const;

	// mutators
public:
	void clear_segments();
	void add_segment( core::Size const startval, core::Size const endval );

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	// protected functions
protected:
	/// @brief determine start and end window, from selector if necessary
	void choose_start_and_end(
		core::Size & start,
		core::Size & end,
		core::Size const segment,
		core::pose::Pose const & pose ) const;

	/// @brief gets aa string, ss string, and abego vector for the area to rebuild
	void get_aa_ss_abego(
		std::string & aa,
		std::string & ss,
		utility::vector1< std::string > & abego,
		core::Size const start,
		core::Size & end,
		core::pose::Pose const & pose ) const;

	/// @brief gets non-const version the pose for the filter to work on
	core::pose::PoseOP generate_pose( core::pose::Pose const & pose ) const;

	/// @brief prepares the pose/segment from start to end for insertion
	void prepare_pose(
		core::pose::Pose & pose,
		core::Size const start,
		core::Size const end ) const;

	/// @brief performs fragment picking/other preparations for building
	protocols::moves::MoverOP
	create_fragment_insertion_mover(
		std::string const & complete_aa,
		std::string const & complete_ss,
		utility::vector1< std::string > const & complete_abego,
		utility::vector1< core::Size > const & chain_endings,
		core::Size const start,
		core::Size const end ) const;

	/// @brief performs fragment insertion and returns number of successful builds. Assumes setup_vlb() has already been called
	core::Size fragment_insertion(
		core::pose::Pose const & pose,
		protocols::moves::Mover & fragment_mover,
		core::Size const end,
		core::conformation::Residue const & end_res ) const;

	/// @brief queries the abego db and calculates the best scoring loop
	core::Real abegodb_score(
		core::pose::Pose const & pose,
		utility::vector1< std::string > const & abego,
		core::Size const start,
		core::Size const stop ) const;

private:   // options
	/// @brief the motif to try to build
	std::string motif_;
	/// @brief try this number of times to build the motif
	core::Size tries_;
	/// @brief residue segments to rebuild
	core::select::residue_selector::ResidueRanges segments_;
	/// @brief "success" is achieved when distance is below this threshold
	core::Real distance_threshold_;
	/// @brief if true, abego values in the input pose will be ignored for the segment we are rebuilding
	bool ignore_pose_abego_;
	/// @brief if true, amino acid identity in the input pose will be used to pick fragments
	bool use_sequence_;
	/// @brief if true, poses will be outputted for each foldability step (default=false)
	bool output_poses_;

private:   // other data
	/// @brief scorefunction to use for folding
	core::scoring::ScoreFunctionOP scorefxn_;
	/// @brief residue selector to identify positions to rebuild
	core::select::residue_selector::ResidueSelectorCOP selector_;
	/// @brief fragment picker
	PickerOP picker_;
	/// @brief the vlb object for testing foldability
	protocols::forge::components::VarLengthBuildOP vlb_;
};


} // filters
} // denovo_design
} // devel

#endif
