// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <protocols/filters/Filter.hh>
#include <protocols/fldsgn/topology/HelixPairing.fwd.hh>
#include <protocols/fldsgn/topology/HSSTriplet.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <protocols/forge/components/VarLengthBuild.fwd.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>

// Core headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>
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
	  core::pose::Pose const & );

  /// @brief Return the name of this mover.
  virtual std::string get_name() const;

  /// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::filters::FilterOP clone() const;

  /// @brief Apply the FoldabilityFilter. Overloaded apply function from filter base class.
	virtual protocols::filters::FilterOP fresh_instance() const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual bool apply( core::pose::Pose const & pose ) const;

	core::Real compute( core::pose::Pose const & pose ) const;

	// mutators
public:
	void set_start_res( core::Size const startval );
	void set_end_res( core::Size const endval );

	// protected functions
protected:
	/// @brief determine start and end window, from selector if necessary
	void choose_start_and_end(
			core::Size & start,
			core::Size & end,
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

	/// @brief deletes the segment from start to end (inclusive) from the pose
	void delete_segment(
			core::pose::Pose & pose,
			core::Size const start,
			core::Size const end ) const;

	/// @brief performs setup on the fragment insertion machinery
	void setup_vlb(
			std::string const & aa,
			std::string const & ss,
			utility::vector1< std::string > const & abego,
			core::Size const start,
			core::Size const end ) const;

	/// @brief performs fragment insertion and returns number of successful builds. Assumes setup_vlb() has already been called
	core::Size fragment_insertion(
			core::pose::Pose const & pose,
			core::Size const end,
			core::conformation::Residue const & end_res ) const;

private:   // options
	/// @brief the motif to try to build
	std::string motif_;
	/// @brief try this number of times to build the motif
	core::Size tries_;
	/// @brief residue number to start building
	core::Size start_res_;
	/// @brief residue number to stop building
	core::Size end_res_;
	/// @brief if true, abego values in the input pose will be ignored for the segment we are rebuilding
	bool ignore_pose_abego_;
	/// @brief if true, poses will be outputted for each foldability step (default=false)
	bool output_poses_;

private:   // other data
	/// @brief residue selector to identify positions to rebuild
	core::pack::task::residue_selector::ResidueSelectorCOP selector_;
	/// @brief the vlb object for testing foldability
	protocols::forge::components::VarLengthBuildOP vlb_;
	mutable std::string cached_aa_;
	mutable std::string cached_ss_;
	mutable core::Size cached_start_;
	mutable core::Size cached_end_;
};


} // filters
} // denovo_design
} // devel

#endif
