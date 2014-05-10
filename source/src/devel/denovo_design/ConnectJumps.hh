// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/ConnectJumps.hh
/// @brief The ConnectJumps Protocol
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_devel_denovo_design_ConnectJumps_hh
#define INCLUDED_devel_denovo_design_ConnectJumps_hh

// Unit headers
#include <devel/denovo_design/ConnectJumps.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers

// Package headers
#include <protocols/forge/components/VarLengthBuild.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <protocols/filters/Filter.hh>

//// C++ headers
#include <string>

#include <core/io/silent/silent.fwd.hh>
#include <utility/vector1.hh>


namespace devel {
namespace denovo_design {

class ConnectJumps : public protocols::moves::Mover {
public:

	/// @brief default constructor
	ConnectJumps();

  /// @brief virtual constructor to allow derivation
	virtual ~ConnectJumps();

  /// @brief Parses the ConnectJumpsTags
	void parse_my_tag(
	  utility::tag::TagCOP tag,
	  basic::datacache::DataMap & data,
	  protocols::filters::Filters_map const &,
	  protocols::moves::Movers_map const &,
	  core::pose::Pose const &
	);

  /// @brief Return the name of this mover.
  virtual std::string get_name() const;

  /// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::moves::MoverOP clone() const;

  /// @brief Apply the ConnectJumps. Overloaded apply function from mover base class.
	virtual void apply( core::pose::Pose & pose );

	/// @brief creates a Ca coordinate constraint for residue resi
	core::scoring::constraints::ConstraintOP
	create_coordinate_cst( core::pose::Pose const & pose,
			core::Size const resi ) const;

private: // options
	std::string motif_;
	core::Size jump1_;
	core::Size jump2_;
	/// @brief the number of residues to be rebuilt before and after the jump -- coordinate constraints will be applied to these.
	core::Size overlap_;
	core::scoring::ScoreFunctionCOP scorefxn_;

private:   // other data
	/// @brief the vlb (so that fragments can be cached)
	protocols::forge::components::VarLengthBuildOP vlb_;
	/// @brief cached data to alert the VLB is something has changed... this will allow it to cache fragments
	std::string cached_ss_;
	std::string cached_aa_;
	core::Size cached_start_;
	core::Size cached_end_;
};




} // denovo_design
} // devel

#endif
