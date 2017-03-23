// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file loopRNA_minimizer.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_LoopCloser_hh
#define INCLUDED_protocols_rna_RNA_LoopCloser_hh

#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>

#include <core/types.hh>

//// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace rna {
namespace movers {

/// @brief The RNA de novo structure modeling protocol
class RNA_LoopCloser: public protocols::moves::Mover {
public:
	/// @brief Construct the protocol object
	RNA_LoopCloser();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const {
		return protocols::moves::MoverOP( new RNA_LoopCloser(*this) );
	}

	/// @brief Apply the loop-rebuild protocol to the input pose
	using protocols::moves::Mover::apply;
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose, std::map< Size, Size > const & connections );

	/// @brief Apply the loop-rebuild protocol to the input pose
	core::Real apply( core::pose::Pose & pose, std::map< Size, Size > const & connections, Size const & cutpoint );

	/// @brief Apply the loop-rebuild protocol to the input pose
	core::Real apply( core::pose::Pose & pose, Size const & cutpoint );

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose, utility::vector1< Size > const & cutpoints );

	//Need get/set functions for options.

	void fast_scan( bool const & setting ) { fast_scan_ = setting; }

	bool
	check_closure( core::pose::Pose const & pose, core::Size const i, core::Real ccd_tolerance = -1.0 /* -1.0 means use absolute_ccd_tolerance_*/ );

	bool
	check_closure( core::pose::Pose const & pose, core::Real ccd_tolerance = -1.0 /* -1.0 means use absolute_ccd_tolerance_*/ );

	core::Real
	get_dist_err( core::pose::Pose const & pose,
		Size const cutpoint
	) const;

	void
	set_verbose( bool const setting ){ verbose_ = setting; }

	void
	set_atom_level_domain_map( toolbox::AtomLevelDomainMapCOP setting ){ atom_level_domain_map_ = setting; }

	utility::vector1< core::Size >
	get_cutpoints_closed( core::pose::Pose const & pose ) const;

private:

	bool
	passes_fast_scan( core::pose::Pose & pose, Size const i ) const;

	// Returns final coordinate error.
	core::Real
	rna_ccd_close( core::pose::Pose & pose, std::map< Size, Size > const & connections, Size const & cutpoint ) const;

	core::Real
	get_chainbreak_xyz( core::pose::Pose const & pose,
		Size const cutpoint,
		utility::vector1< core::Vector > & upstream_xyzs,
		utility::vector1< core::Vector > & downstream_xyzs,
		Size const cutpoint_next_input = 0
	) const;

	core::Real
	get_gap_distance( core::pose::Pose & pose,
		Size const cutpoint
	) const;


private:
	//Make these options:
	bool verbose_;
	Size NUM_ROUNDS_;
	bool check_tolerance_;
	core::Real ccd_tolerance_;
	core::Real absolute_ccd_tolerance_;
	core::Real attempt_closure_cutoff_;
	core::Real gap_distance_cutoff_;
	bool fast_scan_;
	toolbox::AtomLevelDomainMapCOP atom_level_domain_map_;
}; // class RNA_LoopCloser


} //movers
} //rna
} //protocols

#endif
