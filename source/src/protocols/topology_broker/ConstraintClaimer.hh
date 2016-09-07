// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_ConstraintClaimer_hh
#define INCLUDED_protocols_topology_broker_ConstraintClaimer_hh


// Unit Headers
#include <protocols/topology_broker/ConstraintClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh>
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// Project Headers
#include <core/pose/Pose.hh>

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
//#include <basic/Tracer.hh>
//#include <basic/options/option.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>

// option key includes


namespace protocols {
namespace topology_broker {

class ConstraintClaimer : public TopologyClaimer {
	typedef TopologyClaimer Parent;
public:
	ConstraintClaimer(); //for factory
	ConstraintClaimer( std::string cst_file, std::string tag = "" );
	ConstraintClaimer( bool CmdFlag, bool centroid = true, bool fullatom = false );

	TopologyClaimerOP clone() const override {
		return TopologyClaimerOP( new ConstraintClaimer( *this ) );
	}

	void generate_claims( claims::DofClaims& ) override;

	/// @brief type() is specifying the output name of the TopologyClaimer
	std::string type() const override {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "ConstraintClaimer";
	}

	void new_decoy() override;

	bool read_tag( std::string tag, std::istream & ) override;

	void add_constraints( core::pose::Pose& /*pose*/ ) const override;

	// ConstraintSetCOP constraints() const { return constraints_; }
	// ConstraintSetCOP fa_constraints() const { return fa_constraints_; }
	core::Real filter_weight() const { return filter_weight_; }
	std::string const& filter_name() const { return filter_name_; }
	core::Size combine_ratio() const { return combine_ratio_; }
	core::Real drop_random_rate() const { return drop_random_rate_; }
	std::string const& tag() const { return tag_; }

	void set_cst_file( std::string const& );
	void set_fullatom( bool setting );
	void set_centroid( bool setting );
	void set_skip_redundant( core::Size setting );
	void set_combine_ratio( core::Size setting );
	void set_filter_weight( core::Real weight );
	// Undefined, commenting out to fix PyRosetta build  void set_drop_random_rate( core::Real setting );
private:
	std::string filename_, fa_filename_;
	std::string cst_filename_;
	std::string tag_;
	mutable core::scoring::constraints::ConstraintSetOP constraints_;
	mutable core::scoring::constraints::ConstraintSetOP fa_constraints_;

	mutable core::pose::Pose constraint_ref_pose_;
	mutable std::string sequence_;
	mutable utility::vector1< bool > combine_exclude_res_;
	/// @brief true if constraints are active in centroid mode
	bool bCentroid_;

	/// @brief true if constraints are active in full-atom mode
	/// if true use fa_constraints if available... normal constraints otherwise
	bool bFullatom_;

	/// @brief use constraints defined via command line
	bool bCmdFlag_;


	/// @brief combine constraints randomly into Ambiguous Constraints
	core::Size combine_ratio_;  //default 1: no constraint combination

	/// @brief drop restraints with this rate randomly -- default 0 -- no dropping
	core::Real drop_random_rate_;

	/// @brief at most one constraint per residue pair ( does not look at bounds... )
	bool skip_redundant_;

	/// @brief how many residues left and right do we exclude other "redundant" constraints...
	core::Size skip_redundant_width_;

	/// @brief which weight should this constraint set have when used as a filter
	core::Real filter_weight_;
	std::string filter_name_;
}; //class ConstraintClaimer

}
}

#endif
