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


#ifndef INCLUDED_protocols_topology_broker_CoordConstraintClaimer_hh
#define INCLUDED_protocols_topology_broker_CoordConstraintClaimer_hh


// Unit Headers
#include <protocols/topology_broker/CoordConstraintClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh>
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <protocols/loops/Loops.hh>

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

class CoordConstraintClaimer : public TopologyClaimer {
	typedef TopologyClaimer Parent;
public:
	CoordConstraintClaimer(); //for factory
	~CoordConstraintClaimer() override; //for factory
	CoordConstraintClaimer( std::string const & pdb_file );

	TopologyClaimerOP clone() const override {
		return TopologyClaimerOP( new CoordConstraintClaimer( *this ) );
	}

	void generate_claims( claims::DofClaims& ) override;

	/// @brief type() is specifying the output name of the TopologyClaimer
	std::string type() const override {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "CoordConstraintClaimer";
	}


	void add_constraints( core::pose::Pose& /*pose*/ ) const override;

	void new_decoy() override;

	void new_decoy( core::pose::Pose const& ) override;

	/// @brief superimpose xyz coords in constraints_ with pose
	void superimpose( core::pose::Pose const& ) const;

protected:

	/// @brief virtual functions for IO
	bool read_tag( std::string tag, std::istream & ) override;
	void set_defaults() override;
	void init_after_reading() override;

private:
	void read_cst_pose();
	void set_cst_root();
	void generate_constraints( core::pose::Pose const& cst_pose ) const;
	void read_constraints_from_file( core::pose::Pose const& cst_pose ) const;
	std::string filename_;
	std::string cst_filename_;
	loops::Loops rigid_; //if empty ---> all residues, otherwise only on these

	mutable core::scoring::constraints::ConstraintSetOP constraints_;
	mutable std::string sequence_;
	core::scoring::func::FuncOP cst_func_;
	Size root_; //if 0 -- it's ignored. otherwise try to set fold-tree root to this.
	std::string root_from_label_;
	bool bRegenerateFromInputPose_;
	bool bUseXYZ_in_cstfile_;
	core::pose::PoseOP cst_pose_;

	/// @brief true if constraints are active in centroid mode
	bool bCentroid_;

	/// @brief true if constraints are active in full-atom mode
	bool bFullatom_;

	/// @brief tmp for backwards compatibility
	bool bLocal_;

	/// @brief add 0..perturb_ random number to xyz-coords.
	core::Real perturb_;

	/// @brief superimpose xyz coordinates with pose in add_constraints()
	bool bSuperimpose_;

	/// @brief superimpose on these residues
	loops::Loops superimpose_regions_;

}; //class CoordConstraintClaimer

}
}

#endif
