// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody.design/AntibodyDesignModeler.hh
/// @brief Handles modeling of the antibody.  Before, after, and during design
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_AntibodyDesignModeler_hh
#define INCLUDED_protocols_antibody_design_AntibodyDesignModeler_hh

#include <protocols/antibody/design/AntibodyDesignModeler.fwd.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>


#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace antibody {
namespace design {
	
using namespace core::scoring;
using core::pose::Pose;
using std::string;
using core::Real;
using utility::vector1;

///@brief Class that runs the modeling for AntibodyDesign protocols.  Could be moved to antibody.
class AntibodyDesignModeler : public utility::pointer::ReferenceCount {
            
public:
	AntibodyDesignModeler(AntibodyInfoOP & ab_info);
	virtual ~AntibodyDesignModeler();

	///@brief Set to model CDRs.  Default is all of them false.
	void
	set_cdr(CDRNameEnum const cdr, bool setting);
	
	///@brief Set to model only one cdr, or all others but one.
	void
	set_cdr_only(CDRNameEnum const cdr, bool setting);
	
	///@brief Set a range of CDRs.
	void
	set_cdr_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool setting);

	///////////////////////////////////////////////////////////////////////////////////////
	// Grafting
	//
	//
	///@brief Uses (default) grafting to copy cdrs from one antibody to another.
	void
	copy_cdrs(core::pose::Pose& from_pose, core::pose::Pose& to_pose, core::Size overhang);


	///////////////////////////////////////////////////////////////////////////////////////
	// Minimization
	//
	//
	///@brief Relax CDRs only either using FastRelax or CentroidRelax with Cluster constraints.
	/// Constraints should already be set.
	void
	relax_cdrs(Pose & pose, ScoreFunctionOP const & scorefxn, bool centroid_mode);
	
	///@brief Vanilla minimizer.
	void
	minimize_cdrs(Pose & pose, ScoreFunctionOP const & scorefxn);
	
	///@brief Design CDRs - Basic fixedbb design here.
	//void
	//design_cdrs(Pose & pose, ScoreFunctionOP scorefxn);

	///@brief Repack the interface between Antibody and Antigen
	//void
	//repack_interface(Pose & pose, ScoreFunctionOP scorefxn);

	///@brief Repack the CDRs given.  Nothing special.
	//void
	//repack_CDRs(Pose & pose, ScoreFunctionOP const scorefxn);

	///@brief Repack CDRs and their neighbors within distance.
	//void
	//repack_CDRs_and_neighbors(Pose & pose, ScoreFunctionOP scorefxn, core::Real distance=3.0);

	///@brief Repack CDRs and the antigen interface.
	//void
	//repack_CDRs_and_antigen_interface(Pose & pose, ScoreFunctionOP scorefxn);

	///@brief Run Brian's awesome snugdock.
	//void
	//run_snugdock(Pose & pose, ScoreFunctionOP const scorefxn);

private:
            
            
	AntibodyInfoOP ab_info_;
            
            
	vector1 <bool> cdrs_; //These cdrs are cdr's to work on.  Whether copying, refining, designing, etc.etc.

            

            
};
}
}
}



#endif //INCLUDED_protocols_antibody_design_AntibodyDesignMover_hh
