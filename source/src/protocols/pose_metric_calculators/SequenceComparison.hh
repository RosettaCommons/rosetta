// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


//////////////////////////////////////////////////////////////////////
///
/// @brief
/// Compare the sequences between a native and designed protein
///
/// @details
/// This is an implementation taken from Ron Jacak, Douglas Renfrew, Matt O Mera.
/// The main function that is called is the get_sequence_recovery() function. You can
/// pass this function a list of native pdbs and designed pdbs, or just 1 native and 1
/// designed pdb. The sequence recovery will be output in a file called sequencerecovery.txt
/// along with a substitution matrix in a file called submatrix.txt
///
///
/// @references
/// "Native sequences are close to optimal" paper
///
///
/// @author
/// Ron Jacak,
/// Douglas Renfrew (renfrew@nyu.edu) ( added rotamer recovery, cleanup )
/// Steven Combs (moved it into a general use class)
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_SequenceComparison_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_SequenceComparison_hh

#include <utility/vector1.hh>

// C++ headers
#include <set>

#include <core/types.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace pose_metric_calculators {

class SequenceComparison{
public:
	//constructor to assign surface and core cutoff;
	SequenceComparison(core::Size surface, core::Size core) :
		surface_exposure_(surface),
		core_cutoff_(core)
	{

	}


	//default constructor
	SequenceComparison() :
		surface_exposure_(16),
		core_cutoff_(24)
	{

	}


	/// @brief main function that is called. calls measure_sequence_recovery
	void get_sequence_recovery(utility::vector1<core::pose::Pose> & native_poses, utility::vector1<core::pose::Pose> & redesign_poses);

	/// @brief overflowed main function to compare only two proteins
	void get_sequence_recovery(core::pose::Pose & native, core::pose::Pose & designed);

	/// @brief measures the sequence recovery of a list of native proteins and a list of designed proteins. Outputs files to sequencerecovery.txt
	void measure_sequence_recovery( utility::vector1<core::pose::Pose> & native_poses, utility::vector1<core::pose::Pose> & redesign_poses);
	void fill_num_neighbors( core::pose::Pose & pose, utility::vector1< core::Size > & num_nbs );

	/// @brief sets up the task factory used for determing what the neighbor counts...what is on the surface and what is in the core
	core::pack::task::TaskFactoryOP setup_tf( core::pack::task::TaskFactoryOP task_factory_ );
	std::set< core::Size > fill_designable_set( core::pose::Pose & pose, core::pack::task::TaskFactoryOP & tf );

private:

	/// @brief what is the surface cutoff that you are using? default 16
	core::Size surface_exposure_;

	/// @brief what is the core cutoff that you are using? default 24
	core::Size core_cutoff_;

};


}
}

#endif /* SEQUENCECOMPARISON_HH_ */
