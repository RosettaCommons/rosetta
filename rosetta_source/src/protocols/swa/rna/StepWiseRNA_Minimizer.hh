// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_Minimizer_hh
#define INCLUDED_protocols_swa_rna_StepWiseRNA_Minimizer_hh

//#include <numeric/xyzMatrix.hh>
//#include <numeric/xyzVector.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <string>
#include <map>

#include <core/id/AtomID.fwd.hh>


namespace protocols {
namespace swa {
namespace rna {

//	typedef std::map< std::string, core::pose::PoseOP > PoseList;

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class StepWiseRNA_Minimizer: public protocols::moves::Mover {
public:

	//constructor!
	StepWiseRNA_Minimizer(
		utility::vector1 <pose_data_struct2> const & pose_data_list,
		StepWiseRNA_JobParametersCOP & job_parameters
	);

	//destructor -- necessary?
	~StepWiseRNA_Minimizer();

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );
	virtual std::string get_name() const;

	void
	set_silent_file( std::string const & silent_file );

	void
	set_move_map_list(utility::vector1 <core::kinematics::MoveMap> const & move_map_list);

	void
	set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

	core::io::silent::SilentFileDataOP & silent_file_data();

	void
	set_base_centroid_screener( StepWiseRNA_BaseCentroidScreenerOP & screener );

private:

	utility::vector1 <core::kinematics::MoveMap>
	Get_default_movemap( core::pose::Pose const & pose ) const;

	void
	Figure_out_moving_residues( core::kinematics::MoveMap & mm, core::pose::Pose const & pose ) const;

	bool
	Pose_screening(core::pose::Pose & pose, std::string tag, core::io::silent::SilentFileData & silent_file_data) const;

private:

	utility::vector1 <pose_data_struct2> const pose_data_list_;
	StepWiseRNA_JobParametersCOP job_parameters_;

	core::io::silent::SilentFileDataOP sfd_;
	utility::vector1 <core::kinematics::MoveMap> move_map_list_;
	core::scoring::ScoreFunctionOP scorefxn_;
	std::string silent_file_;
	bool const verbose_;
	bool const native_screen_;
	bool const rmsd_cutoff_;
	std::map< core::id::AtomID, core::id::AtomID > pose_to_native_map_;

	utility::vector1< core::Size > fixed_res_;

	StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener_;

};

}
} //swa
} // protocols

#endif
