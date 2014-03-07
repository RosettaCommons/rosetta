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


#ifndef INCLUDED_protocols_stepwise_StepWiseResidueSampler_hh
#define INCLUDED_protocols_stepwise_StepWiseResidueSampler_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>
#include <protocols/stepwise/MainChainTorsionSet.hh> // should make a .fwd.hh probably
#include <string>
#include <map>

//Auto Headers


namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

//	typedef std::map< std::string, core::pose::PoseOP > PoseList;

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class StepWiseResidueSampler: public protocols::moves::Mover {
public:

	//constructor!
	StepWiseResidueSampler(
		utility::vector1< Size > const & moving_residues,
		utility::vector1< MainChainTorsionSetList > const & main_chain_torsion_set_lists );

	//destructor -- necessary?
	~StepWiseResidueSampler();

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;

	void
	set_silent_file( std::string const & setting );

	void
	set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

	core::io::silent::SilentFileDataOP & silent_file_data();

private:

	void
	sample_residues( core::pose::Pose & pose );

	void
	quick_output(
		core::pose::Pose & pose,
		std::string const & tag	 );

	void
	initialize_moving_residues_including_junction( Size const & nres );

	void
	initialize_green_packer( core::Size const & nres );

private:

	utility::vector1< Size > const moving_residues_;
	utility::vector1< MainChainTorsionSetList > const main_chain_torsion_set_lists_;
	//		PoseList pose_list_;
	core::scoring::ScoreFunctionOP scorefxn_;
	protocols::simple_moves::GreenPackerOP green_packer_;
	std::string silent_file_;

	core::io::silent::SilentFileDataOP sfd_;

};

} //protein
} //sampling
} //stepwise
} //protocols

#endif
