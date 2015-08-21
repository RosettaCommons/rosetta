// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RepeatPropagationMover.cc
/// @brief Can repeat both a symmetric, non symmetric & broken pose
/// @author TJ Brunette

// Unit headers
#include <protocols/simple_moves/RepeatPropagationMover.hh>
#include <protocols/simple_moves/RepeatPropagationMoverCreator.hh>


// Project Headers

#include <basic/Tracer.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <protocols/idealize/idealize.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/relax/CentroidRelax.hh>
#include <protocols/relax/cst_util.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>


#include <string>


using namespace core;
using namespace std;
using utility::vector1;

namespace protocols {
namespace simple_moves {
static basic::Tracer TR( "protocols.simple_moves.RepeatPropagationMover" );
std::string RepeatPropagationMoverCreator::keyname() const
{
	return RepeatPropagationMoverCreator::mover_name();
}


protocols::moves::MoverOP
RepeatPropagationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RepeatPropagationMover );
}


std::string RepeatPropagationMoverCreator::mover_name(){
	return "RepeatPropagation";
}

RepeatPropagationMover::RepeatPropagationMover():moves::Mover("RepeatPropagation"){}

void RepeatPropagationMover::apply(core::pose::Pose & pose) {
	if ( pose.total_residue()== last_res_ ) {
		throw utility::excn::EXCN_RosettaScriptsOption("can not handle situation where last_res is being repeated (yet...)");
	}
	if ( repeat_without_replacing_pose_ ) {
		copy_phi_psi_omega(pose,pose);
	} else {
		core::pose::PoseOP repeat_poseOP( new core::pose::Pose() );
		initialize_repeat_pose(pose,*repeat_poseOP);
		copy_phi_psi_omega(pose,*repeat_poseOP);
		pose = *repeat_poseOP;
	}
}

void RepeatPropagationMover::initialize_repeat_pose( core::pose::Pose & pose, core::pose::Pose & repeat_pose){
	duplicate_residues_by_type(pose,repeat_pose);
}

void RepeatPropagationMover::duplicate_residues_by_type( core::pose::Pose & pose, core::pose::Pose & repeat_pose){
	//N-terminal cap-----------------------
	if ( maintain_cap_seq_and_structure_||maintain_cap_sequence_alone_ ) {
		for ( Size res=1; res<=nTerm_cap_size_; ++res ) {
			repeat_pose.append_residue_by_bond(pose.residue(res),true/*ideal bonds*/);
		}
	}
	//Repeat-----------------------
	if ( first_res_ == 1 ) {
		remove_lower_terminus_type_from_pose_residue(pose, 1);
	}
	if ( last_res_ == pose.total_residue() ) {
		remove_upper_terminus_type_from_pose_residue(pose,pose.total_residue());
	}
	Size tmp_numb_repeats = numb_repeats_;
	if ( maintain_cap_sequence_alone_ && nTerm_cap_size_>0 ) {
		tmp_numb_repeats--;
	}
	if ( maintain_cap_sequence_alone_ && cTerm_cap_size_>0 ) {
		tmp_numb_repeats--;
	}
	for ( Size rep=1; rep<=tmp_numb_repeats; ++rep ) {
		for ( Size res=first_res_; res<=last_res_; ++res ) {
			repeat_pose.append_residue_by_bond(pose.residue(res),true/*ideal bonds*/);

		}
	}
	if ( first_res_==1 ) {
		add_lower_terminus_type_to_pose_residue(pose,1);
	}
	if ( last_res_ == pose.total_residue() ) {
		add_upper_terminus_type_to_pose_residue(pose,pose.total_residue());
	}
	//C-terminal cap----------------
	if ( maintain_cap_seq_and_structure_||maintain_cap_sequence_alone_ ) {
		for ( Size res=pose.total_residue()-cTerm_cap_size_+1; res<=pose.total_residue(); ++res ) {
			repeat_pose.append_residue_by_bond(pose.residue(res),true /*ideal bonds*/ );
		}
	}
	//teminal residues
	add_lower_terminus_type_to_pose_residue(repeat_pose,1);
	add_upper_terminus_type_to_pose_residue(repeat_pose,repeat_pose.total_residue());
}

void RepeatPropagationMover::copy_phi_psi_omega(core::pose::Pose & pose, core::pose::Pose & repeat_pose){
	//N-terminal cap-----------------------
	if ( maintain_cap_seq_and_structure_ ) {
		for ( Size res=1; res<=nTerm_cap_size_; ++res ) {
			repeat_pose.set_phi(res , pose.phi(res));
			repeat_pose.set_psi(res, pose.psi(res));
			repeat_pose.set_omega(res, pose.omega(res));
			repeat_pose.set_secstruct(res,pose.secstruct(res));
		}
	}
	//Repeat region
	Real loop_phi = 0;
	Real loop_psi = 0;
	Real loop_omega = 0;
	char loop_secstruct = 'H';
	Size segment_length = last_res_ - first_res_+1;
	for ( Size ii=1; ii<=segment_length; ++ii ) {
		Size pose_res = first_res_-1+ii;
		Size repeat_pose_res = ii;
		if ( maintain_cap_seq_and_structure_ ) {
			repeat_pose_res =nTerm_cap_size_+ii;
		}
		if ( pose_res == 1 ) { //if the first and last positions are involved, loop around
			loop_phi = pose.phi(pose_res+segment_length);
			loop_psi = pose.psi(pose_res+segment_length);
			loop_omega = pose.omega(pose_res+segment_length);
			loop_secstruct = pose.secstruct(pose_res+segment_length);
		} else {
			loop_phi = pose.phi(pose_res);
			loop_psi = pose.psi(pose_res);
			loop_omega =pose.omega(pose_res);
			loop_secstruct = pose.secstruct(pose_res);
		}
		for ( int rep=0; rep<(int)numb_repeats_; ++rep ) {
			repeat_pose.set_phi(repeat_pose_res+( segment_length*rep), loop_phi );
			repeat_pose.set_psi(repeat_pose_res+( segment_length*rep), loop_psi );
			repeat_pose.set_omega( repeat_pose_res+(segment_length*rep),loop_omega );
			repeat_pose.set_secstruct( repeat_pose_res+(segment_length*rep),loop_secstruct );
		}
	}
	//C-terminal cap-----------------------
	if ( maintain_cap_seq_and_structure_ ) {
		for ( Size ii=1; ii<=cTerm_cap_size_; ++ii ) {
			Size pose_res = pose.total_residue()-cTerm_cap_size_+ii;
			Size repeat_pose_res = repeat_pose.total_residue()-cTerm_cap_size_+ii;
			repeat_pose.set_phi(repeat_pose_res , pose.phi(pose_res));
			repeat_pose.set_psi(repeat_pose_res, pose.psi(pose_res));
			repeat_pose.set_omega(repeat_pose_res, pose.omega(pose_res));
			repeat_pose.set_secstruct(repeat_pose_res,pose.secstruct(pose_res));
		}
	}
}


std::string RepeatPropagationMover::get_name() const { return "RepeatPropagationMover"; }

void RepeatPropagationMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & //pose
)
{
	if ( !tag->hasOption("first_template_res") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("repeat mover requires the first residue be entered with first_res tag");
	}
	first_res_ = tag->getOption<Size>("first_template_res");
	if ( !tag->hasOption("last_template_res") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("repeat mover requires the last residue to be entered with last_res tag");
	}
	last_res_ = tag->getOption<Size>("last_template_res");
	if ( !tag->hasOption("numb_repeats") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("repeat mover requires the number of repeats be entered with numb_repeats tag");
	}
	numb_repeats_ = tag->getOption<Size>("numb_repeats");
	repeat_without_replacing_pose_ = tag->getOption<bool>("repeat_without_replacing_pose",false);//for speed.
	maintain_cap_seq_and_structure_ = tag->getOption<bool>("maintain_cap_seq_and_structure",false);
	maintain_cap_sequence_alone_ = tag->getOption<bool>("maintain_cap_sequence_alone",false);
	nTerm_cap_size_=0;
	cTerm_cap_size_=0;
	if ( maintain_cap_sequence_alone_ || maintain_cap_seq_and_structure_ ) {
		if ( (!tag->hasOption("nTerm_cap_size"))||(!tag->hasOption("cTerm_cap_size")) ) {
			throw utility::excn::EXCN_RosettaScriptsOption("repeat mover requires cTerm_cap and nTerm_cap defined if trying to maintain cap sequence or structure");
		}
		nTerm_cap_size_ = tag->getOption<Size>("nTerm_cap_size");
		cTerm_cap_size_ = tag->getOption<Size>("cTerm_cap_size");
	}
	Size repeat_size = last_res_-first_res_+1;
	if ( maintain_cap_sequence_alone_ ) {
		if ( (nTerm_cap_size_ != 0) && (nTerm_cap_size_ != repeat_size) ) {
			throw utility::excn::EXCN_RosettaScriptsOption("if maintaining nTermCap sequence please set nTerm_cap equal to the repeat_size");
		}
	}
	if ( (cTerm_cap_size_ != 0) && (cTerm_cap_size_ != repeat_size) ) {
		throw utility::excn::EXCN_RosettaScriptsOption("if maintaining cTermCap sequence please set cTerm_cap equal to the repeat_size");
	}
}


} // moves
} // protocols
