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
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>
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
#include <protocols/toolbox/superimpose.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


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
		if ( maintain_cap_ ) {
			add_caps(pose,*repeat_poseOP);
		}
		pose = *repeat_poseOP;
	}
}

void RepeatPropagationMover::initialize_repeat_pose( core::pose::Pose & pose, core::pose::Pose & repeat_pose){
	duplicate_residues_by_type(pose,repeat_pose);
}

void RepeatPropagationMover::add_caps(core::pose::Pose & pose, core::pose::Pose & repeat_pose){
	//create c-term subpose first
	Size repeat_size = last_res_-first_res_+1;
	if ( cTerm_cap_size_>0 ) {
		core::pose::Pose cCap_pose;
		vector1<Size> cTerm_pose_positions;
		for ( Size res=pose.total_residue()-cTerm_cap_size_-repeat_size; res<=pose.total_residue(); ++res ) {
			cTerm_pose_positions.push_back(res);
		}
		core::kinematics::FoldTree mod_ft(cTerm_pose_positions.size());
		core::pose::create_subpose(pose,cTerm_pose_positions,mod_ft,cCap_pose);
		Size start_overlap_parent = 0;
		Size end_overlap_parent = 0;
		Size start_overlap_pose = 0;
		Size end_overlap_pose = 0;
		determine_overlap(repeat_pose,cCap_pose,repeat_size,5, "c_term", start_overlap_parent,end_overlap_parent,start_overlap_pose,end_overlap_pose);
		generate_overlap(repeat_pose,cCap_pose, "c_term",start_overlap_parent,end_overlap_parent,start_overlap_pose,end_overlap_pose);
	}
	if ( nTerm_cap_size_>0 ) {
		core::pose::Pose nCap_pose;
		vector1<Size> nTerm_pose_positions;
		for ( Size res=1; res<=nTerm_cap_size_+repeat_size; ++res ) {
			nTerm_pose_positions.push_back(res);
		}
		core::kinematics::FoldTree mod_ft(nTerm_pose_positions.size());
		core::pose::create_subpose(pose,nTerm_pose_positions,mod_ft,nCap_pose);
		Size start_overlap_parent = 0;
		Size end_overlap_parent = 0;
		Size start_overlap_pose = 0;
		Size end_overlap_pose = 0;
		determine_overlap(repeat_pose,nCap_pose,repeat_size,5, "n_term",start_overlap_parent,end_overlap_parent,start_overlap_pose,end_overlap_pose);
		generate_overlap(repeat_pose,nCap_pose, "n_term",start_overlap_parent,end_overlap_parent,start_overlap_pose,end_overlap_pose);
	}
}

void RepeatPropagationMover::duplicate_residues_by_type( core::pose::Pose & pose, core::pose::Pose & repeat_pose){
	//Repeat-----------------------
	if ( first_res_ == 1 ) {
		remove_lower_terminus_type_from_pose_residue(pose, 1);
	}
	if ( last_res_ == pose.total_residue() ) {
		remove_upper_terminus_type_from_pose_residue(pose,pose.total_residue());
	}
	Size tmp_numb_repeats = numb_repeats_;
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
	//teminal residues
	add_lower_terminus_type_to_pose_residue(repeat_pose,1);
	add_upper_terminus_type_to_pose_residue(repeat_pose,repeat_pose.total_residue());
}

void RepeatPropagationMover::copy_phi_psi_omega(core::pose::Pose & pose, core::pose::Pose & repeat_pose){
	//Repeat region
	Real loop_phi = 0;
	Real loop_psi = 0;
	Real loop_omega = 0;
	char loop_secstruct = 'H';
	Size segment_length = last_res_ - first_res_+1;
	for ( Size ii=1; ii<=segment_length; ++ii ) {
		Size pose_res = first_res_-1+ii;
		Size repeat_pose_res = ii;
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
}

void RepeatPropagationMover::determine_overlap(Pose const pose, Pose & parent_pose,Size overlap_max_length,Size overlap_range, std::string overlap_location_pose,Size & start_overlap_parent, Size & end_overlap_parent, Size & start_overlap_pose, Size & end_overlap_pose){
	using namespace core::scoring;
	if ( overlap_location_pose== "c_term" ) {
		Size initial_start_res_parent=1;
		Size end_res_parent=overlap_max_length;
		Size initial_start_res_pose=pose.total_residue()-overlap_max_length+1;
		Size end_res_pose=pose.total_residue();
		Real best_rmsd = 9999;
		for ( Size ii=0; ii<overlap_range; ++ii ) {
			utility::vector1<Size> parent_posepositions;
			utility::vector1<Size> pose_positions;
			Size start_res_parent = initial_start_res_parent+ii;
			Size start_res_pose = initial_start_res_pose+ii;
			for ( Size jj=start_res_parent; jj<=end_res_parent; ++jj ) {
				parent_posepositions.push_back(jj);
			}
			for ( Size jj=start_res_pose; jj<=end_res_pose; ++jj ) {
				pose_positions.push_back(jj);
			}
			pose::Pose parent_poseslice;
			pose::Pose ref_pose_slice;
			pdbslice(parent_poseslice,parent_pose,parent_posepositions);
			pdbslice(ref_pose_slice,pose,pose_positions);
			Real rmsd = CA_rmsd(ref_pose_slice,parent_poseslice);
			if ( rmsd < best_rmsd ) {
				best_rmsd = rmsd;
				start_overlap_parent=start_res_parent;
				end_overlap_parent=end_res_parent;
				start_overlap_pose=start_res_pose;
				end_overlap_pose=end_res_pose;
			}
		}
	}
	if ( overlap_location_pose== "n_term" ) {
		Size start_res_parent=parent_pose.total_residue()-overlap_max_length+1;
		Size initial_end_res_parent=parent_pose.total_residue();
		Size start_res_pose=1;
		Size initial_end_res_pose=overlap_max_length;
		Real best_rmsd = 9999;
		for ( Size ii=0; ii<overlap_range; ++ii ) {
			utility::vector1<Size> parent_posepositions;
			utility::vector1<Size> pose_positions;
			Size end_res_parent = initial_end_res_parent-ii;
			Size end_res_pose = initial_end_res_pose-ii;
			for ( Size jj=start_res_parent; jj<=end_res_parent; ++jj ) {
				parent_posepositions.push_back(jj);
			}
			for ( Size jj=start_res_pose; jj<=end_res_pose; ++jj ) {
				pose_positions.push_back(jj);
			}
			pose::Pose parent_poseslice;
			pose::Pose ref_pose_slice;
			pdbslice(parent_poseslice,parent_pose,parent_posepositions);
			pdbslice(ref_pose_slice,pose,pose_positions);
			Real rmsd = CA_rmsd(ref_pose_slice,parent_poseslice);
			if ( rmsd < best_rmsd ) {
				best_rmsd = rmsd;
				start_overlap_parent=start_res_parent;
				end_overlap_parent=end_res_parent;
				start_overlap_pose=start_res_pose;
				end_overlap_pose=end_res_pose;

			}
		}
	}
}

void RepeatPropagationMover::generate_overlap(Pose & pose, Pose & parent_pose, std::string overlap_location_pose,Size start_overlap_parent, Size end_overlap_parent, Size start_overlap_pose, Size end_overlap_pose){
	using namespace core::id;
	using namespace core::scoring;
	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, BOGUS_ATOM_ID );
	for ( Size ii=0; ii<=end_overlap_pose-start_overlap_pose; ++ii ) {
		core::id::AtomID const id1(pose.residue(start_overlap_pose+ii).atom_index("CA"),start_overlap_pose+ii);
		core::id::AtomID const id2(parent_pose.residue(start_overlap_parent+ii).atom_index("CA"), start_overlap_parent+ii );
		atom_map[id1]=id2;
	}
	superimpose_pose(pose,parent_pose,atom_map);
	//create subpose of pose
	utility::vector1<Size> pose_positions;
	utility::vector1<Size> parent_posepositions;
	if ( overlap_location_pose== "n_term" ) {
		for ( Size ii=end_overlap_pose+1; ii<=pose.total_residue(); ++ii ) {
			pose_positions.push_back(ii);
		}
		for ( Size ii=1; ii<=end_overlap_parent; ++ii ) {
			parent_posepositions.push_back(ii);
		}
	}
	if ( overlap_location_pose== "c_term" ) {
		for ( Size ii=1; ii<start_overlap_pose; ++ii ) {
			pose_positions.push_back(ii);
		}
		for ( Size ii=start_overlap_parent; ii<=parent_pose.total_residue(); ++ii ) {
			parent_posepositions.push_back(ii);
		}
	}
	pose::Pose ref_pose_slice;
	pose::Pose parent_poseslice;

	pdbslice(ref_pose_slice,pose,pose_positions);
	pdbslice(parent_poseslice,parent_pose,parent_posepositions);
	if ( overlap_location_pose== "n_term" ) {
		remove_upper_terminus_type_from_pose_residue(parent_poseslice,parent_poseslice.total_residue());
		remove_lower_terminus_type_from_pose_residue(ref_pose_slice,1);
		append_pose_to_pose(parent_poseslice,ref_pose_slice,false);
		pose = parent_poseslice;
	}
	if ( overlap_location_pose== "c_term" ) {
		remove_upper_terminus_type_from_pose_residue(ref_pose_slice,ref_pose_slice.total_residue());
		remove_lower_terminus_type_from_pose_residue(parent_poseslice,1);
		append_pose_to_pose(ref_pose_slice,parent_poseslice,false);
		pose = ref_pose_slice;
	}
	renumber_pdbinfo_based_on_conf_chains(pose);
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
	maintain_cap_ = tag->getOption<bool>("maintain_cap",false);
	nTerm_cap_size_=0;
	cTerm_cap_size_=0;
	if ( maintain_cap_ ) {
		if ( (!tag->hasOption("nTerm_cap_size"))||(!tag->hasOption("cTerm_cap_size")) ) {
			throw utility::excn::EXCN_RosettaScriptsOption("repeat mover requires cTerm_cap and nTerm_cap defined if trying to maintain cap sequence or structure");
		}
		nTerm_cap_size_ = tag->getOption<Size>("nTerm_cap_size");
		cTerm_cap_size_ = tag->getOption<Size>("cTerm_cap_size");
	}
}


} // moves
} // protocols
