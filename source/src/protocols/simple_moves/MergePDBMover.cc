// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/MergePDBMover.cc
/// @brief This class will allign & combine parts of the pdb.
/// @author TJ Brunette (tjbrunette@gmail.com)
///

// Unit headers
#include <protocols/simple_moves/MergePDBMoverCreator.hh>
#include <protocols/simple_moves/MergePDBMover.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.MergePDBMover" );

#include <utility/tag/Tag.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/types.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Conformation.hh>

#include <core/select/util.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>

#include <core/import_pose/import_pose.hh>
#include <core/scoring/rms_util.hh>
#include <numeric/xyzVector.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

using namespace core;
using core::pose::Pose;

std::string MergePDBMoverCreator::keyname() const
{
	return MergePDBMoverCreator::mover_name();
}

protocols::moves::MoverOP
MergePDBMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MergePDBMover );
}

std::string
MergePDBMoverCreator::mover_name()
{
	return "MergePDB";
}

MergePDBMover::MergePDBMover()
: moves::Mover("MergePDB"),
	xml_input_pose_( /* NULL */ )
{
}

std::string
MergePDBMover::get_name() const {
	return MergePDBMoverCreator::mover_name();
}

moves::MoverOP
MergePDBMover::clone() const
{
	return moves::MoverOP( new MergePDBMover( *this ) );
}

moves::MoverOP
MergePDBMover::fresh_instance() const
{
	return moves::MoverOP( new MergePDBMover );
}

utility::vector1<MergePDBMover::Overlap> MergePDBMover::determine_overlap(Pose const pose){
	vector1<MergePDBMover::Overlap> hits;
	using namespace core::scoring;
	if ( overlap_location_pose_== "c_term" ) {
		Size initial_start_xmlPose=1;
		Size initial_end_xmlPose=initial_start_xmlPose+overlap_length_-1;
		Size initial_start_pose=pose.total_residue()-overlap_length_+1;
		Size initial_end_pose=pose.total_residue();
		for ( Size ii=0; ii<overlap_scan_range_; ++ii ) {
			utility::vector1<Size> xml_input_pose_positions;
			utility::vector1<Size> pose_positions;
			Size start_res_pose = initial_start_pose-ii;
			Size end_res_pose = initial_end_pose-ii;
			for ( Size jj=initial_start_xmlPose; jj<=initial_end_xmlPose; ++jj ) {
				xml_input_pose_positions.push_back(jj);
			}
			for ( Size jj=start_res_pose; jj<=end_res_pose; ++jj ) {
				pose_positions.push_back(jj);
			}
			pose::Pose xml_input_pose_slice;
			pose::Pose ref_pose_slice;
			pdbslice(ref_pose_slice,pose,pose_positions);
			pdbslice(xml_input_pose_slice,*xml_input_pose_,xml_input_pose_positions);
			Real rmsd = CA_rmsd(ref_pose_slice,xml_input_pose_slice);
			//std::cout << "rmsd" << rmsd << std::endl;
			TR << "rmsd" << rmsd << "kinase" << initial_start_xmlPose << ",pdb " <<start_res_pose << std::endl;
			if ( rmsd < overlap_max_rmsd_ ) {
				MergePDBMover::Overlap overlap_tmp(initial_start_xmlPose,initial_end_xmlPose,start_res_pose,end_res_pose);
				hits.push_back(overlap_tmp);
			}
		}
	}
	if ( overlap_location_pose_== "n_term" ) {
		Size initial_start_xmlPose=xml_input_pose_->total_residue()-overlap_length_+1;
		Size initial_end_xmlPose=xml_input_pose_->total_residue();
		Size initial_start_pose=1;
		Size initial_end_pose=1+overlap_length_-1;
		for ( Size ii=0; ii<overlap_scan_range_; ++ii ) {
			utility::vector1<Size> xml_input_pose_positions;
			utility::vector1<Size> pose_positions;
			Size start_res_pose = initial_start_pose+ii;
			Size end_res_pose = initial_end_pose+ii;
			for ( Size jj=initial_start_xmlPose; jj<=initial_end_xmlPose; ++jj ) {
				xml_input_pose_positions.push_back(jj);
			}
			for ( Size jj=start_res_pose; jj<=end_res_pose; ++jj ) {
				pose_positions.push_back(jj);
			}
			pose::Pose xml_input_pose_slice;
			pose::Pose ref_pose_slice;
			pdbslice(xml_input_pose_slice,*xml_input_pose_,xml_input_pose_positions);
			pdbslice(ref_pose_slice,pose,pose_positions);
			Real rmsd = CA_rmsd(ref_pose_slice,xml_input_pose_slice);
			TR << "rmsd" << rmsd << "kinase" << initial_start_xmlPose << ",pdb " <<start_res_pose << std::endl;
			if ( rmsd < overlap_max_rmsd_ ) {
				MergePDBMover::Overlap overlap_tmp(initial_start_xmlPose,initial_end_xmlPose,start_res_pose,end_res_pose);
				hits.push_back(overlap_tmp);
			}
		}
	}
	return(hits);
}

vector1<core::pose::PoseOP> MergePDBMover::generate_overlaps(Pose & pose, vector1<MergePDBMover::Overlap> overlaps) {
	using namespace core::id;
	using namespace core::scoring;
	vector1<core::pose::PoseOP> outputPoses;
	for ( Size kk=1; kk<=overlaps.size(); ++kk ) {
		Size start_overlap_pose = overlaps[kk].start_overlap_pose;
		Size end_overlap_pose = overlaps[kk].end_overlap_pose;
		Size start_overlap_xmlPose = overlaps[kk].start_overlap_xmlPose;
		Size end_overlap_xmlPose = overlaps[kk].end_overlap_xmlPose;
		core::id::AtomID_Map< core::id::AtomID > atom_map;
		core::pose::initialize_atomid_map( atom_map, pose, BOGUS_ATOM_ID );
		for ( Size ii=0; ii<=end_overlap_pose-start_overlap_pose; ++ii ) {
			core::id::AtomID const id1(pose.residue(start_overlap_pose+ii).atom_index("CA"),start_overlap_pose+ii);
			core::id::AtomID const id2(xml_input_pose_->residue(start_overlap_xmlPose+ii).atom_index("CA"), start_overlap_xmlPose+ii );
			atom_map[id1]=id2;
		}
		superimpose_pose(pose,*xml_input_pose_,atom_map);
		//create subpose of pose
		utility::vector1<Size> pose_positions;
		utility::vector1<Size> xml_input_pose_positions;
		if ( overlap_location_pose_== "n_term" ) {
			for ( Size ii=end_overlap_pose+1; ii<=pose.total_residue(); ++ii ) {
				pose_positions.push_back(ii);
			}
			for ( Size ii=1; ii<=end_overlap_xmlPose; ++ii ) {
				xml_input_pose_positions.push_back(ii);
			}
		}
		if ( overlap_location_pose_== "c_term" ) {
			for ( Size ii=1; ii<start_overlap_pose; ++ii ) {
				pose_positions.push_back(ii);
			}
			for ( Size ii=start_overlap_xmlPose; ii<=xml_input_pose_->total_residue(); ++ii ) {
				xml_input_pose_positions.push_back(ii);
			}
		}
		pose::PoseOP ref_pose_slice( new pose::Pose );
		pose::PoseOP xml_input_pose_slice( new pose::Pose );
		pdbslice(*ref_pose_slice,pose,pose_positions);
		pdbslice(*xml_input_pose_slice,*xml_input_pose_,xml_input_pose_positions);
		ref_pose_slice->conformation().detect_disulfides();
		xml_input_pose_slice->conformation().detect_disulfides();
		if ( overlap_location_pose_== "n_term" ) {
			remove_upper_terminus_type_from_pose_residue(*xml_input_pose_slice,xml_input_pose_slice->total_residue());
			remove_lower_terminus_type_from_pose_residue(*ref_pose_slice,1);
			for ( Size ii=1; ii<=ref_pose_slice->total_residue(); ++ii ) {
				xml_input_pose_slice->append_residue_by_bond(ref_pose_slice->residue(ii),false/*ideal bonds*/);
			}
			//append_pose_to_pose(*xml_input_pose_slice,*ref_pose_slice,false);
			renumber_pdbinfo_based_on_conf_chains(*xml_input_pose_slice,true,false,false,false);
			outputPoses.push_back(xml_input_pose_slice);
		}
		if ( overlap_location_pose_== "c_term" ) {
			remove_upper_terminus_type_from_pose_residue(*ref_pose_slice,ref_pose_slice->total_residue());
			remove_lower_terminus_type_from_pose_residue(*xml_input_pose_slice,1);
			for ( Size ii=1; ii<=xml_input_pose_slice->total_residue(); ++ii ) {
				ref_pose_slice->append_residue_by_bond(xml_input_pose_slice->residue(ii),false/*ideal bonds*/);
			}
			//append_pose_to_pose(*ref_pose_slice,*xml_input_pose_slice,false);
			renumber_pdbinfo_based_on_conf_chains(*ref_pose_slice,true,false,false,false);
			outputPoses.push_back(ref_pose_slice);
		}
	}
	return(outputPoses);
}

vector1<core::pose::PoseOP> MergePDBMover::pack_and_minimize(vector1<core::pose::PoseOP> poses, vector1<MergePDBMover::Overlap> overlaps,Real baseline_score) {
	using namespace core::pack::task;
	using namespace core::scoring;
	using namespace core::pack::task::operation;
	vector1<core::pose::PoseOP> outputPoses;
	for ( Size kk=1; kk<=overlaps.size(); ++kk ) {
		poses[kk]->conformation().detect_disulfides();
		sfxn_->score(*poses[kk]);
		Size start_overlap_pose = overlaps[kk].start_overlap_pose;
		Size end_overlap_pose = overlaps[kk].end_overlap_pose;
		Size start_overlap_xmlPose = overlaps[kk].start_overlap_xmlPose;
		Size end_overlap_xmlPose = overlaps[kk].end_overlap_xmlPose;
		vector1<Size> overlap_positions;
		utility::vector1< bool > mutant_resnums( poses[kk]->total_residue(), false);
		if ( overlap_location_pose_== "c_term" ) {
			for ( Size ii=start_overlap_pose; ii<=end_overlap_pose; ++ii ) {
				mutant_resnums[ii]=true;
			}
		}
		if ( overlap_location_pose_== "n_term" ) {
			for ( Size ii=start_overlap_xmlPose; ii<=end_overlap_xmlPose; ++ii ) {
				mutant_resnums[ii]=true;
			}
		}
		core::pack::task::operation::PreventRepacking turn_off_packing;
		core::pack::task::operation::RestrictResidueToRepacking turn_off_design;
		if ( design_range_>0 ) {
			utility::vector1< bool > mutant_resnums_and_neighbors = mutant_resnums;
			core::select::fill_neighbor_residues(*poses[kk], mutant_resnums_and_neighbors, design_range_);
			TR << "designed residues";
			for ( core::Size resnum = 1; resnum <= poses[kk]->total_residue(); ++resnum ) {
				if ( !mutant_resnums_and_neighbors[ resnum ] ) {
					turn_off_design.include_residue( resnum );
				}
				if ( mutant_resnums_and_neighbors[ resnum ] ) {
					TR << resnum <<",";
				}
			}
			TR << std::endl;
		}
		if ( packing_range_>0 ) {
			utility::vector1< bool > mutant_resnums_and_neighbors = mutant_resnums;
			core::select::fill_neighbor_residues(*poses[kk], mutant_resnums_and_neighbors, packing_range_);
			TR << "packed residues";
			for ( core::Size resnum = 1; resnum <= poses[kk]->total_residue(); ++resnum ) {
				if ( !mutant_resnums_and_neighbors[ resnum ] ) {
					turn_off_packing.include_residue( resnum );
				}
				if ( mutant_resnums_and_neighbors[ resnum ] ) {
					TR << resnum <<",";
				}
			}
			TR << std::endl;
		}
		PackerTaskOP task = task_factory_->create_task_and_apply_taskoperations(*poses[kk]);
		turn_off_design.apply(*poses[kk], *task);
		turn_off_packing.apply(*poses[kk], *task);
		Size pack_rounds = 2;
		protocols::simple_moves::PackRotamersMover packer = PackRotamersMover(sfxn_, task, pack_rounds);
		packer.apply(*poses[kk]);
		if ( do_minimize_ ) {
			optimization::MinimizerOptions minopt( "lbfgs_armijo_nonmonotone", 0.02, true, false, false );
			minopt.max_iter( 25 );
			optimization::AtomTreeMinimizer minimizer;
			kinematics::MoveMap mm_loc;
			mm_loc.set_jump( false ); mm_loc.set_bb( false ); mm_loc.set_chi( true );
			minimizer.run( *poses[kk], mm_loc, *sfxn_, minopt );
		}

		Real score = sfxn_->score(*poses[kk]);
		TR << "score after pack and min" << score << "baseline_score" << baseline_score << std::endl;
		if ( score<baseline_score ) {
			outputPoses.push_back(poses[kk]);
		}
	}
	return(outputPoses);
}


core::pose::PoseOP MergePDBMover::get_additional_output(){
	for ( Size ii=1; ii<=outputYet_.size(); ++ii ) {
		if ( outputYet_[ii]==false ) {
			set_last_move_status(protocols::moves::MS_SUCCESS);
			outputYet_[ii]=true;
			return(outputPoses_[ii]);
		}
	}
	set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
	return NULL;
}


void MergePDBMover::apply( Pose & pose )
{
	//bet baseline score------------
	pose.conformation().detect_disulfides();
	xml_input_pose_->conformation().detect_disulfides();
	core::Real baseline_score = sfxn_->score(pose);
	core::Real xml_score = sfxn_->score(*xml_input_pose_);
	if ( baseline_score<xml_score ) {
		baseline_score = xml_score;
	}

	//Step 1 get location
	vector1<Overlap> overlaps = determine_overlap(pose);
	//Step 2 combine
	vector1<core::pose::PoseOP> combined_poses = generate_overlaps(pose, overlaps);
	outputPoses_ = pack_and_minimize(combined_poses,overlaps,baseline_score);

	for ( Size ii=1; ii<=outputPoses_.size(); ++ii ) {
		outputYet_.push_back(false);
	}

	core::pose::PoseOP tmpPoseOP=get_additional_output();
	if ( tmpPoseOP!=NULL ) {
		pose=*tmpPoseOP;
	}
}


void MergePDBMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	overlap_location_pose_ = tag->getOption< std::string >( "attachment_termini" ,"n_term" );
	overlap_length_ = tag->getOption< core::Size >( "overlap_length", 4);
	overlap_max_rmsd_ = tag->getOption< Real >( "overlap_rmsd", 0.3);
	overlap_scan_range_ = tag->getOption< core::Size >( "overlap_scan_range", 20);
	std::string fname( tag->getOption< std::string >( "attach_pdb" ) );
	design_range_ = tag->getOption<Real>("design_range",3);
	packing_range_ = tag->getOption<Real>("packing_range",5);
	do_minimize_ = tag->getOption<bool>("do_minimize_",true);
	xml_input_pose_ = core::pose::PoseOP( new pose::Pose );
	core::import_pose::pose_from_file(*xml_input_pose_, fname , core::import_pose::PDB_file);
	if ( tag->hasOption("scorefxn") ) {
		std::string const scorefxn_key( tag->getOption<std::string>("scorefxn") );
		if ( datamap.has( "scorefxns", scorefxn_key ) ) {
			sfxn_ = datamap.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_key );
		} else {
			throw utility::excn::EXCN_RosettaScriptsOption("ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap.");
		}
	}
	if ( tag->hasOption("task_operations") ) {
		task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, datamap );
	} else {
		task_factory_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory );
	}

}

} // simple_moves
} // protocols

