// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_creation/MergePDBMover.cc
/// @brief This class will allign & combine parts of the pdb.
/// @author TJ Brunette (tjbrunette@gmail.com)
///

// Unit headers
#include <protocols/pose_creation/MergePDBMoverCreator.hh>
#include <protocols/pose_creation/MergePDBMover.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.pose_creation.MergePDBMover" );

#include <utility/tag/Tag.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/OrResidueSelector.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/InterGroupInterfaceByVectorSelector.hh>
#include <core/select/residue_selector/SymmetricalResidueSelector.hh> //need this?
//need this for symmetry interface detection
#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/AddPDBInfoMover.hh>

#include <core/select/residue_selector/util.hh>

#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/types.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/chains_util.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/selection.hh>
#include <core/pose/symmetry/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/select/util.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>

#include <core/import_pose/import_pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <numeric/xyzVector.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/simple_moves/CopyRotamerMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/AddResidueLabelMover.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace pose_creation {

using std::map;
using namespace core;
using core::pose::Pose;
using utility::vector1;
using namespace protocols::simple_moves;

MergePDBMover::MergePDBMover()
: moves::Mover("MergePDB"),
	xml_input_pose_( /* NULL */ )
{
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

void MergePDBMover::determine_overlap(Pose const pose,Size chain_id){
	utility::vector1<core::Size> res_in_chain = get_resnums_for_chain_id(pose,chain_id);
	Size start_res_chain = res_in_chain[1];
	Size end_res_chain = chain_end_res(pose,chain_id);
	using namespace core::scoring;
	//C-TERM
	if ( overlap_location_pose_== "c_term" ) {
		Size initial_start_xmlPose=1;
		Size initial_start_pose=end_res_chain-overlap_length_+1;
		Size initial_end_pose=end_res_chain;
		for ( Size ii=0; ii<overlap_scan_range_cmdLine_; ++ii ) {
			for ( Size kk=0; kk<overlap_scan_range_xml_; ++kk ) {
				Size start_res_xml_pose=initial_start_xmlPose+kk;
				Size end_res_xml_pose=start_res_xml_pose+overlap_length_-1;
				utility::vector1<Size> xml_input_pose_positions;
				utility::vector1<Size> pose_positions;
				Size start_res_pose = initial_start_pose-ii;
				Size end_res_pose = initial_end_pose-ii;
				for ( Size jj=start_res_xml_pose; jj<=end_res_xml_pose; ++jj ) {
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
				TR << "rmsd: " << rmsd << ", xml_pose_position: " << start_res_xml_pose << ", input_pose_position: " << start_res_pose << std::endl;
				if ( rmsd < overlap_max_rmsd_ ) {
					MergePDBMover::Overlap overlap_tmp(start_res_xml_pose,end_res_xml_pose,start_res_pose,end_res_pose);
					overlaps_.push_back(overlap_tmp);
				}
			}
		}
	}
	//N-TERM
	if ( overlap_location_pose_== "n_term" ) {
		Size initial_start_xmlPose=xml_input_pose_->total_residue()-overlap_length_+1;
		Size initial_start_pose=start_res_chain;
		Size initial_end_pose=start_res_chain+overlap_length_-1;
		for ( Size ii=0; ii<overlap_scan_range_cmdLine_; ++ii ) {
			for ( Size kk=0; kk<overlap_scan_range_xml_; ++kk ) {
				utility::vector1<Size> xml_input_pose_positions;
				utility::vector1<Size> pose_positions;
				Size start_res_pose = initial_start_pose+ii;
				Size end_res_pose = initial_end_pose+ii;
				Size start_res_xml_pose = initial_start_xmlPose-kk;
				Size end_res_xml_pose = start_res_xml_pose+overlap_length_-1;
				for ( Size jj=start_res_xml_pose; jj<=end_res_xml_pose; ++jj ) {
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
				TR << "rmsd: " << rmsd << ", xml_pose_position: " << start_res_xml_pose << ", input_pose_position: " << start_res_pose << std::endl;
				if ( rmsd < overlap_max_rmsd_ ) {
					MergePDBMover::Overlap overlap_tmp(start_res_xml_pose,end_res_xml_pose,start_res_pose,end_res_pose);
					overlaps_.push_back(overlap_tmp);
				}
			}
		}
	}
}

bool MergePDBMover::check_duplicate(Pose &pose){
	//setup to be fast because this could be a N^2 runtime.
	Real min_distance=999;
	Size min_position=999;
	for ( Size ii=1; ii<=overlaps_.size(); ++ii ) {
		if ( overlaps_[ii].output_poseOP != NULL ) {
			core::pose::PoseOP output_poseOP = overlaps_[ii].output_poseOP;
			if ( pose.total_residue() == output_poseOP->total_residue() ) {
				Size n_term_res = 1;
				Size c_term_res = pose.total_residue();
				while ( pose.residue(n_term_res).is_virtual_residue() )
						n_term_res++;
				while ( pose.residue(c_term_res).is_virtual_residue() )
						c_term_res--;
				Real n_term_dist = pose.residue(n_term_res).xyz("CA").distance(output_poseOP->residue(n_term_res).xyz("CA"));
				Real c_term_dist = pose.residue(c_term_res).xyz("CA").distance(output_poseOP->residue(c_term_res).xyz("CA"));
				Real dist = n_term_dist+c_term_dist;
				if ( dist<min_distance ) {
					min_distance=dist;
					min_position=ii;
				}
			}
		}
	}
	Real min_rmsd_dist = 999;
	if ( min_distance<3.0 ) { //only calc rmsd if needed
		min_rmsd_dist = core::scoring::CA_rmsd(pose,*(overlaps_[min_position].output_poseOP));
	}
	if ( min_rmsd_dist<duplicate_rmsd_pose_threshold_ ) {
		return true;
	}
	return(false);
}

void MergePDBMover::generate_overlaps(Pose & pose, Size chain_id) {
	using namespace core::id;
	using namespace core::scoring;
	vector1 < std::string > initial_labels;
	protocols::moves::MoverOP tocen = protocols::moves::MoverOP(
		new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
	utility::vector1<core::Size> res_in_chain = get_resnums_for_chain_id(pose,chain_id);
	bool input_pose_symmetric = pose::symmetry::is_symmetric(pose);
	core::scoring::ScoreFunctionOP sfxn0 = scoring::ScoreFunctionFactory::create_score_function( "score0" );
	if ( input_pose_symmetric ) {
		core::scoring::symmetry::symmetrize_scorefunction(*sfxn0);
	}
	Size start_res_chain = res_in_chain[1];
	Size end_res_chain = chain_end_res(pose,chain_id);
	//create a cloned version of xml_pose and xml_pose for rotamer copying.
	for ( Size kk=1; kk<=overlaps_.size(); ++kk ) {
		TR << "Executing overlap #" << kk << std::endl;
		vector1<Size> chunk1_insert_res;
		vector1<Size> chunk2_insert_res;
		Size start_overlap_pose = overlaps_[kk].start_overlap_pose;
		Size end_overlap_pose = overlaps_[kk].end_overlap_pose;
		Size start_overlap_xmlPose = overlaps_[kk].start_overlap_xmlPose;
		Size end_overlap_xmlPose = overlaps_[kk].end_overlap_xmlPose;
		core::id::AtomID_Map< core::id::AtomID > atom_map;
		core::pose::initialize_atomid_map( atom_map, *xml_input_pose_, AtomID::BOGUS_ATOM_ID() );
		for ( Size ii=0; ii<=end_overlap_pose-start_overlap_pose; ++ii ) {
			core::id::AtomID const id2(pose.residue(start_overlap_pose+ii).atom_index("CA"),start_overlap_pose+ii);
			core::id::AtomID const id1(xml_input_pose_->residue(start_overlap_xmlPose+ii).atom_index("CA"), start_overlap_xmlPose+ii );
			atom_map[id1]=id2;
		}
		//remove_lower_terminus_type_from_pose_residue(*xml_input_pose_slice,1);
		superimpose_pose(*xml_input_pose_,pose,atom_map);
		//create subpose of xml_pose
		core::pose::PoseOP output_poseOP( new core::pose::Pose() );
		core::pose::PoseOP tmp_xml_poseOP;
		tmp_xml_poseOP=xml_input_pose_->clone();
		Size end_steal_label_res = pose.size();
		TR.Debug << "pose.size() = " << pose.size() << std::endl;
		if ( input_pose_symmetric ) {
			pose::symmetry::extract_asymmetric_unit(pose,*output_poseOP,false,false);
			end_steal_label_res=end_res_chain;
		} else {
			output_poseOP->detached_copy(pose);
		}
		core::pose::PoseOP backup_input_pose( new core::pose::Pose() );
		core::pose::PoseOP backup_xml_pose( new core::pose::Pose() );
		if ( input_pose_symmetric ) {
			pose::symmetry::extract_asymmetric_unit(pose,*backup_input_pose,false,false);
		} else {
			backup_input_pose = pose.clone();
		}
		backup_xml_pose = xml_input_pose_->clone();

		//remove terminus_type from all backup chains
		vector1< core::Size >backup_chains = get_chains(*backup_input_pose);
		utility::vector1<core::Size> res_in_chain_backup;
		Size start_res_chain_backup;
		Size end_res_chain_backup;
		for ( core::Size cc = 1; cc <= backup_chains.size(); ++cc ) {
			TR << "Removing terminus_type from backup_input_pose chain: " << cc << std::endl;
			res_in_chain_backup = get_resnums_for_chain_id(*backup_input_pose,cc);
			start_res_chain_backup = res_in_chain[1];
			end_res_chain_backup = res_in_chain.back();
			remove_upper_terminus_type_from_pose_residue(*backup_input_pose,end_res_chain_backup);
			remove_lower_terminus_type_from_pose_residue(*backup_input_pose,start_res_chain_backup);
		}
		TR << "Removing terminus_type from backup_xml_pose" << std::endl;
		remove_upper_terminus_type_from_pose_residue(*backup_xml_pose,backup_xml_pose->total_residue());
		remove_lower_terminus_type_from_pose_residue(*backup_xml_pose,1);

		//N-term
		if ( overlap_location_pose_== "n_term" ) {
			TR.Debug << "entering n_term mode" << std::endl;
			TR << "Removing terminus_type from tmp_xml_pose and output_pose" << std::endl;
			remove_upper_terminus_type_from_pose_residue(*tmp_xml_poseOP,tmp_xml_poseOP->total_residue()); //need?
			remove_lower_terminus_type_from_pose_residue(*output_poseOP,start_res_chain); //need? always deleted anyway
			//delete any additional residues
			TR.Debug << "deleting residues: " << start_res_chain << "-" << start_overlap_pose << std::endl;
			output_poseOP->conformation().delete_residue_range_slow(start_res_chain,start_overlap_pose); //always delete the first residue because of the lost phi/psi
			for ( Size ii=start_overlap_xmlPose; ii>=1; --ii ) {
				output_poseOP->conformation().safely_prepend_polymer_residue_before_seqpos(tmp_xml_poseOP->residue(ii),start_res_chain, false); //adds residues from tmp_xml_pose to output_pose
			}
			//if NOT chain A, need to add back initial residue labels of residues in chains BEFORE the merged chain (N-TERM)
			Size chain_offset = 0;
			if ( chain_id != 1 ) { //skips if chain_id = 1
				utility::vector1<core::Size> res_in_chain_before;
				Size start_res_chain_before;
				Size end_res_chain_before;
				for ( core::Size cc = 1; cc < chain_id; ++cc ) { //loop through all chains before chain_id
					TR << "Adding missed initial_interface residue_labels (n_term mode) for chain " << cc << std::endl;
					//need residue numbers for each chain
					res_in_chain_before = get_resnums_for_chain_id(pose,cc);
					start_res_chain_before = res_in_chain_before[1];
					end_res_chain_before = res_in_chain_before.back();
					TR.Debug << "start_res_chain_before: " << start_res_chain_before << ", end_res_chain_before: " << end_res_chain_before << std::endl;
					//loop through each residue in the chain
					for ( core::Size ii = start_res_chain_before; ii <= end_res_chain_before; ++ii ) {
						initial_labels = pose.pdb_info()->get_reslabels(ii);
						//loop through each residue's labels
						for ( core::Size jj=1; jj <= initial_labels.size(); ++jj ) {
							output_poseOP->pdb_info()->add_reslabel(ii, initial_labels[jj]);
						}
					}
					chain_offset += res_in_chain_before.size(); //determine chain_offset (if not chain 1)
				}
			}
			TR << "Total chain_offset: " << chain_offset << std::endl;
			//re-initialize start/end positions
			Size orig_start_overlap_xmlPose = start_overlap_xmlPose; //need this for copy_sequence
			start_overlap_xmlPose += chain_offset;
			end_overlap_xmlPose += chain_offset;
			TR.Debug << "start_overlap_pose: " << start_overlap_pose << ", end_overlap_pose: " << end_overlap_pose << ", start_overlap_xmlPose: " << start_overlap_xmlPose << ", end_overlap_xmlPose: " << end_overlap_xmlPose << std::endl;
			//add initial residue labels
			TR << "Adding initial_interface residue_labels (n_term mode)" << std::endl;
			Size res_loss_n_term = start_overlap_pose-start_res_chain;
			Size res_add = orig_start_overlap_xmlPose;
			Size res_offset = res_add-res_loss_n_term-1;
			TR.Debug << "start_overlap_pose: " << start_overlap_pose << ", end_steal_label_res: " << end_steal_label_res << ", res_offset: " << res_offset << std::endl;
			for ( core::Size ii = start_overlap_pose; ii <= end_steal_label_res; ++ii ) { //chose to do this for the output_pose_size. Should take care of symmetric case
				initial_labels = pose.pdb_info()->get_reslabels(ii);
				for ( core::Size jj=1; jj <= initial_labels.size(); ++jj ) {
					output_poseOP->pdb_info()->add_reslabel(res_offset+ii, initial_labels[jj]);
				}
			}
			//add overlap residue label
			TR << "Adding overlap residue_labels (n_term mode)" << std::endl;
			TR.Debug << "sele overlap, resi " << start_overlap_xmlPose << "-" << end_overlap_xmlPose << std::endl;
			for ( core::Size ii=start_overlap_xmlPose; ii<=end_overlap_xmlPose; ++ii ) {
				output_poseOP->pdb_info()->add_reslabel(ii, "overlap");
			}
			//create residues for residue_selectors
			chunk1_insert_res.clear();
			chunk2_insert_res.clear();
			for ( core::Size ii=chain_offset+1; ii<=end_overlap_xmlPose; ++ii ) {
				chunk2_insert_res.push_back(ii);
			}
			for ( core::Size ii=1; ii<=chain_offset; ++ii ) {
				chunk1_insert_res.push_back(ii);
			}
			for ( core::Size ii=start_overlap_xmlPose; ii<=output_poseOP->total_residue(); ++ii ) {
				chunk1_insert_res.push_back(ii);
			}
			//copy sequence for overlap
			TR << "Copying sequence" << std::endl;
			copy_sequence(start_overlap_xmlPose,end_overlap_xmlPose,start_overlap_pose,orig_start_overlap_xmlPose,*backup_input_pose,*backup_xml_pose,*output_poseOP);
		}

		//C-term
		if ( overlap_location_pose_== "c_term" ) {
			TR.Debug << "entering c_term mode" << std::endl;
			TR << "Removing terminus_type from tmp_xml_pose and output_pose" << std::endl;
			remove_lower_terminus_type_from_pose_residue(*tmp_xml_poseOP,1);
			remove_upper_terminus_type_from_pose_residue(*output_poseOP,end_res_chain);
			//delete any additional residues
			TR.Debug << "deleting residues: " << end_overlap_pose << "-" << end_res_chain << std::endl;
			output_poseOP->conformation().delete_residue_range_slow(end_overlap_pose,end_res_chain);
			for ( Size ii=end_overlap_xmlPose; ii<=tmp_xml_poseOP->total_residue(); ++ii ) {
				Size location_to_insert=end_overlap_pose-1+ii-end_overlap_xmlPose;
				output_poseOP->conformation().safely_append_polymer_residue_after_seqpos(tmp_xml_poseOP->residue(ii),location_to_insert, false); //adds residues from tmp_xml_pose to output_pose
			}
			//add initial residue labels
			TR << "Adding initial_interface residue_labels (c_term mode)" << std::endl;
			TR.Debug << "START: " << 1 << ", end_overlap_pose: " << end_overlap_pose << std::endl;
			for ( core::Size ii = 1; ii <= end_overlap_pose; ++ii ) {
				initial_labels =  pose.pdb_info()->get_reslabels(ii);
				for ( core::Size jj=1; jj <= initial_labels.size(); ++jj ) {
					output_poseOP->pdb_info()->add_reslabel(ii, initial_labels[jj]);
				}
			}
			//need this for adding back missed residues AND for making residue selectors
			Size chain_offset = get_resnums_for_chain_id(*output_poseOP,chain_id).back();
			TR.Debug << "chain_offset: " << chain_offset << std::endl;
			//if NOT last chain, need to add back initial residue labels of residues in chains AFTER the merged chain (C-TERM)
			vector1< core::Size >chains = get_chains(pose);
			core::Size last_chain_in_pose = chains.back();
			if ( ( chain_id != last_chain_in_pose ) && ( ! input_pose_symmetric ) ) { //skips if chain_id = last chain AND if pose is not symmetric
				utility::vector1<core::Size> res_in_chain_after;
				Size start_res_chain_after;
				Size end_res_chain_after;
				for ( core::Size cc = chain_id+1; cc <= last_chain_in_pose; ++cc ) { //loop through all chains AFTER chain_id until last_chain_in_pose
					TR << "Adding missed initial_interface residue_labels (c_term mode) for chain " << cc << std::endl;
					//need residue numbers for each chain
					res_in_chain_after = get_resnums_for_chain_id(pose,cc);
					start_res_chain_after = res_in_chain_after[1];
					end_res_chain_after = res_in_chain_after.back();
					TR.Debug << "start_res_chain_after: " << start_res_chain_after << ", end_res_chain_after: " << end_res_chain_after << std::endl;
					//loop through each residue in the chain
					for ( core::Size ii = start_res_chain_after; ii <= end_res_chain_after; ++ii ) {
						initial_labels =  pose.pdb_info()->get_reslabels(ii);
						Size pose_end_of_chain_id = get_resnums_for_chain_id(pose,chain_id).back();
						//loop through each residue's labels
						for ( core::Size jj=1; jj <= initial_labels.size(); ++jj ) {
							output_poseOP->pdb_info()->add_reslabel(ii+chain_offset-pose_end_of_chain_id, initial_labels[jj]);
						}
					}
				}
			}
			//add overlap residue labels
			TR << "Adding overlap residue_labels (c_term mode)" << std::endl;
			TR.Debug << "start_overlap_pose: " << start_overlap_pose << ", end_overlap_pose: " << end_overlap_pose << std::endl;
			for ( core::Size ii=start_overlap_pose; ii<=end_overlap_pose; ++ii ) {
				output_poseOP->pdb_info()->add_reslabel(ii, "overlap");
			}
			//create residues for residue_selectors
			chunk1_insert_res.clear();
			chunk2_insert_res.clear();
			TR.Debug << "select chunk1, resi ";
			for ( core::Size ii=1; ii<=end_overlap_pose; ++ii ) { //chunk1 residues before fusion
				chunk1_insert_res.push_back(ii);
				TR.Debug << ii << "+";
			}
			for ( core::Size ii=chain_offset+1; ii<=output_poseOP->total_residue(); ++ii ) { //chunk1 residues after fusion
				chunk1_insert_res.push_back(ii);
				TR.Debug << ii << "+";
			} TR.Debug << std::endl;
			TR.Debug << "select chunk2, resi ";
			for ( core::Size ii=start_overlap_pose; ii<=chain_offset; ++ii ) { //chunk2 residues
				chunk2_insert_res.push_back(ii);
				TR.Debug << ii << "+";
			} TR.Debug << std::endl;
			//copy sequence for overlap
			TR << "Copying sequence" << std::endl;
			copy_sequence(start_overlap_pose,end_overlap_pose,start_overlap_pose,start_overlap_xmlPose,*backup_input_pose,*backup_xml_pose,*output_poseOP);
		}
		//seems really weird. When I created a brand new scorefunction it was symmetric
		TR << "Resetting new scorefunction" << std::endl;
		if ( input_pose_symmetric ) {
			asymm_score_->score(*output_poseOP); //needs to be here to calc distances
		} else {
			sfxn_->score(*output_poseOP);
		}

		//discard asymmetric poses with bad clashes
		//if ( ! input_pose_symmetric ) {
		core::pose::PoseOP centroidPoseOP = output_poseOP->clone();
		tocen->apply( *centroidPoseOP );
		Real score0 = sfxn0->score( *centroidPoseOP );
		TR << "ASYMMETRIC score0 for clash_check (threshold: " << clash_threshold_ << " );" << score0 << std::endl;
		bool duplicate = check_duplicate(*output_poseOP);
		if ( score0>clash_threshold_ || duplicate ) {
			//no output
			overlaps_[kk].output_poseOP = NULL;
			TR << "ASYMMETRIC clash_check failed OR duplicate structure, tossing structure #: " << kk << std::endl;
			continue;
		}
		//}

		//step 1 add location where chain A and B now overlap as overlap_2 (other_overlap)
		TR << "Determining other_overlap" << std::endl;
		core::select::residue_selector::ResidueSelectorCOP overlap_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector("overlap"));
		core::select::residue_selector::ResidueSelectorCOP interface_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector(no_design_label_));
		core::select::residue_selector::ResidueSelectorCOP overlap_plus_interface_selector(new core::select::residue_selector::OrResidueSelector(overlap_selector,interface_selector));
		core::select::residue_selector::ResidueSelectorCOP not_overlap_interface_selector(new core::select::residue_selector::NotResidueSelector(overlap_plus_interface_selector));
		core::select::residue_selector::ResidueSelectorCOP partA_selector(new core::select::residue_selector::ResidueIndexSelector(chunk1_insert_res)); //always input_pose
		core::select::residue_selector::ResidueSelectorCOP partB_selector(new core::select::residue_selector::ResidueIndexSelector(chunk2_insert_res)); //always xml_pose
		core::select::residue_selector::InterGroupInterfaceByVectorSelectorOP interface_AB_selector( new core::select::residue_selector::InterGroupInterfaceByVectorSelector );
		interface_AB_selector->group1_selector( partA_selector );
		interface_AB_selector->group2_selector( partB_selector );
		interface_AB_selector->cb_dist_cut( design_range_ ); //defaults to 11
		interface_AB_selector->vector_dist_cut( design_range_*0.8181 ); //defaults to 9 in the selector. This is really arbitrary atm, I just took the default value of 9/11 = 0.8181
		core::select::residue_selector::InterGroupInterfaceByVectorSelectorOP interface_BA_selector( new core::select::residue_selector::InterGroupInterfaceByVectorSelector ); //reverse since the selector has to be downstream
		interface_BA_selector->group1_selector( partB_selector );
		interface_BA_selector->group2_selector( partA_selector );
		interface_BA_selector->cb_dist_cut( design_range_ ); //defaults to 11
		interface_BA_selector->vector_dist_cut( design_range_*0.8181 );
		core::select::residue_selector::ResidueSelectorCOP interface_ABBA_selector(new core::select::residue_selector::OrResidueSelector(interface_AB_selector,interface_BA_selector));
		core::select::residue_selector::ResidueSelectorCOP interface_exclude_overlap_selector(new core::select::residue_selector::AndResidueSelector(interface_ABBA_selector,not_overlap_interface_selector));
		protocols::simple_moves::AddResidueLabelMover add_other_overlap = AddResidueLabelMover(interface_exclude_overlap_selector,"other_overlap");
		TR.Debug << "apply other_overlap to output_pose" << std::endl;
		add_other_overlap.apply(*output_poseOP);


		//store PDBInfoLabels
		std::map<core::Size, vector1<std::string> > res_label_map;
		for ( core::Size ii = 1; ii <=output_poseOP->size(); ++ii ) {
			vector1<std::string> tmp_labels = output_poseOP->pdb_info()->get_reslabels(ii);
			res_label_map.insert(std::pair<Size,vector1<std::string> >(ii,tmp_labels));
		}
		//reset pdbinfos
		renumber_pdbinfo_based_on_conf_chains(*output_poseOP,true,false,false,false);

		//reapply residue_label to output_pose
		std::map<core::Size, vector1<std::string> >::iterator res_label_map_itr;
		vector1<std::string>::iterator res_label_itr;
		for ( res_label_map_itr=res_label_map.begin(); res_label_map_itr!=res_label_map.end(); ++res_label_map_itr ) {
			Size resid = res_label_map_itr->first;
			vector1<std::string> tmp_labels = res_label_map_itr->second;
			for ( res_label_itr=tmp_labels.begin(); res_label_itr!=tmp_labels.end(); ++res_label_itr ) {
				output_poseOP->pdb_info()->add_reslabel(resid, *res_label_itr);
			}
		}

		//store and copy back labels due to symmetry wiping it **symmetry cases only**
		if ( input_pose_symmetric ) {
			//clean and regenerate symmetry
			core::pose::symmetry::make_symmetric_pose(*output_poseOP,symm_file_);
			sfxn_->score(*output_poseOP);

			//discard SYMMETRIC poses with bad clashes
			core::pose::PoseOP centroidPoseOP = output_poseOP->clone();
			tocen->apply( *centroidPoseOP );
			Real score0 = sfxn0->score( *centroidPoseOP );
			TR << "SYMMETRIC score0 for clash_check (threshold: " << clash_threshold_ << " );" << score0 << std::endl;
			bool duplicate = check_duplicate(*output_poseOP);
			if ( score0>clash_threshold_ || duplicate ) {
				//no output
				overlaps_[kk].output_poseOP = NULL;
				TR << "SYMMETRIC clash_check failed OR duplicate structure, tossing structure #: " << kk << std::endl;
				continue;
			}

			//label where the chains hit-each due to symmetry
			TR << "Determining symmetric other_overlap_sym" << std::endl;

			//need to do some gymnastics since building_block_interface is not a residue_selector yet
			//clone output_pose
			core::pose::PoseOP asp_output_poseOP = output_poseOP->clone();
			//extract asymmetric pose, re-apply PDBInfo and sfxn
			protocols::symmetry::ExtractAsymmetricPoseMover extract_asp;
			extract_asp.clear_sym_def( true );
			extract_asp.apply( *asp_output_poseOP );
			protocols::simple_moves::AddPDBInfoMover add_pdb_info; //this thing is KEY.
			add_pdb_info.apply( *asp_output_poseOP );
			asymm_score_->score(*asp_output_poseOP);
			asp_output_poseOP->update_residue_neighbors();
			renumber_pdbinfo_based_on_conf_chains(*asp_output_poseOP,true,false,false,false);

			//reapply residue_label to ASP
			std::map<core::Size, vector1<std::string> >::iterator res_label_map_itr;
			vector1<std::string>::iterator res_label_itr;
			for ( res_label_map_itr=res_label_map.begin(); res_label_map_itr!=res_label_map.end(); ++res_label_map_itr ) {
				Size resid = res_label_map_itr->first;
				vector1<std::string> tmp_labels = res_label_map_itr->second;
				for ( res_label_itr=tmp_labels.begin(); res_label_itr!=tmp_labels.end(); ++res_label_itr ) {
					asp_output_poseOP->pdb_info()->add_reslabel(resid, *res_label_itr);
				}
			}

			//just prints res_id for each chain (for debugging)
			//for ( core::Size cc = 1; cc <= get_chains( *output_poseOP ).back(); ++cc ) {
			// utility::vector1<core::Size> res_in_chainA = get_resnums_for_chain_id(*asp_output_poseOP,cc);
			// std::cout << "output_pose asu chain " << cc << ", residues: " << res_in_chainA[1] << "-" << res_in_chainA.back() << std::endl;
			//}
			//for ( core::Size cc = 1; cc <= get_chains( *asp_output_poseOP ).back(); ++cc ) {
			// utility::vector1<core::Size> res_in_chainA_asp = get_resnums_for_chain_id(*asp_output_poseOP,cc);
			// std::cout << "output_pose_asp asu chain " << cc << ", residues: " << res_in_chainA_asp[1] << "-" << res_in_chainA_asp.back() << std::endl;
			//}

			//residues for chain1, chain1 selector doesn't work for some reason
			utility::vector1<core::Size> res_in_chainA_asp = get_resnums_for_chain_id(*asp_output_poseOP,1);

			//hack to get chainBsym_resis for asp
			core::Size asp_num_chains = get_chains( *asp_output_poseOP ).size()-1; //-1 for the stupid symmetry leftover residues
			vector1<Size> chunk2_insert_res_asp = chunk2_insert_res;
			core::Size chain_len = res_in_chainA_asp.size();
			for ( core::Size cc = 1; cc < asp_num_chains; ++cc ) {
				for ( core::Size ii = 1; ii <= chunk2_insert_res.size(); ++ii ) {
					chunk2_insert_res_asp.push_back( chunk2_insert_res[ii]+(cc*chain_len) );
				}
			}

			//for debugging, pymol selection for chunk2 in asp
			//TR << "select chunk2_insert_res_asp" << ", resi ";
			//for (core::Size ii=1; ii <= chunk2_insert_res_asp.size(); ++ii ){
			// TR << chunk2_insert_res_asp[ii] << "+";
			//} TR << std::endl;

			//add other_overlap_sym
			//select old residue_labels from asym step
			core::select::residue_selector::ResidueSelectorCOP other_overlap_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector("other_overlap"));
			core::select::residue_selector::ResidueSelectorCOP overlap_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector("overlap"));
			core::select::residue_selector::ResidueSelectorCOP interface_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector(no_design_label_));
			core::select::residue_selector::OrResidueSelectorOP old_labels_selector(new core::select::residue_selector::OrResidueSelector);
			old_labels_selector->add_residue_selector(other_overlap_selector);
			old_labels_selector->add_residue_selector(overlap_selector);
			old_labels_selector->add_residue_selector(interface_selector);
			core::select::residue_selector::ResidueSelectorCOP not_old_labels_selector(new core::select::residue_selector::NotResidueSelector(old_labels_selector));
			//select sym_interface
			core::select::residue_selector::ResidueSelectorCOP chain_selector(new core::select::residue_selector::ResidueIndexSelector( res_in_chainA_asp ));
			core::select::residue_selector::ResidueSelectorCOP not_chainA_selector(new core::select::residue_selector::NotResidueSelector(chain_selector));
			core::select::residue_selector::InterGroupInterfaceByVectorSelectorOP interface_sym_AB_selector( new core::select::residue_selector::InterGroupInterfaceByVectorSelector );
			interface_sym_AB_selector->group1_selector( chain_selector );
			interface_sym_AB_selector->group2_selector( not_chainA_selector );
			interface_sym_AB_selector->cb_dist_cut( design_range_ ); //defaults to 11
			interface_sym_AB_selector->vector_dist_cut( design_range_*0.8181 );
			//combine all selections
			core::select::residue_selector::ResidueSelectorCOP chunk2_asp_selector(new core::select::residue_selector::ResidueIndexSelector( chunk2_insert_res_asp ));
			core::select::residue_selector::ResidueSelectorCOP chunk2_around_selector(new core::select::residue_selector::NeighborhoodResidueSelector(chunk2_asp_selector,design_range_,true));
			core::select::residue_selector::AndResidueSelectorOP new_sym_interface_selector(new core::select::residue_selector::AndResidueSelector );
			new_sym_interface_selector->add_residue_selector( chain_selector ); //must be on chain1 (asu)
			new_sym_interface_selector->add_residue_selector( not_old_labels_selector ); //must NOT be overlap or other_overlap or no_design_label
			new_sym_interface_selector->add_residue_selector( interface_sym_AB_selector ); //symmetrical interface residue
			new_sym_interface_selector->add_residue_selector( chunk2_around_selector ); //must be near chunk2
			protocols::simple_moves::AddResidueLabelMover add_other_overlap_v2 = AddResidueLabelMover(new_sym_interface_selector,"other_overlap_sym");

			TR.Debug << "apply other_overlap_sym to sym_output_pose" << std::endl;
			add_other_overlap_v2.apply(*asp_output_poseOP);

			//store asp residue_labels
			TR.Debug << "store symmetric pose residue_labels" << std::endl;
			std::map<core::Size, vector1<std::string> > res_label_map_asp;
			for ( core::Size ii = 1; ii <=get_resnums_for_chain_id(*asp_output_poseOP,1).size(); ++ii ) { //only stores resis for chain1
				vector1<std::string> tmp_labels = asp_output_poseOP->pdb_info()->get_reslabels(ii);
				res_label_map_asp.insert(std::pair<Size,vector1<std::string> >(ii,tmp_labels));
			}
			//reapply asp residue_labels to output_pose
			std::map<core::Size, vector1<std::string> >::iterator res_label_map_itr_asp;
			vector1<std::string>::iterator res_label_itr_asp;
			for ( res_label_map_itr_asp=res_label_map_asp.begin(); res_label_map_itr_asp!=res_label_map_asp.end(); ++res_label_map_itr_asp ) {
				Size resid = res_label_map_itr_asp->first;
				vector1<std::string> tmp_labels = res_label_map_itr_asp->second;
				for ( res_label_itr_asp=tmp_labels.begin(); res_label_itr_asp!=tmp_labels.end(); ++res_label_itr_asp ) {
					output_poseOP->pdb_info()->add_reslabel(resid, *res_label_itr_asp);
				}
			}
		} //end symmetric case

		//accept poses, write output
		overlaps_[kk].output_poseOP = output_poseOP;
		//print pymol selections for labels
		TR << "PDBInfo labels for overlap #" << kk << std::endl;
		TR << "select " << no_design_label_ << ", resi ";
		for ( core::Size ii = 1; ii <= output_poseOP->total_residue(); ++ii ) {
			if ( output_poseOP->pdb_info()->res_haslabel(ii, no_design_label_ ) ) {
				TR << ii << "+";
			}
		} TR << std::endl;
		TR << "select overlap" << ", resi ";
		for ( core::Size ii = 1; ii <= output_poseOP->total_residue(); ++ii ) {
			if ( output_poseOP->pdb_info()->res_haslabel(ii, "overlap" ) ) {
				TR << ii << "+";
			}
		} TR << std::endl;
		TR << "select other_overlap" << ", resi ";
		for ( core::Size ii = 1; ii <= output_poseOP->total_residue(); ++ii ) {
			if ( output_poseOP->pdb_info()->res_haslabel(ii, "other_overlap" ) ) {
				TR << ii << "+";
			}
		} TR << std::endl;
		TR << "select other_overlap_sym" << ", resi ";
		for ( core::Size ii = 1; ii <= output_poseOP->total_residue(); ++ii ) {
			if ( output_poseOP->pdb_info()->res_haslabel(ii, "other_overlap_sym" ) ) {
				TR << ii << "+";
			}
		} TR << std::endl;
	} //end overlap (kk) loop
} //end generate_overlaps

Size MergePDBMover::closest_non_overlap_residue(core::pose::Pose const & pose, core::Size resid, core::Size start_overlap_resid, core::Size end_overlap_resid){
	Real min_dist = 999;
	Size closest_residue = 999999;
	for ( Size ii=1; ii<=pose.total_residue(); ++ii ) {
		if ( ii<start_overlap_resid || ii>end_overlap_resid ) {
			Real tmp_dist = pose.residue(ii).xyz("CA").distance(pose.residue(resid).xyz("CA"));
			if ( tmp_dist < min_dist ) {
				min_dist = tmp_dist;
				closest_residue = ii;
			}
		}
	}
	return(closest_residue);
}

/// @brief increases the range to ignore to include the entire secondary structure element
void MergePDBMover::increase_range_to_ignore_ss_element(core::pose::Pose const & pose, Size init_start, Size init_end, Size & ss_start, Size & ss_end){
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	char ss_type = dssp.get_dssp_secstruct(init_start);
	ss_start = init_start;
	while ( ss_type == dssp.get_dssp_secstruct(ss_start) && ss_start>1 ) {
		ss_start--;
	}
	ss_start++;
	ss_end=init_end;
	ss_type = dssp.get_dssp_secstruct( init_end);
	while ( ss_type == dssp.get_dssp_secstruct(ss_end) && ss_end<pose.total_residue() ) {
		ss_end++;
	}
	ss_end--;
}

void MergePDBMover::copy_sequence(core::Size start_overlap_resid, core::Size end_overlap_resid, core::Size start_overlap_input_pose_resid, core::Size start_overlap_xml_pose_resid, core::pose::Pose const & input_pose, core::pose::Pose const & xml_pose, core::pose::Pose & output_pose){
	protocols::simple_moves::CopyRotamerMover copyRotamerMover = CopyRotamerMover();
	Size numb_res = end_overlap_resid-start_overlap_resid+1;
	TR.Debug << "end_overlap_resid: " << end_overlap_resid << ", start_overlap_resid: " << start_overlap_resid << ", numb_res: " << numb_res << std::endl;
	if ( init_overlap_sequence_=="input_pose" ) {
		for ( Size ii=0; ii<numb_res; ++ii ) {
			TR.Debug << "copying residue " << start_overlap_input_pose_resid+ii << " (input_pose) to " << start_overlap_resid+ii << " (output_pose)" << std::endl;
			copyRotamerMover.apply_from_template_pose(output_pose,input_pose,start_overlap_resid+ii,start_overlap_input_pose_resid+ii);
		}
	}
	if ( init_overlap_sequence_=="xml_pose" ) {
		for ( Size ii=0; ii<numb_res; ++ii ) {
			copyRotamerMover.apply_from_template_pose(output_pose,xml_pose,start_overlap_resid+ii,start_overlap_xml_pose_resid+ii);
		}
	}
	if ( init_overlap_sequence_ =="both" ) {
		//ignores helix that the structure is on for proximity
		Size ss_start_resid=0;
		Size ss_end_resid=0;
		increase_range_to_ignore_ss_element(output_pose,start_overlap_resid,end_overlap_resid, ss_start_resid,ss_end_resid);
		for ( Size ii=0; ii<numb_res; ++ii ) {
			Size closest_resid = closest_non_overlap_residue(output_pose,start_overlap_resid+ii,ss_start_resid,ss_end_resid);
			if ( closest_resid<start_overlap_resid && overlap_location_pose_ == "n_term" ) {
				copyRotamerMover.apply_from_template_pose(output_pose,xml_pose,start_overlap_resid+ii,start_overlap_xml_pose_resid+ii);
			}
			if ( closest_resid<start_overlap_resid && overlap_location_pose_ == "c_term" ) {
				copyRotamerMover.apply_from_template_pose(output_pose,input_pose,start_overlap_resid+ii,start_overlap_input_pose_resid+ii);
			}
			if ( closest_resid>end_overlap_resid && overlap_location_pose_ == "n_term" ) {
				copyRotamerMover.apply_from_template_pose(output_pose,input_pose,start_overlap_resid+ii,start_overlap_input_pose_resid+ii);
			}
			if ( closest_resid>end_overlap_resid && overlap_location_pose_ == "c_term" ) {
				copyRotamerMover.apply_from_template_pose(output_pose,xml_pose,start_overlap_resid+ii,start_overlap_xml_pose_resid+ii);
			}
		}
	}
}

void MergePDBMover::pack_and_minimize(Pose const pose, core::Real baseline_score){
	using namespace core::pack::task;
	using namespace core::scoring;
	using namespace core::pack::task::operation;
	bool pose_symmetric = core::pose::symmetry::is_symmetric(pose);
	bool early_complete = false; //true only if output_only_first_ is set to true and one object has been successfully produced
	for ( Size kk=1; kk<=overlaps_.size(); ++kk ) {
		if ( early_complete ) {
			overlaps_[kk].output_poseOP=NULL; //get rid of all other poses
		}
		if ( overlaps_[kk].output_poseOP != NULL ) { //if there are no overlaps in input
			if ( detect_disulf_before_repack_ ) {
				overlaps_[kk].output_poseOP->conformation().detect_disulfides();
			}
			TaskFactoryOP taskfactoryOP = task_factory_->clone();
			//PackerTaskOP taskOP = task_factory_->create_task_and_apply_taskoperations(*poses[kk]);

			//find residues allowed to design
			core::select::residue_selector::ResidueSelectorCOP overlap_p1_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector("overlap"));
			core::select::residue_selector::ResidueSelectorCOP overlap_p2_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector("other_overlap"));
			core::select::residue_selector::ResidueSelectorCOP overlap_p3_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector("other_overlap_sym"));
			core::select::residue_selector::OrResidueSelectorOP overlap_selector(new core::select::residue_selector::OrResidueSelector);
			overlap_selector->add_residue_selector(overlap_p1_selector);
			overlap_selector->add_residue_selector(overlap_p2_selector);
			overlap_selector->add_residue_selector(overlap_p3_selector);
			core::select::residue_selector::ResidueSelectorCOP no_design_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector(no_design_label_));
			core::select::residue_selector::ResidueSelectorCOP not_no_design_selector(new core::select::residue_selector::NotResidueSelector(no_design_selector));
			core::select::residue_selector::ResidueSelectorCOP design_selector(new core::select::residue_selector::AndResidueSelector(overlap_selector,not_no_design_selector));
			//core::select::residue_selector::ResidueSelectorCOP symmetric_design_selector(new core::select::residue_selector::SymmetricalResidueSelector(design_selector)); //shouldn't need this since packer is sym aware
			utility::vector1< bool > design_res = design_selector->apply( *overlaps_[kk].output_poseOP );

			TR << "select design_res, resi ";
			for ( Size ii=1; ii<=design_res.size(); ++ii ) {
				if ( design_res[ii]==1 ) {
					TR << ii << "+";
				}
			}
			TR << std::endl;

			//selecting residues allowed to pack
			core::select::residue_selector::ResidueSelectorCOP near_overlap_selector(new core::select::residue_selector::NeighborhoodResidueSelector(overlap_selector,packing_range_,false));
			core::select::residue_selector::ResidueSelectorCOP pack_selector(new core::select::residue_selector::AndResidueSelector(near_overlap_selector,not_no_design_selector));
			//core::select::residue_selector::ResidueSelectorCOP symmetric_pack_selector(new core::select::residue_selector::SymmetricalResidueSelector(pack_selector));
			utility::vector1< bool > pack_res = pack_selector->apply( *overlaps_[kk].output_poseOP );

			TR << "select pack_res, resi ";
			for ( Size ii=1; ii<=pack_res.size(); ++ii ) {
				if ( pack_res[ii]==1 ) {
					TR << ii << "+";
				}
			}
			TR << std::endl;

			//selecting residue not allowed to pack/design
			core::select::residue_selector::ResidueSelectorCOP pack_and_design_selector(new core::select::residue_selector::OrResidueSelector(pack_selector,design_selector));
			core::select::residue_selector::ResidueSelectorCOP no_pack_selector(new core::select::residue_selector::NotResidueSelector(pack_and_design_selector));
			//core::select::residue_selector::ResidueSelectorCOP symmetric_no_pack_selector(new core::select::residue_selector::SymmetricalResidueSelector(no_pack_selector));
			utility::vector1< bool > no_pack_res = no_pack_selector->apply( *overlaps_[kk].output_poseOP );

			TR << "select no_pack_res, resi ";
			for ( Size ii=1; ii<=no_pack_res.size(); ++ii ) {
				if ( no_pack_res[ii]==1 ) {
					TR << ii << "+";
				}
			}
			TR << std::endl;

			//set up non-packing
			OperateOnResidueSubsetOP operate_on_non_packing = OperateOnResidueSubsetOP( new OperateOnResidueSubset() );
			PreventRepackingRLTOP prevent_repacking = PreventRepackingRLTOP( new PreventRepackingRLT() );
			operate_on_non_packing->subset( no_pack_res );
			operate_on_non_packing->op( prevent_repacking );
			taskfactoryOP->push_back(operate_on_non_packing);

			//set up non-design
			OperateOnResidueSubsetOP operate_on_packing = OperateOnResidueSubsetOP( new OperateOnResidueSubset() );
			RestrictToRepackingRLTOP restrict_to_repacking = RestrictToRepackingRLTOP( new RestrictToRepackingRLT() );
			operate_on_packing->subset( pack_res );
			operate_on_packing->op( restrict_to_repacking );
			taskfactoryOP->push_back(operate_on_packing);
			Size pack_rounds=5;
			protocols::minimization_packing::PackRotamersMover packer = protocols::minimization_packing::PackRotamersMover(sfxn_);
			packer.nloop( pack_rounds );
			packer.task_factory(taskfactoryOP);

			//pack
			using namespace protocols::minimization_packing;
			if ( pose_symmetric ) {
				protocols::minimization_packing::symmetry::SymPackRotamersMover symm_packer = PackRotamersMover(packer);
				symm_packer.apply(*overlaps_[kk].output_poseOP);
			} else {
				packer.apply(*overlaps_[kk].output_poseOP);
			}
			//minimize
			if ( do_minimize_ ) {
				if ( pose_symmetric ) {
					kinematics::MoveMapOP mm_locOP( new core::kinematics::MoveMap() );
					mm_locOP->set_jump( false ); mm_locOP->set_bb( false ); mm_locOP->set_chi( true );
					protocols::minimization_packing::symmetry::SymMinMover symm_min_mover(mm_locOP,sfxn_,"lbfgs_armijo_nonmonotone",0.02,true);
					symm_min_mover.apply(*overlaps_[kk].output_poseOP);
				} else {
					optimization::MinimizerOptions minopt( "lbfgs_armijo_nonmonotone", 0.02, true, false, false );
					//minopt.max_iter( 25 );
					optimization::AtomTreeMinimizer minimizer;
					kinematics::MoveMap mm_loc;
					mm_loc.set_jump( false ); mm_loc.set_bb( false ); mm_loc.set_chi( true );
					minimizer.run( *overlaps_[kk].output_poseOP, mm_loc, *sfxn_, minopt );
				}
			}

			//checks to see if the pose is at least better than the baseline_score
			//baseline_score = better score between input_pose and xml_pose)
			Real score = sfxn_->score(*overlaps_[kk].output_poseOP);
			TR << "score after pack and min" << score << "baseline_score" << baseline_score << std::endl;
			if ( score>baseline_score ) {
				overlaps_[kk].output_poseOP = NULL;
			}
			if ( score<baseline_score && output_only_first_ ) {
				early_complete=true;
			}
		}
	} //end overlaps (kk)
} //end pack_and_minimize

core::pose::PoseOP MergePDBMover::get_additional_output(){
	for ( Size ii=1; ii<=overlaps_.size(); ++ii ) {
		if ( overlaps_[ii].output_yet==false && overlaps_[ii].output_poseOP != NULL ) {
			set_last_move_status(protocols::moves::MS_SUCCESS);
			overlaps_[ii].output_yet=true;
			if ( output_overlap_positions_ ) {
				std::string design_name = overlaps_[ii].name();
				protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
				job->inner_job_nonconst()->optional_output_name(design_name);
			}
			return(overlaps_[ii].output_poseOP);
		}
	}
	set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
	return nullptr;
}

void MergePDBMover::apply( Pose & pose )
{
	core::Real baseline_score = sfxn_->score(pose);
	core::Real xml_pose_score =0;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::scoring::ScoreFunctionOP sfxn_clone = sfxn_->clone();
		asymm_score_ = core::scoring::symmetry::asymmetrize_scorefunction(*sfxn_clone);
		xml_pose_score = asymm_score_->score(*xml_input_pose_);
	} else {
		xml_pose_score = sfxn_->score(*xml_input_pose_);
	}
	if ( baseline_score<xml_pose_score ) {
		baseline_score = xml_pose_score;
	}
	//check if selected chain exists in the pose
	//vector1< core::Size >chains = get_chains(pose);
	//Size chain_id = 1;
	//if ( chains.size()>1 ) {
	if ( !has_chain(chain_,pose) ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "You have chosen a chain that does not exist");
	}
	Size chain_id = get_chain_id_from_chain(chain_,pose);
	TR << "Operating on chain " << chain_id << std::endl;
	//}
	//Step 1 get location
	TR << "STEP1: Get Location (determine_overlap)" << std::endl;
	determine_overlap(pose,chain_id);
	if ( overlaps_.size()==0 ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "No overlaps_ detected below rmsd cutoff");
	}
	//Step 2 combine
	TR << "STEP2: Combine (generate_overlaps)" << std::endl;
	generate_overlaps(pose,chain_id);

	if ( do_design_ ) {
		//Step 3 combine sequence
		TR << "STEP3: Combine Sequence (pack_and_minimize)" << std::endl;
		pack_and_minimize(pose,baseline_score);
	} else {
		TR << "Skipping STEP3; do_design=" << do_design_ << std::endl;
	}

	core::pose::PoseOP tmpPoseOP=get_additional_output();
	if ( tmpPoseOP!=nullptr ) {
		pose=*tmpPoseOP;
	}
}

void MergePDBMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	overlap_location_pose_ = tag->getOption< std::string >( "attachment_termini" ,"n_term" );
	chain_ = tag->getOption<std::string>("chain","A");
	overlap_length_ = tag->getOption< core::Size >( "overlap_length", 4);
	overlap_max_rmsd_ = tag->getOption< Real >( "overlap_rmsd", 0.3);
	overlap_scan_range_cmdLine_ = tag->getOption< core::Size >( "overlap_scan_range_cmdLine_pose", 20);
	overlap_scan_range_xml_ = tag->getOption< core::Size >( "overlap_scan_range_xml_pose", 20);
	std::string fname( tag->getOption< std::string >( "attach_pdb" ) );
	design_range_ = tag->getOption<Real>("design_range",9);
	packing_range_ = tag->getOption<Real>("packing_range",8);
	do_minimize_ = tag->getOption<bool>("do_minimize",true);
	symm_file_ = tag->getOption<std::string>("symm_file","");
	no_design_label_ = tag->getOption<std::string>("no_design_label","");
	init_overlap_sequence_ = tag->getOption<std::string>("init_overlap_sequence","input_pose");
	detect_disulf_before_repack_ = tag->getOption<bool>("detect_disulf_before_pack",true);
	duplicate_rmsd_pose_threshold_ = tag->getOption<Real>("duplicate_rmsd_pose_threshold",1.0);
	output_only_first_ = tag->getOption< bool >( "output_only_first" ,false );
	output_overlap_positions_ = tag->getOption<bool>("output_overlap_positions",false);
	xml_input_pose_ = core::pose::PoseOP( new pose::Pose );
	do_design_ = tag->getOption<bool>("do_design",true);
	clash_threshold_ = tag->getOption<Real>("clash_threshold",10.0);
	core::import_pose::pose_from_file(*xml_input_pose_, fname , core::import_pose::PDB_file);
	if ( tag->hasOption("scorefxn") ) {
		std::string const scorefxn_key( tag->getOption<std::string>("scorefxn") );
		if ( data.has( "scorefxns", scorefxn_key ) ) {
			sfxn_ = data.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_key );
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap.");
		}
	}
	if ( tag->hasOption("task_operations") ) {
		task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data );
	} else {
		task_factory_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory );
	}

}

std::string MergePDBMover::get_name() const {
	return mover_name();
}

std::string MergePDBMover::mover_name() {
	return "MergePDB";
}

void MergePDBMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	XMLSchemaRestriction attachment_termini_type;
	attachment_termini_type.name( "attachment_termini_type" );
	attachment_termini_type.base_type( xs_string );
	attachment_termini_type.add_restriction( xsr_enumeration, "n_term" );
	attachment_termini_type.add_restriction( xsr_enumeration, "c_term" );
	xsd.add_top_level_element( attachment_termini_type );
	XMLSchemaRestriction init_overlap_sequence_type;
	init_overlap_sequence_type.name( "init_overlap_sequence_type" );
	init_overlap_sequence_type.base_type(xs_string);
	init_overlap_sequence_type.add_restriction( xsr_enumeration, "input_pose" ); //init overlap from input pose
	init_overlap_sequence_type.add_restriction( xsr_enumeration, "xml_pose" ); //init overlap from xml_pose
	init_overlap_sequence_type.add_restriction( xsr_enumeration, "both" ); //init from whichever pose is closest
	xsd.add_top_level_element( init_overlap_sequence_type );
	attlist
		+ XMLSchemaAttribute::attribute_w_default("symm_file", xs_string, "Symmetry definition file if pose symmetric. hack for now","")
		+ XMLSchemaAttribute::attribute_w_default( "attachment_termini", "attachment_termini_type", "termini to  attach to, can be c_term or n_term", "n_term" )
		+ XMLSchemaAttribute::attribute_w_default( "chain", xs_string, "chain", "A" )
		+ XMLSchemaAttribute::attribute_w_default( "overlap_length", xsct_non_negative_integer, "length of overlap between the two structures", "4" )
		+ XMLSchemaAttribute::attribute_w_default( "overlap_rmsd", xsct_real, "How similiar the structures must be", "0.3" )
		+ XMLSchemaAttribute::attribute_w_default( "overlap_scan_range_cmdLine_pose", xsct_non_negative_integer, "how far in from the term to scan on the cmdLine input pose", "20" )
		+ XMLSchemaAttribute::attribute_w_default( "overlap_scan_range_xml_pose", xsct_non_negative_integer, "how far in from the term to scan on the xml input pose", "20" )
		+ XMLSchemaAttribute::required_attribute( "attach_pdb", xs_string, "the pdb to be attached, this is the pdb that moves" )
		+ XMLSchemaAttribute::attribute_w_default( "design_range", xsct_real, "distance from attachment and distance between new interacting residues allowed to design", "9" )
		+ XMLSchemaAttribute::attribute_w_default( "packing_range", xsct_real, "distance from designable_residues allowed to pack", "8" )
		+ XMLSchemaAttribute::attribute_w_default( "detect_disulf_before_pack", xsct_rosetta_bool, "detects disulfides before repacking, be sure to elimate them first with -detect_disulf false flag", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "do_minimize", xsct_rosetta_bool, "Perform energy minimization", "true" )
		+ XMLSchemaAttribute::attribute_w_default("no_design_label", xs_string, "residues to not design or pack","")
		+ XMLSchemaAttribute::attribute_w_default("duplicate_rmsd_pose_threshold", xsct_real, "if poses are the same length they are eliminated if the rmsd check","1.0")
		+ XMLSchemaAttribute::attribute_w_default("init_overlap_sequence","init_overlap_sequence_type","where to take the overlap sequence from input_pose, xml_pose, or both which takes from closest element on a different SS.","input_pose")
		+ XMLSchemaAttribute( "scorefxn", xs_string, "Score function used for packing and design." )
		+ XMLSchemaAttribute::attribute_w_default( "output_only_first", xsct_rosetta_bool, "only does minimization on the first overlap and outputs", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "output_overlap_positions", xsct_rosetta_bool, "outputs overlap positions", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "do_design", xsct_rosetta_bool, "Perform design on sequence", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "clash_threshold", xsct_real, "score0 clash threshold", "10" );
	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Align + combine parts of the pdb", attlist );
}

std::string MergePDBMoverCreator::keyname() const {
	return MergePDBMover::mover_name();
}

protocols::moves::MoverOP
MergePDBMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MergePDBMover );
}

void MergePDBMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MergePDBMover::provide_xml_schema( xsd );
}

} // pose_creation
} // protocols
