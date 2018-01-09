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

static basic::Tracer TR( "protocols.simple_moves.MergePDBMover" );

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
#include <core/select/residue_selector/SymmetricalResidueSelector.hh>

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

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/simple_moves/CopyRotamerMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/AddResidueLabelMover.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace simple_moves {

using std::map;
using namespace core;
using core::pose::Pose;
using utility::vector1;


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

utility::vector1<MergePDBMover::Overlap> MergePDBMover::determine_overlap(Pose const pose,Size chain_id){
	vector1<MergePDBMover::Overlap> hits;
	utility::vector1<core::Size> res_in_chain = get_resnums_for_chain_id(pose,chain_id);
	Size start_res_chain = res_in_chain[1];
	Size end_res_chain = chain_end_res(pose,chain_id);
	using namespace core::scoring;
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
				//std::cout << "rmsd" << rmsd << std::endl;
				TR << "rmsd" << rmsd << "xml_pose_position" << start_res_xml_pose << ",input_pose_position " <<start_res_pose << std::endl;
				if ( rmsd < overlap_max_rmsd_ ) {
					MergePDBMover::Overlap overlap_tmp(start_res_xml_pose,end_res_xml_pose,start_res_pose,end_res_pose);
					hits.push_back(overlap_tmp);
				}
			}
		}
	}
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
				TR << "rmsd" << rmsd << "xml_pose_position" << start_res_xml_pose << ",input_pose_position " <<start_res_pose << std::endl;
				if ( rmsd < overlap_max_rmsd_ ) {
					MergePDBMover::Overlap overlap_tmp(start_res_xml_pose,end_res_xml_pose,start_res_pose,end_res_pose);
					hits.push_back(overlap_tmp);
				}
			}
		}
	}
	return(hits);
}
bool MergePDBMover::check_duplicate(Pose &pose, utility::vector1<core::pose::PoseOP> outputPoses){
	//setup to be fast because this could be a N^2 runtime.
	Real min_distance=999;
	Size min_position=999;
	for ( Size ii=1; ii<=outputPoses.size(); ++ii ) {
		if ( pose.total_residue() == outputPoses[ii]->total_residue() ) {
			Size n_term_res = 1;
			Size c_term_res = pose.total_residue();
			while ( pose.residue(n_term_res).is_virtual_residue() )
					n_term_res++;
			while ( pose.residue(c_term_res).is_virtual_residue() )
					c_term_res--;
			Real n_term_dist = pose.residue(n_term_res).xyz("CA").distance(outputPoses[ii]->residue(n_term_res).xyz("CA"));
			Real c_term_dist = pose.residue(c_term_res).xyz("CA").distance(outputPoses[ii]->residue(c_term_res).xyz("CA"));
			Real dist = n_term_dist+c_term_dist;
			if ( dist<min_distance ) {
				min_distance=dist;
				min_position=ii;
			}
		}
	}
	Real min_rmsd_dist = 999;
	if ( min_distance<3.0 ) { //only calc rmsd if needed
		min_rmsd_dist = core::scoring::CA_rmsd(pose,*(outputPoses[min_position]));
	}
	if ( min_rmsd_dist<duplicate_rmsd_pose_threshold_ ) {
		return true;
	}
	return(false);
}


vector1<core::pose::PoseOP> MergePDBMover::generate_overlaps(Pose & pose, vector1<MergePDBMover::Overlap> overlaps,Size chain_id) {
	using namespace core::id;
	using namespace core::scoring;
	vector1<core::pose::PoseOP> outputPoses;
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
	for ( Size kk=1; kk<=overlaps.size(); ++kk ) {
		vector1<Size> nterm_insert_res;
		vector1<Size> cterm_insert_res;
		Size start_overlap_pose = overlaps[kk].start_overlap_pose;
		Size end_overlap_pose = overlaps[kk].end_overlap_pose;
		Size start_overlap_xmlPose = overlaps[kk].start_overlap_xmlPose;
		Size end_overlap_xmlPose = overlaps[kk].end_overlap_xmlPose;
		core::id::AtomID_Map< core::id::AtomID > atom_map;
		core::pose::initialize_atomid_map( atom_map, *xml_input_pose_, AtomID::BOGUS_ATOM_ID() );
		for ( Size ii=0; ii<=end_overlap_pose-start_overlap_pose; ++ii ) {
			core::id::AtomID const id2(pose.residue(start_overlap_pose+ii).atom_index("CA"),start_overlap_pose+ii);
			core::id::AtomID const id1(xml_input_pose_->residue(start_overlap_xmlPose+ii).atom_index("CA"), start_overlap_xmlPose+ii );
			atom_map[id1]=id2;
		}
		// remove_lower_terminus_type_from_pose_residue(*xml_input_pose_slice,1);
		superimpose_pose(*xml_input_pose_,pose,atom_map);
		//create subpose of xml_pose
		core::pose::PoseOP output_poseOP( new core::pose::Pose() );
		core::pose::PoseOP tmp_xml_poseOP;
		tmp_xml_poseOP=xml_input_pose_->clone();
		Size start_steal_label_res = 1;
		Size end_steal_label_res = pose.size();
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
		remove_upper_terminus_type_from_pose_residue(*backup_input_pose,backup_input_pose->total_residue());
		remove_lower_terminus_type_from_pose_residue(*backup_input_pose,1);
		remove_upper_terminus_type_from_pose_residue(*backup_xml_pose,backup_xml_pose->total_residue());
		remove_lower_terminus_type_from_pose_residue(*backup_xml_pose,1);
		if ( overlap_location_pose_== "n_term" ) {
			remove_upper_terminus_type_from_pose_residue(*tmp_xml_poseOP,tmp_xml_poseOP->total_residue());
			remove_lower_terminus_type_from_pose_residue(*output_poseOP,start_res_chain);
			//delete any additional residues
			output_poseOP->conformation().delete_residue_range_slow(start_res_chain,start_overlap_pose); //always delete the first residue because of the lost phi/psi
			for ( Size ii=start_overlap_xmlPose; ii>=1; --ii ) {
				output_poseOP->conformation().safely_prepend_polymer_residue_before_seqpos(tmp_xml_poseOP->residue(ii),start_res_chain, false);
			}
			//add initial residue labels
			Size res_loss_n_term=start_overlap_pose-start_res_chain;
			Size res_add =start_overlap_xmlPose;
			Size res_offset = res_add-res_loss_n_term-1;
			for ( core::Size ii = start_steal_label_res; ii <= end_steal_label_res; ++ii ) { //chose to do this for the output_pose_size. Should take care of symmetric case
				initial_labels =  pose.pdb_info()->get_reslabels(ii);
				for ( core::Size jj=1; jj <= initial_labels.size(); ++jj ) {
					output_poseOP->pdb_info()->add_reslabel(res_offset+ii, initial_labels[jj]);
				}
			}
			//add overlap residue label
			for ( core::Size ii=start_overlap_xmlPose; ii<=end_overlap_xmlPose; ++ii ) {
				output_poseOP->pdb_info()->add_reslabel(ii, "overlap");
			}
			for ( core::Size ii=1; ii<=end_overlap_xmlPose; ++ii ) {
				nterm_insert_res.push_back(ii);
			}
			for ( core::Size ii=start_overlap_xmlPose; ii<=output_poseOP->total_residue(); ++ii ) {
				cterm_insert_res.push_back(ii);
			}
			copy_sequence(start_overlap_xmlPose,end_overlap_xmlPose,start_overlap_pose,start_overlap_xmlPose,*backup_input_pose,*backup_xml_pose,*output_poseOP);
		}
		if ( overlap_location_pose_== "c_term" ) {
			remove_lower_terminus_type_from_pose_residue(*tmp_xml_poseOP,1);
			remove_upper_terminus_type_from_pose_residue(*output_poseOP,end_res_chain);
			output_poseOP->conformation().delete_residue_range_slow(end_overlap_pose,end_res_chain);
			for ( Size ii=end_overlap_xmlPose; ii<=tmp_xml_poseOP->total_residue(); ++ii ) {
				Size location_to_insert=end_overlap_pose-1+ii-end_overlap_xmlPose;
				output_poseOP->conformation().safely_append_polymer_residue_after_seqpos(tmp_xml_poseOP->residue(ii),location_to_insert, false);
			}
			//add initial residue labels
			for ( core::Size ii = 1; ii <= end_overlap_pose; ++ii ) {
				initial_labels =  pose.pdb_info()->get_reslabels(ii);
				for ( core::Size jj=1; jj <= initial_labels.size(); ++jj ) {
					output_poseOP->pdb_info()->add_reslabel(ii, initial_labels[jj]);
				}
			}
			//add overap residue labels
			for ( core::Size ii=start_overlap_pose; ii<=end_overlap_pose; ++ii ) {
				output_poseOP->pdb_info()->add_reslabel(ii, "overlap");
			}
			for ( core::Size ii=1; ii<=end_overlap_pose; ++ii ) {
				nterm_insert_res.push_back(ii);
			}
			for ( core::Size ii=start_overlap_pose; ii<=output_poseOP->total_residue(); ++ii ) {
				cterm_insert_res.push_back(ii);
			}
			copy_sequence(start_overlap_pose,end_overlap_pose,start_overlap_pose,start_overlap_xmlPose,*backup_input_pose,*backup_xml_pose,*output_poseOP);
		}
		//seems really weird. When I created a brand new scorefunction it was symmetric
		if ( input_pose_symmetric ) {
			asymm_score_->score(*output_poseOP); //needs to be here to calc distances
		} else {
			sfxn_->score(*output_poseOP);
		}
		//step 1 add location where chain A and B now overlap as overlap_2
		core::select::residue_selector::ResidueSelectorCOP overlap_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector("overlap"));
		core::select::residue_selector::ResidueSelectorCOP interface_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector(no_design_label_));
		core::select::residue_selector::ResidueSelectorCOP overlap_plus_interface_selector(new core::select::residue_selector::OrResidueSelector(overlap_selector,interface_selector));
		core::select::residue_selector::ResidueSelectorCOP not_overlap_interface_selector(new core::select::residue_selector::NotResidueSelector(overlap_plus_interface_selector));
		core::select::residue_selector::ResidueSelectorCOP partA_selector(new core::select::residue_selector::ResidueIndexSelector(nterm_insert_res));
		core::select::residue_selector::ResidueSelectorCOP partB_selector(new core::select::residue_selector::ResidueIndexSelector(cterm_insert_res));
		core::select::residue_selector::ResidueSelectorCOP interface_partA_selector(new core::select::residue_selector::NeighborhoodResidueSelector(partB_selector,design_range_,false));
		core::select::residue_selector::ResidueSelectorCOP interface_partB_selector(new core::select::residue_selector::NeighborhoodResidueSelector(partA_selector,design_range_,false));
		core::select::residue_selector::ResidueSelectorCOP interface_AB_selector(new core::select::residue_selector::OrResidueSelector(interface_partA_selector,interface_partB_selector));
		core::select::residue_selector::ResidueSelectorCOP interface_AB_exclude_overlap_selector(new core::select::residue_selector::AndResidueSelector(interface_AB_selector,not_overlap_interface_selector));
		protocols::simple_moves::AddResidueLabelMover add_other_overlap = AddResidueLabelMover(interface_AB_exclude_overlap_selector,"other_overlap");
		//debug code::
		add_other_overlap.apply(*output_poseOP);
		//store res labels
		std::map<core::Size, vector1<std::string> > res_label_map;
		for ( core::Size ii = 1; ii <=output_poseOP->size(); ++ii ) {
			vector1<std::string> tmp_labels = output_poseOP->pdb_info()->get_reslabels(ii);
			res_label_map.insert(std::pair<Size,vector1<std::string> >(ii,tmp_labels));
		}
		renumber_pdbinfo_based_on_conf_chains(*output_poseOP,true,false,false,false);
		if ( input_pose_symmetric ) {
			core::pose::symmetry::make_symmetric_pose(*output_poseOP,symm_file_);
			//lost labels again when applying symmetry
		}

		//reapply residue labels
		std::map<core::Size, vector1<std::string> >::iterator res_label_map_itr;
		vector1<std::string>::iterator res_label_itr;
		for ( res_label_map_itr=res_label_map.begin(); res_label_map_itr!=res_label_map.end(); ++res_label_map_itr ) {
			Size resid = res_label_map_itr->first;
			vector1<std::string> tmp_labels = res_label_map_itr->second;
			for ( res_label_itr=tmp_labels.begin(); res_label_itr!=tmp_labels.end(); ++res_label_itr ) {
				output_poseOP->pdb_info()->add_reslabel(resid, *res_label_itr);
			}
		}
		//label where the chains hit-each but not interface
		//Second detect inter-chain
		vector1< core::Size >chains = get_chains(*output_poseOP);
		if ( chains.size()>1 ) {
			//chainA = shorthand for selected chain
			sfxn_->score(*output_poseOP);
			core::select::residue_selector::ResidueSelectorCOP chain_selector(new core::select::residue_selector::ChainSelector(chain_));
			core::select::residue_selector::ResidueSelectorCOP not_chainA_selector(new core::select::residue_selector::NotResidueSelector(chain_selector));
			core::select::residue_selector::ResidueSelectorCOP interface_chainA_selector(new core::select::residue_selector::NeighborhoodResidueSelector(not_chainA_selector,design_range_,false));
			core::select::residue_selector::ResidueSelectorCOP interface_not_chainA_selector(new core::select::residue_selector::NeighborhoodResidueSelector(not_chainA_selector,design_range_,false));
			core::select::residue_selector::ResidueSelectorCOP interface_btwChain_selector(new core::select::residue_selector::OrResidueSelector(interface_chainA_selector,interface_not_chainA_selector));
			core::select::residue_selector::ResidueSelectorCOP symmetric_interface_btwChain_selector(new core::select::residue_selector::SymmetricalResidueSelector(interface_btwChain_selector));
			core::select::residue_selector::ResidueSelectorCOP symmetric_not_overlap_interface_selector(new core::select::residue_selector::SymmetricalResidueSelector(not_overlap_interface_selector));
			core::select::residue_selector::ResidueSelectorCOP interface_btwChain_exclude_overlap_selector(new core::select::residue_selector::AndResidueSelector(symmetric_interface_btwChain_selector,symmetric_not_overlap_interface_selector));
			protocols::simple_moves::AddResidueLabelMover add_other_overlap_v2 = AddResidueLabelMover(interface_btwChain_exclude_overlap_selector,"other_overlap");
			add_other_overlap_v2.apply(*output_poseOP);
		}
		//discard poses with bad clashes
		core::pose::PoseOP centroidPoseOP = output_poseOP->clone();
		tocen->apply( *centroidPoseOP );
		Real score0 = sfxn0->score( *centroidPoseOP );
		TR << "score0 for clashcheck (cutoff 10):" << score0 << std::endl;
		bool duplicate = check_duplicate(*output_poseOP,outputPoses);
		if ( score0<10.0 && !duplicate ) { //only proceed if few clashes
			outputPoses.push_back(output_poseOP);
		}
	}
	return(outputPoses);
}

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



void MergePDBMover::copy_sequence(core::Size start_overlap_resid, core::Size end_overlap_resid, core::Size start_overlap_input_pose_resid,core::Size start_overlap_xml_pose_resid,core::pose::Pose const & input_pose,core::pose::Pose const & xml_pose,core::pose::Pose & output_pose){
	protocols::simple_moves::CopyRotamerMover copyRotamerMover = CopyRotamerMover();
	Size numb_res = end_overlap_resid-start_overlap_resid+1;
	if ( init_overlap_sequence_=="input_pose" ) {
		for ( Size ii=0; ii<numb_res; ++ii ) {
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


vector1<core::pose::PoseOP> MergePDBMover::pack_and_minimize(vector1<core::pose::PoseOP> poses,Real baseline_score) {
	using namespace core::pack::task;
	using namespace core::scoring;
	using namespace core::pack::task::operation;
	vector1<core::pose::PoseOP> outputPoses;
	bool pose_symmetric = core::pose::symmetry::is_symmetric(*poses[1]);
	for ( Size kk=1; kk<=poses.size(); ++kk ) {
		if ( detect_disulf_before_repack_ ) {
			poses[kk]->conformation().detect_disulfides();
		}
		TaskFactoryOP taskfactorOP = task_factory_->clone();
		//PackerTaskOP taskOP = task_factory_->create_task_and_apply_taskoperations(*poses[kk]);

		//find residues allowed to design
		core::select::residue_selector::ResidueSelectorCOP overlap_p1_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector("overlap"));
		core::select::residue_selector::ResidueSelectorCOP overlap_p2_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector("other_overlap"));
		core::select::residue_selector::ResidueSelectorCOP overlap_selector(new core::select::residue_selector::OrResidueSelector(overlap_p1_selector,overlap_p2_selector));
		//core::select::residue_selector::ResidueSelectorCOP near_overlap_selector(new core::select::residue_selector::NeighborhoodResidueSelector(overlap_selector,design_range_,false));
		core::select::residue_selector::ResidueSelectorCOP interface_selector(new core::select::residue_selector::ResiduePDBInfoHasLabelSelector(no_design_label_));
		//core::select::residue_selector::ResidueSelectorCOP design_plus_interface_selector(new core::select::residue_selector::OrResidueSelector(overlap_selector,near_overlap_selector));
		core::select::residue_selector::ResidueSelectorCOP no_design_selector(new core::select::residue_selector::NotResidueSelector(interface_selector));
		core::select::residue_selector::ResidueSelectorCOP design_selector(new core::select::residue_selector::AndResidueSelector(overlap_selector,no_design_selector));
		core::select::residue_selector::ResidueSelectorCOP symmetric_design_selector(new core::select::residue_selector::SymmetricalResidueSelector(design_selector));
		utility::vector1< bool > design_res = symmetric_design_selector->apply( *poses[kk] );

		TR << std::endl;
		TR << "designed_res" << std::endl;
		for ( Size ii=1; ii<=design_res.size(); ++ii ) {
			if ( design_res[ii]==1 ) {
				TR << ii << "+";
			}
		}
		TR << std::endl;

		//selecting residues allowed to pack
		core::select::residue_selector::ResidueSelectorCOP near_overlap_selector(new core::select::residue_selector::NeighborhoodResidueSelector(overlap_selector,packing_range_,false));
		core::select::residue_selector::ResidueSelectorCOP pack_selector(new core::select::residue_selector::AndResidueSelector(near_overlap_selector,no_design_selector));
		core::select::residue_selector::ResidueSelectorCOP symmetric_pack_selector(new core::select::residue_selector::SymmetricalResidueSelector(pack_selector));
		utility::vector1< bool > pack_res = symmetric_pack_selector->apply( *poses[kk] );
		TR << "pack_res" << std::endl;
		for ( Size ii=1; ii<=pack_res.size(); ++ii ) {
			if ( pack_res[ii]==1 ) {
				TR << ii << "+";
			}
		}
		TR << std::endl;
		//selecting residue not allowed to pack/design
		core::select::residue_selector::ResidueSelectorCOP pack_and_design_selector(new core::select::residue_selector::OrResidueSelector(pack_selector,design_selector));
		core::select::residue_selector::ResidueSelectorCOP no_pack_selector(new core::select::residue_selector::NotResidueSelector(pack_and_design_selector));
		core::select::residue_selector::ResidueSelectorCOP symmetric_no_pack_selector(new core::select::residue_selector::SymmetricalResidueSelector(no_pack_selector));

		utility::vector1< bool > no_pack_res = symmetric_no_pack_selector->apply( *poses[kk] );
		TR << "no_pack_res" << std::endl;
		for ( Size ii=1; ii<=no_pack_res.size(); ++ii ) {
			if ( no_pack_res[ii]==1 ) {
				TR << ii << "+";
			}
		}
		TR << std::endl;
		//find residue not allowed to pack or design

		//set up non-packing
		OperateOnResidueSubsetOP operate_on_non_packing = OperateOnResidueSubsetOP( new OperateOnResidueSubset() );
		PreventRepackingRLTOP prevent_repacking = PreventRepackingRLTOP( new PreventRepackingRLT() );
		operate_on_non_packing->subset( no_pack_res );
		operate_on_non_packing->op( prevent_repacking );
		taskfactorOP->push_back(operate_on_non_packing);

		//set up non-design
		OperateOnResidueSubsetOP operate_on_packing = OperateOnResidueSubsetOP( new OperateOnResidueSubset() );
		RestrictToRepackingRLTOP restrict_to_repacking = RestrictToRepackingRLTOP( new RestrictToRepackingRLT() );
		operate_on_packing->subset( pack_res );
		operate_on_packing->op( restrict_to_repacking );
		taskfactorOP->push_back(operate_on_packing);
		Size pack_rounds=5;
		protocols::simple_moves::PackRotamersMover packer = PackRotamersMover(sfxn_);
		packer.nloop( pack_rounds );
		packer.task_factory(taskfactorOP);
		//pack
		if ( pose_symmetric ) {
			protocols::simple_moves::symmetry::SymPackRotamersMover symm_packer = PackRotamersMover(packer);
			symm_packer.apply(*poses[kk]);
		} else {
			packer.apply(*poses[kk]);
		}
		//minimize
		if ( do_minimize_ ) {
			if ( pose_symmetric ) {
				kinematics::MoveMapOP mm_locOP( new core::kinematics::MoveMap() );
				mm_locOP->set_jump( false ); mm_locOP->set_bb( false ); mm_locOP->set_chi( true );
				protocols::simple_moves::symmetry::SymMinMover symm_min_mover(mm_locOP,sfxn_,"lbfgs_armijo_nonmonotone",0.02,true);
				symm_min_mover.apply(*poses[kk]);
			} else {
				optimization::MinimizerOptions minopt( "lbfgs_armijo_nonmonotone", 0.02, true, false, false );
				//minopt.max_iter( 25 );
				optimization::AtomTreeMinimizer minimizer;
				kinematics::MoveMap mm_loc;
				mm_loc.set_jump( false ); mm_loc.set_bb( false ); mm_loc.set_chi( true );
				minimizer.run( *poses[kk], mm_loc, *sfxn_, minopt );
			}
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
	vector1< core::Size >chains = get_chains(pose);
	Size chain_id = 1;
	if ( chains.size()>1 ) {
		if ( !has_chain(chain_,pose) ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "You have chosen a chain that does not exist");
		}
		chain_id = get_chain_id_from_chain(chain_,pose);
	}
	//Step 1 get location
	vector1<Overlap> overlaps = determine_overlap(pose,chain_id);
	if ( overlaps.size()==0 ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "No overlaps detected below rmsd cutoff");
	}
	//Step 2 combine
	vector1<core::pose::PoseOP> combined_poses = generate_overlaps(pose,overlaps,chain_id);
	if ( combined_poses.size()==0 ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "No overlaps detected after clash check");
	}
	outputPoses_ = pack_and_minimize(combined_poses,baseline_score);

	for ( Size ii=1; ii<=outputPoses_.size(); ++ii ) {
		outputYet_.push_back(false);
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
	design_range_ = tag->getOption<Real>("design_range",3);
	packing_range_ = tag->getOption<Real>("packing_range",5);
	do_minimize_ = tag->getOption<bool>("do_minimize",true);
	symm_file_ = tag->getOption<std::string>("symm_file","");
	no_design_label_ = tag->getOption<std::string>("no_design_label","");
	init_overlap_sequence_ = tag->getOption<std::string>("init_overlap_sequence","input_pose");
	detect_disulf_before_repack_ = tag->getOption<bool>("detect_disulf_before_pack",true);
	duplicate_rmsd_pose_threshold_ = tag->getOption<Real>("duplicate_rmsd_pose_threshold(only works <1.0 rmsd for speed optimization)",1.0);
	xml_input_pose_ = core::pose::PoseOP( new pose::Pose );
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
		+ XMLSchemaAttribute::attribute_w_default( "overlap_length", xsct_non_negative_integer, "XRW TO DO", "4" )
		+ XMLSchemaAttribute::attribute_w_default( "overlap_rmsd", xsct_real, "How similiar the structures must be", "0.3" )
		+ XMLSchemaAttribute::attribute_w_default( "overlap_scan_range_cmdLine_pose", xsct_non_negative_integer, "how far in from the term to scan on the cmdLine input pose", "20" )
		+ XMLSchemaAttribute::attribute_w_default( "overlap_scan_range_xml_pose", xsct_non_negative_integer, "how far in from the term to scan on the xml input pose", "20" )
		+ XMLSchemaAttribute::required_attribute( "attach_pdb", xs_string, "the pdb to be attached, this is the pdb that moves" )
		+ XMLSchemaAttribute::attribute_w_default( "design_range", xsct_real, "dist from attachment allowed to design", "3" )
		+ XMLSchemaAttribute::attribute_w_default( "packing_range", xsct_real, "dist from attachment allowed to pack", "5" )
		+ XMLSchemaAttribute::attribute_w_default( "detect_disulf_before_pack", xsct_rosetta_bool, "detects disulfides before repacking, be sure to elimate them first with -detect_disulf false flag", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "do_minimize", xsct_rosetta_bool, "Perform energy minimization", "true" )
		+ XMLSchemaAttribute::attribute_w_default("no_design_label", xs_string, "residues to not design or pack","")
		+ XMLSchemaAttribute::attribute_w_default("duplicate_rmsd_pose_threshold", xsct_real, "if poses are the same length they are eliminated if the rmsd check","1.0")
		+ XMLSchemaAttribute::attribute_w_default("init_overlap_sequence","init_overlap_sequence_type","where to take the overlap sequence from input_pose, xml_pose, or both which takes from closest element on a different SS.","input_pose")
		+ XMLSchemaAttribute( "scorefxn", xs_string, "Score function used for packing and design." );
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


} // simple_moves
} // protocols
