// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/NearNativeLoopCloser.fwd.hh
/// @brief Uses KIC/CCD and a fast Look back to close loops near native
///
/// @author TJ Brunette tjbrunette@gmail.com
///
// Unit headers
#include <protocols/pose_length_moves/NearNativeLoopCloser.hh>
#include <protocols/pose_length_moves/NearNativeLoopCloserCreator.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/simple_moves/SwitchChainOrderMover.hh>
#include <protocols/jd2/util.hh>

// Core Headers
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/sequence/ABEGOManager.hh>
#include <core/util/SwitchResidueTypeSet.hh>


#include <protocols/indexed_structure_store/SSHashedFragmentStore.hh>
#include <protocols/indexed_structure_store/FragmentStore.hh>
#include <protocols/indexed_structure_store/FragmentStoreManager.hh>


#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>


#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/SSManager.hh>

#include <core/types.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>

#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

#include <numeric/conversions.hh>
#include <numeric/alignment/QCPKernel.hh>
#include <numeric/xyzVector.hh>

#include <Eigen/Core>

#include <vector>
#include <sstream>
#include <map>
//output
#include <utility/io/ozstream.hh>
#include <ObjexxFCL/format.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh> // AUTO IWYU For GeneralizedKIC

//#include <unistd.h>

static basic::Tracer TR( "protocols.pose_length_moves.NearNativeLoopCloser" );

namespace protocols {
namespace pose_length_moves {
using namespace core;
using namespace std;
using utility::vector1;

PossibleLoop::PossibleLoop(int resAdjustmentBeforeLoop, int resAdjustmentAfterLoop,core::Size loopLength,core::Size resBeforeLoop, core::Size resAfterLoop, char resTypeBeforeLoop, char resTypeAfterLoop, core::Size insertedBeforeLoopRes, core::Size insertedAfterLoopRes, core::pose::PoseOP fullLengthPoseOP, core::pose::PoseOP orig_atom_type_fullLengthPoseOP,std::string fragment_store_path,std::string fragment_store_format, std::string fragment_store_compression,core::Size numb_stubs_to_consider){
	TR << "creatingPotentialLoop" <<  resAdjustmentBeforeLoop << "," << resAdjustmentAfterLoop <<"," << loopLength << "," << resBeforeLoop << ","<< resAfterLoop << std::endl;
	resAdjustmentBeforeLoop_ = resAdjustmentBeforeLoop;
	resAdjustmentAfterLoop_ = resAdjustmentAfterLoop;
	loopLength_ = loopLength;
	resBeforeLoop_ = resBeforeLoop;
	fullLength_resBeforeLoop_ = resBeforeLoop;
	resAfterLoop_ = resAfterLoop;
	resTypeBeforeLoop_ = resTypeBeforeLoop;
	resTypeAfterLoop_ = resTypeAfterLoop;
	insertedBeforeLoopRes_ = insertedBeforeLoopRes;
	insertedAfterLoopRes_ = insertedAfterLoopRes;
	original_atom_type_fullLengthPoseOP_ = orig_atom_type_fullLengthPoseOP;
	fullLengthPoseOP_ = fullLengthPoseOP;
	below_distance_threshold_=false;
	stub_rmsd_top_match_= 9999;
	final_rmsd_=9999;
	outputed_ = false;
	fragment_store_path_ = fragment_store_path;
	fragment_store_format_ = fragment_store_format;
	fragment_store_compression_ = fragment_store_compression;
	numb_stubs_to_consider_ = numb_stubs_to_consider;
}

PossibleLoop::~PossibleLoop()= default;

void PossibleLoop::trimRegion(core::pose::PoseOP & poseOP, core::Size resStart, core::Size resStop){
	poseOP->conformation().delete_residue_range_slow(resStart,resStop);
	renumber_pdbinfo_based_on_conf_chains(*poseOP,true,false,false,false);
}


void PossibleLoop::extendRegion(bool towardCTerm, core::Size resStart, core::Size numberAddRes,core::pose::PoseOP & poseOP){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace chemical;
	//could just pick a random phi/psi
	Real tmpPhi = -57.8;
	Real tmpPsi = -47.0;
	Real tmpOmega = 180.0; //I extend with helical parameters. They get replaced during kic.
	core::conformation::ResidueOP new_rsd( nullptr );
	string build_aa_type_one_letter =option[OptionKeys::remodel::generic_aa];
	string build_aa_type = name_from_aa(aa_from_oneletter_code(build_aa_type_one_letter[0]));
	debug_assert( poseOP != nullptr );
	core::chemical::ResidueTypeSetCOP rs( poseOP->residue_type_set_for_pose() );
	kinematics::FoldTree ft;
	if ( towardCTerm == true ) {
		ft = poseOP->fold_tree();
		ft.clear();
		ft.add_edge(1,resStart,core::kinematics::Edge::PEPTIDE);
		ft.add_edge(1,poseOP->size(),1); //edges need to be added in order so the fold tree is completely connected.
		ft.add_edge(poseOP->size(),resStart+1,core::kinematics::Edge::PEPTIDE);
		poseOP->fold_tree(ft);
		for ( core::Size ii=0; ii<numberAddRes; ++ii ) {
			new_rsd = core::conformation::ResidueFactory::create_residue( rs->name_map(build_aa_type) );
			poseOP->conformation().safely_append_polymer_residue_after_seqpos( *new_rsd,resStart+ii, true);
		}
	} else {
		ft = poseOP->fold_tree();
		ft.clear();
		ft.add_edge(1,resStart-1,core::kinematics::Edge::PEPTIDE);
		ft.add_edge(1,poseOP->size(),1); //edges need to be added in order so the fold tree is completely connected.
		ft.add_edge(poseOP->size(),resStart,core::kinematics::Edge::PEPTIDE);
		poseOP->fold_tree(ft);
		for ( core::Size ii=0; ii<numberAddRes; ++ii ) {
			new_rsd = core::conformation::ResidueFactory::create_residue( rs->name_map(build_aa_type) );
			poseOP->conformation().safely_prepend_polymer_residue_before_seqpos( *new_rsd,resStart, true);
		}
	}
	for ( core::Size ii=0; ii<numberAddRes; ++ii ) {
		poseOP->set_phi(resStart+ii, tmpPhi );
		poseOP->set_psi(resStart+ii, tmpPsi );
		poseOP->set_omega(resStart+ii, tmpOmega );
	}
	poseOP->set_phi(resStart+numberAddRes, tmpPhi );
	ft = poseOP->fold_tree();
	ft.clear();
	ft.add_edge(1,poseOP->size(),core::kinematics::Edge::PEPTIDE);
	poseOP->fold_tree(ft);
	renumber_pdbinfo_based_on_conf_chains(*poseOP,true,false,false,false);
}

void PossibleLoop::evaluate_distance_closure(){
	Real maxDist = loopLength_*4.3; //4.3 is the new max. Originally this was 4 but in crystal structure 4tql there is a 3 residue loop
	//(72-74) that when measured from 71-75 measures 12.6.
	core::Size tmpResidueBeforeLoop = resBeforeLoop_+resAdjustmentBeforeLoop_;
	core::Size tmpResidueAfterLoop = resBeforeLoop_+insertedBeforeLoopRes_+insertedAfterLoopRes_-resAdjustmentAfterLoop_+1;
	Real currentDist = fullLengthPoseOP_->residue(tmpResidueBeforeLoop).xyz("CA").distance(fullLengthPoseOP_->residue(tmpResidueAfterLoop).xyz("CA"));
	//std::cout << "XresBeforeLoop_" << tmpResidueBeforeLoop << " resAfterLoop_" << tmpResidueAfterLoop << " loopLength_" << loopLength_ << "resBeforeLoop_" << resBeforeLoop_ << "resAfterLoop_" << resAfterLoop_ << " currentDist" << currentDist << "maxDist" << maxDist << std::endl;
	if ( currentDist < maxDist ) {
		below_distance_threshold_=true;
		stub_rmsd_top_match_=999;//so they get sorted to the front
	} else {
		below_distance_threshold_=false;
	}
}

void PossibleLoop::generate_overlap_range(core::Size & front_overlap, core::Size & back_overlap){
	//The overlap depends on the fragment length. Ideally we would like to have at least 3 residues. But if the loop length is 9 that is impossible
	core::Size fragment_length =  SSHashedFragmentStoreOP_->get_fragment_length();
	front_overlap = 2;
	back_overlap = 2;
	if ( fragment_length>(front_overlap+back_overlap+loopLength_) ) { //within tolerance add a residue to the front
		front_overlap+=1;
	}
	if ( fragment_length>(front_overlap+back_overlap+loopLength_) ) { //still within tolerance add a residue to the back
		back_overlap+=1;
	}
}


void PossibleLoop::generate_stub_rmsd(core::Real stubRmsdThreshold){
	SSHashedFragmentStoreOP_ = protocols::indexed_structure_store::FragmentStoreManager::get_instance()->SSHashedFragmentStore(fragment_store_path_,fragment_store_format_,fragment_store_compression_);
	core::Size tmpResidueBeforeLoop = resBeforeLoop_+resAdjustmentBeforeLoop_;
	core::Size tmpResidueAfterLoop = resBeforeLoop_+insertedBeforeLoopRes_+insertedAfterLoopRes_-resAdjustmentAfterLoop_+1;
	std::vector< numeric::xyzVector<numeric::Real> > coordinates;
	core::Size front_overlap;
	core::Size back_overlap;
	generate_overlap_range(front_overlap,back_overlap);
	vector1<core::Size> residues;
	for ( int ii=front_overlap-1; ii>=0; --ii ) {
		residues.push_back(tmpResidueBeforeLoop-ii);
	}
	for ( core::Size ii=0; ii<back_overlap; ++ii ) {
		residues.push_back(tmpResidueAfterLoop+ii);
	}
	for ( core::Size ii=1; ii<=residues.size(); ++ii ) {
		for ( std::string const & atom_name : SSHashedFragmentStoreOP_->get_fragment_store()->fragment_specification.fragment_atoms ) {
			coordinates.push_back(fullLengthPoseOP_->residue(residues[ii]).xyz(atom_name));
		}
	}
	stub_rmsd_top_match_= 9999;
	SSHashedFragmentStoreOP_->lookback_stub(coordinates,resTypeBeforeLoop_,resTypeAfterLoop_,loopLength_, stub_rmsd_top_match_,stubVector_,stubRmsdThreshold);
}

/* No longer used
void PossibleLoop::generate_uncached_stub_rmsd(){
SSHashedFragmentStoreOP_ = protocols::indexed_structure_store::FragmentStoreManager::get_instance()->SSHashedFragmentStore(fragment_store_path_,fragment_store_format_,fragment_store_compression_);
Real uncached_stub_rmsd_= 9999;
core::Size uncached_stub_index_= 9999;
core::Size tmpResidueBeforeLoop = resBeforeLoop_+resAdjustmentBeforeLoop_;
core::Size tmpResidueAfterLoop = resBeforeLoop_+insertedBeforeLoopRes_+insertedAfterLoopRes_-resAdjustmentAfterLoop_+1;
std::vector< numeric::xyzVector<numeric::Real> > coordinates;
core::Size fragment_length =  SSHashedFragmentStoreOP_->get_fragment_length();
core::Size front_ss_length = 2;
core::Size back_ss_length = 2;
if(fragment_length>(front_ss_length+back_ss_length+loopLength_))//within tolerance add a residue to the front
front_ss_length+=1;
if(fragment_length>(front_ss_length+back_ss_length+loopLength_))//still within tolerance add a residue to the back
back_ss_length+=1;
vector1<core::Size> residues;
for(int ii=front_ss_length-1; ii>=0; --ii){
residues.push_back(tmpResidueBeforeLoop-ii);
}
for(core::Size ii=0; ii<back_ss_length; ++ii){
residues.push_back(tmpResidueAfterLoop+ii);
}
for ( core::Size ii=1; ii<=residues.size(); ++ii ) {
for ( std::string const & atom_name : SSHashedFragmentStoreOP_->get_fragment_store()->fragment_specification.fragment_atoms ) {
coordinates.push_back(fullLengthPoseOP_->residue(residues[ii]).xyz(atom_name));
}
}
SSHashedFragmentStoreOP_->lookback_uncached_stub(coordinates,stub_ss_index_match_,loopLength_,uncached_stub_rmsd_,uncached_stub_index_);
}
*/

void PossibleLoop::setup_finalPose_copy_labels(){
	vector1 < std::string > initial_labels;
	for ( core::Size ii=1; ii<=resBeforeLoop_; ++ii ) {
		initial_labels = original_atom_type_fullLengthPoseOP_->pdb_info()->get_reslabels(ii);
		for ( core::Size jj=1; jj <= initial_labels.size(); ++jj ) {
			finalPoseOP_->pdb_info()->add_reslabel(ii, initial_labels[jj]);
		}
	}
	//fix the residues after the loop
	core::Size residue_after_loop_orig = fullLength_resBeforeLoop_+insertedBeforeLoopRes_+insertedAfterLoopRes_-resAdjustmentAfterLoop_+1;
	//std::cout << "fullLength_resBeforeLoop_" << fullLength_resBeforeLoop_ << "," <<  insertedBeforeLoopRes_ << "," << insertedAfterLoopRes_ << "," << resAdjustmentAfterLoop_ << std::endl;
	core::Size residue_to_change = finalPoseOP_->total_residue()-resAfterLoop_+1;
	for ( core::Size ii=0; ii<residue_to_change; ++ii ) {
		core::Size orig_res_numb =  residue_after_loop_orig+ii;
		core::Size current_res_numb = resAfterLoop_+ii;
		initial_labels = original_atom_type_fullLengthPoseOP_->pdb_info()->get_reslabels(orig_res_numb);
		for ( core::Size jj=1; jj <= initial_labels.size(); ++jj ) {
			finalPoseOP_->pdb_info()->add_reslabel(current_res_numb, initial_labels[jj]);
		}
	}
}

void PossibleLoop::setup_finalPose_copy_rotamers(){
	bool fullatom_original = original_atom_type_fullLengthPoseOP_->is_fullatom();
	bool fullatom_finalPose = finalPoseOP_->is_fullatom();
	if ( !fullatom_original || (fullatom_original && fullatom_finalPose) ) {
		return;
	} else {
		core::util::switch_to_residue_type_set(*finalPoseOP_, core::chemical::FA_STANDARD);
		//fix the residues before the loop
		for ( core::Size ii=1; ii<=resBeforeLoop_; ++ii ) {
			finalPoseOP_->replace_residue(ii, original_atom_type_fullLengthPoseOP_->residue(ii), true );
		}
		//fix the residues after the loop
		core::Size residue_after_loop_orig = fullLength_resBeforeLoop_+insertedBeforeLoopRes_+insertedAfterLoopRes_-resAdjustmentAfterLoop_+1;
		//std::cout << "fullLength_resBeforeLoop_" << fullLength_resBeforeLoop_ << "," <<  insertedBeforeLoopRes_ << "," << insertedAfterLoopRes_ << "," << resAdjustmentAfterLoop_ << std::endl;
		core::Size residue_to_change = finalPoseOP_->total_residue()-resAfterLoop_+1;
		for ( core::Size ii=0; ii<residue_to_change; ++ii ) {
			core::Size orig_res_numb =  residue_after_loop_orig+ii;
			core::Size current_res_numb = resAfterLoop_+ii;
			finalPoseOP_->replace_residue(current_res_numb, original_atom_type_fullLengthPoseOP_->residue(orig_res_numb), true );
		}
	}
}

core::pose::PoseOP PossibleLoop::get_finalPoseOP(){
	return(finalPoseOP_);
}

void PossibleLoop::label_loop(std::string label){
	if ( label != "" ) {
		Size resBeforeLoop_adjusted_for_insert=resBeforeLoop_;
		if ( resAdjustmentBeforeLoop_>0 ) {
			resBeforeLoop_adjusted_for_insert=resBeforeLoop_-resAdjustmentBeforeLoop_;
		}
		Size resAfterLoop_adjusted_for_insert=resAfterLoop_;
		if ( resAdjustmentAfterLoop_>0 ) {
			resAfterLoop_adjusted_for_insert=resAfterLoop_+resAdjustmentAfterLoop_;
		}
		for ( Size ii=resBeforeLoop_adjusted_for_insert+1; ii<=resAfterLoop_adjusted_for_insert-1; ++ii ) {
			finalPoseOP_->pdb_info()->add_reslabel(ii, label);
		}
	}
}

void PossibleLoop::generate_output_pose(bool output_closed,bool ideal_loop, Real rms_threshold,std::string allowed_loop_abegos, std::string closure_type){
	using namespace core::chemical;
	//step 1 : generate trim pose.
	//core::Size tmpResidueBeforeLoop = resBeforeLoop_+resAdjustmentBeforeLoop_;
	//core::Size tmpResidueAfterLoop = resAfterLoop_+insertedBeforeLoopRes_+insertedAfterLoopRes_+resAdjustmentAfterLoop_;
	//trim residues assume the residues after the loop are trimmed first
	core::Size trim_res_start_after_loop =resBeforeLoop_+insertedBeforeLoopRes_+1;
	core::Size trim_res_stop_after_loop =resBeforeLoop_+insertedBeforeLoopRes_+insertedAfterLoopRes_-resAdjustmentAfterLoop_;
	core::Size trim_res_start_before_loop = resBeforeLoop_+resAdjustmentBeforeLoop_+1;//want this there. Should be 1 after this residue
	core::Size trim_res_stop_before_loop = resBeforeLoop_+insertedBeforeLoopRes_;//last residue checks out
	//trim should be 401-401 on right side
	core::pose::PoseOP working_poseOP = fullLengthPoseOP_->clone();
	//remember these could change
	if ( trim_res_start_after_loop <= trim_res_stop_after_loop ) {
		trimRegion(working_poseOP,trim_res_start_after_loop,trim_res_stop_after_loop);
	}
	if ( trim_res_start_before_loop <= trim_res_stop_before_loop ) {
		trimRegion(working_poseOP,trim_res_start_before_loop,trim_res_stop_before_loop);
	}
	resBeforeLoop_=resBeforeLoop_+resAdjustmentBeforeLoop_;
	resAfterLoop_=resBeforeLoop_+1;
	if ( !output_closed ) {
		finalPoseOP_=working_poseOP;
	} else {
		SSHashedFragmentStoreOP_ = protocols::indexed_structure_store::FragmentStoreManager::get_instance()->SSHashedFragmentStore(fragment_store_path_,fragment_store_format_,fragment_store_compression_);
		core::Size fragment_length = SSHashedFragmentStoreOP_->get_fragment_length();
		if ( closure_type == "kic" ) {
			core::scoring::ScoreFunctionOP scorefxn_tmp( core::scoring::ScoreFunctionFactory::create_score_function("score3"));
			extendRegion(true, resBeforeLoop_, loopLength_,working_poseOP);
			resAfterLoop_=resAfterLoop_+loopLength_;
			bool success = kic_closure(scorefxn_tmp,working_poseOP,resBeforeLoop_,resAfterLoop_,200);
			vector1<core::Size> resids;
			for ( int ii=resBeforeLoop_-3; ii<=(int)resBeforeLoop_+3; ii=ii+2 ) {
				core::Size tmp_resid = get_valid_resid(ii,*working_poseOP);
				resids.push_back(tmp_resid);
			}
			Real rmsd = SSHashedFragmentStoreOP_->max_rmsd_in_region(*working_poseOP,resids);
			if ( success ) {
				final_rmsd_=rmsd;
				finalPoseOP_=working_poseOP;
			}
		}
		if ( closure_type == "lookback" ) {
			//increase the length of the loop once before the closure occurs
			resAfterLoop_=resAfterLoop_+loopLength_;
			//Step 1: Sort the stubs
			std::sort(stubVector_.begin(),stubVector_.end(),indexed_structure_store::BackboneStubVectorRMSDComparator());
			Real lowest_loop_rmsd = 99999;
			core::pose::PoseOP ref_poseOP = working_poseOP->clone();
			core::Size MAX_FAILURES=3;
			core::Size failures = 0;
			for ( Size ll=1; (ll<=numb_stubs_to_consider_) && (ll <=stubVector_.size()) && (failures<MAX_FAILURES); ++ll ) {
				if ( ll!=1 ) {
					//workingPose == ref_poseOP the first iteration
					working_poseOP = ref_poseOP->clone();
				}
				Size stub_ss_index_match = stubVector_[ll].ss_index_match;
				Size stub_index_match =stubVector_[ll].index_match;
				//Step 2 : Add loop residues
				extendRegion(true, resBeforeLoop_, loopLength_,working_poseOP);
				//Step 3 : Generate coordinate constraints
				core::Size front_overlap_res_length;
				core::Size back_overlap_res_length;
				generate_overlap_range(front_overlap_res_length, back_overlap_res_length);
				working_poseOP->remove_constraints();
				add_coordinate_csts_from_lookback(stub_ss_index_match,stub_index_match,resBeforeLoop_-front_overlap_res_length+1,true,front_overlap_res_length,working_poseOP);
				add_dihedral_csts_from_lookback(stub_ss_index_match,stub_index_match,resBeforeLoop_-front_overlap_res_length+1,working_poseOP);
				//Step 3 : Generate score function
				core::scoring::ScoreFunctionOP scorefxn( core::scoring::ScoreFunctionFactory::create_score_function("score3") );
				scorefxn->set_weight(core::scoring::coordinate_constraint, 0.7);
				scorefxn->set_weight(core::scoring::dihedral_constraint, 20.0 );
				scorefxn->set_weight(core::scoring::rama,1.0);
				scorefxn->set_weight(core::scoring::cart_bonded, 1.0 );
				scorefxn->set_weight(core::scoring::cart_bonded_length,0.2);
				scorefxn->set_weight(core::scoring::cart_bonded_angle,0.5);
				scorefxn->set_weight(core::scoring::vdw, 2 );
				scorefxn->set_weight(core::scoring::linear_chainbreak,5.0);
				//step 3: assign_phi_psi
				assign_phi_psi_omega_from_lookback(stub_ss_index_match,stub_index_match, working_poseOP,front_overlap_res_length);
				//step 4: initial minimize
				minimize_loop(scorefxn,ideal_loop,working_poseOP);
				//Get rmsd
				vector1<core::Size> resids;
				for ( int ii=resBeforeLoop_-3; ii<=(int)resBeforeLoop_+3; ii=ii+1 ) {
					core::Size tmp_resid = get_valid_resid(ii,*working_poseOP);
					resids.push_back(tmp_resid);
				}
				Real current_rmsd = SSHashedFragmentStoreOP_->max_rmsd_in_region(*working_poseOP,resids);

				//step 5: iterative minimize over multiple lookbacks
				core::sequence::SSManager SM;
				std::string fragStore_ss_string = SM.index2symbolString(stub_ss_index_match,fragment_length);
				core::scoring::dssp::Dssp frag_dssp( *working_poseOP);
				frag_dssp.dssp_reduced();
				std::string tmp_dssp = frag_dssp.get_dssp_secstruct();
				std::string desired_dssp = tmp_dssp.substr(0,resBeforeLoop_-1-1)+fragStore_ss_string+tmp_dssp.substr(resBeforeLoop_+fragment_length-1-1,tmp_dssp.length());
				Real match_rmsd;
				core::Size match_index;
				core::Size match_ss_index;
				scorefxn->set_weight(core::scoring::coordinate_constraint, 0.4);
				//with the improvements the code below results in very marginal improvements. So it has been eliminated
				bool last_round_rmsd_reduction = true;
				for ( core::Size jj=1; jj<=2 && current_rmsd<rms_threshold+0.2 && last_round_rmsd_reduction; ++jj ) {
					core::pose::PoseOP tmp_poseOP = working_poseOP->clone();
					tmp_poseOP->remove_constraints();
					std::string frag_ss = desired_dssp.substr(resBeforeLoop_-front_overlap_res_length+1-1,9);
					SSHashedFragmentStoreOP_->lookback_account_for_dssp_inaccuracy(*tmp_poseOP,resBeforeLoop_-front_overlap_res_length+1,frag_ss,match_rmsd,match_index,match_ss_index);
					add_coordinate_csts_from_lookback(match_ss_index,match_index,resBeforeLoop_-front_overlap_res_length+1,false,front_overlap_res_length,tmp_poseOP);
					minimize_loop(scorefxn,ideal_loop,tmp_poseOP);
					Real redesign_rmsd = SSHashedFragmentStoreOP_->max_rmsd_in_region(*tmp_poseOP,resids);
					if ( redesign_rmsd<current_rmsd ) {
						current_rmsd = redesign_rmsd;
						working_poseOP = tmp_poseOP;
					} else {
						last_round_rmsd_reduction=false;
					}
				}
				//Step 5: Kic to initially close loops if ideal
				if ( ideal_loop ) {
					//uses kic to close the last couple residues in the loop
					core::Size res_from_end_of_loop=3;
					core::Size firstLoopRes= resAfterLoop_-res_from_end_of_loop;
					core::Size lastLoopRes=resAfterLoop_;
					bool success = kic_closure(scorefxn,working_poseOP,firstLoopRes,lastLoopRes,10);
					Real rmsd = SSHashedFragmentStoreOP_->max_rmsd_in_region(*working_poseOP,resids);
					if ( success ) {
						current_rmsd=rmsd;
					} else {
						current_rmsd = 999;
					}
				}
				if ( allowed_loop_abegos!="" ) {
					bool valid_loop_abego = check_loop_abego(working_poseOP,resBeforeLoop_,resAfterLoop_,allowed_loop_abegos,current_rmsd);
					if ( !valid_loop_abego ) {
						current_rmsd=999;
					}
				}
				if ( current_rmsd >rms_threshold ) {
					failures++;
				}
				if ( lowest_loop_rmsd>current_rmsd ) {
					lowest_loop_rmsd = current_rmsd;
					final_rmsd_ = lowest_loop_rmsd;
					finalPoseOP_ = working_poseOP->clone();
				}
			}

		}
	}
}

bool PossibleLoop::check_loop_abego(core::pose::PoseOP & poseOP, core::Size resBeforeLoop, core::Size resAfterLoop, std::string allowed_loop_abegos, Real current_rmsd ){
	utility::vector1< std::string > const abego_v = core::sequence::get_abego( *poseOP );
	core::Size neighboring_residues = 3;
	core::sequence::ABEGOManager am;
	std::string abego = am.get_abego_string(abego_v);
	core::Size loop_cut_length = resAfterLoop+neighboring_residues-(resBeforeLoop-neighboring_residues);
	std::string loop_abego = abego.substr(resBeforeLoop-neighboring_residues-1,loop_cut_length);//1 is for the string conversion
	utility::vector1< std::string > allowed_abegos_v( utility::string_split( allowed_loop_abegos , ',' ) );
	bool found = false;
	for ( auto & allowed_abego : allowed_abegos_v ) {
		core::Size position = loop_abego.find(allowed_abego);
		if ( position != std::string::npos ) {
			found=true;
			TR << "ABEGO loop found:" << allowed_abego << " with RMSD" << current_rmsd << std::endl;
		}
	}
	return(found);
}


void PossibleLoop::assign_phi_psi_omega_from_lookback(core::Size db_index, core::Size fragment_index, core::pose::PoseOP & poseOP,core::Size front_overlap_res_length){
	vector<Real> phi_v = SSHashedFragmentStoreOP_->get_fragment_store(db_index)->realVector_groups["phi"][fragment_index];
	vector<Real> psi_v = SSHashedFragmentStoreOP_->get_fragment_store(db_index)->realVector_groups["psi"][fragment_index];
	vector<Real> omega_v = SSHashedFragmentStoreOP_->get_fragment_store(db_index)->realVector_groups["omega"][fragment_index];
	//to keep the structure ideal I don't assign omega
	kinematics::FoldTree ft;
	ft = poseOP->fold_tree();
	ft.clear();
	ft.add_edge(1,resBeforeLoop_+loopLength_,core::kinematics::Edge::PEPTIDE);
	ft.add_edge(1,poseOP->size(),1); //edges need to be added in order so the fold tree is completely connected.
	ft.add_edge(poseOP->size(),resAfterLoop_,core::kinematics::Edge::PEPTIDE);
	poseOP->fold_tree(ft);
	for ( core::Size ii=0; ii<loopLength_; ++ii ) {//set phi psi for residues after loop
		poseOP->set_phi(resBeforeLoop_+ii,phi_v[ii+front_overlap_res_length-1]);
		poseOP->set_psi(resBeforeLoop_+ii,psi_v[ii+front_overlap_res_length-1]);
		poseOP->set_omega(resBeforeLoop_+ii,omega_v[ii+front_overlap_res_length-1]);
	}
	//at loop end set phi but not psi/omega
	poseOP->set_phi(resBeforeLoop_+loopLength_,phi_v[loopLength_+front_overlap_res_length-1]);
	//at first residue post loop
	poseOP->set_psi(resBeforeLoop_+loopLength_+1,psi_v[loopLength_+front_overlap_res_length-1+1]);
	poseOP->set_omega(resBeforeLoop_+loopLength_+1,omega_v[loopLength_+front_overlap_res_length-1+1]);
}


vector<Real> PossibleLoop::get_center_of_mass(const Real* coordinates, int number_of_atoms){
	vector<Real> center;
	center.push_back(0);
	center.push_back(0);
	center.push_back(0);

	for ( int n = 0; n < number_of_atoms * 3; n += 3 ) {
		center[0] += coordinates[n + 0];
		center[1] += coordinates[n + 1];
		center[2] += coordinates[n + 2];
	}

	center[0] /= number_of_atoms;
	center[1] /= number_of_atoms;
	center[2] /= number_of_atoms;

	return(center);
}

void PossibleLoop::output_fragment_debug(std::vector< numeric::xyzVector<numeric::Real> > coordinates,std::string filename){
	using namespace ObjexxFCL::format;
	utility::io::ozstream out(filename, std::ios_base::app);
	core::Size resid = 1;
	for ( auto & coordinate : coordinates ) {
		out << "ATOM  " << I(5,resid) << "  CA  " <<
			"GLY" << ' ' << 'A' << I(4,resid ) << "    " <<
			F(8,3,coordinate.x()) <<
			F(8,3,coordinate.y()) <<
			F(8,3,coordinate.z()) <<
			F(6,2,1.0) << F(6,2,1.0) << '\n';
		resid++;
	}
	out << "ATOM  " << I(5,resid) << "  CA  " <<
		"GLY" << ' ' << 'A' << I(4,resid) << "    " <<
		F(8,3,0.0) <<
		F(8,3,0.0) <<
		F(8,3,0.0) <<
		F(6,2,1.0) << F(6,2,1.0) << '\n';
	resid++;
	out << "ENDMDL\n";
	out.close();
}



void PossibleLoop::add_coordinate_csts_from_lookback(core::Size stub_ss_index_match, core::Size fragment_index, core::Size pose_residue, bool match_stub_alone, core::Size front_overlap_res_length,core::pose::PoseOP & poseOP){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace chemical;
	using namespace protocols::generalized_kinematic_closure;

	typedef numeric::alignment::QCPKernel< numeric::Real > QCPKernel;
	core::Size fragment_length =  SSHashedFragmentStoreOP_->get_fragment_length();
	std::vector< numeric::xyzVector<numeric::Real> > fragCoordinates = SSHashedFragmentStoreOP_->get_fragment_coordinates(stub_ss_index_match,fragment_index);

	std::vector< numeric::xyzVector<numeric::Real> > fragCoordinates_rot(fragment_length);

	if ( !match_stub_alone ) {
		//full length fragment
		std::vector< numeric::xyzVector<numeric::Real> > coordinates;
		for ( core::Size ii = 0;  ii < fragment_length; ++ii ) {
			coordinates.push_back(poseOP->residue(pose_residue+ii).xyz("CA"));
		}

		Eigen::Transform<numeric::Real, 3, Eigen::Affine> superposition_transform;

		QCPKernel::CoordMap frag_cmap(&fragCoordinates.front().x(), 3, fragCoordinates.size());
		QCPKernel::CoordMap frag_rot_cmap(&fragCoordinates_rot.front().x(), 3, fragCoordinates_rot.size());
		QCPKernel::CoordMap coor_cmap(&coordinates.front().x(), 3, coordinates.size());

		Real rmsd = QCPKernel::calc_coordinate_superposition(frag_cmap, coor_cmap, superposition_transform);
		TR.Debug << "full fragment rmsd" << rmsd << std::endl;

		frag_rot_cmap = superposition_transform * frag_cmap;
	} else { //------Prepare stub match-----------------------------------------------------
		//1. Get coordinates coordintes
		// coordinates
		// fragCoordinates
		// coordinates_stub
		// coordinates_removed_com_stub
		// fragCoordinate_removed_com_stub
		//2. Remove COM from coordinates_removed_com_stub,coordinates_removed_com_stub
		//3. Generate rotMatrix between coordinates_removed_com_stub, coordinates_removed_com_stub
		//4. Generate centers of mass for coordinates_stub call this coordinates_stub_com
		//5. Generate removed_com_vectors using the cooridinates_stub_com called coordinates_removed_stub_com
		//6. apply rotMatrix to fragCoordinates store in fragCoordinates_rot
		//7. Combine data fragCoordinates_rot[ii].x()= coordinates[ii].x()-coordinates_removed_com_stub[ii].x()+fragCoordinates_rot_tmp[ii].x()
		std::vector< numeric::xyzVector<numeric::Real> > coordinates;
		std::vector< numeric::xyzVector<numeric::Real> > coordinates_stub;
		std::vector< numeric::xyzVector<numeric::Real> > fragCoordinates_stub;
		core::Size res_ct = 0;

		for ( core::Size ii = 0;  ii < fragment_length; ++ii ) {
			coordinates.push_back(poseOP->residue(pose_residue+ii).xyz("CA"));

			if ( ii<front_overlap_res_length || ii>=front_overlap_res_length+loopLength_ ) {
				res_ct++;
				coordinates_stub.push_back(poseOP->residue(pose_residue+ii).xyz("CA"));
				fragCoordinates_stub.push_back(fragCoordinates[ii]);
			}
		}

		Eigen::Transform<numeric::Real, 3, Eigen::Affine> superposition_transform;
		Real starting_rmsd = QCPKernel::calc_coordinate_superposition(
			QCPKernel::CoordMap(&fragCoordinates_stub.front().x(), 3, res_ct),
			QCPKernel::CoordMap(&coordinates_stub.front().x(), 3, res_ct),
			superposition_transform
		);
		TR.Debug << "stub fragment rmsd" << starting_rmsd << std::endl;

		QCPKernel::CoordMap frag_cmap(&fragCoordinates.front().x(), 3, fragCoordinates.size());
		QCPKernel::CoordMap frag_rot_cmap(&fragCoordinates_rot.front().x(), 3, fragCoordinates.size());

		frag_rot_cmap = superposition_transform * frag_cmap;

		//output_fragment_debug(fragCoordinates_rot,"X3_fragCoords_rotate.pdb");
		//In test simulations rmsd_v3 was 8.02771 and rmsd_v4 was 8.27703. The only difference is the removing of the COM of the coordintes. But I don't see any errors
		//I'm ignoring this bug for now.  Maybe forever.
		//I think I found this bug. It was due to fragCoordinates_rot not being defined
	}

	ConstraintCOPs csts;
	core::id::AtomID const anchor_atom( core::id::AtomID( poseOP->residue(1).atom_index("CA"), 1) );
	//core::scoring::func::HarmonicFuncOP coord_cst_func = utility::pointer::make_shared< core::scoring::func::HarmonicFunc >( 0.1, 1.0 );
	core::scoring::func::HarmonicFuncOP coord_cst_func = utility::pointer::make_shared< core::scoring::func::HarmonicFunc >( 0.0, 0.2 );
	for ( core::Size ii = 0;  ii < fragment_length; ++ii ) {
		core::Size atomindex =  poseOP->residue( pose_residue ).atom_index( "CA" );
		csts.push_back( utility::pointer::make_shared< CoordinateConstraint >(core::id::AtomID( atomindex, pose_residue+ii),anchor_atom,fragCoordinates_rot[ii],coord_cst_func));
	}
	poseOP->add_constraints( csts );

}

void PossibleLoop::add_dihedral_csts_from_lookback(core::Size stub_ss_index_match,core::Size fragment_index,core::Size pose_residue, core::pose::PoseOP & poseOP){
	using namespace core::scoring::func;
	using numeric::conversions::radians;
	using namespace core::scoring::constraints;
	using namespace core::id;
	core::Size fragment_length = SSHashedFragmentStoreOP_->get_fragment_length();
	vector<Real> phi_v = SSHashedFragmentStoreOP_->get_fragment_store(stub_ss_index_match)->realVector_groups["phi"][fragment_index];
	vector<Real> psi_v = SSHashedFragmentStoreOP_->get_fragment_store(stub_ss_index_match)->realVector_groups["psi"][fragment_index];
	vector<Real> omega_v = SSHashedFragmentStoreOP_->get_fragment_store(stub_ss_index_match)->realVector_groups["omega"][fragment_index];
	core::Real phi_sd_deg = 30; //40
	core::Real psi_sd_deg = 30; //40
	core::Real omega_sd_deg=10; //10
	core::Real phi_sd_rad = numeric::conversions::radians(phi_sd_deg);
	core::Real psi_sd_rad = numeric::conversions::radians(psi_sd_deg);
	core::Real omega_sd_rad = numeric::conversions::radians(omega_sd_deg);
	for ( core::Size ii=pose_residue; ii<pose_residue+fragment_length; ++ii ) {
		AtomID c_0(poseOP->residue(ii-1).atom_index( "C" ),ii-1);
		AtomID n_1(poseOP->residue( ii ).atom_index( "N" ),ii);
		AtomID ca_1( poseOP->residue( ii ).atom_index( "CA" ),ii);
		AtomID c_1( poseOP->residue(ii).atom_index( "C" ), ii );
		AtomID n_2( poseOP->residue(ii+1).atom_index( "N" ),ii+1);
		AtomID ca_2( poseOP->residue(ii+1).atom_index( "CA" ),ii+1);
		// std::cout << "phi" << ii << "phi_v[ii-resBeforeLoop_+1]" << phi_v[ii-resBeforeLoop_+1] << std::endl;
		// std::cout << "psi" << ii << "psi_v[ii-resBeforeLoop_+1]" << psi_v[ii-resBeforeLoop_+1] << std::endl;
		// std::cout << "omega" << ii << "omega_v[ii-resBeforeLoop_+1]" << omega_v[ii-resBeforeLoop_+1] << std::endl;
		Real phi_radians = radians(phi_v[ii-pose_residue]);
		//std::cout << "phi" << phi_v[ii-resBeforeLoop_+1]  <<"," << phi_radians << std::endl;
		CircularHarmonicFuncOP phi_func(utility::pointer::make_shared<CircularHarmonicFunc>( phi_radians, phi_sd_rad ) );
		ConstraintOP phi_cst( utility::pointer::make_shared<DihedralConstraint>(
			c_0,n_1,ca_1,c_1, phi_func ) );
		poseOP->add_constraint( scoring::constraints::ConstraintCOP( phi_cst ) );
		//psi-----
		Real psi_radians = radians(psi_v[ii-pose_residue]);
		//std::cout << "psi" << psi_v[ii-resBeforeLoop_+1]  <<"," << psi_radians << std::endl;
		CircularHarmonicFuncOP psi_func(utility::pointer::make_shared<CircularHarmonicFunc>( psi_radians, psi_sd_rad) );
		ConstraintOP psi_cst( utility::pointer::make_shared<DihedralConstraint>(n_1,ca_1,c_1,n_2, psi_func  ) );
		poseOP->add_constraint( scoring::constraints::ConstraintCOP( psi_cst ) );
		//omega-----
		Real omega_radians = radians(omega_v[ii-pose_residue]);
		//std::cout << "omega" << omega_v[ii-resBeforeLoop_+1]  <<"," << omega_radians << std::endl;
		CircularHarmonicFuncOP omega_func(utility::pointer::make_shared<CircularHarmonicFunc>( omega_radians, omega_sd_rad) );
		ConstraintOP omega_cst( utility::pointer::make_shared<DihedralConstraint>(ca_1,c_1,n_2,ca_2, omega_func ) );
		poseOP->add_constraint( scoring::constraints::ConstraintCOP( omega_cst ) );
	}
}


Size PossibleLoop::get_valid_resid(int resid,core::pose::Pose const pose){
	if ( resid<1 ) {
		TR << "invalid resid encountered n-term" << resid << std::endl;
		resid = 1;
	}
	if ( resid+9-1>(int)pose.size() ) {
		TR << "invalid resid encountered c-term" << resid << std::endl;
		resid = (int)pose.size()-9+1;
	}
	return(core::Size(resid));
}

std::string PossibleLoop::get_description(){
	std::stringstream extraTagInfoString;
	extraTagInfoString << "S" << resBeforeLoop_ << "E" << resAfterLoop_ << "Sa"<< resAdjustmentBeforeLoop_ << "Ea" << resAdjustmentAfterLoop_ << "LL" << loopLength_;
	return(extraTagInfoString.str());
}

std::vector< numeric::xyzVector<numeric::Real> > PossibleLoop::get_coordinates_from_pose(core::pose::PoseOP const poseOP,core::Size resid,core::Size length){
	std::vector< numeric::xyzVector<numeric::Real> > coordinates;
	for ( core::Size ii =resid;  ii < resid+length; ++ii ) {
		coordinates.push_back(poseOP->residue(ii).xyz("CA"));
	}
	return(coordinates);
}

// Real PossibleLoop::get_vdw_change(core::pose::PoseOP poseOP){
//  core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
//  scorefxn->set_weight( core::scoring::vdw, 1.0 );
//  Real loop_vdw_score = scorefxn->score(*poseOP);
//  Real orig_vdw_score = scorefxn->score(*original_poseOP_);
//  return(loop_vdw_score-orig_vdw_score);
// }



Real PossibleLoop::rmsd_between_coordinates(std::vector< numeric::xyzVector<numeric::Real> > fragCoordinates,std::vector< numeric::xyzVector<numeric::Real> > coordinates){
	typedef numeric::alignment::QCPKernel< numeric::Real > QCPKernel;

	return QCPKernel::calc_coordinate_rmsd(
		QCPKernel::CoordMap(&fragCoordinates.front().x(), 3, fragCoordinates.size()),
		QCPKernel::CoordMap(&coordinates.front().x(), 3, coordinates.size()));
}

//closes the loop robustly. However, the phi/psi omega fall into the wrong abego types.
bool PossibleLoop::kic_closure(core::scoring::ScoreFunctionOP scorefxn,core::pose::PoseOP & poseOP,core::Size firstLoopRes, core::Size lastLoopRes,core::Size numb_kic_cycles){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace chemical;
	using namespace protocols::generalized_kinematic_closure;
	GeneralizedKICOP genKIC(utility::pointer::make_shared<GeneralizedKIC>());
	kinematics::FoldTree ft;
	ft = poseOP->fold_tree();
	ft.clear();
	ft.add_edge(1,resBeforeLoop_+loopLength_,core::kinematics::Edge::PEPTIDE);
	ft.add_edge(1,poseOP->size(),1); //edges need to be added in order so the fold tree is completely connected.
	ft.add_edge(poseOP->size(),resAfterLoop_,core::kinematics::Edge::PEPTIDE);
	poseOP->fold_tree(ft);
	genKIC->add_filter( "loop_bump_check" );
	//hold first res to the correct ABEGO
	scorefxn->score(*poseOP);
	for ( core::Size ii=firstLoopRes; ii<=lastLoopRes; ii++ ) {
		genKIC->add_loop_residue( ii );
	}
	genKIC->add_perturber( "randomize_alpha_backbone_by_rama" );
	for ( core::Size ii=firstLoopRes; ii<=lastLoopRes; ii++ ) {
		genKIC->add_residue_to_perturber_residue_list( ii );
	}
	genKIC->add_perturber( "set_dihedral" );
	utility::vector1 < core::id::NamedAtomID > atomset;
	atomset.push_back( core::id::NamedAtomID( "C", lastLoopRes-1 ) );
	atomset.push_back( core::id::NamedAtomID( "N", lastLoopRes ) );
	genKIC->add_atomset_to_perturber_atomset_list( atomset );
	genKIC->add_value_to_perturber_value_list(180.0);
	core::Size midLoop = ((lastLoopRes-firstLoopRes)/2)+firstLoopRes;
	genKIC->set_pivot_atoms( firstLoopRes, "CA", midLoop, "CA" , lastLoopRes, "CA" );
	genKIC->set_closure_attempts(numb_kic_cycles);
	genKIC->set_selector_type("lowest_energy_selector");
	genKIC->set_selector_scorefunction(scorefxn);
	genKIC->close_bond(lastLoopRes-1,"C",lastLoopRes,"N",lastLoopRes-1,"CA",lastLoopRes,"CA",1.32829/*bondlength*/,116.2/*bondangle1*/,121.7/*bondagle2*/,180/*torsion*/,false,false);
	poseOP->conformation().reset_chain_endings();
	genKIC->apply(*poseOP);
	bool kicSuccess = genKIC->last_run_successful();
	if ( kicSuccess ) {
		ft = poseOP->fold_tree();
		ft.clear();
		ft.add_edge(1,poseOP->size(),core::kinematics::Edge::PEPTIDE);
		poseOP->fold_tree(ft);
		//poseOP->conformation().set_polymeric_connection(lastLoopRes-1,lastLoopRes);
		//poseOP->conformation().update_polymeric_connection(lastLoopRes-1,true);
		return true;
	}
	return false;
}


void PossibleLoop::minimize_loop(core::scoring::ScoreFunctionOP scorefxn,bool ideal_loop,core::pose::PoseOP & poseOP){
	using namespace core::optimization;
	using namespace chemical;
	kinematics::FoldTree ft;
	ft = poseOP->fold_tree();
	ft.clear();
	ft.add_edge(1,resBeforeLoop_+loopLength_,core::kinematics::Edge::PEPTIDE);
	ft.add_edge(1,poseOP->size(),1); //edges need to be added in order so the fold tree is completely connected.
	ft.add_edge(poseOP->size(),resAfterLoop_,core::kinematics::Edge::PEPTIDE);
	poseOP->fold_tree(ft);
	core::Size firstLoopRes=resBeforeLoop_;
	core::Size lastLoopRes=resAfterLoop_;
	core::kinematics::MoveMap mm;
	mm.set_bb  ( false );
	mm.set_chi ( false );
	mm.set_jump( false );
	mm.set_bb_true_range(firstLoopRes,lastLoopRes);
	//testing code
	Size cutpoint = resBeforeLoop_+loopLength_;
	Size res_after_cutpoint = resAfterLoop_;
	core::pose::add_variant_type_to_pose_residue( *poseOP, CUTPOINT_LOWER,cutpoint);
	core::pose::add_variant_type_to_pose_residue( *poseOP, CUTPOINT_UPPER,res_after_cutpoint);
	//poseOP->conformation().declare_chemical_bond( cutpoint, poseOP->residue( cutpoint ).atom_name( poseOP->residue( cutpoint ).upper_connect_atom() ),
	// res_after_cutpoint , poseOP->residue( cutpoint ).atom_name( poseOP->residue( cutpoint ).lower_connect_atom() ) );
	//end testing code
	if ( ideal_loop ) {
		MinimizerOptions min_options("lbfgs_armijo_nonmonotone", 0.01, true, false, false);
		min_options.silent(true);
		AtomTreeMinimizer minimizer;

		minimizer.run( *poseOP, mm,*scorefxn, min_options );
		//std::cout << "after atom tree minimize" << std::endl;
		//scorefxn->show(std::cout,*poseOP);
	} else { //use cartesian minimization if it's ok.
		MinimizerOptions min_options("lbfgs_armijo_nonmonotone", 0.01, true, false, false);
		CartesianMinimizer minimizer;
		min_options.silent(true);
		scorefxn->score(*poseOP);
		minimizer.run( *poseOP, mm,*scorefxn, min_options );
	}
	//testing code
	core::pose::remove_variant_type_from_pose_residue( *poseOP, CUTPOINT_LOWER, cutpoint   );
	core::pose::remove_variant_type_from_pose_residue( *poseOP, CUTPOINT_UPPER, res_after_cutpoint );
	poseOP->conformation().declare_chemical_bond( cutpoint, poseOP->residue( cutpoint ).atom_name( poseOP->residue( cutpoint ).upper_connect_atom() ),
		res_after_cutpoint , poseOP->residue( cutpoint ).atom_name( poseOP->residue( cutpoint ).lower_connect_atom() ) );

	ft = poseOP->fold_tree();
	ft.clear();
	ft.add_edge(1,poseOP->size(),core::kinematics::Edge::PEPTIDE);
	poseOP->fold_tree(ft);
}


NearNativeLoopCloser::NearNativeLoopCloser():moves::Mover("NearNativeLoopCloser"){
}






NearNativeLoopCloser::NearNativeLoopCloser(int resAdjustmentStartLow,int resAdjustmentStartHigh,int resAdjustmentStopLow,int resAdjustmentStopHigh,int resAdjustmentStartLow_sheet,int resAdjustmentStartHigh_sheet,int resAdjustmentStopLow_sheet,int resAdjustmentStopHigh_sheet,core::Size loopLengthRangeLow, core::Size loopLengthRangeHigh,core::Size resBeforeLoop,core::Size resAfterLoop,char chainBeforeLoop, char chainAfterLoop,Real rmsThreshold, Real max_vdw_change, bool idealExtension,bool ideal, bool output_closed,std::string closure_type,std::string allowed_loop_abegos,std::string label_loop, std::string fragment_store_path, std::string fragment_store_format, std::string fragment_store_compression,core::Size numb_stubs_to_consider){
	resAdjustmentStartLow_=resAdjustmentStartLow;
	resAdjustmentStartHigh_=resAdjustmentStartHigh;
	resAdjustmentStopLow_=resAdjustmentStopLow;
	resAdjustmentStopHigh_=resAdjustmentStopHigh;
	resAdjustmentStartLow_sheet_=resAdjustmentStartLow_sheet;
	resAdjustmentStartHigh_sheet_=resAdjustmentStartHigh_sheet;
	resAdjustmentStopLow_sheet_=resAdjustmentStopLow_sheet;
	resAdjustmentStopHigh_sheet_=resAdjustmentStopHigh_sheet;
	loopLengthRangeLow_=loopLengthRangeLow;
	loopLengthRangeHigh_=loopLengthRangeHigh;
	resBeforeLoop_=resBeforeLoop;
	resAfterLoop_=resAfterLoop;
	chainBeforeLoop_=chainBeforeLoop;
	chainAfterLoop_=chainAfterLoop;
	rmsThreshold_=rmsThreshold;
	idealExtension_=idealExtension;
	output_closed_ = output_closed;
	numb_stubs_to_consider_ = numb_stubs_to_consider;
	ideal_=ideal;
	output_all_=false;
	top_outputed_ = false;
	numb_outputed_=0;
	max_number_of_results_ = 999999;
	max_vdw_change_ = max_vdw_change;
	closure_type_ = closure_type;
	allowed_loop_abegos_ = allowed_loop_abegos;
	label_loop_=label_loop;
	TR << "native loop closer init"<< resAdjustmentStartLow_ << "," << resAdjustmentStartHigh_ <<"," << resAdjustmentStopLow_ <<"," << resAdjustmentStopHigh_ << "," << loopLengthRangeLow_ << "," << loopLengthRangeHigh_ << std::endl;
	fragment_store_path_= fragment_store_path;
	fragment_store_format_ = fragment_store_format;
	fragment_store_compression_ = fragment_store_compression;
	SSHashedFragmentStoreOP_ = protocols::indexed_structure_store::FragmentStoreManager::get_instance()->SSHashedFragmentStore(fragment_store_path_,fragment_store_format_,fragment_store_compression_);
	SSHashedFragmentStoreOP_->set_threshold_distance(rmsThreshold_);
}


Real NearNativeLoopCloser::close_loop(core::pose::Pose & pose) {
	//time_t start_time = time(NULL);
	//get name of pose for multi-pose mover-------------------
	pose_name_ = "No_Name_Found";
	if ( pose.pdb_info() && ( pose.pdb_info()->name() != "" ) ) {
		pose_name_ = pose.pdb_info()->name();
	} else if ( pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		pose_name_ = static_cast< basic::datacache::CacheableString const & >( pose.data().get( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
	} else {
		pose_name_ = protocols::jd2::current_input_tag();
	}
	//do actual closer
	//-------deal with multiple chains-for easier indexing-----------
	if ( (chainBeforeLoop_!=chainAfterLoop_) || (pose.chain( resBeforeLoop_) != pose.chain( resAfterLoop_ )) ) { //connecting chains
		if ( !(chainBeforeLoop_!=chainAfterLoop_) ) {
			chainBeforeLoop_ = get_chain_from_chain_id(pose.chain( resBeforeLoop_),pose);
			chainAfterLoop_ =  get_chain_from_chain_id(pose.chain( resAfterLoop_),pose);
		}
		combine_chains(pose);
	} else { //internal loop
		vector1<int> chains = get_chains(pose);
		core::Size chain_id;
		if ( chains.size()==1 ) { //chain id need not be specified.
			chain_id = chains[1];
		} else {
			chain_id =  get_chain_id_from_chain(chainBeforeLoop_,pose);
		}
		if ( chain_id != 1 ) { //renumber chains so chain being modified is in front for ease in numbering
			utility::vector1<core::Size> orig_chain_order = get_chains(pose);
			vector1<Size> new_chain_order;
			for ( core::Size ii=1; ii<=orig_chain_order.size(); ++ii ) {
				if ( orig_chain_order[ii]!=chain_id ) {
					new_chain_order.push_back(orig_chain_order[ii]);
				}
			}
			switch_chain_order(pose,new_chain_order);
		}
	}
	//-------make all possible loops-----------
	possibleLoops_ = create_potential_loops(pose);
	if ( !output_closed_ ) {
		for ( core::Size ii=1; ii<=possibleLoops_.size(); ii++ ) {
			possibleLoops_[ii]->generate_output_pose(false,ideal_,rmsThreshold_,allowed_loop_abegos_,closure_type_);
		}
	} else {
		if ( closure_type_ == "kic" ) {
			for ( core::Size ii=1; ii<=possibleLoops_.size(); ii++ ) {
				possibleLoops_[ii]->generate_output_pose(true,ideal_,rmsThreshold_,allowed_loop_abegos_,closure_type_);
			}
		}
		if ( closure_type_ == "lookback" ) {
			core::Real stubRmsdThreshold=rmsThreshold_+0.10;
			//-------check CA-CA distance viability of loop-----------
			for ( core::Size ii=1; ii<=possibleLoops_.size(); ii++ ) {
				possibleLoops_[ii]->evaluate_distance_closure();
				if ( possibleLoops_[ii]->get_below_distance_threshold() && output_closed_ ) {
					possibleLoops_[ii]->generate_stub_rmsd(stubRmsdThreshold);
				}
			}
			//-------sort loops by stub rmsd-----------
			std::sort(possibleLoops_.begin(), possibleLoops_.end(), BackboneStubRMSDComparator());
			//-------get final output------------
			for ( core::Size ii=1; ii<=possibleLoops_.size(); ++ii ) {
				if ( possibleLoops_[ii]->get_stubRMSD() <stubRmsdThreshold ) {
					possibleLoops_[ii]->generate_output_pose(true,ideal_,rmsThreshold_,allowed_loop_abegos_,closure_type_);
				}
			}
			//debug code--------------------
			// for ( core::Size ii=1; ii<possibleLoops_.size(); ++ii ) {
			//  //std::cout <<"rmsds" << ii << "," <<  possibleLoops_[ii]->get_stubRMSD() << "," << possibleLoops_[ii]->get_uncached_stubRMSD() <<  "," << possibleLoops_[ii]->get_final_RMSD() << std::endl;
			//  //std::cout << possibleLoops_[ii]->get_description();
			//  if( possibleLoops_[ii]->get_uncached_stubRMSD() < 0.4){
			//   std::stringstream numbConvert;
			//   numbConvert << "bust_" << ii << ".pdb";
			//   std::string tmp_pdb_name= numbConvert.str();
			//   possibleLoops_[ii]->get_finalPoseOP()->dump_pdb(tmp_pdb_name);
			//  }
			// }
			//
		}
	}
	Real return_rmsd;
	core::pose::PoseOP tmpPoseOP=get_additional_output_with_rmsd(return_rmsd);
	if ( tmpPoseOP!=nullptr ) {
		pose=*tmpPoseOP;
	}

	return(return_rmsd);
	//time_t end_time = time(NULL);
	//std::cout << "total_time" << end_time-start_time << std::endl;
}



void NearNativeLoopCloser::apply(core::pose::Pose & pose) {
	pose.conformation().clear_parameters_set_list();
	close_loop(pose);
}


void NearNativeLoopCloser::combine_chains(core::pose::Pose & pose){
	if ( !(has_chain(chainBeforeLoop_,pose) && has_chain(chainAfterLoop_,pose)) ) {
		utility_exit_with_message("chains not found in pdb");
	}
	//cache labels by chain_id
	std::map< core::Size, vector1 < vector1 < std::string > >  > orig_labels_by_chainid;
	vector1< core::Size > chains = get_chains(pose);
	for ( core::Size ii=1; ii<=chains.size(); ++ii ) {
		vector1<core::Size> res_in_chain = get_resnums_for_chain_id(pose,chains[ii]);
		vector1<vector1 <std::string> > labels_in_chain;
		for ( core::Size jj=1; jj<=res_in_chain.size(); ++jj ) {
			vector1 < std::string > labels_tmp = pose.pdb_info()->get_reslabels(res_in_chain[jj]);
			labels_in_chain.push_back(labels_tmp);
		}
		orig_labels_by_chainid[chains[ii]]=labels_in_chain;
	}


	//extract chains
	core::Size chain1_id =  get_chain_id_from_chain(chainBeforeLoop_,pose);
	core::pose::PoseOP chain1 = pose.split_by_chain(chain1_id);
	core::Size chain2_id =  get_chain_id_from_chain(chainAfterLoop_,pose);
	core::pose::PoseOP chain2 = pose.split_by_chain(chain2_id);
	//assign_residue_number
	resBeforeLoop_=chain1->size();
	resAfterLoop_=chain1->size()+1;
	//Remove terminal residues & attach
	remove_upper_terminus_type_from_pose_residue(*chain1,chain1->size());
	remove_lower_terminus_type_from_pose_residue(*chain2,1);
	for ( core::Size ii=1; ii<=chain2->size(); ++ii ) {
		chain1->append_residue_by_bond(chain2->residue(ii),false,0,chain1->size());
	}
	//append chain to pose
	append_pose_to_pose(pose,*chain1,true);
	//reorder chains with the last one starting at position 3 and becoming position 1 & the original chains removed.
	//combine generates a pose with chain1,2 and the combined as chain3
	//switch then moves chain3 into position 1 and delets chain 1&2
	//In the process the labels need to be cached.

	utility::vector1<core::Size> orig_chain_order = get_chains(pose);
	utility::vector1<core::Size> chain_order_map;
	utility::vector1<core::Size> new_chain_order;
	chain_order_map.push_back(9999);
	core::Size new_chain1 = orig_chain_order[orig_chain_order.size()];
	new_chain_order.push_back(new_chain1);
	for ( core::Size ii=1; ii<orig_chain_order.size(); ++ii ) {//last elment is new pose
		if ( orig_chain_order[ii]!=chain1_id && orig_chain_order[ii]!=chain2_id ) {
			new_chain_order.push_back(orig_chain_order[ii]);
			chain_order_map.push_back(orig_chain_order[ii]);
		}
	}
	switch_chain_order(pose,new_chain_order);

	//put the labels back onto pose.
	//The first chain is a combination of chain1 and chain2 id from above. Put that on first
	vector1< vector1 <std::string> > labels_chain1 = orig_labels_by_chainid[chain1_id];
	vector1< vector1 <std::string> > labels_chain2 = orig_labels_by_chainid[chain2_id];
	//chain 1 begins at position 1
	for ( core::Size ii=1; ii<=labels_chain1.size(); ++ii ) {
		vector1<std::string> pos_labels = labels_chain1[ii];
		for ( core::Size jj=1; jj<=pos_labels.size(); ++jj ) {
			pose.pdb_info()->add_reslabel(ii, pos_labels[jj]);
		}
	}
	//chain 2 begins at position (length of chain1)
	for ( core::Size ii=1; ii<=labels_chain2.size(); ++ii ) {
		vector1<std::string> pos_labels = labels_chain2[ii];
		for ( core::Size jj=1; jj<=pos_labels.size(); ++jj ) {
			pose.pdb_info()->add_reslabel(ii+labels_chain1.size(), pos_labels[jj]);
		}
	}
	//from here on out the resnums for chain_id can be used
	chains = get_chains(pose);
	for ( core::Size ii=2; ii<=chains.size(); ++ii ) {
		vector1<core::Size> res_in_chain = get_resnums_for_chain_id(pose,chains[ii]);
		vector1<vector1 <std::string> > labels_in_chain = orig_labels_by_chainid[chain_order_map[ii]];
		for ( core::Size jj=1; jj<=res_in_chain.size(); ++jj ) {
			vector1<std::string> pos_labels = labels_in_chain[jj];
			for ( core::Size kk=1; kk<=pos_labels.size(); ++kk ) {
				pose.pdb_info()->add_reslabel(res_in_chain[jj], pos_labels[kk]);
			}
		}
	}
}

void NearNativeLoopCloser::switch_chain_order(Pose & pose, utility::vector1<core::Size> new_chain_order){
	//I've had to write this because the switchChainMover loses information from the pose
	core::pose::Pose new_pose=pose;
	new_pose.conformation().clear();
	new_pose.constraint_set(nullptr);//also clears energies
	new_pose.pdb_info( utility::pointer::make_shared< core::pose::PDBInfo >( new_pose, true ) ); //reinitialize the PDBInfo
	for ( Size ii=1; ii<=new_chain_order.size(); ++ii ) {
		core::pose::PoseOP chain = pose.split_by_chain(new_chain_order[ii]);
		append_pose_to_pose(new_pose,*chain,true);
		renumber_pdbinfo_based_on_conf_chains(new_pose,true,false,false,false);
	}
	pose=new_pose;
}


void NearNativeLoopCloser::extendRegion(bool towardCTerm, core::Size resStart, char neighborResType, core::Size numberAddRes,core::pose::PoseOP & poseOP){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace chemical;
	Real tmpPhi =  poseOP->phi(resStart);
	Real tmpPsi =  poseOP->psi(resStart);
	Real tmpOmega = poseOP->omega(resStart);
	if ( idealExtension_ ) {
		if ( neighborResType == 'E' ) {
			tmpPhi = -135;
			tmpPsi = 135;
			tmpOmega = 180;
		}
		if ( neighborResType == 'H' ) {
			tmpPhi = -57.8;
			tmpPsi = -47.0;
			tmpOmega = 180.0;
		}
	}
	core::conformation::ResidueOP new_rsd( nullptr );
	string build_aa_type_one_letter =option[OptionKeys::remodel::generic_aa];
	string build_aa_type = name_from_aa(aa_from_oneletter_code(build_aa_type_one_letter[0]));
	debug_assert( poseOP != nullptr );
	core::chemical::ResidueTypeSetCOP rs( poseOP->residue_type_set_for_pose() );
	kinematics::FoldTree ft;
	if ( towardCTerm == true ) {
		ft = poseOP->fold_tree();
		ft.clear();
		ft.add_edge(1,resStart,core::kinematics::Edge::PEPTIDE);
		ft.add_edge(1,poseOP->size(),1); //edges need to be added in order so the fold tree is completely connected.
		ft.add_edge(poseOP->size(),resStart+1,core::kinematics::Edge::PEPTIDE);
		poseOP->fold_tree(ft);
		for ( core::Size ii=0; ii<numberAddRes; ++ii ) {
			new_rsd = core::conformation::ResidueFactory::create_residue( rs->name_map(build_aa_type) );
			poseOP->conformation().safely_append_polymer_residue_after_seqpos( *new_rsd,resStart+ii, true);
		}
	} else {
		ft = poseOP->fold_tree();
		ft.clear();
		ft.add_edge(1,resStart-1,core::kinematics::Edge::PEPTIDE);
		ft.add_edge(1,poseOP->size(),1); //edges need to be added in order so the fold tree is completely connected.
		ft.add_edge(poseOP->size(),resStart,core::kinematics::Edge::PEPTIDE);
		poseOP->fold_tree(ft);
		for ( core::Size ii=0; ii<numberAddRes; ++ii ) {
			new_rsd = core::conformation::ResidueFactory::create_residue( rs->name_map(build_aa_type) );
			poseOP->conformation().safely_prepend_polymer_residue_before_seqpos( *new_rsd,resStart, true);
		}
	}
	if ( towardCTerm ) {
		for ( core::Size ii=0; ii<numberAddRes; ++ii ) {
			poseOP->set_phi(resStart+ii, tmpPhi );
			poseOP->set_psi(resStart+ii, tmpPsi );
			poseOP->set_omega(resStart+ii, tmpOmega );
		}
		poseOP->set_phi(resStart+numberAddRes,tmpPhi);
	} else {
		poseOP->set_psi(resStart, tmpPsi );
		poseOP->set_omega(resStart, tmpOmega );
		// when extending toward n-term the first psi/omega are undefined
		for ( core::Size ii=1; ii<=numberAddRes; ++ii ) {
			poseOP->set_phi(resStart+ii, tmpPhi );
			poseOP->set_psi(resStart+ii, tmpPsi );
			poseOP->set_omega(resStart+ii, tmpOmega );
		}
	}

	ft = poseOP->fold_tree();
	ft.clear();
	ft.add_edge(1,poseOP->size(),core::kinematics::Edge::PEPTIDE);
	poseOP->fold_tree(ft);
	renumber_pdbinfo_based_on_conf_chains(*poseOP,true,false,false,false);
}

core::pose::PoseOP NearNativeLoopCloser::create_maximum_length_pose(char resTypeBeforeLoop, char resTypeAfterLoop, core::pose::Pose const pose){
	core::pose::PoseOP full_length_poseOP = pose.clone();
	core::Size orig_resAfterLoop_ = resAfterLoop_;
	if ( resBeforeLoop_+1 != resAfterLoop_ ) {
		//trim region:
		kinematics::FoldTree ft;
		ft = full_length_poseOP->fold_tree();
		ft.add_edge(1,resAfterLoop_-1,core::kinematics::Edge::PEPTIDE);
		ft.add_edge(1,full_length_poseOP->size(),1); //edges need to be added in order so the fold tree is completely connected.
		ft.add_edge(full_length_poseOP->size(),resAfterLoop_,core::kinematics::Edge::PEPTIDE);
		full_length_poseOP->fold_tree(ft);
		full_length_poseOP->conformation().delete_residue_range_slow(resBeforeLoop_+1,resAfterLoop_-1);
		//repair fold tree
		ft = full_length_poseOP->fold_tree();
		ft.clear();
		ft.add_edge(1,full_length_poseOP->total_residue(),core::kinematics::Edge::PEPTIDE);
		full_length_poseOP->fold_tree(ft);
		renumber_pdbinfo_based_on_conf_chains(*full_length_poseOP,true,false,false,false);
		resAfterLoop_ = resBeforeLoop_+1;
	}
	int resAdjustmentStopHigh_or_0 = resAdjustmentStopHigh_;//This variable is being used to maintain the labels
	if ( resAdjustmentStopHigh_>0 ) {
		extendRegion(false,resAfterLoop_,resTypeAfterLoop,resAdjustmentStopHigh_,full_length_poseOP);
	} else {
		resAdjustmentStopHigh_or_0=0;
	}
	int resAdjustmentStartHigh_or_0 = resAdjustmentStartHigh_; //This variable is being used to maintain the labels
	if ( resAdjustmentStartHigh_>0 ) {
		extendRegion(true,resBeforeLoop_,resTypeBeforeLoop,resAdjustmentStartHigh_,full_length_poseOP);
	} else {
		resAdjustmentStartHigh_or_0=0;
	}
	for ( core::Size ii=1; ii<=resBeforeLoop_; ++ii ) {
		vector1 < std::string > tmp_labels = pose.pdb_info()->get_reslabels(ii);
		for ( core::Size jj=1; jj <= tmp_labels.size(); ++jj ) {
			full_length_poseOP->pdb_info()->add_reslabel(ii, tmp_labels[jj]);
		}
	}
	core::Size loopLabelAdjustment = orig_resAfterLoop_-resBeforeLoop_-1;
	for ( core::Size ii=orig_resAfterLoop_; ii<=pose.size(); ++ii ) {
		vector1 < std::string > tmp_labels = pose.pdb_info()->get_reslabels(ii);
		Size ii_in_full_length_pose = ( ii - loopLabelAdjustment ) + resAdjustmentStopHigh_or_0 + resAdjustmentStartHigh_or_0;
		for ( core::Size jj=1; jj <= tmp_labels.size(); ++jj ) {
			full_length_poseOP->pdb_info()->add_reslabel(ii_in_full_length_pose, tmp_labels[jj]);
		}
	}
	return(full_length_poseOP);
}

vector1<PossibleLoopOP> NearNativeLoopCloser::create_potential_loops(core::pose::Pose const pose){
	//create vector1 of possible loops
	vector1<PossibleLoopOP> possibleLoops;
	core::scoring::dssp::Dssp pose_dssp( pose );
	pose_dssp.dssp_reduced();
	std::string tmp_dssp = pose_dssp.get_dssp_secstruct();
	core::Size tmpResBeforeLoop = resBeforeLoop_;
	char resTypeBeforeLoop = tmp_dssp[tmpResBeforeLoop];
	//if there is a cut sometimes DSSP assigns the wrong SS
	while ( resTypeBeforeLoop == 'L' ) {
		tmpResBeforeLoop--;
		resTypeBeforeLoop = tmp_dssp[tmpResBeforeLoop];
	}
	core::Size tmpResAfterLoop = resAfterLoop_;
	char resTypeAfterLoop = tmp_dssp[tmpResAfterLoop];
	//if there is a cut sometimes DSSP assigns the wrong SS
	while ( resTypeAfterLoop == 'L' ) {
		tmpResAfterLoop--;
		resTypeAfterLoop = tmp_dssp[tmpResBeforeLoop];
	}
	if ( resTypeBeforeLoop == 'E' && resTypeAfterLoop == 'E' ) {
		//extend/trim equally
		int low_start = resAdjustmentStartLow_sheet_;
		int high_start = resAdjustmentStartHigh_sheet_;
		if ( low_start<resAdjustmentStopLow_sheet_ ) {
			low_start = resAdjustmentStopLow_sheet_;
		}
		if ( high_start>resAdjustmentStopHigh_sheet_ ) {
			high_start = resAdjustmentStopHigh_sheet_;
		}
		resAdjustmentStartHigh_=high_start; //So they are extended an equal length. This value applies to both sheet and helices. Kinda confusing I know.
		resAdjustmentStopHigh_=high_start;
		core::pose::PoseOP max_length_poseOP = create_maximum_length_pose(resTypeBeforeLoop,resTypeAfterLoop,pose);
		core::pose::PoseOP orig_atom_type_max_length_poseOP = max_length_poseOP->clone();
		if ( max_length_poseOP->is_fullatom() ) {
			core::util::switch_to_residue_type_set(*max_length_poseOP, core::chemical::CENTROID);
		}
		for ( int ii=low_start; ii<=high_start; ++ii ) {
			for ( core::Size kk=loopLengthRangeLow_; kk<=loopLengthRangeHigh_; ++kk ) {
				if ( (ii+resBeforeLoop_>=3)&&(ii+resAfterLoop_<=max_length_poseOP->total_residue()-3) ) { //ensures at least a 3 residue SS element next to loop
					PossibleLoopOP tmpLoopOP=utility::pointer::make_shared< PossibleLoop >(ii,ii,kk,resBeforeLoop_,resAfterLoop_,resTypeBeforeLoop,resTypeAfterLoop,resAdjustmentStartHigh_,resAdjustmentStopHigh_,max_length_poseOP,orig_atom_type_max_length_poseOP,fragment_store_path_,fragment_store_format_,fragment_store_compression_,numb_stubs_to_consider_);
					possibleLoops.push_back(tmpLoopOP);
				}
			}
		}
	} else {
		if ( resTypeBeforeLoop=='E' ) {
			TR.Warning << "Not extending 1 side of sheet and only eating in 1 residue" <<std::endl;
			resAdjustmentStartHigh_=0;
			resAdjustmentStartLow_=resAdjustmentStopLow_sheet_;
		}
		if ( resTypeAfterLoop=='E' ) {
			TR.Warning << "Not extending 1 side of sheet and only eating in 1 residue" << std::endl;
			resAdjustmentStopHigh_=0;
			resAdjustmentStopLow_=resAdjustmentStopHigh_sheet_;
		}
		core::pose::PoseOP max_length_poseOP = create_maximum_length_pose(resTypeBeforeLoop,resTypeAfterLoop,pose);
		core::pose::PoseOP orig_atom_type_max_length_poseOP = max_length_poseOP->clone();
		if ( max_length_poseOP->is_fullatom() ) {
			core::util::switch_to_residue_type_set(*max_length_poseOP, core::chemical::CENTROID);
		}
		//added to deal with an issue where the loop range is -1 to -4. In this case the helix is not extended so this should be 0.
		if ( resAdjustmentStartHigh_<0 ) {
			resAdjustmentStartHigh_=0;
		}
		if ( resAdjustmentStopHigh_<0 ) {
			resAdjustmentStopHigh_=0;
		}
		for ( int ii=resAdjustmentStartLow_; ii<=resAdjustmentStartHigh_; ++ii ) {
			for ( int jj=resAdjustmentStopLow_; jj<=resAdjustmentStopHigh_; ++jj ) {
				for ( core::Size kk=loopLengthRangeLow_; kk<=loopLengthRangeHigh_; ++kk ) {
					if ( (ii+resBeforeLoop_>=3)&&(jj+resAfterLoop_<=max_length_poseOP->total_residue()-3) ) {
						PossibleLoopOP tmpLoopOP=utility::pointer::make_shared< PossibleLoop >(ii,jj,kk,resBeforeLoop_,resAfterLoop_,resTypeBeforeLoop,resTypeAfterLoop,resAdjustmentStartHigh_,resAdjustmentStopHigh_,max_length_poseOP,orig_atom_type_max_length_poseOP,fragment_store_path_,fragment_store_format_,fragment_store_compression_,numb_stubs_to_consider_);
						possibleLoops.push_back(tmpLoopOP);
					}
				}
			}
		}
	}
	return(possibleLoops);
}

void
NearNativeLoopCloser::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
){
	std::string loopLengthRange( tag->getOption< std::string >( "loopLengthRange", "1,5") );
	rmsThreshold_ = tag->getOption< core::Real >( "RMSthreshold", 0.4 );
	if ( rmsThreshold_>0.6 ) {
		TR << "********************** rmsThresholds > 0.6 sometimes produces unclosed loops ***********************";
	}
	std::string resAdjustmentRange1( tag->getOption< std::string >( "resAdjustmentRangeSide1", "-3,3") );
	std::string resAdjustmentRange2( tag->getOption< std::string >( "resAdjustmentRangeSide2","-3,3") );
	chainBeforeLoop_ = tag->getOption<char>("chain",'A');
	chainAfterLoop_ = tag->getOption<char>("chain",'A');
	chainBeforeLoop_ = tag->getOption<char>("chainBeforeLoop",'A');
	chainAfterLoop_ = tag->getOption<char>("chainAfterLoop",'A');
	idealExtension_ = tag->getOption<bool>("idealExtension",true);
	label_loop_ = tag->getOption<std::string>("label_loop","");
	allowed_loop_abegos_ = tag->getOption< std::string >( "allowed_loop_abegos","");
	if ( !idealExtension_ ) {
		utility_exit_with_message("the ideal extension flag works but it seems to mess up the pose more than help so delete this line & recompile if you really want to use it");
	}
	max_vdw_change_ = tag->getOption<core::Real>("max_vdw_change",8.0);
	ideal_ = tag->getOption<bool>("ideal",false);
	if ( chainBeforeLoop_==chainAfterLoop_ ) {
		resBeforeLoop_ = tag->getOption<core::Size>("resBeforeLoop");
		resAfterLoop_ = tag->getOption<core::Size>("resAfterLoop");
	}
	output_closed_ = tag->getOption<bool>("close",true);
	closure_type_ = tag->getOption<std::string>("closure_type","lookback");
	if ( output_closed_ ) {
		fragment_store_path_=tag->getOption<std::string>("fragment_store","");
		fragment_store_format_=tag->getOption<std::string>("fragment_store_format","hashed");
		fragment_store_compression_=tag->getOption<std::string>("fragment_store_compression","all");
		numb_stubs_to_consider_ = tag->getOption<core::Size>("numb_stubs_to_consider",1);
		SSHashedFragmentStoreOP_ = protocols::indexed_structure_store::FragmentStoreManager::get_instance()->SSHashedFragmentStore(fragment_store_path_,fragment_store_format_,fragment_store_compression_);
		SSHashedFragmentStoreOP_->set_threshold_distance(rmsThreshold_);
		TR << "database loaded!!" << std::endl;
	}
	output_all_= tag->getOption<bool>("output_all",false);
	max_number_of_results_ = tag->getOption<core::Size>("max_number_of_results",999999);
	numb_outputed_=0;
	top_outputed_ = false;
	utility::vector1< std::string > resAdjustmentRange1_split( utility::string_split( resAdjustmentRange1 , ',' ) );
	utility::vector1< std::string > resAdjustmentRange2_split( utility::string_split( resAdjustmentRange2 , ',' ) );
	utility::vector1< std::string > loopLengthRange_split( utility::string_split( loopLengthRange , ',' ) );
	if ( resAdjustmentRange1_split.size()==2 ) {
		resAdjustmentStartLow_ = atoi(resAdjustmentRange1_split[1].c_str());
		resAdjustmentStartHigh_ = atoi(resAdjustmentRange1_split[2].c_str());
	}
	if ( resAdjustmentRange1_split.size()==1 ) {
		resAdjustmentStartLow_= atoi(resAdjustmentRange1_split[1].c_str());
		resAdjustmentStartHigh_= atoi(resAdjustmentRange1_split[1].c_str());
	}
	if ( resAdjustmentRange2_split.size()==2 ) {
		resAdjustmentStopLow_ = atoi(resAdjustmentRange2_split[1].c_str());
		resAdjustmentStopHigh_ = atoi(resAdjustmentRange2_split[2].c_str());
	}
	if ( resAdjustmentRange2_split.size()==1 ) {
		resAdjustmentStopLow_ = atoi(resAdjustmentRange2_split[1].c_str());
		resAdjustmentStartHigh_ = atoi(resAdjustmentRange2_split[1].c_str());
	}
	if ( loopLengthRange_split.size()==2 ) {
		loopLengthRangeLow_ = atoi(loopLengthRange_split[1].c_str());
		loopLengthRangeHigh_ = atoi(loopLengthRange_split[2].c_str());
	}
	if ( loopLengthRange_split.size()==1 ) {
		loopLengthRangeLow_ = atoi(loopLengthRange_split[1].c_str());
		loopLengthRangeHigh_ = atoi(loopLengthRange_split[1].c_str());
	}
	resAdjustmentStartLow_sheet_ = resAdjustmentStartLow_;
	resAdjustmentStartHigh_sheet_ = resAdjustmentStartHigh_;
	resAdjustmentStopLow_sheet_ = resAdjustmentStopLow_;
	resAdjustmentStopHigh_sheet_ = resAdjustmentStopHigh_;
	//TR << resAdjustmentStartLow_ <<"," << resAdjustmentStartHigh_ << ",:," << resAdjustmentStopLow_ << "," << resAdjustmentStopHigh_ << ",:," << loopLengthRangeLow_ <<"," << loopLengthRangeHigh_ << std::endl;
}

core::pose::PoseOP NearNativeLoopCloser::get_additional_output(){
	Real tmp_rmsd;
	return(get_additional_output_with_rmsd(tmp_rmsd));
}

core::pose::PoseOP NearNativeLoopCloser::get_additional_output_with_rmsd(Real & return_rmsd){
	std::sort(possibleLoops_.begin(), possibleLoops_.end(), FinalRMSDComparator());
	if ( !output_all_&&top_outputed_ ) {
		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
		return nullptr;
	}
	if ( numb_outputed_>=max_number_of_results_ ) {
		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
		return nullptr;
	}
	numb_outputed_+=1;
	for ( core::Size ii=1; ii<=possibleLoops_.size(); ++ii ) {
		if ( !possibleLoops_[ii]->outputed()&&(possibleLoops_[ii]->get_final_RMSD()<rmsThreshold_) ) {
			TR << "Loop outputed with " << possibleLoops_[ii]->get_final_RMSD() << " rmsd" << std::endl;
			possibleLoops_[ii]->setup_finalPose_copy_rotamers();
			possibleLoops_[ii]->setup_finalPose_copy_labels();
			possibleLoops_[ii]->label_loop(label_loop_);
			if ( allowed_loop_abegos_!="" ) {
				utility::vector1< std::string > const abego_v = core::sequence::get_abego( *possibleLoops_[ii]->get_finalPoseOP() );
				core::sequence::ABEGOManager am;
				std::string abego = am.get_abego_string(abego_v);
				TR << "ABEGO of pose" << abego << std::endl;
			}
			possibleLoops_[ii]->outputed(true);
			top_outputed_=true;
			set_last_move_status(protocols::moves::MS_SUCCESS);
			possibleLoops_[ii]->get_finalPoseOP()->data().set(
				core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG,utility::pointer::make_shared< basic::datacache::CacheableString >( pose_name_ ) );
			return_rmsd = possibleLoops_[ii]->get_final_RMSD();
			return(possibleLoops_[ii]->get_finalPoseOP());
		}
		if ( !output_closed_ && (!possibleLoops_[ii]->outputed()) ) {
			possibleLoops_[ii]->outputed(true);
			possibleLoops_[ii]->setup_finalPose_copy_rotamers();
			possibleLoops_[ii]->setup_finalPose_copy_labels();
			possibleLoops_[ii]->label_loop(label_loop_);
			set_last_move_status(protocols::moves::MS_SUCCESS);
			return_rmsd = possibleLoops_[ii]->get_final_RMSD();
			return(possibleLoops_[ii]->get_finalPoseOP());
		}
	}
	// Real low_rmsd=999;
	// for ( core::Size ii=1; ii<possibleLoops_.size(); ++ii ) {
	//  if ( !possibleLoops_[ii]->outputed() ) {
	//   if ( possibleLoops_[ii]->get_final_RMSD()< low_rmsd ) {
	//    low_rmsd = possibleLoops_[ii]->get_final_RMSD();
	//   }
	//  }
	// }
	set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
	return_rmsd = 9999;
	return nullptr;
}

std::string NearNativeLoopCloser::get_name() const {
	return mover_name();
}

std::string NearNativeLoopCloser::mover_name() {
	return "NearNativeLoopCloser";
}

void NearNativeLoopCloser::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"loopLengthRange", xs_string,
		"The range of loops checked. Max is 5", "1,5");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"RMSthreshold", xsct_real,
		"The worst 9-residue RMSD of a loop", "0.4");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"resAdjustmentRangeSide1", xs_string,
		"Trim/extensions of secondary structure elements", "-3,3");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"resAdjustmentRangeSide2", xs_string,
		"Trim/extension of secondary structure elements", "-3,3");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"chain", xs_string,
		"XSD_XRW: TO DO", "A");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"chainBeforeLoop", xs_string,
		"chain id before loop", "A");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"chainAfterLoop", xs_string,
		"chain id after loop ", "A");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"idealExtension", xsct_rosetta_bool,
		"copies use ideal helical coordinates,true, or the last phi/psi before the loop ,false,", "true");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"max_vdw_change", xsct_real,
		"XSD_XRW: TO DO", "8.0");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"ideal", xsct_rosetta_bool,
		"uses kic to make the loop ideal. Used primarily in repeat proteins where the bond lengths needed to be kept ideal", "false");
	attlist + XMLSchemaAttribute(
		"resBeforeLoop", xsct_non_negative_integer,
		"residue before loop");
	attlist + XMLSchemaAttribute(
		"resAfterLoop", xsct_non_negative_integer,
		"residue after loop");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"output_all", xsct_rosetta_bool,
		"output all structures", "false");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"max_number_of_results", xsct_non_negative_integer,
		"for use with output_all. This limits the maximum number of structures", "999999");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"closure_type", xs_string,
		"type of closure, either kic or lookback", "lookback");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"label_loop", xs_string,
		"label loop in pdb_info object", "");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"fragment_store", xs_string,
		"path to fragment store. Note:All fragment stores use the same database", "");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"fragment_store_format", xs_string,
		"Options:hashed,unhashed new format is unhashed", "hashed");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"fragment_store_compression", xs_string,
		"Options:helix_shortLoop,sheet_shortLoop,all", "all");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"numb_stubs_to_consider", xsct_non_negative_integer,
		"turns out the best loop stub doesn't always give the loop with the lowest RMSD. This allows you to consider additional stubs","1");
	attlist + XMLSchemaAttribute::attribute_w_default( "allowed_loop_abegos", xs_string, "comma seperated string of allowed abegos, default=empty all abegos", "" );

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"limit loops to only those that match the given abegos",
		attlist );
}

std::string NearNativeLoopCloserCreator::keyname() const {
	return NearNativeLoopCloser::mover_name();
}

protocols::moves::MoverOP
NearNativeLoopCloserCreator::create_mover() const {
	return utility::pointer::make_shared< NearNativeLoopCloser >();
}

void NearNativeLoopCloserCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	NearNativeLoopCloser::provide_xml_schema( xsd );
}


}//pose_length_moves
}//protocols
