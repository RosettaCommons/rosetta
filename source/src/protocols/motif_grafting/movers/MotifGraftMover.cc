// Project Headers
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/protein_interface_design/movers/MotifGraftMover.cc
/// @brief  Declaration of the MoverCreator class for the MotifGraftCreator
/// @author Daniel-Adriano Silva (dadriano@uw.edu) and Alex Ford (fordas@uw.edu)
/// ToDo:
///  1. TEST in hostile environments. (PARTIALLY DONE)
///  2. Replace TR.Info for TR.Debug in most cases.


//Include C++ classes
#include <sstream>
#include <string>
#include <queue>
//#include <iterator>

//Include Rosetta utilities
#include <utility/vector1.hh>

//Include Rosetta numeric
#include <numeric/xyz.functions.hh>
#include <numeric/MathMatrix_operations.hh>

//Include Rosetta Tracer
#include <basic/Tracer.hh>

//Include Rosetta boost c_ap
#include <boost/algorithm/string.hpp>   
#include <boost/lexical_cast.hpp>   

//Include Rosetta XML tag reader
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

//Include Rosetta Core Stuff
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/VariantType.hh>

//Include Rosetta Scoring functions
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

//Include Rosetta Mover protocols
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/DataMap.fwd.hh>

//Include Rosetta protocols 
#include <protocols/toolbox/pose_manipulation.hh>
#include <protocols/toolbox/superimpose.hh>

//Include Rosetta Job distributor engine
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

//Include Rosetta scripts
#include <protocols/rosetta_scripts/util.hh>

//Include ObjexxFCL
#include <ObjexxFCL/FArray.all.hh>

//Include my headers
#include <protocols/motif_grafting/movers/MotifGraftMover.hh>
#include <protocols/motif_grafting/movers/MotifGraftCreator.hh>


/**@brief This is a protocol... **/
namespace protocols
{
	/**@brief ...that lives inside of motif_grafting... **/
	namespace motif_grafting
	{
		/**@brief ...and belongs to the movers class... **/
		namespace movers
		{
			/**@brief Global class rosetta protocols tracer**/
			static basic::Tracer TR( "protocols.motif_grafting.movers.MotifGraftMover" ); 
			
			
			MotifGraftMover::MotifGraftMover()
			{
			}
			
			//Where does this goes?
			/** :
			gp_p_contextStructure_(NULL),
			gp_p_motif_(NULL),
			gp_r_RMSD_tolerance_(NULL),
			gp_r_clash_score_cutoff_(NULL),
			gp_vp_combinatory_fragment_size_delta_(NULL),
			gp_vp_max_fragment_replacement_size_delta_(NULL),
			gp_s_clash_test_residue_(NULL),
			gp_b_full_motif_bb_alignment_(NULL),
			gp_p_target_pose_(NULL),
			motif_match_results_(NULL)
			 * */
			/**@brief MotifGraftMover parameters and options initializer**/
			void MotifGraftMover::init_parameters(
			std::string const & s_contextStructure,
			std::string const & s_motif,
			core::Real  const & r_RMSD_tolerance,
			core::Real  const & r_clash_score_cutoff,
			std::string const & s_combinatory_fragment_size_delta,
			std::string const & s_max_fragment_replacement_size_delta,
			std::string const & s_clash_test_residue,
			bool        const & b_full_motif_bb_alignment)
			{
				//Parse the arguments in the global space variables
				MotifGraftMover::parse_my_string_arguments_and_cast_to_globalPrivateSpaceVariables(
					s_contextStructure,
					s_motif,
					r_RMSD_tolerance,
					r_clash_score_cutoff,
					s_combinatory_fragment_size_delta,
					s_max_fragment_replacement_size_delta,
					s_clash_test_residue,
					b_full_motif_bb_alignment);
			}
			
			/**@brief MotifGraftMover Destructor**/
			MotifGraftMover::~MotifGraftMover() 
			{
			}
			
			/**@brief Apply mover function**/
			void MotifGraftMover::apply(core::pose::Pose & pose)
			{
				TR.Info << "Executing Motif Graft Mover " << std::endl;
				if (!motif_match_results_.empty() )
				{
					// Error
					TR.Warning << "Warning. Received new input pose with matches remaining." << std::endl;
				}
				
				motif_match_results_ = generate_scaffold_matches(pose, gp_p_motif_, gp_p_contextStructure_);
			
				if (motif_match_results_.size() == 0)
				{
					TR.Info << "Generated no results." << std::endl;
					// Matching failed, set fail_do_not_retry return code.
					set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
					return;
				}
				else if (motif_match_results_.size() > 1)
				{
					TR.Info << "Generated " << motif_match_results_.size() << " results." << std::endl;
					//Create copy of input pose to support result interation
					gp_p_target_pose_ = new Pose(pose);
				}
				
				MotifMatch match = motif_match_results_.top();
				motif_match_results_.pop();
				
				generate_match_pose(pose, match);
			}
			
			/**@brief Iterate over the results to get additional matches in the queue**/
			core::pose::PoseOP MotifGraftMover::get_additional_output()
			{
				core::pose::PoseOP work_pose;
			
				if (motif_match_results_.size() == 0) 
				{
					TR.Debug << "No additional output." << std::endl;
					return NULL;
				}
				else if (motif_match_results_.size() == 1)
				{
					// Last iteration, do not create a new pose reference
					work_pose = gp_p_target_pose_;
				}
				else
				{
					// May produce additional output poses, create a copy of the target pose
					work_pose = new Pose(*gp_p_target_pose_);
				}
			
				TR.Debug << "Returning additional output. " << motif_match_results_.size() << " remaining." << std::endl;
				
				MotifMatch match = motif_match_results_.top();
				motif_match_results_.pop();
				
				generate_match_pose(*work_pose, match);
				
				return work_pose;
			}
			
			/**@brief Identify all potential matches for the given target scaffold (this is where the motif grafting code is called)**/
			std::priority_queue<MotifMatch> MotifGraftMover::generate_scaffold_matches(
				core::pose::Pose const & target_scaffold, 
				core::pose::PoseOP const & target_motif_,
				core::pose::PoseOP const & target_contextStructure_)
			{
				TR.Info << "Target scaffold size: " << target_scaffold.total_residue() << std::endl;
				TR.Info << "Motif size: " << target_motif_->total_residue() << std::endl;
				TR.Info << "contextStructure size: " << target_contextStructure_->total_residue() << std::endl;
				TR.Info << "Target scaffold numChains: " << target_scaffold.conformation().num_chains() << std::endl;
				TR.Info << "Motif numChains: " << target_motif_->conformation().num_chains() << std::endl;
				TR.Info << "contextStructure numChains: " << target_contextStructure_->conformation().num_chains() << std::endl;
				
				//Create the return Priority Queue
				std::priority_queue<MotifMatch> pq_epigraft;
				
				//Just now it cannot work with scaffolds larger than one monomer
				//ToDo: Make it work for multi-chain targets?
				if (target_scaffold.conformation().num_chains() > 1 ){
					set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
					return pq_epigraft;
				}

				//Discriminate_fragments_matching the target_motif_ and compatible with the target_contextStructure_
				//Returns the result by reference to the priority queue (pq_epigraft)
				get_matching_fragments(target_scaffold, target_motif_, target_contextStructure_, gp_r_RMSD_tolerance_, gp_r_clash_score_cutoff_, gp_s_clash_test_residue_, gp_vp_max_fragment_replacement_size_delta_, gp_vp_combinatory_fragment_size_delta_, gp_b_full_motif_bb_alignment_, pq_epigraft);
				
				return pq_epigraft;
			}
			
			/**@brief Return a priority queue with the successful epigrafts **/
			void MotifGraftMover::get_matching_fragments(
				core::pose::Pose const & p_scaffold, 
				core::pose::PoseOP const & p_motif_,
				core::pose::PoseOP const & p_contextStructure_,
				core::Real const & RMSD_tol,
				core::Real const & clash_cutoff,
				std::string const & clash_test_residue,
				utility::vector1 < std::pair< long int, long int > > const & max_fragment_replacement_size_delta,
				utility::vector1 < std::pair< core::Size, core::Size > > const & combinatory_fragment_size_delta,
				bool const & b_full_motif_bb_alignment,
				std::priority_queue<MotifMatch> & pq)
			{
				//Will store all the posible motif-fragments combinations (when variated by combinatory_fragment_size_delta)
				utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > vv_motif_fragments_all_posible_permutations;
				//Generate all the combination of different legths of the motif fragment as requested in combinatory_fragment_size_delta, and iterate the complete epigraft function using it
				generate_combinations_of_motif_fragments_by_delta_variation(p_motif_, combinatory_fragment_size_delta, vv_motif_fragments_all_posible_permutations);
				//Test combinations of fragments (outer loop)
				for ( core::Size numPermutation = 1; numPermutation <= vv_motif_fragments_all_posible_permutations.size(); ++numPermutation )
				{
					TR.Info << "Fragments combination to be tested: ";
					for ( core::Size i = 1; i <= vv_motif_fragments_all_posible_permutations[numPermutation].size(); ++i )
					{
						TR.Info << vv_motif_fragments_all_posible_permutations[numPermutation][i].first << " " << vv_motif_fragments_all_posible_permutations[numPermutation][i].second << " ";
					}
					TR.Info << std::endl;
					
					//Will store the indexes for the [Nter, Cter] of scaffold fragments matching each motif fragment.
					utility::vector1< utility::vector1 < std::pair< core::Size, core::Size > > >  vv_scaffold_fragments_indexes;
					//Will stores the indexes for the [Nter, Cter] of each motif fragment.
					utility::vector1< std::pair< core::Size, core::Size > >  v_motif_fragments_indexes;
					//Get the matching fragments by CA distances
					bool enough_fragments = MotifGraftMover::get_fragments_by_CA_distances(p_scaffold, p_motif_, vv_scaffold_fragments_indexes, v_motif_fragments_indexes, RMSD_tol, max_fragment_replacement_size_delta, vv_motif_fragments_all_posible_permutations[numPermutation]);
					//If we do not have enough fragments to satisfise the restrictions DIE with exception.
					if(!enough_fragments){
						//throw utility::excn::EXCN_Msg_Exception("For this scaffold I couldn't find suitable fragments that match the size of your motifs.");
						TR.Warning << "For this scaffold & fragment/size combination there are not fragments that match the size of your motifs."<< std::endl;
						continue;
					}
					
					//Will store the fragment combinations that pass the combination test
					utility::vector1< motif2scaffold_data > v_m2s_data;
					//Helper. Will store the temporary Buffer for the results of the recursive discontinuous fragments combinations search function
					utility::vector1< std::pair< core::Size, core::Size > > buff_combVec;
					//Automatically generate all the discontinuous fragments combinations that are within the tol restriction
					MotifGraftMover::fragments_permutation_test_by_CA_distances(p_scaffold, p_motif_, vv_scaffold_fragments_indexes, v_motif_fragments_indexes, RMSD_tol, 1, buff_combVec, v_m2s_data);
					TR.Debug << "Num of Fragments so far: " << v_m2s_data.size() << std::endl;
					//If we do not have enough fragments to satisfise the restrictions DIE with exception.
					if(v_m2s_data.size() == 0){
					//throw utility::excn::EXCN_Msg_Exception("For this scaffold I couldn't find any suitable fragment combination that match your motif's geometric restrictions.");
						TR.Warning << "For this scaffold & fragment combination there are not fragments combinations that match the geometric restrictions of your motifs."<< std::endl;
						continue;
					}
				
					//Test each of the fragments based on the RMSD of the extremes of motifs
					//Returns by reference only the fragments that PASS, populates vv_m2s_ndx[i] with the rotation matrix and translation vector
					MotifGraftMover::get_motif_scaffold_superposition_and_RMSD(p_scaffold, p_motif_, RMSD_tol, b_full_motif_bb_alignment, v_m2s_data);
					TR.Debug << "Num of Fragments so far: " << v_m2s_data.size() << std::endl;
					//If we do not have enough fragments to satisfise the restrictions DIE with exception.
					if(v_m2s_data.size() == 0){
						//throw utility::excn::EXCN_Msg_Exception("For this scaffold I couldn't find any suitable fragment combination that pass the RMSD test.");
						TR.Warning << "For this scaffold & fragment combination there are not fragments that pass the RMSD test."<< std::endl;
						continue;
					}
					
					//Determine wich mode of comparision will be used for the clash test (Native or mutated residues, i.e. original vs ... lets say minimal protein)
					//The comparision structure(s) will be allocated on p_scaffold_mono_aa and p_motif_mono_aa.
					//In my point of view the most powerful and reasoneable option is the GLY mutation mode. 
					//If you disagree and want to have a discussion about this drop me an email to:  dadriano at gmail dot com (Daniel Silva)
					core::pose::Pose p_scaffold_mono_aa;
					core::pose::Pose p_motif_mono_aa;
					if(clash_test_residue == "NATIVE"){
						TR.Info << "Using NATIVE mode for clash checking, expect a low number of results" << v_m2s_data.size() << std::endl;
						//Use the native scaffold and pose (using this option is a DUMB, indeed), added for backward/historical consistency
						p_scaffold_mono_aa = p_scaffold;
						p_motif_mono_aa = *p_motif_;
					}else if(clash_test_residue == "GLY" || clash_test_residue == "ALA" || clash_test_residue == "VAL"){
						TR.Info << "Using mutation mode for clash checking, with residue: " << clash_test_residue << std::endl;
						//get mono-/uni- aminoacid copies of the scaffold and motif for the clash test
						p_scaffold_mono_aa = MotifGraftMover::get_mono_aa_pose_copy( p_scaffold , clash_test_residue);
						p_motif_mono_aa = MotifGraftMover::get_mono_aa_pose_copy( *p_motif_ , clash_test_residue);
					}else{
						throw utility::excn::EXCN_Msg_Exception("The aminoacid selected for the clash test is not valid. Valid selections are: \"GLY, ALA, VAL, NATIVE\"");
					}
					
					//Test the remaining fragments based on the clash score
					MotifGraftMover::test_epigraft_and_contextStructure_clashes( p_scaffold_mono_aa, p_motif_mono_aa, *p_contextStructure_, clash_cutoff, v_m2s_data);
					TR.Debug << "Num of Fragments so far: " << v_m2s_data.size() << std::endl;
					if(v_m2s_data.size() == 0){
						//throw utility::excn::EXCN_Msg_Exception("For this scaffold I couldn't find any suitable fragment combination that pass the clash_score test.");
						TR.Warning << "For this scaffold & fragment combination there are not fragments that pass the clash_score test."<< std::endl;
						continue;
					}
					
					//Store the results in the priority queue
					for (core::Size i=1; i <= v_m2s_data.size(); ++i){
						MotifMatch tmpResult = MotifMatch(v_m2s_data[i]);
						pq.push(tmpResult);
					}
					TR.Info << "In the iteration " << numPermutation << ",there are a total of: " << pq.size() << " fragments in the priority queue"<< std::endl;
				}//END Test combinations of fragments (outer loop)
				//return
				if (pq.size() == 0 ){
					throw utility::excn::EXCN_Msg_Exception("For this scaffold there are not suitable scaffold grafts within your constrains");
				}else{
					TR.Info << "After all the iterations there are a total of: " << pq.size() << " fragments in the priority queue"<< std::endl;
				}
			}
			
			/**@brief Generate all the combination of different legths of the motif fragment as requested in combinatory_fragment_size_delta**/
			//ToDo: Check for error: span overlapping
			void MotifGraftMover::generate_combinations_of_motif_fragments_by_delta_variation(
				core::pose::PoseOP const & p_motif_,
				utility::vector1 < std::pair< long int, long int > > const & combinatory_fragment_size_delta,
				utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > & vv_resulting_permutations)
			{
				TR.Debug << "Generating combinations of fragments as requested in combinatory_fragment_size_delta" << std::endl;
				utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > vv_tmp_motif_fragments_indexes;
				//Stores lowbound and highbound of each fragment(s) combination
				std::pair< core::Size, core::Size > tmpMotifNdx;
				//Get the resd IDs for the fragments in the PDB
				//Also get the Ca<->Ca distance between the beggining and end of each fragment
				for ( core::Size i = 1; i <= p_motif_->conformation().num_chains(); ++i )
				{
					utility::vector1< std::pair< core::Size, core::Size > >  v_tmp_motif_fragments_indexes;
					//long int because we need an opperation that can be negative!!!
					for (long int j=p_motif_->conformation().chain_begin(i); j <= (long int)(p_motif_->conformation().chain_begin(i)+combinatory_fragment_size_delta[i].first); ++j )
					{
						//long int because we need an opperation that can be negative!!!
						for (long int k=p_motif_->conformation().chain_end(i); k >= (long int)(p_motif_->conformation().chain_end(i)-combinatory_fragment_size_delta[i].second); --k )
						{
							//The motif fragment cannot be samaler than two aminoacids
							//ToDo: Fix this Bug, core::Size - core::Size is always positive!
							if ( (k-j) < 1 ){
								continue;
							}
							//Store the beggining and the end in motif_extremes vector1
							tmpMotifNdx.first  = (core::Size)(j);
							tmpMotifNdx.second = (core::Size)(k);
							TR.Debug << "Applying deltas to the motif. Variant of motif fragment #" << i << ". Start at: " << tmpMotifNdx.first << ", end: " << tmpMotifNdx.second << std::endl;
							v_tmp_motif_fragments_indexes.push_back(tmpMotifNdx);
						}
					}
					vv_tmp_motif_fragments_indexes.push_back(v_tmp_motif_fragments_indexes);
				}
				utility::vector1< std::pair< core::Size, core::Size > > buff_combVec;
				permutate_n_vv_of_pairs(vv_tmp_motif_fragments_indexes, buff_combVec, 1, vv_resulting_permutations);
				
				/////Finally check that there are not combinations in which the fragments are reversed (ini>end) or smaller than 2 aa
				///utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > >::iterator fragment_permutation = vv_resulting_permutations.begin();
				///while( fragment_permutation != vv_resulting_permutations.end()) 
				///{
				///	bool is_good=true;
				///	for ( core::Size i = 1; i <= (*fragment_permutation).size(); ++i )
				///	{
				///		TR.Info << "F-test: [" << (*fragment_permutation)[i].first << ":" << (*fragment_permutation)[i].second << "]" << std::endl;
				///		if ( ( (*fragment_permutation)[i].second - (*fragment_permutation)[i].first ) < 1 ){
				///			TR.Debug << "Removing a bad combination (a fragment is reversed or Size is minor than 2): [" << (*fragment_permutation)[i].first << ":" << (*fragment_permutation)[i].second << "]" << std::endl;
				///			is_good=false;
				///			break;
				///		}
				///	}
				///	if(is_good){
				///		++fragment_permutation;
				///	}else{
				///		fragment_permutation = vv_resulting_permutations.erase(fragment_permutation);
				///	}
				///}
			}
			
			/** @brief As the name suggests in generates all the permutations of a vector of vectors of pairs (Alex: we should templatize this)**/
			//CAUTION: This loop is self recursive, use extreme caution when calling it
			void MotifGraftMover::permutate_n_vv_of_pairs(
				utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > const & vv_of_pairs,
				utility::vector1< std::pair< core::Size, core::Size > > & buff_combVec,
				core::Size start_index,
				utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > & vv_resulting_permutations)
			{
				//If we reach this condition it means that there is a sucessfull combination of elements
				if ( start_index >= ( vv_of_pairs.size()+1 ) )
				{
					//push the data in the return vector
					vv_resulting_permutations.push_back(buff_combVec);
					//finish with good results
					return;
				}
				for (core::Size i=1; i <= vv_of_pairs[start_index].size(); ++i){
					utility::vector1< std::pair< core::Size, core::Size > > buff_combVec2 = buff_combVec;
					buff_combVec2.push_back(vv_of_pairs[start_index][i]);
					MotifGraftMover::permutate_n_vv_of_pairs(vv_of_pairs, buff_combVec2, start_index+1, vv_resulting_permutations);
				}
				//finish without results (In theory We should never reach this point)
				return;
			}
			
			/**@brief Functions that takes the scaffold, motif, contextStructure and superposition transform data. Deletes from the supperposition data
			 ** those transformations that can't pass the clash score**/
			void MotifGraftMover::test_epigraft_and_contextStructure_clashes(
				core::pose::Pose const & p_scaffold, 
				core::pose::Pose const & p_motif_,
				core::pose::Pose const & p_contextStructure_,
				core::Real const & clash_cutoff,
				utility::vector1< motif2scaffold_data > & v_m2s_data)
			{
				//Create the scoring function
				//ToDo: Add flexibility so the user can control it from Rosetta
				core::scoring::ScoreFunctionOP scorefxn_  = new core::scoring::ScoreFunction();
				scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
				scorefxn_->set_weight( core::scoring::fa_atr, 1.0 );
				
				//Iterate the fragments
				utility::vector1< motif2scaffold_data >::iterator it_fragments = v_m2s_data.begin();
				while( it_fragments != v_m2s_data.end()) 
				{
					//Rotate/translate and stitch the motif and the scaffold
					core::pose::Pose p_frankenstein=MotifGraftMover::stich_motif_in_scaffold_by_indexes_rotation_and_translation(p_scaffold, p_motif_, (*it_fragments), true);
					//Join the contextStructure and epigraft(frankenstein) in a single pose (reuse the container frankenstein)
					//Alex: can I reuse the pose like this? Is/will the size dynamicaly re-allocated?
					p_frankenstein=MotifGraftMover::join_two_poses_by_jump(p_contextStructure_, p_frankenstein);
					//Calculate the clash score
					//ToDo: Implement a more modern way to execute the clash check, 
					//       like the one in "protocols/sic_dock/xyzStripeHashPose.hh" or "numeric/geometry/hashing/xyzStripeHash.hh"
					core::Real clashScore = MotifGraftMover::get_clash_score_from_pose(p_frankenstein, scorefxn_);
					//p_frankenstein.dump_pdb( "test_z"+utility::to_string(clashScore)+".pdb");
					//decide if remove or keep the object/epigraft
					if (clashScore <= clash_cutoff)
					{
						TR.Debug << "Epigraft passed the clash score [fa_rep+fa_atr] test with a value of: " << clashScore << std::endl;
						(*it_fragments).clash_score=clashScore;
						++it_fragments;
					}else{
						//remove from the vector those elements that don't pass the clash test
						TR.Debug << "Epigraft failed the clash score test [fa_rep+fa_atr] with a value of: " << clashScore << std::endl;
						it_fragments = v_m2s_data.erase(it_fragments);
					}
				}
			}
			
			/**@brief returns a pose with two input poses merged (with a jump in-between) and with the PDB info corrected*/
			//Interestingly this is the "slowest part of the code" so far. Why?
			core::pose::Pose MotifGraftMover::join_two_poses_by_jump(
				core::pose::Pose const & p_A, 
				core::pose::Pose const & p_B)
			{
				core::pose::Pose p_result=p_A;
				p_result.conformation().insert_conformation_by_jump(p_B.conformation(), 
											p_result.total_residue() + 1, 
											p_result.fold_tree().num_jump()+1,
											p_result.total_residue());
				//Re-create the PDBinfo (this is really creepy! Please write a proper function to append by jump)
				p_result.pdb_info( new core::pose::PDBInfo(p_result, true) );
				return p_result;
			}
			
			/**@brief score fa_rep and fa_atr in the pose*/
			core::Real MotifGraftMover::get_clash_score_from_pose(
				core::pose::Pose & p_input,
				core::scoring::ScoreFunctionOP const & scorefxn_)
			{
				(*scorefxn_)(p_input);
				core::Real const fa_rep = p_input.energies().total_energies()[ core::scoring::fa_rep ];
				core::Real const fa_atr = p_input.energies().total_energies()[ core::scoring::fa_atr ];
				return (fa_rep+fa_atr);
			}
			
			/**@brief Function that returns by reference a rotated copy of the pose */
			core::pose::Pose MotifGraftMover::get_rotated_and_translated_pose(
				core::pose::Pose const & p_input,
				numeric::xyzMatrix< core::Real > const & RotM,
				numeric::xyzVector< core::Real > const & TvecA,
				numeric::xyzVector< core::Real > const & TvecB)
			{
				//Apply rotation to ALL atoms (analog to Alex Ford code protocols/toolbox/PlaceFragments.cc; thought he denies the code)
				// x_i' <- = R*x_i + com1;
				core::pose::Pose p_result = p_input;
				for ( Size i = 1; i <= p_input.total_residue(); ++i ) {
					for ( Size j = 1; j <= p_input.residue_type(i).natoms(); ++j ) {
						core::id::AtomID id( j, i );
						p_result.set_xyz( id, ( RotM * ( p_input.xyz(id) + TvecB) ) - TvecA );
					}
				}
				return p_result;
			}
			
			/**@brief Helper function to stich (epigraft) two poses given a set of indices in pose A and B stored in a motif2scaffold_data structure**/
			core::pose::Pose MotifGraftMover::stich_motif_in_scaffold_by_indexes_rotation_and_translation(
				core::pose::Pose const & p_scaffold, 
				core::pose::Pose const & p_motif, 
				motif2scaffold_data & m2s_dat,
				bool const & skip_motif_extremes)
			{
				//Will store the stitched pose
				core::pose::Pose p_result;
				//Get a copy of the scaffold rotated to the correct possition
				core::pose::Pose p_scaffold_rotated = get_rotated_and_translated_pose(p_scaffold, m2s_dat.RotM, m2s_dat.TvecA, m2s_dat.TvecB);
				//make a copy of the motif
				core::pose::Pose p_motif_copy = p_motif;
				//Sort the results in vv_m2s_ndx by the scaffold index from the smallest to the largest(needed for avoiding a mess when inserting the fragments)
				//ToDo: add check for fragment overlapping
				std::sort(m2s_dat.v_indexes.begin(), m2s_dat.v_indexes.end(), compare_motif2scaffold_data_by_scaffold_low2high);
				TR.Debug << "Sorted fagments for stitching: ";
				core::Size currScaffoldIndex=1;
				//core::Size currMotifIndex=1;
				//Create the first residue in the pose from the scaffold
				//ToDo: catch (the many) faliure scenarios
				core::conformation::Residue tmpResi = p_scaffold_rotated.residue(1);
				p_result.append_residue_by_jump( tmpResi, currScaffoldIndex );
				++currScaffoldIndex;

				//Loop to remove LOWER and UPPER terminus from motif (should we do this before comming to this point/function?)
				for (core::Size i=1; i <= p_motif_copy.total_residue(); ++i){
					//Remove lower terminus type of the residues in the motif_copy (Needed/obligated to stitch the pieces)
					//ToDo: should we like to move this to a helper function?
					if(p_motif_copy.residue( i ).has_variant_type( core::chemical::LOWER_TERMINUS )){
						TR.Debug << "DASM says: Removing LOWER terminus from residue: "<< i << std::endl;
						core::pose::remove_variant_type_from_pose_residue(p_motif_copy, core::chemical::LOWER_TERMINUS, i);
						//As stupid as it is, if the residue is the first resdue in the chain and you change the type,
						//then you need to return the first atom (hydrogen) coordinates to the original/correct XYZ place by hand
						//ToDo: This is an Uggly Patch, we should improve/correct it.
						if(i == 1)
						{
							//This is an UUUUUggly PATCH, it is Rosetta's fault, somebody that knows how should fix it, please?
							if (p_motif.residue( i ).aa() != core::chemical::aa_pro){
								TR.Debug << "DASM says: Applying patch to 1H->H atom in the LOWER terminus of the first residue!" << std::endl;
								core::id::AtomID id_new(p_motif_copy.residue(1).atom_index("H"),1);
								p_motif_copy.set_xyz(id_new, p_motif.residue(1).xyz( "1H" ));
							}
						}
					}else if(p_motif_copy.residue(i).has_variant_type( core::chemical::UPPER_TERMINUS )){
						//Remove upper terminus type of the residues in the motif_copy (Needed/obligated to stitch the pieces)
						TR.Debug << "DASM says: Removing UPPER terminus from residue: "<< i << std::endl;
						core::pose::remove_variant_type_from_pose_residue(p_motif_copy, core::chemical::UPPER_TERMINUS, i);
					}
				}
				//Iterate the fragments
				for (core::Size i=1; i <= m2s_dat.v_indexes.size(); ++i){
					//Silly output about the fragment that we are working on
					TR.Debug <<" ["<< m2s_dat.v_indexes[i].motifLow << ","  << m2s_dat.v_indexes[i].motifHigh << "]->[" << m2s_dat.v_indexes[i].scaffoldLow << "," << m2s_dat.v_indexes[i].scaffoldHigh << "] " << std::endl; 
					//This is for the scaffold (NOTE: that the loop condition here is < not <=)					
					for (core::Size j=currScaffoldIndex; j < m2s_dat.v_indexes[i].scaffoldLow; ++j){
						//TR.Debug << "Looking for scaffold residueM: " << j << std::endl;
						core::conformation::Residue tmpResi = p_scaffold_rotated.residue(j);
						p_result.append_residue_by_bond( tmpResi, false );
					}
					//This is for the motif (NOTE: that the loop condition here is <= not <)
					for (core::Size j=m2s_dat.v_indexes[i].motifLow; j <= m2s_dat.v_indexes[i].motifHigh; ++j){
						//TR.Debug << "Looking for motif residueM: " << j << std::endl;
						//If skip_motif_extremes do not copy first and last residue in the motif fragment
						if ( skip_motif_extremes && (j==m2s_dat.v_indexes[i].motifLow) )
						{
							continue;
						}else if ( skip_motif_extremes && (j==m2s_dat.v_indexes[i].motifHigh) )
						{
							continue;
						}
						core::conformation::Residue tmpResi = p_motif_copy.residue(j);
						if ( tmpResi.is_lower_terminus() || tmpResi.is_upper_terminus() ){
						}
						p_result.append_residue_by_bond( tmpResi, false );
					}
					//set the next starting point for the scaffold at
					currScaffoldIndex=m2s_dat.v_indexes[i].scaffoldHigh+1;
				}
				//copy the last part of the scaffold (NOTE: that the loop condition here is <= not <)
				for (core::Size i=currScaffoldIndex; i <= p_scaffold_rotated.total_residue(); ++i){
					//TR.Debug << "Looking for scaffold residueF: " << i << std::endl;
					core::conformation::Residue tmpResi = p_scaffold_rotated.residue(i);
					p_result.append_residue_by_bond( tmpResi, false );
				}
				
				TR.Debug << std::endl;
				//TR.Debug << "Stitiching DONE"<< std::endl;
				return p_result;
			}
			
			/** @brief performs soperposition based on the motif and returns fragments within the RMSD_tol and
			 ** also returns the Rotation and translation superposition information in the first [1] vector of each fragment**/
			void MotifGraftMover::get_motif_scaffold_superposition_and_RMSD(
				core::pose::Pose const & p_scaffold, 
				core::pose::PoseOP const & p_motif_,
				core::Real const & RMSD_tol,
				bool b_full_motif_bb_alignment,
				utility::vector1< motif2scaffold_data > & v_m2s_data)
			{
				numeric::xyzMatrix< core::Real >  RotM;
				numeric::xyzVector< core::Real >  TvecA;
				numeric::xyzVector< core::Real >  TvecB;
				//TR.Debug << "Will iterate the fragments for alignment " << std::endl;
				//Iterate the fragments combinations
				utility::vector1< motif2scaffold_data >::iterator it_fragments = v_m2s_data.begin();
				while( it_fragments != v_m2s_data.end()) 
				{
					utility::vector1< core::Size > positions_to_alignA;
					utility::vector1< core::Size > positions_to_alignB;
					//generate a list of the positions to align and get RMSDs. A: motif, B: scaffold
					for ( core::Size i = 1; i<= (*it_fragments).v_indexes.size(); ++i )
					{
						//If we want to make Full BB alignment add all the residues in the fragments to the list
						if(b_full_motif_bb_alignment){
							TR.Debug << "Using full fragment BB alignment" << std::endl;
							if( ((*it_fragments).v_indexes[i].motifHigh-(*it_fragments).v_indexes[i].motifLow) != ((*it_fragments).v_indexes[i].scaffoldHigh - (*it_fragments).v_indexes[i].scaffoldLow ) ){
									throw utility::excn::EXCN_Msg_Exception("Something is WRONG, most likely this is a BUG. You are using full Back Bone alignment but the size of the fragments in the motif and scaffold are different, THIS alignment can't be performed. Please contact the authors of the program and cite this legend.");
								}
							//ToDo: Add a check for the same size of fragments, this alignment can't happen otherwise
							for ( core::Size j = (*it_fragments).v_indexes[i].motifLow; j<= (*it_fragments).v_indexes[i].motifHigh; ++j ){
								positions_to_alignA.push_back(j);
							}
							for ( core::Size j = (*it_fragments).v_indexes[i].scaffoldLow; j<= (*it_fragments).v_indexes[i].scaffoldHigh; ++j ){
								positions_to_alignB.push_back(j);
							}
						}else{//If not FullBB add only the tips
							TR.Debug << "Using only fragment's tips alignment" << std::endl;
							positions_to_alignA.push_back((*it_fragments).v_indexes[i].motifLow);
							positions_to_alignA.push_back((*it_fragments).v_indexes[i].motifHigh);
							positions_to_alignB.push_back((*it_fragments).v_indexes[i].scaffoldLow);
							positions_to_alignB.push_back((*it_fragments).v_indexes[i].scaffoldHigh);
						}
						//TR.Debug << "Positions to align in A(scaffold): " << (*it_fragments)[i].scaffoldLow << " " << (*it_fragments)[i].scaffoldHigh << std::endl;
						//TR.Debug << "Positions to align in B(scaffold): " << (*it_fragments)[i].motifLow << " " << (*it_fragments)[i].motifHigh << std::endl;
					}
					//The actual alignment and RMSD calculation, returns the Rotation Matrix and Traslation Vectors (first the reference frame [motif], second the desired mobile element [scaffold])
					core::Real  RMSD=MotifGraftMover::get_bb_alignment_and_transformation( *p_motif_, positions_to_alignA, p_scaffold, positions_to_alignB, RotM, TvecA, TvecB);
					TR.Debug << "RMSD of fragment = " << RMSD << std::endl;
					if (RMSD <= RMSD_tol)
					{
						TR.Debug << "Fragment passed the RMSD test" << std::endl;
						//ToDo: This [1] index is kind of adhoc, try to remove/improve it
						(*it_fragments).RotM=RotM;
						(*it_fragments).TvecA=TvecA;
						(*it_fragments).TvecB=TvecB;
						(*it_fragments).RMSD=RMSD;
						++it_fragments;
					}else
					{
						//remove from the vector the elements that don't pass the RMSD test
						it_fragments = v_m2s_data.erase(it_fragments);
					}
				}
			}
			
			/** @brief Generates a copy of the pose with the same aminoacid type in all the positions, but keeping the bb 
			 **Note: may crash or give unexpected results if you put a big residue (e.g. VAL or TRP) **/
			core::pose::Pose MotifGraftMover::get_mono_aa_pose_copy(
				core::pose::Pose const & p_input,
				std::string const & aminoacid_code)
			{
				core::pose::Pose p_mono_aa = core::pose::Pose(p_input);
				
				utility::vector1< core::Size > positions_to_replace;
				
				for( core::Size i = 1; i <= p_mono_aa.total_residue() ; ++i) {
					positions_to_replace.push_back( i );
				}
				
				protocols::toolbox::pose_manipulation::construct_poly_uniq_restype_pose( p_mono_aa, positions_to_replace, 
					core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map(aminoacid_code),
					true, true, true );
				return p_mono_aa;
			}
			
			/** @brief Generates all the discontinuous fragments combinations that are within the RMSD_tol restriction.
			 ** Reduce the combinations by matching intra chains distances by pairs.
			 ** The method is/can be exahustive but fast (i.e. iteratively it test the restrains (tree unfolding) and skips to the next combination once one fails (branch removal) ).
			 ** CAUTION, Uses self recursion, so use it wisely **/
			void MotifGraftMover::fragments_permutation_test_by_CA_distances(
				core::pose::Pose const & p_scaffold,
				core::pose::PoseOP const & p_motif_,
				utility::vector1< utility::vector1 < std::pair< core::Size, core::Size > > >  const & vv_scaffold_fragments_indexes, 
				utility::vector1< std::pair< core::Size, core::Size > > const & v_motif_fragments_indexes,
				core::Real const & RMSD_tol,
				core::Size const & start_motif_num, 
				utility::vector1< std::pair< core::Size, core::Size > > & buff_combVec,
				utility::vector1< motif2scaffold_data > & v_m2s_data)
			{
				core::Size numPrevFragments=buff_combVec.size();
				//TR.Debug << "Testing numPrevFragments = " << numPrevFragments << std::endl;
				//TR.Debug << "Target = " << vv_scaffold_fragments_indexes.size()+1 << std::endl;

				//If we reach this condition it means that there is a sucessfull combination of fragments
				if ( start_motif_num >= ( vv_scaffold_fragments_indexes.size()+1 ) )
				{
					//create the storage objects using structures
					motif2scaffold_indexes tmp_m2s_ndx;
					motif2scaffold_data tmp_m2s_data;
					TR.Debug << "Fragment combination that pass all the geometric restraints: ";
					for (core::Size i=1; i <= buff_combVec.size(); ++i){
						TR.Debug << " " << buff_combVec[i].first << "," << buff_combVec[i].second;
						//push the results in a temporary motif2scaffold_indexes structure
						tmp_m2s_ndx.motifLow=v_motif_fragments_indexes[i].first;
						tmp_m2s_ndx.motifHigh=v_motif_fragments_indexes[i].second;
						tmp_m2s_ndx.scaffoldLow=buff_combVec[i].first;
						tmp_m2s_ndx.scaffoldHigh=buff_combVec[i].second;
						//push the first structure in the second
						tmp_m2s_data.v_indexes.push_back(tmp_m2s_ndx);
					}
					TR.Debug << std::endl;
					//push the data in the return vector
					v_m2s_data.push_back(tmp_m2s_data);
					//finish with good results
					return;
				}
				//Precalculate target vs previous pair inter-fragment distance(s) (total 4 distances per pair)
				///This is not ISO safe
				//core::Real distA1B1[numPrevFragments];
				//core::Real distA1B2[numPrevFragments];
				//core::Real distA2B2[numPrevFragments];
				//core::Real distA2B1[numPrevFragments];
				///Hacky solution: use vectors. ToDo: try to improve
				utility::vector1 <core::Real> distA1B1;
				utility::vector1 <core::Real> distA1B2;
				utility::vector1 <core::Real> distA2B2;
				utility::vector1 <core::Real> distA2B1;
				for (core::Size i=1; i <= numPrevFragments; ++i)
				{
					distA1B1.push_back(0.0);
					distA1B2.push_back(0.0);
					distA2B2.push_back(0.0);
					distA2B1.push_back(0.0);
				}
				
				if( start_motif_num > 1 )
				{
					for (core::Size i=1; i <= numPrevFragments; ++i)
					{
						distA1B1[i]=( p_motif_->residue(v_motif_fragments_indexes[i].first).xyz( "CA" ) - p_motif_->residue(v_motif_fragments_indexes[start_motif_num].first).xyz( "CA" ) ).norm();
						//TR.Debug << "TestingA1B1 " << distA1B1[i] << std::endl;
						distA1B2[i]=( p_motif_->residue(v_motif_fragments_indexes[i].first).xyz( "CA" ) - p_motif_->residue(v_motif_fragments_indexes[start_motif_num].second).xyz( "CA" ) ).norm();
						//TR.Debug << "TestingA1B2 " << distA1B2[i] << std::endl;
						distA2B2[i]=( p_motif_->residue(v_motif_fragments_indexes[i].second).xyz( "CA" ) - p_motif_->residue(v_motif_fragments_indexes[start_motif_num].second).xyz( "CA" ) ).norm();
						//TR.Debug << "TestingA2B2 " << distA2B2[i] << std::endl;
						distA2B1[i]=( p_motif_->residue(v_motif_fragments_indexes[i].second).xyz( "CA" ) - p_motif_->residue(v_motif_fragments_indexes[start_motif_num].first).xyz( "CA" ) ).norm();
						//TR.Debug << "TestingA2B1 " << distA2B1[i] << std::endl;
					}
				}
				
				core::Real tmpdist;
				for (core::Size i=1; i <= vv_scaffold_fragments_indexes[start_motif_num].size(); ++i){
					//If we are dealing with fragment number >=2
					if( start_motif_num > 1 ){
						//Test overlaping and intra fragment distances, skip at the first faliure
						bool isGood=true;
						for (core::Size j=1; j <= numPrevFragments; ++j)
						{
							
							//Defacto Remove/skip  overlaping segments
							//ToDo: Should we allow swaping in the order of the fragments?
							if( ( buff_combVec[j].first <= vv_scaffold_fragments_indexes[start_motif_num][i].second ) && ( buff_combVec[j].second >= vv_scaffold_fragments_indexes[start_motif_num][i].first ) )
							{
								isGood=false;
								break;
							}
							//Match the intra distances one by one. Skip if any does not match
							tmpdist=std::abs( distA1B1[j]-( ( p_scaffold.residue(buff_combVec[j].first).xyz( "CA" )- p_scaffold.residue(vv_scaffold_fragments_indexes[start_motif_num][i].first).xyz("CA" ) ).norm() ) );
							if(tmpdist > RMSD_tol)
							{
								isGood=false;
								break;
							}
							tmpdist=std::abs( distA1B2[j]-( ( p_scaffold.residue(buff_combVec[j].first).xyz( "CA" )- p_scaffold.residue(vv_scaffold_fragments_indexes[start_motif_num][i].second).xyz("CA" ) ).norm() ) );
							if(tmpdist > RMSD_tol)
							{
								isGood=false;
								break;
							}
							tmpdist=std::abs( distA2B2[j]-( ( p_scaffold.residue(buff_combVec[j].second).xyz( "CA" )- p_scaffold.residue(vv_scaffold_fragments_indexes[start_motif_num][i].second).xyz("CA" ) ).norm() ) );
							if(tmpdist > RMSD_tol)
							{
								isGood=false;
								break;
							}
							tmpdist=std::abs( distA2B1[j]-( ( p_scaffold.residue(buff_combVec[j].second).xyz( "CA" )- p_scaffold.residue(vv_scaffold_fragments_indexes[start_motif_num][i].first).xyz("CA" ) ).norm() ) );
							if(tmpdist > RMSD_tol)
							{
								isGood=false;
								break;
							}
						}
						if(!isGood) continue;
					}
					utility::vector1< std::pair< core::Size, core::Size > > buff_combVec2 = buff_combVec;
					buff_combVec2.push_back(vv_scaffold_fragments_indexes[start_motif_num][i]);
					MotifGraftMover::fragments_permutation_test_by_CA_distances(p_scaffold, p_motif_, vv_scaffold_fragments_indexes, v_motif_fragments_indexes, RMSD_tol, start_motif_num+1, buff_combVec2, v_m2s_data);
				}
				//finish without results
				return;
			}
			
			/**@brief Generates the combinations of fragments that match the motif(s). Returns by reference:
			 **A) For each motif fragment a vector element with the indices of the motif fragments
			 **B) For each motif fragment a vector element with the indices of the fragments in the Scaffold that are within the CA distance restraint**/
			bool MotifGraftMover::get_fragments_by_CA_distances(
				core::pose::Pose const & p_scaffold, 
				core::pose::PoseOP const & p_motif_,
				utility::vector1< utility::vector1 < std::pair< core::Size, core::Size > > > & vv_scaffold_fragments_indexes,
				utility::vector1< std::pair< core::Size, core::Size > > & v_motif_fragments_indexes,
				core::Real const & RMSD_tol,
				utility::vector1 < std::pair< long int, long int > > const & max_fragment_replacement_size_delta,
				utility::vector1< std::pair< core::Size, core::Size > > const & v_motif_fragments_permutation)
			{
				
				//Stores the Ca<-> distances of each motif fragment
				utility::vector1<core::Real>  motif_distances;
				//Stores lowbound and highbound of each fragment
				std::pair< core::Size, core::Size > tmpMotifNdx;
				//Get the resd IDs for the fragments in the PDB
				//Also get the Ca<->Ca distance between the beggining and end of each fragment
				for ( core::Size i = 1; i <= p_motif_->conformation().num_chains(); ++i )
				{
					//Store the beggining and the end in motif_extremes vector1
					//tmpMotifNdx.first  = p_motif_->conformation().chain_begin(i);
					//tmpMotifNdx.second = p_motif_->conformation().chain_end(i);
					//ToDo: Replace tmpMotifNdx by v_motif_fragments_permutation[i]
					tmpMotifNdx.first=v_motif_fragments_permutation[i].first;
					tmpMotifNdx.second=v_motif_fragments_permutation[i].second;
					TR.Debug << "TEST motif chain begin: " << tmpMotifNdx.first << std::endl;
					TR.Debug << "TEST motif chain end: "   << tmpMotifNdx.second << std::endl;
					v_motif_fragments_indexes.push_back(tmpMotifNdx);
					//Store the distances in the vector motif_distances
					motif_distances.push_back( (p_motif_->residue(v_motif_fragments_indexes[i].first).xyz( "CA" )-p_motif_->residue(v_motif_fragments_indexes[i].second).xyz( "CA" )).norm() );
					TR.Debug << "TEST motif chain end to end Ca<->Ca distance: " << motif_distances[i] << std::endl;
				}
				
				core::Real tmpDist;
				//Iterate all the motif fragments...
				for ( core::Size motifNum = 1; motifNum <= v_motif_fragments_indexes.size(); ++motifNum )
				{
					//The minimum size of the scaffold fragment that we will consider (can be negative so use int)
					int tmpMinDesiredSize=(v_motif_fragments_indexes[motifNum].second-v_motif_fragments_indexes[motifNum].first)+max_fragment_replacement_size_delta[motifNum].first;
					//We cannot replace a fragment smaller than two amino acids, this scenario may happen due to the variation of the lenght in the fragments by the delta-combinations
					if(tmpMinDesiredSize<1){
						tmpMinDesiredSize=1;
					}
					TR.Debug << "For fragment # "<< motifNum << " the minimum size of a fragment in the scaffold to consider will be: " << tmpMinDesiredSize+1 << " residues" << std::endl;
					//The maximum size of the scaffold fragment that we will consider
					core::Size tmpMaxDesiredSize=(v_motif_fragments_indexes[motifNum].second-v_motif_fragments_indexes[motifNum].first)+max_fragment_replacement_size_delta[motifNum].second;
					TR.Debug << "For fragment # "<< motifNum << " the maximum size of the fragment to consider can be:  " << tmpMaxDesiredSize+1 << " residues (if the scaffold is long enough!)" << std::endl;
					//To temporarily store the resulting fragment
					utility::vector1 < std::pair< core::Size, core::Size > > tmp_fragmentWithinCutoff;
					//...to calculate all the distance pairs of CA, then return only those within the +/-tol (Do not consider the first and last residue!) 
					for ( core::Size i = 2; i<= p_scaffold.total_residue()-1; ++i )
					{
						for ( core::Size j = i+tmpMinDesiredSize; j<= p_scaffold.total_residue()-1; ++j )
						{
							//If the fragment comes longer than the size that we request then break and go to the next
							if (j-i > tmpMaxDesiredSize)
							{
								break;
							}
							//Distance of Ca<-> for the current fragment
							tmpDist=(p_scaffold.residue(i).xyz( "CA" )-p_scaffold.residue(j).xyz( "CA" )).norm();
							//TR.Debug << "TEST RMSD:" << tmpDist << std::endl;
							//minus(-) the distance of the test motif. Just to improve readeability
							tmpDist=std::abs(tmpDist-motif_distances[motifNum]);
							if(tmpDist <= RMSD_tol){
								std::pair< core::Size, core::Size > tmpPair;
								tmpPair.first=i;
								tmpPair.second=j;
								tmp_fragmentWithinCutoff.push_back(tmpPair);
								TR.Debug << "For motif #" << motifNum << " & fragment: [" << v_motif_fragments_indexes[motifNum].first << ","<< v_motif_fragments_indexes[motifNum].second << "], the TEST fragment: [" << i << "," << j << "] passed CA check, with dist=" << tmpDist << std::endl;
							}
						}
					}
					//If we do have results (>0) matching this motif fragment...
					if ( tmp_fragmentWithinCutoff.size() > 0 ){
						//...store all in the return scaffold vector
						vv_scaffold_fragments_indexes.push_back(tmp_fragmentWithinCutoff);
					}
				}//end motif fragments iterator
				//Check if we have scaffold fragment matches for all the
				if (vv_scaffold_fragments_indexes.size() != v_motif_fragments_indexes.size())
				{
					return false;
				}
				return true;
			}//end MotifGraftMover::get_fragments_by_CA_distances function
			
			/**@brief Function that performs alignment of the protein BB on the selected aminoacids. 
			 **Returns the RMSD,
			 **Returns by reference the rotation Matrix and Translation Vector,
			 **Will fail if both poses are not protein <-This can be fixed by adding a list of the atoms to align to the function, but I am not doing it now.
			 **Will fail if the number of residues to align is not the same in the two poses. **/
			core::Real MotifGraftMover::get_bb_alignment_and_transformation(
				core::pose::Pose const & poseA,
				utility::vector1< core::Size > const & positions_to_alignA,
				core::pose::Pose const & poseB,
				utility::vector1< core::Size > const & positions_to_alignB,
				numeric::xyzMatrix< core::Real > & RotM,
				numeric::xyzVector< core::Real > & TvecA,
				numeric::xyzVector< core::Real > & TvecB)
			{
				
				core::Size sizeA=positions_to_alignA.size();
				core::Size sizeB=positions_to_alignB.size();
				//Check if the align size vectors are equivalent in size, if not die with error(-1.0)
				if(sizeA != sizeB) return -1.0;
				//To store the RMSD
				core::Real RMSD;
				//To store the positions of the atoms that we want to align
				utility::vector1< numeric::xyzVector< core::Real > > matrixA;
				utility::vector1< numeric::xyzVector< core::Real > > matrixB;
				//Copy the BB atoms XYZ positions for pose A and B
				for( core::Size i = 1; i <= sizeA ; ++i){
					matrixA.push_back(poseA.residue(positions_to_alignA[i]).xyz( "CA" ));
					matrixB.push_back(poseB.residue(positions_to_alignB[i]).xyz( "CA" ));
					matrixA.push_back(poseA.residue(positions_to_alignA[i]).xyz( "C" ));
					matrixB.push_back(poseB.residue(positions_to_alignB[i]).xyz( "C" ));
					matrixA.push_back(poseA.residue(positions_to_alignA[i]).xyz( "O" ));
					matrixB.push_back(poseB.residue(positions_to_alignB[i]).xyz( "O" ));
					matrixA.push_back(poseA.residue(positions_to_alignA[i]).xyz( "N" ));
					matrixB.push_back(poseB.residue(positions_to_alignB[i]).xyz( "N" ));
				}
				//Convert the vectors to Fortran style arrays (as in Alex Ford code)
				ObjexxFCL::FArray2D< numeric::Real > FmatrixA( 3, (sizeA*4) );
				ObjexxFCL::FArray2D< numeric::Real > FmatrixB( 3, (sizeA*4) );
				//??Alex: are this post-increments particulary Roseta-coding-standars correct?
				for (core::Size i = 1; i <= (sizeA*4); i++){
					for (core::Size j = 1; j <= 3; j++){
						FmatrixA(j,i) = matrixA[i](j);
						FmatrixB(j,i) = matrixB[i](j);
						//TR.Debug << "TEST xyzMtxA: " << FmatrixA(j,i) << std::endl;
						//TR.Debug << "TEST xyzMtxB: " << FmatrixB(j,i) << std::endl;
					}
				}
				//Weighted alignment (as in Alex Ford code protocols/seeded_abinitio/util.[cc,hh])
				ObjexxFCL::FArray1D< numeric::Real > weights_fa( (sizeA*4), 1);
				
				//The actual alignment (as in Alex Ford code protocols/seeded_abinitio/util.[cc,hh])
				superposition_transform((sizeA*4), weights_fa, FmatrixA, FmatrixB, RotM, TvecA, TvecB);
				
				RMSD=0.0;
				core::Real tmpDistSqr=0.0;
				//Calculate the RMSD
				for (core::Size i = 1; i <= (sizeA*4); i++){
					//Storage for d^2
					tmpDistSqr=0.0;
					for (int j = 1; j <= 3; j++){
						tmpDistSqr+=(FmatrixA(j,i)-FmatrixB(j,i))*(FmatrixA(j,i)-FmatrixB(j,i));
					}
					RMSD+=tmpDistSqr;
				}
				RMSD=std::sqrt(RMSD/(sizeA*4));
				//Return the RMSD
				return RMSD;
			}
			
			/**@brief Superposition_transform wrapper (as in Alex Ford code protocols/toolbox/superimpose.[cc,hh])**/
			void MotifGraftMover::superposition_transform(
				core::Size natoms,
				ObjexxFCL::FArray1_double const& weights,
				ObjexxFCL::FArray2_double& ref_coords,
				ObjexxFCL::FArray2_double& coords,
				numeric::xyzMatrix< core::Real > &RotM,
				numeric::xyzVector< core::Real > &TvecA,
				numeric::xyzVector< core::Real > &TvecB)
			{
				ObjexxFCL::FArray1D_double ref_transvec( 3 );
				//Get/remove the first COM translation vector
				protocols::toolbox::reset_x( natoms, ref_coords, weights, ref_transvec );
				TvecA = numeric::xyzVector< core::Real >( ref_transvec( 1 ), ref_transvec( 2 ), ref_transvec( 3 ));
				
				//Get/remove the second COM translation vector
				ObjexxFCL::FArray1D_double transvec( 3 );
				protocols::toolbox::reset_x( natoms, coords, weights, transvec );
				TvecB = numeric::xyzVector< core::Real >( transvec( 1 ), transvec( 2 ), transvec( 3 ));
				
				//Fit centered coords, returns the rotation matrix, shifted coordinates (to the center)
				protocols::toolbox::fit_centered_coords( natoms, weights, ref_coords, coords, RotM );
			}
			
			/**@brief Function used generate a match_pose epigraft using the results stored in the MotifMatch database**/
			void MotifGraftMover::generate_match_pose(
				core::pose::Pose & target_pose, 
				MotifMatch motif_match)
			{
				//Rotate/translate and stitch the motif and the scaffold
				TR.Info << "Generating a epigraft structure with RSD score: " << motif_match.get_RMSD() << ", and clash score of: " << motif_match.get_clash_score() << std::endl;
				//get a copy of the graft data
				motif2scaffold_data tmp_fragment_data = motif_match.get_scaffold_fragment_data();
				core::pose::Pose p_frankenstein=MotifGraftMover::stich_motif_in_scaffold_by_indexes_rotation_and_translation(target_pose, *gp_p_motif_, tmp_fragment_data , false);
				
				////Remove Disulfide Cysteins
				////ToDo: Figure out a better way to do this!
				//for (core::Size i = 1; i <= p_frankenstein.total_residue(); ++i){
				//	if ( p_frankenstein.residue( i ).aa() == core::chemical::aa_cys && p_frankenstein.residue( i ).has_variant_type( core::chemical::DISULFIDE ) ){
				//		//core::pose::remove_variant_type_from_pose_residue(p_frankenstein, core::chemical::DISULFIDE, i);
				//		///core::chemical::ResidueType const & simple_cys( p_frankenstein.residue(i).residue_type_set().name_map("CYS") );
				//		///core::pose::replace_pose_residue_copying_existing_coordinates(p_frankenstein, i, simple_cys);
				//	}
				//}
				
				//Recreate disulfide bridges
				TR.Info << "Recreating Disulfide Bridges."<< std::endl;
				core::pose::initialize_disulfide_bonds(p_frankenstein);
				
				//Join the contextStructure and epigraft(frankenstein) in a single pose (reuse the container frankenstein)
				//Alex: can I reuse the pose like this? Is/will the size dynamicaly re-allocated?
				//Alex: contextStructure_ seems a global variable all around, can we encapsulate it?
				p_frankenstein=MotifGraftMover::join_two_poses_by_jump(*gp_p_contextStructure_, p_frankenstein);
				//Copy the results to the return pose
				target_pose=p_frankenstein;
				//Put data in the score table
				protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("graft_RMSD: ", motif_match.get_RMSD());
				protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("graft_clashScore: ", motif_match.get_clash_score());
				protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_string_pair("graft_motif_range: ", motif_match.get_motif_ranges());
				protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_string_pair("graft_scaffold_ranges: ", motif_match.get_scaffold_ranges());
				protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_string_pair("graft_scaffold_size_change: ", motif_match.get_scaffold2motif_size_change());
				return;
			}
			
			/**@brief Fuction to parse the class options and populate their corresponding variables to the global private space**/
			void MotifGraftMover::parse_my_string_arguments_and_cast_to_globalPrivateSpaceVariables(
						std::string const & s_contextStructure,
						std::string const & s_motif,
						core::Real  const & r_RMSD_tolerance,
						core::Real  const & r_clash_score_cutoff,
						std::string const & s_combinatory_fragment_size_delta,
						std::string const & s_max_fragment_replacement_size_delta,
						std::string const & s_clash_test_residue,
						bool        const & b_full_motif_bb_alignment)
			{
				//REQUIRED: context structure
				gp_p_contextStructure_ = core::import_pose::pose_from_pdb( s_contextStructure, false );
				
				//REQUIRED: motif structure
				gp_p_motif_ = core::import_pose::pose_from_pdb( s_motif, false );
				
				//REQUIRED: RMSD_tolerance
				gp_r_RMSD_tolerance_= r_RMSD_tolerance;
				TR.Info << "RMSD tolarance settled to: " << gp_r_RMSD_tolerance_ << std::endl;
				
				//REQUIRED: clash_score_cutoff
				gp_r_clash_score_cutoff_=r_clash_score_cutoff;
				TR.Info << "Clash score cutoff settled to: " << gp_r_clash_score_cutoff_ << std::endl;
				
				//OPTIONAL: combinatory_fragment_size_delta
				//ToDo: add checks for inconsistent inputs
				//ToDo: Maybe this parser from max_fragment_replacement_size_delta_ to vp_max_fragment_replacement_size_delta_ should be outside
				if( s_combinatory_fragment_size_delta.size() > 0){
					utility::vector1< std::string > tmpSplitParser = utility::string_split( s_combinatory_fragment_size_delta, ',' );
					for ( core::Size i = 1; i <= tmpSplitParser.size() ; ++i){
						utility::vector1< std::string > tmpSplitParser2 = utility::string_split( tmpSplitParser[i], ':' );
						std::pair< long int, long int > tmpPairParser;
						tmpPairParser.first=(long int) std::atoi(tmpSplitParser2[1].c_str());
						tmpPairParser.second=(long int) std::atoi(tmpSplitParser2[2].c_str());
						TR.Info << "Parsed combinatory size delta for fragment #"<< i << " = " << tmpPairParser.first <<"," << tmpPairParser.second << std::endl;
						gp_vp_combinatory_fragment_size_delta_.push_back(tmpPairParser);
						if (tmpPairParser.first < 0){
							throw utility::excn::EXCN_RosettaScriptsOption(" The left term [ x: ] of the combinatory fragment size delta for fragment cannot be < 0");
						}else if (tmpPairParser.second < 0){
							throw utility::excn::EXCN_RosettaScriptsOption(" The right term [ :x ] of the combinatory fragment size delta for fragment cannot be < 0");
						}
					}
					//Die if the number of delta deffinitions is different to the number of motif fragments
					if(gp_vp_combinatory_fragment_size_delta_.size() != gp_p_motif_->conformation().num_chains()){
						throw utility::excn::EXCN_RosettaScriptsOption("The number of deffined combinatory_fragment_size_delta pairs (:) must match the number of fragments in your motif");
					}
				}else
				{
					TR.Warning << "No combinatory_fragment_size_delta, assuming that you dont want to test several combinations of your motif's length" << std::endl;
					for ( core::Size i = 1; i <= gp_p_motif_->conformation().num_chains() ; ++i){
						std::pair< long int, long int > tmpPairParser;
						//The minimum number of residues that can be replaced by the code is 2
						tmpPairParser.first=0;
						//The maximum number of residues that can be replaced set to the maximum possible of that the program can manage
						tmpPairParser.second=0;
						TR.Info << "Combinatory size delta for fragment #"<< i << " set by default to: " << tmpPairParser.first <<", " << tmpPairParser.second << std::endl;
						gp_vp_combinatory_fragment_size_delta_.push_back(tmpPairParser);
					}
				}
				
				//OPTIONAL: max_fragment_replacement_size
				//ToDo: add checks for inconsistent inputs
				//ToDo: Maybe this parser from max_fragment_replacement_size_delta_ to vp_max_fragment_replacement_size_delta_ should be outside
				if( s_max_fragment_replacement_size_delta.size() > 0 )
				{
					utility::vector1< std::string > tmpSplitParser = utility::string_split( s_max_fragment_replacement_size_delta, ',' );
					for ( core::Size i = 1; i <= tmpSplitParser.size() ; ++i){
						utility::vector1< std::string > tmpSplitParser2 = utility::string_split( tmpSplitParser[i], ':' );
						std::pair< long int, long int > tmpPairParser;
						tmpPairParser.first=(long int) std::atoi(tmpSplitParser2[1].c_str());
						tmpPairParser.second=(long int) std::atoi(tmpSplitParser2[2].c_str());
						TR.Info << "Parsed replacement size delta for fragment #"<< i << " = " << tmpPairParser.first <<"," << tmpPairParser.second << std::endl;
						if (tmpPairParser.first > 0){
							throw utility::excn::EXCN_RosettaScriptsOption("The left term [ x: ] of the replacement size delta for fragment cannot be > 0");
						}else if (tmpPairParser.second < 0){
							throw utility::excn::EXCN_RosettaScriptsOption("The right term [ :x ] of the replacement size delta for fragment cannot be < 0");
						}
						gp_vp_max_fragment_replacement_size_delta_.push_back(tmpPairParser);
					}
					//Die if the number of delta deffinitions is different to the number of motif fragments
					if(gp_vp_max_fragment_replacement_size_delta_.size() != gp_p_motif_->conformation().num_chains()){
						throw utility::excn::EXCN_RosettaScriptsOption("The number of deffined max_fragment_replacement_size_delta pairs (:) must match the number of fragments in your motif");
					}
				}else
				{
					TR.Warning << "No max_fragment_replacement_size defined, assuming that each motif fragment can replace \
a fragment of any size in the scaffold " << std::endl;
					for ( core::Size i = 1; i <= gp_p_motif_->conformation().num_chains() ; ++i){
						std::pair< long int, long int > tmpPairParser;
						//The minimum number of residues that can be replaced by the code is 2
						  //tmpPairParser.first=( gp_p_motif_->conformation().chain_begin(i) )-( gp_p_motif_->conformation().chain_end(i) )+1;
						tmpPairParser.first=0;
						//The maximum number of residues that can be replaced set to the maximum possible of that the program can manage
						  //tmpPairParser.second=LONG_MAX;
						tmpPairParser.first=0;
						TR.Info << "Replacement size delta for fragment #"<< i << " set by default to: " << tmpPairParser.first <<", " << tmpPairParser.second << std::endl;
						gp_vp_max_fragment_replacement_size_delta_.push_back(tmpPairParser);
					}
					
				}
				
				//OPTIONAL: clash_test_residue
				if( s_clash_test_residue.size() > 0 ){
					gp_s_clash_test_residue_= s_clash_test_residue;
				}else
				{
					TR.Warning << "No clash_test_residue defined, the default \"GLY\" will be used";
					gp_s_clash_test_residue_= "GLY";
				}
				
				//OPTIONAL: b_full_motif_bb_alignment_
				//ToDo: Maybe this parser should be outside
				gp_b_full_motif_bb_alignment_ = b_full_motif_bb_alignment;
				if ( gp_b_full_motif_bb_alignment_ ){
					TR.Warning << "Using Full Back Bone alignment OVERRIDES any max_fragment_replacement_size_delta different to 0:0 .\
This means that using Full Back Bone alignment can only be made between motif and scaffold fragments of EXACTLY the same size. YOU are warned, IF you\
used a max_fragment_replacement_size_delta DIFFERENT to 0:0 I just will overwrite it with 0:0. If this is not what you want disable this flag." << std::endl;
					gp_vp_max_fragment_replacement_size_delta_.clear();
					for ( core::Size i = 1; i <= gp_p_motif_->conformation().num_chains() ; ++i){
						std::pair< long int, long int > tmpPairParser;
						//The minimum number of residues that can be replaced by the code is 2
						tmpPairParser.first=0;
						//The maximum number of residues that can be replaced set to the maximum possible of that the program can manage
						tmpPairParser.second=0;
						TR.Info << "Replacement size delta for fragment #"<< i << " FORCED by full_motif_bb_alignment to: " << tmpPairParser.first <<", " << tmpPairParser.second << std::endl;
						gp_vp_max_fragment_replacement_size_delta_.push_back(tmpPairParser);
					}
				}
				TR.Info << "DONE parsing global arguments" << std::endl;
			}
			
			/**@brief Fuction to parse RosettaScripts XML options**/
			void MotifGraftMover::parse_my_tag(
				utility::tag::TagPtr const tag,
				protocols::moves::DataMap &,
				protocols::filters::Filters_map const &,
				protocols::moves::Movers_map const &,
				core::pose::Pose const &)
			{
				//Generate some temporary variables to store the XML parsed arguments
				std::string s_contextStructure;
				std::string s_motif;
				core::Real  r_RMSD_tolerance;
				core::Real  r_clash_score_cutoff;
				std::string s_combinatory_fragment_size_delta;
				std::string s_max_fragment_replacement_size_delta;
				std::string s_clash_test_residue;
				bool        b_full_motif_bb_alignment;
				
				//Read XML Options
				TR.Info << "Reading XML parameters" << std::endl;
				
				//FOLLOWS: XML options :
				if( tag->hasOption("context_structure") )
				{
					s_contextStructure = tag->getOption< std::string >( "context_structure");
				}else{
					throw utility::excn::EXCN_RosettaScriptsOption("Must specify the context_structure.");
				}
				
				//REQUIRED: motif structure
				if( tag->hasOption("motif_structure") )
				{
					s_motif = tag->getOption< std::string >( "motif_structure" );
				}else
				{
					throw utility::excn::EXCN_RosettaScriptsOption("Must specify the motif_structure.");
				}
				
				//REQUIRED: RMSD_tolerance
				if( tag->hasOption("RMSD_tolerance") ){
					r_RMSD_tolerance = tag->getOption<core::Real>( "RMSD_tolerance", gp_r_RMSD_tolerance_ );
				}else
				{
					throw utility::excn::EXCN_RosettaScriptsOption("Must specify the RMSD_tolerance.");
				}
				
				//REQUIRED: clash_score_cutoff
				if( tag->hasOption("clash_score_cutoff") ){
					r_clash_score_cutoff = tag->getOption<core::Real>( "clash_score_cutoff", gp_r_clash_score_cutoff_ );
				}else
				{
					throw utility::excn::EXCN_RosettaScriptsOption("Must specify the clash_score_cutoff.");
				}
				
				//OPTIONAL: combinatory_fragment_size_delta
				if( tag->hasOption("combinatory_fragment_size_delta") ){
					s_combinatory_fragment_size_delta = tag->getOption< std::string >( "combinatory_fragment_size_delta");
				}else{
					s_combinatory_fragment_size_delta = "";
				}
				
				//OPTIONAL: max_fragment_replacement_size
				if( tag->hasOption("max_fragment_replacement_size_delta") )
				{
					s_max_fragment_replacement_size_delta = tag->getOption< std::string >( "max_fragment_replacement_size_delta");
				}else{
					s_max_fragment_replacement_size_delta = "";
				}
				
				//OPTIONAL: clash_test_residue
				if( tag->hasOption("clash_test_residue") ){
					s_clash_test_residue = tag->getOption< std::string >("clash_test_residue");
				}else
				{
					s_clash_test_residue = "";
				}
				//OPTIONAL: gp_b_full_motif_bb_alignment_
				if( tag->hasOption("full_motif_bb_alignment") ){
					b_full_motif_bb_alignment = tag->getOption< bool >("full_motif_bb_alignment");
				}else
				{
					TR.Warning << "No full_motif_bb_alignment defined, the default \"false\" (only align the tips of fragment[s]) will be used" << std::endl;
					b_full_motif_bb_alignment = false;
				}
				
				TR.Info << "DONE reading XML parmeters" << std::endl;
				//END XML Read
				
				//Parse the arguments in the global space variables
				MotifGraftMover::parse_my_string_arguments_and_cast_to_globalPrivateSpaceVariables(
					s_contextStructure,
					s_motif,
					r_RMSD_tolerance,
					r_clash_score_cutoff,
					s_combinatory_fragment_size_delta,
					s_max_fragment_replacement_size_delta,
					s_clash_test_residue,
					b_full_motif_bb_alignment);
			}
			
			/**@brief Function used by roseta to create clones of movers**/
			protocols::moves::MoverOP MotifGraftMover::clone() const
			{
				return new MotifGraftMover( *this );
			}
			
			/**@brief class instance creator**/
			protocols::moves::MoverOP MotifGraftCreator::create_mover() const
			{
				return new MotifGraftMover;
			}
			
			/**@brief Function that sets the key used to call the mover from RosettaScripts**/
			std::string MotifGraftCreator::keyname() const
			{
				return "MotifGraft";
			}
			
			/**@brief Function that sets the mover name**/
			std::string MotifGraftCreator::mover_name()
			{
				return "MotifGraftMover";
			}
			
		}//END namespace movers
		
	}//END namespace motif_grafting
	
}//END namespace protocols
