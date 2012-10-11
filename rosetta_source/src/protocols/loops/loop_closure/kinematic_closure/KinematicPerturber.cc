// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.cc
/// @brief  implementations for KinematicPerturbers used by the kinematic mover
/// @author Florian Richter, floric@u.washington.edu, march 2009
/// @author Rhiju Das, rhiju@stanford.edu, 2011 -- options of cis/trans prolines, and turn off ca bond geometry variation.
/// @author Amelie Stein, amelie.stein@ucsf.edu, Oct 2012 -- vicinity sampling refactoring & new perturbers

//Unit headers
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>

// Project headers
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>

// Rosetta Headers
#include <core/chemical/AA.hh>

#include <core/pose/Pose.hh>

#include <core/kinematics/AtomTree.hh>

#include <core/scoring/Ramachandran.hh>
#include <core/scoring/Ramachandran2B.hh>
#include <core/scoring/ScoringManager.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh> 
#include <numeric/conversions.hh>

// option key includes
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <core/kinematics/MoveMap.hh>
#include <utility/vector1.hh>

//Auto Headers



namespace protocols {
	namespace loops {
		namespace loop_closure {
			namespace kinematic_closure {
				
				static numeric::random::RandomGenerator RG(43134);
				static basic::Tracer TR("protocols.loops.loop_closure.kinematic_closure.KinematicPerturber");
				
				KinematicPerturber::KinematicPerturber()
				: max_sample_iterations_( basic::options::option[ basic::options::OptionKeys::loops::max_kic_perturber_samples ]() )
				{}
				
				KinematicPerturber::~KinematicPerturber(){}
				
				void KinematicPerturber::set_movemap( core::kinematics::MoveMapCOP mm ) { movemap_ = mm; }
				
				core::kinematics::MoveMapCOP KinematicPerturber::get_movemap() const { return movemap_; }
				
				void
				KinematicPerturber::set_pose_after_closure(
														   core::pose::Pose & pose,
														   utility::vector1<core::Real> const & torsions,
														   utility::vector1<core::Real> const &, //bond_ang,
														   utility::vector1<core::Real> const &, //bond_len,
														   bool //closure_successful
														   ) const
				{
					
					core::Size start( kinmover_->start_res() );
					
					for( core::Size res = 0; res < kinmover_->segment_length(); res++ ){
						
						pose.set_phi( start + res, torsions[ (3*(res+1)) + 1 ] );
						pose.set_psi( start + res, torsions[ (3*(res+1)) + 2 ] );
						if ( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() ) pose.set_omega( start + res, torsions[ (3*(res+1)) + 3 ] );
						
					}
				} //set_pose_after_closure
				
				///////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////TorsionSamplingKinematicPerturber////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////
				
				TorsionSamplingKinematicPerturber::TorsionSamplingKinematicPerturber( KinematicMoverCAP kinmover_in )
				: KinematicPerturber(),
				vary_ca_bond_angles_(false),
				sample_omega_for_pre_prolines_( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() ),
				rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran() )
				{ set_kinmover( kinmover_in ); }
				
				TorsionSamplingKinematicPerturber::~TorsionSamplingKinematicPerturber(){}
				
				///@details randomly varies the torsions (and possibly the bond angles) for the loop.  Respects a MoveMap, if present, for torsions.  Does NOT respect the movemap for angles; does NOT cause any sort of interactions between the MoveMap and the KinematicMover's pivot residues.
				void
				TorsionSamplingKinematicPerturber::perturb_chain(
																 core::pose::Pose const & pose,
																 utility::vector1<core::Real> & torsions,
																 utility::vector1<core::Real> & bond_ang,
																 utility::vector1<core::Real> & //bond_len
																 ) 
				{
					core::kinematics::MoveMapCOP mm(get_movemap());
					
					if( vary_ca_bond_angles_ ){
						
						core::Size pvatom3( (3* (kinmover_->segment_length() + 1)) - 1 );
						
						core::Real bangle_min( kinmover_->BANGLE_MIN() );
						core::Real bangle_sd( kinmover_->BANGLE_SD() );
						
						//what is this iterating over?
						for( Size i = 5; i <= pvatom3; i+=3 ) {
							bond_ang[ i ] = bangle_min + RG.uniform() * bangle_sd;
						}
					}
					
					
					
					core::Size tor_end = torsions.size() - 3;
					
					for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){
						//if(mm) TR << "current residue " << cur_res << "mm reports " << mm->get_bb(cur_res) << std::endl;
						
						if(!mm || mm->get_bb(cur_res)){ //if movemap is null, or (if not null) if movemap says mobile
							
							core::Real rama_phi, rama_psi;
							
							rama_.random_phipsi_from_rama(pose.aa(cur_res), rama_phi, rama_psi);
							
							torsions[i++]=rama_phi; // phi
							torsions[i++]=rama_psi; // psi
							
							i++; // leave omega alone
							
						} else {
							i += 3; //ensure i indexing increases
						}
						
					}
					
					if (  sample_omega_for_pre_prolines_ ) {
						// sample omegas. all omegas before prolines are assumed to be fair game. Note that we're not using move-map, which is phi/psi-specific.
						
						static const core::Real OMEGA_MEAN( 179.8 );
						
						for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){
							
							if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro ) {
								
								i++; //phi
								i++; //psi
								torsions[i++] = ( static_cast<int>( RG.uniform()*2 ) ? OMEGA_MEAN : 0.0 );  // flip a coin -- either 179.8 (trans) or 0.0 (cis)
								
							} else {
								i += 3;
							}
						}
						
					}
					
					
				} //perturb_chain
				
				
				void
				TorsionSamplingKinematicPerturber::set_pose_after_closure(
																		  core::pose::Pose & pose,
																		  utility::vector1<core::Real> const & torsions,
																		  utility::vector1<core::Real> const & bond_ang,
																		  utility::vector1<core::Real> const & bond_len,
																		  bool closure_successful // what is this used for?
																		  ) const
				{
					
					//	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
					parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful );
					
					if( vary_ca_bond_angles_ ){
						
						core::Real offset( 0.0 );
						for (Size res=kinmover_->start_res(), atom=5; res<= kinmover_->end_res(); res++, atom+=3) {
							
							const core::id::AtomID atomid_N (1, res);
							const core::id::AtomID atomid_CA(2, res);
							const core::id::AtomID atomid_C (3, res);
							pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
										 numeric::conversions::radians(180 - bond_ang[atom]));
							
						}
					}
				} //TorsionSamplingKinematicPerturber::set_pose_after_closure(
				
				
				///////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////VicinitySamplingKinematicPerturber////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////
				
				VicinitySamplingKinematicPerturber::VicinitySamplingKinematicPerturber( KinematicMoverCAP kinmover_in )
				: KinematicPerturber(),
				vary_ca_bond_angles_(false),
				degree_vicinity_( basic::options::option[ basic::options::OptionKeys::loops::vicinity_degree ]() ), 
				sample_omega_for_pre_prolines_( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() ) // is this respected at all?
				{ set_kinmover( kinmover_in ); }
				
				VicinitySamplingKinematicPerturber::~VicinitySamplingKinematicPerturber(){}
				
				///@details small variation around the starting phi/psi angles -- order of magnitude is determined by degree_vicinity_
				void
				VicinitySamplingKinematicPerturber::perturb_chain(
																  core::pose::Pose const & pose,
																  utility::vector1<core::Real> & torsions,
																  utility::vector1<core::Real> & bond_ang,
																  utility::vector1<core::Real> & //bond_len
																  ) 
				{
					core::kinematics::MoveMapCOP mm(get_movemap());
					
					if( vary_ca_bond_angles_ ){
						
						core::Size pvatom3( (3* (kinmover_->segment_length() + 1)) - 1 );
						
						core::Real bangle_min( kinmover_->BANGLE_MIN() );
						core::Real bangle_sd( kinmover_->BANGLE_SD() );
						
						//what is this iterating over?
						for( Size i = 5; i <= pvatom3; i+=3 ) {
							bond_ang[ i ] = bangle_min + RG.uniform() * bangle_sd;
						}
					}
					
					
					
					core::Size tor_end = torsions.size() - 3;
					
					for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){
						//if(mm) TR << "current residue " << cur_res << "mm reports " << mm->get_bb(cur_res) << std::endl;
						
						if(!mm || mm->get_bb(cur_res)){ //if movemap is null, or (if not null) if movemap says mobile
							
							core::Real rama_phi, rama_psi;
							
							rama_phi = pose.phi( cur_res ) + degree_vicinity_ * RG.gaussian();
							rama_psi = pose.psi( cur_res ) + degree_vicinity_ * RG.gaussian();
							
							torsions[i++]=rama_phi; // phi
							torsions[i++]=rama_psi; // psi
							
							i++; // leave omega alone
							
						} else {
							i += 3; //ensure i indexing increases
						}
						
					}
					
					/* [currently] no pre-pro-omega sampling in vicinity mode */		
					
				} //perturb_chain
				
				
				void
				VicinitySamplingKinematicPerturber::set_pose_after_closure(
																		   core::pose::Pose & pose,
																		   utility::vector1<core::Real> const & torsions,
																		   utility::vector1<core::Real> const & bond_ang,
																		   utility::vector1<core::Real> const & bond_len,
																		   bool closure_successful // what is this used for?
																		   ) const
				{
					
					//	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
					parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful );
					
					if( vary_ca_bond_angles_ ){
						
						core::Real offset( 0.0 );
						for (Size res=kinmover_->start_res(), atom=5; res<= kinmover_->end_res(); res++, atom+=3) {
							
							const core::id::AtomID atomid_N (1, res);
							const core::id::AtomID atomid_CA(2, res);
							const core::id::AtomID atomid_C (3, res);
							pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
										 numeric::conversions::radians(180 - bond_ang[atom]));
							
						}
					}
				} //VicinitySamplingKinematicPerturber::set_pose_after_closure(
				
				
				////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////TorsionSweepkingKinematicPerturber//////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////
				
				TorsionSweepingKinematicPerturber::TorsionSweepingKinematicPerturber()
				: KinematicPerturber()
				{}
				
				TorsionSweepingKinematicPerturber::~TorsionSweepingKinematicPerturber(){}
				
				void TorsionSweepingKinematicPerturber::set_nonpivot_res_to_sweep( utility::vector1< Size > const & resids )
				{
					nonpivot_res_to_sweep_ = resids;
				}
				
				void TorsionSweepingKinematicPerturber::set_nonpivot_bb_torsion_id( utility::vector1< Size > const & bbtorids )
				{
					assert( nonpivot_res_to_sweep_.size() == bbtorids.size() );
					sweep_torsion_ids_ = bbtorids;
				}
				
				void TorsionSweepingKinematicPerturber::set_sweep_start_angle( utility::vector1< core::Real > const & angles_in_degrees )
				{
					assert( nonpivot_res_to_sweep_.size() == angles_in_degrees.size() );
					sweep_nonpivot_torsion_starts_ = angles_in_degrees;
				}
				
				void TorsionSweepingKinematicPerturber::set_sweep_step_size( utility::vector1< core::Real > const & angle_steps_in_degrees )
				{
					assert( nonpivot_res_to_sweep_.size() == angle_steps_in_degrees.size() );
					sweep_step_sizes_ = angle_steps_in_degrees;
				}
				
				/// @details Initializes the LexicographicalIterator
				void TorsionSweepingKinematicPerturber::set_sweep_nsteps( utility::vector1< core::Size > const & nsteps )
				{
					assert( nonpivot_res_to_sweep_.size() == nsteps.size() );
					sweep_iterator_.set_dimension_sizes( nsteps );
				}
				
				
				void
				TorsionSweepingKinematicPerturber::perturb_chain(
																 core::pose::Pose const &, //pose,
																 utility::vector1<core::Real> & torsions,
																 utility::vector1<core::Real> &,// bond_ang,
																 utility::vector1<core::Real> & //bond_len
																 ) 
				{
					
					if ( sweep_iterator_.at_end() ) {
						utility_exit_with_message("TorsionSweepingKinematicPerturber asked to perturb chain even though sweep iterator is at end."); }
					
					core::Size start( kinmover_->start_res() );
					
					for ( Size ii = 1; ii <= nonpivot_res_to_sweep_.size(); ++ii ) {
						Size torsion_ind = 3 * ( nonpivot_res_to_sweep_[ ii ] - start + 1 ) + sweep_torsion_ids_[ ii ];
						torsions[ torsion_ind ] = sweep_nonpivot_torsion_starts_[ ii ] + sweep_step_sizes_[ ii ] * ( sweep_iterator_[ ii ] - 1 );
						//std::cout << " " <<  nonpivot_res_to_sweep_[ ii ] << " " << sweep_torsion_ids_[ ii ] << " " << torsion_ind << " = " << dt_ang[ torsion_ind ];
					}
					//std::cout << std::endl;
					++sweep_iterator_;
				} //TorsionSweepingKinematicPerturber::perturb_chain(
				
				
				///////////////////////////////////////////////////////////////////////////////////////////
				///////////////////// NeighborDependentTorsionSamplingKinematicPerturber //////////////////
				///////////////////////////////////////////////////////////////////////////////////////////
				
				NeighborDependentTorsionSamplingKinematicPerturber::NeighborDependentTorsionSamplingKinematicPerturber( KinematicMoverCAP kinmover_in )
				: KinematicPerturber(),
				vary_ca_bond_angles_(false),
				sample_omega_for_pre_prolines_( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() ),
				rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran2B() )
				{ set_kinmover( kinmover_in ); }
				
				NeighborDependentTorsionSamplingKinematicPerturber::~NeighborDependentTorsionSamplingKinematicPerturber(){}
				
				///@details randomly varies the torsions (and possibly the bond angles) for the loop, using phi/psi combinations based on rama2b
				void
				NeighborDependentTorsionSamplingKinematicPerturber::perturb_chain(
																				  core::pose::Pose const & pose,
																				  utility::vector1<core::Real> & torsions,
																				  utility::vector1<core::Real> & bond_ang,
																				  utility::vector1<core::Real> & //bond_len
																				  ) 
				{
					core::kinematics::MoveMapCOP mm(get_movemap());
					
					if( vary_ca_bond_angles_ ){
						
						core::Size pvatom3( (3* (kinmover_->segment_length() + 1)) - 1 );
						
						core::Real bangle_min( kinmover_->BANGLE_MIN() );
						core::Real bangle_sd( kinmover_->BANGLE_SD() );
						
						for( Size i = 5; i <= pvatom3; i+=3 ) {
							bond_ang[ i ] = bangle_min + RG.uniform() * bangle_sd;
							//TR << "replacing CA bond angle at " << (kinmover_->start_res()+int((i-4)/3)) << std::endl;
						}
					}
					
					
					
					core::Size tor_end = torsions.size() - 3;
					
					for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){
						
						if(!mm || mm->get_bb(cur_res)){ //if movemap is null, or (if not null) if movemap says mobile
							
							core::Real rama_phi, rama_psi;
							
							// warning: not safe for terminal positions -- however, KIC shouldn't actually include any termini anyway... 
							
							// currently we don't really have data for both neighbors together
							// -- for now, do a coin flip on which one to use, though later we should implement this such that each side has an individual perturber, and we call both with equal likelihood (should have fewer ifs)
							static_cast<int>( RG.uniform()*2 ) ? rama_.random_phipsi_from_rama_left(pose.aa(cur_res-1), pose.aa(cur_res),rama_phi, rama_psi) : rama_.random_phipsi_from_rama_right(pose.aa(cur_res), pose.aa(cur_res+1), rama_phi, rama_psi);
							
							torsions[i++]=rama_phi; // phi
							torsions[i++]=rama_psi; // psi
							
							i++; // leave omega alone
							
						} else {
							i += 3; //ensure i indexing increases
						}
						
					}
					
					if (  sample_omega_for_pre_prolines_ ) {
						// sample omegas. all omegas before prolines are assumed to be fair game. Note that we're not using move-map, which is phi/psi-specific.
						
						static const core::Real OMEGA_MEAN( 179.8 );
						
						for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){
							
							if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro ) {
								
								core::Real rand_omega = ( static_cast<int>( RG.uniform()*2 ) ? OMEGA_MEAN : 0.0 );  // flip a coin -- either 179.8 (trans) or 0.0 (cis)
								
								i++; //phi
								i++; //psi
								torsions[i++] = rand_omega;
								
							} else {
								i += 3;
							}
						}
						
					}
					
					
				} //perturb_chain
				
				
				
				void
				NeighborDependentTorsionSamplingKinematicPerturber::set_pose_after_closure(
																						   core::pose::Pose & pose,
																						   utility::vector1<core::Real> const & torsions,
																						   utility::vector1<core::Real> const & bond_ang,
																						   utility::vector1<core::Real> const & bond_len,
																						   bool closure_successful // what is this used for?
																						   ) const
				{
					
					//	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
					parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful );
					
					if( vary_ca_bond_angles_ ){
						
						core::Real offset( 0.0 );
						for (Size res=kinmover_->start_res(), atom=5; res<= kinmover_->end_res(); res++, atom+=3) {
							
							const core::id::AtomID atomid_N (1, res);
							const core::id::AtomID atomid_CA(2, res);
							const core::id::AtomID atomid_C (3, res);
							pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
										 numeric::conversions::radians(180 - bond_ang[atom]));
							
						}
					}
				} //NeighborDependentTorsionSamplingKinematicPerturber::set_pose_after_closure(
				
				
				
				
				
				
				///////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////TorsionRestrictedKinematicPerturber////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////
				
				TorsionRestrictedKinematicPerturber::TorsionRestrictedKinematicPerturber( KinematicMoverCAP kinmover_in, std::string torsion_bins )
				: KinematicPerturber(),
				vary_ca_bond_angles_(false),
				sample_omega_for_pre_prolines_( true ), // if the torsion string is lowercase we use cis, otherwise trans -- note that this flag is used by the parent::set_pose_after_closure function and thus needs to be set to true, otherwise omega will be ignored -- this should be fixed in the future though
				rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran() )
				{ set_kinmover( kinmover_in ); 
					predefined_torsions_ = torsion_bins; // store torsion string provided by the user or derived from the native structure
					//std::cerr << " predefined torsion string: " << predefined_torsions_ << std::endl;
				}	
				
				TorsionRestrictedKinematicPerturber::~TorsionRestrictedKinematicPerturber(){}
				
				///@details randomly varies the torsions within the given torsion bin 
				void
				TorsionRestrictedKinematicPerturber::perturb_chain(
																   core::pose::Pose const & pose,
																   utility::vector1<core::Real> & torsions,
																   utility::vector1<core::Real> & bond_ang,
																   utility::vector1<core::Real> & //bond_len
																   ) 
				{
					core::kinematics::MoveMapCOP mm(get_movemap());
					
					if( vary_ca_bond_angles_ ){ // is this appropriate here?
						
						core::Size pvatom3( (3* (kinmover_->segment_length() + 1)) - 1 );
						
						core::Real bangle_min( kinmover_->BANGLE_MIN() );
						core::Real bangle_sd( kinmover_->BANGLE_SD() );
						
						for( Size i = 5; i <= pvatom3; i+=3 ) {
							bond_ang[ i ] = bangle_min + RG.uniform() * bangle_sd;
							//TR << "replacing CA bond angle at " << (kinmover_->start_res()+int((i-4)/3)) << std::endl;
						}
					}
					
					
					
					core::Size tor_end = torsions.size() - 3;
					
					core::Size torsion_string_offset = kinmover_->start_res() - kinmover_->loop_begin(); //offset to fetch the torsion bin for the subsegment we're currently sampling -- the string is always for the entire loop -- maybe we should make sure that this doesn't generate a segmentation fault?
					
					for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){
						if(!mm || mm->get_bb(cur_res)){ //if movemap is null, or (if not null) if movemap says mobile
							
							core::Real rama_phi, rama_psi;
							
							rama_.random_phipsi_from_rama_by_torsion_bin(pose.aa(cur_res), rama_phi, rama_psi, toupper(predefined_torsions_[(i-4)/3 + torsion_string_offset])); //here we only care about the phi/psi-based torsion bin, omega is set below
							
							torsions[i++]=rama_phi; // phi
							torsions[i++]=rama_psi; // psi
							
							i++; // leave omega alone
							
						} else {
							i += 3; //ensure i indexing increases
						}
						
					}
					
					// the TorsionRestricted mover will automatically derive cis/trans from the torsion bin string, it [currently] does not care about external flags -- however, note that set_pose_after_closure(...) directly reads the command-line flag, so that it must be set to true in order to accept the cis torsions
					
					static const core::Real OMEGA_MEAN( 179.8 );
					
					for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){
						
						if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro ) {
							
							core::Real rand_omega; // = ( static_cast<int>( RG.uniform()*2 ) ? OMEGA_MEAN : 0.0 );  // flip a coin -- either 179.8 (trans) or 0.0 (cis)
							
							if ( islower(predefined_torsions_[(i-4)/3 + torsion_string_offset]) ) {
								rand_omega = 0; // lowercase is for cis		
								
								//std::cerr << "cis proline at " << cur_res << std::endl; // debug
								
							} else {
								rand_omega = OMEGA_MEAN;
							}
							//}
							
							i++; //phi
							i++; //psi
							torsions[i++] = rand_omega;
							
						} else {
							i += 3;
						}
					}	
					
					//}
					
					
				} //perturb_chain
				
				
				
				void
				TorsionRestrictedKinematicPerturber::set_pose_after_closure(
																			core::pose::Pose & pose,
																			utility::vector1<core::Real> const & torsions,
																			utility::vector1<core::Real> const & bond_ang,
																			utility::vector1<core::Real> const & bond_len,
																			bool closure_successful // what is this used for?
																			) const
				{
					
					//	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
					parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful );
					
					if( vary_ca_bond_angles_ ){
						
						core::Real offset( 0.0 );
						for (Size res=kinmover_->start_res(), atom=5; res<= kinmover_->end_res(); res++, atom+=3) {
							
							const core::id::AtomID atomid_N (1, res);
							const core::id::AtomID atomid_CA(2, res);
							const core::id::AtomID atomid_C (3, res);
							pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
										 numeric::conversions::radians(180 - bond_ang[atom]));
							
						}
					}
				} //TorsionRestrictedKinematicPerturber::set_pose_after_closure(
				
				
				
				
				///////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////// TabooSamplingKinematicPerturber ////////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////
				
				// constructor in case we don't know the sequence (e.g. when remodeling multiple loops)
				// note that you MUST initialize the sequence before calling apply (else segfaults guaranteed)
				TabooSamplingKinematicPerturber::TabooSamplingKinematicPerturber( KinematicMoverCAP kinmover_in )
				: KinematicPerturber(),
				vary_ca_bond_angles_(false),
				sample_omega_for_pre_prolines_( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() ),
				num_strings_(100000), // large numbers could be inefficient here if there are lots of design steps (which invalidate the current torsion strings & taboo map) -- however, generation of the random strings seems to be very fast
				rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran() )
				{ set_kinmover( kinmover_in ); 
				}	
				
				
				TabooSamplingKinematicPerturber::~TabooSamplingKinematicPerturber(){}
				
				void
				TabooSamplingKinematicPerturber::refill_torsion_string_vector() 
				{
					utility::vector1< std::string > torsion_bins_per_position;		// data structure to hold the vector of torsion bins, with frequencies dependent on the amino acid type, for each position in the loop
					utility::vector1< core::chemical::AA > loop_sequence = kinmover_->get_loop_sequence(); 
					torsion_bins_per_position.resize(loop_sequence.size()); // check if we need +1 here somewhere... 
					core::Real ideal_freq, current_freq;
					for (core::Size i = 1; i <= loop_sequence.size(); i++) {
						std::map< char, core::Size > entries_per_torsion_bin;
						rama_.get_entries_per_torsion_bin( loop_sequence[i], entries_per_torsion_bin ); // to be written
						std::string torsion_bins_for_pos;
						// then iterate over the map and fill the string accordingly
						for (std::map< char, core::Size >::const_iterator mcs_i = entries_per_torsion_bin.begin(); 
							 mcs_i != entries_per_torsion_bin.end(); mcs_i++) {
							if (mcs_i->first == 'X') 
								continue; // we only keep X in here for the totals
							ideal_freq = core::Real(mcs_i->second)/entries_per_torsion_bin['X'];
							current_freq = kinmover_->frequency_in_taboo_map( i-1, mcs_i->first ); // calculated on a string, i.e., base is 0
							
							//TR << "for residue " << loop_sequence[i] << ", torsion bin " << mcs_i->first << " has " << mcs_i->second << " entries --> " << ideal_freq << " while current freq is " << current_freq << std::endl;
							
							if (current_freq > 0)
								ideal_freq /= current_freq; // adjust frequency for what has been sampled already -- note that this means that the size of the resulting string isn't necessarily num_strings_ any more
							std::string this_bin = std::string( ceil(ideal_freq*num_strings_), mcs_i->first ); 
							torsion_bins_for_pos += this_bin;
						}
						while (torsion_bins_for_pos.size() < num_strings_) { // make sure they all have the same size -- with adjustments the size may vary (slightly)
							torsion_bins_for_pos += "X"; 
						}
						numeric::random::random_permutation(torsion_bins_for_pos.begin(), torsion_bins_for_pos.end(), numeric::random::RG); 
						torsion_bins_per_position[i] = torsion_bins_for_pos;
						//if (i == 1) 
						//	TR << torsion_bins_for_pos << std::endl; // debug
					}
					
					for (core::Size j = 1; j < num_strings_; j++) {
						
						std::string new_torsion_string; // = new std::string(sequence_.size()); // initialize to correct size
						for (core::Size i = 1; i <= loop_sequence.size(); i++) {
							new_torsion_string += torsion_bins_per_position[i][j]; // make sure this is correct... though it will happily segfault if not
						}
						// check if we have already sampled this string -- if not, add it to the list
						if ( !( kinmover_->is_in_taboo_map(new_torsion_string) ) ) 
							random_torsion_strings_.push_back(new_torsion_string);
						
					}
				}
				
				
				
				std::string
				TabooSamplingKinematicPerturber::next_torsion_string() // core::pose::Pose const & pose ) 
				{
					
					/*
					 // check if the "stack" still contains a string -- if so, return that
					 // if not, generate new torsion strings that should cover the entire space (ha!)
					 
					 -- for each position in the loop, fetch the corresponding residues, and then ask the Rama lookup table for the relative populations of each torsion bin, and generate a string of length n with the frequencies of each letter corresponding to the respective torsion bin's frequency for this residue
					 -- maybe the function to generate a corresponding string can even be part of the Rama function?
					 -- fill up with X if the strings aren't long enough -- they must all have the same length
					 -- for later applications we could think about adjusting the torsion bin frequencies by those that have already been sampled --> we'd need another structure that holds the strings that were already sampled (might be moved there directly the moment they're provided by this function)
					 -- randomly shuffle each string
					 -- generate n random torsion strings for the loop by appending positions 1,2,3... of each of the now randomized strings [for the positions] -- this reflects the respective torsion bin frequencies (because the initial strings do, see above)
					 */
					
					
					if (random_torsion_strings_.size() == 0) {
						random_torsion_strings_.clear();
						refill_torsion_string_vector();
						
						while (random_torsion_strings_.size() == 0) { // when the taboo map gets fuller, it is possible that all 1000 random strings have already been tried, and thus the map is still empty even after "refilling"
							refill_torsion_string_vector();
						}
					}
					
					
					runtime_assert(random_torsion_strings_.size() > 0); // shouldn't happen any more... 
					
					// make sure the torsion string we're returning hasn't been tested yet
					while ( kinmover_->is_in_taboo_map(random_torsion_strings_[random_torsion_strings_.size()]) ) {
						//TR << random_torsion_strings_[random_torsion_strings_.size()] << " has already been tested, next... " << std::endl;
						random_torsion_strings_.pop_back();
						while (random_torsion_strings_.size() == 0) // refill if necessary
							refill_torsion_string_vector();
						
					}	
					
					std::string tb = random_torsion_strings_[random_torsion_strings_.size()];
					//TR << tb << " " << random_torsion_strings_.size() <<  std::endl; // debug	
					
					random_torsion_strings_.pop_back(); 
					
					return tb;
				}
				
				void
				TabooSamplingKinematicPerturber::perturb_chain(
															   core::pose::Pose const & pose,
															   utility::vector1<core::Real> & torsions,
															   utility::vector1<core::Real> & bond_ang,
															   utility::vector1<core::Real> & //bond_len
															   ) 
				{
					core::kinematics::MoveMapCOP mm(get_movemap());
					
					if( vary_ca_bond_angles_ ){ // is this appropriate here?
						
						core::Size pvatom3( (3* (kinmover_->segment_length() + 1)) - 1 );
						
						core::Real bangle_min( kinmover_->BANGLE_MIN() );
						core::Real bangle_sd( kinmover_->BANGLE_SD() );
						
						//what is this iterating over?
						for( Size i = 5; i <= pvatom3; i+=3 ) {
							bond_ang[ i ] = bangle_min + RG.uniform() * bangle_sd;
							//TR << "replacing CA bond angle at " << (kinmover_->start_res()+int((i-4)/3)) << std::endl;
						}
					}
					
					
					
					core::Size tor_end = torsions.size() - 3;
					
					std::string torsion_string = next_torsion_string(); // function to provide the torsion string to be sampled -- to be written
					
					core::Size torsion_string_offset = kinmover_->start_res() - kinmover_->loop_begin();  //offset to fetch the torsion bin for the subsegment we're currently sampling -- the string is always for the entire loop -- maybe we should make sure that this doesn't generate a segmentation fault?
					//std::cerr << "Torsion string is " << torsion_string << ", offset is " << torsion_string_offset << std::endl;
					
					for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){
						if(!mm || mm->get_bb(cur_res)){ //if movemap is null, or (if not null) if movemap says mobile
							
							core::Real rama_phi, rama_psi;
							
							rama_.random_phipsi_from_rama_by_torsion_bin(pose.aa(cur_res), rama_phi, rama_psi, toupper(torsion_string[(i-4)/3 + torsion_string_offset])); //here we only care about the phi/psi-based torsion bin, omega is set below
							
							torsions[i++]=rama_phi; // phi
							torsions[i++]=rama_psi; // psi
							
							i++; // leave omega alone
							
						} else {
							i += 3; //ensure i indexing increases
						}
						
					}
					
					if (  sample_omega_for_pre_prolines_ ) {
						// sample omegas. all omegas before prolines are assumed to be fair game. Note that we're not using move-map, which is phi/psi-specific.
						
						static const core::Real OMEGA_MEAN( 179.8 );
						
						for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){
							
							if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro ) {
								
								core::Real rand_omega = ( static_cast<int>( RG.uniform()*2 ) ? OMEGA_MEAN : 0.0 );  // flip a coin -- either 179.8 (trans) or 0.0 (cis)
								
								/*
								 if ( islower(torsion_string[(i-4)/3 + torsion_string_offset]) ) { // actually at the moment this is impossible... use random value instead?
								 rand_omega = 0; // lowercase is for cis		
								 } else {
								 rand_omega = OMEGA_MEAN;
								 }
								 //}
								 */
								i++; //phi
								i++; //psi
								torsions[i++] = rand_omega;
								
							} else {
								i += 3;
							}
						}
						
					}
					
					
				} //perturb_chain
				
				
				
				void
				TabooSamplingKinematicPerturber::set_pose_after_closure(
																		core::pose::Pose & pose,
																		utility::vector1<core::Real> const & torsions,
																		utility::vector1<core::Real> const & bond_ang,
																		utility::vector1<core::Real> const & bond_len,
																		bool closure_successful // what is this used for?
																		) const
				{
					
					//	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
					parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful );
					
					if( vary_ca_bond_angles_ ){
						
						core::Real offset( 0.0 );
						for (Size res=kinmover_->start_res(), atom=5; res<= kinmover_->end_res(); res++, atom+=3) {
							
							const core::id::AtomID atomid_N (1, res);
							const core::id::AtomID atomid_CA(2, res);
							const core::id::AtomID atomid_C (3, res);
							pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
										 numeric::conversions::radians(180 - bond_ang[atom]));
							
						}
					}
				} //TabooSamplingKinematicPerturber::set_pose_after_closure(
				
				
				
				
				
				
				
				///////////////////////////////////////////////////////////////////////////////////////////
				//////////////////// NeighborDependentTabooSamplingKinematicPerturber /////////////////////
				///////////////////////////////////////////////////////////////////////////////////////////
				
				NeighborDependentTabooSamplingKinematicPerturber::NeighborDependentTabooSamplingKinematicPerturber( KinematicMoverCAP kinmover_in )
				: KinematicPerturber(),
				vary_ca_bond_angles_(false),
				sample_omega_for_pre_prolines_( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() ),
				num_strings_(1000), 
				rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran2B() )
				{ set_kinmover( kinmover_in ); 
				}	
				
				NeighborDependentTabooSamplingKinematicPerturber::~NeighborDependentTabooSamplingKinematicPerturber(){}
				
				void
				NeighborDependentTabooSamplingKinematicPerturber::refill_torsion_string_vector() // the amount of code duplication in here is a shame... can't we have a generic refill_torsion_string_vector() function that takes the respective map as an input? Which class would have to implement this, the parent class? Maybe it'd be good to have a specific parent class for TabooSampling... 
				{
					utility::vector1< std::string > torsion_bins_per_position;		// data structure to hold the vector of torsion bins, with frequencies dependent on the amino acid type, for each position in the loop
					utility::vector1< core::chemical::AA > loop_sequence = kinmover_->get_loop_sequence(); 
					torsion_bins_per_position.resize(loop_sequence.size()); 
					core::Real ideal_freq, current_freq;
					for (core::Size i = 1; i <= loop_sequence.size(); i++) {
						std::map< char, core::Size > entries_per_torsion_bin;
						
						
						if (i == 1) // we don't know the previous position (and this position is solved analytically anyway)
							rama_.get_entries_per_torsion_bin_right( loop_sequence[i], loop_sequence[i+1], entries_per_torsion_bin );
						else if (i == loop_sequence.size()) 
							rama_.get_entries_per_torsion_bin_left( loop_sequence[i-1], loop_sequence[i], entries_per_torsion_bin );
						else {
							// take the minimum count from both sides -- especially if one of them is 0 we cannot ask for a random phi/psi combination as it won't exist
							std::map< char, core::Size > left_tb, right_tb;
							rama_.get_entries_per_torsion_bin_left( loop_sequence[i-1], loop_sequence[i], left_tb );
							rama_.get_entries_per_torsion_bin_right( loop_sequence[i], loop_sequence[i+1], right_tb );	
							// now iterate over both maps and take the minimum	
							for (std::map< char, core::Size >::const_iterator mcs_i = left_tb.begin(); mcs_i != left_tb.end(); mcs_i++) {
								entries_per_torsion_bin[mcs_i->first] = std::min(left_tb[mcs_i->first], right_tb[mcs_i->first]);
							}
						}
						
						
						std::string torsion_bins_for_pos;
						// then iterate over the map and fill the string accordingly
						for (std::map< char, core::Size >::const_iterator mcs_i = entries_per_torsion_bin.begin(); 
							 mcs_i != entries_per_torsion_bin.end(); mcs_i++) {
							if (mcs_i->first == 'X') 
								continue; // we only keep X in here for the totals
							ideal_freq = core::Real(mcs_i->second)/entries_per_torsion_bin['X'];
							current_freq = kinmover_->frequency_in_taboo_map( i-1, mcs_i->first ); // calculated on a string, i.e., base is 0
							
							//TR << "for residue " << core::chemical::AA(loop_sequence[i]) << ", torsion bin " << mcs_i->first << " has " << mcs_i->second << " entries --> " << ideal_freq << " while current freq is " << current_freq << std::endl;
							
							if (current_freq > 0)
								ideal_freq /= current_freq; // adjust frequency for what has been sampled already -- note that this means that the size of the resulting string isn't necessarily num_strings_ any more
							std::string this_bin = std::string( ceil(ideal_freq*num_strings_), mcs_i->first ); 
							torsion_bins_for_pos += this_bin;
						}
						while (torsion_bins_for_pos.size() < num_strings_) { // make sure they all have the same size -- with adjustments the size may vary (slightly)
							torsion_bins_for_pos += "X"; 
						}
						numeric::random::random_permutation(torsion_bins_for_pos.begin(), torsion_bins_for_pos.end(), numeric::random::RG); 
						torsion_bins_per_position[i] = torsion_bins_for_pos;
					}
					
					for (core::Size j = 1; j < num_strings_; j++) {
						
						std::string new_torsion_string;
						for (core::Size i = 1; i <= loop_sequence.size(); i++) {
							new_torsion_string += torsion_bins_per_position[i][j];
						}
						// check if we have already sampled this string -- if not, add it to the list
						if ( !( kinmover_->is_in_taboo_map(new_torsion_string) ) ) 
							random_torsion_strings_.push_back(new_torsion_string);
					}
				}
				
				
				
				
				std::string
				NeighborDependentTabooSamplingKinematicPerturber::next_torsion_string() 
				{
					
					
					if (random_torsion_strings_.size() == 0) {
						random_torsion_strings_.clear();
						refill_torsion_string_vector();
						
						while (random_torsion_strings_.size() == 0) { // when the taboo map gets fuller, it is possible that all 1000 random strings have already been tried, and thus the map is still empty even after "refilling"
							refill_torsion_string_vector();
						}
					}
					
					
					runtime_assert(random_torsion_strings_.size() > 0); 		
					// make sure the torsion string we're returning hasn't been tested yet
					while ( kinmover_->is_in_taboo_map(random_torsion_strings_[random_torsion_strings_.size()]) ) {
						//TR << random_torsion_strings_[random_torsion_strings_.size()] << " has already been tested, next... " << std::endl;
						random_torsion_strings_.pop_back();
						while (random_torsion_strings_.size() == 0) // refill if necessary
							refill_torsion_string_vector();
						
					}	
					
					std::string tb = random_torsion_strings_[random_torsion_strings_.size()];
					
					random_torsion_strings_.pop_back(); 
					
					return tb;
				}
				
				void
				NeighborDependentTabooSamplingKinematicPerturber::perturb_chain(
																				core::pose::Pose const & pose,
																				utility::vector1<core::Real> & torsions,
																				utility::vector1<core::Real> & bond_ang,
																				utility::vector1<core::Real> & //bond_len
																				) 
				{
					core::kinematics::MoveMapCOP mm(get_movemap());
					
					if( vary_ca_bond_angles_ ){ 
						
						core::Size pvatom3( (3* (kinmover_->segment_length() + 1)) - 1 );
						
						core::Real bangle_min( kinmover_->BANGLE_MIN() );
						core::Real bangle_sd( kinmover_->BANGLE_SD() );
						
						for( Size i = 5; i <= pvatom3; i+=3 ) {
							bond_ang[ i ] = bangle_min + RG.uniform() * bangle_sd;
							//TR << "replacing CA bond angle at " << (kinmover_->start_res()+int((i-4)/3)) << std::endl;
						}
					}
					
					
					
					core::Size tor_end = torsions.size() - 3;
					
					std::string torsion_string = next_torsion_string(); 
					
					core::Size torsion_string_offset = kinmover_->start_res() - kinmover_->loop_begin();  //offset to fetch the torsion bin for the subsegment we're currently sampling -- the string is always for the entire loop -- maybe we should make sure that this doesn't generate a segmentation fault?
					//std::cerr << "Torsion string is " << torsion_string << ", offset is " << torsion_string_offset << std::endl;
					
					for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){
						if(!mm || mm->get_bb(cur_res)){ //if movemap is null, or (if not null) if movemap says mobile
							
							core::Real rama_phi, rama_psi;
							
							// coin flip: use information from left or right AA?
							// at either end of the loop, use the "inner" neighbor -- mainly because we haven't checked whether there actually is data for the outer neighbor and the corresponding bin [as we don't have access to those positions during setup]
							// drawback of this approach: we don't consider the residues directly adjacent to the loop
							// -- however, the boundary positions are solved analytically anyway, so this shouldn't matter, practically speaking
							if ( cur_res == kinmover_->loop_begin() ) 
								rama_.random_phipsi_from_rama_by_torsion_bin_right(pose.aa(cur_res), pose.aa(cur_res+1), rama_phi, rama_psi, toupper(torsion_string[(i-4)/3 + torsion_string_offset]));
							else if ( cur_res == kinmover_->loop_end() )
								rama_.random_phipsi_from_rama_by_torsion_bin_left(pose.aa(cur_res-1), pose.aa(cur_res), rama_phi, rama_psi, toupper(torsion_string[(i-4)/3 + torsion_string_offset]));
							else
								static_cast<int>( RG.uniform()*2 ) ? rama_.random_phipsi_from_rama_by_torsion_bin_left(pose.aa(cur_res-1), pose.aa(cur_res), rama_phi, rama_psi, toupper(torsion_string[(i-4)/3 + torsion_string_offset])) : rama_.random_phipsi_from_rama_by_torsion_bin_right(pose.aa(cur_res), pose.aa(cur_res+1), rama_phi, rama_psi, toupper(torsion_string[(i-4)/3 + torsion_string_offset])); //here we only care about the phi/psi-based torsion bin, omega is set below
							
							torsions[i++]=rama_phi; // phi
							torsions[i++]=rama_psi; // psi
							
							i++; // leave omega alone
							
						} else {
							i += 3; //ensure i indexing increases
						}
						
					}
					
					if (  sample_omega_for_pre_prolines_ ) { // this could also be implemented in a parent class -- or is there a specific reason it isn't? 
						// sample omegas. all omegas before prolines are assumed to be fair game. Note that we're not using move-map, which is phi/psi-specific.
						
						static const core::Real OMEGA_MEAN( 179.8 );
						
						for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){
							
							if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro ) {
								
								core::Real rand_omega = ( static_cast<int>( RG.uniform()*2 ) ? OMEGA_MEAN : 0.0 );  // flip a coin -- either 179.8 (trans) or 0.0 (cis)
								
								i++; //phi
								i++; //psi
								torsions[i++] = rand_omega;
								
							} else {
								i += 3;
							}
						}
						
					}
					
					
				} //perturb_chain
				
				
				
				void
				NeighborDependentTabooSamplingKinematicPerturber::set_pose_after_closure(
																						 core::pose::Pose & pose,
																						 utility::vector1<core::Real> const & torsions,
																						 utility::vector1<core::Real> const & bond_ang,
																						 utility::vector1<core::Real> const & bond_len,
																						 bool closure_successful // what is this used for?
																						 ) const
				{
					
					//	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
					parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful );
					
					if( vary_ca_bond_angles_ ){
						
						core::Real offset( 0.0 );
						for (Size res=kinmover_->start_res(), atom=5; res<= kinmover_->end_res(); res++, atom+=3) {
							
							const core::id::AtomID atomid_N (1, res);
							const core::id::AtomID atomid_CA(2, res);
							const core::id::AtomID atomid_C (3, res);
							pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
										 numeric::conversions::radians(180 - bond_ang[atom]));
							
						}
					}
				} //NeighborDependentTabooSamplingKinematicPerturber::set_pose_after_closure(
				
				
				
				
			} // namespace kinematic_closure
		} // namespace loop_closure
	} // namespace loops
} // namespace protocols
