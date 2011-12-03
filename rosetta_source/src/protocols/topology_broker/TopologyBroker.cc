// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @detailed responsibilities:
///           maintains list of ToplogyClaimers
///           maintains DofClaims -- exclusive or non-exclusively markedup dofs like BackboneClaim, IntraResClaim, JumpClaim
///           generates FoldTree, MoveMap, and collects samplers provided by TopologyClaimers
/// @author Oliver Lange

// Unit Headers
#include <protocols/topology_broker/TopologyBroker.hh>

// Package Headers
#include <protocols/topology_broker/Exceptions.hh>
#include <protocols/topology_broker/ClaimerMessage.hh>
#include <protocols/topology_broker/TopologyClaimer.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <protocols/jd2/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/kinematics/Exceptions.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/scoring/dssp/StrandPairing.hh>
#include <protocols/moves/MoverContainer.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// for symmetry
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <basic/Tracer.hh>

// C++ headers
// AUTO-REMOVED #include <iterator>
#include <sstream>
#include <vector>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace ObjexxFCL { } using namespace ObjexxFCL;

static basic::Tracer tr("protocols.topo_broker", basic::t_info);

namespace protocols {
namespace topology_broker {

using namespace utility::excn;
using namespace core;

TopologyBroker::TopologyBroker() :
	nres_( 0 ),
	fold_tree_( NULL ),
	final_fold_tree_( NULL ),
	repack_scorefxn_( NULL ),
	bUseJobPose_( false ),
	current_pose_( NULL )
{}

TopologyBroker::~TopologyBroker() {}

TopologyBroker::TopologyBroker( const TopologyBroker& tp ) :
	ReferenceCount(),
	claimers_( tp.claimers_ )
{
	current_claims_ = tp.current_claims_;
	nres_ = tp.nres_;
	fold_tree_ = tp.fold_tree_;
	final_fold_tree_ = tp.final_fold_tree_;
	repack_scorefxn_ = tp.repack_scorefxn_;
	to_be_closed_cuts_ = tp.to_be_closed_cuts_;
	start_pose_cuts_ = tp.start_pose_cuts_;
	current_pose_ = tp.current_pose_;
}

TopologyBroker const& TopologyBroker::operator = ( TopologyBroker const& src ) {
	if ( this != &src ) {
		claimers_ = src.claimers_;
		current_claims_ = src.current_claims_;
		nres_ = src.nres_;
		fold_tree_ = src.fold_tree_;
		final_fold_tree_ = src.final_fold_tree_;
		repack_scorefxn_ = src.repack_scorefxn_;
		to_be_closed_cuts_ = src.to_be_closed_cuts_;
		start_pose_cuts_ = src.start_pose_cuts_;
		current_pose_ = src.current_pose_;
	}
	return *this;
}

void TopologyBroker::add( TopologyClaimerOP cl ) {
	claimers_.push_back( cl );
	cl->set_broker( this );
}

void TopologyBroker::generate_sequence_claims( DofClaims& all_claims ) {
	for ( TopologyClaimers::iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		(*top)->generate_sequence_claims( all_claims );
	}
	tr.Trace << "All sequence claims: \n" << all_claims << std::endl;
}

void TopologyBroker::relay_message( ClaimerMessage& msg ) const {
	for ( TopologyClaimers::const_iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		if ( msg.matches( (*top)->label() ) ) {
			(*top)->receive_message( msg );
		}
	}
}

void TopologyBroker::generate_round1( DofClaims& all_claims ) {
	for ( TopologyClaimers::iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		tr.Trace << "generate claim for " << (*top)->type() << std::endl;
		(*top)->generate_claims( all_claims );
	}
	tr.Trace << "All round1 claims: \n" << all_claims << std::endl;
}

void TopologyBroker::generate_final_claims( DofClaims& all_claims ) {
	for ( TopologyClaimers::iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		(*top)->finalize_claims( all_claims );
	}
	tr.Trace << "all final claims: \n " << all_claims << std::endl;
}

core::fragment::FragSetCOP TopologyBroker::loop_frags( core::kinematics::MoveMap& movemap ) const {
	fragment::FragSetCOP frags( NULL );
	for ( TopologyClaimers::const_iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		fragment::FragSetCOP new_frags = (*top)->loop_frags( movemap );
		runtime_assert( !(new_frags && frags ) );
		if ( new_frags ) frags = new_frags;
	}
	runtime_assert( frags );
	runtime_assert( !frags->empty() );
	return frags;
}

void TopologyBroker::add_constraints( core::pose::Pose &pose ) {
	pose.constraint_set( NULL );
	for ( TopologyClaimers::const_iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		(*top)->add_constraints( pose );
	}
	if ( tr.Trace.visible() ) {
		tr.Trace << "all constraints\n " << std::endl;
		pose.constraint_set()->show_definition( tr.Trace, pose );
		tr.Trace << std::endl;
	}
}

moves::MoverOP TopologyBroker::mover(core::pose::Pose const& pose, 
									 abinitio::StageID stage_id,
									 core::scoring::ScoreFunction const& scorefxn,
									 core::Real progress ) const {
	
	moves::RandomMoverOP random_mover = new moves::RandomMover;
	for ( TopologyClaimers::const_iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		(*top)->add_mover( *random_mover, pose, stage_id, scorefxn, progress );
	}
	
	//should this be an exception? --- it seems pretty pathological
	if ( !random_mover->size() ) {
		tr.Warning << "[ WARNING ] no mover returned in stage " << stage_id
				   << " progress " << progress << std::endl;
	}
	
	runtime_assert( random_mover->size() ); //seg-fault down the line otherwise
	return random_mover;
}

void TopologyBroker::apply_filter( core::pose::Pose const& pose,
								  abinitio::StageID stage_id,
								  core::Real progress ) const {
	tr.Debug << "apply filter: \n";
	for ( TopologyClaimers::const_iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		std::ostringstream report;
		if ( !(*top)->passes_filter( pose, stage_id, progress, report ) ) {
			tr.Debug.flush();
			throw EXCN_FilterFailed( report.str() );
		}
		if ( report.str().size() ) tr.Debug << "CLAIMER " << (*top)->type() << ":" << (*top)->label() << " FILTER REPORT: \n"
																				<< report.str();
	}
	tr.Debug.flush();
}

void TopologyBroker::accept_claims( DofClaims& claims ) {
	for ( DofClaims::iterator claim=claims.begin();	claim != claims.end(); ++claim ) {
		(*claim)->owner()->claim_accepted( *claim );
	}
}

bool TopologyBroker::broking( DofClaims const& all_claims, DofClaims& pre_accepted ) {
	//	DofClaims pre_accepted;
	tr.Debug << "broking claims..." << std::endl;
	bool fatal( false );
	for ( DofClaims::const_iterator claim = all_claims.begin();
				claim != all_claims.end() && !fatal; ++claim ) 	{
		bool allow( true );
		for ( TopologyClaimers::iterator top = claimers_.begin();
					top != claimers_.end() && allow; ++top ) {
			allow = (*top)->allow_claim( **claim );
		}
		if ( !allow ) {
			tr.Trace << "declined: " << **claim << std::endl;
			fatal = !(*claim)->owner()->accept_declined_claim( **claim );
		} else {
			tr.Trace << "accepted: " << **claim << std::endl;
			pre_accepted.push_back( *claim );
		}
	}
	return !fatal;
}

bool TopologyBroker::has_sequence_claimer() {
	DofClaims seq_claims;
	generate_sequence_claims( seq_claims );

	core::Size nres( 0 );
	for ( DofClaims::iterator claim = seq_claims.begin();	claim != seq_claims.end(); ++claim ) {
		nres += (*claim)->pos( 2 );
	}
	return nres > 0;
}

void TopologyBroker::build_fold_tree( DofClaims& claims, Size nres ) {
	DofClaims exclusive_jumps;
	DofClaims negotiable_jumps;
	DofClaims must_cut;
	DofClaims cut_biases;
	std::vector< int > obligate_cut_points; //yes FoldTree.cc uses still these classes

	Size root( 0 );
	bool excl_root_set( false );

	tr.Debug << "build fold tree ... " << std::endl;
	for ( DofClaims::iterator claim=claims.begin();	claim != claims.end(); ++claim ) {
		if ( (*claim)->type() == DofClaim::JUMP ) {
			if ( (*claim)->exclusive() ) {
				exclusive_jumps.push_back( *claim );
			} else {
				negotiable_jumps.push_back( *claim );
			}
		} else if ( (*claim)->type() == DofClaim::CUT ) {
			must_cut.push_back( *claim );
			tr.Trace << "obligate cut-point requested at " << (*claim)->pos( 1 ) << std::endl;
			obligate_cut_points.push_back( (int) (*claim)->pos( 1 ) );
		} else if ( (*claim)->type() == DofClaim::ROOT ) {
			//we allow only a single non-exclusive setting of root --- this can be overwritten by a single exclusive clain
			if ( ( root && !excl_root_set && (*claim)->exclusive() ) || !root || root == (*claim)->pos( 1 ) ) {
				root = (*claim)->pos( 1 );
				excl_root_set = (*claim)->exclusive();
			} else {
				throw( kinematics::EXCN_InvalidFoldTree( "Shouldn't have two exclusive roots --- ask oliver, throw an exception ?", *fold_tree_ ) );
			}
		}
	}

	if ( !root ) root = 1;

	utility::vector1< core::Real > cut_bias( nres, 1.0 );
	for ( TopologyClaimers::iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		(*top)->manipulate_cut_bias( cut_bias );
	}
	ObjexxFCL::FArray1D_float cut_bias_farray( nres );
	for ( Size i=1; i<=nres; i++ ) cut_bias_farray( i ) = cut_bias[ i ];

	ObjexxFCL::FArray2D_int jumps( 2, exclusive_jumps.size() + negotiable_jumps.size() );
	ObjexxFCL::FArray2D< std::string > jump_atoms( 2, exclusive_jumps.size() + negotiable_jumps.size() );
	ObjexxFCL::FArray2D_int after_loops_jumps( 2, exclusive_jumps.size() + negotiable_jumps.size() );
	ObjexxFCL::FArray2D< std::string > after_loops_jump_atoms( 2, exclusive_jumps.size() + negotiable_jumps.size() );

	Size nexcl( 0 );
	Size n_non_removed( 0 );
	// make list of all exclusive jumps from 1 .. nexecl
	for ( DofClaims::iterator jump=exclusive_jumps.begin();	jump != exclusive_jumps.end(); ++jump ) {
		runtime_assert( (*jump)->size() == 2 );
		JumpClaimOP jump_ptr( dynamic_cast< JumpClaim* >( jump->get() ) );
		runtime_assert( jump_ptr );
		++nexcl;
		jumps( 1, nexcl) = (*jump)->pos( 1 );
		jumps( 2, nexcl) = (*jump)->pos( 2 );
		jump_atoms( 1, nexcl ) = jump_ptr->jump_atom( 1 );
		jump_atoms( 2, nexcl ) = jump_ptr->jump_atom( 2 );
		if ( !jump_ptr->remove() ) {
			++n_non_removed;
			after_loops_jumps( 1, n_non_removed ) = jump_ptr->pos( 1 );
			after_loops_jumps( 2, n_non_removed ) = jump_ptr->pos( 2 );
			after_loops_jump_atoms( 1, n_non_removed ) = jump_ptr->jump_atom( 1 );
			after_loops_jump_atoms( 2, n_non_removed ) = jump_ptr->jump_atom( 2 );
		}
	}
	/// add all negotiable jumps to list ... 1..nexcl+nnegot
	Size nnegot( 0 );
	for ( DofClaims::iterator jump=negotiable_jumps.begin();	jump != negotiable_jumps.end(); ++jump ) {
		runtime_assert( (*jump)->size() == 2 );
		JumpClaimOP jump_ptr( dynamic_cast< JumpClaim* >( jump->get() ) );
		runtime_assert( jump_ptr );
		++nnegot;
		jumps( 1, nexcl+nnegot ) = (*jump)->pos( 1 );
		jumps( 2, nexcl+nnegot ) = (*jump)->pos( 2 );
		jump_atoms( 1, nexcl+nnegot ) = jump_ptr->jump_atom( 1 );
		jump_atoms( 2, nexcl+nnegot ) = jump_ptr->jump_atom( 2 );
	}
	tr.Trace << "negotiable jumps: " << negotiable_jumps << std::endl;
	tr.Trace << "exclusive jumps: " << exclusive_jumps << std::endl;
	bool bValidTree = false;
	bool try_again = true;

	//there might be cutpoints in the starting points these are added for sampling...
	//but not for final fold-tree
	std::vector< int > obligate_sampling_cut_points( obligate_cut_points );
	std::copy( start_pose_cuts_.begin(), start_pose_cuts_.end(), std::back_inserter( obligate_sampling_cut_points ) );

	while ( try_again && !bValidTree ) {
		fold_tree_ = new kinematics::FoldTree;
		Size attempts( 10 );
		while ( !bValidTree && attempts-- > 0 )  {
			bValidTree = fold_tree_->random_tree_from_jump_points( nres, nexcl+nnegot, jumps, obligate_sampling_cut_points, cut_bias_farray, root, true );
		}
		try_again = false;// for now. later we can think of ways to improve by switching negotiable_jumps on and off
	}

	if ( !bValidTree ) {
		std::ostringstream msg;
		for ( Size i = 1; i<=nexcl+nnegot; ++i ) {
		  msg << jumps(1, i ) << " " << jumps(2, i ) << std::endl;
		}
		throw( kinematics::EXCN_InvalidFoldTree( "TopologyBroker failed to make a fold-tree in 10 attempts\n"+msg.str(), *fold_tree_ ) );
	}

	for ( Size i = 1; i <= nexcl+nnegot; ++i ) {
		fold_tree_->set_jump_atoms( i, jumps( 1, i), jump_atoms( 1, i ), jumps( 2, i), jump_atoms( 2, i ) );
	}
	fold_tree_->put_jump_stubs_intra_residue();

	//make final tree: cutbias will be 0 everywhere 1 where we had cuts before ... makes sense?
	for ( Size i=1; i<=nres; i++ ) {
		cut_bias_farray( i ) = 0.0;
		if ( fold_tree_->is_cutpoint( i ) ) cut_bias_farray( i ) = 1.0;
	}
	final_fold_tree_ = new kinematics::FoldTree;
	bool bValidFinalTree =
		final_fold_tree_->random_tree_from_jump_points( nres, n_non_removed, after_loops_jumps, obligate_cut_points, cut_bias_farray, root );
	if ( !bValidFinalTree ) {
		throw( kinematics::EXCN_InvalidFoldTree( "TopologyBroker failed to make a final_fold-tree in 1 attempts ", *final_fold_tree_ ) );
	}
	for ( Size i = 1; i <= n_non_removed; ++i ) {
		final_fold_tree_->set_jump_atoms( i, after_loops_jumps( 1,i), after_loops_jump_atoms(1,i), after_loops_jumps( 2,i ), after_loops_jump_atoms(2,i) );
	}
	final_fold_tree_->put_jump_stubs_intra_residue();
}

void TopologyBroker::initialize_sequence( DofClaims& claims, core::pose::Pose& new_pose ) {
	DofClaims failures;
	sequence_claims_.clear();
	std::set< std::string > labels; //for now sequence-claims define continuous patches. they are moveable if pos(1)==0.

	//first round add those with pos(1) != 0
	for ( DofClaims::iterator claim=claims.begin();	claim != claims.end(); ++claim ) {
		if ( (*claim)->type() == DofClaim::SEQUENCE ) {
			SequenceClaimOP seq_claim = dynamic_cast< SequenceClaim* >( claim->get() );
			if ( seq_claim->pos( 1 ) ) {
				sequence_claims_.push_back( seq_claim );
			}
		}
	}

	//first round add those with pos(1) == 0
	for ( DofClaims::iterator claim=claims.begin();	claim != claims.end(); ++claim ) {
		if ( (*claim)->type() == DofClaim::SEQUENCE ) {
			SequenceClaimOP seq_claim = dynamic_cast< SequenceClaim* >( claim->get() );
			if ( seq_claim->pos( 1 ) == 0 ) {
				sequence_claims_.push_back( seq_claim );
			}
		}
	}

	core::Size nres( 0 );
	for ( SequenceClaims::iterator seq_claim_it =sequence_claims_.begin(); seq_claim_it != sequence_claims_.end(); ++seq_claim_it ) {
		(*seq_claim_it)->set_offset( nres + 1 );
		nres += (*seq_claim_it)->pos( 2 );
		if ( (*seq_claim_it)->pos( 2 ) ) {  //it actually provides sequence and is not a query
			//check label is unique
			tr.Debug << "found SequenceClaim labelled: " << (*seq_claim_it)->label() << std::endl;
			if ( labels.find( (*seq_claim_it)->label() ) == labels.end() ) {
				labels.insert( (*seq_claim_it)->label() );
			} else {
				throw EXCN_Input( "found duplicate sequence label "+(*seq_claim_it)->label() );
			}
		}
	}

	for ( SequenceClaims::iterator claim = sequence_claims_.begin();	claim != sequence_claims_.end(); ++claim ) {
		(*claim)->owner()->initialize_residues( new_pose, *claim, failures );
	}

	runtime_assert( failures.size() == 0 );
}

///@brief get the sequence claim that is consistent with the label,
/// throws EXCN_Unknown_SequenceLabel if not found
SequenceClaim& TopologyBroker::resolve_sequence_label( std::string const& label ) const {
	SequenceClaimOP found(NULL);
	for ( SequenceClaims::const_iterator claim = sequence_claims_.begin();	claim != sequence_claims_.end(); ++claim ) {
		if ( (*claim)->label() == label ) {
			runtime_assert( !found ); //don't allow duplicate labels -- input is checked when this list is made
			found = *claim;
		}
	}
	if ( !found ) throw EXCN_Unknown( "requested SequenceLabel " + label + " not found " );
	return *found;
}

core::Size TopologyBroker::resolve_residue( std::string const& chain_label, core::Size pos ) const {
	return resolve_sequence_label( chain_label ).offset() + pos - 1;
}

void TopologyBroker::initialize_dofs( DofClaims& claims, core::pose::Pose& pose ) {
	DofClaims bb_claims( pose.total_residue(), NULL ); //one claim per position
	DofClaims jumps( pose.num_jump(), NULL );
	for ( DofClaims::iterator claim = claims.begin();	claim != claims.end(); ++claim ) {
		if ( (*claim)->type() == DofClaim::BB ) {
			//BBClaimOP bb_claim = dynamic_cast< BBClaim* >( claim->get() );
			DofClaimOP bb_claim = *claim;
			Size pos = bb_claim->pos( 1 );
			if ( pos > pose.total_residue() ) throw EXCN_BadInput( "attempt to initialize dof "+string_of( pos ) +
				" in "+string_of( pose.total_residue() ) + " pose. ... fragments inconsistent with fasta ? " );
			if ( !bb_claims[ pos ] || bb_claims[ pos ]->right() < bb_claim->right() ) {
				bb_claims[ pos ] = bb_claim;
			} else runtime_assert( bb_claim->right() < DofClaim::EXCLUSIVE );
		} else if ( (*claim)->type() == DofClaim::JUMP ) {
			//JumpClaimOP jump = dynamic_cast< JumpClaim* >( claim->get() );
			DofClaimOP jump = *claim;
			Size jump_nr ( pose.fold_tree().jump_nr( jump->pos( 1 ), jump->pos( 2 ) ) );
			runtime_assert( jump_nr ); //XOR would be even better
			DofClaimOP already ( jumps[ jump_nr ] ); //one of them is 0 -- neutral to addition
			if ( !already || already->right() < jump->right() ) {
				jumps[ jump_nr ] = jump;
			} else runtime_assert( jump->right() < DofClaim::EXCLUSIVE );
		}
	}

	if ( tr.Trace.visible() ) tr.Trace << "init-bb-dofs\n " <<  bb_claims << " init-jump-dofs\n" << jumps << std::endl;

	//check for un-initialized dofs and throw exception with dof_msg if not fully covered
	std::ostringstream dof_msg;
	bool bad( false );
	Size pos( 1 );
	for ( DofClaims::const_iterator it = bb_claims.begin(); it != bb_claims.end(); ++it, ++pos ) {
		if ( !*it ) {
			//throw exception later, but first accumulate all errors in dof_msg
			bad = true;
			dof_msg << "BBTorsion at pos " << pos << "unitialized...unclaimed" << std::endl;
		}
	}

	//check for un-initialized dofs
	for ( DofClaims::const_iterator it = jumps.begin(); it != jumps.end(); ++it ) {
		if ( !*it ) {
			bad = true;
			dof_msg << "Jump unitialized... " << *it << std::endl;
			bool bUnitializedJump = true;
			runtime_assert( !bUnitializedJump );
		}
	}
	if ( bad ) throw( EXCN_Input( dof_msg.str() ) );

	DofClaims cumulated;
	std::copy( bb_claims.begin(), bb_claims.end(), back_inserter( cumulated ) );
	std::copy( jumps.begin(), jumps.end(), back_inserter( cumulated ) );

	DofClaims failures;
	for ( TopologyClaimers::iterator top = claimers_.begin();
				top != claimers_.end(); ++top ) {
		(*top)->initialize_dofs( pose, cumulated, failures );
	}

	if ( failures.size() ) {
		std::ostringstream dof_msg;
		dof_msg << "failed to initialize dofs for these claims: .... " << failures << std::endl;
		throw EXCN_Unknown( dof_msg.str() );
	}
	runtime_assert( failures.size() == 0 ); //should have thrown exception before -- Exception
}

void TopologyBroker::initialize_cuts( DofClaims& claims, core::pose::Pose& pose ) {
	//cuts will contain NULL for automatic cutpoints
	//and the fold-tree cut-point nr for those which we keep since they have been CUT-Claimed
	DofClaims cuts( pose.num_jump(), NULL );
	for ( DofClaims::iterator claim = claims.begin();	claim != claims.end(); ++claim ) {
		if ( (*claim)->type() == DofClaim::CUT ) {
			//DofClaimOP cut = dynamic_cast< JumpClaim* >( claim->get() );
			DofClaimOP cut = *claim;
			Size cut_nr ( pose.fold_tree().cutpoint_map( cut->pos( 1 ) ) );
			// we allow cuts without claim -- random cuts due to jumping ...
			// but if there was a CUT claim, there should be a cut
			runtime_assert( cut_nr );
			cuts[ cut_nr  ] = cut;
			//NOTE in principle we could support mutiple CUT claims at the same position... if they all just want a cut there.
			// maybe need to allow keeping of multiple DofClaims per cutpoint in the list... if the list is actully ever needed.
			// we assume a DofClaim::CUT will never be closed.. of course we could also change that... then however, the rights will play a role.
		}
	}

	//find the unclaimed -- automatic -- cutpoints
	to_be_closed_cuts_.clear();
	for ( Size cut_nr = 1; cut_nr<=cuts.size(); ++cut_nr ) {
		if ( !cuts[ cut_nr ] ) { //automatic cut-point
			to_be_closed_cuts_.push_back( pose.fold_tree().cutpoint( cut_nr ) );
			tr.Debug << "close this cut: " << to_be_closed_cuts_.back() << std::endl;
			// for now we assume that these always should get closed --- add cutpoint variants
			// one could also have Jumpers issue floating CUT claims (without position number)
			// then their pos will be assigned actual cutpoints at this stage

			// could also make a new TopologyClaimer --- CutCloser which will handle these extra cutpoints.
			// but difficult to remove the CutCloser from the list of TopologyClaimers before the next decoy

			//in principle we want to maintain a list of chainbreaks that are to be closed
			//Question: are CUT claims always coding for Cutpoints that shall not be closed ?
			// so far I see it that way...
		}
	}
}

void TopologyBroker::apply( core::pose::Pose& pose ) {
	DofClaims pre_accepted;
	bool ok( true );

	for ( TopologyClaimers::iterator top = claimers_.begin();
				top != claimers_.end(); ++top ) {
		if ( bUseJobPose_ ) {
			(*top)->new_decoy( pose );
		} else {
			(*top)->new_decoy();
		}
	}

	start_pose_cuts_.clear();
	if ( bUseJobPose_ ) {
		for ( Size i(1), cutpoints(pose.fold_tree().num_cutpoint()); i<=cutpoints; i++ ) {
			start_pose_cuts_.push_back( pose.fold_tree().cutpoint( i ) );
		}
	}

	pose.clear(); //yay!
	
	if ( pose::symmetry::is_symmetric( pose ) ) {
		pose::symmetry::make_asymmetric_pose( pose );
	}
	
	tr.Debug << "Initialize Sequence" << std::endl;
	DofClaims fresh_claims;
	
	sequence_claims_.clear();
	generate_sequence_claims( fresh_claims );
	if ( ok ) ok = broking( fresh_claims, pre_accepted );
	initialize_sequence( pre_accepted, pose );
	current_pose_ = new core::pose::Pose( pose );
	
	tr.Debug << "Start Round1-Broking..." << std::endl;
	DofClaims round1_claims;
	generate_round1( round1_claims );
	if ( ok ) ok = broking( round1_claims, pre_accepted );

	tr.Debug << "Start FinalRound-Broking..." << std::endl;
	DofClaims final_claims;
	generate_final_claims( final_claims );
	if ( ok ) ok = broking( final_claims, pre_accepted );

	tr.Debug << "Broking finished" << std::endl;
	//	--> now we know nres
	if ( tr.Debug.visible() )	pose.dump_pdb( "init_seq.pdb" );
	current_pose_ = new core::pose::Pose( pose );

	tr.Debug << "build fold-tree..." << std::endl;
	build_fold_tree( pre_accepted, pose.total_residue() );
	accept_claims( pre_accepted );
	tr.Debug << *fold_tree_ << std::endl;
	pose.fold_tree( *fold_tree_ );

	tr.Debug << "set cuts..." << std::endl;
	initialize_cuts( pre_accepted, pose );

	tr.Debug << "initialize dofs..." << std::endl;
	initialize_dofs( pre_accepted, pose );
	if ( tr.Debug.visible() )	pose.dump_pdb( "init_dofs.pdb" );
	current_pose_ = new core::pose::Pose( pose );

	//we will need this one in switch_to_fullatom
	if ( !repack_scorefxn_ ) repack_scorefxn_ = core::scoring::getScoreFunction();

	// initialize secondary structure from DSSP.
	core::scoring::dssp::Dssp dssp_obj( pose );
	dssp_obj.insert_ss_into_pose( pose );

	// Fix disulfides if a file is given
	if ( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].user() ) {
		io::raw_data::DisulfideFile ds_file( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ]() );
		utility::vector1< std::pair<Size,Size> > disulfides;
		ds_file.disulfides(disulfides, pose);
		pose.conformation().fix_disulfides( disulfides );
	}

	add_chainbreak_variants( pose, 0 /*ignored*/, NULL /*no sequence separation switch*/ );

	//add constraints
	add_constraints( pose );
}

bool TopologyBroker::has_chainbreaks_to_close() const {
	return to_be_closed_cuts_.size();
}

void TopologyBroker::add_chainbreak_variants(
	 pose::Pose &pose,
	 Size max_dist,
	 core::kinematics::ShortestPathInFoldTreeCOP sp
) const {
	pose::Pose init_pose = pose;
	for ( utility::vector1< Size >::const_iterator it = to_be_closed_cuts_.begin();
				it != to_be_closed_cuts_.end(); ++it ) {
		tr.Debug << "consider cut between res " << *it << " and " << *it+1;
		if ( sp ) tr.Debug << " distance is " << sp->dist( *it, *it+1 ) << " of max " << sp->max_dist();
		tr.Debug << " (" << max_dist << ")"<< std::endl;
		if ( sp && max_dist && sp->dist( *it, *it+1 ) > max_dist ) continue;
		if ( !pose.fold_tree().is_cutpoint( *it ) ){
			continue; //maybe we are in full-atom mode, or have some of chainbreaks already closed... ?
		}
		tr.Debug << "add chainbreak variant to residues " << *it << " and " << *it+1 << std::endl;
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, *it );
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, *it+1 );
	}
	pose.constraint_set( init_pose.constraint_set()->remapped_clone( init_pose, pose ) );
}

bool TopologyBroker::check_chainbreak_variants(
	 pose::Pose &pose
) const {
	bool success ( true );
	for ( utility::vector1< Size >::const_iterator it = to_be_closed_cuts_.begin();
				it != to_be_closed_cuts_.end(); ++it ) {
		tr.Debug << "consider cut between res " << *it << " and " << *it+1 << std::endl;
		if ( !pose.fold_tree().is_cutpoint( *it ) ){
			throw( kinematics::EXCN_InvalidFoldTree( "Foldtree missmatch", pose.fold_tree() ) );
		}
		if ( pose.residue( *it ).has_variant_type( chemical::CUTPOINT_LOWER )
				&& pose.residue( *it+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) {
			tr.Debug << "found chainbreak variant at residues " << *it << " and " << *it+1 << std::endl;
		} else {
			tr.Warning << "[WARNING] no chainbreak variant found at residues " << *it << " and " << *it+1 << std::endl;
			tr.Warning << jd2::current_output_name() << std::endl;
			tr.Warning << pose.fold_tree() << std::endl;
			tr.Warning << pose.annotated_sequence() << std::endl;
			success = false;
		}
	}
	return success;
}

///@brief if some claimer wants to influence the movemap for relax he can do it here:
void TopologyBroker::adjust_relax_movemap( core::kinematics::MoveMap& mm) const {
	for ( TopologyClaimers::const_iterator top = claimers_.begin();
				top != claimers_.end(); ++top ) {
		(*top)->adjust_relax_movemap( mm );
	}

}

void TopologyBroker::switch_to_fullatom( core::pose::Pose& pose ) {
	pose::Pose init_pose = pose;
	std::string sequence_old( pose.annotated_sequence() );
	tr.Debug << "switch_to_fullatom... " << std::endl;

	if ( !pose.is_fullatom() ) {
    core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
	}

	tr.Debug << "switched to fullatom... " << std::endl;

	utility::vector1< bool > needToRepack( pose.total_residue(), true );
	for ( TopologyClaimers::const_iterator top = claimers_.begin();
				top != claimers_.end(); ++top ) {
		(*top)->switch_to_fullatom( pose, needToRepack );
	} //this might copy residue-sidechains ... need to have 'disulf' stuff afterwards.
	tr.Debug << "broker switch to fullatom... " << std::endl;

	if ( basic::options::option[ basic::options::OptionKeys::in::detect_disulf ]() ) {
		//disulfide stuff .... do already on the centroid level ? --- that will triger the rebuild disulfides in swtich_to_residye_type_set
		pose.conformation().detect_disulfides();
	}

	// Fix disulfides if a file is given  --- in a claimer? ... a util.cc function() ?
	if ( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].user() ) {
		io::raw_data::DisulfideFile ds_file( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ]() );
		utility::vector1< std::pair<Size,Size> > disulfides;
		ds_file.disulfides(disulfides, pose);
		pose.conformation().fix_disulfides( disulfides );
	}

	pose.conformation().detect_bonds();//apl fix this !

	tr.Debug << "detect bonds... " << std::endl;

	//sanity check
	std::string sequence_new( pose.annotated_sequence() );
	if ( sequence_old != sequence_new ) {
		tr.Warning << "[WARNING] switch_to_fullatom changed sequence/variants:\n " <<
			" before " << sequence_old << "\n after  " << sequence_new << std::endl;
	}

    // repack loop + missing-density residues
    core::pack::task::PackerTaskOP taskstd = core::pack::task::TaskFactory::create_packer_task( pose );
    taskstd->restrict_to_repacking();
    taskstd->or_include_current(true);
	taskstd->restrict_to_residues( needToRepack );

	add_constraints( pose );

	if ( tr.Debug.visible() ) pose.dump_pdb( "before_repack.pdb");
	if ( tr.Debug.visible() ) pose.constraint_set()->show_numbers( tr.Debug );

	protocols::simple_moves::PackRotamersMover pack1( repack_scorefxn_ , taskstd );
    pack1.apply( pose );

    // quick SC minimization
    core::optimization::AtomTreeMinimizer mzr;
    core::optimization::MinimizerOptions options( "dfpmin_armijo_nonmonotone", 1e-5, true, false );
    core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap();
    mm->set_bb( false );
    mm->set_chi( true );

	if ( pose::symmetry::is_symmetric( pose ) ) {
		core::pose::symmetry::make_symmetric_movemap( pose, *mm );
		core::optimization::symmetry::SymAtomTreeMinimizer smzr;
		smzr.run( pose, *mm, *repack_scorefxn_, options );
	} else {
		mzr.run( pose, *mm, *repack_scorefxn_, options );
	}

	tr.Debug << "minimized.. " << std::endl;
}

}
}
