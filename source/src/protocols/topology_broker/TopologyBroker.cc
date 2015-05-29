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
/// @details responsibilities:
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
#include <protocols/topology_broker/SymmetryClaimer.hh>
#include <protocols/topology_broker/SequenceNumberResolver.hh>
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/LegacyRootClaim.hh>
#include <protocols/topology_broker/claims/CutClaim.hh>
#include <protocols/topology_broker/claims/SequenceClaim.hh>
#include <protocols/topology_broker/claims/JumpClaim.hh>
#include <protocols/topology_broker/claims/BBClaim.hh>
#include <protocols/topology_broker/claims/SymmetryClaim.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <protocols/jd2/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
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
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/annotated_sequence.hh>
#include <protocols/moves/MoverContainer.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// for symmetry
#include <core/pose/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <numeric/xyzTriple.hh>

// Utility headers
#include <basic/Tracer.hh>

// C++ headers
#ifdef WIN32
#include <algorithm>
#include <iterator>
#endif

#include <sstream>
#include <vector>
#include <stdexcept>
#include <utility>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace ObjexxFCL { } using namespace ObjexxFCL;

static thread_local basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace utility::excn;
using namespace core;

TopologyBroker::TopologyBroker() :
	nres_( 0 ),
	fold_tree_( /* NULL */ ),
	final_fold_tree_( /* NULL */ ),
	repack_scorefxn_( /* NULL */ ),
	bUseJobPose_( false ),
	use_fold_tree_from_claimer_(false),
	current_pose_( /* NULL */ )
{
	sequence_number_resolver_ = SequenceNumberResolverOP( new SequenceNumberResolver() );
}

TopologyBroker::~TopologyBroker() {}

TopologyBroker::TopologyBroker( const TopologyBroker& tp ) :
	ReferenceCount(),
	utility::pointer::enable_shared_from_this< TopologyBroker >(),
	claimers_( tp.claimers_ )
{
	current_claims_ = tp.current_claims_;
	nres_ = tp.nres_;
	fold_tree_ = tp.fold_tree_;
	final_fold_tree_ = tp.final_fold_tree_;
	repack_scorefxn_ = tp.repack_scorefxn_;
	to_be_closed_cuts_ = tp.to_be_closed_cuts_;
	start_pose_cuts_ = tp.start_pose_cuts_;
	use_fold_tree_from_claimer_ = tp.use_fold_tree_from_claimer_;
	current_pose_ = tp.current_pose_;
	sequence_number_resolver_ = SequenceNumberResolverOP( new SequenceNumberResolver( *tp.sequence_number_resolver_ ) );
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
		use_fold_tree_from_claimer_ = src.use_fold_tree_from_claimer_;
		current_pose_ = src.current_pose_;
	}
	return *this;
}

void TopologyBroker::add( TopologyClaimerOP cl ) {
	claimers_.push_back( cl );
	cl->set_broker( get_self_weak_ptr() );
}

void TopologyBroker::generate_sequence_claims( claims::DofClaims& all_claims ) {
	for ( TopologyClaimers::iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		(*top)->generate_sequence_claims( all_claims );
	}
	tr.Trace << "All sequence claims: \n" << all_claims << std::endl;
}

void TopologyBroker::generate_symmetry_claims( claims::SymmetryClaims& all_claims ) {
	for ( TopologyClaimers::iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		(*top)->generate_symmetry_claims( all_claims );
	}
}

SymmetryClaimerOP TopologyBroker::resolve_symmetry_claims( claims::SymmetryClaims& symm_claims ){
	if ( symm_claims.size() > 1 ){
		tr.Error << "Multiple Symmetry Claims found; Topology Broker does not support multiple symmetries."
                 << std::endl;
		utility_exit();
	}
	else if ( symm_claims.size() == 1 ){
        claims::SymmetryClaimOP symmclaim = symm_claims.at(1);
        SymmetryClaimerOP symmclaimer = utility::pointer::dynamic_pointer_cast< SymmetryClaimer > ( symmclaim->owner().lock() );

		return symmclaimer;
	}

    return NULL;
}

void TopologyBroker::relay_message( ClaimerMessage& msg ) const {
	for ( TopologyClaimers::const_iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		if ( msg.matches( (*top)->label() ) ) {
			(*top)->receive_message( msg );
		}
	}
}

void TopologyBroker::generate_round1( claims::DofClaims& all_claims ) {
	for ( TopologyClaimers::iterator top = claimers_.begin();
					top != claimers_.end(); ++top ) {
		tr.Trace << "generate claim for " << (*top)->type() << " with label "<< (*top)->label() << std::endl;
		(*top)->generate_claims( all_claims );
	}
	tr.Trace << "All round1 claims: \n" << all_claims << std::endl;
}

void TopologyBroker::generate_final_claims( claims::DofClaims& all_claims ) {
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
	runtime_assert( frags != 0 );
	runtime_assert( !frags->empty() );
	return frags;
}

void TopologyBroker::add_constraints( core::pose::Pose &pose ) const {
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

	//tr.Debug << "stage:  " << stage_id << " Progress:  " << progress << std::endl;
	moves::RandomMoverOP random_mover( new moves::RandomMover );
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

//	tr.Debug << "number of movers in random_mover:  " << random_mover->nr_moves() << std::endl;
//	if(tr.Debug.visible())
//	{
//		for ( Size ii = 0; ii < random_mover->size(); ++ii )
//		{
//			tr.Debug << "the " << ii <<"th mover in random_mover is:  " << random_mover->get_mover(ii) << std::endl;
//		}
//		for ( Size i=0; i<random_mover->weights().size(); ++i ) {
//			tr.Debug << "weight at " << i << " is " << random_mover->weights().at(i) << std::endl;
//		}
//	}


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

void TopologyBroker::accept_claims( claims::DofClaims& claims ) {
	for ( claims::DofClaims::iterator claim=claims.begin();	claim != claims.end(); ++claim ) {
		//tr.Debug << "accept_claims: " << (*claim)->owner()->type() << std::endl;
		//if(tr.Trace.visible())
		//{
		//	(*claim)->show(tr.Trace);
		//	tr.Trace << std::endl;
		//}
		(*claim)->owner().lock()->claim_accepted( *claim );
	}
}


// broking is very simple:
// for each claim, ask all claimers if they have a problem with the given game,
// if declined, tell the owner of claim --> accept_declined_claim
// if ok, store in list 'pre_accepted'
bool TopologyBroker::broking( claims::DofClaims const& all_claims, claims::DofClaims& pre_accepted ) {
	//	claims::DofClaims pre_accepted;
	tr.Debug << "broking claims..." << std::endl;
	bool fatal( false );
	for ( claims::DofClaims::const_iterator claim = all_claims.begin();
				claim != all_claims.end() && !fatal; ++claim ) 	{
		bool allow( true );
		for ( TopologyClaimers::iterator top = claimers_.begin();
					top != claimers_.end() && allow; ++top ) {
			allow = (*top)->allow_claim( **claim );
		}
		if ( !allow ) {
			tr.Trace << "declined: " << **claim << std::endl;
			fatal = !(*claim)->owner().lock()->accept_declined_claim( **claim );
		} else {
			tr.Trace << "accepted: " << **claim << std::endl;
			pre_accepted.push_back( *claim );
		}
	}
	return !fatal;
}

bool TopologyBroker::has_sequence_claimer() {
	claims::DofClaims seq_claims;
	generate_sequence_claims( seq_claims );

	core::Size nres( 0 );
	for ( claims::DofClaims::iterator claim = seq_claims.begin();	claim != seq_claims.end(); ++claim ) {
		claims::SequenceClaimOP seq_ptr( utility::pointer::dynamic_pointer_cast< claims::SequenceClaim >( *claim ) );
		runtime_assert( seq_ptr != 0 );
		nres += seq_ptr->length();
	}
	return nres > 0;
}

void TopologyBroker::build_fold_tree( claims::DofClaims& claims, Size nres ) {
	claims::JumpClaims exclusive_jumps;
	claims::JumpClaims negotiable_jumps;
	utility::vector1< std::pair< claims::CutClaimOP , core::Size > > must_cut;
	claims::CutClaims cut_biases;

	Size root( 0 );
	bool excl_root_set( false );

	//building the fold tree from DofClaims and jumps
	tr.Debug << "build fold tree ... " << std::endl;
	for ( claims::DofClaims::iterator claim=claims.begin();	claim != claims.end(); ++claim ) {
		claims::JumpClaimOP jump_ptr( utility::pointer::dynamic_pointer_cast< claims::JumpClaim >( *claim ) );
		claims::CutClaimOP cut_ptr( utility::pointer::dynamic_pointer_cast< claims::CutClaim >( *claim ) );
		claims::LegacyRootClaimOP root_ptr( utility::pointer::dynamic_pointer_cast< claims::LegacyRootClaim >( *claim ) );

		if ( jump_ptr ) { // JumpClaim ------------------------------------
			if ( jump_ptr->exclusive() ) {
				exclusive_jumps.push_back( jump_ptr );
			} else {
				negotiable_jumps.push_back( jump_ptr );
			}
		} else if ( cut_ptr ) { // CutClaim -------------------------------
			Size global_position = sequence_number_resolver_->find_global_pose_number( cut_ptr->get_position() );
			if( global_position >= nres || global_position < 1){
				tr.Debug << "Cut claim requesting cut at (global) position " << global_position << " by claimer '"
								 << cut_ptr->owner().lock()->label() << "' being ignored/erased because it is outside the valid sequence region "
								 << "[1, " << nres << "]. This is expected behavior from one (and only one) sequence claimer." << std::endl;

				claims.erase( claim );
				//Modifying removing an element from the list apparently requires us to go back one in the iterator.
				--claim;
			} else {
				must_cut.push_back( std::make_pair(cut_ptr, global_position) );
				tr.Info << "Obligate cut point requested at (" << cut_ptr->get_position().first << ","
								 << cut_ptr->get_position().second << ")" << std::endl;
			}
		} else if ( root_ptr ) { // LegacyRootClaim -----------------------------
			//we allow only a single non-exclusive setting of root --- this can be overwritten by a single exclusive clain
			if ( ( root && !excl_root_set && root_ptr->exclusive() ) || !root /*|| root == root_ptr->get_position()*/ ) {
				root = sequence_number_resolver_->find_global_pose_number(root_ptr->local_position()); // root_ptr->get_position();
				excl_root_set = root_ptr->exclusive();
			} else {
				throw( kinematics::EXCN_InvalidFoldTree( "Shouldn't have two exclusive roots --- ask oliver, throw an exception ?",
																								 *fold_tree_ ) );
			}
		}
	}

	if ( !root ) root = 1;

	utility::vector1< core::Real > cut_bias( nres, 1.0 );
	for ( TopologyClaimers::iterator top = claimers_.begin(); top != claimers_.end(); ++top ) {
		(*top)->manipulate_cut_bias( cut_bias );
	}
	ObjexxFCL::FArray1D_float cut_bias_farray( nres );
	for ( Size i=1; i<=nres; i++ ) cut_bias_farray( i ) = cut_bias[ i ];

	//Sort jumps by exclusivity, so that they are placed in the ObjexxFCL array correctly
	claims::JumpClaims sorted_jumps;
	std::copy( exclusive_jumps.begin(),  exclusive_jumps.end(), std::back_inserter( sorted_jumps ) );
	std::copy( negotiable_jumps.begin(), negotiable_jumps.end(), std::back_inserter( sorted_jumps ) );

	//Build these dumb ObjexxFCL arrays
	ObjexxFCL::FArray2D_int jumps( 2, sorted_jumps.size() );
	ObjexxFCL::FArray2D< std::string > jump_atoms( 2, sorted_jumps.size() );
	ObjexxFCL::FArray2D_int after_loops_jumps( 2, sorted_jumps.size() );
	ObjexxFCL::FArray2D< std::string > after_loops_jump_atoms( 2, sorted_jumps.size() );

	Size i = 0;
	Size n_non_removed( 0 );
	for ( claims::JumpClaims::iterator jump_it = sorted_jumps.begin();	jump_it != sorted_jumps.end(); ++jump_it ) {
		claims::JumpClaimOP jump_ptr( *jump_it );
		++i;
		jumps( 1, i) = sequence_number_resolver_->find_global_pose_number(jump_ptr->local_pos1());
		jumps( 2, i) = sequence_number_resolver_->find_global_pose_number(jump_ptr->local_pos2());
		jump_atoms( 1, i ) = jump_ptr->jump_atom( 1 );
		jump_atoms( 2, i ) = jump_ptr->jump_atom( 2 );
		if ( !jump_ptr->remove() && jump_ptr->exclusive() ) {
			++n_non_removed;
			after_loops_jumps( 1, n_non_removed ) = sequence_number_resolver_->find_global_pose_number(jump_ptr->local_pos1());
			after_loops_jumps( 2, n_non_removed ) = sequence_number_resolver_->find_global_pose_number(jump_ptr->local_pos2());
			after_loops_jump_atoms( 1, n_non_removed ) = jump_ptr->jump_atom( 1 );
			after_loops_jump_atoms( 2, n_non_removed ) = jump_ptr->jump_atom( 2 );
		}
	}

	bool bValidTree = false;
	bool try_again = true;

	//there might be cutpoints in the starting points these are added for sampling...
	//but not for final fold-tree
	std::vector< int > obligate_sampling_cut_points;
	std::copy( start_pose_cuts_.begin(), start_pose_cuts_.end(), std::back_inserter( obligate_sampling_cut_points ) );
	for( Size i=1; i <= must_cut.size(); ++i ){
		if( std::find( obligate_sampling_cut_points.begin(), obligate_sampling_cut_points.end(),
									 (int) must_cut.at(i).second ) == obligate_sampling_cut_points.end() ){
			obligate_sampling_cut_points.push_back( (int) must_cut.at(i).second );
		}
	}

	while ( try_again && !bValidTree ) {
		fold_tree_ = core::kinematics::FoldTreeOP( new kinematics::FoldTree );
		Size attempts( 10 );
		while ( !bValidTree && attempts-- > 0 )  {
			bValidTree = fold_tree_->random_tree_from_jump_points( nres, sorted_jumps.size(), jumps,
																														 obligate_sampling_cut_points, cut_bias_farray, root, true );
		}
		try_again = false;// for now. later we can think of ways to improve by switching negotiable_jumps on and off
	}

	if ( !bValidTree ) {
		std::ostringstream msg;
		for ( Size i = 1; i<= sorted_jumps.size(); ++i ) {
		  msg << jumps( 1, i ) << " " << jumps( 2, i ) << std::endl;
		}
		throw( kinematics::EXCN_InvalidFoldTree( "TopologyBroker failed to make a fold-tree in 10 attempts\n"+msg.str(),
																						 *fold_tree_ ) );
	}

	//Verify that we got all the cut points we asked for.
	for( Size i=1; i <= must_cut.size(); ++i ) {
		core::Size global_seqpos = must_cut.at(i).second;
		tr.Debug << "Checking CutClaim for global position " << global_seqpos << std::endl;
		if(!fold_tree_->is_cutpoint( (int) global_seqpos  ) ){
			std::ostringstream msg;
			msg << "The requested cutpoint at " << global_seqpos << " (claimed by Claimer with label '"
					<< must_cut.at(i).first->owner().lock()->label() << "') was not found in constructed fold tree. "
					<< "This is sometimes because more cuts are claimed than jumps." << std::endl
					<< "In this case, there are " << must_cut.size() << " CutClaims and " << sorted_jumps.size()
					<< " JumpClaims. A trivial (but scientifically invalid) solution would be to add a "
					<< "BasicJumpClaimer to the topology broker setup file (.tpb) to claim a jump between disconnected chains."
					<< std::endl;
			throw ( utility::excn::EXCN_BadInput( msg.str() ) );
		}
	}


	for ( Size i = 1; i <= sorted_jumps.size(); ++i ) {
		fold_tree_->set_jump_atoms( i, jumps( 1, i), jump_atoms( 1, i ), jumps( 2, i ), jump_atoms( 2, i ) );
	}
	fold_tree_->put_jump_stubs_intra_residue();

	//make final tree: cutbias will be 0 everywhere 1 where we had cuts before ... makes sense?
	for ( Size i=1; i<=nres; i++ ) {
		cut_bias_farray( i ) = 0.0;
		if ( fold_tree_->is_cutpoint( i ) ) cut_bias_farray( i ) = 1.0;
	}

	fold_tree_->show( tr.Info );

	std::vector< int > obligate_cut_points;
	for( Size i=1; i <= must_cut.size(); i++){
		obligate_cut_points.push_back( (int) must_cut.at(i).second );
	}
	final_fold_tree_ = core::kinematics::FoldTreeOP( new kinematics::FoldTree );
	bool bValidFinalTree =
		final_fold_tree_->random_tree_from_jump_points( nres, n_non_removed, after_loops_jumps,
																										obligate_cut_points, cut_bias_farray, root );
	if ( !bValidFinalTree ) {
		throw( kinematics::EXCN_InvalidFoldTree( "TopologyBroker failed to make a final_fold-tree in 1 attempts ",
																						 *final_fold_tree_ ) );
	}
	for ( Size i = 1; i <= n_non_removed; ++i ) {
		final_fold_tree_->set_jump_atoms( i, after_loops_jumps( 1, i ), after_loops_jump_atoms( 1, i ),
																			after_loops_jumps( 2, i ), after_loops_jump_atoms( 2, i ) );
	}
	final_fold_tree_->put_jump_stubs_intra_residue();
	//tr.Debug << "final_fold_tree_ at end of build_fold_tree():";
	//tr.Debug << *final_fold_tree_ << std::endl;
}

//build the fold tree if the claimer builds its on fold tree
void TopologyBroker::build_fold_tree_from_claimer(core::pose::Pose& pose, core::kinematics::FoldTree& fold_tree)
{
	//if(tr.Debug.visible())
	//{
	//	for(core::Size i = 1; i <= num_claimers(); i++)
	//	{
	//		tr.Debug << claimer(i)->type() << std::endl;
	//	}
	//}

	//Loop through all TopologyClaimers and get the final fold tree from it
	for( TopologyClaimers::iterator top = claimers_.begin();top != claimers_.end(); ++top )
	{
		if((*top)->claimer_builds_own_fold_tree())
		{
			tr << "WARNING!!!!  the order of your claimers matter if you build a fold tree in the claimer.  "
					"In the case of the TMHTopologySamplerClaimer, it should come AFTER the MembraneTopologyClaimer" << std::endl;
			(*top)->set_pose_from_broker(pose);
			(*top)->build_fold_tree(pose, fold_tree);
			fold_tree_ = (*top)->get_fold_tree(pose);
			if(tr.Debug.visible()){tr.Debug << *fold_tree_ << std::endl;}
		}
	}
}

//@brief This function is here to enable the std::sort call in initialize_sequence()
bool compare_sequence_claim_priority (claims::SequenceClaimOP const& a, claims::SequenceClaimOP const& b ) {
	return a->priority() > b->priority();
}

void TopologyBroker::initialize_sequence( claims::DofClaims& claims, core::pose::Pose& new_pose ) {

	// type cast all claims into SequenceClaim...
	tr.Debug << "Initializing sequence claims...";
	tr.flush();
	sequence_claims_.clear();
	for ( claims::DofClaims::iterator claim=claims.begin();	claim != claims.end(); ++claim ) {
		claims::SequenceClaimOP seq_ptr( utility::pointer::dynamic_pointer_cast< claims::SequenceClaim >( *claim ) );
		if ( seq_ptr ){
			sequence_claims_.push_back ( seq_ptr );
		}
	}

	// sort SequenceClaims by priority.
	std::sort( sequence_claims_.begin(), sequence_claims_.end(), compare_sequence_claim_priority );

	// build up sequence and store offsets
	core::Size nres( 0 );
	std::string accumulated_sequence;

	for ( claims::SequenceClaims::iterator seq_claim_it =sequence_claims_.begin(); seq_claim_it != sequence_claims_.end();
				++seq_claim_it ) {
		tr.Debug << "Accumulating sequence '" << (*seq_claim_it)->annotated_sequence() << "' from claim labeled '"
						 << (*seq_claim_it)->label() << "'" << std::endl;
		sequence_number_resolver_->register_label_offset( (*seq_claim_it)->label(), nres ); //offset = number of residues in front of claim = nres
		accumulated_sequence+=(*seq_claim_it)->annotated_sequence();
		nres += (*seq_claim_it)->length();
	}

	// make actual sequence
	core::pose::make_pose_from_sequence(
		 new_pose,
		 accumulated_sequence,
		 *( chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ))
	);

	// make extended chain
	for ( Size pos = 1; pos <= new_pose.total_residue(); pos++ ) {
		if ( ! new_pose.residue(pos).is_protein() ) continue;
		new_pose.set_phi( pos, -150 );
		new_pose.set_psi( pos, 150);
		new_pose.set_omega( pos, 180 );
	}

	tr.Debug << "Sequence initialized. Final sequence: " << accumulated_sequence << std::endl;

	//REMEMBER: make cuts at end of each sequence?
}

/// @brief get the sequence claim that is consistent with the label,
/// throws EXCN_Unknown_SequenceLabel if not found
claims::SequenceClaim& TopologyBroker::resolve_sequence_label( std::string const& label ) const {
	claims::SequenceClaimOP found(NULL);
	for ( claims::SequenceClaims::const_iterator claim = sequence_claims_.begin();	claim != sequence_claims_.end(); ++claim ) {
		if ( (*claim)->label() == label ) {
			runtime_assert( !found ); //don't allow duplicate labels -- input is checked when this list is made
			found = *claim;
		}
	}
	if ( !found ) throw EXCN_Unknown( "requested SequenceLabel " + label + " not found " );
	return *found;
}

/*core::Size TopologyBroker::resolve_residue( std::string const& chain_label, core::Size pos ) const {
	return resolve_sequence_label( chain_label ).offset() + pos - 1;
}*/

void TopologyBroker::initialize_dofs( claims::DofClaims& claims, core::pose::Pose& pose ) {
	claims::DofClaims bb_claims( pose.total_residue(), NULL ); //one claim per position
	claims::DofClaims jumps( pose.num_jump(), NULL );
	for ( claims::DofClaims::iterator claim = claims.begin();	claim != claims.end(); ++claim ) {

		claims::BBClaimOP bb_claim ( utility::pointer::dynamic_pointer_cast< claims::BBClaim >( *claim ) );

		if ( bb_claim ){
			claims::LocalPosition pos = bb_claim->local_position();
			core::Size global_position;
			try {
				global_position = sequence_number_resolver_->find_global_pose_number( pos );
			} catch ( std::out_of_range ) {
				std::ostringstream msg;
				msg << "BBClaim at (" << pos.first << ", " << pos.second << ") claimed by " << bb_claim->owner().lock()->type()
						<< " using label '" << bb_claim->owner().lock()->label() << "' is has an invalid label. The global sequence position"
						<< " cannot be globally resolved." << std::endl;
				throw EXCN_BadInput( msg.str() );
			}

			if ( global_position > pose.total_residue() ){
				std::ostringstream msg;
				msg << "BBClaim makes claim to (" << pos.first << ", " << pos.second << " == " << global_position
						<< " global) in pose of " << pose.total_residue() << " residues. This can result when fragments are inconsistent "
						<< "with FASTA, or label are incorrect." << std::endl;
				throw EXCN_BadInput( msg.str() );
			}
			if ( !bb_claims[ global_position ] || bb_claims[ global_position ]->right() < bb_claim->right() ) {
				bb_claims[ global_position ] = bb_claim;
			}
			else runtime_assert( bb_claim->right() < claims::DofClaim::EXCLUSIVE );
		}

		claims::JumpClaimOP jump_ptr( utility::pointer::dynamic_pointer_cast< claims::JumpClaim >( *claim ) );
		if ( jump_ptr ) {
			std::pair< std::string, core::Size > jump_start = jump_ptr->local_pos1();
			std::pair< std::string, core::Size > jump_end   = jump_ptr->local_pos2();

			Size jump_nr = ( pose.fold_tree().jump_nr( sequence_number_resolver_->find_global_pose_number( jump_start ),
																								 sequence_number_resolver_->find_global_pose_number( jump_end ) ) );
			runtime_assert( jump_nr ); //XOR would be even better
			claims::DofClaimOP already ( jumps[ jump_nr ] ); //one of them is 0 -- neutral to addition
			if ( !already || already->right() < jump_ptr->right() ) {
				jumps[ jump_nr ] = jump_ptr;
			} else runtime_assert( jump_ptr->right() < claims::DofClaim::EXCLUSIVE );
		}
	}

	if ( tr.Trace.visible() ) tr.Trace << "init-bb-dofs\n " <<  bb_claims << " init-jump-dofs\n" << jumps << std::endl;

	//check for un-initialized dofs and throw exception with dof_msg if not fully covered
	std::ostringstream dof_msg;
	bool bad( false );
	Size pos( 1 );
	for ( claims::DofClaims::const_iterator it = bb_claims.begin(); it != bb_claims.end(); ++it, ++pos ) {
		if ( !*it ) {
			bad = true;
			dof_msg << "BBTorsion at pos " << pos << "unitialized...unclaimed" << std::endl;
		}
	}

	//check for un-initialized dofs
	for ( claims::DofClaims::const_iterator it = jumps.begin(); it != jumps.end(); ++it ) {
		if ( !*it ) {
			bad = true;
			dof_msg << "Jump unitialized... " << *it << std::endl;
			bool bUnitializedJump = true;
			runtime_assert( !bUnitializedJump );
		}
	}
	if ( bad ) throw( EXCN_Input( dof_msg.str() ) );

	claims::DofClaims cumulated;
	std::copy( bb_claims.begin(), bb_claims.end(), back_inserter( cumulated ) );
	std::copy( jumps.begin(), jumps.end(), back_inserter( cumulated ) );
	//if(tr.Trace.visible())
	//{
	//	tr.Trace << "cumulated Dof Claims: " << std::endl;
	//	for(claims::DofClaims::iterator claim = cumulated.begin(); claim != cumulated.end(); ++claim)
	//	{
	//		(*claim)->show(tr.Trace);
	//		tr.Trace << std::endl;
	//	}
	//}

	claims::DofClaims failures;
	for ( TopologyClaimers::iterator top = claimers_.begin();
				top != claimers_.end(); ++top ) {
	//	if(tr.Debug.visible()) {tr.Debug << "current claimer:  " << (*top)->type() << std::endl;}
		(*top)->initialize_dofs( pose, cumulated, failures );
	}

	if ( failures.size() ) {
		std::ostringstream dof_msg;
		dof_msg << "failed to initialize dofs for these claims: .... " << failures << std::endl;
		throw EXCN_Unknown( dof_msg.str() );
	}
	runtime_assert( failures.size() == 0 ); //should have thrown exception before -- up

}

void TopologyBroker::initialize_cuts( claims::DofClaims& claims, core::pose::Pose& pose ) {
	//cuts will contain NULL for automatic cutpoints
	//and the fold-tree cut-point nr for those which we keep since they have been CUT-Claimed
	claims::DofClaims cuts( pose.num_jump(), NULL );
	for ( claims::DofClaims::iterator claim = claims.begin();	claim != claims.end(); ++claim ) {
		claims::CutClaimOP cut_ptr ( utility::pointer::dynamic_pointer_cast< claims::CutClaim >( *claim ) );
		if ( cut_ptr ){
			Size global_sequence_position = sequence_number_resolver_->find_global_pose_number( cut_ptr->get_position() );
			Size cut_nr ( pose.fold_tree().cutpoint_map( global_sequence_position ) );
			// we allow cuts without claim -- random cuts due to jumping ...
			// but if there was a CUT claim, there should be a cut
			runtime_assert( cut_nr );
			cuts[ cut_nr  ] = cut_ptr;
			//NOTE in principle we could support mutiple CUT claims at the same position... if they all just want a cut there.
			// maybe need to allow keeping of multiple claims::DofClaims per cutpoint in the list... if the list is actully ever needed.
			// we assume a claims::DofClaim::CUT will never be closed.. of course we could also change that... then however, the rights will play a role.
		}
	}

	//find the unclaimed -- automatic -- cutpoints
	to_be_closed_cuts_.clear();
	for ( Size cut_nr = 1; cut_nr<=cuts.size(); ++cut_nr ) {
		if ( !cuts[ cut_nr ] ) { //automatic cut-point
			if(pose.residue(pose.fold_tree().cutpoint(cut_nr)).type().name() == "VRT" || pose.residue(pose.fold_tree().cutpoint(cut_nr)).name3() == "XXX" ||
				pose.residue((pose.fold_tree().cutpoint(cut_nr))-1).type().name() == "VRT" || pose.residue((pose.fold_tree().cutpoint(cut_nr))-1).name3() == "XXX"
				|| pose.residue((pose.fold_tree().cutpoint(cut_nr))+1).type().name() == "VRT" || pose.residue((pose.fold_tree().cutpoint(cut_nr))+1).name3() == "XXX" )
			{
				continue;
			}else{
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

	//for ( utility::vector1< Size >::const_iterator it = to_be_closed_cuts_.begin();
	//			it != to_be_closed_cuts_.end(); ++it ) {
	//	tr.Debug << "consider cut between res " << *it << " and " << *it+1 << std::endl;
	//}
}

void TopologyBroker::apply( core::pose::Pose& pose ) {
	claims::DofClaims pre_accepted;
	bool ok( true );

	for ( TopologyClaimers::iterator top = claimers_.begin();
				top != claimers_.end(); ++top ) {
		(*top)->type();
		if ( bUseJobPose_ ) {
			(*top)->new_decoy( pose );
		} else {
			(*top)->new_decoy();
		}
	}

	//If we are using a tmh topology sampler the fold tree must be constructed differently, and preprocesing needs to occur
	bool tmh_mode = false;
	for ( TopologyClaimers::iterator top = claimers_.begin();top != claimers_.end(); ++top )
	{
		if((*top)->type()  == "TMHTopologySamplerClaimer")
		{
			tmh_mode = true;
			break;
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

	//Checking for Symmetry
	claims::SymmetryClaims symm_claims;
	generate_symmetry_claims( symm_claims );
    SymmetryClaimerOP symm_claimer = resolve_symmetry_claims ( symm_claims );

	tr.Debug << "Initialize Sequence" << std::endl;
	claims::DofClaims fresh_claims;

	sequence_claims_.clear();
	sequence_number_resolver_ = SequenceNumberResolverOP( new SequenceNumberResolver() );

	generate_sequence_claims( fresh_claims );

	if ( ok ) ok = broking( fresh_claims, pre_accepted );
	initialize_sequence( pre_accepted, pose );
    if( symm_claimer ) symm_claimer->symmetry_duplicate( pre_accepted, pose );
	current_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );

	if(!pose.empty() && basic::options::option[basic::options::OptionKeys::abinitio::explicit_pdb_debug] &&
			basic::options::option[basic::options::OptionKeys::run::protocol].value_string()=="broker")
	{
		pose.dump_pdb("pose_broker_apply_begin.pdb");
	}

	//Loop through all claimers and do pre_processing.  Kind of hacky because need to do this after broker::build_fold_tree() but before TMHTopologySampler finalizes fold tree
	if(tmh_mode)
	{
		for ( TopologyClaimers::iterator top = claimers_.begin();top != claimers_.end(); ++top )
		{
			(*top)->pre_process(pose);
		}
	}

	//if(tr.Trace.visible()){tr.Trace << "pose.fold_tree() before generate_round1:  " << pose.fold_tree() << std::endl;}
	current_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );

	tr.Debug << "Start Round1-Broking..." << std::endl;
	claims::DofClaims round1_claims;
	generate_round1( round1_claims );
	if ( ok ) ok = broking( round1_claims, pre_accepted );

	tr.Debug << "Start FinalRound-Broking..." << std::endl;
	claims::DofClaims final_claims;
	generate_final_claims( final_claims );
	if ( ok ) ok = broking( final_claims, pre_accepted );


	tr.Debug << "Broking finished" << std::endl;
	//	--> now we know nres
	current_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );

	core::kinematics::FoldTree fold_tree = pose.fold_tree();
	//if(tr.Debug.visible())
	//{
	//	tr.Debug << "pre_accepted.size():  " << pre_accepted.size() << std::endl;
	//}

	//Loop through all claimers and see if the claimer builds its own fold tree.
	//Also check to see if any
	for ( TopologyClaimers::iterator top = claimers_.begin();top != claimers_.end(); ++top )
	{
		if((*top)->claimer_builds_own_fold_tree())
		{
			use_fold_tree_from_claimer_=true;
			break;
		}
	}


	if(use_fold_tree_from_claimer_==true && tmh_mode)
	{
		build_fold_tree_from_claimer(pose, fold_tree);
	}else{
		//if(tr.Debug.visible()){tr.Debug << "build default Broker FoldTree..." << std::endl;}
		tr.Debug << "build fold-tree..." << std::endl;
		build_fold_tree( pre_accepted, pose.total_residue() );
	}
	accept_claims( pre_accepted);
	tr.Debug << *fold_tree_ << std::endl;

	//if(tr.Debug.visible()){tr.Debug << "pose.fold_tree() after build_fold_tree():  " << fold_tree << std::endl;}

	//commented out because this doesn't seem to do anything anymore after refactor - sld
	//loop through atoms of root
	//core::Size root(0);
	//if(pose.fold_tree().root() == fold_tree_->root())
	//{
	//	root = fold_tree_->root();
	//}

	pose.fold_tree( *fold_tree_ );

	//if(tr.Debug.visible())
	//{
	//	tr.Debug << "TopologyBroker Final FoldTree:  " << pose.fold_tree() << std::endl;
	//	tr.Debug << fold_tree_;
	//}
	tr.Debug << "set cuts..." << std::endl;
	initialize_cuts( pre_accepted, pose );
	tr.Debug << "initialize dofs..." << std::endl;
	initialize_dofs( pre_accepted, pose );

	assert(fold_tree_);

	//	if ( tr.Debug.visible() )	pose.dump_pdb( "init_dofs.pdb" );
	current_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );

	//we will need this one in switch_to_fullatom
	if ( !repack_scorefxn_ ) repack_scorefxn_ = core::scoring::get_score_function();

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

	if(!pose.empty() && basic::options::option[basic::options::OptionKeys::abinitio::explicit_pdb_debug] &&
			basic::options::option[basic::options::OptionKeys::run::protocol].value_string()=="broker")
	{
		pose.dump_pdb("pose_broker_apply_end.pdb");
	}
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
		//tr.Debug << (*it) << " " << pose.residue((*it)).name3() << std::endl;
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

/// @brief if some claimer wants to influence the movemap for relax he can do it here:
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
	core::optimization::MinimizerOptions options( "dfpmin_armijo_nonmonotone", 1e-5, true, false );
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap() );
	mm->set_bb( false );
	mm->set_chi( true );

	add_constraints( pose );

	core::pack::task::PackerTaskOP taskstd = core::pack::task::TaskFactory::create_packer_task( pose );
    taskstd->restrict_to_repacking();
    taskstd->or_include_current(true);
	taskstd->restrict_to_residues( needToRepack );

	if ( tr.Debug.visible() ) pose.constraint_set()->show_numbers( tr.Debug );

 	if ( pose::symmetry::is_symmetric( pose ) ) {
		protocols::simple_moves::symmetry::SymPackRotamersMover pack1( repack_scorefxn_ , taskstd );
		pack1.apply( pose );

		core::pose::symmetry::make_symmetric_movemap( pose, *mm );
		core::optimization::symmetry::SymAtomTreeMinimizer smzr;
		smzr.run( pose, *mm, *repack_scorefxn_, options );
	} else {
		protocols::simple_moves::PackRotamersMover pack1( repack_scorefxn_ , taskstd );
		pack1.apply( pose );

    // quick SC minimization
    core::optimization::AtomTreeMinimizer mzr;
		mzr.run( pose, *mm, *repack_scorefxn_, options );
	}


	tr.Debug << "minimized.. " << std::endl;
}

}
}
