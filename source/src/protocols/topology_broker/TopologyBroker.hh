// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
///           maintains list of ToplogyClaimers
///           maintains DofClaims -- exclusive or non-exclusively markedup dofs like BackboneClaim, IntraResClaim, JumpClaim
///           generates FoldTree, MoveMap, and collects samplers provided by TopologyClaimers
/// @author Oliver Lange

#ifndef INCLUDED_protocols_topology_broker_TopologyBroker_hh
#define INCLUDED_protocols_topology_broker_TopologyBroker_hh

// Unit Headers
#include <protocols/topology_broker/TopologyBroker.fwd.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <protocols/topology_broker/claims/SymmetryClaim.fwd.hh>
#include <protocols/topology_broker/claims/SequenceClaim.fwd.hh>
#include <protocols/topology_broker/SequenceNumberResolver.fwd.hh>
#include <protocols/topology_broker/ClaimerMessage.fwd.hh>
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>
#include <protocols/topology_broker/SymmetryClaimer.fwd.hh>

// Project Headers
#include <protocols/abinitio/FragmentSampler.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/conformation/symmetry/SymmData.hh>


#ifdef __clang__
#include <core/kinematics/ShortestPathInFoldTree.hh>
#endif

#include <core/pose/Pose.fwd.hh> /// FIX THIS
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <utility/assert.hh>
#include <string>

#ifdef WIN32
#include <protocols/topology_broker/TopologyClaimer.hh>
#endif


// option key includes
/* Topology Broker is the central class for the broking mechanism. The broker maintains a list of ToplogyClaimers. These Claimers
will be asked what kind of things they want to do (i.e., they return a list of DofClaims )
then the Broker asks each Claimer if foreign claims are acceptable. A claim can be disallowed. For instance, the RigidChunkClaimer,
can't accept any changes within its region of interest, and thus will disallow any claims to move BB-Torsions or to make jumps within its region.
The Broking process will be handled in two rounds ( generate_round1 and generate_final_claims ) to give Claimers the possibility to react on
the situation presented to them in round1.

The Claimer needs to keep track himself what he can move and what was declined. This is not enforced any further. The idea is that the claimer
being asked to provide a mover will return a Mover class that moves the dofs it had claimed successfully.
-> notify via claim about each change you want to make to the pose
-> don't move anything which was declined.
--> be aware that multiple movers might move the same dof -- unless your claimer takes care of declining foreign claims.

the process is carried out each time the apply is called ( since we might have different choices, e.g., different jumps, different ss-bonds, different chunks)

the life of a job goes thru stages:
first apply is called, but the state is safed for further consulting during the run of the job:

>>>apply

broking ( round1 , final_round )

setting up the pose --> pose is made from scratch!!
whatever you pass into apply will be ignored.(we call pose.clear()  )
1. ) init_sequence ( make residues from sequence )
2. ) init_dofs ( initialize dofs , from fragments (jumps, bbtorsion), from pdb (rigid chunks)

add_constraints

>>>sampling within the protocol
during the run of the samplign protocol the TopologyBroker system can be consulted for
movers ( depending on stage )
adding chainbreaks ( depending on sequence separation )
switch to fullatom ( e.g., rigid chunks provide the full-atom sidechains of their input-pdb if available )
final_fold_tree ( which chainbreaks and jumps should  the loop-closer remove before relax ? )

*/

namespace protocols {
namespace topology_broker {

class TopologyBroker : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< TopologyBroker >
{
	typedef core::Size StageID;
	typedef utility::vector1< TopologyClaimerOP > TopologyClaimers;
public:
	typedef TopologyClaimers::const_iterator const_iterator;
	///constructor
	TopologyBroker();
	virtual ~TopologyBroker();
	TopologyBroker( TopologyBroker const & );
	TopologyBroker & operator = ( TopologyBroker const & );

	/// self pointers
	inline TopologyBrokerCOP get_self_ptr() const { return shared_from_this(); }
	inline TopologyBrokerOP get_self_ptr() { return shared_from_this(); }
	inline TopologyBrokerCAP get_self_weak_ptr() const { return TopologyBrokerCAP( shared_from_this() ); }
	inline TopologyBrokerAP get_self_weak_ptr() { return TopologyBrokerAP( shared_from_this() ); }

	/// ---------------- Application Setup ------------------------------------
	/// @brief add new Claimers to the broker ---- useful before a job is started
	void add( TopologyClaimerOP cl );

	/// @brief use the input pose from the job (i.e., call new_decoy( pose ) )
	void use_job_pose( bool setting ) {
		bUseJobPose_ = setting;
	}

	/// @brief Returns true if we are using the input pose from the job
	/// (i.e. new_decoy(pose)), false otherwise.
	bool use_job_pose() const {
		return bUseJobPose_;
	}

	/// @brief Returns the ith topology claimer if it exists.
	const TopologyClaimerOP & claimer(core::Size i) const {
		assert(i >= 1 && i <= claimers_.size());
		return claimers_[i];
	}

	TopologyClaimers::const_iterator begin() const {
		return claimers_.begin();
	}

	TopologyClaimers::const_iterator end() const {
		return claimers_.end();
	}

	/// @brief Returns the number of claimers associated with the broker
	core::Size num_claimers() const {
		return claimers_.size();
	}

	/// ----------------------- Job Setup ------------------------------------------
	/// @brief at the start of a job this is called, e.g., by the AbrelaxMover
	/// it generates a pose with appropriate foldtree and initializes dofs, adds constraints, etc.
	void apply( core::pose::Pose & );

	////----------------------------------  Consulting  ----------------------------------------------------
	/// the following interface is for the Mover to consult the Broker during the course of a simulation
	/// usually these calls a relayed to the Claimers only some will answer, others will ignore it.
	/// e.g., ConstraintClaimer will not add a Mover, JumpClaimer and RigidChunkClaimer will add chainbreaks
	/// a basic FragmentClaimer will not

	/// @brief return a set of Movers ( RandomMover, i.e. container of movers )
	moves::MoverOP mover( core::pose::Pose const &, abinitio::StageID, core::scoring::ScoreFunction const & scorefxn, core::Real progress ) const;

	/// @brief apply filter (TopologyClaimer::passes_filter() ) and raise exception EXCN_FILTER_FAILED if failed
	void apply_filter( core::pose::Pose const &, abinitio::StageID, core::Real progress ) const;

	/// @brief if some claimer wants to influence the movemap for relax he can do it here:
	void adjust_relax_movemap( core::kinematics::MoveMap & ) const;

	/// @brief the SlidingWindowLoopClosure needs pure fragments, because it changes the the residue number in the ShortLoopClosure part
	/// thus extra hook for this --- > only some Claimers will answer
	core::fragment::FragSetCOP loop_frags( core::kinematics::MoveMap & ) const;

	/// @brief do we need to close loops ( unphysical chainbreaks have been introduced? )
	bool has_chainbreaks_to_close() const;

	/// @brief add chainbreak variant residue to the unphysical chainbreaks
	void add_chainbreak_variants( core::pose::Pose & pose, core::Size max_dist = 0, core::kinematics::ShortestPathInFoldTreeCOP sp = NULL) const;
	/// @brief check that each chainbreak residue has a chainbreak variant
	bool check_chainbreak_variants( core::pose::Pose & pose ) const;

	/// @brief switch to fullatom --- some Claimers might help by providing template based side-chain information
	void switch_to_fullatom( core::pose::Pose & );

	bool does_final_fold_tree_exist() const
	{
		if ( final_fold_tree_ ) {
			return true;
		}
		return false;
	}

	/// @brief access for hacky claimers
	core::kinematics::FoldTree & final_fold_tree() const {
		//std::cout << "Broker FinalFoldTree is:  ";
		//final_fold_tree_->show(std::cout);
		runtime_assert( final_fold_tree_ != 0 );
		return *final_fold_tree_;
	};

	/// @brief get the sequence claim that is consistent with the label,
	/// throws EXCN_Unknown_SequenceLabel if not found
	claims::SequenceClaim & resolve_sequence_label( std::string const & label ) const;

	//core::Size resolve_residue( std::string const & chain_label, core::Size pos ) const;

	const SequenceNumberResolver & sequence_number_resolver() const {
		runtime_assert( sequence_number_resolver_ != 0 );
		return *sequence_number_resolver_;
	}

	void relay_message( ClaimerMessage & msg ) const;
	//// ------------------------------- End Consulting --------------------------------------------

	bool has_sequence_claimer();

	core::pose::Pose const & current_pose() const {
		return *current_pose_;
	}

private:
	/// @brief first round claims are collected
	void generate_sequence_claims( claims::DofClaims & all_claims );

	/// @brief collects symmetry claims
	void generate_symmetry_claims( claims::SymmetryClaims & all_claims );

	/// @brief checks whether only one sequence claim is there, otherwise crashes.
	SymmetryClaimerOP resolve_symmetry_claims( claims::SymmetryClaims & symm_claims );

	void make_sequence_symmetric( claims::DofClaims pre_accepted, core::pose::Pose & pose);


	/// @brief first round claims are collected
	void generate_round1( claims::DofClaims & all_claims );

	/// @brief second round claims are collected
	void generate_final_claims( claims::DofClaims & all_claims );

	/// @brief notify owner of accepted claims
	void accept_claims( claims::DofClaims & claims );

	/// @brief run thru list of claims, ask all claimers if this claims is acceptable --- > returns accepted claims in pre_accepted
	/// throws EXCN_ExclusiveClaimDeclined if the call to the owners TopologyClaimer::accept_declined_claim( declined_claim ) returns false
	bool broking( claims::DofClaims const & all_claims, claims::DofClaims & pre_accepted );

	/// @brief creates a fold-tree from the Jump- and CutClaims
	/// throws EXCN_InvalidFoldTree at failure
	void build_fold_tree( claims::DofClaims & claims, Size nres );

	void build_fold_tree_from_claimer(core::pose::Pose & pose, core::kinematics::FoldTree & fold_tree);

	/// @brief create new pose from SeqClaims
	void initialize_sequence( claims::DofClaims & claims, core::pose::Pose & new_pose );

	/// @brief creates the list "to_be_closed_cuts_" from current fold-tree and CutClaims
	void initialize_cuts( claims::DofClaims & claims, core::pose::Pose & new_pose );

	/// @brief initialize dofs
	void initialize_dofs( claims::DofClaims & claims, core::pose::Pose & new_pose );

	/// @brief add constraints --> referred to Claimers ( e.g., ConstraintClaimer, RigidChunkClaimer )
	void add_constraints( core::pose::Pose & ) const;

private:
	/// @brief vector of Claimers --- RigidChunkClaimer, FragmentClaimer, ConstraintClaimer, etc.
	TopologyClaimers claimers_;

	//=============================================================================
	//all these are derived infos and change each time that we use apply ( generate a new pose with foldtree etc )

	/// @brief list of dof-claims currently active
	claims::DofClaims current_claims_;

	/// @brief current pose has nres total_residues
	core::Size nres_;

	/// @brief the current fold-tree
	core::kinematics::FoldTreeOP fold_tree_;

	//// mutable here is hack to make wrong MembraneTopologyClaimer work --- don't want to clean up that guy right now...
	/// @brief current final-fold-tree --- after removal to_bel_closed_cuts_
	mutable core::kinematics::FoldTreeOP final_fold_tree_;

	/// @brief Scorefunction used in switch_to_fullatom
	core::scoring::ScoreFunctionOP repack_scorefxn_;

	/// @brief these cuts are not physical and should be closed ( i.e., chainbreak energy, loop-closing )
	utility::vector1< Size > to_be_closed_cuts_; //keeps the residue number not the cut-nr --- thats safer.

	/// @brief these cuts are not physical and should be closed ( i.e., chainbreak energy, loop-closing )
	utility::vector1< Size > start_pose_cuts_; //keeps the residue number not the cut-nr --- thats safer.

	claims::SequenceClaims sequence_claims_;

	SequenceNumberResolverOP sequence_number_resolver_;
	/// @brief we restart from the input pose... call steal( pose ) for all claimers
	bool bUseJobPose_;

	bool use_fold_tree_from_claimer_;

	//bool have_fold_tree_;

	core::pose::PoseOP current_pose_;

	//DofClaims pre_accepted_;
}; //class TopologyBroker

}
}

#endif  // INCLUDED_protocols_topology_broker_TopologyBroker_hh
