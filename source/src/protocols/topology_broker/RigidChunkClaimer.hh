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
/// @author Oliver Lange

#ifndef INCLUDED_protocols_topology_broker_RigidChunkClaimer_hh
#define INCLUDED_protocols_topology_broker_RigidChunkClaimer_hh

// Unit Headers
#include <protocols/topology_broker/RigidChunkClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.fwd.hh>
#include <protocols/topology_broker/TopologyClaimer.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <protocols/loops/Loops.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/topology_broker/ClaimerMessage.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace topology_broker {

///@brief defines a rigid part of structure... imagine a loop-relax application core structure is fixed via jumps and loops can move
///@detail the rigid chunk takes a definition of rigid regions in form of an instance of Loops (just taken as bunch of start-end residue numbers ---  here defining the rigid residues and not the loops).
/// the rigid chunk to keep its integrity will need jumps, the claimer will reuse jumps if somebody else claims them,
/// or submit in finalize_claims his own jumps, if not enough jumps are present.
/// in "bExclusive_ mode" the RigidChunk will reclaim any jump claim that is useful and wihin the rigid region. (i.e., foreign claim is dissallowed but own claim with same residues is issued --- in this way the claimer uses e.g., beta-sheet jumps, where they are suggested
/// the input pose is used to initialize the rigid region ( via copying of internal coordinates )
///  e.g., a hole in the structure shouldn't pose a problem, since we basically copy the atom-tree.
class RigidChunkClaimer : public TopologyClaimer {
typedef TopologyClaimer Parent;

public:
	class CM_SuggestFixResidue : public ClaimerMessage {
	public:
		CM_SuggestFixResidue( std::string to ) : ClaimerMessage( to ) {};
		Size good_fix_pos_;
	};

public:
	//c'stor
	RigidChunkClaimer();
	RigidChunkClaimer( core::pose::Pose const& input_pose, loops::Loops rigid );

	//clone
	virtual TopologyClaimerOP clone() const {
		return TopologyClaimerOP( new RigidChunkClaimer( *this ) );
	}

	///@brief type() is specifying the output name of the TopologyClaimer
	virtual std::string type() const {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "RigidChunkClaimer";
	}

	virtual void new_decoy();
	virtual void new_decoy( core::pose::Pose const& );

	///@brief generate DofClaims for BB
	virtual void generate_claims( claims::DofClaims& ); //add to list ( never call clear() on list )

	///@brief has to decline foreign BB claims for rigid regions, reclaim jumps where appropriate
	virtual bool allow_claim( claims::DofClaim const& /*foreign_claim*/ );

	///@brief issue jump-claims for jumps yet missing to keep rigid regions fixed
	virtual void finalize_claims( claims::DofClaims& );

	//	virtual void initialize_residues( core::pose::Pose&, DofClaims const& init_claims, DofClaims& failed_to_init );
	///@brief initialize BB residues and rigid-internal jumps from starting structure --- copying atom-tree dofs
	virtual void initialize_dofs( core::pose::Pose&, claims::DofClaims const& init_claims, claims::DofClaims& failed_to_init );

	///@brief rigid-chunk can probably provide some side-chain info from full-length model
	virtual void switch_to_fullatom( core::pose::Pose&, utility::vector1< bool > bNeedToRepack ) const;

	///@brief will fail if a BB torsion claim of the rigid region has been declined
	virtual bool accept_declined_claim( claims::DofClaim const& was_declined );

	///@brief multiply your bias to this -- if its zero don't change that, i.e., multiply only
	virtual void manipulate_cut_bias( utility::vector1< core::Real >& cut_bias );

	///@brief disallow torsion moves in relax if bRigidInRelax
	virtual void adjust_relax_movemap( core::kinematics::MoveMap& ) const;

	// will be required when we have the option to use coord. csts to fix the rigid chunk.
	//virtual void add_constraints( core::pose::Pose& /*pose*/ );
	//????	virtual void add_score_weights( core::scoring::ScoreFunction& );
	virtual void receive_message( ClaimerMessage& cm );

	///@brief Returns true if we are using loop definitions from ThreadingJob
	bool use_loops_from_threading_job() const {
		return bUseThreadingJobLoops_;
	}

	///@brief Sets whether we should use loop definitions from ThreadingJob
	void use_loops_from_threading_job(bool setting) {
		bUseThreadingJobLoops_ = setting;
	}

protected:
	///@brief select sub-regions from rigid_core_, if skip-rate is specified
	void select_parts();

	virtual bool read_tag( std::string tag, std::istream& is );

	virtual void set_defaults(); //eg before reading starts.

	virtual void init_after_reading();


private:
	///@brief starting pose
	core::pose::Pose input_pose_;

	///@brief starting pose in centroid mode
	core::pose::Pose centroid_input_pose_;

	///@brief regions that can be used for rigid core
	loops::Loops rigid_core_;

	///@brief if skip-rate is given (in loop-definitions) current_rigid_core_ will contain the current "choice" of regions... set in generate_claims()
	loops::Loops current_rigid_core_;

	///@brief jumps used this round --- since generate_claims()
	claims::DofClaims current_jumps_;

	///@brief flag used to specify if the rigid regions should really be treated exclusivity --- i.e., are they really rigid ?
	///@brief changing this flag to false, will allow everything to move, but together with coordinate constraints this might yield
	/// a pose that is easier sampled ?
	bool bExclusive_; //really rigid?

	///@brief jump residues that are just next to rigid region are probably leading to impossibl fold-trees (not always -- but most of the time)
	/// if false we don't allow such jumps
	bool bAllowAdjacentJumps_;

	///@brief use the pose in new_decoy( pose )
	bool bUseInputPose_;

	///@brief use loop-definition from ThreadingJob
	bool bUseThreadingJobLoops_;

	///@brief min_loop_size for Threading-loops
	core::Size min_loop_size_;

	///@brief same effect as OptionKeys::loops::random_grow_loops_by ]() for looprelax
	core::Real random_grow_loops_by_;

	///@brief keep this chunk rigid in relax --- adjust movemap to keep BB-Torsions fixed...
	bool bRigidInRelax_;


	///@brief helper class -- computes if we have all jupms needed to rigidify the chosen chunk and generate more jumps if needed.
	class JumpCalculator : public utility::pointer::ReferenceCount { //helper class do we like this jump, do we need more ?
	public:
		JumpCalculator( loops::Loops const& rigid_,	bool bAllowAdjacentJumps );

		///@brief only called for relevant jumps:
		///*true* if this jump helps us keeping things rigid,
		///*false* if this jump connects rigid regions
		// that are already rigidified via a different jump.
		bool good_jump( core::Size pos1, core::Size pos2 );

		///@brief this jump doesn't help --- it doens't touch two of the rigid regions
		bool irrelevant_jump( core::Size pos1, core::Size pos2 );

		///@brief get the missing jumps
		void generate_rigidity_jumps( RigidChunkClaimer*, claims::DofClaims& extra_jumps, std::string );

	private:
		///@brief what should be rigid
		loops::Loops rigid_;

		///@brief which residues have already been connected via jumps
		utility::vector1< Size > visited_;

		///@brief how many new jumps
		Size new_nr_;

		///@brief use loop-definition from alignment in ThreadingJob
		// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
		// bool bUseThreadingJobLoops_;
		bool bAllowAdjacentJumps_;

	};

	// Types
	typedef utility::pointer::shared_ptr< JumpCalculator >  JumpCalculatorOP;
	JumpCalculatorOP current_jump_calculator_;

}; //class RigidChunkClaimer

}
}

#endif

