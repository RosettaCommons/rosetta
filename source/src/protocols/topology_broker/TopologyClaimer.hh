// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/topology_broker/TopologyClaimer.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_topology_broker_TopologyClaimer_hh
#define INCLUDED_protocols_topology_broker_TopologyClaimer_hh

// Unit headers
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>

// Package headers
#include <protocols/topology_broker/TopologyBroker.fwd.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/claims/SymmetryClaim.hh>
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/SequenceClaim.hh>
#include <protocols/topology_broker/weights/AbinitioMoverWeight.hh>
#include <protocols/topology_broker/weights/ConstAbinitioMoverWeight.hh>
#include <protocols/topology_broker/ClaimerMessage.fwd.hh>

// Project Headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/fragment/FragSet.hh>
#include <protocols/abinitio/FragmentSampler.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1_bool.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace topology_broker {

/// @detail
//    the claim functions are called from the Broker in this sequence
//     generate_sequence_claims()
//    initialize_residues() /// puts in the correct residue-types can copy segments of pose
//    1 x generate_claims()
//    n x allow_claim()
//    [ evtl. n x accept_declined_claim() ]

//    // here we don't call claim_accepted because I want to be able to take individual jumps out if they make fold-trees impossible

//    1 x finalize_claims()
//    n x allow_claim()
//    [ evtl. n x accept_declined_claim() ]
//    n x claim_accepted( ) //hopefully ;-)

//    initialize_dofs( init_claims [ subset of your claims ] ) ///only act on internal dofs -- no structure building

class TopologyClaimer : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< TopologyClaimer >
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~TopologyClaimer(){}

	TopologyClaimer() :
		abinitio_mover_weight_ ( weights::AbinitioMoverWeightOP( new weights::ConstAbinitioMoverWeight( 1.0 ) ) ),
		label_( "NO_LABEL" )
	{};

	/// @brief construct with weight-set ( how important is mover for different abinitio stages ? )
	TopologyClaimer( weights::AbinitioMoverWeightOP weight ) :
		abinitio_mover_weight_ ( weight ),
		label_( "NO_LABEL" )
	{
		if ( !weight ) abinitio_mover_weight_ = weights::AbinitioMoverWeightOP( new weights::ConstAbinitioMoverWeight( 1.0 ) );
	};

	/// self pointers
	inline TopologyClaimerCOP get_self_ptr() const { return shared_from_this(); }
	inline TopologyClaimerOP get_self_ptr() { return shared_from_this(); }
	inline TopologyClaimerCAP get_self_weak_ptr() const { return TopologyClaimerCAP( shared_from_this() ); }
	inline TopologyClaimerAP get_self_weak_ptr() { return TopologyClaimerAP( shared_from_this() ); }

	/// @brief clone it!
	virtual TopologyClaimerOP clone() const = 0;

	/// @brief name of Claimer
	virtual std::string type() const = 0;

	virtual bool claimer_builds_own_fold_tree()
	{
		return false;
	}

	// virtual void set_pose_from_broker(core::pose::Pose& pose) {};

	/// @brief in case a claimer has its own fold_tree.  get_fold_tree() is called by the broker
	virtual core::kinematics::FoldTreeOP get_fold_tree(core::pose::Pose&) {return NULL;}

	//virtual core::pose::PoseOP get_pose() {return NULL;}

	virtual void build_fold_tree(core::pose::Pose&, core::kinematics::FoldTree&) {};

	virtual void set_pose_from_broker(core::pose::Pose&){};

	/// @brief read definition of Claimer from setup file, i.e., a CLAIMER <type> ... END_CLAIMER block
	virtual void read( std::istream & );

	virtual void pre_process(core::pose::Pose&) {};

	/// @brief generate claims that affect the sequence of the pose
	virtual void generate_sequence_claims( claims::DofClaims& ){}; //add to list ( never call clear() on list )

	/// @brief generate claims that affect the sequence of the pose
	virtual void generate_symmetry_claims( claims::SymmetryClaims& ){}; //add to list ( never call clear() on list )

	/// @brief generate first round of DOF claims
	virtual void generate_claims( claims::DofClaims& ) {}; //add to list ( never call clear() on list )

	/// @brief is called after all round1 claims have been approved or retracted -- additional claims can be issued in this round
	virtual void finalize_claims( claims::DofClaims& ) {};

	/// @brief allow a claim from a foreign Claimer
	virtual bool allow_claim( claims::DofClaim const& /*foreign_claim*/ ) { return true; };

	/// @brief initialize sequence ( for approved sequence claims given as init_claim ) Claimer searches init_claims for claims owned by *this
	//void initialize_residues( core::pose::Pose&, claims::SequenceClaimOP init_claim, claims::DofClaims& failed_to_init ){}

	/// @brief initialize dofs -- e.g., torsions, jumps -- Claimer searches init_claims for claims owned by *this
	virtual void initialize_dofs( core::pose::Pose&, claims::DofClaims const& init_claims, claims::DofClaims& failed_to_init );

	/// @brief has this Claimer some side chain conformations to add?
	/// starts with bNeedToRepack true for all residues... if you have a sidechain ---> copy it to pose and set needtoRepack false for this residue
	virtual void switch_to_fullatom( core::pose::Pose&, utility::vector1< bool > /* bNeedToRepack */ ) const {};

	/// @brief multiply your bias to this -- if its zero don't change that, i.e., multiply only
	/// this is used during fold-tree generation to set the cut-points. it starts with 1.0 for all residues.
	/// Fragments can add their loop-fraction
	virtual void manipulate_cut_bias( utility::vector1< core::Real >& /*cut_bias*/ ) {};

	/// @brief notification of declined claims: update your internal representation (e.g., movemap ) to remember this !
	//// return false   -- if you can't live without this claim being accepted. ( e.g., RigidChunks ... )
	virtual bool accept_declined_claim( claims::DofClaim const& /*was_declined*/ ) { return true; }

	/// @brief this claim of yours was accepted.... I so far haven't done anything with this... might go away.
	virtual void claim_accepted( claims::DofClaimOP my_claim ) {
		current_claims_.push_back( my_claim );
	}

	/// @brief add constraints to pose...  might make this stage dependent as with movers...
	virtual void add_constraints( core::pose::Pose& /*pose*/ ) const {};//some constraints are loaded...

	/// @brief return fragments that can be used for loop-sampling... unfortunately some loop-samplers need fragments, rather then fragmovers
	/// (e.g. short-loop closure since it remaps them on a short pose containing only the loop-residues. )
	/// overloaded e.g., by LoopFragmentClaimer.. returns a movemap and fragset good for loop-sampling
	virtual core::fragment::FragSetCOP loop_frags( core::kinematics::MoveMap& /*returned! */ ) const { return NULL; };

	/// @brief claimers can add movers to the RandomMover (Container).
	/// add your moves, make it dependent on stage if you want to. So far this is called only by abinitio...
	/// if you don't want to do anything special --- don't overload this method!
	/// default: adds mover given by virtual call get_mover()  with stage-dependent weight given by abinitio_mover_weight_
	virtual bool passes_filter(
		core::pose::Pose const& /*pose*/,
		abinitio::StageID /*stageID*/, /* abinitio sampler stage */
		core::Real, /* progress */ /* progress within stage */
		std::ostringstream& /*report*/
	) { return true; }

	/// @brief claimers can add movers to the RandomMover (Container).
	/// add your moves, make it dependent on stage if you want to. So far this is called only by abinitio...
	/// if you don't want to do anything special --- don't overload this method!
	/// default: adds mover given by virtual call get_mover()  with stage-dependent weight given by abinitio_mover_weight_
	virtual void add_mover(
		moves::RandomMover& /* random_mover */,
		core::pose::Pose const& /*pose*/,
		abinitio::StageID /*stageID*/, /* abinitio sampler stage */
		core::scoring::ScoreFunction const& /*scorefxn*/,
		core::Real /* progress */ /* progress within stage */
	);

	void set_mover_weight( weights::AbinitioMoverWeightOP wset ) {
		abinitio_mover_weight_ = wset;
	}

	weights::AbinitioMoverWeight& mover_weight() {
		return *abinitio_mover_weight_;
	}

	/// @brief read mover weight from Stream. - so far recognizes: LargeStage, SmallStage, SmoothStage, AllStage
	void read_mover_weight( std::istream& is );

	virtual void adjust_relax_movemap( core::kinematics::MoveMap& ) const {};

	/// @brief don't use this --- it is called by TopologyBroker add_claim only
	void set_broker( TopologyBrokerCAP ptr ) {
		broker_ = ptr;
	}

	/// @brief return the broker we are collaborating with
	TopologyBroker const& broker() const{
		TopologyBrokerCOP broker( broker_ );
		runtime_assert( broker != 0 ); //don't use before claimer is added to its broker
		return *broker; // NEEDS FIXING: returning reference to temporairly locked weak_ptr
	}

	/// @brief a new decoy --- random choices to be made ? make them here
	virtual void new_decoy() {};

	/// @brief an input pose is given, i.e., a starting structure for resampling
	/// don't make random choices, base choices on given pose
	virtual void new_decoy( core::pose::Pose const& ) { new_decoy(); };

	std::string const& label() const {
		return label_;
	}

	void set_label( std::string const& str ) {
		label_ = str;
	}

	virtual void receive_message( ClaimerMessage& ) {};

	//getter function to retrieve this claimer's current_pose (to get after moving spans)
	virtual core::pose::PoseOP get_pose_from_claimer() {return NULL;};

protected:

	/// @brief what is your mover ... called by add_mover --- overload this or add_mover if you have movers too supply
	virtual moves::MoverOP get_mover( core::pose::Pose const& /*pose*/ ) const { return NULL; };

	virtual bool read_tag( std::string tag, std::istream& is );

	virtual void set_defaults(); //eg before reading starts.

	virtual void init_after_reading() {};

private:

	void unknown_tag( std::string tag, std::istream& is ) const;

	/// @brief currently accepted claims... dunno haven't used that yet.
	claims::DofClaims current_claims_;

	/// @brief weight set .. how often shall this claimer's mover be sampled during different stages ?
	weights::AbinitioMoverWeightOP abinitio_mover_weight_;

	/// @brief oh, our master. (use broker() to get this ).
	TopologyBrokerCAP broker_;

	/// @brief a user defined string, can be used to send messages from claimer to claimer
	std::string label_;

	/// @brief in case a claimer has its own fold_tree.  get_fold_tree() is called by the broker
	//core::kinematics::FoldTreeOP fold_tree_;

	//core::pose::Pose pose_;

}; //class TopologyClaimer

inline std::istream& operator >> ( std::istream& is, TopologyClaimer& tc ) {
	tc.read( is );
	return is;
}

}
}

#endif
