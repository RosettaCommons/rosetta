// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/topology_broker/TopologyClaimer.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_topology_broker_TopologyClaimer_hh
#define INCLUDED_protocols_topology_broker_TopologyClaimer_hh

// Unit headers
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>

// Package headers
#include <protocols/topology_broker/TopologyBroker.fwd.hh>
#include <protocols/topology_broker/DofClaim.hh>
#include <protocols/topology_broker/weights/AbinitioMoverWeight.hh>
#include <protocols/topology_broker/weights/ConstAbinitioMoverWeight.hh>
#include <protocols/topology_broker/ClaimerMessage.fwd.hh>

// Project Headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/fragment/FragSet.hh>
#include <protocols/abinitio/FragmentSampler.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
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
// 		 initialize_residues() /// puts in the correct residue-types can copy segments of pose
// 		 1 x generate_claims()
// 		 n x allow_claim()
// 		 [ evtl. n x accept_declined_claim() ]

// 		 // here we don't call claim_accepted because I want to be able to take individual jumps out if they make fold-trees impossible

// 		 1 x finalize_claims()
// 		 n x allow_claim()
// 		 [ evtl. n x accept_declined_claim() ]
// 		 n x claim_accepted( ) //hopefully ;-)

// 		 initialize_dofs( init_claims [ subset of your claims ] ) ///only act on internal dofs -- no structure building

class TopologyClaimer : public utility::pointer::ReferenceCount {
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~TopologyClaimer();
	TopologyClaimer() : abinitio_mover_weight_ ( new weights::ConstAbinitioMoverWeight( 1.0 ) ) {};

	///@brief construct with weight-set ( how important is mover for different abinitio stages ? )
	TopologyClaimer( weights::AbinitioMoverWeightOP weight ) : abinitio_mover_weight_ ( weight ) {
		if ( !weight ) abinitio_mover_weight_ = new weights::ConstAbinitioMoverWeight( 1.0 );
	};

	///@brief clone it!
	virtual TopologyClaimerOP clone() const = 0;

	///@brief name of Claimer
	virtual std::string type() const = 0;


	///@brief read definition of Claimer from setup file, i.e., a CLAIMER <type> ... END_CLAIMER block
	virtual void read( std::istream & );

	///@brief generate claims that affect the sequence of the pose
	virtual void generate_sequence_claims( DofClaims& ){}; //add to list ( never call clear() on list )

	///@brief generate first round of DOF claims
	virtual void generate_claims( DofClaims& ) {}; //add to list ( never call clear() on list )

	///@brief is called after all round1 claims have been approved or retracted -- additional claims can be issued in this round
	virtual void finalize_claims( DofClaims& ) {};

	///@brief allow a claim from a foreign Claimer
	virtual bool allow_claim( DofClaim const& /*foreign_claim*/ ) { return true; };

	///@brief initialize sequence ( for approved sequence claims given as init_claim ) Claimer searches init_claims for claims owned by *this
	virtual void initialize_residues( core::pose::Pose&, SequenceClaimOP init_claim, DofClaims& failed_to_init );

	///@brief initialize dofs -- e.g., torsions, jumps -- Claimer searches init_claims for claims owned by *this
	virtual void initialize_dofs( core::pose::Pose&, DofClaims const& init_claims, DofClaims& failed_to_init );

	///@brief has this Claimer some side chain conformations to add?
	/// starts with bNeedToRepack true for all residues... if you have a sidechain ---> copy it to pose and set needtoRepack false for this residue
	virtual void switch_to_fullatom( core::pose::Pose&, utility::vector1< bool > /* bNeedToRepack */ ) const {};

	///@brief multiply your bias to this -- if its zero don't change that, i.e., multiply only
	/// this is used during fold-tree generation to set the cut-points. it starts with 1.0 for all residues.
	/// Fragments can add their loop-fraction
	virtual void manipulate_cut_bias( utility::vector1< core::Real >& /*cut_bias*/ ) {};

	///@brief notification of declined claims: update your internal representation (e.g., movemap ) to remember this !
	//// return false   -- if you can't live without this claim being accepted. ( e.g., RigidChunks ... )
	virtual bool accept_declined_claim( DofClaim const& /*was_declined*/ ) { return true; }

	///@brief this claim of yours was accepted.... I so far haven't done anything with this... might go away.
	virtual void claim_accepted( DofClaimOP my_claim ) {
		current_claims_.push_back( my_claim );
	}

	///@brief add constraints to pose...  might make this stage dependent as with movers...
	virtual void add_constraints( core::pose::Pose& /*pose*/ ) const {};//some constraints are loaded...

	///@brief return fragments that can be used for loop-sampling... unfortunately some loop-samplers need fragments, rather then fragmovers
	/// (e.g. short-loop closure since it remaps them on a short pose containing only the loop-residues. )
	/// overloaded e.g., by LoopFragmentClaimer.. returns a movemap and fragset good for loop-sampling
	virtual core::fragment::FragSetCOP loop_frags( core::kinematics::MoveMap& /*returned! */ ) const { return NULL; };

	///@brief claimers can add movers to the RandomMover (Container).
	/// add your moves, make it dependent on stage if you want to. So far this is called only by abinitio...
	/// if you don't want to do anything special --- don't overload this method!
	/// default: adds mover given by virtual call get_mover()  with stage-dependent weight given by abinitio_mover_weight_
	virtual bool passes_filter(
		core::pose::Pose const& /*pose*/,
		abinitio::StageID /*stageID*/, /* abinitio sampler stage */
		core::Real, /* progress */ /* progress within stage */
		std::ostringstream& /*report*/
	) { return true; }

	///@brief claimers can add movers to the RandomMover (Container).
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

	///@brief read mover weight from Stream. - so far recognizes: LargeStage, SmallStage, SmoothStage, AllStage
	void read_mover_weight( std::istream& is );

	virtual void adjust_relax_movemap( core::kinematics::MoveMap& ) const {};

	///@brief don't use this --- it is called by TopologyBroker add_claim only
	void set_broker( TopologyBrokerCAP ptr ) {
		broker_ = ptr;
	}

	///@brief return the broker we are collaborating with
	TopologyBroker const& broker() const{
		runtime_assert( broker_ ); //don't use before claimer is added to its broker
		return *broker_;
	}

	///@brief a new decoy --- random choices to be made ? make them here
	virtual void new_decoy() {};

	///@brief an input pose is given, i.e., a starting structure for resampling
	/// don't make random choices, base choices on given pose
	virtual void new_decoy( core::pose::Pose const& ) { new_decoy(); };

	std::string const& label() const {
		return label_;
	}

	void set_label( std::string const& str ) {
		label_ = str;
	}

	virtual void receive_message( ClaimerMessage& ) {};

protected:

	///@brief what is your mover ... called by add_mover --- overload this or add_mover if you have movers too supply
	virtual moves::MoverOP get_mover(	core::pose::Pose const& /*pose*/ ) const { return NULL; };

	virtual bool read_tag( std::string tag, std::istream& is );

	virtual void set_defaults(); //eg before reading starts.

	virtual void init_after_reading() {};

private:

	void unknown_tag( std::string tag, std::istream& is ) const;

	///@brief currently accepted claims... dunno haven't used that yet.
	DofClaims current_claims_;

	///@brief weight set .. how often shall this claimer's mover be sampled during different stages ?
	weights::AbinitioMoverWeightOP abinitio_mover_weight_;

	///@brief oh, our master. (use broker() to get this ).
	TopologyBrokerCAP broker_;

	///@brief a user defined string, can be used to send messages from claimer to claimer
	std::string label_;
}; //class TopologyClaimer

inline std::istream& operator >> ( std::istream& is, TopologyClaimer& tc ) {
	tc.read( is );
	return is;
}

}
}

#endif
