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
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_FragmentClaimer_hh
#define INCLUDED_protocols_topology_broker_FragmentClaimer_hh


// Unit Headers
#include <protocols/topology_broker/FragmentClaimer.fwd.hh>

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/simple_moves/FragmentMover.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <protocols/loops/Loops.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace topology_broker {

class FragmentClaimer : public virtual TopologyClaimer {
	typedef TopologyClaimer Parent;
public:
	FragmentClaimer(); //for factory
	FragmentClaimer( simple_moves::FragmentMoverOP, std::string mover_tag, weights::AbinitioMoverWeightOP weight );
	FragmentClaimer( simple_moves::FragmentMoverOP, std::string mover_tag, weights::AbinitioMoverWeightOP weight, std::string label, core::fragment::FragSetOP fags );
	FragmentClaimer( simple_moves::FragmentMoverOP );
	FragmentClaimer( FragmentClaimer const & src );

	~FragmentClaimer() override;

	FragmentClaimerOP shared_from_this() { return utility::pointer::dynamic_pointer_cast<FragmentClaimer>( TopologyClaimer::shared_from_this() ); }

	TopologyClaimerOP clone() const override {
		return TopologyClaimerOP( new FragmentClaimer( *this ) );
	}

	void generate_claims( claims::DofClaims& ) override;

	/// @brief is called after all round1 claims have been approved or retracted -- additional claims can be issued in this round
	//virtual DofClaims finalize_claims( DofClaims& );

	void initialize_dofs( core::pose::Pose&, claims::DofClaims const& init_claims, claims::DofClaims& failed_to_init ) override;

	bool accept_declined_claim( claims::DofClaim const& was_declined ) override;

	/// @brief type() is specifying the output name of the TopologyClaimer
	std::string type() const override {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "FragmentClaimer";
	}

	void set_mover( simple_moves::FragmentMoverOP mover );

	void set_mover_tag( std::string const& str );

	void set_bInitDofs( bool set ) {
		bInitDofs_ = set;
	}

	std::string const& mover_tag() const {
		return mover_tag_;
	}

	moves::MoverOP get_mover( core::pose::Pose const& /*pose*/ ) const override;

	void get_sequence_region( std::set< core::Size >& start_region ) const;

	/*void current_offset(core::Size offset){
	current_offset_ = offset;
	}

	core::Size current_offset(){
	return current_offset_;
	}*/


protected:
	simple_moves::FragmentMover const & mover() const {
		if ( !mover_ ) throw( utility::excn::EXCN_NullPointer( "mover_ is NULL in FragmentClaimer::mover()" ) );
		return *mover_;
	}

	simple_moves::FragmentMoverOP get_frag_mover_ptr();

	core::fragment::FragSetCOP fragments();

	void set_fragments( core::fragment::FragSetOP );

	void set_claim_right( claims::DofClaim::ClaimRight setting ) {
		claim_right_ = setting;
	}

	bool read_tag( std::string tag, std::istream & ) override;

	core::kinematics::MoveMapOP movemap_;

	void init_after_reading() override;

private:

	simple_moves::FragmentMoverOP mover_;
	std::string mover_tag_;

	//if false the initialize_dofs routine won't do anything -- but also not report the dof as failed-to initialized
	bool bInitDofs_;

	claims::DofClaim::ClaimRight claim_right_; /*default CAN_INIT */

	/// @brief regions that can be used for fragment insertions
	loops::Loops region_;

	/// @brief if non-empty this claimer operates only on the sequences with these labels...
	/// create std::set< Size > with residue numbers by get_sequence_region();
	utility::vector1< std::string > active_sequence_labels_;

	/*
	/// @brief Stores offset of current FragmentClaimer (based on global sequence)
	core::Size current_offset_;
	*/

}; //class FragmentClaimer

}
}

#endif
