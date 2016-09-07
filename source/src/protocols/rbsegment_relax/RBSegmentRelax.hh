// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief RBSegmentRelax protocol
/// @details
///
///
///
/// @author Srivatsan Raman
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_rbsegment_relax_RBSegmentRelax_hh
#define INCLUDED_protocols_rbsegment_relax_RBSegmentRelax_hh

// Package headers
#include <protocols/rbsegment_relax/RBSegment.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>


//C++ headers
#include <map>
#include <list>

#include <core/fragment/FragSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace rbsegment_relax {


//typedef utility::pointer::owning_ptr< RBSegment >  RBSegmentOP;
class RBSegmentRelax : public protocols::moves::Mover {

public:
	RBSegmentRelax();

	RBSegmentRelax( core::scoring::ScoreFunctionOP scorefxn,
		utility::vector1< RBSegment > const &RBSegment_input,
		protocols::loops::Loops const &Loops_input );
	~RBSegmentRelax() override;

	void initialize( utility::vector1< core::fragment::FragSetOP > const &frag_libs , core::Real rnd=0.0);
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	void add_helixMover( RBSegmentMoverOP newMover );
	void add_strandMover( RBSegmentMoverOP newMover );
	void add_genericRBMover( RBSegmentMoverOP newMover );
	void add_compositeSegmentMover( RBSegmentMoverOP newMover );
	void add_wholeStructureMover( protocols::moves::MoverOP newMover );

	void clear_movesets();

	//void reorder_segments();

	void set_temperature( core::Real start, core::Real final );
	void set_helicalMoveStepsize( core::Real onAxisTrans, core::Real onAxisRot,
		core::Real offAxisTrans, core::Real offAxisRot );
	void set_genericRBMoveStepsize( core::Real trans, core::Real rot );
	void set_ncycles( int ncycles );

	void set_cst_weight( core::Real wt );
	void set_cst_width ( core::Real width );
	//void set_cst_type  ( std::string type );


	void set_randomize( core::Size rand ) { rand_ = rand; }
	void set_bootstrap( bool boot ) { bootstrap_ = boot; }
	void set_skip_lr( bool skip ) { no_lr_ = skip; }
	void set_fix_ligands( bool fix ) { fix_ligands_ = fix; }

private:
	void set_default_movemap();

	// scoring stuff
	core::Real cst_weight_, cst_width_, cst_stdev_;
	core::Size cst_seqwidth_;
	//std::string cst_type_;
	core::scoring::ScoreFunctionOP scorefxn_;
	protocols::moves::MonteCarloOP mc_;

	// protocol parameters -- need get and set methods
	core::Size cycles_;
	core::Real init_temp, final_temp; // MC temperature
	core::Real helical_sigT, helical_sigR; // std dev of helical axis rotations and translations
	core::Real helical_sigOffAxisT, helical_sigOffAxisR;  // std dev of helical off-axis rotations and translations
	core::Real strand_sigT, strand_sigR; // std dev of strand axis rotations and translations
	core::Real strand_sigOffAxisT, strand_sigOffAxisR;  // std dev of strand off-axis rotations and translations
	core::Real genericRB_sigT, genericRB_sigR; // stddev of generic RB rotations and translations
	core::Real randomness_;  // Random factor to add to fragment score when choosing RMS-frags

	core::Size rand_;   // number of randomize steps
	bool no_lr_;    // dont do loop relax
	bool bootstrap_;  // initially rebuild pose using fragment insertions only
	bool fix_ligands_;  // [experimental] dont let ligands move during minimize

	// the rigid-body segments
	utility::vector1< RBSegment > rbsegs_input_, rbsegs_remap_;

	// the loops
	protocols::loops::Loops loops_input_;

	// fragments
	utility::vector1< core::fragment::FragSetOP > frag_libs_;

	// sets of allowed moves for each conf. type
	std::vector< protocols::moves::MoverOP > WholeStructureMoveSet_;
	std::vector< RBSegmentMoverOP > HelixMoveSet_;
	std::vector< RBSegmentMoverOP > StrandMoveSet_;
	std::vector< RBSegmentMoverOP > GenericRBMoveSet_;
	std::vector< RBSegmentMoverOP > CompositeSegmentMoveSet_;
};


/*

////////////////////////////////////////////////////////////
/// Sim Anneal Mover
///fd  -- note: THis should probably be moved elsewhere
////////////////////////////////////////////////////////////

class SimAnnealMover : public protocols::moves::Mover {
public:
SimAnnealMover() {}

SimAnnealMover(
protocols::moves::MoverOP mover_in,
protocols::moves::MonteCarloOP mc_in,
core::Real init_temp,
core::Real final_temp,
core::Size cycles
) : moves::Mover("SimAnnealMover"),
mover_( mover_in ),
mc_( mc_in ),
init_temp_( init_temp ),
final_temp_( final_temp ),
cycles_( cycles )
{}


virtual void apply( core::pose::Pose & pose )
{
core::Real gamma = std::pow( ( float( final_temp_/ init_temp_ ) ),
float( 1.0f/cycles_ ) );
core::Real temperature( init_temp_ );
mc_->set_temperature( temperature );

for ( Size i = 1; i <= cycles_; ++i  ) {
mover_->apply( pose );
temperature *= gamma;
//std::cout << "Temperature " << temperature << std::endl;
mc_->set_temperature( temperature );
}
}

protected:
protocols::moves::MoverOP mover_;
protocols::moves::MonteCarloOP mc_;
// MoverStatistics stats_;


private:
core::Real init_temp_, final_temp_;
core::Size cycles_;
};

*/

}
}
#endif
