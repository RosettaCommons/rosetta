// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/loop_mover/LoopMover.cc
/// @brief  loop mover base class
/// @author Mike Tyka

// Unit headers
#include <protocols/loops/loop_mover/LoopMover.hh>

// Package Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>
#include <protocols/loops/loops_main.hh>

#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh> // tracer output

#include <core/fragment/FragSet.hh>
#include <core/conformation/ppo_torsion_bin.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/kinematics/MoveMap.hh>
#include <protocols/checkpoint/CheckPointer.hh>

#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

//Utility Headers


/// ObjexxFCL headers
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <iostream>
#include <map>
#include <string>
#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif
// option key includes
// AUTO-REMOVED #include <basic/options/keys/loops.OptionKeys.gen.hh>

// AUTO-REMOVED #include <basic/options/option.hh>
#include <numeric/random/random.hh>

//Auto Headers

// AUTO-REMOVED #include <ObjexxFCL/format.hh>

// Auto-header: duplicate removed #include <utility/vector1.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end





namespace protocols {
namespace loops {
namespace loop_mover {

///////////////////////////////////////////////////////////////////////////////
using namespace core;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

LoopMover::LoopMover() :
	Mover(),
	guarded_loops_( GuardedLoopsFromFileOP( new GuardedLoopsFromFile ) )
{
   init();
}

/// @details GuardedLoops constructed with a pointer to a loops object:
/// sets the GuardedLoops object into it's "not in charge" mode; i.e.
/// the external source of the LoopsOP is responsible for resolving loop indices
LoopMover::LoopMover( protocols::loops::LoopsOP loops_in )
:
	Mover(),
	guarded_loops_( GuardedLoopsFromFileOP( new GuardedLoopsFromFile( loops_in ) ))
{
    init();
}

/// @details GuardedLoops constructed with an unresolve-loop-indices object:
/// sets the GuardedLoops object into it's "in charge" mode; i.e.
/// before the loops object may be used, the loop indices must be resolved.
LoopMover::LoopMover( protocols::loops::LoopsFileData const & loops_from_file ) :
	Mover(),
	guarded_loops_( GuardedLoopsFromFileOP( new GuardedLoopsFromFile( loops_from_file ) ))
{
	init();
}

LoopMover::LoopMover( protocols::loops::GuardedLoopsFromFileOP guarded_loops ) :
	Mover(),
	guarded_loops_( guarded_loops )
{
	init();
}


void LoopMover::init()
{
	Mover::type( "LoopMover" );
	loops_from_observer_cache_ = false;
	checkpoints_ = checkpoint::CheckPointerOP( new checkpoint::CheckPointer( "LoopMover" ) );
	false_movemap_ = MoveMapOP( new core::kinematics::MoveMap() );
}

// destructor
LoopMover::~LoopMover(){}

void
LoopMover::set_guarded_loops_not_in_charge()
{
	guarded_loops_->in_charge( false );
}


//////////////////////////////////////////////////////////////////////////////////

// accessors //

void LoopMover::set_scorefxn( const core::scoring::ScoreFunctionOP score_in )
{
	scorefxn_ = score_in;
}

const core::scoring::ScoreFunctionOP & LoopMover::scorefxn() const
{
	return scorefxn_;
}

void LoopMover::loops( protocols::loops::LoopsOP lptr )
{
	// pointer assignment
	guarded_loops_->set_loops_pointer( lptr );
}

void LoopMover::loops( protocols::loops::LoopsFileData const & loops )
{
	// pointer assignment
	guarded_loops_->loops( loops );
}

/// @brief Set the guarded_loops pointer
void LoopMover::loops( protocols::loops::GuardedLoopsFromFileOP guarded_loops )
{
	guarded_loops_ = guarded_loops;
}

protocols::loops::LoopsCOP LoopMover::loops() const
{
	return guarded_loops_->loops();
}

protocols::loops::LoopsOP LoopMover::loops()
{
   return guarded_loops_->loops();
}



const utility::vector1< core::fragment::FragSetOP > & LoopMover::frag_libs() const
{
    return frag_libs_;
}

void LoopMover::set_use_loops_from_observer_cache( bool const loops_from_observer_cache )
{
    loops_from_observer_cache_ = loops_from_observer_cache;
}

bool LoopMover::use_loops_from_observer_cache() const
{
    return loops_from_observer_cache_;
}

checkpoint::CheckPointerOP & LoopMover::get_checkpoints()
{
    return checkpoints_;
}

void LoopMover::false_movemap( MoveMapOP const & mm )
{
    false_movemap_ = mm;
}

LoopMover::MoveMapOP const & LoopMover::false_movemap() const
{
    return false_movemap_;
}

// end of accessors

// Additional MoveMap method
Size LoopMover::enforce_false_movemap( MoveMapOP & mm ) const
{
    return mm->import_false( *false_movemap() );
}


/// @details  Set a loop to extended torsion angles.
void LoopMover::set_extended_torsions(
	core::pose::Pose & pose,
	Loop const & loop
)
{
	Real const init_phi  ( -150.0 );
	Real const init_psi  (  150.0 );
	Real const init_omega(  180.0 );

	static int counter = 0;

	tr().Debug << "Extending loop torsions" << loop.start() << " " << loop.stop()
		<< std::endl;

	idealize_loop(pose, loop );

	Size start_extended = std::max((Size)1,loop.start());
	Size end_extended   = std::min(pose.total_residue(),loop.stop());
	for ( Size i = start_extended; i <= end_extended; ++i ) {
		if ( i != start_extended )	pose.set_phi( i, init_phi );
		if ( i != end_extended ) pose.set_psi( i, init_psi );
		if ( ( i != start_extended ) && ( i != end_extended ) ) pose.set_omega( i, init_omega );
	}

	counter++;
}

std::string
LoopMover::get_name() const {
	return "LoopMover";
}

void LoopMover::add_fragments(
	core::fragment::FragSetOP fragset
) {
	if ( fragset->size() > 0 ) {
		frag_libs_.push_back( fragset->clone() );
	}
}

void LoopMover::clear_fragments() {
	frag_libs_.clear();
}

void
LoopMover::set_loops_from_pose_observer_cache( core::pose::Pose const & pose ){

	if( pose.observer_cache().has( core::pose::datacache::CacheableObserverType::SPECIAL_SEGMENTS_OBSERVER) ){
		LoopsOP loops( new Loops() );
		utility::vector1< std::pair< core::Size, core::Size > > const & segments = utility::pointer::static_pointer_cast< core::pose::datacache::SpecialSegmentsObserver const >(pose.observer_cache().get_const_ptr( core::pose::datacache::CacheableObserverType::SPECIAL_SEGMENTS_OBSERVER ) )->segments();
		for( core::Size i = 1; i <= segments.size(); ++i ){
			core::Size loop_end = segments[i].second - 1; //segment convention
			if( loop_end <= segments[i].first ) continue; //safeguard against faulty or segments of length 1
			tr() << "Setting loop from observer cache between seqpos " << segments[i].first << " and " << loop_end << "." << std::endl;
			loops->add_loop( segments[i].first, loop_end, numeric::random::random_range( int(segments[i].first), int(loop_end)  ) );
		}
		guarded_loops_->loops( *loops ); // <--- deep copy into the GuardedLoops existing Loops object
	}
	else{
		utility_exit_with_message("trying to set loops from observer cache even though no cache was detected in the pose");
	}
}


// Used by both LoopMover_Refine_KIC and LoopMover_Perturb_KIC - maybe it should be there, then?
void
loops_set_chainbreak_weight( core::scoring::ScoreFunctionOP scorefxn, Size const round ){
	bool const use_linear_chainbreak( basic::options::option[basic::options::OptionKeys::loops::kic_use_linear_chainbreak]() );
	if ( use_linear_chainbreak ){
		scorefxn->set_weight( core::scoring::linear_chainbreak, float( round ) * 150.0 / 3.0 );
		scorefxn->set_weight( core::scoring::chainbreak, 0.0 );
	} else { // default behavior.
		scorefxn->set_weight( core::scoring::chainbreak, float(round) * 10.0 / 3.0 );
	}
}

///@brief copy ctor
LoopMover::LoopMover( LoopMover const & rhs ) :
	Mover(rhs)
{
	initForEqualOperatorAndCopyConstructor(*this, rhs);
}

///@brief assignment operator
LoopMover & LoopMover::operator=( LoopMover const & rhs ){
	//abort self-assignment
	if (this == &rhs) return *this;

	Mover::operator=(rhs);
	initForEqualOperatorAndCopyConstructor(*this, rhs);
	return *this;
}

void LoopMover::initForEqualOperatorAndCopyConstructor(LoopMover & lhs, LoopMover const & rhs){
	//going through all of member data and assigning it
	lhs.guarded_loops_ = rhs.guarded_loops_; // shallow copy of the guarded_loops_ pointer.
	lhs.scorefxn_ = rhs.scorefxn_;
	lhs.frag_libs_ = rhs.frag_libs_;
	lhs.checkpoints_ = rhs.checkpoints_;
	lhs.loops_from_observer_cache_ = rhs.loops_from_observer_cache_;
	lhs.false_movemap_ = rhs.false_movemap_;
}

void LoopMover::resolve_loop_indices( core::pose::Pose const & p )
{
	guarded_loops_->resolve_loop_indices( p );
}

/// @brief bin torsion angles as described in http://www.ncbi.nlm.nih.gov/pubmed/19646450
/// @details generates a string with the torsion angle bins, using uppercase letters as in the
/// publication above for omega ~ 180, and lowercase letters for omega ~ 0; to be used in
/// loop sampling analysis
///@author Amelie Stein
///@date April 26, 2012
core::conformation::torsion_bin_string
LoopMover::torsion_features_string( core::pose::Pose const & pose ) const
{
	core::conformation::torsion_bin_string torsion_bins, pos_bin;

	// currently this generates one string for all loops, so the user has to map the positions --
	// alternatively one could return a vector of torsion feature strings or some other more complex construct
	// note: won't work properly for multiple loops at the moment -- TODO
	for ( Loops::const_iterator it=loops()->begin(), it_end=loops()->end(); it != it_end; ++it ) {
		if (torsion_bins.size() != 0 ) {
			torsion_bins.push_back( core::conformation::ppo_torbin_U ); // next loop
		}
		for ( core::Size i = it->start(), l_end = it->stop(); i <= l_end; i++ ) {
			torsion_bins.push_back( core::conformation::get_torsion_bin(pose.phi(i), pose.psi(i), pose.omega(i) ));
			//std::cerr << i << " " << pose.phi(i) << " " << pose.psi(i) << " " <<  pose.omega(i) << " " << torsion_bins << std::endl; // debugging
		}
	}

	return torsion_bins;

} // torsion_features_string

} // namespace loop_mover
} // namespace loops
} // namespace protocols
