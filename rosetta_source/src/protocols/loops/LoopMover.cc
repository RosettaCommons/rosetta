// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/IndependentLoopMover.cc
/// @brief  loop mover base class
/// @author Mike Tyka

#include <protocols/loops/LoopMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh> // tracer output

#include <core/fragment/FragSet.hh>

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
// AUTO-REMOVED #include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end





namespace protocols {
namespace loops {

///////////////////////////////////////////////////////////////////////////////
using namespace core;
using namespace ObjexxFCL;
using namespace ObjexxFCL::fmt;

basic::Tracer tr("protocol.loops.LoopMover");

LoopMover::LoopMover() : Mover()
{
   init( new Loops() ); 
}

LoopMover::LoopMover( protocols::loops::LoopsOP loops_in ) : Mover()
{
    init( loops_in );
}
    
void LoopMover::init( protocols::loops::LoopsOP loops_in )
{
    Mover::type( "LoopMover" );
    loops_from_observer_cache_ = false;
    checkpoints_ = new checkpoint::CheckPointer( "LoopMover" );
    false_movemap_ = new core::kinematics::MoveMap();
    loops_ = loops_in;
    
}

// destructor
LoopMover::~LoopMover(){}

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

void LoopMover::loops( protocols::loops::LoopsOP const l )
{ 
    loops_ = l; 
}

const protocols::loops::LoopsOP LoopMover::loops() const 
{
    return loops_;
}


const utility::vector1< core::fragment::FragSetOP > & LoopMover::frag_libs() const 
{
    return frag_libs_;
}

void LoopMover::set_use_loops_from_observer_cache( bool const loops_from_observer_cache )
{
    loops_from_observer_cache_ = loops_from_observer_cache;
}

bool const LoopMover::use_loops_from_observer_cache() const
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

	tr.Debug << "Extending loop torsions" << loop.start() << " " << loop.stop()
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
		loops_->clear();
		utility::vector1< std::pair< core::Size, core::Size > > const & segments = utility::pointer::static_pointer_cast< core::pose::datacache::SpecialSegmentsObserver const >(pose.observer_cache().get_const_ptr( core::pose::datacache::CacheableObserverType::SPECIAL_SEGMENTS_OBSERVER ) )->segments();
		for( core::Size i = 1; i <= segments.size(); ++i ){
			core::Size loop_end = segments[i].second - 1; //segment convention
			if( loop_end <= segments[i].first ) continue; //safeguard against faulty or segments of length 1
			tr << "Setting loop from observer cache between seqpos " << segments[i].first << " and " << loop_end << "." << std::endl;
			loops_->add_loop( segments[i].first, loop_end, numeric::random::random_range( int(segments[i].first), int(loop_end)  ) );
		}
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
	lhs.loops_ = rhs.loops_;
    lhs.scorefxn_ = rhs.scorefxn_;
    lhs.frag_libs_ = rhs.frag_libs_;
    lhs.checkpoints_ = rhs.checkpoints_;
    lhs.loops_from_observer_cache_ = rhs.loops_from_observer_cache_;
    lhs.false_movemap_ = rhs.false_movemap_;
}

} // namespace loops
} // namespace protocols
