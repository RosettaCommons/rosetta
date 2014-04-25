// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @details
/// @author Oliver Lange


// Unit Headers
#include <protocols/loops/loop_closure/ccd/CcdLoopClosureMover.hh>
#include <protocols/loops/loop_closure/ccd/ccd_closure.hh>
#include <core/kinematics/MoveMap.hh>
// Package Headers

// Project Headers
// AUTO-REMOVED #include <core/pose/Pose.hh>

// Utility Headers
#include <utility/PyAssert.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

//Auto Headers



//last boolean causes this tracer to be muted by default
static basic::Tracer TR_loop("protocols.loops.loop_closure.ccd.CcdLoopClosureMover", basic::t_info, true);

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

using namespace core;

CcdLoopClosureMover::CcdLoopClosureMover(
     Loop const& loop_def,
     core::kinematics::MoveMapCOP mm
) :
	loop_( loop_def ),
	movemap_( mm ),
	max_rama_score_increase_( 2.0 ), //defaults copied from rosetta++, map_squence.cc scored_frag_close
	max_total_delta_helix_( 10 ),
	max_total_delta_strand_( 50 ),
	max_total_delta_loop_( 75 ),
	tolerance_( 0.01 ),
	ccd_cycles_( 100 ),
	bRama_check_( true )
{}

void
CcdLoopClosureMover::set_max_rama_score_increase(
	Real input_max_rama_score_increase
) {
	max_rama_score_increase_ = input_max_rama_score_increase;
	return;
}

void
CcdLoopClosureMover::set_max_total_delta(
	Real input_max_delta_helix,
	Real input_max_delta_strand,
	Real input_max_delta_loop
) {
	max_total_delta_helix_ = input_max_delta_helix;
	max_total_delta_strand_ = input_max_delta_strand;
	max_total_delta_loop_ = input_max_delta_loop;
	return;
}

void
CcdLoopClosureMover::set_tolerance( Real input_tolerance ) {
	tolerance_ =  input_tolerance;
	return;
}

void
CcdLoopClosureMover::set_ccd_cycles( Size input_ccd_cycles ) {
	ccd_cycles_ = input_ccd_cycles;
	return;
}

void
CcdLoopClosureMover::set_bRama_check( bool input_bRama_check ) {
	bRama_check_ = input_bRama_check;
	return;
}

CcdLoopClosureMover::~CcdLoopClosureMover(){}

void
CcdLoopClosureMover::apply( pose::Pose &pose ) {
	TR_loop << "Closing loop " << loop_.start() << " " << loop_.stop() << std::endl;
  actual_cycles_ = fast_ccd_loop_closure(
    pose,
		*movemap_,
    loop_.start(),
    loop_.stop(),
    loop_.cut(),
    ccd_cycles_,
    tolerance_,
    bRama_check_,
    max_rama_score_increase_,
    max_total_delta_helix_,
    max_total_delta_strand_,
    max_total_delta_loop_,
    forward_deviation_,
    backward_deviation_,
    torsion_delta_,
    rama_delta_
  );
}

std::string
CcdLoopClosureMover::get_name() const {
	return "CcdLoopClosureMover";
}

void
CcdLoopClosureMover::show(std::ostream & output) const
{
	Mover::show(output);
	output << "LOOP start: " << get_loop_start() << "  stop: " << get_loop_stop()<< "  cut: " << get_loop_cut() <<
			"  size: " << get_loop_size() << "  skip rate: " << get_loop_skip_rate() << "  extended?: " <<
			get_loop_extended() <<  std::endl;
	output <<   "Max rama score increase: " << get_max_rama_score_increase() <<
				"\nNumber of CCD cycles:    " << get_ccd_cycles() <<
				"\nTolerance:               " << get_tolerance() <<
				"\nbRama check:             " << bRama_check() <<
				"\nMax total delta helix:   " << get_max_total_delta("H") <<
				"\nMax total delta strand:  " << get_max_total_delta("E") <<
				"\nMax total delta loop:    " << get_max_total_delta("L") << "\nMovemap: " << std::endl;
	movemap().show(output);
}

core::Real
CcdLoopClosureMover::get_max_rama_score_increase() const {
	return max_rama_score_increase_;
}

core::Real
CcdLoopClosureMover::get_tolerance() const {
	return tolerance_;
}

core::Real
CcdLoopClosureMover::get_max_total_delta( std::string secstr ) const {
	assert((secstr == "H") || (secstr == "E") || (secstr == "L"));
	PyAssert((secstr == "H") || (secstr == "E") || (secstr == "L"),
			"CcdLoopClosureMover::get_max_total_delta( std::string secstr ): secstr must be 'H', 'E', or 'L'.");
	if (secstr == "H")
		return max_total_delta_helix_;
	else if (secstr == "E")
		return max_total_delta_strand_;
	else /*secstr == "L"*/
		return max_total_delta_loop_;
}

core::kinematics::MoveMap
CcdLoopClosureMover::movemap() const { 
	return *movemap_; 
}

std::ostream &operator<< ( std::ostream &os, CcdLoopClosureMover const &mover )
{
	mover.show(os);
	return os;
}

CcdMover::CcdMover(
	Loop const& loop_def,
	core::kinematics::MoveMapCOP mm
) :
	total_moves_( loop_def.size() ),
  loop_( loop_def ),
  movemap_( mm )
{}

CcdMover::~CcdMover() {}

void
CcdMover::apply( pose::Pose &pose ) {
  ccd_moves(
		total_moves_,
		pose,
		*movemap_,
		loop_.start(),
		loop_.stop(),
		loop_.cut()
  );
}

std::string
CcdMover::get_name() const {
	return "CcdMover";
}

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols
