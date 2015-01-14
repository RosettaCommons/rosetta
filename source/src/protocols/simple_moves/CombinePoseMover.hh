// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/CombinePoseMover.hh
/// @brief
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_simple_moves_CombinePoseMover_hh
#define INCLUDED_protocols_simple_moves_CombinePoseMover_hh

#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/simple_moves/CombinePoseMover.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/io/silent/SilentStruct.hh>
//#include <core/io/silent/SilentStruct.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {

using namespace core;

class CombinePoseMover : public protocols::moves::Mover

{
public:
	//constructor
  CombinePoseMover( core::scoring::ScoreFunctionCOP sfxn,
										pose::Pose const &pose_ref );

  virtual ~CombinePoseMover();

	virtual 
	void apply( pose::Pose &pose2 );

	virtual 
	std::string get_name() const { return "CombinePoseMover"; }

	void set_default();

	void set_max_struct( Size const n ){ max_struct_ = n; }
	void set_max_try( Size const n ){ max_struct_try_ = n; }
	void set_minfrac_crossover( Real const value ){ minfrac_crossover_ = value; }
	void set_maxfrac_crossover( Real const value ){ maxfrac_crossover_ = value; }
	void set_nonideal( bool const setting ){ nonideal_ = setting; }
	void set_minimize( bool const setting ){ do_minimize_ = setting; }
	void set_cartesian( bool const setting ){ cartesian_crossover_ = setting; }
	void set_rmsdcut( Real const setting ){ rmsdcut_ = setting; }

	Size max_struct() const { return max_struct_; }
	Size max_struct_try() const { return max_struct_try_; }
 
	pose::Pose const pose_ref() const { return pose_ref_; }

	std::string pose_tag() const { return pose_tag_; }
	void set_pose_tag( std::string const value ) { pose_tag_ = value; }

	void set_store_silents( bool const value ){ store_silents_ = value; }

	// Structures / history
	std::vector< io::silent::SilentStructOP >
	return_silent() const { return sampled_structures_; }
	
	void clear_combine_history() { combine_history_.resize( 0 ); }
	void append_combine_history( std::vector< Size > const v )
	{ combine_history_.push_back( v ); }

	std::vector< std::vector< Size > > 
	return_combine_history() const { return combine_history_; }

private:
	pose::Pose pose_ref_;

	// Parameters/options
	Size max_struct_;
	Size max_struct_try_;
	Real minfrac_crossover_;
	Real maxfrac_crossover_;
	bool cartesian_crossover_;
	bool nonideal_;
	bool store_silents_;
	bool do_minimize_;

	// Tag
	std::string pose_tag_;

	// Minimization score
	core::scoring::ScoreFunctionCOP sfxn_;

	// Simple filter for clash
	core::scoring::ScoreFunctionCOP sfxn0_;
	Real vdwcut_;
	Real rmsdcut_;

	// structure storage
	std::vector< std::vector< Size > > combine_history_;
	std::vector< io::silent::SilentStructOP > sampled_structures_;

};

} // namespace simple_moves
} // namespace protocols

#endif
