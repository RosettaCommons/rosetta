// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/Loops.hh
/// @brief
/// @author Chu Wang
/// @author Mike Tyka
/// @author James Thompson

#ifndef INCLUDED_protocols_loops_Loops_hh
#define INCLUDED_protocols_loops_Loops_hh

// Unit header
#include <protocols/loops/Loops.fwd.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <numeric/xyzVector.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/id/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Package headers
#include <protocols/loops/Loop.hh>

namespace protocols {
namespace loops {

///////////////////////////////////////////////////////////////////////////
// a list of loops
class Loops : public utility::pointer::ReferenceCount {

public:
	typedef utility::vector1< Loop > LoopList;
	typedef LoopList::iterator iterator;
	typedef LoopList::const_iterator const_iterator;


public:
	inline core::Size num_loop() const { return loops_.size(); }
	inline const_iterator begin() const { return loops_.begin(); }
	inline const_iterator end() const { return loops_.end(); }
	inline iterator v_begin() { return loops_.begin(); }
	inline iterator v_end() { return loops_.end(); }

	//constructor
	Loops(){};
	//copy constructor
	Loops( const Loops & src ):
		utility::pointer::ReferenceCount(),
		loops_(src.loops_) {}
	//operator
	Loops & operator =( Loops const & src ) {
		loops_ = src.loops_;
		return *this;
	}

	friend std::ostream & operator<<( std::ostream & os, const Loops & loops );

	void read_loop_file(
		std::string filename,
		bool strict_looprelax_checks = true,
		std::string token = "LOOP"
	);
	void read_stream_to_END(
    std::istream &is,
		bool strict_looprelax_checks,
		std::string filename, /*for error msg */
		std::string token = "LOOP"
	);
	void
	write_loops_to_file(
		std::string const & filename,
		std::string token = "LOOP"
	) const;

	void
	write_loops_to_stream(
  	std::ostream& data,
		std::string token
	) const;

	void
	add_loop( Loop loop, core::Size minimal_gap = 0 );

	void
	add_loop(
		core::Size const start,
		core::Size const stop,
		core::Size const cut = 0,
		core::Real skip_rate = 0.0,
		bool const extended = false
	);

	void
	add_loop(
		const Loops::const_iterator & it
	);

	void
	add_loop(
		const Loops::iterator & it
	);

	void push_back( Loop loop ) {
		add_loop( loop );
	}

	void
	push_back(
		core::Size const start,
		core::Size const stop,
		core::Size const cut = 0,
		core::Real skip_rate = 0.0,
		bool const extended = false
	) {
		add_loop( start, stop, cut, skip_rate, extended );
	}

	void
	add_overlap_loop(
		Loops loops
	);

	void
	add_overlap_loop(
		const Loop loop
	);

	void
	delete_loop(
		core::Size const start,
		core::Size const stop
	);

	const_iterator one_random_loop() const;

	core::Size
	loop_size(
		core::Size const loop_num
	) const;

	///@brief return number of residues in all loops of this definition -- sum_i( loop_size( i ) )
	core::Size
	loop_size() const;

	core::Size size() const {
		return loops_.size();
	}

	core::Size nr_residues() const {
		return loop_size();
	}

	void sequential_order();

	void clear();

	LoopList const& loops() const { return loops_; }

	/// @brief  Is seqpos contained in any of my loops?
	bool
	is_loop_residue( core::Size const seqpos, int const offset = 0 ) const;

	/// @brief is seqpos a residue in this Loops container ?
	bool
	has( core::Size const seqpos, int const offset = 0 ) const {
		return is_loop_residue( seqpos, offset );
	}

	void set_extended( bool input );

 	void make_sequence_shift( int shift );

	/// @brief yield the Loop which contains the residue seqpos, returns false if seqpos is not in any residue.
	bool loop_of_residue( core::Size const seqpos, Loop& loop ) const;

	///@brief yields the regions that are NOT loops in this loop-object.
	Loops invert( core::Size total_res ) const;

	//@brief switch DOF_Type for residues in loop. id::CHI, id::BB --- don't use with id::JUMP
	void switch_movemap( core::kinematics::MoveMap& movemap, core::id::TorsionType, bool allow_moves = true ) const;

	//@brief return index in list of the loop, 0 if not found
	core::Size loop_index_of_residue( core::Size const seqpos ) const;

	//@brief Autochoose a cutpoint using the secondary structure of the pose unless cutpoint is already set
	void auto_choose_cutpoints( core::pose::Pose const & pose );

	//@brief Autochoose a cutpoint using the secondary structure of the pose unless cutpoint is already set
	void choose_cutpoints( core::pose::Pose const & pose );

	void verify_against(	core::pose::Pose const & pose ) const;

	void remove_terminal_loops( core::pose::Pose const & pose );

	/// @brief Extend a loop
	void grow_all_loops(
		core::Size nres,
		core::Real magnitude
	);

	/// @brief Extend a loop (don't extend across cutpoints)
	void grow_all_loops(
		core::pose::Pose const &pose,
		core::Real magnitude
	);

	/// @brief Extend a loop (don't extend across cutpoints)
	void grow_loop(
		core::pose::Pose const &pose,
		Loop & loop,
		core::Real magnitude
	);

	/// @brief Extend a loop
	void grow_loop(
		core::Size nres,
		Loop & loop,
		core::Real magnitude
	);

	/// @brief Extend a loop unequally in both dirs
	void grow_loop(
		core::Size nres,
		Loop & loop,
		core::Real mag_left,
		core::Real mag_right
	);

  /// @brief Computes the center of mass of the Ca atoms specified by this
  /// instance, writing the result to <center>. Assumes there is no missing
  /// backbone density.
  ///
  /// Note: if this method is called on an instance without any Loop's, returns (0,0,0).
  void center_of_mass(const core::pose::Pose& pose, numeric::xyzVector<core::Real>* center) const;

	///@brief set each loop-residue in the vector to val.
	/// input vector of nres length ( if shorter last residues of loop are ignored )
	template< class T >
	void transfer_to_residue_vector( utility::vector1< T >&, T val );

	///@brief add all residues within this loop definition into selection
	void get_residues( utility::vector1< Size>& selection ) const;
	// i know this encourages old style for-loops (i.e. without iterators) but so much of the code
	// already used such loops, i opted for safety.
	inline
	const Loop&
	operator[] ( core::Size const i ) const
	{
		return loops_[i];
	}

	inline
	Loop&
	operator[] ( core::Size const i )
	{
		return loops_[i];
	}

	bool operator==( Loops const& other ) const;

	bool operator!=( Loops const& other ) const {
		return !( *this == other);
	}

private:
	LoopList loops_;
}; // Loops

Loops get_loops_from_file();
std::string get_loop_file_name();

template< class T >
void Loops::transfer_to_residue_vector( utility::vector1< T > & vector, T val ) {
	core::Size nres = vector.size();
	for ( const_iterator it = begin(); it  != end(); ++it ) {
		if ( it->start() <= nres ) {
			for ( core::Size pos = it->start(); pos <= std::min( it->stop(), nres ); pos++ ) {
				vector[ pos ] = val;
			}
		}
	}
}

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_Loops_HH
