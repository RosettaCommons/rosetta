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

#ifndef INCLUDED_protocols_loops_Loops_HH
#define INCLUDED_protocols_loops_Loops_HH

// Unit header
#include <protocols/loops/Loops.fwd.hh>

// Package headers
#ifdef WIN32
#include <protocols/loops/Loop.hh>
#else
#include <protocols/loops/Loop.fwd.hh>
#endif

#include <protocols/loops/LoopsFileIO.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <numeric/xyzVector.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

// C/C++ headers
#include <ostream>
#include <string>

namespace protocols {
namespace loops {

///////////////////////////////////////////////////////////////////////////
// a list of loops
class Loops : public utility::pointer::ReferenceCount {

public:
  typedef utility::vector1< Loop > LoopList;
  typedef utility::vector1< SerializedLoop > SerializedLoopList;
  typedef LoopList::iterator iterator;
  typedef LoopList::const_iterator const_iterator;


public:
  bool empty() const;
  core::Size num_loop() const;
  const_iterator begin() const;
  const_iterator end() const;
  iterator v_begin();
  iterator v_end();

  //constructor
  Loops();

    //copy constructor
  Loops( const Loops & src );

  Loops( SerializedLoopList const & src );
  Loops( std::string const & loop_file_name );
  Loops( bool setup_loops_from_options_system );

  Loops( utility::vector1< bool > const& selection,
         bool randomize_cutpoint = true );

  // assignment operator
  Loops & operator =( Loops const & src );

  // destructor
  virtual ~Loops();

  friend std::ostream & operator<<( std::ostream & os, const Loops & loops );

    void read_loops_options();

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
    const const_iterator & it
  );

  void
  add_loop(
    const iterator & it
  );

  void push_back( Loop loop );

  void
  push_back(
    core::Size const start,
    core::Size const stop,
    core::Size const cut = 0,
    core::Real skip_rate = 0.0,
    bool const extended = false
  );

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
  core::Size loop_size() const;

  core::Size size() const;

  core::Size nr_residues() const;

  void sequential_order();

  void clear();

  LoopList const & loops() const;

  LoopsFileIOOP get_loop_file_reader() const;

  /// @brief  Is seqpos contained in any of my loops?
  bool
  is_loop_residue( core::Size const seqpos, int const offset = 0 ) const;

  /// @brief is seqpos a residue in this Loops container ?
  bool has( core::Size const seqpos, int const offset = 0 ) const;

  void set_extended( bool input );

   void make_sequence_shift( int shift );

  /// @brief yield the Loop which contains the residue seqpos, returns false if seqpos is not in any residue.
  bool loop_of_residue( core::Size const seqpos, Loop& loop ) const;

  /// Given the total number of residues, returns the inverse of this selection.
  Loops invert(core::Size num_residues) const;

  //@brief switch DOF_Type for residues in loop. id::CHI, id::BB --- don't use with id::JUMP
  void switch_movemap( core::kinematics::MoveMap& movemap, core::id::TorsionType, bool allow_moves = true ) const;

  //@brief return index in list of the loop, 0 if not found
  core::Size loop_index_of_residue( core::Size const seqpos ) const;

  //@brief Autochoose a cutpoint using the secondary structure of the pose unless cutpoint is already set
  void auto_choose_cutpoints( core::pose::Pose const & pose );

  //@brief Autochoose a cutpoint using the secondary structure of the pose unless cutpoint is already set
  void choose_cutpoints( core::pose::Pose const & pose );

  void verify_against(  core::pose::Pose const & pose ) const;

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

  /// @brief if possible grows loop will not cross cutpoints or if possible into sheets
  void grow_loop_away_from_sheets(
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

    std::string const & loop_file_name();
    void set_loop_file_name_and_reset( std::string const & loop_filename );

    bool strict_looprelax_checks();
    void set_strict_looprelax_checks( bool const check );

    std::string const & file_reading_token();
    void set_file_reading_token( std::string const & token );

    /// @brief Computes the center of mass of the Ca atoms specified by this
    /// instance, writing the result to <center>. Assumes there is no missing
    /// backbone density.
    ///
    /// Note: if this method is called on an instance without any Loop's, returns (0,0,0).
    void center_of_mass(const core::pose::Pose& pose, numeric::xyzVector<core::Real>* center) const;

  ///@brief set each loop-residue in the vector to val.
  /// input vector of nres length ( if shorter last residues of loop are ignored )
  template< class T >
  void transfer_to_residue_vector( utility::vector1< T >&, T val ) const;

  ///@brief add all residues within this loop definition into selection
  void get_residues( utility::vector1< Size>& selection ) const;

  // i know this encourages old style for-loops (i.e. without iterators) but so much of the code
  // already used such loops, i opted for safety.
  const Loop & operator[] ( core::Size const i ) const;

  Loop & operator[] ( core::Size const i );

  bool operator==( Loops const& other ) const;

  bool operator!=( Loops const& other ) const;

private:
  void init(
    LoopList const & loops_in,
    bool const read_loop_file_from_options = false,
    std::string const & passed_in_filename = ""
  );

  LoopList setup_loops_from_data( SerializedLoopList const & loop_data );

  void read_loop_file();

private:

  LoopList loops_;
  mutable LoopsFileIOOP loop_file_reader_;

  bool strict_looprelax_checks_on_file_reads_;

  std::string loop_filename_;
  std::string file_reading_token_;

}; // Loops

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_Loops_HH
