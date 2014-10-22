// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/ShortBackrubMover.hh
/// @brief definition of ShortBackrubMover class and functions
/// @author Noah Ollikainen (nollikai@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_ShortBackrubMover_hh
#define INCLUDED_protocols_simple_moves_ShortBackrubMover_hh

// Unit headers
#include <protocols/simple_moves/ShortBackrubMover.fwd.hh>

// Project headers
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/backrub/BackrubMover.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

class ShortBackrubMover : public protocols::moves::Mover {

public:

	typedef protocols::moves::MoverOP MoverOP;
	
	// default constructor
	ShortBackrubMover();

	/// @brief constructor that sets input pose
	ShortBackrubMover( core::pose::PoseOP pose );
	
	/// @brief copy constructor
	ShortBackrubMover( ShortBackrubMover const & rval );

	/// @brief destructor
	virtual ~ShortBackrubMover();

	/// @brief clone this object
	virtual	protocols::moves::MoverOP clone() const;

	/// @brief create this type of object
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	
	// setters
	void set_resnum( core::Size resnum );
	void set_rotation_std_dev( core::Real rotation_std_dev );
	void set_randomize_resnum( bool randomize_resnum );
	void set_uniform_backrub( bool uniform_backrub );
	virtual void set_input_pose( core::pose::PoseCOP pose );
	
	// getters
	core::Size get_resnum() const;
	core::Real get_rotation_std_dev() const;
	bool get_randomize_resnum() const;
	bool get_uniform_backrub() const;
	protocols::backrub::BackrubMoverOP get_backrubmover() const;

	void parse_my_tag(
  	TagCOP tag,
  	basic::datacache::DataMap & data,
	  Filters_map const &,
	  protocols::moves::Movers_map const &,
	  Pose const & );

private:

	/// @brief mover used for Backrub moves
	protocols::backrub::BackrubMoverOP backrubmover_;
	
	/// @brief residue number specifying the center residue for the next Backrub move
	core::Size resnum_;
	
	/// @brief standard deviation of rotation angle (degrees) used for Backrub moves
	core::Real rotation_std_dev_;
	
	/// @brief if true, choose a random residue for the next move
	bool randomize_resnum_;
	
	/// @brief if true, sample rotation angle from a uniform distribution from -20 to 20
	bool uniform_backrub_;
	
};  // class ShortBackrubMover

} // moves
} // protocols

#endif
