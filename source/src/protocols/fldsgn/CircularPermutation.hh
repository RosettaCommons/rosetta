// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fldsgn/CircularPermutation.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_CircularPermutation_hh
#define INCLUDED_protocols_fldsgn_CircularPermutation_hh

// unit headers
#include <protocols/fldsgn/CircularPermutation.fwd.hh>

// type headers
#include <core/types.hh>

// package headers

// project headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>


// C++ headers
#include <string>

// parser headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>


#include <utility/vector1.hh>

namespace protocols {
namespace fldsgn {


class CircularPermutation : public protocols::moves::Mover {


private: // typedefs


	typedef protocols::moves::Mover Super;


public: // typedefs


	//typedef std::string String;
	//typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef core::kinematics::MoveMap MoveMap;
	typedef core::kinematics::MoveMapOP MoveMapOP;
	typedef protocols::moves::MoverOP MoverOP;

	// for parser
	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public: // construct/destruct


	/// @brief default constructor
	CircularPermutation();

	/// @brief copy constructor
	CircularPermutation( CircularPermutation const & rval );

	/// @brief default destructor
	virtual ~CircularPermutation();


private: // disallow assignment


	/// @brief copy assignment
	/// @remarks Mover base class prevents this from working properly...
	CircularPermutation & operator =( CircularPermutation const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	MoverOP clone() const;

	/// @brief create this type of object
	virtual
	MoverOP fresh_instance() const;


public: // accessors


	/// @brief
	Size new_terminal_pos() const;


public: // mutators


	/// @brief
	void new_terminal_pos( Size const s );


public: // helper functions


	/// @brief
	Size which_chain( Size const s, Pose const & pose ) const;

	/// @brief
	void split_chains( Pose & pose, utility::vector1< Size > const & pos );


public: // virtual main methods


	/// @brief apply defined moves to given Pose
	virtual
	void apply( Pose & pose );

	virtual std::string get_name() const;

public: //parser


	/// @brief parse xml file
	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		Movers_map const &,
		Pose const & );


private: // data


	/// @brief new N- & C- terminal position
	Size new_terminal_pos_;

	/// @brief
	bool ignore_chain_;

	Size split_;


};


} // namespace fldsgn
} // namespace protocols


#endif /* INCLUDED_protocols_fldsgn_CircularPermutation_HH */
