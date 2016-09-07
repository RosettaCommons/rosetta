// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/matdes/GetRBDOFValues.hh
/// @brief  header file for GetRBDOFValues class
/// @author Jacob Bale (balej@u.washington.edu)


#ifndef INCLUDED_protocols_matdes_GetRBDOFValues_hh
#define INCLUDED_protocols_matdes_GetRBDOFValues_hh

// Unit Headers
#include <protocols/matdes/GetRBDOFValues.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


//// C++ headers

namespace protocols {
namespace matdes {

/// @brief WARNING WARNING WARNING THIS CLASS USES THE THREAD UNSAFE SymDofMoverSampler
/// AND MAKES ANY PROTOCOL THAT USES THIS CLASS THREAD UNSAFE.
class GetRBDOFValues : public protocols::filters::Filter {
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	// @brief default constructor
	GetRBDOFValues();

	// @brief constructor with arguments
	GetRBDOFValues( int jump, std::string dof_name, bool verb, char ax, bool disp, bool ang, core::Real init_d, core::Real init_a, bool get_init );

	// @brief copy constructor
	GetRBDOFValues( GetRBDOFValues const & rval );

	~GetRBDOFValues() override;


public:// virtual constructor


	// @brief make clone
	protocols::filters::FilterOP clone() const override;

	// @brief make fresh instance
	protocols::filters::FilterOP fresh_instance() const override;


public:// accessor

	// @brief get name of this filter
	std::string name() const override { return "GetRBDOFValues"; }

public:// setters

	void jump_id( core::Size const jump );
	void sym_dof_name( std::string const & dof_name );
	void verbose( bool const verb );
	void axis( char const ax );
	void radial_disp( bool const disp );
	void angle( bool const ang );
	void init_disp( core::Real const init_d );
	void init_angle( core::Real const init_a );
	void get_init_value( bool const get_init );

public:// getters
	core::Size jump_id() const;
	std::string sym_dof_name() const;
	bool verbose() const;
	char axis() const;
	bool radial_disp() const;
	bool angle() const;
	core::Real init_disp() const;
	core::Real init_angle() const;
	bool get_init_value() const;

public:// parser

	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		Movers_map const &,
		Pose const & ) override;

public:// virtual main operation

	// @brief always returns true.
	bool apply( Pose const & pose ) const override;

	/// @brief
	core::Real report_sm( Pose const & pose ) const override;
	void report( std::ostream & out, Pose const & pose ) const override;

	/// @brief get the translation and rotation for a user specified jump
	core::Real compute( Pose const & pose, bool const & verb, std::string const & dof_name, int const & jump, char const & ax, bool const & disp, bool const & ang, core::Real const & init_d, core::Real const & init_a, bool const & get_init ) const;

private:

	int jump_id_;
	std::string sym_dof_name_;
	bool verbose_, radial_disp_, angle_, get_init_value_;
	char axis_;
	core::Real init_disp_, init_angle_;

};

} // matdes
} // protocols

#endif
