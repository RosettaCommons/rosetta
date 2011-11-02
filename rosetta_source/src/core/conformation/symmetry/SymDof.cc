// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief  Class to store ingformation about symmetrical dofs
/// @file   core/conformation/symmetry/SymDof.cc
/// @author Ingemar Andre

// Unit headers
#include <core/conformation/symmetry/SymDof.hh>

// Utility header
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <utility/string_util.hh>

#include <utility/vector1.hh>


namespace core {
namespace conformation {
namespace symmetry {

SymDof::SymDof()
{
	for ( Size i=1; i<=6; ++i ) {
    allowed_dof_jumps_.push_back(false);
		lower_range_dof_jumps1_.push_back(0.0);
		upper_range_dof_jumps1_.push_back(0.0);
		lower_range_dof_jumps2_.push_back(0.0);
    upper_range_dof_jumps2_.push_back(0.0);
		has_range1_lower_.push_back(false);
	  has_range1_upper_.push_back(false);
	  has_range2_lower_.push_back(false);
	  has_range2_upper_.push_back(false);
		jump_dir_.push_back( 1 );
  }

}

SymDof::SymDof( SymDof const & src )
{
	allowed_dof_jumps_ = src.allowed_dof_jumps_;
	lower_range_dof_jumps1_ = src.lower_range_dof_jumps1_;
	upper_range_dof_jumps1_ = src.upper_range_dof_jumps1_;
	lower_range_dof_jumps2_ = src.lower_range_dof_jumps2_;
  upper_range_dof_jumps2_ = src.upper_range_dof_jumps2_;
	has_range1_lower_ = src.has_range1_lower_;
  has_range1_upper_ = src.has_range1_upper_;
  has_range2_lower_ = src.has_range2_lower_;
  has_range2_upper_ = src.has_range2_upper_;
	jump_dir_ = src.jump_dir_;
}

	SymDof &
  SymDof::operator=( SymDof const & src ) {
		allowed_dof_jumps_ = src.allowed_dof_jumps_;
		lower_range_dof_jumps1_ = src.lower_range_dof_jumps1_;
		upper_range_dof_jumps1_ = src.upper_range_dof_jumps1_;
		lower_range_dof_jumps2_ = src.lower_range_dof_jumps2_;
		upper_range_dof_jumps2_ = src.upper_range_dof_jumps2_;
		has_range1_lower_ = src.has_range1_lower_;
		has_range1_upper_ = src.has_range1_upper_;
		has_range2_lower_ = src.has_range2_lower_;
		has_range2_upper_ = src.has_range2_upper_;
		jump_dir_ = src.jump_dir_;
	return *this;
	}

	SymDof::~SymDof() {}

	// @details is df allowed to move?
	bool
	SymDof::allow_dof( int df ) const
	{
		assert( df >= X_DOF && df <= Z_ANGLE_DOF );
		assert( allowed_dof_jumps_.size() == 6 );
		return allowed_dof_jumps_[df];

	}

	bool
	SymDof::has_dof()
	{
		for (int i = X_DOF; i <= Z_ANGLE_DOF; ++i ) {
			if ( allow_dof(i) ) return true;
		}
		return false;
	}

	// @details the lower boundary of range1
	core::Real
	SymDof::range1_lower( int df ) const
	{
		assert( df >= X_DOF && df <= Z_ANGLE_DOF );
		assert( lower_range_dof_jumps1_.size() == 6 );
		return lower_range_dof_jumps1_[df];
	}

	// @details the upper boundary of range1
	core::Real
  SymDof::range1_upper( int df ) const
  {
    assert( df >= X_DOF && df <= Z_ANGLE_DOF );
    assert( upper_range_dof_jumps1_.size() == 6 );
    return upper_range_dof_jumps1_[df];
  }

	// @details the lower boundary of range2
	core::Real
	SymDof::range2_lower( int df ) const
	{
		assert( df >= X_DOF && df <= Z_ANGLE_DOF );
		assert( lower_range_dof_jumps2_.size() == 6 );
		return lower_range_dof_jumps2_[df];
	}

	// @details the upper boundary of range1
	core::Real
  SymDof::range2_upper( int df ) const
  {
    assert( df >= X_DOF && df <= Z_ANGLE_DOF );
    assert( upper_range_dof_jumps2_.size() == 6 );
    return upper_range_dof_jumps2_[df];
  }

	// details Have a range1 been specified?
	bool
	SymDof::has_range1( int df ) const
	{
		if ( has_range1_lower_[df]
					&& has_range1_upper_[df] ) return true;
		else return false;
	}

	// details Have a range2 been specified?
	bool
	SymDof::has_range2( int df ) const
	{
		if ( has_range2_lower_[df]
          && has_range2_upper_[df] ) return true;
    else return false;
	}

	// @details has a lower boundary of range1 been specified?
	bool
	SymDof::has_range1_lower( int df ) const
	{
		assert( df >= X_DOF && df <= Z_ANGLE_DOF );
		return has_range1_lower_[df];
	}

	// @details has a upper boundary of range1 been specified?
	bool
  SymDof::has_range1_upper( int df ) const
  {
		assert( df >= X_DOF && df <= Z_ANGLE_DOF );
    return has_range1_upper_[df];
  }

	// @details has a lower boundary of range2 been specified?
	bool
  SymDof::has_range2_lower( int df ) const
  {
		assert( df >= X_DOF && df <= Z_ANGLE_DOF );
    return has_range2_lower_[df];
  }

	// @details has a upper boundary of range2 been specified?
	bool
  SymDof::has_range2_upper( int df ) const
  {
		assert( df >= X_DOF && df <= Z_ANGLE_DOF );
    return has_range2_upper_[df];
  }

	// @detail return the direction( upstream or downstream )
	// of the jump for a dof
	int
	SymDof::jump_direction( int df ) const
	{
		assert( df >= X_DOF && df <= Z_ANGLE_DOF );
		return jump_dir_[df];
	}

// @details function to parse a string describing dofs and parameters associated with them
// This function is used in reading symmetry_definition files. A typical line would look like this:
// set_dof BASEJUMP x(50) angle_x angle_y angle_z
// x, y, z are translations along the cartesian axis. angle_x, angle_y, angle_z are rotations around
// the cartesian exis.
// There are two ranges that can be specified enclosed by parenthesises ie. x(0-50:2-3)
void SymDof::read( std::string dof_line )
{
	// replace parenthesis with space for easier parsing
	std::replace( dof_line.begin(), dof_line.end(), ')', ' ' );
	std::istringstream l( dof_line );
	while ( true ) {
		std::string j;
		l >> j;
		if ( l.fail() ) break;
		// first read dof_type
		int dof_type(0);
		//  Split for parsing
		std::vector< std::string> split ( utility::string_split( j, '(' ) );
		//std::cout<< dof_line <<" "<<j<<" "<<split[0]<<std::endl;
		if ( split[0] == "x" ) dof_type = X_DOF;
		if ( split[0] == "y" ) dof_type = Y_DOF;
		if ( split[0] == "z" ) dof_type = Z_DOF;
		if ( split[0] == "angle_x" ) dof_type = X_ANGLE_DOF;
		if ( split[0] == "angle_y" ) dof_type = Y_ANGLE_DOF;
		if ( split[0] == "angle_z" ) dof_type = Z_ANGLE_DOF;
		if ( dof_type == 0 ) utility_exit_with_message("Dof type must be x,y,z,x_angle,y_angle,z_angle...");

		// Allow dof is true
		allowed_dof_jumps_[dof_type] = true;

		// Parse the rest
		if ( split.size() == 2 ) {
			std::vector< std::string> direction_split ( utility::string_split( split[1], ';' ) );
			// Parse the range1
			if ( direction_split.size() >= 1 ) {
				std::vector< std::string> range_split ( utility::string_split( direction_split[ 0 ], ':' ) );
				if ( range_split.size() >= 1 ) {
					lower_range_dof_jumps1_[dof_type] = std::atof( range_split[0].c_str() );
					has_range1_lower_[dof_type] = true;
					if ( range_split.size() == 2  && range_split[0] != range_split[1] ) {
						upper_range_dof_jumps1_[dof_type] = std::atof( range_split[1].c_str() );
						has_range1_upper_[dof_type] = true;
					} else {
						upper_range_dof_jumps1_[dof_type] = lower_range_dof_jumps1_[dof_type];
					}
				}
			}
			// Parse the range2
			if ( direction_split.size() >= 2 ) {
				std::vector< std::string> range_split ( utility::string_split( direction_split[ 1 ], ':' ) );
				if ( range_split.size() >= 1 ) {
					lower_range_dof_jumps2_[dof_type] = std::atof( range_split[0].c_str() );
					has_range2_lower_[dof_type] = true;
					if ( range_split.size() == 2 && range_split[0] != range_split[1] ) {
						upper_range_dof_jumps2_[dof_type] = std::atof( range_split[1].c_str() );
						has_range2_upper_[dof_type] = true;
					} else {
						upper_range_dof_jumps2_[dof_type] = lower_range_dof_jumps2_[dof_type];
					}
				}
			}
			// Parse the jump direction. Either dof_type(n2c) or dof_type(c2n)
			if ( direction_split.size() == 3 ) {
				if ( direction_split[2] == "n2c" ) {
					jump_dir_[dof_type] = 1;
				} else if ( direction_split[2] == "c2n" ) {
					jump_dir_[dof_type] = -1;
				} else {
					utility_exit_with_message("Unknown jump direction in Dof parsing...");
				}
			}
		}
	}
}

std::ostream& operator<< ( std::ostream & s, const SymDof & dof )
{

	 for ( Size i=1; i<=6; ++i ) {
		 if( dof.allow_dof(i) ) {
			 std::string dir ( "n2c" );
			 if ( dof.jump_direction(i) == -1 ) dir = "c2n";
			 if ( i == 1 ) s << "x";
			 if ( i == 2 ) s << "y";
			 if ( i == 3 ) s << "z";
			 if ( i == 4 ) s << "angle_x";
			 if ( i == 5 ) s << "angle_y";
			 if ( i == 6 ) s << "angle_z";
			 s << "(" << dof.range1_lower(i) << ":" << dof.range1_upper(i) << ";"
				 << dof.range2_lower(i) << ":" << dof.range2_upper(i) << ";" << dir << ")";
		 }
	 }

	 return s;
}

void
SymDof::add_dof_from_string(
										std::vector< std::string > dof_line )
{
	assert( dof_line.size() >= 3 );

	for ( Size i = 3; i <= dof_line.size(); ++i ) {
		read(dof_line[i-1]);
	}
}

} // symmetry
} // conformation
} // core
