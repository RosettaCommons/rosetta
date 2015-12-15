// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  Class to store ingformation about symmetrical dofs
/// @file   core/conformation/symmetry/SymDof.cc
/// @author Ingemar Andre

// Unit headers
#include <core/conformation/symmetry/SymDof.hh>

// Utility header
#include <utility/exit.hh>
#include <utility/string_util.hh>

#include <utility/vector1.hh>
#include <algorithm>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

#endif // SERIALIZATION


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
		range2_is_bound_.push_back(false);
		jump_dir_.push_back( 1 );
	}

}

SymDof::SymDof( SymDof const & src )
: range2_is_bound_( src.range2_is_bound_ ),
	allowed_dof_jumps_( src.allowed_dof_jumps_ ),
	lower_range_dof_jumps1_( src.lower_range_dof_jumps1_ ),
	upper_range_dof_jumps1_( src.upper_range_dof_jumps1_ ),
	lower_range_dof_jumps2_( src.lower_range_dof_jumps2_ ),
	upper_range_dof_jumps2_( src.upper_range_dof_jumps2_ ),
	has_range1_lower_( src.has_range1_lower_ ),
	has_range1_upper_( src.has_range1_upper_ ),
	has_range2_lower_( src.has_range2_lower_ ),
	has_range2_upper_( src.has_range2_upper_ ),
	jump_dir_( src.jump_dir_ )
{
}

SymDof &
SymDof::operator=( SymDof const & src ) {
	if ( this != &src ) {
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
		range2_is_bound_ = src.range2_is_bound_;
	}
	return *this;
}

SymDof::~SymDof() {}

// @details is df allowed to move?
bool
SymDof::allow_dof( int df ) const
{
	debug_assert( df >= X_DOF && df <= Z_ANGLE_DOF );
	debug_assert( allowed_dof_jumps_.size() == 6 );
	return allowed_dof_jumps_[df];

}

// @details is df allowed to move?
void
SymDof::set_allow_dof( int df, bool newval )
{
	debug_assert( df >= X_DOF && df <= Z_ANGLE_DOF );
	debug_assert( allowed_dof_jumps_.size() == 6 );
	allowed_dof_jumps_[df] = newval;
}

bool
SymDof::has_dof()
{
	for ( int i = X_DOF; i <= Z_ANGLE_DOF; ++i ) {
		if ( allow_dof(i) ) return true;
	}
	return false;
}

// @details the lower boundary of range1
core::Real
SymDof::range1_lower( int df ) const
{
	debug_assert( df >= X_DOF && df <= Z_ANGLE_DOF );
	debug_assert( lower_range_dof_jumps1_.size() == 6 );
	return lower_range_dof_jumps1_[df];
}

// @details the upper boundary of range1
core::Real
SymDof::range1_upper( int df ) const
{
	debug_assert( df >= X_DOF && df <= Z_ANGLE_DOF );
	debug_assert( upper_range_dof_jumps1_.size() == 6 );
	return upper_range_dof_jumps1_[df];
}

// @details the lower boundary of range2
core::Real
SymDof::range2_lower( int df ) const
{
	debug_assert( df >= X_DOF && df <= Z_ANGLE_DOF );
	debug_assert( lower_range_dof_jumps2_.size() == 6 );
	return lower_range_dof_jumps2_[df];
}

// @details the upper boundary of range1
core::Real
SymDof::range2_upper( int df ) const
{
	debug_assert( df >= X_DOF && df <= Z_ANGLE_DOF );
	debug_assert( upper_range_dof_jumps2_.size() == 6 );
	return upper_range_dof_jumps2_[df];
}

// details Have a range1 been specified?
bool
SymDof::has_range1( int df ) const
{
	if ( has_range1_lower_[df] && has_range1_upper_[df] ) return true;
	else return false;
}

// details Have a range2 been specified?
bool
SymDof::has_range2( int df ) const
{
	if ( has_range2_lower_[df] && has_range2_upper_[df] ) return true;
	else return false;
}

// @details has a lower boundary of range1 been specified?
bool
SymDof::has_range1_lower( int df ) const
{
	debug_assert( df >= X_DOF && df <= Z_ANGLE_DOF );
	return has_range1_lower_[df];
}

// @details has a upper boundary of range1 been specified?
bool
SymDof::has_range1_upper( int df ) const
{
	debug_assert( df >= X_DOF && df <= Z_ANGLE_DOF );
	return has_range1_upper_[df];
}

// @details has a lower boundary of range2 been specified?
bool
SymDof::has_range2_lower( int df ) const
{
	debug_assert( df >= X_DOF && df <= Z_ANGLE_DOF );
	return has_range2_lower_[df];
}

// @details has a upper boundary of range2 been specified?
bool
SymDof::has_range2_upper( int df ) const
{
	debug_assert( df >= X_DOF && df <= Z_ANGLE_DOF );
	return has_range2_upper_[df];
}

// @details has a upper boundary of range2 been specified?
bool
SymDof::range2_is_bound( int df ) const
{
	debug_assert( df >= X_DOF && df <= Z_ANGLE_DOF );
	return range2_is_bound_[df];
}

// @detail return the direction( upstream or downstream )
// of the jump for a dof
int
SymDof::jump_direction( int df ) const
{
	debug_assert( df >= X_DOF && df <= Z_ANGLE_DOF );
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
	// replace parenthesis/bracket with space for easier parsing
	std::replace( dof_line.begin(), dof_line.end(), ')', ' ' );
	std::replace( dof_line.begin(), dof_line.end(), ']', ' ' );

	std::istringstream l( dof_line );
	while ( true ) {
		std::string j;
		l >> j;
		if ( l.fail() ) break;
		// first read dof_type
		int dof_type(0);

		//  Split for parsing
		bool has_bracket = (j.find('[') != j.npos);
		utility::vector1< std::string> split ( utility::string_split( j, has_bracket?'[':'(' ) );

		if ( split[1] == "x" ) dof_type = X_DOF;
		if ( split[1] == "y" ) dof_type = Y_DOF;
		if ( split[1] == "z" ) dof_type = Z_DOF;
		if ( split[1] == "angle_x" ) dof_type = X_ANGLE_DOF;
		if ( split[1] == "angle_y" ) dof_type = Y_ANGLE_DOF;
		if ( split[1] == "angle_z" ) dof_type = Z_ANGLE_DOF;
		if ( dof_type == 0 ) utility_exit_with_message("Dof type must be x,y,z,x_angle,y_angle,z_angle...");

		// bracket implies that range2 is an absolute bound
		if ( has_bracket ) range2_is_bound_[dof_type] = true;

		// range2_is_bound_ unsupported for rotations
		runtime_assert ( !has_bracket || dof_type==X_DOF || dof_type==Y_DOF || dof_type==Z_DOF );

		// Allow dof is true
		allowed_dof_jumps_[dof_type] = true;

		// Parse the rest
		if ( split.size() == 2 ) {
			utility::vector1< std::string> direction_split ( utility::string_split( split[2], ';' ) );
			// Parse the range1
			if ( direction_split.size() >= 1 ) {
				utility::vector1< std::string> range_split ( utility::string_split( direction_split[ 1 ], ':' ) );
				if ( range_split.size() >= 1 ) {
					lower_range_dof_jumps1_[dof_type] = std::atof( range_split[1].c_str() );
					has_range1_lower_[dof_type] = true;
					if ( range_split.size() == 2  && range_split[1] != range_split[2] ) {
						upper_range_dof_jumps1_[dof_type] = std::atof( range_split[2].c_str() );
						has_range1_upper_[dof_type] = true;
					} else {
						upper_range_dof_jumps1_[dof_type] = lower_range_dof_jumps1_[dof_type];
					}
				}
			}
			// Parse the range2
			if ( direction_split.size() >= 2 ) {
				utility::vector1< std::string> range_split ( utility::string_split( direction_split[ 2 ], ':' ) );
				if ( range_split.size() >= 1 ) {
					lower_range_dof_jumps2_[dof_type] = std::atof( range_split[1].c_str() );
					has_range2_lower_[dof_type] = true;
					if ( range_split.size() == 2 && range_split[1] != range_split[2] ) {
						upper_range_dof_jumps2_[dof_type] = std::atof( range_split[2].c_str() );
						has_range2_upper_[dof_type] = true;
					} else {
						upper_range_dof_jumps2_[dof_type] = lower_range_dof_jumps2_[dof_type];
					}
				}

				// quick sanity check
				if ( range2_is_bound_[dof_type] ) {
					runtime_assert ( lower_range_dof_jumps2_[dof_type] <= lower_range_dof_jumps1_[dof_type] &&
						upper_range_dof_jumps2_[dof_type] >= upper_range_dof_jumps1_[dof_type] );
				}
			}
			// Parse the jump direction. Either dof_type(n2c) or dof_type(c2n)
			if ( direction_split.size() == 3 ) {
				if ( direction_split[3] == "n2c" ) {
					jump_dir_[dof_type] = 1;
				} else if ( direction_split[3] == "c2n" ) {
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
		if ( dof.allow_dof(i) ) {
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
SymDof::add_dof_from_string( utility::vector1< std::string > dof_line )
{
	debug_assert( dof_line.size() >= 3 );

	for ( Size i = 3; i <= dof_line.size(); ++i ) {
		read(dof_line[i]);
	}
}

bool
operator==(
	SymDof const & a,
	SymDof const & b
) {
	return
		std::equal(a.allowed_dof_jumps_.begin(), a.allowed_dof_jumps_.end(), b.allowed_dof_jumps_.begin()) &&
		std::equal(a.lower_range_dof_jumps1_.begin(), a.lower_range_dof_jumps1_.end(), b.lower_range_dof_jumps1_.begin()) &&
		std::equal(a.upper_range_dof_jumps1_.begin(), a.upper_range_dof_jumps1_.end(), b.upper_range_dof_jumps1_.begin()) &&
		std::equal(a.lower_range_dof_jumps2_.begin(), a.lower_range_dof_jumps2_.end(), b.lower_range_dof_jumps2_.begin()) &&
		std::equal(a.upper_range_dof_jumps2_.begin(), a.upper_range_dof_jumps2_.end(), b.upper_range_dof_jumps2_.begin()) &&
		std::equal(a.has_range1_lower_.begin(), a.has_range1_lower_.end(), b.has_range1_lower_.begin()) &&
		std::equal(a.has_range1_upper_.begin(), a.has_range1_upper_.end(), b.has_range1_upper_.begin()) &&
		std::equal(a.has_range2_lower_.begin(), a.has_range2_lower_.end(), b.has_range2_lower_.begin()) &&
		std::equal(a.has_range2_upper_.begin(), a.has_range2_upper_.end(), b.has_range2_upper_.begin()) &&
		std::equal(a.jump_dir_.begin(), a.jump_dir_.end(), b.jump_dir_.begin());
}

bool
operator!=(
	SymDof const & a,
	SymDof const & b
) {
	return !(a == b);
}

} // symmetry
} // conformation
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::symmetry::SymDof::save( Archive & arc ) const {
	arc( CEREAL_NVP( range2_is_bound_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( allowed_dof_jumps_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( lower_range_dof_jumps1_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( upper_range_dof_jumps1_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( lower_range_dof_jumps2_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( upper_range_dof_jumps2_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( has_range1_lower_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( has_range1_upper_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( has_range2_lower_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( has_range2_upper_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( jump_dir_ ) ); // utility::vector1<int>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::symmetry::SymDof::load( Archive & arc ) {
	arc( range2_is_bound_ ); // utility::vector1<_Bool>
	arc( allowed_dof_jumps_ ); // utility::vector1<_Bool>
	arc( lower_range_dof_jumps1_ ); // utility::vector1<Real>
	arc( upper_range_dof_jumps1_ ); // utility::vector1<Real>
	arc( lower_range_dof_jumps2_ ); // utility::vector1<Real>
	arc( upper_range_dof_jumps2_ ); // utility::vector1<Real>
	arc( has_range1_lower_ ); // utility::vector1<_Bool>
	arc( has_range1_upper_ ); // utility::vector1<_Bool>
	arc( has_range2_lower_ ); // utility::vector1<_Bool>
	arc( has_range2_upper_ ); // utility::vector1<_Bool>
	arc( jump_dir_ ); // utility::vector1<int>
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::symmetry::SymDof );
#endif // SERIALIZATION
