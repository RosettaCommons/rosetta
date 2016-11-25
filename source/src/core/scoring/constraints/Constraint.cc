// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

// Unit Headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/pose/Pose.hh>

// Package Headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/XYZ_Func.hh>

// Project Headers

#include <core/id/AtomID.hh>

// Utility header
#include <utility/vector1.hh>

// C++ Headers

#include <algorithm>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

Constraint::Constraint( ScoreType const & t ): score_type_(t) {}

Constraint::~Constraint() {}

ConstraintOP Constraint::clone( core::scoring::func::FuncOP ) const {
	unimplemented_method_error( std::string("clone" ) );
	return NULL;
}

ConstraintOP Constraint::remapped_clone(
	pose::Pose const& /*src*/,
	pose::Pose const& /*dest*/,
	id::SequenceMappingCOP map=NULL ) const {
	unimplemented_method_error( std::string("remapped_clone" ) );
	if ( !map ) return NULL; // to make compile happy
	return NULL;
}

utility::vector1< core::Size >
Constraint::residues() const {
	utility::vector1< int > pos_list;
	for ( Size i=1; i<= natoms(); ++i ) {
		int const seqpos( atom(i).rsd() );
		// seqpos already in list?
		if ( std::find( pos_list.begin(), pos_list.end(), seqpos ) == pos_list.end() ) {
			pos_list.push_back( seqpos );
		}
	}
	return pos_list;
}

void Constraint::read_constraint( std::istream & /*in*/, core::pose::Pose const & /*pose*/) {
	unimplemented_method_error( std::string( "read_constraint" ) );
}

void Constraint::read_data( std::istream & ) {}

ConstraintOP
Constraint::remap_resid( core::id::SequenceMapping const &/*seqmap*/ ) const
{
	unimplemented_method_error( std::string( "remap_resid" ) );
	return NULL;
}

/// @brief Returns the unweighted score of this constraint computed over the given pose.
Real
Constraint::score( pose::Pose const& pose ) const {
	// Just set every weight to 1.0, to capture all results
	EnergyMap weights;
	for ( Size st(1); st <= Size(n_score_types); ++st ) {
		weights[ ScoreType(st) ] = 1.0;
	}

	return score( pose, weights );
}

/// @brief Returns the weighted score of this constraint computed over the given pose.
Real
Constraint::score( pose::Pose const& pose,  EnergyMap const & weights ) const {
	func::ConformationXYZ xyz_func( pose.conformation() );
	EnergyMap emap;
	score( xyz_func, weights, emap );
	return emap.dot( weights );
}

/// @brief return the raw "distance" before that distance is handed to the FUNC object
Real
Constraint::dist( core::pose::Pose const & pose ) const {
	func::ConformationXYZ xyz_func( pose.conformation() );
	return dist( xyz_func );
}

std::string Constraint::type() const {
	return "UNKNOWN_TYPE";
}

void Constraint::setup_for_scoring( core::scoring::func::XYZ_Func const &, ScoreFunction const & ) const {}

void Constraint::setup_for_derivatives( core::scoring::func::XYZ_Func const &, ScoreFunction const & ) const {}

void Constraint::show( std::ostream & /*out*/ ) const {
	unimplemented_method_error( std::string( "show" ) );
}

void Constraint::show_def( std::ostream & /*out*/, pose::Pose const & ) const {
	unimplemented_method_error( std::string( "show_def" ) );
}

void Constraint::read_def( std::istream &, pose::Pose const &, core::scoring::func::FuncFactory const & ) {
	unimplemented_method_error( std::string( "read_def" ) );
}

void Constraint::steal_def( pose::Pose const& ) {
	unimplemented_method_error( std::string( "steal_def" ) );
}

std::string Constraint::to_string() const {
	std::ostringstream out;
	show(out);
	return out.str();
}

Size
Constraint::show_violations(
	std::ostream & out,
	pose::Pose const &,
	Size,
	Real /*threshold*/
) const {
	out << "Constraint_show_violation stubbed out!\n" ;
	//threshold = 1; //to make compile happy // set but never used ~Labonte
	return 0;
}

core::scoring::func::Func const &
Constraint::get_func() const {
	unimplemented_method_error( std::string( "get_func" ) );
	static core::scoring::func::HarmonicFunc dummy_func( 0.0, 0.0);
	return dummy_func; // satisfy compiler
}

bool Constraint::operator != ( Constraint const & other ) const{
	return !(*this == other);
}

core::Size Constraint::choose_effective_sequence_separation(
	core::kinematics::ShortestPathInFoldTree const& sp,
	numeric::random::RandomGenerator&
) {
	return effective_sequence_separation( sp );
}

core::Size Constraint::effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& ) const {
	return 0;
}

void Constraint::unimplemented_method_error( std::string const & method_name ) const {
	utility_exit_with_message(
		"Called Constraint::" + method_name + " method from derived class " +
		type() + "," + "ended up in Constraint::" + method_name + "\n"
	);
}

std::ostream& operator<< ( std::ostream & out, Constraint const & cst ) {
	cst.show( out );
	return out;
}


}
}
}


#ifdef    SERIALIZATION
/// @details Default constructor that's needed by Cereal in order to deserialize a constraint.
/// This should ONLY be used by derived classes in their default constructors and then only
/// for the sake of deserializing a constraint.  The default constructor initializes the score_type_
/// to atom_pair_constraint -- a dummy value -- and then will const-cast the score_type_ to
/// deserialize the actual score_type_ value in its load method.
core::scoring::constraints::Constraint::Constraint() :
	score_type_( atom_pair_constraint )
{}

template< class Archive >
void
core::scoring::constraints::Constraint::save( Archive & arc ) const {
	arc( CEREAL_NVP( score_type_ ) );
}

template< class Archive >
void
core::scoring::constraints::Constraint::load( Archive & arc ) {
	arc( const_cast< ScoreType & > ( score_type_ ) ); // const enum core::scoring::ScoreType
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::Constraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::Constraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_Constraint )
#endif // SERIALIZATION
