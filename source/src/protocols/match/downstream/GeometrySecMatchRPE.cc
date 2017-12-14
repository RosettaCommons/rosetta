// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/GeometrySecMatchRPE.cc
/// @brief
/// @author Florian Richter, floric@u.washington.edu, june 09

// Unit headers
#include <protocols/match/downstream/GeometrySecMatchRPE.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>

// Package headers

// Project headers
#include <core/conformation/Residue.hh>
#include <basic/basic.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/string_util.hh>

// C++ headers
#include <list>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

AtomGeometrySecMatchRPE::AtomGeometrySecMatchRPE(
	protocols::toolbox::match_enzdes_util::GeomSampleInfo const & gsi
) :
	SecMatchResiduePairEvaluator(),
	lowval_(gsi.ideal_val() - gsi.tolerance() ),
	highval_( gsi.ideal_val() + gsi.tolerance() )
{
	clear_at_inds();
}

AtomGeometrySecMatchRPE::~AtomGeometrySecMatchRPE() = default;


bool
AtomGeometrySecMatchRPE::require_all_target_residue_atom_coordinates() const
{
	return false;
}

bool
AtomGeometrySecMatchRPE::require_target_atom_coordinate( Size target_atom_id ) const
{
	for ( Size ii = 1; ii <= at_inds_.size(); ++ii ) {
		if ( at_inds_[ ii ].first == 2 && at_inds_[ ii ].second == target_atom_id ) {
			return true;
		}
	}
	return false;
}


bool
AtomGeometrySecMatchRPE::check_value(
	core::Real value ) const
{
	if ( value > highval_ ) return false;
	else if ( value < lowval_ ) return false;
	else return true;
}

void
AtomGeometrySecMatchRPE::add_at_ind(
	core::Size which_cst_res,
	core::Size atom_ind_in_res
){
	at_inds_.push_back( std::make_pair( which_cst_res, atom_ind_in_res ) );
}


void
AtomGeometrySecMatchRPE::clear_at_inds(){
	at_inds_.clear();
}

void
AtomGeometrySecMatchRPE::set_lowval(
	core::Real lowval
)
{
	lowval_ = lowval;
}

void
AtomGeometrySecMatchRPE::set_highval(
	core::Real highval
)
{
	highval_ = highval;
}


AtomDistanceSecMatchRPE::AtomDistanceSecMatchRPE(
	protocols::toolbox::match_enzdes_util::GeomSampleInfo const & gsi
) :
	AtomGeometrySecMatchRPE( gsi )
{
	set_lowval( lowval() * lowval() );
	set_highval( highval() * highval() );
}


bool
AtomDistanceSecMatchRPE::evaluate_residues(
	core::conformation::Residue const & candidate_res,
	core::conformation::Residue const & target_res
) const
{

	debug_assert( at_inds().size() == 2 );

	core::Real distance_squared( candidate_res.xyz( at_inds()[1].second ).distance_squared( target_res.xyz( at_inds()[2].second ) ) );

	/*if ( ! check_value( distance_squared ) ) {
	std::cout << "AtomDistanceSecMatchRPE::evaluate_residues fail: " << distance_squared << std::endl;
	}*/
	return check_value( distance_squared );

}

std::string
AtomDistanceSecMatchRPE::print(
	core::chemical::ResidueTypeCOP candidate_restype,
	core::chemical::ResidueTypeCOP target_restype
) const
{
	return "AtomDistance range " + utility::to_string( std::sqrt(lowval()) ) +
		" A to " +  utility::to_string( std::sqrt(highval()) ) + " A between " +
		utility::trim( at_inds()[ 1 ].first == 1 ? candidate_restype->name() : target_restype->name() ) +
		" atom " +
		utility::trim( at_inds()[ 1 ].first == 1 ? candidate_restype->atom_name( at_inds()[ 1 ].second ) : target_restype->atom_name( at_inds()[ 1 ].second )) +
		" and " +
		utility::trim( at_inds()[ 2 ].first == 1 ? candidate_restype->name() : target_restype->name() ) +
		" atom " +
		utility::trim( at_inds()[ 2 ].first == 1 ? candidate_restype->atom_name( at_inds()[ 2 ].second ) : target_restype->atom_name( at_inds()[ 2 ].second ));
}


bool
AtomDistanceSecMatchRPE::require_candidate_residue_atoms_to_lie_near_target_atom(
	Size target_atom_id
) const
{
	debug_assert( at_inds().size() == 2 );
	debug_assert( at_inds()[ 2 ].first == 2 ); // second atom is to the target residue

	return at_inds()[ 2 ].second == target_atom_id;
}

utility::vector1< AtomDistanceSecMatchRPE::Size >
AtomDistanceSecMatchRPE::candidate_res_atoms_reqd_near_target_atom(
	Size target_atom_id
) const
{
	if ( at_inds()[ 2 ].second == target_atom_id ) {
		utility::vector1< Size > other_atom( 1, at_inds()[ 1 ].second );
		return other_atom;
	} else {
		utility::vector1< Size > empty;
		return empty;
	}
}

AtomDistanceSecMatchRPE::Real
AtomDistanceSecMatchRPE::max_separation_dist_to_target_atom( Size target_atom_id ) const
{
	if ( at_inds()[ 2 ].second == target_atom_id ) {
		return std::sqrt( highval() );
	} else {
		return -1.0;
	}
}


AtomAngleSecMatchRPE::AtomAngleSecMatchRPE(
	protocols::toolbox::match_enzdes_util::GeomSampleInfo const & gsi
) :
	AtomGeometrySecMatchRPE( gsi )
{
	set_lowval( lowval() * numeric::constants::d::degrees_to_radians  );
	set_highval(  highval() * numeric::constants::d::degrees_to_radians  );

}

bool
AtomAngleSecMatchRPE::evaluate_residues(
	core::conformation::Residue const & candidate_res,
	core::conformation::Residue const & target_res
) const
{
	debug_assert( at_inds().size() == 3 );

	//first figure out which of the residues the atoms belong to
	core::PointPosition p1( at_inds()[1].first == 1 ? candidate_res.xyz( at_inds()[1].second ) : target_res.xyz( at_inds()[1].second ) );
	core::PointPosition p2( at_inds()[2].first == 1 ? candidate_res.xyz( at_inds()[2].second ) : target_res.xyz( at_inds()[2].second ) );
	core::PointPosition p3( at_inds()[3].first == 1 ? candidate_res.xyz( at_inds()[3].second ) : target_res.xyz( at_inds()[3].second ) );

	//same approach as angle constraint
	//core/scoring/constraints/AngleConstraint.cc
	numeric::xyzVector< core::Real > u1( p1 - p2 );
	numeric::xyzVector< core::Real > u2( p3 - p2 );
	core::Real const n1( u1.length() );
	core::Real const n2( u2.length() );

	core::Real angle( numeric::arccos( dot( u1,u2 ) / ( n1 * n2 ) ) );


	/*if ( ! check_value( angle ) ) {
	std::cout << "AtomAngleSecMatchRPE::evaluate_residues fail: " <<
	numeric::constants::d::radians_to_degrees * angle << " low " <<
	numeric::constants::d::radians_to_degrees * lowval() << " high " <<
	numeric::constants::d::radians_to_degrees * highval() << std::endl;
	}*/


	return check_value( angle );

}

std::string
AtomAngleSecMatchRPE::print(
	core::chemical::ResidueTypeCOP candidate_restype,
	core::chemical::ResidueTypeCOP target_restype
) const
{
	return "AtomAngle range " + utility::to_string( numeric::constants::d::radians_to_degrees * lowval() ) +
		" degrees to " +  utility::to_string( numeric::constants::d::radians_to_degrees * highval() ) + " degrees between " +
		utility::trim( at_inds()[ 1 ].first == 1 ? candidate_restype->name() : target_restype->name() ) +
		" atom " +
		utility::trim( at_inds()[ 1 ].first == 1 ? candidate_restype->atom_name( at_inds()[ 1 ].second ) : target_restype->atom_name( at_inds()[ 1 ].second )) +
		", "  +
		utility::trim( at_inds()[ 2 ].first == 1 ? candidate_restype->name() : target_restype->name() ) +
		" atom " +
		utility::trim( at_inds()[ 2 ].first == 1 ? candidate_restype->atom_name( at_inds()[ 2 ].second ) : target_restype->atom_name( at_inds()[ 2 ].second )) +
		", and " +
		utility::trim( at_inds()[ 3 ].first == 1 ? candidate_restype->name() : target_restype->name() ) +
		" atom " +
		utility::trim( at_inds()[ 3 ].first == 1 ? candidate_restype->atom_name( at_inds()[ 3 ].second ) : target_restype->atom_name( at_inds()[ 3 ].second ));
}

AtomDihedralSecMatchRPE::AtomDihedralSecMatchRPE(
	protocols::toolbox::match_enzdes_util::GeomSampleInfo const & gsi
) : AtomGeometrySecMatchRPE( gsi ),
	check_periodicity_(gsi.periodicity() != 360.0 ),
	periodicity_(gsi.periodicity() * numeric::constants::d::degrees_to_radians ),
	offset_( basic::periodic_range( gsi.ideal_val() * numeric::constants::d::degrees_to_radians,  periodicity_ ) )
{
	set_lowval( -( gsi.tolerance() * numeric::constants::d::degrees_to_radians ) );
	set_highval(  gsi.tolerance() * numeric::constants::d::degrees_to_radians  );
}

bool
AtomDihedralSecMatchRPE::evaluate_residues(
	core::conformation::Residue const & candidate_res,
	core::conformation::Residue const & target_res
) const
{
	debug_assert( at_inds().size() == 4 );

	//first figure out which of the residues the atoms belong to
	core::PointPosition p1( at_inds()[1].first == 1 ? candidate_res.xyz( at_inds()[1].second ) : target_res.xyz( at_inds()[1].second ) );
	core::PointPosition p2( at_inds()[2].first == 1 ? candidate_res.xyz( at_inds()[2].second ) : target_res.xyz( at_inds()[2].second ) );
	core::PointPosition p3( at_inds()[3].first == 1 ? candidate_res.xyz( at_inds()[3].second ) : target_res.xyz( at_inds()[3].second ) );
	core::PointPosition p4( at_inds()[4].first == 1 ? candidate_res.xyz( at_inds()[4].second ) : target_res.xyz( at_inds()[4].second ) );

	Real value( basic::periodic_range( numeric::dihedral_radians( p1, p2, p3, p4 ) - offset_, periodicity_ ) );

	/*if ( ! check_value( value )) {
	std::cout << "AtomDihedralSecMatchRPE::evaluate_residues fail: " <<
	numeric::constants::d::radians_to_degrees * numeric::dihedral_radians( p1, p2, p3, p4 ) << " offsetper corrected " <<
	numeric::constants::d::radians_to_degrees * value << " tol: " << numeric::constants::d::radians_to_degrees * highval() << std::endl;
	}*/

	return check_value( value );

}

std::string
AtomDihedralSecMatchRPE::print(
	core::chemical::ResidueTypeCOP candidate_restype,
	core::chemical::ResidueTypeCOP target_restype
) const
{
	std::string prefix = "AtomDihedral range" + std::string( std::abs( periodicity_ - numeric::constants::d::pi_2 ) > 1e-6 ? "s:" : ":" );
	auto const n_periods = static_cast< Size > ( numeric::constants::d::pi_2 / periodicity_ );
	for ( Size ii = 0; ii < n_periods; ++ii ) {
		Real lo_deg = numeric::constants::d::radians_to_degrees * ( ii * periodicity_ + offset_ + lowval() );
		Real hi_deg = numeric::constants::d::radians_to_degrees * ( ii * periodicity_ + offset_ + highval() );
		prefix += " [ " + utility::to_string( lo_deg ) + ", " + utility::to_string( hi_deg ) + " ]";
	}
	prefix += " degrees between " +
		utility::trim( at_inds()[ 1 ].first == 1 ? candidate_restype->name() : target_restype->name() ) +
		" atom " +
		utility::trim( at_inds()[ 1 ].first == 1 ? candidate_restype->atom_name( at_inds()[ 1 ].second ) : target_restype->atom_name( at_inds()[ 1 ].second )) +
		", " +
		utility::trim( at_inds()[ 2 ].first == 1 ? candidate_restype->name() : target_restype->name() ) +
		" atom " +
		utility::trim( at_inds()[ 2 ].first == 1 ? candidate_restype->atom_name( at_inds()[ 2 ].second ) : target_restype->atom_name( at_inds()[ 2 ].second )) +
		", " +
		utility::trim( at_inds()[ 3 ].first == 1 ? candidate_restype->name() : target_restype->name() ) +
		" atom " +
		utility::trim( at_inds()[ 3 ].first == 1 ? candidate_restype->atom_name( at_inds()[ 3 ].second ) : target_restype->atom_name( at_inds()[ 3 ].second )) +
		", and " +
		utility::trim( at_inds()[ 4 ].first == 1 ? candidate_restype->name() : target_restype->name() ) +
		" atom " +
		utility::trim( at_inds()[ 4 ].first == 1 ? candidate_restype->atom_name( at_inds()[ 4 ].second ) : target_restype->atom_name( at_inds()[ 4 ].second ));
	return prefix;
}

void
GeometrySecMatchRPE::add_atomgeom_evaluator(
	AtomGeometrySecMatchRPECOP evaluator
)
{
	atom_geom_rpes_.push_back( evaluator );
}


GeometrySecMatchRPE::GeometrySecMatchRPE(
	protocols::toolbox::match_enzdes_util::MatchConstraintFileInfo const & mcfi,
	utility::vector1< core::Size > const & downstream_inds,
	utility::vector1< core::Size > const & upstream_inds
){


	// TODO: Figure out who uses this to determine the best way to handle an mcfi with multiple constraints.
	// For now, let's just grab the first one:
	using namespace protocols::toolbox::match_enzdes_util;
	utility::vector1< SingleConstraint > constraints( mcfi.constraints() );
	assert( constraints.size() );
	SingleConstraint const & constraint = constraints[ 1 ];

	if ( constraint.dis_U1D1 ) {
		AtomDistanceSecMatchRPEOP adist( new AtomDistanceSecMatchRPE( *(constraint.dis_U1D1 ) ) );
		adist->add_at_ind( 1, upstream_inds[1] );
		adist->add_at_ind( 2, downstream_inds[1] );
		atom_geom_rpes_.push_back( adist );
	}

	if ( constraint.ang_U2D1 ) {
		AtomAngleSecMatchRPEOP aang1( new AtomAngleSecMatchRPE( *(constraint.ang_U2D1 ) ) );
		aang1->add_at_ind( 1, upstream_inds[2] );
		aang1->add_at_ind( 1, upstream_inds[1] );
		aang1->add_at_ind( 2, downstream_inds[1] );
		atom_geom_rpes_.push_back( aang1 );
	}

	if ( constraint.ang_U1D2 ) {
		AtomAngleSecMatchRPEOP aang2( new AtomAngleSecMatchRPE( *(constraint.ang_U1D2 ) ) );
		aang2->add_at_ind( 1, upstream_inds[1] );
		aang2->add_at_ind( 2, downstream_inds[1] );
		aang2->add_at_ind( 2, downstream_inds[2] );
		atom_geom_rpes_.push_back( aang2 );
	}

	if ( constraint.tor_U3D1 ) {
		AtomDihedralSecMatchRPEOP adih1( new AtomDihedralSecMatchRPE( *(constraint.tor_U3D1 ) ) );
		adih1->add_at_ind( 1, upstream_inds[3] );
		adih1->add_at_ind( 1, upstream_inds[2] );
		adih1->add_at_ind( 1, upstream_inds[1] );
		adih1->add_at_ind( 2, downstream_inds[1] );
		atom_geom_rpes_.push_back( adih1 );
	}

	if ( constraint.tor_U2D2 ) {
		AtomDihedralSecMatchRPEOP adih2( new AtomDihedralSecMatchRPE( *(constraint.tor_U2D2 ) ) );
		adih2->add_at_ind( 1, upstream_inds[2] );
		adih2->add_at_ind( 1, upstream_inds[1] );
		adih2->add_at_ind( 2, downstream_inds[1] );
		adih2->add_at_ind( 2, downstream_inds[2] );
		atom_geom_rpes_.push_back( adih2 );
	}

	if ( constraint.tor_U1D3 ) {
		AtomDihedralSecMatchRPEOP adih3( new AtomDihedralSecMatchRPE( *(constraint.tor_U1D3 ) ) );
		adih3->add_at_ind( 1, upstream_inds[1] );
		adih3->add_at_ind( 2, downstream_inds[1] );
		adih3->add_at_ind( 2, downstream_inds[2] );
		adih3->add_at_ind( 2, downstream_inds[3] );
		atom_geom_rpes_.push_back( adih3 );
	}
}

bool
GeometrySecMatchRPE::evaluate_residues(
	core::conformation::Residue const & candidate_res,
	core::conformation::Residue const & target_res
) const
{
	for ( auto const & atom_geom_rpe : atom_geom_rpes_ ) {
		if ( ! atom_geom_rpe->evaluate_residues( candidate_res, target_res ) ) return false;
	}
	//if we've made it to here, that means all AtomGeomRPEs returned true
	return true;
}


bool
GeometrySecMatchRPE::require_all_target_residue_atom_coordinates() const
{
	return false;
}

bool
GeometrySecMatchRPE::require_target_atom_coordinate( Size target_atom_id ) const
{
	for ( auto const & atom_geom_rpe : atom_geom_rpes_ ) {
		if ( atom_geom_rpe->require_target_atom_coordinate( target_atom_id ) ) return true;
	}
	return false;
}

bool
GeometrySecMatchRPE::require_candidate_residue_atoms_to_lie_near_target_atom( Size target_atom_id ) const
{
	for ( auto const & atom_geom_rpe : atom_geom_rpes_ ) {
		if ( atom_geom_rpe->require_candidate_residue_atoms_to_lie_near_target_atom( target_atom_id ) ) return true;
	}
	return false;
}

/// @details Aggregate the sets of atoms that are required to be near a given
/// target atom from the various AtomGeometry evaluators
utility::vector1< GeometrySecMatchRPE::Size >
GeometrySecMatchRPE::candidate_res_atoms_reqd_near_target_atom(
	Size target_atom_id
) const
{
	std::list< Size > atoms;
	for ( auto const & atom_geom_rpe : atom_geom_rpes_ ) {
		utility::vector1< Size > reqd_atoms = atom_geom_rpe->candidate_res_atoms_reqd_near_target_atom( target_atom_id );
		if ( reqd_atoms.size() != 0 ) {
			for ( Size ii = 1; ii <= reqd_atoms.size(); ++ii ) {
				atoms.push_back( reqd_atoms[ ii ] );
				atoms.sort();
				atoms.unique();
			}
		}
	}

	utility::vector1< Size > atom_vect( atoms.size() );
	Size count( 1 );
	for ( std::list< Size >::const_iterator iter = atoms.begin(),
			iter_end = atoms.end(); iter != iter_end; ++iter ) {
		atom_vect[ count ] = *iter;
		++count;
	}
	return atom_vect;
}

/// @details Return the shortest of the distance cutoffs from the AtomGeometry that
/// do describe a distance cutoff to a particular target atom.
GeometrySecMatchRPE::Real
GeometrySecMatchRPE::max_separation_dist_to_target_atom( Size target_atom_id ) const
{
	Real min_max_dis( -1.0 );
	for ( auto const & atom_geom_rpe : atom_geom_rpes_ ) {
		Real max_dis = atom_geom_rpe->max_separation_dist_to_target_atom( target_atom_id );
		if ( max_dis > 0.0 && min_max_dis > max_dis ) {
			min_max_dis = max_dis;
		}
	}

	return min_max_dis;

}


}
}
}
