// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
// This file is based off of code from the Biochemistry Library (BCL).
// The BCL is copyright Vanderbilt University (Meiler Lab), a RosettaCommons member

/// @file   core/chemical/modifications/ValenceHandler.cc
/// @brief  A bunch of functions to determine where to place atoms based on hybridization of atoms and number of bonds
/// It should be noted that a vector of coordinates are returned. This is because it creates the angles for all
/// hydrogens that can be placed. If you are combining a fragment, you really only need the first index of coords
/// @author Rosetta conversion: Steven Combs / Sandeep Kothiwale

#include <core/chemical/modifications/ValenceHandler.hh>
#include <numeric/constants.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <numeric/xyzVector.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>

#include <basic/Tracer.hh>

namespace core {
namespace chemical {
namespace modifications {

static basic::Tracer TR("core.chemical.modifications.ValenceHandler");

utility::vector1<numeric::xyzVector<core::Real> > determine_coordinates(core::chemical::MutableResidueType const & res, core::chemical::VD const & atom) {
	utility::vector1<numeric::xyzVector<core::Real> > positions;

	core::chemical::gasteiger::GasteigerAtomTypeDataCOP gasteiger_type( res.atom(atom).gasteiger_atom_type() );
	if ( ! gasteiger_type ) {
		// Skip "unknown" atoms
		TR << "Skipping atom " << res.atom_name( atom ) << " when determining missing hydrogen coordinates because it does not have an assigned Gasteiger type." << std::endl;
		return positions;
	}

	// get the hybrid orbital for the atom type
	core::chemical::gasteiger::GasteigerAtomTypeData::HybridOrbitalType orbital_type( gasteiger_type->get_hybrid_orbital_type() );

	if ( orbital_type == core::chemical::gasteiger::GasteigerAtomTypeData::Unhybridized ) {
		TR << "Skipping atom " << res.atom_name( atom ) << " when determining missing hydrogen coordinates because its Gasteiger type (" << gasteiger_type->get_name() <<") is unhybridized." << std::endl;
		return positions;
	}

	// Can't use index-based bonded_neighbor(), as this ResidueType might not be finalized yet.
	utility::vector1< core::chemical::VD > bonded_atoms;
	for ( core::chemical::AdjacentIterPair bonded_iter( res.bonded_neighbor_iterators(atom) );
			bonded_iter.first != bonded_iter.second;
			++bonded_iter.first ) {
		bonded_atoms.push_back( *bonded_iter.first );
	}

	core::Size const nr_known_bonds( bonded_atoms.size() );

	//We assume here that each bonding partner uses one hybrid bond.
	core::Size const nr_possible_bonds( gasteiger_type->get_number_hybrid_bonds() );

	if ( nr_known_bonds > nr_possible_bonds ) {
		TR.Warning << "For atom " << res.atom_name(atom) << " with type " << gasteiger_type->get_name()
			<< " number of known bonds (" << nr_known_bonds << ") exceeds number of possible bonds (" << nr_possible_bonds
			<< ") - refusing to add valences." << std::endl;
		return positions;
	}

	core::Size const nr_missing_bonds( nr_possible_bonds - nr_known_bonds );
	TR << "Missing " << nr_missing_bonds << " bonding partners on " << res.atom_name(atom) << " (" << nr_known_bonds << " known bonds.)" << std::endl;
	if ( nr_missing_bonds == 0 ) {
		//Exit early - nothing to do
		return positions;
	}

	// store the expected bond length, which is the covalent radius of H plus the covalent radius of ATOM
	// TODO: Does this need to be calculated, or is a constant "good enough"?
	const double bond_length( 1.0);

	// get the position of this atom
	const numeric::xyzVector<core::Real> position( res.atom(atom).ideal_xyz());

	//maybe we will need this functionality in the future, but for now, I can not imagine a case where the number of knwon bonds would be 0.

	/* if( nr_known_bonds == 0) {
	positions = GetIdealizedGeometry( orbital_type, bond_length);

	// translate the points such that they are relative to position
	for
	(
	storage::Vector< linal::Vector3D>::iterator itr( positions.Begin()), itr_end( positions.End());
	itr != itr_end;
	++itr
	)
	{
	*itr += position;
	}
	}*/
	if ( nr_known_bonds == 1 ) {
		if ( orbital_type == core::chemical::gasteiger::GasteigerAtomTypeData::SP ) {
			TR << res.atom_name(atom) << ": sp1 atom with 1 known bond to "<< res.atom(bonded_atoms[1]).name() << " and " << nr_missing_bonds << " missing. " << std::endl;
			assert( nr_missing_bonds == 1 );

			positions.push_back( linear_coordinates( position, res.atom(bonded_atoms[1]).ideal_xyz(), bond_length));
		} else if ( orbital_type == core::chemical::gasteiger::GasteigerAtomTypeData::SP2 ) {
			TR << res.atom_name(atom) << ": sp2 atom with 1 known bond to "<< res.atom(bonded_atoms[1]).name() << " and " << nr_missing_bonds << " missing. " << std::endl;
			assert( nr_missing_bonds <= 2);

			utility::vector1< core::chemical::VD > bonds_neighbor;
			for ( core::chemical::AdjacentIterPair bonded_iter( res.bonded_neighbor_iterators( bonded_atoms[1] ) );
					bonded_iter.first != bonded_iter.second;
					++bonded_iter.first ) {
				bonds_neighbor.push_back( *bonded_iter.first );
			}

			utility::vector1<numeric::xyzVector<core::Real> > neighboring_positions;
			for ( core::Size i=1; i <= bonds_neighbor.size(); ++i ) {
				if ( ( bonds_neighbor[i] != bonded_atoms[1] ) || ( bonds_neighbor[i] != atom ) ) {
					neighboring_positions.push_back( res.atom(bonds_neighbor[i]).ideal_xyz() );
					break;
				}
			}

			if ( !neighboring_positions.empty() ) {
				//Calculate the trans position from the first neighbor
				TR << "Calculating sp2 one-bonded from first neighbor." << std::endl;
				positions.push_back(position - (neighboring_positions[1] - res.atom(bonded_atoms[1]).ideal_xyz() ).normalize() * bond_length );
			} else {
				TR << "Calculating sp2 one-bonded from no-neighbors." << std::endl;
				// coordinate 1
				positions.push_back(
					angle_coordinates(
					position,
					res.atom(bonded_atoms[1]).ideal_xyz(),
					numeric::xyzVector<core::Real>( bond_length),
					bond_length,
					120.0 / 180.0 * numeric::constants::d::pi,
					120.0 / 180.0 * numeric::constants::d::pi,
					false,
					numeric::xyzVector<core::Real>( 0.0)
					)
				);
			}

			if ( nr_missing_bonds == 2 ) {
				TR << "Adding second position to sp2 one-bonded." << std::endl;
				// coordinate 2
				positions.push_back(
					triganol_coordinates(
					position,
					res.atom(bonded_atoms[1]).ideal_xyz(),
					positions[1],
					bond_length
					)
				);
			}
		} else { // SP3
			TR << res.atom_name(atom) << ": sp3 atom with 1 known bond to "<< res.atom(bonded_atoms[1]).name() << " and " << nr_missing_bonds << " missing. " << std::endl;
			assert( nr_missing_bonds <= 3 );

			// coordinate 1
			positions.push_back(
				angle_coordinates(
				position,
				res.atom(bonded_atoms[1]).ideal_xyz(),
				numeric::xyzVector<core::Real>( bond_length),
				bond_length,
				109.5 / 180.0 * numeric::constants::d::pi,
				109.5 / 180.0 * numeric::constants::d::pi,
				false,
				numeric::xyzVector<core::Real>( 0.0)
				)
			);

			if ( nr_missing_bonds >= 2 ) {
				// helper coordinates
				numeric::xyzVector<core::Real> foot_point(
					triganol_coordinates(
					position,
					res.atom(bonded_atoms[1]).ideal_xyz(),
					positions[1],
					bond_length * std::cos( 54.75 / 180 * numeric::constants::d::pi)
					)
				);
				numeric::xyzVector<core::Real> offset(
					bond_length * std::sin( 54.75 / 180 * numeric::constants::d::pi) *
					cross_product(
					res.atom(bonded_atoms[1]).ideal_xyz() - position,
					positions[1] - position
					).normalize()
				);
				// coordinate 2
				positions.push_back( foot_point + offset);

				if ( nr_missing_bonds >= 3 ) {
					// coordinate 3
					positions.push_back( foot_point - offset);
				}
			}
		}
	} else if ( nr_known_bonds == 2 ) { // 2 known bonds
		if ( orbital_type == core::chemical::gasteiger::GasteigerAtomTypeData::SP2 ) { // trigonal
			TR << res.atom_name(atom) << ": sp2 atom with 2 known bonds to "<< res.atom(bonded_atoms[1]).name() << " and " << res.atom(bonded_atoms[2]).name()
				<< " and " << nr_missing_bonds << " missing. " << std::endl;
			assert( nr_missing_bonds == 1 );
			positions.push_back(
				triganol_coordinates(
				position,
				res.atom(bonded_atoms[1]).ideal_xyz(),
				res.atom(bonded_atoms[2]).ideal_xyz(),
				bond_length
				)
			);
		} else { // tetrahedral/SP3 (Shouldn't get here with SP1)
			TR << res.atom_name(atom) << ": sp3 atom with 2 known bonds to "<< res.atom(bonded_atoms[1]).name() << " and " << res.atom(bonded_atoms[2]).name()
				<< " and " << nr_missing_bonds << " missing. " << std::endl;
			assert( nr_missing_bonds <= 2 );
			// helper coordinates
			numeric::xyzVector<core::Real> foot_point(
				triganol_coordinates(
				position,
				res.atom(bonded_atoms[1]).ideal_xyz(),
				res.atom(bonded_atoms[2]).ideal_xyz(),
				bond_length * std::cos( 54.75 / 180 * numeric::constants::d::pi)
				)
			);

			numeric::xyzVector<core::Real> offset
				(
				bond_length * std::sin( 54.75 / 180 * numeric::constants::d::pi) *
				cross_product(
				res.atom(bonded_atoms[1]).ideal_xyz() - position,
				res.atom(bonded_atoms[2]).ideal_xyz() - position
				).normalize()
			);

			// coordinate 1
			positions.push_back( foot_point + offset);

			if ( nr_missing_bonds >= 2 ) {
				// coordinate 2
				positions.push_back( foot_point - offset);
			}
		}
	} else if ( nr_known_bonds == 3 ) { // tetrahedral geometry, only missing one bond
		TR << res.atom_name(atom) << ": sp3 atom with 2 known bonds to "<< res.atom(bonded_atoms[1]).name() << ", " << res.atom(bonded_atoms[2]).name()
			<< " and " << res.atom(bonded_atoms[3]).name() << " and " << nr_missing_bonds << " missing. " << std::endl;
		assert( nr_missing_bonds == 1 );
		positions.push_back
			(
			tetrahedral_coordinates(
			position,
			res.atom(bonded_atoms[1]).ideal_xyz(),
			res.atom(bonded_atoms[2]).ideal_xyz(),
			res.atom(bonded_atoms[3]).ideal_xyz(),
			bond_length
			)
		);
	} else {
		TR.Warning << "In core::chemical::modifications::determine_coordinates(), don't know how to handle atom " << res.atom_name(atom) << " with "
			<< nr_known_bonds << " known bonds (and " << nr_missing_bonds << " missing)" << std::endl;
	}

	return positions;
}
/*

//! @brief get the idealized geometry for a particular hybrid orbital type and bond length
//! @param ORBITAL_TYPE the type of orbital associated with the geometry
//! @param BOND_LENGTH the length of the bonds in the geometry
storage::Vector< linal::Vector3D> ValenceHandler::GetIdealizedGeometry
(
const HybridOrbitalType &ORBITAL_TYPE,
const double &BOND_LENGTH
)
{
storage::Vector< linal::Vector3D> positions;

if( ORBITAL_TYPE == GetHybridOrbitalTypes().e_SP3)
{
positions.Resize( 4);
// put the H's at the vertices of a tetrahedron
// Thus, the bond lengths should have coordinates of bond length / sqrt( 3) is
const double bond_length_norm( BOND_LENGTH / std::sqrt( 3.0));
positions( 0) = linal::Vector3D( bond_length_norm, -bond_length_norm, -bond_length_norm);
positions( 1) = linal::Vector3D( -bond_length_norm, bond_length_norm, -bond_length_norm);
positions( 2) = linal::Vector3D( -bond_length_norm, -bond_length_norm, bond_length_norm);
positions( 3) = linal::Vector3D( -bond_length_norm, -bond_length_norm, -bond_length_norm);
}
else if( ORBITAL_TYPE == GetHybridOrbitalTypes().e_SP2)
{
positions.Resize( 3);
// put the vertices at the corners of an equalateral triangle centered at position
positions( 0) = linal::Vector3D( BOND_LENGTH, 0.0, 0.0);
positions( 1) = linal::Vector3D( -BOND_LENGTH / 2.0, BOND_LENGTH * std::sqrt( 3.0) / 2.0, 0.0);
positions( 2) = linal::Vector3D( -BOND_LENGTH / 2.0, -BOND_LENGTH * std::sqrt( 3.0) / 2.0, 0.0);
}
else if( ORBITAL_TYPE == GetHybridOrbitalTypes().e_SP)
{
positions.Resize( 2);
positions( 0) = linal::Vector3D( BOND_LENGTH, 0.0, 0.0);
positions( 1) = linal::Vector3D( -BOND_LENGTH, 0.0, 0.0);
}

return positions;

}
*/





//! point X in B->A->X where A, B and X are on the same line
//! @brief calculates point X in B->A->X where A, B and X are on the same line
//! @param VECTOR_A first point
//! @param VECTOR_B second point
//! @param DISTANCE_XA distance between X and VECTOR_A
//! @return point X in B->A->X where A, B and X are on the same line
numeric::xyzVector<core::Real> linear_coordinates(
	numeric::xyzVector<core::Real> const & vector_a, numeric::xyzVector<core::Real> const & vector_b, core::Real const distance_xa){
	// get the unit vector directed from B to A
	numeric::xyzVector<core::Real> x = (vector_a - vector_b).normalize();
	// extend it the desired distance
	x *= distance_xa;
	// add it to point A
	x += vector_a;
	// return the calculated point
	return x;
}

//! @brief calculates coordinates using angle information (point X in B,C->A->X)
//! @param VECTOR_A first point
//! @param VECTOR_B second point
//! @param VECTOR_C third point
//! @param DISTANCE_XA distance between X and VECTOR_A
//! @param ANGLE_XAB angle between X, VECTOR_A, and VECTOR_B
//! @param ANGLE_XAC angle between X, VECTOR_A, and VECTOR_C
//! @param SIDE true if on a side
//! @param VECTOR_SIDE vector of the side
//! @return point X in B,C->A->X
numeric::xyzVector<core::Real> angle_coordinates(
	const numeric::xyzVector<core::Real> &VECTOR_A,
	const numeric::xyzVector<core::Real> &VECTOR_B,
	const numeric::xyzVector<core::Real> &VECTOR_C,
	const core::Real DISTANCE_XA,
	const core::Real ANGLE_XAB,
	const core::Real ANGLE_XAC,
	const bool SIDE,
	const numeric::xyzVector<core::Real> &VECTOR_SIDE)
{

	// in plane components and direction
	const numeric::xyzVector<core::Real> a( ( VECTOR_A - VECTOR_B).normalize() );
	const numeric::xyzVector<core::Real> b( ( VECTOR_A - VECTOR_C).normalize() );
	const numeric::xyzVector<core::Real> c( std::cos( numeric::constants::d::pi - ANGLE_XAB) * a);
	const numeric::xyzVector<core::Real> d( std::cos( numeric::constants::d::pi - ANGLE_XAC) * b);

	const double bac_angl( projection_angle( c, d));
	const double a_dist( c.norm());
	const double b_dist( d.norm());
	const double c_dist
		(
		std::sqrt
		(
		std::max
		(
		0.0,
		( ( a_dist * a_dist ) + ( b_dist * b_dist) - 2 * a_dist * b_dist * std::cos( bac_angl))
		/ ( std::sin( bac_angl) * std::sin( bac_angl) ) -
		( a_dist * a_dist)
		)
		)
	);

	//in a,b plane
	const numeric::xyzVector<core::Real> e( cross_product( c, d).normalize());
	const numeric::xyzVector<core::Real> f( cross_product( e, c).normalize() * c_dist);
	const numeric::xyzVector<core::Real> g( c + f);

	//orthogonal to a,b plane
	numeric::xyzVector<core::Real> h( e * std::sqrt( std::max( 0.0, 1.0 - g.dot_product(g))));

	//side
	const numeric::xyzVector<core::Real> i( VECTOR_SIDE - VECTOR_A);
	const double dihe_1( dihedral_coordinates( i, numeric::xyzVector<core::Real>(0.0), a, b));
	const double dihe_2( dihedral_coordinates( e, numeric::xyzVector<core::Real>(0.0), a, b));
	if ( dihe_2 && ( ( dihe_1 / dihe_2 > 0.0) != SIDE) ) {
		h *= ( -1.0);
	}

	//sum and length
	const numeric::xyzVector<core::Real> x( numeric::xyzVector<core::Real>( g + h).normalize() * DISTANCE_XA);

	return VECTOR_A + x;
}

//! @brief calculates the unit vector starting from one linal::Vector3D to another
//! @param ORIGIN vector of origin
//! @param TARGET target vector
//! @return the unit vector between ORIGIN and TARGET
numeric::xyzVector<core::Real> triganol_coordinates
(
	const numeric::xyzVector<core::Real> &VECTOR_A,
	const numeric::xyzVector<core::Real> &VECTOR_B,
	const numeric::xyzVector<core::Real> &VECTOR_C,
	const double DISTANCE_XA
)
{
	// get the average unit vector between B->A and C->A
	numeric::xyzVector<core::Real> x( ( VECTOR_A - VECTOR_B).normalize());
	x += ( VECTOR_A -  VECTOR_C).normalize();
	x.normalize();

	// extend the unit vector to the desired length
	x *= DISTANCE_XA;

	// add it to vector A ( the offset)
	x += VECTOR_A;
	return x;
}


//! point X in B,C,D->A->X
//! @brief calculates point X in B,C,D->A->X
//! @param VECTOR_A first point
//! @param VECTOR_B second point
//! @param VECTOR_C third point
//! @param VECTOR_D fourth point
//! @param DISTANCE_XA distance between X and VECTOR_A
//! @return point X in B,C,D->A->X
numeric::xyzVector<core::Real> tetrahedral_coordinates
(
	const numeric::xyzVector<core::Real> &VECTOR_A,
	const numeric::xyzVector<core::Real> &VECTOR_B,
	const numeric::xyzVector<core::Real> &VECTOR_C,
	const numeric::xyzVector<core::Real> &VECTOR_D,
	const double DISTANCE_XA
)
{
	// get the average unit vector between B->A and C->A, D->A
	numeric::xyzVector<core::Real> x( ( VECTOR_A - VECTOR_B).normalize());
	x += ( VECTOR_A - VECTOR_C).normalize();
	x += ( VECTOR_A - VECTOR_D).normalize();
	x.normalize();

	// extend out to the desired distance from A
	x *= DISTANCE_XA;

	// make position relative to A
	x += VECTOR_A;

	return x;
}


//! @brief calculates projection angle between two linal::Vector3D
//! @param VECTOR_A first vector (point)
//! @param VECTOR_B second vector (point)
//! @return projection angle between two linal::Vector3D
double projection_angle( const numeric::xyzVector<core::Real> &VECTOR_A, const numeric::xyzVector<core::Real> &VECTOR_B)
{
	const double projection_angle_cosinus( projection_angle_cosin( VECTOR_A, VECTOR_B));

	// through numerical drift it could be possible that the value is slightly higher or lower than -1 or 1
	// (VECTOR_A * VECTOR_B) / ( Norm( VECTOR_A) * Norm( VECTOR_B)) is the actual cos angle between vectors ab and cd
	if ( projection_angle_cosinus >= double( 1.0) ) {
		// acos( 1) == 0.0
		return 0.0;
	} else if ( projection_angle_cosinus <= double( -1.0) ) {
		// acos( -1) == pi
		return numeric::constants::d::pi;
	}

	// -1 < projection_angle_cosinus < 1, so calling acos( projection_angle_cosinus) yields the proper angle
	return std::acos( projection_angle_cosinus);
}


double projection_angle_cosin( const numeric::xyzVector<core::Real> &VECTOR_A, const numeric::xyzVector<core::Real> &VECTOR_B)
{
	const double projection_angle_cosinus(  VECTOR_A.dot_product( VECTOR_B )  / ( VECTOR_A.norm() * VECTOR_B.norm()  ) );
	return projection_angle_cosinus;
}

//! @brief dihedral angle between four points (A->B -x-> C->D)
//! @brief see http://en.wikipedia.org/wiki/Dihedral_angle for reference
//! @param VECTOR_A first vector (point)
//! @param VECTOR_B second vector (point)
//! @param VECTOR_C third vector (point)
//! @param VECTOR_D fourth vector (point)
//! @return dihedral angle between four points
double dihedral_coordinates
(
	const numeric::xyzVector<core::Real> &VECTOR_A,
	const numeric::xyzVector<core::Real> &VECTOR_B,
	const numeric::xyzVector<core::Real> &VECTOR_C,
	const numeric::xyzVector<core::Real> &VECTOR_D
)
{
	// calculate the two cross products (b1xb2 and b2xb3)
	const numeric::xyzVector<core::Real> cross_b1_b2( cross_product( (VECTOR_B-VECTOR_A), (VECTOR_C - VECTOR_B)));
	const numeric::xyzVector<core::Real> cross_b2_b3( cross_product( (VECTOR_C- VECTOR_B),(VECTOR_D - VECTOR_C) ));

	// calculate the vectors b1 and b2
	const numeric::xyzVector<core::Real> b1( VECTOR_B - VECTOR_A);

	// get the distance b -c

	const double distance( VECTOR_C.distance(VECTOR_B));

	// calculate dihedral angle
	const double dihedral_angle( std::atan2( distance * b1.dot_product(cross_b2_b3), cross_b1_b2.dot_product( cross_b2_b3 ) ));

	return dihedral_angle;
}


}
}
}
