// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class for reading in the orbital type properties
///
/// @details
/// This class contains the ORBITALS INTERNAL_ICOOR data that is read in from residue_io.cc. Actually,
/// the data is set when residue_io.cc calls the command from residuetype.cc set_orbital_icoor. The data
/// is set and chills out in memory until needed. The actual xyz coordinates are not stored at this point.
/// xyz coordinates are constructed when conformation/Residue.hh asks for the build function in this class.
/// At that point, the coordinates are built and returned.
///
/// But wait, you say, why do you store the names of the atoms instead of the index of the atoms!? Well, the
/// problem occurs when residuetype reorders the indices of the atoms. When this occurrs, the indices for the icoor
/// are not reordered here. Another problem ocurs because orbital indices are not reordered in this process because
/// orbital indices are seperate from the atom indices. Regardless, when you build the xyz coords, this step is transparent
/// because the function orbital_xyz() in residue.hh takes care of this conversion of indice to string.
///
/// @note NOTE!!!!!!!!!!! The internal coordinates cannot contain an orbital as the stub1, stub2, or stub3 atom.
/// This is because the xyz coordinates are not updated when the conformation changes. The stub1, stub2, stub2 atoms
/// must be actual atoms and not orbitals!!!!! (design feature or flaw? you decide)
///
///
/// @author
/// Steven Combs
///
///
/////////////////////////////////////////////////////////////////////////


#include <core/chemical/orbitals/ICoorOrbitalData.hh>

//#include <core/chemical/ResidueType.hh>
//#include <core/conformation/Residue.hh>

// Utility headers
//#include <utility/exit.hh>
//#include <numeric/xyz.functions.hh>
#include <core/kinematics/Stub.hh>
//#include <core/kinematics/Jump.hh>
//#include <utility/vector1.hh>

#ifdef    SERIALIZATION
#include <core/chemical/ResidueGraphTypes.srlz.hh>
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#endif // SERIALIZATION

namespace core {
namespace chemical {
namespace orbitals {

/// @brief construct ICoorOrbitalData.
ICoorOrbitalData::ICoorOrbitalData():
	phi_(0.0),
	theta_(0.0),
	distance_(0.0),
	stub1_(0),
	stub2_(0),
	stub3_(0),
	vertex1_(0),
	vertex2_(0),
	vertex3_(0)
{}


//Testing size orbitals
ICoorOrbitalData::ICoorOrbitalData(
	Real phi,
	Real theta,
	Real distance,
	Size stub1,
	Size stub2,
	Size stub3,
	VD vertex1,
	VD vertex2,
	VD vertex3
):
	phi_(phi),
	theta_(theta),
	distance_(distance),
	stub1_(stub1),
	stub2_(stub2),
	stub3_(stub3),
	vertex1_(vertex1),
	vertex2_(vertex2),
	vertex3_(vertex3)
{

}

void
ICoorOrbitalData::remap_atom_vds( std::map< VD, VD > const & old_to_new ) {
	if ( old_to_new.count( vertex1_ ) == 1 ) {
		vertex1_ = old_to_new.at( vertex1_ );
	}
	if ( old_to_new.count( vertex2_ ) == 1 ) {
		vertex2_ = old_to_new.at( vertex2_ );
	}
	if ( old_to_new.count( vertex3_ ) == 1 ) {
		vertex3_ = old_to_new.at( vertex3_ );
	}
}


/// @brief return the phi for a given orbital
Real ICoorOrbitalData::phi() const
{
	return phi_;
}

/// @brief return the theta for a given orbital
Real ICoorOrbitalData::theta() const
{
	return theta_;
}

/// @brief return the distance for a given orbital
Real ICoorOrbitalData::distance() const
{
	return distance_;
}


/// @brief build the xyz coordinates for an orbital based upon the stub1, stub2, stub3 xyz coordinates.
/// @note NOTE!!!!!!!!!!! The internal coordinates cannot contain an orbital as the stub1, stub2, or stub3 atom.
/// This is because the xyz coordinates are not updated when the conformation changes. The stub1, stub2, stub2 atoms
/// must be actual atoms and not orbitals!!!!11111!!!!!!!11111!(design feature or flaw? you decide)
Vector
ICoorOrbitalData::build(
	Vector stub1_xyz,
	Vector stub2_xyz,
	Vector stub3_xyz) const
{
	debug_assert( kinematics::Stub( stub1_xyz,
		stub2_xyz,
		stub3_xyz).is_orthogonal( 0.001 ) );

	return kinematics::Stub(stub1_xyz, stub2_xyz, stub3_xyz).spherical(phi_, theta_, distance_);
}


}
}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::orbitals::ICoorOrbitalData::save( Archive & arc ) const {
	arc( CEREAL_NVP( phi_ ) ); // Real
	arc( CEREAL_NVP( theta_ ) ); // Real
	arc( CEREAL_NVP( distance_ ) ); // Real
	arc( CEREAL_NVP( stub1_ ) ); // Size
	arc( CEREAL_NVP( stub2_ ) ); // Size
	arc( CEREAL_NVP( stub3_ ) ); // Size
	SERIALIZE_VD( arc, vertex1_ ); // EXEMPT vertex1_
	SERIALIZE_VD( arc, vertex2_ ); // EXEMPT vertex2_
	SERIALIZE_VD( arc, vertex3_ ); // EXEMPT vertex3_
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::orbitals::ICoorOrbitalData::load( Archive & arc ) {
	arc( phi_ ); // Real
	arc( theta_ ); // Real
	arc( distance_ ); // Real
	arc( stub1_ ); // Size
	arc( stub2_ ); // Size
	arc( stub3_ ); // Size
	DESERIALIZE_VD( arc, vertex1_ ); // EXEMPT vertex1_
	DESERIALIZE_VD( arc, vertex2_ ); // EXEMPT vertex2_
	DESERIALIZE_VD( arc, vertex3_ ); // EXEMPT vertex3_
}
SAVE_AND_LOAD_SERIALIZABLE( core::chemical::orbitals::ICoorOrbitalData );
#endif // SERIALIZATION
