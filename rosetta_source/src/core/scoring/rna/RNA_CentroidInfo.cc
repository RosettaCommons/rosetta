// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_BaseBasePotential.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/rna/RNA_CentroidInfo.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_CentroidInfo.fwd.hh>

// Package headers

// Project headers
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <numeric/xyzMatrix.hh>

#include <utility/vector1.hh>


// AUTO-REMOVED #include <numeric/xyz.functions.hh>

// C++

///////////////////////////////////////////////////////
// Keep track of some base geometry that is
// useful for RNA scoring.
///////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace rna {

/// @details Copy constructors must copy all data, not just some...
RNA_CentroidInfo::RNA_CentroidInfo( RNA_CentroidInfo const & src ) :
	CacheableData()
{
  base_centroids_ = src.base_centroids_;
  base_stubs_ = src.base_stubs_;
  calculated_ = src.calculated_;
}

typedef  numeric::xyzMatrix< Real > Matrix;

///////////////////////////////////////////////////////////////////////////////
// Note that, in principle, the centroid could be a virtual atom,
// and I could let the atom-tree folding machinery figure out where the hell it is.
//
// Also, this copies some code from Phil's dna/base_geometry.cc
//
//Comments (Parin Sep 23 ,2009)...possible problem if every atoms in the nucleotide is virtual...in that case numatoms=0....will this crash the code??

Vector
RNA_CentroidInfo::get_base_centroid( conformation::Residue const & rsd ) const
{
  assert( rsd.is_RNA() );

  Vector centroid( 0.0 );
  Size numatoms = 0;

//	std::cout << "In RNA_CentroidInfo::get_base_centroid: Base atoms" << std::endl; //Parin Sep 23, 2009
  for ( Size i=rsd.first_sidechain_atom()+1; i<= rsd.nheavyatoms(); ++i ) {

/*		/////////////////////////////////////////////////////////////////////Parin Sep 23, 2009/////////////////////////////////////////////////////
		std::cout << "atom " << i  << " " << 	"name= " << rsd.type().atom_name(i) << " type= " << rsd.atom_type(i).name()  << " " << rsd.atom_type_index(i) << " " << rsd.atomic_charge(i);

		if(rsd.is_virtual(i)){
			std::cout << "  Virtual type: Ignore! " << std::endl;
			continue;
		}
		std::cout << std::endl;
*/		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    centroid += rsd.xyz(i);
    numatoms++;
  }
  centroid /= static_cast< Real >( numatoms );

  return centroid;
}



//My local version Sep 23, 2009
/*

Vector
get_base_centroid( conformation::Residue const & rsd , bool verbose=false)
{
  assert( rsd.is_RNA() );

  Vector centroid( 0.0 );
  Size numatoms = 0;

	if(verbose)  std::cout << "Base atoms" << std::endl;
  for ( Size i=rsd.first_sidechain_atom()+1; i<= rsd.nheavyatoms(); ++i ) { //rsd.first_sidechain_atom()+1 to not include the O2star oxygen.

		if(verbose) std::cout << "atom " << i  << " " << 	"name= " << rsd.type().atom_name(i) << " type= " << rsd.atom_type(i).name()  << " " << rsd.atom_type_index(i) << " " << rsd.atomic_charge(i);

		if(rsd.is_virtual(i)){
			if(verbose) std::cout << "  Virtual type: Ignore! " << std::endl;
			continue;
		}
		if(verbose) std::cout << std::endl;

    centroid += rsd.xyz(i);
    numatoms++;
  }
  centroid /= static_cast< Real >( numatoms );
  return centroid;
}
*/
///////////////////////////////////////////////////////////////////////////////
kinematics::Stub
RNA_CentroidInfo::get_base_coordinate_system( conformation::Residue const & rsd, Vector const & centroid ) const
{
  using namespace chemical;
  Size res_type = rsd.aa();

  assert( rsd.is_RNA() );

  Vector x,y,z;

	if ( rsd.has( " CEN" ) ) {
		x = rsd.xyz(" X  ") - rsd.xyz(" CEN" );
		y = rsd.xyz(" Y  ") - rsd.xyz(" CEN" );
		z = cross( x, y );
	} else {

		// Make an axis pointing from base centroid to Watson-Crick edge.
		std::string WC_atom;
		if ( res_type == na_rad ) WC_atom = " N1 ";
		if ( res_type == na_rcy ) WC_atom = " N3 ";
		if ( res_type == na_rgu ) WC_atom = " N1 ";
		if ( res_type == na_ura ) WC_atom = " N3 ";

		Vector const WC_coord (rsd.xyz( WC_atom ) );
		x = WC_coord - centroid;
		x.normalize();

		// Make a perpendicular axis pointing from centroid towards
		// Hoogstein edge (e.g., major groove in a double helix).
		std::string H_atom;
		if ( res_type == na_rad ) H_atom = "N7";
		if ( res_type == na_rcy ) H_atom = "C5";
		if ( res_type == na_rgu ) H_atom = "N7";
		if ( res_type == na_ura ) H_atom = "C5";

		Vector const H_coord (rsd.xyz( H_atom ) );
		y = H_coord - centroid; //not orthonormal yet...
		z = cross(x, y);
		z.normalize(); // Should point roughly 5' to 3' if in a double helix.

		y = cross(z, x);
		y.normalize(); //not necessary but doesn't hurt.

		//  std::cout << "WC : " << WC_coord << "   H : " << H_coord << "    centroid: " << centroid << std::endl;
	}

  return kinematics::Stub( Matrix::cols( x, y, z ), centroid );
}


//////////////////////////////////////////////////////////////////////////////////////
void
RNA_CentroidInfo::initialize_base_centroids_and_stubs( pose::Pose const & pose )
{

  base_centroids_.clear();
  base_stubs_.clear();

  for (Size i = 1; i <= pose.total_residue(); i++ ){
    conformation::Residue const & res_i  ( pose.residue(i) );

    Vector centroid_i( 0.0  );
    kinematics::Stub stub_i;

    if ( res_i.is_RNA() ) {
      centroid_i = get_base_centroid( res_i );
      stub_i     = get_base_coordinate_system(res_i, centroid_i );
    }

    base_centroids_.push_back( centroid_i );
    base_stubs_.push_back( stub_i );

  }

}

//////////////////////////////////////////////////////////////////////////////////////
void
RNA_CentroidInfo::update( pose::Pose const & pose )
{
	initialize_base_centroids_and_stubs( pose );
}


}
}
}
