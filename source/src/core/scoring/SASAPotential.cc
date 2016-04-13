// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/SASAPotential.cc
/// @brief
/// @author Jim Havranek

// Project headers
#include <core/scoring/SASAPotential.hh>
#include <core/scoring/power_diagram/PowerDiagram.hh>

#include <core/scoring/Energies.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
//#include <core/scoring/etable/count_pair/CountPairIntraResC3.hh>
//#include <core/scoring/etable/count_pair/CountPair1BC3.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/AtomTree.hh>
// AUTO-REMOVED #include <core/kinematics/DomainMap.hh>
//#include <core/pack/task/PackerTask.hh>
#include <basic/prof.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>

#include <numeric/constants.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.io.hh>
#include <numeric/xyz.functions.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <boost/lexical_cast.hpp>

#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <map>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.SASAPotential" );

typedef numeric::xyzVector< core::Real > Vector;
typedef numeric::xyzMatrix< core::Real > Matrix;

using namespace core::scoring::power_diagram;

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {

bool
SASAShouldItCount(
	conformation::Residue const & rsd,
	Size const & atm
)
{
	//  if( rsd.atom_name( atm ) == " CB " ||  rsd.atom_name( atm ) == " O  " )  {
	//  if( rsd.atom_name( atm ) == " CB " )  {
	//   return true;
	//  }


	if ( rsd.seqpos() == 2 ) return false;

	//  if( rsd.atom_name( atm ) == " C  " ||  rsd.atom_name( atm ) == " O  " ||  rsd.atom_name( atm ) == " CB " ) {
	if ( rsd.atom_name( atm ) == " CB " ) {
		return true;
	} else {
		return false;
	}

}














////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
SASAPotential::setup_for_scoring(
	pose::Pose & pose
) const
{
	// TR << "in setup_for_scoring use_nblist is " << pose.energies().use_nblist() << std::endl;

	pd_->construct_from_pose( pose );

	// TR << "Exiting setup_for_scoring" << std::endl;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
SASAPotential::get_res_res_sasa(
	Residue const & rsd1,
	Residue const & rsd2
) const
{
	using namespace etable::count_pair;

	Real sasaE( 0.0 );
	Size natoms1 = rsd1.natoms();
	// Size natoms2 = rsd2.natoms();

	Size const res1( rsd1.seqpos() );

	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );

	// As a convention, surface are is strictly considered intrares
	if ( !same_res ) { return 0.0; }

	for ( Size atm1 = 1 ; atm1 <= natoms1 ; ++atm1 ) {
		if ( rsd1.is_virtual( atm1 ) ) continue;

		// Get the surface area for this atom
		Real const atom_atomE( pd_->extract_sasa_for_atom( res1, atm1 ) );
		sasaE += atom_atomE;
	}


	return sasaE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
SASAPotential::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & , // provides context
	Real const & factor,
	utility::vector1< DerivVectorPair > & r1_at_derivs,
	utility::vector1< DerivVectorPair > & r2_at_derivs
) const
{
	using namespace etable::count_pair;

	Size natoms1 = rsd1.natoms();
	Size natoms2 = rsd2.natoms();

	Size res1( rsd1.seqpos() );
	Size res2( rsd2.seqpos() );
	bool const same_res( res1 == res2 );

	// assert( pose.energies().use_nblist() );
	assert( r1_at_derivs.size() >= rsd1.natoms() );
	assert( r2_at_derivs.size() >= rsd2.natoms() );

	// Real SASA_check( 0.0 );

	for ( Size atm1 = 1 ; atm1 <= natoms1 ; ++atm1 ) {
		if ( rsd1.is_virtual( atm1 ) ) continue;

		PDsphereOP patom1( pd_->sphere_lookup( res1, atm1 ) );
		//  std::list< power_diagram::PDinterOP > intersections1( pd_->get_intersections_for_atom( res1, atm1 ) );

		//  for( std::list< power_diagram::PDinterOP >::iterator itr = intersections1.begin() ;
		//   itr != intersections1.end() ; ++itr ) {
		//   find_common_intersection_atoms( (*itr) );
		//  }

		//  utility::vector1< utility::vector1< power_diagram::SAnode > >
		//    cycles1( get_cycles_from_intersections( intersections1, patom1 ) );

		//  if( !SASAShouldItCount( rsd1, atm1 ) ) continue;

		Size const start_atom( same_res ? atm1 + 1 : 1 );
		for ( Size atm2 = start_atom ; atm2 <= natoms2 ; ++atm2 ) {
			if ( rsd2.is_virtual( atm2 ) ) continue;

			//   if( same_res && ( atm1 == atm2 ) ) { continue; }

			PDsphereOP patom2( pd_->sphere_lookup( res2, atm2 ) );
			//   std::list< power_diagram::PDinterOP > intersections2( pd_->get_intersections_for_atom( res2, atm2 ) );

			//   for( std::list< power_diagram::PDinterOP >::iterator itr = intersections2.begin() ;
			//    itr != intersections2.end() ; ++itr ) {
			//    find_common_intersection_atoms( (*itr) );
			//   }

			//   utility::vector1< utility::vector1< power_diagram::SAnode > >
			//     cycles2( get_cycles_from_intersections( intersections2, patom2 ) );


			//   if( !SASAShouldItCount( rsd2, atm2 ) ) continue;

			// Figure out if atom 2 is involved with any of the arcs for atom 1

			Vector f1( 0.0 );
			Vector f2( 0.0 );

			//    TR << "Getting derivs for res1 " << res1 << " atom1 " << atm1 << " and res2 " << res2 << " atom2 " << atm2 << std::endl;

			get_derivs_from_cycles( patom1->cycles(), patom1, patom2, f1, f2 );
			r1_at_derivs[ atm1 ].f1() += factor*f1;
			r1_at_derivs[ atm1 ].f2() += factor*f2;
			r2_at_derivs[ atm2 ].f1() -= factor*f1;
			r2_at_derivs[ atm2 ].f2() -= factor*f2;

			f1 = 0.0;
			f2 = 0.0;

			get_derivs_from_cycles( patom2->cycles(), patom2, patom1, f1, f2 );
			r1_at_derivs[ atm1 ].f1() -= factor*f1;
			r1_at_derivs[ atm1 ].f2() -= factor*f2;
			r2_at_derivs[ atm2 ].f1() += factor*f1;
			r2_at_derivs[ atm2 ].f2() += factor*f2;


		}


		//  Real const atom_atomE( get_sasa_from_cycles( cycles, patom1 ) );
		//  SASA_check += atom_atomE;

	}
	// TR << "SASA for residue " << res1 << " is " << SASA_check << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace scoring
} // namespace core
