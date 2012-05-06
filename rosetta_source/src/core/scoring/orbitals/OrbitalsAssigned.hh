// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/orbitals/OrbitalsAssigned.hh
/// @brief  Created on: Jun 2, 2010
/// @author combss

#ifndef INCLUDED_core_scoring_orbitals_ORBITALSASSIGNED_HH
#define INCLUDED_core_scoring_orbitals_ORBITALSASSIGNED_HH

#include <utility/vector1.hh>
// AUTO-REMOVED #include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/types.hh>

//Auto Headers
#include <core/conformation/Residue.fwd.hh>
#include <string>


namespace core{
namespace scoring{
namespace orbitals{

class OrbitalsAssigned {

public:

	utility::vector1< numeric::xyzVector< core::Real > > get_lp_xyz(
			core::conformation::Residue const & residue
	);

	utility::vector1< numeric::xyzVector<core::Real> >  CoordinatesDihedral(
			std::string atom1,
			std::string atom2,
			std::string atom3,
			core::Real distance_xa,
			core::conformation::Residue const & residue
	);

	utility::vector1< numeric::xyzVector<core::Real> >  CoordinatesTetrahedral(
			std::string atom1,
			std::string atom2,
			std::string atom3,
			core::Real distance,
			core::conformation::Residue const & residue

	);

	utility::vector1< numeric::xyzVector<core::Real> > aromatic_ring_center(
			numeric::xyzVector<core::Real> vector_d,
			numeric::xyzVector<core::Real> vector_f,
			core::conformation::Residue const & residue,
			core::Real dist

);

	utility::vector1< numeric::xyzVector<core::Real> > cp_function(
			std::string atomtype,
			numeric::xyzVector<core::Real> vector_d,
			numeric::xyzVector<core::Real> vector_f,
			core::conformation::Residue const & residue,
			core::Real dist

	);

	utility::vector1< std::pair< numeric::xyzVector<core::Real>,  std::string > >  get_hydrogens(
			core::conformation::Residue const & resid1
	);

private:


};//OrbitalsAssigned



}//namespace orbitals
}//namespace scoring
}//namespace core



#endif /* ORBITALSASSIGNED_HH_ */
