// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packstat/SimplePDB.hh
///
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_packstat_SimplePDB_hh
#define INCLUDED_core_scoring_packstat_SimplePDB_hh

// Project forward headers
#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/SimplePDB.fwd.hh>
#include <core/scoring/packstat/SimplePDB_Atom.hh>


#include <core/scoring/packstat/AtomRadiusMap.fwd.hh>


namespace core {
namespace scoring {
namespace packstat {

typedef utility::vector1<SimplePDB_Atom> SPAtoms;
typedef utility::vector1<SimplePDB_Atom>::iterator SPAtomIter;
typedef utility::vector1<SimplePDB_Atom>::const_iterator SPAtomCIter;

/// @brief
class SimplePDB
{

	friend std::istream & operator>> ( std::istream & in , SimplePDB       & pdb  );
	friend std::ostream & operator<< ( std::ostream & out, SimplePDB const & pdb  );

public: // Creation

	SimplePDB() {}
	~SimplePDB() {}

	Spheres get_spheres( AtomRadiusMap const & arm ) const;

	utility::vector1< numeric::xyzVector<PackstatReal> > get_res_centers() const;

	void remove_surface_waters();

	core::Size num_water() const;

	utility::vector1<std::string> & res_labels() const { return res_labels_; }

	PosePackDataOP get_pose_pack_data() const;

private: // fields

	utility::vector1<SimplePDB_Atom> atoms_;

public: // accessors

	std::size_t natom() const { return atoms_.size(); }
	utility::vector1<SimplePDB_Atom>       & atoms()       { return atoms_; }
	utility::vector1<SimplePDB_Atom> const & atoms() const { return atoms_; }
	SimplePDB_Atom       & atom( std::size_t const i )       { return atoms_[i]; }
	SimplePDB_Atom const & atom( std::size_t const i ) const { return atoms_[i]; }
	mutable utility::vector1<std::string> res_labels_;


}; // SimplePDB


} // namespace packstat
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_packstat_SimplePDB_HH
