// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packstat/SimplePDB.cc
///
/// @brief
/// @author will sheffler


// Unit header or inline function header
#include <core/scoring/packstat/SimplePDB.hh>
#include <basic/Tracer.hh>
#include <core/scoring/packstat/compute_sasa.hh>

#include <set>

#include <core/scoring/packstat/AtomRadiusMap.hh>
#include <utility/vector1.hh>

static thread_local basic::Tracer TR( "protocols.packstat.SimplePDB" );

namespace core {
namespace scoring {
namespace packstat {

Spheres
SimplePDB::get_spheres(
  AtomRadiusMap const & arm
) const {
  using namespace std;
  set<pair<string,string> > errors;
  Spheres spheres;
  for( SPAtomCIter i = atoms_.begin(); i != atoms_.end(); ++i ) {
    SimplePDB_Atom const & atom( *i );
    numeric::xyzVector<PackstatReal> xyz( atom.x, atom.y, atom.z );
    PackstatReal radius = arm.get_radius( atom.type, atom.res );
    if( radius > 0 ) {
      spheres.push_back( Sphere( xyz, radius ) );
    } else if( radius < 0 ) {
      errors.insert( pair<string,string>(atom.type,atom.res) );
    }
  }
  if( errors.size() ) {
    TR << "ignoring unknown atom types:" << std::endl;
    for(set<pair<string,string> >::iterator i = errors.begin(); i != errors.end(); ++i) {
			TR << "(" << i->first << "," << i->second << "), ";
    }
    TR << std::endl;
    // spheres.clear();
  }
  return spheres;
}


PosePackDataOP
SimplePDB::get_pose_pack_data() const
{
	PosePackDataOP pd( new PosePackData );
	AtomRadiusMap arm;
	pd->spheres = get_spheres(arm);
	pd->centers = get_res_centers();
	pd->labels = res_labels();
	return pd;
}


utility::vector1< numeric::xyzVector<PackstatReal> >
SimplePDB::get_res_centers() const
{
	Size MIN_RES_ATOM_COUNT = 1;
	using namespace std;
	using namespace numeric;
	using namespace utility;
	Size num_res = 0;
	int last_res_num = -12345;
	char last_chain = '*';
	char last_icode = ' ';
	for( SPAtomCIter i = atoms_.begin(); i != atoms_.end(); ++i ) {
		if( i->resnum != last_res_num || i->chain != last_chain  || i->icode != last_icode) {
			num_res++;
			last_res_num = i->resnum;
			last_chain = i->chain;
			last_icode = i->icode;
		}
	}
	vector1< xyzVector<PackstatReal> > centers;
	res_labels_.clear();
	xyzVector<PackstatReal> center(0,0,0);
	last_res_num = -12345;
	last_chain = '*';
	Size res_atom_count = 0;
	for( SPAtomCIter i = atoms_.begin(); i != atoms_.end(); ++i ) {
		if( i->resnum != last_res_num || i->chain != last_chain || i->icode != last_icode) {
			if( num_res > 0 ) {
				if( res_atom_count >= MIN_RES_ATOM_COUNT ){
					centers.push_back(center/res_atom_count);
					res_labels_.push_back(i->res);
				}
				res_atom_count = 0;
			}
			center = xyzVector<PackstatReal>(0,0,0);
			last_res_num = i->resnum;
			last_chain = i->chain;
		}
		center += xyzVector<PackstatReal>( i->x, i->y, i->z );
		res_atom_count++;
	}
	if( res_atom_count >= MIN_RES_ATOM_COUNT ) centers.push_back( center/res_atom_count );
	return centers;
}


core::Size SimplePDB::num_water() const {
	core::Size count = 0;
	for( core::Size i = 1; i <= atoms_.size(); ++i ) {
		if( atoms_[i].res == "HOH" ) {
			count++;
		}
	}
	return count;
}


void SimplePDB::remove_surface_waters() {
	// std::cerr << "remove_surface_waters()" << std::endl;
	using namespace core;
	for( Size i = 1; i <= atoms_.size(); ++i ) {
		atoms_[i].xyz.x(atoms_[i].x);
		atoms_[i].xyz.y(atoms_[i].y);
		atoms_[i].xyz.z(atoms_[i].z);
		atoms_[i].radius = 2; // huge hack, but avoids needing AtomRadiusMap
	}
	int count = 0;
	while( true ) {
		// std::cerr << "REMOVE_SURFACE_WATERS " << count << std::endl;
		count++;
		compute_sasa_generic<SimplePDB_Atom>( atoms_, 3.0 );
		vector1<SimplePDB_Atom> newatoms;
		bool removed = false;
		for( Size i = 1; i <= atoms_.size(); ++i ) {
			// std::cerr << "ATOM " << atoms_[i].res << " " << atoms_[i].sasa << std::endl;
			if( "HOH" != atoms_[i].res || atoms_[i].sasa == 0 ) { // if not water or buried
				newatoms.push_back( atoms_[i] );
			} else {
				removed = true;
			}
		}
		if( !removed ) break;
		atoms_ = newatoms;
	}
}


} // namespace packstat
} // namespace scoring
} // namespace core
