// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/GridHash3D.cc
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#ifndef INCLUDED_protocols_ligand_docking_GALigandDock_GridHash3D_hh
#define INCLUDED_protocols_ligand_docking_GALigandDock_GridHash3D_hh

#include <core/types.hh>
#include <core/conformation/Atom.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/lkball/LK_BallInfo.hh>

#include <unordered_map>


#if defined(B0)  // work around B0 define in termios.h on Mac
	#undef B0
#endif

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {


struct atmInfo {
	numeric::xyzVector<core::Real>
	nbr_atom() { return atom.xyz(); }

	core::conformation::Atom atom;
	core::Real atomicCharge;

	core::Size nAtomWaters;
	core::scoring::lkball::WaterCoords atomWaters;
};

struct hbDon {
	numeric::xyzVector<core::Real>
	nbr_atom() { return D; }

	numeric::xyzVector< core::Real > D,H;
	core::scoring::hbonds::HBDonChemType dontype;
};

struct hbAcc {
	numeric::xyzVector<core::Real>
	nbr_atom() { return A; }

	numeric::xyzVector< core::Real > A,B,B_0;
	core::scoring::hbonds::HBAccChemType acctype;
};


template <class T>
class GridHash3D {
public:
	GridHash3D() {
		resolution_ = 5.0;
	}

	void
	set_resolution( core::Real reso ) {
		resolution_ = reso;
	}

	void
	add_point(
		T info
	);

	void
	get_neighbors(
		numeric::xyzVector<core::Real> const &queryPoint,
		utility::vector1<T> &neighbors,
		core::Real maxDis=7.0
	);

	int
	npoints() { return nodes_.size(); }

private:
	std::unordered_multimap< int, T > nodes_;
	core::Real resolution_;
};

template <class T>
void
GridHash3D<T>::add_point( T info ) {

	int x,y,z;
	x = (int) std::floor( (info.nbr_atom()[0]+0.5)/resolution_ );
	x = x%1024; if ( x<0 ) x+=1024;
	y = (int) std::floor( (info.nbr_atom()[1]+0.5)/resolution_ );
	y = y%1024; if ( y<0 ) y+=1024;
	z = (int) std::floor( (info.nbr_atom()[2]+0.5)/resolution_ );
	z = z%1024; if ( z<0 ) z+=1024;

	int hashid = (x<<20)+(y<<10)+z;

	nodes_.insert(std::make_pair(hashid, info));
}


template <class T>
void
GridHash3D<T>::get_neighbors(
	numeric::xyzVector<core::Real> const &queryPoint,
	utility::vector1<T> &neighbors,
	core::Real maxDis/*=7.0*/ )
{
	int x,y,z;
	x = (int) std::floor( (queryPoint[0]+0.5)/resolution_ );
	y = (int) std::floor( (queryPoint[1]+0.5)/resolution_ );
	z = (int) std::floor( (queryPoint[2]+0.5)/resolution_ );

	int boundcheck = std::ceil( maxDis/resolution_ );

	neighbors.clear();

	core::Real maxDis2 = maxDis*maxDis;
	for ( int ix=x-boundcheck; ix<=x+boundcheck; ++ix ) {
		int xx = ix%1024; if ( xx<0 ) xx+=1024;
		for ( int iy=y-boundcheck; iy<=y+boundcheck; ++iy ) {
			int yy = iy%1024; if ( yy<0 ) yy+=1024;
			for ( int iz=z-boundcheck; iz<=z+boundcheck; ++iz ) {
				int zz = iz%1024; if ( zz<0 ) zz+=1024;

				int hashid = (xx<<20)+(yy<<10)+zz;
				auto it_cell = nodes_.equal_range( hashid );
				for ( auto it=it_cell.first; it!=it_cell.second; it++ ) {
					core::Real distAct2 = (queryPoint - it->second.nbr_atom()).length_squared();
					if ( distAct2 <= maxDis2 ) {
						neighbors.push_back( it->second );
					}
				}
			}
		}
	}
}
}
}
}

#endif
