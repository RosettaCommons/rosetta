// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packing/surf_vol.cc
/// @brief  Packing Score
/// @author Will Sheffler

//Unit headers
#include <core/scoring/packing/surf_vol.hh>

//Package headers

#include <basic/options/option.hh>

#include <basic/prof.hh>

//numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/xyzMatrix.hh>

//utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#ifdef  __native_client__
#define WIN32
#endif

#ifndef WIN32
#include <pstream.h>
#endif

//C++ headers
#include <iostream>
#include <iomanip>

// option key includes

#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/packing/PoseBalls.hh>
#include <core/scoring/packing/PoseBallsLite.hh>


namespace core {
namespace scoring {
namespace packing {


Real
get_surf_tot(
	pose::Pose const & pose,
	core::Real const   probe_radius
) {
	using namespace core;

	Real tot_surf = 0;

#ifndef WIN32
	PROF_START( basic::DALPHABALL );
	PoseBallsLite pb(pose);

	std::string cmd = basic::options::option[ basic::options::OptionKeys::holes::dalphaball ]();
	redi::pstream proc( cmd + " totalsurf" );
	// proc.precision(20);
	proc << "NPOINTS" << std::endl << pb.nballs() << std::endl << "COORDS" << std::endl;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		Ball b(pb.ball(i));
		proc << b.x() << " " << b.y() << " " << b.z() << " " << b.r() + probe_radius << " " << std::endl;
	}
	// no weights in surf_vol mode
	proc << "END" << std::endl << redi::peof;

	proc >> tot_surf;
	// std::string s;
	// proc >> s; std::cerr << s << std::endl;
	// proc >> s; std::cerr << s << std::endl;
	// proc >> s; std::cerr << s << std::endl;
	// proc >> s; std::cerr << s << std::endl;
	// Size index;
	// Real s,v;
	// for( Size i = 1; i <= pb.nballs(); i++ ) {
	// 	proc >> index >> s >> v;
	// 	if( i != index ) {
	// 		std::cerr << "DALPHABALL output indicies not matching! " << i << "!=" << index << std::endl;
	// 		std::exit(-1);
	// 	}
	// 	tot_surf += s;
	// }

	PROF_STOP( basic::DALPHABALL );

#endif
	return tot_surf;
}


SurfVol
get_surf_vol(
	pose::Pose const & pose,
	core::Real const   probe_radius
) {
	core::id::AtomID_Mask atoms;
	core::pose::initialize_atomid_map_heavy_only(atoms,pose,true);
	return get_surf_vol(pose,atoms,probe_radius);
}

SurfVol
get_surf_vol(
	pose::Pose const & pose,
	core::id::AtomID_Mask const & whichatoms,
	core::Real const   probe_radius
) {
	using namespace core;

	SurfVol result;

#ifndef WIN32

	PoseBallsLite pb(pose,whichatoms);

	initialize_AtomID_Map<Real>(result.surf,pb);
	initialize_AtomID_Map<Real>(result.vol ,pb);

	PROF_START( basic::DALPHABALL );

	std::string cmd = basic::options::option[ basic::options::OptionKeys::holes::dalphaball ]();
	redi::pstream proc( cmd + " surf_vol" );
	// proc.precision(20);
	proc << "NPOINTS" << std::endl << pb.nballs() << std::endl << "COORDS" << std::endl;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		Ball b(pb.ball(i));
		proc << b.x() << " " << b.y() << " " << b.z() << " " << b.r() + probe_radius << " " << std::endl;
	}
	// no weights in surf_vol mode
	proc << "END" << std::endl << redi::peof;

	result.tot_surf = 0.0;
	result.tot_vol  = 0.0;
	Size index;
	Real s,v;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		proc >> index >> s >> v;
		if( i != index ) {
			std::cerr << "DALPHABALL output indicies not matching! " << i << "!=" << index << std::endl;
			std::exit(-1);
		}
		result.surf[ pb.index_to_id(i) ] = s;
		result.vol [ pb.index_to_id(i) ] = v;
		result.tot_surf += s;
		result.tot_vol  += v;
	}

	PROF_STOP( basic::DALPHABALL );

#endif

	return result;
}

SurfVolDeriv
get_surf_vol_deriv(
	pose::Pose const & pose,
	core::Real const   probe_radius
) {
	core::id::AtomID_Mask atoms;
	core::pose::initialize_atomid_map_heavy_only(atoms,pose,true);
	return get_surf_vol_deriv(pose,atoms,probe_radius);
}

SurfVolDeriv
get_surf_vol_deriv(
	pose::Pose const & pose,
	core::id::AtomID_Mask const & whichatoms,
	core::Real const   probe_radius
) {
	using namespace core;

	SurfVolDeriv result;

#ifndef WIN32

	PROF_START( basic::DALPHABALL );

	PoseBalls pb(pose,whichatoms);

	initialize_AtomID_Map<Real>(result.surf,pb);
	initialize_AtomID_Map<Real>(result.vol ,pb);
	initialize_AtomID_Map<numeric::xyzVector<Real> >(result.dsurf,pb);
	initialize_AtomID_Map<numeric::xyzVector<Real> >(result.dvol ,pb);

	std::string cmd = basic::options::option[ basic::options::OptionKeys::holes::dalphaball ]();
	redi::pstream proc( cmd + " surf_vol_deriv" );
	proc.precision(20);
	proc << "NPOINTS" << std::endl << pb.nballs() << std::endl << "COORDS" << std::endl;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		Ball b(pb.ball(i));
		proc << b.x() << " " << b.y() << " " << b.z() << " " << b.r() + probe_radius << " " << std::endl;
	}
	// no weights in surf_vol mode
	proc << "END" << std::endl << redi::peof;

	result.tot_surf = 0.0;
	result.tot_vol  = 0.0;
	Size index;
	Real s,v,dsx,dsy,dsz,dvx,dvy,dvz;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		proc >> index >> s >> v >> dsx >> dsy >> dsz >> dvx >> dvy >> dvz;
		if( i != index ) {
			std::cerr << "DALPHABALL output indicies not matching! " << i << "!=" << index << std::endl;
			std::exit(-1);
		}
		result.surf[ pb.index_to_id(i) ] = s;
		result.vol [ pb.index_to_id(i) ] = v;
		result.dsurf[ pb.index_to_id(i) ] = numeric::xyzVector<Real>(dsx,dsy,dsz);
		result.dvol [ pb.index_to_id(i) ] = numeric::xyzVector<Real>(dvx,dvy,dvz);
		result.tot_surf += s;
		result.tot_vol  += v;
	}

	PROF_STOP( basic::DALPHABALL );

#endif

	return result;
}


}
}
}
