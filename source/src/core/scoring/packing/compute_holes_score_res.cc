// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packing/compute_holes_score_res.cc
/// @brief  Packing Score
/// @author Will Sheffler

//Unit headers

#include <core/scoring/packing/PoseBalls.hh>

//Package headers
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>

#include <core/scoring/packing/compute_holes_score_res.hh>
#include <utility/exit.hh>

#include <basic/prof.hh>

//numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

//utility headers
#include <utility/vector1.hh>

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

#include <core/conformation/Residue.hh>


namespace core {
namespace scoring {
namespace packing {


Real
compute_holes_score_res(
	pose::Pose const & pose,
	PoseBalls const & pb,
	HolesParamsRes const & params
) {
	Real raw_score = 0;

#ifndef WIN32

	using namespace core;

	PROF_START( basic::DALPHABALL );

	std::string cmd = basic::options::option[ basic::options::OptionKeys::holes::dalphaball ]();
	redi::pstream proc( cmd );
	proc.precision(20);
	proc << "NPOINTS" << std::endl << pb.nballs() << std::endl << "COORDS" << std::endl;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		Ball b(pb.ball(i));
		proc << b.x() << " " << b.y() << " " << b.z() << " " << b.r() << " " << std::endl;
	}
	proc << "WEIGHTS" << std::endl;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		core::Size res_num = pb.res_num(i);
		std::string res_name = pb.res_name(i);
		core::Size atom_num = pb.atom_num(i);
		bool skip = false;
		if( res_num <= pose.total_residue() && pose.residue(res_num).is_upper_terminus() )	skip = true;
		if( res_num <= pose.total_residue() && pose.residue(res_num).is_lower_terminus() )	skip = true;
		if( !params.have_params(res_name) ) skip = true;
		if( skip ) {
			proc << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0" << std::endl;
		} else {
			for( Size j = (atom_num-1)*24+5; j <= atom_num*24; j++ ) {
				proc << params.param(res_name)[j] / 12.56637 << " ";
			}
			proc << std::endl;
		}
	}
	proc << "END" << std::endl << redi::peof;

	Size index,ialpha;
	Real val;
	// Real total_weight = 0;
	for( Size a = 1; a <= 20; a++ ) {
		for( Size i = 1; i <= pb.nballs(); i++ ) {
			proc >> ialpha >> index >> val;
			if( a != ialpha || i != index ) {
				std::cerr << "DALPHABALL output indicies not matching! " << i << "!=" << index << " " << a << "!=" << ialpha << std::endl;
				utility_exit_with_message( "DALPHABALL output indicies not matching!" );
			}
			core::Size  res_num = pb. res_num(i);
			core::Size atom_num = pb.atom_num(i);
			std::string res_name = pb.res_name(i);
			bool skip = false;
			if( res_num <= pose.total_residue() && pose.residue(res_num).is_upper_terminus() )	skip = true;
			if( res_num <= pose.total_residue() && pose.residue(res_num).is_lower_terminus() )	skip = true;
			if( !params.have_params(res_name) ) skip = true;
			if( skip ) continue;

			raw_score += val;
			if( a==1 && atom_num==1 ) {
				raw_score -= params.rho( pb.res_name(i) );
				raw_score -= params.rho();
			}
		}
	}

	PROF_STOP( basic::DALPHABALL );

	for( Size i = 1; i <= pb.nballs(); i++ ) {
		std::string res_name = pb.res_name(i);
		Size res_num  = pb.res_num(i);
		Size atom_num = pb.atom_num(i);
		bool skip = false;
		if( res_num <= pose.total_residue() && pose.residue(res_num).is_upper_terminus() )	skip = true;
		if( res_num <= pose.total_residue() && pose.residue(res_num).is_lower_terminus() )	skip = true;
		if( !params.have_params(res_name) ) skip = true;
		if( skip ) continue;

		Real anb5=0,anb10=0,anb15=0,anb20=0;
		for( Size j = 1; j < i; j++ ) {
			Real dis2 = pb.ball(i).xyz().distance_squared(pb.ball(j).xyz());
			if( 400 >= dis2 ) {
				anb20++;
				if( 225 >= dis2 ) anb15++;
				if( 100 >= dis2 ) anb10++;
				if(  25 >= dis2 ) anb5 ++;
			}
		}
		raw_score += params.param(res_name)[24*(atom_num-1)+1] * anb5  / 14.61869 ;
		raw_score += params.param(res_name)[24*(atom_num-1)+2] * anb10 / 116.94952;
		raw_score += params.param(res_name)[24*(atom_num-1)+3] * anb15 / 394.70462;
		raw_score += params.param(res_name)[24*(atom_num-1)+4] * anb20 / 935.59614;
	}

#endif

	return raw_score;
}

Real
compute_holes_deriv_res(
	pose::Pose const & pose,
	PoseBalls const & pb,
	HolesParamsRes const & params,
	core::id::AtomID_Map< numeric::xyzVector<core::Real> > & derivs
) {

	Real raw_score = 0;

#ifndef WIN32

	using namespace core;
	using namespace numeric;
	using id::AtomID;

	PROF_START( basic::DALPHABALL_DERIV );

	for( Size ir = 1; ir <= derivs.size(); ir++ ) {
		for( Size ia = 1; ia <= derivs.n_atom(ir); ia++ ) {
			derivs[AtomID(ia,ir)] = xyzVector<Real>(0.0,0.0,0.0);
		}
	}
	redi::pstream proc( std::string(basic::options::option[ basic::options::OptionKeys::holes::dalphaball ]()) + " DERIV" );
	proc.precision(20);
	proc << "NPOINTS" << std::endl << pb.nballs() << std::endl << "COORDS" << std::endl;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		Ball b(pb.ball(i));
		proc << b.x() << " " << b.y() << " " << b.z() << " " << b.r() << " " << std::endl;
	}
	proc << "WEIGHTS" << std::endl;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		core::Size res_num = pb.res_num(i);
		std::string res_name = pb.res_name(i);
		core::Size atom_num = pb.atom_num(i);
		bool skip = false;
		if( res_num <= pose.total_residue() && pose.residue(res_num).is_upper_terminus() )	skip = true;
		if( res_num <= pose.total_residue() && pose.residue(res_num).is_lower_terminus() )	skip = true;
		if( !params.have_params(res_name) ) skip = true;
		if( skip ) {
			proc << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0" << std::endl;
		} else {
			for( Size j = (atom_num-1)*24+5; j <= atom_num*24; j++ ) {
				proc << params.param(res_name)[j] / 12.56637 << " ";
			}
			proc << std::endl;
		}
	}
	proc << "END" << std::endl << redi::peof;

	Real dx,dy,dz,val;
	for( Size a = 1; a <= 20; a++ ) {
		for( Size i = 1; i <= pb.nballs(); i++ ) {
			proc /*>> ialpha >> index*/ >> val >> dx >> dy >> dz;
			core::Size  res_num = pb. res_num(i);
			core::Size atom_num = pb.atom_num(i);
			std::string res_name = pb.res_name(i);
			bool skip = false;
			if( res_num <= pose.total_residue() && pose.residue(res_num).is_upper_terminus() )	skip = true;
			if( res_num <= pose.total_residue() && pose.residue(res_num).is_lower_terminus() )	skip = true;
			if( !params.have_params(res_name) ) skip = true;
			if( skip ) continue;

			derivs[ id::AtomID(atom_num,res_num) ] += xyzVector<Real>(dx,dy,dz);
			raw_score += val;
			if( a==1 && atom_num==1 ) {
				raw_score -= params.rho( pb.res_name(i) );
				raw_score -= params.rho();
			}
		}
	}

	// xyzVector<Real> mean(0,0,0);
	// Real meanlen = 0.0;
	// Real count = 0.0;
	// for( Size ir = 1; ir <= derivs.size(); ir++ ) {
	// 	for( Size ia = 1; ia <= derivs.n_atom(ir); ia++ ) {
	// 		mean    += derivs[AtomID(ia,ir)];
	// 		meanlen += derivs[AtomID(ia,ir)].length();
	// 		count++;
	// 	}
	// }
	// mean /= count;
	// meanlen /= count;
	// // std::cerr << "meanlen " << meanlen << " mean grad " << mean.x() << " " << mean.y() << " " << mean.z() << std::endl;
	//
	// for( Size ir = 1; ir <= derivs.size(); ir++ ) {
	// 	for( Size ia = 1; ia <= derivs.n_atom(ir); ia++ ) {
	// 		derivs[AtomID(ia,ir)] -= mean;
	// 	}
	// }

	PROF_STOP( basic::DALPHABALL );

#endif

	return raw_score;
}

/// wrappers
Real
compute_holes_score_res(
	pose::Pose const & pose,
	HolesParamsRes const & params
) {
	return compute_holes_score_res(pose,PoseBalls(pose),params);
}

Real
compute_holes_deriv_res(
	pose::Pose const & pose,
	HolesParamsRes const & params,
	core::id::AtomID_Map< numeric::xyzVector<core::Real> > & derivs
) {
	return compute_holes_deriv_res( pose, PoseBalls(pose), params, derivs );
}


}
}
}
