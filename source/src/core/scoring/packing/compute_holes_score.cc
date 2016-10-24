// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packing/compute_holes_score.cc
/// @brief  Packing Score
/// @author Will Sheffler

#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>  // REQUIRED FOR WINDOWS
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/PoseBalls.hh>
#include <basic/database/open.hh>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <numeric/numeric.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>
#include <utility/numbers.hh>

#ifndef WIN32
#ifndef  __native_client__
#include <pstream.h>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <boost/math/special_functions/fpclassify.hpp>


#endif
#endif

static THREAD_LOCAL basic::Tracer TR( "core.scoring.packing.compute_holes_score" );

namespace core {
namespace scoring {
namespace packing {

/// @details Auto-generated virtual destructor
HolesResult::~HolesResult() {}


void compute_holes_surfs(PoseBalls & pb, std::string ) {

#ifndef WIN32
#ifndef  __native_client__

	using namespace std;
	using namespace basic::options;

	for ( Size tries = 1; tries <= 10; ++tries ) {
		TR << "compute_holes_surfs try: " << tries << std::endl;
		pb.reset_surf();

		std::string cmd = basic::options::option[ OptionKeys::holes::dalphaball ]();

		// build DAlphaBall input commands
		std::ostringstream oss;
		oss << "NPOINTS" << endl << pb.nballs() << endl << "COORDS" << endl;
		for ( Size i = 1; i <= pb.nballs(); i++ ) {
			Ball const & b(pb.ball(i));
			oss << b.x() << " " << b.y() << " " << b.z() << " " << b.r() << " " << endl;
		}
		oss << "END" << endl;

		if ( option[OptionKeys::holes::debug]() ) {
			utility::io::ozstream out("debug_holes_DalphaBall_command.daball");
			out << oss.str();
			out.close();
			utility_exit_with_message("dumped DAlphaBall debug commands");
		}

		// run DAlphaBall and check
		redi::pstream proc( cmd + " alpha20_surf" );
		proc << oss.str() << redi::peof;
		Size nlines = 0;
		string buf, tmp;
		while ( getline(proc,tmp) ) {
			buf += tmp+"\n";
			++nlines;
		}
		proc.close();
		if ( !proc.eof() ) {
			utility_exit_with_message("too much DAlphaBall output!");
		}
		if ( nlines != 20*pb.nballs() ) {
			utility::io::ozstream out("debug_holes_DalphaBall_command.daball");
			out << oss.str();
			out.close();
			utility::io::ozstream out2("debug_holes_DalphaBall_output.log");
			out2 << buf;
			out2.close();
			TR << "nlines: " << nlines << "\npb nballs: " << pb.nballs() << std::endl;
			utility_exit_with_message("incorrect DAlphaBall output! see debug_holes_DalphaBall_command.log");
		}

		std::istringstream iss(buf);
		bool fail = false;
		for ( Size a = 1; a <= 20; a++ ) {
			for ( Size i = 1; i <= pb.nballs(); i++ ) {
				Size index,ialpha;
				std::string val;
				iss >> ialpha >> index >> val;
				Real rval = atof(val.c_str());
				if ( utility::isnan(rval) ) rval = 0.0;
				//TR << "DAlphaBall output index " << ialpha << " " << a << "   " << index << " " << i << std::endl;
				if ( i != index || a != ialpha ) {
					TR << "DAlphaBall output index mismatch " << ialpha << " " << a << "   " << index << " " << i << std::endl;
					fail = true;
					break;
				}
				pb.set_surf(i,a,rval);
			}
			if ( fail ) break;
		}

		if ( !fail ) {
			TR << "compute_holes_surfs completed successfully" << std::endl;
			return;
		} else {
			TR << "failed to run DAlphaBall" << std::endl;
		}
		//  utility_exit_with_message("dab test");
	}
	std::cerr << "Too many compute_holes_surfs failures" << std::endl;
	std::exit(-1);

#endif
#endif
}


HolesResult
compute_rosettaholes_score(
	pose::Pose  const & pose,
	PoseBalls         & pb,
	HolesParams const & resl_params,
	HolesParams const & dec_params,
	HolesParams const & dec15_params,
	bool use_cached_surfs,
	std::string cmd
) {
	if ( cmd == "" ) {
		cmd = basic::options::option[ basic::options::OptionKeys::holes::dalphaball ]();
	}

	HolesResult result;
	if ( !use_cached_surfs ) compute_holes_surfs(pb,cmd);

	core::pose::initialize_atomid_map(result.atom_scores,pose);
	for ( Size i = 1; i <= pb.nballs(); ++i ) {
		Real resl_score = 0.0, dec_score = 0.0, dec15_score = 0.0;
		Size at = pb.atom_type(i);
		char ss = pb.secstruct(i);
		for ( Size a = 1; a <= 20; ++a ) {
			resl_score   += resl_params  .sa_weight(at,ss,a) * pb.surf(i,a);// now in weights / 12.56637;
			dec_score   += dec_params   .sa_weight(at,ss,a) * pb.surf(i,a);// now in weights / 12.56637;
			dec15_score += dec15_params .sa_weight(at,ss,a) * pb.surf(i,a);// now in weights / 12.56637;
		}
		resl_score   += resl_params.nb_weight(at,ss) * pb.smooth_nb(i);// now in weights / 150.0;
		dec_score   += dec_params .nb_weight(at,ss) * pb.smooth_nb(i);// now in weights / 150.0;
		dec15_score += dec15_params .nb_weight(at,ss) * pb.smooth_nb(i);// now in weights / 150.0;
		resl_score   -= resl_params.intercept(at,ss) + resl_params.intercept();
		dec_score   -= dec_params .intercept(at,ss) + dec_params .intercept();
		dec15_score -= dec15_params .intercept(at,ss) + dec15_params .intercept();
		result.resl_score  += resl_score;
		result.decoy_score +=  dec_score;
		result.dec15_score +=  dec15_score;
		Real tmp = 1.0 - (1.0 / (1.0 + exp( 3.768941 * dec_score - 0.5842765 ) ));
		result.atom_scores.set( pb.index_to_id(i), resl_score + 3*tmp );
	}
	result.resl_score  /= pb.nballs();
	result.decoy_score /= pb.nballs();
	result.dec15_score /= pb.nballs();
	result.natom = pb.nballs();

	result.score = 1.0 - (1.0 / (1.0 + exp( 3.768941 * result.decoy_score - 0.5842765 ) ));
	result.score = result.resl_score + 3*result.score;

	TR << "compute_rosettaholes_score done: " << result.score << std::endl;

	return result;
}

Real
compute_dec15_score(
	pose::Pose  const & pose
) {
	std::string cmd = basic::options::option[ basic::options::OptionKeys::holes::dalphaball ]();

	HolesParams dec15_params;
	dec15_params.read_data_file(basic::database::full_name("scoring/rosettaholes/decoy15.params"));
	PoseBalls pb(pose);

	Real dec15_score = 0.0;
	compute_holes_surfs(pb,cmd);

	for ( Size i = 1; i <= pb.nballs(); ++i ) {
		Size at = pb.atom_type(i);
		char ss = pb.secstruct(i);
		for ( Size a = 1; a <= 20; ++a ) {
			dec15_score += dec15_params .sa_weight(at,ss,a) * pb.surf(i,a);// now in weights / 12.56637;
		}
		dec15_score += dec15_params .nb_weight(at,ss) * pb.smooth_nb(i);// now in weights / 150.0;
		dec15_score -= dec15_params .intercept(at,ss) + dec15_params .intercept();
	}
	dec15_score /= pb.nballs();

	TR << "compute_dec15_score done: " << dec15_score << std::endl;

	return dec15_score;
}


HolesResult
compute_holes_score(
	pose::Pose  const & pose,
	PoseBalls         & pb,
	HolesParams const & params,
	bool use_cached_surfs,
	std::string cmd
) {
	if ( cmd == "" ) {
		cmd = basic::options::option[ basic::options::OptionKeys::holes::dalphaball ]();
	}

	HolesResult result;
	if ( !use_cached_surfs ) compute_holes_surfs(pb,cmd);

	core::pose::initialize_atomid_map(result.atom_scores,pose);
	for ( Size i = 1; i <= pb.nballs(); ++i ) {
		Size at = pb.atom_type(i);
		char ss = pb.secstruct(i);
		Real tmp = 0.0;
		for ( Size a = 1; a <= 20; ++a ) {
			tmp += params.sa_weight(at,ss,a) * pb.surf(i,a);// now in weights / 12.56637;
		}
		tmp += params.nb_weight(at,ss) * pb.smooth_nb(i);// now in weights / 150.0;
		tmp -= params.intercept(at,ss) + params.intercept();
		result.score += tmp;
		result.atom_scores.set( pb.index_to_id(i), tmp );
	}
	TR << "compute_holes_score done: " << result.score << std::endl;
	return result;
}

inline core::Real sqr( core::Real x ) {
	return x*x;
}

// smoothed neighbor is between 9 and 11 A
inline core::Real sigmoidish_neighbor( core::Real sqdist ) {
	if ( sqdist >= 121.0 ) {
		return 0.0;
	} else if ( sqdist <= 81.0 ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
		return sqr(1.0 - sqr( (dist - 9.0) / (11.0 - 9.0) ) );
	}
}

// smoothed neighbor is between 9 and 11 A
inline core::Real sigmoidish_neighbor_deriv( core::Real sqdist ) {
	if ( sqdist >= 121.0 ) {
		return 0.0;
	} else if ( sqdist <= 81.0 ) {
		return 0.0;
	} else {
		Real dist = sqrt( sqdist );
		Real x = (dist - 9.0)/ (11.0 - 9.0);
		return -2*x * 2*(1-x*x) / (11.0 - 9.0);
	}
}

core::Real compute_smooth_nb_deriv(
	PoseBalls         & pb,
	HolesParams const & params,
	core::id::AtomID_Map< numeric::xyzVector<core::Real> > & deriv
) {
	using namespace numeric;
	Real tot_snb = 0.0;
	for ( core::Size i = 1; i <= pb.nballs(); i++ ) {
		if ( !pb.is_heavy(i) ) continue;
		deriv.set(pb.index_to_id(i),deriv.get(pb.index_to_id(i)));
		xyzVector<Real> & ixyz( pb.ball(i).xyz() );
		Size at1 = pb.atom_type(i);
		char ss1 = pb.secstruct(i);
		tot_snb += params.nb_weight(at1,ss1);
		for ( core::Size j = 1; j < i; j++ ) {
			if ( !pb.is_heavy(j) ) continue;
			xyzVector<Real> & jxyz( pb.ball(j).xyz() );
			Real d2( ixyz.distance_squared(jxyz) );
			if ( d2 >= 121.0 ) continue;

			Real  sn = sigmoidish_neighbor(d2);
			Real dsn = sigmoidish_neighbor_deriv(d2);
			Size at2 = pb.atom_type(j);
			char ss2 = pb.secstruct(j);
			xyzVector<core::Real> dxyz = ixyz-jxyz;
			dxyz.normalize();
			Real w = 0.0;
			if ( params.have_params(at1,ss1) ) w += params.nb_weight(at1,ss1);
			if ( params.have_params(at2,ss2) ) w += params.nb_weight(at2,ss2);
			tot_snb += sn * w;
			if ( d2 > 81 ) {
				// // TODO: sheffler why 4*?
				deriv[pb.index_to_id(i)] += dxyz*dsn*w;
				deriv[pb.index_to_id(j)] -= dxyz*dsn*w;
			}
		}
	}
	return tot_snb;
}


HolesResult
compute_holes_deriv(
	pose::Pose  const & /*pose*/,
	PoseBalls         & pb,
	HolesParams const & params,
	core::id::AtomID_Map< numeric::xyzVector<core::Real> > & deriv
) {


	HolesResult result;

#ifndef WIN32
#ifndef  __native_client__
	using namespace std;
	using namespace basic::options;

	for ( Size tries = 1; tries <= 10; ++tries ) {
		TR << "compute_holes_deriv try:" << tries << std::endl;
		deriv.clear(numeric::xyzVector<Real>(0.0,0.0,0.0));
		/*Real orig_tot_snb = */compute_smooth_nb_deriv( pb, params, deriv );

		// TESTING SMOOTH NB DERIV
		// Real test_snb = 0.0;
		// for( Size i = 1; i <= pb.nballs(); i++ ) {
		//    Size at = pb.atom_type(i);
		//    char ss = pb.secstruct(i);
		//    // std::cerr << "test " << pb.smooth_nb(i) << " " << params.nb_weight(at,ss) << std::endl;
		//    test_snb = test_snb + pb.smooth_nb(i) * params.nb_weight(at,ss);
		// }
		// std::cerr << "snb: " << test_snb << " " << orig_tot_snb << std::endl;
		//
		// Real DELTA = 0.0001;
		// for( Size i = 1; i < pb.nballs(); i++ ) {
		//    pb.ball(i).x() += DELTA;
		//    core::id::AtomID_Map< numeric::xyzVector<core::Real> > dummy;
		//    Real tot_snb = compute_smooth_nb_deriv( pb, params, dummy );
		//    pb.ball(i).x() -= DELTA;
		//    std::cerr << "dsnb " << (tot_snb-orig_tot_snb)/DELTA << " " << deriv[pb.index_to_id(i)].x() << std::endl;
		// }

		std::string cmd = basic::options::option[ OptionKeys::holes::dalphaball ]();
		redi::pstream proc( cmd + " alpha20_deriv_surf" );
		proc << "NPOINTS" << endl << pb.nballs() << endl << "COORDS" << endl;
		for ( Size i = 1; i <= pb.nballs(); i++ ) {
			Ball const & b(pb.ball(i));
			proc << b.x() << " " << b.y() << " " << b.z() << " " << b.r() << " " << endl;
		}
		proc << "WEIGHTS" << endl;
		for ( Size i = 1; i <= pb.nballs(); i++ ) {
			Size at = pb.atom_type(i);
			char ss = pb.secstruct(i);
			if ( params.have_params(at,ss) ) {
				for ( Size j = 1; j <=20; ++j ) {
					proc << params.sa_weight(at,ss,j) << " ";
				}
				proc << endl;
			} else {
				proc << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0" << endl;
			}
		}
		proc << "END" << endl << redi::peof;

		bool fail = false;
		result.score = 0.0;
		for ( Size i = 1; i <= pb.nballs(); i++ ) {
			Size index;
			Real val,dx,dy,dz;
			proc >> index >> val >> dx >> dy >> dz;
			if ( i != index ) {
				TR << "compute_holes_deriv DAlphaBall output index mismatch " << i << " " << index << std::endl;
				fail = true;
				break;
			}
			Size at = pb.atom_type(i);
			char ss = pb.secstruct(i);
			if ( params.have_params(at,ss) ) {
				val += params.nb_weight(at,ss) * pb.smooth_nb(i);
				val -= params.intercept(at,ss) + params.intercept();
			} else if ( std::abs(val) > 9e-6 ) {
				TR << "no holes params but dalphaball score for atom is not 0!!!! " << val << std::endl;
			}
			result.score += val;
			result.atom_scores.set( pb.index_to_id(i), val );
			deriv[pb.index_to_id(i)] += numeric::xyzVector<core::Real>(dx,dy,dz);
		}
		if ( !fail ) {
			TR << "compute_holes_deriv done: " << result.score << std::endl;
			return result;
		}
		// result.score /= pb.nballs();

		// test that score from alpha20 and alpha20_deriv match
		// HolesResult result2 = compute_holes_score(pb,pose,params);
		// TR << "result1: " << result.score << " result2: " << result2.score << std::endl;
		// for( Size i = 1; i <= pose.size(); ++i ) {
		//    TR << i << " " << result2.atom_scores[core::id::AtomID(1,i)] << " " << result.atom_scores[core::id::AtomID(1,i)] << std::endl;
		// }

		//    pb.compute_smooth_nb();
		// Real base = compute_holes_score(pose,pb,params).score;
		//    Real delta = 0.0001;
		//    for( Size i = 1; i <= 10; ++i ) {
		//       pb.ball(i).x() += delta;
		//       pb.compute_smooth_nb();
		//    Real delta_x = (compute_holes_score(pose,pb,params).score - base) / delta;
		//       pb.ball(i).x() -= delta;
		//       pb.ball(i).y() += delta;
		//       pb.compute_smooth_nb();
		//    Real delta_y = (compute_holes_score(pose,pb,params).score - base) / delta;
		//       pb.ball(i).y() -= delta;
		//       pb.ball(i).z() += delta;
		//       pb.compute_smooth_nb();
		//    Real delta_z = (compute_holes_score(pose,pb,params).score - base) / delta;
		//       pb.ball(i).z() -= delta;
		//
		//       std::cerr << params.intercept()
		//          << " " << deriv[pb.index_to_id(i)].x() << " " << delta_x << " "
		//          << " " << deriv[pb.index_to_id(i)].y() << " " << delta_y << " "
		//          << " " << deriv[pb.index_to_id(i)].z() << " " << delta_z << endl;
		//    }
	}
	std::cerr << "too many dalphaball deriv fails!" << std::endl;
	std::exit(-1);

#endif
#endif

	return result;

}


HolesResult
compute_rosettaholes_score(
	pose::Pose const & pose,
	PoseBalls        & pb
) {
	HolesParams hp_resl, hp_dec, hp_dec15;
	hp_resl .read_data_file(basic::database::full_name("scoring/rosettaholes/resl.params"));
	hp_dec  .read_data_file(basic::database::full_name("scoring/rosettaholes/decoy25.params"));
	hp_dec15.read_data_file(basic::database::full_name("scoring/rosettaholes/decoy15.params"));
	return compute_rosettaholes_score(pose,pb,hp_resl,hp_dec,hp_dec15);
}


HolesResult
compute_holes_deriv(
	pose::Pose  const & pose,
	HolesParams const & params,
	core::id::AtomID_Map< numeric::xyzVector<core::Real> > & deriv
) {
	PoseBalls pb(pose);
	return compute_holes_deriv( pose, pb, params, deriv );
}

HolesResult
compute_rosettaholes_score(
	pose::Pose  const & pose,
	HolesParams const & resl_params,
	HolesParams const & dec_params,
	HolesParams const & dec15_params
) {
	PoseBalls pb(pose);
	return compute_rosettaholes_score( pose, pb, resl_params, dec_params, dec15_params );
}

HolesResult
compute_rosettaholes_score(
	pose::Pose const & pose
) {
	PoseBalls pb(pose);
	return compute_rosettaholes_score( pose, pb );
}

HolesResult
compute_holes_score(
	pose::Pose  const & pose,
	HolesParams const & params
) {
	PoseBalls pb(pose);
	return compute_holes_score(pose,pb,params);
}

HolesResult
compute_holes_score( pose::Pose const & pose, std::string const & cmd )
{
	PoseBalls pb(pose);
	HolesParams hp_resl, hp_dec, hp_dec15;
	hp_resl .read_data_file(basic::database::full_name("scoring/rosettaholes/resl.params"));
	hp_dec  .read_data_file(basic::database::full_name("scoring/rosettaholes/decoy25.params"));
	hp_dec15.read_data_file(basic::database::full_name("scoring/rosettaholes/decoy15.params"));
	bool use_cached_surfs( false );
	return compute_rosettaholes_score( pose, pb, hp_resl, hp_dec, hp_dec15, use_cached_surfs, cmd );
}


} // packing
} // scoring
} // core
