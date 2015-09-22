// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
//#include <basic/options/option.hh>

#include <protocols/packstat/types.hh>
#include <protocols/packstat/SimplePDB_Atom.hh>
#include <protocols/packstat/SimplePDB.hh>
#include <protocols/packstat/io.hh>
#include <protocols/packstat/AtomRadiusMap.hh>
#include <protocols/packstat/compute_sasa.hh>
#include <protocols/packstat/sasa_dot_locations.hh>
#include <protocols/packstat/packing_score_params.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/options/option.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/sasa.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctype.h>

#include <time.h>

// option key includes

#include <basic/options/keys/packstat.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <numeric/NumericTraits.hh>


using core::Real;

static THREAD_LOCAL basic::Tracer TRps( "packstat" );

using namespace core::scoring::packstat;

// void test_io() {
//   using namespace core::scoring::packstat;
//   using namespace std;
//
//   istringstream iss("ATOM     12  N   GLU A   2     -13.565  31.875  -5.182  1.00 51.33           N\n");
//   SimplePDB_Atom atom;
//   iss >> atom;
//   cout << atom << endl;
//
//   SimplePDB pdb;
//   ifstream in("test_in.pdb");
//   in >> pdb;
//
//   cout << pdb;
//
//
// }
//
// void test_spheres(std::string fname) {
//   using namespace core::scoring::packstat;
//   using namespace std;
//
//   AtomRadiusMap arm;
//   SimplePDB pdb;
//   utility::io::izstream in(fname.c_str());
//   in >> pdb;
//
// 	// for( SPAtomIter i = pdb.atoms().begin(); i != pdb.atoms().end(); ++i ) {
// 	// 	TRps << *i << " " << arm.get_radius(i->type,i->res) << std::endl;
// 	// }
//
// 	// TRps << fname << std::endl;
//   Spheres spheres( pdb.get_spheres(arm) );
//   for( SphereIter i = spheres.begin(); i != spheres.end(); ++i ) {
//      TRps << *i << std::endl;
//   }
//
// }
//
// void test_sasa(std::string fname) {
//   using namespace core::scoring::packstat;
//   using namespace std;
// 	using namespace core;
//
// 	pose::Pose pose;
// 	core::import_pose::pose_from_pdb(pose,fname);
//
// 	utility::vector1< Real > radii;
// 	{
// 		chemical::AtomTypeSet const & atom_set( pose.residue(1).atom_type_set() );
// 		Size const SASA_RADIUS_INDEX( atom_set.extra_parameter_index( "SASA_RADIUS_LEGACY" ) );
// 		radii.resize( atom_set.n_atomtypes() );
// 		for ( Size i=1; i<= radii.size(); ++i ) {
// 			chemical::AtomType const & at( atom_set[i] );
// 			radii[i] = atom_set[i].extra_parameter( SASA_RADIUS_INDEX );
// 		}
// 	}
//
// 	Spheres pose_spheres;
// 	{
// 		id::AtomID_Map<Real> atom_sasa;
// 		utility::vector1<Real> rsd_sasa;
// 		scoring::calc_per_atom_sasa(pose,atom_sasa,rsd_sasa,3.0);
// 		for( size_t i = 1; i <= pose.total_residue(); ++i ) {
// 			for( size_t j = 1; j <= pose.residue(i).natoms(); ++j ) {
// 				TRps << "3.0 " << i << " " << j << " " << atom_sasa[i][j] << std::endl;
// 				conformation::Atom const & a( pose.residue(i).atom(j) );
// 				pose_spheres.push_back( Sphere( a.xyz(), radii[a.type()] ) );
// 			}
// 		}
// 		TRps << std::endl;
// 	}
// 	{
// 		id::AtomID_Map<Real> atom_sasa;
// 		utility::vector1<Real> rsd_sasa;
// 		scoring::calc_per_atom_sasa(pose,atom_sasa,rsd_sasa,1.4);
// 		for( size_t i = 1; i <= pose.total_residue(); ++i ) {
// 			for( size_t j = 1; j <= pose.residue(i).natoms(); ++j ) {
// 				TRps << "1.4 " << i << " " << j << " " << atom_sasa[i][j] << std::endl;
// 			}
// 		}
// 		TRps << std::endl;
// 	}
//
// 	AtomRadiusMap arm;
//   SimplePDB pdb;
//   utility::io::izstream in(fname.c_str());
//   in >> pdb;
//
// 	// TRps << fname << std::endl;
//   Spheres spheres( pdb.get_spheres(arm) );
//   // for( SphereIter i = spheres.begin(); i != spheres.end(); ++i ) {
//   //   TRps << *i << std::endl;
//   // }
// 	spheres = pose_spheres;
//
// 	SasaOptions opts;
// 	opts.probe_radii.push_back(3.0);
// 	opts.probe_radii.push_back(1.4);
// 	SasaResultOP sr = compute_sasa( spheres, opts );
// 	for( size_t p = 1; p <= 2; ++p ) {
// 		for( size_t i = 1; i <= spheres.size(); ++i ) {
// 			TRps << opts.probe_radii[p] << " " << i << " " << sr->sphere_sasa(i,p) << std::endl;
// 		}
// 	}
//
// }
//
// struct OrderSphereOnX {
//   bool operator()( core::scoring::packstat::Sphere const & a, core::scoring::packstat::Sphere const & b ) {
// 		return a.xyz.x() < b.xyz.x();
//   }
// };
// struct RealAscending {
//   bool operator()( core::Real const a, core::Real const b ) {
// 		return a < b;
//   }
// };
//
// core::Real median( utility::vector1<core::Real> l) {
// 	std::sort(l.begin(),l.end(),RealAscending());
// 	return( l[ std::max((size_t)1,l.size()/2) ] );
// }
//
// Spheres pose_to_spheres( core::pose::Pose & pose ) {
// 	using namespace std;
// 	ostringstream oss;
// 	core::io::pdb::dump_pdb(pose,oss);
// 	istringstream iss( oss.str() );
// 	AtomRadiusMap arm;
//   SimplePDB pdb;
//   iss >> pdb;
// 	return pdb.get_spheres(arm);
// }
//
// void write_spheres_pdb( Spheres & spheres, std::string fname ) {
// 	using namespace core::scoring::packstat;
// 	std::ofstream out(fname.c_str());
// 	for( SphereIter i = spheres.begin(); i != spheres.end(); ++i ) {
// 		int rnum = 0, anum = 0;
// 		PackstatReal occ = 0.0f;
// 	   out << "ATOM  " + I( 5, ( anum ) ) + "  V   PRT Z"
// 				+ I( 4, rnum ) + "    "
// 				+ F( 8, 3, i->xyz.x() ) + F( 8, 3, i->xyz.y() ) + F( 8, 3, i->xyz.z() )
// 				+ F( 6, 2, occ ) + ' ' + F( 5, 2, i->radius ) << std::endl;
// 	}
// 	out.close();
// }

inline std::string base_name(const std::string& str) {
  size_t begin = 0;
  size_t end = str.length();

  for (size_t i=0; i<str.length(); ++i) {
    if (str[i] == '/') begin = i+1;
  }

  return str.substr(begin,end);
}

std::string get_out_tag(std::string fname) {
	std::string base = base_name(fname);
	std::transform( base.begin(), base.end(), base.begin(), tolower );
	system( ("mkdir -p out/" + base.substr(1,2)).c_str() );
	std::string OUT_TAG = "out/" + base.substr(1,2) + "/" + base;
	return OUT_TAG;
}

void output_packstat_pdb( std::string fname, utility::vector1<core::Real> const & res_scores ) {
	using namespace core::scoring::packstat;
  using namespace std;
	using namespace core;
	using namespace basic::options;
	using namespace ObjexxFCL::format;
	using namespace numeric;
	using namespace utility;

	bool surface_accessibility = option[ OptionKeys::packstat::surface_accessibility ]();
	core::Real burial_radius   = option[ OptionKeys::packstat::cavity_burial_probe_radius ]();

	string OUT_TAG = get_out_tag(fname);

	AtomRadiusMap arm;
  SimplePDB pdb;
  utility::io::izstream in(fname.c_str());
  in >> pdb;

	Spheres spheres;
	TRps << fname << std::endl;
	spheres = pdb.get_spheres(arm);
	std::cerr << "spheres len: " << spheres.size() << std::endl;
	vector1< xyzVector<PackstatReal> > centers( pdb.get_res_centers() );
	std::cerr << "centers len: " << centers.size() << std::endl;

	SasaOptions opts;
	opts.prune_max_iters = 999;
	opts.prune_max_delta = 0;
	opts.num_cav_ball_layers = 10;
	opts.frac_cav_ball_required_exposed = 0.00;
	opts.area_cav_ball_required_exposed = 0.00;
	opts.surrounding_sasa_smoothing_window = 3;
	opts.min_cav_ball_radius = 0.9;
	for( PackstatReal pr = 3.0; pr >= 0.4; pr -= 0.1 ) opts.probe_radii.push_back(pr);
	opts.prune_cavity_burial_probe_radii.push_back( burial_radius );
	if( surface_accessibility ) {
		for( PackstatReal pr = burial_radius-0.1; pr >= 0.1; pr -= 0.1 ) {
			opts.prune_cavity_burial_probe_radii.push_back(pr);
		}
	}

	TRps << "compute MSAs" << std::endl;
	SasaResultOP sr = compute_sasa( spheres, opts );

	////////////////////////////////////////////////////////////////////////////////////////////////
	CavBalls cavballs = sr->cavballs;
	TRps << "pruning hidden cav balls " << cavballs.size() << std::endl;
	cavballs = prune_hidden_cavity_balls( cavballs, opts );

	TRps << "pruning exposed cav balls " << cavballs.size() << std::endl;
	cavballs = prune_cavity_balls( spheres, cavballs, opts );

	TRps << "compute cav ball volumes	" << cavballs.size() << std::endl;
	compute_cav_ball_volumes( cavballs, opts );

	///////////////////////////////////////////////////////////////////////////////////////////////
	TRps << "writting stupid pdb to "+OUT_TAG+".pdb" << std::endl;
	ostringstream out;
	// if( pdb.atoms().size() != spheres.size() ) { // hydrogens!
	// 	TRps << "WARNING: output_packstat_pdb: size of pdb.atoms() != size of spheres" << std::endl;
	// 	TRps << "         spheres: " << spheres.size() << ", pdb.atoms(): " << pdb.atoms().size() << std::endl;
	// }
	int prev_resnum = -12345;
	int resnum = 0;
	int res_atom_count = 999;
	for( int i = 1; i <= (int)pdb.atoms().size(); ++i ) {
		SimplePDB_Atom & atom( pdb.atoms()[i] );
		// std::cerr << "FOO " << resnum << " " << res_atom_count << " " << prev_resnum << " " << atom.resnum << std::endl;
		if( prev_resnum != atom.resnum ) {
			if( res_atom_count > 3 ) {
				resnum++;
			}
			prev_resnum = atom.resnum;
			res_atom_count = 0;
		}
		res_atom_count++;
		core::Real bfac = 1.0;
		if( resnum > (int)res_scores.size() ) {
			TRps << "ERROR: more residues than res scores " << resnum << " " << res_scores.size() << std::endl;
		} else {
			bfac = res_scores[resnum];
		}
		// std::cerr << "bfac " << resnum << " " << res_scores.size() << " " << res_scores[resnum] << std::endl;
		out << atom.whole_line.substr(0,61) << F(4,2,bfac) << atom.whole_line.substr(65) << std::endl;
	  // out << LJ(6,atom.ATOM) + I( 5, ( atom.num ) ) + " "
	  // 				+ LJ(4,atom.type) + " " + LJ(3,atom.res) + " " + atom.chain
	  // 				+ I( 4, atom.resnum ) + "    "
	  // 				+ F( 8, 3, atom.x ) + F( 8, 3, atom.y ) + F( 8, 3, atom.z )
	  // 				+ F( 6, 2, atom.occ ) + ' ' + F( 5, 2, atom.bfac ) << std::endl;
	}
	for( size_t cb = 1; cb <= cavballs.size(); ++cb ) {
		if( cavballs[cb].radius() > 0.9 )
			out << cavballs[cb].hetero_atom_line() << std::endl;
	}
	// for( size_t cb = 1; cb <= sel_cbs.size(); ++cb) {
	// 	out << sel_cbs[cb].hetero_atom_line(7) << std::endl;
	// }
	utility::io::ozstream outz((OUT_TAG+".pdb").c_str());
	outz << out.str();
	outz.close();
	out.clear();

}

void output_packstat( std::string fname ) {

	using namespace core::scoring::packstat;
  using namespace std;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace ObjexxFCL::format;
	using namespace numeric;
	using namespace utility;

	int oversample      = option[ OptionKeys::packstat::oversample ]();
	bool include_water  = option[ OptionKeys::packstat::include_water ]();
	bool residue_scores = option[ OptionKeys::packstat::residue_scores ]();
	bool packstat_pdb   = option[ OptionKeys::packstat::packstat_pdb ]();

	AtomRadiusMap arm;
  SimplePDB pdb;
  utility::io::izstream in(fname.c_str());
  in >> pdb;

	Spheres spheres;
	TRps << fname << std::endl;
	spheres = pdb.get_spheres(arm);
	// write_spheres_pdb(spheres,"spheres_from_pdb.pdb");
	// std::sort( spheres.begin(), spheres.end(), OrderSphereOnX() );
	// std::cerr << "spheres len: " << spheres.size() << std::endl;
	vector1< xyzVector<PackstatReal> > centers( pdb.get_res_centers() );
	// std::cerr << "centers len: " << centers.size() << std::endl;

	PosePackData pd;
	pd.spheres = spheres;
	pd.centers = centers;
	core::Real packing_score = compute_packing_score( pd, oversample );

	TRps << "packing score: " << fname << " " << packing_score;
	if( include_water ) {
		TRps << " ( with " << pdb.num_water() << " waters )";
	}
	TRps << std::endl;

	utility::vector1<core::Real> res_scores; // needed if output res scores or pdb
	if( packstat_pdb || residue_scores ) {
		res_scores = compute_residue_packing_scores( pd, oversample );
	}

	if( packstat_pdb ) {
		output_packstat_pdb( fname, res_scores );
	}

	if( residue_scores ) {
		for( int i = 1; i <= (int)res_scores.size(); ++i ) {
			TRps << "packing score: residue " << i << " " << res_scores[i] << std::endl;
		}
	}

}


int main (int argc, char *argv[])
{

	try {


	devel::init( argc, argv );

  using namespace core::scoring::packstat;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace utility;

  const Real PI = numeric::NumericTraits<Real>::pi();

  // test_io();

	// test_sasa_dots();

	if( option[ in::file::s ].user() ) {
  	vector1<file::FileName> files( option[ in::file::s ]() );
  	for( size_t i = 1; i <= files.size(); ++i ) {
			output_packstat( files[i] );
  	}
	} else if( option[ in::file::l ].user() ) {
  	vector1<file::FileName> files( option[ in::file::l ]() );
  	for( size_t i = 1; i <= files.size(); ++i ) {
			utility::io::izstream list( files[i] );
			std::string fname;
			while( list >> fname ) {
				// std::cerr << "'" << fname << "'" << std::endl;
				output_packstat( fname );
			}
  	}
	}
	return 0;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


// void print_packing_score( std::string fname ) {
// 	using namespace core::scoring::packstat;
//   using namespace std;
// 	using namespace core;
// 	using namespace ObjexxFCL::format;
// 	using namespace numeric;
// 	using namespace utility;
//
// 	AtomRadiusMap arm;
//   SimplePDB pdb;
//   utility::io::izstream in(fname.c_str());
//   in >> pdb;
//
// 	Spheres spheres;
// 	TRps << fname << std::endl;
// 	spheres = pdb.get_spheres(arm);
// 	// write_spheres_pdb(spheres,"spheres_from_pdb.pdb");
// 	std::sort( spheres.begin(), spheres.end(), OrderSphereOnX() );
// 	// std::cerr << "spheres len: " << spheres.size() << std::endl;
// 	vector1< xyzVector<PackstatReal> > centers( pdb.get_res_centers() );
// 	// std::cerr << "centers len: " << centers.size() << std::endl;
//
// 	PosePackData pd;
// 	pd.spheres = spheres;
// 	pd.centers = centers;
// 	core::Real packing_score = compute_packing_score( pd, 0 );
//
// 	std::cerr << "packing score: " << fname << " " << packing_score;
// 	if( basic::options::option[ basic::options::OptionKeys::packstat::include_water ]() ) {
// 		std::cerr << " (" << pdb.num_water() << " waters)";
// 	}
// 	std::cerr << std::endl;
//
// }
//
// void test_cav_balls( std::string fname ) {
// 	using namespace core::scoring::packstat;
//   using namespace std;
// 	using namespace core;
// 	using namespace ObjexxFCL::format;
// 	using namespace numeric;
// 	using namespace utility;
//
// 	time_t before_comp_ps = clock();
// 	if(true) {
//
// 		AtomRadiusMap arm;
// 	  SimplePDB pdb;
// 	  utility::io::izstream in(fname.c_str());
// 	  in >> pdb;
// 		Spheres spheres;
// 		TRps << fname << std::endl;
// 		spheres = pdb.get_spheres(arm);
// 		std::sort( spheres.begin(), spheres.end(), OrderSphereOnX() );
// 		vector1< xyzVector<PackstatReal> > centers( pdb.get_res_centers() );
// 		PosePackData pd;
// 		pd.centers = centers;
// 		pd.spheres = spheres;
// 		core::Real discrim = 0, respred = 0, discrim2 = 0, respred2 = 0;
// 		std::cerr << "PD_PACKSCORE " << fname << " "
// 						  << compute_packing_score(pd,1)
// 							<< std::endl;
//   	return;
// 		// for( int i = 1; i <= 10; ++i ) {
// 		// 	pair<Real,Real> ps = compute_packing_scores( pd, 0 );
// 		// 	discrim += ps.first;
// 		// 	respred += ps.second;
// 		// 	discrim2 += ps.first * ps.first;
// 		// 	respred2 += ps.second * ps.second;
// 		// 	std::cout << fname << " " << spheres.size() << " " << centers.size() << " " << ps.first << " " << ps.second << std::endl;
// 		// }
// 		// discrim /= 10.0; respred /= 10.0; discrim2 /= 10.0; respred2 /= 10.0;
// 		// std::cout << fname << " mean " << discrim << " " << respred << " sd " << sqrt(discrim2-discrim*discrim) << " " << sqrt(respred2-respred*respred) << " " << std::endl;
// 		// TRps << "comp_ps time: " << ((float)(clock()-before_comp_ps))/1000000.0 << std::endl;
//
//
// 		// Pose pose;
// 		// core::import_pose::pose_from_pdb( pose, fname );
// 		//
// 		// std::cerr << "PACKSCORE " << fname << " " << compute_packing_score(pose) << std::endl;
// 		// return;
// 		//
// 		// // core::Real discrim, respred;
// 		// // compute_packing_scores( pose, discrim, respred );
// 		// // std::cerr << "SCORES: " << discrim << " " << respred << std::endl;
// 		// // std::cerr << "comp_ps time: " << ((float)(clock()-before_comp_ps))/1000000.0 << std::endl;
// 		// core::id::AtomID_Map<Real> atom_scores = compute_atom_packing_scores(pose);
// 		// ofstream out((fname+".byatom.pdb").c_str());
// 		// core::io::pdb::dump_bfactor_pdb(pose,atom_scores,out);
// 		// out.close();
// 		//
// 		// return;
// 	}
//
// 	time_t start_t = clock();
//
// 	string OUT_TAG = get_out_tag(fname);
// 	TRps << "OUT_TAG " << OUT_TAG << std::endl;
//
// 	AtomRadiusMap arm;
//   SimplePDB pdb;
//   utility::io::izstream in(fname.c_str());
//   in >> pdb;
//
// 	Spheres spheres;
// 	TRps << fname << std::endl;
// 	spheres = pdb.get_spheres(arm);
// 	write_spheres_pdb(spheres,"spheres_from_pdb.pdb");
// 	std::sort( spheres.begin(), spheres.end(), OrderSphereOnX() );
// 	std::cerr << "spheres len: " << spheres.size() << std::endl;
// 	vector1< xyzVector<PackstatReal> > centers( pdb.get_res_centers() );
// 	std::cerr << "centers len: " << centers.size() << std::endl;
//
// 	SasaOptions opts;
// 	opts.prune_max_iters = 999;
// 	opts.prune_max_delta = 0;
// 	opts.num_cav_ball_layers = 10;
// 	opts.frac_cav_ball_required_exposed = 0.00;
// 	opts.area_cav_ball_required_exposed = 0.00;
// 	opts.surrounding_sasa_smoothing_window = 3;
// 	opts.min_cav_ball_radius = 0.9;
// 	for( PackstatReal pr = 3.0; pr >= 0.4; pr -= 0.0333 ) opts.probe_radii.push_back(pr);
// 	// opts.prune_cavity_burial_probe_radii.push_back( 1.6 );
// 	for( PackstatReal pr = 1.4; pr >= 0.1; pr -= 0.1    ) opts.prune_cavity_burial_probe_radii.push_back(pr);
//
// 	time_t io_t = clock();
//
// 	TRps << "compute MSAs" << std::endl;
// 	SasaResultOP sr = compute_sasa( spheres, opts );
//
// 	time_t compute_sasa_t = clock();
//
// 	////////////////////////////////////////////////////////////////////////////////////////////////
// 	CavBalls cavballs = sr->cavballs;
// 	TRps << "pruning hidden cav balls " << cavballs.size() << std::endl;
// 	cavballs = prune_hidden_cavity_balls( cavballs, opts );
// 	time_t prune_1_t = clock();
//
// 	TRps << "pruning exposed cav balls " << cavballs.size() << std::endl;
// 	cavballs = prune_cavity_balls( spheres, cavballs, opts );
// 	time_t prune_2_t = clock();
//
// 	TRps << "compute cav ball volumes	" << cavballs.size() << std::endl;
// 	compute_cav_ball_volumes( cavballs, opts );
// 	time_t compute_vols_t = clock();
//
// 	TRps << "selecting representatives" << std::endl;
// 	CavBalls sel_cbs( select_cav_balls(cavballs,4.0) );
// 	time_t sel_cb_t = clock();
//
// 	std::cerr << "times: "  << std::endl;
// 	std::cerr << "comp_sasa " << ((float)(compute_sasa_t - io_t          ))/1000000.0 << std::endl;
// 	std::cerr << "prune_1   " << ((float)(prune_1_t      - compute_sasa_t))/1000000.0 << std::endl;
// 	std::cerr << "prune_2   " << ((float)(prune_2_t      - prune_1_t     ))/1000000.0 << std::endl;
// 	std::cerr << "sa & vol  " << ((float)(compute_vols_t - prune_2_t     ))/1000000.0 << std::endl;
//
// 	ostringstream cav_info;
// 	for( CavBallIter i = cavballs.begin(); i != cavballs.end(); ++i ) {
// 		// PackstatReal maxa = i->radius() * i->radius() * 4.0 * PI * 1.00;
// 		// PackstatReal maxv = i->radius() * i->radius() * i->radius() * 4.0/3.0 * PI * 1.00;
// 		cav_info << F( 7, 3, i->radius() ) << " "
// 		         << F( 7, 3, i->exposed_radius ) << " "
// 		         << F( 7, 3, i->area     ) << " "
// 		         << F( 7, 3, i->vol      ) << " "
// 		         << std::endl;
// 		// if( i->area > maxa * 1.01 ) TRps << "area too high " << i->area << " " << maxa << std::endl;
// 		// if( i->vol  > maxv * 1.01 ) TRps << "vol. too high " << i->vol  << " " << maxv << std::endl;
// 	}
// 	utility::io::ozstream cav_info_file(((OUT_TAG+".cavities").c_str()));
// 	cav_info_file << cav_info.str();
// 	cav_info_file.close();
// 	cav_info.clear();
//
//
// 	////////////////////////////////////////////////////////////////////////////////////////////////
// 	// TRps << "SASA" << std::endl;
// 	// for( size_t i = 1; i <= spheres.size(); ++i ) {
// 	// 	TRps << I(5,i) << " ";
// 	// 	for( size_t p = opts.probe_radii.size(); p >= 1; --p ) {
// 	// 		TRps << F( 5, 2, sr->sphere_sasa(i,p) ) << " ";
// 	// 	}
// 	// 	TRps << std::endl;
// 	// }
//
// 	///////////////////////////////////////////////////////////////////////////////////////////////
// 	TRps << "writting stupid pdb to "+OUT_TAG+".pdb" << std::endl;
// 	ostringstream out;
// 	for( SphereIter i = spheres.begin(); i != spheres.end(); ++i ) {
// 		int rnum = 0, anum = 0;
// 		PackstatReal occ = 0.0f;
// 	  out << "ATOM  " + I( 5, ( anum ) ) + "  V   PRT Z"
// 				+ I( 4, rnum ) + "    "
// 				+ F( 8, 3, i->xyz.x() ) + F( 8, 3, i->xyz.y() ) + F( 8, 3, i->xyz.z() )
// 				+ F( 6, 2, occ ) + ' ' + F( 5, 2, i->radius ) << std::endl;
// 	}
// 	for( size_t cb = 1; cb <= cavballs.size(); ++cb ) {
// 		if( cavballs[cb].radius() > 0.9 )
// 			out << cavballs[cb].hetero_atom_line() << std::endl;
// 	}
// 	// for( size_t cb = 1; cb <= sel_cbs.size(); ++cb) {
// 	// 	out << sel_cbs[cb].hetero_atom_line(7) << std::endl;
// 	// }
// 	utility::io::ozstream outz((OUT_TAG+".pdb").c_str());
// 	outz << out.str();
// 	outz.close();
// 	out.clear();
//
// 	//////////////////////////////////////////////////////////////////////////////////////////////////
// 	TRps << "output_res_surrounding_sasa " << centers.size() << std::endl;
// 	ostringstream outss2;
// 	for( size_t i = 1; i <= centers.size(); ++i ) {
// 		outss2 << i << " ";
// 		PackingScoreResDataOP d= compute_surrounding_sasa( centers[i], spheres, sr, opts );
// 		for( Size i = 1; i <= d->nrad(); ++i ) {
// 			for( Size j = 1; j <= d->npr(); ++j ) {
// 				outss2 << d->msa(i,j) << " ";
// 			}
// 		}
// 		outss2 << std::endl;
// 	}
// 	utility::io::ozstream outss2z((OUT_TAG+".res_sur_sasa").c_str());
// 	outss2z << outss2.str();
// 	outss2z.close();
// 	outss2.clear();
//
//
// 	// TRps << "output_surrounding_sasa " << sel_cbs.size() << std::endl;
// 	// ostringstream outss;
// 	// for( size_t i = 1; i <= sel_cbs.size(); ++i ) {
// 	// 	CavityBall & cb( sel_cbs[i] );
// 	// 	outss << cb.radius_ << " "
// 	// 				<< cb.exposed_radius << " "
// 	// 				<< cb.area << " "
// 	// 				<< cb.vol << " ";
// 	// 	output_surrounding_sasa( outss, sel_cbs[i].xyz(), spheres, sr, opts );
// 	// 	outss << std::endl;
// 	// }
// 	// utility::io::ozstream outssz((OUT_TAG+".cav_sur_sasa").c_str());
// 	// outssz << outss.str();
// 	// outssz.close();
// 	// outss.clear();
// }
//
// numeric::xyzMatrix<PackstatReal> rand_rot() {
// 	using namespace numeric;
// 	using namespace numeric::random;
// 	xyzVector<PackstatReal> axis(uniform(),uniform(),uniform());
// 	while( axis.length() > 1 ) axis = xyzVector<PackstatReal>(uniform(),uniform(),uniform());
// 	return rotation_matrix<PackstatReal>( axis, uniform() * 2 * PI );
// }
//
// void test_sasa_dots() {
// 	using namespace utility;
// 	using namespace numeric;
// 	using namespace core::scoring::packstat;
// 	using namespace std;
//
// 	vector1<xyzVector<PackstatReal> > sasa_dots = get_sasa_dot_locations();
//
// 	{
// 		ofstream out("sasa_dots_randomized.pdb");
// 		for( PackstatReal pr = 350.0; pr <= 600.0; pr += 10.0 ) {
// 			xyzMatrix<PackstatReal> rot = rand_rot();
// 			for( size_t i = 1; i <= sasa_dots.size(); ++i ) {
// 				xyzVector<PackstatReal> dot = rot * sasa_dots[i];
// 				int rnum = 0, anum = 0;
// 				PackstatReal occ = 0.0f;
// 			    out << "ATOM  " + I( 5, ( anum ) ) + "  V   PRT Z"
// 				+ I( 4, rnum ) + "    "
// 				+ F( 8, 3, dot.x()*pr ) + F( 8, 3, dot.y()*pr ) + F( 8, 3, dot.z()*pr )
// 				+ F( 6, 2, occ ) + ' ' + F( 5, 2, 5.0 ) << std::endl;
// 			}
// 		}
// 		out.close();
// 	}
// 	{
// 		ofstream out("sasa_dots.pdb");
// 		for( PackstatReal pr = 350.0; pr <= 600.0; pr += 10.0 ) {
// 			for( size_t i = 1; i <= sasa_dots.size(); ++i ) {
// 				xyzVector<PackstatReal> dot = sasa_dots[i];
// 				int rnum = 0, anum = 0;
// 				PackstatReal occ = 0.0f;
// 			    out << "ATOM  " + I( 5, ( anum ) ) + "  V   PRT Z"
// 				+ I( 4, rnum ) + "    "
// 				+ F( 8, 3, dot.x()*pr ) + F( 8, 3, dot.y()*pr ) + F( 8, 3, dot.z()*pr )
// 				+ F( 6, 2, occ ) + ' ' + F( 5, 2, 5.0 ) << std::endl;
// 			}
// 		}
// 		out.close();
// 	}
//
// }

