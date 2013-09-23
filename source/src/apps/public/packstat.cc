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

#include "core/scoring/packstat/types.hh"
#include "core/scoring/packstat/SimplePDB_Atom.hh"
#include "core/scoring/packstat/SimplePDB.hh"
// AUTO-REMOVED #include "core/scoring/packstat/io.hh"
#include "core/scoring/packstat/AtomRadiusMap.hh"
#include "core/scoring/packstat/compute_sasa.hh"
// AUTO-REMOVED #include "core/scoring/packstat/sasa_dot_locations.hh"
// AUTO-REMOVED #include "core/scoring/packstat/packing_score_params.hh"

#include <protocols/analysis/PackStatMover.hh>

#include <devel/init.hh>
#include "core/types.hh"
// AUTO-REMOVED #include "core/id/AtomID_Map.hh"
#include <basic/options/option.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>

// AUTO-REMOVED #include "core/pose/Pose.hh"
// AUTO-REMOVED #include "core/io/pdb/pose_io.hh"
// AUTO-REMOVED #include "core/scoring/sasa.hh"

#include "basic/Tracer.hh"

#include "utility/vector1.hh"
#include "utility/file/FileName.hh"
#include "utility/io/izstream.hh"
#include "utility/io/ozstream.hh"

// AUTO-REMOVED #include "numeric/xyz.functions.hh"
// AUTO-REMOVED #include "numeric/random/random.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
// AUTO-REMOVED #include <ctype.h>

// AUTO-REMOVED #include <time.h>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packstat.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>

#include <protocols/jd2/JobDistributor.hh>

#include <utility/excn/Exceptions.hh>

using core::Real;

basic::Tracer TRps("packstat");

using namespace core::scoring::packstat;

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
	if (system( ("mkdir -p out/" + base.substr(1,2)).c_str() ) == -1) {
		TRps.Error << "Unable to make directory!" << std::endl;
	}
	std::string OUT_TAG = "out/" + base.substr(1,2) + "/" + base;
	return OUT_TAG;
}

void output_packstat_pdb( std::string fname, utility::vector1<core::Real> const & res_scores ) {
	using core::Size;
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
	opts.min_cav_ball_radius = option[ OptionKeys::packstat::min_cav_ball_radius ]();
	opts.min_cluster_overlap = option[ OptionKeys::packstat::min_cluster_overlap ]();
	opts.cluster_min_volume = option[ OptionKeys::packstat::cluster_min_volume ]();
	for( PackstatReal pr = 3.0; pr >= 0.4; pr -= 0.1 ) opts.probe_radii.push_back(pr);
	opts.prune_cavity_burial_probe_radii.push_back( burial_radius );
	if( surface_accessibility ) {
		for( PackstatReal pr = burial_radius-0.1; pr >= 0.1; pr -= 0.1 ) {
			opts.prune_cavity_burial_probe_radii.push_back(pr);
		}
	}

	//TRps << "compute MSAs" << std::endl;
	SasaResultOP sr = compute_sasa( spheres, opts );

	////////////////////////////////////////////////////////////////////////////////////////////////
	CavBalls cavballs = sr->cavballs;
	//TRps << "pruning hidden cav balls " << cavballs.size() << std::endl;
	cavballs = prune_hidden_cavity_balls( cavballs, opts );

	//TRps << "pruning exposed cav balls " << cavballs.size() << std::endl;
	cavballs = prune_cavity_balls( spheres, cavballs, opts );

	//TRps << "compute cav ball volumes	" << cavballs.size() << std::endl;
	compute_cav_ball_volumes( cavballs, opts );

	vector1< CavityBallCluster > clusters =
	compute_cav_ball_clusters( cavballs, opts );

	///////////////////////////////////////////////////////////////////////////////////////////////
	//TRps << "writting stupid pdb to "+OUT_TAG+".pdb" << std::endl;
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
	// for( size_t cb = 1; cb <= cavballs.size(); ++cb ) {
	// 	// if( cavballs[cb].radius() > 0.9 )
	// 		out << cavballs[cb].hetero_atom_line() << std::endl;
	// }
	// for( size_t cb = 1; cb <= sel_cbs.size(); ++cb) {
	// 	out << sel_cbs[cb].hetero_atom_line(7) << std::endl;
	// }

	 std::cout << "Add Cavities to PDB:" << std::endl;
	Size count = 1;
	for( Size i = 1; i <= clusters.size(); i++ ) {
		// if( clusters[i].volume < 150.0 ) continue;
		if( clusters[i].surface_accessibility < option[ OptionKeys::packstat::min_surface_accessibility ]() ) continue;
		std::cout << "Cavity: " << count
			<< " volume "       << clusters[i].volume
			<< " surf "         << clusters[i].surface_area
			<< " surf. acc. "   << clusters[i].surface_accessibility
  			 << std::endl;
		for( Size j = 1; j <= clusters[i].cavballs.size(); j++ ) {
			// if( clusters[i].cavballs[j].radius() > 0.7 )
				out << clusters[i].cavballs[j].hetero_atom_line( centers.size()+count, count, 0.0 ) << std::endl;
		}
		count++;
	}

	// for( Size i = 1; i <= clusters.size(); i++ ) {
	// 	numeric::xyzVector<core::Real> xyz_ = clusters[i].center;
	// 	out << "HETATM" + I( 5, 1 ) + "  V   CTR Z"
	// 		+ I( 4, 1 ) + "    "
	// 		+ F( 8, 3, xyz_.x() ) + F( 8, 3, xyz_.y() ) + F( 8, 3, xyz_.z() )
	// 		+ F( 6, 2, 0.0 ) + ' ' + F( 5, 2, 0.0 ) << std::endl;
	// }

	utility::io::ozstream outz((OUT_TAG+".pdb").c_str());
	outz << out.str();
	outz.close();
	out.clear();

}

void output_packstat( std::string fname ) {

	using core::Size;
	using namespace core::scoring::packstat;
  	using namespace std;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace ObjexxFCL::format;
	using namespace numeric;
	using namespace utility;

	core::Size oversample = option[ OptionKeys::packstat::oversample ]();
	bool include_water  = option[ OptionKeys::packstat::include_water ]();
	bool residue_scores = option[ OptionKeys::packstat::residue_scores ]();
	bool packstat_pdb   = option[ OptionKeys::packstat::packstat_pdb ]();
	bool raw_stats      = option[ OptionKeys::packstat::raw_stats ]();

	AtomRadiusMap arm;
	SimplePDB pdb;
	utility::io::izstream in(fname.c_str());
	in >> pdb;

	Spheres spheres;
	//TRps << fname << std::endl;
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

	TRps << "packing score: " << fname << " " << centers.size() << " " << spheres.size() << " " << packing_score;
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

	if( raw_stats ) {	// stupid duplicate code....

		assert( pd.spheres.size() > 0 );
		assert( pd.centers.size() > 0 );

		SasaOptions opts;
		opts.surrounding_sasa_smoothing_window = 1+2*oversample;
		opts.num_surrounding_sasa_bins = 7;
		for( core::Size ipr = 1; ipr <= 31; ++ipr ) {
			PackstatReal pr = 3.0 - ((double)(ipr-1))/10.0;
			PackstatReal ostep = 0.1 / (oversample*2.0+1.0);
			for( core::Size i = 1; i <= oversample; ++i )	opts.probe_radii.push_back( pr + i*ostep );
			opts.probe_radii.push_back( pr );
			for( core::Size i = 1; i <= oversample; ++i )	opts.probe_radii.push_back( pr - i*ostep );
		}

		SasaResultOP result = compute_sasa( pd.spheres, opts );
		for( core::Size i = 1; i <= pd.centers.size(); ++i ) {
			// std::cout << i << " ";
			PackingScoreResDataCOP dat( compute_surrounding_sasa( pd.centers[i], pd.spheres, result, opts ) );
			TRps << "RAW_STATS " << i << " ";
			for( Size i =1; i <= dat->nrad(); ++i ) {
				for( Size j =1; j <= dat->npr(); ++j ) {
					TRps << dat->msa(i,j) << " ";
				}
			}
			TRps << std::endl;

		} // end raw_stats

	}

}

int main (int argc, char *argv[]) {
	try {

	devel::init( argc, argv );

  using namespace core::scoring::packstat;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace utility;

  // test_io();

	// test_sasa_dots();
	if( !option[ out::output ].user() ){
		option[ out::nooutput ].value( true );
	}

	if( option[ in::file::silent ].user() ) {

		protocols::analysis::PackStatMoverOP pack_mover;
		//protocols::jobdist::universal_main(m);
		protocols::jd2::JobDistributor::get_instance()->go(pack_mover);

	} else {

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

	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;

}


