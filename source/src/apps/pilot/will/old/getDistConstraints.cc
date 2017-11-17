// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#include <core/id/AtomID_Map.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <core/scoring/packing/PoseBalls.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <basic/prof.hh>

#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/FArray2D.hh>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <fstream>
#include <map>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using std::string;
using namespace core;
using utility::vector1;
using core::id::AtomID;
using numeric::xyzVector;
using namespace core::scoring::packing;

// HolesParams params("/Users/sheffler/project/holes/holes_params.txt");


void test( std::string fname ) {

	// Real m,a,d = 0.00001;
	// while( d < 1.0 ) {
	// 	test_deriv(d,m,a);
	// 	std::cout << d << " " << m << " " << a << std::endl;
	// 	d *= 10;
	// }
	// return;

	using namespace std;
	using namespace core;
	using namespace io::pdb;
	using namespace pose;
	using namespace scoring;
	using namespace packing;
	using namespace core::scoring::packstat;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Pose pose;
	core::import_pose::pose_from_file(pose,fname, core::import_pose::PDB_file);

	utility::vector1<Size> cst_rsd( option[ in::target_residues ]() );
	if( cst_rsd.size() == 0 ) {
		utility_exit_with_message( "You need to specify -target_residue stupid!" );
	}

	PoseOP native;
	if( option[ in::file::native ].user() ) {
		native = new Pose;
		core::import_pose::pose_from_file(*native,option[ in::file::native ](), core::import_pose::PDB_file);
		core::scoring::calpha_superimpose_pose( *native, pose );
	}

	utility::vector1<std::map<id::AtomID,Real> > cst_dat = cavity_distance_constraint( pose, cst_rsd, native );
	for( Size j = 1; j <= cst_dat.size(); ++j ) {
		std::ofstream out((fname+".res_"+string_of(cst_rsd[j])+"_cav_dist.cst").c_str());
		std::map<id::AtomID,Real>::iterator i;
		for( i = cst_dat[j].begin(); i != cst_dat[j].end(); ++i ) {
			id::AtomID id = i->first;
			Real dist     = i->second;
			string aname1 = pose.residue(id.rsd()).atom_name(id.atomno());
			out << "AtomPair CB " << cst_rsd[j] << " " << aname1 << " " << id.rsd() << " GAUSSIANFUNC " << dist << " 2.0" << std::endl;
		}
		out.close();
	}

}

int
main (int argc, char *argv[])
{

	try {


	devel::init( argc, argv );

  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace utility;

	// test_gradient();
	// return 0;

	if( option[ in::file::s ].user() ) {
  	vector1<file::FileName> files( option[ in::file::s ]() );
  	for( size_t i = 1; i <= files.size(); ++i ) {
    	test( files[i] );
  	}
	} else if( option[ in::file::l ].user() ) {
  		vector1<file::FileName> files( option[ in::file::l ]() );
  		for( size_t i = 1; i <= files.size(); ++i ) {
				utility::io::izstream list( files[i] );
				std::string fname;
				while( list >> fname ) {
					// std::cerr << "'" << fname << "'" << std::endl;
					test( fname );
				}
  		}
	}

	return 0;


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
