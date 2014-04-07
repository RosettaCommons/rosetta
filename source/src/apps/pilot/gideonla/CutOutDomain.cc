// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file   CutOutDomain.cc
//
/// @brief Created to extract domain
/// @author Gideon Lapidoth (glapidoth@gmail.com)
/// @date 08 May, 2013


// Mini-Rosetta headers
#include <time.h>
#include <numeric/constants.hh>
#include <numeric/random/random.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/util.hh>//option.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/rosetta_scripts/util.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/pose/util.hh>
#include <core/import_pose/pose_stream/util.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cutoutdomain.OptionKeys.gen.hh>




using basic::T;
using basic::Error;
using basic::Warning;
clock_t clk1; //for timing the process
clock_t clk2;



core::Size
find_nearest_res( core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size const chain/*=0*/ ){
  core::Real min_dist( 100000 ); core::Size nearest_res( 0 );
  core::Size i;
  for( i = 1; i < target.total_residue(); ++i ){
		if( target.residue( i ).is_ligand() ) continue;
		if( chain && target.residue( i ).chain() != chain ) continue;
    core::Real const dist( source.residue( res ).xyz( "CA" ).distance( target.residue( i ).xyz( "CA" ) ) );
    if( dist <= min_dist ){
      min_dist = dist;
      nearest_res = i;
    }
  }
  static basic::Tracer TR("This is nearest_res");

      TR<<nearest_res<<std::endl;
  if( min_dist <= 10.0 ) return nearest_res;
  else return 0;
}
///////////////////////////////////////////////////////////////////////////////
using namespace basic::options;

int
main( int argc, char * argv [] )
{
    try {
	  clk1 = clock();

	devel::init(argc, argv);
	using namespace protocols::jd2;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::cutoutdomain;
	static basic::Tracer TR("CutOutDomain:");

   option.add_relevant( OptionKeys::cutoutdomain::start                       );

   if ( !option[ OptionKeys::in::file::s ].user() ||
    !option[ OptionKeys::cutoutdomain::start ].user() ||
    !option[ OptionKeys::cutoutdomain::end ].user()
  ) {
    TR  << "Usage: " <<std:: endl <<
      "-s <Template file> <aligned PDB files> " <<
      "  -cutoutdomain:start <residue number> -cutoutdomain:end <residue number> "
      "  -database <minirosetta_db>"  << std::endl;
    exit(-1);
  }


	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
	core::Size count = 0;
	core::pose::Pose pose; //Contains the first pdb (template pdb) in the command line
	input.fill_pose( pose);
	//std::string name(pose.pdb_info()->name());
	//TR<<name<<std::endl;
	core::Size start(
		option[ basic::options::OptionKeys::cutoutdomain::start ]()
	);

	core::Size end(
		option[ basic::options::OptionKeys::cutoutdomain::end ]()
	);

	while( input.has_another_pose() && (count < 1000 ) ) {
		try{
		core::pose::Pose Temp_pose;
		core::pose::Pose new_pose;
		input.fill_pose( Temp_pose);
		std::string name(Temp_pose.pdb_info()->name());
		TR<<name<<", total number of residues: "<<Temp_pose.total_residue()<<std::endl;
		/*if (pdb.empty()){
			continue;
		}*/
		core::Size from = find_nearest_res(pose,Temp_pose,start, 1/*chain*/ );
	//	TR<<from<<std::endl;
		core::Size to  = find_nearest_res(pose,Temp_pose,end, 1/*chain*/ );
		if (to==0||from==0){
			continue;
		}
		//TR<<to<<std::endl;
		TR<<"First resdiue in the target pose "<<name<<" : "<<Temp_pose.residue(from).name1()<<from<<std::endl;
		TR<<"End resdiue in the target pose "<<name<<" : "<<Temp_pose.residue(to).name1()<<to<<std::endl;
		Temp_pose.conformation().delete_residue_range_slow( to+1,Temp_pose.total_residue() );
		Temp_pose.conformation().delete_residue_range_slow( 1,from-1 );

		TR<<"First resdiue in the target pose "<<name<<" : "<<to<<std::endl;
		std::string output_name = Temp_pose.pdb_info()->name()+"_loop.pdb";
		Temp_pose.dump_pdb(output_name);


		count ++;
	}
	catch (int k) {
		continue;
	}
	}

	clk2 = clock();
	TR<<"The entire process took: "<<(clk2-clk1)/CLOCKS_PER_SEC<<std::endl;
    } catch ( utility::excn::EXCN_Base const & e ) {
			std::cerr << "caught exception " << e.msg() << std::endl;
			return -1;
    }
	return 0;
}
