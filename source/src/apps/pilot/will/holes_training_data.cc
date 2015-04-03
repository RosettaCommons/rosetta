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

#include <core/types.hh>
#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/packing/PoseBalls.hh>

#include <basic/options/option.hh>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <pstream.h>

#include <time.h>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


std::map<std::string,utility::io::ozstream*> outs;

void test( std::string fname ) {
	using namespace std;
	using namespace core;
	using namespace io::pdb;
	using namespace pose;
	using namespace scoring;
	using namespace packing;

	Pose pose;
	core::import_pose::pose_from_pdb(pose,fname);

	PoseBalls pb( pose );

	// neighbor counts
	std::cout << fname << " neighbor counts" << std::endl;
	utility::vector1<Size> anb5(pb.nballs()),anb10(pb.nballs()),anb15(pb.nballs()),anb20(pb.nballs());
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		for( Size j = 1; j < i; j++ ) {
			Real dis2 = numeric::distance_squared(pb.ball(i).xyz(),pb.ball(j).xyz());
			if(  25 >= dis2 ) anb5 [i]++;
			if( 100 >= dis2 ) anb10[i]++;
			if( 225 >= dis2 ) anb15[i]++;
			if( 400 >= dis2 ) anb20[i]++;
		}
	}
	// int infp,outfp,status;
	// char buf[9999999];

	time_t t = clock();
	std::cerr << "starting surf calculation: " << std::endl;

	string cmd("/Users/sheffler/project/koehl/DAlphaBall_new/DAlphaBall_surf.mactel");
	redi::pstream proc(cmd);
	proc << "NPOINTS" << endl << pb.nballs() << endl << "COORDS" << endl;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		Ball b(pb.ball(i));
	// 	// cout << "PoseBalls " << i << " " << pb.index_to_id(i) << " " << b.x() << " " << b.y() << " " << b.z() << " " << b.r() <<  endl;
		proc << b.x() << " " << b.y() << " " << b.z() << " " << b.r() << " " << endl;
	}
	proc << "WEIGHTS" << endl;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		proc << "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1" << endl;
	}
	proc << "END" << endl << redi::peof;

	int index,ialpha;
	Real val;
	// while( proc >> buf ) result += buf + "\n";
	// std::cout << result << std::endl;

	for( Size a = 1; a <= 20; a++ ) {
		for( Size i = 1; i <= pb.nballs(); i++ ) {
			proc >> ialpha >> index >> val;
			pb.set_surf(i,a,val);
		}
	}


	// output
	map<string,string> lines;
	std::string prefix = basic::options::option[ basic::options::OptionKeys::out::prefix ]();
	system(("mkdir -p "+prefix).c_str());

	Size old_resnum = 0;
	string tag;
	for( Size i = 1; i <= pb.nballs(); i++ ) {
		Size resnum = pb.res_num(i);
		tag = pb.res_name(i);
		if( resnum <= pose.total_residue() && pose.residue(resnum).is_upper_terminus() )	tag += "_Cterm";
		if( resnum <= pose.total_residue() && pose.residue(resnum).is_lower_terminus() )	tag += "_Nterm";
		if( pb.res_num(i) > pose.total_residue() ) {
			tag = "UNK";
			if( resnum == old_resnum ) lines[tag] += "\n" + pb.res_name(i) +" "+ fname +" "+ string_of(resnum) + " ";
		}
		if( old_resnum != resnum ) {
			if( lines[tag].size() ) lines[tag] += "\n";
			old_resnum = resnum;
			lines[tag] += pb.res_name(i) +" "+ fname +" "+ string_of(resnum) + " ";
		}
		lines[tag] += pb.atom_name(i) + " ";
		for( Size a = 1; a <=20; a++ ) lines[tag] += string_of(pb.surf(i,a)) + " ";
	}
	// lines[tag] += "\n";
	std::cerr << "time: " << Real(clock()-t)/Real(CLOCKS_PER_SEC) << std::endl;

	// for( map<string,string>::iterator i = lines.begin(); i != lines.end(); i++ ) {
	// 	cout << "writing: " << i->first << endl;
	// 	utility::io::ozstream out( prefix+'/'+i->first+".lr", std::ios_base::app );
	// 	out << i->second << endl;
	// 	out.close();
	// }

}

int
main (int argc, char *argv[])
{

	try {


	devel::init( argc, argv );

  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace utility;

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

	for( std::map<std::string,utility::io::ozstream*>::iterator i = outs.begin(); i != outs.end(); ++i ) {
		i->second->close();
		delete i->second;
	}

	return 0;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
