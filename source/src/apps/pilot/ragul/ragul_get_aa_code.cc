// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ragul Gowthaman

//GPU enabling is not default
//To test how many threads are fastest for your computer,
//use -gpu:threads 1024 (or other number) on the command line

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>

// Protocol Headers
#include <devel/init.hh>
#include <core/scoring/rms_util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <fstream>
#include <iostream>
#include <string>

using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace basic::options::OptionKeys;
using namespace utility;

OPT_KEY( StringVector, decoy_protein )
OPT_KEY( IntegerVector, list_of_residues )

int main( int argc, char * argv [] ) {

	try{

		NEW_OPT( list_of_residues, "comma seperated list of residue number", 1 );

		devel::init(argc, argv);

		utility::vector1 < core::Size > res_pos = option[list_of_residues]();

		utility::vector1<std::string> decoy_files;

		if ( basic::options::option[in::file::s].user() ) {
			for ( core::Size h=1; h<=option[in::file::s]().size(); h++ ) {
				decoy_files.push_back( option[in::file::s]()[h] );
			}
		} else if ( option[in::file::l].user() ) {
			for ( core::Size h=1; h<=option[in::file::l]().size(); h++ ) {
				utility::io::izstream pdbs(option[in::file::l]()[h]);
				std::string fname;
				while ( pdbs >> fname ) {
					decoy_files.push_back(fname);
				}
			}
		}

		for ( core::Size f=1; f <= decoy_files.size(); f++ ) {
			std::string const input_decoy_name = decoy_files[f];
			std::cout<<"Reading decoy " << input_decoy_name <<std::endl;
			core::pose::Pose decoy_pose;
			core::import_pose::pose_from_pdb(decoy_pose, decoy_files[f]);
			std::cout<<"SEQUENCE "<<input_decoy_name<<" ";
			for ( core::Size g=1; g <= res_pos.size(); g++ ) {
				core::Size h = option[list_of_residues][g];
				std::cout<<decoy_pose.residue(h).name1();
			}
			std::cout<<std::endl;
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
