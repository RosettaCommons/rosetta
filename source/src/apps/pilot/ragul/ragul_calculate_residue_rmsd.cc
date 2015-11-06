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
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>

#include <basic/Tracer.hh>



using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace basic::options::OptionKeys;

OPT_KEY( StringVector, decoy_protein )
OPT_KEY( String, reference_protein )
OPT_KEY( String, res_num )

int main( int argc, char * argv [] ) {

  try{

		NEW_OPT( decoy_protein, "decoy protein name", "" );
		NEW_OPT( reference_protein, "reference protein name", "reference_protien.pdb" );
		NEW_OPT( res_num, "reference residue number", "1" );

		devel::init(argc, argv);

  std::string const reference_pose_name = option[ reference_protein ];
	std::string const reference_res_num = option[ res_num ];

  core::pose::Pose reference_pose;
  core::import_pose::pose_from_pdb( reference_pose, reference_pose_name );
	core::Size reference_res = core::pose::parse_resnum( reference_res_num, reference_pose);
	core::conformation::Residue const res_res1( reference_pose.conformation().residue( reference_res ) );

	utility::vector1<string> input_decoy_list = option[ decoy_protein ]();
	for (core::Size f=1; f <= input_decoy_list.size(); f++) {
		std::string const input_decoy_name = input_decoy_list[f];
		std::cout<<"Reading decoy " << input_decoy_name <<std::endl;
		core::pose::Pose decoy_pose;
		core::import_pose::pose_from_pdb( decoy_pose, input_decoy_name );
		core::conformation::Residue const res_res2( decoy_pose.conformation().residue( reference_res ) );
		core::Real rmsd (0.0);
		// make sure we're comparing the same amino acid type
		runtime_assert( res_res1.aa() == res_res2.aa() );
		rmsd = core::scoring::automorphic_rmsd( res_res1, res_res2, false /*superimpose*/ );
		std::cout<<"RMSD "<<input_decoy_name<<" "<<rmsd<<std::endl;
	}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
