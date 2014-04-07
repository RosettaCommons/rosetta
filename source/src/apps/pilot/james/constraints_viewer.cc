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
/// @author James Thompson

// libRosetta headers


#include <core/types.hh>


#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/func/MixtureFunc.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

#include <core/scoring/constraints/util.hh>
#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <basic/Tracer.hh>

#include <core/scoring/constraints/ConstraintIO.hh>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>



int
main( int argc, char * argv [] ) {
	try {

	// options, random initialization
	devel::init( argc, argv );
	using namespace core::scoring::constraints;
	using namespace basic::options::OptionKeys;
	using namespace basic::options;

  // setup residue types
  core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	// read in a dummy pose
  core::pose::Pose native_pose;
  if ( option[ in::file::native ].user() ) {
    core::import_pose::pose_from_pdb(
			native_pose,
			*rsd_set,
			option[ in::file::native ]()
		);
  }

  // read in constraints
  ConstraintSetOP cstset_;
  std::string cstfile = core::scoring::constraints::get_cst_file_option();
  cstset_ = ConstraintIO::read_constraints(
		cstfile,
		new ConstraintSet,
		native_pose
	);

	// print the constraints out to a file.
	std::string filename = option[ out::file::silent ]();
	if ( !utility::file::file_exists( filename ) ) {
		utility::io::ozstream output;
		output.open( filename );
		if ( option[ james::debug ]() ) {
			cstset_->show_definition( output, native_pose );
		} else {
			cstset_->show(output);
		}
		output.close();
	} else {
		utility_exit_with_message( "Error: file exists: " + filename  );
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // main
