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
/// @author James Thompson

// libRosetta headers


#include <core/types.hh>


#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/func/Func.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <basic/Tracer.hh>

#include <core/scoring/constraints/ConstraintIO.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] )
{
	try {

	// options, random initialization
	devel::init( argc, argv );
	using namespace core::scoring::constraints;
	using namespace basic::options::OptionKeys;
	using namespace basic::options;

  // setup residue types
  core::chemical::ResidueTypeSetCAP rsd_set;
  if ( option[ in::file::fullatom ]() ) {
    rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
  } else {
    rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
  }

	// read in a dummy pose
  core::pose::Pose native_pose;
  if ( option[ in::file::native ].user() ) {
    core::import_pose::pose_from_file( native_pose, *rsd_set, option[ in::file::native ]() , core::import_pose::PDB_file);
  }

  // read in constraints
  ConstraintSetOP cstset;
  std::string cstfile = option[ basic::options::OptionKeys::constraints::cst_file ]();
  cstset = ConstraintIO::get_instance()->read( cstfile, new ConstraintSet, native_pose );
	utility::vector1< ConstraintCOP > csts = cstset->get_all_constraints();

	for ( utility::vector1< ConstraintCOP >::iterator it = csts.begin(), end = csts.end(); it != end; ++it ) {
		FuncOP f = it->get_cst_func();
		f->show_definition( std::cout );
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
