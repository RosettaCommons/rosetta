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

#include <utility/excn/Exceptions.hh>


#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/Func.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/Tracer.hh>



#include <utility/file/FileName.hh>


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



using namespace std;

int
main( int argc, char* argv [] )
{
    try {
	// options, random initialization
	devel::init( argc, argv );

	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::file;

	// setup residue types
	core::chemical::ResidueTypeSetCAP rsd_set;
	if ( option[ in::file::fullatom ]() ) {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	} else {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	}

	// configure score functions (energy and constraint score functions)
  core::scoring::ScoreFunctionOP scorefxn( core::scoring::getScoreFunction() );
	core::scoring::ScoreFunctionOP constraintscorefxn( new core::scoring::ScoreFunction );
	constraintscorefxn->set_weight( atom_pair_constraint, 1.0 );
	ConstraintSetOP cstset_ = NULL;


	// initialize ConstraintSet (required!)
	std::string cstfile = option[ basic::options::OptionKeys::constraints::cst_file ]();
	if ( ! utility::file::file_exists( cstfile ) ) {
		utility_exit_with_message( "Error: can't read file " + cstfile );
	}

	// initialize native pose (if user specifies to do so)
	core::pose::Pose native_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( native_pose, *rsd_set, option[ in::file::native ]() );
	}

	utility::vector1< FileName > pdb_file_names;
	if ( option[ in::file::s ].user() ) {
		pdb_file_names = option[ in::file::s ]().vector();
	}

	for ( utility::vector1< FileName >::const_iterator it = pdb_file_names.begin(), end = pdb_file_names.end();
				it != end;
				++it ) {
		core::pose::Pose target_pose;
		core::import_pose::pose_from_pdb( target_pose, *rsd_set, (std::string) *it );
		if ( !cstset_ ) {
			cstset_ = ConstraintIO::read( cstfile, new ConstraintSet, target_pose );
		} else {
			target_pose.constraint_set( cstset_ );
		}

		core::Real constraint_score = (*constraintscorefxn)(target_pose);
		core::Real score = (*scorefxn)(target_pose);
		cout << *it << '\t' << constraint_score << '\t' << score << '\n';
	}
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
                                  }
	return 0;
}


