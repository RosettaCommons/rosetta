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

#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/sequtil.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>

#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>

#include <utility/file/FileName.hh>

#include <protocols/viewer/viewers.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>
using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

using namespace core;

#include "homolog_cst.hh"


// option key includes

#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <utility/excn/Exceptions.hh>


///////////////////////////////////////////////////////////////////////////////
void*
build_template( void* )
{
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::chemical::ResidueTypeSetCAP rsd_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	core::pose::Pose native_pose, fold_pose;
	if ( option[ in::file::native ].user() ) {
		std::string native_fn = option[ in::file::native  ]();
		core::import_pose::pose_from_pdb( native_pose, native_fn );
	}

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	scorefxn->set_weight( atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ]() );

	std::string sequence = core::sequence::read_fasta_file( option[ in::file::fasta ]() )[1];

	core::pose::make_pose_from_sequence( fold_pose, sequence, *rsd_set );
	create_starting_template( fold_pose, option[ OptionKeys::constraints::cst_file ]() );

	fold_pose.dump_pdb( "starting_template.pdb" );

	using namespace optimization;
	kinematics::MoveMap mm;

	core::scoring::ScoreFunctionOP bump_scorefxn( new core::scoring::ScoreFunction() );
	bump_scorefxn->set_weight( vdw, 1.0 );
	bump_scorefxn->set_weight( atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ]() );

	mm.set_bb ( true );

	fold_pose.dump_pdb("before_min.pdb");
	AtomTreeMinimizer().run(
		fold_pose, mm, (*bump_scorefxn), MinimizerOptions( "dfpmin_armijo_nonmonotone", 0.001, true )
	);
	fold_pose.dump_pdb("after_min.pdb");

	return 0;
} // build_template

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {

	// options, random initialization
	devel::init( argc, argv );
	protocols::viewer::viewer_main( build_template );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
