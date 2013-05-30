// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief test app for packer design that is biased toward a multiple sequence alignment/PSSM
/// @author ashworth

#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/io/pdb/pose_io.hh> // pose_from_pdb
#include <basic/options/option.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/dna/setup.hh>

#include <protocols/viewer/viewers.hh>

#include <core/sequence/SequenceProfile.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>

#include <basic/Tracer.hh>
static basic::Tracer TR("apps.pilot.MSA_design");

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh> // file_exists
#include <utility/file/FileName.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
using utility::vector1;

// c++ headers
#include <iostream>

// option key includes
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>


using namespace core;
	using namespace chemical;
	using namespace pack;
		using namespace task;
			using namespace operation;
	using namespace scoring;
		using namespace constraints;
	using namespace sequence;

using namespace protocols;

////////////////////////////////////////////////////////////////////////////////
void
MSA_design(
	pose::PoseOP pose
)
{
	using namespace basic::options;

	TR << '\n' << "*** Starting MSA design ***" << std::endl;

	// register SequenceProfileConstraint with the ConstraintFactory so that it can be constructed from a constraint file
	//ConstraintIO::get_cst_factory().add_type(
	//	new core::scoring::constraints::SequenceProfileConstraint( Size(), utility::vector1< id::AtomID >(), NULL ) );

	// set up scoring
	// this sets BasePartner in pose cacheable data
	// DnaChains is redundant with BasePartner for the time being
	// (BasePartner is simpler, and is used for DNA scoring methods in core)
	scoring::dna::set_base_partner( *pose );

	ScoreFunctionOP scorefxn( getScoreFunction() );

	if ( ! option[ OptionKeys::constraints::cst_file ].user() ) {
		if ( ! option[ OptionKeys::in::file::pssm ].user() ) {
			utility_exit_with_message("a cst file must be specified by the option -constraints::cst_file, or a sequence profile must be specified by the option -in::file::pssm.");
		}
		// add constraints to bias design toward a sequence profile
		SequenceProfileOP profile = new SequenceProfile;
		utility::file::FileName filename( option[ OptionKeys::in::file::pssm ]().front() );
		profile->read_from_checkpoint( filename );

		for ( Size seqpos(1), end( pose->total_residue() ); seqpos <= end; ++seqpos ) {
			// add individual profile constraint for each residue position
			// because of the underlying constraint implementation, this enures that the constraint is a context-independent 1-body energy, or (intra)residue constraint
			pose->add_constraint( new core::scoring::constraints::SequenceProfileConstraint( *pose, seqpos, profile ) );
		}
	} else {
	  std::string cst_file( option[ OptionKeys::constraints::cst_file ]().front() );
	  ConstraintSetOP cst_set =
			ConstraintIO::get_instance()->read_constraints_new( cst_file, new ConstraintSet, *pose );
		pose->constraint_set( cst_set );
	}


	// the PackerTask tells the Packer how to pack sidechains
	TaskFactory tf;
	tf.push_back( new InitializeFromCommandline );
	if ( option[ OptionKeys::packing::resfile ].user() ) {
		tf.push_back( new ReadResfile );
	}
	tf.push_back( new OperateOnCertainResidues( new PreventRepackingRLT, new ResidueHasProperty("DNA") ) );
	PackerTaskCOP ptask = tf.create_task_and_apply_taskoperations( *pose );

	// run the "Packer" (design sidechains)
	Size const packer_runs( option[ OptionKeys::packing::ndruns ]() );
	// containers in which to store the results of packing
	utility::vector1< std::pair< Real, std::string > > results;
	utility::vector1< pose::PoseOP > pose_list;

	pack_rotamers_loop( *pose, *scorefxn, ptask, packer_runs, results, pose_list );

	// output structures
	Size counter(0);
	for ( utility::vector1< pose::PoseOP >::iterator result( pose_list.begin() ),
				end( pose_list.end() ); result != end; ++result ) {
		std::ostringstream ss;
		ss << "MSA_design" << '_' << counter++ << ".pdb";
		(*result)->dump_pdb( ss.str() );
	}
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void*
MSA_design_main( void* )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	TR << "Getting input filename(s)" << '\n';

	typedef vector1< utility::file::FileName > Filenames;
	Filenames pdbnames;

	if ( option[ in::file::l ].user() ) {
		Filenames listnames( option[ in::file::l ]().vector() );
		for ( Filenames::const_iterator filename( listnames.begin() );
		      filename != listnames.end(); ++filename ) {
			utility::io::izstream list( (*filename).name().c_str() );
			while ( list ) {
				std::string pdbname;
				list >> pdbname;
				pdbnames.push_back( pdbname );
			}
		}

	} else if ( option[ in::file::s ].user() ) {
		pdbnames = option[ in::file::s ]().vector();

	} else {
		std::cerr << "No files given: Use either -file:s or -file:l "
		          << "to designate a single pdb or a list of pdbs"
		          << std::endl;
	}

	for ( Filenames::const_iterator filename( pdbnames.begin() );
	      filename != pdbnames.end(); ++filename ) {
		std::cout << "Input pdb file " << *filename;
		if ( !utility::file::file_exists( *filename ) ) {
			std::cout << " not found, skipping" << std::endl;
			continue;
		}
		std::cout << std::endl;
		pose::PoseOP pose = new pose::Pose;
		protocols::viewer::add_conformation_viewer( pose->conformation(), "Main Pose" );
		core::import_pose::pose_from_pdb( *pose, *filename );
		MSA_design( pose );
	}

	return 0;
}

int
main( int argc, char * argv[] )
{
	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	devel::init( argc, argv );
	protocols::viewer::viewer_main( MSA_design_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
}
