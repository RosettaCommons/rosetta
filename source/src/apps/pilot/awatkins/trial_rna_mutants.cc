// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Take a vector of name3s and positions. Do all combinations.


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/rna/RNA_ResidueLevelTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <protocols/rna/denovo/movers/RNA_Minimizer.hh>
#include <core/pack/pack_rotamers.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/init/init.hh>

#include <core/pose/rna/util.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>


#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

//RNA stuff.
#include <protocols/rna/denovo/util.hh>
//#include <protocols/simple_moves/MutateResidue.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( IntegerVector, seqpos )
OPT_KEY( String, residue_types )

///////////////////////////////////////////////////////////////////////////////
void
rna_design_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	ScoreFunctionOP scorefxn = get_score_function();

	pose::Pose pose;
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );
	std::string pdb_file  = option[ in::file::s ][1];
	core::import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
	protocols::rna::denovo::ensure_phosphate_nomenclature_matches_mini( pose );

	// ignore -- unused
	utility::vector1< std::pair< Real, std::string > > results;
	utility::vector1< pose::PoseOP > pose_list;

	dump_pdb( pose, "start.pdb");
	pose::Pose save_pose( pose );
	utility::vector1< Size> sequence_positions = option[ seqpos ].value();
	utility::vector1< std::string > restypes = utility::string_split( option[ residue_types ].value(), ',' );

	// scorefxn->energy_method_options().exclude_DNA_DNA( exclude_DNA_DNA );
	methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
	options.exclude_DNA_DNA( false );
	scorefxn->set_energy_method_options( options );
	scorefxn->show( std::cout,pose );
	pose.dump_pdb( "start.pdb" );

	for ( Size const pos : sequence_positions ) {
		for ( std::string const & restype : restypes ) {
			pose = save_pose;

			// First, mutate residue.
			core::pose::rna::mutate_position( pose, pos, restype );
			std::stringstream ss;
			ss << pos << "_" << restype << ".pdb";
			//pose.dump_scored_pdb( ss.str(), *scorefxn );

			// Then, pack
			pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));

			for ( Size ii = 1; ii <= pose.size(); ++ii ) {
				if ( ii == pos ) continue;
				auto const & r1 = pose.residue( pos );
				auto const & r2 = pose.residue( ii );
				task->nonconst_residue_task( ii ).or_ex1( true );
				task->nonconst_residue_task( ii ).or_ex4( true );
				if ( r1.xyz( r1.nbr_atom() ).distance_squared( r2.xyz( r2.nbr_atom() ) ) < 15 * 15 ) continue;

				task->nonconst_residue_task( ii ).prevent_repacking();
			}

			task->nonconst_residue_task( pos ).allow_noncanonical_aa( restype );
			task->nonconst_residue_task( pos ).and_extrachi_cutoff( 0 );
			task->nonconst_residue_task( pos ).nonconst_rna_task().set_sample_rna_chi( true );
			task->nonconst_residue_task( pos ).or_ex4( true );
			task->nonconst_residue_task( pos ).or_ex1( true );

			pack::pack_rotamers( pose, *scorefxn, task );
			//pose.dump_scored_pdb( ss.str(), *scorefxn );

			//pack::pack_rotamers_loop( pose, *scorefxn, task, 1, results, pose_list);

			// Finally, minimize (input probably was)
			protocols::rna::denovo::movers::RNA_Minimizer min;
			min.set_score_function(scorefxn);
			min.apply( pose );
			pose.dump_scored_pdb( ss.str(), *scorefxn );


			std::cout << pos << " " << restype << " " << ( *scorefxn )( pose ) << std::endl;

		}
	}
}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	rna_design_test();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -s <pdb file> [ -resfile <resfile>] " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		//Uh, options? MOVE THESE TO OPTIONS NAMESPACE INSIDE CORE/OPTIONS.
		NEW_OPT( seqpos, "sequence positions to examine", 0);
		NEW_OPT( residue_types, "string of residue type name3s", "");

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
