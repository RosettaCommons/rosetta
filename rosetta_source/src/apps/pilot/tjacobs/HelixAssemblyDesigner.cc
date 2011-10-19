// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file /rosetta/rosetta_source/src/apps/pilot/tjacobs/HelixAssemblyDesigner.cc
/// @brief 
/// @author Tim Jacobs

// Core
#include<core/pack/task/TaskFactory.hh>
#include<core/pack/task/operation/TaskOperations.hh>
#include<core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>

// Devel
#include<devel/init.hh>
#include<devel/helixAssembly/AddResiduesRotamerSetOperation.hh>
#include<devel/helixAssembly/DeleteAllRotamerSetOperation.hh>
#include<devel/helixAssembly/NativeResidueReader.hh>

// Protocols
#include <protocols/moves/MinPackMover.hh>
#include <protocols/moves/PackRotamersMover.hh>
#include <protocols/moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/TaskAwareMinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>

#include <protocols/jd2/JobDistributor.hh>

// option key includes
#include <basic/options/util.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

//local options
namespace basic{ namespace options{ namespace OptionKeys{
basic::options::FileOptionKey const native_residue_files("native_residue_files");
}}}//basic::options::OptionKeys

int
main( int argc, char * argv [] )
{
	using namespace std;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	option.add( native_residue_files, "File with all native residue files").def("");
	devel::init(argc, argv);

	// make symmetric pose if necessary
	if ( !option[ native_residue_files ].user() )  {
		utility::exit("Must provide a native residue file for helix design executable!", 1);
	}
	NativeResidueReader native_res_reader;
	utility::file::FileName native_res_file( option[ native_residue_files ]() );

	std::map<core::Size, utility::vector1<core::conformation::ResidueOP> > nat_ro_map =
			native_res_reader.generateResiduesFromFile(native_res_file.name());

	cout << "Finished populating native rotamers map" << endl;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// end of setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//create a task factory: this will create a new PackerTask for each input pose
	core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;

	//Delete all rotamers
	DeleteAllRotamerSetOperation delete_all_rotamers_operation;
	core::pack::task::operation::AppendRotamerSet delete_rotamers(delete_all_rotamers_operation.clone());
	main_task_factory->push_back( delete_rotamers.clone() );

	//Add rotamers defined in the native residue file
	for(std::map<core::Size, utility::vector1<core::conformation::ResidueOP> >::const_iterator map_it = nat_ro_map.begin();
			map_it != nat_ro_map.end(); ++map_it){

//		cout << "Resnum: " << map_it->first << endl;
//		for(core::Size p=1; p<=map_it->second.size();p++){
//			cout << "Residue: " << map_it->second[p]->name3() << endl;
//			cout << "Num atoms: " << map_it->second[p]->natoms() << endl;
//		}

		AddResiduesRotamerSetOperation nat_ro_set(map_it->second);
		core::pack::task::operation::AppendResidueRotamerSet append_res(map_it->first, nat_ro_set.clone());
		main_task_factory->push_back( append_res.clone() );
	}

	main_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( new core::pack::task::operation::ReadResfile );
	}

	//create a ScoreFunction from commandline options (default is score12)
	core::scoring::ScoreFunctionOP score_fxn = core::scoring::getScoreFunction();

	//create the PackRotamersMover which will do the packing
	protocols::moves::PackRotamersMoverOP pack_mover = new protocols::moves::PackRotamersMover;

//	// Use the symmetric packer if necessary
//	if ( option[ symmetry::symmetry_definition ].user() ) {
//		pack_mover = new protocols::moves::symmetry::SymPackRotamersMover;
//	}

	pack_mover->task_factory( main_task_factory );
	pack_mover->score_function( score_fxn );

//	//This sequence mover will contain packing for sure, and may contain minimization
//	protocols::moves::SequenceMoverOP seq_mover = new protocols::moves::SequenceMover;

//	// make symmetric pose if necessary
//	if ( option[ symmetry::symmetry_definition ].user() )  {
//	    seq_mover->add_mover( new protocols::moves::symmetry::SetupForSymmetryMover );
//	}

//	seq_mover->add_mover( pack_mover );

	cout << "GO GO GO GO GO!!!!" << endl;

//	protocols::jd2::JobDistributor::get_instance()->go(seq_mover);
	protocols::jd2::JobDistributor::get_instance()->go(pack_mover);

	cout << "-------------DONE-------------" << endl;
}
