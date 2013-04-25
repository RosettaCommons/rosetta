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
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSetCollection.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>

// Devel
#include<devel/init.hh>
#include<devel/helixAssembly/AddResiduesRotamerSetOperation.hh>
#include<devel/helixAssembly/DeleteAllRotamerSetOperation.hh>
#include<devel/helixAssembly/NativeResidueReader.hh>
#include<devel/helixAssembly/BridgeFragmentMover.hh>

// Protocols
#include <protocols/simple_moves/MinPackMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <protocols/jd2/JobDistributor.hh>

// option key includes
#include <basic/options/util.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

// utility
#include <utility/io/izstream.hh>

//local options
namespace basic{ namespace options{ namespace OptionKeys{
basic::options::FileOptionKey const native_residue_files("native_residue_files");
basic::options::FileOptionKey const bridge_fragments("bridge_fragments");
}}}//basic::options::OptionKeys

int
main( int argc, char * argv [] )
{

	try {

	using namespace std;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add( native_residue_files, "File with all native residue files");
	option.add( bridge_fragments, "File containing bridge fragments");
	devel::init(argc, argv);

	if ( !option[ native_residue_files ].user() )  {
		utility::exit("Must provide a native residue file for helix design executable!", 1);
	}
	utility::file::FileName native_res_file( option[ native_residue_files ]() );

	if ( !option[ bridge_fragments ].user() )  {
		utility::exit("Must provide a bridge fragments file for helix design executable!", 1);
	}
	utility::file::FileName bridge_fragment_list_file( option[ bridge_fragments ]() );

	//////////////////////////////////////////////////////////////////////////////////
	///////////////////Close helical bundles using bridge fragments///////////////////
	utility::vector1<utility::file::FileName> bridge_fragment_files;
	utility::io::izstream fragments_stream( bridge_fragment_list_file );
	if ( !fragments_stream.good() ) {
		utility_exit_with_message("unable to open input file file: "+bridge_fragment_list_file.name()+"\n");
	}
	while ( fragments_stream.good() ) {
		std::string name;
		fragments_stream.getline(name);
		if ( fragments_stream.good() ) bridge_fragment_files.push_back( utility::file::FileName(name) );
	}

	core::Size total_bridge_frags(0);
	utility::vector1<core::fragment::FragSetOP> frag_sets;
	for(core::Size i=1; i<=bridge_fragment_files.size(); ++i){
		core::fragment::FragSetOP frag_set(core::fragment::FragmentIO().read_data(bridge_fragment_files[i].name()));
		total_bridge_frags+=frag_set->size();
		frag_sets.push_back(frag_set);
	}
	cout << "Total number of bridge fragments: " << total_bridge_frags << endl;

	protocols::moves::SequenceMoverOP seq_mover = new protocols::moves::SequenceMover;
	BridgeFragmentMoverOP bridge_frag_mover = new BridgeFragmentMover(frag_sets);
	seq_mover->add_mover( bridge_frag_mover );
	///////////////////Done adding bridge fragments///////////////////
	//////////////////////////////////////////////////////////////////




	///////////////////Setup task factory to design using native residue files///////////////////
	NativeResidueReader native_res_reader;
	std::map<core::Size, utility::vector1<core::conformation::ResidueOP> > nat_ro_map =
			native_res_reader.generateResiduesFromFile(native_res_file.name());

	cout << "Finished populating native rotamers map" << endl;

	//create a task factory: this will create a new PackerTask for each input pose
	core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;

	//Delete all rotamers
//	DeleteAllRotamerSetOperation delete_all_rotamers_operation;
//	core::pack::task::operation::AppendRotamerSet delete_rotamers(delete_all_rotamers_operation.clone());
//	main_task_factory->push_back( delete_rotamers.clone() );

	//Add rotamers defined in the native residue file
	for(std::map<core::Size, utility::vector1<core::conformation::ResidueOP> >::const_iterator map_it = nat_ro_map.begin();
			map_it != nat_ro_map.end(); ++map_it){

		AddResiduesRotamerSetOperation nat_ro_set(map_it->second);
		core::pack::task::operation::AppendResidueRotamerSet append_res(map_it->first, nat_ro_set.clone());
		main_task_factory->push_back( append_res.clone() );
	}

	main_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( new core::pack::task::operation::ReadResfile );
	}
	///////////////////Finished task factory setup///////////////////


	//create a ScoreFunction from commandline options (default is score12)
	core::scoring::ScoreFunctionOP score_fxn = core::scoring::getScoreFunction();

	//create the PackRotamersMover which will do the packing
	protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;

	pack_mover->task_factory( main_task_factory );
	pack_mover->score_function( score_fxn );
	seq_mover->add_mover( pack_mover );

	//////MINIMIZATION///////
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_bb(true);
	movemap->set_chi(true);
	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover(
			movemap,
			score_fxn,
			basic::options::option[ basic::options::OptionKeys::run::min_type ].value(),
			0.01,
			true
	);
	seq_mover->add_mover( min_mover );
	//////END MINIMIZATION///////
	cout << "GO GO GO GO GO!!!!" << endl;

	protocols::jd2::JobDistributor::get_instance()->go(seq_mover);

	cout << "-------------DONE-------------" << endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
