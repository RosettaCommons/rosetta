// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /rosetta/rosetta_source/src/apps/pilot/tjacobs/HelixAssemblyDesigner.cc
/// @brief
/// @author Tim Jacobs

// Core
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/AddResiduesRotamerSetOperation.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSetCollection.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pose/Pose.hh>

// Devel
#include <devel/init.hh>
#include <protocols/sewing/DeleteAllRotamerSetOperation.hh>
#include <protocols/sewing/NativeResidueReader.hh>
#include <protocols/sewing/BridgeFragmentMover.hh>
#include <protocols/sewing/LoopCreationMover.hh>

// Protocols
#include <protocols/simple_moves/MinPackMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/jd2/JobDistributor.hh>

// option key includes
#include <basic/options/util.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// utility
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

namespace HelixAssemblyDesigner {
	basic::options::IntegerOptionKey const num_residues_to_match( "num_residues_to_match" ); // Number of residues to match before and after the jump
	basic::options::IntegerOptionKey const num_helices_in_repeat( "num_helices_in_repeat" ); // Number of residues to match before and after the jump
	basic::options::BooleanOptionKey const autodetect_loops( "autodetect_loops" ); //Try to find loops based on gaps in the protein
	basic::options::FileOptionKey const native_residue_file( "native_residue_file" ); //Try to find loops based on gaps in the protein
}

int
main( int argc, char * argv [] )
{

	try {

	using namespace std;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task;
	using namespace protocols::sewing;

	option.add( HelixAssemblyDesigner::num_residues_to_match, "Number of residues to match before and after the jump");
	option.add( HelixAssemblyDesigner::num_helices_in_repeat, "Number of helices in the repeating unit");
	option.add( HelixAssemblyDesigner::autodetect_loops, "Try to find loops based on gaps in the protein");
	option.add( HelixAssemblyDesigner::native_residue_file, "File of rotamers from native residues");
	devel::init(argc, argv);


	///////////////////Setup task factory to design using native residue files///////////////////
	if ( !option[ HelixAssemblyDesigner::native_residue_file ].user() )  {
		utility_exit_with_message("Must provide the native_residue_file option for helix design executable!");
	}
	utility::file::FileName native_res_file( option[ HelixAssemblyDesigner::native_residue_file ]() );

	NativeResidueReader native_res_reader;
	std::map<core::Size, utility::vector1<core::conformation::ResidueOP> > nat_ro_map =
			native_res_reader.generateResiduesFromFile(native_res_file.name());

	cout << "Finished populating native rotamers map" << endl;

	//create a task factory: this will create a new PackerTask for each input pose
	TaskFactoryOP main_task_factory = new TaskFactory;

	//Delete all rotamers
//	DeleteAllRotamerSetOperation delete_all_rotamers_operation;
//	operation::AppendRotamerSet delete_rotamers(delete_all_rotamers_operation.clone());
//	main_task_factory->push_back( delete_rotamers.clone() );

	//Add rotamers defined in the native residue file
	for(std::map<core::Size, utility::vector1<core::conformation::ResidueOP> >::const_iterator map_it = nat_ro_map.begin();
			map_it != nat_ro_map.end(); ++map_it){

		//Create rotamer set operation from a list of residues
		core::pack::rotamer_set::AddResiduesRotamerSetOperation nat_ro_set(map_it->second);

		//Add AppendResidueRotamerSet task operation to the task factory. This task operation
		//adds the rotamer set to the residue-level task for the given residue
		main_task_factory->push_back(
			new operation::AppendResidueRotamerSet(map_it->first, nat_ro_set.clone()) );
	}

	main_task_factory->push_back( new operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( new operation::ReadResfile );
	}

	///////////////////Setup pack rotamers mover with task factory/////////////////

	//create a ScoreFunction from commandline options
	core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

	protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;

	pack_mover->task_factory( main_task_factory );
	pack_mover->score_function( score_fxn );


	////////////////Setup minmover//////////////////
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_bb(true);
	movemap->set_chi(true);
	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover(
		movemap,
		score_fxn,
		basic::options::option [ basic::options::OptionKeys::run::min_type ].value(),
		0.01,
		true
	);


	///////////////////Setup loop creation mover/////////////////
	core::Real min_rms = option[lh::min_rms];
	core::Real max_rms = option[lh::max_rms];
	core::Size max_radius = option[lh::max_radius];
	utility::vector1<core::Size> loop_sizes = option[lh::loopsizes]();
	core::Size num_residues_to_match = option[HelixAssemblyDesigner::num_residues_to_match].def(3);
	core::Size num_helices_in_repeat = option[HelixAssemblyDesigner::num_helices_in_repeat].def(0);
	bool autodetect_loops = option[HelixAssemblyDesigner::autodetect_loops].def(false);

	protocols::sewing::LoopCreationMoverOP loop_creation_mover;
	if(option[loops::loop_file].user())
	{
		if(autodetect_loops){
			utility_exit_with_message("You cannot provide a loops file and use the autodetect loops features");
		}
		std::string loops_file = option[loops::loop_file]()[1];
		protocols::loops::Loops loops_to_close = protocols::loops::Loops(loops_file);
		loop_creation_mover =
			new protocols::sewing::LoopCreationMover(
				loops_to_close, loop_sizes, max_radius, max_rms, min_rms, num_residues_to_match, num_helices_in_repeat, nat_ro_map);
	}
	else
	{
		if(!autodetect_loops){
			utility_exit_with_message("You must either specify the autodetect_loops flag or provide a loops file with loops::loop_file");
		}
		loop_creation_mover =
			new protocols::sewing::LoopCreationMover(
				loop_sizes, max_radius, max_rms, min_rms, num_residues_to_match, num_helices_in_repeat, nat_ro_map);
	}

	///////////Add all movers and RUN//////////////
	protocols::moves::SequenceMoverOP seq_mover = new protocols::moves::SequenceMover;
	seq_mover->add_mover( pack_mover );
//	seq_mover->add_mover( min_mover );
	seq_mover->add_mover( loop_creation_mover );
//	seq_mover->add_mover( min_mover );

	protocols::jd2::JobDistributor::get_instance()->go(seq_mover);
	cout << "-------------DONE-------------" << endl;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


//	if ( !option[ bridge_fragments ].user() )  {
//		utility::exit("Must provide a bridge fragments file for helix design executable!", 1);
//	}
//	utility::file::FileName bridge_fragment_list_file( option[ bridge_fragments ]() );
//
//	//////////////////////////////////////////////////////////////////////////////////
//	///////////////////Close helical bundles using bridge fragments///////////////////
//	utility::vector1<utility::file::FileName> bridge_fragment_files;
//	utility::io::izstream fragments_stream( bridge_fragment_list_file );
//	if ( !fragments_stream.good() ) {
//		utility_exit_with_message("unable to open input file file: "+bridge_fragment_list_file.name()+"\n");
//	}
//	while ( fragments_stream.good() ) {
//		std::string name;
//		fragments_stream.getline(name);
//		if ( fragments_stream.good() ) bridge_fragment_files.push_back( utility::file::FileName(name) );
//	}
//
//	core::Size total_bridge_frags(0);
//	utility::vector1<core::fragment::FragSetOP> frag_sets;
//	for(core::Size i=1; i<=bridge_fragment_files.size(); ++i){
//		core::fragment::FragSetOP frag_set(core::fragment::FragmentIO().read_data(bridge_fragment_files[i].name()));
//		total_bridge_frags+=frag_set->size();
//		frag_sets.push_back(frag_set);
//	}
//	cout << "Total number of bridge fragments: " << total_bridge_frags << endl;
//
//	BridgeFragmentMoverOP bridge_frag_mover = new BridgeFragmentMover(frag_sets);
//	seq_mover->add_mover( bridge_frag_mover );
//	///////////////////Done adding bridge fragments///////////////////
//	//////////////////////////////////////////////////////////////////
//
//
//
//
//
//	//////END MINIMIZATION///////
//	cout << "GO GO GO GO GO!!!!" << endl;
