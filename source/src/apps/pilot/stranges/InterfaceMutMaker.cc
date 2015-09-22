// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/stranges/InterfaceMutMaker.cc
/// @brief Makes a mutation from a command line flag and outputs wt and mut complex and separated structures
/// @author Ben Stranges

//Format of mutation string is as follows Chain,WT,#,Mut
//  -InterfaceMutMaker::mutation A,E,45,T  B,N,22,W

#include <utility/io/izstream.hh>
// Unit headers
#include <devel/init.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/conformation/Conformation.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/SetReturningPackRotamersMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/format.hh>

// C++ headers
#include <sstream>
#include <iostream>
#include <string>
#include <stdlib.h>

#include <core/import_pose/import_pose.hh>

static THREAD_LOCAL basic::Tracer TR( "InterfaceMutMaker" );

using namespace core;
using namespace utility;
using namespace protocols;
using namespace core::pack::task;
using namespace core::pack::task::operation;
using namespace basic::options;
using namespace basic::options::OptionKeys;

//  Index for location of tokens in a defined mutation
typedef enum Index{
  chain_index = 1,
  wt_index = 2,
  res_index = 3,
  mut_index = 4
};

// application specific options
namespace InterfaceMutMaker {
	BooleanOptionKey const allow_rbmin( "InterfaceMutMaker::allow_rbmin" );
	BooleanOptionKey const output_all_packruns ( "InterfaceMutMaker::output_all_packruns" );
	StringVectorOptionKey const mutation( "InterfaceMutMaker::mutation" );
	StringOptionKey const mutation_file( "InterfaceMutMaker::mutation_file" );
	//IntegerOptionKey const rbjump( "InterfaceMutMaker::rbjump" );
	StringVectorOptionKey const movechains( "InterfaceMutMaker::movechains" );
}

// //calc rmsd helper function
// core::Real calc_rmds( core::pose::Pose & pose1, core::pose::Pose & pose2){
// 	Real rms_value(0.0);
// }

//helper function to read in a mutation string from a vector of strings
void tokenize_string( const vector1< std::string > & str_vector, utility::vector1< utility::vector1< std::string > > & tokens, const std::string & delimiters = " " ) {

  tokens.resize( str_vector.size() );
  for( Size ii = 1; ii <= str_vector.size(); ++ii){
    std::string str ( str_vector[ii] );

    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);

    // Find first "non-delimiter".
    std::string::size_type pos = str.find_first_of(delimiters, lastPos);

    while ( std::string::npos != pos || std::string::npos != lastPos ) {
      // Found a token, add it to the vector.
      tokens[ii].push_back( str.substr( lastPos, pos - lastPos ) );

      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of( delimiters, pos );

      // Find next "non-delimiter"
      pos = str.find_first_of( delimiters, lastPos );
    }
  }
}

//reorder the fold tree and set the new jump
pose::Pose reorder_foldtree(core::pose::Pose & pose, utility::vector1<std::string> & mobile_chains_string,
														core::Size & mobile_jump){
	using core::kinematics::Edge;
	//parse the mobile chains to figure out chain nums
	std::set<Size> mobile_chains;
	for(Size j = 1; j <= mobile_chains_string.size(); ++j){
		char this_chain (mobile_chains_string[ j ][0]);
	 	for (Size i = 1; i<=pose.total_residue(); ++i){
			if (pose.pdb_info()->chain( i ) == this_chain){
				mobile_chains.insert( pose.chain(i) );
				break;
			}
		}
		TR << "Mobile chains are: "  << this_chain << ", ";
	}
	TR << "these will be moved together." << std::endl;
	if(mobile_chains.empty())
		utility_exit_with_message_status( "Can't find given movechain.  Exiting...\n", 1 );
	//now get info about chains
	Size numchains ( pose.conformation().num_chains() );
	vector1<Size> chain_starts;
	vector1<Size> chain_ends;
	for(Size ii = 1; ii<= numchains; ++ii){
		chain_starts.push_back( pose.conformation().chain_begin( ii ) );
		chain_ends.push_back( pose.conformation().chain_end( ii ) );
	}

	//find a non mobile chain
	Size anchor_chain (1); //switch later if needed
	for(Size ii = 1; ii<= numchains; ++ii){
		if ( !mobile_chains.count( ii ) ){
			anchor_chain = ii;
			break;
		}
	}

	//now build our fold tree
	//this will move the mobile chains together
	kinematics::FoldTree foldtree (pose.fold_tree());
	Size this_jump(1);
	foldtree.clear();
	bool previous_mobile (false);
	for( Size ii=1; ii<= numchains; ++ii){
		foldtree.add_edge( Edge(chain_starts[ii], chain_ends[ii], kinematics::Edge::PEPTIDE) );
		if (ii != anchor_chain) {
			if ( mobile_chains.count(ii) && previous_mobile  ){ //not in the first mobile chain
				foldtree.add_edge( Edge(chain_starts[ii], chain_ends[*(mobile_chains.begin())], this_jump) );
			}
			else if ( mobile_chains.count(ii) ){ //ie in the first mobile chain
				foldtree.add_edge( Edge(chain_starts[ii], chain_ends[anchor_chain], this_jump) );
				mobile_jump = this_jump; //sets this jump to be the new mobile one
				previous_mobile=true;
			}
			else foldtree.add_edge( Edge(chain_starts[ii], chain_ends[anchor_chain], this_jump) );
			++this_jump;
		}
	}
	foldtree.reorder(chain_starts[anchor_chain]);
	pose.fold_tree(foldtree);
	TR << "Remade foldtree:\n"<< foldtree << std::endl;

	return pose;

}//end reorder foldtree


// move the pose apart by the jump and return the new pose
pose::Pose split_chains( core::pose::Pose & pose, scoring::ScoreFunctionOP & scorefxn, Size & jump ){

	pose::Pose copy(pose);

	protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( copy, jump ) );
	translate->step_size( 300.0 );
	translate->apply( copy );

  ((*scorefxn)(copy));

  return copy;
}


//begin main
int
main( int argc, char* argv[] ) {

	try {


  option.add( InterfaceMutMaker::mutation, "The mutation to make to the wild-type pose. Format to use is chain,old-residue-type,PDB-res-num,new-residue-type." );
  option.add( InterfaceMutMaker::mutation_file, "The mutation to make to the wild-type pose. Format is a white space or return separated text file. chain,old-residue-type,PDB-res-num,new-residue-type." );
	//option.add( InterfaceMutMaker::rbjump, "Defines the rigid body jump(s) to split the interface apart by.");
	//option.add( InterfaceMutMaker::n_outputs, "Defines the number of structures to output for each mutation.  Also defines the number of final minimization runs.");
	option.add( InterfaceMutMaker::allow_rbmin, "Allow rigid body minimization." );
	option.add( InterfaceMutMaker::output_all_packruns, "Output all runs of SetReturningPackRotamersMover." );
	option.add( InterfaceMutMaker::movechains, "Which chain(s) is/are moved away from the others." );

	//initialize
  devel::init( argc, argv );

	//inputs
	Size pack_cycles = option[packing::ndruns].def(1);
	bool allow_rbmin = option[ InterfaceMutMaker::allow_rbmin].def(true);
	bool output_all = option[ InterfaceMutMaker::output_all_packruns].def(false);

	TR << "Number of packing runs: "<< pack_cycles << std::endl;
	TR << "Allow rigid body minimization: " << allow_rbmin << std::endl;
	TR << "Output all packing runs: " << output_all << std::endl;

  // read in mutation
  vector1< vector1< std::string > > parsed_tokens;
  vector1< std::string > mutation_string_vector;  //fill later

	// read file if there is one...
	std::string mutation_filename;  //fill later
	vector1< std::string > mutation_file_vector; //fill below
	if (option[ InterfaceMutMaker::mutation_file ].active() ){
		mutation_filename = option[ InterfaceMutMaker::mutation_file ];
#ifndef NDEBUG
	TR << "Mutation file:  " << mutation_filename << std::endl;
#endif
		utility::io::izstream mutation_file ( mutation_filename );
		std::string line, mut_name;
		while ( getline( mutation_file , line) ){
			std::istringstream stream (line);
			while( stream >> mut_name ){
				//push back mut_name into a vector
				mutation_file_vector.push_back( mut_name );
			}
		}
	} //end read file

	//check to make sure a mutation is defined and fill the vector
	if (  option[ InterfaceMutMaker::mutation ].active() && option[ InterfaceMutMaker::mutation_file ].active() ) {
		mutation_string_vector = option[ InterfaceMutMaker::mutation ];
		for (Size ii = 1; ii <= mutation_file_vector.size(); ++ii){
			mutation_string_vector.push_back( mutation_file_vector[ii] );
		}
	}
	else if (option[ InterfaceMutMaker::mutation_file ].active() ) {
		mutation_string_vector = mutation_file_vector;
	}
	else if (option[ InterfaceMutMaker::mutation ].active() ){
		mutation_string_vector = option[ InterfaceMutMaker::mutation ];
	}
	else utility_exit_with_message_status( "No mutation given!  Exiting...\n", 1 );

	// tokenize the mutation string vector
	tokenize_string( mutation_string_vector, parsed_tokens, "," );

  vector1< std::string > chains_name;
  vector1< std::string > wt_residues;
  vector1< std::string > mut_residues;
  vector1< int > resid;

  for ( Size ii = 1; ii <= parsed_tokens.size(); ++ii ){
    vector1< std::string > this_token ( parsed_tokens[ii] );
    chains_name.push_back( this_token[ chain_index ] );
    wt_residues.push_back( this_token[ wt_index ] );
    mut_residues.push_back( this_token[ mut_index ] );
    int resnum;
    std::istringstream ss( this_token[ res_index ] );
    ss >> resnum;
    resid.push_back( resnum );

    TR << "Finished parsing mutation string. Found mutation of chain: "<< this_token[ chain_index ] <<" residue: '" << this_token[ wt_index ] << "' to '" <<this_token[ mut_index ] << "' at PDB resid: " << this_token[ res_index ] << std::endl;
  }

	//set up stringstream of the mutation(s)
	std::stringstream mut_stringstream_output;
	for(Size jj = 1; jj <= resid.size(); ++jj){
		mut_stringstream_output << chains_name[jj] << wt_residues[jj]
														<< resid[jj] << mut_residues[jj];
		if( jj != resid.size())
			mut_stringstream_output << ".";
	}
	std::string mut_string_output ( mut_stringstream_output.str() );


  //start getting score and pose inputs

  utility::file::FileName pdb_file_name;
  if ( option[ in::file::s ].active() ) {
    pdb_file_name = utility::file::FileName( basic::options::start_file() );
  }
	else utility_exit_with_message_status( "If you want to make a mutant structure it's probably a good idea to include a structure in input.  Use -s flag ONLY.  Exiting...\n", 1 );//add error checking later

	//define a native pose if given one for RMSD calc later...
	pose::Pose native_pose;

	if (basic::options::option[ in::file::native ].user()){
			core::import_pose::pose_from_pdb( native_pose, basic::options::option[ in::file::native ]());
	}

  scoring::ScoreFunctionOP scorefxn;
  scorefxn = scoring::get_score_function();
  scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn->energy_method_options() );
  energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
  scorefxn->set_energy_method_options( energymethodoptions );

  pose::Pose pose;
  core::import_pose::pose_from_pdb( pose, pdb_file_name.name() );

  // score the input pose for kicks and giggles
  (*scorefxn)( pose );
	if(pose.conformation().num_chains() < 2)
		utility_exit_with_message_status( "This protocol will only work on pdbs with >= 2 chains. Exiting\n", 1 );

	//debug ounly output
#ifndef NDEBUG
	kinematics::FoldTree foldtree (pose.fold_tree());
	TR << "Number of jumps: "<< pose.num_jump() << " Foldtree for input pose: \n"
			 << foldtree << std::endl;
	for(Size ii = 1; ii <= pose.num_jump(); ++ii ){
		TR << "Jump: " << ii << " Downstream: "<< pose.pdb_info()->chain( foldtree.downstream_jump_residue(ii) )
			 << " Upstream: " << pose.pdb_info()->chain( foldtree.upstream_jump_residue(ii) ) << std::endl;
	}
	for(Size jj=1; jj <= pose.conformation().num_chains(); ++jj){
		Size chainstart(pose.conformation().chain_begin( jj ) );
		Size chainend(pose.conformation().chain_end( jj ) );
		TR << "Chain: "<< jj << " start: " << chainstart << " (" <<pose.pdb_info()->pose2pdb(chainstart) << ")"
			 << " end: " << chainend << " (" <<pose.pdb_info()->pose2pdb(chainend) << ")" << std::endl;
	}
#endif

	//reorder the fold tree and set the new rbjump
	Size jump_num;
	vector1<std::string> move_chains;
	if(pose.conformation().num_chains() > 2){
		if ( option[InterfaceMutMaker::movechains].active() ){
			move_chains = option[InterfaceMutMaker::movechains].value();
			TR << "Movable chains are: ";
			for( Size ii = 1; ii<=move_chains.size(); ++ii){
				TR<< move_chains[ii] << ", ";
			}
			TR << std::endl;
			//actual function call
			pose = reorder_foldtree(pose, move_chains, jump_num);
			TR << "More than 2 chains present, setting jumpnum to: " << jump_num << std::endl;
		}
		else
			utility_exit_with_message_status( "More than 2 chains found.  Need to define which chains are mobile with -movechains flag.  Exiting...\n", 1 );
	}
	else {
		jump_num = 1;
		TR << "Only 2 chains present, setting jump to: " << jump_num << std::endl;
	}

	//more debugging output
#ifndef NDEBUG
 foldtree = pose.fold_tree();
	TR << "New fold_tree: \n"
		 << foldtree << "New Jump: "<< jump_num << std::endl;
	for(Size ii = 1; ii <= pose.num_jump(); ++ii ){
		TR << "Jump: " << ii << " Downstream: "<< pose.pdb_info()->chain( foldtree.downstream_jump_residue(ii) )
			 << " Upstream: " << pose.pdb_info()->chain( foldtree.upstream_jump_residue(ii) ) << std::endl;
	}
	for(Size jj=1; jj <= pose.conformation().num_chains(); ++jj){
		Size chainstart(pose.conformation().chain_begin( jj ) );
		Size chainend(pose.conformation().chain_end( jj ) );
		TR << "Chain: "<< jj << " start: " << chainstart << " (" <<pose.pdb_info()->pose2pdb(chainstart) << ")"
			 << " end: " << chainend << " (" <<pose.pdb_info()->pose2pdb(chainend) << ")" << std::endl;
	}
#endif

  //rewrite resid vector in terms of pose numbering
	//also make a set that contains the mutation residues
	std::set< Size > mutation_set;
  vector1< Size > pose_resid;
  for ( Size ii = 1; ii <= chains_name.size(); ++ii ) {
    if ( chains_name[ii].length() != 1 )
      	utility_exit_with_message( "Chain name is too big! \n" );
    //should return only the chain char  hopefully no error check needed
		std::string this_chain ( chains_name[ii] );
    char chain_char (this_chain[0]);
    Size pose_resnum ( pose.pdb_info()->pdb2pose(chain_char, resid[ii]) );
    if ( pose_resnum == 0 )
			utility_exit_with_message( "Input mutation chain/residue not found! \n" );
    pose_resid.push_back( pose_resnum );
		mutation_set.insert( pose_resnum );

		//check to make sure the wt residue is actually right, uses enum compare
		char this_aa (wt_residues[ ii ][0]);
		if( pose.residue( pose_resnum ).aa() != core::chemical::aa_from_oneletter_code( this_aa ) ){
			TR << "ERROR:  For input: "<< pdb_file_name.base() << " (residue chain) #: "
				 << pose.pdb_info()->pose2pdb( pose_resnum ) << " AA: "
				 <<core::chemical::name_from_aa( pose.residue( pose_resnum ).aa() ) << ".  You entered AA: "
				 <<core::chemical::name_from_aa( core::chemical::aa_from_oneletter_code( this_aa ) ) << std::endl;
			utility_exit_with_message( "Incorrect residue type found.\n" );
		}
  }


  TR << "Input mutations are correct!" << std::endl;


	////////////////////////////////////////////////
	// Now set up task factories for WT and Mut
	////////////////////////////////////////////////
	pose::Pose wt_pose (pose), mut_pose(pose);

	TaskFactoryOP wt_tf = new TaskFactory();

	InitializeFromCommandlineOP init_op = new InitializeFromCommandline();
	//IncludeCurrentOP ic_op = new IncludeCurrent();
	//ReadResfileOP resfile_op = new ReadResfile();

	//set common taskfactory options
	wt_tf->push_back( init_op );
	//	wt_tf->push_back( resfile_op );
	wt_tf->push_back( new operation::IncludeCurrent );

	//calculator to figure out set of neighbors around mutation
	std::string const nb_calc("complex_neighboorhood_calculator");
	pose::metrics::CalculatorFactory::Instance().register_calculator( nb_calc, new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( mutation_set ) );
	basic::MetricValue< std::set< Size > > neighbor_mv;
	mut_pose.metric( nb_calc, "neighbors", neighbor_mv);
	std::set< Size > const neighbor_set ( neighbor_mv.value() );

	//set up different tasks for mutation
	PreventRepackingOP wt_prevent_op = new PreventRepacking();
	RestrictResidueToRepackingOP wt_repack_op = new RestrictResidueToRepacking();

	// prevent packing on all residues not in neighbor set
	// restrict to repacking on all in neighbor set
	for( Size ii = 1; ii<= mut_pose.n_residue(); ++ii) {
		if( neighbor_set.count( ii ) )
			wt_repack_op->include_residue( ii );
		else
			wt_prevent_op->include_residue( ii );
	}

	//fill task factory with these restrictions
	wt_tf->push_back( wt_prevent_op ) ;
	wt_tf->push_back( wt_repack_op )  ;


#ifndef NDEBUG
	TR<< "Neighbors to mutations are: \n";
	for( std::set< core::Size >::const_iterator it(neighbor_set.begin()), end(neighbor_set.end());
			 it != end; ++it){
		TR << *it << ", ";
	}
	TR << std::endl;
#endif

	//need a special PackerTask to make the mutation before minimization
	//doesn't move any other residues, just changes the designated one
	//stolen from Steven
	core::pack::task::PackerTaskOP muttask(core::pack::task::TaskFactory::create_packer_task((mut_pose)));
	utility::vector1_bool packable(mut_pose.total_residue(), false); //false = nobody is packable

	for( Size jj = 1; jj<= pose_resid.size(); ++jj){
		//pose_resid and mut_residues have have refs to residue and aa of mutation
		utility::vector1< bool > allowed_aa( chemical::num_canonical_aas, false );
		allowed_aa[ core::chemical::aa_from_oneletter_code( mut_residues[jj][0] ) ] = true;
		muttask->nonconst_residue_task( pose_resid[jj] ).restrict_absent_canonical_aas(allowed_aa);
		//sample rotamers here to avoid conflicts with MinMover later on
		muttask->nonconst_residue_task( pose_resid[jj] ).or_ex1( true );
		muttask->nonconst_residue_task( pose_resid[jj] ).or_ex2( true );
		muttask->nonconst_residue_task( pose_resid[jj] ).or_ex3( true );
		muttask->nonconst_residue_task( pose_resid[jj] ).or_ex4( true );
		packable[ pose_resid[jj] ] = true;
	}
	muttask->restrict_to_residues(packable);  // prevents non mutant residues from moving

#ifndef NDEBUG
	TR << "Mutation task to apply: \n"<< *muttask << std::endl;
#endif

	//apply mutation
	protocols::simple_moves::PackRotamersMoverOP packrot_mover( new protocols::simple_moves::PackRotamersMover (scorefxn, muttask) );
	packrot_mover->apply ( mut_pose );  //makes the passed in mutation


#ifndef NDEBUG
	TR<< "WT Packer Task: " << *(wt_tf->create_task_and_apply_taskoperations(wt_pose));
	TR<< "MUT Packer Task: " << *(wt_tf->create_task_and_apply_taskoperations(mut_pose));

#endif

	TR << "Finished creating all TaskOperations and TaskFactories and making the mutations." << std::endl;


	/////////////////////////////////////////////////////////
	// Now set up the move map for minimization
	/////////////////////////////////////////////////////////
	//need two MoveMaps because RB minimization of the split apart complexs is stupid.
	kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
	kinematics::MoveMapOP movemap_sc = new core::kinematics::MoveMap();
	kinematics::MoveMapOP movemap_split = new core::kinematics::MoveMap();
	//allow bb and sc minimization for the neighbors of the mutations
	for( std::set< core::Size >::const_iterator it(neighbor_set.begin()), end(neighbor_set.end());
			 it != end; ++it){
		movemap_sc->set_chi(*it, true);
		movemap->set_bb( *it, true );
		movemap_split->set_bb( *it, true );
		movemap->set_chi( *it, true );
		movemap_split->set_chi( *it, true );
	}
	// allow RB minimization for the defined jump in the complex only
	// under option control for mutations that may push the complex away in minimization

	if(allow_rbmin)
		movemap->set_jump( jump_num, true);
	else movemap->set_jump( jump_num, false);

	// Now that we have a mutation, split the poses up
	pose::Pose sep_wt_pose (wt_pose);
	pose::Pose sep_mut_pose (mut_pose);
	sep_wt_pose = split_chains( sep_wt_pose, scorefxn, jump_num );
	sep_mut_pose = split_chains( sep_mut_pose, scorefxn, jump_num );


	/////////////////////////////////////////////////////////
	// Make MinMovers and PackRotsMover
	////////////////////////////////////////////////////////

	//first do minimization mover making
	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( movemap, scorefxn, option[ OptionKeys::run::min_type ].value(), 0.01, true /*use_nblist*/ );
	protocols::simple_moves::MinMoverOP min_mover_sc = new protocols::simple_moves::MinMover( movemap_sc, scorefxn, option[ OptionKeys::run::min_type ].value(), 0.01, true /*use_nblist*/ );
	protocols::simple_moves::MinMoverOP min_mover_split = new protocols::simple_moves::MinMover( movemap_split, scorefxn, option[ OptionKeys::run::min_type ].value(), 0.01, true /*use_nblist*/ );

	TR << "Minimizing with: " << option[ OptionKeys::run::min_type ].value() << "." << std::endl;
	//make vectors to hold repacked poses
	vector1< pose::Pose > packed_wt_complex_poses = vector1< pose::Pose >( pack_cycles );
	vector1< pose::Pose > packed_wt_separate_poses = vector1< pose::Pose >( pack_cycles );
	vector1< pose::Pose > packed_mut_complex_poses = vector1< pose::Pose >( pack_cycles );
	vector1< pose::Pose > packed_mut_separate_poses = vector1< pose::Pose >( pack_cycles );


	// Make PackRots movers
	protocols::simple_moves::SetReturningPackRotamersMoverOP wt_repacker = new protocols::simple_moves::SetReturningPackRotamersMover( pack_cycles );
	wt_repacker->task_factory( wt_tf );
	wt_repacker->score_function( scorefxn );


	///////////////////////////////////////////////////////////////////////////////
	//  Now do actual applys  (wahoo!!!)
	///////////////////////////////////////////////////////////////////////////////

	//do we need monte carlo on these?

	//Print out initial scores
	TR << "Initial scores for inputs: \n"
		 << pdb_file_name.base() << "_complex_wt_" << mut_string_output  <<": "<< (*scorefxn)(wt_pose) << "\n"
		 << pdb_file_name.base() << "_separate_wt_" << mut_string_output  <<": " <<(*scorefxn)(sep_wt_pose) << "\n"
		 << pdb_file_name.base() << "_complex_mut_" << mut_string_output  <<": " <<(*scorefxn)(mut_pose) << "\n"
		 << pdb_file_name.base() << "_separate_mut_" << mut_string_output  <<": " <<(*scorefxn)(sep_mut_pose)
		 << std::endl;


	// Step 1 side chain only minimization...
	min_mover_sc->apply(wt_pose);
	min_mover_sc->apply(mut_pose);
	min_mover_sc->apply(sep_wt_pose);
	min_mover_sc->apply(sep_mut_pose);

	//Print out scores after minimization
	TR << "Scores after Minimization \n"
		 << pdb_file_name.base() << "_complex_wt_" << mut_string_output  <<": " << (*scorefxn)(wt_pose) << "\n"
		 << pdb_file_name.base() << "_separate_wt_" << mut_string_output  <<": " << (*scorefxn)(sep_wt_pose) << "\n"
		 << pdb_file_name.base() << "_complex_mut_" << mut_string_output  <<": " << (*scorefxn)(mut_pose) << "\n"
		 << pdb_file_name.base() << "_separate_mut_" << mut_string_output  <<": " << (*scorefxn)(sep_mut_pose)
		 << std::endl;


	//Step 2 packing in case minimization not enough in some cases (big clashes)
	//set up vector of repacked poses
	TR << "Packing results for: "<< pdb_file_name.base() << "_complex_wt_" << mut_string_output   << std::endl;
	wt_repacker->apply(wt_pose);
	wt_repacker->get_repacked_poses( packed_wt_complex_poses );
	TR << "Packing results for: "<< pdb_file_name.base() << "_separate_wt_" << mut_string_output   << std::endl;
	wt_repacker->apply(sep_wt_pose);
	wt_repacker->get_repacked_poses( packed_wt_separate_poses );
	TR << "Packing results for: "<< pdb_file_name.base() << "_complex_mut_" << mut_string_output   << std::endl;
	wt_repacker->apply(mut_pose);
	wt_repacker->get_repacked_poses(  packed_mut_complex_poses );
	TR << "Packing results for: "<< pdb_file_name.base() << "_separate_mut_" << mut_string_output   << std::endl;
	wt_repacker->apply(sep_mut_pose);
	wt_repacker->get_repacked_poses(  packed_mut_separate_poses );


//Print out scores after packing
	TR << "Best scores after Packing \n"
		 << pdb_file_name.base() << "_complex_wt_" << mut_string_output  <<": " <<(*scorefxn)(wt_pose) << "\n"
		 << pdb_file_name.base() << "_separate_wt_" << mut_string_output  <<": " <<(*scorefxn)(sep_wt_pose) << "\n"
		 << pdb_file_name.base() << "_complex_mut_" << mut_string_output  <<": " <<(*scorefxn)(mut_pose) << "\n"
		 << pdb_file_name.base() << "_separate_mut_" << mut_string_output  <<": " <<(*scorefxn)(sep_mut_pose)
		 << std::endl;

	//Step 3 //////////////////////////////////////////////////////////////////
	// One more round of minimization and output.
	// Just outputs the best result from packing unless otherwise specificed by output_all_packruns
	std::string wt_complex_filename = pdb_file_name.base() + "_complex_wt_" + mut_string_output  + ".pdb";
	std::string mut_complex_filename = pdb_file_name.base() + "_complex_mut_" + mut_string_output  + ".pdb";
	std::string wt_split_filename = pdb_file_name.base() + "_separate_wt_" + mut_string_output + ".pdb";
	std::string mut_split_filename = pdb_file_name.base() + "_separate_mut_" + mut_string_output +  ".pdb";
	//apply minimization
	min_mover->apply( wt_pose );
	min_mover->apply( mut_pose );
	min_mover_split->apply( sep_wt_pose );
	min_mover_split->apply( sep_mut_pose );

	//figure out an rmsd if the native structure is given
	core::Real  wt_rms (0.0);
	core::Real  mut_rms (0.0);
	if (basic::options::option[ in::file::native ].user()){
		// only calculate on complexes
		// allow superposition because RB min is allowed
		wt_rms = scoring::rmsd_with_super( native_pose, wt_pose , scoring::is_protein_CA ) ;
		mut_rms = scoring::rmsd_with_super( native_pose, mut_pose, scoring::is_protein_CA ) ;
	}

	//print out tracer output
	TR << "Scores of best after minimization: "<<" \n"
		 << wt_complex_filename <<": "<<(*scorefxn)( wt_pose ) <<  "  rms: "<< wt_rms<<"\n"
		 << wt_split_filename <<": " <<(*scorefxn)( sep_wt_pose ) <<  "\n"
		 << mut_complex_filename  <<": " <<(*scorefxn)( mut_pose ) <<  "  rms: "<< mut_rms<<"\n"
		 << mut_split_filename << ": " <<(*scorefxn)( sep_mut_pose ) << std::endl;
	//calc ddG mutant - wt
	TR << "DDGbind " <<pdb_file_name.base()<<" "<< mut_string_output << " = "
		 << ( (*scorefxn)( mut_pose ) - (*scorefxn)( sep_mut_pose ) ) - ((*scorefxn)( wt_pose ) -(*scorefxn)( sep_wt_pose ) ) << std::endl;

	//output pdbs!
	wt_pose.dump_scored_pdb( wt_complex_filename, *(scorefxn()) );
	mut_pose.dump_scored_pdb( mut_complex_filename, *(scorefxn()) );
	sep_wt_pose.dump_scored_pdb( wt_split_filename, *(scorefxn()) );
	sep_mut_pose.dump_scored_pdb( mut_split_filename, *(scorefxn()) );

	// Dump all the packing results if desired.
	if ( output_all ){
		TR << "Outputting all poses from packing."<< std::endl;
		// minimize all the poses from packing
		for(Size ii=1; ii<= packed_wt_complex_poses.size(); ++ii){
			min_mover->apply( packed_wt_complex_poses[ ii ] );
			min_mover->apply( packed_mut_complex_poses[ ii ]);
			min_mover_split->apply( packed_wt_separate_poses[ ii ]);
			min_mover_split->apply( packed_mut_separate_poses[ ii ]);

			//figure out an rmsd if the native structure is given
			core::Real  wt_rms (0.0);
			core::Real  mut_rms (0.0);
			if (basic::options::option[ in::file::native ].user()){
				// only calculate on complexes
				// allow superposition because RB min is allowed
				wt_rms = scoring::rmsd_with_super( native_pose, packed_wt_complex_poses[ ii ], scoring::is_protein_CA ) ;
				mut_rms = scoring::rmsd_with_super( native_pose, packed_mut_complex_poses[ ii ], scoring::is_protein_CA ) ;
			}

			//now dump the pdbs (all 4 of them for each run)
			//first set up naming system...
			std::string wt_complex_filename = pdb_file_name.base() + "_complex_wt_" + mut_string_output + "_" + ObjexxFCL::format::I( 3, 3, ii ) + ".pdb";
			std::string mut_complex_filename = pdb_file_name.base() + "_complex_mut_" + mut_string_output + "_" + ObjexxFCL::format::I( 3, 3, ii ) + ".pdb";
			std::string wt_split_filename = pdb_file_name.base() + "_separate_wt_" + mut_string_output+ "_" + ObjexxFCL::format::I( 3, 3, ii ) + ".pdb";
			std::string mut_split_filename = pdb_file_name.base() + "_separate_mut_" + mut_string_output+ "_" + ObjexxFCL::format::I( 3, 3, ii ) + ".pdb";

			packed_wt_complex_poses[ ii ].dump_scored_pdb( wt_complex_filename, *(scorefxn()) );
			packed_mut_complex_poses[ ii ].dump_scored_pdb( mut_complex_filename, *(scorefxn()) );
			packed_wt_separate_poses[ ii ].dump_scored_pdb( wt_split_filename, *(scorefxn()) );
			packed_mut_separate_poses[ ii ].dump_scored_pdb( mut_split_filename, *(scorefxn()) );

			TR << "Scores after Final Minimization round: "<< ii<<" \n"
				 << wt_complex_filename <<": "<<(*scorefxn)(packed_wt_complex_poses[ ii ]) <<  "  rms: "<< wt_rms<<"\n"
				 << wt_split_filename <<": " <<(*scorefxn)(packed_wt_separate_poses[ ii ]) <<  "\n"
				 << mut_complex_filename  <<": " <<(*scorefxn)(packed_mut_complex_poses[ ii ]) <<  "  rms: "<< mut_rms<<"\n"
				 << mut_split_filename << ": " <<(*scorefxn)(packed_mut_separate_poses[ ii ]) << std::endl;

		} //end for all packed poses
	}// end output all packed poses

	std::cout << "Complete..." << std::endl;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

} //end main


