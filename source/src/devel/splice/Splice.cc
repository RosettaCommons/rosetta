// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/splice/Splice.cc
/// @brief
/// @author Sarel Fleishman (sarel@weizmann.ac.il)
/// @modified by gideonla (glapidoth@gmail.com)

///Clarifications about Variables:
/*
 * Template_pose = The pose that is used as refrence for all the start and end positions of loop segments
 * Souce_pose = The
 * startn = the position on the pose where the loop starts (N-ter)
 * startc = The position on the pose where the loop ends (C-ter)
 * from_res = The user supplies loop start residues accoding to Template PDB. Splice.cc updates this to be the correct residue accrding to the current pose (both
 * 			  Structures are lianged)
 * to_res = Same as from_res but apllies to the C-ter of the segment
 *

 */


// Unit headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/Job.hh>
#include <devel/splice/Splice.hh>
#include <devel/splice/SpliceSegment.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <devel/splice/SpliceCreator.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/AA.hh>
#include <protocols/protein_interface_design/filters/TorsionFilter.hh>
#include <protocols/protein_interface_design/util.hh>
#include <boost/foreach.hpp>
#include <protocols/toolbox/task_operations/RestrictChainToRepackingOperation.hh>
#include <protocols/rigid/RB_geometry.hh>
#define foreach BOOST_FOREACH
// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/DataMapObj.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <protocols/protein_interface_design/movers/AddChainBreak.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <iostream>
#include <sstream>
#include <algorithm>
//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <protocols/toolbox/task_operations/ProteinInterfaceDesignOperation.hh>
#include <protocols/toolbox/task_operations/ThreadSequenceOperation.hh>
#include <protocols/toolbox/task_operations/SeqprofConsensusOperation.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <numeric/xyzVector.hh>
#include <protocols/loops/FoldTreeFromLoopsWrapper.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/protein_interface_design/movers/LoopLengthChange.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/CircularHarmonicFunc.hh>
#include <numeric/constants.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/scoring/Energies.hh>
#include <numeric/xyz.functions.hh>
//////////////////////////////////////////////////
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
#include <core/pose/util.hh>
///////////////////////////////////////////////////
#include <fstream>
#include <ctime>

namespace devel {
namespace splice {

static basic::Tracer TR( "devel.splice.Splice" );
static basic::Tracer TR_ccd( "devel.splice.Splice_ccd" );
		
static numeric::random::RandomGenerator RG( 78289 );
std::string
SpliceCreator::keyname() const
{
	return SpliceCreator::mover_name();
}

protocols::moves::MoverOP
SpliceCreator::create_mover() const {
	return new Splice;
}

std::string
SpliceCreator::mover_name()
{
	return "Splice";
}

Splice::Splice() :
			Mover( SpliceCreator::mover_name() ),
			from_res_( 0 ),
			to_res_( 0 ),
			saved_from_res_( 0 ),
			saved_to_res_( 0 ),
			source_pdb_( "" ),
			ccd_( true ),
			rms_cutoff_( 999999 ),
			res_move_( 4 ),
			randomize_cut_( false ),
			cut_secondarystruc_( false ),
			task_factory_( NULL ),
			design_task_factory_( NULL ),
			torsion_database_fname_( "" ),
			database_entry_( 0 ),
			database_pdb_entry_( "" ),
			template_file_( "" ),
			poly_ala_( true ),
			equal_length_( false ),
			template_pose_( NULL ),
			start_pose_( NULL ),
			saved_fold_tree_( NULL ),
			design_( false ),
			dbase_iterate_( false ),
			first_pass_( true ),
			locked_res_( NULL ),
			locked_res_id_( ' ' ),
			checkpointing_file_ ( "" ),
			loop_dbase_file_name_( "" ),
			loop_pdb_source_( "" ),
			mover_tag_( NULL ),
			splice_filter_( NULL ),
			use_sequence_profiles_( false ),
			segment_type_( "" ),
			profile_weight_away_from_interface_( 1.0 ),
			restrict_to_repacking_chain2_( true ),
			design_shell_(6.0),
			repack_shell_(8.0),
			scorefxn_(NULL),
			Pdb4LetName_("")
{
	add_sequence_constraints_only_ = false;
	torsion_database_.clear();
	delta_lengths_.clear();
	dbase_subset_.clear();
	splice_segments_.clear();
	pdb_segments_.clear();
	end_dbase_subset_ = new protocols::moves::DataMapObj< bool >;
	end_dbase_subset_->obj = false;
	basic::options::option[ basic::options::OptionKeys::out::file::pdb_comments ].value(true);
}

//Tell the JD to output comments



utility::vector1< std::string > segment_names_ordered_;//This vector will hold the segment names by order so when the segments are concatented into a single profile it is done by user defined order
std::string dofs_pdb_name; //This variable hold the name of the pdb in the torsion db

Splice::~Splice() {}

/// @brief copy a stretch of aligned phi-psi dofs from source to target. No repacking no nothing.
/// The core function, copy_segment, copies residues from the source to the target without aligning the residues, thereby delivering all of their dofs
void
copy_stretch( core::pose::Pose & target, core::pose::Pose const & source, core::Size const from_res, core::Size const to_res ){
	using namespace core::pose;
	using namespace protocols::rosetta_scripts;
	using namespace core::chemical;

	core::Size const host_chain( 1 ); /// in certain cases, when the partner protein sterically overlaps with the designed protein, there are amibguities about which chain to search. The use of host_chain removes these ambiguities. Here, ugly hardwired
	core::Size const from_nearest_on_source( find_nearest_res( source, target, from_res, host_chain ) );
	core::Size const to_nearest_on_source( find_nearest_res( source, target, to_res, host_chain ) );
	TR<<"target: "<<from_res<<" "<<to_res<<" source: "<<from_nearest_on_source<<" "<<to_nearest_on_source<<std::endl;
	runtime_assert( from_nearest_on_source && to_nearest_on_source );
	// change loop length:
	core::Size const residue_diff( to_nearest_on_source - from_nearest_on_source - (to_res - from_res ));
	//	if( residue_diff == 0 ){
	//		TR<<"skipping copy_stretch since loop lengths are identical"<<std::endl;
	//		return;
	//	}
	core::kinematics::FoldTree const saved_ft( target.fold_tree() );
	TR<<"DEBUG: copy_stretch foldtree: "<<saved_ft<<std::endl;
	protocols::protein_interface_design::movers::LoopLengthChange llc;
	llc.loop_start( from_res );
	llc.loop_end( to_res );
	llc.delta( residue_diff );
	//	target.dump_pdb( "before_copy_stretch_llc_test.pdb" );
	llc.apply( target );
	//	target.dump_pdb( "after_copy_stretch_llc_test.pdb" );

	target.copy_segment( to_nearest_on_source - from_nearest_on_source + 1, source, from_res, from_nearest_on_source );
	//target.dump_pdb( "after_copy_stretch_test.pdb" );
}

/// The checkpointing file has the following structure: the first line contains an ordered list of the dbase_subset_ for splice to iterate over the loop database. The second line contains the last element tested (the loop-entry number in the database; not the iterator to it!) and the third line contains the best element tested (again, the loop number from the database, not the iterator!).
/// To recover from a checkpoint the following reads the dbase_subset_ then, if this is a first_pass_ the best entry becomes current, and if it is not a first_pass then the current entry is current.
void
Splice::load_from_checkpoint()
{
	using namespace std;

	if( checkpointing_file_ == "" ) return;
	utility::io::izstream data( checkpointing_file_ );
	if ( !data ) return;
	TR<<"Loading from checkpoint"<<std::endl;
	/// first read the dbase_subset from the checkpointing file
	{
		string line;
		getline( data, line );
		if( line.length() == 0 ){
			TR<<"Checkpointing file empty or corrupted. Not loading."<<std::endl;
			return;
		}
		istringstream line_stream( line );
		dbase_subset_.clear();
		while( !line_stream.eof() ){
			core::Size entry;
			line_stream >> entry;
			dbase_subset_.push_back( entry );
		}
	}
	TR<<"dbase subset order loaded from checkpoint is: ";
	foreach( core::Size const i, dbase_subset_ ){
		TR<<i<<' ';
	}

	{
		std::string line;
		getline( data, line );
		istringstream line_stream( line );
		core::Size entry;
		line_stream >> entry;
		current_dbase_entry_ = std::find( dbase_subset_.begin(), dbase_subset_.end(), entry );
	}
	TR << "current dbase entry loaded from checkpoint is: " << *current_dbase_entry_ << std::endl;
}

void
Splice::save_to_checkpoint() const{
	if( checkpointing_file_ == "" )
		return;
	TR<<"Splice checkpointing to file: "<<checkpointing_file_<<std::endl;
	std::ofstream data;
	data.open( checkpointing_file_.c_str(), std::ios::out );
	if( !data.good() )
		utility_exit_with_message( "Unable to open splice checkpointing file for writing: " + checkpointing_file_ + "\n" );
	foreach( core::Size const dbase_entry, dbase_subset_ ){
		TR<<' '<<dbase_entry;
		data << ' ' << dbase_entry;
	}
	if( current_dbase_entry_ == dbase_subset_.end() )
		data<<'\n'<<99999<<std::endl;
	else
		data<<'\n'<<*current_dbase_entry_<<std::endl;
	data.close();
}

///@brief controls which dbase entry will be used. Three options: 1. specific one according to user instruction; 2. randomized out of a subset of the dbase with fitting sequence lengths (if user specified 0); 3. iterating over dbase subset
core::Size
Splice::find_dbase_entry( core::pose::Pose const & pose )
{
	core::Size dbase_entry( database_entry() );
	if( first_pass_ ){/// setup the dbase subset where loop lengths fit the selection criteria
		for( core::Size i = 1; i <= torsion_database_.size(); ++i ){// find entries that fit the length criteria
			using namespace protocols::rosetta_scripts;

			ResidueBBDofs const & dofs( torsion_database_[ i ] );
			core::Size const nearest_to_entry_start_on_pose( find_nearest_res( pose, *template_pose_, dofs.start_loop(), 1/*chain*/ ) );
			core::Size const nearest_to_entry_stop_on_pose( find_nearest_res( pose, *template_pose_, dofs.stop_loop(), 1/*chain*/ ) );
			core::Size const pose_residues = nearest_to_entry_stop_on_pose - nearest_to_entry_start_on_pose + 1;
			int const delta( dofs.size() - pose_residues );
			if( locked_res() >= nearest_to_entry_start_on_pose && locked_res() <= nearest_to_entry_stop_on_pose ){
				/// if locked_res is within the loop, don't select different loop lengths
				if( delta != 0 )
					continue;
			}
			bool const fit = std::find( delta_lengths_.begin(), delta_lengths_.end(), delta ) != delta_lengths_.end();
			if( fit || database_pdb_entry_ != "" || dbase_entry != 0 )
				dbase_subset_.push_back( i );
		}
		if( dbase_subset_.empty() ){
			TR<<"Loop of appropriate length not found in database. Returning"<<std::endl;
			retrieve_values();
			return 0;
		}
		TR<<"Found "<<dbase_subset_.size()<<" entries in the torsion dbase that match the length criteria"<<std::endl;
		numeric::random::random_permutation( dbase_subset_.begin(), dbase_subset_.end(), RG );
		current_dbase_entry_ = dbase_subset_.begin();
		load_from_checkpoint();
		first_pass_ = false;
	} // fi first_pass
	if( dbase_iterate() ){
		load_from_checkpoint();
		if( current_dbase_entry_ == dbase_end() ){
			TR<<"Request to read past end of dbase. Splice returns without doing anything."<<std::endl;
			return 0;
		}
		dbase_entry = *current_dbase_entry_;
		if( !first_pass_ )
			current_dbase_entry_++;
		if( current_dbase_entry_ == dbase_end() ){
			TR<<"Reached last dbase entry"<<std::endl;
			end_dbase_subset_->obj = true;
		}
	} // fi dbase_iterate
	else if( dbase_entry == 0 ){
		if( database_pdb_entry_ == "" )//randomize dbase entry
			dbase_entry = ( core::Size )( RG.uniform() * dbase_subset_.size() + 1 );
		else{ // look for the pdb_entry name
			for( core::Size count = 1; count <= dbase_subset_.size(); ++count ){
				if( torsion_database_[ dbase_subset_[ count ] ].source_pdb() == database_pdb_entry_ ){
					TR<<"Found entry for "<<database_pdb_entry_<<" at number "<<dbase_subset_[ count ]<<std::endl;
					dbase_entry = dbase_subset_[ count ];
					break;
				}
			}
			runtime_assert( dbase_entry <= dbase_subset_.size() );
		}
	}//fi dbase_entry==0

	return dbase_entry;
}

void
Splice::apply( core::pose::Pose & pose )
{
	using namespace protocols::rosetta_scripts;
	using core::chemical::DISULFIDE;

	if( add_sequence_constraints_only() ){
		TR<<"Only adding sequence constraints!!! Not doing any splice!!!"<<std::endl;
		add_sequence_constraints( pose );
		return;
	}

	set_last_move_status( protocols::moves::MS_SUCCESS );
	TR<<"Starting splice apply"<<std::endl;
	save_values();
	if( locked_res() ){
		locked_res_id( pose.residue( locked_res() ).name1() );
		TR<<"locked residue/locked_residue_id set to: "<<locked_res()<<','<<locked_res_id()<<std::endl;
	}

	/// from_res() and to_res() can be determined directly on the tag, through a taskfactory, or through a template file. If through a template file,
	/// we start by translating from_res/to_res from the template file to the in coming pose as in the following paragraph
	if( template_file_ != "" ){ /// using a template file to determine from_res() to_res()
		core::Size template_from_res( 0 ), template_to_res( 0 );

		if( from_res() && to_res() ){
			template_from_res = find_nearest_res( pose, *template_pose_, from_res(), 1/*chain*/);
			template_to_res   = find_nearest_res( pose, *template_pose_, to_res(), 1/*chain*/  );
			runtime_assert( template_from_res );
			runtime_assert( template_to_res );
		}

		from_res( template_from_res );
		to_res( template_to_res );
	}// fi template_file != ""

	core::pose::Pose const in_pose_copy( pose );
	pose.conformation().detect_disulfides(); // just in case; but I think it's unnecessary

	/// from_res/to_res can also be determined through task factory, by identifying the first and last residues that are allowed to design in this tf
	if( torsion_database_fname_ == "" && from_res() == 0 && to_res() == 0 ){/// set the splice site dynamically according to the task factory
		utility::vector1< core::Size > designable( protocols::rosetta_scripts::residue_packer_states( pose, task_factory(), true/*designable*/, false/*packable*/ ) );
		std::sort( designable.begin(), designable.end() );
		from_res( designable[ 1 ] );
		to_res( designable[ designable.size() ] );
	}
	core::pose::Pose source_pose;
	core::Size nearest_to_from( 0 ), nearest_to_to( 0 ), residue_diff( 0 ); // residues on source_pose that are nearest to from_res and to_res; what is the difference in residue numbers between incoming pose and source pose
	ResidueBBDofs dofs; /// used to store the torsion/resid dofs from any of the input files
	dofs.clear();
	core::Size cut_site( 0 );
	if( torsion_database_fname_ == "" ){ // read dofs from source pose rather than database
		core::import_pose::pose_from_pdb( source_pose, source_pdb_ );
		nearest_to_from = find_nearest_res( source_pose, pose, from_res(), 1/*chain*/ );
		nearest_to_to = find_nearest_res( source_pose, pose, to_res(), 1/*chain*/ );
		residue_diff = nearest_to_to - nearest_to_from - ( to_res() - from_res() );
		if( nearest_to_from == 0 || nearest_to_to == 0 ){
			std::ostringstream os;
			os<<"nearest_to_from: "<<nearest_to_from<<" nearest_to_to: "<<nearest_to_to<<". Failing"<<std::endl;
			utility_exit_with_message(os.str());
		/*	TR<<"nearest_to_from: "<<nearest_to_from<<" nearest_to_to: "<<nearest_to_to<<". Failing"<<std::endl;
			set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
			retrieve_values();
			return;
		*/

		}
		for( core::Size i = nearest_to_from; i <= nearest_to_to; ++i ){
			if( source_pose.residue( i ).has_variant_type( DISULFIDE ) ){/// in future, using disulfides would be a great boon as it rigidifies loops.
				std::ostringstream os;
				os<<"Residue "<<i<<" is a disulfide. Failing"<<std::endl;
				utility_exit_with_message(os.str());
				/*TR<<"Residue "<<i<<" is a disulfide. Failing"<<std::endl;
				set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
				retrieve_values();
				return;
				*/

			}
			/// Feed the source_pose dofs into the BBDofs array
			BBDofs residue_dofs;
			residue_dofs.resid( i ); /// resid is probably never used
			residue_dofs.phi( source_pose.phi( i ) );
			residue_dofs.psi( source_pose.psi( i ) );
			residue_dofs.omega( source_pose.omega( i ) );

			//core::Size const nearest_on_target( find_nearest_res( pose, source_pose, i ) );

			/// convert 3let residue code to 1let code
			std::stringstream ss; std::string s;
			ss << source_pose.residue( i ).name1();
			ss >> s;
			residue_dofs.resn( s );

			dofs.push_back( residue_dofs );
		}// for i nearest_to_from..nearest_to_to
		cut_site = dofs.cut_site() ? dofs.cut_site() + from_res() - 1: to_res(); // isn't this always going to be to_res()? I think so...
	}// fi torsion_database_fname==NULL
	else{/// read from dbase
		core::Size const dbase_entry( find_dbase_entry( pose ) );
		if( dbase_entry == 0 )// failed to read entry
			return;
		dofs = torsion_database_[ dbase_entry ];
		std::string const source_pdb_name( dofs.source_pdb() );
		dofs_pdb_name = source_pdb_name;
		if( use_sequence_profiles_ ){
			load_pdb_segments_from_pose_comments( pose );
			modify_pdb_segments_with_current_segment( source_pdb_name );
		}
		TR<<"Taking loop from source pdb "<<source_pdb_name<<std::endl;
		if( mover_tag_() != NULL )
			mover_tag_->obj = "segment_" + source_pdb_name;
		foreach( BBDofs & resdofs, dofs ){/// transform 3-letter code to 1-letter code
			using namespace core::chemical;
			if( resdofs.resn() == "CYD" ){// at one point it would be a good idea to use disfulfides rather than bail out on them...; I think disulfided cysteins wouldn't be written as CYD. This requires something more clever...
				TR<<"Residue "<<resdofs.resid()<<" is a disulfide. Failing"<<std::endl;
				set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
				retrieve_values();
				return;
			}
			std::stringstream ss; std::string s;
			ss << oneletter_code_from_aa( aa_from_name( resdofs.resn() ) );
			ss >> s;
			resdofs.resn( s );
		}/// foreach resdof
		nearest_to_to = dofs.size(); /// nearest_to_to and nearest_to_from are used below to compute the difference in residue numbers...
		nearest_to_from = 1;
		/// set from_res/to_res/cut_site on the incoming pose
		if( template_file_ != "" ){/// according to the template pose
			from_res( find_nearest_res( pose, *template_pose_, dofs.start_loop(), 1/*chain*/ ) );
			to_res( find_nearest_res( pose, *template_pose_, dofs.stop_loop(), 1/*chain*/ ) );
			//			to_res( from_res() + dofs.size() -1);
			runtime_assert( from_res() );
			runtime_assert( to_res() );
			cut_site = dofs.cut_site() - dofs.start_loop() + from_res();
		} // fi template_file != ""
		else{/// according to the dofs array (taken from the dbase)
			from_res( dofs.start_loop() );
			to_res( dofs.stop_loop() );
			cut_site = dofs.cut_site();
			runtime_assert( from_res() && to_res() && cut_site );
		}
		residue_diff = dofs.size() - ( dofs.stop_loop() - dofs.start_loop()  + 1 );
	}// read from dbase
	TR<<"From res: "<<from_res()<<" to_res: "<<to_res()<<std::endl;
	runtime_assert( to_res() > from_res() );
	//	if( saved_fold_tree_ )/// is saved_fold_tree_ being used?
	//		pose.fold_tree( *saved_fold_tree_ );

	/// The database is computed with respect to the template pose, so before applying dofs from the dbase it's important to make that stretch identical to
	/// the template. from_res() and to_res() were previously computed to be with respect to the incoming pose, so within this subroutine the refer to pose rather
	/// than template_pose (this is a bit confusing, but it works!)
	copy_stretch( pose, *template_pose_, from_res(), to_res() );
	//	( *scorefxn() ) ( pose );

	using namespace utility;
	/// randomize_cut() should not be invoked with a database entry, b/c the dbase already specified the cut sites.
	/// this is important b/c nearest_to_from/nearest_to_to are degenerate if the dbase is used.
	if( randomize_cut() ){
		/// choose cutsite randomly within loop residues on the loop (no 2ary structure)
		core::scoring::dssp::Dssp dssp( source_pose );
		dssp.dssp_reduced(); // switch to simplified H E L notation
		std::vector< core::Size > loop_positions_in_source;
		loop_positions_in_source.clear();
		TR<<"DSSP of source segment: ";
		for( core::Size i = nearest_to_from; i <= std::min( nearest_to_to, to_res() - from_res() + nearest_to_from ); ++i ){
			if( dssp.get_dssp_secstruct( i ) == 'L' || cut_secondarystruc() ) // allow site for cutting if it's either in a loop or if cutting secondary structure is allowed
				loop_positions_in_source.push_back( i );
			TR<<dssp.get_dssp_secstruct( i );
		}
		TR<<std::endl;
		//New test to see what is the sequence of the new loop
		TR.Debug<<"The sequence of the source loop is: ";
		for( core::Size i = nearest_to_from; i <= std::min( nearest_to_to, to_res() - from_res() + nearest_to_from ); ++i ){
			TR.Debug<<source_pose.residue( i ).name1()<< " ";
		}
		TR<<std::endl;
		cut_site = loop_positions_in_source[ (core::Size) ( RG.uniform() * loop_positions_in_source.size()) ] - nearest_to_from + from_res();
		TR<<"Cut placed at: "<<cut_site<<std::endl;
	}// fi randomize_cut
	//	pose.dump_pdb( "before_ft_test.pdb" ); //this is the strucutre before changing the loop
	fold_tree( pose, from_res(), pose.total_residue()/*to_res() SJF DEBUGGING 7Dec12*/, cut_site );/// the fold_tree routine will actually set the fold tree to surround the loop
	//	pose.dump_pdb( "after_ft_test.pdb" );
	/// change the loop length
	TR<<"Foldtree before loop length change: "<<pose.fold_tree()<<std::endl;
	protocols::protein_interface_design::movers::LoopLengthChange llc;
	llc.loop_start( from_res() );
	llc.loop_end( cut_site + residue_diff < from_res() ? to_res() : cut_site );
	llc.delta( residue_diff );
	llc.apply( pose );
	TR<<"Foldtree after loop length change: "<<pose.fold_tree()<<std::endl;

	//pose.dump_pdb( "after_2ndllc_test.pdb" );
	/// set torsions
	core::Size const total_residue_new( dofs.size() );
	TR<<"Changing dofs\n";
	for( core::Size i = 0; i < total_residue_new; ++i ){
		core::Size const pose_resi( from_res() + i );
		//		TR<<"Previous phi/psi/omega at resi: "<<pose_resi<<" "<<pose.phi( pose_resi )<<'/'<<pose.psi( pose_resi )<<'/'<<pose.omega( pose_resi )<<'\n';
		pose.set_phi( pose_resi, dofs[ i + 1 ].phi() );
		pose.set_psi(  pose_resi, dofs[ i + 1 ].psi() );
		pose.set_omega( pose_resi, dofs[ i + 1 ].omega() );
		//		pose.dump_pdb( "dump"+ utility::to_string( i ) + ".pdb" );
		TR<<"resi, phi/psi/omega: "<< pose_resi<<' '<<pose.phi( pose_resi )<<'/'<<pose.psi( pose_resi )<<'/'<<pose.omega( pose_resi )<<std::endl;
		//		TR<<"requested phi/psi/omega: "<<dofs[ i + 1 ].phi()<<'/'<<dofs[i+1].psi()<<'/'<<dofs[i+1].omega()<<std::endl;
	}
	//	pose.dump_pdb( "after_changedofs_test.pdb" );
	TR<<std::endl;
	std::string threaded_seq( "" );/// will be all ALA except for Pro/Gly on source pose and matching identities on source pose
	/// Now decide on residue identities: Alanine throughout except when the template pose has Gly, Pro or a residue that is the same as that in the original pose
	utility::vector1< core::Size > pro_gly_res; //keeping track of where pro/gly residues are placed
	pro_gly_res.clear();
	//pose.dump_pdb("before_threading.pdb");
	for( core::Size i = 0; i < total_residue_new; ++i ){
		core::Size const pose_resi( from_res() + i );
		std::string const dofs_resn( dofs[ i + 1 ].resn() );
		runtime_assert( dofs_resn.length() == 1 );
		if( pose_resi == locked_res() ){
			threaded_seq += locked_res_id();
			continue;
		}
		if( design() ){ // all non pro/gly residues in template are allowed to design
			if( dofs_resn == "G" || dofs_resn == "P" ){
				threaded_seq += dofs_resn;
				pro_gly_res.push_back( pose_resi );
				TR<<"Pro/Gly will be allowed at: "<<pose_resi<<std::endl;
			}
			else
				threaded_seq += 'x';
			//pose.replace_residue(pose_resi,source_pose.residue(nearest_to_from+i),0);
			continue;
		}
		core::Size const host_chain( 1 );
		core::Size const nearest_in_copy( find_nearest_res( in_pose_copy, pose, pose_resi, host_chain ) );
		if( ( nearest_in_copy > 0 && dofs_resn[ 0 ] == in_pose_copy.residue( nearest_in_copy ).name1() )  || dofs_resn == "G" || dofs_resn == "P" )
			threaded_seq += dofs_resn;
		else{
			if( poly_ala() )
				threaded_seq += "A";
			else{
				char orig_residue( 0 );
				if( nearest_in_copy )
					orig_residue = in_pose_copy.residue( nearest_in_copy ).name1();
				if( orig_residue == 0 || orig_residue == 'G' || orig_residue == 'P' )
					threaded_seq += 'x'; // residues that were originally Gly/Pro can be designed now
				else
					threaded_seq += ' '; // only repack
			}
		}
	}
	//pose.dump_pdb( "after_sequence_thread.pdb" );
	using namespace protocols::toolbox::task_operations;
	using namespace core::pack::task;
	ThreadSequenceOperationOP tso = new ThreadSequenceOperation;
	tso->target_sequence( threaded_seq );
	tso->start_res( from_res() );
	tso->allow_design_around( true ); // 21Sep12: from now on the design shell is determined downstream //false );
	TR<<"Threading sequence: "<<threaded_seq<<" starting from "<<from_res()<<std::endl;
	TaskFactoryOP tf;
	if( design_task_factory()() == NULL )
		tf = new TaskFactory;
	else
		tf = new TaskFactory( *design_task_factory() );

	if( restrict_to_repacking_chain2() ){
		for( core::Size i = 2; i <= pose.conformation().num_chains(); ++i ){
			TR<<"Restricting chain "<<i<<" to repacking only"<<std::endl;
			tf->push_back( new protocols::toolbox::task_operations::RestrictChainToRepackingOperation( i ) );
		}
	}

	tf->push_back( new operation::InitializeFromCommandline );
	tf->push_back( new operation::NoRepackDisulfides );
	tf->push_back( tso );
	DesignAroundOperationOP dao = new DesignAroundOperation;
	dao->design_shell( (design_task_factory()() == NULL ? 0.0 : design_shell()) ); // threaded sequence operation needs to design, and will restrict design to the loop, unless design_task_factory is defined, in which case a larger shell can be defined
	dao->repack_shell( repack_shell() );
	for( core::Size i = from_res(); i <= from_res() + total_residue_new - 1; ++i ){
		if( !pose.residue( i ).has_variant_type( DISULFIDE ) )
			dao->include_residue( i );
	}
	tf->push_back( dao );
	TR<<"allowing pro/gly only at positions (29Mar13, given sequence profiles, now allowing pro/gly/his at all designed positions. The following is kept for benchmarking): ";
	for(core::Size res_num=1; res_num <= pose.total_residue(); res_num++ ){
		if( std::find( pro_gly_res.begin(), pro_gly_res.end(), res_num ) == pro_gly_res.end() ){
			operation::RestrictAbsentCanonicalAASOP racaas = new operation::RestrictAbsentCanonicalAAS;
			racaas->keep_aas( "ADEFGHIKLMNPQRSTVWY" ); /// disallow pro/gly/cys/his /// 29Mar13 now allowing all residues other than Cys. Expecting sequence profiles to take care of gly/pro/his
			racaas->include_residue( res_num );
			tf->push_back( racaas);
		}
		else
			TR<<res_num<<", ";
	}
	TR<<std::endl;
	//	if( locked_res() ){
	//		operation::PreventRepackingOP pr = new operation::PreventRepacking;
	//		pr->include_residue( locked_res() );
	//		tf->push_back( pr );
	//		TR<<"preventing locked residue "<<locked_res()<<" from repacking"<<std::endl;
	//	}


	protocols::protein_interface_design::movers::AddChainBreak acb;
	acb.resnum( utility::to_string( cut_site + residue_diff ) );
	acb.find_automatically( false );
	acb.change_foldtree( false );
	acb.apply( pose );
	TR<<"Adding chainbreak at: "<<cut_site + residue_diff <<std::endl;
	//SJF debug	pose.conformation().detect_disulfides();
	//	( *scorefxn() ) ( pose );
	//	pose.update_residue_neighbors();
	if( use_sequence_profiles_ )
		TR<<"NOW ADDING SEQUENCE CONSTRAINTS"<<std::endl;
	add_sequence_constraints( pose );

	if( ccd() ){
		using namespace protocols::loops;

		core::Size const startn(from_res());
		core::Size const startc(from_res() + total_residue_new - 1 );

///		Loop loop( std::max( (core::Size) 2, from_res() - 6 )/*start*/, std::min( pose.total_residue()-1, to_res() + 6 )/*stop*/, cut_site/*cut*/ );
		Loop loop( startn, startc, cut_site ); /// Gideon & Sarel (8Jul13): we're now respecting the user's choice of from_res to_res and not melting the framework
		LoopsOP loops = new Loops();
		loops->push_back( loop );

		/// Gideon & Sarel (8Jul13): the comment below is no longer true, see comments above.
		/// Set ccd to minimize 4 residues at each loop terminus including the first residue of the loop. This way,
		/// the torsion in the loop are maintained. Allow repacking around the loop.
		/// If disulfide occurs in the range that is allowed to minimize, adjust that region to not include disulf
		core::scoring::ScoreFunctionOP scorefxn_local( scorefxn()->clone() );/// in case you want to modify the scorefxn. Currently not used
		protocols::loops::loop_mover::refine::LoopMover_Refine_CCD ccd_mover( loops, scorefxn_local );
		ccd_mover.temp_initial( 1.5 );
		ccd_mover.temp_final( 0.5 );
		core::kinematics::MoveMapOP mm;
		mm = new core::kinematics::MoveMap;
		mm->set_chi( false ); mm->set_bb( false ); mm->set_jump( false );
		mm->set_jump( 2, true ); /// 1Feb13 for cases in which the we're splicing in the presence of a ligand
		/// First look for disulfides. Those should never be moved.
///8Jul13: Gideon & Sarel: we previously melted the framework down from the from_res to_res points. Now, we respect the
/// user's choice of from_res to_res, but the user has to make sure that the span doesn't contain a disulfide.
/// a disulfide within the span will result in an assertion
/*		core::Size disulfn( 0 ), disulfc( 0 );
		for( core::Size i = from_res() - 3; i <= from_res(); ++i ){
			if( pose.residue( i ).has_variant_type( DISULFIDE ) ){
				disulfn = i;
			}
		}
		for( core::Size i = from_res() + total_residue_new - 1; i <= from_res() + total_residue_new + 2; ++i ){
			if( pose.residue( i ).has_variant_type( DISULFIDE ) ){
				disulfc = i;
				break;
			}
		}*/
		//TR<<"startn="<<startn<<std::endl;
		//TR<<"Residues allowed to minimize"<<std::endl;
		for( core::Size i = startn; i <= std::min( startn + res_move()-1, startc ); ++i ){
		//TR<<i<<"+";
			mm->set_chi( i, true );
			mm->set_bb( i, true );

		}
		//TR<<""<<std::endl;;
			//	TR<<"Residues allowed to minimize"<<std::endl;
		//TR<<"starcn="<<startc<<std::endl;

		for( core::Size i = std::max( startn, startc - res_move() + 1); i <= startc; ++i ){
			//Debuging stuff:
		//TR<<i<<"+";
			mm->set_chi( i, true ); //allowing chi angle movement
			mm->set_bb( i, true ); //allowing bb movement
		}
		//TR<<""<<std::endl;;

		ccd_mover.set_task_factory( tf );
		ccd_mover.move_map( mm );


		//pose.dump_pdb("before_ccd.pdb");
		add_dihedral_constraints(pose, source_pose,nearest_to_from,nearest_to_to,cut_site  );//add dihedral constraints loop
		add_coordinate_constraints(pose, source_pose,nearest_to_from,nearest_to_to);//add coordiante constraints to loop

		//as control I want to see which residues are allowed to be design according to design shell.Gideonla may13
		utility::vector1< core::Size > designable_residues = residue_packer_states (pose,tf, true, false);
		TR<<"Residues Allowed to Design:"<<std::endl;
		for (utility::vector1<core::Size>::const_iterator i (designable_residues.begin()); i != designable_residues.end(); ++i) {
			TR<<pose.residue(*i).name1()<<*i<<",";
		}
		TR<<std::endl;

		utility::vector1< core::Size > packable_residues = residue_packer_states (pose,tf, false, true);
		TR<<"Residues Allowed to Repack:"<<std::endl;
		for (utility::vector1<core::Size>::const_iterator i (packable_residues.begin()); i != packable_residues.end(); ++i) {
			TR<<pose.residue(*i).name1()<<*i<<",";
		}
		TR<<std::endl;
		TR<<"Weighted score function before ccd:"<<std::endl;
		scorefxn()->show(pose);//before ccd starts make sure we have all constratins in place, gidoenla Jul13
		ccd_mover.apply( pose );
		TR<<"Weighted score function after ccd:"<<std::endl;
		scorefxn()->show(pose);//before ccd starts make sure we have all constratins in place, gidoenla Jul13
		//pose.dump_pdb("after_ccd.pdb");
		/// following ccd, compute rmsd to source loop to ensure that you haven't moved too much. This is a pretty decent filter
		if( torsion_database_fname_ == "" ){ // no use computing rms if coming from a database (no coordinates)
			core::Real rms( 0 );
			for( core::Size i = 0; i <= total_residue_new - 1; ++i ){
				core::Real const dist( pose.residue( from_res() + i ).xyz( "CA" ).distance( source_pose.residue( nearest_to_from+ i ).xyz("CA" ) ) );
				rms += dist;
			}
			core::Real const average_rms( rms / total_residue_new );
			TR<<"Average distance of spliced segment to original: "<< average_rms<<std::endl;
			if( average_rms >= rms_cutoff() ){
				TR<<"Failing because rmsd = "<<average_rms<<std::endl;
				set_last_move_status( protocols::moves::FAIL_RETRY );
				retrieve_values();
				return;
			}

			//print rms value to output pdb structure
				std::string Result;
				std::string Result_filter;        
				std::ostringstream convert; 
				std::ostringstream convert_filter; 
				convert << average_rms;     
				Result = convert.str(); 

				core::pose::add_comment(pose,"RMSD to source loop",Result);//change correct association between current loop and pdb file
				
				//print ChainBreak value to output pdb structure
				convert_filter << splice_filter()->score( pose );
				Result_filter = convert_filter.str(); 
				core::pose::add_comment(pose,"Chainbreak Val:",Result_filter);
			if( !splice_filter()->apply( pose ) ){
				TR<<"Failing because filter fails"<<std::endl;
				set_last_move_status( protocols::moves::FAIL_RETRY );
				retrieve_values();
				return;
			}
		}
		/// tell us what the torsions of the new (closed) loop are. This is used for dbase construction. At one point, might be a good idea to make the mover
		/// output the dofs directly to a dbase file rather than to a log file.
		TaskFactoryOP tf_dofs = new TaskFactory;
		DesignAroundOperationOP dao_dofs = new DesignAroundOperation;
		for( core::Size i = startn; i <= std::min(startc + res_move() - 1,startc); ++i )
			dao_dofs->include_residue( i );
		dao_dofs->design_shell( 0 );/// only include the loop residues
		tf_dofs->push_back( dao_dofs );
		protocols::protein_interface_design::filters::Torsion torsion;
		torsion.task_factory( tf_dofs );
		torsion.task_factory_set( true );
		torsion.apply( pose );
		core::Size const stop_on_template( to_res());
		TR_ccd << "start, stop, cut: "<<startn<<" "<<stop_on_template<<" "<<cut_site<<std::endl; /// used for the dbase

		/// Now write to dbase disk file
		if( loop_dbase_file_name_ != "" ){
			std::ofstream dbase_file;
			dbase_file.open( loop_dbase_file_name_.c_str(), std::ios::app );
			for( core::Size i = startn; i <= std::min( startc + res_move() - 1, startc ); ++i )
				dbase_file << pose.phi( i )<<' '<<pose.psi( i )<<' '<<pose.omega( i )<<' '<<pose.residue( i ).name3()<<' ';
			dbase_file << startn<<' '<<stop_on_template<<' '<<cut_site<<' ';
			if( loop_pdb_source_ != "" )
				dbase_file <<Pdb4LetName_<<std::endl;
			else
				dbase_file << "cut" << std::endl; // the word cut is used as a placeholder. It is advised to use instead the source pdb file in this field so as to keep track of the origin of dbase loops
			dbase_file.close();
		}
	}// fi ccd
	else{ // if no ccd, still need to thread sequence
		//Debugging, remove after, gideonla aug13
		//TR<<"NOT DOING CCD, DOING REPACKING INSTEAD"<<std::endl;
		PackerTaskOP ptask = tf()->create_task_and_apply_taskoperations( pose );
		protocols::simple_moves::PackRotamersMover prm( scorefxn(), ptask );
		//		pose.conformation().detect_disulfides();
		//		pose.update_residue_neighbors();
		//		(*scorefxn())(pose);
		prm.apply( pose );
		//pose.dump_pdb("before_rtmin.pdb");
		//After Re-packing we add RotamerTrialMover to resolve any left over clashes, gideonla Aug13
		TaskFactoryOP tf_rtmin  = new TaskFactory(*tf);//this taskfactory (tf_rttmin) is only used here. I don't want to affect other places in splice, gideonla aug13
		tf_rtmin->push_back( new operation::RestrictToRepacking()); //W don't rtmin to do design
		ptask = tf_rtmin()->create_task_and_apply_taskoperations( pose );
		protocols::simple_moves::RotamerTrialsMinMover rtmin( scorefxn(), *ptask );
		rtmin.apply(pose);
		//pose.dump_pdb("after_rtmin.pdb");
		
	}
	saved_fold_tree_ = new core::kinematics::FoldTree( pose.fold_tree() );
	retrieve_values();
}

/// splice apply might change the from_res/to_res internals since they sometimes refer to the template file. If that happens, we want the values to
/// revert to their original values before the end of the apply function (so retrieve_values) below must be called before return.
void
Splice::save_values(){
	saved_from_res_ = from_res();
	saved_to_res_ = to_res();
}

void
Splice::retrieve_values(){
	from_res( saved_from_res_ );
	to_res( saved_to_res_ );
	first_pass_ = false;
	save_to_checkpoint();
}

std::string
Splice::get_name() const {
	return SpliceCreator::mover_name();
}

void
Splice::parse_my_tag( TagPtr const tag, protocols::moves::DataMap &data, protocols::filters::Filters_map const & filters, protocols::moves::Movers_map const &, core::pose::Pose const & pose )
{
	utility::vector1< TagPtr > const sub_tags( tag->getTags() );

	typedef utility::vector1< std::string > StringVec;
	///Adding Checker. If the "Current Segment" does not appear in the list of segements than exit with error
	segment_names_ordered_.clear(); //This string vector hold all the segment names inserted bythe user to ensure that the sequence profile is built according to tthe user
	bool check_segment = true;
	foreach( TagPtr const sub_tag, sub_tags ){
		if( sub_tag->getName() == "Segments" ){
			check_segment = false;//This is set to false unless current segment appears in the segment list
			use_sequence_profiles_ = true;
			profile_weight_away_from_interface( tag->getOption< core::Real >( "profile_weight_away_from_interface", 1.0 ) );
			segment_type_ = sub_tag->getOption< std::string >( "current_segment" );
			TR<<"reading segments in splice "<<tag->getName()<<std::endl;
			/* e.g.,
<Splice name=splice_L2...
  <Segments current_segment=L1>
				<L1 pdb_profile_match="pdb_profile_match.L1" profiles="L1:L1.pssm"/>
				<L2 pdb_profile_match="pdb_profile_match.L2" profiles="L2:L2.pssm"/>
				<L3 pdb_profile_match="pdb_profile_match.L3" profiles="L3:L3.pssm"/>
				<Frm1 pdb_profile_match="pdb_profile_match.Frm1" profiles="Frm1:frm1.pssm"/>
				<Frm2 pdb_profile_match="pdb_profile_match.Frm2" profiles="Frm2:frm2.pssm"/>
				<Frm3 pdb_profile_match="pdb_profile_match.Frm3" profiles="Frm3:frm3.pssm"/>
                <Frm4 pdb_profile_match="pdb_profile_match.Frm4" profiles="Frm4:frm4.pssm"/>
			</Segments>
</Splice>
			 */
			utility::vector1< TagPtr > const segment_tags( sub_tag->getTags() );
			foreach( TagPtr const segment_tag, segment_tags ){
				std::string const segment_name( segment_tag->getName() );//get name of segment from xml
				std::string const pdb_profile_match( segment_tag->getOption< std::string >( "pdb_profile_match" ) ); // get name of pdb profile match, this file contains all the matching between pdb name and sub segment name, i.e L1.1,L1.2 etc
				std::string const profiles_str( segment_tag->getOption< std::string >( "profiles" ) );
				StringVec const profile_name_pairs( utility::string_split( profiles_str, ',' ) );
				SpliceSegmentOP splice_segment( new SpliceSegment );
				//TR<<"Now working on segment:"<<segment_name<<"\n";
				foreach( std::string const s, profile_name_pairs ){
					StringVec const profile_name_file_name( utility::string_split( s, ':' ) );
					//TR<<"pssm file:"<<profile_name_file_name[ 2 ]<<",segment name:"<<profile_name_file_name[ 1 ]<<std::endl;

					splice_segment->read_profile( profile_name_file_name[ 2 ], profile_name_file_name[ 1 ] );

				}
				splice_segment->read_pdb_profile( pdb_profile_match );

				//TR<<"the segment name is: "<<segment_name<<std::endl;
				if (segment_name.compare(segment_type_) == 0){
						check_segment=true;
				}
				splice_segments_.insert( std::pair< std::string, SpliceSegmentOP >( segment_name, splice_segment ) );
				segment_names_ordered_.push_back(segment_name);
			}//foreach segment_tag
		}// fi Segments
	}//foreach sub_tag
	if (!check_segment){
			utility_exit_with_message( "Segment "+segment_type_+" was not found in the list of segemnts. Check XML file\n");
	}
	 scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	add_sequence_constraints_only( tag->getOption< bool >( "add_sequence_constraints_only", false ) );
	if( add_sequence_constraints_only() ){
		TR<<"add_sequence_constraints only set to true. Therefore not parsing any of the other Splice flags."<<std::endl;
		return;
	}

	start_pose_ = new core::pose::Pose( pose );
	runtime_assert( tag->hasOption( "task_operations" ) != (tag->hasOption( "from_res" ) || tag->hasOption( "to_res" ) ) || tag->hasOption( "torsion_database" ) ); // it makes no sense to activate both taskoperations and from_res/to_res.
	runtime_assert( tag->hasOption( "torsion_database" ) != tag->hasOption( "source_pdb" ) );
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	if( !tag->hasOption( "task_operations" ) ){
		from_res( core::pose::parse_resnum( tag->getOption< std::string >( "from_res", "0" ), pose ) );
		to_res( core::pose::parse_resnum( tag->getOption< std::string >( "to_res", "0" ), pose ) );
	}
	if( tag->hasOption( "design_task_operations" ) ){
		TR<<"Defined design_task_factory, which will be used during splice design"<<std::endl;
		design_task_factory( protocols::rosetta_scripts::parse_task_operations( tag->getOption< std::string >( "design_task_operations" ), data ) );
	}
	if( tag->hasOption( "residue_numbers_setter" ) ){
		runtime_assert( !tag->hasOption( "locked_res" ) );
		locked_res_ = protocols::moves::get_set_from_datamap< protocols::moves::DataMapObj< utility::vector1< core::Size > > >( "residue_numbers", tag->getOption< std::string >( "residue_numbers_setter" ), data );
	}
	if( tag->hasOption( "torsion_database" ) ){
		torsion_database_fname( tag->getOption< std::string >( "torsion_database" ) );
		database_entry( tag->getOption< core::Size >( "database_entry", 0 ) );
		database_pdb_entry( tag->getOption< std::string >( "database_pdb_entry", "" ) );
		runtime_assert( !( tag->hasOption( "database_entry" ) && tag->hasOption( "database_pdb_entry" ) ) );
		runtime_assert( !( tag->hasOption( "delta_lengths" ) && (tag->hasOption( "database_pdb_entry" ) || tag->hasOption( "database_entry" ) ) ) );
		read_torsion_database();
		TR<<"torsion_database: "<<torsion_database_fname()<<" ";
		if( database_entry() == 0 ){
			if( database_pdb_entry_ == "" )
				TR<<" database entry will be randomly picked at run time. ";
			else
				TR<<" picking database entry "<<database_pdb_entry()<<std::endl;
		}
		else{
			TR<<" database_entry: "<<database_entry()<<" ";
			runtime_assert( database_entry() <= torsion_database_.size() );
		}
	}
	else
		source_pdb( tag->getOption< std::string >( "source_pdb" ) );

	ccd( tag->getOption< bool >( "ccd", 1 ) );
	//dihedral_const(tag->getOption< core::Real >( "dihedral_const", 0 ) );//Added by gideonla Apr13, set here any real posiive value to impose dihedral constraints on loop
	//coor_const(tag->getOption< core::Real >( "coor_const", 0 ) );//Added by gideonla May13, set here any real to impose coordinate constraint on loop
	design_shell(tag->getOption< core::Real >( "design_shell", 6.0 ) );//Added by gideonla May13,
	repack_shell(tag->getOption< core::Real >( "repack_shell", 8.0 ) );//Added by gideonla May13,
	rms_cutoff( tag->getOption< core::Real >( "rms_cutoff", 999999 ) );
	runtime_assert( !(tag->hasOption( "torsion_database" ) && tag->hasOption( "rms_cutoff" )) ); // torsion database doesn't specify coordinates so no point in computing rms
	res_move( tag->getOption< core::Size >( "res_move", 4 ) );
	randomize_cut( tag->getOption< bool >( "randomize_cut", false ) );
	runtime_assert( ( tag->hasOption( "randomize_cut" ) && tag->hasOption( "source_pose" ) ) || !tag->hasOption( "source_pose" ) );
	cut_secondarystruc( tag->getOption< bool >( "cut_secondarystruc", false ) );
	//	runtime_assert( (tag->hasOption( "cut_secondarystruc ") && tag->hasOption( "randomize_cut" )) || !tag->hasOption( "cut_secondarystruc" ) );
	template_file( tag->getOption< std::string >( "template_file", "" ) );
	equal_length( tag->getOption< bool >( "equal_length", false ) );
	poly_ala( tag->getOption< bool >( "thread_ala", true ) );

	std::string delta;
	if( tag->hasOption( "delta_lengths" ) ){
		delta = tag->getOption< std::string >( "delta_lengths" );
		StringVec const lengths_keys( utility::string_split( delta, ',' ) );
		foreach( std::string const delta, lengths_keys ){
			if( delta == "" ) continue;
			int const delta_i( 1 * atoi( delta.c_str() ) );
			delta_lengths_.push_back( delta_i );
		}
	}
	else
		delta_lengths_.push_back( 0 );
	std::sort( delta_lengths_.begin(), delta_lengths_.end() );
	std::unique( delta_lengths_.begin(), delta_lengths_.end() );

	if( template_file_ != "" ){ /// using a template file to determine from_res() to_res()
		if( data.has( "poses", template_file_ ) ){
			template_pose_ = data.get< core::pose::Pose * >( "poses", template_file_ );
			TR<<"using template pdb from datamap"<<std::endl;
		}
		else if( tag->hasOption( "template_file" ) ){
			template_pose_ = new core::pose::Pose;
			core::import_pose::pose_from_pdb( *template_pose_, template_file_ );
			data.add( "poses", template_file_, template_pose_ );
			TR<<"loading template_pose from "<<template_file_<<std::endl;
		}
	}
	else
		template_pose_ = new core::pose::Pose( pose );

	design( tag->getOption< bool >( "design", false ) );
	dbase_iterate( tag->getOption< bool >( "dbase_iterate", false ) );
	if( dbase_iterate() ){ /// put the end_dbase_subset_ variable on the datamap for LoopOver & MC to be sensitive to it
		std::string const curr_mover_name( tag->getOption< std::string >( "name" ) );
		data.add( "stopping_condition", curr_mover_name, end_dbase_subset_ );
		TR<<"Placed stopping_condition "<<curr_mover_name<<" on the DataMap"<<std::endl;
	}
	if( tag->hasOption( "locked_residue" ) ){
		locked_res( core::pose::parse_resnum( tag->getOption< std::string >( "locked_residue" ), pose ) );
		locked_res_id( pose.residue( locked_res() ).name1() );
		TR<<"locking residue "<<locked_res()<<" of identity "<<locked_res_id()<<std::endl;
	}
	checkpointing_file( tag->getOption< std::string > ( "checkpointing_file", "" ) );
	loop_dbase_file_name( tag->getOption< std::string > ( "loop_dbase_file_name", "" ) );
	if( tag->hasOption( "splice_filter" ))
		splice_filter( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "splice_filter" ), filters ) );
	if( tag->hasOption( "mover_tag" ) )
		mover_tag_ = protocols::moves::get_set_from_datamap< protocols::moves::DataMapObj< std::string > >( "tags", tag->getOption< std::string >( "mover_tag" ), data );
	loop_pdb_source( tag->getOption< std::string >( "loop_pdb_source", "" ) );

	restrict_to_repacking_chain2( tag->getOption< bool >( "restrict_to_repacking_chain2", true ) );

	TR<<"from_res: "<<from_res()<<" to_res: "<<to_res()<<" dbase_iterate: "<<dbase_iterate()<<" randomize_cut: "<<randomize_cut()<<" cut_secondarystruc: "<<cut_secondarystruc()<<" source_pdb: "<<source_pdb()<<" ccd: "<<ccd()<<" rms_cutoff: "<<rms_cutoff()<<" res_move: "<<res_move()<<" template_file: "<<template_file()<<" checkpointing_file: "<<checkpointing_file_<<" loop_dbase_file_name: "<<loop_dbase_file_name_<<" loop_pdb_source: "<<loop_pdb_source()<<" mover_tag: "<<mover_tag_<<" torsion_database: "<<torsion_database_fname_<<" restrict_to_repacking_chain2: "<<restrict_to_repacking_chain2()<<std::endl;
}

protocols::moves::MoverOP
Splice::clone() const {
	return( protocols::moves::MoverOP( new Splice( *this ) ));
}

void
Splice::scorefxn( core::scoring::ScoreFunctionOP sf ){
	scorefxn_ = sf;
}

core::scoring::ScoreFunctionOP
Splice::scorefxn() const{
	return scorefxn_;
}

core::pack::task::TaskFactoryOP
Splice::task_factory() const{ return task_factory_; }

void
Splice::task_factory( core::pack::task::TaskFactoryOP tf ){ task_factory_ = tf; }

core::pack::task::TaskFactoryOP
Splice::design_task_factory() const{ return design_task_factory_; }

void
Splice::design_task_factory( core::pack::task::TaskFactoryOP tf ){ design_task_factory_ = tf; }


/// the torsion dbase should have the following structure:
/// each line represents a single loop. Each four values represent <phi> <psi> <omega> <3-let resid>; the last entry in a line represents <loop start> <loop stop> <cut site> cut; where cut signifies that this is the loop designator
void
Splice::read_torsion_database(){
	using namespace std;

	TR<<"Reading torsion database"<<std::endl;
	utility::io::izstream data( torsion_database_fname_ );
	if ( !data ) {
		TR << "cannot open torsion database " << torsion_database_fname_ << std::endl;
		utility_exit();
	}
	std::string line;
	while( getline( data, line ) ) {
		utility::vector1< std::string > const elements_in_line( utility::string_split( line, ' ' ) );
		if( elements_in_line.size() % 4 != 0 )
			utility_exit_with_message( "While reading torsion database "+torsion_database_fname_+" found a line where the number of elements is not divisible by 4. This likely stems from an error in the database:\n" + line );
		std::istringstream line_stream( line );
		ResidueBBDofs bbdof_entry;
		bbdof_entry.clear();
		while( !line_stream.eof() ){
			core::Real phi, psi, omega;
			std::string resn;
			line_stream >> phi >> psi >> omega >> resn;
			if( line_stream.eof() ){// the end of the line signifies that we're reading the start, stop, cut, source_pdb fields
				bbdof_entry.start_loop( (core::Size ) phi );
				bbdof_entry.stop_loop( (core::Size ) psi );
				bbdof_entry.cut_site( (core::Size ) omega );
				bbdof_entry.source_pdb( resn );
			}
			else
				bbdof_entry.push_back( BBDofs( 0/*resid*/, phi, psi, omega, resn ) ); /// resid may one day be used. Currently it isn't
		}
		torsion_database_.push_back( bbdof_entry );
	}
	TR<<"Finished reading torsion database with "<<torsion_database_.size()<<" entries"<<std::endl;
}

///@brief set the fold tree around start/stop/cut sites.
/// presently makes a simple fold tree, but at one point may be a more complicated function to include two poses
void
Splice::fold_tree( core::pose::Pose & pose, core::Size const start, core::Size const stop, core::Size const cut ) const{
	using namespace protocols::loops;
	core::conformation::Conformation const & conf( pose.conformation() );
	core::Size const s1 = std::max( (core::Size) 2, start - 6 );
	core::Size const s2 = std::min( conf.chain_end( 1 )/* - 1*/, stop + 6 );
	core::kinematics::FoldTree ft;
	ft.clear();
	if( conf.num_chains() == 1 ){/// build simple ft for the cut
		ft.add_edge( 1, s1, -1 );
		ft.add_edge( s1, s2, 1 );
		ft.add_edge( s1, cut, -1 );
		ft.add_edge( s2, cut + 1, -1 );
		ft.add_edge( s2, pose.total_residue(), -1 );
		ft.delete_self_edges();
		TR<<"single chain ft: "<<ft<<std::endl;
		pose.fold_tree( ft );
		return;
	}
	//core::Size from_res( 0 );
	for( core::Size resi = conf.chain_begin( 1 ); resi <= conf.chain_end( 1 ); ++resi ){
		if( pose.residue( resi ).has_variant_type( core::chemical::DISULFIDE ) ){
			//from_res = resi;  // set but never used ~Labonte
			break;
		}
	}
	ft.add_edge( 1, s1, -1 );
	ft.add_edge( s1, s2, 1 );
	ft.add_edge( s2, conf.chain_end( 1 ), -1 );
	if( locked_res() > 0 && (locked_res() <= s2 && locked_res() >= s1 )){
		TR<<"s1,s2,locked_res: "<<s1<<','<<s2<<','<<locked_res()<<std::endl;
		if( locked_res() < cut ){
			ft.add_edge( s1, locked_res(), -1 );
			ft.add_edge( locked_res(), cut, -1 );
			ft.add_edge( s2, cut+1, -1 );
		}
		if( locked_res() > cut ){
			ft.add_edge( s1, cut, -1 );
			ft.add_edge( s2, locked_res(), -1 );
			ft.add_edge( locked_res(), cut + 1, -1 );
		}
		if( locked_res() == cut ){
			ft.add_edge( s1, cut, -1 );
			ft.add_edge( s2, cut+1, -1 );
		}
		using namespace protocols::protein_interface_design;
		std::string const from_atom( optimal_connection_point( pose.residue( locked_res() ).name3() ) );
		core::Real min_dist( 100000 );
		core::Size nearest_res( 0 );
		core::Size nearest_atom( 0 );
		for( core::Size resi = conf.chain_begin( 2 ); resi <= conf.chain_end( 2 ); ++resi ){
			core::conformation::Residue const residue( conf.residue( resi ) );
			if( residue.is_ligand() ) continue;
			for( core::Size atomi = 1; atomi <= residue.natoms(); ++atomi ){
				core::Real const dist( conf.residue( locked_res() ).xyz( from_atom ).distance( residue.xyz( atomi ) ) );
				if( dist <= min_dist ){
					nearest_res = resi;
					nearest_atom = atomi;
					min_dist = dist;
				}
			}
		}
		runtime_assert( nearest_res );
		ft.add_edge( locked_res(), nearest_res, 2 );
		ft.add_edge( nearest_res, conf.chain_begin( 2 ), -1 );
		ft.add_edge( nearest_res, conf.chain_end( 2 ), -1 );
		ft.set_jump_atoms( 2, from_atom, conf.residue( nearest_res ).atom_name( nearest_atom ) );
	}
	else{
		if(locked_res() > 0 && ! ( locked_res() > s1 && locked_res() < s2 ) ){
			TR<<"locked_res "<<locked_res()<<" is outside loop scope so ignoring"<<std::endl;
		}
		ft.add_edge( s1, cut, -1 );
		ft.add_edge( s2, cut + 1, -1 );
		ft.add_edge( 1, conf.chain_begin( 2 ), 2 );
	}
	if( (!locked_res() || ( locked_res() <= s1 || locked_res() >= s2 ) ) && !pose.residue( conf.chain_begin( 2 ) ).is_ligand() )
		ft.add_edge( conf.chain_begin( 2 ), conf.chain_end( 2 ), -1 );
	ft.reorder(1);
	TR<<"Previous ft: "<<pose.fold_tree()<<std::endl;
	//	pose.dump_pdb( "before_ft.pdb" );
	pose.fold_tree( ft );
	//	pose.dump_pdb( "after_ft.pdb" );
	TR<<"Current ft: "<<pose.fold_tree()<<std::endl;
}

BBDofs::~BBDofs() {}

ResidueBBDofs::~ResidueBBDofs() {}

utility::vector1< core::Size >::const_iterator
Splice::dbase_begin() const{ return dbase_subset_.begin(); }

utility::vector1< core::Size >::const_iterator
Splice::dbase_end() const{ return dbase_subset_.end(); }

core::Size
Splice::locked_res() const {
	if( locked_res_ )
		return locked_res_->obj[ 1 ];
	else
		return 0;
}

void
Splice::locked_res( core::Size const r ) { locked_res_->obj[ 1 ] = r; }

void
Splice::locked_res_id( char const c ){ locked_res_id_ = c ; }

char
Splice::locked_res_id() const{ return locked_res_id_; }

std::string
Splice::checkpointing_file() const { return checkpointing_file_; }

void
Splice::checkpointing_file( std::string const cf ){ checkpointing_file_ = cf; }

void
Splice::loop_dbase_file_name( std::string const s ){
	loop_dbase_file_name_ = s;
}

std::string
Splice::loop_dbase_file_name() const{
	return loop_dbase_file_name_;
}

void
Splice::loop_pdb_source( std::string const s ){
	loop_pdb_source_ = s;
}

std::string
Splice::loop_pdb_source() const{
	return loop_pdb_source_;
}

protocols::filters::FilterOP
Splice::splice_filter() const{
	return splice_filter_;
}

void
Splice::splice_filter( protocols::filters::FilterOP f ){
	splice_filter_ = f;
}

void
Splice::read_splice_segments( std::string const segment_type, std::string const segment_name, std::string const file_name ){
	if(use_sequence_profiles_){
		splice_segments_[ segment_type ]->read_profile( file_name, segment_name );
		TR<<"In segment_type "<<segment_type_<<": reading profile for segment "<<segment_name<<" from file "<<file_name<<std::endl;
	}
}

core::sequence::SequenceProfileOP
Splice::generate_sequence_profile(core::pose::Pose & pose)
{
	if (use_sequence_profiles_){
		using namespace core::sequence;
		using namespace std;

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using protocols::jd2::JobDistributor;

		std::string pdb_tag= JobDistributor::get_instance()->current_job()->inner_job()->input_tag();

	//	std::string temp_pdb_name = pdb_tag.erase(pdb_tag.size() - 5); //JD adds "_0001" to input name, we need to erase it
		//pdb_tag = temp_pdb_name +".pdb";
		TR<<" The scaffold file name is :"<<pdb_tag<<std::endl;//file name of -s pdb file
//		core::pose::read_comment_pdb(pdb_tag,pose); //read comments from pdb file
/*		std::string pdb_dump_fname_("test2");
		std::ofstream out( pdb_dump_fname_.c_str() );
		pose.dump_pdb(out); //Testing out comment pdb, comment this out after test (GDL) */
		map< string, string > const comments = core::pose::get_all_comments( pose );
		if (comments.size()<3 ){
			utility_exit_with_message("Please check comments field in the pdb file (header= ##begin comments##), could not find any comments");
		}
		///This code will cut the source pdb file name and extract the four letter code
		if (torsion_database_fname_!=""){
			TR<<"Torsion data base filename is: "<<torsion_database_fname_<<std::endl;
			Pdb4LetName_ = dofs_pdb_name;
			//TR<<"Pdb4LetName_: "<<Pdb4LetName_<<std::endl;
/*	    std::string pdb_dump_fname_("before_splice.pdb");
		std::ofstream out( pdb_dump_fname_.c_str() );
		pose.dump_pdb(out); //Testing out comment pdb, comment this out after test (GDL) */

		}
		else {
			Pdb4LetName_ = source_pdb_;//
		}
		// Remove directory if present.
		// Do this before extension removal in case directory has a period character.
		const size_t last_slash_idx = Pdb4LetName_.find_last_of("\\/");
		if (std::string::npos != last_slash_idx)
		{
			Pdb4LetName_.erase(0, last_slash_idx + 1);
		}
		// Remove extension if present.
		const size_t period_idx = Pdb4LetName_.rfind('.');
		if (std::string::npos != period_idx) {
			Pdb4LetName_.erase(period_idx);
		}
		if (!add_sequence_constraints_only_){///If only doing sequence constraints then don't add to pose comments source name
		TR<<"The currnet segment is: "<<segment_type_<<" and the source pdb is "<<Pdb4LetName_<<std::endl;
		core::pose::add_comment(pose,"segment_"+segment_type_,Pdb4LetName_);//change correct association between current loop and pdb file
		}

		load_pdb_segments_from_pose_comments(pose); // get segment name and pdb accosiation from comments in pdb file
		TR<<"There are "<<pdb_segments_.size()<<" PSSM segments"<<std::endl;


		runtime_assert( pdb_segments_.size() ); //This assert is in place to make sure that the pdb file has the correct comments, otherwise this function will fail
//		pose.dump_pdb("test"); //Testing out comment pdb, comment this out after test (GDL)

		// test that all PDB_segments are present
		//for( map< string, string >::const_iterator i = pdb_segments_.begin(); i != pdb_segments_.end(); ++i ){
			//TR<<i->first<<std::endl;
		//}

		//	std::string segment_name_ordered [14] = {"Frm1.light", "L1","Frm2.light", "L2","Frm3.light", "L3","Frm4.light", "Frm1.heavy","H1","Frm2.heavy","H2","Frm3.heavy","H3","Frm4.heavy"}; //This string array will be used as refernce to arrange the segment profiles correctly
		utility::vector1< SequenceProfileOP > profile_vector;

		profile_vector.clear(); //this vector holds all the pdb segment profiless

		foreach( std::string const segment_type, segment_names_ordered_ ){ //<- Start of PDB segment iterator
			if (splice_segments_[ segment_type ]->pdb_profile(pdb_segments_[segment_type])==0){
				utility_exit_with_message(" could not find the source pdb name:: "+ pdb_segments_[segment_type]+ ", in pdb_profile_match file."+segment_type+"\n");
			}
			profile_vector.push_back( splice_segments_[ segment_type ]->pdb_profile( pdb_segments_[segment_type] ));
		} // <- End of PDB segment iterator
		TR<<"The size of the profile vector is: "<<profile_vector.size()<<std::endl;

		///Before upwweighting constraint we check that the PSSMs are properly aligned by making sure
		///that the PSSM score of the falnking segments of the current designed segment agree with the identity of the aa (i.e if
		/// we are designing L1 then we would expect that segments Frm.light1 and Frm.light2 have concensus aa identities)


		if (ccd_==0){//if CDD if false we are only doing splice_in, then we should check whether PSSM and pose are aligned correctly, gideonla Aug13
		core::Size aapos=0;
		//TR<<"TESTING PSSMs"<<std::endl;
		for(core::Size seg=1; seg <=profile_vector.size() ; seg++ ){//go over all the PSSM sements provided by the user
			for( core::Size pos /*go over profile ids*/ = 1; pos <= profile_vector[seg]->size(); ++pos ){
				++aapos;//go over pose residue
				//TR<<pose.residue(aapos).name1()<<aapos<<" CYS score is: "<<profile_vector[seg]->prof_row(pos)[2]<<std::endl;
				if ((profile_vector[seg]->prof_row(pos)[2])>8){//If the profile vector holds a disulfide Cys it will have a pssm score over 8
					std::stringstream ss; std::string s;
					ss << pose.residue(aapos).name1();
					ss >> s;
					//TR<<"found a dis cys="<<s<<std::endl;
					if (s.compare("C") != 0){
						std::string seqpos;         
						std::ostringstream convert;
						convert << aapos;      // insert the textual representation of 'Number' in the characters in the stream
						seqpos = convert.str(); // set 'Result' to the contents of the stream
						utility_exit_with_message(" PSSM and pose might be misaligned, position "+ s+seqpos+ " should be a CYS\n");
						//utility_exit_with_message(" could not find the source pdb name:: "+ pdb_segments_[segment_type]+ ", in pdb_profile_match file."+segment_type+"\n");
					} //fi
				}//fi
			}//end inner segment for
		}//end pssm segment for
		}




		return concatenate_profiles( profile_vector,segment_names_ordered_ );
	}
	return NULL;  // Control can reach end of non-void function without this line. ~ Labonte
}

void
Splice::load_pdb_segments_from_pose_comments( core::pose::Pose const & pose ){
	if(use_sequence_profiles_){
		//If we are using sequence profiles then the condition is true and function can run
		using namespace std;
		map< string, string > const comments = core::pose::get_all_comments( pose );
		TR<<"The size of comments is: "<<comments.size()<<std::endl;
		core::Size j = 1; //for testing
		for( std::map< string, string >::const_iterator i = comments.begin(); i != comments.end(); ++i ){
			//TR<<"the size of j is: "<<j<<std::endl;
			std::string const key( i->first );
			//TR<<"the size of j after i->first is: "<<j<<std::endl;
			std::string const val( i->second );
			//TR<<"the size of j after i->second is: "<<j<<std::endl;
			if( key.substr( 0, 7 ) != "segment" )/// the expected format is segment_??, where we're interested in ??
				continue;
			std::string const short_key( key.substr(8, 1000 ) );
			pdb_segments_[ short_key ] = val;
			TR<<"recording segment/pdb pair: "<<short_key<<'/'<<val<<std::endl;
			++j;
		}
	}
}

void
Splice::modify_pdb_segments_with_current_segment( std::string const pdb_name ){
	pdb_segments_[ segment_type_ ] = pdb_name;
}

// @brief utility function for computing which residues on chain1 are away from the interface
utility::vector1< core::Size >
find_residues_on_chain1_inside_interface( core::pose::Pose const & pose ){
	using namespace protocols::toolbox::task_operations;
	ProteinInterfaceDesignOperationOP pido = new ProteinInterfaceDesignOperation;
	//TR<<"test pido @line 1163 "<<pose<<pose<<"\n";
	pido->repack_chain1( true );
	pido->design_chain1( true );
	pido->repack_chain2( false );
	pido->design_chain2( false );
	pido->interface_distance_cutoff( 8.0 );
	core::pack::task::TaskFactoryOP tf_interface( new core::pack::task::TaskFactory );
	tf_interface->push_back( pido );
	///// FIND COMPLEMENT ////////
	utility::vector1< core::Size > const chain1_interface( protocols::rosetta_scripts::residue_packer_states( pose, tf_interface, true, true ) ); /// find packable but not designable residues; according to pido specifications above these will be on chain1 outside an 8A shell around chain2

	return chain1_interface;
}

void
Splice::add_sequence_constraints( core::pose::Pose & pose){
	if(use_sequence_profiles_){
		using namespace core::scoring::constraints;

		/// first remove existing sequence constraints
		TR<<"Removing existing sequence profile constraints from pose"<<std::endl;
		ConstraintCOPs constraints( pose.constraint_set()->get_all_constraints() );
		TR<<"Total number of constraints at start: "<<constraints.size()<<std::endl;
		core::Size cst_num( 0 );
		foreach( ConstraintCOP const c, constraints ){
			if( c->type() == "SequenceProfile" ){//only remove profile sequence constraints
				pose.remove_constraint( c );
				cst_num++;
			}
		}
		TR<<"Removed a total of "<<cst_num<<" sequence constraints."<<std::endl;
		TR<<"After removal the total number of constraints is: "<<pose.constraint_set()->get_all_constraints().size()<<std::endl;
		/// then impose new sequence constraints
		core::sequence::SequenceProfileOP seqprof( generate_sequence_profile(pose) );
		TR<<"Chain length/seqprof size: "<<pose.conformation().chain_end( 1 ) - pose.conformation().chain_begin( 1 ) + 1<<", "<<seqprof->size()-1<<std::endl;
/*		std::string pdb_dump_fname_("after_splice.pdb");
		std::ofstream out( pdb_dump_fname_.c_str() );
		pose.dump_pdb(out); //Testi*/
		runtime_assert( seqprof->size()-1 == pose.conformation().chain_end( 1 ) - pose.conformation().chain_begin( 1 ) + 1 ); //Please note that the minus 1 after seqprof size is because seqprof size is always +1 to the actual size. Do not chnage this!!
		cst_num = 0;
		TR<<"Upweighting sequence constraint for residues: ";
		if (pose.conformation().num_chains() == 1){//If pose has only one chain (no ligand) than all residues are weighted the same
			for( core::Size seqpos = pose.conformation().chain_begin( 1 ); seqpos <= pose.conformation().chain_end( 1 ); ++seqpos ) {
				TR.Debug<<"Now adding constraint to aa: "<<seqpos<<pose.aa(seqpos)<<std::endl;
				TR.Debug<<"The sequence profile row for that residue is: "<<seqprof->prof_row(seqpos)<<std::endl;
				SequenceProfileConstraintOP spc( new SequenceProfileConstraint( pose, seqpos, seqprof ) );
				//spc->weight( 1000 );
				pose.add_constraint( spc );
			}
			ConstraintCOPs constraints( pose.constraint_set()->get_all_constraints() );
			TR<<"Total number of constraints at End: "<<constraints.size()<<std::endl;
		}
		else{ //if pose has two chains than there is also a ligand therefore we weight antibody rediues according to distance from ligand
			utility::vector1< core::Size > const non_upweighted_residues( find_residues_on_chain1_inside_interface( pose ) );
			for( core::Size seqpos = pose.conformation().chain_begin( 1 ); seqpos <= pose.conformation().chain_end( 1 ); ++seqpos ){
				using namespace core::scoring::constraints;
				SequenceProfileConstraintOP spc( new SequenceProfileConstraint( pose, seqpos, seqprof ) );
				if( std::find( non_upweighted_residues.begin(), non_upweighted_residues.end(), seqpos ) == non_upweighted_residues.end() ){//seqpos not in interface so upweight
					spc->weight( profile_weight_away_from_interface());
					TR<<seqpos<<",";
				}
				TR<<std::endl;
				pose.add_constraint( spc );
				cst_num++;
			}
			TR<<"Added a total of "<<cst_num<<" sequence constraints."<<std::endl;
			TR<<"Now the pose has a total of "<<pose.constraint_set()->get_all_constraints().size()<<" constraints"<<std::endl;
		}



		/// just checking that the scorefxn has upweighted res_type_constraint
		core::Real const score_weight( scorefxn()->get_weight( core::scoring::res_type_constraint ) );
		TR<<"res_type_constraint weight is set to "<<score_weight<<std::endl;
		if( score_weight <= 0.001 )
			TR<<"Warning! res_type_constraint weight is low, even though I've just added sequence constraints to the pose! These sequence constraints will have no effect. This could be an ERROR"<<std::endl;
	}
}


void
Splice::add_coordinate_constraints( core::pose::Pose & pose, core::pose::Pose const & source_pose,core::Size nearest_to_from,core::Size nearest_to_to )
{

	//pose.dump_pdb("during_coor_constraint.pdb");
	core::scoring::constraints::ConstraintOPs cst;
	core::Size from = protocols::rosetta_scripts::find_nearest_res(pose, source_pose,nearest_to_from, 1/*chain*/ ); //The following for loop itterates over the source pose residues
	core::Size const fixed_res( protocols::rosetta_scripts::find_nearest_res(pose, source_pose,nearest_to_from, 1/*chain*/ )-3); //The fixed res is set to be 3 residues up from start loop (should be in the stem region, gideonla 2may13
	TR<<"Anchor residue for the coordinate constraint is "<<fixed_res<<std::endl;
	TR<<"Current pose CA xyz coordinate/source pdb CA xyz coordinate:"<<std::endl;
	core::id::AtomID const anchor_atom( core::id::AtomID( pose.residue( fixed_res ).atom_index( "CA" ), fixed_res ) );
		for( core::Size i = nearest_to_from; i <= nearest_to_to; ++i ){
			core::scoring::constraints::FuncOP coor_cont_fun = new core::scoring::constraints::HarmonicFunc(0.0,1);
			cst.push_back( new core::scoring::constraints::CoordinateConstraint( core::id::AtomID(pose.residue(from).atom_index("CA"), from), anchor_atom, source_pose.residue(i).atom("CA").xyz(),coor_cont_fun));
			//Print xyz coor of current pose CA atoms vs. source pose
			TR<<from<<pose.aa(from)<<" "<<pose.residue(from).atom("CA").xyz()[0]<<","<<pose.residue(from).atom("CA").xyz()[1]<<","<<pose.residue(from).atom("CA").xyz()[2]<<" / "<<
			i<<source_pose.aa(i)<<" "<<source_pose.residue(i).atom("CA").xyz()[0]<<","<<pose.residue(i).atom("CA").xyz()[1]<<","<<pose.residue(i).atom("CA").xyz()[2]<<std::endl;
			from++;
			pose.add_constraints(cst);
		}
		//scorefxn()->show(pose);
}

void
Splice::add_dihedral_constraints( core::pose::Pose & pose, core::pose::Pose const & source_pose,core::Size nearest_to_from,core::Size nearest_to_to, core::Size/* cut_site */){
	core::Size from = protocols::rosetta_scripts::find_nearest_res(pose, source_pose,nearest_to_from, 1/*chain*/ ); //The following for loop itterates over the source pose residues
	//inorder to compare the angles we keep track of the corresponding residues in the template pose
	TR<<"Apllying Dihedral constraints to pose, internal weight = "<< dihedral_const_<<", Pose/Source PDB:"<<std::endl;
	for( core::Size i = nearest_to_from; i <= nearest_to_to; ++i ){
		core::scoring::constraints::ConstraintOPs csts; //will hold Dihedral constraints
		//Set up constraints for the phi angle
		core::id::AtomID phi_resi_n( source_pose.residue_type( i ).atom_index( "N" ), i );
		numeric::xyzVector< core::Real > xyz_Ni = source_pose.residue(i).atom("N").xyz();

		core::id::AtomID phi_resj_c( source_pose.residue_type( i-1 ).atom_index( "C" ), i-1 );
		numeric::xyzVector< core::Real > xyz_Cj = source_pose.residue(i-1).atom( "C" ).xyz();

		core::id::AtomID phi_resi_co( source_pose.residue_type( i ).atom_index( "C" ), i );
		numeric::xyzVector< core::Real > xyz_Ci = source_pose.residue(i).atom( "C" ).xyz();

		core::id::AtomID phi_resi_ca( source_pose.residue_type( i ).atom_index( "CA" ), i );
		numeric::xyzVector< core::Real > xyz_Cai = source_pose.residue(i).atom( "CA" ).xyz();



		TR<<"Phi: "<<from<<pose.aa(from)<<":"<<pose.phi(from)<<" / "<<i<<source_pose.aa(i)<<":"<<numeric::dihedral_degrees(xyz_Cj,xyz_Ni,xyz_Cai,xyz_Ci) <<std::endl;
		core::scoring::constraints::FuncOP di_const_func_phi = new core::scoring::constraints::CircularHarmonicFunc((source_pose.phi(i)*numeric::constants::d::pi_2)/360,1);
		csts.push_back( new core::scoring::constraints::DihedralConstraint(phi_resj_c,phi_resi_n,phi_resi_ca,phi_resi_co, di_const_func_phi ) );
		//for debuggin comment this out

		//Set up constraints for the psi angle
		core::id::AtomID psi_resi_n( source_pose.residue_type( i ).atom_index( "N" ), i );
		xyz_Ni = source_pose.residue(i).atom("N").xyz();

		core::id::AtomID psi_resj_n( source_pose.residue_type( i+1 ).atom_index( "N" ), i+1 );
		numeric::xyzVector< core::Real > xyz_Nj = source_pose.residue(i+1).atom("N").xyz();

		core::id::AtomID psi_resi_co( source_pose.residue_type( i ).atom_index( "C" ), i );
		xyz_Ci = source_pose.residue(i).atom("C").xyz();

		core::id::AtomID psi_resi_ca( source_pose.residue_type( i ).atom_index( "CA" ), i );
		xyz_Cai = source_pose.residue(i).atom("CA").xyz();

		//for each residue the ideal angle is taken from the "donor" pdb
		core::scoring::constraints::FuncOP di_const_func_psi = new core::scoring::constraints::CircularHarmonicFunc((source_pose.psi(i)*numeric::constants::d::pi_2)/360,1);
		csts.push_back( new core::scoring::constraints::DihedralConstraint(psi_resi_n,psi_resi_ca,psi_resi_co,psi_resj_n, di_const_func_psi ) );
		TR<<"Psi: "<<from<<pose.aa(from)<<":"<<pose.psi(from)<<" / "<<i<<source_pose.aa(i)<<":"<<numeric::dihedral_degrees(xyz_Ni,xyz_Cai,xyz_Ci,xyz_Nj) <<std::endl;
		//Set up constraints for the omega angle
		core::id::AtomID omega_resj_n( source_pose.residue_type( i+1 ).atom_index( "N" ), i+1);
		xyz_Ni = source_pose.residue(i+1).atom("N").xyz();

		core::id::AtomID omega_resi_ca( source_pose.residue_type( i ).atom_index( "CA" ), i );
		xyz_Cai = source_pose.residue(i).atom("CA").xyz();

		core::id::AtomID omega_resi_co( source_pose.residue_type( i ).atom_index( "C" ), i );
		xyz_Ci = source_pose.residue(i).atom("C").xyz();

		core::id::AtomID omega_resj_ca( source_pose.residue_type( i+1 ).atom_index( "CA" ), i+1);
		numeric::xyzVector< core::Real > xyz_Caj = source_pose.residue(i+1).atom("CA").xyz();
		TR<<"omega: "<<from<<pose.aa(from)<<":"<<pose.omega(from)<<" / "<<i<<source_pose.aa(i)<<":"<<numeric::dihedral_degrees(xyz_Cai,xyz_Ci,xyz_Nj,xyz_Caj) <<std::endl;
		//for each residue the ideal angle is taken from the "donor" pdb
		core::scoring::constraints::FuncOP di_const_func_omega = new core::scoring::constraints::CircularHarmonicFunc((source_pose.omega(i)*numeric::constants::d::pi_2)/360,1);
		csts.push_back( new core::scoring::constraints::DihedralConstraint(omega_resi_ca,omega_resi_co, omega_resj_n, omega_resj_ca, di_const_func_omega ) );

		pose.add_constraints(csts);
		from++;//every itteration we must increment "from"
		}
		core::Real const score_weight( scorefxn()->get_weight( core::scoring::dihedral_constraint ) );
		TR<<"dihedral_constraint weight is set to "<<score_weight<<std::endl;
		//scorefxn()->show(pose);
		//pose.dump_pdb("at_end_of_dihedral_const.pdb");
}

core::Real
Splice::profile_weight_away_from_interface() const{
	return profile_weight_away_from_interface_;
}

void
Splice::profile_weight_away_from_interface( core::Real const p ){
	profile_weight_away_from_interface_ = p;
}


//This function is not called anymore, consider removing; gideonla Aug13
bool
Splice::check_aa(std::string curraa, utility::vector1<core::Real > profRow)
 {
	// order of amino acids in the .pssm file
	static std::map< std::string, core::Size > order;
	order["A"] = 1;
	order["C"] = 2;
	order["D"] = 3;
	order["E"] = 4;
	order["F"] = 5;
	order["G"] = 6;
	order["H"] = 7;
	order["I"] = 8;
	order["K"] = 9;
	order["L"] = 10;
	order["M"] = 11;
	order["N"] = 12;
	order["P"] = 13;
	order["Q"] = 14;
	order["R"] = 15;
	order["S"] = 16;
	order["T"] = 17;
	order["V"] = 18;
	order["W"] = 19;
	order["Y"] = 20;
			if ((profRow[order[curraa]])<0){
				return true;
			}
			else{
			return false;
		}

}

} //splice
} //devel

