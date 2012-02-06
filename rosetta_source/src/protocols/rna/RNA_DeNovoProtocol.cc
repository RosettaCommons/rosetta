// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA de novo fragment assembly
/// @brief protocols that are specific to RNA_DeNovoProtocol
/// @detailed
/// @author Rhiju Das


// Unit headers
#include <protocols/rna/RNA_DeNovoProtocol.hh>

// Package headers
#include <protocols/toolbox/AllowInsert.hh>
#include <protocols/rna/RNA_BasePairClassifier.hh>
#include <protocols/rna/RNA_DataReader.hh>
#include <protocols/rna/RNA_DataReader.fwd.hh>
#include <protocols/rna/FullAtomRNA_Fragments.hh>
#include <protocols/rna/RNA_LoopCloser.hh>
#include <protocols/rna/RNA_LoopCloser.fwd.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <protocols/rna/RNA_Minimizer.fwd.hh>
#include <protocols/rna/RNA_Minimizer.hh>
#include <protocols/rna/RNA_Relaxer.fwd.hh>
#include <protocols/rna/RNA_Relaxer.hh>
#include <protocols/rna/RNA_StructureParameters.hh>
#include <protocols/rna/RNA_ChunkLibrary.hh>
#include <protocols/rna/RNA_ChunkLibrary.fwd.hh>
#include <protocols/swa/StepWiseUtil.hh> //move this to toolbox/

// Project headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.hh>
#include <core/pose/Pose.hh>
#include <basic/database/open.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <protocols/viewer/viewers.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/scoring/rna/RNA_BaseDoubletClasses.hh>

#include <utility/file/file_sys_util.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>

// External library headers

//C++ headers
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <fstream>
#ifdef WIN32
#include <ctime>
#endif

//Auto Headers
#include <protocols/viewer/GraphicsState.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end


using namespace core;

namespace protocols {
namespace rna {

static numeric::random::RandomGenerator RG(12320);  // <- Magic number, do not change it!

static basic::Tracer TR( "protocols.rna.rna_denovo_protocol" ) ;

RNA_DeNovoProtocol::RNA_DeNovoProtocol(
	 Size const nstruct,
	 Size const monte_carlo_cycles,
	 std::string const silent_file,
	 bool const heat_structure /*= true*/,
	 bool const minimize_structure /*= false*/,
	 bool const relax_structure /*=false*/):
    Mover(),
		nstruct_( nstruct ),
		monte_carlo_cycles_( monte_carlo_cycles ),
		all_rna_fragments_file_( basic::database::full_name("chemical/rna/1jj2.torsions") ),
		silent_file_( silent_file ),
		lores_silent_file_( "" ),
		heat_structure_( heat_structure ),
		dump_pdb_( false ), //RHIJU DO NOT CHECK THIS IN AS TRUE!
		minimize_structure_( minimize_structure ),
		relax_structure_( relax_structure ),
		ignore_secstruct_( false ),
		close_loops_( false ),
		close_loops_after_each_move_( false ),
		simple_rmsd_cutoff_relax_( false ),
		m_Temperature_( 2.0 ),
		frag_size_( 3 ),
		rna_params_file_( "" ),
		rna_data_file_( "" ),
		jump_library_file_( basic::database::full_name("chemical/rna/1jj2_RNA_jump_library.dat" ) ),
		rna_structure_parameters_( RNA_StructureParametersOP( new RNA_StructureParameters ) ),
		rna_data_reader_( RNA_DataReaderOP( new RNA_DataReader ) ),
		output_lores_silent_file_( false ),
		filter_lores_base_pairs_( false ),
		vary_bond_geometry_( false ),
		binary_rna_output_( false ),
		jump_change_frequency_( 0.1 ),
		lores_scorefxn_( "rna_lores.wts" ),
		chunk_coverage_( 0.0 ),
		staged_constraints_( false )
{
	Mover::type("RNA_DeNovoProtocol");
	rna_loop_closer_ = protocols::rna::RNA_LoopCloserOP( new protocols::rna::RNA_LoopCloser );
	//	rna_loop_closer_->fast_scan( true );
	local_rna_low_resolution_potential_.more_precise_base_pair_classification( true );
}

/// @brief Clone this object
protocols::moves::MoverOP RNA_DeNovoProtocol::clone() const {
	return new RNA_DeNovoProtocol(*this);
}

//////////////////////////////////////////////////
RNA_DeNovoProtocol::~RNA_DeNovoProtocol()
{
}

/// @details  Apply the RNA de novo modeling protocol to a pose.
///
void RNA_DeNovoProtocol::apply( core::pose::Pose & pose	) {

	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::io::pdb;
	using namespace core::io::silent;
	using namespace protocols::rna;

	///////////////////////////////////////////////////////////////////////////
	// A bunch of initialization
	///////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////
	// Some useful movers...
	////////////////////////////////////////
	initialize_movers( pose );
	if (dump_pdb_) pose.dump_pdb( "init.pdb" );

	// RNA low resolution score function.
	denovo_scorefxn_ = ScoreFunctionFactory::create_score_function( lores_scorefxn_ );
	initialize_constraints( pose );
	initial_denovo_scorefxn_ = denovo_scorefxn_->clone();

	//Keep a copy for resetting after each decoy.
	Pose start_pose = pose;

	monte_carlo_ = protocols::moves::MonteCarloOP( new protocols::moves::MonteCarlo( pose, *denovo_scorefxn_, 2.0 ) );

	// Some other silent file setup
	initialize_lores_silent_file();
	initialize_tag_is_done();


	Size max_tries( 1 );
	if (filter_lores_base_pairs_)  max_tries = 10;

	///////////////////////////////////////////////////////////////////////////
	// Main Loop.
	///////////////////////////////////////////////////////////////////////////
	for (Size n = 1; n <= nstruct_; n++ ) {

		std::string const out_file_tag = "S_"+lead_zero_string_of( n, 6 );
		if (tag_is_done_[ out_file_tag ] ) continue;

		Size ntries( 0 );
		bool found_good_decoy( false );
		while( ++ntries <= max_tries && !found_good_decoy ) {

			time_t pdb_start_time = time(NULL);

			pose = start_pose;
			rna_structure_parameters_->setup_fold_tree_and_jumps_and_variants( pose );
			rna_structure_parameters_->setup_base_pair_constraints( pose ); // needs to happen after setting cutpoint variants, etc.
			rna_chunk_library_->initialize_random_chunks( pose ); //actually not random if only one chunk in each region.

			if (dump_pdb_) dump_pdb( pose, "start.pdb" );

			if (heat_structure_ ) do_random_fragment_insertions( pose );

			if (dump_pdb_) dump_pdb( pose, "random.pdb" );
			monte_carlo_->reset( pose );

			TR << "Beginning main loop... " << std::endl;

			Size const rounds = 10; //for now.
			frag_size_ = 3;

			for (Size r = 1; r <= rounds; r++ ) {

				//Keep score function coarse for early rounds.
				update_denovo_scorefxn_weights( r, rounds );

				monte_carlo_->score_function( *denovo_scorefxn_ );

				pose = monte_carlo_->lowest_score_pose();

				// Introduce constraints in stages.
				update_pose_constraints( r, rounds, pose );
				monte_carlo_->reset( pose );

				// Finer and finer fragments
				update_frag_size( r, rounds );

				//////////////////////
				// This is it ... do the loop.
				//////////////////////
				for( Size i=1; i <= monte_carlo_cycles_/rounds ; ++i ) {
					// Make this generic fragment/jump multimover next?
					RNA_move_trial( pose );
				}

				if ( get_native_pose() ) {
					Real const rmsd = all_atom_rmsd( *get_native_pose(), pose );
					TR << "All atom rmsd: " << rmsd << std::endl;
				}

				monte_carlo_->recover_low( pose );
				monte_carlo_->show_counters();
				monte_carlo_->reset_counters();
			}

			pose = monte_carlo_->lowest_score_pose();
			denovo_scorefxn_->show( std::cout, pose );

			time_t pdb_end_time = time(NULL);
			TR << "Finished fragment assembly of " << out_file_tag << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;

			if (filter_lores_base_pairs_) found_good_decoy = rna_structure_parameters_->check_base_pairs( pose );
		}

		if (output_lores_silent_file_ ) align_and_output_to_silent_file( pose, lores_silent_file_, out_file_tag );

		if (close_loops_) {
			//rna_loop_closer_->close_loops_carefully( pose, rna_structure_parameters_->connections() );
			rna_loop_closer_->apply( pose, rna_structure_parameters_->connections() );
			denovo_scorefxn_->show( std::cout, pose );
		}

		if (minimize_structure_){
			rna_minimizer_->apply( pose );

			if (close_loops_) {
				rna_loop_closer_->apply( pose, rna_structure_parameters_->connections() );
			}
		}

		if (relax_structure_)	rna_relaxer_->apply( pose );

		std::string const out_file_name = out_file_tag + ".pdb";
		if (dump_pdb_)	 dump_pdb( pose,  out_file_name );

		align_and_output_to_silent_file( pose, silent_file_, out_file_tag );

	} //nstruct

}


///////////////////////////////////////////////////////////////////////////////////////////////////////
std::string
RNA_DeNovoProtocol::get_name() const {
	return "RNA_DeNovoProtocol";
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::initialize_constraints( core::pose::Pose & pose ) {

	using namespace core::scoring;

	if (pose.constraint_set()->has_constraints() )	{
		denovo_scorefxn_->set_weight( atom_pair_constraint, 1.0 );
		constraint_set_ = pose.constraint_set()->clone();
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::initialize_movers( core::pose::Pose & pose ){

	// all jumping, secondary structure, base pair constraint, allow_insert information
	// will be stored in a .prm file.
	rna_structure_parameters_->initialize( pose, rna_params_file_, jump_library_file_, ignore_secstruct_ );

	// reads in any data on, e.g., exposure of different bases --> saves inside the pose's rna_data_info.
	rna_data_reader_->initialize( pose, rna_data_file_ );

	all_rna_fragments_ = RNA_FragmentsOP( new FullAtomRNA_Fragments( all_rna_fragments_file_ ) );

	if ( chunk_res_.size() > 0 ){
		rna_chunk_library_ = RNA_ChunkLibraryOP( new RNA_ChunkLibrary( chunk_silent_files_, pose, chunk_res_ ) );
	} else {
		rna_chunk_library_ = RNA_ChunkLibraryOP( new RNA_ChunkLibrary( chunk_silent_files_, pose, rna_structure_parameters_->connections() ) );
	}


	chunk_coverage_ = rna_chunk_library_->chunk_coverage();
	TR << "CHUNK_COVERAGE: " << chunk_coverage_ << std::endl;
	rna_structure_parameters_->allow_insert()->and_allow_insert( rna_chunk_library_->allow_insert() );

	//	rna_structure_parameters_->allow_insert()->show();


	// Do this in the main loop to ensure diverse fold trees
	//	rna_structure_parameters_->setup_fold_tree_and_jumps_and_variants( pose );

	//	std::cout << "allow insert: after changes to fold tree and variants! " << std::endl;
	//	rna_structure_parameters_->allow_insert()->show();

	rna_chunk_library_->set_allow_insert( rna_structure_parameters_->allow_insert() );

	// do this in main loop to ensure diversity
	//	rna_chunk_library_->initialize_random_chunks( pose, dump_pdb_ );

	rna_fragment_mover_ = RNA_FragmentMoverOP( new RNA_FragmentMover( all_rna_fragments_, rna_structure_parameters_->allow_insert() ) );

	rna_minimizer_ = RNA_MinimizerOP( new RNA_Minimizer );
	rna_minimizer_->set_allow_insert( rna_structure_parameters_->allow_insert() );
	rna_minimizer_->vary_bond_geometry( vary_bond_geometry_ );

	rna_relaxer_ = RNA_RelaxerOP( new RNA_Relaxer( rna_fragment_mover_, rna_minimizer_) );
	rna_relaxer_->simple_rmsd_cutoff_relax( simple_rmsd_cutoff_relax_ );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::initialize_tag_is_done()
{

	using namespace core::io::silent;

	tag_is_done_.clear();

	utility::vector1< std::string > tags_done;

	SilentFileData silent_file_data;
	if ( utility::file::file_exists( silent_file_ ) ) {
		tags_done = silent_file_data.read_tags_fast( silent_file_ );
		for ( utility::vector1< std::string >::const_iterator iter = tags_done.begin(); iter != tags_done.end(); iter++ ) {
			std::cout << "Already done: " << *iter << std::endl;
			tag_is_done_[ *iter ] = true;
		}
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::initialize_lores_silent_file() {

	if ( !output_lores_silent_file_ ) return;

	static std::string const new_prefix( "_LORES.out" );

	std::string::size_type pos = silent_file_.find( ".out", 0 );
	if (pos == std::string::npos ){
		utility_exit_with_message(  "If you want to output a lores silent file, better use .out suffix ==> " + silent_file_ );
	}
	lores_silent_file_ = silent_file_;
	lores_silent_file_.replace( pos, new_prefix.length(), new_prefix );
}

//////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::calc_rmsds( core::io::silent::SilentStruct & s, core::pose::Pose & pose, std::string const & out_file_tag ) const
{
	using namespace core::scoring;

	Real const rmsd = all_atom_rmsd( *get_native_pose(), pose );
	TR << "All atom rmsd: " << rmsd << " for " << out_file_tag << std::endl;
	s.add_energy( "rms", rmsd );

	Real rmsd_stems = 0.0;
	std::list< Size > stem_residues( rna_structure_parameters_->get_stem_residues( pose ) );

	if ( stem_residues.size() > 0 ) {
		rmsd_stems = all_atom_rmsd( *get_native_pose(), pose, stem_residues );
		TR << "All atom rmsd over stems: " << rmsd_stems << " for " << out_file_tag << std::endl;
	}
	s.add_energy( "rms_stem", rmsd_stems );

}

///////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::output_silent_struct(
										core::io::silent::SilentStruct & s, core::io::silent::SilentFileData & silent_file_data,
										std::string const & silent_file, pose::Pose & pose, std::string const out_file_tag,
										bool const score_only /* = false */) const
{

	using namespace core::io::silent;
	using namespace core::scoring;

	if ( get_native_pose() ) calc_rmsds( s, pose, out_file_tag  );

	add_number_base_pairs( pose, s );
	if ( get_native_pose() ) add_number_native_base_pairs( pose, s );

	TR << "Outputting to silent file: " << silent_file << std::endl;
	silent_file_data.write_silent_struct( s, silent_file, score_only );

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::output_to_silent_file( core::pose::Pose & pose, std::string const & silent_file, std::string const & out_file_tag, bool const score_only /* = false */) const
{

	using namespace core::io::silent;
	using namespace core::scoring;

	// Silent file setup?
	//static SilentFileData silent_file_data;
	SilentFileData silent_file_data;

	// What is all this rigamarole, making the silent struct data?
	// Why do I need to supply the damn file name? That seems silly.
	TR << "Making silent struct for " << out_file_tag << std::endl;

	if ( binary_rna_output_ ) {
		BinaryRNASilentStruct s( pose, out_file_tag );
		output_silent_struct( s, silent_file_data, silent_file, pose, out_file_tag, score_only );
	} else {
		RNA_SilentStruct s( pose, out_file_tag );
		output_silent_struct( s, silent_file_data, silent_file, pose, out_file_tag, score_only );
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::align_and_output_to_silent_file( core::pose::Pose & pose, std::string const & silent_file, std::string const & out_file_tag ) const
{

	if ( get_native_pose() ){
		Pose const & native_pose = *get_native_pose();

		//realign to native for ease of viewing.
		// check for any fixed domains.
		utility::vector1< Size > superimpose_res;

		protocols::toolbox::AllowInsertOP const & allow_insert = rna_structure_parameters_->allow_insert();
		for( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( !allow_insert->get( n ) ) superimpose_res.push_back( n );
		}

		// if no fixed domains, just superimpose over all residues.
		if ( superimpose_res.size() == 0 ){
			for( Size n = 1; n <= pose.total_residue(); n++ )  superimpose_res.push_back( n );
		}

		id::AtomID_Map< id::AtomID > const & alignment_atom_id_map_native =
			protocols::swa::create_alignment_id_map( pose, native_pose, superimpose_res ); // perhaps this should move to toolbox.

		std::cout << "Aligning pose to native." << std::endl;

		//pose.dump_pdb( "before_align.pdb");
		//		native_pose.dump_pdb( "native.pdb" );
		core::scoring::superimpose_pose( pose, native_pose, alignment_atom_id_map_native );
		//		pose.dump_pdb( "after_align.pdb");

	}

	output_to_silent_file( pose, silent_file, out_file_tag );
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::do_random_fragment_insertions( core::pose::Pose & pose ) {

	rna_chunk_library_->initialize_random_chunks( pose );

	if (dump_pdb_) pose.dump_pdb( "add_chunks.pdb" );

	Size const heat_cycles = 3 * pose.total_residue();
	TR << "Heating up... " << std::endl;

	for (Size i = 1; i <= heat_cycles; i++ ){
		rna_fragment_mover_->random_fragment_insertion( pose, 1 /*frag_size*/ );
	}

	if (dump_pdb_) 	pose.dump_pdb( "random_frag1.pdb" );

	rna_chunk_library_->initialize_random_chunks( pose );

	if (dump_pdb_) 	pose.dump_pdb( "random_frag2.pdb" );

}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::update_denovo_scorefxn_weights( Size const & r, Size const & rounds )
{
	using namespace core::scoring;

	Real const rna_base_axis_final_weight        = initial_denovo_scorefxn_->get_weight( rna_base_axis );
	Real const rna_base_stagger_final_weight     = initial_denovo_scorefxn_->get_weight( rna_base_stagger );
	Real const rna_base_stack_axis_final_weight  = initial_denovo_scorefxn_->get_weight( rna_base_stack_axis );
	Real const linear_chainbreak_final_weight    = initial_denovo_scorefxn_->get_weight( linear_chainbreak );
	Real const atom_pair_constraint_final_weight = initial_denovo_scorefxn_->get_weight( atom_pair_constraint );

	//Keep score function coarse for early rounds.
	Real const suppress  = (r - 1.0)/(rounds - 1.0);

	denovo_scorefxn_->set_weight( rna_base_axis,      suppress*rna_base_axis_final_weight  );
	denovo_scorefxn_->set_weight( rna_base_stagger,   suppress*rna_base_stagger_final_weight  );
	denovo_scorefxn_->set_weight( rna_base_stack_axis,suppress*rna_base_stack_axis_final_weight  );
	denovo_scorefxn_->set_weight( linear_chainbreak,  suppress*linear_chainbreak_final_weight  );
	denovo_scorefxn_->set_weight( atom_pair_constraint,  suppress*atom_pair_constraint_final_weight  );
}


////////////////////////////////////////////////////////////////////////////////////////
Size
RNA_DeNovoProtocol::figure_out_constraint_separation_cutoff( Size const & r, Size const & rounds, Size const & max_dist )
{

	//Keep score function coarse for early rounds.
	Real const suppress  = ( r )/(rounds - 4.0);

	Size separation_cutoff = static_cast< Size > ( suppress * max_dist ) + 2;
	if ( separation_cutoff > max_dist ) separation_cutoff = max_dist;

	return separation_cutoff;

}


////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::update_pose_constraints( Size const & r, Size const & rounds, core::pose::Pose & pose )
{
	using namespace core::scoring::constraints;

	if ( !staged_constraints_) return;

	if ( !constraint_set_ ) return;

	ConstraintSetOP cst_set_new( new scoring::constraints::ConstraintSet );

	static core::kinematics::ShortestPathInFoldTree shortest_path_in_fold_tree( pose.fold_tree() );
	Size const separation_cutoff = figure_out_constraint_separation_cutoff( r, rounds, shortest_path_in_fold_tree.max_dist() );
	TR << "ROUND " << r << " out of " << rounds << std::endl;
	TR << "FOLD_TREE CURRENT SEPARATION CUTOFF " << separation_cutoff << " out of " << shortest_path_in_fold_tree.max_dist() << std::endl;

	ConstraintCOPs csts( constraint_set_->get_all_constraints() );

	for ( Size n = 1; n <= csts.size(); n++ ) {

		ConstraintCOP const & cst( csts[n] );

		if ( cst->natoms() == 2 )  { // currently only defined for pairwise distance constraints.
			Size const i = cst->atom( 1 ).rsd();
			Size const j = cst->atom( 2 ).rsd();
			Size const dist( shortest_path_in_fold_tree.dist( i , j ) );
			if ( dist  > separation_cutoff ) continue;
		}

		cst_set_new->add_constraint( cst );
	}

	pose.constraint_set( cst_set_new );

	TR << "NUM CONSTRAINTS " << pose.constraint_set()->get_all_constraints().size() << " out of " <<
		csts.size() << std::endl;

}



////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::update_frag_size( Size const & r, Size const & rounds )
{
	frag_size_ = 3;
	if ( r > 1.0 * (rounds/3.0) ) frag_size_ = 2;
	if ( r > 2.0 * (rounds/3.0) ) frag_size_ = 1;
	TR << "Fragment size: " << frag_size_ << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::RNA_move_trial( pose::Pose & pose ) {

	//There are now two kinds of insertions --
	// (1) fragment insertions for, e.g., contiguous 3-mers
	//   and
	// (2) "chunk insertions", which change out whole loops, motifs, or
	//     junctions based on previous models stored in silent files
	//

	if ( RG.uniform() < chunk_coverage_ ) {
		random_chunk_trial( pose );
	} else {
		random_fragment_trial( pose );
	}


	//Following returns early if there are no jumps.
	if  ( RG.uniform() < jump_change_frequency_ )  random_jump_trial( pose );

}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::random_jump_trial( pose::Pose & pose ) {

	bool success( false );
	success = rna_structure_parameters_->random_jump_change( pose );
	if (!success) return;

	if ( close_loops_after_each_move_ ) rna_loop_closer_->apply( pose );

	monte_carlo_->boltzmann( pose, "jump_change" );

}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::random_fragment_trial( pose::Pose & pose ) {

	rna_fragment_mover_->random_fragment_insertion( pose, frag_size_ );
	if ( close_loops_after_each_move_ ) rna_loop_closer_->apply( pose );

	monte_carlo_->boltzmann( pose, "frag" + SS(frag_size_) );

}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::random_chunk_trial( pose::Pose & pose ) {

	//	if ( frag_size_ == 2 ) {
	//		pose.dump_pdb( "before_chunk.pdb" );
	//		std::cout << "BEFORE: " << (*denovo_scorefxn_)( pose ) << std::endl;
	//		denovo_scorefxn_->show( std::cout, pose  );
	//	}

	rna_chunk_library_->random_chunk_insertion( pose );

	//	if ( frag_size_ == 2 ) {
	//		pose.dump_pdb( "after_chunk.pdb" );
	//		std::cout << "AFTER: " << (*denovo_scorefxn_)( pose ) << std::endl;
	//		denovo_scorefxn_->show( std::cout, pose );
	// 	}

	if ( close_loops_after_each_move_ ) rna_loop_closer_->apply( pose );

	monte_carlo_->boltzmann( pose, "chunk" );

	//	if ( frag_size_ == 2 ) {
	//		utility_exit_with_message(  "After chunk." );
	//	}

}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
// Following may better fit in a util.cc ,or pose_metrics...
void
RNA_DeNovoProtocol::add_number_base_pairs( pose::Pose const & pose, io::silent::SilentStruct & s ) const
{
	using namespace scoring::rna;
	using namespace conformation;

	//	pose::Pose pose = pose_input;
	//	(*denovo_scorefxn_)( pose );
	//	local_rna_low_resolution_potential_.update_rna_base_pair_list( pose );

	//	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	//	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	//	Energy_base_pair_list const & scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	utility::vector1< core::scoring::rna::Base_pair > base_pair_list;
	utility::vector1< bool > is_bulged;
	classify_base_pairs( pose, base_pair_list, is_bulged );

	Size N_WC( 0 ), N_NWC( 0 );

	//	for ( Energy_base_pair_list::const_iterator it = scored_base_pair_list.begin();
	//				it != scored_base_pair_list.end(); ++it ){
	for ( Size n = 1; n <= base_pair_list.size(); n++ ) {

		Base_pair const base_pair = base_pair_list[ n ];

		Size const i = base_pair.res1;
		Size const j = base_pair.res2;

		Size const k = base_pair.edge1;
		Size const m = base_pair.edge2;

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		if ( ( k == WATSON_CRICK && m == WATSON_CRICK
					 && base_pair.orientation == 1 )  &&
				 possibly_canonical( rsd_i.aa(), rsd_j.aa() ) )		{
			N_WC++;
		} else {
			N_NWC++;
		}
	}

 	s.add_string_value( "N_WC",  ObjexxFCL::fmt::I( 9, N_WC) );
	s.add_string_value( "N_NWC", ObjexxFCL::fmt::I( 9, N_NWC ) );
	s.add_string_value( "N_BS",  ObjexxFCL::fmt::I( 9, get_number_base_stacks( pose ) ) );

 	//s.add_energy( "N_WC",  N_WC );
	//	s.add_energy( "N_NWC", N_NWC );
	//	s.add_energy( "N_BS",  get_number_base_stacks( pose ) );
}

/////////////////////////////////////////////////////////////////////
bool
check_in_base_pair_list( scoring::rna::Base_pair const & base_pair /*from native*/,
												 utility::vector1< core::scoring::rna::Base_pair > const & base_pair_list /*for decoy*/)
{
	using namespace scoring::rna;

	bool in_list( false );

	for ( Size n = 1; n <= base_pair_list.size(); n++ ) {

		Base_pair const base_pair2 = base_pair_list[ n ];

		if ( ( base_pair.res1 == base_pair2.res1 && base_pair.res2 == base_pair2.res2 )  &&
				 ( base_pair.edge1 == base_pair2.edge1 && base_pair.edge2 == base_pair2.edge2 )  &&
				 base_pair.orientation == base_pair2.orientation ) {
			in_list = true;
			break;
		}

		if ( ( base_pair.res2 == base_pair2.res1 && base_pair.res1 == base_pair2.res2 )  &&
				 ( base_pair.edge2 == base_pair2.edge1 && base_pair.edge1 == base_pair2.edge2 )  &&
				 base_pair.orientation == base_pair2.orientation ) {
			in_list = true;
			break;
		}

	}

	return in_list;

}

/////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::add_number_native_base_pairs(pose::Pose & pose, io::silent::SilentStruct & s ) const
{
	if ( !get_native_pose() ) return;

	using namespace scoring::rna;
	using namespace conformation;

	pose::Pose native_pose = *get_native_pose();

	utility::vector1< core::scoring::rna::Base_pair > base_pair_list;
	utility::vector1< bool > is_bulged;
	classify_base_pairs( pose, base_pair_list, is_bulged );

	utility::vector1< core::scoring::rna::Base_pair > base_pair_list_native;
	utility::vector1< bool > is_bulged_native;
	classify_base_pairs( native_pose, base_pair_list_native, is_bulged_native );


	//(*denovo_scorefxn_)( pose );
	//	(*denovo_scorefxn_)( native_pose );
	//local_rna_low_resolution_potential_.update_rna_base_pair_list( native_pose );
	//local_rna_low_resolution_potential_.update_rna_base_pair_list( pose );

	//	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	//	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	//	Energy_base_pair_list const & scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	//	RNA_ScoringInfo const & rna_scoring_info_native( rna_scoring_info_from_pose( native_pose ) );
	//	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info_native( rna_scoring_info_native.rna_filtered_base_base_info() );
	//	Energy_base_pair_list const & scored_base_pair_list_native( rna_filtered_base_base_info_native.scored_base_pair_list() );

	Size N_WC_NATIVE( 0 ), N_NWC_NATIVE( 0 );
	Size N_WC( 0 ), N_NWC( 0 );

	//std::cout << "BASE PAIR LIST " << std::endl;
	for ( Size n = 1; n <= base_pair_list_native.size(); n++ ) {

		//Real const score = it->first;
		//		Real const SCORE_CUTOFF( -1.0 );
		//		if (score > SCORE_CUTOFF) continue;

		core::scoring::rna::Base_pair const base_pair = base_pair_list_native[ n ];

		Size const i = base_pair.res1;
		Size const j = base_pair.res2;

		Size const k = base_pair.edge1;
		Size const m = base_pair.edge2;

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		//std::cout << " NATIVE BASE PAIR " << i << " " << j << " " << k << " " << m << " " << it->first << std::endl;

		if ( ( k == WATSON_CRICK && m == WATSON_CRICK
					 && base_pair.orientation == 1 )  &&
				 possibly_canonical( rsd_i.aa(), rsd_j.aa() ) )		{
			N_WC_NATIVE++;
			if ( check_in_base_pair_list( base_pair /*from native*/, base_pair_list /*for decoy*/) ) N_WC++;
		} else {
			N_NWC_NATIVE++;
			if ( check_in_base_pair_list( base_pair /*from native*/, base_pair_list /*for decoy*/) ){
				N_NWC++;
			} else {
				std::cout << "Missing native base pair " << pose.residue( i ).name1() << i << " " << pose.residue(j).name1() << j << "  " << get_edge_from_num( k ) << " " << get_edge_from_num( m ) << " " << std::endl;
			}
		}
	}

	Real f_natWC( 0.0 ), f_natNWC( 0.0 ), f_natBP( 0.0 );
	if (N_WC_NATIVE > 0 ) f_natWC = ( N_WC / (1.0 * N_WC_NATIVE) );
	if (N_NWC_NATIVE > 0 ) f_natNWC = ( N_NWC / (1.0 * N_NWC_NATIVE) );
	if ( (N_WC_NATIVE + N_NWC_NATIVE) > 0 ) f_natBP = ( (N_WC+N_NWC) / (1.0 * (N_WC_NATIVE + N_NWC_NATIVE) ));

	s.add_energy( "f_natWC" , f_natWC );
	s.add_energy( "f_natNWC", f_natNWC );
	s.add_energy( "f_natBP" , f_natBP );

}


} // namespace rna
} // namespace protocols
