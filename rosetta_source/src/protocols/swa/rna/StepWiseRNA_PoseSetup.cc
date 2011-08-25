// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_PoseSetup
/// @brief Sets up pose and job parameters for RNA stepwise building.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_PoseSetup.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
//#include <protocols/rna/RNA_StructureParameters.hh>
//#include <protocols/swa/rna/StepWiseRNA_Util.hh>

//////////////////////////////////

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/database/open.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <utility/exit.hh>
#include <time.h>

#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>


using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise sampling of proteins (and probably other
// biopolymers soon). Take a starting pose and a list of residues to sample,
//  and comprehensively sample all backbone torsion angles by recursion.
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.stepwise_residue_sampler" ) ;

//typedef std::map< core::Size, core::Size > ResMap;

namespace protocols {
namespace swa {
namespace rna {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseRNA_PoseSetup::StepWiseRNA_PoseSetup( Size const & moving_res,
																								std::string const & desired_sequence,
																								utility::vector1< std::string > const & input_tags,
																								utility::vector1< std::string > const & silent_files_in,
																								utility::vector1< core::Size > const & input_res,
																								utility::vector1< core::Size > const & input_res2,
																								utility::vector1< core::Size > const & cutpoint_open,
																								Size const & cutpoint_closed ):
		moving_res_( moving_res ),
		desired_sequence_( desired_sequence ),
		rsd_set_( core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ) ),
		input_tags_( input_tags ),
		silent_files_in_( silent_files_in ),
		cutpoint_open_( cutpoint_open ),
		cutpoint_closed_( cutpoint_closed ),
		is_cutpoint_( desired_sequence_.size(), false ),
		job_parameters_( new StepWiseRNA_JobParameters ),
		virtualize_5prime_phosphates_( true )
  {
		///////////////////////////////////////////////////////
		// Cutpoint setup
		for ( Size n = 1; n <= cutpoint_open_.size();   n++ ) {
			is_cutpoint_( cutpoint_open_[ n ] ) = true;
			if ( cutpoint_open_[ n ] == cutpoint_closed_ ) utility_exit_with_message( "Position cannot be both cutpoint_open and cutpoint_closed" );
		}
		if ( cutpoint_closed > 0 ) is_cutpoint_( cutpoint_closed_ ) = true;

		///////////////////////////////////////////////////////
		assert( input_tags_.size() >= 1 );
		assert( input_tags_.size() <= 2 );

		input_res_vectors_.push_back( input_res );
		input_res_vectors_.push_back( input_res2 );
	}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseRNA_PoseSetup::~StepWiseRNA_PoseSetup()
  {}

	//////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::apply( core::pose::Pose & pose ) {

		using namespace core::pose;
		using namespace core::kinematics;
		using namespace core::chemical;

		// Define chain boundaries, mapping to working pose, etc.
		figure_out_working_sequence_and_mapping();
		figure_out_jump_partners();
		figure_out_cuts();

		// actually make the pose, set fold tree, copy in starting
		//  templates from disk.
		make_pose( pose );
		read_input_pose_and_copy_dofs( pose );

		/////////////////////////////////////////////////
		// useful in job definition. Need to carry out
		// following in exactly the right order!
		figure_out_Prepend_Internal( pose );
		figure_out_partition_definition( pose ); // this also re-roots the fold tree.
		figure_out_gap_size_and_five_prime_chain_break_res();

		/////////////////////////////////
		reroot_fold_tree( pose );
		apply_cutpoint_variants( pose );
		check_close_chain_break( pose );
		if ( virtualize_5prime_phosphates_ ) apply_virtual_phosphate_variants( pose );
		apply_bulge_variants( pose );

		// slice up native pose.
		slice_native();

		//Misc. parameters for job
		job_parameters_->set_working_fixed_res(  apply_full_to_sub_mapping( fixed_res_ ) );
		job_parameters_->set_working_terminal_res(  apply_full_to_sub_mapping( terminal_res_ ) );

		pose.dump_pdb( "start.pdb" );
	}

// 	////////////////////////////////////////////////////////////////
// 	void
// 	StepWiseRNA_PoseSetup::check_moving_res_in_chain( Size const & start_chain, Size const & end_chain,
// 																										Size const & num_chains, Size & which_chain_has_moving_res ) {

// 		if ( start_chain <= moving_res_ &&
// 				 end_chain >= moving_res_ ) which_chain_has_moving_res = num_chains;

// 		if ( which_chain_has_moving_res == 0 ) {
// 			// moving_res may be built from scratch. But it still should be adjacent to a chain!
// 			if ( moving_res_ == (start_chain - 1)  && !is_cutpoint_( moving_res_) ) which_chain_has_moving_res = num_chains;
// 			if ( moving_res_ == (end_chain + 1)  && !is_cutpoint_( end_chain ) )    which_chain_has_moving_res = num_chains;
// 		}

// }

std::string
StepWiseRNA_PoseSetup::get_name() const {
	return "StepWiseRNA_PoseSetup";
}

	/////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::figure_out_working_sequence_and_mapping(){

		utility::vector1< std::pair< core::Size, core::Size > > chain_boundaries;
		core::Size num_chains;
		std::map< core::Size, core::Size > full_to_sub;
		std::string working_sequence;

		// There are up to two input poses. Need to merge their input residues,
		// and figure out chain boundaries.
		Size const nres( desired_sequence_.size() );
		ObjexxFCL::FArray1D< Size > is_working_res( nres );
		is_working_res = 0;
		for ( Size i = 1; i <= input_res_vectors_.size(); i++ ) {
			for ( Size n = 1; n <= input_res_vectors_[i].size(); n++ )  {
				is_working_res( input_res_vectors_[i][ n ] ) = i;
			}
		}
		is_working_res( moving_res_ ) = 999;

		Size start_chain( 0 );
		Size end_chain( 0 );

		Size n( 0 );
		for ( Size pos = 1; pos <= nres; pos++ ) {

			if ( !is_working_res( pos ) ) continue;
			n++;

			if ( n == 1 ) start_chain = pos;

			if (n > 1 &&
					( pos > end_chain + 1 ||
						is_cutpoint_( end_chain ) ) ) {

				chain_boundaries.push_back( std::make_pair( start_chain, end_chain ) );
				//				std::cout << "FOUND CHAIN " << start_chain << " " << end_chain << std::endl;
				//check_moving_res_in_chain( start_chain, end_chain, chain_boundaries.size(), which_chain_has_moving_res );

				start_chain = pos;
			}

			end_chain = pos;
		}

		// For now, need to have at least one chain defined in the input!
		assert( start_chain > 0 );
		chain_boundaries.push_back( std::make_pair( start_chain, end_chain ) );
		//		check_moving_res_in_chain( start_chain, end_chain, chain_boundaries.size(), which_chain_has_moving_res );

		num_chains = chain_boundaries.size();

		//		if  ( which_chain_has_moving_res== 0 ) {
		//			utility_exit_with_message( "Could not figure out which chain to attach moving_res to!" );
		//		}

		working_sequence= "";
		full_to_sub.clear();
		Size count( 0 );
		for ( Size i = 1; i <= desired_sequence_.size(); i++ ) {
			if ( is_working_res( i ) ) {
				working_sequence += desired_sequence_[ i-1 ];
				count++;
				full_to_sub[ i ] = count;
			}
		}

		Size const working_moving_res = full_to_sub[ moving_res_ ];

		job_parameters_->set_sequence( desired_sequence_ );
		job_parameters_->set_working_sequence( working_sequence );
		job_parameters_->set_moving_res( moving_res_ );
		job_parameters_->set_working_moving_res( working_moving_res );
		job_parameters_->set_is_working_res( is_working_res );
		job_parameters_->set_full_to_sub( full_to_sub );
		job_parameters_->set_chain_boundaries( chain_boundaries );
		//job_parameters_->set_which_chain_has_moving_res( which_chain_has_moving_res );

	}

	///////////////////////////////////////////////////////////////////////////////////
	// Changed this to match Parin's favorite convention --
	// jumps far away from rebuilt residue.
	void
	StepWiseRNA_PoseSetup::figure_out_jump_partners() {

		jump_partners_.clear();

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries(	 job_parameters_->chain_boundaries() );
		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
		ObjexxFCL::FArray1D< Size > const & is_working_res = job_parameters_->is_working_res();

		Size const num_chains( chain_boundaries.size() );

		for ( Size n = 1 ; n < num_chains; n++ ) {
			Size const jump_partner1 = chain_boundaries[n].second;
			Size const jump_partner2 = chain_boundaries[n+1].first;

			// If there are two input poses, make sure their
			// rigid body arrangment is sampled --> no jumps between these
			// different regions!
 			if ( is_working_res( jump_partner1 ) == 1 &&
					 is_working_res( jump_partner2 ) == 2 ) continue;
			if ( is_working_res( jump_partner2 ) == 1 &&
					 is_working_res( jump_partner1 ) == 2 ) continue;

			if (moving_res_ == jump_partner1 ) continue;
			if (moving_res_ == jump_partner2 ) continue;

			jump_partners_.push_back( std::make_pair( full_to_sub[ jump_partner1 ], full_to_sub[ jump_partner2 ] ) );
		}

		if ( jump_partners_.size() < (num_chains-1) ) {
			// Need another jump ... "circularize".
			Size const jump_partner1 = chain_boundaries[ 1 ].first;
			Size const jump_partner2 = chain_boundaries[ num_chains ].second;
			jump_partners_.push_back( std::make_pair( full_to_sub[ jump_partner1 ], full_to_sub[ jump_partner2 ] ) );
		}


		if ( jump_partners_.size() < (num_chains-1) ) {
			// Still not enough jumps. This can happen if there are two input poses
			// and there are internal cutpoints (yea weird). Following is not too elegant
			// and may break!!!
			for ( Size n = 1 ; n < num_chains; n++ ) {
				Size const jump_partner1 = chain_boundaries[n].second;
				Size const jump_partner2 = chain_boundaries[n+1].first;
				if ( is_working_res( jump_partner1 ) == 1 &&
						 is_working_res( jump_partner2 ) == 2 ) continue;
				if ( is_working_res( jump_partner2 ) == 1 &&
						 is_working_res( jump_partner1 ) == 2 ) continue;

				if (moving_res_ == jump_partner1 || moving_res_ == jump_partner2 ) {
					// passed over this jump connection in the first pass. Try it now.
					jump_partners_.push_back( std::make_pair( full_to_sub[ jump_partner1 ], full_to_sub[ jump_partner2 ] ) );
					// are we done yet?
					if ( jump_partners_.size() == (num_chains-1) ) break;
				}

			}
		}

	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::figure_out_cuts() {

		cuts_.clear();

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries(	 job_parameters_->chain_boundaries() );
		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
		Size const num_chains = chain_boundaries.size();

		for ( Size n = 1; n < num_chains; n++ ) {
			cuts_.push_back( full_to_sub[ chain_boundaries[n].second ] );
		}

	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::make_pose( pose::Pose & pose ){

		using namespace core::conformation;

		std::string const working_sequence = job_parameters_->working_sequence();
		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );

		core::pose::make_pose_from_sequence( pose, working_sequence, *rsd_set_, false /*auto_termini*/);

		Size const nres( pose.total_residue() );
		core::kinematics::FoldTree f( nres );
		assert( cuts_.size() == jump_partners_.size() );
		Size const num_cuts( cuts_.size() );

		ObjexxFCL::FArray2D< int > jump_point( 2, num_cuts, 0 );
		ObjexxFCL::FArray1D< int > cuts( num_cuts, 0 );

		for ( Size i = 1; i <= num_cuts; i++ ) {
			jump_point( 1, i ) = jump_partners_[i].first;
			jump_point( 2, i ) = jump_partners_[i].second;
			cuts( i ) = cuts_[ i ];
			//			std::cout << " HELLO: " << jump_point( 1, i ) << " " << jump_point( 2, i ) << " " << cuts( i ) << std::endl;
		}

		Size const root_res = ( full_to_sub[ moving_res_ ] == 1 ) ? nres : 1;

		f.tree_from_jumps_and_cuts( nres, num_cuts, jump_point, cuts, root_res );

		for ( Size i = 1; i <= num_cuts; i++ ) {
			Size const k = f.upstream_jump_residue( i );
			Size const m = f.downstream_jump_residue( i );

			Residue const & rsd1( pose.residue( k ) );
			Residue const & rsd2( pose.residue( m ) );

			f.set_jump_atoms( i, rsd1.atom_name( rsd1.chi_atoms(1)[4] ), rsd2.atom_name( rsd2.chi_atoms(1)[4] ) );

		}

		pose.fold_tree( f );

	}


	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::read_input_pose_and_copy_dofs( pose::Pose & pose )
	{

		using namespace core::chemical;

		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );

		Pose start_pose;
		for ( Size i = 1; i <= input_tags_.size(); i++ ){
			Import_pose( i, start_pose );

			utility::vector1< Size > const & input_res = input_res_vectors_[i];

			// No virtual anything!
			for ( Size n = 1; n <= start_pose.total_residue(); n++  ) {
				core::pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", n );
				core::pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_O2STAR_HYDROGEN", n );
				core::pose::remove_variant_type_from_pose_residue( pose, "CUTPOINT_LOWER", n );
				core::pose::remove_variant_type_from_pose_residue( pose, "CUTPOINT_UPPER", n );
			}

			start_pose.dump_pdb( "import.pdb" );

			/////////////////////////////////////////////////////
			//Now actually copy into the pose.
			// Need to know correspondence of residues between imported pose and the working pose...
			if ( input_res.size() != start_pose.total_residue() ) {
				utility_exit_with_message( "Need to specify -already_built_residues, and the number of residues needs to match pdb file inputted after -s" );
			}

			std::map< core::Size, core::Size > res_map;
			for ( Size n = 1; n <= input_res.size(); n++ ) {
				res_map[ full_to_sub[ input_res[n] ] ] = n;
				//std::cout << full_to_sub_[ input_res_[ n ] ] << " " << n << std::endl;
			}
			//		std::cout <<  pose.annotated_sequence( true ) << " " << start_pose.annotated_sequence( true ) << std::endl;

			//copy_dofs( pose, start_pose, res_map, true /*copy_dofs_for_junction_residues*/ );
			copy_dofs( pose, start_pose, res_map );

			//		pose.dump_pdb( "copy_dof.pdb" );
			std::cout << pose.fold_tree() << std::endl;
		}


	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::figure_out_partition_definition( pose::Pose const & pose ){

		// trick to figure out which residues are upstream vs. downstream of the moving suite --
		// there's already a fold_tree function to do this, but it partitions based on a JUMP.
		//  So put in a fake jump between the moving_residue and the neighbor it is connected to.
		Size const & moving_suite( job_parameters_->working_moving_suite() );
		Pose pose_with_cut_at_moving_suite = pose;
		Size const jump_at_moving_suite = make_cut_at_moving_suite( pose_with_cut_at_moving_suite, moving_suite );
		core::kinematics::FoldTree const & f( pose_with_cut_at_moving_suite.fold_tree() );

		Size const nres( pose.total_residue() );
		ObjexxFCL::FArray1D_bool partition_definition( nres, false );
		f.partition_by_jump( jump_at_moving_suite, partition_definition );
		job_parameters_->set_partition_definition( partition_definition ); //this is a useful decomposition.
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Figure out a good root residue -- which partition of the pose has the fewest residues to move around?
	void
	StepWiseRNA_PoseSetup::reroot_fold_tree( pose::Pose & pose ){

		ObjexxFCL::FArray1D< bool > const & partition_definition = job_parameters_->partition_definition();
		Size const nres = pose.total_residue();

		Size num_partition_0( 0 ), num_partition_1( 0 );
		Size some_residue_in_partition_0( 0 ), some_residue_in_partition_1( 0 ), root_res( 0 );

		for ( Size n = 1; n <= nres; n++ ) {
			if( partition_definition( n ) ) {
				num_partition_1 += 1;
				some_residue_in_partition_1 = n;
			} else {
				num_partition_0 += 1;
				some_residue_in_partition_0 = n;
			}
		}

		assert( num_partition_0 > 0 );
		assert( num_partition_1 > 0 );

		if ( num_partition_1 >= num_partition_0 ){
			// best to put the root in partition 1 -- it is bigger, and will stay anchored.
			if ( partition_definition( 1 ) ) {
				root_res = 1;
			} else if ( partition_definition( nres ) ) {
				root_res = nres;
			} else {
				root_res = some_residue_in_partition_1;
			}
		} else {
			if ( !partition_definition( 1 ) ) {
				root_res = 1;
			} else if ( !partition_definition( nres ) ) {
				root_res = nres;
			} else {
				root_res = some_residue_in_partition_0;
			}
		}

		Size const moving_res( job_parameters_->working_moving_res() );
		std::cout << "Num. res in partition 0:  " << num_partition_0 << ".   Num. res in partition 1: " << num_partition_1 <<
			". Moving_res is in  partition " << partition_definition( moving_res ) <<
			".  New root residue " << root_res << std::endl;

		assert( root_res > 0 );
		core::kinematics::FoldTree f = pose.fold_tree();
		f.reorder( root_res );
		pose.fold_tree( f );

		// moving positions
		utility::vector1< Size > moving_positions;
		bool const root_partition = partition_definition( pose.fold_tree().root() );
		for (Size seq_num=1; seq_num<=pose.total_residue(); seq_num++){
			if ( partition_definition( seq_num ) != root_partition ) 	moving_positions.push_back( seq_num );
		}
		job_parameters_->set_moving_pos( moving_positions );

		//////////////////////////////////////////////////////////////////////////////////////////////////
		// For "internal" suites, we can also better define whether we are prepending or appending.
		//  This actually does not affect very much -- only how we *cluster* (suite_rmsd
		//  needs to pick a side_chain to calculate rmsds over!).
		Size const moving_suite( job_parameters_->working_moving_suite() );
		if ( job_parameters_->Is_internal() ){
			if ( partition_definition( moving_suite ) == partition_definition( root_res ) ){
				job_parameters_->set_Is_prepend( false );
			} else {
				job_parameters_->set_Is_prepend( true );
			}
		} else {
			/// consistency check!
			bool const should_be_prepend = ( partition_definition( moving_suite ) != partition_definition( root_res ) );
			if ( should_be_prepend != job_parameters_->Is_prepend() ) {
				utility_exit_with_message( "Possible problem with prepend/append assignment!!" );
			}
		}

	}


	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::figure_out_gap_size_and_five_prime_chain_break_res(){

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries(	 job_parameters_->chain_boundaries() );

 		Size gap_size = 999; // junk value... totally "free" end.
 		job_parameters_->set_gap_size( gap_size /*DUMMY*/ );
 		job_parameters_->set_five_prime_chain_break_res( 0 );

		/////////////////////////////////////////////////////////////////////////////////////////////
		// Need to look for a chainbreak whose ends actually will move relative to each other if
		// we change degrees of freedom in the "moving residues".
		ObjexxFCL::FArray1D< bool > const & partition_definition = job_parameters_->partition_definition();

		/////////////////////////////////////////////////////////////////////////////////////////////
 		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
		Size found_moving_gap( 0 );
		Size const num_chains = chain_boundaries.size();
		for ( Size n = 1; n < num_chains; n++ ) {
			Size chain_end = chain_boundaries[ n ].second;
			Size next_chain_start = chain_boundaries[ n+1 ].first;
			if ( partition_definition( full_to_sub[ chain_end ] ) !=
					 partition_definition( full_to_sub[ next_chain_start ] ) ){

				bool found_cutpoint_open( false );
				for (Size i = 1; i <= cutpoint_open_.size(); i++ ){
					if ( cutpoint_open_[i] >= chain_end && cutpoint_open_[i] < next_chain_start ) {
						found_cutpoint_open = true;
						break;
					}
				}
				if ( found_cutpoint_open ) continue;

				job_parameters_->set_gap_size(  next_chain_start - chain_end - 1 );
				job_parameters_->set_five_prime_chain_break_res( full_to_sub[ chain_end ] );
				found_moving_gap++;

			}
		}

		if ( found_moving_gap > 1 ){
			utility_exit_with_message( "Had trouble figure out which gap might be the one to close! Try to renumber input poses sequentially.");
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::apply_cutpoint_variants( pose::Pose & pose ) const {

		using namespace core::id;

		pose::Pose pose_copy = pose;

		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );

		if ( cutpoint_closed_ == 0 ) return;

		if ( full_to_sub.find( cutpoint_closed_ )   != full_to_sub.end() &&
				 full_to_sub.find( cutpoint_closed_+1 ) != full_to_sub.end() ) {

			std::cout << "Applying cutpoint variants to " << cutpoint_closed_ << std::endl;

			Size const cutpos = full_to_sub[ cutpoint_closed_ ];

			// Taken from Parin's code. Need to make sure virtual atoms are correctly positioned
			// next to O1P, O2P.
			Correctly_position_cutpoint_phosphate_torsions( pose, cutpos, false /*verbose*/ );

			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpos   );
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpos+1 );


			for (Size i = cutpos; i <= cutpos + 1; i++ ){
				for (Size j = 1; j <= scoring::rna::NUM_RNA_MAINCHAIN_TORSIONS; j++ ) {
					id::TorsionID torsion_id( i, id::BB, j );
					pose.set_torsion( torsion_id, pose_copy.torsion( torsion_id ) ) ;
				} // j
			} // i


		} else {
			utility_exit_with_message( "User provided cutpoint_closed not in working pose?" );
		}

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::apply_bulge_variants( pose::Pose & pose ) const {
		using namespace core::conformation;
		using namespace core::pose;

		ObjexxFCL::FArray1D< Size > const & is_working_res( job_parameters_->is_working_res() );
		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );

		for ( Size n = 1; n <= bulge_res_.size(); n++ ) {
			if ( !is_working_res( bulge_res_[ n ] ) ) continue;
			core::pose::add_variant_type_to_pose_residue( pose, "BULGE", 	full_to_sub[ bulge_res_[ n ] ] );
		}


	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::check_close_chain_break( core::pose::Pose const & pose ) const{

		if ( !Is_close_chain_break( pose ) && job_parameters_->gap_size() == 0 ) {
			utility_exit_with_message( "mismatch --> gap_size = 0, but no cutpoint variants defined?" );
		}
		if ( Is_close_chain_break( pose ) && job_parameters_->gap_size() != 0 ) {
			utility_exit_with_message( "mismatch --> gap_size != 0, but cutpoint variants defined?" );
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::Import_pose( Size const & i, core::pose::Pose & import_pose) const {

		using namespace core::conformation;

		assert( i <= input_tags_.size() );

		if ( i > silent_files_in_.size() ) {
			// not a silent file, read in from pdb text file.

			std::string pose_name = input_tags_[ i ];
			std::size_t found=pose_name.find(".pdb");
			if (found==std::string::npos) {
				pose_name.append(".pdb");
			}

			//  if(verbose) std::cout << "	The following pose will be imported :" << pose_name << std::endl;
			core::import_pose::pose_from_pdb( import_pose, *rsd_set_, pose_name );

		} else {

			core::io::silent::SilentFileData silent_file_data;
			silent_file_data.read_file( silent_files_in_[ i ] );

			std::string const & input_tag = input_tags_[ i ];
			bool found_tag( false );
			for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(),
							end = silent_file_data.end(); iter != end; ++iter ) {
				if ( iter->decoy_tag() != input_tag ) continue;
				found_tag = true;
				iter->fill_pose( import_pose, *rsd_set_ );
				break;
			}

			if ( !found_tag ) utility_exit_with_message( "Could not find specified tag " + input_tag + " in silent file!" );

		}

		///////////////////////////////
		// Check for sequence match.
		utility::vector1< Size > const & input_res = input_res_vectors_[ i ];
		if ( import_pose.total_residue() != input_res.size() ) utility_exit_with_message( "input pose does not have same # residues as input res" );
		bool match( true );
		for( Size n = 1; n <= import_pose.total_residue(); n++ ) {
			if (  import_pose.sequence()[ n-1 ]  !=
						desired_sequence_[ input_res[n] - 1 ] ) {
				match = false; break;
			}
		}
		if (!match) utility_exit_with_message( "mismatch in sequence between input pose and desired sequence, given input_res " );
	}

	////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::slice_native(){

		using namespace core::conformation;
		using namespace core::pose;

		ObjexxFCL::FArray1D< Size > const & is_working_res( job_parameters_->is_working_res() );
		std::string const & working_sequence( job_parameters_->working_sequence() );

		if ( !get_native_pose() ) return;
		Pose const & native_pose( *get_native_pose() );

		PoseOP working_native_pose( new Pose );

		for ( Size i = 1; i <= native_pose.total_residue(); i++ ) {

			if ( !is_working_res( i ) ) continue;
			ResidueOP residue_to_add = native_pose.residue( i ).clone() ;
			working_native_pose->append_residue_by_bond(  *residue_to_add  ) ;

		}

		assert( working_native_pose->sequence() == working_sequence );

		job_parameters_->set_working_native_pose( working_native_pose );

		//		working_native_pose_->dump_pdb( "working1.pdb" );

	}

	//////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size > const
	StepWiseRNA_PoseSetup::apply_full_to_sub_mapping( utility::vector1< Size > & res_vector) const{

		using namespace core::conformation;
		using namespace core::pose;

		ObjexxFCL::FArray1D< Size > const & is_working_res( job_parameters_->is_working_res() );
		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );

		utility::vector1< core::Size > working_res_vector;
		for ( Size n = 1; n <= res_vector.size(); n++ ) {
			if ( !is_working_res( res_vector[ n ] ) ) continue;
			working_res_vector.push_back( full_to_sub[ res_vector[ n ] ]);
		}

		return working_res_vector;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	StepWiseRNA_JobParametersOP & StepWiseRNA_PoseSetup::job_parameters(){
		return job_parameters_;
	}

	/////////////////////////////////////////////////////////////////
	void StepWiseRNA_PoseSetup::figure_out_Prepend_Internal( pose::Pose const & pose ){

		Size const working_moving_res( job_parameters_->working_moving_res() );

		bool Is_prepend( false ), Is_internal( false );
		Size moving_suite( 0 );
		if ( working_moving_res == 1 ||
				 pose.fold_tree().is_cutpoint( working_moving_res - 1 ) ) {
			Is_prepend = true;
			Is_internal = false;
			moving_suite = working_moving_res;
		} else if ( pose.fold_tree().is_cutpoint( working_moving_res ) ){
			// append
			Is_prepend = false;
			Is_internal = false;
			moving_suite = working_moving_res - 1;
		} else {
			Is_prepend = false;
			Is_internal = true;
			moving_suite = working_moving_res;
		}

		job_parameters_->set_Is_prepend( Is_prepend );
		job_parameters_->set_Is_internal( Is_internal );
		job_parameters_->set_working_moving_suite( moving_suite );

	}

	//////////////////////////////////////////////////////////////////////////////
	// Assume we have done a crappy job of placing 5' phosphates.
	//  this would be true if, for example, the 5' phosphate was built
	//  by prepending in a previous rebuild-from-scratch effort.
	// The only exception is if the 5'-phosphate is involved in *chain closure*.
	void
	StepWiseRNA_PoseSetup::apply_virtual_phosphate_variants( pose::Pose & pose ) const{

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries(	 job_parameters_->chain_boundaries() );
		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );

		Size const num_chains = chain_boundaries.size();

		for ( Size n = 1; n <= num_chains; n++ ) {
			Size const chain_start = chain_boundaries[ n ].first;
			if ( cutpoint_closed_ > 0   &&  full_to_sub[ chain_start ] == full_to_sub[ cutpoint_closed_ ]+1 ) continue;

			if ( pose.residue( full_to_sub[ chain_start ] ).type().has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
				utility_exit_with_message( "Should not be trying to virtualize phosphate on close cutpoint residue!" );
			}

			core::pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", full_to_sub[ chain_start ] );
		}
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::set_fixed_res( utility::vector1 < core::Size > const & fixed_res ){
		fixed_res_ = fixed_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::set_terminal_res( utility::vector1 < core::Size > const & terminal_res ){
		terminal_res_ = terminal_res;
	}


  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_PoseSetup::set_bulge_res( utility::vector1 < core::Size > const & bulge_res ){
		bulge_res_ = bulge_res;
	}

}
}
}
