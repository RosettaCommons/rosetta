// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinPoseSetup
/// @brief Sets up pose and job parameters for protein or RNA stepwise building.
/// @detailed
/// @author Rhiju Das
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/stepwise/legacy/modeler/protein/StepWiseProteinPoseSetup.hh>
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/legacy/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/constraint_util.hh>

//////////////////////////////////
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
//////////////////////////////////////////////
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/import_pose/import_pose.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <basic/Tracer.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>

// RNA stuff
#include <protocols/farna/util.hh>
#include <core/chemical/rna/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <utility/stream_util.hh>
#include <utility/exit.hh>

#include <string>

//Auto Headers
#include <core/scoring/rms_util.hh>
#include <utility/vector1.hh>

using namespace core;
using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::modeler::protein;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise modeler of proteins (and probably other
// biopolymers soon). Take a starting pose and a list of residues to sample,
//  and comprehensively sample all backbone torsion angles by recursion.
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.stepwise.legacy.modeler.protein.StepWiseProteinPoseSetup" ) ;

//typedef std::map< core::Size, core::Size > ResMap;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace protein {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
	StepWiseProteinPoseSetup::StepWiseProteinPoseSetup( utility::vector1< core::Size > const & moving_res_list,
																											std::string const & desired_sequence,
																											utility::vector1< InputStreamWithResidueInfoOP > & input_streams_with_residue_info,
																											utility::vector1< core::Size > const & cutpoint_open,
																											utility::vector1< core::Size > const & cutpoint_closed ):
		// do we still need this variable, moving_res_, since we have moving_res_list? --- Rhiju, feb. 2010
		moving_res_list_( moving_res_list ),
		desired_sequence_( desired_sequence ), //This is the full/global sequence Jan 2, 2009 Parin S.
		rsd_set_( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ),
		input_streams_with_residue_info_( input_streams_with_residue_info ),
		cutpoint_open_( cutpoint_open ),
		cutpoint_closed_( cutpoint_closed ),
		is_cutpoint_( desired_sequence_.size(), false ),
		secstruct_( "" ),
		working_parameters_( new stepwise::modeler::working_parameters::StepWiseWorkingParameters ),
		virtualize_5prime_phosphates_( true ),
		add_peptide_plane_variants_( false ),
		remove_nterminus_variant_( false ),
		remove_cterminus_variant_( false ),
		parin_favorite_output_( false ),
		add_virt_res_(false),
		cst_file_( "" ),
		disulfide_file_( "" ),
		align_file_( "" ),
		ready_to_align_( false ),
		dump_( false )
	{
		///////////////////////////////////////////////////////
		// Cutpoint setup
		for ( Size n = 1; n <= cutpoint_closed_.size();   n++ ) {
			is_cutpoint_( cutpoint_closed_[ n ] ) = true;
		}

		for ( Size n = 1; n <= cutpoint_open_.size();   n++ ) {
			is_cutpoint_( cutpoint_open_[ n ] ) = true;
			for ( Size m = 1; m <= cutpoint_closed_.size();   m++ ) {
				if ( cutpoint_open_[ n ] == cutpoint_closed_[ m ] ) utility_exit_with_message( "Position cannot be both cutpoint_open and cutpoint_closed" );
			}
		}

	}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseProteinPoseSetup::~StepWiseProteinPoseSetup()
  {}

	/////////////////////
	std::string
	StepWiseProteinPoseSetup::get_name() const {
		return "StepWiseProteinPoseSetup";
	}


	//////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::apply( core::pose::Pose & pose ) {

		using namespace core::pose;
		using namespace core::kinematics;
		using namespace core::chemical;

		// Define chain boundaries, mapping to working pose, etc.
		figure_out_working_sequence_and_mapping();
		figure_out_jump_partners();
		figure_out_cuts();

		// actually make the pose, set fold tree, copy in starting
		//  templates from disk.
		make_pose( pose ); //Create pose with random torsions + setup the fold_tree

		initialize_pose_from_streams( pose ); //Grab info from input poses/silent-files. Also initializes full_to_sub in the streams.

		/////////////////////////////////////////////////
		// useful in job definition. Need to carry out
		// following in exactly the right order!
		figure_out_Prepend_Internal( pose );
		figure_out_partition_definition( pose ); // this DOES NOT re-roots the fold tree. Parin S. Jan 30, 2010
		figure_out_gap_size_and_first_chain_break_res();

		////////////////////////////////
		//Misc. parameters for job --
		// fixed res, terminal res will be useful in finalizing fold tree.
		// superimpose_res, rmsd_res will be necessary for clustering across subset of residues, and rmsds to native.
		working_parameters_->set_working_fixed_res(  apply_full_to_sub_mapping( fixed_res_ ) );
		working_parameters_->set_working_superimpose_res(  apply_full_to_sub_mapping( superimpose_res_ ) );
		working_parameters_->set_working_calc_rms_res(  apply_full_to_sub_mapping( calc_rms_res_ ) );
		working_parameters_->set_working_bridge_res(  apply_full_to_sub_mapping( bridge_res_ ) );

		/////////////////////////////////
		// More fold tree stuff.
		reroot_fold_tree( pose );
		apply_cutpoint_variants( pose );
		check_close_chain_break( pose );

		/////////////////////////////////
		// Residue variants.
		if ( virtualize_5prime_phosphates_ ) apply_virtual_phosphate_variants( pose );
		if ( add_peptide_plane_variants_ )   apply_peptide_plane_variants( pose );
		apply_bulge_variants( pose );
		apply_virtual_res_variant( pose);

		setup_constraints( pose );
		setup_disulfides( pose );
		if (add_virt_res_)
			add_aa_virt_rsd_as_root(pose);

		//////////////////////////////////////////////////
		// Final cleanup. This is manual and quite silly, but the pose machinery currently
		// does a poor job with phi/psi's at the ends of chains.
		// ACTUALLY IMPROVEMENTS TO COPY_DOFS() SEEMS TO FIX THIS ISSUE -- AND
		//  NOW FIX_PHI_PSI_OFFSETS() IS ACTUALLY CAUSING A PROBLEM!!
		//		if ( dump_ ) pose.dump_pdb( "before_fix_phi_psi.pdb" );
		//		fix_phi_psi_offsets( pose );
		//		if ( dump_ ) pose.dump_pdb( "after_fix_phi_psi.pdb" );

		//Would it be possible to setup the fold_tree, chain_break terminus and virtual phosphate of the native_pose as well so that it can be minimized? Parin Jan 29, 2010
		//This will be easy to do it we create a new class to create working_parameters_ independent of pose_setup
		//		check_superimpose_res( pose );

		setup_working_native_pose(); //slice up native pose
		align_poses( pose ); // read in align pose, then align pose and native_pose to it.
		working_parameters_->set_working_native_pose( working_native_pose );

		add_terminal_res_repulsion( pose );

		setup_secstruct( pose );

		// new way to pass information about what is fixed & minimized - -necessary for StepWiseMonteCarlo.
		setup_full_model_info( pose ); // fixed_domain [fixed_res], extra_minimize_res.

		//		pose.dump_pdb( "start.pdb" );
		//		TR.Debug << "FOLD_TREE " << pose.fold_tree() << std::endl;
		//		if ( pose.fold_tree().num_jump() > 0 ) TR.Debug << "RT " << pose.jump( 1 ).rt() << std::endl;
	}

	// silly util
	utility::vector1< Size >
	convert_to_vector1( ObjexxFCL::FArray1D< Size > farray ){
		utility::vector1< Size > vec1;
		for ( Size n = 1; n <= farray.size(); n++ ) vec1.push_back( farray(n) );
		return vec1;
	}

	/////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::figure_out_working_sequence_and_mapping(){

		utility::vector1< std::pair< core::Size, core::Size > > chain_boundaries;
		//core::Size num_chains;
		std::map< core::Size, core::Size > full_to_sub;
		std::string working_sequence;

		// There are up to two input poses. Need to merge their input residues,
		// and figure out chain boundaries.
		Size const nres( desired_sequence_.size() );
		ObjexxFCL::FArray1D< Size > is_working_res( nres );
		ObjexxFCL::FArray1D< Size > is_moving_res( nres );
		is_working_res = 0;
		is_moving_res = 0;

		// Then input residues. Note that this can overwrite domain definitions in loop
		for ( Size i = 1; i <= input_streams_with_residue_info_.size(); i++ ) {
			utility::vector1< Size > const & input_res_vector =  input_streams_with_residue_info_[i]->input_res();
			for ( Size n = 1; n <= input_res_vector.size(); n++ )  {
				is_working_res( input_res_vector[ n ] ) = i; //indicate whether belong to pose 1 or 2, what about the case where there are overlap residues between pose1 and pose2? Jan 29 2010, Parin S.
			}
		}

		//First any residues in a loop to be closed.
		for ( Size i = 1; i <= bridge_res_.size(); i++ ) {
			is_working_res( bridge_res_[ i ] ) = BRIDGE_RES; /*Some number*/
		}


		for( Size i = 1; i <= moving_res_list_.size(); i++){
			if ( !is_working_res( moving_res_list_[ i ] ) ) is_working_res( moving_res_list_[i] ) = MOVING_RES;
			is_moving_res( moving_res_list_[i] ) = MOVING_RES;
		}

		Size start_chain( 0 );
		Size end_chain( 0 );

		Size n( 0 );
		for ( Size pos = 1; pos <= nres; pos++ ) {

			if ( !is_working_res( pos ) ) continue;
			n++;

			if ( n == 1 ) start_chain = pos;

			if (n > 1 &&
					( pos > end_chain + 1 || //In what situation does pos > end_chain + 1??
						is_cutpoint_( end_chain ) ) ) {

				chain_boundaries.push_back( std::make_pair( start_chain, end_chain ) ); //The last chain...
				//				TR.Debug << "FOUND CHAIN " << start_chain << " " << end_chain << std::endl;
				//check_moving_res_in_chain( start_chain, end_chain, chain_boundaries.size(), which_chain_has_moving_res );

				start_chain = pos;
			}

			end_chain = pos;
		}

		// For now, need to have at least one chain defined in the input!
		assert( start_chain > 0 );
		chain_boundaries.push_back( std::make_pair( start_chain, end_chain ) );
		//		TR.Debug << "FOUND CHAIN " << start_chain << " " << end_chain << std::endl;
		//		check_moving_res_in_chain( start_chain, end_chain, chain_boundaries.size(), which_chain_has_moving_res );

		//num_chains = chain_boundaries.size();  // set but never used ~Labonte

		//		if  ( which_chain_has_moving_res== 0 ) {
		//			utility_exit_with_message( "Could not figure out which chain to attach moving_res to!" );
		//		}

		working_sequence= "";
		full_to_sub.clear();
		Size count( 0 );
		utility::vector1< core::Size > working_res_list;
		for ( Size i = 1; i <= desired_sequence_.size(); i++ ) {
			if ( is_working_res( i ) ) {
				working_sequence += desired_sequence_[ i-1 ]; //i-1 because std::string elements starts at 0...
				count++;
				full_to_sub[ i ] = count;
				working_res_list.push_back( i );
			}
		}

		TR.Debug << "desired_sequence= " << desired_sequence_ << std::endl;
		TR.Debug << "working_sequence= " << working_sequence << std::endl;

		utility::vector1< core::Size > working_moving_res_list;
		for(Size i=1; i <= desired_sequence_.size(); i++){
			if ( is_working_res( i ) && is_moving_res( i ) ){
				working_moving_res_list.push_back(full_to_sub[ i ]);
			}
		}

		working_parameters_->set_full_to_sub( full_to_sub );  //res_map
		working_parameters_->set_sequence( desired_sequence_ ); //full_sequence
		//		working_parameters_->set_working_sequence( working_sequence ); //partial_sequence, auto-determined based on full_to_sub
		//		working_parameters_->set_working_res_list( working_res_list ); // auto-determined from full_to_sub
		working_parameters_->set_working_moving_res_list(  working_moving_res_list );
		working_parameters_->set_is_working_res( convert_to_vector1( is_working_res ) ); //A vector to indicate if res exist in the working pose and if it belong to input_pose 1 or input_pose 2 or working_res_list.
		working_parameters_->set_is_moving_res( is_moving_res );
		working_parameters_->set_chain_boundaries( chain_boundaries );
		//working_parameters_->set_which_chain_has_moving_res( which_chain_has_moving_res );

	}


	///////////////////////////////////////////////////////////////////////////////////
	// Changed this to match Parin's favorite convention --
	// jumps far away from rebuilt residue.
	void
	StepWiseProteinPoseSetup::figure_out_jump_partners() {

		jump_partners_.clear();

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries(	 working_parameters_->chain_boundaries() );
		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );
		utility::vector1< Size > const & is_working_res = working_parameters_->is_working_res();

		Size const num_chains( chain_boundaries.size() );

		//		TR.Debug << "NUM_CHAINS: " << num_chains << std::endl;

		////////////////////////////////////////////////
		// Make a list of potential jumps.
		utility::vector1< std::pair< Size, Size > > potential_jump_partners, potential_jump_partners_low_priority;

		// start with user-defined jumps (top priority list)
		utility::vector1< Size > working_jump_res = working_parameters_->apply_full_to_sub_mapping( jump_res_ );
		for ( Size n = 1; n <= working_jump_res.size() /2; n++ ){
			Size const i = working_jump_res[ 2*(n-1) + 1 ];
			Size const j = working_jump_res[ 2*(n-1) + 2 ];
			potential_jump_partners.push_back( std::make_pair( i, j) );
		}

		// Then, jumps that go from the end of one chain to the beginning of the next --
		//  note that we don't like this if moving residues are involved, so those go
		//  to a third priority list.
		for ( Size n = 1 ; n <= num_chains; n++ ) {

			Size const i = chain_boundaries[n].second;
			Size const j = n < num_chains ? chain_boundaries[n+1].first : chain_boundaries[1].first;

			// If there are two input poses, make sure their
			// rigid body arrangment is sampled --> no jumps between these
			// different regions!
			TR.Debug << "CHAIN " << n << ". Considering jump:  " << i << "  to  " << j << ". WORKING_RES: " << is_working_res[i] << " " << is_working_res[j] << std::endl;
 			if ( is_working_res[ i ] == 1 &&
					 is_working_res[ j ] == 2 ) continue;
			if ( is_working_res[ j ] == 1 &&
					 is_working_res[ i ] == 2 ) continue;
			if ( is_working_res[ i ] == BRIDGE_RES ) continue;
			if ( is_working_res[ j ] == BRIDGE_RES ) continue;

			TR.Debug << "Adding potential jump! " << std::endl;

			std::pair< Size, Size > jump_pair = std::make_pair( full_to_sub[ i ], full_to_sub[ j ] );
			if ( !moving_res_list_.has_value( i) && !moving_res_list_.has_value( j) ){
				potential_jump_partners.push_back( jump_pair );
			} else {
				potential_jump_partners_low_priority.push_back( jump_pair );
			}
		}

		// Tack on the low priority jump pairs at the end.
		for ( Size n = 1; n <= potential_jump_partners_low_priority.size(); n++ ) {
			potential_jump_partners.push_back(  potential_jump_partners_low_priority[ n ]  );
		}

		// Might as well figure out (ahead of time) the chain assignments.
		utility::vector1< std::pair< Size, Size > > potential_chain_partners;
		for ( Size n = 1; n <= potential_jump_partners.size(); n++ ) {
			potential_chain_partners.push_back(  std::make_pair( which_chain( potential_jump_partners[n].first ),
																													 which_chain( potential_jump_partners[n].second ) ) );
			TR.Debug << "POTENTIAL JUMPS: " << potential_jump_partners[n].first << " "
								<< potential_jump_partners[n].second << "   chain: "
								<< which_chain( potential_jump_partners[n].first ) << ' '
								<< which_chain( potential_jump_partners[n].second ) << std::endl;
		}


		////////////////////////////////////////////////
		// Then grab enough jumps so that we have all the chains
		//  somehow connected to each other.
		jump_partners_.clear();
		utility::vector1< std::pair< Size, Size > > chain_partners;
		Size ntries( 0 );

		while ( chain_partners.size() < (num_chains-1) && ntries++ < 10000 ){
			for ( Size n = 1; n <= potential_jump_partners.size(); n++ ){
				TR.Debug << "Going to test jump: " << potential_jump_partners[ n ].first << " "
									<< potential_jump_partners[ n ].second << "? " << already_connected( potential_chain_partners[ n ], chain_partners ) << std::endl;

				if (  already_connected( potential_chain_partners[ n ], chain_partners ) ) continue;
				chain_partners.push_back( potential_chain_partners[ n ] );
				jump_partners_.push_back( potential_jump_partners[ n ] );
				TR.Debug << "ADDING JUMP: " << potential_jump_partners[n].first << ' ' << potential_jump_partners[n].second << std::endl;
				break;
			}
		}

		if ( chain_partners.size() != (num_chains-1) ) utility_exit_with_message( "Problem setting up jump partners!" );

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	StepWiseProteinPoseSetup::which_chain( Size const & i ){

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries(	 working_parameters_->chain_boundaries() );
		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );

		for ( Size n = 1; n <= chain_boundaries.size(); n++ ){
			if ( full_to_sub[ chain_boundaries[n].first ] <= i &&
					 full_to_sub[ chain_boundaries[n].second ] >= i ) return n;
		}
		return 0;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseProteinPoseSetup::already_connected( std::pair< Size, Size > const & potential_chain_partner,
																				utility::vector1< std::pair< Size, Size > > const & chain_partners ) const {
		utility::vector1< bool > already_checked( chain_partners.size(), false );
		// traverse our way from chain to chain until we get from the first partner to the second partner -- or to the end of the list.
		return already_connected(  potential_chain_partner.first, potential_chain_partner.second, chain_partners, already_checked );
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseProteinPoseSetup::already_connected( Size const start_chain,
																				Size const stop_chain,
																				utility::vector1< std::pair< Size, Size > > const & chain_partners,
																				utility::vector1< bool > already_checked ) const {

		bool made_the_connection( false );
		for ( Size n = 1; n <= chain_partners.size(); n++ ) {

			if ( already_checked[ n ] ) continue;

			Size other_chain( 0 );
			if ( chain_partners[n].first == start_chain ){
				other_chain = chain_partners[ n ].second;
			}
			if ( chain_partners[n].second == start_chain ){
				other_chain = chain_partners[ n ].first;
			}
			if ( other_chain < 1 ) continue;

			if ( other_chain == stop_chain ){
				made_the_connection = true;
			} else {
				utility::vector1< bool > already_checked_new = already_checked;
				already_checked_new[ n ] = true;
				made_the_connection = already_connected( other_chain, stop_chain, chain_partners, already_checked_new );
			}

			if ( made_the_connection ) return true;
		}

		// end of the road
		return false;

	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::figure_out_cuts() {

		cuts_.clear();

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries(	 working_parameters_->chain_boundaries() );
		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );
		Size const num_chains = chain_boundaries.size();

		for ( Size n = 1; n < num_chains; n++ ) {
			cuts_.push_back( full_to_sub[ chain_boundaries[n].second ] );
		}

	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::make_full_pose( pose::Pose & pose ){

		make_pose_from_sequence( pose, desired_sequence_, *rsd_set_); //, false /*auto_termini*/);

		if ( remove_nterminus_variant_ ) 	pose::remove_lower_terminus_type_from_pose_residue( pose, 1 );
		if ( remove_cterminus_variant_ ) 	pose::remove_upper_terminus_type_from_pose_residue( pose, pose.total_residue() );

		TR.Debug << remove_nterminus_variant_ << ' ' << remove_cterminus_variant_ << "  FULL POSE ANNOTATED SEQUENCE: " << pose.annotated_sequence( true ) << std::endl;

	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::make_pose( pose::Pose & pose ){

		using namespace core::conformation;
		using namespace core::scoring::constraints;

		//std::string const working_sequence = working_parameters_->working_sequence();
		//		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );

		make_full_pose( pose );

		// Extend pose. (protein)
		for ( Size n = 1; n <= pose.total_residue();n++ ) {
			if ( !pose.residue( n ).is_protein() ) continue;
			pose.set_phi( n, -150 );
			pose.set_psi( n, 150);
			pose.set_omega( n, 180 );
		}

		// slice it down.
		utility::vector1< Size > const & is_working_res = working_parameters_->is_working_res();
		utility::vector1< Size > working_res_list;
		for ( Size n = 1; n <= pose.total_residue(); n++ ) if ( is_working_res[ n ] ) working_res_list.push_back( n );
		pdbslice( pose, working_res_list );


		Size const nres( pose.total_residue() );
		core::kinematics::FoldTree f( nres );
		assert( cuts_.size() == jump_partners_.size() );
		Size const num_cuts( cuts_.size() );

		ObjexxFCL::FArray2D< int > jump_point( 2, num_cuts, 0 );
		ObjexxFCL::FArray1D< int > cuts( num_cuts, 0 );

		for ( Size i = 1; i <= num_cuts; i++ ) {
			TR.Debug << "Upstream   jump_point= " << jump_partners_[i].first;
			TR.Debug << "  Downstream jump_point= " << jump_partners_[i].second;
			TR.Debug << "  cut_point= " << cuts_[ i ] << std::endl;

			jump_point( 1, i ) = jump_partners_[i].first;
			jump_point( 2, i ) = jump_partners_[i].second;
			cuts( i ) = cuts_[ i ];
			//			TR.Debug << " HELLO: " << jump_point( 1, i ) << " " << jump_point( 2, i ) << " " << cuts( i ) << std::endl;
		}

		//		if ( moving_res_list_.size() < 1 ) return;
		//		Size const root_res = ( full_to_sub[ moving_res_list_[1] ] == 1 ) ? nres : 1;

		f.tree_from_jumps_and_cuts( nres, num_cuts, jump_point, cuts ); //order of element in jump_point and cuts does not have to match. Jan 29, 2010 Parin S.

		for ( Size i = 1; i <= num_cuts; i++ ) {
			Size const k = f.upstream_jump_residue( i );
			Size const m = f.downstream_jump_residue( i );

			Residue const & rsd1( pose.residue( k ) );
			Residue const & rsd2( pose.residue( m ) );

			bool const KeepStubInResidue( true ); //rsd1.is_protein() & rsd2.is_protein() );
			f.set_jump_atoms( i, get_stepwise_jump_atom( rsd1 ), get_stepwise_jump_atom( rsd2 ),  KeepStubInResidue );

		}

		f.reassign_atoms_for_intra_residue_stubs(); // it seems silly that we need to do this separately.

		pose.fold_tree( f );

	}


	////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::setup_disulfides( pose::Pose & pose ){

		if ( disulfide_file_.size() == 0 ) return;

		Pose full_pose;
		make_full_pose( full_pose );

		// move to its own function?
		utility::vector1< std::pair<Size,Size> > disulf_bonds;
		core::io::raw_data::DisulfideFile ds_file( disulfide_file_ );
		ds_file.disulfides(disulf_bonds, full_pose);

		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );
		utility::vector1< Size > const & is_working_res = working_parameters_->is_working_res();

		utility::vector1< std::pair<Size,Size> > working_disulf_bonds;

		for ( Size n = 1; n <= disulf_bonds.size(); n++ ){
			//			TR.Debug <<  disulf_bonds[n].first  << " " << is_working_res[ disulf_bonds[n].first  ] << "    " << disulf_bonds[n].second << " " << is_working_res[ disulf_bonds[n].second  ] << std::endl;
			if ( is_working_res[ disulf_bonds[n].first  ]>0 &&
					 is_working_res[ disulf_bonds[n].second  ]>0  ){

				TR.Debug << "FOUND PAIR: " << disulf_bonds[n].first << "--" << disulf_bonds[n].second << "[in subpose: " << full_to_sub[ disulf_bonds[n].first ] << "--" << full_to_sub[ disulf_bonds[n].second ] << "]" << std::endl;

				working_disulf_bonds.push_back(  std::make_pair( full_to_sub[ disulf_bonds[n].first ], full_to_sub[ disulf_bonds[n].second ] ) );
			}
		}

		pose.conformation().fix_disulfides( working_disulf_bonds );
		for ( Size n = 1; n<= pose.total_residue();n++ ) TR.Debug << pose.residue_type( n ).has_variant_type( chemical::DISULFIDE );
		TR.Debug << std::endl;

		// this is not showing up right...
		TR.Debug << "AFTER DISULF: " << pose.annotated_sequence( true ) << std::endl;

	}

	////////////////////////////////////////////////////////////////////////////////////
	std::string
	StepWiseProteinPoseSetup::get_stepwise_jump_atom( core::conformation::Residue const & rsd ){

		//if (rsd.is_RNA() ) return rsd.atom_name( rsd.chi_atoms(1)[4] );
		if (rsd.is_RNA() ) return " C4'";

		if (rsd.is_protein() ) return " CA ";

		utility_exit_with_message( "Cannot deal with pose that is not RNA or protein yet." );
		return "";
	}

	////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::setup_constraints( pose::Pose & pose ){

		using namespace core::scoring::constraints;

		if ( cst_file_.size() == 0 ) return;

		utility::vector1< Size > const & working_res_list = working_parameters_->working_res_list();

		Pose full_pose;
		make_full_pose( full_pose );

		// Constraints!
		// To read in constraints, need to use full desired_sequence pose,
		// since there's a careful check of atom names...
		// Hey, this has to happen after the variants are figured out!!!
		cst_set_ = ConstraintIO::get_instance()->read_constraints( cst_file_, new ConstraintSet, full_pose );
		core::scoring::ScoreFunctionOP scorefxn = new core::scoring::ScoreFunction;
		scorefxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
		scorefxn->set_weight( core::scoring::coordinate_constraint, 1.0 );
		full_pose.constraint_set( cst_set_ );
		(*scorefxn)( full_pose );

		cst_set_ = constraint_set_slice( cst_set_, working_res_list, pose, full_pose );
		pose.constraint_set( cst_set_ );

		//(*scorefxn)( pose );
	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::initialize_phi_psi_offsets(  pose::Pose const & pose ){

		Size const & nres = pose.total_residue();

		phi_offsets_.dimension( nres, 0.0 );
		psi_offsets_.dimension( nres, 0.0 );

		for ( Size n = 1; n <= nres; n++ ) {
			if (!pose.residue( n ).is_protein()) continue;
			// The "pretend" phi and psi are torsion angles based on:
      //     C-CA-N-H (rather than C-CA-N-C)   for phi, and
      //     CA-N-C-O (rather than CA-N-C-CA)  for psi.
			phi_offsets_( n )  = protocols::stepwise::legacy::modeler::protein::get_pretend_phi_explicit( pose, n );
			psi_offsets_( n )  = protocols::stepwise::legacy::modeler::protein::get_pretend_psi_explicit( pose, n );
			//			TR.Debug << " HELLO! " << protocols::stepwise::legacy::modeler::protein::get_pretend_phi_explicit( pose, n )<< " " << phi_offsets_( n ) << " " << psi_offsets_( n ) << std::endl;
		}
	}


	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::save_phi_psi_offsets(
      pose::Pose const & start_pose,
			utility::vector1< core::Size > const & input_res,
			utility::vector1< core::Size > const & slice_res ){

		assert( slice_res.size() == input_res.size() );
		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );

		for ( Size n = 1; n <= input_res.size(); n++ ) {
			if (!start_pose.residue( slice_res[ n ] ).is_protein()) continue;
			phi_offsets_( full_to_sub[ input_res[ n ] ] )  = protocols::stepwise::legacy::modeler::protein::get_pretend_phi_explicit( start_pose, slice_res[ n ] ) ;
			psi_offsets_( full_to_sub[ input_res[ n ] ] )  = protocols::stepwise::legacy::modeler::protein::get_pretend_psi_explicit( start_pose, slice_res[ n ] ) ;
		}

	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::fix_phi_psi_offsets( pose::Pose & pose ) const{

		for ( Size n = 1; n <= pose.total_residue(); n++ ) {

			if (!pose.residue( n ).is_protein()) continue;
			Real const phi_current = protocols::stepwise::legacy::modeler::protein::get_pretend_phi_explicit( pose, n );
			Real const psi_current = protocols::stepwise::legacy::modeler::protein::get_pretend_psi_explicit( pose, n );

			pose.set_phi( n, pose.phi( n ) - phi_current + phi_offsets_( n ) ) ;
			pose.set_psi( n, pose.psi( n ) - psi_current + psi_offsets_( n ) ) ;
		}

	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::copy_rna_chi( pose::Pose & pose,
																	 pose::Pose const import_pose,
																	 utility::vector1< core::Size > const & input_res,
																	 utility::vector1< core::Size > const & slice_res ){

		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );

		for ( Size n = 1; n <= input_res.size(); n++ ) {
			if (!import_pose.residue( slice_res[ n ] ).is_RNA()) continue;
			pose.set_chi( 1, full_to_sub[ input_res[ n ] ],  import_pose.chi( 1, slice_res[n] ) );
		}

	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::initialize_pose_from_streams( pose::Pose & pose )
	{

		using namespace core::chemical;

		initialize_phi_psi_offsets( pose );

		Pose import_pose;

		if (dump_) pose.dump_pdb( "before_copy_dofs.pdb" );

 		// go ahead and align pose to this first input. If native or align_pose is specified, we will realign again later...
 		bool align_pose_to_import_pose( true );
 		if (get_native_pose() || align_file_.size() > 0 ) align_pose_to_import_pose = false; // will happen later in align_poses()

		for ( Size i = 1; i <= input_streams_with_residue_info_.size(); i++ ){

 			if ( i > 1) align_pose_to_import_pose = false;

			InputStreamWithResidueInfoOP & stream = input_streams_with_residue_info_[ i ];
			stream->set_full_to_sub( working_parameters_->full_to_sub() );
			stream->set_rsd_set( rsd_set_ );

 			stream->copy_next_pose_segment( pose, import_pose, true /*check_sequence_matches*/, align_pose_to_import_pose );

			save_phi_psi_offsets( import_pose,  stream->input_res(), stream->slice_res() );
			copy_rna_chi( pose, import_pose,  stream->input_res(), stream->slice_res() ); //hope this is OK.

			if (dump_) pose.dump_pdb( "after_copy_dofs"+ObjexxFCL::string_of(i)+".pdb" );
			if (dump_) import_pose.dump_pdb( "import"+ObjexxFCL::string_of(i)+".pdb" );

		}

		//		protocols::farna::print_internal_coords( pose );

	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::figure_out_partition_definition( pose::Pose const & pose ){

		ObjexxFCL::FArray1D_bool partition_definition( pose.total_residue(), false );

		utility::vector1< Size > const & moving_suite_list( working_parameters_->working_moving_suite_list() );
		utility::vector1< Size > const & moving_res_list( working_parameters_->working_moving_res_list() );

		// THIS SUCKS!  SHOULD GET RID OF IT AND SPECIFY FROM COMMAND LINE MOVING_SUITE!!!

		if ( moving_res_list.size() > 0 ){ //may not be filled if we are setting up for a "prepack"
			// trick to figure out which residues are upstream vs. downstream of the moving suite --
			Size partition_res( 0 );
			if ( (moving_suite_list.size() > 0 &&
						moving_suite_list[ 1 ] == (moving_res_list[ 1 ]-1) ) /*append*/ ||
					 ( moving_res_list[1] > 1 &&
						 moving_res_list[ moving_res_list.size() ] == pose.total_residue() ) ){
				partition_res = moving_res_list[ 1 ] - 1;
			} else {
				partition_res = moving_res_list[ moving_res_list.size() ];
			}
			TR.Debug << "PARTITION_RES ==> " << partition_res << std::endl;
			pose.fold_tree().partition_by_residue( partition_res, partition_definition );
		}

		TR.Debug << "PARTITION_DEFINITION ==> ";
		for ( Size i = 1; i <= pose.total_residue(); i++ ) TR.Debug << partition_definition( i );
		TR.Debug << std::endl;
		TR.Debug << "FOLD_TR.DebugEE " << pose.fold_tree() << std::endl;
		//		TR.Debug << "NUM_JUMPS " << pose.fold_tree().num_jump() << std::endl;

		working_parameters_->set_partition_definition( partition_definition ); //this is a useful decomposition.

		if ( moving_res_list.size() == 0 ) return;

		bool const start_partition = partition_definition( moving_res_list[ 1 ] );
		working_parameters_->set_working_moving_partition_res( figure_out_moving_partition_res( pose, moving_res_list ) );

		ObjexxFCL::FArray1D< Size > const & is_moving_res( working_parameters_->is_moving_res() );
 		std::map< core::Size, core::Size > & sub_to_full( working_parameters_->sub_to_full() );
		for ( Size i = 1; i <= pose.total_residue(); i++ ){
			if ( is_moving_res( sub_to_full[ i ] ) ) continue;
			if ( partition_definition( i )  == start_partition  ){
				working_parameters_->set_is_internal( true );
				break;
			}
		}


	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Figure out a good root residue -- which partition of the pose has the fewest residues to move around?
	void
	StepWiseProteinPoseSetup::reroot_fold_tree( pose::Pose & pose ){

		utility::vector1< Size > const & moving_res_list = working_parameters_->working_moving_res_list();
		if ( moving_res_list.size() < 1 ) return; // might be setting up for a prepack -- root is irrelevant.

		ObjexxFCL::FArray1D< bool > const & partition_definition = working_parameters_->partition_definition();
		// 		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );
		Size const nres = pose.total_residue();

		Size num_partition_0( 0 ), num_partition_1( 0 );
		Size possible_new_root_residue_in_partition_0( 0 ), possible_new_root_residue_in_partition_1( 0 ), root_res( 0 );

		for ( Size n = 1; n <= nres; n++ ) {
			if( partition_definition( n ) ) {
				num_partition_1 += 1;
				if ( pose.fold_tree().possible_root( n ) ) possible_new_root_residue_in_partition_1 = n;
			} else {
				num_partition_0 += 1;
				if ( pose.fold_tree().possible_root( n ) ) possible_new_root_residue_in_partition_0 = n;
			}
		}

		//If moving_res=1 or nres, then it is best to put root_res at one of the fixed residue in the correct partition. Parin Jan 18, 2010
		//CHANGE FROM fix_res to terminal_res!!! Parin Feb 7, 2010
		utility::vector1< core::Size > working_fixed_res( working_parameters_->working_fixed_res() );
		for(Size i=1; i<=working_fixed_res.size(); i++){
			Size const seq_num = working_fixed_res[i];
			if(partition_definition( seq_num ) &&  pose.fold_tree().possible_root( seq_num ) ){
				possible_new_root_residue_in_partition_1= seq_num;
				break;
			}
		}

		for(Size i=1; i<=working_fixed_res.size(); i++){
			Size const seq_num = working_fixed_res[i];
			if( !partition_definition( seq_num ) && pose.fold_tree().possible_root( seq_num ) ){
				possible_new_root_residue_in_partition_0= seq_num;
				break;
			}
		}

		//assert( num_partition_0 > 0 );
		//assert( num_partition_1 > 0 );

		//Made some changes here .... Parin Jan 18, 2009
		Size const moving_res( moving_res_list[1] );
		if ( working_parameters_->is_internal() ){
			//			TR.Debug << "is_internal" << std::endl;
			if ( num_partition_1 >= num_partition_0 ){
				// best to put the root in partition 1 -- it is bigger, and will stay anchored.
				if ( partition_definition( 1 ) && partition_definition( nres ) ) {
					root_res = 1;
				} else {
					root_res = possible_new_root_residue_in_partition_1;
				}
			} else {
				if ( !partition_definition( 1 ) && !partition_definition( nres )) {
					root_res = 1;
				} else {
					root_res = possible_new_root_residue_in_partition_0;
				}
			}
		} else {
			//			TR.Debug << "Is not internal" << std::endl;
			// Put root in the partition that does *not* contain the moving residue.
			if ( partition_definition( moving_res ) ){
				root_res = possible_new_root_residue_in_partition_0;
				//				TR.Debug << "part1 " << root_res << std::endl;
			} else {
				root_res = possible_new_root_residue_in_partition_1;
				// 			TR.Debug << "part0 " << root_res << std::endl;
			}
		}

		if ( root_res == 0 ) root_res = 1; // happens if all the residues are moving.

		TR.Debug << "Num. res in partition 0:  " << num_partition_0 << ".   Num. res in partition 1: " << num_partition_1 <<
			". Moving_res is in  partition " << partition_definition( moving_res ) <<
			".  New root residue " << root_res << std::endl;

		assert( root_res > 0 );
		core::kinematics::FoldTree f = pose.fold_tree();
		bool reorder_went_OK = f.reorder( root_res );
		runtime_assert( reorder_went_OK );
		pose.fold_tree( f );


		// moving positions // deprecated.
		//utility::vector1< Size > moving_positions;
		//		bool const root_partition = partition_definition( pose.fold_tree().root() );
		//		for (Size seq_num=1; seq_num<=pose.total_residue(); seq_num++){
		//			if ( partition_definition( seq_num ) != root_partition ) 	moving_positions.push_back( seq_num );
		//		}
		//		working_parameters_->set_moving_pos( moving_positions ); // deprecated.

		//////////////////////////////////////////////////////////////////////////////////////////////////
		// For "internal" suites, we can also better define whether we are prepending or appending.
		//  This actually does not affect very much -- only how we *cluster* (suite_rmsd
		//  needs to pick a side_chain to calculate rmsds over!). This way the suite_rmsd is calculated over the base in the smaller/moving paritition.

		utility::vector1< Size > working_moving_suite_list( working_parameters_->working_moving_suite_list() );
		if ( working_moving_suite_list.size() < 1 ) return;
		Size const first_moving_suite = working_moving_suite_list[ 1 ];

		if ( working_parameters_->is_internal() ){
			if ( partition_definition( first_moving_suite ) == partition_definition( root_res ) ){
				working_parameters_->set_is_prepend( false );
			} else {
				working_parameters_->set_is_prepend( true );
			}
		} else {
			/// consistency check!
			bool const should_be_prepend = ( partition_definition( first_moving_suite ) != partition_definition( root_res ) );
			if ( should_be_prepend != working_parameters_->is_prepend() ) {

				// Wait there's one way the root res is a moving res -- if we're building from scratch.
				utility::vector1< Size > const & is_working_res( working_parameters_->is_working_res() );
				ObjexxFCL::FArray1D< Size > const & is_moving_res( working_parameters_->is_moving_res() );
				bool ok( true );
				for ( Size i = 1; i <= desired_sequence_.size(); i++ ) {
					if ( is_working_res[ i ] && !is_moving_res( i ) ) {
						ok = false; break;
					}
				}

				if (!ok) {
					TR.Debug << " is_prepend: " << working_parameters_->is_prepend() << std::endl;
					TR.Debug << " should_be_prepend: " << should_be_prepend << std::endl;
					TR.Debug << " first_moving_suite: " << first_moving_suite << "  partition: " << partition_definition( first_moving_suite ) << std::endl;
					TR.Debug << " root_res: " << root_res << "   partition: " << partition_definition( root_res ) << std::endl;
					utility_exit_with_message( "Possible problem with prepend/append assignment!!" );
				}

			}
		}

	}


	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::figure_out_gap_size_and_first_chain_break_res(){

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries(	 working_parameters_->chain_boundaries() );

 		Size gap_size = 999; // junk value... totally "free" end.
 		working_parameters_->set_gap_size( gap_size /*DUMMY*/ );

		/////////////////////////////////////////////////////////////////////////////////////////////
		// Need to look for a chainbreak whose ends actually will move relative to each other if
		// we change degrees of freedom in the "moving residues".
		ObjexxFCL::FArray1D< bool > const & partition_definition = working_parameters_->partition_definition();
		utility::vector1< Size > const & is_working_res =  working_parameters_->is_working_res();

		/////////////////////////////////////////////////////////////////////////////////////////////
 		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );
		Size found_moving_gap( 0 );
		Size const num_chains = chain_boundaries.size();
		for ( Size n = 1; n < num_chains; n++ ) {
			Size chain_end = chain_boundaries[ n ].second;
			Size next_chain_start = chain_boundaries[ n+1 ].first;

			TR.Debug << "CHECK_CHAIN --> "  << chain_end << ' ' << partition_definition( full_to_sub[ chain_end ] ) << ' ' <<  next_chain_start << ' ' << partition_definition( full_to_sub[ next_chain_start ] ) << std::endl;

			if ( partition_definition( full_to_sub[ chain_end ] ) !=
					 partition_definition( full_to_sub[ next_chain_start ] ) ||
					 ( is_working_res[ chain_end ] == BRIDGE_RES ) ||
					 ( is_working_res[ next_chain_start ] == BRIDGE_RES ) ){

				bool found_cutpoint_open( false );
				for (Size i = 1; i <= cutpoint_open_.size(); i++ ){
					if ( cutpoint_open_[i] >= chain_end && cutpoint_open_[i] < next_chain_start ) {
						found_cutpoint_open = true;
						break;
					}
				}
				if ( found_cutpoint_open ) continue;

				TR.Debug << "CUTPOINT CLOSED --> "  << chain_end << ' ' << next_chain_start << std::endl;
				working_parameters_->set_gap_size(  next_chain_start - chain_end - 1 );
				found_moving_gap++;

			}
		}

		if ( found_moving_gap > 1 ){
			TR.Debug <<  "WARNING: Had trouble figure out which gap might be the one to close! Try to renumber input poses sequentially." << std::endl;
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::apply_cutpoint_variants( pose::Pose & pose ) const {

		using namespace core::id;

		pose::Pose pose_copy = pose;
		kinematics::FoldTree f_simple( pose.total_residue() );
		pose_copy.fold_tree( f_simple );

		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );

		for ( Size n = 1; n <= cutpoint_closed_.size(); n++ ){

			Size const cutpoint = cutpoint_closed_[ n ];

			if ( full_to_sub.find( cutpoint )   != full_to_sub.end() &&
					 full_to_sub.find( cutpoint+1 ) != full_to_sub.end() ) {

				Size const cutpos = full_to_sub[ cutpoint ];

				pose::correctly_add_cutpoint_variants( pose, cutpos ); //perhaps should not be in rna namespace...
				TR.Debug << "Applied cutpoint variants to " << cutpoint << std::endl;

				for (Size i = cutpos; i <= cutpos + 1; i++ ){
					utility::vector1< Real > const mainchain_torsions = pose_copy.residue( i ).mainchain_torsions();
					for (Size j = 1; j <= mainchain_torsions.size(); j++ ) {
						id::TorsionID torsion_id( i, id::BB, j );
						//TR.Debug << "TORSION AT CHAINBREAK: " << torsion_id << mainchain_torsions[ j ] << std::endl;
						pose.set_torsion( torsion_id, mainchain_torsions[ j ] ); //This makes sure that the chain_break torsions have the correct value
					} // j
				} // i

			}
		}

		//		pose.dump_pdb( "AFTER_CUTPOINTS.pdb" );

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::apply_bulge_variants( pose::Pose & pose ) const {
		using namespace core::conformation;
		using namespace core::pose;

		utility::vector1< Size > const & is_working_res( working_parameters_->is_working_res() );
		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );

		for ( Size n = 1; n <= bulge_res_.size(); ++n ) {
			if ( !is_working_res[ bulge_res_[ n ] ] ) continue;
			pose::add_variant_type_to_pose_residue( pose, core::chemical::BULGE, full_to_sub[ bulge_res_[ n ] ] );
		}


	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//applies terminus variants to any internal residues. Note that
	void
	StepWiseProteinPoseSetup::apply_terminus_variants_at_protein_rna_boundaries( pose::Pose & pose ) const {
		using namespace core::conformation;
		using namespace core::pose;

		for ( Size n = 1; n <= pose.total_residue()-1; n++ ) {
			if ( (pose.residue( n ).is_RNA() && pose.residue( n+1).is_protein())  ||
					 (pose.residue( n ).is_protein() && pose.residue( n+1).is_RNA()) ){
				pose::add_variant_type_to_pose_residue( pose, core::chemical::UPPER_TERMINUS_VARIANT, n  );
				pose::add_variant_type_to_pose_residue( pose, core::chemical::LOWER_TERMINUS_VARIANT, n+1 );
			}
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::check_close_chain_break( core::pose::Pose const & pose ) const{

		if ( !is_close_chain_break( pose ) && working_parameters_->gap_size() == 0 ) {
			utility_exit_with_message( "mismatch --> gap_size = 0, but no cutpoint variants defined?" );
		}
		//		if ( is_close_chain_break( pose ) && working_parameters_->gap_size() != 0 ) {
			//utility_exit_with_message( "mismatch --> gap_size != 0, but cutpoint variants defined?" );
		//		}

	}


	////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::check_superimpose_res( core::pose::Pose const & pose ){

		utility::vector1< core::Size > const & working_superimpose_res( working_parameters_->working_superimpose_res() );
		utility::vector1< core::Size > const & working_fixed_res( working_parameters_->working_fixed_res() );
		ObjexxFCL::FArray1D< bool > const & partition_definition = working_parameters_->partition_definition();

		bool const root_partition = partition_definition( pose.fold_tree().root() );

		for ( Size i = 1; i <= working_superimpose_res.size(); i++ ) {
			if ( !working_fixed_res.has_value( working_superimpose_res[ i ]) ){
				TR.Debug <<  "working superimpose_res " << working_superimpose_res[ i] << " not in fixed_res? " << std::endl;
				utility_exit_with_message( "You will get weird results if superimpose_res residues are allowed to move!" );
			}
			if ( !partition_definition( working_superimpose_res[ i ] ) == root_partition ) {
				TR.Debug << "working_superimpose_res " << working_superimpose_res[ i ] << std::endl;
				utility_exit_with_message( "Superimpose_res residues need to end up in root paritition and stay fixed!" );
			}
		}
	}


	////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::get_working_pose( pose::Pose const & pose, pose::Pose & working_pose ){

		using namespace core::pose;
		using namespace core::conformation;

		working_pose = Pose(); //clear this out.

		utility::vector1< Size > const & is_working_res( working_parameters_->is_working_res() );

		// Need a simple fold tree for following to work...
		Pose full_pose_copy = pose;
		full_pose_copy.fold_tree(  core::kinematics::FoldTree(  full_pose_copy.total_residue() ) );

		for ( Size i = 1; i <= full_pose_copy.total_residue(); i++ ) {
			if ( !is_working_res[ i ] ) continue;
			//TR.Debug << "About to append: " << i << std::endl;
			ResidueOP residue_to_add = full_pose_copy.residue( i ).clone() ;

			if ( i > 1 && residue_to_add->is_lower_terminus() ) {
				working_pose.append_residue_by_jump(  *residue_to_add, working_pose.total_residue()  ) ;
			} else {
				working_pose.append_residue_by_bond(  *residue_to_add  ) ;
			}
		}

		assert( working_pose.sequence() == working_parameters_->working_sequence() );

		// RNA thing.
		protocols::farna::make_phosphate_nomenclature_matches_mini( working_pose );


		//also carry over disulfides?


	}

	////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::setup_working_native_pose(){
		using namespace core::conformation;
		using namespace core::pose;
		using namespace protocols::farna;

		if ( get_native_pose() == 0 ) return;

		TR.Debug << "NATIVE sequence: " << get_native_pose()->annotated_sequence( true ) << std::endl;

		working_native_pose = new Pose;
		get_working_pose( *get_native_pose(), *working_native_pose );
		if (add_virt_res_)	add_aa_virt_rsd_as_root(*working_native_pose);

		setup_full_model_info( *working_native_pose ); // fixed_domain [fixed_res], extra_minimize_res.

	}

	//////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::align_poses( pose::Pose & pose ){

		if ( align_file_.size() > 0 ){
			Pose full_align_pose;
			import_pose::pose_from_pdb( full_align_pose, *rsd_set_, align_file_ );

			working_align_pose_ = new pose::Pose;
			get_working_pose( full_align_pose, *working_align_pose_ );

			TR.Debug << "PoseSetup: Align to file: " << align_file_ << std::endl;
		} else {
			if ( get_native_pose() ) {

				working_align_pose_ = new pose::Pose;
				*working_align_pose_ = *working_native_pose;

				TR.Debug << "PoseSetup: Align to NATIVE" << std::endl;

			} else {

				TR.Debug << "PoseSetup: Do not do any alignment." << std::endl;
				return;

			}
		}

		////////////////////////////////////////////////////////////////////////////////////////////////
		// Following is needed to help calculate rmsds to native... there are some runs,
		// for example where we just model a little loop. We only want to calculate rmsd over that loop, with
		// the rest of the structure pre-aligned. User can specify "-superimpose_res" from the command line.
		utility::vector1< core::Size > superimpose_res( working_parameters_->working_superimpose_res() );

		if ( superimpose_res.size() == 0 ){
			for ( Size i = 1; i <= pose.total_residue(); i++ ) superimpose_res.push_back( i );
		}

		alignment_atom_id_map_ = align::create_aligment_id_map_legacy( pose, *working_align_pose_, superimpose_res );
		ready_to_align_ = true;
		align_pose( pose );

		if ( get_native_pose() ) {
			id::AtomID_Map< id::AtomID > const & alignment_atom_id_map_native =
				align::create_aligment_id_map_legacy( *working_native_pose, *working_align_pose_, superimpose_res );
			core::scoring::superimpose_pose( *working_native_pose, *working_align_pose_, alignment_atom_id_map_native );
		}

		//		working_native_pose->dump_pdb( "WORKING_NATIVE.pdb" );
		//		working_align_pose_->dump_pdb( "WORKING_ALIGN.pdb" );
		//		pose.dump_pdb( "WORKING_POSE.pdb" );


	}

	//////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::align_pose( core::pose::Pose & pose ) const {
		if (!ready_to_align_ || !working_align_pose_ ){
			utility_exit_with_message( "Called align pose, but PoseSetup wasn't ready to do it..." );
		}
		core::scoring::superimpose_pose( pose, *working_align_pose_, alignment_atom_id_map_);
	}

	//////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size > const
	StepWiseProteinPoseSetup::apply_full_to_sub_mapping( utility::vector1< Size > & res_vector) const{
		return working_parameters_->apply_full_to_sub_mapping( res_vector );
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	stepwise::modeler::working_parameters::StepWiseWorkingParametersOP & StepWiseProteinPoseSetup::working_parameters(){
		return working_parameters_;
	}

	/////////////////////////////////////////////////////////////////
	void StepWiseProteinPoseSetup::figure_out_Prepend_Internal( pose::Pose const & pose ){

		utility::vector1< core::Size > const & working_moving_res_list( working_parameters_->working_moving_res_list() );
		if ( working_moving_res_list.size() < 1 ) return; //happens if we are just setting up the pose to, e.g., prepack.

		Size const first_working_moving_res  = working_moving_res_list[ 1 ];
		Size const last_working_moving_res   = working_moving_res_list[ working_moving_res_list.size() ];

		//		Size moving_suite( 0 );
		utility::vector1< core::Size > working_moving_suite_list;

		bool moving_res_attached_at_start = true;
		bool moving_res_attached_at_end = true;

		//////////////////////////////////////////////////////
		if ( first_working_moving_res == 1 ||
			pose.fold_tree().is_cutpoint( first_working_moving_res - 1 ) ) {
			moving_res_attached_at_start = false;
		}

		if ( last_working_moving_res == pose.total_residue() ||
				 pose.fold_tree().is_cutpoint( last_working_moving_res ) ){
			moving_res_attached_at_end = false;
		}

		// It might actually make more sense to replace is_prepend and is_internal
		//  with moving_res_attached_at_start  and   moving_res_attached_at_end...
		//
		// Or perhaps we should just get rid of these booleans -- what about the more complex
		// case in which the N-terminus and C-terminus are being sampled?
		//
		//		bool const is_prepend = moving_res_attached_at_end;
		bool const is_internal = ( moving_res_attached_at_end && moving_res_attached_at_start );

		//////////////////////////////////////////////////////
		// Fill out working_moving_suite.
		if ( moving_res_attached_at_start && !is_internal ){
			working_moving_suite_list.push_back( first_working_moving_res - 1 );
		}

		for ( Size n = 1; n < working_moving_res_list.size(); n++ ){
			working_moving_suite_list.push_back( working_moving_res_list[ n ] );
		}

		if ( moving_res_attached_at_end ){
			working_moving_suite_list.push_back( last_working_moving_res );
		}

		//		TR.Debug << "WORKING_MOVING_SUITE_LIST: " << std::endl;
		//		for ( Size i = 1; i <= working_moving_suite_list.size(); i++ )  TR.Debug << ' ' << working_moving_suite_list[ i ];
		//		TR.Debug << std::endl;

		//working_parameters_->set_moving_res_attached_at_start( moving_res_attached_start  );
		//working_parameters_->set_moving_res_attached_at_end(   moving_res_attached_at_end );

		//		TR.Debug << "ATTACH AT START: " << moving_res_attached_at_start << std::endl;
		//		TR.Debug << "ATTACH AT END:   " << moving_res_attached_at_end << std::endl;

		working_parameters_->set_is_prepend(  moving_res_attached_at_end  );

		// Watch out! is_internal gets replaced later based on "partition_definition" --
		//  occurs in complex fold_trees.
		working_parameters_->set_is_internal( moving_res_attached_at_end && moving_res_attached_at_start  );

		//		working_parameters_->set_working_moving_suite_list( working_moving_suite_list ); // auto updated...

	}

	//////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseProteinPoseSetup::is_working_cutpoint_closed( Size const res, std::map< Size, Size > & full_to_sub ) const {

		for ( Size m = 1; m <= cutpoint_closed_.size();   m++ ) {
			if ( res > 0 &&
					 res == full_to_sub[ cutpoint_closed_[m] ]   &&
					 (res+1) == full_to_sub[ cutpoint_closed_[m]+1 ] ) return true;
		}
		return false;
	}

	//////////////////////////////////////////////////////////////////////////////
	// Assume we have done a crappy job of placing 5' phosphates.
	//  this would be true if, for example, the 5' phosphate was built
	//  by prepending in a previous rebuild-from-scratch effort.
	// The only exception is if the 5'-phosphate is involved in *chain closure*.
	//(5' phosphate actually refer to the phosphate group of the residue 3' of the chain_break!!! Jan 29, 2010 Parin S.)
	void
	StepWiseProteinPoseSetup::apply_virtual_phosphate_variants( pose::Pose & pose ) const{

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries(	 working_parameters_->chain_boundaries() );
		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );

		Size const num_chains = chain_boundaries.size();

		for ( Size n = 1; n <= num_chains; n++ ) {
			Size const chain_start = chain_boundaries[ n ].first;
			Size const potential_five_prime_res = full_to_sub[ chain_start ];

			if (is_working_cutpoint_closed( potential_five_prime_res-1, full_to_sub ) ) continue;

			if ( !pose.residue( potential_five_prime_res ).is_RNA() ) continue;

			if ( pose.residue( potential_five_prime_res ).type().has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
				utility_exit_with_message( "Should not be trying to virtualize phosphate on close cutpoint residue!" );
			}

			pose::add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, full_to_sub[ chain_start ] );
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// remove?
	void
	StepWiseProteinPoseSetup::apply_peptide_plane_variants_OLD( pose::Pose & /*pose*/ ) const{

		if ( !add_peptide_plane_variants_ ) return;

// 		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );
// 		utility::vector1< core::Size > working_moving_res_list = working_parameters_->working_moving_res_list();

// 		for ( Size n = 1; n <= working_moving_res_list.size(); n++ ) {

// 			Size const moving_res( working_moving_res_list[ n ] );

// 			if (!pose.residue_type( moving_res ).is_protein() ) continue;

// 			if ( ( moving_res == pose.total_residue() || pose.fold_tree().is_cutpoint( moving_res ) ) &&
// 					 ( moving_res != full_to_sub[ cutpoint_closed_ ] ) &&
// 					 !pose.residue_type( moving_res ).has_variant_type( "UPPER_TERMINUS" ) ) {
// 				pose::add_variant_type_to_pose_residue( pose, "C_METHYLAMIDATION", moving_res );
// 			}

// 			if ( ( moving_res == 1 || pose.fold_tree().is_cutpoint( moving_res-1 ) ) &&
// 					 ( moving_res != full_to_sub[ cutpoint_closed_ ] + 1) &&
// 					 !pose.residue_type( moving_res ).has_variant_type( "LOWER_TERMINUS" ) ) {
// 				pose::add_variant_type_to_pose_residue( pose, "N_ACETYLATION", moving_res );
// 			}

// 		}
	}

	//////////////////////////////////////////////////////////////////////////////
	// Add peptide plane "caps" to moving ends of protein... lets phi/psi be sampled
	// and moves around hydrogen bond donor (NH) or acceptor (backbone C=O).
	void
	StepWiseProteinPoseSetup::apply_peptide_plane_variants( pose::Pose & pose ) const{

		if ( !add_peptide_plane_variants_ ) return;

		utility::vector1< std::pair< core::Size, core::Size > > const & chain_boundaries(	 working_parameters_->chain_boundaries() );
		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );

		Size const num_chains = chain_boundaries.size();

		for ( Size n = 1; n <= num_chains; n++ ) {

			Size const chain_start = chain_boundaries[ n ].first;
			Size const potential_Nterm_res = full_to_sub[ chain_start ];
			if ( !is_working_cutpoint_closed( potential_Nterm_res-1, full_to_sub ) &&
					 pose.residue( potential_Nterm_res ).is_protein() &&
					 !pose.residue( potential_Nterm_res ).has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) {
				pose::add_variant_type_to_pose_residue( pose, core::chemical::N_ACETYLATION, potential_Nterm_res );
			}

			Size const chain_end = chain_boundaries[ n ].second;
			Size const potential_Cterm_res = full_to_sub[ chain_end ];
			if ( !is_working_cutpoint_closed( potential_Cterm_res, full_to_sub ) &&
					 pose.residue( potential_Cterm_res ).is_protein() &&
					 !pose.residue( potential_Cterm_res ).has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT ) ) {
				pose::add_variant_type_to_pose_residue( pose, core::chemical::C_METHYLAMIDATION, potential_Cterm_res );
			}


		}

	}


	////////////////////////////////////////////////////////////////////////////////////////
 	void
 	StepWiseProteinPoseSetup::apply_virtual_res_variant(pose::Pose & pose ) const {
		//Parin Jan 17, 2010
		// For RNA only, for now.
		// Rhiju added proteins May 5, 2010
		using namespace core::id;

		pose::Pose pose_copy = pose;

		std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );

		for(Size i=1; i<=virtual_res_list_.size(); i++){
			Size const seq_num=virtual_res_list_[i];
			if (full_to_sub.find( seq_num ) != full_to_sub.end() ) {
				pose::add_variant_type_to_pose_residue( pose,
						core::chemical::VIRTUAL_RESIDUE_VARIANT, full_to_sub[ seq_num] );
			}
		}
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::add_terminal_res_repulsion( core::pose::Pose & pose ) const
	{
		using namespace core::conformation;
		using namespace core::id;
		using namespace core::scoring::constraints;
		using namespace core::chemical::rna;

		ConstraintSetOP cst_set( pose.constraint_set()->clone() );
		assert( cst_set );

		utility::vector1< core::Size > const & working_terminal_res = working_parameters_->working_terminal_res();

		/////////////////////////////////////////////////
		Size const nres( pose.total_residue() );

		ObjexxFCL::FArray1D<bool> is_moving_res( nres, false );
		ObjexxFCL::FArray1D<bool> is_fixed_res( nres, false );

		ObjexxFCL::FArray1D< bool > const & partition_definition = working_parameters_->partition_definition();
		bool const root_partition = partition_definition( pose.fold_tree().root() );

		for (Size seq_num=1; seq_num <= nres; seq_num++){
			if ( partition_definition( seq_num ) == root_partition ) {
				is_fixed_res( seq_num ) = true;
			} else {
				is_moving_res( seq_num ) = true;
			}
		}


		/////////////////////////////////////////////////
		Distance const DIST_CUTOFF = 8.0;
		core::scoring::func::FuncOP const repulsion_func( new core::scoring::func::FadeFunc( -2.0 /*min*/, DIST_CUTOFF /*max*/, 1.0 /*fade zone width*/, 100.0 /*penalty*/ ) );

		for ( Size i = 1; i <= working_terminal_res.size(); i++ ) {

			Size const k = working_terminal_res[ i ];

			for ( Size m = 1; m <= nres; m++ ) {

				Residue const & rsd1( pose.residue( k ) );
				Residue const & rsd2( pose.residue( m ) );

				AtomID const atom_id1( first_base_atom_index( rsd1 ), rsd1.seqpos() );
				AtomID const atom_id2( first_base_atom_index( rsd2 ), rsd2.seqpos() );

				// the one exception -- close contacts WITHIN a partition
				if ( ( ( is_moving_res( k )  && is_moving_res( m ) ) ||
							 ( is_fixed_res(  k )  && is_fixed_res(  m ) ) ) &&
						 ( pose.xyz( atom_id1 ) - pose.xyz( atom_id2 ) ).length() < DIST_CUTOFF ) {
					//TR.Debug << "Not adding repulsive constraint between " << k << " and " << m << " already closeby in same partition" << std::endl;
					continue;
				}

				// distance from O3' to P
				cst_set->add_constraint( new AtomPairConstraint( atom_id1, atom_id2, repulsion_func ) );

			}
		}

		pose.constraint_set( cst_set );

		/////////////////////////////////////////////////////////////
		// New virtualize terminal res.
		//   keep the hydrogens in there though -- they supply sterics.
		//		for ( Size i = 1; i <= working_terminal_res.size(); i++ ) {
		//			pose::add_variant_type_to_pose_residue( pose, chemical::VIRTUAL_BASE_HEAVY_ATOM, working_terminal_res[ i ] );
		//		}


	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::setup_secstruct( core::pose::Pose & pose ) const {

		utility::vector1< Size > const & is_working_res( working_parameters_->is_working_res() );

		if ( secstruct_.size() > 0 ){

			if ( secstruct_.size() != desired_sequence_.size() ) utility_exit_with_message( "Input secstruct does not have same size as full length sequence" );

			Size count = 0;
			for( Size n = 1; n <= desired_sequence_.size(); n++ ){
				if( !is_working_res[ n ] ) continue;
				count++;
				pose.set_secstruct( count, secstruct_[n-1] );
			}

		} else {

			for( Size n = 1; n <= pose.total_residue(); n++ ) {
				pose.set_secstruct( n, 'L' );
			}

		}


	}


	/////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::setup_full_model_info( core::pose::Pose & pose ) const {
		FullModelInfoOP full_model_info = new FullModelInfo( desired_sequence_, cutpoint_open_,
																												 working_parameters_->working_res_list() );

		FullModelParametersOP full_model_parameters = full_model_info->full_model_parameters()->clone();
		utility::vector1< Size > const & is_working_res = working_parameters_->is_working_res();
		utility::vector1< Size > extra_minimize_res = extra_minimize_res_; // may be updated.
		utility::vector1< Size > fixed_domain; // convert from old convention to new one stored in full_model info.
		for ( Size n = 1; n <= is_working_res.size(); n++ ){
			if ( is_working_res[n] == MOVING_RES || is_working_res[n] == BRIDGE_RES ){
				fixed_domain.push_back( 0 ); // moving and sample-able.
			} else if ( fixed_res_.has_value( n ) ) {
				fixed_domain.push_back( 1 ); // freeze this.
			} else {
				fixed_domain.push_back( is_working_res[n] ); // 1 or 2, for input domains 1 or 2.
				if ( is_working_res[n] > 0 && !extra_minimize_res_.has_value( n ) ) extra_minimize_res.push_back( n );
			}
		}
		full_model_parameters->set_parameter( FIXED_DOMAIN, fixed_domain );
		full_model_parameters->set_parameter_as_res_list( EXTRA_MINIMIZE, extra_minimize_res );
		full_model_parameters->set_parameter_as_res_list( CALC_RMS, calc_rms_res_ );
		full_model_info->set_full_model_parameters( full_model_parameters );
		set_full_model_info( pose, full_model_info );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	//Adding Virtual res as root
	void
	StepWiseProteinPoseSetup::add_aa_virt_rsd_as_root( core::pose::Pose & pose){  //Fang's electron density code

		Size const nres = pose.total_residue();
		//if already rooted on virtual residue , return
		if ( pose.residue( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) {
			TR.Warning << "addVirtualResAsRoot() called but pose is already rooted on a VRT residue ... continuing." << std::endl;
			return;
		}

		core::chemical::ResidueTypeSet const & residue_set = pose.residue_type(1).residue_type_set();
		core::chemical::ResidueTypeCOPs const & rsd_type_list( residue_set.name3_map("VRT") );
		core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *rsd_type_list[1] ) );
		pose.append_residue_by_jump( *new_res , 1 );

		// make the virt atom the root
		kinematics::FoldTree newF( pose.fold_tree() );
		newF.reorder( nres+1 );
		pose.fold_tree( newF );

		utility::vector1< core::Size > working_fixed_res = working_parameters_->working_fixed_res();
		working_fixed_res.push_back( nres+1 );
		working_parameters_->set_working_fixed_res( working_fixed_res );

	}
	////////////////////////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_fixed_res( utility::vector1 < core::Size > const & fixed_res ){
		fixed_res_ = fixed_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_jump_res( utility::vector1 < core::Size > const & jump_res ){
		jump_res_ = jump_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_virtual_res( utility::vector1 < core::Size > const & set_virtual_res_list){
		virtual_res_list_ = set_virtual_res_list;
	}
  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_terminal_res( utility::vector1 < core::Size > const & terminal_res ){
		terminal_res_ = terminal_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_bulge_res( utility::vector1 < core::Size > const & bulge_res ){
		bulge_res_ = bulge_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_bridge_res( utility::vector1 < core::Size > const & bridge_res ){
		bridge_res_ = bridge_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_rsd_set( core::chemical::ResidueTypeSetCAP & rsd_set ){
		rsd_set_ = rsd_set;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_cst_file( std::string const cst_file ){
		cst_file_ = cst_file;
	}


  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_disulfide_file( std::string const disulfide_file ){
		disulfide_file_ = disulfide_file;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_align_file( std::string const align_file ){
		align_file_ = align_file;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_add_peptide_plane_variants( bool const & setting ){
		add_peptide_plane_variants_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res ){
		calc_rms_res_ = calc_rms_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_superimpose_res( utility::vector1< core::Size > const & superimpose_res ){
		superimpose_res_ = superimpose_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_dump( bool const dump ){
		dump_ = dump;
	}

  //////////////////////////////////////////////////////////////////////////
	bool
	StepWiseProteinPoseSetup::ready_to_align() const{
		return ready_to_align_;
	}


  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::set_secstruct( std::string const secstruct ){
		secstruct_ = secstruct;
	}

	void
	StepWiseProteinPoseSetup::set_add_virt_res( bool const setting ){
		add_virt_res_ = setting;
	}

} //protein
} //modeler
} //legacy
} //stepwise
} //protocols
