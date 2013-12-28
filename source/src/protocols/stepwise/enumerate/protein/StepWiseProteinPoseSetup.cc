// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinPoseSetup
/// @brief Create starting pose for rebuilding, from input pose
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/enumerate/protein/StepWiseProteinPoseSetup.hh>
#include <protocols/stepwise/enumerate/protein/StepWiseProteinUtil.hh>
#include <protocols/stepwise/StepWiseUtil.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/import_pose/import_pose.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>


#include <utility/vector1.hh>

#include <list>

//Auto Headers
#include <core/id/AtomID.hh>

using namespace core;
using core::Real;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
	StepWiseProteinPoseSetup::StepWiseProteinPoseSetup( std::string const & desired_sequence, utility::vector1< std::string > const & start_tags,
																				utility::vector1< std::string > const & silent_files_in ):
		Mover(),
		desired_sequence_( desired_sequence ),
		start_tags_( start_tags ),
		silent_files_in_( silent_files_in ),
		rsd_set_( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ),
		nres_( desired_sequence_.size() ),
		input_residue_array_( nres_, false ),
		junction_residue_array_( nres_, false ),
		moving_residue_array_( nres_, false ),
		sample_junction_( true ),
		n_terminus_( false ),
		c_terminus_( false ),
		add_peptide_plane_( false ),
		verbose_( false )
  {
		Mover::type("StepWiseProteinPoseSetup");

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


	///////////////////////////////////////////////////////////////////////////////////
	void
	fix_end_phi_psi( pose::Pose & pose, pose::Pose & start_pose, Size const & start_res, Size const & end_res ){

		if ( end_res < pose.total_residue() ){
			Real const psi_current = get_pretend_psi_explicit( pose, end_res );
			Real const psi_scratch = get_pretend_psi_explicit( start_pose, start_pose.total_residue() );
			pose.set_psi( end_res, pose.psi( end_res ) + psi_scratch - psi_current );
		}

		if ( start_res > 1 ){
			Real const phi_current = get_pretend_phi_explicit( pose, start_res );
			Real const phi_scratch = get_pretend_phi_explicit( start_pose, 1 );
			pose.set_phi( start_res, pose.phi( start_res ) + phi_scratch - phi_current );
		}

	}


	////////////////////////////////////////////////////
	////////////////////////////////////////////////////
	// Need to be very careful about TERMINI.
	// Note, in principle, we could inherit these if the target pose
	// desires termini variants at end.

  //////////////////////////////////////////////////////////////////////////
	void
  StepWiseProteinPoseSetup::apply( Pose & pose )
	{

		using namespace core::pose;

		//		if ( start_tags_.size() == 0 ){
		//			initialize_from_scratch( pose );
		//			return;
		//		}
		make_pose_from_sequence( pose, desired_sequence_, *rsd_set_, false /*auto_termini*/);

		if ( verbose_ )		pose.dump_pdb( "START1.pdb" );
		add_end_variants( pose );
		if ( verbose_ )		pose.dump_pdb( "START2.pdb" );

		// make extended chain
		for ( Size pos = 1; pos <= pose.total_residue(); pos++ ) {
			if ( ! pose.residue(pos).is_protein() ) continue;
			pose.set_phi( pos, -150 );
			pose.set_psi( pos, 150);
			pose.set_omega( pos, 180 );
		}

		if ( verbose_ )		pose.dump_pdb( "START2x.pdb" );

		// Later this will become more sophisticated, allowing for several pose moving_elements to be read in,
		//  aligned to a sequence alignment, and in-between residues samples (stored in "moving_residues_" ).

		utility::vector1< Size > start_res_all, end_res_all;

		if ( verbose_ )		pose.dump_pdb( "START_BEFORE_COPY_DOFS.pdb" );

		for ( Size i = 1; i <= start_tags_.size(); i++ ) {

			std::string const & start_tag( start_tags_[ i ] );

			//std::cout << "--- ABOUT TO READ IN INITIAL POSE" << i << std::endl;

			Pose start_pose;
			if ( silent_files_in_.size() > 0 ) {
				// Read in from silent file
				core::io::silent::SilentFileData silent_file_data;
				silent_file_data.read_file( silent_files_in_[ i ] );
				bool found_tag( false );
				for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(),
								end = silent_file_data.end(); iter != end; ++iter ) {
					if ( iter->decoy_tag() != start_tag ) continue;
					found_tag = true;
					iter->fill_pose( start_pose, *rsd_set_ );
					break;
				}
				if ( !found_tag ) utility_exit_with_message( "Could not find specified tag in silent file!" );
			} else {
				// Read in from PDB
				import_pose::pose_from_pdb( start_pose, *rsd_set_, start_tag );
			}

			//			remove_end_variants( start_pose );

			Size const start_res = figure_out_nested_positions( start_pose.sequence(), desired_sequence_ );
			Size const end_res = start_res +  start_pose.sequence().size() - 1;

			for ( Size k = start_res; k <= end_res; k++ )	 input_residue_array_( k ) = true;
			junction_residue_array_( start_res ) = true;
			junction_residue_array_( end_res ) = true;

			//Now actually copy into the pose.
			ResMap res_map;
			Size count( 0 );
			for ( Size n = start_res; n <= end_res; n++ ) {
				count++;
				res_map[ n ] = count;
			}

			if ( verbose_ )		start_pose.dump_pdb( "INPUT.pdb" );

			match_end_variants( pose, start_pose, start_res, end_res );

			if ( verbose_ )		start_pose.dump_pdb( "MATCH_END.pdb" );

			copy_dofs( pose, start_pose, res_map );

			fix_end_phi_psi( pose, start_pose, start_res, end_res );

		}

		junction_residue_array_( 1 ) = false;
		junction_residue_array_( nres_ ) = false;

		// List of moving residues -- useful for outside routines.
		moving_residues_.clear();
		moving_residue_array_ = false;
		for (Size n = 1; n <= nres_; n++ ) {
			if ( ! input_residue_array_( n ) ) moving_residue_array_( n ) = true;
			if (sample_junction_ && junction_residue_array_( n ) ) moving_residue_array_( n ) = true;
			if ( moving_residue_array_( n ) ) {
				moving_residues_.push_back( n );
				//std::cout << "MOVING: " << n << std::endl;
			}
		}

		if ( verbose_ )			pose.dump_pdb( "START3.pdb" );

		//		prepend_residues( pose );
		//		append_residues ( pose, start_pose );

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::add_end_variants( core::pose::Pose & pose ){

		if ( n_terminus_ ) {
			pose::add_variant_type_to_pose_residue( pose, "LOWER_TERMINUS", 1  );
		} else {
			if (add_peptide_plane_) pose::add_variant_type_to_pose_residue( pose, "N_ACETYLATION", 1  );
		}

		if ( c_terminus_ ) {
			pose::add_variant_type_to_pose_residue( pose, "UPPER_TERMINUS", pose.total_residue()  );
		} else {
			if (add_peptide_plane_) pose::add_variant_type_to_pose_residue( pose, "C_METHYLAMIDATION", pose.total_residue()  );
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::match_specific_variants( pose::Pose const & pose, Size const & pose_res,
																							pose::Pose & scratch_pose, Size const & scratch_pose_res,
																							utility::vector1< chemical::VariantType > const & variant_types ) const
	{

		using namespace core::chemical;
		utility::vector1< VariantType > variant_types_to_add, variant_types_to_remove;

		for (Size n = 1; n <= variant_types.size(); n++ ) {

			VariantType const & variant_type = variant_types[ n ];
			if ( pose.residue( pose_res ).has_variant_type( variant_type ) &&
					 !scratch_pose.residue( scratch_pose_res ).has_variant_type( variant_type ) ){
				variant_types_to_add.push_back( variant_type );
			}

			if ( !pose.residue( pose_res ).has_variant_type( variant_type ) &&
					 scratch_pose.residue( scratch_pose_res ).has_variant_type( variant_type ) ){
				variant_types_to_remove.push_back( variant_type );
			}

		}

		core::pose::Pose scratch_pose_save = scratch_pose;

		// Not exactly sure if I need to do this -- all the removes, then all the adds...
		for ( Size n = 1; n <= variant_types_to_remove.size(); n++ )  {
			pose::remove_variant_type_from_pose_residue( scratch_pose, variant_types_to_remove[n], scratch_pose_res   );

			// This is really silly.
			if ( variant_types_to_remove[n] == "LOWER_TERMINUS" ) {
							scratch_pose.set_xyz( core::id::AtomID( scratch_pose.residue( scratch_pose_res).atom_index("H"),
																											scratch_pose_res ),
																		scratch_pose_save.xyz( core::id::AtomID( scratch_pose_save.residue( scratch_pose_res).atom_index( "1H" ), scratch_pose_res ) ) );
			}

		}


		for ( Size n = 1; n <= variant_types_to_add.size(); n++ )  {
			pose::add_variant_type_to_pose_residue( scratch_pose, variant_types_to_add[n], scratch_pose_res  );
		}

	}

	///////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseSetup::match_end_variants( pose::Pose const & pose, pose::Pose & scratch_pose, Size const & start_res, Size const & end_res ) const
	{
		using namespace core::chemical;

		utility::vector1< VariantType >  N_variants;
		N_variants.push_back( "LOWER_TERMINUS" );
		N_variants.push_back( "N_ACETYLATION" );

		match_specific_variants( pose, start_res,   scratch_pose, 1, N_variants );

		utility::vector1< VariantType >  C_variants;
		C_variants.push_back( "UPPER_TERMINUS" );
		C_variants.push_back( "C_METHYLAMIDATION" );

		match_specific_variants( pose, end_res,   scratch_pose, scratch_pose.total_residue(), C_variants );

	}


	////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	StepWiseProteinPoseSetup::figure_out_nested_positions(
   std::string const & inside_sequence,
	 std::string const & desired_sequence ) const
	{

		Size const max_start_res = desired_sequence.size() - inside_sequence.size()+1;

		for (Size potential_start_res = 1; potential_start_res <= max_start_res; potential_start_res++ ) {
			bool found_match = true;

			//Position may already be "called for" by a previous assignment.
			if ( input_residue_array_( potential_start_res )  ) continue;

			//Look for exact sequence match.
			for ( Size i = 0; i < inside_sequence.size(); i++ ) {
				if ( inside_sequence[i] != desired_sequence[i - 1 + potential_start_res] ) {
					found_match = false;
					break;
				}
			}
			if ( found_match ) {
				return potential_start_res;
			}
		}

		return 0;

	}

	//////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size > const &
	StepWiseProteinPoseSetup::moving_residues() const
	{
		return moving_residues_;
	}

// 	//////////////////////////////////////////////////////////////////////////////
// 	void
// 	StepWiseProteinPoseSetup::prepend_residues( pose::Pose & pose ) {

// 		using namespace core::chemical;
// 		using namespace core::conformation;

// 		if ( start_res_ <= 1 )  return;

// 		pose::remove_lower_terminus_type_from_pose_residue( pose, 1   );
// 		for ( Size i = start_res_ - 1; i > 0; i-- ) {
// 			ResidueTypeCOP new_rsd_type( ResidueSelector().set_name1( desired_sequence_[i-1] ).exclude_variants().select( *rsd_set_ )[1] );
// 			//			ResidueTypeCOP new_rsd_type( rsd_set_->aa_map( aa_from_oneletter_code( desired_sequence[i-1] ) )[1] );
// 			ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
// 			pose.conformation().safely_prepend_polymer_residue_before_seqpos( *new_rsd, 1, true );
// 			pose.set_omega( 1, 180.0 ); //NOTE SIDE CHAINS ARE MESSED UP!
// 		}
// 		for (Size i = 1; i < start_res_ ; i++ ) {
// 			moving_residues_.push_back( i );
// 		}
// 	}


// 	//////////////////////////////////////////////////////////////////////////////
// 	void
// 	StepWiseProteinPoseSetup::append_residues( pose::Pose & pose, pose::Pose const & start_pose ) {

// 		using namespace core::chemical;
// 		using namespace core::conformation;

// 		if ( end_res_ >= desired_sequence_.size() ) return;

// 		pose::remove_upper_terminus_type_from_pose_residue( pose, pose.total_residue()   );
// 		for ( Size i = end_res_ + 1; i <= desired_sequence_.size(); i++ ) {
// 			ResidueTypeCOP new_rsd_type( ResidueSelector().set_name1( desired_sequence_[i-1] ).exclude_variants().select( *rsd_set_ )[1] );
// 			ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
// 			pose.conformation().append_polymer_residue_after_seqpos( *new_rsd, pose.total_residue(), true );
// 			pose.set_omega( pose.total_residue(), 180.0 );
// 			moving_residues_.push_back( i );
// 		}
// 		//		pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_CTERM_O", pose.total_residue() );
// 		pose.set_psi( end_res_, start_pose.psi( start_pose.total_residue() ) );
// 	}

// 	//////////////////////////////////////////////////////////////////////////////
// 	void
// 	StepWiseProteinPoseSetup::initialize_from_scratch( pose::Pose & pose )
// 	{
// 		make_pose_from_sequence( pose, desired_sequence_, *rsd_set_, false /*auto_termini*/);

// 		for ( Size n = 1; n <= pose.total_residue(); n++ ) {
// 			moving_residues_.push_back( n );
// 		}
// 	}

/////////////////////////////////////////////////////////////////

	void
	StepWiseProteinPoseSetup::set_add_peptide_plane( bool const & setting ) {
		add_peptide_plane_ = setting;
	}

} //protein
} //enumerate
} //stepwise
} //protocols
