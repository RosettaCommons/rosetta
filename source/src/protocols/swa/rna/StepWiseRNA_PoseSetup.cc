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
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/StepWiseUtil.hh>

#include <protocols/rna/RNA_ProtocolUtil.hh>
//////////////////////////////////
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/FadeFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <core/io/silent/BinaryRNASilentStruct.hh> //Feb 24, 2011, FARFAR start_silent_file.

#include <utility/exit.hh>
#include <time.h>

#include <string>

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

static basic::Tracer TR( "protocols.swa.rna.stepwise_rna_pose_setup" ) ;

namespace protocols {
namespace swa {
namespace rna {

//////////////////////////////////////////////////////////////////////////
//constructor!
StepWiseRNA_PoseSetup::StepWiseRNA_PoseSetup( StepWiseRNA_JobParametersOP & job_parameters ):
	rsd_set_( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA ) ),
	job_parameters_( job_parameters ),
	copy_DOF_( false ),
	verbose_( true ),
	rebuild_bulge_mode_( false ), //Nov 26, 2010
	output_pdb_( false ), //Sept 24, 2011
	apply_virtual_res_variant_at_dinucleotide_( true ), //Jul 12, 2012
	align_to_native_( false ), //Jul 12, 2012
	use_phenix_geo_( false )
{}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//destructor
StepWiseRNA_PoseSetup::~StepWiseRNA_PoseSetup()
{}

/////////////////////
std::string
StepWiseRNA_PoseSetup::get_name() const {
	return "StepWiseRNA_PoseSetup";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
StepWiseRNA_PoseSetup::apply( core::pose::Pose & pose ) {

	using namespace ObjexxFCL;

	if ( verbose_ ) Output_title_text( "Enter StepWiseRNA_PoseSetup::apply", TR.Debug );

	// actually make the pose, set fold tree, copy in starting templates from disk.

	pose::Pose pose_without_cutpoints;

	if ( copy_DOF_ ){

		make_pose( pose_without_cutpoints ); //Create pose with random torsions
		TR.Debug << "read_input_pose_and_copy_dofs( pose_without_cutpoints )" << std::endl;
		read_input_pose_and_copy_dofs( pose_without_cutpoints );

		make_pose( pose ); //Create pose with random torsions
		pose.fold_tree( job_parameters_->fold_tree() );
		TR.Debug << "read_input_pose_and_copy_dofs( pose )" << std::endl;
		read_input_pose_and_copy_dofs( pose );

	} else{
		pose_without_cutpoints = pose; //Hard copy, slow
		pose.fold_tree( job_parameters_->fold_tree() );
	}

	//Fang: apply ideal coord to residue being modeled if use_phenix_geo
	if ( use_phenix_geo_ ) {
		Size const moving_res = job_parameters_->working_moving_res();
		pose::rna::apply_pucker( pose, moving_res );
	}

	if ( output_pdb_ ) pose.dump_pdb( "test.pdb" );

	//WARNING STILL NEED TO IMPLEMENT harmonic_chainbreak HERE!
	apply_cutpoint_variants( pose, pose_without_cutpoints );
	apply_virtual_phosphate_variants( pose );
	add_protonated_H1_adenosine_variants( pose );
	apply_bulge_variants( pose );
	apply_virtual_res_variant( pose ); //User specified virtaul_res and bulge_res of it being built
	verify_protonated_H1_adenosine_variants( pose );
	add_terminal_res_repulsion( pose );

	if ( job_parameters_->add_virt_res_as_root() ) add_aa_virt_rsd_as_root( pose ); //Fang's electron density code.

	setup_pdb_info_with_working_residue_numbers( pose );

	if ( output_pdb_ ) pose.dump_pdb( "start.pdb" );

	if ( verbose_ ) Output_title_text( "Exit StepWiseRNA_PoseSetup::apply", TR.Debug );
}

////////////////////////////////////////////////////////


void
StepWiseRNA_PoseSetup::setup_native_pose( core::pose::Pose & pose ){
	using namespace core::conformation;
	using namespace core::pose;
	using namespace protocols::rna;

	if ( verbose_ ) Output_title_text( "Enter StepWiseRNA_PoseSetup::setup_native_pose", TR.Debug );

	utility::vector1< core::Size > const & is_working_res( job_parameters_->is_working_res() );
	std::string const & working_sequence( job_parameters_->working_sequence() );
	utility::vector1< core::Size > const & working_best_alignment( job_parameters_->working_best_alignment() );
	utility::vector1< core::Size > const & working_moving_res_list( job_parameters_->working_moving_res_list() );

	utility::vector1< core::Size > const working_native_virtual_res_list_ = apply_full_to_sub_mapping( native_virtual_res_list_, StepWiseRNA_JobParametersCOP( job_parameters_ ) );

	utility::vector1 < core::Size > const & working_moving_partition_pos = job_parameters_->working_moving_partition_pos(); //Sept 14, 2011.

	if ( !get_native_pose() ) return;
	// Need a simple fold tree for following to work...
	Pose native_pose_copy = *get_native_pose();
	native_pose_copy.fold_tree(  core::kinematics::FoldTree(  native_pose_copy.total_residue() ) );

	PoseOP working_native_pose = new Pose;

	//First option is to pass in the full length native
	if ( native_pose_copy.sequence() ==  job_parameters_->full_sequence() ){
		TR.Debug << "User passed in full length native pose" << std::endl;

		( *working_native_pose ) = native_pose_copy;

		for ( Size seq_num = native_pose_copy.total_residue(); seq_num >= 1 ; seq_num-- ){
			if ( !is_working_res[ seq_num ] ){
				working_native_pose->conformation().delete_residue_slow( seq_num );
			}
		}

	} else if ( native_pose_copy.sequence() ==  working_sequence ){ //Could also pass in the working_native_pose directly
		TR.Debug << "User passed in working native pose" << std::endl;
		( *working_native_pose ) = native_pose_copy;
	} else{
		TR.Debug <<  std::setw( 50 ) << "  native_pose_copy.sequence() = " << native_pose_copy.sequence()  << std::endl;
		TR.Debug <<  std::setw( 50 ) << "  job_parameters_->full_sequence() = " << job_parameters_->full_sequence() << std::endl;
		TR.Debug <<  std::setw( 50 ) << "  working_sequence = " << working_sequence << std::endl;
		utility_exit_with_message( "The native pose passed in by the User does not match both the full_sequence and the working sequence of the inputted fasta_file" );
	}

	if ( working_native_pose->sequence() !=  working_sequence ){
		TR.Debug <<  std::setw( 50 ) << "working_native_pose->sequence() = " << working_native_pose->sequence();
		TR.Debug <<  std::setw( 50 ) << "working_sequence = " << working_sequence << std::endl;
		utility_exit_with_message( "working_native_pose->sequence() !=  working_sequence" );
	}


	protocols::rna::make_phosphate_nomenclature_matches_mini( ( *working_native_pose ) );


	utility::vector1< core::Size > act_working_alignment;
	act_working_alignment.clear();

	utility::vector1< core::Size > const & working_native_alignment = job_parameters_->working_native_alignment();

	if ( working_native_alignment.size() != 0 ){ //User specified the native alignment res.
		TR.Debug << "working_native_alignment.size() != 0!, align native_pose with working_native_alignment!" << std::endl;

		for ( Size n = 1; n <= working_native_alignment.size(); n++ ){
			Size const seq_num = working_native_alignment[n];
			if ( working_moving_res_list.has_value( seq_num ) ) continue;
			if ( working_moving_partition_pos.has_value( seq_num ) ) continue; //Sept 14, 2011.
			act_working_alignment.push_back( seq_num );
		}


	} else{ //Use default alignment res.
		TR.Debug << "working_native_alignment.size() == 0!, align native_pose with working_best_alignment!" << std::endl;

		for ( Size n = 1; n <= working_best_alignment.size(); n++ ){
			Size const seq_num = working_best_alignment[n];
			if ( working_moving_res_list.has_value( seq_num ) ) continue;
			if ( working_moving_partition_pos.has_value( seq_num ) ) continue; //Sept 14, 2011.
			act_working_alignment.push_back( seq_num );
		}
	}


	Output_seq_num_list( "act_working_alignment = ", act_working_alignment, TR.Debug );
	Output_seq_num_list( "working_moving_res_list = ", working_moving_res_list, TR.Debug );
	Output_seq_num_list( "working_moving_partition_pos = ", working_moving_partition_pos, TR.Debug );
	Output_seq_num_list( "working_native_alignment = ", working_native_alignment, TR.Debug );
	Output_seq_num_list( "working_best_alignment = ", working_best_alignment, TR.Debug );

	if ( act_working_alignment.size() == 0 ) utility_exit_with_message( "act_working_alignment.size() == 0" );


	if ( job_parameters_ -> add_virt_res_as_root() ) { //Fang's electron density code

		//Fang why do you need this?////
		pose::Pose dummy_pose = *working_native_pose;
		apply_cutpoint_variants( *working_native_pose, dummy_pose );
		add_protonated_H1_adenosine_variants( *working_native_pose );
		verify_protonated_H1_adenosine_variants( *working_native_pose );
		////////////////////////////////

		add_aa_virt_rsd_as_root( *working_native_pose );
		( *working_native_pose ).fold_tree( pose.fold_tree() );
		align_poses( pose, "working_pose", ( *working_native_pose ), "working_native_pose", act_working_alignment );

	} else if ( align_to_native_ ){
		// Rhiju -- prefer to align pose to native pose...
		align_poses( pose, "working_pose", ( *working_native_pose ), "working_native_pose", act_working_alignment );

	} else{ //March 17, 2012 Standard. This ensures that adding native_pdb does not change the final output silent_struct coordinates.
		align_poses( ( *working_native_pose ), "working_native_pose", pose, "working_pose", act_working_alignment );

	}




	job_parameters_->set_working_native_pose( working_native_pose );

	//Setup cutpoints, chain_breaks, variant types here!

	if ( output_pdb_ ) working_native_pose->dump_pdb( "working_native_pose.pdb" );
	if ( output_pdb_ ) pose.dump_pdb( "sampling_pose.pdb" );


	for ( Size n = 1; n <= working_native_virtual_res_list_.size(); n++ ){
//			pose::add_variant_type_to_pose_residue( (*working_native_pose), "VIRTUAL_RNA_RESID#UE", working_native_virtual_res_list_[n]);
		apply_virtual_rna_residue_variant_type( ( *working_native_pose ), working_native_virtual_res_list_[n], false /*apply_check*/ ) ;
	}


	if ( output_pdb_ ) working_native_pose->dump_pdb( "working_native_pose_with_virtual_res_variant_type.pdb" );



	if ( verbose_ ) Output_title_text( "Exit StepWiseRNA_PoseSetup::setup_native_pose", TR.Debug );

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::Import_pose( Size const & i, core::pose::Pose & import_pose ) const {

	using namespace core::conformation;


	if ( i > input_tags_.size() ) utility_exit_with_message( "i > input_tags_.size()!" );
	if ( copy_DOF_ == false ) utility_exit_with_message( "copy_DOF_ == false!" );

	if ( i > silent_files_in_.size() ) {
		// not a silent file, read in from pdb text file.

		std::string pose_name = input_tags_[ i ];
		std::size_t found = pose_name.find( ".pdb" );
		if ( found == std::string::npos ) {
			pose_name.append( ".pdb" );
		}

		//  if(verbose) TR.Debug << "	The following pose will be imported :" << pose_name << std::endl;
		import_pose::pose_from_pdb( import_pose, *rsd_set_, pose_name );
		protocols::rna::make_phosphate_nomenclature_matches_mini( import_pose );

	} else {

		import_pose_from_silent_file( import_pose, silent_files_in_[ i ], input_tags_[i] );

	}

	///////////////////////////////
	// Check for sequence match.
	utility::vector1< Size > const & input_res = job_parameters_->input_res_vectors()[i];
	std::string const & full_sequence = job_parameters_->full_sequence();

	if ( import_pose.total_residue() != input_res.size() ){
		 TR.Debug << "import_pose.total_residue() = " << import_pose.total_residue() << " input_res.size() = " << input_res.size() << std::endl;
		 utility_exit_with_message( "input pose does not have same # residues as input res" );
	}
	bool match( true );
	for ( Size n = 1; n <= import_pose.total_residue(); n++ ) {
		if (  import_pose.sequence()[ n - 1 ]  != full_sequence[ input_res[n] - 1 ] ) {
			match = false; break;
		}
	}
	if ( !match ) {
		TR.Error << "IMPORT_POSE SEQUENCE " << import_pose.sequence() << std::endl;
		TR.Error << "FULL SEQUENCE " << full_sequence << std::endl;
		utility_exit_with_message( "mismatch in sequence between input pose and desired sequence, given input_res " );
	}
}

////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::make_pose( pose::Pose & pose ){

	using namespace core::conformation;

	std::string const working_sequence = job_parameters_->working_sequence();

	//if ( copy_DOF_==false) utility_exit_with_message( "copy_DOF_==false!" );


	make_pose_from_sequence( pose, working_sequence, *rsd_set_, false /*auto_termini*/ );

	if ( output_pdb_ ){
		TR.Debug << "outputting extended_chain.pdb" << std::endl;
		pose.dump_pdb( "extended_chain.pdb" );
	}

}

///////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::read_input_pose_and_copy_dofs( pose::Pose & pose )
{

	std::string const working_sequence = job_parameters_->working_sequence();
	utility::vector1< utility::vector1< Size > > const & input_res_vectors = job_parameters_->input_res_vectors();

	if ( input_tags_.size() < 1 ) utility_exit_with_message( "input_tags_.size() < 1" );

	if ( copy_DOF_ == false ) utility_exit_with_message( "copy_DOF_ == false!" );

	using namespace core::pose;
	using namespace ObjexxFCL;

	std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );

	utility::vector1< pose::Pose > start_pose_list;

	Pose start_pose;
	Pose start_pose_with_variant;
	for ( Size i = 1; i <= input_tags_.size(); i++ ){

		TR.Debug << "import_pose " << i << std::endl;

		//Ok account for special case of build loop outward from scratch...
		if ( input_tags_[i] == "build_from_scratch" ){
			if ( input_res_vectors[i].size() != 1 ) utility_exit_with_message( "input_tags_[i] == \"build_from_scratch\" but input_res_vectors[i].size() != 1" );
			continue;
		}

		Import_pose( i, start_pose );
		start_pose_with_variant = start_pose;
		start_pose_list.push_back( start_pose );

		utility::vector1< Size > const & input_res = job_parameters_->input_res_vectors()[i]; //This is from input_pose numbering to full_pose numbering

		// Remove all variant types
		// DO TO LIST: Need to remove atom constraint and remove angle constaint as well
		//NOTES: June 16, 2011
		//Should LOWER_TERMINUS and UPPER_TERMINUS be removed as well? LOWER_TERMINUS does determine the position of  OP2 and OP1?
		//Also should then check that pose.residue_type(i).variant_types() is the empty?
		//Alternatively could convert to FARFAR way and use the NEW_copy_dof that match atom names (MORE ROBUST!). This way doesn't need to remove any variant type from the chunk_pose?

		utility::vector1< std::string > variant_type_list;
		variant_type_list.push_back( "VIRTUAL_PHOSPHATE" );
		variant_type_list.push_back( "VIRTUAL_O2STAR_HYDROGEN" );
		variant_type_list.push_back( "CUTPOINT_LOWER" );
		variant_type_list.push_back( "CUTPOINT_UPPER" );
		variant_type_list.push_back( "VIRTUAL_RNA_RESIDUE" );
		variant_type_list.push_back( "VIRTUAL_RNA_RESIDUE_UPPER" );
		variant_type_list.push_back( "BULGE" );
		variant_type_list.push_back( "VIRTUAL_RIBOSE" );
		variant_type_list.push_back( "PROTONATED_H1_ADENOSINE" );
		variant_type_list.push_back( "3PRIME_END_OH" ); 				//Fang's electron density code
		variant_type_list.push_back( "5PRIME_END_PHOSPHATE" ); //Fang's electron density code
		variant_type_list.push_back( "5PRIME_END_OH" ); 				//Fang's electron density code

		for ( Size seq_num = 1; seq_num <= start_pose.total_residue(); seq_num++  ) {
			for ( Size k = 1; k <= variant_type_list.size(); ++k ) {
				std::string const & variant_type = variant_type_list[k];

				if ( start_pose.residue( seq_num ).has_variant_type( variant_type ) ) {
					remove_variant_type_from_pose_residue( start_pose, variant_type, seq_num );
				}

			}
		}

		if ( output_pdb_ ) {
			start_pose.dump_pdb( "import_" + string_of( i ) + ".pdb" );
			start_pose_with_variant.dump_pdb( "import_orig_" + string_of( i ) + ".pdb" );
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Now actually copy into the pose.
		// Need to know correspondence of residues between imported pose and the working pose...
		if ( input_res.size() != start_pose.total_residue() ) {
			utility_exit_with_message( "Need to specify -already_built_residues, and the number of residues needs to match pdb file inputted after -s" );
		}

		std::map< core::Size, core::Size > res_map;  //This is map from sub numbering to input_res numbering..
		for ( Size n = 1; n <= input_res.size(); n++ ) {
			res_map[ full_to_sub[ input_res[n] ] ] = n;
			//TR.Debug << full_to_sub_[ input_res_[ n ] ] << " " << n << std::endl;
		}

		//Does this work for the "overlap residue" case?? If there is a overlap residue, then order of input_res will manner...Parin Jan 2, 2010.

		//Dec 24, 2011 Parin S.:Convert to Rhiju's NEW version
		//copy_dofs( pose, start_pose, res_map, true /*copy_dofs_for_junction_residues*/ );
		copy_dofs_match_atom_names( pose, start_pose, res_map, false /*backbone_only*/, false /*ignore_virtual*/ );

		//////////////////////Add virtual_rna_residue, virtual_rna_residue variant_upper and VIRTUAL_RIBOSE variant type back//////////////////////
		//////////////////////For standard dincleotide move, there should be virtual_rna_residue and virtual_rna_residue variant_upper //////////
		//////////////////////However if floating base sampling, will contain the virtual_ribose as well..


		utility::vector1< Size > const & cutpoint_closed_list = job_parameters_->cutpoint_closed_list();
		utility::vector1< Size > const & working_cutpoint_closed_list = apply_full_to_sub_mapping( cutpoint_closed_list, StepWiseRNA_JobParametersCOP( job_parameters_ ) );

		for ( Size n = 1; n <= start_pose_with_variant.total_residue(); n++  ) {
			if ( has_virtual_rna_residue_variant_type( start_pose_with_variant, n ) ){
				apply_virtual_rna_residue_variant_type( pose, full_to_sub[ input_res[n] ], working_cutpoint_closed_list ) ;
			}

			if ( start_pose_with_variant.residue( n ).has_variant_type( "VIRTUAL_RIBOSE" ) ){
				add_variant_type_to_pose_residue( pose, "VIRTUAL_RIBOSE", full_to_sub[ input_res[n] ] );
			}

			if ( start_pose_with_variant.residue( n ).has_variant_type( "3PRIME_END_OH" ) ) {
				add_variant_type_to_pose_residue( pose, "3PRIME_END_OH", full_to_sub[ input_res[n] ] );
			} else if ( start_pose_with_variant.residue( n ).has_variant_type( "5PRIME_END_PHOSPHATE" ) ) {
				add_variant_type_to_pose_residue( pose, "5PRIME_END_PHOSPHATE", full_to_sub[ input_res[n] ] );
			} else if ( start_pose_with_variant.residue( n ).has_variant_type( "5PRIME_END_OH" ) ) {
				add_variant_type_to_pose_residue( pose, "5PRIME_END_OH", full_to_sub[ input_res[n] ] );
			}

			if ( i > silent_files_in_.size() ) { // not a silent file, read in from pdb text file. May 04, 2011

				//start_pose_with_variant does not have PROTONATED_H1_ADENOSINE variant type since the input_pdb does not have the variant type or loses the variant when imported into Rosetta.

				if ( start_pose_with_variant.residue( n ).has_variant_type( "PROTONATED_H1_ADENOSINE" ) ) { //May 03, 2011
					Output_seq_num_list( "protonate_H1_adenosine_list = ", job_parameters_->protonated_H1_adenosine_list(), TR.Debug );
					utility_exit_with_message( "start_pose have PROTONATED_H1_ADENOSINE variant type at full_seq_num = " + ObjexxFCL::string_of( input_res[n] ) + " even though it was read in from PDB file!" );
				}

				if ( ( job_parameters_->protonated_H1_adenosine_list() ).has_value( input_res[n] ) ){
					apply_protonated_H1_adenosine_variant_type( pose, full_to_sub[ input_res[n] ] );
				}

			} else{

				if ( start_pose_with_variant.residue( n ).has_variant_type( "PROTONATED_H1_ADENOSINE" ) ) { //May 03, 2011

					if ( start_pose.residue( n ).aa() != core::chemical::na_rad ){
						utility_exit_with_message( "start_pose have PROTONATED_H1_ADENOSINE variant type at full_seq_num = " + ObjexxFCL::string_of( input_res[n] ) + " but rsd.aa() != core::chemical::na_rad!" );
					}

					if ( ( job_parameters_->protonated_H1_adenosine_list() ).has_value( input_res[n] ) == false ){
						Output_seq_num_list( "protonate_H1_adenosine_list = ", job_parameters_->protonated_H1_adenosine_list(), TR.Debug );
						utility_exit_with_message( "full_seq_num = " + ObjexxFCL::string_of( input_res[n] ) + " is have a PROTONATED_H1_ADENOSINE variant type in the start_pose but is not in the protonate_H1_adenosine_list" );
					}

					apply_protonated_H1_adenosine_variant_type( pose, full_to_sub[ input_res[n] ] );
				}
			}

		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if ( output_pdb_ ) pose.dump_pdb( "copy_dof_" + string_of( i ) + ".pdb" );
		TR.Debug << pose.fold_tree() << std::endl;
	}

	protocols::rna::assert_phosphate_nomenclature_matches_mini( pose ); //Just to be safe, Jun 11, 2010

	correctly_copy_HO2star_positions( pose, start_pose_list );


}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::correctly_copy_HO2star_positions( pose::Pose & working_pose, utility::vector1< pose::Pose > const & start_pose_list ){

	using namespace ObjexxFCL;
	using namespace core::conformation;
	using namespace core::id;

	if ( verbose_ ) Output_title_text( "Enter StepWiseRNA_PoseSetup::correctly_copy_HO2star_position", TR.Debug );

	if ( output_pdb_ ) working_pose.dump_pdb( "copy_dof_pose_BEFORE_correctly_copy_HO2star_positions.pdb" );

	std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );


	//Problem aside in the long_loop mode, where we combine two silent_file pose. The two pose both contain the fixed background residues
	//In current implementation of COPY_DOFS, the positioning of the background residues in SECOND silent file is used.
	//However, some background HO2star interacts with the first pose and hence need to correctly copys these HO2star position.
	//When we implemenet protien/RNA interactions, will extend this to include all "packed" side chains.

	utility::vector1< utility::vector1< Size > > const & input_res_vectors = job_parameters_->input_res_vectors();

	std::map< core::Size, core::Size > full_to_input_res_map_ONE = create_full_to_input_res_map( input_res_vectors[1] );
	std::map< core::Size, core::Size > full_to_input_res_map_TWO = create_full_to_input_res_map( input_res_vectors[2] );

	if ( input_res_vectors.size() != 2 ) return; //This problem occur only when copying two input pose.

	//first find residues that are common between the two input pose1

	utility::vector1< Size > common_res_list;

	for ( Size n = 1; n <= input_res_vectors[1].size(); n++ ){
		Size const seq_num = input_res_vectors[1][n];
		if ( input_res_vectors[2].has_value( seq_num ) == true ){
			common_res_list.push_back( seq_num );
		}
	}

	if ( verbose_ ){
		Output_seq_num_list( "input_res_vectors[1] = ", input_res_vectors[1], TR, 30 );
		Output_seq_num_list( "input_res_vectors[2] = ", input_res_vectors[2], TR, 30 );
		Output_seq_num_list( "common_res_list = ", common_res_list, TR, 30 );
	}

	if ( common_res_list.size() == 0 ) {
		if ( verbose_ ) Output_title_text( "Exit StepWiseRNA_PoseSetup::correctly_copy_HO2star_position", TR.Debug );
		return; //No common/background residues...
	}

	for ( Size n = 1; n <= common_res_list.size(); n++ ){
		Size const full_seq_num = common_res_list[n];
		Size const working_seq_num = full_to_sub[full_seq_num];

		Size const input_ONE_seq_num = full_to_input_res_map_ONE.find( full_seq_num )->second;
		Real const nearest_dist_ONE = get_nearest_dist_to_O2star( input_ONE_seq_num, start_pose_list[1], input_res_vectors[1], common_res_list );

		//SML PHENIX conference cleanup
		if ( basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value() ){
			if ( !start_pose_list[1].residue_type( input_ONE_seq_num ).is_RNA() ) continue;
		}

		Size const input_TWO_seq_num = full_to_input_res_map_TWO.find( full_seq_num )->second;
		Real const nearest_dist_TWO = get_nearest_dist_to_O2star( input_TWO_seq_num, start_pose_list[2], input_res_vectors[2], common_res_list );

		//SML PHENIX conference cleanup
		if ( basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value() ){
			if ( !start_pose_list[2].residue_type( input_TWO_seq_num ).is_RNA() ) continue;
		}

		//pose.set_torsion( TorsionID( moving_res, id::CHI, 4 ), 0 );  //This torsion is not sampled. Arbitary set to zero to prevent randomness
		id::TorsionID const torsion_id( working_seq_num, id::CHI, 4 );

		if ( verbose_ ) TR.Debug << "full_seq_num = " << full_seq_num << " nearest_dist_ONE = " << nearest_dist_ONE << " nearest_dist_TWO = " << nearest_dist_TWO;

		if ( nearest_dist_ONE < nearest_dist_TWO ){
			working_pose.set_torsion( torsion_id, start_pose_list[1].torsion( TorsionID( input_ONE_seq_num, id::CHI, 4 ) ) );
			TR.Debug << " NEARER TO INPUT_POSE_ONE " ;
		} else{
			working_pose.set_torsion( torsion_id, start_pose_list[2].torsion( TorsionID( input_TWO_seq_num, id::CHI, 4 ) ) );
			TR.Debug << " NEARER TO INPUT_POSE_TWO " ;
		}

		TR.Debug << std::endl;

	}

	if ( output_pdb_ ) working_pose.dump_pdb( "copy_dof_pose_AFTER_correctly_copy_HO2star_positions.pdb" );
	if ( verbose_ ) Output_title_text( "Exit StepWiseRNA_PoseSetup::correctly_copy_HO2star_position", TR.Debug );

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

core::Real
StepWiseRNA_PoseSetup::get_nearest_dist_to_O2star( Size const O2star_seq_num, pose::Pose const & input_pose, utility::vector1< Size > const input_res_list, utility::vector1< Size > const & common_res_list ){

	using namespace core::conformation;
	using namespace ObjexxFCL;


	//Real nearest_dist_SQ=9999999999999.99; //Feb 12, 2012 This might lead to server-test error at R47200
	Real nearest_dist_SQ = 99999999.9; //Feb 12, 2012

	conformation::Residue const & input_pose_O2star_rsd = input_pose.residue( O2star_seq_num );

	if ( !input_pose_O2star_rsd.has( " O2'" ) ) utility_exit_with_message( "rsd at input_seq_num = " + string_of( O2star_seq_num ) + " doesn't have O2' atom! " );

	numeric::xyzVector< core::Real > const O2star_xyz = input_pose_O2star_rsd.xyz( input_pose_O2star_rsd.atom_index( " O2'" ) );


	for ( Size input_pose_seq_num = 1; input_pose_seq_num <= input_res_list.size(); input_pose_seq_num++ ){
		Size const full_seq_num = input_res_list[input_pose_seq_num];

		if ( common_res_list.has_value( full_seq_num ) == true ) continue; //A common/background res..not a residue built by SWA.

		conformation::Residue const & input_pose_rsd = input_pose.residue( input_pose_seq_num );

		for ( Size at = 1; at <= input_pose_rsd.natoms(); at++ ){

			Real const dist_to_o2star_SQ = ( input_pose_rsd.xyz( at ) - O2star_xyz ).length_squared();
			if ( dist_to_o2star_SQ < nearest_dist_SQ ) nearest_dist_SQ = dist_to_o2star_SQ;

		}


	}

	Real const nearest_dist = sqrt( nearest_dist_SQ );

	return nearest_dist;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::apply_cutpoint_variants( pose::Pose & pose, pose::Pose & pose_without_cutpoints ){

	using namespace core::id;

	std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
	utility::vector1< Size > const & cutpoint_closed_list = job_parameters_->cutpoint_closed_list();

	for ( Size n = 1; n <= cutpoint_closed_list.size(); n++ ){

		Size const cutpoint_closed = cutpoint_closed_list[n];

		if ( cutpoint_closed == 0 ) continue;

		if ( full_to_sub.find( cutpoint_closed ) != full_to_sub.end() && full_to_sub.find( cutpoint_closed + 1 ) != full_to_sub.end() ) {

			TR.Debug << "Applying cutpoint variants to " << cutpoint_closed << std::endl;

			Size const cutpos = full_to_sub[ cutpoint_closed];

			correctly_add_cutpoint_variants( pose, cutpos, false /*check fold tree*/);

			TR.Debug << "pose ( before copy ): " << std::endl;
			print_backbone_torsions( pose, cutpos );
			print_backbone_torsions( pose, cutpos + 1 );

			for ( Size i = cutpos; i <= cutpos + 1; i++ ){
				for ( Size j = 1; j <= chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; j++ ) {
					id::TorsionID torsion_id( i, id::BB, j );
					pose.set_torsion( torsion_id, pose_without_cutpoints.torsion( torsion_id ) ); //This makes sure that the chain_break torsions have the correct value
				} // j
			} // i

			if ( verbose_ ) {
				TR.Debug << "pose_without_cutpoints " << std::endl;
				print_backbone_torsions( pose_without_cutpoints, cutpos );
				print_backbone_torsions( pose_without_cutpoints, cutpos + 1 );
				TR.Debug << "pose: " << std::endl;
				print_backbone_torsions( pose, cutpos );
				print_backbone_torsions( pose, cutpos + 1 );
			}

		} else {
			utility_exit_with_message( "User provided cutpoint_closed not in working pose?" );
		}
	}

}

//////////////////////////////////////////////////////////////////////////////
// Assume we have done a crappy job of placing 5' phosphates.
//  this would be true if, for example, the 5' phosphate was built
//  by prepending in a previous rebuild-from-scratch effort.
// The only exception is if the 5'-phosphate is involved in *chain closure*.
//(5' phosphate actually refer to the phosphate group of the residue 3' of the chain_break!!! Jan 29, 2010 Parin S.)
void
StepWiseRNA_PoseSetup::apply_virtual_phosphate_variants( pose::Pose & pose ) const{

	using namespace ObjexxFCL;

	utility::vector1< std::pair < core::Size, core::Size > > const & chain_boundaries(	 job_parameters_->chain_boundaries() );
	std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
	utility::vector1< Size > const & cutpoint_closed_list = job_parameters_->cutpoint_closed_list();

	Size const num_chains = chain_boundaries.size();

	for ( Size n = 1; n <= num_chains; n++ ) {
		Size const chain_start = chain_boundaries[ n ].first;

		bool chain_start_is_cutpoint_closed = false;

		for ( Size i = 1; i <= cutpoint_closed_list.size(); i++ ){
			Size const cutpoint_closed = cutpoint_closed_list[i];

			if ( cutpoint_closed > 0   &&  full_to_sub[ chain_start ] == full_to_sub[ cutpoint_closed ] + 1 ) {
				chain_start_is_cutpoint_closed = true;
				break;
			}
		}
		if ( chain_start_is_cutpoint_closed ) continue;

		if ( pose.residue( full_to_sub[ chain_start ] ).type().has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
			utility_exit_with_message( "Should not be trying to virtualize phosphate on close cutpoint residue ( chain_start = " + string_of( chain_start )  + " )" );
		}

		//SML PHENIX conference cleanup
		if ( basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value() ){
			if ( !pose.residue_type( full_to_sub[ chain_start ] ).is_RNA() ) continue;
		}

		pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", full_to_sub[ chain_start ] );
	}
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::add_terminal_res_repulsion( core::pose::Pose & pose ) const
{
	using namespace core::conformation;
	using namespace core::id;
	using namespace core::scoring::constraints;
	using namespace core::chemical::rna;

	ConstraintSetOP cst_set( pose.constraint_set()->clone() );
	assert( cst_set );

	utility::vector1< core::Size > const & working_terminal_res = job_parameters_->working_terminal_res();

	if ( working_terminal_res.size() == 0 ) return;


	/////////////////////////////////////////////////
	Size const nres( pose.total_residue() );

	ObjexxFCL::FArray1D < bool > is_moving_res( nres, false );
	ObjexxFCL::FArray1D < bool > is_fixed_res( nres, false );

	ObjexxFCL::FArray1D < bool > const & partition_definition = job_parameters_->partition_definition();
	bool const root_partition = partition_definition( pose.fold_tree().root() );

	for ( Size seq_num = 1; seq_num <= nres; seq_num++ ){
		if ( partition_definition( seq_num ) == root_partition ) {
			is_fixed_res( seq_num ) = true;
		} else {
			is_moving_res( seq_num ) = true;
		}
	}

	Output_seq_num_list( "working_terminal_res_list = ", working_terminal_res, TR.Debug );
	/////////////////////////////////////////////////
	Distance const DIST_CUTOFF = 8.0;
	FuncOP const repulsion_func( new FadeFunc( -2.0 /*min*/, DIST_CUTOFF /*max*/, 1.0 /*fade zone width*/, 100.0 /*penalty*/ ) );

	for ( Size i = 1; i <= working_terminal_res.size(); i++ ) {

		Size const k = working_terminal_res[ i ];
		Residue const & rsd1( pose.residue( k ) );
		if ( rsd1.has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
			TR.Debug << "rsd1.has_variant_type( \"VIRTUAL_RNA_RESIDUE\" ), seq_num = " << k << " Ignore terminal_res_repulsion distance constraint " << std::endl;
			continue;
		}
		for ( Size m = 1; m <= nres; m++ ) {

			Residue const & rsd2( pose.residue( m ) );
			if ( rsd2.has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
				 TR.Debug << "rsd2.has_variant_type( \"VIRTUAL_RNA_RESIDUE\" ), seq_num = " << m << " Ignore terminal_res_repulsion distance constraint " << std::endl;
				 continue;
			}

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

	/* For debugging....
	TR.Debug << "constraints " << std::endl;
	core::scoring::constraints::ConstraintSetOP cst_set(  pose.constraint_set()->clone() );
	cst_set->show( TR );
	pose.remove_constraints();
	*/
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::add_protonated_H1_adenosine_variants( pose::Pose & pose ) const {

	utility::vector1< core::Size > const & working_protonated_H1_adenosine_list = job_parameters_->working_protonated_H1_adenosine_list();
	Size const working_moving_res =  job_parameters_->working_moving_res(); // corresponds to user input.

	bool apply_check = true;

	if ( rebuild_bulge_mode_ ){
		//OK as long as we will definitely remove the virtual variant type from this res before sampling and minimizing!
		apply_check = false;
		TR.Debug << "rebuild_bulge_mode_ = true, setting apply_check for apply_protonated_H1_adenosine_variant_type to false" << std::endl;
	}

	if ( working_protonated_H1_adenosine_list.has_value( working_moving_res ) ){
		apply_protonated_H1_adenosine_variant_type( pose, working_moving_res, apply_check );
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::verify_protonated_H1_adenosine_variants( pose::Pose & pose ) const {

	using namespace core::conformation;
	using namespace core::pose;
	using namespace ObjexxFCL;


	utility::vector1< core::Size > const & working_protonated_H1_adenosine_list = job_parameters_->working_protonated_H1_adenosine_list();

//Check that all protonated_H1_adenosine exist in the pose!
	for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

		if ( working_protonated_H1_adenosine_list.has_value( seq_num ) ){

			if ( pose.residue( seq_num ).aa() != core::chemical::na_rad ){
				print_JobParameters_info( StepWiseRNA_JobParametersCOP( job_parameters_ ), "DEBUG job_parameters", TR.Debug );
				utility_exit_with_message( "working_protonated_H1_adenosine_list.has_value( seq_num ) == true but pose.residue( seq_num ).aa() != core::chemical::na_rad, seq_num = " + string_of( seq_num ) );
			}

			if ( pose.residue( seq_num ).has_variant_type( "PROTONATED_H1_ADENOSINE" ) == false && pose.residue( seq_num ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) == false ){
				print_JobParameters_info( StepWiseRNA_JobParametersCOP( job_parameters_ ), "DEBUG job_parameters", TR.Debug );
				utility_exit_with_message( "working_protonated_H1_adenosine_list.has_value( seq_num ) == true but residue doesn't either PROTONATED_H1_ADENOSINE or VIRTUAL_RNA_RESIDUE variant type, seq_num = " + string_of( seq_num ) );
			}
		} else{
			if ( pose.residue( seq_num ).has_variant_type( "PROTONATED_H1_ADENOSINE" ) ){

				print_JobParameters_info( StepWiseRNA_JobParametersCOP( job_parameters_ ), "DEBUG job_parameters", TR.Debug );
				TR.Debug << "ERROR: seq_num = " << seq_num << std::endl;
				TR.Debug << "ERROR: start_pose.residue( n ).aa() = " << name_from_aa( pose.residue( seq_num ).aa() ) << std::endl;
				utility_exit_with_message( "working_protonated_H1_adenosine_list.has_value( seq_num ) == false but pose.residue( seq_num ).has_variant_type( \"PROTONATED_H1_ADENOSINE\" ) ) == false" );
			}

		}

	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::apply_bulge_variants( pose::Pose & pose ) const {

	using namespace core::conformation;
	using namespace core::pose;
	using namespace ObjexxFCL;


	std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
	utility::vector1< Size > const & working_terminal_res = job_parameters_->working_terminal_res();

	utility::vector1< Size > terminal_res = apply_sub_to_full_mapping( working_terminal_res, job_parameters_ );

	for ( Size i = 1; i <= bulge_res_.size(); i++ ) {
		Size const seq_num = bulge_res_[i];

		if ( terminal_res.has_value( seq_num ) == true ) utility_exit_with_message( "seq_num: " + string_of( seq_num ) + " cannot be both both a bulge_res and a terminal res!" );

		if ( full_to_sub.find( seq_num ) == full_to_sub.end() ) continue;

		pose::add_variant_type_to_pose_residue( pose, "BULGE", 	full_to_sub[ seq_num] );
	}


}

////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::apply_virtual_res_variant( pose::Pose & pose ) const {
	//Parin Jan 17, 2009

	using namespace core::id;
	using namespace ObjexxFCL;


	std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
	utility::vector1< Size > const & working_terminal_res = job_parameters_->working_terminal_res();

	utility::vector1< Size > terminal_res = apply_sub_to_full_mapping( working_terminal_res, job_parameters_ );

	///////////////////////////////////////////////////////////////////////////////////////
	Size const working_moving_res(  job_parameters_->working_moving_res() ); // corresponds to user input.
	bool const Is_prepend(  job_parameters_->Is_prepend() );
	Size const working_bulge_moving_res = ( Is_prepend ) ? working_moving_res + 1 : working_moving_res - 1;

	bool const Is_dinucleotide = ( job_parameters_->working_moving_res_list().size() == 2 );

	if ( apply_virtual_res_variant_at_dinucleotide_ &&  Is_dinucleotide ){
		if ( working_bulge_moving_res != job_parameters_->working_moving_res_list()[2] ){
			Output_boolean( "Is_prepend = ", job_parameters_->Is_prepend(), TR.Debug );
			TR.Debug << " working_moving_res = " << working_moving_res << std::endl;
			TR.Debug << "working_bulge_moving_res = " << working_bulge_moving_res << " working_moving_res_list()[2] = " << job_parameters_->working_moving_res_list()[2] << std::endl;
			utility_exit_with_message( "working_bulge_moving_res != working_moving_res_list[2]" );
		}

		if ( working_terminal_res.has_value( working_bulge_moving_res ) == true ){
			 utility_exit_with_message( "working_bulge_moving_res cannot be both both a virtual_res and a terminal res!" );
		}
		apply_virtual_rna_residue_variant_type( pose, working_bulge_moving_res );

//			pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESID#UE", working_bulge_moving_res);
	}
	/////////////////////////////////////////////////////////////////////////////////////////

	for ( Size i = 1; i <= virtual_res_list_.size(); i++ ){
		Size const seq_num = virtual_res_list_[i];

		if ( terminal_res.has_value( seq_num ) == true ) utility_exit_with_message( "seq_num: " + string_of( seq_num ) + " cannot be both both a virtual_res and a terminal res!" );

		if ( full_to_sub.find( seq_num ) == full_to_sub.end() ) continue;

		apply_virtual_rna_residue_variant_type( pose, full_to_sub[ seq_num] );

//			pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESID#UE", full_to_sub[ seq_num] );

	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
//Adding Virtual res as root
void
StepWiseRNA_PoseSetup::add_aa_virt_rsd_as_root( core::pose::Pose & pose ){  //Fang's electron density code

	Size const nres = pose.total_residue();
	Size const working_moving_res(  job_parameters_->working_moving_res() );
	//if already rooted on virtual residue, return
	if ( pose.residue( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) {
		TR.Warning << "addVirtualResAsRoot() called but pose is already rooted on a VRT residue ... continuing." << std::endl;
		return;
	}

	core::chemical::ResidueTypeSet const & residue_set = pose.residue_type( 1 ).residue_type_set();
	core::chemical::ResidueTypeCOPs const & rsd_type_list( residue_set.name3_map( "VRT" ) );
	core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *rsd_type_list[1] ) );
	if ( working_moving_res == 1 ) {
		pose.append_residue_by_jump( *new_res, nres );
	} else {
		pose.append_residue_by_jump( *new_res, 1 );
	}

	// make the virt atom the root
	kinematics::FoldTree newF( pose.fold_tree() );
	newF.reorder( nres + 1 );
	pose.fold_tree( newF );
}
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::setup_pdb_info_with_working_residue_numbers( pose::Pose & pose ) const {

	using namespace core::pose;

	utility::vector1< core::Size > const & is_working_res( job_parameters_->is_working_res() );
	std::string const & full_sequence = job_parameters_->full_sequence();

	utility::vector1< Size > working_res;
	for ( Size i = 1; i <= full_sequence.size(); i++ ){
		if ( is_working_res[ i ] ){
			working_res.push_back( i );
		}
	}

	if ( job_parameters_->add_virt_res_as_root() ){
		working_res.push_back( full_sequence.size() + 1  );
	}

	PDBInfoOP pdb_info = new PDBInfo( pose );
	pdb_info->set_numbering( working_res );
	pose.pdb_info( pdb_info );
}


}
}
}
