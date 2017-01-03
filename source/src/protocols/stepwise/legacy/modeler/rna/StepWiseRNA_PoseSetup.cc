// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StepWiseRNA_PoseSetup
/// @brief Sets up pose and job parameters for RNA stepwise building.
/// @details
/// @author Rhiju Das, Parin Sripakdeevong


//////////////////////////////////
#include <protocols/stepwise/legacy/modeler/rna/StepWiseRNA_PoseSetup.hh>
#include <protocols/stepwise/legacy/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/rna/sugar/util.hh>
#include <protocols/stepwise/modeler/output_util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/farna/util.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/copydofs/util.hh>
#include <core/pose/rna/util.hh>
#include <protocols/scoring/VDW_CachedRepScreenInfo.hh>
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <core/io/silent/BinarySilentStruct.hh> //Feb 24, 2011, FARFAR start_silent_file.

#include <utility/exit.hh>
#include <time.h>

#include <string>

using namespace core;
using core::Real;
using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::modeler::rna;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//
// Setup pose for StepWiseRNA Assembly, by Parin Sripakdeevong.
// Some code cleanup by Rhiju (2013).
//
// To do: clean up & generalize protonated A variant.
//        unification with new FullModelInfo?
//
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.legacy.modeler.rna.StepWiseRNA_PoseSetup" );

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace rna {

//////////////////////////////////////////////////////////////////////////
//constructor!
StepWiseRNA_PoseSetup::StepWiseRNA_PoseSetup( stepwise::modeler::working_parameters::StepWiseWorkingParametersOP & working_parameters ):
	rsd_set_( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ),
	working_parameters_( working_parameters ),
	copy_DOF_( false ),
	verbose_( true ),
	output_pdb_( false ), //Sept 24, 2011
	apply_virtual_res_variant_at_dinucleotide_( true ), //Jul 12, 2012
	align_to_native_( true ),
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

	::output_title_text( "Enter StepWiseRNA_PoseSetup::apply", TR.Debug );
	create_starting_pose( pose );
	if ( use_phenix_geo_ ) pose::rna::apply_pucker( pose, working_parameters_->working_moving_res() );
	apply_cutpoint_variants( pose );
	apply_virtual_phosphate_variants( pose );
	add_protonated_H1_adenosine_variants( pose );
	apply_bulge_variants( pose );
	apply_virtual_res_variant( pose );
	verify_protonated_H1_adenosine_variants( pose );
	add_terminal_res_repulsion( pose );
	additional_setup_for_floating_base( pose );
	instantiate_residue_if_rebuilding_bulge( pose );
	update_fold_tree_at_virtual_sugars( pose );

	if ( working_parameters_->add_virt_res_as_root() ) add_aa_virt_rsd_as_root( pose ); //Fang's electron density code.

	// These should be called explicitly by the object (i.e., stepremovewise_rna_pose_setup->setup_full_model_info( pose ) )
	setup_full_model_info( pose );
	setup_vdw_cached_rep_screen_info( pose );

	if ( output_pdb_ ) pose.dump_pdb( "start.pdb" );
	if ( verbose_ ) output_title_text( "Exit StepWiseRNA_PoseSetup::apply", TR.Debug );
}

////////////////////////////////////////////////////////


void
StepWiseRNA_PoseSetup::setup_native_pose( core::pose::Pose & pose ){
	using namespace core::conformation;
	using namespace core::pose;
	using namespace protocols::farna;

	if ( verbose_ ) output_title_text( "Enter StepWiseRNA_PoseSetup::setup_native_pose", TR.Debug );

	utility::vector1< core::Size > const & is_working_res( working_parameters_->is_working_res() );
	std::string const & working_sequence( working_parameters_->working_sequence() );
	utility::vector1< core::Size > const & working_best_alignment( working_parameters_->working_best_alignment() );
	utility::vector1< core::Size > const & working_moving_res_list( working_parameters_->working_moving_res_list() );

	utility::vector1< core::Size > const working_native_virtual_res_list_ = apply_full_to_sub_mapping( native_virtual_res_list_, stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP( working_parameters_ ) );

	utility::vector1 < core::Size > const & working_moving_partition_res = working_parameters_->working_moving_partition_res(); //Sept 14, 2011.

	if ( !get_native_pose() ) return;
	// Need a simple fold tree for following to work...
	Pose native_pose_copy = *get_native_pose();
	native_pose_copy.fold_tree(  core::kinematics::FoldTree(  native_pose_copy.size() ) );

	PoseOP working_native_pose( new Pose );

	//First option is to pass in the full length native
	if ( native_pose_copy.sequence() == core::pose::rna::remove_bracketed( working_parameters_->full_sequence() ) ) {
		TR.Debug << "User passed in full length native pose" << std::endl;

		( *working_native_pose ) = native_pose_copy;

		for ( Size seq_num = native_pose_copy.size(); seq_num >= 1 ; seq_num-- ) {
			if ( !is_working_res[ seq_num ] ) {
				working_native_pose->conformation().delete_residue_slow( seq_num );
			}
		}

	} else if ( native_pose_copy.sequence() == core::pose::rna::remove_bracketed( working_sequence ) ) { //Could also pass in the working_native_pose directly
		TR.Debug << "User passed in working native pose" << std::endl;
		( *working_native_pose ) = native_pose_copy;
	} else {
		TR.Debug <<  std::setw( 50 ) << "  native_pose_copy.sequence() = " << native_pose_copy.sequence()  << std::endl;
		TR.Debug <<  std::setw( 50 ) << "  working_parameters_->full_sequence() = " << working_parameters_->full_sequence() << std::endl;
		TR.Debug <<  std::setw( 50 ) << "  working_sequence = " << core::pose::rna::remove_bracketed( working_sequence ) << std::endl;
		utility_exit_with_message( "The native pose passed in by the User does not match both the full_sequence and the working sequence of the inputted fasta_file" );
	}

	if ( working_native_pose->sequence() != core::pose::rna::remove_bracketed( working_sequence ) ) {
		TR.Debug <<  std::setw( 50 ) << "working_native_pose->sequence() = " << working_native_pose->sequence();
		TR.Debug <<  std::setw( 50 ) << "working_sequence = " << core::pose::rna::remove_bracketed( working_sequence ) << std::endl;
		utility_exit_with_message( "working_native_pose->sequence() !=  working_sequence" );
	}

	core::pose::rna::make_phosphate_nomenclature_matches_mini( ( *working_native_pose ) );

	utility::vector1< core::Size > act_working_alignment;
	act_working_alignment.clear();

	utility::vector1< core::Size > const & working_native_alignment = working_parameters_->working_native_alignment();

	if ( working_native_alignment.size() != 0 ) { //User specified the native alignment res.
		TR.Debug << "working_native_alignment.size() != 0!, align native_pose with working_native_alignment!" << std::endl;

		for ( Size n = 1; n <= working_native_alignment.size(); n++ ) {
			Size const seq_num = working_native_alignment[n];
			if ( working_moving_res_list.has_value( seq_num ) ) continue;
			if ( working_moving_partition_res.has_value( seq_num ) ) continue; //Sept 14, 2011.
			act_working_alignment.push_back( seq_num );
		}


	} else { //Use default alignment res.
		TR.Debug << "working_native_alignment.size() == 0!, align native_pose with working_best_alignment!" << std::endl;

		for ( Size n = 1; n <= working_best_alignment.size(); n++ ) {
			Size const seq_num = working_best_alignment[n];
			if ( working_moving_res_list.has_value( seq_num ) ) continue;
			if ( working_moving_partition_res.has_value( seq_num ) ) continue; //Sept 14, 2011.
			act_working_alignment.push_back( seq_num );
		}
	}

	output_seq_num_list( "act_working_alignment = ", act_working_alignment, TR.Debug );
	output_seq_num_list( "working_moving_res_list = ", working_moving_res_list, TR.Debug );
	output_seq_num_list( "working_moving_partition_res = ", working_moving_partition_res, TR.Debug );
	output_seq_num_list( "working_native_alignment = ", working_native_alignment, TR.Debug );
	output_seq_num_list( "working_best_alignment = ", working_best_alignment, TR.Debug );

	if ( act_working_alignment.size() == 0 ) utility_exit_with_message( "act_working_alignment.size() == 0" );

	if ( working_parameters_->add_virt_res_as_root() ) { //Fang's electron density code

		//Fang why do you need this?////
		pose::Pose dummy_pose = *working_native_pose;
		apply_cutpoint_variants( *working_native_pose );
		add_protonated_H1_adenosine_variants( *working_native_pose );
		verify_protonated_H1_adenosine_variants( *working_native_pose );
		////////////////////////////////

		add_aa_virt_rsd_as_root( *working_native_pose );
		( *working_native_pose ).fold_tree( pose.fold_tree() );
		align_poses( pose, "working_pose", ( *working_native_pose ), "working_native_pose", act_working_alignment );

	} else if ( align_to_native_ ) {
		// Rhiju -- prefer to align pose to native pose...
		align_poses( pose, "working_pose", ( *working_native_pose ), "working_native_pose", act_working_alignment );

	} else { //March 17, 2012 (Previous) Classic. This ensures that adding native_pdb does not change the final output silent_struct coordinates.
		align_poses( ( *working_native_pose ), "working_native_pose", pose, "working_pose", act_working_alignment );

	}
	setup_full_model_info( *working_native_pose );
	working_parameters_->set_working_native_pose( working_native_pose );

	//Setup cutpoints, chain_breaks, variant types here!

	if ( output_pdb_ ) working_native_pose->dump_pdb( "working_native_pose.pdb" );
	if ( output_pdb_ ) pose.dump_pdb( "modeler_pose.pdb" );

	for ( Size n = 1; n <= working_native_virtual_res_list_.size(); n++ ) {
		//   pose::add_variant_type_to_pose_residue( (*working_native_pose), chemical::VIRTUAL_RNA_RESIDUE, working_native_virtual_res_list_[n]);
		core::pose::rna::apply_virtual_rna_residue_variant_type( ( *working_native_pose ), working_native_virtual_res_list_[n], false /*apply_check*/ ) ;
	}

	if ( output_pdb_ ) working_native_pose->dump_pdb( "working_native_pose_with_virtual_res_variant_type.pdb" );
	if ( verbose_ ) output_title_text( "Exit StepWiseRNA_PoseSetup::setup_native_pose", TR.Debug );
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
		if ( found == std::string::npos ) pose_name.append( ".pdb" );

		//  if(verbose) TR.Debug << " The following pose will be imported :" << pose_name << std::endl;
		core::chemical::ResidueTypeSetCOP rsd_set( rsd_set_ );
		import_pose::pose_from_file( import_pose, *rsd_set, pose_name, core::import_pose::PDB_file );
		core::pose::rna::make_phosphate_nomenclature_matches_mini( import_pose );
	} else {
		import_pose_from_silent_file( import_pose, silent_files_in_[ i ], input_tags_[i] );
	}

	///////////////////////////////
	// Check for sequence match.
	utility::vector1< Size > const & input_res = working_parameters_->input_res_vectors()[i];
	std::string const & full_sequence = working_parameters_->full_sequence();
	std::string stripped_full_sequence = core::pose::rna::remove_bracketed( full_sequence );
	runtime_assert ( import_pose.size() == input_res.size() );
	bool match( true );
	for ( Size n = 1; n <= import_pose.size(); n++ ) {
		if (  import_pose.sequence()[ n - 1 ]  != stripped_full_sequence[ input_res[n] - 1 ] ) {
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
StepWiseRNA_PoseSetup::make_extended_pose( pose::Pose & pose ){
	core::chemical::ResidueTypeSetCOP rsd_set( rsd_set_ );
	make_pose_from_sequence( pose, working_parameters_->working_sequence(),
		*rsd_set, false /*auto_termini*/ );
	if ( output_pdb_ ) {
		TR.Debug << "outputting extended_chain.pdb" << std::endl;
		pose.dump_pdb( "extended_chain.pdb" );
	}
}


///////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::create_starting_pose( pose::Pose & pose ){
	if ( copy_DOF_ ) {
		create_pose_from_input_poses( pose );
	} else {
		pose.fold_tree( working_parameters_->fold_tree() );
	}
}

///////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::create_pose_from_input_poses( pose::Pose & pose )
{
	make_extended_pose( pose );
	pose.fold_tree( working_parameters_->fold_tree() );
	TR.Debug << "read_input_pose_and_copy_dofs( pose )" << std::endl;
	read_input_pose_and_copy_dofs( pose );
}

///////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::read_input_pose_and_copy_dofs( pose::Pose & pose )
{
	using namespace core::pose;
	using namespace core::pose::rna;
	using namespace ObjexxFCL;

	runtime_assert( input_tags_.size() >= 1 );
	runtime_assert( copy_DOF_ );

	utility::vector1< utility::vector1< Size > > const & input_res_vectors = working_parameters_->input_res_vectors();
	std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );
	utility::vector1< pose::Pose > start_pose_list;

	Pose start_pose;
	Pose start_pose_with_variant;
	for ( Size i = 1; i <= input_tags_.size(); i++ ) {

		TR.Debug << "import_pose " << i << std::endl;
		//Ok account for special case of build loop outward from scratch...
		if ( input_tags_[i] == "build_from_scratch" ) {
			if ( input_res_vectors[i].size() != 1 ) utility_exit_with_message( "input_tags_[i] == \"build_from_scratch\" but input_res_vectors[i].size() != 1" );
			continue;
		}

		Import_pose( i, start_pose );
		start_pose_with_variant = start_pose;
		start_pose_list.push_back( start_pose );

		utility::vector1< Size > const & input_res = working_parameters_->input_res_vectors()[i]; //This is from input_pose numbering to full_pose numbering

		// Remove all variant types
		// TO DO LIST: Need to remove atom constraint and remove angle constaint as well
		//NOTES: June 16, 2011
		//Should LOWER_TERMINUS and UPPER_TERMINUS be removed as well? LOWER_TERMINUS does determine the position of  OP2 and OP1?
		//Also should then check that pose.residue_type(i).variant_types() is the empty?
		utility::vector1< core::chemical::VariantType > variant_type_list;
		variant_type_list.push_back( core::chemical::VIRTUAL_PHOSPHATE );
		variant_type_list.push_back( core::chemical::VIRTUAL_O2PRIME_HYDROGEN );
		variant_type_list.push_back( core::chemical::CUTPOINT_LOWER );
		variant_type_list.push_back( core::chemical::CUTPOINT_UPPER );
		variant_type_list.push_back( core::chemical::VIRTUAL_RNA_RESIDUE );
		variant_type_list.push_back( core::chemical::BULGE );
		variant_type_list.push_back( core::chemical::VIRTUAL_RIBOSE );
		variant_type_list.push_back( core::chemical::PROTONATED_N1_ADENOSINE );
		variant_type_list.push_back( core::chemical::THREE_PRIME_END_OH );     //Fang's electron density code
		variant_type_list.push_back( core::chemical::FIVE_PRIME_END_PHOSPHATE ); //Fang's electron density code
		variant_type_list.push_back( core::chemical::FIVE_PRIME_END_OH );     //Fang's electron density code
		for ( Size seq_num = 1; seq_num <= start_pose.size(); seq_num++  ) {
			for ( Size k = 1; k <= variant_type_list.size(); ++k ) {
				core::chemical::VariantType const variant_type = variant_type_list[k];
				if ( start_pose.residue( seq_num ).has_variant_type( variant_type ) ) {
					remove_variant_type_from_pose_residue( start_pose, variant_type, seq_num );
				}
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Now actually copy into the pose.
		// Need to know correspondence of residues between imported pose and the working pose...
		runtime_assert ( input_res.size() == start_pose.size() );
		std::map< core::Size, core::Size > res_map;  //This is map from sub numbering to input_res numbering..
		for ( Size n = 1; n <= input_res.size(); n++ )  res_map[ full_to_sub[ input_res[n] ] ] = n;

		//Does this work for the "overlap residue" case?? If there is a overlap residue, then order of input_res will manner...Parin Jan 2, 2010.
		core::pose::copydofs::copy_dofs_match_atom_names( pose, start_pose, res_map, false /*backbone_only*/, false /*side_chain_only*/, false /*ignore_virtual*/ );

		//////////////////////Add virtual_rna_residue, virtual_rna_residue variant_upper and VIRTUAL_RIBOSE variant type back//////////////////////
		//////////////////////For standard dinucleotide move, there should be virtual_rna_residue and virtual_rna_residue variant_upper //////////
		//////////////////////And if floating base modeler, will contain the virtual_sugar as well...
		utility::vector1< Size > const & cutpoint_closed_list = working_parameters_->cutpoint_closed_list();
		utility::vector1< Size > const & working_cutpoint_closed_list = apply_full_to_sub_mapping( cutpoint_closed_list, stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP( working_parameters_ ) );

		for ( Size n = 1; n <= start_pose_with_variant.size(); n++  ) {
			if ( has_virtual_rna_residue_variant_type( start_pose_with_variant, n ) ) {
				apply_virtual_rna_residue_variant_type( pose, full_to_sub[ input_res[n] ], working_cutpoint_closed_list ) ;
			}
			if ( start_pose_with_variant.residue( n ).has_variant_type( core::chemical::VIRTUAL_RIBOSE ) ) {
				add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_RIBOSE, full_to_sub[ input_res[n] ] );
			}
			if ( start_pose_with_variant.residue( n ).has_variant_type( core::chemical::THREE_PRIME_END_OH ) ) { // are these still in use?
				add_variant_type_to_pose_residue( pose, core::chemical::THREE_PRIME_END_OH, full_to_sub[ input_res[n] ] );
			} else if ( start_pose_with_variant.residue( n ).has_variant_type( core::chemical::FIVE_PRIME_END_PHOSPHATE ) ) {
				add_variant_type_to_pose_residue( pose, core::chemical::FIVE_PRIME_END_PHOSPHATE, full_to_sub[ input_res[n] ] );
			} else if ( start_pose_with_variant.residue( n ).has_variant_type( core::chemical::FIVE_PRIME_END_OH ) ) {
				add_variant_type_to_pose_residue( pose, core::chemical::FIVE_PRIME_END_OH, full_to_sub[ input_res[n] ] );
			}
			do_checks_and_apply_protonated_H1_adenosine_variant( pose, start_pose_with_variant, n, i, input_res, full_to_sub );
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if ( output_pdb_ ) pose.dump_pdb( "copy_dof_" + string_of( i ) + ".pdb" );
		TR.Debug << pose.fold_tree() << std::endl;
	}

	protocols::farna::assert_phosphate_nomenclature_matches_mini( pose ); //Just to be safe, Jun 11, 2010
	correctly_copy_HO2prime_positions( pose, start_pose_list );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::do_checks_and_apply_protonated_H1_adenosine_variant( pose::Pose & pose,
	pose::Pose const & start_pose_with_variant,
	Size const & n /*res num*/,
	Size const & i /*input pose num*/,
	utility::vector1< Size > const & input_res,
	std::map< Size, Size > & full_to_sub ){

	if ( i > silent_files_in_.size() ) { // not a silent file, read in from pdb text file. May 04, 2011
		//start_pose_with_variant does not have PROTONATED_N1_ADENOSINE variant type since the input_pdb does not have the variant type or loses the variant when imported into Rosetta.
		runtime_assert( !start_pose_with_variant.residue( n ).has_variant_type( core::chemical::PROTONATED_N1_ADENOSINE ) ); //May 03, 2011
	} else { // from silent file -- may have information on variants.
		if ( start_pose_with_variant.residue( n ).has_variant_type( core::chemical::PROTONATED_N1_ADENOSINE ) ) { //May 03, 2011
			runtime_assert( start_pose_with_variant.residue( n ).aa() == core::chemical::na_rad );
			runtime_assert( working_parameters_->protonated_H1_adenosine_list().has_value( input_res[n] ) );
		}
	}

	if ( ( working_parameters_->protonated_H1_adenosine_list() ).has_value( input_res[n] ) ) {
		apply_protonated_H1_adenosine_variant_type( pose, full_to_sub[ input_res[n] ] );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// For 2' hydroxyls shared between input poses, they might be held in place by hydrogen bonds defined in one input pose,
//  but not the other. So set those 2' hydroxyl torsions based on which (non-shared) residues they are closest to.
//  [This is of course a heuristic.]
//
void
StepWiseRNA_PoseSetup::correctly_copy_HO2prime_positions( pose::Pose & working_pose, utility::vector1< pose::Pose > const & start_pose_list ){

	using namespace ObjexxFCL;
	using namespace core::conformation;
	using namespace core::id;

	if ( verbose_ ) output_title_text( "Enter StepWiseRNA_PoseSetup::correctly_copy_HO2prime_position", TR.Debug );
	if ( output_pdb_ ) working_pose.dump_pdb( "copy_dof_pose_BEFORE_correctly_copy_HO2prime_positions.pdb" );

	std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );

	//Problem aside in the long_loop mode, where we combine two silent_file pose. The two pose both contain the fixed background residues
	//In current implementation of COPY_DOFS, the positioning of the background residues in SECOND silent file is used.
	//However, some background HO2prime interacts with the first pose and hence need to correctly copys these HO2prime position.
	//When we implemenet protien/RNA interactions, will extend this to include all "packed" side chains.
	utility::vector1< utility::vector1< Size > > const & input_res_vectors = working_parameters_->input_res_vectors();
	std::map< core::Size, core::Size > full_to_input_res_map_ONE = create_full_to_input_res_map( input_res_vectors[1] );
	std::map< core::Size, core::Size > full_to_input_res_map_TWO = create_full_to_input_res_map( input_res_vectors[2] );
	if ( input_res_vectors.size() != 2 ) return; //This problem occurs only when copying two input pose.

	//first find residues that are common between the two input pose1
	utility::vector1< Size > common_res_list;
	for ( Size n = 1; n <= input_res_vectors[1].size(); n++ ) {
		Size const seq_num = input_res_vectors[1][n];
		if ( input_res_vectors[2].has_value( seq_num ) == true ) {
			common_res_list.push_back( seq_num );
		}
	}

	if ( verbose_ ) {
		output_seq_num_list( "input_res_vectors[1] = ", input_res_vectors[1], TR, 30 );
		output_seq_num_list( "input_res_vectors[2] = ", input_res_vectors[2], TR, 30 );
		output_seq_num_list( "common_res_list = ", common_res_list, TR, 30 );
	}
	if ( common_res_list.size() == 0 ) {
		if ( verbose_ ) output_title_text( "Exit StepWiseRNA_PoseSetup::correctly_copy_HO2prime_position", TR.Debug );
		return; //No common/background residues...
	}

	for ( Size const full_seq_num : common_res_list ) {
		Size const working_seq_num = full_to_sub[full_seq_num];

		Size const input_ONE_seq_num = full_to_input_res_map_ONE.find( full_seq_num )->second;
		Real const nearest_dist_ONE = get_nearest_dist_to_O2prime( input_ONE_seq_num, start_pose_list[1], input_res_vectors[1], common_res_list );
		if ( !start_pose_list[1].residue_type( input_ONE_seq_num ).is_RNA() ) continue;

		Size const input_TWO_seq_num = full_to_input_res_map_TWO.find( full_seq_num )->second;
		Real const nearest_dist_TWO = get_nearest_dist_to_O2prime( input_TWO_seq_num, start_pose_list[2], input_res_vectors[2], common_res_list );
		if ( !start_pose_list[2].residue_type( input_TWO_seq_num ).is_RNA() ) continue;

		id::TorsionID const torsion_id( working_seq_num, id::CHI, 4 );
		if ( verbose_ ) TR.Debug << "full_seq_num = " << full_seq_num << " nearest_dist_ONE = " << nearest_dist_ONE << " nearest_dist_TWO = " << nearest_dist_TWO;
		if ( nearest_dist_ONE < nearest_dist_TWO ) {
			working_pose.set_torsion( torsion_id, start_pose_list[1].torsion( TorsionID( input_ONE_seq_num, id::CHI, 4 ) ) );
			TR.Debug << " NEARER TO INPUT_POSE_ONE " ;
		} else {
			working_pose.set_torsion( torsion_id, start_pose_list[2].torsion( TorsionID( input_TWO_seq_num, id::CHI, 4 ) ) );
			TR.Debug << " NEARER TO INPUT_POSE_TWO " ;
		}

		TR.Debug << std::endl;
	}

	if ( output_pdb_ ) working_pose.dump_pdb( "copy_dof_pose_AFTER_correctly_copy_HO2prime_positions.pdb" );
	if ( verbose_ ) output_title_text( "Exit StepWiseRNA_PoseSetup::correctly_copy_HO2prime_position", TR.Debug );

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
core::Real
StepWiseRNA_PoseSetup::get_nearest_dist_to_O2prime( Size const O2prime_seq_num,
	pose::Pose const & input_pose,
	utility::vector1< Size > const & input_res_list,
	utility::vector1< Size > const & common_res_list ){

	using namespace core::conformation;
	using namespace ObjexxFCL;

	Real nearest_dist_SQ = 99999999.9; //Feb 12, 2012
	conformation::Residue const & input_pose_O2prime_rsd = input_pose.residue( O2prime_seq_num );
	runtime_assert( input_pose_O2prime_rsd.has( " O2'" ) );
	numeric::xyzVector< core::Real > const O2prime_xyz = input_pose_O2prime_rsd.xyz( input_pose_O2prime_rsd.atom_index( " O2'" ) );

	for ( Size input_pose_seq_num = 1; input_pose_seq_num <= input_res_list.size(); input_pose_seq_num++ ) {
		Size const full_seq_num = input_res_list[input_pose_seq_num];

		if ( common_res_list.has_value( full_seq_num ) ) continue; //A common/background res..not a residue built by SWA.

		conformation::Residue const & input_pose_rsd = input_pose.residue( input_pose_seq_num );

		for ( Size at = 1; at <= input_pose_rsd.natoms(); at++ ) {
			Real const dist_to_o2prime_SQ = ( input_pose_rsd.xyz( at ) - O2prime_xyz ).length_squared();
			if ( dist_to_o2prime_SQ < nearest_dist_SQ ) nearest_dist_SQ = dist_to_o2prime_SQ;
		}
	}

	Real const nearest_dist = sqrt( nearest_dist_SQ );

	return nearest_dist;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::apply_cutpoint_variants( pose::Pose & pose ){

	using namespace core::id;
	std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );
	utility::vector1< Size > const & cutpoint_closed_list = working_parameters_->cutpoint_closed_list();

	// this copy contains original torsions across chainbreak [zeta, alpha, beta].
	Pose pose_without_cutpoints = pose;
	kinematics::FoldTree simple_fold_tree( pose.size() );
	pose_without_cutpoints.fold_tree( simple_fold_tree ); // permits read out of original torsions across cutpoints. I think.

	for ( Size const cutpoint_closed : cutpoint_closed_list ) {
		if ( cutpoint_closed == 0 ) continue;

		if ( full_to_sub.find( cutpoint_closed )     != full_to_sub.end() &&
				full_to_sub.find( cutpoint_closed + 1 ) != full_to_sub.end() ) {

			TR.Debug << "Applying cutpoint variants to " << cutpoint_closed << std::endl;
			Size const cutpos = full_to_sub[ cutpoint_closed];
			core::pose::correctly_add_cutpoint_variants( pose, cutpos, false /*check fold tree*/);

			TR.Debug << "pose ( before copy ): " << std::endl;
			print_backbone_torsions( pose, cutpos );
			print_backbone_torsions( pose, cutpos + 1 );
			for ( Size i = cutpos; i <= cutpos + 1; i++ ) {
				for ( Size j = 1; j <= chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; j++ ) {
					id::TorsionID torsion_id( i, id::BB, j );
					pose.set_torsion( torsion_id, pose_without_cutpoints.torsion( torsion_id ) ); //This makes sure that the chain_break torsions has the correct value
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

	utility::vector1< std::pair < core::Size, core::Size > > const & chain_boundaries(  working_parameters_->chain_boundaries() );
	std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );
	utility::vector1< Size > const & cutpoint_closed_list = working_parameters_->cutpoint_closed_list();

	Size const num_chains = chain_boundaries.size();

	for ( Size n = 1; n <= num_chains; n++ ) {
		Size const chain_start = chain_boundaries[ n ].first;

		bool chain_start_is_cutpoint_closed = false;

		for ( Size i = 1; i <= cutpoint_closed_list.size(); i++ ) {
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
		if ( !pose.residue_type( full_to_sub[ chain_start ] ).is_RNA() ) continue;

		pose::add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, full_to_sub[ chain_start ] );
	}
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::add_terminal_res_repulsion( core::pose::Pose & pose ) const
{
	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::id;
	using namespace core::scoring::constraints;
	using namespace core::chemical::rna;

	ConstraintSetOP cst_set( pose.constraint_set()->clone() );
	assert( cst_set );

	utility::vector1< core::Size > const & working_terminal_res = working_parameters_->working_terminal_res();

	if ( working_terminal_res.size() == 0 ) return;


	/////////////////////////////////////////////////
	Size const nres( pose.size() );
	ObjexxFCL::FArray1D < bool > is_moving_res( nres, false );
	ObjexxFCL::FArray1D < bool > is_fixed_res( nres, false );

	ObjexxFCL::FArray1D < bool > const & partition_definition = working_parameters_->partition_definition();
	bool const root_partition = partition_definition( pose.fold_tree().root() );

	for ( Size seq_num = 1; seq_num <= nres; seq_num++ ) {
		if ( partition_definition( seq_num ) == root_partition ) {
			is_fixed_res( seq_num ) = true;
		} else {
			is_moving_res( seq_num ) = true;
		}
	}

	output_seq_num_list( "working_terminal_res_list = ", working_terminal_res, TR.Debug );
	/////////////////////////////////////////////////
	Distance const DIST_CUTOFF = 8.0;
	core::scoring::func::FuncOP const repulsion_func( new core::scoring::func::FadeFunc( -2.0 /*min*/, DIST_CUTOFF /*max*/, 1.0 /*fade zone width*/, 100.0 /*penalty*/ ) );

	for ( Size const k : working_terminal_res ) {
		ResidueType const & rsd1( pose.residue_type( k ) );
		if ( rsd1.has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
			TR.Debug << "rsd1.has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ), seq_num = " << k << " Ignore terminal_res_repulsion distance constraint " << std::endl;
			continue;
		}
		for ( Size m = 1; m <= nres; m++ ) {

			ResidueType const & rsd2( pose.residue_type( m ) );
			if ( rsd2.has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
				TR.Debug << "rsd2.has_variant_type( chemical::VIRTUAL_RNA_RESIDUE ), seq_num = " << m << " Ignore terminal_res_repulsion distance constraint " << std::endl;
				continue;
			}

			AtomID const atom_id1( first_base_atom_index( rsd1 ), k );
			AtomID const atom_id2( first_base_atom_index( rsd2 ), m );

			// the one exception -- close contacts WITHIN a partition
			if ( ( ( is_moving_res( k )  && is_moving_res( m ) ) ||
					( is_fixed_res( k )  && is_fixed_res(  m ) ) ) &&
					( pose.xyz( atom_id1 ) - pose.xyz( atom_id2 ) ).length() < DIST_CUTOFF ) {
				//TR.Debug << "Not adding repulsive constraint between " << k << " and " << m << " already closeby in same partition" << std::endl;
				continue;
			}

			// distance from O3' to P
			cst_set->add_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint( atom_id1, atom_id2, repulsion_func ) ) ) );
		}
	}

	pose.constraint_set( cst_set );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::additional_setup_for_floating_base( pose::Pose & pose ) const {

	Size const num_nucleotides = working_parameters_->working_moving_res_list().size();
	bool const is_dinucleotide = ( num_nucleotides == 2 );

	if ( !working_parameters_->floating_base() ) return;

	virtualize_sugar_and_backbone_at_moving_res( pose );

	if ( !is_dinucleotide ) return;

	Size const moving_res_ =  working_parameters_->working_moving_res(); // Might not correspond to user input.
	Size const reference_res_ =  working_parameters_->working_reference_res(); //the last static_residues that this attach to the moving residues
	if ( num_nucleotides == 1 ) {
		setup_chain_break_jump_point( pose, moving_res_, reference_res_ );
	}
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::instantiate_residue_if_rebuilding_bulge( pose::Pose & pose ){
	if ( working_parameters_->rebuild_bulge_mode() ) core::pose::rna::remove_virtual_rna_residue_variant_type( pose, working_parameters_->working_moving_res() );
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::virtualize_sugar_and_backbone_at_moving_res( pose::Pose & pose ) const {
	////////Only base is instantiated in 'floating base' at moving res /////////////////////////////
	Size const moving_res = working_parameters_->working_moving_res();
	pose::add_variant_type_to_pose_residue( pose,
		core::chemical::VIRTUAL_PHOSPHATE, moving_res ); //This is unique to floating_base_mode
	pose::add_variant_type_to_pose_residue( pose,
		core::chemical::VIRTUAL_RIBOSE, moving_res); //This is unique to floating_base_mode
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::add_protonated_H1_adenosine_variants( pose::Pose & pose ) const {

	utility::vector1< core::Size > const & working_protonated_H1_adenosine_list = working_parameters_->working_protonated_H1_adenosine_list();
	Size const working_moving_res =  working_parameters_->working_moving_res(); // corresponds to user input.

	bool apply_check = true;
	if ( working_parameters_->rebuild_bulge_mode() ) {
		//OK as long as we will definitely remove the virtual variant type from this res before modeler and minimizing!
		apply_check = false;
		TR.Debug << "rebuild_bulge_mode_ = true, setting apply_check for apply_protonated_H1_adenosine_variant_type to false" << std::endl;
	}
	if ( working_protonated_H1_adenosine_list.has_value( working_moving_res ) ) {
		apply_protonated_H1_adenosine_variant_type( pose, working_moving_res, apply_check );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::verify_protonated_H1_adenosine_variants( pose::Pose & pose ) const {

	using namespace core::conformation;
	using namespace core::pose;
	using namespace ObjexxFCL;


	utility::vector1< core::Size > const & working_protonated_H1_adenosine_list = working_parameters_->working_protonated_H1_adenosine_list();

	//Check that all protonated_H1_adenosine exist in the pose!
	for ( Size seq_num = 1; seq_num <= pose.size(); seq_num++ ) {
		if ( working_protonated_H1_adenosine_list.has_value( seq_num ) ) {

			if ( pose.residue_type( seq_num ).aa() != core::chemical::na_rad ) {
				print_WorkingParameters_info( stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP( working_parameters_ ), "DEBUG working_parameters", TR.Debug );
				utility_exit_with_message( "working_protonated_H1_adenosine_list.has_value( seq_num ) == true but pose.residue( seq_num ).aa() != core::chemical::na_rad, seq_num = " + string_of( seq_num ) );
			}

			if ( ! pose.residue_type( seq_num ).has_variant_type( core::chemical::PROTONATED_N1_ADENOSINE ) &&
					! pose.residue_type( seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
				print_WorkingParameters_info( stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP(
					working_parameters_ ), "DEBUG working_parameters", TR.Debug );
				utility_exit_with_message( "working_protonated_H1_adenosine_list.has_value( seq_num ) == true but "
					"residue doesn't either PROTONATED_N1_ADENOSINE or VIRTUAL_RNA_RESIDUE variant type, seq_num = " + string_of( seq_num ) );
			}
		} else if ( pose.residue_type( seq_num ).has_variant_type( core::chemical::PROTONATED_N1_ADENOSINE ) ) {
			print_WorkingParameters_info( stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP( working_parameters_ ), "DEBUG working_parameters", TR.Debug );
			TR.Debug << "ERROR: seq_num = " << seq_num << std::endl;
			TR.Debug << "ERROR: start_pose.residue( n ).aa() = " << name_from_aa( pose.residue_type( seq_num ).aa() ) << std::endl;
			utility_exit_with_message( "working_protonated_H1_adenosine_list.has_value( seq_num ) == false but pose.residue( seq_num ).has_variant_type( \"PROTONATED_N1_ADENOSINE\" ) ) == false" );
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Virtual sugars (VIRTUAL_RIBOSE variant type) should only occur in SWA if the residue was created via a
// 'floating base' move, with a jump to an anchor (reference) residue in the rest of the pose.
// Let's update the fold tree, to include these jumps in the workign pose.
// This information can be used in, e.g., VirtualSugarSampler.
void
StepWiseRNA_PoseSetup::update_fold_tree_at_virtual_sugars( pose::Pose & pose ){

	std::map< Size, Size > const & reference_res_for_each_virtual_sugar = sugar::get_reference_res_for_each_virtual_sugar_without_fold_tree( pose,  working_parameters_->working_moving_suite() );

	TR.Debug << "BEFORE VIRTUAL SUGAR UPDATE " << pose.fold_tree() << std::endl;
	for ( auto const & sugar_and_ref : reference_res_for_each_virtual_sugar ) {
		Size const virtual_sugar_res = sugar_and_ref.first;
		Size const reference_res = sugar_and_ref.second;

		if ( pose.fold_tree().jump_exists( virtual_sugar_res, reference_res ) ) continue;
		setup_chain_break_jump_point( pose, virtual_sugar_res, reference_res );
	}
	TR.Debug << "AFTER VIRTUAL SUGAR UPDATE " << pose.fold_tree() << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::apply_bulge_variants( pose::Pose & pose ) const {

	using namespace core::conformation;
	using namespace core::pose;
	using namespace ObjexxFCL;

	std::map< core::Size, core::Size > const & full_to_sub( working_parameters_->full_to_sub() );
	utility::vector1< Size > const & working_terminal_res = working_parameters_->working_terminal_res();
	utility::vector1< Size > terminal_res = apply_sub_to_full_mapping( working_terminal_res, working_parameters_ );

	for ( Size const seq_num : bulge_res_ ) {
		runtime_assert ( !terminal_res.has_value( seq_num ) );
		if ( full_to_sub.find( seq_num ) == full_to_sub.end() ) continue;
		pose::add_variant_type_to_pose_residue( pose, core::chemical::BULGE, full_to_sub.find( seq_num )->second );
	}
}

////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::apply_virtual_res_variant( pose::Pose & pose ) const {
	//Parin Jan 17, 2009

	using namespace core::id;
	using namespace ObjexxFCL;

	std::map< core::Size, core::Size > & full_to_sub( working_parameters_->full_to_sub() );
	utility::vector1< Size > const & working_terminal_res = working_parameters_->working_terminal_res();

	utility::vector1< Size > terminal_res = apply_sub_to_full_mapping( working_terminal_res, working_parameters_ );

	///////////////////////////////////////////////////////////////////////////////////////
	Size const working_moving_res(  working_parameters_->working_moving_res() ); // corresponds to user input.
	bool const is_prepend(  working_parameters_->is_prepend() );
	Size const working_bulge_moving_res = ( is_prepend ) ? working_moving_res + 1 : working_moving_res - 1;

	bool const is_dinucleotide = ( working_parameters_->working_moving_res_list().size() == 2 );

	if ( apply_virtual_res_variant_at_dinucleotide_ &&  is_dinucleotide ) {
		if ( working_bulge_moving_res != working_parameters_->working_moving_res_list()[2] ) {
			output_boolean( "is_prepend = ", working_parameters_->is_prepend(), TR.Debug );
			TR.Debug << " working_moving_res = " << working_moving_res << std::endl;
			TR.Debug << "working_bulge_moving_res = " << working_bulge_moving_res << " working_moving_res_list()[2] = " << working_parameters_->working_moving_res_list()[2] << std::endl;
			utility_exit_with_message( "working_bulge_moving_res != working_moving_res_list[2]" );
		}

		if ( working_terminal_res.has_value( working_bulge_moving_res ) == true ) {
			utility_exit_with_message( "working_bulge_moving_res cannot be both both a virtual_res and a terminal res!" );
		}
		core::pose::rna::apply_virtual_rna_residue_variant_type( pose, working_bulge_moving_res );

		//   pose::add_variant_type_to_pose_residue( pose, chemical::VIRTUAL_RNA_RESIDUE, working_bulge_moving_res);
	}
	/////////////////////////////////////////////////////////////////////////////////////////

	for ( Size const seq_num : virtual_res_list_ ) {
		if ( terminal_res.has_value( seq_num ) ) utility_exit_with_message( "seq_num: " + string_of( seq_num ) + " cannot be both both a virtual_res and a terminal res!" );
		if ( full_to_sub.find( seq_num ) == full_to_sub.end() ) continue;

		core::pose::rna::apply_virtual_rna_residue_variant_type( pose, full_to_sub[ seq_num] );
		//   pose::add_variant_type_to_pose_residue( pose, chemical::VIRTUAL_RNA_RESIDUE, full_to_sub[ seq_num] );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
//Adding Virtual res as root
void
StepWiseRNA_PoseSetup::add_aa_virt_rsd_as_root( core::pose::Pose & pose ) {  //Fang's electron density code

	Size const nres = pose.size();
	Size const working_moving_res( working_parameters_->working_moving_res() );
	//if already rooted on virtual residue, return
	if ( pose.residue( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) {
		TR.Warning << "add_aa_virt_rsd_as_root() called but pose is already rooted on a VRT residue ... continuing." << std::endl;
		return;
	}

	core::chemical::ResidueTypeCOP rsd_type( core::pose::virtual_type_for_pose(pose) );
	core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *rsd_type ) );
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
// how about just setting full_model_info here?
void
StepWiseRNA_PoseSetup::setup_full_model_info( pose::Pose & pose ) const {

	using namespace core::pose;
	using namespace core::pose::full_model_info;

	utility::vector1< core::Size > const & is_working_res( working_parameters_->is_working_res() );
	utility::vector1< core::Size > const & cutpoint_open_list( working_parameters_->cutpoint_open_list() );
	utility::vector1< core::Size > const & fixed_res( working_parameters_->fixed_res() );
	std::string full_sequence = working_parameters_->full_sequence();
	Size const clean_length = core::pose::rna::remove_bracketed( full_sequence ).size();
	utility::vector1< Size > working_res;
	for ( Size i = 1; i <= clean_length; ++i ) {
		if ( is_working_res[ i ] ) {
			working_res.push_back( i );
		}
	}
	if ( working_parameters_->add_virt_res_as_root() ) {
		working_res.push_back( clean_length + 1  );
		full_sequence += 'X';
	}

	FullModelInfoOP full_model_info( new FullModelInfo( full_sequence,  cutpoint_open_list, working_res ) );
	FullModelParametersOP full_model_parameters = full_model_info->full_model_parameters()->clone();
	utility::vector1< core::Size > extra_minimize_res;
	utility::vector1< Size > fixed_domain; // convert from old convention to new one stored in full_model info.
	for ( Size n = 1; n <= is_working_res.size(); n++ ) {
		if ( is_working_res[n] == MOVING_RES ) {
			fixed_domain.push_back( 0 ); // moving and sample-able.
		} else if ( fixed_res.has_value( n ) ) {
			fixed_domain.push_back( 1 ); // freeze this.
		} else {
			fixed_domain.push_back( is_working_res[n] ); // 1 or 2, for input domains 1 or 2.
			if ( is_working_res[n] > 0 ) extra_minimize_res.push_back( n );
		}
	}
	if ( working_parameters_->add_virt_res_as_root() ) fixed_domain.push_back( 1 );
	full_model_parameters->set_parameter( FIXED_DOMAIN, fixed_domain );
	full_model_parameters->set_parameter( INPUT_DOMAIN, fixed_domain );

	// extra_minimize_res was never in use in RNA modeling -- everything was set via -fixed_res or -minimize_res flags, which
	// get loaded in StepWiseRNA_PoseSetupFromCommandLine.
	full_model_parameters->set_parameter_as_res_list( EXTRA_MINIMIZE, extra_minimize_res );
	full_model_parameters->set_parameter_as_res_list( CALC_RMS,     working_parameters_->calc_rms_res() );
	full_model_parameters->set_parameter_as_res_list( RNA_TERMINAL, working_parameters_->terminal_res() );
	full_model_parameters->set_parameter_as_res_list( RNA_SYN_CHI,  working_parameters_->force_syn_chi_res_list() );
	full_model_parameters->set_parameter_as_res_list( RNA_ANTI_CHI,  working_parameters_->force_anti_chi_res_list() );
	full_model_parameters->set_parameter_as_res_list( RNA_NORTH_SUGAR, working_parameters_->force_north_sugar_list() );
	full_model_parameters->set_parameter_as_res_list( RNA_SOUTH_SUGAR, working_parameters_->force_south_sugar_list() );
	full_model_info->set_full_model_parameters( full_model_parameters );
	set_full_model_info( pose, full_model_info );
}

////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseRNA_PoseSetup::setup_vdw_cached_rep_screen_info( pose::Pose & pose ) const {
	//using namespace checker;
	protocols::scoring::fill_vdw_cached_rep_screen_info_from_command_line( pose );
}



} //rna
} //modeler
} //legacy
} //stepwise
} //protocols
