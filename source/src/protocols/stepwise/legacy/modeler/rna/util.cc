// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/legacy/modeler/rna/util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/legacy/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/output_util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/chemical/rna/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>

#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.legacy.modeler.rna.util" );

using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::modeler::rna;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace rna {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	check_can_prepend( utility::vector1< core::Size > const & seq_num_list ){
		for ( Size n = 1; n <= seq_num_list.size() - 1; n++ ){ //[11, 12, 13]
			if ( ( seq_num_list[n] + 1 ) != seq_num_list[n + 1] ) return false;
		}
		return true;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	check_can_append( utility::vector1< core::Size > const & seq_num_list ){
		for ( Size n = 1; n <= seq_num_list.size() - 1; n++ ){ //[14, 13, 12]
			if ( ( seq_num_list[n] - 1 ) != seq_num_list[n + 1] ) return false;
		}
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	//Sort by the first element. Low number on the top of the list
	bool
	pair_sort_criterion( std::pair < Size, Size > pair_one, std::pair < Size, Size > pair_two ){
		return ( pair_one.first < pair_two.first );
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	void
	sort_pair_list( utility::vector1< std::pair < Size, Size > > pair_list ){
		sort( pair_list.begin(), pair_list.end(), pair_sort_criterion );
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	void
	output_pair_size( std::pair < Size, Size > const & pair_size, std::ostream & outstream /* = std::cout */ ){
		outstream << "( " << pair_size.first << ", " << pair_size.second << " ) ";
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	void
	output_pair_size( utility::vector1 < std::pair < Size, Size > > const & pair_size_vector, std::string const & output_string, std::ostream & outstream /* = std::cout */, core::Size const spacing ){
		outstream << std::setw( spacing ) << std::left << output_string << " :";
		for ( Size n = 1; n <= pair_size_vector.size(); n++ ){
			output_pair_size( pair_size_vector[n], outstream );
		}
		outstream << std::endl;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	void output_is_prepend_map( std::string const tag, std::map< core::Size, bool > const & my_map, core::Size const max_seq_num, std::ostream & outstream /* = std::cout */, core::Size const tag_spacing ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;

		outstream << std::setw( tag_spacing ) << tag;

		Size spacing = 4;
//		outstream << std::setw(30) << "is_residue_prepend:";
		for ( Size seq_num = 1; seq_num <= max_seq_num; seq_num++ ){
			char prepend_char;
			if ( my_map.find( seq_num ) != my_map.end() ){
				prepend_char = ( my_map.find( seq_num )->second ) ? 'P' : 'A';
			} else{
				prepend_char = '-';
			}
			outstream << std::setw( spacing ) << prepend_char;
		}
		outstream << std::endl;

	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	void
	output_bool_list( std::string const tag, utility::vector1< Size > const & size_list, std::ostream & outstream /* = std::cout */, core::Size const spacing ){
		utility::vector1< bool > bool_list;

		for ( Size n = 1; n <= size_list.size(); n++ ){
			bool_list.push_back( size_list[n] );
		}
		output_bool_list( tag, bool_list, outstream, spacing );
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	void
	output_bool_list( std::string const tag, utility::vector1< bool > const & bool_list, std::ostream & outstream /* = std::cout */, core::Size const spacing ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;

		outstream <<  std::setw( spacing ) << tag;

		for ( Size seq_num = 1; seq_num <= bool_list.size(); seq_num++ ){
			output_boolean( bool_list[seq_num], outstream );
		}
		outstream << std::endl;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	void
	output_size_list( std::string const tag, utility::vector1< Size > const & size_list, std::ostream & outstream /* = std::cout */, core::Size const spacing ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;

		outstream <<  std::setw( spacing ) << tag;

		for ( Size seq_num = 1; seq_num <= size_list.size(); seq_num++ ){
			outstream << I( 4, size_list[seq_num] );
		}
		outstream << std::endl;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	void
	output_fold_tree_info( kinematics::FoldTree const & fold_tree, std::string const pose_name, std::ostream & outstream /* = std::cout */ ){

		outstream << "fold tree of " << pose_name << ": " << std::endl;
		for ( int i = 1; i <= fold_tree.num_cutpoint(); i++ ){
			outstream << std::setw( 30 ) << "jump_point_num = " << i;
			outstream << "   cutpoint = " << fold_tree.cutpoint( i );
			outstream << "   5' jump_point = " << fold_tree.jump_point( 1, i ) << ", " << fold_tree.upstream_atom( i );
			outstream << "   3' jump_point = " << fold_tree.jump_point( 2, i ) << ", " << fold_tree.downstream_atom( i ) << std::endl;
		}
	}

	void
	output_fold_tree_info( pose::Pose const & pose, std::string pose_name, std::ostream & outstream /* = std::cout */ ){
		output_fold_tree_info( pose.fold_tree(), pose_name, outstream );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Ok this function should only be called if pose contain full sequence.
	//This one include edge phosphates
	core::Real
	full_length_rmsd_over_residue_list( pose::Pose const & pose1, pose::Pose const & pose2, utility::vector1 < Size > const & residue_list, std::string const & full_sequence, bool const verbose, bool const ignore_virtual_atom ){

		using namespace ObjexxFCL;


		if ( pose1.sequence() != full_sequence ){
			TR << "pose1.sequence() = " << pose1.sequence() << std::endl;
			TR << "pose2.sequence() = " << pose2.sequence() << std::endl;
			TR << "full_sequence = " << full_sequence << std::endl;
			utility_exit_with_message( "pose1.sequence() != full_sequence" );
		}

		if ( pose2.sequence() != full_sequence ){
			TR << "pose1.sequence() = " << pose1.sequence() << std::endl;
			TR << "pose2.sequence() = " << pose2.sequence() << std::endl;
			TR << "full_sequence = " << full_sequence << std::endl;
			utility_exit_with_message( "pose2.sequence() != full_sequence" );
		}

		Size const total_res = pose1.total_residue();

		if ( verbose ){
			output_title_text( "Enter full_length_rmsd_over_residue_list function", TR );
			output_boolean( "ignore_virtual_atom = ", ignore_virtual_atom, TR ); TR << std::endl;
			output_seq_num_list( "residue_list = ", residue_list, TR, 30 );
		}

		Size atom_count = 0;
		Real sum_sd = 0;

		for ( Size i = 1; i <= residue_list.size(); i++ ){

			Size const full_seq_num = residue_list[i];

			bool is_prepend = false;
			bool both_pose_res_is_virtual = false;

			if ( pose1.residue( full_seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) &&
					pose2.residue( full_seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ){
				both_pose_res_is_virtual = true;
			}

			if ( ( full_seq_num + 1 ) <= total_res ){
				if ( pose1.residue( full_seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ){
					runtime_assert ( pose1.residue( full_seq_num + 1 ).has_variant_type(
							core::chemical::VIRTUAL_RNA_RESIDUE_UPPER ) );
				}

				if ( pose2.residue( full_seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ){
					runtime_assert ( pose2.residue( full_seq_num + 1 ).has_variant_type(
							core::chemical::VIRTUAL_RNA_RESIDUE_UPPER ) );
				}
			}

			if ( verbose ){
				TR << "full_seq_num = " << full_seq_num;
				output_boolean( " is_prepend = ", is_prepend, TR );
				output_boolean( " both_pose_res_is_virtual = ", both_pose_res_is_virtual, TR ); TR << std::endl;
			}

			if ( both_pose_res_is_virtual ) continue;

			//add atom in the suites to atom_count
			//add sd of each atom to sum_sd
			suite_square_deviation( pose1, pose2, is_prepend, full_seq_num, full_seq_num, atom_count, sum_sd, verbose, ignore_virtual_atom );

			if ( ( ( full_seq_num + 1 ) <= total_res ) && residue_list.has_value( full_seq_num + 1 ) == false ){

				if ( verbose ) TR << "Phosphate_edge_res_( full_seq_num + 1 ) = " << full_seq_num + 1 << std::endl;

				phosphate_square_deviation( pose1, pose2, full_seq_num + 1, full_seq_num + 1, atom_count, sum_sd, verbose, ignore_virtual_atom );
			}

		}


		sum_sd = sum_sd/( atom_count );
		Real rmsd = sqrt( sum_sd );

		if ( atom_count == 0 ) rmsd = 0.0; //special case...implement this on May 5, 2010

		if ( verbose ){
			TR << "sum_sd = " << sum_sd << " atom_count = " << atom_count << " rmsd = " << rmsd << std::endl;
			output_title_text( "Exit In full_length_rmsd_over_residue_list function", TR );
		}

		return ( std::max( 0.01, rmsd ) );

	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	print_backbone_torsions( pose::Pose const & pose, Size const five_prime_chainbreak ){

		using namespace core::id;

		conformation::Residue const & suite_lower_res = pose.residue( five_prime_chainbreak );
		TR << std::setw( 5 ) << " ep = " << std::setw( 15 ) << suite_lower_res.mainchain_torsion( 5 );
		TR << std::setw( 5 ) << " z = "  << std::setw( 15 ) << suite_lower_res.mainchain_torsion( 6 );


		Size const three_prime_chainbreak = five_prime_chainbreak + 1;

		if ( three_prime_chainbreak <= pose.total_residue() ) {
			conformation::Residue const & suite_upper_res = pose.residue( three_prime_chainbreak );
			TR << std::setw( 5 ) << " a = "  << std::setw( 15 ) << suite_upper_res.mainchain_torsion( 1 );
			TR << std::setw( 5 ) << " b = "  << std::setw( 15 ) << suite_upper_res.mainchain_torsion( 2 );
			TR << std::setw( 5 ) << " g = "  << std::setw( 15 ) << suite_upper_res.mainchain_torsion( 3 );
		}

		TR << std::endl;

	}

	//////////////////////////////////////////////////////////////////////////
 	core::Size
	setup_chain_break_jump_point( core::pose::Pose & pose,
																core::Size const moving_res,
																core::Size const reference_res ){
		return setup_bulge_jump_point( pose, moving_res, reference_res, false /*verbose*/ );
	}

	//////////////////////////////////////////////////////////////////////////
	void
	remove_chain_break_jump_point( core::pose::Pose & pose,
																 core::Size const moving_res,
																 core::Size const reference_res ){
		kinematics::FoldTree f = pose.fold_tree();
		Size const n = f.jump_nr( moving_res, reference_res );
		runtime_assert( n > 0 );
		f.delete_jump_and_intervening_cutpoint( n );
		pose.fold_tree( f );
	}

	//////////////////////////////////////////////////////////////////////////
	core::Size
	setup_bulge_jump_point( pose::Pose & pose,
													Size const & moving_base,
													Size const & reference_base,
													bool const verbose ){

		using namespace core::conformation;

		runtime_assert ( moving_base != reference_base );
		int i, j;
		Size cutpoint;
		if ( moving_base > reference_base ){
			i = reference_base;
			j = moving_base;
			cutpoint = moving_base - 1;
		} else{
			i = moving_base;
			j = reference_base;
			cutpoint = moving_base;
		}

		core::kinematics::FoldTree fold_tree = pose.fold_tree();	//HARD COPY?

		if( verbose ) output_fold_tree_info( fold_tree, "Before add bulge jump point", TR.Debug );
		fold_tree.new_jump( reference_base, moving_base, cutpoint ); //Choose the residue five_prime of the actual cutpoint position

		if ( verbose ) TR.Debug << "after add new jump point" << std::endl;

		Size jump_num = fold_tree.jump_nr( i, j );
		runtime_assert( jump_num > 0 );

		Residue const & rsd1( pose.residue( i ) );
		Residue const & rsd2( pose.residue( j ) );
		fold_tree.set_jump_atoms( jump_num, rsd1.atom_name( rsd1.chi_atoms( 1 )[4] ), rsd2.atom_name( rsd2.chi_atoms( 1 )[4] ) ); //Base atoms...
		if ( verbose ) output_fold_tree_info( fold_tree, "New fold_tree with bulge jump point", TR.Debug );
		pose.fold_tree( fold_tree );

		return cutpoint;
	}

	//////////////////////////////////////////////////////////////////////////////////////
	core::Size
	get_residue_base_state( core::pose::Pose const & pose, Size const seq_num ){

		using namespace core::scoring;
		using namespace core::chemical::rna;

		Real const CHI_CUTOFF = 15.0; //Kinda RANDOM..ROUGH average between north chi_anti (~79) and north chi_syn (~-50)

		conformation::Residue const & rsd = pose.residue( seq_num );
		Real const chi = numeric::principal_angle_degrees( rsd.chi( CHI - NUM_RNA_MAINCHAIN_TORSIONS ) );

		if ( chi <= CHI_CUTOFF ){
			return SYN;
		} else{
			return ANTI;
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////

	core::Size
	get_residue_pucker_state( core::pose::Pose const & pose, Size const seq_num, bool const verbose ){

		using namespace core::scoring;
		using namespace core::chemical::rna;

		static RNA_FittedTorsionInfo const rna_fitted_torsion_info;
		Real const DELTA_CUTOFF( rna_fitted_torsion_info.delta_cutoff() );

		if ( verbose ) TR.Debug << "  DELTA_CUTOFF angle = " << DELTA_CUTOFF;

		conformation::Residue const & rsd( pose.residue( seq_num ) );
		Real delta = numeric::principal_angle_degrees( rsd.mainchain_torsion( DELTA ) );

		if ( ( delta > 1.0 && delta < 179.00 ) == false ){
			TR.Debug << " seq_num = " << seq_num << " delta angle = " << delta << std::endl;

/////////////////////////
			if ( pose.residue( seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ){
				TR.Debug << "Warning: delta angle is out of range for virtual_residue at seq_num " << seq_num << "!" << std::endl;
			} else{
				//This part is now obsolete .... Apr 30, 2010...
				//A possibility for a out of range delta is for imported pdbs (upper element of 1q93 for example).
				Real principal_delta = numeric::principal_angle_degrees( delta );

				//Consistency check
				if ( delta <  - 180 ){
					if ( ( delta + 359 < principal_delta ) == false || ( delta + 361 > principal_delta ) == false ){
						utility_exit_with_message( "delta <  - 180 but delta + 359 < principal_delta ) == false || ( delta + 361 > principal_delta ) == false!" );
					}
				}

				if ( delta > 180 ){
					if ( ( delta - 359 > principal_delta ) == false || ( delta - 361 < principal_delta ) == false ){
						utility_exit_with_message( "delta > 180 but delta - 359 > principal_delta ) == false || ( delta - 361 < principal_delta ) == false!" );
					}
				}

				delta = principal_delta;

				//Check again
				if ( ( delta > 1.0 && delta < 179.00 ) == false ) utility_exit_with_message( "principal delta angle out of range!" );
//////////////////////////
			}
		}

		if ( verbose ) TR.Debug << "  delta angle = " << delta << std::endl;

		if ( delta <= DELTA_CUTOFF ) {
			return NORTH;
		} else {
			return SOUTH;
		}
	}


	////////////////////////////////////////////////////////////////////////
	void
	apply_rotamer( pose::Pose & pose, utility::vector1< Torsion_Info > const & rotamer_list ){
		for ( Size i = 1; i <= rotamer_list.size(); i++ ) {
			pose.set_torsion( rotamer_list[ i ].id, rotamer_list[i].value );
		}
	}


	bool
	is_same_sugar_pucker( core::pose::Pose const & current_pose, core::pose::Pose const & cluster_center_pose, Size const seq_num ){
		if ( get_residue_pucker_state( current_pose, seq_num ) == get_residue_pucker_state( cluster_center_pose, seq_num ) ){
			return true;
		} else{
			return false;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	setup_simple_fold_tree( core::pose::Pose & pose ){

//		using namespace core::chemical;

		Size const nres = pose.total_residue();

		kinematics::FoldTree simple_fold_tree( nres ); //create a simple fold tree

		simple_fold_tree.simple_tree( nres ); //Just to make sure.

		pose.fold_tree( simple_fold_tree );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void
	import_pose_from_silent_file( core::pose::Pose & import_pose, std::string const & silent_file , std::string const & input_tag ){

		using namespace core::chemical;
		using namespace core::conformation;

		static const ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance() ->
			residue_type_set( core::chemical::FA_RNA );

		core::io::silent::SilentFileData silent_file_data;
		silent_file_data.read_file( silent_file );

		Size num_matching_tag = 0;

		for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(), end = silent_file_data.end(); iter != end; ++iter ){
			if ( iter->decoy_tag() != input_tag ) continue;
			num_matching_tag += 1;
			iter->fill_pose( import_pose, *rsd_set );
		}

		if ( num_matching_tag != 1 ){
			utility_exit_with_message( "num_matching_tag = ( " + ObjexxFCL::string_of( num_matching_tag ) + " ) != 1 for tag " + input_tag + " in silent file ( " + silent_file + " )!" );
		}

		if ( check_for_messed_up_structure( import_pose, input_tag ) == true ){
		 	utility_exit_with_message( "import_pose " + input_tag + " from silent_file " + silent_file + " is a messed up pose!" );
		}

	}

	/////////////////New function on Nov 11, 2010///////////////
	std::string
	get_tag_from_pdb_filename( std::string const pdb_filename ){

		std::string tag;

		size_t found = pdb_filename.rfind( '/' );

		if ( found != std::string::npos ){
			tag = pdb_filename.substr( found + 1 );
		} else {
			tag = pdb_filename;
		}

		size_t found_2 = tag.rfind( ".pdb" );

		if ( found_2 != std::string::npos ){
			tag = tag.substr( 0, tag.size() - 4 );
		}

		return tag;
	}

	//DUPLICATE OF CODE IN StepWiseWorkingParametersSetup.cc
	void
	move_jump_atom_to_base( core::kinematics::FoldTree & fold_tree, std::string const & working_sequence ){

		Size const num_cutpoint = fold_tree.num_cutpoint();

		for ( Size i = 1; i <= num_cutpoint; i++ ) {
			Size const k = fold_tree.upstream_jump_residue( i );
			Size const m = fold_tree.downstream_jump_residue( i );

			char upstream_res = working_sequence[k - 1];
			char downstream_res = working_sequence[m - 1];

			//Base atoms...
			//chi_atoms(1)[4] )=  C2 if URA or RCY
			//chi_atoms(1)[4] )=  C4 if RGU or RAD
			std::string upstream_jump_atom;
			std::string downstream_jump_atom;

			if ( upstream_res == 'u' || upstream_res == 'c' ){
				upstream_jump_atom = " C2 ";
			} else if ( upstream_res == 'a' || upstream_res == 'g' ){
				upstream_jump_atom = " C4 ";
			} else{
				utility_exit_with_message( "Invalid upstream_res!!" );
			}

			if ( downstream_res == 'u' || downstream_res == 'c' ){
				downstream_jump_atom = " C2 ";
			} else if ( downstream_res == 'a' || downstream_res == 'g' ){
				downstream_jump_atom = " C4 ";
			} else{
				utility_exit_with_message( "Invalid downstream_res!!" );
			}

			TR << "upstream_res = " << k << upstream_res << " upstream_jump_atom = " << upstream_jump_atom;
			TR << " downstream_res = " << k << downstream_res << " downstream_jump_atom = " << downstream_jump_atom << std::endl;

			fold_tree.set_jump_atoms( i, downstream_jump_atom, upstream_jump_atom );

		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	print_WorkingParameters_info( working_parameters::StepWiseWorkingParametersCOP const & const_WP, std::string const WP_name, std::ostream & outstream /* = std::cout */, bool const is_simple_full_length_WP  ){

		working_parameters::StepWiseWorkingParametersOP WP = new working_parameters::StepWiseWorkingParameters;

		( *WP ) = ( *const_WP );

		print_WorkingParameters_info( WP, WP_name, outstream, is_simple_full_length_WP );

	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	print_WorkingParameters_info( working_parameters::StepWiseWorkingParametersOP const & WP, std::string const WP_name, std::ostream & outstream /* = std::cout */, bool const is_simple_full_length_WP ){

		using namespace ObjexxFCL;

		output_title_text( "printing WorkingParameters Information for " + WP_name, outstream );

		utility::vector1< Size > empty_seq_num_list;
		empty_seq_num_list.clear();

		outstream << "full_sequence = " <<  WP->full_sequence();
		outstream << " working_sequence = " <<  WP->working_sequence();
		outstream << " moving_res = " <<  WP->moving_res();
		outstream << " working_moving_res = " << WP->working_moving_res();
		outstream << " working_moving_suite = " <<  WP->working_moving_suite() << std::endl;


		outstream << "gap_size = " << WP->gap_size();
		outstream << " five_prime_chain_break_res = " << WP->five_prime_chain_break_res();
		output_boolean( " is_prepend = ",  WP->is_prepend(), outstream ) ;
		output_boolean( " is_internal = ", WP->is_internal(), outstream );
		output_boolean( " output_extra_RMSDs = ", WP->	output_extra_RMSDs(), outstream ); outstream << std::endl;
		output_boolean( " floating_base = ", WP->floating_base(), outstream ); outstream << std::endl;


		//std::map< core::Size, core::Size > full_to_sub_;
		//std::map< core::Size, core::Size > sub_to_full_;

		//utility::vector1< std::pair< core::Size, core::Size > > chain_boundaries_;

		//	ObjexxFCL::FArray1D< bool > partition_definition_;

		//core::pose::PoseOP working_native_pose_;

		outstream << "------------full_stuff------------" << std::endl;
		output_bool_list( "is_working_res = ", WP->is_working_res(), outstream );

		for ( Size n = 1; n <= WP->input_res_vectors().size(); n++ ){
			output_seq_num_list( "input_res_vectors[" + string_of( n ) + "]", WP->input_res_vectors()[n], outstream );
		}
		output_seq_num_list( "global_sample_res_list = ", WP->global_sample_res_list(), outstream );

		if ( WP->force_syn_chi_res_list().size() > 0 )  output_seq_num_list( "force_syn_chi_res_list = ", WP->force_syn_chi_res_list(), outstream );
		if ( WP->force_north_sugar_list().size() > 0 ) output_seq_num_list( "force_north_sugar_list = ", WP->force_north_sugar_list(), outstream );
		if ( WP->force_south_sugar_list().size() > 0 ) output_seq_num_list( "force_south_sugar_list = ", WP->force_south_sugar_list(), outstream );
		if ( WP->protonated_H1_adenosine_list().size() > 0 ) output_seq_num_list( "protonated_H1_adenosine_list ", WP->protonated_H1_adenosine_list(), outstream );





		output_is_prepend_map( "is_prepend_map = ", WP->is_prepend_map(), WP->full_sequence().size(), outstream );
		output_seq_num_list( "calc_rms_res = ", 								WP->calc_rms_res() 								, outstream );

		output_seq_num_list( "native_alignment = ",  							WP->native_alignment(), outstream );
		output_seq_num_list( "cutpoint_closed_list = ", 					WP->cutpoint_closed_list(), outstream );


		outstream << "------------working_stuff------------" << std::endl;

		output_seq_num_list( "working_global_sample_res_list = ", WP->working_global_sample_res_list(), outstream );
		if ( WP->force_syn_chi_res_list().size() > 0 )  output_seq_num_list( "working_force_syn_chi_res_list = ", WP->working_force_syn_chi_res_list(), outstream );
		if ( WP->force_north_sugar_list().size() > 0 ) output_seq_num_list( "working_force_north_sugar_list = ", WP->working_force_north_sugar_list(), outstream );
		if ( WP->force_south_sugar_list().size() > 0 ) output_seq_num_list( "working_force_south_sugar_list = ", WP->working_force_south_sugar_list(), outstream );
		if ( WP->protonated_H1_adenosine_list().size() > 0 ) output_seq_num_list( "working_protonated_H1_adenosine_list = ", WP->working_protonated_H1_adenosine_list(), outstream );


		output_seq_num_list( "working_fixed_res = ",						WP->working_fixed_res() 						, outstream );

		if ( is_simple_full_length_WP == false ){
			output_seq_num_list( "working_moving_res_list = ", WP->working_moving_res_list(), outstream );
			output_seq_num_list( "working_moving_suite_list = ", WP->working_moving_suite_list(), outstream );
		} else{
			output_seq_num_list( "line_filler ", empty_seq_num_list, outstream );
			output_seq_num_list( "line_filler ", empty_seq_num_list, outstream );
		}

		output_seq_num_list( "working_terminal_res = ", 					WP->working_terminal_res() 					, outstream );
		output_seq_num_list( "working_moving_partition_res = ",  WP->working_moving_partition_res(), outstream );

		output_seq_num_list( "working_best_alignment = ", 				WP->working_best_alignment(), outstream );
		output_seq_num_list( "working_native_alignment = ", 			WP->working_native_alignment(), outstream );

		if ( is_simple_full_length_WP == false ){
			utility::vector1< bool > vector1_partition_definition;

			for ( Size n = 1; n <= WP->partition_definition().size(); n++ ){
				vector1_partition_definition.push_back( WP->partition_definition()( n ) );
			}

			output_bool_list( "partition_definition = ", 	vector1_partition_definition, outstream );

			Size const root_res = (WP->fold_tree().size() > 0) ? WP->fold_tree().root() : 0;
			outstream << "root_res = " << root_res << std::endl;
			outstream << "working_reference_res = " << WP->working_reference_res() << std::endl;
			output_fold_tree_info( WP->fold_tree(), "fold_tree", outstream );
		}

		output_title_text( "", outstream );


	}



/////////////////////////////////////////////////////////////////////////////////////////////
	void
	set_nucleotide_to_A_form( pose::Pose & pose, Size const seq_num ){
		//Torsion value extracted from 3DNA (website) (A-U BP repeating) idealized A-form helix. Note that bond angle and bond length of idealized Rosetta doesn't exactly match the values in 3DNA

		using namespace core::id;

		pose.set_torsion( TorsionID( seq_num, id::BB,  1 ), -68.9 ); //alpha
		pose.set_torsion( TorsionID( seq_num, id::BB,  2 ), 179.5 ); //beta
		pose.set_torsion( TorsionID( seq_num, id::BB,  3 ), 54.5 ); //gamma
		pose.set_torsion( TorsionID( seq_num, id::BB,  5 ), -154.0 ); //epsilon
		pose.set_torsion( TorsionID( seq_num, id::BB,  6 ), -70.8 ); //zeta

		pose.set_torsion( TorsionID( seq_num, id::BB,  4 ), 82.2 ); //delta
		pose.set_torsion( TorsionID( seq_num, id::CHI, 1 ), 79.2 ); //chi
		pose.set_torsion( TorsionID( seq_num, id::CHI, 2 ), 36.9 ); //nu2
		pose.set_torsion( TorsionID( seq_num, id::CHI, 3 ), 94.7 ); //nu1


	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::string
	path_basename( std::string const full_path ){

		size_t found = full_path.rfind( '/' );

		std::string basename;

		if ( found != std::string::npos ){
			basename = full_path.substr( found + 1 );
		} else {
			basename = full_path;
		}

		return basename;
 	}


} //rna
} //modeler
} //legacy
} //stepwise
} //protocols
