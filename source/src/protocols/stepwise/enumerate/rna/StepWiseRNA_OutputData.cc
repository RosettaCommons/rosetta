// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_OutputData.cc (Created on Sept 26, 2011)
/// @brief Output silent_file_data functions for Stepwise Assembly RNA.
/// @detailed
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_OutputData.hh> //Oct 22, 2011...Not sure why the code worked without this!
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Util.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_FloatingBaseSamplerUtil.hh> //Sept 26, 2011
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_ResidueInfo.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/rotamer_sampler/rigid_body/EulerAngles.hh>
#include <protocols/farna/RNA_BasePairClassifier.hh>
#include <core/scoring/rna/RNA_BaseDoubletClasses.hh>

#include <core/scoring/ScoreType.hh> //Parin Sept 20, 2011.
//////////////////////////////////

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <numeric/conversions.hh>
#include <numeric/NumericTraits.hh>

#include <basic/Tracer.hh>

#include <iostream>
#include <fstream>
#include <sstream>
#include <ObjexxFCL/format.hh>
#include <set>
#include <time.h>
#include <map>

#include <stdio.h> //Sept 26, 2011

//for process_mem_usage:
#include <ios>


using namespace core;

static basic::Tracer TR( "protocols.stepwise.rna.StepWiseRNA_OutputData" );

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	core::io::silent::BinaryRNASilentStruct
	get_binary_rna_silent_struct_safe( pose::Pose const & const_pose, std::string const & tag, std::string const & silent_file ){

		// What's the deal with this creation of silent struct and deletion? -- rhiju
		//    do you remember the scenario in which you really needed all this?
		//
		//        What I encountered was that sometime a pose structure that have entirely
		//        normal coordinates inside Rosetta becomes messed up when written as a
		//        silent_struct. I also noticed that if the pose is then slightly rotated
		//        (euler angles), then this problem will disappear.
		//
		//        So essentially, the get_binary_rna_silent_struct_safe() function writes a
		//        pose as a silent_struct to file and then checks whether the silent_struct is
		//        messed up. If it is messed up, the function rotates the pose and rewrite the
		//        silent_struct. This process continues until the silent_struct is no longer
		//        messed-up or maximum trial (10) is reached.
		//
		//        For most cases, the silent_struct is successfully written on the first
		//        trial since a very small fraction (perhaps 0.01%) of pose experiences
		//        this problem. -- parin (2013)

		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::chemical;
		using namespace core::conformation;


		std::string const debug_silent_file = silent_file + "_CONVERSION_DEBUG";
		std::string const debug_tag = tag + "_CONVERSION_DEBUG";
                const Real RADS_PER_DEG = numeric::NumericTraits < Real > ::pi() / 180.;

		SilentFileData silent_file_data;

		static const ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance() ->
			residue_type_set( core::chemical::RNA );

		Size NUM_trials = 10;
		Real const local_angle_bin_size = 20;
		Real const local_z_bin_size = 0.05;

		pose::Pose first_trial_pose_from_silent_file;
		pose::Pose first_trial_pose;

		for ( Size trial_num = 1; trial_num <= NUM_trials; trial_num++ ){ //Found that just rigid problem rotation of the pose solves the silent_file conversion problem

			pose::Pose pose = const_pose;

			if ( trial_num != 1 ){

				///////////////////////////get centroid of the structure/////////////////////////////

				numeric::xyzVector< core::Real > centroid = Vector( 0.0, 0.0, 0.0 );
				Size numatoms = 0;

				for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

					conformation::Residue const & rsd( pose.residue( seq_num ) );

					for ( Size at = 1; at <= rsd.natoms(); at++ ){

						if ( rsd.is_virtual( at ) ) continue;

			    		centroid += rsd.xyz( at );
 					  	numatoms++;
					}
				}

				if ( numatoms == 0 ) utility_exit_with_message( "numatoms == 0" );

				centroid = centroid/numatoms;

				////////////////////////////////////////////////////////////////////////////////
				rotamer_sampler::rigid_body::EulerAngles euler_angles;
				Matrix rotation_matrix;
				euler_angles.set_alpha( ( 0.25*trial_num )*local_angle_bin_size*( RADS_PER_DEG ) );
				euler_angles.set_gamma( ( 0.25*trial_num )*local_angle_bin_size*( RADS_PER_DEG ) );
				euler_angles.set_z( ( 0.25*trial_num )*local_z_bin_size );	//MAKE SURE THIS DOESN'T GET OUT OF BOUND!
				euler_angles.convert_to_rotation_matrix( rotation_matrix );

				for ( Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++ ){

					conformation::Residue const & rsd( pose.residue( seq_num ) );

					for ( Size at = 1; at <= rsd.natoms(); at++ ){

						id::AtomID const id( at, seq_num );

						pose.set_xyz( id, pose.xyz( id ) - centroid ); //This should minimize the error introduced by the rigid body rotation!
						pose.set_xyz( id, rotation_matrix * pose.xyz( id ) );
					}
				}
			}

			BinaryRNASilentStruct DEBUG_silent_struct( pose, debug_tag );
			BinaryRNASilentStruct const silent_struct( pose, tag );

			if ( file_exists( debug_silent_file ) ) remove_file( debug_silent_file );

			silent_file_data.write_silent_struct( DEBUG_silent_struct, debug_silent_file, false );

			///////////////////////////////////////////////////////////////////////////////

			core::io::silent::SilentFileData import_silent_file_data;
			import_silent_file_data.read_file( debug_silent_file );
			pose::Pose pose_from_silent_file;

			bool found_tag = false;
			Size num_struct = 0;

			for ( core::io::silent::SilentFileData::iterator iter = import_silent_file_data.begin(), end = import_silent_file_data.end(); iter != end; ++iter ){
				num_struct += 1;
				if ( iter->decoy_tag() != debug_tag ) continue;
				found_tag = true;
				iter->fill_pose( pose_from_silent_file, *rsd_set );
			}

			if ( num_struct != 1 ) utility_exit_with_message( "num_struct = ( " + string_of( num_struct ) + " ) != 1" );
			if ( found_tag == false ) utility_exit_with_message( "Could not find specified tag ( " + debug_tag + " ) in silent file ( " + debug_silent_file + " )!" );

			if ( file_exists( debug_silent_file ) == false ){
				utility_exit_with_message( "debug_silent_file ( " + debug_silent_file + " ) SHOULD exist!" );
			}

			remove_file( debug_silent_file );

			if ( trial_num == 1 ){
				first_trial_pose_from_silent_file = pose_from_silent_file;
				first_trial_pose = pose;
			}

			if ( check_for_messed_up_structure( pose_from_silent_file, debug_tag ) == false ){
				return silent_struct;

			} else{
			 	TR << "WARNING: Problem with writing pose ( " << debug_tag << " ) to silent_file [Attempt #" << trial_num << "]" << std::endl;
			}

		}

		first_trial_pose_from_silent_file.dump_pdb( "SILENT_FILE_CONVERSION_PROBLEM_" + tag + "_pose_from_silent_file.pdb" );
		first_trial_pose.dump_pdb( "SILENT_FILE_CONVERSION_PROBLEM_" + tag + ".pdb" );
		BinaryRNASilentStruct ERROR_silent_struct( first_trial_pose, debug_tag );
		std::string const ERROR_silent_file = "SILENT_FILE_CONVERSION_PROBLEM_" + tag + ".out";
		silent_file_data.write_silent_struct( ERROR_silent_struct, ERROR_silent_file, false );

		utility_exit_with_message( "Fail to write pose ( " + debug_tag + " ) to silent_file after " + string_of( NUM_trials ) + " trials " );

		////////////This is just to prevent compiler WARNING MESSAGES/////////
		BinaryRNASilentStruct EMPTY_silent_struct;
		return EMPTY_silent_struct;
		//////////////////////////////////////////////////////////////////////

	}


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	core::io::silent::BinaryRNASilentStruct
	get_binary_rna_silent_struct_safe_wrapper( pose::Pose const & const_pose, std::string const & tag, std::string const & silent_file, bool const write_score_only ){

		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace core::pose;

		if ( write_score_only ){

			BinaryRNASilentStruct s( const_pose, tag ); //If write score only, don't have to safe about pose to silent_struct conversion!
			return s;

		} else{

			return ( get_binary_rna_silent_struct_safe( const_pose, tag, silent_file ) );

		}

		////////////This is just to prevent compiler WARNING MESSAGES/////////
		BinaryRNASilentStruct EMPTY_silent_struct;
		return EMPTY_silent_struct;
		//////////////////////////////////////////////////////////////////////

	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	output_data( std::string const & silent_file, std::string const & tag, bool const write_score_only, pose::Pose const & pose, core::pose::PoseCOP native_poseCOP, StepWiseRNA_JobParametersCOP job_parameters_ ){
		static core::io::silent::SilentFileData silent_file_data;
		output_data( silent_file_data, silent_file, tag, write_score_only, pose, native_poseCOP, job_parameters_ );
	}


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Accept the job_parameter instead.
	void
	output_data( core::io::silent::SilentFileData& silent_file_data, std::string const & silent_file, std::string const & tag, bool const write_score_only, pose::Pose const & pose, core::pose::PoseCOP native_poseCOP, StepWiseRNA_JobParametersCOP job_parameters_ ){

		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace core::pose;

		utility::vector1 < core::Size > const & rmsd_res_list = job_parameters_->rmsd_res_list();
		std::map< core::Size, core::Size > const & full_to_sub = job_parameters_->const_full_to_sub();
		std::map< core::Size, bool > const & is_prepend_map = job_parameters_->is_prepend_map();
		bool const is_prepend(  job_parameters_->is_prepend() ); // if true, moving_suite+1 is fixed. Otherwise, moving_suite is fixed.
		Size const moving_base_residue( job_parameters_->actually_moving_res() );

		BinaryRNASilentStruct s = get_binary_rna_silent_struct_safe_wrapper( pose, tag, silent_file, write_score_only );

		bool const output_extra_RMSDs = job_parameters_->output_extra_RMSDs();

		if ( native_poseCOP ){

			if ( write_score_only ){ //Basically the optimal alignment, align working_res as well if it is part of the alignment res list.

				s.add_energy( "all_rms", rms_at_corresponding_heavy_atoms( pose, *native_poseCOP ) );

				// This assumes that pose and native_pose are correctly syperimposed.
				// I added a function in Pose_Setup to make sure this happens. Parin Jan 28, 2010

				s.add_energy( "rmsd", suite_rmsd( pose, *native_poseCOP, moving_base_residue, is_prepend, false ) );
				s.add_energy( "loop_rmsd", rmsd_over_residue_list( pose, *native_poseCOP, rmsd_res_list, full_to_sub, is_prepend_map, false, false ) );

				s.add_energy( "V_rms", suite_rmsd( pose, *native_poseCOP, moving_base_residue, is_prepend, true ) );
				s.add_energy( "V_loop_rms", rmsd_over_residue_list( pose, *native_poseCOP, rmsd_res_list, full_to_sub, is_prepend_map, false, true ) );

				if ( job_parameters_->gap_size() == 0 ){
					s.add_energy( "PBP_rmsd", phosphate_base_phosphate_rmsd( pose, *native_poseCOP, moving_base_residue,  false ) );
				} else{
					s.add_energy( "PBP_rmsd", 0.0 );
				}

			} else{

				utility::vector1< core::Size > const & working_native_alignment = job_parameters_->working_native_alignment();
				utility::vector1< core::Size > const & working_best_alignment = job_parameters_->working_best_alignment();

				if ( output_extra_RMSDs ){

					s.add_energy( "all_rms", rms_at_corresponding_heavy_atoms( pose, *native_poseCOP ) );

					pose::Pose current_pose = pose; //hard copy, computationally expensive

					if ( working_native_alignment.size() != 0 ){ //user specify which residue to align with native.
						align_poses( current_pose, tag, *native_poseCOP, "native", working_native_alignment );
					} else{ //default
						align_poses( current_pose, tag, *native_poseCOP, "native", working_best_alignment );
					}
					s.add_energy( "O_rmsd", suite_rmsd( current_pose, *native_poseCOP, moving_base_residue, is_prepend, false ) );
					s.add_energy( "O_loop_rmsd", rmsd_over_residue_list( current_pose, *native_poseCOP, rmsd_res_list, full_to_sub, is_prepend_map, false, false ) );

					s.add_energy( "O_V_rms", suite_rmsd( current_pose, *native_poseCOP, moving_base_residue, is_prepend, true ) );
					s.add_energy( "O_V_loop_rms", rmsd_over_residue_list( current_pose, *native_poseCOP, rmsd_res_list, full_to_sub, is_prepend_map, false, true ) );

					if ( job_parameters_->gap_size() == 0 ){
						s.add_energy( "O_PBP_rmsd", phosphate_base_phosphate_rmsd( current_pose, *native_poseCOP, moving_base_residue,  false ) );
					} else{
						s.add_energy( "O_PBP_rmsd", 0.0 );
					}
				}

				////////Simple loop RMSD exclude only virtual atoms in native_pdb (mostly just the native virtual_res)///////
				core::pose::Pose curr_pose_no_variants = pose;

				// rhiju, 2013 -- I want the rmsd over non-virtual atoms!! Can't strip off virtuals.
				//remove_all_variant_types( curr_pose_no_variants ); //This removes all virtual_atoms

				if ( working_native_alignment.size() != 0 ){ //user specify which residue to align with native.
					align_poses( curr_pose_no_variants, tag + "_no_variants", ( *native_poseCOP ), "native",  working_native_alignment );
				} else{ //default
					align_poses( curr_pose_no_variants, tag + "_no_variants", ( *native_poseCOP ), "native",  working_best_alignment );
				}

				s.add_energy( "NAT_rmsd", rmsd_over_residue_list( curr_pose_no_variants,
																													*native_poseCOP,
																													rmsd_res_list,
																													full_to_sub,
																													is_prepend_map,
																													false /*verbose*/,
																													true /*ignore_virtual_atom*/ ) );

			}
		}

		silent_file_data.write_silent_struct( s, silent_file, write_score_only );

	}

} //rna
} //enumerate
} //stepwise
} //protocols
