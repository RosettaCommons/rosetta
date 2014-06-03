// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNAResidueSampler
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/sampling/util.hh>
#include <protocols/stepwise/sampling/protein/util.hh>
#include <protocols/stepwise/sampling/rna/util.hh>
#include <protocols/stepwise/sampling/align/StepWisePoseAligner.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/full_model_info/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/farna/util.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/rna/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/copydofs/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/SequenceMapping.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.hh>
#include <utility/stream_util.hh>
#include <utility/file/file_sys_util.hh>

#include <iostream>
#include <fstream>
#include <sstream>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

using namespace core;
using numeric::conversions::degrees;
using numeric::conversions::radians;
using ObjexxFCL::string_of;
using utility::operator<<;

static numeric::random::RandomGenerator RG(539155021);  // <- Magic number, do not change it!
static basic::Tracer TR( "protocols.stepwise.util" );

//////////////////////////////////////////////////////////////////////////
//
// Currently a grab-bag of helper funcitons -- should eventually
// separate out into separate util.hh functions, or move into
// core/util.
//////////////////////////////////////////////////////////////////////////


namespace protocols {
namespace stepwise {
namespace sampling {

	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	Size
	make_cut_at_moving_suite( core::kinematics::FoldTree & fold_tree, Size const & moving_suite ){

		fold_tree.new_jump( moving_suite, moving_suite + 1, moving_suite );
		return find_jump_number_at_suite( fold_tree, moving_suite );

	}

	///////////////////////////////////////////////////////////////////////
	Size
	make_cut_at_moving_suite( pose::Pose & pose, Size const & moving_suite ){

		core::kinematics::FoldTree fold_tree = pose.fold_tree();

		Size jump_number( 0 );
		if ( fold_tree.is_cutpoint( moving_suite ) ){ // already a cutpoint there
			return find_jump_number_at_suite( fold_tree, moving_suite );
		} else {
			jump_number = make_cut_at_moving_suite( fold_tree, moving_suite );
		}

		pose.fold_tree( fold_tree );

		return jump_number;
	}

	///////////////////////////////////////////////////////////////////////
	Size
	find_jump_number_at_suite( kinematics::FoldTree const & fold_tree, Size const & moving_suite ){

		int const i( moving_suite ), j( moving_suite + 1 );
		for ( Size n = 1; n <= fold_tree.num_jump(); n++ ) {
			if ( fold_tree.upstream_jump_residue( n ) == i && fold_tree.downstream_jump_residue( n ) == j ) return n;
			if ( fold_tree.upstream_jump_residue( n ) == j && fold_tree.downstream_jump_residue( n ) == i ) return n;
		}

		utility_exit_with_message( "Problem with jump number" );
		return 0;
	}


	///////////////////////////////////////////////////////////////////////////////////////
	Size
	look_for_unique_jump_to_moving_res( kinematics::FoldTree const & fold_tree, Size const & i ){

		for ( Size n = 1; n <= fold_tree.num_jump(); n++ ) {
			if ( fold_tree.downstream_jump_residue( n ) == int( i )  ) {
				return n;
			}
		}


		Size num_jump( 0 ), jump_idx( 0 );
		for ( Size n = 1; n <= fold_tree.num_jump(); n++ ) {
			if ( fold_tree.upstream_jump_residue( n ) == int( i ) || fold_tree.downstream_jump_residue( n ) == int( i ) ) {
				jump_idx = n;
				num_jump++;
			}
		}
		runtime_assert( num_jump == 1 );
		return jump_idx;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	bool
	is_cutpoint_closed( pose::Pose const & pose, Size const seq_num ){
		runtime_assert( seq_num > 0 );
		runtime_assert( seq_num <= pose.total_residue() );
		if ( pose.residue( seq_num  ).has_variant_type( chemical::CUTPOINT_LOWER )  ){
			runtime_assert ( pose.residue( seq_num+1 ).has_variant_type( chemical::CUTPOINT_UPPER )  );
			return true;
		}
		return false;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	get_cutpoint_closed( pose::Pose const & pose ){
		utility::vector1< Size > cutpoint_closed;
		for (Size seq_num = 1; seq_num < pose.total_residue(); seq_num++) {
			if ( is_cutpoint_closed( pose, seq_num ) ) cutpoint_closed.push_back( seq_num );
		}
		return cutpoint_closed;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	bool
	is_close_chain_break(pose::Pose const & pose){
		for (Size seq_num = 1; seq_num < pose.total_residue(); seq_num++) {
			if ( is_cutpoint_closed( pose, seq_num ) ) return true;
		}
		return false;
	}

	///////////////////////////////////////////////////////////////////////////////
	void
	pdbslice( core::pose::Pose & new_pose,
						core::pose::Pose const & pose,
						utility::vector1< core::Size > const & slice_res )
	{
		using namespace core::pose;
		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::chemical::rna;

		new_pose.clear();

		for ( Size i = 1; i <= slice_res.size(); i++ ) {

			ResidueOP residue_to_add = pose.residue( slice_res[ i ] ).clone() ;

			if ( (i > 1 &&  ( slice_res[i] != slice_res[i-1] + 1 )) /*new segment*/ ||
					 residue_to_add->is_lower_terminus() ||
					 residue_to_add->has_variant_type( "N_ACETYLATION") ||
					 (i>1 && pose.fold_tree().is_cutpoint( slice_res[i-1] ) ) ){
				if( residue_to_add->is_RNA() && (i>1) && new_pose.residue_type(i-1).is_RNA() ){

					new_pose.append_residue_by_jump(  *residue_to_add, i-1,
																						chi1_torsion_atom( new_pose.residue(i-1) ),
																						chi1_torsion_atom( *residue_to_add ), true /*new chain*/ );
				} else {

					new_pose.append_residue_by_jump(  *residue_to_add, i-1, "", "", true /*new chain*/ );
				}
			} else {

				new_pose.append_residue_by_bond(  *residue_to_add  ) ;
			}
		}

		if ( full_model_info_defined( pose ) ){
			FullModelInfoOP full_model_info = const_full_model_info( pose ).clone_info();
			utility::vector1< Size > const & res_list = full_model_info->res_list();
			utility::vector1< Size > new_res_list;
			for ( Size n = 1; n <= new_pose.total_residue(); n++ ){
				new_res_list.push_back( res_list[ slice_res[ n ] ] );
			}
			full_model_info->set_res_list( new_res_list );
			set_full_model_info( new_pose, full_model_info );
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	pdbslice( core::pose::Pose & pose,
						utility::vector1< core::Size > const & slice_res ){

		pose::Pose mini_pose;
		pdbslice( mini_pose, pose, slice_res );
		pose = mini_pose;

	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/// @brief  Superimpose mod_pose onto ref_pose using the mapping of residues from
	/// mod_pose to ref_pose given by res_map. Simple wrapper around superimpose_pose using IDs.
	Real
	superimpose_pose(
									 pose::Pose & mod_pose,
									 pose::Pose const & ref_pose,
									 std::map< Size, Size > const & res_map )
	{
		id::AtomID_Map< id::AtomID > atom_ID_map = create_alignment_id_map( mod_pose, ref_pose, res_map );
		return scoring::superimpose_pose( mod_pose, ref_pose, atom_ID_map );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 	id::AtomID_Map< id::AtomID >
 	create_alignment_id_map(	pose::Pose const & mod_pose,
													pose::Pose const & ref_pose,
													utility::vector1< core::Size > const & superimpose_res ){

		std::map< core::Size, core::Size > res_map;

 		for ( Size seq_num = 1; seq_num <= mod_pose.total_residue(); ++seq_num ) {
			if ( !superimpose_res.has_value(seq_num) ) continue;
			res_map[ seq_num ] = seq_num;
		}

		return create_alignment_id_map( mod_pose, ref_pose, res_map );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 	id::AtomID_Map< id::AtomID >
 	create_alignment_id_map(	pose::Pose const & mod_pose,
													pose::Pose const & ref_pose,
													std::map< core::Size, core::Size > res_map ){

 		using namespace chemical;
 		using namespace protocols::stepwise::sampling::protein;
 		using namespace protocols::stepwise::sampling::rna;
 		using namespace core::id;

 		AtomID_Map< AtomID > atom_ID_map;
 		pose::initialize_atomid_map( atom_ID_map, mod_pose, BOGUS_ATOM_ID );

 		for ( Size seq_num = 1; seq_num <= mod_pose.total_residue(); ++seq_num ) {

			if ( mod_pose.residue( seq_num ).is_RNA() && res_map.find( seq_num ) != res_map.end() && res_map[ seq_num ] > 0) {
				// Parin please update this function!!! Can't we just superimpose over C4'?
				setup_suite_atom_id_map( mod_pose, ref_pose, seq_num,  res_map[ seq_num ], atom_ID_map);

			} else if ( mod_pose.residue( seq_num ).is_protein() ){ // superimpose over CA.
				setup_protein_backbone_atom_id_map( mod_pose, ref_pose, seq_num, res_map[ seq_num ], atom_ID_map); // This will superimpose over N, C-alpha, C
			}

 		}

 		return atom_ID_map;

 	}

	//////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size > const
	convert_to_working_res( utility::vector1< Size > const & res_vector,
													utility::vector1< Size > const & working_res
													) {

		if ( res_vector.size() == 0 ) return res_vector;

		if ( working_res.size() == 0 ) return res_vector;

		std::map< Size, Size > full_to_sub;
		for ( Size i = 1; i <= working_res.size(); i++ ) {
			full_to_sub[ working_res[ i ] ] = i;
		}

		utility::vector1< Size > convert_res_vector;

		for ( Size i = 1; i <= res_vector.size(); i++ ) {
			if ( full_to_sub.find( res_vector[ i ] ) == full_to_sub.end() ) continue;
			convert_res_vector.push_back( full_to_sub[ res_vector[ i ] ] );
		}

		return convert_res_vector;

	}


	///////////////////////////////////////////////////////////////////////////////
	// Currently only handles atom pair constraints.
	///////////////////////////////////////////////////////////////////////////////
	core::scoring::constraints::ConstraintSetOP
	constraint_set_slice( core::scoring::constraints::ConstraintSetOP & cst_set,
												utility::vector1< core::Size > const & slice_res,
												pose::Pose const & pose,
												pose::Pose const & full_pose )
	{

		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::id;

		ConstraintSetOP cst_set_new( new scoring::constraints::ConstraintSet );

		ConstraintCOPs csts( cst_set->get_all_constraints() );

		//		std::map< Size, Size > slice_map;
		//		for (Size i = 1; i <= slice_res.size(); i++) slice_map[ slice_res[ i ] ] = i;
		utility::vector1< Size > mapping( full_pose.total_residue(), 0);
		for (Size i = 1; i <= slice_res.size(); i++) mapping[ slice_res[ i ] ] = i;
		SequenceMappingOP smap = new SequenceMapping( mapping );

		for ( Size n = 1; n <= csts.size(); n++ ) {

			ConstraintCOP const & cst( csts[n] );
			ConstraintOP cst_new = cst->remapped_clone( full_pose, pose, smap );
			if ( cst_new ) {
				cst_set_new->add_constraint( cst_new );
				//				std::cout << "HEY CONSTRAINTS!!! "
				//									<< cst_new->atom(1).rsd() << " " << pose.residue_type( cst_new->atom(1).rsd() ).atom_name( cst_new->atom(1).atomno() )
				//									<< " to "
				//									<< cst_new->atom(2).rsd() << " " << pose.residue_type( cst_new->atom(2).rsd() ).atom_name( cst_new->atom(2).atomno() ) << std::endl;
			}

			// currently only defined for pairwise distance constraints,
			//  and coordinate constraints
			//			if ( cst->score_type() == atom_pair_constraint)  {

// 				Size const i = cst->atom( 1 ).rsd();
// 				Size const j = cst->atom( 2 ).rsd();
// 				//			Size const dist( shortest_path_in_fold_tree.dist( i , j ) );
// 				//			if ( dist  > separation_cutoff ) continue;

// 				if ( slice_map.find( i ) == slice_map.end()  ) continue;
// 				if ( slice_map.find( j ) == slice_map.end()  ) continue;

// 				std::cout << "CST MAP: " << i << " " << slice_map[ i] << "          " << j << " " << slice_map[ j ] << std::endl;

// 				std::string const & atom_name1 = full_pose.residue_type( i ).atom_name( cst->atom(1).atomno() );
// 				std::string const & atom_name2 = full_pose.residue_type( j ).atom_name( cst->atom(2).atomno() );

// 				AtomID atom1_new( named_atom_id_to_atom_id( NamedAtomID( atom_name1, slice_map[ i ] ), pose );
// 				AtomID atom2_new( named_atom_id_to_atom_id( NamedAtomID( atom_name2, slice_map[ j ] ), pose );

// 				ConstraintOP cst_new = new AtomPairConstraint( atom1_new, atom2_new,
// 																											 cst->get_func().clone() /*is this defined?*/, cst->score_type() );

//			if ( cst_new ) cst_set_new->add_constraint( cst_new );


				//			} else if ( cst->score_type() == coordinate_constraint)  {




		}


		std::cout << "NUM CONSTRAINTS " << cst_set_new->get_all_constraints().size() << " out of " <<
			csts.size() << std::endl;

		return cst_set_new;
	}

	///////////////////////////////////////////////////////////////////////
	std::string
	get_file_name( std::string const & silent_file, std::string const & tag )
	{
		int pos( silent_file.find( ".out" ) );
		runtime_assert( pos > -1 );
		std::string silent_file_sample( silent_file );
		silent_file_sample.replace( pos, 4, tag+".out" );
		return silent_file_sample;

	}

	///////////////////////////////////////////////////////////////////////
	void
	check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose::Pose const & pose,
																															 core::scoring::ScoreFunctionOP & scorefxn ){

		using namespace core::scoring;

		if ( pose.constraint_set()->has_constraints() )	{
			if ( scorefxn->has_zero_weight( atom_pair_constraint ) ||
					 scorefxn->has_zero_weight( coordinate_constraint ) ) {
				utility_exit_with_message( "Since we want constraints, need to use a scorefunction with non-zero atom_pair_constraint and coordinate_constraint weight");
			}
		}

	}

	///////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	merge_vectors( utility::vector1< Size > const & vec1,
								 utility::vector1< Size > const & vec2 ){

		std::map< Size, bool > silly_map;
		for ( Size n = 1; n <= vec1.size(); n++ ) silly_map[ vec1[n] ] = true;
		for ( Size n = 1; n <= vec2.size(); n++ ) silly_map[ vec2[n] ] = true;

		utility::vector1< Size > merged_vec;
		for ( std::map<Size,bool>::iterator it = silly_map.begin(); it != silly_map.end(); it++ ){
			merged_vec.push_back( it->first );
		}
		return merged_vec;

	}

	///////////////////////////////////////////////////////////////////////
	void
	create_euler_rotation(
												Matrix & M,
												Real const & alpha,
												Real const & beta,
												Real const & gamma,
												Vector const & /* axis1 not actually used*/,
												Vector const & axis2,
												Vector const & axis3
												)
	{
		// Z-axis assumed to be long axis.
		Matrix M1 = numeric::rotation_matrix( axis3, Real( radians( alpha ) ) );
		Matrix M2 = numeric::rotation_matrix( axis2, Real( radians( beta ) ) );
		Matrix M3 = numeric::rotation_matrix( axis3, Real( radians( gamma ) ) );

		M = M3 * M2 * M1;
	}

	///////////////////////////////////////////////////////////////////////
	void
	create_euler_rotation(
												Matrix & M,
												Real const & alpha,
												Real const & beta,
												Real const & gamma )
	{
		static Vector const axis1( 1.0, 0.0, 0.0 );
		static Vector const axis2( 0.0, 1.0, 0.0 );
		static Vector const axis3( 0.0, 0.0, 1.0 );
		create_euler_rotation( M, alpha, beta, gamma, axis1, axis2, axis3 );
	}

	//////////////////////////////////////////////////////////////////////////////////
	void
	get_euler_angles( Real & alpha, Real & beta, Real & gamma, Matrix M1, Matrix M2, bool const verbose /*=true*/ ){

		Matrix M_test = M2;

		// Figure out what axis system2 looks like in axis system1.
		M2 = M1.transposed() * M2;

		// First figure out how to backrotate z rotation.
		Vector z_vec = M2.col_z();
		Real const gamma_radians = std::atan2( z_vec(2), z_vec(1) );
		gamma = degrees( gamma_radians );

		M2 = rotation_matrix( Vector( 0.0, 0.0, 1.0 ), -1.0 * gamma_radians ) * M2;
		z_vec = M2.col_z();
		if ( verbose ) std::cout << "This better have a zero in y position " << z_vec(1) << ' ' << z_vec(2) << ' ' << z_vec(3) << std::endl;

		// Then figure out how to backrotate y rotation.
		Real const beta_radians = std::atan2( z_vec(1), z_vec(3) );
		beta = degrees( beta_radians );

		M2 = rotation_matrix( Vector( 0.0, 1.0, 0.0 ), -1.0 * beta_radians ) * M2;
		z_vec = M2.col_z();
		if ( verbose ) std::cout << "This better have a zero in x and y position " << z_vec(1) << ' ' << z_vec(2) << ' ' << z_vec(3) << std::endl;

		// Finally, backrotate z rotation.
		Vector x_vec = M2.col_x();
		Real const alpha_radians = std::atan2( x_vec(2), x_vec(1) );
		alpha = degrees( alpha_radians );

		M2 = rotation_matrix( Vector( 0.0, 0.0, 1.0 ), -1.0 * alpha_radians ) * M2;
		x_vec = M2.col_x();
		if ( verbose ) std::cout << "This better have a zero in y and z position " << x_vec(1) << ' ' << x_vec(2) << ' ' << x_vec(3) << std::endl;

		if ( verbose ){
			Matrix M;

			Vector const xaxis1 = M1.col_x();
			Vector const yaxis1 = M1.col_y();
			Vector const zaxis1 = M1.col_z();
			create_euler_rotation( M, alpha, beta, gamma, xaxis1, yaxis1, zaxis1 );

			//	Matrix M_test;
			//	create_euler_rotation( M_test, alpha, beta, gamma, Vector( 1.0,0.0,0.0), Vector( 0.0,1.0,0.0), Vector(0.0,0.0,1.0) );
			//	M_test = M1 * M_test * M1.transposed();
			//	M_test = M1.transposed() * M_test * M1; // Can we rotate M1 into M2?
			M = M * M1;

			std::cout << "These better match:  " << std::endl;
			std::cout << M(1,1) <<  ' ' << M(1,2)  << ' ' << M(1,3) << std::endl;
			std::cout << M(2,1) <<  ' ' << M(2,2)  << ' ' << M(2,3) << std::endl;
			std::cout << M(3,1) <<  ' ' << M(3,2)  << ' ' << M(3,3) << std::endl;
			std::cout << std::endl;
			std::cout << M_test(1,1) <<  ' ' << M_test(1,2) <<  ' ' << M_test(1,3) << std::endl;
			std::cout << M_test(2,1) <<  ' ' << M_test(2,2)  << ' ' << M_test(2,3) << std::endl;
			std::cout << M_test(3,1) <<  ' ' << M_test(3,2)  << ' ' << M_test(3,3) << std::endl;
		}

	}

///////////////////////////////////////////////////////////////////////
void
translate( pose::Pose & pose, Vector const shift,
					 pose::Pose const & ref_pose,
					 utility::vector1< Size > const & moving_res ){

	using namespace core::id;

	for ( Size n = 1; n <= moving_res.size(); n++ ) {
		Size const i = moving_res[ n ];

		for ( Size m = 1; m <= pose.residue_type( i ).natoms(); m++ ) {
			pose.set_xyz( AtomID(m,i),   ref_pose.xyz( AtomID(m,i) ) + shift );
		}

	}

}


///////////////////////////////////////////////////////////////////////
void
rotate( pose::Pose & pose, Matrix const M,
				pose::Pose const & ref_pose,
				utility::vector1< Size > const & moving_res,
				Vector const & centroid ){

	using namespace core::id;

	for ( Size n = 1; n <= moving_res.size(); n++ ) {
		Size const i = moving_res[ n ];

		for ( Size m = 1; m <= pose.residue_type( i ).natoms(); m++ ) {
			pose.set_xyz( AtomID(m,i),  M * ( ref_pose.xyz( AtomID(m,i) ) - centroid ) + centroid );
		}
	}

}

//////////////////////////////////////////////////////////////////
void
rotate( pose::Pose & pose, Matrix const M,
				pose::Pose const & ref_pose,
				utility::vector1< Size > const & moving_res ){

	Vector centroid( 0.0, 0.0, 0.0 );
	rotate( pose, M, ref_pose, moving_res, centroid );
}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_base_centroid_and_rotation_matrix( pose::Pose const & pose, Size const i, Vector & centroid, Matrix & M ){

		using namespace scoring::rna;
		using namespace kinematics;
		static RNA_CentroidInfo rna_centroid_info;

		centroid = rna_centroid_info.get_base_centroid( pose.residue( i ) );
		Stub s = rna_centroid_info.get_base_coordinate_system( pose.residue( i ), centroid );
		M = s.M;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	translate_and_rotate_residue_to_origin( pose::Pose & pose, Size const i, utility::vector1< Size > const moving_res, bool const do_not_rotate  ){
		using namespace protocols::stepwise;

		Vector centroid;
		Matrix M;
		get_base_centroid_and_rotation_matrix( pose, i, centroid, M);

		translate( pose, -centroid, pose, moving_res);
		if ( !do_not_rotate )	rotate( pose, M.transposed(), pose, moving_res);

	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	translate_and_rotate_residue_to_origin( pose::Pose & pose, Size const i ){
		translate_and_rotate_residue_to_origin( pose, i, utility::tools::make_vector1( i ) );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// minimal wrapper around PoseAligner.
	void
	superimpose_at_fixed_res( pose::Pose & pose, pose::Pose const & native_pose,
														Real & rmsd, Size & natoms_rmsd,
														bool skip_bulges ){

		sampling::align::StepWisePoseAligner pose_aligner( native_pose );
		pose_aligner.set_skip_bulges( skip_bulges );
		pose_aligner.apply( pose );
		rmsd = pose_aligner.rmsd();
		natoms_rmsd = pose_aligner.natoms_rmsd();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Real
	superimpose_at_fixed_res_and_get_all_atom_rmsd( pose::Pose & pose, pose::Pose const & native_pose,
																									bool skip_bulges /* = false */) {
		Real rmsd;
		Size natoms_rmsd;
		superimpose_at_fixed_res( pose, native_pose, rmsd, natoms_rmsd, skip_bulges );
		return rmsd;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// calculate rmsd for pose and any 'other poses', and add up in quadrature.
	void
	superimpose_recursively( pose::Pose & pose, pose::Pose const & native_pose,
													 Real & rmsd, Size & natoms, bool skip_bulges = false ){

		using namespace core::pose;
		using namespace core::pose::full_model_info;

		Real rmsd_pose;
		Size natoms_pose;
		superimpose_at_fixed_res( pose, native_pose, rmsd_pose, natoms_pose, skip_bulges );

		Real const total_sd = ( rmsd * rmsd * natoms) + (rmsd_pose * rmsd_pose * natoms_pose );
		natoms += natoms_pose;
		if ( natoms > 0 ) {
			rmsd = std::sqrt( total_sd / Real( natoms ) );
		} else {
			runtime_assert( std::abs( rmsd ) < 1e-5 );
		}

		utility::vector1< PoseOP > const & other_pose_list = nonconst_full_model_info( pose ).other_pose_list();
		for ( Size n = 1; n <= other_pose_list.size(); n++ ){
			superimpose_recursively( *( other_pose_list[ n ] ), native_pose, rmsd, natoms, skip_bulges );
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// calculate rmsd for pose and any 'other poses', and add up in quadrature.
	Real
	superimpose_recursively( pose::Pose & pose, pose::Pose const & native_pose ){
		Real rmsd( 0.0 );
		Size natoms( 0 );
		superimpose_recursively( pose, native_pose, rmsd, natoms, false /*skip bulges*/ );
		return rmsd;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	get_number_missing_residue_connections( pose::Pose & pose ) {

		using namespace core::pose::full_model_info;
		utility::vector1< Size > const pose_domain_map = figure_out_pose_domain_map( pose );
		utility::vector1< Size > const & cutpoint_open = const_full_model_info( pose ).cutpoint_open_in_full_model();
		utility::vector1< Size > const & fixed_domain_map = const_full_model_info( pose ).fixed_domain_map();

		Size nmissing( 0 );
		Size const nres = pose_domain_map.size();
		for ( Size n = 1; n <= nres; n++ ){
			if ( fixed_domain_map[ n ] == 0 ){
				if ( pose_domain_map[ n ] == 0 ) nmissing++;
			} else if ( n < nres &&
									!cutpoint_open.has_value( n ) &&
									fixed_domain_map[ n+1 ] > 0 &&
									fixed_domain_map[ n+1 ] != fixed_domain_map[ n ] ) {
				if ( pose_domain_map[ n ] == 0 ||
						 pose_domain_map[ n+1 ] == 0 ||
						 ( pose_domain_map[ n ] != pose_domain_map[ n+1 ] ) ) {
					nmissing++;
				}
			}
		}

		return nmissing;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	is_at_terminus( core::pose::Pose const & pose, Size const i ){

		if ( i == pose.total_residue() || pose.fold_tree().is_cutpoint( i ) ){ // could be a 3' chain terminus
			return true;
		} else if ( i == 1 || pose.fold_tree().is_cutpoint( i-1 ) ) {
			return true;
		}

		return false;
	}


	////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_jump_partners_from_pose( utility::vector1< Size > & jump_partners1,
															 utility::vector1< Size > & jump_partners2,
															 utility::vector1< std::string > & jump_atoms1,
															 utility::vector1< std::string > & jump_atoms2,
															 pose::Pose const & pose,
															 utility::vector1< Size > const & working_res ){

		using namespace core::kinematics;

		FoldTree const & f = pose.fold_tree();
		for ( Size n = 1; n <= f.num_jump(); n++ ){
			Size const j1 = f.upstream_jump_residue( n );
			Size const j2 = f.downstream_jump_residue( n );
			jump_partners1.push_back(  working_res[ j1 ] );
			jump_partners2.push_back(  working_res[ j2 ] );
			jump_atoms1.push_back( f.upstream_atom( n ) );
			jump_atoms2.push_back( f.downstream_atom( n ) );
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	void
	get_endpoints_from_pose( utility::vector1< Size > & endpoints,
													 pose::Pose const & pose,
													 utility::vector1< Size > const & working_res ){

		for (Size i = 1; i < pose.total_residue(); i++ ){
			if ( pose.fold_tree().is_cutpoint( i ) ){
				endpoints.push_back( working_res[ i ] );
			}
		}
		endpoints.push_back( working_res[ pose.total_residue() ] );

	}


	////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	map_to_local_numbering( utility::vector1< Size > const & vec,
													utility::vector1< Size > const & working_res ){
		utility::vector1< Size > vec_new;
		for ( Size n = 1; n <= vec.size(); n++ ){
			vec_new.push_back( working_res.index( vec[n] ) );
		}
		return vec_new;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	std::map< Size, Size >
	get_res_map( utility::vector1< Size > const & working_res,
							 utility::vector1< Size > const & source_res ){

		std::map< Size, Size > res_map;

		for ( Size n = 1; n <= working_res.size(); n++ ){
			if ( !source_res.has_value( working_res[n] ) ) continue;
			res_map[ n ] = source_res.index( working_res[n] );
		}

		return res_map;
	}


	////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	merge_disjoint_vectors( utility::vector1< Size > const & res_vector1,
													utility::vector1< Size > const & res_vector2 ){

		utility::vector1< Size > res_vector = res_vector1;

		for ( Size n = 1; n <= res_vector2.size(); n++ ) {
			runtime_assert( ! res_vector.has_value( res_vector2[n] ) );
			res_vector.push_back( res_vector2[ n ] );
		}

		std::sort( res_vector.begin(), res_vector.end() ); // hope this works.

		return res_vector;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// this could go into Residue.hh itself
	void
	remove_upper_terminus( core::conformation::ResidueOP & rsd ){
		using namespace chemical;
		using namespace conformation;
		if ( rsd->has_variant_type( UPPER_TERMINUS ) ){
			ResidueTypeSet const & rsd_set( rsd->residue_type_set() );
			ResidueType const & new_rsd_type( rsd_set.get_residue_type_with_variant_removed( rsd->type(), UPPER_TERMINUS ) );
			rsd = new Residue( new_rsd_type, true );
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	void
	remove_lower_terminus( core::conformation::ResidueOP & rsd ){
		using namespace chemical;
		using namespace conformation;
		if ( rsd->has_variant_type( LOWER_TERMINUS ) ){
			ResidueTypeSet const & rsd_set( rsd->residue_type_set() );
			ResidueType const & new_rsd_type( rsd_set.get_residue_type_with_variant_removed( rsd->type(), LOWER_TERMINUS ) );
			rsd = new Residue( new_rsd_type, true );
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	find_root_without_virtual_ribose( kinematics::FoldTree const & f, pose::Pose const & pose ){
		for ( Size n = 1; n <= f.nres(); n++ ){
			if ( f.possible_root( n ) && !pose.residue_type( n ).has_variant_type( "VIRTUAL_RIBOSE" ) ){
				return n;
			}
		}
		utility_exit_with_message( "Fail: find_root_without_virtual_ribose" );
		return 0;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	definite_terminal_root( pose::Pose const & pose, Size const i ){
		using namespace core::pose::full_model_info;
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & res_list = full_model_info.res_list();
		utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
		if ( res_list[ i ] == 1 ||
				 cutpoint_open_in_full_model.has_value( res_list[ i ] - 1 ) ){ // great, nothing will ever get prepended here.
			return true;
		}
		if ( res_list[ i ] == full_model_info.size() ||
				 cutpoint_open_in_full_model.has_value( res_list[ i ] ) ){ // great, nothing will ever get appended here.
			return true;
		}
		return false;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	void
	try_reroot_at_fixed_domain( pose::Pose & pose ){

		using namespace core::pose::full_model_info;
		kinematics::FoldTree f = pose.fold_tree();
		Size new_root( 0 );
		for ( Size n = 1; n <= f.nres(); n++ ){
			if ( definite_terminal_root( pose, n ) ){
				new_root = n; break;
			}
		}
		if ( new_root > 0 ){
			f.reorder( new_root );
			pose.fold_tree( f );
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	find_first_root_residue( kinematics::FoldTree const & f,
													 utility::vector1< Size > const & working_res_subset,
													 utility::vector1< Size > const & working_res ){
		runtime_assert( f.nres() == working_res.size() );
		bool found_root( false );
		for ( Size n = 1; n <= f.nres(); n++ ){
			if ( working_res_subset.has_value( working_res[ n ] ) &&
					 f.possible_root( n ) ){
				found_root = true;
				return n;
			}
		}
		runtime_assert( found_root );
		return 0;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	void
	reroot_based_on_full_model_info( pose::Pose & pose,
																	 utility::vector1< Size > const & root_partition_res ){
		using namespace core::kinematics;
		Size new_root( 0 ), possible_root( 0 );
		for ( Size n = 1; n <= root_partition_res.size(); n++ ){
			Size const i = root_partition_res[ n ];
			if ( !pose.fold_tree().possible_root( i ) ) continue;
			if ( definite_terminal_root( pose, i ) ) {
				new_root = i; break;
			}
			if ( possible_root == 0 ) possible_root = i; // not as desirable, but sometimes necessary
		}
		if ( new_root == 0 ){
			if ( possible_root == 0 ){
				std::cerr << pose.fold_tree() << std::endl;
			}
			runtime_assert( possible_root > 0 );
			new_root = possible_root;
		}

		FoldTree f = pose.fold_tree();
		if ( static_cast<int>(new_root) == f.root() ) return;
		f.reorder( new_root );
		pose.fold_tree( f );
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	void
	merge_in_other_pose_by_jump( pose::Pose & pose, pose::Pose const & pose2,
															 Size const lower_merge_res, Size const upper_merge_res  ){
		runtime_assert( lower_merge_res < upper_merge_res );
		merge_in_other_pose( pose, pose2,
												 lower_merge_res, upper_merge_res,
												 false /*connect_residues_by_bond*/ );
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	void
	merge_in_other_pose_by_bond( pose::Pose & pose, pose::Pose const & pose2, Size const merge_res ){
		Size const lower_merge_res( merge_res ), upper_merge_res( merge_res + 1 );
		merge_in_other_pose( pose, pose2,
												 lower_merge_res, upper_merge_res,
												 true /*connect_residues_by_bond*/ );
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	void
	merge_in_other_pose( pose::Pose & pose, pose::Pose const & pose2,
											 Size const lower_merge_res, Size const upper_merge_res,
											 bool const connect_residues_by_bond ) {

		using namespace core::pose::datacache;
		using namespace core::pose::full_model_info;

		if ( connect_residues_by_bond ) runtime_assert( upper_merge_res == lower_merge_res + 1 );

		FullModelInfo & full_model_info = nonconst_full_model_info( pose );

		core::pose::Pose pose_scratch;
		utility::vector1< Size > const new_res_list = merge_two_poses_using_full_model_info( pose_scratch, pose, pose2,
																																												 lower_merge_res, upper_merge_res,
																																												 connect_residues_by_bond );
		pose.conformation() = pose_scratch.conformation();

		full_model_info.set_res_list( new_res_list );
		//		pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info );
		update_pdb_info_from_full_model_info( pose ); // for output pdb or silent file -- residue numbering.

	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	merge_two_poses_using_full_model_info( pose::Pose & pose,
																				 pose::Pose const & pose1,
																				 pose::Pose const & pose2,
																				 Size const lower_merge_res,
																				 Size const upper_merge_res,
																				 bool const connect_residues_by_bond ) {

		using namespace core::pose::full_model_info;
		using namespace core::pose::datacache;
		using namespace basic::datacache;

		// get working_residue information from each pose
		utility::vector1< Size > const & working_res1 = get_res_list_from_full_model_info_const( pose1 );
		utility::vector1< Size > const & working_res2 = get_res_list_from_full_model_info_const( pose2 );

		if ( connect_residues_by_bond ) runtime_assert( upper_merge_res == lower_merge_res + 1 );

		return merge_two_poses( pose, pose1, pose2, working_res1, working_res2,
														lower_merge_res, upper_merge_res, connect_residues_by_bond );

	}


	////////////////////////////////////////////////////////////////////////////////////////////////
	// Note that following is really general, and potentially very valuable, for all stepwise approaches
	//  and even in pose setup for fragment assembly of RNA or proteins.
	// Consider including in a GeneralPoseSetup class.
	//
	// Trying to get correct fold-tree handling, and variants.
	// Also trying to properly handle sequence reorderings (which can get pretty complicated in general )
	//
	utility::vector1< Size >
	merge_two_poses( pose::Pose & pose,
									 pose::Pose const & pose1,
									 pose::Pose const & pose2,
									 utility::vector1< Size > const & working_res1,
									 utility::vector1< Size > const & working_res2,
									 Size const lower_merge_res,
									 Size const upper_merge_res,
									 bool const connect_residues_by_bond,
									 bool const fix_first_pose /* = true */) {

		using namespace core::kinematics;
		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::pose::copydofs;

		if ( connect_residues_by_bond ) runtime_assert( upper_merge_res == lower_merge_res + 1 );

		if ( working_res1.has_value( lower_merge_res ) ){
			runtime_assert( working_res2.has_value( upper_merge_res ) );
		} else {
			runtime_assert( working_res2.has_value( lower_merge_res   ) );
			runtime_assert( working_res1.has_value( upper_merge_res ) );
			TR.Debug << "merge_two_poses: order is switched.  " << std::endl;
			// order is switched. No problem. Well... this could probably be avoided if we instead
			// specify merge_res1 and merge_res2, and allow them to be in arbitrary sequential order.
			return merge_two_poses( pose, pose2, pose1, working_res2, working_res1,
															lower_merge_res, upper_merge_res, connect_residues_by_bond,
															! fix_first_pose );
		}

		// define full working res -- union of the working res of the individual poses.
		utility::vector1< Size > const working_res = merge_disjoint_vectors( working_res1, working_res2 );

		for ( Size n = 1; n <= working_res.size(); n++ ){
			Size const k = working_res[ n ];
			ResidueOP rsd;
			bool after_cutpoint( false );
			if ( working_res1.has_value( k ) ) {
				Size const j = working_res1.index( k );
				if ( j == 1 || pose1.fold_tree().is_cutpoint( j-1 ) ) after_cutpoint = true;
				rsd = pose1.residue( j ).clone();
			} else {
				runtime_assert( working_res2.has_value( k ) );
				Size const j = working_res2.index( k );
				if ( j == 1 || pose2.fold_tree().is_cutpoint( j-1 ) ) after_cutpoint = true;
				rsd = pose2.residue( j ).clone();
			}

			if ( connect_residues_by_bond && k == lower_merge_res ) {
				remove_upper_terminus( rsd );
				rsd = remove_variant_type_from_residue( *rsd, "THREE_PRIME_PHOSPHATE", pose ); // got to be safe.
				rsd = remove_variant_type_from_residue( *rsd, "C_METHYLAMIDATION", pose ); // got to be safe.
			}
			if ( connect_residues_by_bond && k == upper_merge_res ) {
				runtime_assert( after_cutpoint );
				remove_lower_terminus( rsd );
				rsd = remove_variant_type_from_residue( *rsd, "FIVE_PRIME_PHOSPHATE", pose ); // got to be safe.
				rsd = remove_variant_type_from_residue( *rsd, "N_ACETYLATION", pose ); // got to be safe.
				after_cutpoint = false; // we're merging after all.
			}
			if ( n == 1 || !after_cutpoint ){
				pose.append_residue_by_bond( *rsd, true /* build_ideal_geometry */ );
			} else {
				pose.append_residue_by_jump( *rsd, pose.total_residue() );
			}
		}

		//////////////////////////////////////////////////////////////////////////////////
		// figure out fold tree for this merged pose --> move this to its own function.

		// figure out jumps
		utility::vector1< Size > jump_partners1, jump_partners2, endpoints, cuts;
		utility::vector1< std::string > jump_atoms1, jump_atoms2;
		get_jump_partners_from_pose( jump_partners1, jump_partners2, jump_atoms1, jump_atoms2, pose1, working_res1 );
		get_jump_partners_from_pose( jump_partners1, jump_partners2, jump_atoms1, jump_atoms2, pose2, working_res2 );
		if ( !connect_residues_by_bond ){
			jump_partners1.push_back( lower_merge_res );
			jump_partners2.push_back( upper_merge_res );
			jump_atoms1.push_back( chemical::rna::default_jump_atom( pose1.residue( working_res1.index( lower_merge_res ) ) ) ); // rna specific?
			jump_atoms2.push_back( chemical::rna::default_jump_atom( pose2.residue( working_res2.index( upper_merge_res ) ) ) ); // rna specific?
		}
		runtime_assert( jump_partners1.size() == jump_partners2.size() );

		// figure out cuts
		get_endpoints_from_pose( endpoints, pose1, working_res1 );
		get_endpoints_from_pose( endpoints, pose2, working_res2 );
		if ( connect_residues_by_bond ) runtime_assert( endpoints.has_value( lower_merge_res ) );
		Size const last_res = working_res[ working_res.size() ];
		runtime_assert( endpoints.has_value( last_res ) );
		for ( Size n = 1; n <= endpoints.size(); n++ ){
			if ( ( !connect_residues_by_bond || endpoints[ n ] != lower_merge_res ) &&
					 endpoints[ n ] != last_res ) cuts.push_back( endpoints[ n ] );
		}

		Size num_cuts = cuts.size();
		runtime_assert( num_cuts == jump_partners1.size() );

		jump_partners1 = map_to_local_numbering( jump_partners1, working_res );
		jump_partners2 = map_to_local_numbering( jump_partners2, working_res );
		cuts = map_to_local_numbering( cuts, working_res );

		FoldTree f = get_tree( pose.total_residue(), cuts, jump_partners1, jump_partners2, jump_atoms1, jump_atoms2 );

		Size root( 0 );
		if ( fix_first_pose ) {
			root = find_first_root_residue( f, working_res1, working_res );
		} else {
			root = find_first_root_residue( f, working_res2, working_res );
		}
		runtime_assert( f.reorder( root ) );

		pose.fold_tree( f );

		// map (internal) coordinates from separate poses into merged one.
		// this is potentially dangerous, if we accumulate floating point errors, and scores change.
		std::map< Size, Size > res_map1 = get_res_map( working_res, working_res1 );
		std::map< Size, Size > res_map2 = get_res_map( working_res, working_res2 );

		copy_dofs_match_atom_names( pose, pose1, res_map1, false /*backbone_only*/, false /*side_chain_only*/, false /*ignore_virtual*/ );
		copy_dofs_match_atom_names( pose, pose2, res_map2, false /*backbone_only*/, false /*side_chain_only*/, false /*ignore_virtual*/ );

		if ( fix_first_pose ){
			superimpose_pose( pose, pose1, res_map1 );
		} else {
			superimpose_pose( pose, pose2, res_map2 );
		}

		return working_res;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	apply_numbering( utility::vector1< Size > const & res,
									 utility::vector1< Size > const & numbering ) {
		utility::vector1< Size > res_renumbered;
		for ( Size n = 1; n <= res.size(); n++ ) {
			runtime_assert( res[ n ] > 0 && res[ n ] <= numbering.size() );
			res_renumbered.push_back( numbering[ res[ n ] ] );
		}
		return res_renumbered;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	get_other_residues( utility::vector1< Size > const & res, Size const & nres ){
		utility::vector1< Size > other_res;
		for ( Size n = 1; n <= nres; n++ ){
			if ( !res.has_value( n ) ) other_res.push_back( n );
		}
		return other_res;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// could/should also make a version slice_out_pose_and_update_full_model_info...
	void
	slice_out_pose( pose::Pose & pose,
									pose::Pose & sliced_out_pose,
									utility::vector1< Size > const & residues_to_delete ) {

		using namespace core::pose;
		using namespace core::pose::full_model_info;

		// need this for later.
		FullModelInfo & full_model_info = nonconst_full_model_info( pose );
		utility::vector1< Size > const original_res_list = full_model_info.res_list();

		// piece sliced out
		slice( sliced_out_pose, pose, residues_to_delete );

		FullModelInfoOP sliced_out_full_model_info = full_model_info.clone_info();
		sliced_out_full_model_info->set_res_list( apply_numbering( residues_to_delete, original_res_list ) );
		sliced_out_full_model_info->clear_other_pose_list();
		sliced_out_pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, sliced_out_full_model_info );
		update_pdb_info_from_full_model_info( sliced_out_pose ); // for output pdb or silent file -- residue numbering.

		// remainder piece.
		utility::vector1< Size > const residues_to_retain = get_other_residues( residues_to_delete, pose.total_residue() );

		Pose pose_scratch;
		slice( pose_scratch, pose, residues_to_retain );
		pose.conformation() = pose_scratch.conformation();

		full_model_info.set_res_list( apply_numbering( residues_to_retain, original_res_list )  );
		update_pdb_info_from_full_model_info( pose ); // for output pdb or silent file -- residue numbering.

		//		try_reroot_at_fixed_domain( pose );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// this may falter if the sliced out pose is not contiguous. That's what the boolean return value is for.
	// Could also return MovingResidueCase
	void
	slice( pose::Pose & sliced_out_pose,
				 pose::Pose const & pose,
				 utility::vector1< Size > const & slice_res ){

		using namespace core::kinematics;
		using namespace core::conformation;
		using namespace core::pose::copydofs;

		Size num_five_prime_connections( 0 ), num_three_prime_connections( 0 );
		Size num_jumps_to_previous( 0 ), num_jumps_to_next( 0 );

		sliced_out_pose.conformation().clear();

		for ( Size n = 1; n <= slice_res.size(); n++ ){

			Size const j = slice_res[ n ];
			ResidueOP rsd = pose.residue( j ).clone();

			bool const after_cutpoint = ( j == 1 ||	( pose.fold_tree().is_cutpoint( j - 1 ) ) );
			if ( after_cutpoint ) remove_lower_terminus( rsd );
			if ( !after_cutpoint && !slice_res.has_value( j - 1 ) ) num_five_prime_connections++;

			bool const jump_to_previous = check_jump_to_previous_residue_in_chain( pose, j, slice_res );
			if ( jump_to_previous ) num_jumps_to_previous++;

			bool const before_cutpoint = ( j == pose.total_residue() || ( pose.fold_tree().is_cutpoint( j ) ) );
			if ( before_cutpoint ) remove_upper_terminus( rsd );
			if ( !before_cutpoint && !slice_res.has_value( j + 1 ) ) num_three_prime_connections++;

			bool const jump_to_next = check_jump_to_next_residue_in_chain( pose, j, slice_res );
			if ( jump_to_next ) num_jumps_to_next++;

			if ( n == 1 || ( !after_cutpoint && slice_res.has_value( j - 1 ) ) ){
				sliced_out_pose.append_residue_by_bond( *rsd, true /* build_ideal_geometry */ );
			} else {
				sliced_out_pose.append_residue_by_jump( *rsd, sliced_out_pose.total_residue() );
			}
		}

		runtime_assert ( num_five_prime_connections <= 1 );
		runtime_assert ( num_three_prime_connections <= 1 );
		runtime_assert ( num_jumps_to_previous <= 1 );
		runtime_assert ( num_jumps_to_next <= 1 );
		//		TR << num_five_prime_connections << " " <<  num_three_prime_connections << " " << num_jumps_to_previous << " " <<  num_jumps_to_next << std::endl;
		// requirement for a clean slice:
		runtime_assert( (num_five_prime_connections + num_three_prime_connections + num_jumps_to_previous + num_jumps_to_next) == 1 );

		// fold tree!
		// figure out jumps
		utility::vector1< Size > jump_partners1, jump_partners2, endpoints, cuts;
		utility::vector1< std::string > jump_atoms1, jump_atoms2;
		FoldTree const & f = pose.fold_tree();
		for ( Size n = 1; n <= f.num_jump(); n++ ){
			Size const j1 = f.upstream_jump_residue( n );
			Size const j2 = f.downstream_jump_residue( n );
			if ( slice_res.has_value( j1 ) && slice_res.has_value( j2 ) ){
				jump_partners1.push_back( slice_res.index( j1 ) );
				jump_partners2.push_back( slice_res.index( j2 ) );
				jump_atoms1.push_back( f.upstream_atom( n ) );
				jump_atoms2.push_back( f.downstream_atom( n ) );
			}
		}

		for (Size n = 1; n < slice_res.size(); n++ ){
			Size const & i = slice_res[ n ];
			if ( !slice_res.has_value( i+1 ) || // boundary to non-sliced segment.
					 pose.fold_tree().is_cutpoint( i ) ) {  // cut internal to sliced segment.
				cuts.push_back( n );
			}
		}

		Size const num_cuts = cuts.size();
		runtime_assert( jump_partners1.size() == jump_partners2.size() );
		runtime_assert( num_cuts == jump_partners1.size() );

		FoldTree f_slice;
		ObjexxFCL::FArray2D< int > jump_point_( 2, num_cuts, 0 );
		ObjexxFCL::FArray1D< int > cuts_( num_cuts, 0 );
		for ( Size i = 1; i <= num_cuts; i++ ) {
			jump_point_( 1, i ) = std::min( jump_partners1[ i ], jump_partners2[ i ] );
			jump_point_( 2, i ) = std::max( jump_partners1[ i ], jump_partners2[ i ] );
			cuts_( i ) = cuts[ i ];
		}
		f_slice.tree_from_jumps_and_cuts( sliced_out_pose.total_residue(), num_cuts, jump_point_, cuts_ );

		// fix jump atoms.
		bool const KeepStubInResidue( true );
		for ( Size i = 1; i <= num_cuts; i++ ){
			Size const n = f_slice.jump_nr( jump_partners1[ i ], jump_partners2[ i ] );
			f_slice.set_jump_atoms( n,
												jump_partners1[ i ], jump_atoms1[ i ],
												jump_partners2[ i ], jump_atoms2[ i ], KeepStubInResidue );
		}
		if ( slice_res.has_value( pose.fold_tree().root() ) ) f_slice.reorder( slice_res.index( pose.fold_tree().root() ) );

		sliced_out_pose.fold_tree( f_slice );

		// map (internal) coordinates from separate poses into merged on.
		// this is potentially dangerous, if we accumulate floating point errors, and scores change.
		std::map< Size, Size > res_map;
		for ( Size n = 1; n <= slice_res.size(); n++ ) res_map[ n ] =  slice_res[ n ];
		copy_dofs_match_atom_names( sliced_out_pose, pose, res_map, false /*backbone_only*/, false /*side_chain_only*/, false /*ignore_virtual*/ );
		superimpose_pose( sliced_out_pose, pose, res_map );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	check_jump_to_previous_residue_in_chain( pose::Pose const & pose, Size const i,
																					 utility::vector1< Size > const & current_element,
																					 FullModelInfo const & full_model_info ){
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		utility::vector1< Size > const & chains_in_full_model = full_model_info.chains_in_full_model();
		return 	check_jump_to_previous_residue_in_chain( pose, i, current_element, res_list, chains_in_full_model );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	check_jump_to_previous_residue_in_chain( pose::Pose const & pose, Size const i,
																					 utility::vector1< Size > const & current_element ){
		utility::vector1< Size > res_list, chains_in_full_model;
		for ( Size n = 1; n <= pose.total_residue(); n++ ) {
			res_list.push_back( n );
			chains_in_full_model.push_back( 1 ); // could easily fix this for virtual residue.
		}
		return 	check_jump_to_previous_residue_in_chain( pose, i, current_element, res_list, chains_in_full_model );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	check_jump_to_previous_residue_in_chain( pose::Pose const & pose, Size const i,
																					 utility::vector1< Size > const & current_element,
																					 utility::vector1< Size > const & res_list,
																					 utility::vector1< Size > const & chains_in_full_model ){
		Size previous_res = i - 1;
		Size const current_chain = chains_in_full_model[ res_list[ i ] ];
		while ( previous_res >= 1 &&
						chains_in_full_model[ res_list[ previous_res ] ] == current_chain ){
			if ( !current_element.has_value( res_list[ previous_res ] ) &&
					 pose.fold_tree().jump_nr( previous_res, i ) > 0 ) return previous_res;
			previous_res--;
		}
		return 0;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	check_jump_to_next_residue_in_chain( pose::Pose const & pose, Size const i,
																					 utility::vector1< Size > const & current_element,
																					 FullModelInfo const & full_model_info ){
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		utility::vector1< Size > const & chains_in_full_model = full_model_info.chains_in_full_model();
		return 	check_jump_to_next_residue_in_chain( pose, i, current_element, res_list, chains_in_full_model );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	check_jump_to_next_residue_in_chain( pose::Pose const & pose, Size const i,
																					 utility::vector1< Size > const & current_element ){
		utility::vector1< Size > res_list, chains_in_full_model;
		for ( Size n = 1; n <= pose.total_residue(); n++ ) {
			res_list.push_back( n );
			chains_in_full_model.push_back( 1 ); // could easily fix this for virtual residue.
		}
		return 	check_jump_to_next_residue_in_chain( pose, i, current_element, res_list, chains_in_full_model );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	check_jump_to_next_residue_in_chain( pose::Pose const & pose, Size const i,
																						 utility::vector1< Size > const & current_element,
																						 utility::vector1< Size > const & res_list,
																						 utility::vector1< Size > const & chains_in_full_model ){
		Size subsequent_res = i + 1;
		Size const current_chain = chains_in_full_model[ res_list[ i ] ];
		while ( subsequent_res <= pose.total_residue() &&
						chains_in_full_model[ res_list[ subsequent_res ] ] == current_chain ){
			if ( !current_element.has_value( res_list[ subsequent_res ] ) &&
					 pose.fold_tree().jump_nr( subsequent_res, i ) > 0 ) return subsequent_res;
			subsequent_res++;
		}
		return 0;
	}


	////////////////////////////////////////////////////////////////////
	void
	fix_up_residue_type_variants_at_strand_end( pose::Pose & pose, Size const res ) {

		using namespace core::chemical;
		using namespace core::pose::full_model_info;

		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
		utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
		utility::vector1< Size > const chains_full = get_chains_full( pose );

		// Could this be a chainbreak (cutpoint_closed )?
		TR.Debug << "checking for cutpoint after append: " << res << " " << res_list[ res ]  << " " << cutpoint_open_in_full_model.size() << std::endl;

		if ( res < pose.total_residue() &&
				 res_list[ res ] + 1 == res_list[ res + 1 ] &&
				 ! cutpoint_open_in_full_model.has_value( res_list[ res ]) ){

			// can happen after additions
			core::pose::correctly_add_cutpoint_variants( pose, res );

			// leave virtual riboses in this should actually get instantiated by the modeler
			//	remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", res );

		} else{

			// can happen after deletions
			remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, res );

			// proteins...
			if ( pose.residue_type( res ).is_protein() ){
				if ( res_list[ res ] < full_model_info.size() &&
						 chains_full[ res_list[ res ] + 1 ] == chains_full[ res_list[ res ] ] &&
						 ( res == pose.total_residue()  || res_list[ res ] + 1 < res_list[ res + 1 ] ) ){
					remove_variant_type_from_pose_residue( pose, UPPER_TERMINUS, res );
					add_variant_type_to_pose_residue( pose, "C_METHYLAMIDATION", res );
				} else {
					remove_variant_type_from_pose_residue( pose, "C_METHYLAMIDATION", res );
					add_variant_type_to_pose_residue( pose, UPPER_TERMINUS, res );
				}
			}

		}

	}

	////////////////////////////////////////////////////////////////////
	void
	fix_up_residue_type_variants_at_strand_beginning( pose::Pose & pose, Size const res ) {

		using namespace core::chemical;
		using namespace core::pose::full_model_info;

		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
		utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
		utility::vector1< Size > const chains_full = get_chains_full( pose );

		// Could this be a chainbreak (cutpoint_closed )?

		TR.Debug << "checking for cutpoint after prepend: " << res << " " << res_list[ res ] << " " << cutpoint_open_in_full_model.size() << std::endl;

		if ( res > 1 &&
				 res_list[ res ] - 1 == res_list[ res - 1 ] &&
				 ! cutpoint_open_in_full_model.has_value( res_list[ res - 1 ])  ){

			// can happen after additions
			core::pose::correctly_add_cutpoint_variants( pose, res - 1 );


		} else {

			// can happen after additions
			if ( pose.residue_type( res ).is_RNA() &&
					 !pose.residue_type( res ).has_variant_type( "FIVE_PRIME_PHOSPHATE" ) ){
				add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", res );
			}

			// can happen after deletions
			remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, res );

			// proteins...
			if ( pose.residue_type( res ).is_protein() ){
				if ( res_list[ res ] > 1 &&
						 chains_full[ res_list[ res ] - 1 ] == chains_full[ res_list[ res ] ] &&
						 ( res == 1 || res_list[ res ] - 1 > res_list[ res - 1 ] ) ){
					remove_variant_type_from_pose_residue( pose, LOWER_TERMINUS, res );
					add_variant_type_to_pose_residue( pose, "N_ACETYLATION", res );
				} else {
					remove_variant_type_from_pose_residue( pose, "N_ACETYLATION", res );
					add_variant_type_to_pose_residue( pose, LOWER_TERMINUS, res );
				}
			}

		}
	}

	////////////////////////////////////////////////////////////////////
	void
	fix_up_residue_type_variants_at_floating_base( pose::Pose & pose, Size const res ) {

		if ( !pose.residue(res ).is_RNA() ) return;
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
		utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
		utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();

		if ( fixed_domain_map[ res_list[ res ] ] != 0 ) return;

		if ( res > 1 &&
				 res_list[ res ] - 1 == res_list[ res - 1 ] &&
				 ! cutpoint_open_in_full_model.has_value( res_list[ res - 1 ])  ) return;

		if ( res < pose.total_residue() &&
				 res_list[ res ] + 1 == res_list[ res + 1 ] &&
				 ! cutpoint_open_in_full_model.has_value( res_list[ res ]) ) return;

		if ( pose.residue_type( res ).has_variant_type( "FIVE_PRIME_PHOSPHATE" ) )  return;
		if ( pose.residue_type( res ).has_variant_type( "THREE_PRIME_PHOSPHATE" ) ) return;

		add_variant_type_to_pose_residue( pose, "VIRTUAL_RIBOSE", res );

	}


	////////////////////////////////////////////////////////////////////
	void
	fix_up_residue_type_variants( pose::Pose & pose_to_fix ) {

		using namespace core::chemical;
		pose::Pose pose = pose_to_fix; // costly, but prevents seg fault with graphics.

		for ( Size n = 1; n <= pose.total_residue(); n++ ){

			// Are we at a strand beginning?
			bool const at_strand_beginning = ( n == 1 || pose.fold_tree().is_cutpoint( n-1 ) );
			if ( at_strand_beginning ) {
				fix_up_residue_type_variants_at_strand_beginning( pose, n );
			} else { // make sure there is nothing crazy here
				remove_variant_type_from_pose_residue( pose, LOWER_TERMINUS, n );
				remove_variant_type_from_pose_residue( pose, VIRTUAL_PHOSPHATE, n );
				remove_variant_type_from_pose_residue( pose, "FIVE_PRIME_PHOSPHATE", n );
				remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", n );
				runtime_assert( !pose.residue_type( n ).has_variant_type( CUTPOINT_UPPER ) );
			}

			// Look for strand ends.
			bool const at_strand_end = ( n == pose.total_residue() || pose.fold_tree().is_cutpoint( n ) );
			if ( at_strand_end ) {
				fix_up_residue_type_variants_at_strand_end( pose, n );
			} else {
				remove_variant_type_from_pose_residue( pose, UPPER_TERMINUS, n );
				remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", n );
				remove_variant_type_from_pose_residue( pose, "THREE_PRIME_PHOSPHATE", n );
				runtime_assert( !pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) );
			}

			 // check for floating_base
			if ( at_strand_end && at_strand_beginning ) fix_up_residue_type_variants_at_floating_base( pose,  n );
		}

		// Just copying the conformation() makes sure that other objects (such as other_pose_list) don't get cloned --
		//  can be important if external functions are holding OPs to those objects.
		pose_to_fix.conformation() = pose.conformation();
		pose_to_fix.pdb_info( pose.pdb_info() ); // silly -- ensures that PDBInfo is not flagged as 'obsolete'.
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	fix_up_jump_atoms( pose::Pose & pose ){
		// appears necessary for protein
		core::kinematics::FoldTree f = pose.fold_tree();
		f.reassign_atoms_for_intra_residue_stubs(); // it seems silly that we need to do this separately.
		pose.fold_tree( f );
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	void
	fix_up_jump_atoms_and_residue_type_variants( pose::Pose & pose_to_fix ) {
	 	fix_up_jump_atoms( pose_to_fix );
	 	fix_up_residue_type_variants( pose_to_fix );
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	void
	switch_focus_to_other_pose( pose::Pose & pose, Size const & focus_pose_idx ){

		using namespace core::pose;
		using namespace core::pose::full_model_info;

		if ( focus_pose_idx == 0 ) return;

		utility::vector1< PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();

		// the original pose (or its clone, actually) is relieved of the
		// duty of holding information on other poses
		PoseOP original_pose_clone = pose.clone();
		utility::vector1< PoseOP > blank_pose_list;
		nonconst_full_model_info( *original_pose_clone ).set_other_pose_list( blank_pose_list );

		// need to shift focus to the other pose. It will now be responsible for holding
		// pointers to the other poses.
		PoseCOP other_pose = other_pose_list[ focus_pose_idx ];
		FullModelInfoOP new_full_model_info = const_full_model_info( *other_pose ).clone_info();

		// copy in pose list from original pose.
		utility::vector1< PoseOP > new_other_pose_list = new_full_model_info->other_pose_list();
		// for now -- but later could have pose 'trees', in which other_poses themselves have other_poses.
		runtime_assert( new_other_pose_list.size() == 0 );

		new_other_pose_list.push_back( original_pose_clone );
		for ( Size i = 1; i <= other_pose_list.size(); i++ ){
			if ( i == focus_pose_idx ) continue;
			new_other_pose_list.push_back( other_pose_list[ i ]->clone() );
		}
		new_full_model_info->set_other_pose_list( new_other_pose_list );

		// OK, now shift focus! Hope this works.
		pose = ( *other_pose ); // makes a copy.
		pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, new_full_model_info );

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	switch_focus_among_poses_randomly( pose::Pose & pose, scoring::ScoreFunctionOP scorefxn,
																		 bool force_switch /* = false */ ) {

		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace protocols::stepwise;

		Size const num_other_poses = const_full_model_info( pose ).other_pose_list().size();
		if ( force_switch ) runtime_assert( num_other_poses > 0 );

		Size const focus_pose_idx = force_switch ? RG.random_range( 1, num_other_poses ) : RG.random_range( 0, num_other_poses );
		if ( focus_pose_idx == 0 ) return false;

		Real const score_before_switch_focus = ( scorefxn != 0 ) ? (*scorefxn)( pose ) : 0.0;
		TR.Debug << TR.Green << "SWITCHING FOCUS! SWITCHING FOCUS! SWITCHING FOCUS! SWITCHING FOCUS! to: " << focus_pose_idx << TR.Reset << std::endl;
		switch_focus_to_other_pose( pose, focus_pose_idx );
		Real const score_after_switch_focus = ( scorefxn != 0 ) ? (*scorefxn)( pose ) : 0.0;

		// originally set threshold at 0.001, but triggered rare errors. At some point worth tracking down...
		if (  std::abs( score_before_switch_focus - score_after_switch_focus ) > 0.10 ){
			utility_exit_with_message( "Energy change after switching pose focus: " + string_of( score_before_switch_focus ) + " to " +string_of( score_after_switch_focus ) );
		}

		return true;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// specific for two helix test case. Sort of a unit test. Maybe I should make it a unit test.
	void
	test_merge_and_slice_with_two_helix_test_case( 	utility::vector1< core::pose::PoseOP > const & input_poses,
																									core::scoring::ScoreFunctionOP scorefxn ){

		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace protocols::stepwise;

		scorefxn->show( *( input_poses[ 1 ] ) );
		scorefxn->show( *( input_poses[ 2 ] ) );

		input_poses[1]->dump_pdb( "input_pose1.pdb" );
		input_poses[2]->dump_pdb( "input_pose2.pdb" );

		pose::Pose merged_pose;
		merge_two_poses_using_full_model_info( merged_pose, *(input_poses[ 2 ]), *( input_poses[ 1 ] ),
																					 12 /*lower_merge_res*/, 13 /*upper_merge_res*/, true /*connect_by_bond*/ );
		merged_pose.dump_pdb( "merged_pose.pdb" );
		scorefxn->show( merged_pose );

		pose::Pose sliced_out_pose;
		utility::vector1< Size > slice_res = const_full_model_info( *(input_poses[2]) ).res_list();
		slice_out_pose( merged_pose, sliced_out_pose, slice_res );
		merged_pose.dump_pdb( "pose_after_slice.pdb" );
		sliced_out_pose.dump_pdb( "sliced_out_pose.pdb" );
		scorefxn->show( merged_pose );
		scorefxn->show( sliced_out_pose );
	}


	///////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_moving_chain_break_res( pose::Pose const & pose, kinematics::MoveMap const & mm ) {

		utility::vector1< Size > moving_chainbreak_res;
		utility::vector1< Size > const cutpoint_closed = get_cutpoint_closed( pose );

		// kind of brute-force check over all moving residues. Can optimize if needed.
		for ( Size n = 1; n < pose.total_residue(); n++ ){

			if ( cutpoint_closed.has_value( n ) ){
				if ( moving_chainbreak_res.has_value( n ) ) continue;
				if (mm.get( id::TorsionID( n, id::BB, 5 ) ) ||
						mm.get( id::TorsionID( n, id::BB, 6 ) ) ||
						mm.get( id::TorsionID( n+1, id::BB, 1 ) ) ||
						mm.get( id::TorsionID( n+1, id::BB, 2 ) ) ||
						mm.get( id::TorsionID( n+1, id::BB, 3 ) ) ) moving_chainbreak_res.push_back( n );
			} else {

				if ( pose.fold_tree().is_cutpoint( n ) ) continue;

				if ( ! mm.get( id::TorsionID( n+1, id::BB, 1 ) ) ) continue;

				ObjexxFCL::FArray1D_bool partner1( pose.total_residue(), false );
				pose.fold_tree().partition_by_residue( n, partner1 );

				for ( Size m = 1; m <= cutpoint_closed.size(); m++ ){
					Size const k = cutpoint_closed[ m ];
					if ( moving_chainbreak_res.has_value( k ) ) continue;
					if ( partner1( k ) != partner1( k+1 ) ) moving_chainbreak_res.push_back( k );
				}
			}
		}

		return moving_chainbreak_res;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	check_for_fixed_domain( pose::Pose const & pose,
													utility::vector1< Size> const & partition_res ){
		utility::vector1< Size > const domain_map = core::pose::full_model_info::get_fixed_domain_from_full_model_info_const( pose );
		for ( Size i = 1; i <= partition_res.size(); i++ ){
			if ( domain_map[ partition_res[ i ] ] > 0 ) return true;
		}
		return false;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	primary_fixed_domain( pose::Pose const & pose,
											utility::vector1< Size> const & partition_res ){
		utility::vector1< Size > const domain_map = core::pose::full_model_info::get_fixed_domain_from_full_model_info_const( pose );
		Size d_min( 0 );
		for ( Size i = 1; i <= partition_res.size(); i++ ){
			Size const d = domain_map[ partition_res[ i ] ];
			if ( d > 0 && ( d_min == 0 || d < d_min ) ) d_min = d;
		}
		if ( d_min == 0 ) return 999; // big number, to signal no domain found.
		return d_min;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	check_for_fixed_domain( pose::Pose const & pose ){
		utility::vector1< Size > partition_res;
		for ( Size i = 1; i <= pose.total_residue(); i++ ) partition_res.push_back( i );
		return check_for_fixed_domain( pose, partition_res );
	}

	///////////////////////////////////////////////////////////////
	void
	make_variants_match( pose::Pose & pose,
											 pose::Pose const & reference_pose,
											 Size const & n,
											 chemical::VariantType const variant_type ){
		if ( reference_pose.residue( n ).has_variant_type( variant_type ) ){
			add_variant_type_to_pose_residue( pose, variant_type, n );
		} else {
			remove_variant_type_from_pose_residue( pose, variant_type, n );
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_moving_cutpoints_closed_from_moving_res( pose::Pose const & pose, Size const moving_res ){
		return figure_out_moving_cutpoints_closed( pose, figure_out_moving_partition_res( pose, moving_res ) );
	}

	//////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_moving_cutpoints_closed_from_moving_res( pose::Pose const & pose, utility::vector1< Size > const & moving_res_list ){
		utility::vector1< Size > moving_cutpoints;
		for ( Size n = 1; n <= moving_res_list.size(); n++ ){
			moving_cutpoints = merge_vectors( moving_cutpoints,  figure_out_moving_cutpoints_closed( pose, figure_out_moving_partition_res( pose, moving_res_list[n] ) ) );
		}
		return moving_cutpoints;
	}

	//////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_moving_cutpoints_closed( pose::Pose const & pose,
																			utility::vector1< Size > const & moving_partition_res ){
		utility::vector1< Size > cutpoints_closed, five_prime_chain_breaks, three_prime_chain_breaks, chain_break_gap_sizes;
		figure_out_moving_chain_breaks( pose, moving_partition_res, cutpoints_closed, five_prime_chain_breaks, three_prime_chain_breaks, chain_break_gap_sizes );
		return cutpoints_closed;
	}

	//////////////////////////////////////////////////////////////////////////////////
	void
	figure_out_moving_chain_breaks( pose::Pose const & pose,
																	utility::vector1< Size > const & moving_partition_res,
																	utility::vector1< Size > & cutpoints_closed,
																	utility::vector1< Size > & five_prime_chain_breaks,
																	utility::vector1< Size > & three_prime_chain_breaks,
																	utility::vector1< Size > & chain_break_gap_sizes ){

		utility::vector1< Size > const chains =	figure_out_chains_from_full_model_info_const( pose );
		utility::vector1< Size > const & res_list =	get_res_list_from_full_model_info_const( pose );
		Size five_prime_chain_break( 0 ), three_prime_chain_break( 0 );
		for ( Size n = 1; n < pose.total_residue(); n++ ){

			if ( !pose.fold_tree().is_cutpoint( n ) ) continue;

			// Skip virtual anchors
			if ( pose.residue( n ).aa() == core::chemical::aa_vrt ) continue;
			if ( pose.residue( n+1 ).aa() == core::chemical::aa_vrt ) continue;

			// must be in different partitions to qualify as 'moving'
			if ( moving_partition_res.has_value( n ) == moving_partition_res.has_value( n+1 ) ) continue;

			// must be in same chain to qualify as a chain break
			if ( chains[ n ] != chains[ n+1 ] ) continue;

			// rewind to non-virtual residue closest to chainbreak
			for ( five_prime_chain_break = n; five_prime_chain_break >= 1; five_prime_chain_break-- ){
				if ( !pose.residue_type( five_prime_chain_break ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ) break;
			}

			// fast-forward to non-virtual residue closest to chainbreak
			for ( three_prime_chain_break = n+1; three_prime_chain_break <= pose.total_residue(); three_prime_chain_break++ ){
				if ( !pose.residue_type( three_prime_chain_break ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ) break;
			}

			five_prime_chain_breaks.push_back( five_prime_chain_break );
			three_prime_chain_breaks.push_back( three_prime_chain_break );
			runtime_assert( res_list[ three_prime_chain_break ] > res_list[ five_prime_chain_break ] );
			Size const gap_size = res_list[ three_prime_chain_break ] - res_list[ five_prime_chain_break ] - 1;
			chain_break_gap_sizes.push_back( gap_size );
			if ( gap_size == 0 ){
				runtime_assert( five_prime_chain_break == n ); // no rewind past bulges
				runtime_assert( three_prime_chain_break == n + 1 ); // no fast forward past bulges
				runtime_assert( pose.residue_type( five_prime_chain_break ).has_variant_type( "CUTPOINT_LOWER" ) );
				runtime_assert( pose.residue_type( three_prime_chain_break ).has_variant_type( "CUTPOINT_UPPER" ) );
				cutpoints_closed.push_back( n );
			}

		}

	}

	///////////////////////////////////////////////////////////////////////
	utility::vector1< bool >
	get_partition_definition_by_jump( pose::Pose const & pose, Size const & jump_nr /*jump_number*/ ){
		ObjexxFCL::FArray1D<bool> partition_definition( pose.total_residue(), false );

		pose.fold_tree().partition_by_jump( jump_nr, partition_definition );

		//silly conversion. There may be a faster way to do this actually.
		utility::vector1< bool > partition_definition_vector1;
		for ( Size n = 1; n <= pose.total_residue(); n++ )	partition_definition_vector1.push_back( partition_definition(n) );

		return partition_definition_vector1;
	}

	///////////////////////////////////////////////////////////////////////
	utility::vector1< bool >
	get_partition_definition( pose::Pose const & pose, Size const & moving_suite ){

		ObjexxFCL::FArray1D<bool> partition_definition( pose.total_residue(), false );
		if ( moving_suite > 0 ) {
			pose.fold_tree().partition_by_residue( moving_suite, partition_definition );
		} else {
			partition_definition.dimension( pose.total_residue(), true );
		}

		//silly conversion. There may be a faster way to do this actually.
		utility::vector1< bool > partition_definition_vector1;
		for ( Size n = 1; n <= pose.total_residue(); n++ )	partition_definition_vector1.push_back( partition_definition(n) );

		return partition_definition_vector1;

	}


	///////////////////////////////////////////////////////////////////////////
	Size
	figure_out_reference_res_for_suite( pose::Pose const & pose, Size const moving_res ){
		bool connected_by_jump( false );
		Size reference_res =  pose.fold_tree().get_parent_residue( moving_res, connected_by_jump );
		runtime_assert( !connected_by_jump );
		return reference_res;
	}

	///////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_moving_partition_res_for_suite( pose::Pose const & pose,
																						 Size const moving_res,
																						 Size const reference_res ){
		Size const moving_suite = ( moving_res < reference_res ) ? moving_res : reference_res;
		utility::vector1< bool > partition_definition = get_partition_definition( pose, moving_suite );
		utility::vector1< Size > moving_partition_res;
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( partition_definition[ n ] == partition_definition[ moving_res ] ) moving_partition_res.push_back( n );
		}
		return moving_partition_res;
	}

	///////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_moving_partition_res_for_jump( pose::Pose const & pose,
																						Size const jump_nr ){
		utility::vector1< bool > partition_definition = get_partition_definition_by_jump( pose, jump_nr );
		utility::vector1< Size > moving_partition_res;
		Size const moving_res = pose.fold_tree().downstream_jump_residue( jump_nr );
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( partition_definition[ n ] == partition_definition[ moving_res ] ) moving_partition_res.push_back( n );
		}
		return moving_partition_res;
	}

	///////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_moving_partition_res( pose::Pose const & pose,
																	 Size const moving_res ){
		utility::vector1< Size > root_partition_res, moving_partition_res;
		figure_out_root_and_moving_partition_res( pose, moving_res, root_partition_res, moving_partition_res );
		return moving_partition_res;
	}

	///////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_moving_partition_res( pose::Pose const & pose,
																	 utility::vector1< Size > const & moving_res_list ){

		utility::vector1< Size > moving_partition_res_all;
		for( Size n = 1; n <= moving_res_list.size(); n++ ){
			if ( pose.residue_type( moving_res_list[n] ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ) continue; // legacy insanity from SWA 'dinucleotide' move.
			utility::vector1< Size > moving_partition_res = figure_out_moving_partition_res( pose, moving_res_list[n] );
			for ( Size i = 1; i <= moving_partition_res.size(); i++ ) {
				if ( !moving_partition_res_all.has_value( moving_partition_res[i] ) ) moving_partition_res_all.push_back( moving_partition_res[i] );
			}
		}
		std::sort( moving_partition_res_all.begin(), moving_partition_res_all.end() );
		return moving_partition_res_all;
	}

	///////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	figure_out_root_partition_res( pose::Pose const & pose,
																 utility::vector1< Size > const & moving_res_list ){
		utility::vector1< Size > root_partition_res;
		utility::vector1< Size > const moving_partition_res = figure_out_moving_partition_res( pose, moving_res_list );
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( !moving_partition_res.has_value( n ) )  root_partition_res.push_back( n ) ;
		}
		return root_partition_res;
	}

	/////////////////////////////////////////////////////////////////////////////
	void
	figure_out_root_and_moving_partition_res( pose::Pose const & pose, Size const moving_res,
																						utility::vector1< Size > & root_partition_res,
																						utility::vector1< Size > & moving_partition_res ) {

		FoldTree const & f = pose.fold_tree();
		Size const reference_res = f.get_parent_residue( moving_res );
		Size const jump_nr = f.jump_nr( moving_res, reference_res );
		if ( !jump_nr && reference_res > 0 ) runtime_assert( moving_res == reference_res + 1 || reference_res == moving_res + 1 );
		utility::vector1< bool > partition_definition = jump_nr ? get_partition_definition_by_jump( pose, jump_nr ) : get_partition_definition( pose, std::min( moving_res, reference_res ) );
		for ( Size n = 1; n <= pose.total_residue(); n++ ) {
			if ( partition_definition[ n ] == partition_definition[ moving_res ] ){
				moving_partition_res.push_back( n );
			} else {
				root_partition_res.push_back( n );
			}
		}
	}

	///////////////////////////////////////////////////////////////////
	bool
	revise_root_and_moving_res( pose::Pose & pose, Size & moving_res /* note that this can change too*/ ){
		if ( moving_res == 0 ) return false;

		utility::vector1< Size > root_partition_res, moving_partition_res;
		figure_out_root_and_moving_partition_res( pose, moving_res, root_partition_res, moving_partition_res );

		bool switch_moving_and_root_partitions = ( root_partition_res.size() < moving_partition_res.size() );
		if ( root_partition_res.size() == moving_partition_res.size() ){
			switch_moving_and_root_partitions = primary_fixed_domain( pose, moving_partition_res ) < primary_fixed_domain( pose, root_partition_res );
		}
		if ( switch_moving_and_root_partitions ){ // the way things should be:
			Size const moving_res_original = moving_res; // to check switching went OK.
			Size const reference_res = pose.fold_tree().get_parent_residue( moving_res );
			moving_res = reference_res;
			reroot_based_on_full_model_info( pose, moving_partition_res /* new root_parition_res*/  );
			runtime_assert( pose.fold_tree().get_parent_residue( moving_res ) == static_cast<int>( moving_res_original ) );
		} else {
			reroot_based_on_full_model_info( pose, root_partition_res );
		}

		return switch_moving_and_root_partitions;
	}



	/////////////////////////////////////////////////////////////////////////////////////
	bool
	revise_root_and_moving_res_list( pose::Pose & pose,
																	 utility::vector1< Size > & moving_res_list /* note that this can change too*/ ){

		if ( moving_res_list.size() == 0 ) return false; // maybe after a delete -- just minimize.
		if ( pose.residue_type( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) return false; // for ERRASER.

		// weird case, happens in some protein early moves -- get the root out of the moving res!
		if ( moving_res_list.has_value( pose.fold_tree().root() ) ){
			utility::vector1< Size > other_res = get_other_residues( moving_res_list, pose.total_residue() );
			if ( other_res.size() == 0 ) return false; // from scratch.
			reroot_based_on_full_model_info( pose, other_res );
		}

		// find connection point to 'fixed res'
		Size moving_res_at_connection( 0 ), reference_res( 0 );
		for ( Size n = 1; n <= moving_res_list.size(); n++ ){
			Size const & moving_res = moving_res_list[n];
			Size const parent_res = pose.fold_tree().get_parent_residue( moving_res );
			if ( moving_res_list.has_value( parent_res ) ) continue;
			runtime_assert( moving_res_at_connection == 0 ); // moving_res_list must be contiguous in fold tree.
			moving_res_at_connection = moving_res;
			reference_res            = parent_res;
			// break; // don't break, actually -- let's make sure there's only one connection point.
		}
		runtime_assert( moving_res_at_connection != 0 );

		if ( reference_res == 0 ) return false; // special case -- happens with 'from-scratch' moves.

		bool const is_jump = pose.fold_tree().jump_nr( reference_res, moving_res_at_connection ) > 0 ;
		bool switched_moving_and_root_partitions = revise_root_and_moving_res( pose, moving_res_at_connection );
		if ( is_jump && moving_res_list.size() != 1 ) runtime_assert( !switched_moving_and_root_partitions );

		// revise_root_and_moving_res() handled revision of single residue only... need to translate this
		// switch to the whole list of residues.
		if (  switched_moving_and_root_partitions ){
			utility::vector1< Size > const moving_res_list_original = moving_res_list;
			moving_res_list = utility::tools::make_vector1( moving_res_at_connection );
			if ( is_jump ){ // easy switch-a-roo. jump connected single residues, and switch them.
				runtime_assert( moving_res_at_connection /*was switched*/  == reference_res );
			} else {
				// might be more residues. Example:
				//
				//              M   M   M   <-- moving_res_list (original)
				// ROOT...- 3 - 4 - 5 - 6 - 7 - ....
				//            M   M   M     <-- defines these moving connections
				//
				// after re-rooting, now should be:
				//
				//          M   M   M       <-- moving_res_list (new)
				//     ...- 3 - 4 - 5 - 6 - 7 - .... NEW ROOT
				//            M   M   M     <-- defines the same moving connections
				//
				for ( Size n = 1; n <= moving_res_list_original.size(); n++ ){
					Size const & moving_res = moving_res_list_original[n];
					if ( moving_res_list_original.has_value( pose.fold_tree().get_parent_residue( moving_res ) ) ) moving_res_list.push_back( moving_res );
				}
			}
			runtime_assert( moving_res_list.size() == moving_res_list_original.size() );
		}

		return switched_moving_and_root_partitions;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	Size
	split_pose( pose::Pose & pose, Size const moving_res, Size const reference_res ){

		Size jump_at_moving_suite = pose.fold_tree().jump_nr( moving_res, reference_res );
		if ( !jump_at_moving_suite ){
			runtime_assert( (moving_res == reference_res + 1) || (moving_res == reference_res - 1) );
			Size const moving_suite = ( moving_res < reference_res ) ? moving_res : reference_res;
			runtime_assert( !pose.fold_tree().is_cutpoint( moving_suite ) );
			jump_at_moving_suite = make_cut_at_moving_suite( pose, moving_suite );
			if ( pose.residue_type( moving_suite + 1 ).is_RNA() ) add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, moving_suite+1 );
		}

		kinematics::Jump j = pose.jump( jump_at_moving_suite );
		j.set_translation( Vector( 1.0e4, 0.0, 0.0 ) );
		pose.set_jump( jump_at_moving_suite, j );
		return jump_at_moving_suite;
	}

	////////////////////////////////////////////////////////////////////////////////
	void
	split_pose( pose::Pose & pose, utility::vector1< Size > const & moving_res_list ){
		for ( Size n = 1; n <= moving_res_list.size(); n++ ){
			Size const moving_res = moving_res_list[n];
			Size const reference_res = pose.fold_tree().get_parent_residue( moving_res );
			if ( reference_res == 0 ) continue; // happens when root atom is 'moving', e.g. in moves from_scratch.
			split_pose( pose, moving_res, reference_res );
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	void
	fix_protein_jump_atom( pose::Pose & pose, Size const res, std::string const atom_name ){
		using namespace core::kinematics;
		if ( !pose.residue_type( res ).is_protein() ) return;
		FoldTree f = pose.fold_tree();
		for ( Size n = 1; n <= f.num_jump(); n++ ){
			if ( f.upstream_jump_residue( n ) == static_cast<int>( res ) ){
				f.set_jump_atoms( n, atom_name, f.downstream_atom( n ), false /*intra_residue_stubs*/ );
			}
			if ( f.downstream_jump_residue( n ) == static_cast<int>( res ) ){
				f.set_jump_atoms( n, f.upstream_atom( n ), atom_name, false /*intra_residue_stubs*/ );
			}
		}
		pose.fold_tree( f );
	}


	////////////////////////////////////////////////////////////////////////////////
	void
	align_pose_and_add_rmsd_constraints( pose::Pose & pose,
																			 pose::PoseCOP native_pose,
																			 utility::vector1< Size > const & moving_res_list,
																			 Real const rmsd_screen ) {

	if ( native_pose == 0 ) return;

	utility::vector1< Size > root_partition_res = figure_out_root_partition_res( pose, moving_res_list );
	if ( root_partition_res.size() == 0 ) root_partition_res.push_back( pose.fold_tree().root() );

	// can later generalize to use 'reference_pose', not necessarily native_pose.
	sampling::align::StepWisePoseAligner pose_aligner( *native_pose );
	pose_aligner.set_root_partition_res( root_partition_res );
	Pose pose_save = pose;
	pose_aligner.apply( pose );
	TR.Debug << "SUPERIMPOSE RMSD: " <<  pose_aligner.rmsd_over_alignment_atoms() << std::endl;
	if ( pose_aligner.rmsd_over_alignment_atoms() < 1.0e-5 ) pose = pose_save; // to avoid floating point deviations.
	if ( rmsd_screen > 0.0 ) pose_aligner.create_coordinate_constraints( pose, rmsd_screen );
}

///////////////////////////////////////////////////////////////////////////////
void
add_to_pose_list( utility::vector1< core::pose::PoseOP > & pose_list, pose::Pose const & pose, std::string const pose_tag ) {
	core::pose::PoseOP pose_op = pose.clone();
	tag_into_pose( *pose_op, pose_tag );
	pose_list.push_back( pose_op );
}

//////////////////////////////////////////////////////////////////////////////
// will deprecate this soon.
bool
is_protein( pose::Pose const & pose, utility::vector1< Size > const & moving_res_list ) {
	Size const example_res = ( moving_res_list.size() > 0 ) ? moving_res_list[1] : 1;
	if ( pose.residue_type( example_res ).is_protein() ){
		return true;
	} else {
		runtime_assert( pose.residue_type( example_res ).is_RNA() );
		return false;
	}
}


/////////////////////////////////////////////////////////
core::scoring::ScoreFunctionOP
get_minimize_scorefxn( core::pose::Pose const & pose,
											 core::scoring::ScoreFunctionCOP scorefxn,
											 StepWiseModelerOptionsCOP options ){
	using namespace core::scoring;
	ScoreFunctionOP minimize_scorefxn = scorefxn->clone();
	if (minimize_scorefxn->get_weight( atom_pair_constraint ) == 0.0) minimize_scorefxn->set_weight( atom_pair_constraint, 1.0 ); //go ahead and turn these on
	if (minimize_scorefxn->get_weight( coordinate_constraint) == 0.0) minimize_scorefxn->set_weight( coordinate_constraint, 1.0 ); // go ahead and turn these on
	check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose, minimize_scorefxn );
	//	minimize_scorefxn->set_weight( linear_chainbreak, 150.0 ); // original SWA protein value.
	minimize_scorefxn->set_weight( linear_chainbreak, 5.0 ); // unify with RNA.
	if ( options->cart_min() && ( minimize_scorefxn->get_weight( cart_bonded ) == 0.0 ) ) minimize_scorefxn->set_weight( cart_bonded, 1.0 );
	if ( options->mapfile_activated() && minimize_scorefxn->get_weight( elec_dens_atomwise ) == 0.0 ) minimize_scorefxn->set_weight( elec_dens_atomwise, 10.0 );
	return minimize_scorefxn;
}


///////////////////////////////////////////////////////////////////////////////////
// suites that connect different domains, e.g., different input structures.
utility::vector1< core::Size >
get_domain_boundary_suites( pose::Pose const & pose ){
	utility::vector1< Size > domain_boundary_suites;
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & res_list = full_model_info.res_list();
	utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();
	utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();
	for ( Size n = 1; n < pose.total_residue(); n++ ){
		if ( !cutpoint_open_in_full_model.has_value( res_list[ n ] ) &&
				 !pose.fold_tree().is_cutpoint( n ) &&
				( res_list[ n + 1 ]  == res_list[ n ] + 1 ) &&
				 ( fixed_domain_map[ res_list[ n + 1 ] ] !=  fixed_domain_map[ res_list[ n ] ] ) &&
				 !domain_boundary_suites.has_value( n ) ){
			TR.Debug << "ADDING NEW SUITE TO BE MINIMIZED BASED ON LOCATION AT DOMAIN BOUNDARY: " << n << std::endl;
			domain_boundary_suites.push_back( n );
		}
	}
	return domain_boundary_suites;
}

/////////////////////////////////////////////////////////////////////////////////
// suites that connect different domains, e.g., different input structures. This
// function returns downstream residues within each suite connection.
utility::vector1< core::Size >
get_domain_boundary_res( pose::Pose const & pose ){

	// convert from domain boundaries suites to domain boundary residues.
	utility::vector1< core::Size > domain_boundary_suites =	get_domain_boundary_suites( pose );

	utility::vector1< core::Size > domain_boundary_res;
	for ( Size n = 1; n <= domain_boundary_suites.size(); n++ ){
		Size const & i = domain_boundary_suites[ n ];
		if ( pose.fold_tree().get_parent_residue( i ) == int(i+1) ){
			domain_boundary_res.push_back( i );
		} else {
			runtime_assert( pose.fold_tree().get_parent_residue( i+1 ) == int(i) );
			domain_boundary_res.push_back( i+1 );
		}
	}
	return domain_boundary_res;
}

////////////////////////////////////////////////////////////////////////////////
utility::vector1< core::Size >
get_moving_res_including_domain_boundaries( pose::Pose const & pose, utility::vector1< core::Size > const & moving_res_list ){
	return merge_vectors( moving_res_list, get_domain_boundary_res( pose ) );
}

////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_all_working_moving_res( working_parameters::StepWiseWorkingParametersCOP working_parameters) {
	return merge_vectors( working_parameters->working_moving_res_list(),
												working_parameters->working_bridge_res() );
}

////////////////////////////////////////////////////////////////////////////////
void
virtualize_side_chains( pose::Pose & pose ) {
	for ( Size n = 1; n <= pose.total_residue(); n++ ){
		if ( pose.residue_type( n ).is_RNA() ){
			add_variant_type_to_pose_residue( pose, "VIRTUAL_O2PRIME_HYDROGEN", n );
		} else if ( pose.residue_type( n ).is_protein() ){
			add_variant_type_to_pose_residue( pose, "VIRTUAL_SIDE_CHAIN", n );
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
remove_silent_file_if_it_exists( std::string const & silent_file){
	if ( utility::file::file_exists( silent_file ) ) {
		TR << "WARNING: silent_file " << silent_file << " already exists! removing..." << std::endl;
		runtime_assert( std::remove( silent_file.c_str() ) );
	}
}

////////////////////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionCOP
initialize_sample_scorefxn( core::scoring::ScoreFunctionCOP scorefxn,
														pose::Pose const & pose,
														StepWiseModelerOptionsCOP modeler_options ){

	using namespace core::scoring;

	ScoreFunctionOP sample_scorefxn;
	std::string const & sample_weights = modeler_options->pack_weights();

	if ( sample_weights.size() == 0 ) {
		sample_scorefxn = scorefxn->clone();
		/////////////////////////////////////////////////////
		sample_scorefxn->set_weight( fa_rep, 0.12 ); // from RNA.
		sample_scorefxn->set_weight( linear_chainbreak, 0.0 ); // from RNA.
		sample_scorefxn->set_weight( chainbreak, 0.0 ); // from RNA.
		/////////////////////////////////////////////////////
	} else { // this used to happen for proteins -- may decide to deprecate.
		sample_scorefxn = ScoreFunctionFactory::create_score_function( sample_weights );
		check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose, sample_scorefxn );
		sample_scorefxn->set_weight( linear_chainbreak, 0.2 /*arbitrary*/ ); // will cause problem with RNA?!
		if ( modeler_options->mapfile_activated()  && sample_scorefxn->get_weight( elec_dens_atomwise ) == 0.0 ){
			sample_scorefxn->set_weight( elec_dens_atomwise, 10.0 );
		}
	}

	return sample_scorefxn;
}

////////////////////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionCOP
initialize_pack_scorefxn( core::scoring::ScoreFunctionCOP sample_scorefxn,
													pose::Pose const & /* put back in if we want to try contains_protein hack*/ ){

	using namespace core::scoring;

	//	if ( !contains_protein( pose ) ) return initialize_o2prime_pack_scorefxn( sample_scorefxn );

	ScoreFunctionOP pack_scorefxn = sample_scorefxn->clone();

	// hack for speed -- geom_sol & lk_nonpolar are too slow right now.
	// [see also: O2PrimePacker]
	if ( sample_scorefxn->has_nonzero_weight( geom_sol ) ||
			 sample_scorefxn->has_nonzero_weight( geom_sol_fast ) ){

		runtime_assert( sample_scorefxn->has_nonzero_weight( lk_nonpolar ) );
		Real const lk_weight = sample_scorefxn->get_weight( lk_nonpolar );

		pack_scorefxn->set_weight( geom_sol,      0.0 );
		pack_scorefxn->set_weight( geom_sol_fast, 0.0 );
		pack_scorefxn->set_weight( lk_nonpolar,   0.0 );
		pack_scorefxn->set_weight( fa_sol,        lk_weight );

	} else {
		runtime_assert( !pack_scorefxn->has_nonzero_weight( lk_nonpolar ) ); // only allow lk_nonpolar if with geom_sol.
	}

	// these also take too long to compute.
	pack_scorefxn->set_weight( fa_stack, 0.0 );
	pack_scorefxn->set_weight( ch_bond, 0.0 );

	return pack_scorefxn;
}

////////////////////////////////////////////////////////////////////////
// may be deprecated in favor of initialize_pack_scorefxn above.
core::scoring::ScoreFunctionCOP
initialize_o2prime_pack_scorefxn( core::scoring::ScoreFunctionCOP const & scorefxn ){

	using namespace core::scoring;

	ScoreFunctionOP o2prime_pack_scorefxn = new ScoreFunction;

	// Each of the following terms have been pretty optimized for the packer (trie, etc.)
	o2prime_pack_scorefxn->set_weight( fa_atr, scorefxn->get_weight( fa_atr ) );
	o2prime_pack_scorefxn->set_weight( fa_rep, scorefxn->get_weight( fa_rep ) );
	o2prime_pack_scorefxn->set_weight( hbond_lr_bb_sc, scorefxn->get_weight( hbond_lr_bb_sc ) );
	o2prime_pack_scorefxn->set_weight( hbond_sr_bb_sc, scorefxn->get_weight( hbond_sr_bb_sc ) );
	o2prime_pack_scorefxn->set_weight( hbond_sc, scorefxn->get_weight( hbond_sc ) );
	o2prime_pack_scorefxn->set_weight( free_2HOprime, scorefxn->get_weight( free_2HOprime ) );
	///Warning, don't include hbond_intra, since hbond_intra HAS NOT been been optimized for packing!


	if ( scorefxn->has_nonzero_weight( lk_nonpolar ) ) {
		//// note that geom_sol is not optimized well --> replace with lk_sol for now.
		o2prime_pack_scorefxn->set_weight( fa_sol, scorefxn->get_weight( lk_nonpolar ) );
	} else {
		runtime_assert( scorefxn->has_nonzero_weight( fa_sol ) );
		o2prime_pack_scorefxn->set_weight( fa_sol, scorefxn->get_weight( fa_sol ) );
	}
	//This sets NO_HB_ENV_DEP, INCLUDE_INTRA_RES_RNA_HB and etcs.
	o2prime_pack_scorefxn->set_energy_method_options( scorefxn->energy_method_options() );

	return o2prime_pack_scorefxn;
}

//////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_all_residues( pose::Pose const & pose ){
	utility::vector1< Size > all_res;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) all_res.push_back( n );
	return all_res;
}

//////////////////////////////////////////////////////////////////////
Size
find_downstream_connection_res( pose::Pose const & pose,
																utility::vector1< Size > const & moving_partition_res ){
	Size downstream_connection_res( 0 );
	for ( Size n = 1; n <= moving_partition_res.size(); n++ ) {
		Size const parent_res = pose.fold_tree().get_parent_residue( moving_partition_res[n] );
		if ( parent_res > 0 && !moving_partition_res.has_value( parent_res ) ){
			runtime_assert( downstream_connection_res == 0 );
			downstream_connection_res = moving_partition_res[n];
		}
	}
	return downstream_connection_res;
}


//////////////////////////////////////////////////////////////////////
Size
get_unique_connection_res( pose::Pose const & pose,
													 utility::vector1< Size > const & moving_partition_res ) {
	Size moving_res1 = find_downstream_connection_res( pose, moving_partition_res );

	utility::vector1< Size > const other_partition_res = get_other_residues( moving_partition_res, pose.total_residue() );
	Size moving_res2 = find_downstream_connection_res( pose, other_partition_res );

	// can't have 'downstream connections' both ways -- something wrong with partition.
	runtime_assert( (moving_res1 == 0) || ( moving_res2 == 0 ) );
	if ( moving_res1 > 0 ) return moving_res1;
	if ( moving_res2 > 0 ) return moving_res2;
	return pose.fold_tree().root();
}


} //sampling
} //stepwise
} //protocols
