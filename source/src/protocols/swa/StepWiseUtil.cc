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
#include <protocols/swa/StepWiseUtil.hh>
#include <protocols/swa/protein/StepWiseProteinUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>

#include <protocols/rna/RNA_ProtocolUtil.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/SequenceMapping.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/izstream.hh>

#include <iostream>
#include <fstream>
#include <sstream>
#include <ObjexxFCL/format.hh>
#include <set>


//Auto Headers
#include <core/kinematics/MoveMap.hh>
#include <utility/vector1.hh>
using namespace core;
using numeric::conversions::degrees;
using numeric::conversions::radians;

static basic::Tracer TR( "protocols.swa.StepWiseUtil" );

namespace protocols {
namespace swa {

	//////////////////////////////////////////////////////////////////////////
	Size
	make_cut_at_moving_suite( pose::Pose & pose, Size const & moving_suite ){

		core::kinematics::FoldTree f( pose.fold_tree() );
		//		Size jump_at_moving_suite( 0 );
		f.new_jump( moving_suite, moving_suite+1, moving_suite );
		pose.fold_tree( f );

		int const i( moving_suite ), j( moving_suite+1 );
		for ( Size n = 1; n <= f.num_jump(); n++ ) {
			if ( f.upstream_jump_residue(n) == i && f.downstream_jump_residue(n) == j ) return n;
			if ( f.upstream_jump_residue(n) == j && f.downstream_jump_residue(n) == i ) return n;
		}

		utility_exit_with_message( "Problem with jump number" );

		return 0; // we never get here.
	}

	///////////////////////////////////////////////////////////////////////////////////////
	bool
	Is_close_chain_break(pose::Pose const & pose){

		for(Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++) {

			if ( !pose.residue( seq_num  ).has_variant_type( chemical::CUTPOINT_LOWER )  ) continue;
			if ( !pose.residue( seq_num+1 ).has_variant_type( chemical::CUTPOINT_UPPER )  ) continue;

			return true;
		}
		return false;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	Output_boolean(std::string const & tag, bool boolean, std::ostream & TR ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;
		TR << tag;

		if(boolean==true){
			TR << A(4,"T");
		} else {
			TR << A(4,"F");
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	Output_boolean(bool boolean, std::ostream & TR ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;

		if(boolean==true){
			TR << A(4,"T");
		} else {
			TR << A(4,"F");
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	Output_movemap(kinematics::MoveMap const & mm, Size const total_residue, std::ostream & TR){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;
		using namespace core::kinematics;
		using namespace core::id;
		Size spacing=10;

		std::cout << "Movemap (in term of partial_pose seq_num): " << std::endl;
		std::cout << A(spacing,"res_num") << A(spacing,"alpha") << A(spacing,"beta") << A(8,"gamma") << A(8,"delta") <<A(8,"eplison") <<A(8,"zeta");
		std::cout << A(spacing,"chi_1") << A(spacing,"nu_2") << A(spacing,"nu_1") << A(8,"chi_O2") << std::endl;

		for(Size n=1; n<= total_residue; n++){

			std::cout << I(spacing, 3 , n);
			Output_boolean(mm.get(TorsionID( n , id::BB,  1 )), TR ); A(spacing-4, "");
			Output_boolean(mm.get(TorsionID( n , id::BB,  2 )), TR ); A(spacing-4, "");
			Output_boolean(mm.get(TorsionID( n , id::BB,  3 )), TR ); A(spacing-4, "");
			Output_boolean(mm.get(TorsionID( n , id::BB,  4 )), TR ); A(spacing-4, "");
			Output_boolean(mm.get(TorsionID( n , id::BB,  5 )), TR ); A(spacing-4, "");
			Output_boolean(mm.get(TorsionID( n , id::BB,  6 )), TR ); A(spacing-4, "");
			Output_boolean(mm.get(TorsionID( n , id::CHI, 1 )), TR ); A(spacing-4, "");
			Output_boolean(mm.get(TorsionID( n , id::CHI, 2 )), TR ); A(spacing-4, "");
			Output_boolean(mm.get(TorsionID( n , id::CHI, 3 )), TR ); A(spacing-4, "");
			Output_boolean(mm.get(TorsionID( n , id::CHI, 4 )), TR ); A(spacing-4, "");
			std::cout << std::endl;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////
	// This is similar to code in RNA_Minimizer.cc
	void
	Figure_out_moving_residues( core::kinematics::MoveMap & mm, core::pose::Pose const & pose,
															utility::vector1< core::Size > const & fixed_res,
															bool const move_takeoff_torsions,
															bool const move_jumps_between_chains
															)
  {

		using namespace core::id;

		Size const nres( pose.total_residue() );

		ObjexxFCL::FArray1D< bool > allow_insert( nres, true );
		for (Size i = 1; i <= fixed_res.size(); i++ ) allow_insert( fixed_res[ i ] ) = false;
		for (Size n = 1; n <= pose.total_residue(); n++ ){
			if ( pose.residue( n ).has_variant_type( "VIRTUAL_RESIDUE" ) ) allow_insert( n ) = false;
		}

		mm.set_bb( false );
		mm.set_chi( false );
		mm.set_jump( false );

		std::cout << "ALLOWING BB TO MOVE: ";
		for  (Size i = 1; i <= nres; i++ )  {

			//std::cout << "ALLOW INSERT " << i << " " << allow_insert(i) << std::endl;
			if ( !move_takeoff_torsions && !allow_insert(i) ) continue; // don't allow, e.g., psi/omega of residue before loop to move.

			utility::vector1< TorsionID > torsion_ids;

			for ( Size torsion_number = 1; torsion_number <= pose.residue( i ).mainchain_torsions().size(); torsion_number++ ) {
				torsion_ids.push_back( TorsionID( i, id::BB, torsion_number ) );
			}
			for ( Size torsion_number = 1; torsion_number <= pose.residue( i ).nchi(); torsion_number++ ) {
				torsion_ids.push_back( TorsionID( i, id::CHI, torsion_number ) );
			}


			for ( Size n = 1; n <= torsion_ids.size(); n++ ) {

				TorsionID const & torsion_id  = torsion_ids[ n ];

				id::AtomID id1,id2,id3,id4;
				bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
				if (fail) continue;

				// If there's any atom that is in a moving residue by this torsion, let the torsion move.
				//  should we handle a special case for cutpoint atoms? I kind of want those all to move.
				if ( !allow_insert( id1.rsd() ) && !allow_insert( id2.rsd() ) && !allow_insert( id3.rsd() )  && !allow_insert( id4.rsd() ) ) continue;
				mm.set(  torsion_id, true );

				if ( n == 1 )	std::cout << ' ' <<  i;

			}

		}
		std::cout << std::endl;

		utility::vector1< Size > chain_index;
		Size chain_number( 0 );
		for (Size n = 1; n <= pose.total_residue(); n++ ){
			if ( pose.residue_type( n ).is_lower_terminus()  && !pose.residue_type( n ).has_variant_type( chemical::N_ACETYLATION ) ) chain_number++;
			chain_index.push_back( chain_number );
		}

		//		for (Size n = 1; n <= pose.total_residue(); n++ ) std::cout << "CHAIN "<< n << ' ' << chain_index[ n ] << std::endl;

		for (Size n = 1; n <= pose.fold_tree().num_jump(); n++ ){
			Size const jump_pos1( pose.fold_tree().upstream_jump_residue( n ) );
			Size const jump_pos2( pose.fold_tree().downstream_jump_residue( n ) );

			if ( allow_insert( jump_pos1 ) || allow_insert( jump_pos2 ) ) 	 {
				mm.set_jump( n, true );
				std::cout << "allow_insert ALLOWING JUMP " << n << " to move. It connects " << jump_pos1 << " and " << jump_pos2 << "." << std::endl;
			}

			if ( move_jumps_between_chains ){
				if ( chain_index[ jump_pos1 ] != chain_index[ jump_pos2 ] ){
					std::cout << "move_jumps_between_chains ALLOWING JUMP " << n << " to move. It connects " << jump_pos1 << " and " << jump_pos2 << "." << std::endl;
				}
			}

		}

	}

	///////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////
	void
	pdbslice( core::pose::Pose & new_pose,
						core::pose::Pose const & pose,
						utility::vector1< core::Size > const & slice_res )
	{
		using namespace core::pose;
		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::scoring::rna;

		new_pose.clear();

		for ( Size i = 1; i <= slice_res.size(); i++ ) {

			ResidueOP residue_to_add = pose.residue( slice_res[ i ] ).clone() ;

			if ( (i > 1 &&  ( slice_res[i] != slice_res[i-1] + 1 )) /*new segment*/ || residue_to_add->is_lower_terminus() || residue_to_add->has_variant_type( "N_ACETYLATION") || (i>1 && pose.fold_tree().is_cutpoint( slice_res[i-1] ) )  ){
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

		// //how about a new fold tree?
		// core::kinematics::FoldTree f( new_pose.total_residue() );
		// for ( Size i = 1; i < slice_res.size(); i++ ) {
		// 	if ( slice_res[i+1] > (slice_res[i]+1) ){
		// 		f.new_jump( i, i+1, i);
		// 	}
		// }
		// new_pose.fold_tree( f );

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
 	id::AtomID_Map< id::AtomID >
 	create_alignment_id_map(	pose::Pose const & mod_pose,
													pose::Pose const & ref_pose,
													utility::vector1< core::Size > const & superimpose_res ){

		std::map< core::Size, core::Size > res_map;

		// 		if( ref_pose.sequence()!=mod_pose.sequence() ){
		// 			utility_exit_with_message( "ref_pose.sequence()!=mod_pose.sequence()");
		// 		}

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
 		using namespace protocols::swa::protein;
 		using namespace protocols::swa::rna;
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


	////////////////////////////////////////////////////////
	// Move ths out of stepwise_protein_test.cc!!
	////////////////////////////////////////////////////////
	/////////////////////////////////////////////////
	utility::vector1< std::string > load_s_and_l()
	{
		using basic::options::option;
		using utility::vector1;
		using namespace basic::options::OptionKeys;

		// concatenate -s and -l flags together to get total list of PDB files
		vector1< std::string > pdb_file_names;
		if ( option[ in::file::s ].active() ) {
			pdb_file_names = option[ in::file::s ]().vector(); // make a copy (-s)
		}

		vector1< std::string > list_file_names;
		if ( option[ in::file::l ].active() )
			list_file_names = option[ in::file::l ]().vector(); // make a copy (-l)
		if ( option[ in::file::list ].active() ){
			vector1< std::string > better_list_file_names;
			better_list_file_names= option[in::file::list ]().vector(); // make a copy (-list)
			for(vector1< std::string >::iterator i = better_list_file_names.begin(), i_end = better_list_file_names.end(); i != i_end; ++i) {
				list_file_names.push_back(*i); // make a copy (-l)
			}
		}

		for(vector1< std::string >::iterator i = list_file_names.begin(), i_end = list_file_names.end(); i != i_end; ++i) {
			std::string filename( *i );
			utility::io::izstream data( filename.c_str() );
			if ( !data.good() ) {
				utility_exit_with_message( "Unable to open file: " + filename + '\n' );
			}
			std::string line;
			while( getline(data, line) ) {
				pdb_file_names.push_back( std::string(line) );
			}
			data.close();
		}

		return pdb_file_names;
	}


	///////////////////////////////////////////////////////////////////////
	std::string
	get_file_name( std::string const & silent_file, std::string const & tag )
	{
		Size pos( silent_file.find( ".out" ) );
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

		// Oh come on, testing basic matrix algebra here.

		// std::cout << "M1:  " << std::endl;
		// std::cout << M1(1,1) <<  ' ' << M1(1,2)  << ' ' << M1(1,3) << std::endl;
		// std::cout << M1(2,1) <<  ' ' << M1(2,2)  << ' ' << M1(2,3) << std::endl;
		// std::cout << M1(3,1) <<  ' ' << M1(3,2)  << ' ' << M1(3,3) << std::endl;
		// std::cout << "x: " << xaxis1(1) << " " << xaxis1(2) << " " << xaxis1(3) << std::endl;
		// std::cout << "y: " << yaxis1(1) << " " << yaxis1(2) << " " << yaxis1(3) << std::endl;
		// std::cout << "z: " << zaxis1(1) << " " << zaxis1(2) << " " << zaxis1(3) << std::endl;

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
			//			std::cout << "ROTATE " << m << " " << i << " " << pose.xyz( AtomID(m,i) )(1)  << " " << ref_pose.xyz( AtomID(m,i) )(1) << std::endl;
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
		using namespace protocols::swa;

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
	void
	add_to_atom_id_map_after_checks( std::map< id::AtomID, id::AtomID> & atom_id_map,
																	 std::string const & atom_name,
																	 Size const & n1, Size const & n2,
																	 pose::Pose const & pose1, pose::Pose const & pose2 ){

		using namespace core::id;

		runtime_assert ( n1 >= 1 && n1 <= pose1.total_residue() );
		runtime_assert ( n2 >= 1 && n2 <= pose2.total_residue() );
		runtime_assert( pose1.residue_type( n1 ).aa() == pose2.residue_type( n2 ).aa() );

		if ( ! pose1.residue_type( n1 ).has( atom_name ) ) return;
		if ( ! pose2.residue_type( n2 ).has( atom_name ) ) return;

		Size const idx1 = pose1.residue_type( n1 ).atom_index( atom_name );
		Size const idx2 = pose2.residue_type( n2 ).atom_index( atom_name );

		if ( pose1.residue_type( n1 ).is_virtual( idx1 ) ) return;
		if ( pose2.residue_type( n2 ).is_virtual( idx2 ) ) return;

		atom_id_map[ AtomID( idx1, n1 ) ] = AtomID( idx2, n2 );

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Real
	get_all_atom_rmsd( pose::Pose const & pose, pose::Pose const & native_pose, utility::vector1< Size > const & rmsd_res ){

		using namespace core::chemical;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace protocols::rna;

		Pose native_pose_local = native_pose; // local working copy, mutated in cases where nucleotides have been designed ('n')

		// first need to slice up native_pose to match residues in actual pose.
		// define atoms over which to compute RMSD, using rmsd_res.
		FullModelInfo const & full_model_info = const_full_model_info_from_pose( pose );
		utility::vector1< Size > const sub_to_full = full_model_info.sub_to_full();
		std::string const full_sequence = full_model_info.full_sequence();

		utility::vector1< Size > working_rmsd_res;
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( rmsd_res.has_value( sub_to_full[ n ] ) ) {
				working_rmsd_res.push_back( n );

				char const pose_nt = pose.sequence()[ n-1 ];
				if ( full_sequence[ sub_to_full[ n ] - 1 ] == 'n' ){
					mutate_position( native_pose_local, sub_to_full[ n ], pose_nt );
				} else {
					runtime_assert( full_sequence[ sub_to_full[ n ] - 1 ] == pose_nt);
				}
				runtime_assert( native_pose_local.sequence()[ sub_to_full[ n ] - 1] == pose_nt );
			}
		}

		std::map< AtomID, AtomID > atom_id_map;
		for ( Size k = 1; k <= working_rmsd_res.size(); k++ ){

			Size const n = working_rmsd_res[ k ];

			for ( Size q = 1; q <= pose.residue_type( n ).nheavyatoms(); q++ ){
				add_to_atom_id_map_after_checks( atom_id_map,
																				 pose.residue_type( n ).atom_name( q ),
																				 n, sub_to_full[ n ],
																				 pose, native_pose_local );
			}

			if ( ! working_rmsd_res.has_value( n + 1 ) &&
					 ( n + 1 ) <= pose.total_residue() &&
					 ( !pose.fold_tree().is_cutpoint( n ) || pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) ) ) {
				// RNA-specific. Would be trivial to expand to proteins.
				runtime_assert( pose.residue_type( n + 1 ).is_RNA() );
				add_to_atom_id_map_after_checks( atom_id_map, " P  ", n + 1, sub_to_full[ n + 1 ], pose, native_pose_local );
				add_to_atom_id_map_after_checks( atom_id_map, " OP1", n + 1, sub_to_full[ n + 1 ], pose, native_pose_local );
				add_to_atom_id_map_after_checks( atom_id_map, " OP2", n + 1, sub_to_full[ n + 1 ], pose, native_pose_local );
				add_to_atom_id_map_after_checks( atom_id_map, " O5'", n + 1, sub_to_full[ n + 1 ], pose, native_pose_local );
			}
		}

		TR << "Calculating RMSD over " << atom_id_map.size() << " atoms." << std::endl;
		Real rmsd( 0.0 );
		if ( atom_id_map.size() > 0 ) rmsd = scoring::rms_at_all_corresponding_atoms( pose, native_pose, atom_id_map );

		return rmsd;

	}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	get_number_missing_residues( pose::Pose const & pose ) {

		using namespace core::pose::full_model_info;

		FullModelInfo const & full_model_info = const_full_model_info_from_pose( pose );
		utility::vector1< Size > const sub_to_full = full_model_info.sub_to_full();
		std::string const full_sequence = full_model_info.full_sequence();

		return ( full_sequence.size() - sub_to_full.size() );

	}


}
}
