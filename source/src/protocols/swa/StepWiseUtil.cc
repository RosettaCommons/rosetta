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
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/SequenceMapping.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/izstream.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/scoring/Energies.hh>

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
		Size num_jump( 0 ), jump_idx( 0 );
		for ( Size n = 1; n <= fold_tree.num_jump(); n++ ) {
			if ( fold_tree.upstream_jump_residue( n ) == i || fold_tree.downstream_jump_residue( n ) == i ) {
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


	//////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	output_boolean(std::string const & tag, bool boolean, std::ostream & TR ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;
		TR << tag;

		if(boolean==true){
			TR << A(4,"T");
		} else {
			TR << A(4,"F");
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	output_boolean(bool boolean, std::ostream & TR ){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;

		if(boolean==true){
			TR << A(4,"T");
		} else {
			TR << A(4,"F");
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	output_movemap(kinematics::MoveMap const & mm, Size const total_residue, std::ostream & TR){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::format;
		using namespace core::kinematics;
		using namespace core::id;
		Size spacing=10;

		std::cout << "Movemap (in term of partial_pose seq_num): " << std::endl;
		std::cout << A(spacing,"res_num") << A(spacing,"alpha") << A(spacing,"beta") << A(8,"gamma") << A(8,"delta") <<A(8,"eplison") <<A(8,"zeta");
		std::cout << A(spacing,"chi_1") << A(spacing,"nu_2") << A(spacing,"nu_1") << A(8,"chi_O2") << std::endl;

		for(Size n=1; n<= total_residue; n++){

			std::cout << I(spacing, 3 , n);
			output_boolean(mm.get(TorsionID( n , id::BB,  1 )), TR ); A(spacing-4, "");
			output_boolean(mm.get(TorsionID( n , id::BB,  2 )), TR ); A(spacing-4, "");
			output_boolean(mm.get(TorsionID( n , id::BB,  3 )), TR ); A(spacing-4, "");
			output_boolean(mm.get(TorsionID( n , id::BB,  4 )), TR ); A(spacing-4, "");
			output_boolean(mm.get(TorsionID( n , id::BB,  5 )), TR ); A(spacing-4, "");
			output_boolean(mm.get(TorsionID( n , id::BB,  6 )), TR ); A(spacing-4, "");
			output_boolean(mm.get(TorsionID( n , id::CHI, 1 )), TR ); A(spacing-4, "");
			output_boolean(mm.get(TorsionID( n , id::CHI, 2 )), TR ); A(spacing-4, "");
			output_boolean(mm.get(TorsionID( n , id::CHI, 3 )), TR ); A(spacing-4, "");
			output_boolean(mm.get(TorsionID( n , id::CHI, 4 )), TR ); A(spacing-4, "");
			std::cout << std::endl;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////
	// This is similar to code in RNA_Minimizer.cc
	void
	figure_out_moving_residues( core::kinematics::MoveMap & mm, core::pose::Pose const & pose,
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

		for  (Size i = 1; i <= nres; i++ )  {

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

			}

		}

		utility::vector1< Size > chain_index;
		Size chain_number( 0 );
		for (Size n = 1; n <= pose.total_residue(); n++ ){
			if ( pose.residue_type( n ).is_lower_terminus()  && !pose.residue_type( n ).has_variant_type( chemical::N_ACETYLATION ) ) chain_number++;
			chain_index.push_back( chain_number );
		}

		for (Size n = 1; n <= pose.fold_tree().num_jump(); n++ ){
			Size const jump_pos1( pose.fold_tree().upstream_jump_residue( n ) );
			Size const jump_pos2( pose.fold_tree().downstream_jump_residue( n ) );

			if ( allow_insert( jump_pos1 ) || allow_insert( jump_pos2 ) ) 	 {
				mm.set_jump( n, true );
				//				std::cout << "allow_insert ALLOWING JUMP " << n << " to move. It connects " << jump_pos1 << " and " << jump_pos2 << "." << std::endl;
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
		using namespace core::chemical::rna;

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
	bool
	residue_is_bulged( pose::Pose const & pose, Size const & resid ) {
		boost::unordered_map < core::Size , core::Size > num_stacks = pose.get_stacking_map();

		if ( num_stacks[ resid ] < 2 ) {
			return true;
		}

		return false;
	}

	//////////////////
	void superimpose_at_fixed_res( pose::Pose & pose, pose::Pose const & native_pose, Real & rmsd, Size & natoms_rmsd, core::pose::full_model_info::FullModelInfo const & full_model_info, bool skip_bulges = false ) {

		using namespace core::chemical;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace core::scoring;
		using namespace protocols::rna;

		Pose native_pose_local = native_pose; // local working copy, mutated in cases where nucleotides have been designed ('n')

		// first need to slice up native_pose to match residues in actual pose.
		// define atoms over which to compute RMSD, using rmsd_res.
		utility::vector1< Size > const & res_list = full_model_info.res_list();
		utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();
		std::string const full_sequence = full_model_info.full_sequence();

		// following needs to be updated.
		utility::vector1< Size > const rmsd_res = full_model_info.moving_res_in_full_model();

		utility::vector1< Size > calc_rms_res;
		utility::vector1< Size > skipped_residues;

		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( rmsd_res.has_value( res_list[ n ] ) ) {
				if ( skip_bulges && residue_is_bulged( pose, n ) && residue_is_bulged( native_pose_local, res_list[ n ] ) ) {
					skipped_residues.push_back( n );
					continue;
				}

				calc_rms_res.push_back( n );

				char const pose_nt = pose.sequence()[ n-1 ];
				if ( full_sequence[ res_list[ n ] - 1 ] == 'n' ){
					mutate_position( native_pose_local, res_list[ n ], pose_nt );
				} else {
					runtime_assert( full_sequence[ res_list[ n ] - 1 ] == pose_nt);
				}
				runtime_assert( native_pose_local.sequence()[ res_list[ n ] - 1] == pose_nt );
			}
		}

		std::map< AtomID, AtomID > calc_rms_atom_id_map;

		for ( Size k = 1; k <= calc_rms_res.size(); k++ ){
			Size const n = calc_rms_res[ k ];
			for ( Size q = 1; q <= pose.residue_type( n ).nheavyatoms(); q++ ){
				add_to_atom_id_map_after_checks( calc_rms_atom_id_map,
												pose.residue_type( n ).atom_name( q ),
												n, res_list[ n ],
												pose, native_pose_local );
			}
		}

		utility::vector1< Size > calc_rms_suites;
		// additional RNA suites over which to calculate RMSD
		for ( Size n = 1; n < pose.total_residue(); n++ ){

			if ( !pose.residue_type( n ).is_RNA() || !pose.residue_type( n + 1 ).is_RNA() ) continue;
			if ( calc_rms_res.has_value( n+1 ) ) continue;

			// Atoms at ends of rebuilt loops:
			if ( calc_rms_res.has_value( n ) &&
				( !pose.fold_tree().is_cutpoint( n ) || pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) ) ) {
				calc_rms_suites.push_back( n ); continue;
			}

			// Domain boundaries:
			if ( (res_list[ n+1 ] == res_list[ n ] + 1) &&
				fixed_domain_map[ res_list[ n ] ] != 0 &&
				fixed_domain_map[ res_list[ n+1 ] ] != 0 &&
				fixed_domain_map[ res_list[ n ] ] != fixed_domain_map[ res_list[ n+1 ] ] ){
				calc_rms_suites.push_back( n );
			}
		}

		utility::vector1< std::string > const extra_suite_atoms = utility::tools::make_vector1( " P  ", " OP1", " OP2", " O5'" );
		for ( Size k = 1; k <= calc_rms_suites.size(); k++ ){
			Size const n = calc_rms_suites[ k ];
			for ( Size q = 1; q <= extra_suite_atoms.size(); q++ ){
				add_to_atom_id_map_after_checks( calc_rms_atom_id_map, extra_suite_atoms[ q ],
												n+1, res_list[ n+1 ],
												pose, native_pose_local );
			}
		}

		//		for ( std::map < AtomID, AtomID >::const_iterator it = calc_rms_atom_id_map.begin();
		//					it != calc_rms_atom_id_map.end(); it++ ){
		//			TR << it->first << " mapped to " << it->second << std::endl;
		//		}

		// define superposition atoms. Should be over atoms in any fixed domains. This should be
		// the 'inverse' of calc_rms atoms.

		std::map< AtomID, AtomID > superimpose_atom_id_map;
		for ( Size n = 1; n < pose.total_residue(); n++ ){
			if ( skipped_residues.has_value( n ) ) continue;

			for ( Size q = 1; q <= pose.residue_type( n ).nheavyatoms(); q++ ){
				if ( calc_rms_atom_id_map.find( AtomID( q, n ) ) == calc_rms_atom_id_map.end() ){
					add_to_atom_id_map_after_checks( superimpose_atom_id_map,
													pose.residue_type( n ).atom_name( q ),
													n, res_list[ n ],
													pose, native_pose_local );
				}
			}
		}

		// What if there weren't any fixed atoms? superimpose over everything.
		if ( superimpose_atom_id_map.size() == 0 ) superimpose_atom_id_map = calc_rms_atom_id_map;

		rmsd = 0.0;
		natoms_rmsd = calc_rms_atom_id_map.size();
		if ( natoms_rmsd > 0 && superimpose_atom_id_map.size() > 0 ) {
			//		Real const rmsd0 = rms_at_corresponding_atoms( pose, native_pose, atom_id_map );
			scoring::superimpose_pose( pose, native_pose, superimpose_atom_id_map );
			rmsd = rms_at_corresponding_atoms_no_super( pose, native_pose, calc_rms_atom_id_map );
		}
		TR << "Pose " << make_tag_with_dashes(res_list) << ": RMSD " << rmsd << " over " << natoms_rmsd << " atoms, superimposing on " << superimpose_atom_id_map.size() << " atoms. " << std::endl;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	superimpose_at_fixed_res( pose::Pose & pose, pose::Pose const & native_pose, Real & rmsd, Size & natoms_rmsd, core::pose::full_model_info::FullModelInfoOP full_model_pointer, bool skip_bulges = false ){

		using namespace core::chemical;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace core::scoring;
		using namespace protocols::rna;

		if ( full_model_pointer ) {
			superimpose_at_fixed_res( pose, native_pose, rmsd, natoms_rmsd, *full_model_pointer, skip_bulges );
		} else {
			FullModelInfo const & full_model_info = const_full_model_info( pose );
			superimpose_at_fixed_res( pose, native_pose, rmsd, natoms_rmsd, full_model_info, skip_bulges );
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	superimpose_recursively( pose::Pose & pose, pose::Pose const & native_pose, Real & rmsd, Size & natoms, core::pose::full_model_info::FullModelInfoOP full_model_pointer, bool skip_bulges = false ){

		using namespace core::pose;
		using namespace core::pose::full_model_info;

		Real rmsd_pose;
		Size natoms_pose;
		superimpose_at_fixed_res( pose, native_pose, rmsd_pose, natoms_pose, full_model_pointer, skip_bulges );

		Real const total_sd = ( rmsd * rmsd * natoms) + (rmsd_pose * rmsd_pose * natoms_pose );
		natoms += natoms_pose;
		if ( natoms > 0 ) {
			rmsd = std::sqrt( total_sd / Real( natoms ) );
		} else {
			runtime_assert( std::abs( rmsd ) < 1e-5 );
		}

		utility::vector1< PoseOP > const & other_pose_list = nonconst_full_model_info( pose ).other_pose_list();
		for ( Size n = 1; n <= other_pose_list.size(); n++ ){
			superimpose_recursively( *( other_pose_list[ n ] ), native_pose, rmsd, natoms, full_model_pointer, skip_bulges );
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Real
	superimpose_at_fixed_res_and_get_all_atom_rmsd( pose::Pose & pose, pose::Pose const & native_pose, core::pose::full_model_info::FullModelInfoOP full_model_pointer, bool skip_bulges ){
		Real rmsd( 0.0 );
		Size natoms( 0 );
		superimpose_recursively( pose, native_pose, rmsd, natoms, full_model_pointer, skip_bulges );
		return rmsd;
	}

	Real
	superimpose_at_fixed_res_and_get_all_atom_rmsd( pose::Pose & pose, pose::Pose const & native_pose, bool skip_bulges ){
		core::pose::full_model_info::FullModelInfoOP dummy_pointer;
		return superimpose_at_fixed_res_and_get_all_atom_rmsd( pose, native_pose, dummy_pointer, skip_bulges );
	}


	////////////////////////
	///////////////////////////////
	void
	clear_constraints_recursively( pose::Pose & pose ) {
		using namespace core::pose;
		using namespace core::pose::full_model_info;

		core::scoring::constraints::ConstraintSetOP cst_set = pose.constraint_set()->clone();
		cst_set->clear();
		pose.constraint_set( cst_set );

		utility::vector1< PoseOP > const & other_pose_list = nonconst_full_model_info( pose ).other_pose_list();
		for ( Size n = 1; n <= other_pose_list.size(); n++ ){
			clear_constraints_recursively( *( other_pose_list[ n ] ) );
		}
	}

	/////////////////////
	void
	add_coordinate_constraints_from_map( pose::Pose & pose, pose::Pose const & native_pose, std::map< id::AtomID, id::AtomID > const & superimpose_atom_id_map, core::Real const & constraint_x0, core::Real const & constraint_tol ) {

		Size const my_anchor( pose.fold_tree().root() ); //Change to use the root of the current foldtree as done by Rocco in AtomCoordinateCstMover - JAB.

		core::scoring::constraints::ConstraintSetOP cst_set = pose.constraint_set()->clone();

		for ( std::map< id::AtomID, id::AtomID >::const_iterator
			 it=superimpose_atom_id_map.begin(), it_end = superimpose_atom_id_map.end(); it != it_end; ++it ) {
			id::AtomID const mapped_atom = it->second;
			cst_set->add_constraint( new core::scoring::constraints::CoordinateConstraint ( it->first, id::AtomID(1, my_anchor), native_pose.residue(mapped_atom.rsd()).xyz(mapped_atom.atomno()), new core::scoring::func::FlatHarmonicFunc( constraint_x0, 1.0, constraint_tol )) );
			//cst_set->add_constraint( new core::scoring::constraints::CoordinateConstraint ( it->first, id::AtomID(1, my_anchor), native_pose.residue(mapped_atom.rsd()).xyz(mapped_atom.atomno()), new core::scoring::constraints::FadeFunc( -0.7, 1.5, 0.8, -1.0, 0.0)) );
		}

		pose.constraint_set( cst_set );

	}

	////////////////////////////////////
	void
	superimpose_at_fixed_res_and_add_constraints( pose::Pose & pose, pose::Pose const & native_pose, core::Real const & constraint_x0, core::Real const & constraint_tol ) {

		using namespace core::chemical;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace core::scoring;
		using namespace protocols::rna;

		Pose native_pose_local = native_pose; // local working copy, mutated in cases where nucleotides have been designed ('n')

		// first need to slice up native_pose to match residues in actual pose.
		// define atoms over which to compute RMSD, using rmsd_res.
		FullModelInfo const & full_model_info = const_full_model_info( pose );
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();
		std::string const full_sequence = full_model_info.full_sequence();

		// following needs to be updated.
		utility::vector1< Size > const rmsd_res = full_model_info.moving_res_in_full_model();

		utility::vector1< Size > calc_rms_res;
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( rmsd_res.has_value( res_list[ n ] ) ) {
				calc_rms_res.push_back( n );
				char const pose_nt = pose.sequence()[ n-1 ];
				if ( full_sequence[ res_list[ n ] - 1 ] == 'n' ){
					mutate_position( native_pose_local, res_list[ n ], pose_nt );
				} else {
					runtime_assert( full_sequence[ res_list[ n ] - 1 ] == pose_nt);
				}
				runtime_assert( native_pose_local.sequence()[ res_list[ n ] - 1] == pose_nt );
			}
		}

		std::map< AtomID, AtomID > calc_rms_atom_id_map;

		for ( Size k = 1; k <= calc_rms_res.size(); k++ ){
			Size const n = calc_rms_res[ k ];
			for ( Size q = 1; q <= pose.residue_type( n ).nheavyatoms(); q++ ){
				add_to_atom_id_map_after_checks( calc_rms_atom_id_map,
												pose.residue_type( n ).atom_name( q ),
												n, res_list[ n ],
												pose, native_pose_local );
			}
		}

		utility::vector1< Size > calc_rms_suites;
		// additional RNA suites over which to calculate RMSD
		for ( Size n = 1; n < pose.total_residue(); n++ ){

			if ( !pose.residue_type( n ).is_RNA() || !pose.residue_type( n + 1 ).is_RNA() ) continue;
			if ( calc_rms_res.has_value( n+1 ) ) continue;

			// Atoms at ends of rebuilt loops:
			if ( calc_rms_res.has_value( n ) &&
				( !pose.fold_tree().is_cutpoint( n ) || pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) ) ) {
				calc_rms_suites.push_back( n ); continue;
			}

			// Domain boundaries:
			if ( (res_list[ n+1 ] == res_list[ n ] + 1) &&
				fixed_domain_map[ res_list[ n ] ] != 0 &&
				fixed_domain_map[ res_list[ n+1 ] ] != 0 &&
				fixed_domain_map[ res_list[ n ] ] != fixed_domain_map[ res_list[ n+1 ] ] ){
				calc_rms_suites.push_back( n );
			}
		}

		utility::vector1< std::string > const extra_suite_atoms = utility::tools::make_vector1( " P  ", " OP1", " OP2", " O5'" );
		for ( Size k = 1; k <= calc_rms_suites.size(); k++ ){
			Size const n = calc_rms_suites[ k ];
			for ( Size q = 1; q <= extra_suite_atoms.size(); q++ ){
				add_to_atom_id_map_after_checks( calc_rms_atom_id_map, extra_suite_atoms[ q ],
												n+1, res_list[ n+1 ],
												pose, native_pose_local );
			}
		}

		add_coordinate_constraints_from_map( pose, native_pose_local, calc_rms_atom_id_map, constraint_x0, constraint_tol );
	}

	////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	superimpose_recursively_and_add_constraints( pose::Pose & pose, pose::Pose const & native_pose,
																							 core::Real const & constraint_x0, core::Real const & constraint_tol ) {

		using namespace core::pose;
		using namespace core::pose::full_model_info;

		superimpose_at_fixed_res_and_add_constraints( pose, native_pose, constraint_x0, constraint_tol );

		utility::vector1< PoseOP > const & other_pose_list = nonconst_full_model_info( pose ).other_pose_list();
		for ( Size n = 1; n <= other_pose_list.size(); n++ ){
			superimpose_recursively_and_add_constraints( *( other_pose_list[ n ] ), native_pose, constraint_x0, constraint_tol );
		}

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
	void
	try_reroot_at_fixed_domain( pose::Pose & pose ){

		using namespace core::pose::full_model_info;
		kinematics::FoldTree f = pose.fold_tree();
		utility::vector1< Size > const pose_domain_map = figure_out_pose_domain_map( pose );
		Size new_root( 0 );
		for ( Size n = 1; n <= f.nres(); n++ ){
			if ( pose_domain_map[ n ] > 0 &&
					 f.possible_root( n ) &&
					 ( n == 1 || f.is_cutpoint( n - 1 ) ||
						 pose_domain_map[ n - 1 ] != pose_domain_map[ n ] ) ){
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
	merge_in_other_pose( pose::Pose & pose, pose::Pose const & pose2, Size const & merge_res ){

		using namespace core::pose::datacache;
		using namespace core::pose::full_model_info;

		FullModelInfo & full_model_info = nonconst_full_model_info( pose );

		core::pose::Pose pose_scratch;
		utility::vector1< Size > const new_res_list = merge_two_poses_using_full_model_info( pose_scratch, pose, pose2, merge_res );
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
																				 Size const & merge_res ) {

		using namespace core::pose::full_model_info;
		using namespace core::pose::datacache;
		using namespace basic::datacache;

		// following assert may not be absolutely necessary.
		//		runtime_assert( & const_full_model_info( pose1 ) == & const_full_model_info( pose2 ) );

		// get working_residue information from each pose
		utility::vector1< Size > const & working_res1 = get_res_list_from_full_model_info_const( pose1 );
		utility::vector1< Size > const & working_res2 = get_res_list_from_full_model_info_const( pose2 );

		return merge_two_poses( pose, pose1, pose2, working_res1, working_res2, merge_res );

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
									 Size const & merge_res,
									 bool const fix_first_pose /* = true */) {

		using namespace core::kinematics;
		using namespace core::chemical;
		using namespace core::conformation;

		if ( working_res1.has_value( merge_res ) ){
			runtime_assert( working_res2.has_value( merge_res+1 ) );
		} else {
			runtime_assert( working_res2.has_value( merge_res   ) );
			runtime_assert( working_res1.has_value( merge_res+1 ) );
			TR.Debug << "merge_two_poses: order is switched.  " << std::endl;
			// order is switched. No problem.
			return merge_two_poses( pose, pose2, pose1, working_res2, working_res1, merge_res, ! fix_first_pose );
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

			if ( k == merge_res ) remove_upper_terminus( rsd );
			if ( k == ( merge_res + 1)  ) {
				runtime_assert( after_cutpoint );
				remove_lower_terminus( rsd );
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
		runtime_assert( jump_partners1.size() == jump_partners2.size() );

		// figure out cuts
		get_endpoints_from_pose( endpoints, pose1, working_res1 );
		get_endpoints_from_pose( endpoints, pose2, working_res2 );
		runtime_assert( endpoints.has_value( merge_res ) );
		Size const last_res = working_res[ working_res.size() ];
		runtime_assert( endpoints.has_value( last_res ) );
		for ( Size n = 1; n <= endpoints.size(); n++ ){
			if ( endpoints[ n ] != merge_res &&
					 endpoints[ n ] != last_res ) cuts.push_back( endpoints[ n ] );
		}

		Size num_cuts = cuts.size();
		runtime_assert( num_cuts == jump_partners1.size() );

		jump_partners1 = map_to_local_numbering( jump_partners1, working_res );
		jump_partners2 = map_to_local_numbering( jump_partners2, working_res );
		cuts = map_to_local_numbering( cuts, working_res );

		FoldTree f;
		ObjexxFCL::FArray2D< int > jump_point_( 2, num_cuts, 0 );
		ObjexxFCL::FArray1D< int > cuts_( num_cuts, 0 );
		for ( Size i = 1; i <= num_cuts; i++ ) {
			jump_point_( 1, i ) = jump_partners1[ i ];
			jump_point_( 2, i ) = jump_partners2[ i ];
			cuts_( i ) = cuts[ i ];
		}
		f.tree_from_jumps_and_cuts( pose.total_residue(), num_cuts, jump_point_, cuts_ );

		// fix jump atoms.
		bool const KeepStubInResidue( true );
		for ( Size i = 1; i <= num_cuts; i++ ){
			Size const n = f.jump_nr( jump_partners1[ i ], jump_partners2[ i ] );
			f.set_jump_atoms( n,
												jump_partners1[ i ], jump_atoms1[ i ],
												jump_partners2[ i ], jump_atoms2[ i ], KeepStubInResidue );
		}

		Size root( 0 );
		if ( fix_first_pose ) {
			root = find_first_root_residue( f, working_res1, working_res );
		} else {
			root = find_first_root_residue( f, working_res2, working_res );
		}
		runtime_assert( f.reorder( root ) );

		pose.fold_tree( f );

		// map (internal) coordinates from separate poses into merged one.
		std::map< Size, Size > res_map1 = get_res_map( working_res, working_res1 );
		std::map< Size, Size > res_map2 = get_res_map( working_res, working_res2 );

		copy_dofs_match_atom_names( pose, pose1, res_map1, false /*backbone_only*/, false /*ignore_virtual*/ );
		copy_dofs_match_atom_names( pose, pose2, res_map2, false /*backbone_only*/, false /*ignore_virtual*/ );

		if ( fix_first_pose ){
			superimpose_pose( pose, pose1, res_map1 );
		} else {
			superimpose_pose( pose, pose2, res_map2 );
		}

		return working_res;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// try to unify all cutpoint addition into this function.
	void
	correctly_add_cutpoint_variants( core::pose::Pose & pose,
																	 Size const res_to_add,
																	 bool const check_fold_tree /* = true*/){

		using namespace core::chemical;

		runtime_assert( res_to_add < pose.total_residue() );
		if ( check_fold_tree ) runtime_assert( pose.fold_tree().is_cutpoint( res_to_add ) );

		if ( pose.residue_type( res_to_add ).has_variant_type( UPPER_TERMINUS ) ) {
			remove_variant_type_from_pose_residue( pose, UPPER_TERMINUS, res_to_add );
		}
		if ( pose.residue_type( res_to_add + 1 ).has_variant_type( LOWER_TERMINUS ) ) {
			remove_variant_type_from_pose_residue( pose, LOWER_TERMINUS, res_to_add + 1 );
		}

		if ( pose.residue_type( res_to_add ).is_RNA() ){
			// could also keep track of alpha, beta, etc.
			runtime_assert( pose.residue_type( res_to_add + 1 ).is_RNA() );
			protocols::swa::rna::correctly_position_cutpoint_phosphate_torsions( pose, res_to_add, false /*verbose*/ );
		}
		add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, res_to_add   );
		add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, res_to_add + 1 );
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

		Size num_five_prime_connections( 0 ), num_three_prime_connections( 0 );
		Size num_jumps_to_previous( 0 ), num_jumps_to_subsequent( 0 );

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

			bool const jump_to_subsequent = check_jump_to_subsequent_residue_in_chain( pose, j, slice_res );
			if ( jump_to_subsequent ) num_jumps_to_subsequent++;

			if ( n == 1 || !after_cutpoint ){
				sliced_out_pose.append_residue_by_bond( *rsd, true /* build_ideal_geometry */ );
			} else {
				sliced_out_pose.append_residue_by_jump( *rsd, sliced_out_pose.total_residue() );
			}
		}

		runtime_assert ( num_five_prime_connections <= 1 );
		runtime_assert ( num_three_prime_connections <= 1 );
		runtime_assert ( num_jumps_to_previous <= 1 );
		runtime_assert ( num_jumps_to_subsequent <= 1 );
		TR.Debug << num_five_prime_connections << " " <<  num_three_prime_connections << " " << num_jumps_to_previous << " " <<  num_jumps_to_subsequent << std::endl;
		// requirement for a clean slice:
		runtime_assert( (num_five_prime_connections + num_three_prime_connections + num_jumps_to_previous + num_jumps_to_subsequent) == 1 );

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
			jump_point_( 1, i ) = jump_partners1[ i ];
			jump_point_( 2, i ) = jump_partners2[ i ];
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

		sliced_out_pose.fold_tree( f_slice );

		// map (internal) coordinates from separate poses into merged on.
		std::map< Size, Size > res_map;
		for ( Size n = 1; n <= slice_res.size(); n++ ) res_map[ n ] =  slice_res[ n ];
		copy_dofs_match_atom_names( sliced_out_pose, pose, res_map, false /*backbone_only*/, false /*ignore_virtual*/ );
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
	check_jump_to_subsequent_residue_in_chain( pose::Pose const & pose, Size const i,
																					 utility::vector1< Size > const & current_element,
																					 FullModelInfo const & full_model_info ){
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		utility::vector1< Size > const & chains_in_full_model = full_model_info.chains_in_full_model();
		return 	check_jump_to_subsequent_residue_in_chain( pose, i, current_element, res_list, chains_in_full_model );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	check_jump_to_subsequent_residue_in_chain( pose::Pose const & pose, Size const i,
																					 utility::vector1< Size > const & current_element ){
		utility::vector1< Size > res_list, chains_in_full_model;
		for ( Size n = 1; n <= pose.total_residue(); n++ ) {
			res_list.push_back( n );
			chains_in_full_model.push_back( 1 ); // could easily fix this for virtual residue.
		}
		return 	check_jump_to_subsequent_residue_in_chain( pose, i, current_element, res_list, chains_in_full_model );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	check_jump_to_subsequent_residue_in_chain( pose::Pose const & pose, Size const i,
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

		// Could this be a chainbreak (cutpoint_closed )?
		TR.Debug << "checking for cutpoint after append: " << res << " " << res_list[ res ]  << " " << cutpoint_open_in_full_model.size() << std::endl;

		if ( res < pose.total_residue() &&
				 res_list[ res ] + 1 == res_list[ res + 1 ] &&
				 ! cutpoint_open_in_full_model.has_value( res_list[ res ]) ){

			// can happen after additions
			if ( pose.residue_type( res + 1 ).has_variant_type( VIRTUAL_PHOSPHATE) ) {
				remove_variant_type_from_pose_residue( pose, VIRTUAL_PHOSPHATE, res + 1 );
			}

			// can happen after additions
			correctly_add_cutpoint_variants( pose, res );

			// leave virtual riboses in this should actually get instantiated by the modeler
			//if ( pose.residue_type( res ).has_variant_type( "VIRTUAL_RIBOSE" ) )	remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", res );

		} else{

			// can happen after deletions
			if ( pose.residue_type( res ).has_variant_type( CUTPOINT_LOWER ) ){
				remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, res );
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

		// Could this be a chainbreak (cutpoint_closed )?

		TR.Debug << "checking for cutpoint after prepend: " << res << " " << res_list[ res ] << " " << cutpoint_open_in_full_model.size() << std::endl;

		if ( res > 1 &&
				 res_list[ res ] - 1 == res_list[ res - 1 ] &&
				 ! cutpoint_open_in_full_model.has_value( res_list[ res - 1 ])  ){

			// can happen after additions
			correctly_add_cutpoint_variants( pose, res - 1 );

			// This should actually be instantiated by the modeler.
			//if ( pose.residue_type( res ).has_variant_type( "VIRTUAL_RIBOSE" ) )	remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", res );

		} else {

			// can happen after additions
			if ( pose.residue_type( res ).is_RNA() && !pose.residue_type( res ).has_variant_type( "VIRTUAL_PHOSPHATE" ) ) {
				add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", res );
			}

			// can happen after deletions
			if ( pose.residue_type( res ).has_variant_type( CUTPOINT_UPPER ) ){
				remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, res );
			}

		}
	}

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

		if ( !pose.residue_type( res ).has_variant_type( "VIRTUAL_RIBOSE" ) )	{
			add_variant_type_to_pose_residue( pose, "VIRTUAL_RIBOSE", res );
		}

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
				if ( pose.residue_type( n ).has_variant_type( VIRTUAL_PHOSPHATE ) )	remove_variant_type_from_pose_residue( pose, VIRTUAL_PHOSPHATE, n );
				// don't remove virtual sugars -- those will be instantiated inside the residue sampler
				//				if ( pose.residue_type( n ).has_variant_type( "VIRTUAL_RIBOSE" ) )	remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", n );
				runtime_assert( !pose.residue_type( n ).has_variant_type( CUTPOINT_UPPER ) );
			}

			// Look for strand ends.
			bool const at_strand_end = ( n == pose.total_residue() || pose.fold_tree().is_cutpoint( n ) );
			if ( at_strand_end ) {
				fix_up_residue_type_variants_at_strand_end( pose, n );
			} else {
				// don't remove virtual sugars -- those will be instantiated inside the residue sampler
				//				if ( pose.residue_type( n ).has_variant_type( "VIRTUAL_RIBOSE" ) )	remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", n );
				runtime_assert( !pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) );
			}

			 // check for floating_base
			if ( at_strand_end && at_strand_beginning ) fix_up_residue_type_variants_at_floating_base( pose,  n );
		}

		pose_to_fix = pose;
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
	// specific for two helix test case. Sort of a unit test. Maybe I should make it a unit test.
	void
	test_merge_and_slice_with_two_helix_test_case( 	utility::vector1< core::pose::PoseOP > const & input_poses,
																									core::scoring::ScoreFunctionOP scorefxn ){

		using namespace core::pose;
		using namespace core::pose::full_model_info;
		using namespace protocols::swa;

		scorefxn->show( *( input_poses[ 1 ] ) );
		scorefxn->show( *( input_poses[ 2 ] ) );

		input_poses[1]->dump_pdb( "input_pose1.pdb" );
		input_poses[2]->dump_pdb( "input_pose2.pdb" );

		pose::Pose merged_pose;
		merge_two_poses_using_full_model_info( merged_pose, *(input_poses[ 2 ]), *( input_poses[ 1 ] ), 12 /*merge_res*/ );
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

	//////////////////////////////////////////////////////////////////////////////////////
	// might be better to move these into core (e.g., core::pose::full_model_info ),
	// or into a new protocols/full_model_setup/ directory.
	core::pose::PoseOP
	get_pdb_and_cleanup( std::string const input_file,
											 core::chemical::ResidueTypeSetCAP rsd_set )
	{
		using namespace core::pose;
		PoseOP input_pose = new Pose;
		import_pose::pose_from_pdb( *input_pose, *rsd_set, input_file );
		cleanup( *input_pose );
		return input_pose;
	}

	//////////////////////////////////////////////////////////////////////////////////////
	// might be better to move these into core (e.g., core::pose::full_model_info ),
	// or into a new protocols/full_model_setup/ directory.
	void
	cleanup( pose::Pose & pose ){
		protocols::rna::figure_out_reasonable_rna_fold_tree( pose );
		protocols::rna::virtualize_5prime_phosphates( pose );
	}

	//////////////////////////////////////////////////////////////////////////////////////
	// might be better to move these into core (e.g., core::pose::full_model_info ),
	// or into a new protocols/full_model_setup/ directory.
	void
	get_other_poses( utility::vector1< pose::PoseOP > & other_poses,
									 utility::vector1< std::string > const & other_files,
									core::chemical::ResidueTypeSetCAP rsd_set ){

		for ( Size n = 1; n <= other_files.size(); n++ ){
			other_poses.push_back( get_pdb_and_cleanup( other_files[ n ], rsd_set ) );
		}
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
		utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
		for ( Size i = 1; i <= partition_res.size(); i++ ){
			if ( domain_map[ partition_res[ i ] ] > 0 ) return true;
		}
		return false;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	check_for_fixed_domain( pose::Pose const & pose ){
		utility::vector1< Size > partition_res;
		for ( Size i = 1; i <= pose.total_residue(); i++ ) partition_res.push_back( i );
		return check_for_fixed_domain( pose, partition_res );
	}

}
}
