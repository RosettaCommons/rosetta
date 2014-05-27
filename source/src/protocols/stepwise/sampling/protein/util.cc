// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinUtil
/// @brief a few functions used by several StepWiseProteinAnsatz classes
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/sampling/protein/util.hh>
#include <protocols/stepwise/sampling/util.hh>

//////////////////////////////////
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/types.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>

#include <numeric/angle.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>
#include <ObjexxFCL/format.hh>

#include <string>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>


#include <utility/numbers.hh>

using core::Real;
using core::Size;
using core::pose::Pose;
using utility::tools::make_vector1;
using ObjexxFCL::string_of;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

/////////////////////////////////////////////////////////////
	Real get_rotamer_angle( core::Size const & i, core::Size const & N_SAMPLE ){
		return  ( -180.0 + ( 360.0 / N_SAMPLE ) * i + 0.001 );
	}



	/////////////////////////////////////////////////////////
	void
	output_pose_list( utility::vector1< core::pose::PoseOP > pose_list,
										core::pose::PoseCOP native_pose,
										std::string const & silent_file,
										utility::vector1< Size > const & working_calc_rms_res	){
		core::io::silent::SilentFileDataOP sfd = new core::io::silent::SilentFileData; // silly
		for ( Size n = 1; n <= pose_list.size(); n++ ){
			Pose & pose = *pose_list[ n ];
			std::string const tag = "S_"+ ObjexxFCL::string_of( n-1 /* start with zero */);
			output_silent_struct( pose, native_pose,
														silent_file,
														tag, sfd, working_calc_rms_res );
		}
	}


  /////////////////////////////////////////////////////////////////////////////////////////////////
  void
  output_silent_struct( core::pose::Pose const & pose, core::pose::PoseCOP const & native_pose_op,
												std::string const & silent_file, std::string const & tag,
												core::io::silent::SilentFileDataOP sfd_in /* = 0 */ ){

		utility::vector1< core::Size > calc_rms_res; //blank vector
		output_silent_struct( pose, native_pose_op, silent_file, tag, sfd_in, calc_rms_res );
	}

  /////////////////////////////////////////////////////////////////////////////////////////////////
  void
  output_silent_struct( core::pose::Pose const & pose, core::pose::PoseCOP const & native_pose_op,
												std::string const & silent_file, std::string const & tag,
												core::io::silent::SilentFileDataOP sfd_in,
												utility::vector1< core::Size > const & calc_rms_res	){

    using namespace core::io::silent;
    using namespace core::scoring;

    BinarySilentStruct s( pose, tag );

		if ( native_pose_op != 0 ){

			core::pose::Pose pose_for_superpose = pose; // is this necessary? quite a bit of overhead.
			// core::pose::Pose const & native_pose_for_superpose( *native_pose_op ); // Unused variable causes warning.

			Real rmsd( -1.0), backbone_rmsd( -1.0), all_rmsd( -1.0 );

			std::map< core::id::AtomID, core::id::AtomID > CA_map, backbone_heavy_map, heavy_map;
			setup_matching_CA_atoms(                     pose, *native_pose_op, CA_map );
			setup_matching_protein_backbone_heavy_atoms( pose, *native_pose_op, backbone_heavy_map );
			setup_matching_heavy_atoms(                  pose, *native_pose_op, heavy_map );

			if ( calc_rms_res.size() == 0 ) { // superimpose over everything
				rmsd = rms_at_corresponding_atoms( pose, *native_pose_op, CA_map );
				backbone_rmsd = rms_at_corresponding_atoms( pose, *native_pose_op, backbone_heavy_map );
				all_rmsd = rms_at_corresponding_atoms( pose, *native_pose_op, heavy_map );
			} else {
				rmsd = rms_at_corresponding_atoms_no_super( pose, *native_pose_op, CA_map, calc_rms_res );
				backbone_rmsd = rms_at_corresponding_atoms_no_super( pose, *native_pose_op, backbone_heavy_map, calc_rms_res );
				all_rmsd = rms_at_corresponding_atoms_no_super( pose, *native_pose_op, heavy_map, calc_rms_res );
			}

			if ( !utility::is_nan( rmsd ) ) s.add_energy( "rms", rmsd );
			s.add_energy( "all_rms", all_rmsd );
			if ( !utility::is_nan( backbone_rmsd ) ) s.add_energy( "backbone_rms", backbone_rmsd );
		}

    static const SilentFileData silent_file_data;
    if ( silent_file.size() > 0 ) silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );
		if ( sfd_in != 0 ) sfd_in->add_structure( s );

  }

	////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	remove_end_variants( core::pose::Pose & pose ){

		using namespace core::chemical;
		core::pose::remove_variant_type_from_pose_residue( pose, "LOWER_TERMINUS", 1  );
		core::pose::remove_variant_type_from_pose_residue( pose, "N_ACETYLATION", 1  );

		core::pose::remove_variant_type_from_pose_residue( pose, "UPPER_TERMINUS", pose.total_residue()  );
		core::pose::remove_variant_type_from_pose_residue( pose, "C_METHYLAMIDATION", pose.total_residue()  );

	}


	///////////////////////////////////////////////////////////////////
	// sigh. bookkeeping to keep junctions inherited properly...
	Real get_pretend_psi_explicit( core::pose::Pose const & pose, Size const & res )
	{
		// Psi may not be defined in pose properly because there's no "next residue".
		// But... we kind of know phi anyway, because the carbonyl oxygen is in the
		// current residue. Would have been better to have phi, psi, etc. defined *internally*
		// for each residue, using the O for psi, H for phi, etc., and/or virtual atoms for
		// other crazy torsions.
		using namespace core::id;
		Real const psi = numeric::dihedral_radians(
																							 pose.xyz( NamedAtomID( " O  ", res ) ),
																							 pose.xyz( NamedAtomID( " C  ", res ) ),
																							 pose.xyz( NamedAtomID( " CA ", res ) ),
																							 pose.xyz( NamedAtomID( " N  ", res ) ) );
		// This is phi, up to an offset.
		return numeric::principal_angle_degrees( 180.0 + numeric::conversions::degrees( psi ) );
	}

	///////////////////////////////////////////////////////////////////
	// sigh. bookkeeping to keep junctions inherited properly...
	Real get_pretend_phi_explicit( core::pose::Pose const & pose, Size const & res )
	{
		// Psi may not be defined in pose properly because there's no "next residue".
		// But... we kind of know phi anyway, because the carbonyl oxygen is in the
		// current residue. Would have been better to have phi, psi, etc. defined *internally*
		// for each residue, using the O for psi, H for phi, etc., and/or virtual atoms for
		// other crazy torsions.
		using namespace core::id;
		Real phi( 0.0 );

		if ( pose.residue(res).has( " H  " ) ) {
			phi = numeric::dihedral_radians(
																			pose.xyz( NamedAtomID( " C  ", res ) ),
																			pose.xyz( NamedAtomID( " CA ", res ) ),
																			pose.xyz( NamedAtomID( " N  ", res ) ),
																			pose.xyz( NamedAtomID( " H  ", res ) ) );
		} else if ( pose.residue(res).has( "1H  " ) ) {
			phi = numeric::dihedral_radians(
																			pose.xyz( NamedAtomID( " C  ", res ) ),
																			pose.xyz( NamedAtomID( " CA ", res ) ),
																			pose.xyz( NamedAtomID( " N  ", res ) ),
																			pose.xyz( NamedAtomID( "1H  ", res ) ) );
		}

		// This is phi, up to an offset.
		return numeric::principal_angle_degrees( 180.0 + numeric::conversions::degrees( phi ) );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////
 	void
 	setup_protein_CA_atom_id_map( core::pose::Pose const & pose_1,
															  core::pose::Pose const & pose_2,
																Size const base_res,
																core::id::AtomID_Map< core::id::AtomID > & atom_ID_map){

		using namespace core::id;

		if ( !pose_1.residue_type( base_res ).is_protein() ) return;
		if ( !pose_2.residue_type( base_res ).is_protein() ) return;

		AtomID atom1(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " CA ", base_res ), pose_1 ));
		AtomID atom2(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " CA ", base_res ), pose_2 ));
		atom_ID_map.set( atom1, atom2 );

 	}

	//////////////////////////////////////////////////////////////////////////////////////////////////
 	void
  	setup_protein_backbone_atom_id_map( core::pose::Pose const & pose_1,
  																			core::pose::Pose const & pose_2,
  																			Size const base_res,
  																			core::id::AtomID_Map< core::id::AtomID > & atom_ID_map){
 		setup_protein_backbone_atom_id_map( pose_1, pose_2, base_res, atom_ID_map );
 	}

 	//////////////////////////////////////////////////////////////////////////////////////////////////
  	void
  	setup_protein_backbone_atom_id_map( core::pose::Pose const & pose_1,
 																			core::pose::Pose const & pose_2,
 																			Size const base_res,
 																			Size const base_res2,
 																			core::id::AtomID_Map< core::id::AtomID > & atom_ID_map){

		using namespace core::id;

 		if( base_res == 0 ) return;
 		if( base_res2 == 0 ) return;

		if ( !pose_1.residue_type( base_res ).is_protein() ) return;
		if ( !pose_2.residue_type( base_res2 ).is_protein() ) return;

		if ( pose_1.residue_type( base_res ).has_variant_type( "VIRTUAL_RESIDUE" ) ) return;
		if ( pose_2.residue_type( base_res2 ).has_variant_type( "VIRTUAL_RESIDUE" ) ) return;

		{
			AtomID atom1(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " CA ", base_res ), pose_1 ));
			AtomID atom2(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " CA ", base_res2 ), pose_2 ));
			atom_ID_map.set( atom1, atom2 );
		}

		{
			AtomID atom1(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " C  ", base_res ), pose_1 ));
			AtomID atom2(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " C  ", base_res2 ), pose_2 ));
			atom_ID_map.set( atom1, atom2 );
		}

		{
			AtomID atom1(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " N  ", base_res ), pose_1 ));
			AtomID atom2(  core::pose::named_atom_id_to_atom_id( NamedAtomID( " N  ", base_res2 ), pose_2 ));
			atom_ID_map.set( atom1, atom2 );
		}

 	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	figure_out_protein_modeling_info( pose::Pose const & pose,
																		Size const moving_res,
																		utility::vector1< Size > & moving_res_list ){

		if ( !pose.residue( moving_res ).is_protein() ) return;

		// go back another residue -- this was the default choice in protein SWA.
		utility::vector1< Size > const & fixed_domain_map = core::pose::full_model_info::get_fixed_domain_from_full_model_info_const( pose );
		Size const upstream_res = pose.fold_tree().get_parent_residue( moving_res );
		if ( pose.residue_type( upstream_res ).is_protein() &&
				 fixed_domain_map[ upstream_res ] == 0 ) moving_res_list.push_back( upstream_res );
	}


	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	get_bridge_res( pose::Pose const & pose,
									utility::vector1< Size > const & moving_res_list /*working*/ ){

		// check if CCD-closure  is necessary.
		utility::vector1< Size > bridge_res;
		if ( moving_res_list.size() == 0 ) return bridge_res;

		Size const moving_res = moving_res_list[ 1 ];
		utility::vector1< Size > const cutpoints_closed = figure_out_moving_cutpoints_closed_from_moving_res( pose, moving_res );
		if ( cutpoints_closed.size() == 0 ) return bridge_res;

		for ( Size n = 1; n <= cutpoints_closed.size(); n++ ) {

			int cutpoint_closed = static_cast<int>( cutpoints_closed[ n ] );
			if ( !pose.residue_type( cutpoint_closed ).is_protein() ) continue;

			// as in protein SWA, choose two bridge residues for CCD closure on the 'other side' of sampled residue.
			utility::vector1< int > offsets;
			utility::vector1< Size > const & fixed_domain_map = core::pose::full_model_info::get_fixed_domain_from_full_model_info_const( pose );
			if ( moving_res_list.has_value( cutpoint_closed ) ){
				offsets = make_vector1( +1, +2 );
			} else if ( moving_res_list.has_value( cutpoint_closed+1 ) ){
				offsets = make_vector1( -1, 0 );
			} else {
				offsets = make_vector1( 0, +1 ); // bracket CCD closure point, since moving_res does not.
			}
			utility::vector1< Size > working_bridge_res;
			for ( Size n = 1; n <= offsets.size(); n++ ) {
				int const bridge_res = cutpoint_closed + offsets[n];
				if( bridge_res < 1 ) continue;
				runtime_assert( !moving_res_list.has_value( bridge_res ) );
				if ( fixed_domain_map[ bridge_res ] == 0 ) working_bridge_res.push_back( bridge_res );
			}
			bridge_res = merge_vectors( bridge_res, const_full_model_info( pose ).sub_to_full( working_bridge_res ) );
		}
		return bridge_res;

	}


//////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
just_protein( utility::vector1< Size > const & res_list, pose::Pose const & pose ){
	utility::vector1< Size > protein_res_list;
	for ( Size n = 1; n <= res_list.size(); n++ ) { if ( pose.residue_type( res_list[n] ).is_protein() ) protein_res_list.push_back( res_list[n] ); }
	return protein_res_list;
}

//////////////////////////////////////////////////////////////////////////
bool
contains_protein( core::pose::Pose const & pose ){
	utility::vector1< Size > protein_res_list;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) { if ( pose.residue_type( n ).is_protein() ) return true; }
	return false;
}


} //protein
} //sampling
} //stepwise
} //protocols
