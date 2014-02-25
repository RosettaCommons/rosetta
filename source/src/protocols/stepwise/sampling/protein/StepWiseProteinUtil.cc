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
#include <protocols/stepwise/sampling/protein/StepWiseProteinUtil.hh>

//////////////////////////////////
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/types.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>

#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>

#include <numeric/angle.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>

#include <string>

//Auto Headers
#include <utility/vector1.hh>


using core::Real;
using core::Size;
using core::pose::Pose;



namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

/////////////////////////////////////////////////////////////
	Real get_rotamer_angle( core::Size const & i, core::Size const & N_SAMPLE ){
		return  ( -180.0 + ( 360.0 / N_SAMPLE ) * i + 0.001 );
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

    BinaryProteinSilentStruct s( pose, tag );

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

				//rmsd = CA_rmsd( pose_for_superpose, native_pose_for_superpose );
				//backbone_rmsd = rmsd_with_super(  pose_for_superpose, native_pose_for_superpose, is_protein_backbone_including_O);
				//all_rmsd = rms_at_corresponding_heavy_atoms( pose_for_superpose, native_pose_for_superpose );
			} else {
				rmsd = rms_at_corresponding_atoms_no_super( pose, *native_pose_op, CA_map, calc_rms_res );
				backbone_rmsd = rms_at_corresponding_atoms_no_super( pose, *native_pose_op, backbone_heavy_map, calc_rms_res );
				all_rmsd = rms_at_corresponding_atoms_no_super( pose, *native_pose_op, heavy_map, calc_rms_res );
			}

			s.add_energy( "rms", rmsd );
			s.add_energy( "all_rms", all_rmsd );
			s.add_energy( "backbone_rms", backbone_rmsd );
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
	fragment_set_slice( core::fragment::ConstantLengthFragSetOP & fragset,
											utility::vector1< core::Size > const & slice_res ){

		using namespace core::fragment;

		Size const len( fragset->max_frag_length() );

		ConstantLengthFragSetOP fragset_new = new ConstantLengthFragSet;

		for ( Size n = 1; n <= (slice_res.size() - len + 1); n++ ) {

			Size const & pos = slice_res[ n ];

			FrameList frames;

			if ( pos > (fragset->max_pos()-len+1) ) {
				std::cout << "WARNING: NO FRAGS FOR POSITION " << pos << std::endl;
				continue;
			}

			fragset->frames( pos, frames );

			// CURRENTLY ONLY WORKS FOR CONST FRAG LENGTH SETS!!!! ASSUMES ONE FRAME!!!
			assert( frames.size() == 1 );

			FrameOP & frame( frames[1] );
			FrameOP frame_new = new Frame( n, len );

			for ( Size n = 1; n <= frame->nr_frags(); n++ ) {
				frame_new->add_fragment( frame->fragment_ptr( n ) );
			}

			fragset_new->add( frame_new );

		}

		fragset = fragset_new;

	}




} //protein
} //sampling
} //stepwise
} //protocols
