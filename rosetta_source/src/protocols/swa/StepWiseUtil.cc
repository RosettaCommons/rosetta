// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseUtil
/// @brief a few functions used by several StepWiseAnsatz classes
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/StepWiseUtil.hh>

//////////////////////////////////
#include <core/types.hh>

#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

#include <numeric/angle.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>

#include <string>

#include <core/pose/util.hh>
#include <utility/vector1.hh>



using core::Real;
using core::Size;
using core::pose::Pose;



namespace protocols {
namespace swa {

/////////////////////////////////////////////////////////////
	Real get_rotamer_angle( core::Size const & i, core::Size const & N_SAMPLE ){
		return  ( -180.0 + ( 360.0 / N_SAMPLE ) * i + 0.001 );
	}


  /////////////////////////////////////////////////////////////////////////////////////////////////
  void
  output_silent_struct( core::pose::Pose const & pose, core::pose::PoseCOP const & native_pose_op,
												std::string const & silent_file, std::string const & tag,
												core::io::silent::SilentFileDataOP sfd_in /* = 0 */ ){

    using namespace core::io::silent;
    using namespace core::scoring;

    BinaryProteinSilentStruct s( pose, tag );

		if ( native_pose_op != 0 ){

			core::pose::Pose pose_for_superpose = pose;
			core::pose::Pose const & native_pose_for_superpose( *native_pose_op );
			//remove_end_variants( pose_for_superpose );
			//remove_end_variants( native_pose_for_superpose );

			Real const rmsd = CA_rmsd( pose_for_superpose, native_pose_for_superpose );
			Real const backbone_rmsd = rmsd_with_super(  pose_for_superpose, native_pose_for_superpose, is_protein_backbone_including_O);
			//			Real const all_rmsd_old = all_atom_rmsd( pose_for_superpose, native_pose_for_superpose );
			Real const all_rmsd = rms_at_corresponding_heavy_atoms( pose_for_superpose, native_pose_for_superpose );

			s.add_energy( "rms", rmsd );
			//			s.add_energy( "all_rms_old", all_rmsd_old );
			s.add_energy( "all_rms", all_rmsd );
			s.add_energy( "backbone_rms", backbone_rmsd );
		}
    //	s.add_energy( "phi", pose.phi( n1 ) );
    //	s.add_energy( "psi", pose.psi( n1 ) );
    //	s.add_energy( "rama1", pose.energies().onebody_energies( n1 )[ rama ]);
    //	s.add_energy( "centroid_score", centroid_score );

    //SilentFileData silent_file_data;

		//    silent_file_data.write_silent_struct( s, silent_file+".sc", true /*write score only*/ );
    static const SilentFileData silent_file_data;
    silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );

		if ( sfd_in != 0 ) sfd_in->add_structure( s );

    //	pose.dump_pdb( tag+".pdb" );
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

		if ( pose.residue(res).has( " H  ") ) {
			phi = numeric::dihedral_radians(
																			pose.xyz( NamedAtomID( " C  ", res ) ),
																			pose.xyz( NamedAtomID( " CA ", res ) ),
																			pose.xyz( NamedAtomID( " N  ", res ) ),
																			pose.xyz( NamedAtomID( " H  ", res ) ) );
		} else if ( pose.residue(res).has("1H  ") ) {
			phi = numeric::dihedral_radians(
																			pose.xyz( NamedAtomID( " C  ", res ) ),
																			pose.xyz( NamedAtomID( " CA ", res ) ),
																			pose.xyz( NamedAtomID( " N  ", res ) ),
																			pose.xyz( NamedAtomID( "1H  ", res ) ) );
		}

		// This is phi, up to an offset.
		return numeric::principal_angle_degrees( 180.0 + numeric::conversions::degrees( phi ) );
	}

}
}
