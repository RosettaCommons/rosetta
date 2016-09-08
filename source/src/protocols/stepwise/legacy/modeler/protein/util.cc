// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StepWiseProteinUtil
/// @brief a few functions used by several StepWiseProteinAnsatz classes
/// @details
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/legacy/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/file_util.hh>
#include <protocols/stepwise/modeler/util.hh>

//////////////////////////////////
#include <core/fragment/Frame.hh>
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
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.legacy.modeler.protein.util" );

using core::Real;
using core::Size;
using core::pose::Pose;
using utility::tools::make_vector1;
using ObjexxFCL::string_of;
using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::modeler::protein;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace protein {

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

/////////////////////////////////////////////////////////
void
output_pose_list( utility::vector1< core::pose::PoseOP > pose_list,
	core::pose::PoseCOP native_pose,
	std::string const & silent_file,
	utility::vector1< Size > const & working_calc_rms_res,
	bool const overwrite /* = true */ ){

	core::io::silent::SilentFileDataOP sfd( new core::io::silent::SilentFileData ); // silly
	if ( overwrite ) remove_silent_file_if_it_exists( silent_file );

	for ( Size n = 1; n <= pose_list.size(); n++ ) {
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
	utility::vector1< core::Size > const & calc_rms_res ){

	using namespace core::io::silent;
	using namespace core::scoring;

	BinarySilentStruct s( pose, tag );

	if ( native_pose_op != 0 ) {

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

///////////////////////////////////////////////////////////////////////////////////////
bool
is_close_chain_break(pose::Pose const & pose){
	for ( Size seq_num = 1; seq_num < pose.size(); seq_num++ ) {
		if ( is_cutpoint_closed( pose, seq_num ) ) return true;
	}
	return false;
}


} //protein
} //modeler
} //legacy
} //stepwise
} //protocols
