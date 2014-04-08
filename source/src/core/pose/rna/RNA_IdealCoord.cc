// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file   core/pose/rna/RNA_IdealCoord.cc
/// @brief  Apply ideal RNA geometry to a residue or a pose
/// @author  Fang-Chieh Chou

// Unit headers
#include <core/pose/rna/RNA_IdealCoord.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/io/pdb/file_data.hh>
#include <basic/database/open.hh>
#include <core/id/DOF_ID_Map.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/constants.hh>

//// C++ headers
#include <string>
#include <cmath>

using namespace core::chemical::rna;

namespace core {
namespace pose {
namespace rna {

//////////////////////////////////////////////////////
RNA_IdealCoord::RNA_IdealCoord():
	utility::pointer::ReferenceCount(),
	path_( basic::database::full_name("chemical/residue_type_sets/rna_phenix/ideal_geometry/") )
{
	init();
}

RNA_IdealCoord::~RNA_IdealCoord() {}

/////////////////////////////////////////////////////
bool RNA_IdealCoord::is_torsion_exists(
	Pose const & pose,
	id::TorsionID const & torsion_id
) const {
	using namespace id;
	Size res_index = torsion_id.rsd();
	if ( res_index < 1 || res_index > pose.total_residue() ) return false;

	AtomID id1,id2,id3,id4;
	bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
	return (! fail);
}
/////////////////////////////////////////////////////
void RNA_IdealCoord::init() {
	RNA_FittedTorsionInfo const torsion_info;
	delta_cutoff_ = torsion_info.delta_cutoff();

	//Names of the pdb files
	utility::vector1 < std::string > pdb_file_list;
	pdb_file_list.push_back( path_ + "/A_n_std.pdb" );
	pdb_file_list.push_back( path_ + "/A_s_std.pdb" );
	pdb_file_list.push_back( path_ + "/G_n_std.pdb" );
	pdb_file_list.push_back( path_ + "/G_s_std.pdb" );
	pdb_file_list.push_back( path_ + "/C_n_std.pdb" );
	pdb_file_list.push_back( path_ + "/C_s_std.pdb" );
	pdb_file_list.push_back( path_ + "/U_n_std.pdb" );
	pdb_file_list.push_back( path_ + "/U_s_std.pdb" );

	//Initialize the reference poses
	chemical::ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set(chemical::RNA);
	for ( Size i = 1; i <= pdb_file_list.size(); ++i ) {
		PoseOP ref_pose = new Pose();
		io::pdb::build_pose_from_pdb_as_is( *ref_pose, *rsd_set, pdb_file_list[i] );
		MiniPoseOP ref_mini_pose = new MiniPose( *ref_pose );
		ref_mini_pose_list_.push_back( ref_mini_pose );
	}
}

/////////////////////////////////////////////////////
//Apply ideal coords to a residue in pose.
//pucker_conformations: 0 for maintaining current, 1 for North, 2 for South
//std::map < id::DOF_ID , Real > RNA_IdealCoord::apply_and_return(
//	Pose & pose,
//	Size const seqpos,
//	Size pucker,
//	bool const keep_backbone_torsion
//) const {
//	using namespace id;
//	using namespace chemical;
//	using namespace conformation;
//
//	std::map < id::DOF_ID , Real > result;
//	assert( pucker <= 2 );
//
//	Residue const & res = pose.residue( seqpos );
//	if ( !res.is_RNA() ) return result;
//
//	if ( pucker == WHATEVER ) {
//		Real const delta  = pose.torsion( TorsionID(seqpos, id::BB, DELTA) );
//		if ( delta > delta_cutoff_ ) {
//			pucker = SOUTH;
//		} else {
//			pucker = NORTH;
//		}
//	}
//
//	//Figure out the residue_type.
//	Size res_class = 0;
//	switch ( res.aa() ) {
//		case na_rad:
//			res_class = 1; break;
//		case na_rgu:
//			res_class = 3; break;
//		case na_rcy:
//			res_class = 5; break;
//		case na_ura:
//			res_class = 7; break;
//		default:
//			utility_exit_with_message( "Invalid res.aa()!" );
//	}
//
//	if ( pucker == SOUTH ) ++res_class;
//
//	//Record the torsions in starting pose
//	utility::vector1 < TorsionID > saved_torsion_id;
//	utility::vector1 < Real > saved_torsions;
//	if ( keep_backbone_torsion ) {
//		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  ALPHA   ) );
//		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  BETA    ) );
//		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  GAMMA   ) );
//		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  EPSILON ) );
//		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  ZETA    ) );
//		saved_torsion_id.push_back( TorsionID( seqpos,   id::CHI, 1       ) ); //CHI
//		saved_torsion_id.push_back( TorsionID( seqpos,   id::CHI, 4       ) ); //O2H
//		saved_torsion_id.push_back( TorsionID( seqpos-1, id::BB,  ZETA    ) );
//		saved_torsion_id.push_back( TorsionID( seqpos+1, id::BB,  ALPHA   ) );
//		for ( Size i = 1; i <= saved_torsion_id.size(); ++i ) {
//			bool const is_exists = is_torsion_exists( pose, saved_torsion_id[i] );
//			if (is_exists) {
//				saved_torsions.push_back( pose.torsion( saved_torsion_id[i] ) );
//			} else {
//				saved_torsions.push_back( -9999 );
//			}
//		}
//	}
//
//	//Apply ideal dofs
//	std::map <Size, Size> res_map;
//	res_map.insert( std::pair <Size, Size> (seqpos, 2) ); //Only the center res (#2) matters in ref_pose
//	MiniPoseOP const ref_mini_pose = ref_mini_pose_list_[res_class];
//	result = copy_dofs_match_atom_names_and_return( pose, *ref_mini_pose, res_map );
//
//	//Copy back the original torsions
//	if ( keep_backbone_torsion ) {
//		for ( Size i = 1; i <= saved_torsion_id.size(); ++i ) {
//			if ( saved_torsions[i] > -1000 ) pose.set_torsion( saved_torsion_id[i], saved_torsions[i] );
//		}
//	}
//
//	return result;
//}
//
//////////////////////////
void RNA_IdealCoord::apply_coords(
		Pose & pose,
		Size const seqpos,
		Size const res_class,
		bool const ignore_base,
		bool const keep_backbone_torsion
) const {
	using namespace id;
	using namespace chemical;
	using namespace conformation;

	//Record the torsions in starting pose
	utility::vector1 < TorsionID > saved_torsion_id;
	utility::vector1 < Real > saved_torsions;
	if ( keep_backbone_torsion ) {
		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  ALPHA   ) );
		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  BETA    ) );
		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  GAMMA   ) );
		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  EPSILON ) );
		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  ZETA    ) );
		saved_torsion_id.push_back( TorsionID( seqpos,   id::CHI, 1       ) ); //CHI
		saved_torsion_id.push_back( TorsionID( seqpos,   id::CHI, 4       ) ); //O2H
		saved_torsion_id.push_back( TorsionID( seqpos-1, id::BB,  ZETA    ) );
		saved_torsion_id.push_back( TorsionID( seqpos+1, id::BB,  ALPHA   ) );
		for ( Size i = 1; i <= saved_torsion_id.size(); ++i ) {
			bool const is_exists = is_torsion_exists( pose, saved_torsion_id[i] );
			if (is_exists) {
				saved_torsions.push_back( pose.torsion( saved_torsion_id[i] ) );
			} else {
				saved_torsions.push_back( -9999 );
			}
		}
	}

	MiniPoseOP const ref_mini_pose = ref_mini_pose_list_[res_class];
	//Apply ideal dofs
	if ( ignore_base ) {
		std::map < core::id::AtomID , core::id::AtomID > atom_id_map;
		utility::vector1< std::string > const & ref_atom_names(
				ref_mini_pose->atom_names_list()[2] );
		utility::vector1< std::string > const & non_base_atoms(
				chemical::rna::non_base_atoms );
		chemical::ResidueType const & rsd_type1( pose.residue_type( seqpos ) );
		for ( Size i = 1; i <= non_base_atoms.size(); ++i ) {
			std::string const & atom_name( non_base_atoms[i] );
			Size const index2 =
					std::find( ref_atom_names.begin(), ref_atom_names.end(), atom_name ) -
					ref_atom_names.begin() + 1;
			if ( index2 <= ref_atom_names.size() ) {
				Size index1( 0 );
				for ( Size j = 1; j <= rsd_type1.natoms(); ++j ) {
					if ( rsd_type1.atom_name( j ) == atom_name ) {
						index1 = j;
						break;
					}
				}
				if ( index1 != 1 )
						atom_id_map[AtomID( index1, seqpos )] = AtomID( index2, 2 );
			}
		}
		copy_dofs( pose, *ref_mini_pose, atom_id_map );
	} else {
		std::map <Size, Size> res_map;
		res_map.insert( std::pair <Size, Size> (seqpos, 2) ); //Only the center res (#2) matters in ref_pose
		copy_dofs_match_atom_names( pose, *ref_mini_pose, res_map );
	}

	//Copy back the original torsions
	if ( keep_backbone_torsion ) {
		for ( Size i = 1; i <= saved_torsion_id.size(); ++i ) {
			if ( saved_torsions[i] > -1000 ) pose.set_torsion( saved_torsion_id[i], saved_torsions[i] );
		}
	}
}
//////////////////////////
void RNA_IdealCoord::apply(
		Pose & pose,
		Size const seqpos,
		PuckerState pucker,
		bool const keep_backbone_torsion
) const {
	using namespace id;
	using namespace chemical;
	using namespace conformation;

	Residue const & res = pose.residue( seqpos );
	if ( !res.is_RNA() ) return;
	runtime_assert( pucker <= 2 );

	if ( pucker == ANY_PUCKER ) {
		Real const delta  = pose.torsion( TorsionID(seqpos, id::BB, DELTA) );
		if ( delta > delta_cutoff_ ) {
			pucker = SOUTH;
		} else {
			pucker = NORTH;
		}
	}

	//Figure out the residue_type.
	Size res_class = 0;
	switch ( res.aa() ) {
		case na_rad:
			res_class = 1; break;
		case na_rgu:
			res_class = 3; break;
		case na_rcy:
			res_class = 5; break;
		case na_ura:
			res_class = 7; break;
		default:
			utility_exit_with_message( "Invalid res.aa()!" );
	}

	if ( pucker == SOUTH ) ++res_class;
	apply_coords( pose, seqpos,
			res_class, false /*ignore_base*/, keep_backbone_torsion );
}
/////////////////////////////////////////////////////
void RNA_IdealCoord::apply_pucker(
		Pose & pose,
		Size const seqpos,
		PuckerState pucker,
		bool const keep_backbone_torsion
) const {
	using namespace id;
	using namespace chemical;
	using namespace conformation;

	assert( pucker <= 2 );

	Residue const & res = pose.residue( seqpos );
	if ( !res.is_RNA() ) return;

	if ( pucker == ANY_PUCKER ) {
		Real const delta  = pose.torsion( TorsionID(seqpos, id::BB, DELTA) );
		if ( delta > delta_cutoff_ ) {
			pucker = SOUTH;
		} else {
			pucker = NORTH;
		}
	}

	// Assume the pucker coord of A for all bases
	Size res_class = 1;
	if ( pucker == SOUTH ) ++res_class;
	apply_coords( pose, seqpos,
			res_class, true /*ignore_base*/, keep_backbone_torsion );
}
/////////////////////////////////////////////////////
//Apply ideal coords to whole pose.
//pucker_conformations: 0 for maintaining current, 1 for North, 2 for South
void RNA_IdealCoord::apply(
	Pose & pose,
	utility::vector1 < PuckerState > const & puckers,
	bool const keep_backbone_torsion
) const {
	assert ( pose.total_residue() == puckers.size() );
	for ( Size i = 1; i <= pose.total_residue(); ++i )
			apply( pose, i, puckers[i], keep_backbone_torsion );
}
/////////////////////////////////////////////////////
//Apply ideal coords to whole pose, maintain current pucker state.
void RNA_IdealCoord::apply(
	Pose & pose,
	bool const keep_backbone_torsion
) const {
	for ( Size i = 1; i <= pose.total_residue(); ++i )
			apply( pose, i, ANY_PUCKER, keep_backbone_torsion );
}
/////////////////////////////////////////////////

}
}
}

