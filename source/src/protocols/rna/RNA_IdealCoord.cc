// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file   protocols/rna/RNA_IdealCoord.cc
/// @brief  Apply ideal RNA geometry to a residue or a pose
/// @author  Fang-Chieh Chou

// Unit headers
#include <protocols/rna/RNA_IdealCoord.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/database/open.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/constants.hh>

//// C++ headers
#include <string>
#include <cmath>

using namespace core;
using namespace core::pose;
using namespace core::chemical::rna;

namespace protocols {
namespace rna {

//////////////////////////////////////////////////////
RNA_IdealCoord::RNA_IdealCoord():
	path_( basic::database::full_name("chemical/residue_type_sets/rna_phenix/ideal_geometry/") )
{
	init();
}

RNA_IdealCoord::~RNA_IdealCoord() {}

/////////////////////////////////////////////////////
bool RNA_IdealCoord::is_torsion_exists(Pose const & pose, id::TorsionID const & torsion_id) const {
	using namespace core::id;
	Size res_index = torsion_id.rsd();
	if ( res_index < 1 || res_index > pose.total_residue() ) return false;

	AtomID id1,id2,id3,id4;
	bool fail = pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
	return (! fail);
}
/////////////////////////////////////////////////////
void RNA_IdealCoord::init() {
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
	chemical::ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set ( "rna" );
	Pose ref_pose;
	for (Size i = 1; i <= pdb_file_list.size(); ++i) {
		import_pose::pose_from_pdb ( ref_pose, *rsd_set, pdb_file_list[i] );
		ref_pose_list_.push_back(ref_pose);
	}
}

/////////////////////////////////////////////////////
void RNA_IdealCoord::apply( Pose & pose, Size const seqpos, bool const is_north, bool const keep_backbone_torsion ) const {

	using namespace core::id;
	using namespace core::chemical;
	using namespace core::conformation;

	Residue const & res = pose.residue( seqpos );
	if ( !res.is_RNA() ) return;

	//Figure out the residue_type.
	Size res_class = 0;
	if ( res.aa() == na_rad ) {
		res_class = 1;
	} else if ( res.aa() == na_rgu ) {
		res_class = 3;
	} else if ( res.aa() == na_rcy ) {
		res_class = 5;
	} else if ( res.aa() == na_ura ) {
		res_class = 7;
	}
	if (! is_north) ++res_class;

	//Record the torsions in starting pose
	utility::vector1 < TorsionID > saved_torsion_id;
	utility::vector1 < Real > saved_torsions;
	if (keep_backbone_torsion) {
		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  ALPHA   ) );
		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  BETA    ) );
		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  GAMMA   ) );
		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  EPSILON ) );
		saved_torsion_id.push_back( TorsionID( seqpos,   id::BB,  ZETA    ) );
		saved_torsion_id.push_back( TorsionID( seqpos,   id::CHI, 1       ) );
		saved_torsion_id.push_back( TorsionID( seqpos,   id::CHI, 4       ) );
		saved_torsion_id.push_back( TorsionID( seqpos-1, id::BB,  ZETA    ) );
		saved_torsion_id.push_back( TorsionID( seqpos+1, id::BB,  ALPHA   ) );
		for (Size i = 1; i <= saved_torsion_id.size(); ++i) {
			bool const is_exists = is_torsion_exists( pose, saved_torsion_id[i] );
			if (is_exists) {
				saved_torsions.push_back( pose.torsion( saved_torsion_id[i] ) );
			} else {
				saved_torsions.push_back( -9999 );
			}
		}
	}

	//Apply ideal dofs
	std::map <Size, Size> res_map;
	res_map.insert( std::pair <Size, Size> (seqpos, 2) ); //Only the center res (#2) matters in ref_pose
	Pose const & ref_pose = ref_pose_list_[res_class];
	copy_dofs_match_atom_names(pose, ref_pose, res_map);

	//Copy back the original torsions
	if (keep_backbone_torsion) {
		for (Size i = 1; i <= saved_torsion_id.size(); ++i) {
			if (saved_torsions[i] > -1000) {
				pose.set_torsion( saved_torsion_id[i], saved_torsions[i] );
			}
		}
	}
}
/////////////////////////////////////////////////////
//Apply ideal coords to whole pose.
//pucker_conformations: 0 for skipping, 1 for North, 2 for South
void RNA_IdealCoord::apply( Pose & pose, utility::vector1 < Size > const & pucker_conformations, bool const keep_backbone_torsion ) const {
	if ( pose.total_residue() != pucker_conformations.size() ) {
		utility_exit_with_message("RNA_IdealCoord::apply--pose.total_resdiue() != pucker_conformations.size() !!!!");
	}
	for (Size i = 1; i <= pose.total_residue(); ++i) {
		switch ( pucker_conformations[i] ) {
			case 0: //Skip
				break;
			case 1: //North
				apply(pose, i, true, keep_backbone_torsion);
				break;
			case 2: //South
				apply(pose, i, false, keep_backbone_torsion);
				break;
			default :
				utility_exit_with_message("RNA_IdealCoord::apply--Invalid value for pucker_conformations!!!!");
		}
	}
}
/////////////////////////////////////////////////

}
}
