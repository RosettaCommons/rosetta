// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock.cc
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#include <protocols/ligand_docking/GALigandDock/LigandConformer.hh>
#include <protocols/ligand_docking/GALigandDock/TorsionSampler.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/util.hh>
#include <core/pose/init_id_map.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <numeric/Quaternion.hh>
#include <numeric/angle.functions.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_xyz.hh>
#include <numeric/conversions.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID.hh>
#include <core/id/types.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray3D.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <numeric/fourier/SHT.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

static basic::Tracer TR( "protocols.ligand_docking.GALigandDock.LigandConformer" );

LigandConformer::LigandConformer() {
	init_params();
}

LigandConformer::~LigandConformer() {}

LigandConformer::LigandConformer(
	core::pose::PoseCOP pose,
	utility::vector1< core::Size > const &ligids,
	utility::vector1< core::Size > movingscs,
	bool freeze_ligand_backbone,
	bool freeze_ligand
) {
	init_params();
	initialize( pose, ligids, movingscs, freeze_ligand_backbone, freeze_ligand );
}

void
LigandConformer::init_params() {
	torsmutationRate_ = 0.1;
	rtmutationRate_ = 0.5;
	transmutWidth_ = 2.0;
	rotmutWidth_ = 60.0;
	ligchimutWidth_ = 120.0;
	generation_tag_ = "";
	rg_ = 0.0;
	sample_ring_conformers_ = true;
	jumpid_ = 0;
	negTdS_ = 0;
	freeze_ligand_backbone_ = false;
	freeze_ligand_ = false;
}

void
LigandConformer::initialize(
	core::pose::PoseCOP pose,
	utility::vector1< core::Size > const &ligids,
	utility::vector1< core::Size > movingscs,
	bool freeze_ligand_backbone,
	bool freeze_ligand
) {
	ligids_ = ligids;
	ref_pose_ = pose;
	movingscs_ = movingscs;
	freeze_ligand_backbone_ = ( freeze_ligand_backbone || freeze_ligand );
	freeze_ligand_ = freeze_ligand;

	// 1) get jumpid
	utility::vector1< core::kinematics::Edge > jumps = pose->fold_tree().get_jump_edges();
	// skip if ligand-only case
	if ( jumps.size() >= 1 ) {
		for ( auto j : jumps ) {
			core::Size up = j.start(), down = j.stop();
			if ( std::find( ligids_.begin(), ligids_.end(), up) == ligids_.end()
					&& std::find( ligids_.begin(), ligids_.end(), down) != ligids_.end()
					) {
				runtime_assert(jumpid_==0);
				jumpid_ = j.label();
			}
		}
		runtime_assert(jumpid_!=0);
	}

	// 2) set ligand name
	ligand_typenames_.clear();
	for ( auto ligid : ligids ) {
		ligand_typenames_.push_back( pose->residue_type(ligid).name() );
	}
	update_conf( pose );

	// get ligand chi dependence... this is a stupid way, should have a better way ...
        core::conformation::Residue ligand( pose->residue( ligids[1] ) );
        core::Size nligchi( ligand.nchi());
        ligandchi_downstream_.resize( nligchi );

        utility::vector1< core::Size > atmindex_defining_chi( nligchi, 0 );

        for ( core::Size ichi = 1; ichi <= nligchi; ++ichi ) {
                utility::vector1< core::Size > const &chiatms = ligand.chi_atoms( ichi );
                atmindex_defining_chi[ichi] = (ligand.atom_base(chiatms[2]) == chiatms[3])? chiatms[1] : chiatms[4];
        }

        // WARNING! This assumes nbr atom is atom tree root.
        // THIS IS NOT NECESSARILY TRUE!
        // if not true this will hang
        for ( core::Size ichi = 1; ichi <= nligchi; ++ichi ) {
                core::Size iatm( atmindex_defining_chi[ichi] );
                core::Size ibase( ligand.atom_base(iatm) );
                while ( ibase != ligand.nbr_atom() ) { // recurrsive until reaches to nbr atom
                        if ( atmindex_defining_chi.contains(iatm) ) {
                                core::Size chi_parent_of_ichi = atmindex_defining_chi.index_of(iatm);
                                if ( chi_parent_of_ichi != ichi ) ligandchi_downstream_[chi_parent_of_ichi].push_back( ichi );
                        }
                        iatm = ibase;
                        ibase = ligand.atom_base(iatm);
                }
        }

	for ( core::Size ichi = 1; ichi <= ligandchi_downstream_.size(); ++ichi ) {
		utility::vector1 <core::Size > downstream_chis = ligandchi_downstream_[ichi];
		utility::vector1< core::Size > const &chiatms = ligand.chi_atoms( ichi );
		for ( core::Size jchi = 1; jchi <= downstream_chis.size(); ++jchi ) {
			utility::vector1< core::Size > const &jchiatms = ligand.chi_atoms( downstream_chis[jchi] );
		}
	}

	core::pose::PoseOP ref_pose_ligand ( new core::pose::Pose );
	for ( core::Size i=1; i<=ligids.size(); ++i ) {
		if ( i==1 ) {
			ref_pose_ligand->append_residue_by_jump( pose->residue( ligids[i] ), ref_pose_ligand->total_residue() );
		} else {
			ref_pose_ligand->append_residue_by_bond( pose->residue( ligids[i] ) );
		}
	}
	core::pose::addVirtualResAsRoot( *ref_pose_ligand );

	ref_pose_ligand_ = ref_pose_ligand;
}

void
LigandConformer::update_ligchi_types( core::conformation::Residue const& ligres ){
	ligandchi_types_.resize(ligres.nchi());
	TorsionType torsion_type;
	for ( core::Size i=1; i<=ligres.nchi(); ++i ) {
		torsion_type.at1 = ligres.atom_type_index(ligres.chi_atoms(i)[1]);
		torsion_type.at2 = ligres.atom_type_index(ligres.chi_atoms(i)[2]);
		torsion_type.at3 = ligres.atom_type_index(ligres.chi_atoms(i)[3]);
		torsion_type.at4 = ligres.atom_type_index(ligres.chi_atoms(i)[4]);
		torsion_type.bn = ligres.type().bond_type(ligres.chi_atoms(i)[2], ligres.chi_atoms(i)[3]);
		torsion_type.br = ligres.type().bond_ringness(ligres.chi_atoms(i)[2], ligres.chi_atoms(i)[3]);
		ligandchi_types_[i] = torsion_type;
	}
}

// pose -> gene
// update internal representation based on this conformation
void
LigandConformer::update_conf( core::pose::PoseCOP pose ) {
	using namespace core::id;

	rb_.resize(7);
	proteinchis_.resize(movingscs_.size());
	proteinrestypes_.resize(movingscs_.size());

	// skip if ligand-only case
	if ( jumpid_ != 0 ) {
		core::kinematics::Jump const &ligjump = pose->jump( jumpid_ );
		numeric::xyzVector< core::Real > T = ligjump.get_translation();
		numeric::Quaternion< core::Real > Q;
		numeric::R2quat( ligjump.get_rotation(), Q );
		rb_[1] = Q.w(); rb_[2] = Q.x(); rb_[3] = Q.y(); rb_[4] = Q.z();
		rb_[5] = T.x(); rb_[6] = T.y(); rb_[7] = T.z();
	}

	// movable sidechains
	for ( core::Size i=1; i<= movingscs_.size(); ++i ) {
		core::Size nchi = pose->residue( movingscs_[i] ).nchi();
		utility::vector1< core::Real > chis_i(nchi,0.0);
		for ( core::Size j=1; j<=nchi; ++j ) {
			chis_i[j] = pose->chi( j, movingscs_[i] );
		}
		proteinchis_[i] = chis_i;
		proteinrestypes_[i] = pose->residue(movingscs_[i]).type_ptr();
	}

	// internal torsions
	ligandchis_.clear();
	ligandchi_types_.clear();
	ligandtorsionids_.clear();
	TorsionType torsion_type;
	AtomID atid1, atid2, atid3, atid4;
	for ( core::Size i : ligids_ ) {
		if ( freeze_ligand_ ) continue;
		core::conformation::Residue const & ligres = pose->residue(i);
		ligid_restype_map_[i] = ligres.type_ptr();
		for ( core::Size r=1; r<=2; ++r ) {
			if ( freeze_ligand_backbone_ && r == 1 ) continue;
			core::Size const n_torsions( r==1 ?ligres.mainchain_atoms().size() : ligres.nchi() );
			for  ( core::Size j=1; j <= n_torsions; ++j ) {
				core::id::TorsionType const id_tor_type( r == 1 ? core::id::BB : core::id::CHI );
				TorsionID const tor_id(i, id_tor_type, j);
				bool const fail( pose->conformation().get_torsion_angle_atom_ids(tor_id, atid1, atid2, atid3, atid4) );
				if ( fail ) {
					if ( ( r == 1 ) &&
							( ( j == 1 && pose->fold_tree().is_cutpoint( i-1 ) ) ||
							( j >= n_torsions-1 && pose->fold_tree().is_cutpoint( i ) ) ) ) continue;
					if ( TR.Debug.visible() ) TR.Debug << " missed torsion: " << tor_id << std::endl;
					continue;
				}
				torsion_type.at1 = ligres.atom_type_index( atid1.atomno() );
				torsion_type.at2 = ligres.atom_type_index( atid2.atomno() );
				torsion_type.at3 = ligres.atom_type_index( atid3.atomno() );
				torsion_type.at4 = ligres.atom_type_index( atid4.atomno() );
				torsion_type.bn = ligres.type().bond_type( atid2.atomno(), atid3.atomno() );
				torsion_type.br = ligres.type().bond_ringness( atid2.atomno(), atid3.atomno() );
				ligandchi_types_.push_back( torsion_type );
				ligandchis_.push_back( r == 1 ? ligres.mainchain_torsion(j): ligres.chi(j) );
				ligandtorsionids_.push_back(tor_id);
			}
		}
	}

	// internal ring torsions
	if ( !freeze_ligand_ && sample_ring_conformers_ && ligids_.size() == 1 ) {
		// *** for now assume this only used in single-ligand case
		//   it might affect some non-canonicals
		core::chemical::ResidueType const &ligrt = pose->residue(ligids_[1]).type();
		core::Size nnu = 0;
		core::Size nligring( ligrt.n_rings() );
		for ( core::Size j=1; j<=nligring; ++j ) {
			nnu += ligrt.ring_atoms( j ).size() - 1;
		}
		ligandnus_.resize(nnu);
		ligandtaus_.resize(nnu+nligring);
		core::Size offset = 0;
		for ( core::Size j=1; j<=nligring; ++j ) {
			nnu = ligrt.ring_atoms( j ).size() - 1;
			for ( core::Size k=1; k<=nnu; ++k ) {
				ligandnus_[offset+k] = pose->torsion( core::id::TorsionID( ligids_[1], core::id::NU, offset + k ) );
				core::id::AtomID a1( ligrt.nu_atoms( offset + k )[ 1 ], ligids_[1] );
				core::id::AtomID a2( ligrt.nu_atoms( offset + k )[ 2 ], ligids_[1] );
				core::id::AtomID a3( ligrt.nu_atoms( offset + k )[ 3 ], ligids_[1] );
				ligandtaus_[offset+(j-1)+k] = numeric::conversions::degrees( pose->conformation().bond_angle( a1, a2, a3 ) );
			}
			core::id::AtomID a1( ligrt.nu_atoms( offset + nnu )[ 2 ], ligids_[1] );
			core::id::AtomID a2( ligrt.nu_atoms( offset + nnu )[ 3 ], ligids_[1] );
			core::id::AtomID a3( ligrt.nu_atoms( offset + nnu )[ 4 ], ligids_[1] );
			ligandtaus_[offset+j+nnu] = numeric::conversions::degrees( pose->conformation().bond_angle( a1, a2, a3 ) );
			offset += nnu;
		}
	}

	// copy ligandxyz
	if ( !freeze_ligand_ || rg_ == 0.0 ) {
		ligandxyz_.clear();
		core::Vector com( 0.0 );
		for ( core::Size i : ligids_ ) {
			core::conformation::Residue const &lig( pose->residue(i) );
			for ( core::Size iatm = 1; iatm <= lig.nheavyatoms(); ++iatm ) {
				ligandxyz_.push_back(lig.xyz( iatm ));
				com += lig.xyz( iatm );
			}
		}
		com /= ligandxyz_.size();

		// compute radius of gyration
		rg_ = 0.0;
		for ( auto atm : ligandxyz_ ) {
			rg_ += com.distance_squared( atm );
		}
		rg_ /= ligandxyz_.size();
		rg_ = std::sqrt(rg_);
	}

	ligandxyz_synced_ = true;
}

void
LigandConformer::update_ligand_conf( core::pose::PoseCOP pose ) {
	using namespace core::id;

	if ( pose->size() != ligids_.size() ) {
		utility_exit_with_message("The size of pose used in update_ligand_conf is different than the ligand size.");
	}

	// internal torsions
	ligandchis_.clear();
	ligandchi_types_.clear();
	ligandtorsionids_.clear();
	TorsionType torsion_type;
	AtomID atid1, atid2, atid3, atid4;
	for ( core::Size i=1;  i<= ligids_.size(); ++i ) {
		if ( freeze_ligand_ ) continue;
		core::conformation::Residue const & ligres = pose->residue(i);
		ligid_restype_map_[ ligids_[i] ] = ligres.type_ptr();
		for ( core::Size r=1; r<=2; ++r ) {
			if ( freeze_ligand_backbone_ && r == 1 ) continue;
			core::Size const n_torsions( r==1 ?ligres.mainchain_atoms().size() : ligres.nchi() );
			for  ( core::Size j=1; j <= n_torsions; ++j ) {
				core::id::TorsionType const id_tor_type( r == 1 ? core::id::BB : core::id::CHI );
				TorsionID const tor_id_curr(i, id_tor_type, j);
				TorsionID const tor_id_ref(ligids_[i], id_tor_type, j);
				bool const fail( pose->conformation().get_torsion_angle_atom_ids(tor_id_curr, atid1, atid2, atid3, atid4) );
				if ( fail ) {
					if ( ( r == 1 ) &&
							( ( j == 1 && pose->fold_tree().is_cutpoint( i-1 ) ) ||
							( j >= n_torsions-1 && pose->fold_tree().is_cutpoint( i ) ) ) ) continue;
					if ( TR.Debug.visible() ) TR.Debug << " missed torsion: " << tor_id_curr << std::endl;
					continue;
				}
				torsion_type.at1 = ligres.atom_type_index( atid1.atomno() );
				torsion_type.at2 = ligres.atom_type_index( atid2.atomno() );
				torsion_type.at3 = ligres.atom_type_index( atid3.atomno() );
				torsion_type.at4 = ligres.atom_type_index( atid4.atomno() );
				torsion_type.bn = ligres.type().bond_type( atid2.atomno(), atid3.atomno() );
				torsion_type.br = ligres.type().bond_ringness( atid2.atomno(), atid3.atomno() );
				ligandchi_types_.push_back( torsion_type );
				ligandchis_.push_back( r == 1 ? ligres.mainchain_torsion(j): ligres.chi(j) );
				ligandtorsionids_.push_back(tor_id_ref);
			}
		}
	}

	// internal ring torsions
	if ( !freeze_ligand_ && sample_ring_conformers_ && ligids_.size() == 1 ) {
		// *** for now assume this only used in single-ligand case
		//   it might affect some non-canonicals
		core::chemical::ResidueType const &ligrt = pose->residue(1).type();
		core::Size nnu = 0;
		core::Size nligring( ligrt.n_rings() );
		for ( core::Size j=1; j<=nligring; ++j ) {
			nnu += ligrt.ring_atoms( j ).size() - 1;
		}
		ligandnus_.resize(nnu);
		ligandtaus_.resize(nnu+nligring);
		core::Size offset = 0;
		for ( core::Size j=1; j<=nligring; ++j ) {
			nnu = ligrt.ring_atoms( j ).size() - 1;
			for ( core::Size k=1; k<=nnu; ++k ) {
				ligandnus_[offset+k] = pose->torsion( core::id::TorsionID( 1, core::id::NU, offset + k ) );
				core::id::AtomID a1( ligrt.nu_atoms( offset + k )[ 1 ], 1 );
				core::id::AtomID a2( ligrt.nu_atoms( offset + k )[ 2 ], 1 );
				core::id::AtomID a3( ligrt.nu_atoms( offset + k )[ 3 ], 1 );
				ligandtaus_[offset+(j-1)+k] = numeric::conversions::degrees( pose->conformation().bond_angle( a1, a2, a3 ) );
			}
			core::id::AtomID a1( ligrt.nu_atoms( offset + nnu )[ 2 ], 1 );
			core::id::AtomID a2( ligrt.nu_atoms( offset + nnu )[ 3 ], 1 );
			core::id::AtomID a3( ligrt.nu_atoms( offset + nnu )[ 4 ], 1 );
			ligandtaus_[offset+j+nnu] = numeric::conversions::degrees( pose->conformation().bond_angle( a1, a2, a3 ) );
			offset += nnu;
		}
	}

	// copy ligandxyz
	if ( !freeze_ligand_ || rg_ == 0.0 ) {
		ligandxyz_.clear();
		core::Vector com( 0.0 );
		for ( core::Size i=1;  i<= ligids_.size(); ++i ) {
			core::conformation::Residue const &lig( pose->residue(i) );
			for ( core::Size iatm = 1; iatm <= lig.nheavyatoms(); ++iatm ) {
				ligandxyz_.push_back(lig.xyz( iatm ));
				com += lig.xyz( iatm );
			}
		}
		com /= ligandxyz_.size();

		// compute radius of gyration
		rg_ = 0.0;
		for ( auto atm : ligandxyz_ ) {
			rg_ += com.distance_squared( atm );
		}
		rg_ /= ligandxyz_.size();
		rg_ = std::sqrt(rg_);
	}

	ligandxyz_synced_ = true;
}

// gene -> pose
// generate a pose based on this conformation
void
LigandConformer::to_pose( core::pose::PoseOP pose ) const {
	if ( pose->total_residue() == 0 ) *pose = *ref_pose_;

	// jump
	numeric::Quaternion< core::Real > Q( rb_[1], rb_[2], rb_[3], rb_[4] );
	numeric::xyzVector< core::Real >  T( rb_[5], rb_[6], rb_[7] );
	numeric::xyzMatrix< core::Real >  R;
	numeric::quat2R( Q, R );

	if ( jumpid_ != 0 ) {
		core::kinematics::Jump ligjump = pose->jump( jumpid_ );
		ligjump.set_translation( T );
		ligjump.set_rotation( R );
		pose->set_jump( jumpid_, ligjump );
	}

	// movable sidechains
	for ( core::Size i=1; i<=movingscs_.size(); ++i ) {
		core::Size resid_i = movingscs_[i];

		if ( pose->residue_type(resid_i).name() != proteinrestypes_[i]->name() ) {   // is there a better way to compare restypes?
			core::conformation::Residue newres(proteinrestypes_[i], true);
			pose->replace_residue( resid_i, newres, true );
		}
		for ( core::Size j=1; j<=proteinchis_[i].size(); ++j ) {
			pose->set_chi( j, resid_i, proteinchis_[i][j] );
		}
	}

	// internal torsions
	for ( core::Size i=1; i <= ligandtorsionids_.size(); ++i ) {
		pose->conformation().set_torsion(ligandtorsionids_[i], ligandchis_[i]);
	}

	// internal ring torsions
	if ( sample_ring_conformers_ && ligids_.size() == 1 ) {
		core::chemical::ResidueType const &ligrt = pose->residue(ligids_[1]).type();
		core::Size nligring( ligrt.n_rings() );
		core::Size offset = 0;
		for ( core::Size j=1; j<=nligring; ++j ) {
			core::Size nnu = ligrt.ring_atoms( j ).size() - 1;
			for ( core::Size k=1; k<=nnu; ++k ) {
				pose->set_torsion( core::id::TorsionID( ligids_[1], core::id::NU, offset + k ), ligandnus_[offset+k] );
				core::id::AtomID a1( ligrt.nu_atoms( offset + k )[ 1 ], ligids_[1] );
				core::id::AtomID a2( ligrt.nu_atoms( offset + k )[ 2 ], ligids_[1] );
				core::id::AtomID a3( ligrt.nu_atoms( offset + k )[ 3 ], ligids_[1] );
				pose->conformation().set_bond_angle( a1, a2, a3, numeric::conversions::radians(ligandtaus_[offset+(j-1)+k]) );
			}
			core::id::AtomID a1( ligrt.nu_atoms( offset + nnu )[ 2 ], ligids_[1] );
			core::id::AtomID a2( ligrt.nu_atoms( offset + nnu )[ 3 ], ligids_[1] );
			core::id::AtomID a3( ligrt.nu_atoms( offset + nnu )[ 4 ], ligids_[1] );
			pose->conformation().set_bond_angle( a1, a2, a3, numeric::conversions::radians(ligandtaus_[offset+j+nnu]) );
			offset += nnu;
		}
	}
}


///
/// create a reduced pose representation (to be used in minimization)
void
LigandConformer::to_minipose( core::pose::PoseOP pose, LigandConformer &minilig ) const {
	core::pose::PoseOP fullpose( new core::pose::Pose );
	to_pose( fullpose );

	minilig = *this;
	pose->clear();

	// sidechains
	for ( core::Size i=1; i<=movingscs_.size(); ++i ) {
		if ( i == 1 ||  (movingscs_[i-1] == movingscs_[i]-1 && !fullpose->residue( movingscs_[i-1] ).is_upper_terminus() ) ) {
			pose->append_residue_by_bond( fullpose->residue( movingscs_[i] ) );
		} else {
			pose->append_residue_by_jump( fullpose->residue( movingscs_[i] ), pose->total_residue() );
		}
		minilig.movingscs_[i] = pose->total_residue();
	}

	// ligand
	minilig.ligids_.clear();
	for ( core::Size i=1; i<=ligids_.size(); ++i ) {
		if ( i==1 ) {
			pose->append_residue_by_jump( fullpose->residue( ligids_[i] ), pose->total_residue() );
			minilig.jumpid_ = pose->fold_tree().num_jump();
		} else {
			pose->append_residue_by_bond( fullpose->residue( ligids_[i] ) );
		}
		minilig.ligids_.push_back(pose->total_residue());
	}

	// root ligand on virtual if no residues to anchor jump
	// (so we have a jump to minimize)
	if ( minilig.jumpid_ == 0 ) {
		core::pose::addVirtualResAsRoot(*pose);
		minilig.jumpid_ = 1;
	}
}

///
/// update internal information from the reduced pose representation
void
LigandConformer::update_conf_from_minipose( core::pose::PoseCOP pose ) {
	core::pose::PoseOP fullpose( new core::pose::Pose );
	to_pose( fullpose );

	int resctr = 1;
	for ( core::Size i=1; i<=movingscs_.size(); ++i ) {
		fullpose->replace_residue( movingscs_[i], pose->residue(resctr++), false );
	}
	for ( core::Size i=1; i<=ligids_.size(); ++i ) {
		fullpose->replace_residue( ligids_[i], pose->residue(resctr++), false );
	}

	update_conf( fullpose );
}

std::string
LigandConformer::to_string() const {
	using namespace ObjexxFCL::format;
	std::ostringstream oss;
	oss << "[[ rt: ";
	for ( core::Size i=1; i<=4; ++i ) {
		oss << F( 6, 3, rb_[i]) << " ";
	}
	for ( core::Size i=5; i<=7; ++i ) {
		oss << F( 6, 1, rb_[i]) << " ";
	}
	oss << "c: ";
	for ( core::Size i=1; i<=ligandchis_.size(); ++i ) {
		oss << F( 6, 1, ligandchis_[i]) << " ";
	}
	oss << " ]]";
	return oss.str();
}



void
LigandConformer::dump_pose( std::string pdbname ) const {
	core::pose::PoseOP pose( new core::pose::Pose() );
	to_pose( pose );
	pose->dump_pdb( pdbname );
}


// generate only the ligand residue based on this conformation
core::conformation::Residue
LigandConformer::ligand_residue( core::Size ires ) const {
	core::pose::PoseOP pose (new core::pose::Pose( *ref_pose_ ));
	to_pose( pose );
	return (pose->residue(ligids_[ires]));
}

// generate only the protein residue based on this conformation
core::conformation::Residue
LigandConformer::protein_residue( core::Size ires ) const {
	core::pose::PoseOP pose (new core::pose::Pose( *ref_pose_ ));
	to_pose( pose );
	return (pose->residue( ires ));
}

// generate a pose containing the _input_ receptor structure
core::pose::PoseOP
LigandConformer::receptor( ) const {
	core::pose::PoseOP pose (new core::pose::Pose( *ref_pose_ ));
	for ( auto ligid = ligids_.rbegin(); ligid != ligids_.rend(); ++ligid ) {
		pose->delete_residue_slow( *ligid );
	}
	return (pose);
}


// const access for distance calculation;
// ref_pose_ should be synced with rb,ligchi,and so on through update_conf
utility::vector1< core::Vector > const &
LigandConformer::ligand_xyz( ) {
	if ( !ligandxyz_synced_ ) {
		TR.Debug << "WARNING! Expensive ligand resynch" << std::endl; // it's not that expensive
		core::pose::PoseOP pose (new core::pose::Pose( *ref_pose_ ));
		to_pose( pose );

		// copy ligandxyz
		ligandxyz_.clear();
		core::Vector com( 0.0 );
		for ( core::Size i : ligids_ ) {
			core::conformation::Residue const &lig( pose->residue(i) );
			for ( core::Size iatm = 1; iatm <= lig.nheavyatoms(); ++iatm ) {
				ligandxyz_.push_back(lig.xyz( iatm ));
				com += lig.xyz( iatm );
			}
		}
		com /= ligandxyz_.size();

		// compute radius of gyration
		rg_ = 0.0;
		for ( auto atm : ligandxyz_ ) {
			rg_ += com.distance_squared( atm );
		}
		rg_ /= ligandxyz_.size();
		rg_ = std::sqrt(rg_);

		ligandxyz_synced_ = true;
	}
	return ligandxyz_;
}

// the initial perturbation randomizes ligand conf
void
LigandConformer::randomize( core::Real transmax ) {
	ligand_xyz(); // force synch

	numeric::Quaternion< core::Real > Q;
	numeric::xyzMatrix< core::Real > R=numeric::random::random_rotation();
	numeric::R2quat( R, Q );
	rb_[1] = Q.w();
	rb_[2] = Q.x();
	rb_[3] = Q.y();
	rb_[4] = Q.z();

	core::Vector Taxis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
	core::Real Tlen = numeric::random::rg().uniform()*transmax;
	for ( core::Size k = 5; k <= 7; ++k ) {
		rb_[k] += Tlen*Taxis[k-5];
	}

	// ligand chis
	core::scoring::RamaPrePro const& rama_prepro( core::scoring::ScoringManager::get_instance()->get_RamaPrePro() );
	core::Size nligchi = ligandchis_.size();

	if ( !freeze_ligand_ ) {
		for ( core::Size j=1; j<=nligchi; ++j ) {
			core::id::TorsionID const& torid=ligandtorsionids_[j];
			core::Real angle_j = 360.0 * numeric::random::rg().uniform();
			core::Size ligid = torid.rsd();
			core::chemical::ResidueTypeCOP restype( ligid_restype_map_[ligid] );
			if ( torid.type() == core::id::BB && torid.torsion() == core::id::phi_torsion && ligid_restype_map_.count(ligid+1) > 0 ) {
				utility::vector1 < core::Real > torsions;
				rama_prepro.random_mainchain_torsions(ref_pose_->conformation(), restype, ligid_restype_map_[ligid+1], torsions);
				ligandchis_[j] = torsions[1];
				j++;
				ligandchis_[j] = torsions[2];
				continue;
			}
			// a quick fix for omega torsion, set them to trans
			if ( torid.type() == core::id::BB && torid.torsion() == core::id::omega_torsion ) {
				ligandchis_[j] = 180.0;
				continue;
			}
			ligandchis_[j] = angle_j;
		}
	}

	ligandxyz_synced_ = false;

	// ligand nus
	// load ring confs
	if ( !freeze_ligand_ && sample_ring_conformers_ && ligids_.size() == 1 ) {
		core::chemical::ResidueType const &ligrt = ref_pose_->residue(ligids_[1]).type();
		core::Size const n_rings( ligrt.n_rings() );
		core::Size offset = 0;
		for ( core::Size j=1; j<=n_rings; ++j ) {
			// pick a random conformation
			utility::vector1< core::chemical::rings::RingConformer > const & ringconfs
				= ligrt.ring_conformer_set( j )->get_all_nondegenerate_conformers();
			core::Size nringconfs = ringconfs.size();
			core::Size pickedconf = numeric::random::random_range( 1, nringconfs );
			TR << "RING " << j << " picked conf " << pickedconf << " of " << nringconfs << std::endl;
			core::Size nnu = ligrt.ring_atoms( j ).size() - 1;
			for ( core::Size k=1; k<=nnu; ++k ) {
				ligandnus_[offset+k] = ringconfs[pickedconf].nu_angles[k];
				ligandtaus_[offset+(j-1)+k] = ringconfs[pickedconf].tau_angles[k];
			}
			ligandtaus_[offset+j+nnu] = ringconfs[pickedconf].tau_angles[nnu+1];
			offset += nnu;
		}
	}
}

void
LigandConformer::sample_conformation(
	core::Real transmax,
	TorsionSamplerCOP const & sampler )
{
	ligand_xyz( ); // force synch

	numeric::Quaternion< core::Real > Q;
	numeric::xyzMatrix< core::Real > R=numeric::random::random_rotation();
	numeric::R2quat( R, Q );
	rb_[1] = Q.w();
	rb_[2] = Q.x();
	rb_[3] = Q.y();
	rb_[4] = Q.z();

	core::Vector Taxis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
	core::Real Tlen = numeric::random::rg().uniform()*transmax;
	for ( core::Size k = 5; k <= 7; ++k ) {
		rb_[k] += Tlen*Taxis[k-5];
	}

	// ligand chis
	core::Size nligchi = ligandchis_.size();
	for ( core::Size j=1; j<=nligchi; ++j ) {
		TorsionType const& ttype = ligandchi_types_[j];
		core::Real degree = sampler->sample(ttype.bn, ttype.br, ttype.at1,
			ttype.at2,ttype.at3,ttype.at4);
		ligandchis_[j] = degree;
	}

	ligandxyz_synced_ = false;

	// ligand nus
	// load ring confs
	if ( sample_ring_conformers_ && ligids_.size() == 1 ) {
		core::chemical::ResidueType const &ligrt = ref_pose_->residue(ligids_[1]).type();
		core::Size const n_rings( ligrt.n_rings() );
		core::Size offset = 0;
		for ( core::Size j=1; j<=n_rings; ++j ) {
			// pick a random conformation
			utility::vector1< core::chemical::rings::RingConformer > const & ringconfs
				= ligrt.ring_conformer_set( j )->get_all_nondegenerate_conformers();
			core::Size nringconfs = ringconfs.size();
			core::Size pickedconf = numeric::random::random_range( 1, nringconfs );
			TR << "RING " << j << " picked conf " << pickedconf << " of " << nringconfs << std::endl;
			core::Size nnu = ligrt.ring_atoms( j ).size() - 1;
			for ( core::Size k=1; k<=nnu; ++k ) {
				ligandnus_[offset+k] = ringconfs[pickedconf].nu_angles[k];
				ligandtaus_[offset+(j-1)+k] = ringconfs[pickedconf].tau_angles[k];
			}
			ligandtaus_[offset+j+nnu] = ringconfs[pickedconf].tau_angles[nnu+1];
			offset += nnu;
		}
	}
}

void
LigandConformer::assign_ligand_trans( core::Vector transv ) {
	for ( core::Size k = 5; k <= 7; ++k ) {
		rb_[k] = transv[k-5];
	}
}


void
LigandConformer::superimpose_to_ref_pose( utility::vector1< core::id::AtomID > const & ids ) {
	// init fullpose
	core::pose::PoseOP fullpose ( new core::pose::Pose );
	to_pose( fullpose );

	// init pose_ligand
	core::pose::PoseOP pose_ligand ( new core::pose::Pose );
	for ( core::Size i=1; i<=ligids_.size(); ++i ) {
		if ( i==1 ) {
			pose_ligand->append_residue_by_jump( fullpose->residue( ligids_[i] ), pose_ligand->total_residue() );
		} else {
			pose_ligand->append_residue_by_bond( fullpose->residue( ligids_[i] ) );
		}
	}
	core::pose::addVirtualResAsRoot( *pose_ligand );

	// init atom_map
	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map(atom_map, *pose_ligand, core::id::AtomID::BOGUS_ATOM_ID());
	for ( core::id::AtomID const & id : ids ) {
		atom_map[id] = id;
	}

	// align pose_ligands
	core::scoring::superimpose_pose( *pose_ligand, *ref_pose_ligand_, atom_map );

	// update fullpose
	for ( core::Size i=1; i<=ligids_.size(); ++i ) {
		fullpose->replace_residue( ligids_.at(i), pose_ligand->residue(i), false );
	}
	update_conf( fullpose );
}

// mutation
LigandConformer
mutate(LigandConformer const &l ) {
	using namespace ObjexxFCL::format;

	LigandConformer retval( l );
	std::string tag;

	// rigid body
	bool randomize_rt = (numeric::random::rg().uniform() < l.rtmutationRate_ );

	tag += " ";
	if ( randomize_rt ) {
		numeric::Quaternion< core::Real > Q;
		if ( l.rotmutWidth_ > 180.0 ) { // complete randomization
			numeric::xyzMatrix< core::Real > R=numeric::random::random_rotation();
			numeric::R2quat( R, Q );
		} else {  // perturbation
			core::Vector Raxis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
			core::Real angle = numeric::NumericTraits<core::Real>::deg2rad() * l.rotmutWidth_ * numeric::random::rg().gaussian();
			core::Real ca = cos(angle), sa=sin(angle);
			numeric::Quaternion< core::Real > Q1( ca, sa*Raxis[0], sa*Raxis[1], sa*Raxis[2]);
			Q = Q1*numeric::Quaternion< core::Real >(l.rb_[1],l.rb_[2],l.rb_[3],l.rb_[4]);
		}
		retval.rb_[1] = Q.w(); retval.rb_[2] = Q.x(); retval.rb_[3] = Q.y(); retval.rb_[4] = Q.z();

		core::Vector Taxis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
		core::Real len = numeric::random::rg().uniform()*l.transmutWidth_;
		for ( core::Size k = 5; k <= 7; ++k ) {
			retval.rb_[k] += len*Taxis[k-5];
		}
		tag += "rt ";
	} else {
		// ligand chis
		core::Size nligchi = l.ligandchis_.size();
		core::Size ntors_changed(0);
		for ( core::Size j=1; j<=nligchi; ++j ) {
			if ( numeric::random::rg().uniform() < l.torsmutationRate_ ) {
				core::Real angle_j;
				if ( l.ligchimutWidth_ > 180.0 ) { // complete randomization
					angle_j = 360.0 * numeric::random::rg().uniform();
				} else {
					core::Real angleDel = l.ligchimutWidth_*(1.0 - 2.0*numeric::random::rg().uniform());
					angle_j = l.ligandchis_[j] + angleDel;
				}
				retval.ligandchis_[j] = angle_j;
				tag += "c"+utility::to_string(j)+" ";
				ntors_changed++;
			}
		}

		// make sure at least one torsion is mutated
		if ( ntors_changed == 0 && !randomize_rt && nligchi > 0 ) {
			core::Size j( numeric::random::rg().random_range(1,nligchi) );
			core::Real angle_j;
			if ( l.ligchimutWidth_ > 180.0 ) { // complete randomization
				angle_j = 360.0 * numeric::random::rg().uniform();
			} else {
				core::Real angleDel = l.ligchimutWidth_*(1.0 - 2.0*numeric::random::rg().uniform());
				angle_j = l.ligandchis_[j] + angleDel;
			}
			retval.ligandchis_[j] = angle_j;
			tag += "c"+utility::to_string(j)+" ";
		}

		// ligand nus
		//   mutate entire ring conf. with rate=torsmutationRate_
		if ( l.sample_ring_conformers() && l.ligids_.size() == 1 ) {
			core::chemical::ResidueType const &ligrt = l.ref_pose_->residue(l.ligids_[1]).type();
			core::Size const n_rings( ligrt.n_rings() );
			core::Size offset = 0;
			for ( core::Size j=1; j<=n_rings; ++j ) {
				utility::vector1< core::chemical::rings::RingConformer > const & ringconfs
					= ligrt.ring_conformer_set( j )->get_all_nondegenerate_conformers();
				core::Size nringconfs = ringconfs.size();
				core::Size pickedconf = numeric::random::random_range( 1, nringconfs );
				core::Size nnu = ligrt.ring_atoms( j ).size() - 1;
				if ( numeric::random::rg().uniform() < l.torsmutationRate_ ) {
					// pick a random conformation
					for ( core::Size k=1; k<=nnu; ++k ) {
						retval.ligandnus_[offset+k] = ringconfs[pickedconf].nu_angles[k];
						retval.ligandtaus_[offset+(j-1)+k] = ringconfs[pickedconf].tau_angles[k];
					}
					retval.ligandtaus_[offset+j+nnu] = ringconfs[pickedconf].tau_angles[nnu+1];
					tag += "r"+utility::to_string(j)+" ";
				}
				offset += nnu;
			}
		}
	}

	retval.set_generation_tag( tag );
	retval.ligandxyz_synced_ = false;

	return retval;
}

// mutation
LigandConformer
mutate_ft( LigandConformer const &l, bool single_mutation ) {
	using namespace ObjexxFCL::format;

	LigandConformer retval( l );
	std::string tag;

	// rigid body
	bool randomize_rt = (numeric::random::rg().uniform() < l.rtmutationRate_ );

	tag += " ";
	if ( randomize_rt ) {
		numeric::Quaternion< core::Real > Q;
		if ( l.rotmutWidth_ > 180.0 ) { // complete randomization
			numeric::xyzMatrix< core::Real > R=numeric::random::random_rotation();
			numeric::R2quat( R, Q );
		} else {  // perturbation
			core::Vector Raxis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
			core::Real angle = numeric::NumericTraits<core::Real>::deg2rad() * l.rotmutWidth_ * numeric::random::rg().gaussian();
			core::Real ca = cos(angle), sa=sin(angle);
			numeric::Quaternion< core::Real > Q1( ca, sa*Raxis[0], sa*Raxis[1], sa*Raxis[2]);
			Q = Q1*numeric::Quaternion< core::Real >(l.rb_[1],l.rb_[2],l.rb_[3],l.rb_[4]);
		}
		retval.rb_[1] = Q.w(); retval.rb_[2] = Q.x(); retval.rb_[3] = Q.y(); retval.rb_[4] = Q.z();

		core::Vector Taxis( numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() ) );
		core::Real len = numeric::random::rg().uniform()*l.transmutWidth_;
		for ( core::Size k = 5; k <= 7; ++k ) {
			retval.rb_[k] += len*Taxis[k-5];
		}
		tag += "rt ";
	} else {
		// ligand chis
		core::Size nligchi = l.ligandchis_.size();

		// make sure at least one torsion is mutated
		if ( nligchi > 0 ) {
			if ( single_mutation ) {
				core::Size j( numeric::random::rg().random_range(1,nligchi) );
				core::Real angle_j;
				if ( l.ligchimutWidth_ > 180.0 ) { // complete randomization
					angle_j = 360.0 * numeric::random::rg().uniform();
				} else {
					core::Real angleDel = l.ligchimutWidth_*(1.0 - 2.0*numeric::random::rg().uniform());
					angle_j = l.ligandchis_[j] + angleDel;
				}
				retval.ligandchis_[j] = angle_j;
				tag += "c"+utility::to_string(j)+" ";
			}
		

			else {
				core::Size nligchi = l.ligandchis_.size();
				core::Size ichi_to_mutate( numeric::random::rg().random_range(1,nligchi) );
				utility::vector1< core::Size > const &chis_down = l.ligandchi_downstream_[ichi_to_mutate];
				
				core::Real angleDel = l.ligchimutWidth_*(1.0 - 2.0*numeric::random::rg().uniform());
				retval.ligandchis_[ichi_to_mutate] = l.ligandchis_[ichi_to_mutate] + angleDel;
				tag += "c"+utility::to_string(ichi_to_mutate)+" ";
				for ( core::Size j = 1; j <= chis_down.size(); ++j ) {
					if ( numeric::random::rg().uniform() < 0.5 ) {
						core::Size jchi = chis_down[j];
						tag += "c"+utility::to_string(jchi)+" ";
						angleDel = /*l.ligchimutWidth_*/180*(1.0 - 2.0*numeric::random::rg().uniform());
						retval.ligandchis_[jchi] = l.ligandchis_[jchi] + angleDel;
					}
				}
			}
		}
	}

	retval.set_generation_tag( tag );
	retval.ligandxyz_synced_ = false;

	return retval;
}

// crossover
LigandConformer
crossover(LigandConformer const &l1, LigandConformer const &l2) {
	using namespace ObjexxFCL::format;

	LigandConformer retval;
	std::string tag;

	// rigid body
	//   a) take rb as a single unit
	//   b) take sidechains from same gene we inherit rb (will be repacked though)
	bool takeRBfromIn1 = (numeric::random::rg().uniform() < 0.5 );
	if ( takeRBfromIn1 ) {
		retval = l1;
		tag += " rt:1 ";
	} else {
		retval = l2;
		tag += " rt:2 ";
	}
	// TODO: superimpose random subset of atoms instead of "stealing" jump

	// torsional
	core::Size nligchi = l1.ligandchis_.size();
	core::Size ntors_changed(0);
	tag += "c:";
	for ( core::Size j=1; j<=nligchi; ++j ) {
		bool takeChiJfromIn1 = (numeric::random::rg().uniform() < 0.5 );
		if ( takeChiJfromIn1 ) {
			retval.ligandchis_[j] = l1.ligandchis_[j];
			tag += "1";
		} else {
			ntors_changed++;
			retval.ligandchis_[j] = l2.ligandchis_[j];
			tag += "2";
		}
	}
	tag += " ";

	// make sure at least one torsion is mutated
	if ( ntors_changed == 0 && nligchi > 0 ) {
		core::Size j = core::Size(nligchi*numeric::random::rg().uniform())+1;
		retval.ligandchis_[j] = l2.ligandchis_[j];
		tag += "x"+utility::to_string(j)+" ";
	}

	// ligand nus
	//   take ring confs as a single unit
	if ( l1.sample_ring_conformers() && l1.ligids_.size() == 1 ) {
		core::chemical::ResidueType const &ligrt = l1.ref_pose_->residue(l1.ligids_[1]).type();
		core::Size const n_rings( ligrt.n_rings() );
		core::Size offset = 0;
		if ( n_rings>0 ) tag += "r:";
		for ( core::Size j=1; j<=n_rings; ++j ) {
			//utility::vector1< core::chemical::rings::RingConformer > const & ringconfs
			//= ligrt.ring_conformer_set( j )->get_all_nondegenerate_conformers();
			core::Size nnu = ligrt.ring_atoms( j ).size() - 1;

			bool takeNuJfromIn1 = (numeric::random::rg().uniform() < 0.5 );
			if ( takeNuJfromIn1 ) {
				for ( core::Size k=1; k<=nnu; ++k ) {
					retval.ligandnus_[offset+k] = l1.ligandnus_[offset+k];
					retval.ligandtaus_[offset+(j-1)+k] = l1.ligandtaus_[offset+(j-1)+k];
				}
				retval.ligandtaus_[offset+j+nnu] = l1.ligandtaus_[offset+j+nnu];
				tag += "1";
			} else {
				for ( core::Size k=1; k<=nnu; ++k ) {
					retval.ligandnus_[offset+k] = l2.ligandnus_[offset+k];
					retval.ligandtaus_[offset+(j-1)+k] = l2.ligandtaus_[offset+(j-1)+k];
				}
				retval.ligandtaus_[offset+j+nnu] = l2.ligandtaus_[offset+j+nnu];
				tag += "2";
			}
			offset += nnu;
		}
	}

	retval.set_generation_tag( tag );

	retval.ligandxyz_synced_ = false;

	return retval;
}

LigandConformer
crossover_ft(LigandConformer const &l1, LigandConformer const &l2 ){
        // pick (randomly) one to be "parent"
        //fd let's give every pose a shot to be parent
        bool l1_is_parent = true; // (numeric::random::rg().uniform() <= 0.5 );
        LigandConformer const &l_parent = l1_is_parent? l1 : l2;
        LigandConformer const &l_child = l1_is_parent? l2 : l1;

        LigandConformer retval(l_parent);
        std::string tag;
        if ( l1_is_parent ) tag = "parent1 ";
        else tag = "parent2 ";

        // pick ONE chi and set all downstream to l_child
        core::Size nligchi = l1.ligandchis_.size();
        core::Size ichi_to_cross( numeric::random::rg().random_range(1,nligchi) );
        utility::vector1< core::Size > const &chis_down = l1.ligandchi_downstream_[ichi_to_cross];
        retval.ligandchis_[ichi_to_cross] = l_child.ligandchis_[ichi_to_cross];
        tag += "flip:c"+utility::to_string(ichi_to_cross)+" ";
        for ( core::Size j = 1; j <= chis_down.size(); ++j ) {
                core::Size jchi = chis_down[j];
                tag += "down:c"+utility::to_string(jchi)+" ";
                retval.ligandchis_[jchi] = l_child.ligandchis_[jchi];
        }

        retval.set_generation_tag( tag );
        retval.ligandxyz_synced_ = false;

        return retval;
}

core::Real
distance_slow( LigandConformer &gene1, LigandConformer &gene2 ){
	if ( gene1.ligids_.size() > 1 || gene2.ligids_.size() > 1 ) {
		return (distance_fast(gene1,gene2));
	}

	return core::scoring::automorphic_rmsd(
		gene1.ligand_residue(1),
		gene2.ligand_residue(1),
		false
	);
}

core::Real
distance_fast( LigandConformer &gene1, LigandConformer &gene2 ) {
	core::Real d( 0.0 );
	core::Size n( 0 );
	utility::vector1< core::Vector >const &xyz1 = gene1.ligand_xyz();
	utility::vector1< core::Vector >const &xyz2 = gene2.ligand_xyz();
	for ( core::Size iatm = 1; iatm <= xyz1.size(); ++iatm ) {
		d += xyz1[iatm].distance_squared( xyz2[iatm] );
		n++;
	}

	return std::sqrt(d/n);
}

// return distance based on chi difference
// FD: NOTE THIS WILL NEED TO BE MODIFIED FOR DESIGN! (e.g. if nchi(gene1) != nchi(gene2) at some position)
std::pair< core::Real, core::Real >
distance_internal( LigandConformer const &gene1, LigandConformer const &gene2 ){
	std::pair< core::Real, core::Real > d;

	numeric::Quaternion< core::Real > q1 = gene1.quat();
	numeric::Quaternion< core::Real > q2 = gene2.quat();
	numeric::xyzVector< core::Real > t1 = gene1.trans();
	numeric::xyzVector< core::Real > t2 = gene2.trans();

	// more weight on rotation as ligand gets bigger
	// perhaps better estimation using radius of generation?
	//core::Real const r = gene1.ligand_residue().nbr_radius();
	core::Real const r = 0.5*(gene1.ligand_rg() + gene2.ligand_rg());

	numeric::Quaternion< core::Real > qdiff = q1.invert()*q2;
	d.first = t1.distance(t2) + r*qdiff.angle(); // angle is radian

	core::Size nligchi = gene1.get_ligandchis().size();
	core::Size ndof( nligchi );
	for ( core::Size ires = 1; ires <= gene1.moving_scs().size(); ++ires ) {
		ndof += gene1.get_protein_chis(ires).size();
	}
	//chidiff.resize(ndof);
	utility::vector1< core::Real > chidiff( ndof, 0.0 );

	for ( core::Size ichi = 1; ichi <= gene1.get_ligandchis().size(); ++ichi ) {
		core::Real diff( gene1.get_ligandchi(ichi) - gene2.get_ligandchi(ichi) );
		chidiff[ichi] = diff;
	}

	for ( core::Size ires = 1; ires <= gene1.moving_scs().size(); ++ires ) {
		utility::vector1< core::Real > const &reschi1 = gene1.get_protein_chis(ires);
		utility::vector1< core::Real > const &reschi2 = gene2.get_protein_chis(ires);
		for ( core::Size ichi = 1; ichi <= reschi1.size(); ++ichi ) {
			core::Real diff( reschi1[ichi] - reschi2[ichi] );
			chidiff[nligchi+ichi] = diff;
		}
	}

	for ( core::Size idof = 1; idof <= ndof; ++idof ) {
		d.second += std::abs(chidiff[idof]); // use Hamming instead of MSD or RMSD
	}

	return d;
}


}
}
}
