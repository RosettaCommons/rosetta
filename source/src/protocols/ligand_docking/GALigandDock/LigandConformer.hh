// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/LigandConformer.hh
///
/// @brief  Compactly represent a docked pose by storing jump + torsion dofs + pose sidechain dofs
/// @author Hahnbeom Park and Frank DiMaio

#ifndef INCLUDED_protocols_ligand_docking_GALigandDock_LigandConformer_hh
#define INCLUDED_protocols_ligand_docking_GALigandDock_LigandConformer_hh

#include <protocols/ligand_docking/GALigandDock/LigandConformer.fwd.hh>

#include <utility/VirtualBase.hh>
#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/id/TorsionID.hh>
#include <numeric/Quaternion.hh>
#include <utility/vector1.hh>
#include <protocols/ligand_docking/GALigandDock/TorsionSampler.fwd.hh>
#include <core/id/AtomID.hh>
#include <map>
#include <ObjexxFCL/FArray3D.hh>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

/// @brief
/// Gene representation of ligand & flexible sidechains in receptor
/// @details
/// Gene is preresented by rigid body (rb_), ligandchis, and receptorchis
/// Also has functions to tranform back and forth to pose object
/// Uses friend functions to perform mutation / crossovers with others within gene representation

struct TorsionType{
	core::chemical::BondName bn;
	core::chemical::BondRingness br;
	core::Size at1, at2, at3, at4;
};

class PointScoreComparator{
public:
	bool operator()(std::pair< numeric::xyzVector<core::Real>, core::Real > Pair1, std::pair< numeric::xyzVector<core::Real>, core::Real > Pair2){
		return (Pair1.second < Pair2.second);
	}

};

class LigandConformer : public utility::VirtualBase {
public:
	friend LigandConformer
	mutate(LigandConformer const &l );

	// alternative mutation for many torsions
	friend LigandConformer
	mutate_ft(LigandConformer const &l, bool single_mutation );

	// crossover
	friend LigandConformer
	crossover(LigandConformer const &l1, LigandConformer const &l2);

	// alternative mutation for many torsions
	friend LigandConformer
	crossover_ft(LigandConformer const &l1, LigandConformer const &l2);

	// distance without equivalent atom substitution
	friend core::Real
	distance_fast( LigandConformer &gene1, LigandConformer &gene2 );

	// distance _with_ equivalent atom substitution
	friend core::Real
	distance_slow( LigandConformer &gene1, LigandConformer &gene2 );

	// internal coordinate distance
	friend std::pair< core::Real, core::Real >
	distance_internal( LigandConformer const &gene1, LigandConformer const &gene2 );


public:
	LigandConformer();
	~LigandConformer() override;

	LigandConformer(
		core::pose::PoseCOP pose,
		utility::vector1 <core::Size > const &ligids,
		utility::vector1< core::Size > movingscs,
		bool freeze_ligand_backbone=false,
		bool freeze_ligand=false

	);

	// common pathway for parameter initialization
	void
	init_params();

	// set the reference pose and initial state
	void
	initialize(
		core::pose::PoseCOP pose,
		utility::vector1 <core::Size > const &ligids,
		utility::vector1< core::Size > movingscs,
		bool freeze_ligand_backbone,
		bool freeze_ligand
	);

	bool
	defined() const { return (rb_.size() > 0); }

	// update internal representation based on this conformation
	void
	update_conf( core::pose::PoseCOP pose );

	// update internal representation of the ligand based on this ligand conformation
	void
	update_ligand_conf( core::pose::PoseCOP pose );

	// generate a pose based on this conformation
	void
	to_pose( core::pose::PoseOP pose) const;

	// generate a pose based on this conformation
	void
	to_minipose( core::pose::PoseOP pose, LigandConformer &minilig) const;

	// update internal information from the reduced pose representation
	void
	update_conf_from_minipose( core::pose::PoseCOP pose );

	numeric::Quaternion< core::Real > const
	quat() const { return numeric::Quaternion< core::Real >( rb_[1], rb_[2], rb_[3], rb_[4] ); }

	numeric::xyzVector< core::Real > const
	trans() const { return numeric::xyzVector< core::Real >( rb_[5], rb_[6], rb_[7] ); }

	// generate only the ligand residue based on this conformation
	utility::vector1< core::Vector > const &
	ligand_xyz();

	core::Real
	ligand_rg() const { return rg_; };

	// generate only the ligand residue based on this conformation
	core::conformation::Residue
	ligand_residue( core::Size ires ) const;

	// generate only the protein residue based on this conformation
	core::conformation::Residue
	protein_residue( core::Size ires ) const;

	// generate only the protein residue based on this conformation
	core::pose::PoseOP
	receptor(  ) const;

	// get the id of the ligand
	utility::vector1< core::Size >
	ligand_ids() const { return ligids_; }

	utility::vector1< core::Real > const &
	get_ligandchis() const { return ligandchis_; }

	core::Real
	get_ligandchi( core::Size ichi ) const { return ligandchis_[ichi]; }

	void
	set_ligandchi( core::Size ichi, core::Real value ) {
		ligandchis_[ichi] = value;
	}

	core::Size
	n_ligandchis() const { return ligandchis_.size(); }

	// get the moving sidechains
	utility::vector1< core::Size > const &
	moving_scs() const { return movingscs_; }

	void set_moving_scs( utility::vector1< core::Size > setting ){ movingscs_ = setting; }

	// get sidechain chi info; ires should be the movingsc index (not protein seqpos)
	void
	set_protein_restype( core::Size ires, core::chemical::ResidueTypeCOP restype ) {
		proteinrestypes_[ires] = restype;
		proteinchis_[ires].resize( restype->nchi(), 0.0 );
		ligandxyz_synced_ = false;
	}

	void
	set_protein_chis( core::Size ires, utility::vector1< core::Real > const & newchis ) {
		proteinchis_[ires] = newchis;
		ligandxyz_synced_ = false;
	}

	void
	set_sample_ring_conformers( bool setting ) { sample_ring_conformers_ = setting; }

	bool
	sample_ring_conformers() const { return sample_ring_conformers_; }

	core::chemical::ResidueTypeCOP
	get_protein_restype( core::Size ires ) const {
		return proteinrestypes_[ires];
	}

	utility::vector1< core::Real >
	get_protein_chis( core::Size ires ) const {
		return proteinchis_[ires];
	}

	// initial perturbation
	void
	randomize( core::Real transmax );

	void
	sample_conformation( core::Real transmax, TorsionSamplerCOP const & sampler);

	// Not defined yet
	// void
	// superimpose_to_alternative_frame( LigandConformer const &refconf );

	void set_rotwidth   ( core::Real setting ){ rotmutWidth_ = setting; }
	void set_transwidth ( core::Real setting ){ transmutWidth_ = setting; }
	void set_chiwidth   ( core::Real setting ){ ligchimutWidth_ = setting; }
	void set_torsmutrate( core::Real setting ){ torsmutationRate_ = setting; }
	void set_rtmutrate  ( core::Real setting ){ rtmutationRate_ = setting; }

	// assign_translation vector; used for building "apo" structure
	void
	assign_ligand_trans( core::Vector transv );

	//utility::vector1< numeric::xyzVector< core::Real > > points_to_search() const { return points_to_search_; }

	void score( core::Real scorein ) { score_ = scorein; }
	core::Real score() const { return score_; }

	void density_score( core::Real scorein ) { density_score_ = scorein; }
	core::Real density_score() const { return density_score_; }

	void rms( core::Real rmsin ) { rms_ = rmsin; }
	core::Real rms() const { return rms_; }

	void
	dump_pose( std::string pdbname ) const;

	void
	set_generation_tag( std::string tag ){ generation_tag_ = tag; }

	std::string
	generation_tag() const { return generation_tag_; }

	std::string
	to_string() const;

	core::kinematics::FoldTree const &
	get_reference_ft() const { return ref_pose_->fold_tree(); }

	core::Size
	get_jumpid() const { return jumpid_; }

	// fd returns true if the ligand is the last residue in the input pose
	bool
	is_ligand_terminal() const {
		return (
			std::find( ligids_.begin(), ligids_.end(), ref_pose_->total_residue() ) != ligids_.end()
		);
	}

	std::string
	ligand_typename() const {
		std::string retval = ligand_typenames_[1];
		for ( core::Size i = 2; i<=ligand_typenames_.size(); ++i ) {
			retval += "-"+ligand_typenames_[i];
		}
		return retval;
	}

	std::string
	ligand_typename(core::Size i) const {
		return ligand_typenames_[i];
	}

	void
	set_negTds( core::Real inval ) { negTdS_ = inval; }

	core::Real
	neg_Tds() const { return negTdS_; }

	core::pose::PoseCOP
	get_ref_pose() const { return ref_pose_; }

	// update the ligand chi torsion types
	void
	update_ligchi_types( core::conformation::Residue const& ligres);

	TorsionType const&
	get_ligchi_type( core::Size ndx) const {
		assert( ndx<=ligandchi_types_.size() );
		return ligandchi_types_[ndx];
	}
	utility::vector1< TorsionType > const&
	get_ligchi_types() const { return ligandchi_types_; }

	// would be nice to generalize this method
	void
	superimpose_to_ref_pose( utility::vector1< core::id::AtomID > const & ids  );

	void
	set_has_density_map( bool setting ) { has_density_map_ = setting; }

	bool
	has_density_map() const { return has_density_map_; }

	bool
	is_ligand_frozen() const { return freeze_ligand_; }

	bool
	is_ligand_bb_frozen() const { return ( freeze_ligand_ || freeze_ligand_backbone_ ); }


private:

	// reference pose
	core::pose::PoseCOP ref_pose_;

	// reference pose ligand
	core::pose::PoseCOP ref_pose_ligand_;

	// the resid of the ligand
	utility::vector1< core::Size > ligids_;

	// the sidechains that are allowed to move
	utility::vector1< core::Size > movingscs_;
	utility::vector1< core::chemical::ResidueTypeCOP > proteinrestypes_;
	utility::vector1< utility::vector1< core::Real > > proteinchis_;
	//initialized in init function
	bool sample_ring_conformers_;
	bool freeze_ligand_backbone_;
	bool freeze_ligand_;

	// the score
	core::Real score_;
	core::Real density_score_;

	// the rms to native (if given)
	core::Real rms_;

	// the internal representation of the pose
	utility::vector1< core::Real > rb_;
	utility::vector1< core::Real > ligandchis_;
	utility::vector1< TorsionType > ligandchi_types_;
	utility::vector1< core::id::TorsionID > ligandtorsionids_;
	std::map< core::Size, core::chemical::ResidueTypeCOP > ligid_restype_map_;
	utility::vector1< core::Real > ligandnus_;    // ring torsions
	utility::vector1< core::Real > ligandtaus_; // ring angles

	utility::vector1< utility::vector1< core::Size > > ligandchi_downstream_; // directional connections of torsions from root for alt crossover

	// radius of gyration of ligand
	core::Real rg_;

	// ligand jumpid
	core::Size jumpid_;

	//
	bool ligandxyz_synced_;
	utility::vector1< core::Vector > ligandxyz_;

	// mutation rate parameters
	core::Real torsmutationRate_, rtmutationRate_, transmutWidth_, rotmutWidth_, ligchimutWidth_;

	// history of conformation generation
	std::string generation_tag_;

	// ligand type name
	utility::vector1<std::string> ligand_typenames_;

	// -Tds
	core::Real negTdS_;

	bool has_density_map_ = false;
};

typedef utility::vector1< LigandConformer > LigandConformers;

// mutate a gene
LigandConformer
mutate(LigandConformer const &l );

// alternative mutation for many torsions
LigandConformer
mutate_ft( LigandConformer const &l, bool single_mutation );

// crossover
LigandConformer
crossover(LigandConformer const &l1, LigandConformer const &l2);

// alternative crossover for many torsions
LigandConformer
crossover_ft(LigandConformer const &l1, LigandConformer const &l2);

core::Real
distance_fast( LigandConformer &gene1, LigandConformer &gene2 ); // non-const because of sync check

std::pair< core::Real, core::Real >
distance_internal( LigandConformer const &gene1, LigandConformer const &gene2 );

core::Real
distance_slow( LigandConformer &gene1, LigandConformer &gene2 );

}
}
}

#endif
