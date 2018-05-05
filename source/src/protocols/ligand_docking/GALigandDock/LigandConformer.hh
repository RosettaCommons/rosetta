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

#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <numeric/Quaternion.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <map>


namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {


class LigandConformer : public utility::pointer::ReferenceCount {
public:
	friend LigandConformer
	mutate(LigandConformer const &l );

	// crossover
	friend LigandConformer
	crossover(LigandConformer const &l1, LigandConformer const &l2);

	// crossover2
	friend LigandConformer
	crossover_ft(LigandConformer const &l1, LigandConformer const &l2);

	friend core::Real
	distance_fast( LigandConformer &gene1, LigandConformer &gene2 );

	friend std::pair< core::Real, core::Real >
	distance_internal( LigandConformer const &gene1, LigandConformer const &gene2 );

	friend core::Real
	distance_slow( LigandConformer const &gene1, LigandConformer const &gene2 );

public:
	LigandConformer();
	~LigandConformer();

	LigandConformer(
		core::pose::PoseCOP pose,
		core::Size ligid,
		utility::vector1< core::Size > movingscs
	);

	// common pathway for parameter initialization
	void
	init_params();

	// set the reference pose and initial state
	void
	initialize(
		core::pose::PoseCOP pose,
		core::Size ligid,
		utility::vector1< core::Size > movingscs
	);

	bool
	defined() const { return (rb_.size() > 0); }

	// update internal representation based on this conformation
	void
	update_conf( core::pose::PoseCOP pose );

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
	ligand_residue() const;

	// generate only the protein residue based on this conformation
	core::conformation::Residue
	protein_residue( core::Size ires ) const;

	// get the id of the ligand
	core::Size
	ligand_id() const { return ligid_; }

	utility::vector1< core::Real > const &
	get_ligandchis() const { return ligandchis_; }

	core::Real
	get_ligandchi( core::Size ichi ) const { return ligandchis_[ichi]; }

	// get the moving sidechains
	utility::vector1< core::Size > const &
	moving_scs() const { return movingscs_; }

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

	void set_rotwidth   ( core::Real setting ){ rotmutWidth_ = setting; }
	void set_transwidth ( core::Real setting ){ transmutWidth_ = setting; }
	void set_chiwidth   ( core::Real setting ){ ligchimutWidth_ = setting; }
	void set_torsmutrate( core::Real setting ){ torsmutationRate_ = setting; }
	void set_rtmutrate  ( core::Real setting ){ rtmutationRate_ = setting; }

	// assign_translation vector; used for building "apo" structure
	void
	assign_ligand_trans( core::Vector transv );

	void score( core::Real scorein ) { score_ = scorein; }
	core::Real score() const { return score_; }

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

private:
	// reference pose
	core::pose::PoseCOP ref_pose_;

	// the resid of the ligand
	core::Size ligid_;

	// the sidechains that are allowed to move
	utility::vector1< core::Size > movingscs_;
	utility::vector1< core::chemical::ResidueTypeCOP > proteinrestypes_;
	utility::vector1< utility::vector1< core::Real > > proteinchis_;

	// the score
	core::Real score_;

	// the rms to native (if given)
	core::Real rms_;

	// the internal representation of the pose
	utility::vector1< core::Real > rb_;
	utility::vector1< core::Real > ligandchis_;
	utility::vector1< utility::vector1< core::Size > > ligandchi_downstream_;

	// radius of gyration of ligand
	core::Real rg_;

	//
	bool ligandxyz_synced_;
	utility::vector1< core::Vector > ligandxyz_;

	// mutation rate parameters
	core::Real torsmutationRate_, rtmutationRate_, transmutWidth_, rotmutWidth_, ligchimutWidth_, protchimutWidth_;

	// history of conformation generation
	std::string generation_tag_;
};

typedef utility::vector1< LigandConformer > LigandConformers;

// mutate a gene
LigandConformer
mutate(LigandConformer const &l );

// crossover
LigandConformer
crossover(LigandConformer const &l1, LigandConformer const &l2);

// crossover using some foldtree knowledge
LigandConformer
crossover_ft(LigandConformer const &l1, LigandConformer const &l2);

core::Real
distance_fast( LigandConformer &gene1, LigandConformer &gene2 ); // non-const because of sync check

std::pair< core::Real, core::Real >
distance_internal( LigandConformer const &gene1, LigandConformer const &gene2 );

core::Real
distance_slow( LigandConformer const &gene1, LigandConformer const &gene2 );

}
}
}

#endif
