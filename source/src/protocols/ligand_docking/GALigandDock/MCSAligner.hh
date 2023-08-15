// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/MCSAligner.hh
///
/// @brief  Use RDKit maximum common substructure to align a ligand to a reference ligand
/// @author Guangfeng Zhou and Frank DiMaio

#ifndef INCLUDED_protocols_ligand_docking_GALigandDock_MCSAligner_hh
#define INCLUDED_protocols_ligand_docking_GALigandDock_MCSAligner_hh

#include <protocols/ligand_docking/GALigandDock/LigandConformer.fwd.hh>
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <core/id/AtomID.hh> // AUTO IWYU For AtomID


namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

/// @brief
/// Aligns ligand using maximum commom substructure (MCS)
/// @details
/// Performs ligand alignments using MCS bewteen
/// - reference ligand pose
/// - ligand to dock

struct MCSAlignerOptions{
	core::chemical::rdkit::RestypeToRDMolOptions restype_to_rdmol_options;
	bool perturb_rb = false;
	bool perturb_torsion = false;
	core::Real perturb_rb_translation = 2.0; //Angstrom
	core::Real perturb_rb_rotation = 15.0; //Degree
	core::Real perturb_torsion_rotation = 15.0; //Degree
	MCSAlignerOptions(){
		restype_to_rdmol_options.neutralize = false;
		restype_to_rdmol_options.keep_hydro = true;
		restype_to_rdmol_options.sanitize = false;
		restype_to_rdmol_options.noImplicitHs = true;
		restype_to_rdmol_options.skipHs = true;
		restype_to_rdmol_options.aro2double = true;
	}
};

class MCSAligner {
public:

	MCSAligner( core::pose::Pose const& reference_pose, core::Size reference_ligres_idx, MCSAlignerOptions & options );

	/// @brief set MCSAligner options
	void set_options(MCSAlignerOptions const& options) { options_ = options; }

	/// @brief main apply function
	void
	apply( LigandConformer & lig );

	void
	align_pose(core::pose::Pose const& template_pose, core::pose::Pose &ligand_pose,
		std::map<core::Size, core::Size> const& pair_indices_map, core::Size const& template_idx, core::Size const& ligand_idx);

	void
	set_torsion_and_align(core::pose::Pose const& template_pose, core::pose::Pose &ligand_pose,
		std::map<core::Size, core::Size> const& pair_indices_map, core::Size const& template_idx, core::Size const& ligand_idx);

	utility::vector1<bool> const& torsion_in_align() const { return torsion_in_align_; }

private:
	/// @brief set constraints to target
	void set_constraints(
		core::pose::Pose & pose,
		utility::vector1<core::Size> ligids,
		utility::vector1< std::pair< core::Size, core::Size > > &marked_pairs,
		core::Real const w_prior = 1.0, // default no upper limit
		utility::vector1< core::Size > const &SrcPriorIDs = utility::vector1< core::Size >(),
		utility::vector1< core::Size > const &TgtPriorIDs = utility::vector1< core::Size >()
	);


private:
	core::pose::Pose const& reference_pose_;
	core::Size reference_ligres_idx_;
	MCSAlignerOptions options_;
	bool perturb_rb_, perturb_torsion_;
	utility::vector1<bool> torsion_in_align_;

};

typedef utility::vector1<std::pair<core::Size, core::Size>> SizePairVec;
typedef utility::keys::Key4Tuple< core::Size, core::Size, core::Size, core::Size > DihedralAtomTuple;
typedef std::map< utility::keys::Key4Tuple< core::Size, core::Size, core::Size, core::Size >, core::Size> DihedralAtomTuple2ChiIdxMap;

typedef utility::pointer::shared_ptr< MCSAligner > MCSAlignerOP;
typedef utility::pointer::shared_ptr< MCSAligner const > MCSAlignerCOP;

}
}
}

#endif

