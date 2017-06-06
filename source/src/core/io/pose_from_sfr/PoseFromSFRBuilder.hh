// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/pose_builder/PoseFromSFRBuilder.hh
/// @brief  Declarations for PoseFromSFRBuilder and related classes.
/// @author Sergey Lyskov
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_io_pose_builder_PoseFromSFRBuilder_HH
#define INCLUDED_core_io_pose_builder_PoseFromSFRBuilder_HH

// Package headers
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/ResidueInformation.hh>
#include <core/io/pose_from_sfr/chirality_resolution.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/NamedAtomID_Map.hh>
#include <core/chemical/MergeBehaviorManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/PDBInfo.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>

// C++ headers
#include <map>
#include <string>


namespace core {
namespace io {
namespace pose_from_sfr {

/// @brief The %PoseFromSFRBuilder is responsible for taking a partially-processed representation of a
/// structure file, a structure-file repersentation or StructFileRep, and it constructs a Pose from it.
/// Its primary output is a Pose, but it also keeps track of the "missing atoms:" the set of atoms in
/// the Pose whose coordinates were not given in the original file.
///
/// @details The construction involves renaming some atoms, for example, the O2* atom in RNA would be
/// renamed O2', and resolving the chain- and branch-connectivity.  It uses the ResidueTypeFinder to
/// select the best-fit ResidueType for each residue.
///
/// This code handles n-terminal acetylation, and some other similar chemical modifications, by merging
/// the atoms in the "ACE" residue with the next residue.  Instructions of what residue types to merge
/// and how to rename the atoms to avoid naming collisions are read in from the database by the
/// MergeBehaviorManager.  After the atoms are merged and renamed, the appropriate AcetylatedNtermProteinFull
/// is applied to that next residue.
///
/// The process of building a pose is handled in four passes.
class PoseFromSFRBuilder
{
public:

	PoseFromSFRBuilder( chemical::ResidueTypeSetCOP rts, StructFileRepOptions const & options );
	~PoseFromSFRBuilder();

	/// @brief Build a Pose from an SFR using the ResidueTypeSet and the StructFileReaderOptions
	/// object that were passed in to the %Builder in it's constructor.
	void build_pose( StructFileRep const & sfr, pose::Pose & pose );

	// @brief Returns a record of which atoms in the most-recently-built Pose were
	/// not present in the SFR
	id::AtomID_Mask const &
	missing_atoms() const;



private:
	void setup( StructFileRep const & sfr );
	void pass_1_merge_residues_as_necessary();
	void pass_2_resolve_residue_types();
	void pass_3_verify_sufficient_backbone_atoms();
	void pass_4_redo_termini();
	void pass_5_note_discarded_atoms();
	void build_initial_pose( pose::Pose & pose );
	void refine_pose( pose::Pose & pose );
	void build_pdb_info_1_everything_but_temps( pose::Pose & pose );
	void build_pdb_info_2_temps( pose::Pose & pose );

	bool update_atom_information_based_on_occupancy( AtomInformation & ai ) const;

	/// @brief randomize missing density
	void randomize_missing_coords( AtomInformation & ai ) const;

	/// @brief  This function uses linkage information to determine main-chain
	/// and branch polymer connectivity.
	void
	determine_residue_branching_info(
		Size seqpos,
		std::string const & name3,
		bool & same_chain_prev,
		bool & same_chain_next,
		utility::vector1< std::string > & branch_points_on_this_residue,
		bool & is_branch_point,
		bool & is_branch_lower_terminus );

	bool
	is_residue_type_recognized(
		Size const pdb_residue_index,
		std::string const & rosetta_residue_name3,
		core::chemical::ResidueTypeCOPs const & rsd_type_list
	);

	bool
	is_residue_type_recognized(
		Size const pdb_residue_index,
		std::string const & rosetta_residue_name3
	);

	chemical::ResidueTypeCOP
	get_rsd_type(
		std::string const & name3,
		Size seqpos,
		utility::vector1< std::string > const &  branch_points_on_this_residue,
		std::string const & resid,
		bool const is_lower_terminus,
		bool const is_upper_terminus,
		bool const is_branch_point,
		bool const is_branch_lower_terminus,
		bool const is_d_aa,
		bool const is_l_aa
	);

	void convert_nucleic_acid_residue_info_to_standard();

	Size prev_residue_skipping_merges( Size resid ) const;
	Size prev_residue_skipping_null_residue_types( Size resid ) const;
	Size next_residue_skipping_merges( Size resid ) const;
	Size next_residue_skipping_null_residue_types( Size resid ) const;

	bool determine_separate_chemical_entity( char chainID ) const;
	bool determine_same_chain_prev( Size resid, bool separate_chemical_entity ) const;
	bool determine_same_chain_next( Size resid, bool separate_chemical_entity ) const;
	bool determine_check_Ntermini_for_this_chain( char chainID ) const;
	bool determine_check_Ctermini_for_this_chain( char chainID ) const;

	void fill_name_map( Size seqpos );

	Size find_atom_tree_root_for_metal_ion( core::pose::Pose const & pose, conformation::ResidueCOP metal_rsd );
	bool last_residue_was_recognized( Size seqpos ) const;
	void output_ignore_water_warning_once();

private:
	chemical::ResidueTypeSetCOP residue_type_set_;
	StructFileRep sfr_;
	StructFileRepOptions options_;

	// Records which atoms in the output Pose were not present in the SFR
	id::AtomID_Mask missing_;

	utility::vector1< core::io::ResidueInformation > rinfos_;
	id::NamedAtomID_Mask coordinates_assigned_;
	utility::vector1< Size > pose_to_rinfo_;
	utility::vector1< StructFileRep::ResidueTemps> pose_temps_;
	utility::vector1< core::pose::UnrecognizedAtomRecord > unrecognized_atoms_;
	StructFileRep::Strings branch_lower_termini_;
	StructFileRep::Strings branch_lower_termini_extra_;

	utility::vector1< chemical::merge_residue_behavior > merge_behaviors_;
	utility::vector1< std::map< std::string, std::string > > merge_atom_maps_;

	utility::vector1< std::string > rosetta_residue_name3s_;
	utility::vector1< core::Size > glycan_positions_;
	utility::vector1< std::string > glycan_tree_roots_;
	utility::vector1< chemical::ResidueTypeCOP > residue_types_;
	utility::vector1< NameBimap > remapped_atom_names_;
	utility::vector1< bool > is_lower_terminus_;
	utility::vector1< bool > is_branch_lower_terminus_;
	utility::vector1< bool > same_chain_prev_;
	utility::vector1< bool > residue_was_recognized_;

	bool outputted_ignored_water_warning_;
};

std::string convert_atom_name( std::string const & res_name, std::string const & atom_name );

std::string convert_res_name( std::string const & name );

void
create_working_data(
	StructFileRepOptions const & options,
	StructFileRep const & sfr,
	utility::vector1< core::io::ResidueInformation > & rinfos
);

bool
update_atom_information_based_on_occupancy(
	StructFileRepOptions const & options,
	AtomInformation & ai
);

/// @brief randomize missing density
void randomize_missing_coords( AtomInformation & ai );

} // namespace pose_from_sfr
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_pdb_file_data_HH
