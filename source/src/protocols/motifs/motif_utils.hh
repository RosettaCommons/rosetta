// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file motif_utils.hh
/// @brief Header for motif helper/conversion/io functions
/// @author havranek, sthyme (sthyme@gmail.com)

#ifndef INCLUDED_protocols_motifs_motif_utils_hh
#define INCLUDED_protocols_motifs_motif_utils_hh

// Unit Headers

// Package Headers
#include <protocols/motifs/BuildPosition.fwd.hh>
#include <protocols/motifs/Motif.fwd.hh>
#include <protocols/motifs/MotifLibrary.fwd.hh>
#include <protocols/motifs/SingleMotif.fwd.hh>

// Project Headers (protocols)
#include <protocols/dna/DnaDesignDef.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/file/FileName.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <map>
#include <set>

#include <utility/vector1.hh>
#include <utility/io/izstream.fwd.hh>


namespace protocols {
namespace motifs {

// If ever there are more than two states, use enum
core::Size const NO_LIGAND = 0;
core::Size const LIGAND = 1;


// This group of functions loads the motifs from external files
SingleMotifOP
single_motif_from_filename(
	utility::file::FileName const & motif_filename
);// uses a constructor that accepts two residues, not input atom names
// so motif atoms are set in the static map in Motif.cc

SingleMotifOP
single_motif_from_stream(
	utility::io::izstream & motif_info
);

SingleMotifOP
single_motif_from_stream(
	std::istream & motif_info
);

SingleMotifOP
single_ligand_motif_from_stream(
	std::istream & motif_info
);

// Dot product of the  z-axis of nucleotides (1.0 = parallel, 0.0 = not)
// This test is only applicable in the case of motifs that involve nucleic acids
core::Real
parallel_base_test(
	core::conformation::Residue const & pose_dna,
	core::conformation::Residue const & motif_dna
);

// The following group of functions are for havranek motif search
core::Real
backbone_stub_match(
	core::conformation::Residue const & r1,
	core::conformation::Residue const & r2
);

void
add_motif_bb_constraints(
	core::scoring::constraints::ConstraintSetOP cst_set,
	core::pose::Pose & pose,
	core::Size this_pos,
	core::conformation::Residue const & inv_rotamer
);

void add_motif_sc_constraints(
	core::scoring::constraints::ConstraintSetOP cst_set,
	core::pose::Pose & pose,
	core::Size this_pos,
	core::conformation::Residue const & inv_rotamer ,
	MotifCOP this_motif,
	bool const is_it_forward
);

void
mutate_loops_for_search(
	core::pose::Pose & pose,
	protocols::loops::Loops & flex_regions
);

void
mutate_position_vector_for_search(
	core::pose::Pose & pose,
	utility::vector1< core::Size > & trim_positions
);

// Functions for sthyme loading external data and setting up DNA mutations
MotifLibrary const
get_MotifLibrary_user();

MotifLibrary const
get_LigandMotifLibrary_user();

utility::vector1< core::conformation::ResidueOP >  const
get_targetconformers_user();

std::map< std::string, utility::vector1< core::conformation::ResidueOP > > const
setup_conformer_map(
	utility::vector1< core::conformation::ResidueOP > const & conformerOPs
);

utility::vector1< core::Size >
get_target_positions_make_dna_mutations(
	core::pose::Pose & pose
);

std::map< core::Size, std::set< std::string > >
get_target_position_map_make_dna_mutations(
	core::pose::Pose & pose
);

void
make_dna_mutations(
	core::pose::Pose & pose
);

void
make_dna_mutations(
	core::pose::Pose & pose,
	protocols::dna::DnaDesignDefOPs const & targeted_dna
);

utility::vector1< core::Size >
defs2vector(
	core::pose::Pose const & pose,
	protocols::dna::DnaDesignDefOPs const & targeted_dna
);

utility::vector1< std::pair< core::Size, utility::vector1< std::string > > >
defs2allowedtypes(
	core::pose::Pose const & pose,
	protocols::dna::DnaDesignDefOPs const & targeted_dna
);

std::map< core::Size, std::set< std::string > >
defs2map(
	core::pose::Pose const & pose,
	protocols::dna::DnaDesignDefOPs const & targets
);

std::map< core::Size, std::set< std::string > >
bpdefs2map(
	core::pose::Pose const & pose,
	protocols::dna::DnaDesignDefOPs const & targets
);

std::string
name3_from_oneletter(
	std::string const & oneletter
);

utility::vector1< core::Size >
get_motif_build_positions_user(
	core::pose::Pose const & pose
);

protocols::dna::DnaDesignDefOPs
get_motif_build_position_defs_user();

void
load_build_position_data(
	BuildPosition & bp,
	std::string const & filename,
	core::pose::Pose & pose,
	core::Size const ligand_marker = NO_LIGAND
);

utility::vector1< utility::file::FileName >
get_filenames(
	utility::vector1< utility::file::FileName > const & listnames
);

core::conformation::ResidueOP
single_residue_from_stream(
	utility::io::izstream & residueinfo
);

utility::vector1< bool >
bools_from_sizes(
	core::Size const nres,
	utility::vector1< core::Size > const & v
);

// Makes a base pair mutation, taken from devel for include reasons
void
make_base_pair_mutation(
	core::pose::Pose & pose,
	core::Size const seqpos,
	core::chemical::AA const & na
);

core::Real
atom_specific_rms(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	utility::vector1< std::string > const & atoms
);

core::Real
atom_specific_rms(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	utility::vector1< core::Size > const & atoms
);

core::pack::rotamer_set::RotamerSetOP
build_rotamers_lite(
	core::pose::Pose & pose,
	core::Size const rotamer_build_position,
	utility::vector1< bool > aa_info,
	core::Size const ex_,
	bool bump_check = true
);

} // namespace motifs
} // namespace protocols

#endif // INCLUDED_protocols_motifs_motif_utils
