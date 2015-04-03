// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_dna_util_hh
#define INCLUDED_protocols_dna_util_hh

#include <protocols/dna/DnaChains.fwd.hh>
#include <protocols/dna/typedefs.hh>
#include <protocols/dna/DnaDesignDef.fwd.hh>

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector0.fwd.hh>
#include <utility/vector1.hh> // there is no forward declaration possible for const_iterator(?)

// C++ headers
#include <list>
#include <iosfwd>


namespace protocols {
namespace dna {

/// @brief basic struct for remembering position/type information before/during/after design
class PositionType {
public:
	PositionType(
		core::Size pos = 0,
		core::chemical::ResidueTypeCOP rt = 0,
		bool design = false
	) : position( pos ), type( rt ), designable( design ) {}

	core::Size position;
	core::chemical::ResidueTypeCOP type;
	bool designable;
};

typedef utility::vector1< PositionType > PositionTypes;


/// @details checks c-beta (except glycine) to base atom-atom distances, not including ribose or phosphate backbone.
/// @author ashworth
bool
close_to_dna(
	core::conformation::Residue const & pres,
	core::conformation::Residue const & dres,
	core::Real cut2,
	bool base_only = false
);

/// @details arginine rotamer sweep at a protein residue to see if it should be considered a (potentially) 'dna-contacting' residue
/// @author ashworth
core::Real
argrot_dna_dis2(
	core::pose::Pose const & pose,
	core::Size presid,
	core::conformation::Residue const & pres,
	core::conformation::Residue const & dres,
	core::Real threshold,
	bool base_only = false
);

/// @details distance check for contact between two sets of atoms
/// @author ashworth
core::Real
contact_distance2(
	core::conformation::Atoms::const_iterator a_begin,
	core::conformation::Atoms::const_iterator a_end,
	core::conformation::Atoms::const_iterator b_begin,
	core::conformation::Atoms::const_iterator b_end,
	core::Real threshold = 0.0
);

/// @details A sanity check for the arginine rotamer screen. Can prevent the design of positions that are best left alone because they are too far away along the helical axis ('laterally').
/// @author ashworth
core::Real
z_axis_dist(
	core::conformation::Residue const & pres,
	core::conformation::Residue const & dres
);

/// @brief also consider using the dna_base_partner function below
/// @author ashworth
std::string dna_comp_name_str( std::string const & dna );

/// @brief intended to convert any DNA "threeletter code" into the full three-letter code. Note that this does not (necessarily) return the same thing as residue_type::name3 (which returns "  N" format as of Dec 2008)
std::string dna_full_name3( std::string const & name3 );

core::chemical::AA
dna_base_partner( core::chemical::AA const & na );

/// @details DnaChains version, adapted from pbradley's code.  More paranoid geometry checks, in order to allow highly distorted basepairs without making mistakes
/// @author ashworth
void
find_basepairs(
	core::pose::Pose const & pose,
	DnaChains & dna_chains,
	bool include_unpaired = true
);

/* void
make_sequence_combinations(
	ResTypeSequences & sequences,
	utility::vector1< core::Size > const & seq_indices,
	core::pack::task::PackerTaskCOP ptask
); */

/// @brief make a list of all single mutants from a base sequence
/// @author ashworth
void
make_single_mutants(
	ResTypeSequence const & sequence,
	core::pack::task::PackerTaskCOP ptask,
	ResTypeSequences & sequences
);

void
design_residues_list(
	std::list< PositionType > & design_residues,
	core::pose::Pose const & pose,
	core::pack::task::PackerTask const & ptask
);

void
make_sequence_combinations(
	utility::vector1< core::Size >::const_iterator seqset_iter,
	utility::vector1< core::Size > const & seq_indices,
	core::pack::task::PackerTaskCOP ptask,
	ResTypeSequence & sequence,
	ResTypeSequences & sequences
);

std::ostream & operator << ( std::ostream & os, ResTypeSequence const & seq );

std::string seq_to_str( ResTypeSequence const & seq );

std::ostream & operator << ( std::ostream & os, ResTypeSequences const & seqs );

std::string seq_pdb_str(
	ResTypeSequence const &,
	core::pose::Pose const &
);

void print_sequence_pdb_nums(
	ResTypeSequence const &,
	core::pose::Pose const &,
	std::ostream &
);

void print_sequences_pdb_nums(
	ResTypeSequences const &,
	core::pose::Pose const &,
	std::ostream &
);

/// @details for packing a single DNA sequence out of a multi-DNA-sequence RotamerSet
/// @author ashworth
void
restrict_dna_rotamers(
	core::pack::rotamer_set::RotamerSetsCOP rotsets,
	ResTypeSequence const & seq,
	utility::vector0<int> & rot_to_pack
);

/// @details for packing a single sequence out of a RotamerSets that (potentially) represents sequence variability
/// @author ashworth
void
restrict_to_single_sequence(
	core::pack::rotamer_set::RotamerSetsCOP rotamer_sets,
	utility::vector1< core::chemical::ResidueTypeCOP > const & single_sequence,
	utility::vector0< int > & rot_to_pack
);

void
substitute_residue(
	core::pose::Pose & pose,
	core::Size index,
	core::chemical::ResidueType const & new_type
);

void write_checkpoint( core::pose::Pose & pose, core::Size iter );
void load_checkpoint( core::pose::Pose & pose, core::Size & iter );
void checkpoint_cleanup();

/// @brief loads command-line dna design definitions (shorthand alternative to using resfile)
/// option value is string vector
///   i.e. -dna_defs C.-6 C.-5
///   or   -dna_defs C.-6.GUA C.-5.CYT
/// @author ashworth
void
load_dna_design_defs_from_strings(
	DnaDesignDefOPs & defs,
	utility::vector1< std::string > const & str_defs
);

void
load_dna_design_defs_from_file(
	DnaDesignDefOPs & defs,
	std::string const & filename,
	std::string const & pdb_prefix
);

void
load_dna_design_defs_from_options(
	DnaDesignDefOPs & defs,
	std::string pdb_prefix = std::string()
);

void
add_constraints_from_file(
	core::pose::Pose & pose
);

core::kinematics::FoldTree
make_base_pair_aware_fold_tree(
	core::pose::Pose const & pose
);

bool
not_already_connected(
	core::pose::Pose const & pose,
	core::Size const num_jumps,
	char const this_chain,
	char const other_chain,
	ObjexxFCL::FArray2D_int & jump_pairs
);


void
set_base_segment_chainbreak_constraints(
	core::pose::Pose & pose,
	core::Size const this_base,
	core::Size const end_base
);


} // namespace dna
} // namespace protocols

#endif
