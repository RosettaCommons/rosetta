// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/RotamerSet_.hh
/// @brief  rotamer set implementation class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_rotamer_set_rotamer_building_functions_hh
#define INCLUDED_core_pack_rotamer_set_rotamer_building_functions_hh

//Package headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>

// //Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/graph/Graph.fwd.hh>

// // Utility headers
#include <utility/io/izstream.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace rotamer_set {

typedef utility::vector1< Real > DihedralSet;

void
read_DNA_rotlib(
	utility::io::izstream & lib_stream,
	utility::vector1< DihedralSet* > & library
);


void
build_lib_dna_rotamers(
	utility::vector1< DihedralSet* > const & library,
	Size const resid,
	pose::Pose const & pose,
	chemical::ResidueTypeCOP concrete_residue,
	utility::vector1< conformation::ResidueOP > & rotamers
);


void
build_random_dna_rotamers(
	Size const resid,
	pose::Pose const & pose,
	chemical::ResidueTypeCOP concrete_residue,
	pack::task::ExtraRotSample const & level,
	utility::vector1< conformation::ResidueOP > & rotamers
);


void
build_dna_rotamers(
	Size const resid,
	pose::Pose const & pose,
	chemical::ResidueTypeCOP concrete_residue,
	pack::task::PackerTask const & task,
	utility::vector1< conformation::ResidueOP > & rotamers
);


/// @brief Make a rotamer (Residue) for every combination of torsion angle in the rotamer bins listed in the params
/// file for the given Residue.
void build_rotamers_from_rotamer_bins(conformation::Residue const & residue,
	utility::vector1<conformation::ResidueOP> & rotamers,
	uint current_chi_index=1,
	utility::vector1<uint> *current_bin_indices=NULL);


conformation::ResidueOP
create_oriented_water_rotamer(
	chemical::ResidueType const & h2o_type,
	Vector const & xyz_atom1,
	Vector const & xyz_O,
	Vector const & xyz_atom2,
	std::string const & name1,
	std::string const & name2,
	conformation::Residue const & tp5 // for the approx geometry of the lone pairs
);


void
build_independent_water_rotamers(
	Size const seqpos_water,
	chemical::ResidueType const & h2o_type,
	pack::task::PackerTask const & task,
	pose::Pose const & pose,
	graph::GraphCOP packer_neighbor_graph,
	utility::vector1< conformation::ResidueOP > & new_rotamers
);

void
debug_dump_rotamers(
	utility::vector1< conformation::ResidueOP > & rotamers
);

void
build_dependent_water_rotamers(
	RotamerSets const & rotsets,
	Size const seqpos_water,
	chemical::ResidueType const & h2o_type,
	pack::task::PackerTask const & task,
	pose::Pose const & pose,
	graph::GraphCOP packer_neighbor_graph,
	utility::vector1< conformation::ResidueOP > & new_rotamers
);


} // namespace rotamer_set
} // namespace pack
} // namespace core

#endif // INCLUDED_core_pack_RotamerSet_RotamerSet__HH

