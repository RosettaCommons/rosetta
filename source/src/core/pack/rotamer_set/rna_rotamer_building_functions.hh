// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/RotamerSet_.hh
/// @brief  rotamer set implementation class
/// @author Rhiju Das (rhiju@stanford.edu)

#ifndef INCLUDED_core_pack_rotamer_set_rna_rotamer_building_functions_hh
#define INCLUDED_core_pack_rotamer_set_rna_rotamer_building_functions_hh

//Package headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>

// //Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/graph/Graph.fwd.hh>

// // Utility headers
#include <utility/io/izstream.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace rotamer_set {

void
build_rna_rotamers(
	Size const resid,
	pose::Pose const & pose,
	chemical::ResidueTypeCOP concrete_residue,
	pack::task::PackerTask const & task,
	utility::vector1< conformation::ResidueOP > & rotamers,
	Size & id_for_current_rotamer
);


} // namespace rotamer_set
} // namespace pack
} // namespace core

#endif // INCLUDED_core_pack_RotamerSet_RotamerSet__HH

