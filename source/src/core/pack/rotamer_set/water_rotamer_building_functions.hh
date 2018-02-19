// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/water_rotamer_building_functions.hh
/// @brief  a few water rotamer building functions
/// @author Frank DiMaio & Ryan Pavlovicz


#ifndef INCLUDED_core_pack_rotamer_set_water_rotamer_building_functions_hh
#define INCLUDED_core_pack_rotamer_set_water_rotamer_building_functions_hh

//Package headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>

// //Project headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/graph/Graph.fwd.hh>

// // Utility headers
#include <utility/io/izstream.fwd.hh>
#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>


namespace core {
namespace pack {
namespace rotamer_set {

struct WaterRot {
	WaterRot( numeric::xyzVector< core::Real > const &coords, std::string aa, std::string aatm, std::string abase1, std::string abase2) :
		coords_(coords), aa_(aa), aatm_(aatm), abase1_(abase1), abase2_ (abase2)
	{ }

	numeric::xyzVector< core::Real > coords_;
	std::string aa_, aatm_, abase1_, abase2_;
};

class WaterRotsDB {
private:
	utility::vector1< WaterRot > protein_rots_,ligand_rots_;
	bool init_;

public:
	WaterRotsDB() : init_(false) {};

	bool
	is_initialized() {
		return init_;
	}

	core::Size
	n_protein_rots() {
		return protein_rots_.size();
	}

	core::Size
	n_ligand_rots() {
		return ligand_rots_.size();
	}

	void
	initialize();

	WaterRot const &
	protein_rot(core::Size i) {
		return protein_rots_[i];
	}

	WaterRot const &
	ligand_rot(core::Size i) {
		return ligand_rots_[i];
	}
};

core::PackerEnergy
bump_check(
	core::conformation::ResidueCOP rotamer,
	core::Size resid,
	scoring::ScoreFunction const & sf,
	pose::Pose const & pose,
	task::PackerTask const & task,
	utility::graph::GraphCOP packer_neighbor_graph
);

void
build_rotated_water_rotamers(
	Size const seqpos_water,
	pack::task::PackerTask const & task,
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	utility::graph::GraphCOP packer_neighbor_graph,
	utility::vector1< conformation::ResidueOP > & new_rotamers,
	bool incl_vrt
);

void
build_backbone_point_water_rotamers(
	Size const seqpos_water,
	pack::task::PackerTask const & task,
	pose::Pose const & pose,
	utility::graph::GraphCOP , //packer_neighbor_graph,
	utility::vector1< conformation::ResidueOP > & new_rotamers,
	bool incl_vrt
);


} // namespace rotamer_set
} // namespace pack
} // namespace core

#endif // INCLUDED_core_pack_RotamerSet_RotamerSet__HH

