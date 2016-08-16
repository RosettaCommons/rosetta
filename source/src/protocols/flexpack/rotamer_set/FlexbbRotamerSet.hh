// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/rotamer_set/FlexbbRotamerSet.hh
/// @brief  Declaration for a class to hold rotamers for a single backbone conformation in
///  a flexible packing run
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), Florian Richter (floric@u.washington.edu), sep 08

#ifndef INCLUDED_protocols_flexpack_rotamer_set_FlexbbRotamerSet_hh
#define INCLUDED_protocols_flexpack_rotamer_set_FlexbbRotamerSet_hh

//Unit headers
#include <protocols/flexpack/rotamer_set/FlexbbRotamerSet.fwd.hh>

#include <core/pack/rotamer_set/RotamerSet_.hh>

#include <protocols/flexpack/rotamer_set/FlexbbRotamerSets.fwd.hh>

#include <utility/vector1.hh>


// Package Headers

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace flexpack {
namespace rotamer_set {

class FlexbbRotamerSet : public core::pack::rotamer_set::RotamerSet_
{

public:
	typedef core::pack::rotamer_set::RotamerSet_ parent;

	FlexbbRotamerSet();
	virtual ~FlexbbRotamerSet();

	void
	set_owner( FlexbbRotamerSetsCAP owner );

	/// @brief  Build rotamers that depend on positions of rotamers built in a previous pass
	/// This function won't work...
	virtual
	void build_dependent_rotamers(
		core::pack::rotamer_set::RotamerSets const & rotamer_sets,
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & scorefxn,
		core::pack::task::PackerTask const & task,
		core::graph::GraphCOP packer_neighbor_graph
	);

	void
	set_existing_residue( core::conformation::ResidueCOP residue );

	//core::conformation::Residue const &
	//existing_residue() const;


protected:

	/// @brief Creates a set of rotamers for a particular residue type
	/// (the concrete residue type) while relying on the rotamer-
	/// building instructions within the PackerTask.
	/// Use the residue in the input pose at position resid_ as the existing residue.
	virtual
	void build_rotamers_for_concrete_virt(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & scorefxn,
		core::pack::task::PackerTask const & task,
		core::chemical::ResidueTypeCOP concrete_residue,
		core::graph::GraphCOP packer_neighbor_graph,
		bool use_neighbor_context = true
	);

	/// @brief Computes the "bump energy" of a rotamer: the bump energy is the
	/// sum of rotamer's interactions with 1) the backbone-and-side chains of
	/// neighboring residues that are held fixed during this repacking optimization
	/// and 2) the backbones of neighboring residues that are changable during this
	/// repacking optimization.
	virtual
	core::PackerEnergy
	bump_check(
		core::conformation::ResidueCOP rotamer,
		core::scoring::ScoreFunction const & sf,
		core::pose::Pose const & pose,
		core::pack::task::PackerTask const & task,
		core::graph::GraphCOP packer_neighbor_graph
	) const;

	virtual
	void
	initialize_pose_for_rotsets_creation(
		core::pose::Pose & /*pose*/
	) const {}

private:
	ResidueCOP existing_residue_;
	FlexbbRotamerSetsCAP owner_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_flexpack_rotamer_set_FlexbbRotamerSet )
#endif // SERIALIZATION


#endif
