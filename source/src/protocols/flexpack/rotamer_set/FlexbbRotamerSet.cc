// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/flexpack/rotamer_set/FlexbbRotamerSet.hh
/// @brief  Declaration for a class to hold rotamers for a single backbone conformation in
///  a flexible packing run
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), Florian Richter (floric@u.washington.edu), sep 08

#include <protocols/flexpack/rotamer_set/FlexbbRotamerSet.hh>
#include <protocols/flexpack/rotamer_set/FlexbbRotamerSets.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace flexpack {
namespace rotamer_set {

FlexbbRotamerSet::FlexbbRotamerSet() :
	parent(),
	existing_residue_( /* 0 */ ),
	owner_( /* 0 */ )
{}

FlexbbRotamerSet::~FlexbbRotamerSet()
{}

void
FlexbbRotamerSet::set_owner( FlexbbRotamerSetsCAP owner )
{
	owner_ = owner;
}


void
FlexbbRotamerSet::build_dependent_rotamers(
	core::pack::rotamer_set::RotamerSets const & /*rotamer_sets*/,
	core::pose::Pose const & /*pose*/,
	core::scoring::ScoreFunction const & /*scorefxn*/,
	core::pack::task::PackerTask const & /*task*/,
	core::graph::GraphCOP /*packer_neighbor_graph*/
)
{
	utility_exit_with_message("ERROR: FlexbbRotamerSet does not support build_dependent_rotamers");
}

void
FlexbbRotamerSet::set_existing_residue( core::conformation::ResidueCOP residue )
{
	existing_residue_ = residue;
}


void
FlexbbRotamerSet::build_rotamers_for_concrete_virt(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & scorefxn,
	core::pack::task::PackerTask const & task,
	core::chemical::ResidueTypeCOP concrete_residue,
	core::graph::GraphCOP packer_neighbor_graph,
	bool use_neighbor_context
)
{
	build_rotamers_for_concrete( pose, scorefxn, task,
		concrete_residue, *existing_residue_,
		packer_neighbor_graph, use_neighbor_context );
}

/// @brief Computes the "bump energy" of a rotamer: the bump energy is the
/// sum of rotamer's interactions with 1) the backbone-and-side chains of
/// neighboring residues that are held fixed during this repacking optimization
/// and 2) the backbones of neighboring residues that are changable during this
/// repacking optimization.
core::PackerEnergy
FlexbbRotamerSet::bump_check(
	core::conformation::ResidueCOP rotamer,
	core::scoring::ScoreFunction const & sf,
	core::pose::Pose const & pose,
	core::pack::task::PackerTask const & task,
	core::graph::GraphCOP packer_neighbor_graph
) const
{
	/// iterate across neighbors and do shit... that's right

	using namespace core::scoring;
	using namespace core::conformation;

	core::Real bumpE( 0.0 );
	protocols::flexpack::rotamer_set::FlexbbRotamerSetsCOP owner( owner_ );

	for ( core::graph::Graph::EdgeListConstIter
			ir  = packer_neighbor_graph->get_node( resid() )->const_edge_list_begin(),
			ire = packer_neighbor_graph->get_node( resid() )->const_edge_list_end();
			ir != ire; ++ir ) {

		core::Real smallest_neighbor_bumpE( 0.0 );

		int const neighbor_id( (*ir)->get_other_ind( resid() ) );

		core::Size bbconfs_this_neighbor( owner->nbbconfs_for_res( neighbor_id ) );

		utility::vector1< Residue const * > check_residues;

		if ( bbconfs_this_neighbor == 1 ) check_residues.push_back( & pose.residue( neighbor_id ) );

		else {
			for ( core::Size bbconf = 1; bbconf <= bbconfs_this_neighbor ; ++bbconf ) {
				check_residues.push_back( & owner->backbone_for_resid_bbconf( neighbor_id, bbconf ) );
			}
		}

		for ( core::Size check_res = 1; check_res <= check_residues.size();  ++check_res ) {

			EnergyMap emap;

			if ( ! task.pack_residue( neighbor_id ) ) {
				sf.bump_check_full( *rotamer, *check_residues[ check_res ], pose, emap);
			} else {
				sf.bump_check_backbone( *rotamer, *check_residues[ check_res ], pose, emap);
			}
			core::Real cur_bumpE = sf.weights().dot( emap );

			if ( check_res == 1 || cur_bumpE < smallest_neighbor_bumpE ) smallest_neighbor_bumpE = cur_bumpE;

		} //iterator over flexible states of this neighbor

		bumpE += smallest_neighbor_bumpE;

	} // iterator over neighbors

	return static_cast< core::PackerEnergy >( bumpE );
} //bump_check


}
}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::flexpack::rotamer_set::FlexbbRotamerSet::save( Archive & arc ) const {
	arc( cereal::base_class< core::pack::rotamer_set::RotamerSet_ >( this ) );
	arc( CEREAL_NVP( existing_residue_ ) ); // ResidueCOP
	arc( CEREAL_NVP( owner_ ) ); // FlexbbRotamerSetsCAP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::flexpack::rotamer_set::FlexbbRotamerSet::load( Archive & arc ) {
	arc( cereal::base_class< core::pack::rotamer_set::RotamerSet_ >( this ) );

	std::shared_ptr< core::conformation::Residue > local_existing_residue;
	arc( local_existing_residue ); // ResidueCOP
	existing_residue_ = local_existing_residue; // copy the non-const pointer(s) into the const pointer(s)

	std::weak_ptr< FlexbbRotamerSets > local_owner;
	arc( local_owner ); // FlexbbRotamerSetsCAP
	owner_ = local_owner;
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::flexpack::rotamer_set::FlexbbRotamerSet );
CEREAL_REGISTER_TYPE( protocols::flexpack::rotamer_set::FlexbbRotamerSet )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_flexpack_rotamer_set_FlexbbRotamerSet )
#endif // SERIALIZATION
