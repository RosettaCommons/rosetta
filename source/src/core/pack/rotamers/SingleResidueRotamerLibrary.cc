// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamers/SingleResidueRotamerLibrary.cc
/// @brief  SingleResidueRotamerLibrary class
/// @author Andrew Leaver-Fay
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit Headers
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/BumpSelector.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/graph/Graph.hh>

#include <basic/Tracer.hh>

namespace core {
namespace pack {
namespace rotamers {

static THREAD_LOCAL basic::Tracer TR( "core.pack.rotamers.SingleResidueRotamerLibrary" );

SingleResidueRotamerLibrary::~SingleResidueRotamerLibrary()
{}

void
SingleResidueRotamerLibrary::bump_filter(
	RotamerVector & rotamers,
	core::Size resid,
	scoring::ScoreFunction const & scorefxn,
	pose::Pose const & pose,
	task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph
) const
{
	using namespace core::pack::rotamer_set; // Temporary for BumpSelectors

	if ( rotamers.size() == 0 || ! task.bump_check() ) { return; } //Nothing to do
	BumpSelector bump_selector( task.max_rotbump_energy() );

	RotamerVector passing;

	TR.Debug << "Filtering rotamers on bump energy for residue " << resid << " " << rotamers[1]->name() << std::endl;
	for ( Size ii = 1; ii <= rotamers.size(); ++ii ) {
		core::conformation::ResidueOP rot = rotamers[ ii ];
		TR.Debug << " Suggested rotamer " << ii << ": ";
		for ( Size jj = 1; jj <= rot->nchi(); ++jj ) {
			TR.Debug << rot->chi()[jj] << ' ';
		}

		core::PackerEnergy bumpenergy = bump_check( rot, resid, scorefxn, pose, task, packer_neighbor_graph );
		TR.Debug << " Bump energy: " << bumpenergy;
		BumpSelectorDecision decision =  bump_selector.iterate_bump_selector( bumpenergy );
		switch ( decision ) {
		case KEEP_ROTAMER :
			TR.Debug << " ... added" << std::endl;
			passing.push_back( rot );
			break;
		case DELETE_PREVIOUS_ROTAMER :
			TR.Debug << " ... replace previous" << std::endl;
			if ( passing.size() == 0 ) {
				utility_exit_with_message("Internal consistency error: cannot replace non-existant previous residue.");
			}
			passing[ passing.size() ] = rot;
			break;
		case DELETE_ROTAMER : // Do nothing.
			TR.Debug << " ... deleted" << std::endl;
			break;
		}
	}
	// Swap the in/out value with the contents of the temporary vector to pass back the value
	TR.Debug << " N rotamers before: " << rotamers.size() << " after: " << passing.size() << std::endl;
	rotamers.swap(passing);
}

/// @details Bump check does not include long range energies,
/// though, maybe this should change.
core::PackerEnergy
SingleResidueRotamerLibrary::bump_check(
	core::conformation::ResidueCOP rotamer,
	core::Size resid,
	scoring::ScoreFunction const & sf,
	pose::Pose const & pose,
	task::PackerTask const & task,
	graph::GraphCOP packer_neighbor_graph
) const
{
	using namespace scoring;
	using namespace conformation;

	EnergyMap emap;

	for ( graph::Graph::EdgeListConstIter
			ir  = packer_neighbor_graph->get_node( resid )->const_edge_list_begin(),
			ire = packer_neighbor_graph->get_node( resid )->const_edge_list_end();
			ir != ire; ++ir ) {
		int const neighbor_id( (*ir)->get_other_ind( resid ) );
		Residue const & neighbor( pose.residue( neighbor_id ) );

		if ( ! task.pack_residue( neighbor_id ) ) {
			sf.bump_check_full( *rotamer, neighbor, pose, emap);
		} else {
			sf.bump_check_backbone( *rotamer, neighbor, pose, emap);
		}
	}
	return static_cast< core::PackerEnergy > (sf.weights().dot( emap ));
}

core::Size
SingleResidueRotamerLibrary::current_rotamer(
	RotamerVector & rotamers,
	core::Size resid,
	task::PackerTask const & task,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const & existing_residue
) const
{
	if ( task.include_current( resid ) && existing_residue.name() == concrete_residue->name() ) {
		TR.Debug << " Including current rotamer: " <<
			existing_residue.name() << ' ' << existing_residue.seqpos() << std::endl;
		conformation::ResidueOP rot = existing_residue.create_rotamer();
		rotamers.push_back( rot );
		return rotamers.size();
	} else {
		return 0;
	}
}

void
SingleResidueRotamerLibrary::emergency_rotamer(
	RotamerVector & rotamers,
	core::Size /*resid*/,
	pose::Pose const & pose,
	task::PackerTask const & /*task*/,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const & existing_residue
) const
{
	if ( rotamers.size() == 0 && concrete_residue->nchi() == 0 ) {
		TR.Debug << " Using emergency rotamer" << std::endl;
		conformation::ResidueOP rot = conformation::ResidueFactory::create_residue( *concrete_residue, existing_residue, pose.conformation() );
		rotamers.push_back( rot );
	}
}

RotamerVector
SingleResidueRotamerLibrary::virtual_sidechain(
	RotamerVector const & rotamers,
	core::Size resid,
	pose::Pose const & pose,
	task::PackerTask const & task,
	chemical::ResidueTypeCOP /*concrete_residue*/,
	conformation::Residue const & existing_residue
) const {
	RotamerVector retval;
	if ( task.residue_task( resid ).include_virtual_side_chain() ) {
		if ( existing_residue.nchi() > 0 &&
				existing_residue.aa() != chemical::aa_pro &&
				existing_residue.n_non_polymeric_residue_connections() == 0 ) {
			dunbrack::RotamerLibraryScratchSpace scratch;
			Size n_min( 0 );
			Real fa_dun_min( 0.0 );
			for ( Size n = 1; n <= rotamers.size(); n++ ) {
				Real const fa_dun = this->rotamer_energy( *rotamers[ n ], scratch );
				if ( n_min == 0 || fa_dun < fa_dun_min ) {
					n_min = n; fa_dun_min = fa_dun;
				}
				// hack for now. PackerTask must be instantiated with a pose without virtual sidechains.
				runtime_assert( !rotamers[n]->has_variant_type( chemical::VIRTUAL_SIDE_CHAIN ) );
			}
			conformation::ResidueOP rot = core::pose::add_variant_type_to_residue( *rotamers[ n_min ], chemical::VIRTUAL_SIDE_CHAIN, pose);
			retval.push_back( rot );
		}
	}
	return retval;
}

/// @brief Equality test for equivalence.
/// Two SingleResidueRotamerLibraries test equal if and only if they represent the exact same behavior
bool
SingleResidueRotamerLibrary::operator ==( SingleResidueRotamerLibrary const & ) const {
	// If you're comparing arbitrary SingleResidueRotamerLibrary, chances are they aren't equal.
	// (Override your subclass if this doesn't work for you.)
	TR.Warning << "[ WARNING ] Program is trying to compare two arbitrary SingleResidueRotamerLibraries - this is probably a bug." << std::endl;
	return false;
}

} // namespace rotamers
} // namespace pack
} // namespace core
