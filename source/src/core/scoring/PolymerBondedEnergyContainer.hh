// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/PolymerBondedEnergyContainer.hh
/// @brief  A container interface long range energies for polymer-bonded residue interactions only.
/// @author Frank DiMaio
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Modified 21 February 2016 so that this no longer just scores n->n+1 interactions, but includes
/// anything that's polymer-bonded (whether or not it's adjacent in linear sequence).  This includes N-to-C cyclic peptide bonds.

#ifndef INCLUDED_core_scoring_PolymerBondedEnergyContainer_hh
#define INCLUDED_core_scoring_PolymerBondedEnergyContainer_hh

// Unit headers
#include <core/scoring/PolymerBondedEnergyContainer.fwd.hh>

// Package headers
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

// STL headers:
#include <map>
#include <set>


namespace core {
namespace scoring {

///////////////////////////////////////////////////////

class PolymerBondedNeighborIterator : public ResidueNeighborIterator
{
	PolymerBondedNeighborIterator & operator = (PolymerBondedNeighborIterator const & src );

public:
	~PolymerBondedNeighborIterator() override;

	// Moves pos_in, so by-value
	PolymerBondedNeighborIterator(
		Size const base_in,
		utility::vector1< Size > const & positions_in,
		PolymerBondedEnergyContainer & parent
	);

	ResidueNeighborIterator & operator = ( ResidueNeighborIterator const & src ) override;

	ResidueNeighborIterator const & operator ++ () override;

	bool operator == ( ResidueNeighborIterator const & other ) const override;

	bool operator != ( ResidueNeighborIterator const & other ) const override;

	Size upper_neighbor_id() const override;

	Size lower_neighbor_id() const override;

	Size residue_iterated_on() const override;

	Size neighbor_id() const override;

	void save_energy( EnergyMap const & emap ) override;

	void retrieve_energy( EnergyMap & emap ) const override;

	void accumulate_energy( EnergyMap & emap ) const override;

	void mark_energy_computed() override;

	void mark_energy_uncomputed() override;

	bool energy_computed() const override;

private:
	Size base_, curr_idx_;
	utility::vector1< Size > pos_; // positions to loop over (bonded neighbors of base)
	PolymerBondedEnergyContainer *parent_;
};


///////////////////////////////////////////////////////

class PolymerBondedNeighborConstIterator : public ResidueNeighborConstIterator
{
	PolymerBondedNeighborConstIterator & operator = (PolymerBondedNeighborConstIterator const & src );
public:
	~PolymerBondedNeighborConstIterator() override;

	// Moves pos_in, so by value
	PolymerBondedNeighborConstIterator(
		Size const base_in,
		utility::vector1< Size > const & positions_in,
		PolymerBondedEnergyContainer const & parent
	);

	ResidueNeighborConstIterator & operator = ( ResidueNeighborConstIterator const & src ) override;

	ResidueNeighborConstIterator const & operator ++ () override;

	bool operator == ( ResidueNeighborConstIterator const & other ) const override;

	bool operator != ( ResidueNeighborConstIterator const & other ) const override;

	Size upper_neighbor_id() const override;

	Size lower_neighbor_id() const override;

	Size residue_iterated_on() const override;

	Size neighbor_id() const override;

	void retrieve_energy( EnergyMap & emap ) const override;

	void accumulate_energy( EnergyMap & emap ) const override;

	bool energy_computed() const override;

private:
	Size base_, curr_idx_;
	utility::vector1< Size > pos_; // positions to loop over (bonded neighbors of base)
	PolymerBondedEnergyContainer const *parent_;
};

///////////////////////////////////////////////////////////////////////////

class PolymerBondedEnergyContainer : public LREnergyContainer {
public:

	~PolymerBondedEnergyContainer() override;


	LREnergyContainerOP clone() const override;

	/// @brief Pose constructor.
	/// @details Initializes PolymerBondedEnergyContainer from a pose, facilitating calculations involving non-canonical connections
	/// (e.g. terminal peptide bonds).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	PolymerBondedEnergyContainer( core::pose::Pose const & pose, utility::vector1< ScoreType > const & score_type_in );


	bool empty() const override;


	bool
	any_neighbors_for_residue( int /*resid*/ ) const override;


	bool
	any_upper_neighbors_for_residue( int /*resid*/ ) const override;

	Size
	size() const;


	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_begin( int resid ) const override;


	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_end( int resid ) const override;


	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_begin( int resid ) const override;


	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_end( int resid ) const override;

	//////////////////// non-const versions

	ResidueNeighborIteratorOP
	neighbor_iterator_begin( int resid ) override;


	ResidueNeighborIteratorOP
	neighbor_iterator_end( int resid ) override;


	ResidueNeighborIteratorOP
	upper_neighbor_iterator_begin( int resid ) override;


	ResidueNeighborIteratorOP
	upper_neighbor_iterator_end( int resid ) override;

	/// @brief Is this PolymerBondedEnergyContainer properly set up for the pose?
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	bool is_valid( core::pose::Pose const &pose ) const;

	utility::vector1< ScoreType > const &
	score_types() const {
		return score_types_;
	}

	bool
	get_computed(core::Size res1, core::Size res2) const {
		std::pair< core::Size, core::Size > edge = std::make_pair(std::min(res1,res2),std::max(res1,res2));

		auto it = computed_.find(edge);
		if ( it == computed_.end() ) { return true; }
		return it->second;
	}

	void
	set_computed(core::Size res1, core::Size res2, bool val) {
		std::pair< core::Size, core::Size > edge = std::make_pair(std::min(res1,res2),std::max(res1,res2));
		computed_[edge] = val;
	}

	core::Real
	get_energy(core::Size res1, core::Size res2, core::Size scoreterm) const {
		std::pair< core::Size, core::Size > edge = std::make_pair(std::min(res1,res2),std::max(res1,res2));

		auto
			it = tables_.find(edge);
		if ( it == tables_.end() ) { return 0.0; }

		return (it->second[scoreterm]);
	}

	void
	set_energy(core::Size res1, core::Size res2, core::Size scoreterm, core::Real E) {
		std::pair< core::Size, core::Size > edge = std::make_pair(std::min(res1,res2),std::max(res1,res2));
		tables_[edge][scoreterm] = E;
	}


private:
	/// @brief Given a pose, set up a list of pairs of residues.
	/// @details The first index is the lower residue (the residue contributing the carboxyl/UPPER_CONNECT) and the second is the
	/// upper residue (the residue contributing the amide/LOWER_CONNECT).
	/// @note fpd: this now includes _all_ chemical edges
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void initialize_peptide_bonded_pair_indices( core::pose::Pose const &pose );

private:
	/// @brief number of chemical connections in the pose
	Size size_;

	/// @brief A map of all chemical edges in the pose
	std::multimap< core::Size, core::Size > chemical_edges_;
	std::set< std::pair< core::Size, core::Size > > chemical_edge_set_;

	/// @brief The vector of score types that this PolymerBondedEnergyContainer will be used to calculate.
	utility::vector1< ScoreType > score_types_;

	/// @brief energy tables. maps each edge to a vector (one for each type)
	std::map< std::pair<core::Size, core::Size> , utility::vector1< core::Real > > tables_;

	/// @brief computed tables. maps each edge to a bool
	std::map< std::pair<core::Size, core::Size> , bool > computed_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	PolymerBondedEnergyContainer();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_PolymerBondedEnergyContainer )
#endif // SERIALIZATION


#endif
