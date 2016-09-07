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


namespace core {
namespace scoring {

///////////////////////////////////////////////////////

class PolymerBondedNeighborIterator : public ResidueNeighborIterator
{
	PolymerBondedNeighborIterator & operator = (PolymerBondedNeighborIterator const & src );

public:
	~PolymerBondedNeighborIterator() override;

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
	get_computed(core::Size resnum) const {
		return computed_[resnum];
	}

	void
	set_computed(core::Size resnum, bool val) {
		computed_[resnum] = val;
	}

	core::Real
	get_energy(core::Size resnum, core::Size scoreterm) const {
		return (tables_[resnum][scoreterm]);
	}

	void
	set_energy(core::Size resnum, core::Size scoreterm, core::Real E) {
		tables_[resnum][scoreterm] = E;
	}


private:
	/// @brief Given a pose, set up a list of pairs of peptide-bonded residues.
	/// @details The first index is the lower residue (the residue contributing the carboxyl/UPPER_CONNECT) and the second is the
	/// upper residue (the residue contributing the amide/LOWER_CONNECT).
	/// @note Since the PolymerBondedEnergyContainer is also apparently used for RNA, I'm going to base this on UPPER_CONNECT/LOWER_CONNECT,
	/// not on any peptide-specific geometry.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void initialize_peptide_bonded_pair_indices( core::pose::Pose const &pose );

	/// @brief Logic to get the residue that this residue is connected to at its upper_connect, and which is connected to this residue at its lower_connect.
	/// @details Returns 0 if this residue is not polymeric, if it lacks an upper_connect, if it's not connected to anything at its upper connect, if the
	/// thing that it's connected to is not polymeric, if the thing that it's connected to lacks a lower_connect, or if the thing that it's connected to is not
	/// connected to it at its lower connect.  Phew, what a mouthful!
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::Size other_res_index( core::pose::Pose const &pose, core::Size const this_res_index ) const;

private:
	/// @brief Number of residues in pose (asymmetric case) or asymmetric unit (symmetric case).
	///
	Size size_;

	/// @brief The vector of score types that this PolymerBondedEnergyContainer will be used to calculate.
	///
	utility::vector1< ScoreType > score_types_;

	/// @brief
	///
	utility::vector1< utility::vector1< Real > > tables_;

	/// @brief Has the energy been computed for this residue?
	/// @details Defaults to a vector of "false".  Switched to "true" as energies are computed.
	utility::vector1< bool > computed_;

	/// @brief A map of (residue providing UPPER_CONNECT -> residue providing LOWER_CONNECT).  In the
	/// typical case, this is a map of (lower res -> upper res), though this may not be true for cyclic
	/// geometry.
	std::map <core::Size, core::Size > peptide_bonded_pair_indices_;

	/// @brief Inverse of the map above
	std::map <core::Size, core::Size > inv_peptide_bonded_pair_indices_;

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
