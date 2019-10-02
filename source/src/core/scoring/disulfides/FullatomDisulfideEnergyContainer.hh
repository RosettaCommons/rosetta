// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/disulfides/FullatomDisulfideEnergyContainer.hh
/// @brief  Disulfide Energy Container class declaration
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_disulfides_FullatomDisulfideEnergyContainer_hh
#define INCLUDED_core_scoring_disulfides_FullatomDisulfideEnergyContainer_hh

// Unit headers
#include <core/scoring/disulfides/FullatomDisulfideEnergyContainer.fwd.hh>

// Package headers
#include <core/scoring/LREnergyContainer.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>

#include <core/scoring/EnergyMap.fwd.hh>


#include <core/pose/Pose.fwd.hh>

#include <core/scoring/disulfides/DisulfideAtomIndices.fwd.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace disulfides {

class DisulfResNeighbIterator : public ResidueNeighborIterator {
	DisulfResNeighbIterator & operator = (DisulfResNeighbIterator const & );
public:

	DisulfResNeighbIterator(
		FullatomDisulfideEnergyContainer * owner,
		Size focused_node,
		Size disulfide_index
	);

	DisulfResNeighbIterator(
		FullatomDisulfideEnergyContainer * owner
	);

	~DisulfResNeighbIterator() override;

	ResidueNeighborIterator & operator = ( ResidueNeighborIterator const & ) override;
	ResidueNeighborIterator const & operator ++ () override;
	bool operator == ( ResidueNeighborIterator const & ) const override;
	bool operator != ( ResidueNeighborIterator const & ) const override;

	Size upper_neighbor_id() const override;
	Size lower_neighbor_id() const override;

	Size residue_iterated_on() const override;
	Size neighbor_id() const override;

	void save_energy( EnergyMap const & ) override;
	void retrieve_energy( EnergyMap & ) const override;
	void accumulate_energy( EnergyMap & ) const override;

	void mark_energy_computed() override;
	void mark_energy_uncomputed() override;

	bool energy_computed() const override;

private:
	FullatomDisulfideEnergyContainer * owner_;
	Size focused_residue_;
	Size disulfide_index_;

};

class DisulfResNeighbConstIterator : public ResidueNeighborConstIterator {
	DisulfResNeighbConstIterator & operator = (DisulfResNeighbConstIterator const & );
public:
	DisulfResNeighbConstIterator(
		FullatomDisulfideEnergyContainer const * owner,
		Size focused_node,
		Size disulfide_index
	);

	DisulfResNeighbConstIterator( FullatomDisulfideEnergyContainer const * owner );

	~DisulfResNeighbConstIterator() override;

	ResidueNeighborConstIterator & operator = ( ResidueNeighborConstIterator const & ) override;
	ResidueNeighborConstIterator const & operator ++ () override;
	bool operator == ( ResidueNeighborConstIterator const & ) const override;
	bool operator != ( ResidueNeighborConstIterator const & ) const override;

	Size upper_neighbor_id() const override;
	Size lower_neighbor_id() const override;

	Size residue_iterated_on() const override;
	Size neighbor_id() const override;

	void retrieve_energy( EnergyMap & ) const override;
	void accumulate_energy( EnergyMap & ) const override;

	bool energy_computed() const override;

private:
	FullatomDisulfideEnergyContainer const * owner_;
	Size focused_residue_;
	Size disulfide_index_;

};

class FullatomDisulfideEnergyComponents
{
public:
	FullatomDisulfideEnergyComponents() :
		dslf_ss_dst_( 0.0 ),
		dslf_cs_ang_( 0.0 ),
		dslf_ss_dih_( 0.0 ),
		dslf_ca_dih_( 0.0 ),
		dslf_cbs_ds_( 0.0 ),
		dslf_fa13_( 0.0 )
	{}

	Energy dslf_ss_dst() const { return dslf_ss_dst_;}
	Energy dslf_cs_ang() const { return dslf_cs_ang_;}
	Energy dslf_ss_dih() const { return dslf_ss_dih_;}
	Energy dslf_ca_dih() const { return dslf_ca_dih_;}
	Energy dslf_cbs_ds() const { return dslf_cbs_ds_;}
	Energy dslf_fa13() const { return dslf_fa13_;}

	Energy & dslf_ss_dst() { return dslf_ss_dst_;}
	Energy & dslf_cs_ang() { return dslf_cs_ang_;}
	Energy & dslf_ss_dih() { return dslf_ss_dih_;}
	Energy & dslf_ca_dih() { return dslf_ca_dih_;}
	Energy & dslf_cbs_ds() { return dslf_cbs_ds_;}
	Energy & dslf_fa13() { return dslf_fa13_;}

private:

	Energy dslf_ss_dst_;
	Energy dslf_cs_ang_;
	Energy dslf_ss_dih_;
	Energy dslf_ca_dih_;
	Energy dslf_cbs_ds_;
	Energy dslf_fa13_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class FullatomDisulfideEnergyContainer : public LREnergyContainer {
public:
	static Size const NO_DISULFIDE;

public:
	FullatomDisulfideEnergyContainer();

	FullatomDisulfideEnergyContainer( pose::Pose const & );

	void
	update( pose::Pose const & );

	~FullatomDisulfideEnergyContainer() override;

	bool empty() const override;

	LREnergyContainerOP clone() const override;

	void
	set_num_nodes( Size ) override;

	bool
	any_neighbors_for_residue( int resid ) const override;

	bool
	any_upper_neighbors_for_residue( int resid ) const override;

	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_begin( int resid ) const override;

	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_end( int resid ) const override;

	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_begin( int resid ) const override;

	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_end( int resid ) const override;

	ResidueNeighborIteratorOP
	neighbor_iterator_begin( int resid ) override;

	ResidueNeighborIteratorOP
	neighbor_iterator_end( int resid ) override;

	ResidueNeighborIteratorOP
	upper_neighbor_iterator_begin( int resid ) override;

	ResidueNeighborIteratorOP
	upper_neighbor_iterator_end( int resid ) override;

	bool
	disulfide_bonded( Size res1id, Size res2id ) const;

	bool residue_forms_disulfide( Size resid ) const;

	// What residue is resid forming a disulfide with?
	Size other_neighbor_id( Size resid ) const;

	DisulfideAtomIndices const &
	disulfide_atom_indices( Size resid ) const;

	DisulfideAtomIndices const &
	other_neighbor_atom_indices( Size resid ) const;

	// Read and write access granted to Disulf Iterators. Other classes should not use these methods.
public:

	// Mutators
	void save_energy( Size disulfide_index, EnergyMap const & emap );
	void mark_energy_computed( Size disulfide_index );
	void mark_energy_uncomputed( Size disulfide_index );

	// Accessors
	Size lower_neighbor_id( Size disulfide_index ) const;
	Size upper_neighbor_id( Size disulfide_index ) const;
	Size other_neighbor_id( Size disulfide_index, Size resid ) const;

	void accumulate_energy( Size disulfide_index, EnergyMap & emap ) const;
	void retrieve_energy( Size disulfide_index, EnergyMap & emap ) const;
	bool energy_computed( Size disulfide_index ) const;
	Size num_residues() const;


private:
	void find_disulfides( pose::Pose const & pose );
	bool disulfides_changed( pose::Pose const & pose );
	Size num_disulfides() const;

private:
	utility::vector1< Size > resid_2_disulfide_index_;
	utility::vector1< chemical::ResidueTypeCOP > disulfide_residue_types_;
	utility::vector1< std::pair< Size, Size > > disulfide_partners_;
	utility::vector1< std::pair< DisulfideAtomIndices, DisulfideAtomIndices > > disulfide_atom_indices_;
	utility::vector1< std::pair< FullatomDisulfideEnergyComponents, bool > > disulfide_info_;

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
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_disulfides_FullatomDisulfideEnergyContainer )
#endif // SERIALIZATION


#endif
