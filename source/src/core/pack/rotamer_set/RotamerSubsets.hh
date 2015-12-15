// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/RotamerSet/RotamerSets.hh
/// @brief  RotamerSets class declaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_rotamer_set_RotamerSubsets_hh
#define INCLUDED_core_pack_rotamer_set_RotamerSubsets_hh

// Unit Headers
#include <core/pack/rotamer_set/RotamerSubsets.fwd.hh>

// Package Headers

#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/task/PackerTask.fwd.hh>

#ifdef WIN32 //VC++ needs full class declaration
#include <core/pack/rotamer_set/RotamerSet.hh> // WIN32 INCLUDE
#endif

// Project Headers
#include <core/conformation/Residue.fwd.hh>

// Utility Headers
#include <utility/vector0.hh>
#include <core/types.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {

//typedef utility::vector1< RotamerSetOP > RotamerSetVector;

class RotamerSubsets : public FixbbRotamerSets
{
public:
	typedef task::PackerTaskCOP PackerTaskCOP;

public:
	RotamerSubsets( FixbbRotamerSets & source, utility::vector0< int > const & rotamer_subset  );
	~RotamerSubsets();

	virtual uint nrotamers() const;
	virtual uint nrotamers_for_moltenres( uint ) const;

	virtual uint total_residue() const;

	virtual uint nmoltenres() const;

	virtual
	uint
	moltenres_2_resid( uint ) const;

	virtual
	uint
	resid_2_moltenres( uint ) const;

	virtual
	uint
	moltenres_for_rotamer( uint ) const;

	virtual
	uint
	res_for_rotamer( uint ) const;

	virtual
	core::conformation::ResidueCOP
	rotamer( uint ) const;

	virtual
	core::conformation::ResidueCOP
	rotamer_for_moltenres( uint moltenres_id, uint rotamerid ) const;

	virtual
	uint
	nrotamer_offset_for_moltenres( uint ) const;

	virtual
	RotamerSetCOP
	rotamer_set_for_residue( uint resid ) const;

	virtual
	RotamerSetOP
	rotamer_set_for_residue( uint resid );

	virtual
	RotamerSetCOP
	rotamer_set_for_moltenresidue( uint moltenresid ) const;

	virtual
	RotamerSetOP
	rotamer_set_for_moltenresidue( uint moltenresid );

	RotamerSetVector::const_iterator begin()
	{ return set_of_rotamer_sets_.begin(); }

	RotamerSetVector::const_iterator end()
	{ return set_of_rotamer_sets_.end(); }

	/// convert rotid in full rotamer enumeration into rotamer id on its source residue
	virtual
	uint
	rotid_on_moltenresidue( uint rotid ) const;

	/// convert moltenres rotid to id in full rotamer enumeration
	virtual
	uint
	moltenres_rotid_2_rotid( uint moltenres, uint moltenresrotid ) const;

	virtual
	void
	show( std::ostream & out ) const;

private:
	void update_offset_data();


protected:

	virtual
	utility::vector1< uint > const &
	resid_2_moltenres_vector() const {
		return resid_2_moltenres_;
	}

	virtual
	utility::vector1< uint > const &
	moltenres_2_resid_vector() const {
		return moltenres_2_resid_;
	}

private:
	uint nmoltenres_;
	uint total_residue_;

	uint nrotamers_;

	RotamerSetVector set_of_rotamer_sets_;
	utility::vector1< uint > resid_2_moltenres_;
	utility::vector1< uint > moltenres_2_resid_;
	utility::vector1< uint > nrotamer_offsets_;

	// originating moltenres for a particular rotamer in the enumeration of all rotamers
	utility::vector1< uint > moltenres_for_rotamer_;
	utility::vector1< uint > nrotamers_for_moltenres_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	RotamerSubsets();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace rotamer_set
} // namespace pack
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_rotamer_set_RotamerSubsets )
#endif // SERIALIZATION


#endif // INCLUDED_core_pack_RotamerSet_RotamerSets_HH
