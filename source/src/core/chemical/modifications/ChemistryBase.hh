// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/modifications/ChemistryBase.hh
/// @brief  Base class for chemical modifications such as add atoms, delete atoms, replace atom, hydrogen manipulationsds
/// @author Steven Combs

#ifndef INCLUDED_core_chemical_modifications_ChemistryBase_hh
#define INCLUDED_core_chemical_modifications_ChemistryBase_hh


#include <utility/VirtualBase.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/modifications/ChemistryBase.fwd.hh>
#include <core/chemical/AtomRefMapping.hh>

#include <core/types.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <string>
#include <map>

namespace core {
namespace chemical {
namespace modifications {

/// @brief What was the status of the last call to apply (or get_additional_output)
enum ChemistryStatus {
	SUCCESS = 0, // Successfully modified the ResidueType
	FAIL_RETRY, // Did not succeed, but could possibly do so if apply() was called again.
	FAIL_DO_NOT_RETRY, // Something is off with the Chemistry/ResidueType pairing, but a different ResidueType might work
	FAIL_BAD_INPUT, // Something is wrong with the input ResidueType - future calls will always fail
	FAIL_BAD_SETTINGS // Something is wrong with the Chemistry object itself - it will never work on any ResidueType
	// FAIL_BAD_SETTINGS should probably be a utility_exit instead ...
};

class ChemistryBase : public utility::VirtualBase
{
public:

	ChemistryBase(std::string const & name):
		name_(name),
		last_status_( SUCCESS )
	{}

	/// @brief Return the name of this Chemistry object
	std::string name() const { return name_; }

	/// @brief Modify the passed ResidueType
	virtual void apply(MutableResidueType & ) = 0;

	/// @brief Are there alternate ResidueTypes which are availible from the last time we called apply?
	/// (That is, will get_addtional_output() return non-null?)
	virtual
	bool
	has_additional_output() const { return false; }

	/// @brief Get additional generated ResidueTypes, if any.
	/// This allows for 1-to-many Chemistries
	virtual
	core::chemical::MutableResidueTypeOP
	get_additional_output() {
		if ( get_last_status() == SUCCESS ) {
			set_last_status( FAIL_RETRY );
		}
		return nullptr;
	}

	/// @brief Get the vertex mapping that was used for the last apply() or get_additional_output()
	/// This is a mapping FROM the vds in the BEFORE MutableResidueType TO the vds in the AFTER MutableResidueType.
	/// The base class implementation defaults to an identity mapping.
	virtual
	VDVDMapping
	get_mapping() const;

	/// @brief What was the status of the last call to apply()/get_additional_output()
	ChemistryStatus
	get_last_status() const { return last_status_; }

	/// @brief Set the status of the chemistry object
	void
	set_last_status( ChemistryStatus setting ) { last_status_ = setting; }

private:

	std::string name_;

	ChemistryStatus last_status_;
};


}
}
}

#endif
