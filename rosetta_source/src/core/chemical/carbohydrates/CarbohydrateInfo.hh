// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    CarbohydrateInfo.hh
/// @brief   Declarations and simple accessor/mutator definitions for CarbohydrateInfo.
/// @author  labonte

#ifndef INCLUDED_core_chemical_carbohydrates_CarbohydrateInfo_HH
#define INCLUDED_core_chemical_carbohydrates_CarbohydrateInfo_HH

// Unit header
#include <core/chemical/carbohydrates/CarbohydrateInfo.fwd.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <iostream>

namespace core {
namespace chemical {
namespace carbohydrates {

/// @remarks  For now, this class functions in a similar manner as RNA_ResidueType, which is not a ResidueType at all,
/// but rather a container of properties specific to RNA.  In the future, I may copy over the functionality in this
/// class to a CarbohydrateResidueType class that derives from ResidueType.  This seems a smarter way to do things.
/// (Thanks, Brian, for the idea!) ~Labonte
class CarbohydrateInfo : public utility::pointer::ReferenceCount {
public:
	// Standard methods ////////////////////////////////////////////////////////
	/// @brief  Empty constructor
	CarbohydrateInfo();

	/// @brief  Standard Constructor
	CarbohydrateInfo(core::chemical::ResidueTypeCOP residue_type);

	/// @brief  Copy constructor
	CarbohydrateInfo(CarbohydrateInfo const & object_to_copy);

	// Assignment operator
	CarbohydrateInfo & operator=(CarbohydrateInfo const & object_to_copy);

	// Destructor
	~CarbohydrateInfo();


	// Standard Rosetta methods ////////////////////////////////////////////////
	/// @brief  Generate string representation of CarbohydrateInfo for debugging purposes.
	void show(std::ostream & output=std::cout) const;

	// Insertion operator (overloaded so that CarbohydrateInfo can be "printed" in PyRosetta).
	friend std::ostream & operator<<(std::ostream & output, CarbohydrateInfo const & object_to_output);


	// Accessors/Mutators
	// Nomenclature
	/// @brief  Return the full IUPAC name of the monosaccharide.
	std::string
	full_name() const {
		return full_name_;
	}


	// Oxidation type
	/// @brief  Return true if the monosaccharide is an aldose.
	bool
	is_aldose() const {
		return is_aldose_;
	}

	/// @brief  Return true if the monosaccharide is a ketose.
	bool
	is_ketose() const {
		return !is_aldose_;
	}


	// Size
	/// @brief  Get the number of carbons in the monosaccharide.
	core::Size
	n_carbons() const {
		return n_carbons_;
	}

	/// @brief  Return true if the monosaccharide is a triose.
	bool
	is_triose() const {
		if (n_carbons_ == 3) {
			return true;
		}
		return false;
	}

	/// @brief  Return true if the monosaccharide is a tetrose.
	bool
	is_tetrose() const {
		if (n_carbons_ == 4) {
			return true;
		}
		return false;
	}

	/// @brief  Return true if the monosaccharide is a pentose.
	bool
	is_pentose() const {
		if (n_carbons_ == 5) {
			return true;
		}
		return false;
	}

	/// @brief  Return true if the monosaccharide is a hexose.
	bool
	is_hexose() const {
		if (n_carbons_ == 6) {
			return true;
		}
		return false;
	}

	/// @brief  Return true if the monosaccharide is a heptose.
	bool
	is_heptose() const {
		if (n_carbons_ == 7) {
			return true;
		}
		return false;
	}


	// Stereochemistry
	/// @brief   Get the stereochemical designation for the monosaccharide.
	/// @return  'L' or 'D'
	char
	stereochem() const {
		return stereochem_;
	}

	/// @brief  Return true if the monosaccharide is an L-sugar.
	bool
	is_L_sugar() const {
		if (stereochem_ == 'L') {
			return true;
		}
		return false;
	}

	/// @brief  Return true if the monosaccharide is a D-sugar.
	bool
	is_D_sugar() const {
		if (stereochem_ == 'D') {
			return true;
		}
		return false;
	}


	// Ring properties
	/// @brief    Get the size of the carbohydrate ring.
	/// @details  A linear monosaccharide has a ring size of zero.
	core::Size
	ring_size() const {
		return ring_size_;
	}

	/// @brief  Return true if the monosaccharide is linear.
	bool
	is_acyclic() const {
		if (ring_size_ == 0) {
			return true;
		}
		return false;
	}

	/// @brief  Return true if the monosaccharide is a ring.
	bool
	is_cyclic() const {
		return !is_acyclic();
	}

	/// @brief  Return true if the monosaccharide is a furanose.
	bool
	is_furanose() const {
		if (ring_size_ == 5) {
			return true;
		}
		return false;
	}

	/// @brief  Return true if the monosaccharide is a pyranose.
	bool
	is_pyranose() const {
		if (ring_size_ == 6) {
			return true;
		}
		return false;
	}

	/// @brief  Return true if the monosaccharide is a septanose.
	bool
	is_septanose() const {
		if (ring_size_ == 7) {
			return true;
		}
		return false;
	}


	// Anomeric form
	/// @brief   Get the anomeric form for the monosaccharide.
	/// @return  "alpha", "beta", or ""
	std::string
	anomer() const {
		return anomer_;
	}

	/// @brief  Return true if the cyclic monosaccharide is an alpha sugar.
	bool
	is_alpha_sugar() const {
		if (anomer_ == "alpha") {
			return true;
		}
		return false;
	}

	/// @brief  Return true if the cyclic monosaccharide is a beta sugar.
	bool
	is_beta_sugar() const {
		if (anomer_ == "beta") {
			return true;
		}
		return false;
	}


	// Polymer
	/// @brief    Return true if the monosaccharide is attached to something at the anomeric carbon.
	/// @remarks  not yet implemented
	bool
	is_glycoside() const {
		return is_glycoside_;
	}

	// TODO: Determine a good way to store connectivity and branching.
	// maybe pairs: (1->4), (1->2), where first designates the anomeric carbon and second the position on the
	// previous residue to which it attaches?
	// But we'll need a forward way also so that a residue knows what is connected to it.


	// Side-chain modifications
	// TODO: Determine a good way to track modifications at each carbon.
	// For now, use booleans.
	/// @brief  Return true if the primary hydroxyl group is oxidized to the acid.
	bool
	is_uronic_acid() const {
		return is_uronic_acid_;
	}

	/// @brief    Return the CHI identifier for the requested nu (internal ring torsion) angle.
	std::pair<core::id::TorsionType, core::Size> nu_id(core::Size subscript) const;

private:
	// Private methods /////////////////////////////////////////////////////////
	// Initialize data members from properties.
	void init(core::chemical::ResidueTypeCOP residue_type);

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data(CarbohydrateInfo object_to_copy_to, CarbohydrateInfo object_to_copy_from);

	// Return the number of carbon atoms (not counting R groups) in the ResidueType.
	core::Size get_n_carbons();

	// Read through all the properties.  Check for impossible cases.  If any property type is not set, the default
	// value will be maintained.
	void read_and_set_properties();

	// If cyclic, define nu angles in terms of CHI ids.
	void define_nu_ids();

	// Private data ////////////////////////////////////////////////////////////
	core::chemical::ResidueTypeCOP residue_type_;
	std::string full_name_;
	bool is_aldose_;
	core::Size n_carbons_;
	char stereochem_;  // L or D
	core::Size ring_size_;  // 0 indicates linear sugar
	std::string anomer_;  // alpha, beta, or null
	bool is_glycoside_;
	bool is_uronic_acid_;

	// Definitions of nu angles in terms of Rosetta 3 CHI angles for this particular sugar.
	// The nu angles should always be the last CHI angles defined in the param file or by a patch file.
	utility::vector1<std::pair<core::id::TorsionType, core::Size> > nu_id_;

};  // class CarbohydrateInfo

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_carbohydrates_CarbohydrateInfo_HH
