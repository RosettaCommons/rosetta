// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    CarbohydrateInfo.cc
/// @brief   Method definitions for CarbohydrateInfo.
/// @author  labonte

// Unit header
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>

// Package headers
#include <core/chemical/ResidueType.hh>

// Utility headers
#include <utility/PyAssert.hh>
#include <utility/exit.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>


// Construct tracer.
static basic::Tracer TR("core.chemical.carbohydrates.CarbohydrateInfo");


namespace core {
namespace chemical {
namespace carbohydrates {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
// Empty constructor
CarbohydrateInfo::CarbohydrateInfo()
{
	chemical::ResidueTypeCOP residue_type;

	init(residue_type);
}

// Standard constructor
/// @param    <residue_type>: the ResidueType object containing this CarbohydrateInfo
CarbohydrateInfo::CarbohydrateInfo(core::chemical::ResidueTypeCOP residue_type)
{
	init(residue_type);
}

// Copy constructor
CarbohydrateInfo::CarbohydrateInfo(CarbohydrateInfo const & object_to_copy)
{
	copy_data(*this, object_to_copy);
}

// Assignment operator
CarbohydrateInfo &
CarbohydrateInfo::operator=(CarbohydrateInfo const & object_to_copy)
{
	// Abort self-assignment.
	if (this == &object_to_copy) {
		return *this;
	}

	copy_data(*this, object_to_copy);
	return *this;
}

// Destructor
CarbohydrateInfo::~CarbohydrateInfo() {}


// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods
void
CarbohydrateInfo::show(std::ostream & output) const
{
	using namespace std;

	// Parse properties.
	string prefix, suffix, ring_form, modifications;
	if (is_aldose_) {
		prefix = "aldo";
	} else {
		// TODO: Differentiate between 2-ketoses/3-ketoses, etc.
		prefix = "keto";
	}
	switch (n_carbons_) {
		case 3:
			suffix = "triose";
			break;
		case 4:
			suffix = "tetrose";
			break;
		case 5:
			suffix = "pentose";
			break;
		case 6:
			suffix = "hexose";
			break;
		case 7:
			suffix = "heptose";
			break;
	}
	switch (ring_size_) {
		case 5:
			ring_form = "furanose";
			break;
		case 6:
			ring_form = "pyranose";
			break;
		case 7:
			ring_form = "septanose";
			break;
	}
	if (is_uronic_acid_) {
		modifications += string("  uronic acid\n");
	}
	// TODO: Add more modifications.
	if (modifications == "") {
		modifications = "  none";
	}

	// Produce output.
	output << "Carbohydrate Properties for this Residue:" << endl;
	output << " IUPAC Name: " << full_name_ << endl;
	output << " Classification: " << prefix << suffix << endl;
	output << " Stereochemistry: " << stereochem_ << endl;
	if (ring_size_ != 0) {
		output << " Ring Form: " << ring_form << endl;
		output << " Anomeric Form: " << anomer_ << endl;
	}
	output << " Modifications: " << endl << modifications << endl;
}
	

// Accessors/Mutators
// Return the CHI identifier for the requested nu (internal ring torsion) angle.
/// @param    <subscript>: the subscript for nu, which must be between 1 and 2 less than the ring size, inclusive
/// @return	  a pair of values corresponding to the atom tree torsion definitions, in which the first element is
/// either the TorsionID BB or CHI and the second element is an integer
/// @remarks  The atom tree in Rosetta 3 does not allow for rings, so cyclic carbohydrates are implemented as
/// linear residues.  Because of this, the atom tree assigns backbone (BB) torsions to what it considers the main-
/// chain.  Thus, only one side of the ring is considered backbone.  Side-chain (CHI) angles must be defined in
/// the .params file for the residue; they are not automatically assigned.  nu angles, which are the torsion
/// angles defining the ring, not considered BB by the atom tree must therefore be defined as CHI angles in the
/// .params file, even though they are not in reality side-chain torsions.  Since a ring also has a multiplicity
/// of actual side-chains, the indices for those CHI angles that are actually nu angles will vary.
std::pair<core::id::TorsionType, core::Size>
CarbohydrateInfo::nu_id(core::Size subscript) const
{
	assert((subscript > 0) && (subscript <= ring_size_ - 2));
	PyAssert((subscript > 0) && (subscript <= ring_size_ - 2), "CarbohydrateInfo::nu_id(core::Size subscript): nu(subscript) does not have a CHI identifier.");

	return nu_id_[subscript];
}


// Private methods /////////////////////////////////////////////////////////////
// Initialize data members from properties.
void
CarbohydrateInfo::init(core::chemical::ResidueTypeCOP residue_type)
{
	// Set default values.
	residue_type_ = residue_type;
	full_name_ = "";  // TEMP: not yet implemented
	is_aldose_ = true;  // assumes that most sugars will be aldoses if not specified by .params file
	n_carbons_ = get_n_carbons();
	stereochem_ = 'D';  // assumes that most sugars will have D stereochemistry
	ring_size_ = 0;  // assumes linear
	anomer_ = "";  // assumes linear
	is_glycoside_ = true;  // TEMP: not yet implemented
	is_uronic_acid_ = false;
	
	read_and_set_properties();

	define_nu_ids();
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
CarbohydrateInfo::copy_data(
		CarbohydrateInfo object_to_copy_to,
		CarbohydrateInfo object_to_copy_from)
{
	object_to_copy_to.residue_type_ = object_to_copy_from.residue_type_;
	object_to_copy_to.full_name_ = object_to_copy_from.full_name_;
	object_to_copy_to.is_aldose_ = object_to_copy_from.is_aldose_;
	object_to_copy_to.n_carbons_ = object_to_copy_from.n_carbons_;
	object_to_copy_to.stereochem_ = object_to_copy_from.stereochem_;
	object_to_copy_to.ring_size_ = object_to_copy_from.ring_size_;
	object_to_copy_to.anomer_ = object_to_copy_from.anomer_;
	object_to_copy_to.is_glycoside_ = object_to_copy_from.is_glycoside_;
	object_to_copy_to.is_uronic_acid_ = object_to_copy_from.is_uronic_acid_;
	object_to_copy_to.nu_id_ = object_to_copy_from.nu_id_;
}

// Return the number of carbon atoms (not counting R groups) in the ResidueType.
core::Size
CarbohydrateInfo::get_n_carbons()
{
	using namespace std;

	// If we ever add rare sugars larger than 7 carbons, increase the initial i value.
	for (Size carbon_num = 7; carbon_num >= 3; --carbon_num) {
		char carbon_num_char = '0' + carbon_num;  // quick way to convert int to char
		if (residue_type_->has(string(1, 'C') + string(1, carbon_num_char) /*convert chars to strings to concatenate*/)) {
			return carbon_num;
		}
	}
	utility_exit_with_message(
			"This residue is not a sugar or else there is an error in C atom labeling in the .params file.");
	return 0;  // will never be reached
}

// Read through all the properties.  Check for impossible cases.  If any property type is not set, the default
// value will be maintained.
void
CarbohydrateInfo::read_and_set_properties()
{
	using namespace std;
	using namespace utility;

	vector1<string> properties = residue_type_->properties();

	bool aldose_or_ketose_set = false;
	bool stereochem_set = false;
	bool ring_size_set = false;
	bool anomer_set = false;

	for (Size i = 1, n_properties = properties.size(); i <= n_properties; ++i) {
		if (properties[i] == "ALDOSE") {
			if (!is_aldose_) {
				utility_exit_with_message("A sugar cannot be both an aldose and a ketose; check the .param file.");
			} else {
				is_aldose_ = true;
				aldose_or_ketose_set = true;
			}
		} else if (properties[i] == "KETOSE") {
			if (aldose_or_ketose_set && is_aldose_) {
				utility_exit_with_message("A sugar cannot be both an aldose and a ketose; check the .param file.");
			} else {
				is_aldose_ = false;
				aldose_or_ketose_set = true;
			}
		} else if (properties[i] == "L_SUGAR") {
			if (stereochem_set && (stereochem_ == 'D')) {
				utility_exit_with_message("A sugar cannot have both L and D stereochem.; check the .param file.");
			} else {
				stereochem_ = 'L';
				stereochem_set = true;
			}
		} else if (properties[i] == "D_SUGAR") {
			if (stereochem_ == 'L') {
				utility_exit_with_message("A sugar cannot have both L and D stereochem.; check the .param file.");
			} else {
				stereochem_ = 'D';
				stereochem_set = true;
			}
		} else if (properties[i] == "FURANOSE") {
			if (ring_size_set && (ring_size_ != 5)) {
				utility_exit_with_message("A sugar cannot have multiple ring sizes; check the .param file.");
			} else {
				ring_size_ = 5;
				ring_size_set = true;
			}
		} else if (properties[i] == "PYRANOSE") {
			if (ring_size_set && (ring_size_ != 6)) {
				utility_exit_with_message("A sugar cannot have multiple ring sizes; check the .param file.");
			} else {
				ring_size_ = 6;
				ring_size_set = true;
			}
		} else if (properties[i] == "SEPTANOSE") {
			if (ring_size_set && (ring_size_ != 7)) {
				utility_exit_with_message("A sugar cannot have multiple ring sizes; check the .param file.");
			} else {
				ring_size_ = 7;
				ring_size_set = true;
			}
		} else if (properties[i] == "ALPHA_SUGAR") {
			if (anomer_set && (anomer_ == "beta")) {
				utility_exit_with_message("A sugar cannot be both alpha and beta; check the .param file.");
			} else {
				anomer_ = "alpha";
				anomer_set = true;
			}
		} else if (properties[i] == "BETA_SUGAR") {
			if (anomer_set && (anomer_ == "alpha")) {
				utility_exit_with_message("A sugar cannot be both alpha and beta; check the .param file.");
			} else {
				anomer_ = "beta";
				anomer_set = true;
			}
		} else if (properties[i] == "URONIC_ACID") {
			is_uronic_acid_ = true;
		}
	}

	if ((ring_size_ != 0) && (anomer_ == "")) {
		utility_exit_with_message("A cyclic sugar must have its anomeric property declared; check the .param file.");
	}
	if ((ring_size_ == 0) && (anomer_ != "")) {
		utility_exit_with_message("An acyclic sugar cannot be alpha or beta; check the .param file.");
	}
}

// If cyclic, define nu angles in terms of CHI ids.
void
CarbohydrateInfo::define_nu_ids()
{
	using namespace std;
	using namespace utility;
	using namespace id;

	if (ring_size_ != 0) {
		// Get the number of torsions need to define a ring conformation.
		// The two remaining ring torsions (e.g., for a six-membered ring, nu(0) and nu(5)) will have to be determined
		// using vector calculus, because of the cut bond required by the atom tree.
		Size n_torsions_needed = ring_size_ - 2;
		Size n_CHIs = residue_type_->nchi();

		// The final CHIs in the .params file define the (needed) ring torsions.
		Size first_CHI = n_CHIs - n_torsions_needed + 1;

		for (Size i = first_CHI; i <= n_CHIs; ++i) {
			nu_id_.push_back(make_pair(CHI, i));
		}
	}
}


// Friend methods //////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that CarbohydrateInfo can be "printed" in PyRosetta).
std::ostream &
operator<<(std::ostream & output, CarbohydrateInfo const & object_to_output)
{
	object_to_output.show(output);
	return output;
}

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core
