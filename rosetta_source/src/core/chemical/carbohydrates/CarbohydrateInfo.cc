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
#include <sstream>


// Construct tracer.
static basic::Tracer TR("core.chemical.carbohydrates.CarbohydrateInfo");


namespace core {
namespace chemical {
namespace carbohydrates {

// Define static data.
// If we ever add rare sugars larger than 7 carbons, increase the value.
const core::Size CarbohydrateInfo::MAX_C_SIZE_LIMIT = 7;
const core::Size CarbohydrateInfo::MIN_C_SIZE_LIMIT = 3;


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
	if (is_aldose()) {
		prefix = "aldo";
	} else /*is ketose*/ {
		char num = '0' + anomeric_carbon_;
		prefix = string(1, num) + string("-keto");
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
		modifications = "  none\n";
	}

	// Produce output.
	output << "Carbohydrate Properties for this Residue:" << endl;
	output << " Basic Name: " << base_name() << endl;
	output << " IUPAC Name: " << full_name_ << endl;
	output << " Classification: " << prefix << suffix << endl;
	output << " Stereochemistry: " << stereochem_ << endl;
	if (ring_size_ != 0) {
		output << " Ring Form: " << ring_form << endl;
		output << " Anomeric Form: " << anomer_ << endl;
	}
	output << " Modifications: " << endl << modifications << endl;
	output << " Polymeric Information:" << endl;
	if (mainchain_glycosidic_bond_acceptor_) {
		output << "  Main chain connection: (_->" << mainchain_glycosidic_bond_acceptor_ << ')' << endl;
	} else {
		output << "  Main chain connection: N/A" << endl;
	}
	output << "  Branch connections: " << "branches not yet implemented" << endl;
}


// Accessors/Mutators
// Return the standard/common, non-residue, short name of the monosaccharide.
std::string
CarbohydrateInfo::base_name() const
{
	return root_from_code(residue_type_->name3()) + "ose";
}

// Return the attachment point of the downstream saccharide residue attached to ith branch off of this residue.
/// @param    <i>: the branch point index
/// @return   an integer n of (1->n) of polysaccharide nomenclature, where n specifies the attachment point on the
/// upstream monosaccharide residue; e.g., 4 specifies O4
/// @details  A monosaccharide with a group linked to it at one position is a distinct residue type from the same
/// monosaccharide with the same group linked to it at another position.  For example, Rosetta treats (1->4)-beta-
/// D-glucopyranose as an entirely distinct residue type from (1->3)-beta-D-glucopyranose, with separate .params
/// files for each.\n
/// \n
/// See also:\n
///  CarbohydrateInfo.mainchain_glycosidic_bond_acceptor()\n
///  CarbohydrateInfo.n_branches()
/// @remarks  Branches are not yet implemented.
core::Size
CarbohydrateInfo::branch_point(core::Size i) const
{
	assert((i > 0) && (i <= n_branches()));
	PyAssert((i > 0) && (i <= n_branches()),
			"CarbohydrateInfo::branch_point(core::Size i): "
			"There is no ith branch point on this carbohydrate residue.");

	return branch_points_[i];
}

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
	PyAssert((subscript > 0) && (subscript <= ring_size_ - 2),
			"CarbohydrateInfo::nu_id(core::Size subscript): "
			"nu(subscript) does not have a CHI identifier.");

	return nu_id_[subscript];
}

// Return the BB or CHI identifier for the requested glycosidic linkage torsion angle.
/// @param    <torsion_index>: an integer corresponding to phi (1), psi (2), or omega (3)
/// @return	  a pair of values corresponding to the atom tree torsion definitions, in which the first element is
/// either the TorsionID BB or CHI and the second element is an integer
/// @details  It is crucial to note that this data structure stores information to identify:\n
///  phi(n)\n
///  psi(n+1), NOT psi(n)\n
///  omega(n+1), NOT omega(n)\n
/// \n
/// See Also:\n
///  Pose.phi()\n
///  Pose.set_phi()\n
///  Pose.psi()\n
///  Pose.set_psi()\n
///  Pose.omega()\n
///  Pose.set_omega()
/// @remarks  An enum would be better than an integer for input to this function; however, static constants
/// phi_torsion, psi_torsion, and omega_torsion were already defined in core/id/types.hh.\n
/// The atom tree in Rosetta 3 does not allow for rings, so cyclic carbohydrates are implemented as
/// linear residues.  Because of this, the atom tree assigns backbone (BB) torsions to what it considers the main-
/// chain.  Thus, only one side of the ring is considered backbone.  Side-chain (CHI) angles must be defined in
/// the .params file for the residue; they are not automatically assigned.  Glycosidic linkage torsions are not
/// necessarily defined as main chain torsions by the atom tree, so they must be designated here, in some cases with
/// the use of CHI angles.
std::pair<core::id::TorsionType, core::Size>
CarbohydrateInfo::glycosidic_linkage_id(core::Size torsion_index) const
{
	Size upper_bound = 2;  // for phi and psi
	if (has_exocyclic_linkage_) {
		upper_bound = 3;  // for omega
	}
	assert((torsion_index >= 1) && (torsion_index <= upper_bound));
	PyAssert((torsion_index >= 1) && (torsion_index <= upper_bound),
			"CarbohydrateInfo::glycosidic_linkage_id(core::Size torsion_index): "
			"no defined torsion angle for this index.");

	return glycosidic_linkage_id_[torsion_index];
}


// Private methods /////////////////////////////////////////////////////////////
// Initialize data members from properties.
void
CarbohydrateInfo::init(core::chemical::ResidueTypeCOP residue_type)
{
	// Set default values.
	residue_type_ = residue_type;
	anomeric_carbon_ = 1;  // assumes that most sugars will be aldoses if not specified by .params file
	n_carbons_ = get_n_carbons();
	stereochem_ = 'D';  // assumes that most sugars will have D stereochemistry
	ring_size_ = 0;  // assumes linear
	anomer_ = "";  // assumes linear
	if (residue_type_->is_lower_terminus()){
		is_glycoside_ = false;
	} else {
		is_glycoside_ = true;
	}
	is_uronic_acid_ = false;

	read_and_set_properties();

	determine_polymer_connections();

	determine_IUPAC_names();

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
	object_to_copy_to.short_name_ = object_to_copy_from.short_name_;
	object_to_copy_to.anomeric_carbon_ = object_to_copy_from.anomeric_carbon_;
	object_to_copy_to.n_carbons_ = object_to_copy_from.n_carbons_;
	object_to_copy_to.stereochem_ = object_to_copy_from.stereochem_;
	object_to_copy_to.ring_size_ = object_to_copy_from.ring_size_;
	object_to_copy_to.anomer_ = object_to_copy_from.anomer_;
	object_to_copy_to.is_glycoside_ = object_to_copy_from.is_glycoside_;
	object_to_copy_to.is_uronic_acid_ = object_to_copy_from.is_uronic_acid_;
	object_to_copy_to.nu_id_ = object_to_copy_from.nu_id_;
	object_to_copy_to.mainchain_glycosidic_bond_acceptor_ = object_to_copy_from.mainchain_glycosidic_bond_acceptor_;
	object_to_copy_to.branch_points_ = object_to_copy_from.branch_points_;
	object_to_copy_to.has_exocyclic_linkage_ = object_to_copy_from.has_exocyclic_linkage_;
	object_to_copy_to.glycosidic_linkage_id_ = object_to_copy_from.glycosidic_linkage_id_;
}

// Return the number of carbon atoms (not counting R groups) in the ResidueType.
core::Size
CarbohydrateInfo::get_n_carbons() const
{
	using namespace std;

	for (Size carbon_num = MAX_C_SIZE_LIMIT; carbon_num >= MIN_C_SIZE_LIMIT; --carbon_num) {
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
			if (anomeric_carbon_ != 1) {
				utility_exit_with_message("A sugar cannot be both an aldose and a ketose; check the .param file.");
			} else {
				anomeric_carbon_ = 1;
				aldose_or_ketose_set = true;
			}
		} else if (properties[i] == "KETOSE") {
			if (aldose_or_ketose_set && (anomeric_carbon_ == 1)) {
				utility_exit_with_message("A sugar cannot be both an aldose and a ketose; check the .param file.");
			} else {
				anomeric_carbon_ = 2;  // TODO: Provide method for dealing with non-ulose ketoses.
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

// Get connection data from the residue type.
void
CarbohydrateInfo::determine_polymer_connections()
{
	using namespace std;
	using namespace id;

	if (!residue_type_->is_upper_terminus()) {
		Size upper_atom_index = residue_type_->upper_connect_atom();
		string atom_name = residue_type_->atom_name(upper_atom_index);
		//char atom_number = atom_name[2];
		mainchain_glycosidic_bond_acceptor_ = atoi(&atom_name[2]);
		//Size position = atom_number - '0';
	} else {
		mainchain_glycosidic_bond_acceptor_ = 0;
	}

	// TODO: Implement branching.

	// Exocyclic linkage?
	Size carbons_in_ring = ring_size_ - 1 /*oxygen*/;
	Size last_carbon_in_ring = carbons_in_ring + anomeric_carbon_ - 1;
	if (mainchain_glycosidic_bond_acceptor_ > last_carbon_in_ring) {
		has_exocyclic_linkage_ = true;
	} else {
		has_exocyclic_linkage_ = false;
	}

	// Define phi (phi_torsion = 1 in core/id/types.hh).
	// For aldopyranoses, phi(n) is defined as: O5(n)-C1(n)-OX(n-1)-CX(n-1)
	// BB X+1 is: CX-OX-UPPER1-UPPER2
	// However, CHI 1 is O5-C1-O1-HO1, which for an internal residue with virtual atoms for O1 and HO1, and is
	// the same as phi(n), provided the virtual atoms are made to move with any rotation of BB X+1.
	// The same concept holds for aldofuranoses; however, ketoses are more complicated.  the cyclic oxygen must
	// be the reference for phi, yet CHI 2 at the anomeric position is defined with C1 as the reference atom,
	// not the cyclic oxygen (O5 for furanoses, O6 for pyranoses).  So a unique CHI angle must be defined in
	// .params file for phi.
	// TODO: Correct these. Two virtual atoms in a row in a CHI gives NAN.  Will need to use vector calculus for phi.
	if (is_aldose()) {
		glycosidic_linkage_id_.push_back(make_pair(CHI, 1));
	} else {
		// TODO: Correct this. This is the correct bond but the wrong angle definition.  I need to decide where to
		// this CHI in the .params file.
		glycosidic_linkage_id_.push_back(make_pair(CHI, anomeric_carbon_));
	}

	// Define psi (psi_torsion = 2 in core/id/types.hh).
	// psi(n) is defined as: C(anomeric)(n)-OX(n-1)-CX(n-1)-CX-1(n-1)
	// BB X is: CX-1-CX-OX-UPPER
	// Thus, this is actually the psi angle of the NEXT residue!
	glycosidic_linkage_id_.push_back(make_pair(BB, mainchain_glycosidic_bond_acceptor_));

	// Define omega (omega_torsion = 3 core/id/types.hh).
	if (has_exocyclic_linkage_) {
		glycosidic_linkage_id_.push_back(make_pair(CHI, mainchain_glycosidic_bond_acceptor_ - 1));
	}
}

// Determine and set the full and abbreviated IUPAC names.
// The NAME property in the .params file is actually the standard IUPAC abbreviation (of an internal/unpatched
// residue), not the full name.  It, combined with any patches, is the Rosetta name for the ResidueType.  The IUPAC
// names will change depending on the residue's place in the sequence and/or any patches.
void
CarbohydrateInfo::determine_IUPAC_names()
{
	using namespace std;

	// Determine prefixes.
	stringstream prefixes(stringstream::out);
	if (!residue_type_->is_upper_terminus()) {
		prefixes << "(_->" << mainchain_glycosidic_bond_acceptor_ << ")-";
	}
	if (!residue_type_->is_lower_terminus()) {
		prefixes << anomer_ << '-';
	}
	prefixes << stereochem_ << '-';

	// Determine root.
	string code = residue_type_->name3();
	string root = root_from_code(code);

	// Determine suffix.
	stringstream long_suffix(stringstream::out);
	stringstream short_suffix(stringstream::out);
	switch (ring_size_) {
		case 5:
			long_suffix << "ofuran";
			short_suffix << 'f';
			break;
		case 6:
			long_suffix << "opyran";
			short_suffix << 'p';
			break;
		case 7:
			long_suffix << "oseptan";
			short_suffix << 's';
			break;
	}
	if (residue_type_->is_lower_terminus()) {
		if (is_glycoside_) {
			if (is_uronic_acid_) {
				long_suffix << "uronoside";
				short_suffix << "A";
			} else {
				long_suffix << "oside";
			}
		} else {
			if (is_uronic_acid_) {
				long_suffix << "uronate";
				short_suffix << "A";
			} else {
				long_suffix << "ose";
			}
		}
	} else {
		if (is_uronic_acid_) {
			long_suffix << "uronoyl";
			short_suffix << "A";
		} else {
			long_suffix << "osyl";
		}
		short_suffix << '-';
	}

	full_name_ = prefixes.str() + root + long_suffix.str();
	short_name_ = prefixes.str() + code + short_suffix.str();
}

// If cyclic, define nu angles in terms of CHI ids.
void
CarbohydrateInfo::define_nu_ids()
{
	using namespace std;
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
