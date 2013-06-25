// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/carbohydrates/CarbohydrateInfo.hh
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
#include <map>
#include <string>
#include <iostream>


namespace core {
namespace chemical {
namespace carbohydrates {

class CarbohydrateInfo : public utility::pointer::ReferenceCount {
public:
	// Standard methods ////////////////////////////////////////////////////////
	/// @brief  Empty constructor
	CarbohydrateInfo();

	/// @brief  Standard constructor
	CarbohydrateInfo(core::chemical::ResidueTypeCAP residue_type);

	/// @brief  Copy constructor
	CarbohydrateInfo(CarbohydrateInfo const & object_to_copy);

	// Assignment operator
	CarbohydrateInfo & operator=(CarbohydrateInfo const & object_to_copy);

	// Destructor
	~CarbohydrateInfo();


	// Standard Rosetta methods ////////////////////////////////////////////////
	/// @brief  Generate string representation of CarbohydrateInfo for debugging purposes.
	virtual void show(std::ostream & output=std::cout) const;


	// Static constant data access
	/// @brief A list of allowed properties in carbohydrate .params and patch files.
	static utility::vector1<std::string> const & sugar_properties();

	/// @brief A list of Rosetta PDB 3-letter codes for saccharide residues mapped to the corresponding root.
	static std::map<std::string, std::string> const & code_to_root_map();

	//static std::map<std::string, std::string> const CODE_TO_ROOT_MAP;


	// Accessors/Mutators
	// Nomenclature
	/// @brief  Return the full IUPAC name of the monosaccharide.
	std::string
	full_name() const
	{
		return full_name_;
	}

	/// @brief  Return the abbreviated IUPAC name of the monosaccharide (for use in polysaccharide sequences).
	std::string
	short_name() const
	{
		return short_name_;
	}

	/// @brief  Return the standard/common, non-residue, short name of the monosaccharide.
	std::string	base_name() const;


	// Oxidation type
	/// @brief    Return true if the monosaccharide is an aldose.
	/// @details  An aldose sugar is an aldehyde derivative.
	/// See also:\n
	///  CarbohydrateInfo.is_ketose()
	///  CarbohydrateInfo.is_ulose()
	///  CarbohydrateInfo.anomeric_carbon()
	bool
	is_aldose() const
	{
		return (anomeric_carbon_ == 1);
	}

	/// @brief    Return true if the monosaccharide is a ketose.
	/// @details  A ketose sugar is a ketone derivative.\n
	/// Does not distinguish between 2-ketoses (uloses) and 3-ketoses.\n
	/// \n
	/// See also:\n
	///  CarbohydrateInfo.is_aldose()
	///  CarbohydrateInfo.is_ulose()
	///  CarbohydrateInfo.anomeric_carbon()
	bool
	is_ketose() const
	{
		return (anomeric_carbon_ != 1);
	}

	/// @brief    Return true if the monosaccharide is an n-ketose.
	/// @details  A ketose sugar is a ketone derivative.\n
	/// See also:\n
	///  CarbohydrateInfo.is_aldose()
	///  CarbohydrateInfo.is_ulose()
	///  CarbohydrateInfo.anomeric_carbon()
	bool
	is_ketose(core::uint const n) const
	{
		return (anomeric_carbon_ == n);
	}

	/// @brief    Return true if the monosaccharide is a 2-ketose.
	/// @details  See also:\n
	///  CarbohydrateInfo.is_aldose()
	///  CarbohydrateInfo.is_ketose()
	///  CarbohydrateInfo.anomeric_carbon()
	bool
	is_ulose() const
	{
		return is_ketose(2);
	}

	/// @brief    Return the anomeric carbon number.
	/// @details  For linear monosaccharides, this number corresponds to the carbon that is oxidized to the aldehyde
	/// or ketone.
	/// See also:\n
	///  CarbohydrateInfo.anomeric_carbon_name()
	///  CarbohydrateInfo.is_aldose()
	///  CarbohydrateInfo.is_ketose()
	///  CarbohydrateInfo.is_ulose()
	core::uint
	anomeric_carbon() const
	{
		return anomeric_carbon_;
	}

	/// @brief    Return the atom name of the anomeric carbon.
	/// @details  See also:\n
	///  CarbohydrateInfo.anomeric_carbon()
	///  CarbohydrateInfo.is_aldose()
	///  CarbohydrateInfo.is_ketose()
	///  CarbohydrateInfo.is_ulose()
	std::string
	anomeric_carbon_name() const
	{
		return anomeric_carbon_name_;
	}


	/// @brief    Return the cyclic oxygen number or 0, if linear.
	/// @details  This atom is used as a reference atom for certain torsion angles.\n
	/// See also:\n
	///  CarbohydrateInfo.anomeric_carbon()
	///  CarbohydrateInfo.cyclic_oxygen_name()
	core::uint
	cyclic_oxygen() const
	{
		return cyclic_oxygen_;
	}

	/// @brief    Return the atom name of the cyclic oxygen.
	/// @details  This atom is used as a reference atom for certain torsion angles.\n
	/// See also:\n
	///  CarbohydrateInfo.cyclic_oxygen()
	///  CarbohydrateInfo.anomeric_carbon_name()
	std::string
	cyclic_oxygen_name() const
	{
		return cyclic_oxygen_name_;
	}


	// Size
	/// @brief    Get the number of carbons in the monosaccharide.
	/// @details  This ignores carbons found in modifications to the base sugar.
	core::Size
	n_carbons() const
	{
		return n_carbons_;
	}

	/// @brief  Return true if the monosaccharide is a triose.
	bool
	is_triose() const
	{
		return (n_carbons_ == 3);
	}

	/// @brief  Return true if the monosaccharide is a tetrose.
	bool
	is_tetrose() const
	{
		return (n_carbons_ == 4);
	}

	/// @brief  Return true if the monosaccharide is a pentose.
	bool
	is_pentose() const
	{
		return (n_carbons_ == 5);
	}

	/// @brief  Return true if the monosaccharide is a hexose.
	bool
	is_hexose() const
	{
		return (n_carbons_ == 6);
	}

	/// @brief  Return true if the monosaccharide is a heptose.
	bool
	is_heptose() const
	{
		return (n_carbons_ == 7);
	}

	/// @brief    Return true if the monosaccharide is an octose.
	bool
	is_octose() const
	{
		return (n_carbons_ == 8);
	}

	/// @brief    Return true if the monosaccharide is a nonose.
	bool
	is_nonose() const
	{
		return (n_carbons_ == 9);
	}


	// Stereochemistry
	/// @brief   Get the stereochemical designation for the monosaccharide.
	/// @return  'L' or 'D'
	char
	stereochem() const
	{
		return stereochem_;
	}

	/// @brief  Return true if the monosaccharide is an L-sugar.
	bool
	is_L_sugar() const
	{
		return (stereochem_ == 'L');
	}

	/// @brief  Return true if the monosaccharide is a D-sugar.
	bool
	is_D_sugar() const
	{
		return (stereochem_ == 'D');
	}


	// Ring properties
	/// @brief    Get the size of the carbohydrate ring.
	/// @details  A linear monosaccharide has a ring size of zero.
	core::Size
	ring_size() const
	{
		return ring_size_;
	}

	/// @brief  Return true if the monosaccharide is linear.
	bool
	is_acyclic() const
	{
		return (ring_size_ == 0);
	}

	/// @brief  Return true if the monosaccharide is a ring.
	bool
	is_cyclic() const
	{
		return !is_acyclic();
	}

	/// @brief    Return true if the monosaccharide is a furanose.
	/// @details  A furanose has a five-membered ring (like furan).
	bool
	is_furanose() const
	{
		return (ring_size_ == 5);
	}

	/// @brief    Return true if the monosaccharide is a pyranose.
	/// @details  A pyranose has a six-membered ring (like pyran).
	bool
	is_pyranose() const
	{
		return (ring_size_ == 6);
	}

	/// @brief    Return true if the monosaccharide is a septanose.
	/// @details  A septanose has a seven-membered ring.
	bool
	is_septanose() const
	{
		return (ring_size_ == 7);
	}


	// Anomeric form
	/// @brief    Get the anomeric form for the monosaccharide.
	/// @return   "alpha", "beta", or ""
	/// @details  "alpha" and "beta" designate the stereochemistry at the anomeric carbon of a cyclic sugar.
	std::string
	anomer() const
	{
		return anomer_;
	}

	/// @brief    Return true if the cyclic monosaccharide is an alpha sugar.
	/// @details  "alpha" and "beta" designate the stereochemistry at the anomeric carbon of a cyclic sugar.
	bool
	is_alpha_sugar() const
	{
		return (anomer_ == "alpha");
	}

	/// @brief    Return true if the cyclic monosaccharide is a beta sugar.
	/// @details  "alpha" and "beta" designate the stereochemistry at the anomeric carbon of a cyclic sugar.
	bool
	is_beta_sugar() const
	{
		return (anomer_ == "beta");
	}


	// Polymer info
	/// @brief    Return true if the monosaccharide is attached to something at the anomeric carbon.
	/// @remarks  not yet fully implemented
	bool
	is_glycoside() const
	{
		return is_glycoside_;
	}

	/// @brief    Return the attachment point of the downstream saccharide residue of the main chain.
	/// @return   an integer n of (1->n) of polysaccharide nomenclature, where n specifies the attachment point on the
	/// upstream monosaccharide residue; e.g., 4 specifies O4; n = 0 specifies an upper terminus
	/// @details  A monosaccharide with a group linked to it at one position is a distinct residue type from the same
	/// monosaccharide with the same group linked to it at another position.  For example, Rosetta treats (1->4)-beta-
	/// D-glucopyranose as an entirely distinct residue type from (1->3)-beta-D-glucopyranose, with separate .params
	/// files for each.\n
	/// \n
	/// See also:\n
	///  CarbohydrateInfo.n_branches()\n
	///  CarbohydrateInfo.branch_point()
	core::uint
	mainchain_glycosidic_bond_acceptor() const
	{
		return mainchain_glycosidic_bond_acceptor_;
	}

	/// @brief    Return the number of branches off of this residue.
	/// @details  A monosaccharide with a group linked to it at one position is a distinct residue type from the same
	/// monosaccharide with the same group linked to it at another position.  For example, Rosetta treats (1->4)-beta-
	/// D-glucopyranose as an entirely distinct residue type from (1->3)-beta-D-glucopyranose, with separate .params
	/// files for each.\n
	/// \n
	/// See also:\n
	///  CarbohydrateInfo.mainchain_glycosidic_bond_acceptor()\n
	///  CarbohydrateInfo.branch_point()
	core::Size
	n_branches() const
	{
		return branch_points_.size();
	}

	/// @brief  Return the attachment point of the downstream saccharide residue attached to ith branch off of this
	/// residue.
	core::uint branch_point(core::uint i) const;

	/// @brief  Return true if the attachment point of the downstream saccharide is on an exocyclic carbon.
	bool
	has_exocyclic_linkage() const
	{
		return has_exocyclic_linkage_;
	}


	// Side-chain modifications
	// (Be sure to update database/chemical/carbohydrates/sugar_properties.list if adding functions to this section.)
	/// @brief  Return true if any hydroxyl group has been modified to an acetylated amino group.
	bool
	is_N_acetylated() const {
		return modifications_.contains("acetylamino sugar");
	}

	/// @brief  Return true if any hydroxyl group has been modified by acetylation.
	bool
	is_O_acetylated() const {
		return modifications_.contains("acetyl sugar");
	}

	/// @brief  Return true if the sugar has been acetylated at any position.
	bool
	is_acetylated() const {
		return is_N_acetylated() || is_O_acetylated();
	}

	/// @brief  Return true if any hydroxyl group has been modified to an amino group or an acetylated amino group.
	bool
	is_amino_sugar() const {
		return modifications_.contains("amino sugar") || modifications_.contains("acetylamino sugar");
	}

	/// @brief  Return true if the primary hydroxyl group is oxidized to the acid.
	bool
	is_uronic_acid() const {
		return modifications_.contains("uronic acid");
	}


	// Torsion angle mappings
	/// @brief  Return the BB or CHI identifier for the requested glycosidic linkage torsion angle.
	std::pair<core::id::TorsionType, core::uint> glycosidic_linkage_id(core::uint torsion_index) const;

	/// @brief  Return the CHI identifier for the requested nu (internal ring torsion) angle.
	std::pair<core::id::TorsionType, core::uint> nu_id(core::uint subscript) const;

private:
	// Private methods /////////////////////////////////////////////////////////
	// Initialize data members from properties.
	void init(core::chemical::ResidueTypeCAP residue_type);

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data(CarbohydrateInfo object_to_copy_to, CarbohydrateInfo object_to_copy_from);

	// Return the number of carbon atoms (not counting R groups) in the ResidueType.
	core::Size get_n_carbons() const;

	// Read through all the properties.  Check for impossible cases.  If any property type is not set, the default
	// value will be maintained.
	void read_and_set_properties();

	// Get connection data from the residue type.
	void determine_polymer_connections();

	// Determine and set the full and abbreviated IUPAC names.
	void determine_IUPAC_names();

	// Get monosaccharide root name from the 3-letter code.
	std::string
	root_from_code(std::string code) const
	{
		return code_to_root_map().find(code)->second;  // operator[] is not overloaded for const maps.
	}

	// If cyclic, define nu angles in terms of CHI ids.
	void define_nu_ids();

	// Private data ////////////////////////////////////////////////////////////
	core::chemical::ResidueTypeCAP residue_type_;
	std::string full_name_;
	std::string short_name_;
	core::uint anomeric_carbon_;  // also indicative of location of aldehyde/ketone oxidation
	std::string anomeric_carbon_name_;  // string for quick reference
	core::uint cyclic_oxygen_;  // 0 if linear
	std::string cyclic_oxygen_name_;  // string for quick reference
	core::Size n_carbons_;
	char stereochem_;  // L or D
	core::Size ring_size_;  // 0 indicates linear
	std::string anomer_;  // alpha, beta, or null
	bool is_glycoside_;
	utility::vector1<std::string> modifications_;  // indexed by position

	// Glycosidic bond attachment points, i.e., the second integer in (1->n) notations.
	core::uint mainchain_glycosidic_bond_acceptor_;  // 0 if N/A, i.e., if residue type is an upper terminus
	utility::vector1<core::uint> branch_points_;
	bool has_exocyclic_linkage_;

	// Torsion angle mappings
	// Definitions of phi, psi, omega, and nu angles in terms of Rosetta 3 BB and CHI angles for this particular
	// sugar.  chi angles directly correspond to CHI torsions.  (Other torsion angles further up the side chains can
	// be accessed and set with CHI torsions but do not have an official designation, so they are not mapped.)
	utility::vector1<std::pair<core::id::TorsionType, core::uint> > glycosidic_linkage_id_;

	// The nu angles should always be the last CHI angles defined in the param file or by a patch file.
	utility::vector1<std::pair<core::id::TorsionType, core::uint> > nu_id_;

	// Constants.
	static core::Size const MAX_C_SIZE_LIMIT;  // maximum size of a carbohydrate carbon chain in Rosetta
	static core::Size const MIN_C_SIZE_LIMIT;
};  // class CarbohydrateInfo

// Insertion operator (overloaded so that CarbohydrateInfo can be "printed" in PyRosetta).
std::ostream & operator<<(std::ostream & output, CarbohydrateInfo const & object_to_output);



}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_carbohydrates_CarbohydrateInfo_HH
