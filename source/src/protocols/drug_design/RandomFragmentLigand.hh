// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/RandomFragmentLigand.hh
/// @brief Randomly frament a ligand to make it smaller
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_RandomFragmentLigand_hh
#define INCLUDED_protocols_drug_design_RandomFragmentLigand_hh

#include <protocols/drug_design/RandomFragmentLigand.fwd.hh>
#include <protocols/chemistries/Chemistry.hh>
#include <core/chemical/AtomRefMapping.hh>

#include <core/chemical/MutableResidueType.fwd.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <string>

namespace protocols {
namespace drug_design {

class RandomFragmentLigand  : public protocols::chemistries::Chemistry {
public:
	RandomFragmentLigand();

	void keep_bigger( bool keep_bigger ) { keep_bigger_ = keep_bigger; }
	void keep_atom( std::string const & keep_atom );
	void ccbond( bool setting ) { ccbond_ = setting; }

	bool keep_bigger() const { return keep_bigger_; }
	std::string keep_atom() const { return keep_atom_; }
	bool ccbond() const { return ccbond_; }

	void apply( core::chemical::MutableResidueType & ) override;

	/// @brief Initialize any data members of this instance from an input tag
	/// and a DataMap object
	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &
	) override;

	core::chemical::VDVDMapping
	get_mapping() const override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	/// @brief Keep the fragment which is bigger.
	bool keep_bigger_;

	/// @brief Keep the fragment which contains the given atom.
	std::string keep_atom_;

	/// @brief Only break C-C bonds
	bool ccbond_;

	core::chemical::VDVDMapping mapping_;
};

} // namespace drug_design
} // namespace protocols

#endif
