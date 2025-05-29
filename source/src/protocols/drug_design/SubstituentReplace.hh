// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/SubstituentReplace.hh
/// @brief use RDKit to replace a subsitutent of a substructure in a ResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_SubstituentReplace_hh
#define INCLUDED_protocols_drug_design_SubstituentReplace_hh

#include <protocols/drug_design/SubstituentReplace.fwd.hh>

#include <protocols/drug_design/substitution_support.fwd.hh>

#include <core/chemical/rdkit/RDKit.fwd.hh>

#include <protocols/chemistries/Chemistry.hh>
#include <core/chemical/AtomRefMapping.hh>

#include <core/chemical/MutableResidueType.fwd.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <string>

namespace protocols {
namespace drug_design {

class SubstituentReplace  : public protocols::chemistries::Chemistry {
public:
	SubstituentReplace();

	void apply( core::chemical::MutableResidueType & ) override;

	/// @brief The file which contains the templates to search for to add to input residue type.
	void template_database( std::string filename, bool append=false );

	/// @brief The file which contains the substituents to add to input residue type.
	void substituents_database( std::string filename, bool append=false );

	/// @brief Will hydrogens in the input be converted to dummy stubs?
	void H_as_dummy( core::Real setting ) { H_as_dummy_ = setting; }
	bool H_as_dummy() const { return H_as_dummy_; }

	/// @brief Will V atoms (Vanadium, but used commonly in Rosetta for "virtual" designations)
	/// in the input be converted to dummy stubs?
	void V_as_dummy( core::Real setting ) { V_as_dummy_ = setting; }
	bool V_as_dummy() const { return V_as_dummy_; }

	/// @brief If not empty, use property weighting based on the given property.
	void weight_by_property( std::string const & setting ) { property_name_ = setting; }
	std::string const & weight_by_property() const { return property_name_; }

	/// @brief Initialize any data members of this instance from an input tag
	/// and a DataMap object
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datacache
	) override;

	core::chemical::VDVDMapping
	get_mapping() const override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	MoleculeSubstitutionOP
	pick_replacement(MoleculeSubstitutionOP template_molsub, unsigned int & tdx_out) const;

private:

	/// @brief The templates to possibly use
	utility::vector1< ::RDKit::ROMolOP > templates_;
	/// @brief The fragments to apply
	utility::vector1< ::RDKit::ROMolOP > substituents_;

	/// @brief Do we consider Hydrogens to be dummy atoms?
	bool H_as_dummy_;
	/// @brief Do we consider V atoms to be dummy atoms?
	bool V_as_dummy_;

	/// @brief If not empty, pick fragments based on the weighting by the given property.
	std::string property_name_;

	core::chemical::VDVDMapping mapping_;

};

} // namespace drug_design
} // namespace protocols

#endif
