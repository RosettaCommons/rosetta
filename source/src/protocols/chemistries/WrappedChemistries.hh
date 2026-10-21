// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/WrappedChemistries.hh
/// @brief  XML-able Wrappers for the ChemistryBase-derived chemistries in core
/// @author Rocco Moretti

#ifndef INCLUDED_protocols_chemistries_WrappedChemistries_hh
#define INCLUDED_protocols_chemistries_WrappedChemistries_hh

#include <protocols/chemistries/Chemistry.hh>

#include <core/chemical/modifications/ChemistryBase.fwd.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace chemistries {

class WrappedBaseChemistry : public protocols::chemistries::Chemistry {
public:

	WrappedBaseChemistry( std::string const & name,
		core::chemical::modifications::ChemistryBaseOP subchem = nullptr ):
		Chemistry( name ),
		sub_chemistry_( subchem )
	{}

	void apply(core::chemical::MutableResidueType & ) override;

	void apply(core::chemical::MutableResidueType &, core::pose::Pose const & ) override;

	bool
	has_additional_output() const override;

	core::chemical::MutableResidueTypeOP
	get_additional_output() override;

	core::chemical::VDVDMapping
	get_mapping() const override;

	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &) override = 0;

	void set_subchemistry( core::chemical::modifications::ChemistryBaseOP setting ) { sub_chemistry_ = setting; }

private:

	core::chemical::modifications::ChemistryBaseOP sub_chemistry_;

};

///////////////////////////////////////////////////////////////
// Actual usage

class ReprotonateChemistry : public  WrappedBaseChemistry {
public:
	ReprotonateChemistry();

	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &) override;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string class_name();
};

}
}

#endif

