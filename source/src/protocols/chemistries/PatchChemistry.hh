// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/PatchChemistry.hh
/// @brief  A chemistry which applies a patch
/// @author Rocco Moretti

#ifndef INCLUDED_protocols_chemistries_PatchChemistry_hh
#define INCLUDED_protocols_chemistries_PatchChemistry_hh

#include <protocols/chemistries/PatchChemistry.fwd.hh>
#include <protocols/chemistries/Chemistry.hh>

#include <core/chemical/Patch.fwd.hh>
#include <core/chemical/PatchOperation.fwd.hh>

#include <utility/vector1.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace chemistries {

class PatchChemistry : public protocols::chemistries::Chemistry {
public:

	PatchChemistry():
		Chemistry( class_name() )
	{}

	void apply(core::chemical::MutableResidueType & ) override;

	void apply(core::chemical::MutableResidueType & restype, core::pose::Pose const & pose ) override;

	core::chemical::VDVDMapping
	get_mapping() const override;

	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &) override;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string class_name();


	void patch_name( std::string const & setting ) { patch_name_ = setting; }
	void patch_file( std::string const & setting ) { patch_file_ = setting; }
	void add_patch_operation_line( std::string const & line );

private:

	core::chemical::MutableResidueTypeOP
	apply_patches( core::chemical::MutableResidueType const & restype, utility::vector1< core::chemical::PatchCOP > const & patches ) const;

	core::chemical::MutableResidueTypeOP
	apply_patch( core::chemical::MutableResidueType const & restype, core::chemical::Patch const & patch ) const;

	void
	apply_operations( core::chemical::MutableResidueType & restype ) const;

	void
	update_type( core::chemical::MutableResidueType & restype, core::chemical::MutableResidueType const & new_restype );

private:

	/// @brief The name of the patch from the ResidueTypeSet
	std::string patch_name_;

	/// @brief The name of a file from which to get the patch
	std::string patch_file_;

	/// @brief PatchOperations
	utility::vector1< core::chemical::PatchOperationOP > operations_;

	/// @brief The cached mapping of before to after.
	core::chemical::VDVDMapping mapping_;
};

}
}

#endif

