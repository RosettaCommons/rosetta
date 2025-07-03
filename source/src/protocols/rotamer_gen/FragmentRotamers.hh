// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rotamer_gen/FragmentRotamers.hh
/// @brief  generate rotamers for a residuetype
/// @author Steven Combs

#ifndef INCLUDED_protocols_rotamer_gen_FragmentRotamers_hh
#define INCLUDED_protocols_rotamer_gen_FragmentRotamers_hh

// Unit Headers
#include <protocols/rotamer_gen/FragmentRotamers.fwd.hh>
#include <protocols/chemistries/Chemistry.hh>

//// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

#ifdef USEBCL
#include <chemistry/bcl_chemistry_sample_conformations.h>
#endif

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace rotamer_gen {


class FragmentRotamers: public protocols::chemistries::Chemistry
{
public:

	FragmentRotamers():
		Chemistry(class_name())
	{}

	void
	parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &) override;

	void
	apply(core::chemical::MutableResidueType & restype) override;

	// We can use the identity mapping for get_mapping()

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

#ifdef USEBCL
::bcl::util::ShPtr< ::bcl::chemistry::SampleConformations> conformation_sampler_;
#endif

}; // class FragmentRotamers

} //namespace rotamer_gen
} //namespace protocols

#endif
