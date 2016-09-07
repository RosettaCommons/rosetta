// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/ResidueLoaderOptions.hh
/// @brief forward headers for the Residue Loader Options
/// @author Sam DeLuca

#ifndef INCLUDED_core_chemical_ResidueLoaderOptions_hh
#define INCLUDED_core_chemical_ResidueLoaderOptions_hh

#include <core/chemical/ResidueLoaderOptions.fwd.hh>

//project headers
#include <basic/resource_manager/ResourceOptions.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>


namespace core {
namespace chemical {

class ResidueLoaderOptions : public basic::resource_manager::ResourceOptions
{
public:
	ResidueLoaderOptions();
	~ResidueLoaderOptions() override= default;

	
	void
	parse_my_tag(
		utility::tag::TagCOP tag
	) override;

	
	std::string
	type() const override;

	std::string atom_type_set_tag() const;
	std::string mm_atom_type_set_tag() const;
	std::string element_set_tag() const;
	std::string residue_type_set_tag() const;
	std::string orbital_set_tag() const;

private:
	std::string atom_type_set_tag_;
	std::string mm_atom_type_set_tag_;
	std::string element_set_tag_;
	std::string residue_type_set_tag_;
	std::string orbital_set_tag_;

};


}
}

#endif

