// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/ResidueLoaderOptions.cc
/// @brief forward headers for the Residue Loader Options
/// @author Sam DeLuca

#include <core/chemical/ResidueLoaderOptions.hh>
#include <core/chemical/ResidueLoaderOptionsCreator.hh>

//utility headers
#include <utility/tag/Tag.hh>

namespace core {
namespace chemical {

std::string ResidueLoaderOptionsCreator::options_type() const
{
	return "ResidueLoaderOptions";
}

basic::resource_manager::ResourceOptionsOP ResidueLoaderOptionsCreator::create_options() const
{
	return new ResidueLoaderOptions;
}


ResidueLoaderOptions::ResidueLoaderOptions() :
		atom_type_set_tag_("fa_standard"),
		mm_atom_type_set_tag_("fa_standard"),
		element_set_tag_("default"),
		residue_type_set_tag_("fa_standard"),
		orbital_set_tag_("fa_standard")
{

}

void ResidueLoaderOptions::parse_my_tag(
	utility::tag::TagCOP tag
)
{

	atom_type_set_tag_ = tag->getOption< std::string >( "atom_type_set_tag", "fa_standard");
	mm_atom_type_set_tag_ = tag->getOption<std::string>("mm_atom_type_set_tag", "fa_standard");
	element_set_tag_ = tag->getOption<std::string>("element_set_tag", "default");
	residue_type_set_tag_ = tag->getOption<std::string>("residue_type_set_tag", "fa_standard");
	orbital_set_tag_ = tag->getOption<std::string>("residue_type_set_tag", "fa_standard");
}

std::string ResidueLoaderOptions::type() const
{
	return "ResidueLoaderOptions";
}

std::string ResidueLoaderOptions::atom_type_set_tag() const
{
	return atom_type_set_tag_;
}

std::string ResidueLoaderOptions::mm_atom_type_set_tag() const
{
	return mm_atom_type_set_tag_;
}

std::string ResidueLoaderOptions::element_set_tag() const
{
	return element_set_tag_;
}

std::string ResidueLoaderOptions::residue_type_set_tag() const
{
	return residue_type_set_tag_;
}

std::string ResidueLoaderOptions::orbital_set_tag() const
{
	return orbital_set_tag_;
}

}
}
