// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   devel/metal_interface/bder/ParseMetalSite.hh
/// @brief  Searches pose for a zinc residue, then fills a vector of MetalSiteResidue objects with info including sequence position of coordinating sidechains, ligand atom xyz, ligand atom name, and atom ids to provide a convenient way for protocols to add metalsite constraints (ligand refers to protein sidechains)
/// @author Bryan Der

#ifndef INCLUDED_devel_metal_interface_ParseMetalSite_HH
#define INCLUDED_devel_metal_interface_ParseMetalSite_HH

#include <devel/metal_interface/ParseMetalSite.fwd.hh>
#include <devel/metal_interface/MetalSiteResidue.fwd.hh> // abbrev. msr
#include <core/pose/Pose.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>


namespace devel {
namespace metal_interface {

class ParseMetalSite : public utility::pointer::ReferenceCount {

public:

	ParseMetalSite();
	ParseMetalSite( core::Size zinc_res );

	virtual ~ParseMetalSite();

	utility::vector1< devel::metal_interface::MetalSiteResidueOP > parse_metalsite( core::pose::Pose const & pose );

	virtual void set_expecting_n_ligands ( Size n );
	virtual bool check_for_parse_error();

private:
	core::Size n_ligands_;
	core::Size zinc_res_;
	bool parse_error_;
	utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr_;

};//end ParseMetalSite

}//namespace metal_interface
}//namespace devel

#endif // INCLUDED_devel_metal_interface_ParseMetalSite_HH
