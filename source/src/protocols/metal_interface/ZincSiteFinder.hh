// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/metal_interface/bder/ZincSiteFinder.hh
/// @brief  Searches pose for a zinc residue, then fills a vector of MetalSiteResidue objects with info including sequence position of coordinating sidechains, ligand atom xyz, ligand atom name, and atom ids to provide a convenient way for protocols to add metalsite constraints (ligand refers to protein sidechains)
/// @author Bryan Der

#ifndef INCLUDED_protocols_metal_interface_ZincSiteFinder_HH
#define INCLUDED_protocols_metal_interface_ZincSiteFinder_HH

#include <protocols/metal_interface/ZincSiteFinder.fwd.hh>
#include <protocols/metal_interface/MetalSiteResidue.fwd.hh> // abbrev. msr
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace metal_interface {

class ZincSiteFinder : public utility::pointer::ReferenceCount {

public:

	ZincSiteFinder();
	ZincSiteFinder( core::Size zinc_res );

	virtual ~ZincSiteFinder();

	utility::vector1< protocols::metal_interface::MetalSiteResidueOP > find_zinc_site( core::pose::Pose const & pose );

	virtual void set_expecting_n_ligands ( Size n );
	virtual bool check_for_parse_error();

private:
	core::Size n_ligands_;
	core::Size zinc_res_;
	bool parse_error_;
	utility::vector1< protocols::metal_interface::MetalSiteResidueOP > msr_;

};//end ZincSiteFinder

}//namespace metal_interface
}//namespace protocols

#endif // INCLUDED_protocols_metal_interface_ZincSiteFinder_HH
