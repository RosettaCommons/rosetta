// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/StructureDataCreator.hh
/// @brief Creator for StructureData WriteableCacheableData
/// @details
/// @author Tom Linsky

#ifndef INCLUDED_protocols_denovo_design_components_StructureDataCreator_hh
#define INCLUDED_protocols_denovo_design_components_StructureDataCreator_hh

// Unit headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <basic/datacache/WriteableCacheableDataCreator.hh>

namespace protocols {
namespace denovo_design {
namespace components {

class StructureDataCreator : public basic::datacache::WriteableCacheableDataCreator {
public:
	virtual basic::datacache::WriteableCacheableDataOP
	create_data( std::istream & in ) const;

	virtual std::string
	keyname() const;
};

} // namespace components
} // namespace denovo_design
} // namespace protocols

#endif // header guard
