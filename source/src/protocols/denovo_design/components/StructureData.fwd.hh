// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/denovo_design/components/StructureData.fwd.hh
/// @brief  StructureData forward header
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_components_StructureData_fwd_hh
#define INCLUDED_protocols_denovo_design_components_StructureData_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace denovo_design {
namespace components {

// Forward
class StructureData;

// Pointer Types
typedef utility::pointer::shared_ptr< StructureData > StructureDataOP;
typedef utility::pointer::shared_ptr< StructureData const > StructureDataCOP;

typedef utility::pointer::weak_ptr< StructureData > StructureDataAP;
typedef utility::pointer::weak_ptr< StructureData const > StructureDataCAP;

// Forward
class SingleChainStructureData;

// Pointer Types
typedef utility::pointer::shared_ptr< SingleChainStructureData > SingleChainStructureDataOP;
typedef utility::pointer::shared_ptr< SingleChainStructureData const > SingleChainStructureDataCOP;

typedef utility::pointer::weak_ptr< SingleChainStructureData > SingleChainStructureDataAP;
typedef utility::pointer::weak_ptr< SingleChainStructureData const > SingleChainStructureDataCAP;

// Forward
class MultiChainStructureData;

// Pointer Types
typedef utility::pointer::shared_ptr< MultiChainStructureData > MultiChainStructureDataOP;
typedef utility::pointer::shared_ptr< MultiChainStructureData const > MultiChainStructureDataCOP;

typedef utility::pointer::weak_ptr< MultiChainStructureData > MultiChainStructureDataAP;
typedef utility::pointer::weak_ptr< MultiChainStructureData const > MultiChainStructureDataCAP;

} // namespace components
} // namespace denovo_design
} // namespace protocols

#endif
