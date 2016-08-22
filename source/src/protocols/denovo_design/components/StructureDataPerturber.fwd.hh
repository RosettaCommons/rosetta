// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/components/StructureDataPerturber.fwd.hh
/// @brief Classes for altering StructureData objects on the fly
/// @author Tom Linsky (tlinsky@uw.edu)


#ifndef INCLUDED_protocols_denovo_design_components_StructureDataPerturber_fwd_hh
#define INCLUDED_protocols_denovo_design_components_StructureDataPerturber_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

// Forward
namespace protocols {
namespace denovo_design {
namespace components {

class StructureDataPerturber;

typedef utility::pointer::shared_ptr< StructureDataPerturber > StructureDataPerturberOP;
typedef utility::pointer::shared_ptr< StructureDataPerturber const > StructureDataPerturberCOP;
typedef utility::vector1< StructureDataPerturber > StructureDataPerturbers;
typedef utility::vector1< StructureDataPerturberOP > StructureDataPerturberOPs;
typedef utility::vector1< StructureDataPerturberCOP > StructureDataPerturberCOPs;

class ConnectionPerturber;
typedef utility::pointer::shared_ptr< ConnectionPerturber > ConnectionPerturberOP;
typedef utility::pointer::shared_ptr< ConnectionPerturber const > ConnectionPerturberCOP;

class CompoundPerturber;
typedef utility::pointer::shared_ptr< CompoundPerturber > CompoundPerturberOP;
typedef utility::pointer::shared_ptr< CompoundPerturber const > CompoundPerturberCOP;

class HelixPerturber;
typedef utility::pointer::shared_ptr< HelixPerturber > HelixPerturberOP;
typedef utility::pointer::shared_ptr< HelixPerturber const > HelixPerturberCOP;

} //protocols
} //denovo_design
} //components


#endif //INCLUDED_protocols_denovo_design_components_StructureDataPerturber_fwd_hh

