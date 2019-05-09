// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/membrane/AddMembraneMover.hh
/// @brief   Initialize the RosettaMP framework by adding representations of the membrane
///          environemnt to the conformation data held by the pose
///
/// @details Given a pose, configure RosettaMP by adding the following information:
///             1. Create and append a membrane residue (MEM) to the pose
///             2. Create and store a SpanningTopology object
///             3. Setup the initial membrane coordinates (typically centered at the origin)
///             4. (Optional) Initialize per-atom lipid accessibility data
///             5. (Optional) Initialize dimensions of the aqueous pore
///          Upon completion, the call pose.conformation().is_membrane() will return true
///
/// @note    If you add a new step, please document and ensure all data is properly
///          initialized by constructors, parse_my_tag, init_from_cmd, serialization
///          routines, and xsd routines. This class is a data loading mammoth
///
/// @note    Docs last updated: 8/22/17
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @author  JKLeman (julia.koehler.leman@gmail.com)


#ifndef INCLUDED_protocols_membrane_AddMembraneMoverCreator_hh
#define INCLUDED_protocols_membrane_AddMembraneMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {

/// @brief Mover Creator
class AddMembraneMoverCreator : public protocols::moves::MoverCreator {

public:

	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneMoverCreator_hh
