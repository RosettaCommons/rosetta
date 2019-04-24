// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzymatic_movers/KinaseMoverCreator.hh
/// @brief  Method declarations for KinaseMoverCreator.
/// @author Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_protocols_enzymatic_movers_KinaseMoverCreator_HH
#define INCLUDED_protocols_enzymatic_movers_KinaseMoverCreator_HH

// Project header
#include <protocols/moves/MoverCreator.hh>


namespace protocols {
namespace enzymatic_movers {

/// @brief  MoverCreator allowing the MoverFactory to create a KinaseMover
class KinaseMoverCreator : public moves::MoverCreator {
public:
	/// @brief  Return an up-casted owning pointer (MoverOP) to the mover.
	protocols::moves::MoverOP create_mover() const override;

	/// @brief  Return the key name for the associated Mover (KinaseMover).
	std::string keyname() const override;

	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

}  // namespace enzymatic_movers
}  // namespace protocols

#endif  // INCLUDED_protocols_enzymatic_movers_KinaseMover_HH
