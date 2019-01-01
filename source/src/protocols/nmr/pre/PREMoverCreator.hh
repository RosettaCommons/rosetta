// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pre/PREMoverCreator.hh
/// @brief   This class will create instances of PREMover for the MoverFactory
/// @details last Modified: 05/30/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_pre_PREMoverCreator_HH
#define INCLUDED_protocols_nmr_pre_PREMoverCreator_HH

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace nmr {
namespace pre {

class PREMoverCreator : public protocols::moves::MoverCreator {
public:
	virtual protocols::moves::MoverOP create_mover() const override;
	virtual std::string keyname() const override;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} // pre
} // nmr
} // protocols

#endif // INCLUDED_protocols_nmr_pre_PREMoverCreator_HH
