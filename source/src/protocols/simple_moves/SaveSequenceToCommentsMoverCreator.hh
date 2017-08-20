// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/SaveSequenceToCommentsMoverCreator.hh
/// @brief  Declaration of the MoverCreator class for the SaveSequenceToCommentsMover
/// @author Christoffer Norn (chnorn@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_SaveSequenceToCommentsMoverCreator_hh
#define INCLUDED_protocols_simple_moves_SaveSequenceToCommentsMoverCreator_hh

// Project headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace simple_moves {

class SaveSequenceToCommentsMoverCreator : public moves::MoverCreator
{
public:
	// XRW TEMP  virtual moves::MoverOP create_mover() const;
	// XRW TEMP  virtual std::string keyname() const;
	// XRW TEMP  static  std::string mover_name();
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

}
}

#endif
