// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/qsar/scoring_grid/PCSMultiGridCreator.hh
/// @brief   creator class for nested PCS scoring grid
/// @details last Modified: 05/24/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_qsar_scoring_grid_PCSMultiGridCreator_HH
#define INCLUDED_protocols_qsar_scoring_grid_PCSMultiGridCreator_HH

#include <protocols/qsar/scoring_grid/GridCreator.hh>
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class PCSMultiGridCreator : public GridCreator
{
public:
	virtual GridBaseOP create_grid(utility::tag::TagCOP tag) const;
	virtual GridBaseOP create_grid() const;

	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
	//static std::string grid_name();
};

} // namespace scoring_grid
} // namespace qsar
} // namespace protocols


#endif // INCLUDED_protocols_qsar_scoring_grid_PCSMultiGridCreator_HH
