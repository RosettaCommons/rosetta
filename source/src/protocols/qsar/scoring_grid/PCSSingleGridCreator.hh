// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/qsar/scoring_grid/PCSSingleGridCreator.hh
/// @brief   creator class for simple PCS scoring grid
/// @details last Modified: 05/17/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_qsar_scoring_grid_PCSSingleGridCreator_HH
#define INCLUDED_protocols_qsar_scoring_grid_PCSSingleGridCreator_HH

#include <protocols/qsar/scoring_grid/GridCreator.hh>
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class PCSSingleGridCreator : public GridCreator
{
public:
	GridBaseOP create_grid(utility::tag::TagCOP tag) const override;
	GridBaseOP create_grid() const override;

	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
	//static std::string grid_name();
};

} // namespace scoring_grid
} // namespace qsar
} // namespace protocols


#endif // INCLUDED_protocols_qsar_scoring_grid_PCSSingleGridCreator_HH
