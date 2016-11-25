// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/GridCreator.hh
/// @brief Base class for GridCreators for the Grid load time factory registration scheme
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_GridCreator_HH
#define INCLUDED_protocols_qsar_scoring_grid_GridCreator_HH

#include <core/types.hh>
#include <protocols/qsar/scoring_grid/GridBase.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace qsar {
namespace scoring_grid {

/// @brief Abstract class fora  mover factory.  The creator class is responsible
/// for creating a particular mover class
class GridCreator : public utility::pointer::ReferenceCount
{
public:
	GridCreator();
	virtual ~GridCreator();

	virtual GridBaseOP create_grid(utility::tag::TagCOP tag) const = 0;
	virtual GridBaseOP create_grid() const = 0;
	virtual std::string keyname() const = 0;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const = 0;
private:
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Real weight_;
};

typedef utility::pointer::shared_ptr<GridCreator> GridCreatorOP;
typedef utility::pointer::shared_ptr<GridCreator const> GridCreatorCOP;

}
}
}

#endif /* GRIDCREATOR_HH_ */
