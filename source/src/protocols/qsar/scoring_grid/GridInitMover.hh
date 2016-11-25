// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/GridInitMover.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_GridInitMover_HH_
#define INCLUDED_protocols_qsar_scoring_grid_GridInitMover_HH_


#include <protocols/moves/Mover.hh>
#include <protocols/qsar/scoring_grid/GridInitMover.fwd.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class GridInitMover : public protocols::moves::Mover
{
public:
	GridInitMover();

	virtual ~GridInitMover();
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	// XRW TEMP  virtual std::string get_name() const;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	void apply(core::pose::Pose & pose) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string chain_;

};

}
}
}


#endif /* GRIDINITMOVER_HH_ */
