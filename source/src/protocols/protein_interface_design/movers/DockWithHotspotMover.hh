// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Lei Shi <shilei@u.washington.edu>
/// @date 3/19/2013

#ifndef INCLUDED_protocols_protein_interface_design_movers_DockWithHotspotMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_DockWithHotspotMover_hh

#include <protocols/protein_interface_design/movers/DockWithHotspotMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

//parsing
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map

#include <utility/vector1.hh>


// Utility headers

// C++ headers

// Unit headers

namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief A mover to mutate a single residue
class DockWithHotspotMover : public protocols::moves::Mover
{
private:
	typedef protocols::moves::Mover parent;
public:
	/// @brief default ctor
	DockWithHotspotMover();
	virtual ~DockWithHotspotMover();

	virtual void apply( core::pose::Pose & pose ) override;
	virtual protocols::moves::MoverOP clone() const override {
		return (protocols::moves::MoverOP( new protocols::protein_interface_design::movers::DockWithHotspotMover( *this ) ) );
	}
	virtual protocols::moves::MoverOP fresh_instance() const override {
		return protocols::moves::MoverOP( new DockWithHotspotMover );
	}

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	utility::vector1< std::string > hotspot_filenames_;
	utility::vector1< core::Real > hotspot_distcb_weight_;
	core::Real hotspot_score_weight_;
	core::Real centroidscore_filter_;
	core::Real hotspotcst_filter_;
};

} //movers
} //protein_interface_design
} //protocols

#endif //INCLUDED_protocols_protein_interface_design_movers_DockWithHotspotMover_HH_
