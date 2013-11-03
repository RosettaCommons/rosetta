// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ContingentFilter.hh
/// @brief A filter that is contingent on some other mover to set its pass/fail value
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_StubScoreLoopsFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_StubScoreLoopsFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <protocols/protein_interface_design/filters/StubScoreLoopsFilter.fwd.hh>

#include <protocols/simple_filters/ConstraintScoreCutoffFilter.hh>

#include <protocols/hotspot_hashing/HotspotStub.fwd.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.fwd.hh>
#include <utility/vector1.hh>


// Unit headers

namespace protocols {
namespace protein_interface_design{
namespace filters {

class StubScoreLoopsFilter : public protocols::simple_filters::ConstraintScoreCutoffFilter
{
private:
	typedef protocols::simple_filters::ConstraintScoreCutoffFilter Parent;
	typedef std::pair< protocols::hotspot_hashing::HotspotStubSetOP, std::pair< protocols::hotspot_hashing::HotspotStubOP, core::Size > > StubSetStubPos;
public:
	/// @brief default ctor
	StubScoreLoopsFilter();

	virtual ~StubScoreLoopsFilter();

	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
private:
	hotspot_hashing::HotspotStubSetOP stub_set_;
	std::string resfile_;
	core::Real cb_force_;
	core::Size loop_start_;
	core::Size loop_stop_;
};

} // filters
} //protein_interface_design
} // protocols

#endif //INCLUDED_protocols_Filters_StubScoreLoopsFilter_HH_
