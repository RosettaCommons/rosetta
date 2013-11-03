// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/HotspotHasherMover.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_HotspotHasherMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_HotspotHasherMover_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
// AUTO-REMOVED #include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class HotspotHasherMover : public protocols::moves::Mover
{
public:
	HotspotHasherMover();
	HotspotHasherMover( std::vector<std::string> const resnames,
			core::scoring::ScoreFunctionCOP scorefxn,
			core::Size const n_stubs,
			core::Size const target_resnum,
			protocols::filters::FilterOP hotspot_filter,
			core::Real const target_distance,
			std::string const hashin_fname,
			std::string const hashout_fname );
	protocols::moves::MoverOP clone() const;
	virtual ~HotspotHasherMover();
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new HotspotHasherMover ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	// SET VARIABLES BASED ON THE COMMAND LINE
	// Residues to use for hashing (defaults to all, sans Gly, Cys, or Pro)
	core::scoring::ScoreFunctionCOP scorefxn_;
	std::vector< std::string > resnames_;
	core::Size n_stubs_, target_resnum_;
	core::Real target_distance_, score_threshold_;
	std::string hashin_fname_, hashout_fname_;
	protocols::filters::FilterOP hotspot_filter_; // a filter for each hotspot. defaults to TrueFilter
};


} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_HotspotHasherMover_HH*/
