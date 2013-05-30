// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/HotspotDisjointedFoldTreeMover.hh
/// author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_HotspotDisjointedFoldTreeMover_HH
#define INCLUDED_protocols_protein_interface_design_movers_HotspotDisjointedFoldTreeMover_HH

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <set>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class HotspotDisjointedFoldTreeMover : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;

public:
	HotspotDisjointedFoldTreeMover();
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	void parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new HotspotDisjointedFoldTreeMover ); }
	virtual ~HotspotDisjointedFoldTreeMover();

	void add_residue( core::Size const r );
	std::set< core::Size > get_residues() const;
	void chain( core::Size const c );
	core::Size chain() const;
	core::Real ddG_threshold() const;
	void ddG_threshold( core::Real const );
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	void interface_radius( core::Real const rad );
	core::Real interface_radius() const;
	core::kinematics::FoldTreeOP make_disjointed_foldtree( core::pose::Pose const & pose ) const;
private:
	core::Real ddG_threshold_; //dflt 1.0; ala-scan energy above which residues will be considered for the disjointed foldtree
	std::set< core::Size > residues_; // the list of residues to make disjointed
	core::Size chain_; //deflt 2; what is the host chain
	core::Real interface_radius_; //dflt 8 ; value used to look for hotspot residues
	core::scoring::ScoreFunctionOP scorefxn_; //dflt NULL
};


} //movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_HotspotDisjointedFoldTreeMover_HH*/

