// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/AddSidechainConstraintsToHotspots.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_AddSidechainConstraintsToHotspots_hh
#define INCLUDED_protocols_protein_interface_design_movers_AddSidechainConstraintsToHotspots_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class AddSidechainConstraintsToHotspots : public protocols::moves::Mover
{
public:
	AddSidechainConstraintsToHotspots();
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	void parse_my_tag( utility::tag::TagCOP const tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	protocols::moves::MoverOP clone() const { return( protocols::moves::MoverOP( new AddSidechainConstraintsToHotspots( *this ) ) ); }
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new AddSidechainConstraintsToHotspots ); }
	virtual ~AddSidechainConstraintsToHotspots();
	core::Size chain() const;
	void chain( core::Size const c );
	core::Real coord_sdev() const;
	void coord_sdev( core::Real const sdev );
	void add_residue( core::Size const res );
	std::set< core::Size > const & residues() const;
private:
	core::Size chain_; //dflt 2
	core::Real coord_sdev_; //dflt 1.0
	std::set< core::Size > residues_; //dflt empty; used to decide which residues to add constraints to.
};


} //movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_AddSidechainConstraintsToHotspots_HH*/

