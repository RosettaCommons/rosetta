// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/ShoveResidueMover.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_ShoveResidueMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_ShoveResidueMover_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class ShoveResidueMover : public protocols::moves::Mover {
public:
	ShoveResidueMover();
	ShoveResidueMover( core::Size resnum );

	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	void parse_my_tag( utility::tag::TagCOP const tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	protocols::moves::MoverOP clone() const { return new ShoveResidueMover( *this ); }
	protocols::moves::MoverOP fresh_instance() const { return new ShoveResidueMover; }

private:
	bool remove_shove_variant_;
	core::Size resnum_;
	utility::vector1< core::Size > shove_residues_; // a list of residues for which to use the shove_bb atom type, so that backbone atoms might clash.
};


} //movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_ShoveResidueMover_HH*/

