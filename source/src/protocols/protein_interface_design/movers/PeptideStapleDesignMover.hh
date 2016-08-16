// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/PeptideStapleDesignMover.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_PeptideStapleDesignMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_PeptideStapleDesignMover_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <protocols/simple_moves/PeptideStapleMover.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief Introduces a peptide staple (ala G. Verdine) to the pose.
class PeptideStapleDesignMover : public protocols::moves::Mover
{
public:
	PeptideStapleDesignMover();
	PeptideStapleDesignMover( core::Size const seqpos, core::Size const staple_gap );
	PeptideStapleDesignMover( PeptideStapleDesignMover const & init );
	virtual ~PeptideStapleDesignMover();

	protocols::moves::MoverOP clone() const;
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new PeptideStapleDesignMover ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	protocols::simple_moves::PeptideStapleMoverOP stapler_;
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_PeptideStapleDesignMover_HH*/
