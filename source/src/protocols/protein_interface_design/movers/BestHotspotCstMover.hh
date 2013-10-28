// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/BestHotspotCstMover.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_BestHotspotCstMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_BestHotspotCstMover_hh

#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/hotspot_hashing/HotspotStubSet.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief remove all HotspotCst's from the pose except the best X
class BestHotspotCstMover : public protocols::moves::Mover
{
public:
	BestHotspotCstMover();
	BestHotspotCstMover(
		protocols::hotspot_hashing::HotspotStubSetOP stub_set,
		core::Size const host_chain,
		core::Size const n_resi
	);
	BestHotspotCstMover(
		BestHotspotCstMover const & init
	);
	protocols::moves::MoverOP clone() const {
		return( protocols::moves::MoverOP( new BestHotspotCstMover( *this ) ) );
	}
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new BestHotspotCstMover ); }
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	void parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	virtual ~BestHotspotCstMover();
private:
	core::Size host_chain_; //where is the stub to be placed
	protocols::hotspot_hashing::HotspotStubSetOP stub_set_;
	core::Real cb_force_constant_;
	core::Size n_resi_; // number of best residues for which to find cst's
};
} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_BestHotspotCstMover_HH*/
