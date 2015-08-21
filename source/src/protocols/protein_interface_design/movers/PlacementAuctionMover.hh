// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/PlacementAuctionMover.hh
/// @brief definition of a class for making the placement auction used by PlaceSimultaneouslyMover
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_PlacementAuctionMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_PlacementAuctionMover_hh

// Project Headers
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/protein_interface_design/movers/PlacementAuctionMover.fwd.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <map>
#include <string>

#include <protocols/hotspot_hashing/HotspotStub.fwd.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.fwd.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class PlacementAuctionMover : public simple_moves::DesignRepackMover
{
public:
	typedef std::pair< protocols::hotspot_hashing::HotspotStubSetOP, std::pair< protocols::hotspot_hashing::HotspotStubOP, core::Size > > StubSetStubPos;
	typedef std::pair< protocols::hotspot_hashing::HotspotStubSetOP, protocols::hotspot_hashing::HotspotStubOP > StubsetStubPair;
	typedef std::pair< core::Real, std::pair< core::Size, StubsetStubPair > > ResidueAuctionItem;
	typedef std::multimap< core::Real, std::pair< core::Size, StubsetStubPair > > ResidueAuction;
	typedef ResidueAuction::iterator iterator;
	typedef ResidueAuction::const_iterator const_iterator;

public:
	PlacementAuctionMover();
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	ResidueAuction auction_results() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const {
		return protocols::moves::MoverOP( new PlacementAuctionMover );
	}
	void insert( ResidueAuctionItem const & item );
	core::Size size() const;
	void erase( iterator it );
	void clear();
	iterator begin();
	iterator end();
	const_iterator begin() const;
	const_iterator end() const;
	virtual void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	//// mutators for Placement movers to copy their internals onto auctionMover
	void host_chain( core::Size const hc );
	void max_cb_cb_dist( core::Real const mccd );
	void cb_force( core::Real const cf );
	void stub_sets( utility::vector1< StubSetStubPos > const & vec );
	std::string get_stub_scorefxn() const;
	utility::vector1< StubSetStubPos > const & stub_sets() const;
	utility::vector1< StubSetStubPos > & stub_sets();
	~PlacementAuctionMover();

private:
	ResidueAuction auction_results_;
	core::Size host_chain_;
	core::Real max_cb_cb_dist_, cb_force_;
	std::string stub_energy_fxn_;
	utility::vector1< StubSetStubPos > stub_sets_;
};

} //movers
} //protein_interface_design
} //protocols

#endif /*INCLUDED_protocols_protein_interface_design_movers_PlacementAuction_HH*/
