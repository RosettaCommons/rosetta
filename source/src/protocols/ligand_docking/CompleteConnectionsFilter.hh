// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/ligand_docking/CompleteConnectionsFilter.hh
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_ligand_docking_CompleteConnectionsFilter_hh
#define INCLUDED_protocols_ligand_docking_CompleteConnectionsFilter_hh


#include <core/types.hh>
#include <protocols/filters/Filter.hh>

#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

class CompleteConnectionsFilter : public protocols::filters::Filter
{
public:
	CompleteConnectionsFilter() :
		//utility::pointer::ReferenceCount(),
		protocols::filters::Filter( "CompleteConnections" ) {}

	CompleteConnectionsFilter(std::string chain) :
		//utility::pointer::ReferenceCount(),
		protocols::filters::Filter( "CompleteConnections" ),
		chain_(chain)
	{}

	CompleteConnectionsFilter( CompleteConnectionsFilter const & init ) :
		//utility::pointer::ReferenceCount(),
		protocols::filters::Filter( init ),
		chain_(init.chain_)
	{};

	virtual ~CompleteConnectionsFilter(){};

	bool apply( core::pose::Pose const & pose ) const;
	protocols::filters::FilterOP clone() const {
		return protocols::filters::FilterOP( new CompleteConnectionsFilter( *this ) );
	}

	protocols::filters::FilterOP fresh_instance() const{
		return protocols::filters::FilterOP( new CompleteConnectionsFilter() );
	}

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & reference_pose );

private:
	std::string chain_;
};
} // ligand_docking
} // protocols

#endif //INCLUDED_protocols_ProteinInterfaceDesign_ligand_docking_CompleteConnectionsFilter_HH_
