// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ligand_docking/ChainExistsFilter.cc
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/ligand_docking/ChainExistsFilter.hh>
#include <protocols/ligand_docking/ChainExistsFilterCreator.hh>


#include <protocols/filters/Filter.hh>
// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <utility/exit.hh>

#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

static THREAD_LOCAL basic::Tracer atom_tracer( "protocols.ligand_docking.ChainExistsFilter" );

bool
ChainExistsFilter::apply( core::pose::Pose const & pose ) const {
	assert(chain_.size()==1 );
	utility::vector1<core::Size> chain_ids= core::pose::get_chain_ids_from_chain(chain_, pose);

	if ( chain_ids.empty() ) return false;

	return true;
}

void
ChainExistsFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{

	if ( tag->getName() != "ChainExists" ) {
		atom_tracer << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( ! tag->hasOption("chain") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("ChainExists filter needs a 'chain' option");
	}
	chain_ = tag->getOption<std::string>("chain");
}

protocols::filters::FilterOP
ChainExistsFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ChainExistsFilter ); }

std::string
ChainExistsFilterCreator::keyname() const { return "ChainExists"; }


} // ligand_docking
} // protocols
