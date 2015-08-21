// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/ligand_docking/HeavyAtomFilter.cc
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/ligand_docking/HeavyAtomFilter.hh>
#include <protocols/ligand_docking/HeavyAtomFilterCreator.hh>


#include <protocols/filters/Filter.hh>
// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/pose/util.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer heavy_atom_tracer( "protocols.ligand_docking.HeavyAtomFilter" );

bool
HeavyAtomFilter::apply( core::pose::Pose const & pose ) const {
	assert(chain_.size()==1 );
	assert(heavy_atom_limit_ >0 );
	core::Size const chain_id= core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size const start = pose.conformation().chain_begin(chain_id);
	core::Size const end = pose.conformation().chain_end(chain_id);

	if ( core::pose::num_heavy_atoms(start,end,pose) > heavy_atom_limit_ ) {
		heavy_atom_tracer<< "Reached heavy atom limit"<< std::endl;
		return false;
	}
	return true;
}

void
HeavyAtomFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{

	if ( tag->getName() != "HeavyAtom" ) {
		heavy_atom_tracer << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( ! (tag->hasOption("chain") && tag->hasOption("heavy_atom_limit") ) ) {
		throw utility::excn::EXCN_RosettaScriptsOption("HeavyAtom filter needs a 'chain' and a 'heavy_atom_limit' option");
	}
	chain_ = tag->getOption<std::string>("chain");
	heavy_atom_limit_ = tag->getOption<core::Size>("heavy_atom_limit");

}

protocols::filters::FilterOP
HeavyAtomFilterCreator::create_filter() const { return protocols::filters::FilterOP( new HeavyAtomFilter ); }

std::string
HeavyAtomFilterCreator::keyname() const { return "HeavyAtom"; }


} // ligand_docking
} // protocols
