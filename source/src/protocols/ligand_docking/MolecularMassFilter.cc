// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/ligand_docking/MolecularMassFilter.cc
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/ligand_docking/MolecularMassFilter.hh>
#include <protocols/ligand_docking/MolecularMassFilterCreator.hh>


#include <protocols/filters/Filter.hh>
// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer mass_tracer( "protocols.ligand_docking.MolecularMassFilter" );

bool
MolecularMassFilter::apply( core::pose::Pose const & pose ) const {
	assert(chain_.size()==1 );
	assert(mass_limit_ >0 );
	core::Size const chain_id= core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size const start = pose.conformation().chain_begin(chain_id);
	core::Size const end = pose.conformation().chain_end(chain_id);

	//core::Real mass=0;


	if(	core::pose::mass(start,end,pose) > mass_limit_ ){
		mass_tracer<< "Reached atom limit"<< std::endl;
		return false;
	}
	return true;
}

void
MolecularMassFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{

	if ( tag->getName() != "MolecularMass" ) {
		mass_tracer << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( ! (tag->hasOption("chain") && tag->hasOption("mass_limit") ) ){
		throw utility::excn::EXCN_RosettaScriptsOption("MolecularMass filter needs a 'chain' and an 'mass_limit' option");
	}
	chain_ = tag->getOption<std::string>("chain");
	mass_limit_ = tag->getOption<core::Size>("mass_limit");


}

protocols::filters::FilterOP
MolecularMassFilterCreator::create_filter() const { return protocols::filters::FilterOP( new MolecularMassFilter ); }

std::string
MolecularMassFilterCreator::keyname() const { return "MolecularMass"; }



} // ligand_docking
} // protocols
