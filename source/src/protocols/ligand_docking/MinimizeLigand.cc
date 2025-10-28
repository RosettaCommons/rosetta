// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com), adapted from the ResfileReader code
/// by Steven Lewis (smlewi@gmail.com) and Andrew Leaver-Fay

// Project Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

// Unit Headers
#include <protocols/ligand_docking/MinimizeLigand.hh>
#include <protocols/ligand_docking/ResidueTorsionRestraints.hh>
#include <core/pose/chains_util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <utility/vector1.hh>



//STL headers


namespace protocols {
namespace ligand_docking {

static basic::Tracer MinimizeLigand_tracer( "protocols.ligand_docking.MinimizeLigand", basic::t_debug );

MinimizeLigand::MinimizeLigand():
	//utility::VirtualBase(),
	protocols::moves::Mover("MinimizeLigand")
{
	ligand_torsion_restraints_.clear();
}

MinimizeLigand::MinimizeLigand(std::string const & chain, core::Real degrees):
	chain_(chain), degrees_(degrees)
{
	ligand_torsion_restraints_.clear();
}

MinimizeLigand::MinimizeLigand(MinimizeLigand const & that):
	//utility::VirtualBase(),
	protocols::moves::Mover( that ),
	chain_(that.chain_),
	degrees_(that.degrees_)
{}

MinimizeLigand::~MinimizeLigand() = default;

std::string MinimizeLigand::get_name() const{
	return "MinimizeLigand";
}

void
MinimizeLigand::apply( core::pose::Pose & pose ){
	core::Size chain_id= core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size begin = pose.conformation().chain_begin(chain_id);
	core::Size const end = pose.conformation().chain_end(chain_id);
	for ( ; begin <= end; ++begin ) {
		ligand_torsion_restraints_.push_back(
			utility::pointer::make_shared< protocols::ligand_docking::ResidueTorsionRestraints >(pose, begin, degrees_));
	}
}

/// @brief Remove the constraints added by this mover to the pose in apply, if any.
void
MinimizeLigand::remove_constraints( core::pose::Pose & pose ) {
	for ( core::Size ii(1); ii<= ligand_torsion_restraints_.size(); ++ii ) {
		ligand_torsion_restraints_[ii]->disable(pose);
	}
}

utility::vector1<protocols::ligand_docking::ResidueTorsionRestraintsOP>::iterator
MinimizeLigand::begin(){
	return ligand_torsion_restraints_.begin();
}
utility::vector1<protocols::ligand_docking::ResidueTorsionRestraintsOP>::iterator
MinimizeLigand::end(){
	return ligand_torsion_restraints_.end();
}


} //namespace ligand_docking
} //namespace protocols
