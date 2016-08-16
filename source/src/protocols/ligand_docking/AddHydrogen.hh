// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/AddHydrogen.hh
///
/// @brief
/// @author Gordon Lemmon


#ifndef INCLUDED_protocols_ligand_docking_AddHydrogen_hh
#define INCLUDED_protocols_ligand_docking_AddHydrogen_hh

#include <protocols/ligand_docking/AddHydrogen.fwd.hh>

#include <core/pose/Pose.fwd.hh>
//#include <core/conformation/Residue.hh>
#include <protocols/moves/Mover.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

static THREAD_LOCAL basic::Tracer add_hydrogen_tracer( "protocols.ligand_docking.AddHydrogen", basic::t_debug );

/// @brief
///
/// @details
///
class AddHydrogen : public protocols::moves::Mover{

public:
	AddHydrogen();
	AddHydrogen(core::Size const residue, core::Size const connection_id);
	virtual ~AddHydrogen();
	AddHydrogen(AddHydrogen const & that);
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	core::Size residue_index_;
	core::Size connection_id_;

}; // class AddHydrogen

std::string generate_unique_name(std::string input_name="");

} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_AddHydrogen_HH
