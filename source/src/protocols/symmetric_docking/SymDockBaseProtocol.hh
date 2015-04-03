// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/symmetric_docking/SymDockBaseProtocol.hh
///
/// @brief
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_symmetric_docking_SymDockBaseProtocol_hh
#define INCLUDED_protocols_symmetric_docking_SymDockBaseProtocol_hh

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/scoring/ScoreFunction.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace symmetric_docking {

class SymDockBaseProtocol : public protocols::moves::Mover
{
public:

	SymDockBaseProtocol();

	virtual ~SymDockBaseProtocol();

	core::scoring::ScoreFunctionOP get_lowres_scorefxn() { return scorefxn_lowres_; }
	core::scoring::ScoreFunctionOP get_highres_scorefxn() { return scorefxn_hires_; }

	virtual void apply( core::pose::Pose & /*pose*/ ) { utility_exit_with_message("Not intended to actually be used!"); }
	virtual std::string get_name() const;

protected:

	core::scoring::ScoreFunctionOP scorefxn_lowres_, scorefxn_hires_;
	std::string symm_definition_file_;

};

} // symmetric_docking
} // protocols
#endif // INCLUDED_protocols_symmetric_docking_SymDockBaseProtocol_HH
