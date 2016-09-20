// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/coupled_moves/CoupledMovesProtocol.hh
/// @brief Mover that implements the CoupledMovesProtocol
/// @author Noah <nollikai@gmail.com>, refactored into header by Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_coupled_moves_CoupledMovesProtocol_hh
#define INCLUDED_protocols_coupled_moves_CoupledMovesProtocol_hh

#include <protocols/coupled_moves/CoupledMovesProtocol.fwd.hh>

#include <protocols/moves/Mover.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace coupled_moves {

class CoupledMovesProtocol : public protocols::moves::Mover {
public:
	CoupledMovesProtocol();
	CoupledMovesProtocol(CoupledMovesProtocol const & cmp);

	virtual void apply( core::pose::Pose& pose );
	std::string get_name() const { return "CoupledMovesProtocol"; }

	virtual protocols::moves::MoverOP clone() const {
		return protocols::moves::MoverOP( new CoupledMovesProtocol( *this ) );
	}

	virtual protocols::moves::MoverOP fresh_instance() const {
		return protocols::moves::MoverOP( new CoupledMovesProtocol );
	}

	core::Real compute_ligand_score_bonus(
		core::pose::PoseOP pose,
		utility::vector1<core::Size> ligand_resnums,
		core::Real ligand_weight);

private:
	core::scoring::ScoreFunctionOP score_fxn_;
	core::pack::task::TaskFactoryOP main_task_factory_;

}; //CoupledMovesProtocol

} //coupled_moves
} //protocols

#endif //INCLUDED_protocols_coupled_moves_CoupledMovesProtocol_HH
