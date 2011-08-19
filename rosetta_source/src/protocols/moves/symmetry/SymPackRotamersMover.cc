// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author Ingemar Andre

// Unit headers
#include <protocols/moves/symmetry/SymPackRotamersMover.hh>

// Project headers
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

// AUTO-REMOVED #include <core/util/prof.hh>
#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.moves.symmetry.SymPackRotamersMover");

// Utility Headers
#include <utility/exit.hh>
#include <utility/vector0.hh>

namespace protocols {
namespace moves {
namespace symmetry {

using namespace core;
using namespace scoring;
using namespace pack;

/// PackRotamersMover

SymPackRotamersMover::SymPackRotamersMover()
	: PackRotamersMover(),
		sym_rotamer_sets_( new rotamer_set::symmetry::SymmetricRotamerSets() ),
		symmetric_ig_(0)
{}

	// constructors with arguments
SymPackRotamersMover::SymPackRotamersMover(
	ScoreFunctionCOP scorefxn,
	task::PackerTaskCOP task,
	Size nloop
) : PackRotamersMover( scorefxn, task, nloop ),
		sym_rotamer_sets_( new rotamer_set::symmetry::SymmetricRotamerSets() ),
		symmetric_ig_(0)
{}

SymPackRotamersMover::~SymPackRotamersMover(){}

SymPackRotamersMover::SymPackRotamersMover( PackRotamersMover const & other )
	: PackRotamersMover( other )
{
	sym_rotamer_sets_ = new rotamer_set::symmetry::SymmetricRotamerSets();

}

/*void
SymPackRotamersMover::apply( pose::Pose & pose )
{
	// jec update_residue_neighbors() required to update EnergyGraph (ensures graph_state == GOOD) when calling Interface.cc
	pose.update_residue_neighbors();
	// guarantee of valid ScoreFunction and PackerTask postponed until now

//	else assert( task_is_valid( pose ) );

	// get rotamers, energies
	this->setup( pose );

	this->run( pose );

}*/

std::string
SymPackRotamersMover::get_name() const {
	return "SymPackRotamersMover";
}

void SymPackRotamersMover::setup( pose::Pose & pose )
{

  // jec update_residue_neighbors() required to update EnergyGraph (ensures graph_state == GOOD) when calling Interface.cc
  pose.update_residue_neighbors();
  // guarantee of valid ScoreFunction and PackerTask postponed until now
  if ( score_function() == 0 ) {
    Warning() << "undefined ScoreFunction -- creating a default one" << std::endl;
    ScoreFunctionCOP scfx ( ScoreFunctionFactory::create_score_function( STANDARD_WTS ) );
		score_function( scfx );
  }

   // if present, task_factory_ always overrides/regenerates task_
		if ( task() != 0 ) {
			symmetric_task_ = (task())->clone();
		}
	  if ( task_factory() != 0 ) {
			symmetric_task_ = task_factory()->create_task_and_apply_taskoperations( pose );
		} else if ( task() == 0 ) {
    Warning() << "undefined PackerTask -- creating a default one" << std::endl;
    symmetric_task_ = core::pack::task::TaskFactory::create_packer_task( pose );
  }
  // in case PackerTask was not generated locally, verify compatibility with pose
  //else runtime_assert( task_is_valid( pose ) );
	make_symmetric_task( pose, symmetric_task_ );

	symmetric_pack_rotamers_setup( pose, *( score_function() ), symmetric_task_, sym_rotamer_sets_, symmetric_ig_ );

	setup_IG_res_res_weights( pose, symmetric_task_, sym_rotamer_sets_, symmetric_ig_ );
}

core::PackerEnergy SymPackRotamersMover::run( pose::Pose & pose, utility::vector0< int > rot_to_pack ) const
{
	return symmetric_pack_rotamers_run( pose, symmetric_task_, sym_rotamer_sets_, symmetric_ig_, rot_to_pack );
}

void
SymPackRotamersMover::make_symmetric_task(
	pose::Pose & pose,
	task::PackerTaskOP task
)
{
	assert( pose::symmetry::is_symmetric( pose ) );
	SymmetricConformation & SymmConf (
    dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
  core::conformation::symmetry::SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	utility::vector1<bool> allow_repacked( pose.total_residue(), false );
	for (Size res=1; res <= pose.total_residue(); ++res )
    {
        if ( pose.residue(res).aa() != core::chemical::aa_vrt && symm_info->fa_is_independent(res) )
					allow_repacked.at(res) = true;
    }
    task->restrict_to_residues( allow_repacked );

}

}
} // moves
} // protocols
