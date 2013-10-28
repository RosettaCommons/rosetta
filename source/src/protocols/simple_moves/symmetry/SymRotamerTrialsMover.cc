// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author  Ingemar Andre

// Unit headers
#include <protocols/simple_moves/symmetry/SymRotamerTrialsMoverCreator.hh>
#include <protocols/simple_moves/symmetry/SymRotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_trials.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
// AUTO-REMOVED #include <protocols/rosetta_scripts/util.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {
namespace symmetry {

// creator
std::string
SymRotamerTrialsMoverCreator::keyname() const {
	return SymRotamerTrialsMoverCreator::mover_name();
}

protocols::moves::MoverOP
SymRotamerTrialsMoverCreator::create_mover() const {
	return new SymRotamerTrialsMover;
}

std::string
SymRotamerTrialsMoverCreator::mover_name() {
	return "SymRotamerTrialsMover";
}

//////////////////////////
// default constructor
SymRotamerTrialsMover::SymRotamerTrialsMover() : protocols::simple_moves::RotamerTrialsMover()
{
	protocols::moves::Mover::type( "SymRotamerTrials" );
}

// constructor with arguments
SymRotamerTrialsMover::SymRotamerTrialsMover(
	ScoreFunctionCOP scorefxn_in,
	PackerTask & task_in
) : protocols::simple_moves::RotamerTrialsMover(scorefxn_in, task_in )
{
	protocols::moves::Mover::type( "SymRotamerTrials" );
}

// constructor with arguments
SymRotamerTrialsMover::SymRotamerTrialsMover(
	ScoreFunctionCOP scorefxn_in,
	TaskFactoryCOP factory_in
) : protocols::simple_moves::RotamerTrialsMover(scorefxn_in, factory_in )
{
	protocols::moves::Mover::type( "SymRotamerTrials" );
}

SymRotamerTrialsMover::~SymRotamerTrialsMover() {}

void
SymRotamerTrialsMover::apply( core::pose::Pose & pose )
{
	core::pack::task::PackerTaskOP symmetric_task( task(pose)->clone() );
	make_symmetric_task( pose, symmetric_task );
	core::pack::symmetric_rotamer_trials( pose, *scorefxn(), symmetric_task );
}

std::string
SymEnergyCutRotamerTrialsMover::get_name() const {
	return "SymEnergyCutRotamerTrialsMover";
}

std::string
SymRotamerTrialsMover::get_name() const {
	return "SymRotamerTrialsMover";
}

void
SymRotamerTrialsMover::make_symmetric_task(
  core::pose::Pose & pose,
  core::pack::task::PackerTaskOP task
)
{
  assert( core::pose::symmetry::is_symmetric( pose ) );
  SymmetricConformation & SymmConf (
    dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
  core::conformation::symmetry::SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

  utility::vector1<bool> allow_repacked( pose.total_residue(), false );
  for (Size res=1; res <= pose.total_residue(); ++res )
    {
        if ( symm_info->fa_is_independent(res) ) allow_repacked.at(res) = true;
    }
    task->restrict_to_residues( allow_repacked );
}

/// @brief parse xml
void
SymRotamerTrialsMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &fm,
	protocols::moves::Movers_map const &mm,
	Pose const &pose ) 
{
	RotamerTrialsMover::parse_my_tag( tag,data,fm,mm,pose );
}

/////////////////////////
// default constructor
SymEnergyCutRotamerTrialsMover::SymEnergyCutRotamerTrialsMover() :
  SymRotamerTrialsMover()
{
  protocols::moves::Mover::type( "SymEnergyCutRotamerTrials" );
}

// constructor with arguments
SymEnergyCutRotamerTrialsMover::SymEnergyCutRotamerTrialsMover(
  ScoreFunctionCOP scorefxn_in,
  PackerTask & task_in,
  protocols::moves::MonteCarloOP mc_in,
  core::Real energycut_in
) : SymRotamerTrialsMover(scorefxn_in, task_in), mc_( mc_in ), energycut_( energycut_in )
{
  protocols::moves::Mover::type( "SymEnergyCutRotamerTrials" );
}

// constructor with arguments
SymEnergyCutRotamerTrialsMover::SymEnergyCutRotamerTrialsMover(
  ScoreFunctionCOP scorefxn_in,
  TaskFactoryCOP factory_in,
  protocols::moves::MonteCarloOP mc_in,
  core::Real energycut_in
) : SymRotamerTrialsMover(scorefxn_in, factory_in), mc_( mc_in ), energycut_( energycut_in )
{
  protocols::moves::Mover::type( "SymEnergyCutRotamerTrials" );
}

SymEnergyCutRotamerTrialsMover::~SymEnergyCutRotamerTrialsMover() {}

void
SymEnergyCutRotamerTrialsMover::apply( core::pose::Pose & pose )
{
	  PackerTaskOP rottrial_task( task(pose)->clone() );
  ( *scorefxn() )(pose);
  /// Now handled automatically.  scorefxn()->accumulate_residue_total_energies( pose );
  setup_energycut_task( pose, *mc(), *rottrial_task );
  /// This call is dangerous.  If sequence or length has changed since task was created, it will crash.
  /// Not a problem if you used a TaskFactory
  core::pack::symmetric_rotamer_trials( pose, *scorefxn(), rottrial_task );
}

/// @details starting from a fresh task, it reduces the number of residues to be repacked to only
/// those whose energy has increased by energycut_ since the application of the last move.
void
SymEnergyCutRotamerTrialsMover::setup_energycut_task(
	core::pose::Pose const & pose,
	protocols::moves::MonteCarlo const & mc,
	core::pack::task::PackerTask & task_in
) const
{
	using namespace core;
	using core::scoring::total_score;

	//Size count_fixed( 0 ), count_repacked( 0 );

	task_in.restrict_to_repacking();

	for ( int i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
		core::Real const resE ( pose.energies().residue_total_energy(i) );
		core::Real const lowest_resE( mc.lowest_score_pose().energies().residue_total_energy(i) );
		core::Real const deltaE ( resE - lowest_resE );
		if ( deltaE < energycut_ ) {
			task_in.nonconst_residue_task(i).prevent_repacking();
			//++count_fixed;
		} else {
			// let this residue be repacked
			//++count_repacked;
		}
	}
}
// protected accessor function for derived mover
protocols::moves::MonteCarloOP
SymEnergyCutRotamerTrialsMover::mc()
{
	return mc_;
}

void
SymEnergyCutRotamerTrialsMover::make_symmetric_task(
  core::pose::Pose & pose,
  core::pack::task::PackerTaskOP task
)
{
  assert( core::pose::symmetry::is_symmetric( pose ) );
  if( task->symmetrize_by_union() || task->symmetrize_by_intersection() ) return; // new machinery
  SymmetricConformation & SymmConf (
    dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
  core::conformation::symmetry::SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

  utility::vector1<bool> allow_repacked( pose.total_residue(), false );
  for (Size res=1; res <= pose.total_residue(); ++res )
    {
        if ( symm_info->fa_is_independent(res) ) allow_repacked.at(res) = true;
    }
    task->restrict_to_residues( allow_repacked );

}


} // symmetry
} // moves
} // protocols
