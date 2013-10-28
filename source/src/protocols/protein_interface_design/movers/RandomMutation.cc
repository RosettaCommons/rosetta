// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/RandomMutation.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/RandomMutation.hh>
#include <protocols/protein_interface_design/movers/RandomMutationCreator.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
// Package headers
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperations.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <boost/foreach.hpp>

#include <utility/vector0.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>


#define foreach BOOST_FOREACH

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace std;
using namespace core::scoring;

static basic::Tracer TR( "protocols.protein_interface_design.movers.RandomMutation" );
static numeric::random::RandomGenerator RG( 2111918 );

std::string
RandomMutationCreator::keyname() const
{
	return RandomMutationCreator::mover_name();
}

protocols::moves::MoverOP
RandomMutationCreator::create_mover() const {
	return new RandomMutation;
}

std::string
RandomMutationCreator::mover_name()
{
	return "RandomMutation";
}

RandomMutation::RandomMutation() :
	Mover( RandomMutationCreator::mover_name() ),
	task_factory_( NULL ),
	scorefxn_( NULL )
{
}


RandomMutation::~RandomMutation() {}

void
RandomMutation::apply( core::pose::Pose & pose )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

	PackerTaskCOP task;
	if ( cache_task_ && task_ ) {
		if( pose.total_residue() == task_->total_residue() ) {
			task = task_;
		} else {
			task_ = 0; // Invalidate cached task.
		}
	}
	if ( ! task ) {
		task = task_factory()->create_task_and_apply_taskoperations( pose );
	}
	if ( cache_task_ && !task_ ) {
		task_ = task;
	}

	utility::vector1< core::Size > being_designed;
	being_designed.clear();

	for( core::Size resi = 1; resi <= pose.total_residue(); ++resi ){
		if( task->residue_task( resi ).being_designed() && pose.residue(resi).is_protein() )
			being_designed.push_back( resi );
	}
	if( being_designed.empty() ) {
		TR.Warning << "WARNING: No residues are listed as designable." << std::endl;
		return;
	}
	core::Size const random_entry = being_designed[ (core::Size) floor( RG.uniform() * being_designed.size() )+1 ];
  typedef list< ResidueTypeCOP > ResidueTypeCOPList;
  ResidueTypeCOPList const & allowed( task->residue_task( random_entry ).allowed_residue_types() );
  utility::vector1< AA > allow_temp;
  allow_temp.clear();
  foreach( ResidueTypeCOP const t, allowed ){
		if( t->aa() != pose.residue( random_entry ).aa() )
    	allow_temp.push_back( t->aa() );
	}

  AA const target_aa( allow_temp[ (core::Size) floor( RG.uniform() * allow_temp.size() ) + 1 ] );
  utility::vector1< bool > allowed_aas;
  allowed_aas.clear();
  allowed_aas.assign( num_canonical_aas, false );
  allowed_aas[ target_aa ] = true;
  //PackerTaskOP mutate_residue = task_factory()->create_task_and_apply_taskoperations( pose );
	PackerTaskOP mutate_residue( task->clone() );
	mutate_residue->initialize_from_command_line().or_include_current( true );
  for( core::Size resi = 1; resi <= pose.total_residue(); ++resi ){
    if( resi != random_entry )
      mutate_residue->nonconst_residue_task( resi ).restrict_to_repacking();
    else
      mutate_residue->nonconst_residue_task( resi ).restrict_absent_canonical_aas( allowed_aas );
  }
  TR<<"Mutating residue "<<pose.residue( random_entry ).name3()<<random_entry<<" to ";
	protocols::simple_moves::PackRotamersMoverOP pack;
  if( core::pose::symmetry::is_symmetric( pose ) )
  	pack =  new protocols::simple_moves::symmetry::SymPackRotamersMover( scorefxn(), mutate_residue );
  else
    pack = new protocols::simple_moves::PackRotamersMover( scorefxn(), mutate_residue );
  pack->apply( pose );
  TR<<pose.residue( random_entry ).name3()<<std::endl;
  (*scorefxn())(pose);
}

std::string
RandomMutation::get_name() const {
	return RandomMutationCreator::mover_name();
}

void
RandomMutation::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
  task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
  scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	cache_task_ = tag->getOption< bool >( "cache_task", false );
}

protocols::moves::MoverOP
RandomMutation::clone() const {
    return( protocols::moves::MoverOP( new RandomMutation( *this ) ));
}

core::scoring::ScoreFunctionOP
RandomMutation::scorefxn() const{
	return scorefxn_;
}

void
RandomMutation::scorefxn( core::scoring::ScoreFunctionOP scorefxn )
{
	scorefxn_ = scorefxn;
}

core::pack::task::TaskFactoryOP
RandomMutation::task_factory() const{
	return( task_factory_ );
}

void
RandomMutation::task_factory( core::pack::task::TaskFactoryOP task_factory){
	task_factory_ = task_factory;
}

bool RandomMutation::cache_task() const {
	return( cache_task_ );
}

void RandomMutation::cache_task( bool cache ) {
	cache_task_ = cache;
}

} //movers
} //protein_interface_design
} //protocols
