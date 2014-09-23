// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   ForceDisulfidesMover.cc
///
/// @brief
/// @author Sarel Fleishman

// unit headers
#include <protocols/simple_moves/ForceDisulfidesMover.hh>
#include <protocols/simple_moves/ForceDisulfidesMoverCreator.hh>
#include <boost/foreach.hpp>

// type headers
#include <core/types.hh>
#include <core/id/types.hh>

// project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>

// utility header
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/scoring/ScoreFunction.hh>

namespace protocols {
namespace simple_moves {

static thread_local basic::Tracer TR( "protocols.simple_moves.ForceDisulfidesMover" );

std::string
ForceDisulfidesMoverCreator::keyname() const
{
	return ForceDisulfidesMoverCreator::mover_name();
}

protocols::moves::MoverOP
ForceDisulfidesMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ForceDisulfidesMover );
}

std::string
ForceDisulfidesMoverCreator::mover_name()
{
	return "ForceDisulfides";
}

ForceDisulfidesMover::ForceDisulfidesMover() :
	protocols::moves::Mover("ForceDisulfidesMover"),
	scorefxn_( /* NULL */ )
{ disulfides_.clear(); }

ForceDisulfidesMover::~ForceDisulfidesMover() {}

protocols::moves::MoverOP
ForceDisulfidesMover::clone() const
{
	return protocols::moves::MoverOP( new ForceDisulfidesMover( *this ) );
}

protocols::moves::MoverOP
ForceDisulfidesMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new ForceDisulfidesMover() );
}

void
ForceDisulfidesMover::apply( Pose & pose ) {
	TR<<"Fixing disulfides"<<std::endl;
	pose.conformation().fix_disulfides( disulfides_ );

	using namespace core::pack::task;
	using namespace protocols::toolbox::task_operations;
	using core::pack::task::operation::TaskOperationCOP;
	
	DesignAroundOperationOP dao( new DesignAroundOperation );
	dao->design_shell( 0.0 );
	dao->repack_shell( 6.0 );
	for( utility::vector1< std::pair< core::Size, core::Size > > ::const_iterator pair = disulfides_.begin(); pair != disulfides_.end(); ++pair ){
		dao->include_residue( pair->first );
		dao->include_residue( pair->second );
	}
	TaskFactoryOP tf( new TaskFactory );
	tf->push_back( dao );
	tf->push_back( TaskOperationCOP( new operation::InitializeFromCommandline ) );
	PackerTaskOP ptask = tf->create_task_and_apply_taskoperations( pose );
	PackRotamersMover prm( scorefxn(), ptask );
	TR<<"repacking disulfide surroundings"<<std::endl;
	prm.apply( pose );
}


void
ForceDisulfidesMover::scorefxn( core::scoring::ScoreFunctionOP sf ){
  scorefxn_ = sf;
}

core::scoring::ScoreFunctionOP
ForceDisulfidesMover::scorefxn() const{
  return scorefxn_;
}

std::string
ForceDisulfidesMover::get_name() const {
	return "ForceDisulfidesMover";
}

void
ForceDisulfidesMover::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose )
{
  scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	utility::vector1< std::string > const residue_pairs( utility::string_split( tag->getOption< std::string >( "disulfides" ), ',' ) );
	TR<<"Setting fix disulfides on residues: ";
	BOOST_FOREACH( std::string const residue_pair, residue_pairs ){
		utility::vector1< std::string > const residues( utility::string_split( residue_pair, ':' ));
		runtime_assert( residues.size() == 2);
		core::Size const res1( core::pose::parse_resnum( residues[ 1 ], pose ) );
		core::Size const res2( core::pose::parse_resnum( residues[ 2 ], pose ) );
		disulfides_.push_back( std::pair< core::Size, core::Size >( res1, res2 ) );
		TR<<res1<<':'<<res2<<',';
	}
	TR<<std::endl;
}

void
ForceDisulfidesMover::disulfides( utility::vector1< std::pair < core::Size, core::Size >  > d ){ disulfides_ = d; }

utility::vector1< std::pair< core::Size, core::Size > >
ForceDisulfidesMover::disulfides() const { return disulfides_; }

} // moves
} // protocols
