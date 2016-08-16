// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/FindConsensusSequence.cc
/// @brief Takes in multiple poses from the MSDJobDistributor and finds the consensus sequence that optimizes energy of all input poses
/// @author Alex Sevy (alex.sevy@gmail.com)

#include <protocols/simple_moves/FindConsensusSequence.hh>
#include <protocols/simple_moves/FindConsensusSequenceCreator.hh>

// type headers
#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/VectorPoseMover.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/string_util.hh>
#include <protocols/simple_moves/MutateResidue.hh>

namespace protocols {
namespace simple_moves {

static basic::Tracer TR("protocols.simple_moves.FindConsensusSequence");

std::string
FindConsensusSequenceCreator::keyname() const
{
	return FindConsensusSequenceCreator::mover_name();
}

moves::MoverOP
FindConsensusSequenceCreator::create_mover() const {
	return FindConsensusSequenceOP( new FindConsensusSequence );
}

std::string
FindConsensusSequenceCreator::mover_name()
{
	return "FindConsensusSequence";
}

FindConsensusSequence::FindConsensusSequence() :
	moves::VectorPoseMover( FindConsensusSequenceCreator::mover_name() ),
	sfxn_( core::scoring::get_score_function() )
{}

FindConsensusSequence::~FindConsensusSequence() {}

moves::MoverOP FindConsensusSequence::clone() const {
	return FindConsensusSequenceOP( new FindConsensusSequence( *this ) );
}

moves::MoverOP FindConsensusSequence::fresh_instance() const {
	return FindConsensusSequenceOP ( new FindConsensusSequence );
}

void FindConsensusSequence::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & datamap,
	filters::Filters_map const &,
	moves::Movers_map const & ,//movers,
	core::pose::Pose const &
)  {

	if ( tag->hasOption("scorefxn") ) {
		std::string const scorefxn_key( tag->getOption<std::string>("scorefxn") );
		if ( datamap.has( "scorefxns", scorefxn_key ) ) {
			sfxn_ = datamap.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_key );
		} else {
			throw utility::excn::EXCN_RosettaScriptsOption("ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap.");
		}
	}

	if ( tag->hasOption("task_operations") ) {
		task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, datamap );
	} else {
		task_factory_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory );
	}

	if ( tag->hasOption( "resfiles" ) ) {
		resfiles_ = utility::string_split( tag->getOption<std::string>( "resfiles" ), ',' );
	}

}

void
FindConsensusSequence::parse_resfiles ()
{
	res_links_.clear();
	if ( resfiles_.size() > 0 ) {
		utility::vector1< core::pack::task::PackerTaskOP > design_tasks ( poses_.size() );
		for ( core::Size ii = 1; ii <= poses_.size(); ++ii ) {
			design_tasks[ ii ] = core::pack::task::TaskFactory::create_packer_task( *poses_[ ii ] );
			core::pack::task::parse_resfile( *poses_[ ii ], *design_tasks[ ii ], resfile_at( ii ) );

		}
		core::Size designable_residues = -1;
		for ( core::Size i = 1; i <= design_tasks.size(); ++i ) {
			res_links_.push_back( parse_resfile( design_tasks[ i ]) );
			if ( i == 1 )  designable_residues = res_links_[ 1 ].size();
			else {
				if ( res_links_[ i ].size() != designable_residues ) {
					utility_exit_with_message( "All resfiles must have the same number of designable residues");
				}
			}
		}
	} else {
		for ( core::Size i = 1; i <= poses_.size(); ++i ) {
			utility::vector1< core::Size > positions;
			for ( core::Size j = 1; j <= poses_[ i ]->total_residue(); ++j ) {
				positions.push_back( j );
			}
			res_links_.push_back( positions );
			if ( positions.size() != res_links_[ 1 ].size() ) {
				utility_exit_with_message( "All poses must have the same number of residues, otherwise provide a resfile");
			}
		}
	}

}

utility::vector1< core::Size >
FindConsensusSequence::parse_resfile ( core::pack::task::PackerTaskCOP design_task )
{
	utility::vector1< core::Size > vector;
	utility::vector1<bool> designing = design_task->designing_residues();
	for ( core::Size i = 1; i <= designing.size(); ++i ) {
		if ( designing[ i ] ) {
			vector.push_back( i );
		}
	}
	return vector;
}

std::string
FindConsensusSequence::resfile_at ( core::Size index ) {
	if ( resfiles_.size() == 1 ) {
		return resfiles_[ 1 ];
	} else if ( resfiles_.size() == 0 ) {
		utility_exit_with_message("No resfiles initialized");
	} else {
		if ( index > resfiles_.size() ) {
			utility_exit_with_message("Number of resfiles does not match number of structures");
		} else {
			return resfiles_[ index ];
		}
	}
}

void FindConsensusSequence::apply( core::pose::Pose & /*pose*/ ) {
	if ( poses_.size() == 0 ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Error: no poses initialized. Did you try to run this from Rosetta Scripts? That's not supported right now");
	}
	using namespace core::pack::task;
	simple_moves::PackRotamersMoverOP packer ( new simple_moves::PackRotamersMover );
	packer->score_function( sfxn_ );
	operation::RestrictToRepackingOP rtr ( new operation::RestrictToRepacking );
	if ( !task_factory_ ) {
		task_factory_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory );
	}
	task_factory_->push_back( rtr );
	packer->task_factory( task_factory_ );
	parse_resfiles();
	simple_moves::MutateResidueOP mutation_mover;
	//bool diff;
	for ( core::Size i = 1; i <= res_links_[ 1 ].size(); ++i ) {
		core::Size seqpos1 = res_links_[ 1 ][ i ];
		//check to see if they are all the same
		bool diff = false;
		for ( core::Size j = 2; j <= poses_.size(); ++j ) {
			core::Size seqpos2 = res_links_[ j ][ i ];
			if ( poses_[ j ]->residue( seqpos2 ).aa() != poses_[ 1 ]->residue( seqpos1 ).aa() ) {
				diff = true;
				break;
			}
		}
		if ( diff ) {
#ifndef NDEBUG
			TR << "picking the best residue at position " << seqpos1 << std::endl;
#endif


			core::Size min_index = 0;
			core::Real min_score = 0;
			utility::vector1< core::Real > scores ( poses_.size(), 0 );
			for ( core::Size ref = 1; ref <= poses_.size(); ++ref ) {
				for ( core::Size comp = 1; comp <= poses_.size(); ++comp ) {
					//if residues at this position are same
					if ( poses_[ ref ]->residue( res_links_[ ref ][ i ] ).aa() ==
							poses_[ comp ]->residue( res_links_[ comp ][ i ] ).aa() ) {
						scores[ ref ] += poses_[ comp ]->energies().total_energy();
						continue;
					} else { //if they're different

						core::pose::PoseOP mutpose = poses_[ comp ]->clone();
						mutation_mover = simple_moves::MutateResidueOP ( new simple_moves::MutateResidue (
							res_links_[ comp ][ i ], //position
							poses_[ ref ]->residue( res_links_[ ref ][ i ] ).name3() //residue
							) );
						mutation_mover->apply( *mutpose );
						packer->apply( *mutpose );
						scores[ ref ] += mutpose->energies().total_energy();
					}

				}
				if ( scores[ ref ] < min_score || ref == 1 ) {
					min_index = ref;
					min_score = scores[ ref ];
				}
			}
#ifndef NDEBUG
			TR << "best residue is " << poses_[ min_index ]->residue( res_links_[ min_index ][ i ] ).name3() << std::endl;
#endif
			for ( core::Size pose = 1; pose <= poses_.size(); ++pose ) {
				if ( poses_[ pose ]->residue( res_links_[ pose ][ i ] ).aa() !=
						poses_[ min_index ]->residue( res_links_[ min_index ][ i ] ).aa() ) {
					mutation_mover = simple_moves::MutateResidueOP ( new simple_moves::MutateResidue (
						res_links_[ pose ][ i ], //position
						poses_[ min_index ]->residue( res_links_[ min_index ][ i ] ).name3() //residue
						) );
					mutation_mover->apply( *poses_[ pose ] );
					packer->apply( *poses_[ pose ] );
				}
			}
		}
	}
}

std::string FindConsensusSequence::get_name() const {
	return "FindConsensusSequence";
}

core::scoring::ScoreFunctionOP FindConsensusSequence::score_function() const {
	return sfxn_;
}

void FindConsensusSequence::score_function( core::scoring::ScoreFunctionOP sfxn) {
	sfxn_ = sfxn;
}

core::pack::task::TaskFactoryOP FindConsensusSequence::task_factory() const {
	return task_factory_;
}

void FindConsensusSequence::task_factory( core::pack::task::TaskFactoryOP task_factory ) {
	task_factory_ = task_factory;
}

void FindConsensusSequence::resfiles( utility::vector1< std::string > & resfiles ) {
	resfiles_ = resfiles;
}

utility::vector1< utility::vector1< core::Size > > FindConsensusSequence::res_links () { return res_links_; }

} //simple_moves
} //protocols
