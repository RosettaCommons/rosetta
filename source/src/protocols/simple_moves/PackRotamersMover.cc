// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov

// Unit headers
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/PackRotamersMoverCreator.hh>

// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/elscripts/util.hh>

#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
// AUTO-REMOVED #include <utility/string_util.hh> // string_split

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

using namespace core;
	using namespace basic::options;
	using namespace pack;
		using namespace task;
			using namespace operation;
	using namespace scoring;

using basic::Warning;
using basic::t_warning;
static thread_local basic::Tracer TR( "protocols.simple_moves.PackRotamersMover" );

// PackRotamersMover

std::string
PackRotamersMoverCreator::keyname() const
{
	return PackRotamersMoverCreator::mover_name();
}

protocols::moves::MoverOP
PackRotamersMoverCreator::create_mover() const {
	return new PackRotamersMover;
}

std::string
PackRotamersMoverCreator::mover_name()
{
	return "PackRotamersMover";
}

PackRotamersMover::PackRotamersMover() :
	protocols::moves::Mover("PackRotamersMover"),
	scorefxn_(0),
	task_(0),
	nloop_( option[ OptionKeys::packing::ndruns ].value() ),
	task_factory_(0),
	rotamer_sets_( new rotamer_set::RotamerSets ),
	ig_(0)
{}

PackRotamersMover::PackRotamersMover( std::string const & type_name ) :
	protocols::moves::Mover( type_name ),
	scorefxn_(0),
	task_(0),
	nloop_( option[ OptionKeys::packing::ndruns ].value() ),
	task_factory_(0),
	rotamer_sets_( new rotamer_set::RotamerSets ),
	ig_(0)
{}

// constructors with arguments
PackRotamersMover::PackRotamersMover(
	ScoreFunctionCOP scorefxn,
	PackerTaskCOP task,
	Size nloop
) :
	protocols::moves::Mover("PackRotamersMover"),
	scorefxn_( scorefxn ),
	task_( task ),
	nloop_( nloop ),
	task_factory_(0),
	rotamer_sets_( new rotamer_set::RotamerSets ),
	ig_(0)
{}

PackRotamersMover::~PackRotamersMover(){}

PackRotamersMover::PackRotamersMover( PackRotamersMover const & other ) :
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( other )
{
	scorefxn_ = other.score_function();
	task_ = other.task();
	nloop_ = other.nloop();
	task_factory_ = other.task_factory();
	rotamer_sets_ = new rotamer_set::RotamerSets;
	ig_ = 0;
}

void
PackRotamersMover::apply( Pose & pose )
{
	// get rotamers, energies
	// also performs lazy initialization of ScoreFunction, PackerTask
	this->setup( pose );

	core::PackerEnergy best_energy(0.);
	Pose best_pose;
	best_pose = pose;
	for ( Size run(1); run <= nloop_; ++run ) {
		// run SimAnnealer
		core::PackerEnergy packer_energy( this->run( pose ) );
		// Real const score( scorefxn_( pose ) ); another option for deciding which is the 'best' result
		if ( run == 1 || packer_energy < best_energy ) {
			best_pose = pose;
			best_energy = packer_energy;
		}
	}
	if ( nloop_ > 1 ) pose = best_pose;

	//guaruntees proper scoring if this mover is used as a protocol (as in fixbb)
	(*scorefxn_)(pose);
	/// Now handled automatically.  scorefxn_->accumulate_residue_total_energies(pose);
}

std::string
PackRotamersMover::get_name() const {
	return PackRotamersMoverCreator::mover_name();
}

void
PackRotamersMover::show(std::ostream & output) const
{
	Mover::show(output);
	if ( score_function() != 0 ) {
		output << "Score function: " << score_function()->get_name() << std::endl;
	}
	else { output << "Score function: none" << std::endl; }
}

///@brief when the PackerTask was not generated locally, verify compatibility with pose
///@details the pose residue types must be equivalent to the ones used to generate the ResidueLevelTasks, because of the way that prevent_repacking and its associated flags work
bool
PackRotamersMover::task_is_valid( Pose const & pose ) const
{
	if ( task_->total_residue() != pose.total_residue() ) return false;
	for ( Size i(1); i <= pose.total_residue(); ++i ) {
		if ( ! task_->residue_task(i).is_original_type( &pose.residue_type(i) ) ) return false;
	}
	return true;
}

void PackRotamersMover::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & score_fxns,
				utility::lua::LuaObject const & tasks,
				protocols::moves::MoverCacheSP /*cache*/ ) {
	if( def["nloop"] ) {
		nloop_ = def["nloop"].to<Size>();
		runtime_assert( nloop_ > 0 );
	}

	if( def["scorefxn"] ) {
		score_function( protocols::elscripts::parse_scoredef( def["scorefxn"], score_fxns ) );
	} else {
		score_function( score_fxns["score12"].to<ScoreFunctionSP>()->clone()  );
	}
	if( def["tasks"] ) {
		TaskFactoryOP new_task_factory( protocols::elscripts::parse_taskdef( def["tasks"], tasks ));
		if ( new_task_factory == 0) return;
		task_factory( new_task_factory );
	}
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
PackRotamersMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & pose
)
{
	//fpd (commenting out below) classes derived from PackRotamers may call this function
	//if ( tag->getName() != "PackRotamersMover" ) {
	//	TR(t_warning) << " received incompatible Tag " << tag << std::endl;
	//	assert(false);
	//	return;
	//}
	if ( tag->hasOption("nloop") ) {
		nloop_ = tag->getOption<Size>("nloop",1);
		runtime_assert( nloop_ > 0 );
	}
	parse_score_function( tag, datamap, filters, movers, pose );
	parse_task_operations( tag, datamap, filters, movers, pose );
}

///@brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
void
PackRotamersMover::parse_score_function(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	ScoreFunctionOP new_score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
	if ( new_score_function == 0 ) return;
	score_function( new_score_function );
}

///@brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
void
PackRotamersMover::parse_task_operations(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == 0) return;
	task_factory( new_task_factory );
}

///@brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PackRotamersMover::fresh_instance() const
{
	return new PackRotamersMover;
}

///@brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PackRotamersMover::clone() const
{
	return new protocols::simple_moves::PackRotamersMover( *this );
}

///@brief get rotamers, energies. Also performs lazy initialization of ScoreFunction, PackerTask.
void PackRotamersMover::setup( Pose & pose )
{
	// jec update_residue_neighbors() required to update EnergyGraph (ensures graph_state == GOOD) when calling Interface.cc
	pose.update_residue_neighbors();
	// guarantee of valid ScoreFunction and PackerTask postponed until now
	if ( scorefxn_ == 0 ) {
		Warning() << "undefined ScoreFunction -- creating a default one" << std::endl;
		scorefxn_ = get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	}

	// if present, task_factory_ always overrides/regenerates task_
	if ( task_factory_ != 0 ) {
		task_ = task_factory_->create_task_and_apply_taskoperations( pose );
	} else if ( task_ == 0 ) {
		Warning() << "undefined PackerTask -- creating a default one" << std::endl;
		task_ = TaskFactory::create_packer_task( pose );
	}
	// in case PackerTask was not generated locally, verify compatibility with pose
	else runtime_assert( task_is_valid( pose ) );

	note_packertask_settings( pose );

	pack_rotamers_setup( pose, *scorefxn_, task_, rotamer_sets_, ig_ );

	setup_IG_res_res_weights( pose, task_, rotamer_sets_, ig_ );
}

core::PackerEnergy PackRotamersMover::run( Pose & pose, utility::vector0< int > rot_to_pack ) const
{
	return pack_rotamers_run( pose, task_, rotamer_sets_, ig_, rot_to_pack );
}

///@brief note PackerTask's packable and designable residues as string info
void PackRotamersMover::note_packertask_settings( Pose const & pose )
{
	std::ostringstream packable, designable;
	packable << "REMARK PackingRes";
	designable << "REMARK DesignRes";

	for ( Size i(1), end( task_->total_residue() ); i <= end; ++i ) {
		ResidueLevelTask const & rtask( task_->residue_task(i) );
		if ( rtask.being_designed() ) {
			designable << ", ";
			if ( pose.pdb_info() ) {
				designable << pose.pdb_info()->number(i) << " " << pose.pdb_info()->chain(i);
			} else {
				designable << i;
			}
		} else if ( rtask.being_packed() ) {
			packable << ", ";
			if ( pose.pdb_info() ) {
				packable << pose.pdb_info()->number(i) << " " << pose.pdb_info()->chain(i);
			} else {
				packable << i;
			}
		}
	}
	info().clear();
	info().push_back( packable.str() );
	info().push_back( designable.str() );
}

// setters
void PackRotamersMover::score_function( ScoreFunctionCOP sf )
{
	runtime_assert( sf );
	scorefxn_ = sf;
}

void PackRotamersMover::task( task::PackerTaskCOP t ) { task_ = t; }

void PackRotamersMover::task_factory( TaskFactoryCOP tf )
{
	runtime_assert( tf );
	task_factory_ = tf;
}

void PackRotamersMover::nloop( Size nloop_in ) { nloop_ = nloop_in; }

// accessors
ScoreFunctionCOP PackRotamersMover::score_function() const { return scorefxn_; }
PackerTaskCOP PackRotamersMover::task() const { return task_; }
TaskFactoryCOP PackRotamersMover::task_factory() const { return task_factory_; }
rotamer_set::RotamerSetsCOP PackRotamersMover::rotamer_sets() const { return rotamer_sets_; }
interaction_graph::InteractionGraphBaseCOP PackRotamersMover::ig() const { return ig_; }

std::ostream &operator<< (std::ostream &os, PackRotamersMover const &mover)
{
	mover.show(os);
	return os;
}

} // moves
} // protocols

