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
/// @author Chris King

// Unit headers
#include <protocols/simple_moves/LabelPDBInfoMover.hh>
#include <protocols/simple_moves/LabelPDBInfoMoverCreator.hh>

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
static thread_local basic::Tracer TR( "protocols.simple_moves.LabelPDBInfoMover" );

// LabelPDBInfoMover

std::string
LabelPDBInfoMoverCreator::keyname() const
{
	return LabelPDBInfoMoverCreator::mover_name();
}

protocols::moves::MoverOP
LabelPDBInfoMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new protocols::simple_moves::LabelPDBInfoMover::LabelPDBInfoMover );
}

std::string
LabelPDBInfoMoverCreator::mover_name()
{
	return "LabelPDBInfo";
}

LabelPDBInfoMover::LabelPDBInfoMover() :
	protocols::moves::Mover("LabelPDBInfo")
{}

LabelPDBInfoMover::LabelPDBInfoMover( std::string const & type_name ) :
	protocols::moves::Mover( type_name )
{}

// constructors with arguments
LabelPDBInfoMover::LabelPDBInfoMover(
	PackerTaskCOP task,
	std::string label
) :
	protocols::moves::Mover("LabelPDBInfo")
{
	task_ = task;
	label_ = label;
}

LabelPDBInfoMover::~LabelPDBInfoMover(){}

LabelPDBInfoMover::LabelPDBInfoMover( LabelPDBInfoMover const & other ) :
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( other )
{
	task_ = other.task();
	task_factory_ = other.task_factory();
	label_= other.label();
}

void
LabelPDBInfoMover::apply( Pose & pose )
{
//TODO write apply task ops and iter thru designable res and add label
//TODO option for setting pose as a diff pdb
//TODO diff labels for design/repack/prevent
	// if present, task_factory_ always overrides/regenerates task_
	if ( task_factory_ != 0 ) {
		task_ = task_factory_->create_task_and_apply_taskoperations( pose );
	} else if ( task_ == 0 ) {
		Warning() << "undefined PackerTask -- creating a default one" << std::endl;
		task_ = TaskFactory::create_packer_task( pose );
	}
	// in case PackerTask was not generated locally, verify compatibility with pose
	else runtime_assert( task_is_valid( pose ) );

	for ( Size i(1), end( task_->total_residue() ); i <= end; ++i ) {
		//seqpos label will be PDB res index if we have PDB info in the pose
		Size seqpos( i );
		if ( pose.pdb_info() ) seqpos = pose.pdb_info()->number( i );

		if ( task_->residue_task( i ).being_designed() && !label().empty() ) {
			pose.pdb_info()->add_reslabel( seqpos, label() );
		} else if ( task_->residue_task( i ).being_designed() ) {
		}
	}


}

std::string
LabelPDBInfoMover::get_name() const {
	return LabelPDBInfoMoverCreator::mover_name();
}

void
LabelPDBInfoMover::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Adding label " << std::endl; //output label name
}

///@brief when the PackerTask was not generated locally, verify compatibility with pose
///@details the pose residue types must be equivalent to the ones used to generate the
///@details ResidueLevelTasks, because of the way that prevent_repacking and its associated flags work
bool
LabelPDBInfoMover::task_is_valid( Pose const & pose ) const
{
	if ( task_->total_residue() != pose.total_residue() ) return false;
	for ( Size i(1); i <= pose.total_residue(); ++i ) {
		chemical::ResidueTypeCOP r = pose.residue_type(i).get_self_ptr();
		if ( ! task_->residue_task(i).is_original_type( r ) ) return false;
	}
	return true;
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
LabelPDBInfoMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & pose
)
{
	label( tag->getOption< std::string >( "label", "" ) );
	TR << "Adding PDB REMARK \"" << label() << "\" to designable residues" << std::endl;
	parse_task_operations( tag, datamap, filters, movers, pose );
}

///@brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
void
LabelPDBInfoMover::parse_task_operations(
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
LabelPDBInfoMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new protocols::simple_moves::LabelPDBInfoMover::LabelPDBInfoMover );
}

///@brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
LabelPDBInfoMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::simple_moves::LabelPDBInfoMover( *this ) );
}

// setters
void LabelPDBInfoMover::task( task::PackerTaskCOP t ) { task_ = t; }
void LabelPDBInfoMover::label( std::string l ) { label_ = l; }

void LabelPDBInfoMover::task_factory( TaskFactoryCOP tf )
{
	runtime_assert( tf != 0 );
	task_factory_ = tf;
}

// accessors
PackerTaskCOP LabelPDBInfoMover::task() const { return task_; }
TaskFactoryCOP LabelPDBInfoMover::task_factory() const { return task_factory_; }
std::string LabelPDBInfoMover::label() const { return label_; }

std::ostream &operator<< (std::ostream &os, LabelPDBInfoMover const &mover)
{
	mover.show(os);
	return os;
}

} // moves
} // protocols

