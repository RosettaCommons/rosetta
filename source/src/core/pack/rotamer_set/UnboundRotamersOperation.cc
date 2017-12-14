// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/UnboundRotamersOperation.cc
///
/// @brief
/// @author Ian W. Davis


#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>

#include <core/conformation/Residue.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <basic/options/option.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

#ifdef WIN32
#include <utility/graph/Graph.hh>
#endif


// option key includes

#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/vector1.hh>

//Auto Headers
//#include <core/import_pose/import_pose.hh>


namespace core {
namespace pack {
namespace rotamer_set {


static basic::Tracer TR( "core.pack.rotamer_set.UnboundRotamersOperation" );


UnboundRotamersOperation::UnboundRotamersOperation():
	RotamerSetOperation(),
	total_rot_(0),
	poses_()
{
}


UnboundRotamersOperation::~UnboundRotamersOperation() = default;


core::pack::rotamer_set::RotamerSetOperationOP
UnboundRotamersOperation::clone() const
{
	return core::pack::rotamer_set::RotamerSetOperationOP( new UnboundRotamersOperation( *this ) );
}


void UnboundRotamersOperation::add_pose(core::pose::PoseCOP pose)
{
	poses_.push_back(pose);
	if ( pose->size() > total_rot_ ) total_rot_ = pose->size();
}


Size UnboundRotamersOperation::total_residue()
{
	return total_rot_;
}


void UnboundRotamersOperation::initialize_from_command_line()
{
	using namespace basic::options;
	if ( !option[ OptionKeys::packing::unboundrot ].active() ) return;
	for ( Size i = 1; i <= option[ OptionKeys::packing::unboundrot ]().size(); ++i ) {
		std::string filename = option[ OptionKeys::packing::unboundrot ]()[i].name();
		TR << "Adding 'unbound' rotamers from " << filename << std::endl;
		core::pose::PoseOP pose( new core::pose::Pose() );
		//core::import_pose::pose_from_file( *pose, filename , core::import_pose::PDB_file);
		core::io::pdb::build_pose_from_pdb_as_is( *pose, filename );
		this->add_pose( pose );
	}
}


/// @brief Helper function, combines existing's metadata with conformer's conformation.
conformation::ResidueOP
dup_residue(
	conformation::Residue const & existing,
	conformation::Residue const & conformer
)
{
	// The above is also bad:  existing may not be the same residue type as conformer!
	conformation::ResidueOP newrsd = conformer.clone();
	newrsd->chain( existing.chain() );
	newrsd->seqpos( existing.seqpos() );
	newrsd->copy_residue_connections_from( existing ); // this is probably not good enough if residue types diverge more than protonation state...

	return newrsd;
}


void
UnboundRotamersOperation::alter_rotamer_set(
	pose::Pose const & pose,
	scoring::ScoreFunction const & /*sfxn*/,
	task::PackerTask const & ptask,
	utility::graph::GraphCOP /*packer_neighbor_graph*/,
	core::pack::rotamer_set::RotamerSet & rotamer_set
)
{
	auto const seqnum = (Size) rotamer_set.resid();
	debug_assert( seqnum <= ptask.total_residue() );
	core::pack::task::ResidueLevelTask const & rtask = ptask.residue_task(seqnum);
	for ( Size i = 1; i <= poses_.size(); ++i ) {
		core::pose::Pose const & ubr_pose = *(poses_[i]);
		if ( seqnum > ubr_pose.size() ) continue;
		core::chemical::ResidueType const & restype = ubr_pose.residue_type(seqnum);
		bool type_is_allowed = false;
		for ( auto j = rtask.allowed_residue_types_begin(),
				j_end = rtask.allowed_residue_types_end(); j != j_end; ++j ) {
			if ( restype.name() == (**j).name() ) {
				type_is_allowed = true;
				break;
			}
		}
		if ( type_is_allowed ) {
			TR.Debug << "Adding 'unbound' rotamer at position " << seqnum << std::endl;
			conformation::ResidueOP newrsd = dup_residue( pose.residue(seqnum), ubr_pose.residue(seqnum) );
			newrsd->place( pose.residue(seqnum), pose.conformation() );
			rotamer_set.add_rotamer_into_existing_group( *newrsd );
		} else {
			TR.Debug << "Residue names do not match. Skipping 'unbound' rotamer at position " << seqnum << std::endl;
		}
	}
}


} // namespace rotamer_set
} // namespace pack
} // namespace core
