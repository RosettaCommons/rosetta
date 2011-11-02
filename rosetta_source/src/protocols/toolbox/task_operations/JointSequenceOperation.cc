// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/task_operations/JointSequenceOperation.cc
/// @brief set designable residues to those observed in a set of structures
/// @author Rocco Moretti, rmoretti@u.washington.edu


// Unit Headers
#include <protocols/toolbox/task_operations/JointSequenceOperation.hh>
#include <protocols/toolbox/task_operations/JointSequenceOperationCreator.hh>

// Project Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/sequence/Sequence.hh>

#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <string>

#include <utility/vector0.hh>


static basic::Tracer TR("protocols.toolbox.tas_operations.JointSequenceOperation");

namespace protocols{
namespace toolbox{
namespace task_operations{

core::pack::task::operation::TaskOperationOP
JointSequenceOperationCreator::create_task_operation() const
{
	return new JointSequenceOperation;
}

/// @brief default constructor
JointSequenceOperation::JointSequenceOperation():
	TaskOperation(),
	use_current_pose_(true),
	use_natro_(false),
	ubr_(0)
{
}

/// @brief destructor
JointSequenceOperation::~JointSequenceOperation() {}

/// @brief clone
core::pack::task::operation::TaskOperationOP
JointSequenceOperation::clone() const {
	return new JointSequenceOperation( *this );
}

/// @brief all AA that have a higher probability in the seqprofile
/// than the native residue are allowed. probability also
/// needs to be higher than min_aa_probability_
void
JointSequenceOperation::apply( Pose const & pose, PackerTask & task ) const
{
	for( std::vector<core::sequence::SequenceOP>::const_iterator iter(sequences_.begin()); iter != sequences_.end(); iter++ ) {
		if( (*iter)->length() != pose.total_residue() ) {
				std::string name("current pdb");
				if(! pose.pdb_info() ) {
					name = pose.pdb_info()->name();
				}
				TR.Warning << "WARNING: Pose " << (*iter)->id() << " contains a different number of residues than " << name << std::endl;
		}
	}
	for( core::Size ii = 1; ii <= pose.total_residue(); ++ii){
		if( !pose.residue_type( ii ).is_protein() ) continue;

		utility::vector1< bool > allowed(core::chemical::num_canonical_aas, false);

		if(use_current_pose_) {
			allowed[ pose.aa(ii) ] = true;
		}
		for( std::vector<core::sequence::SequenceOP>::const_iterator iter(sequences_.begin()); iter != sequences_.end(); iter++ ) {
			if ( ii > (*iter)->length() ) continue; // ignore short references
			char aa( (*(*iter))[ ii ] );
			if( core::chemical::oneletter_code_specifies_aa(aa) ) {
				allowed[ core::chemical::aa_from_oneletter_code(aa)  ] = true;
			}
		}

		task.nonconst_residue_task(ii).restrict_absent_canonical_aas( allowed );
	} //loop over all residues

	if( use_natro_ && ubr_) {
		task.append_rotamerset_operation( ubr_ );
	}
} // apply

void
JointSequenceOperation::parse_tag( TagPtr tag )
{
	use_current_pose( tag->getOption< bool >( "use_current", true ) );
	use_natro( tag->getOption< bool >( "use_natro", false ) );
	if( tag->getOption< bool >( "use_native", false )) {
		if( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			add_native_pdb( basic::options::option[ basic::options::OptionKeys::in::file::native ] );
		}
		else {
			TR.Warning << "WARNING: Native PDB not specified on command line." << std::endl;
		}
	}
	if( tag->hasOption("filename") ){
		add_pdb( tag->getOption< String >( "filename" ) );
	}
	if( tag->hasOption("native") ){
		add_native_pdb( tag->getOption< String >( "native" ) );
	}
}

/// @brief Add the sequence from the given filename to the set of allowed aas.
void
JointSequenceOperation::add_pdb( std::string filename )
{
	core::pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, filename );
	add_pose( pose );
}

/// @brief Add the sequence from the given pose to the set of allowed aas.
void
JointSequenceOperation::add_pose( Pose const & pose )
{
	std::string name("unknown");
	if(! pose.pdb_info() ) {
		name = pose.pdb_info()->name();
	}
	sequences_.push_back( new core::sequence::Sequence(pose.sequence(), name) );
}

/// @brief Add the sequence from the given filename to the set of allowed aas
/// and add the rotamers to the set of possible rotamers (if use_natro_ is set)
void
JointSequenceOperation::add_native_pdb( std::string filename ) {
	core::pose::PoseOP poseop(new core::pose::Pose);
	core::import_pose::pose_from_pdb( *poseop, filename );
	add_native_pose( poseop );
}

/// @brief Add the sequence from the given pose to the set of allowed aas
/// and add the rotamers to the set of possible rotamers
void
JointSequenceOperation::add_native_pose( core::pose::PoseCOP posecop ){ // PoseCOP needed by UnboundRot, unfortunately
	if( use_natro_ ) { // Deliberate check now to avoid keeping native poses around if we're never going to need them.
		ubr_->add_pose(posecop);
	}
	add_pose(*posecop);
}

/// @brief Should the current pose (pose supplied to apply) be used in addition to the other ones?
void
JointSequenceOperation::use_current_pose( bool ucp )
{
	use_current_pose_ = ucp;
}

/// @brief Should the rotamers for the native poses be used?
void
JointSequenceOperation::use_natro( bool unr ) {
	if( !use_natro_ && unr ) {
		ubr_ = new core::pack::rotamer_set::UnboundRotamersOperation;
	}
	if( use_natro_ && !unr ) {
		ubr_ = 0; // Allow owning pointer to garbage collect as necessary.
	}
	use_natro_ = unr;
}

} // TaskOperations
} // toolbox
} // protocols

