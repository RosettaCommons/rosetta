// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/task_operations/RestrictNativeResiduesOperation.cc
/// @brief Restrict every residue in the current pose that is native to repacking. ie, only allow mutated positions to be designed.
/// @author Jacob Bale, balej@u.washington.edu


// Unit Headers
#include <protocols/toolbox/task_operations/RestrictNativeResiduesOperation.hh>
#include <protocols/toolbox/task_operations/RestrictNativeResiduesOperationCreator.hh>

// Project Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Conformation.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <string>
#include <ObjexxFCL/format.hh>
#include <utility/vector0.hh>


static thread_local basic::Tracer TR( "protocols.toolbox.task_operations.RestrictNativeResiduesOperation" );

namespace protocols{
namespace toolbox{
namespace task_operations{

core::pack::task::operation::TaskOperationOP
RestrictNativeResiduesOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictNativeResiduesOperation );
}

/// @brief default constructor
RestrictNativeResiduesOperation::RestrictNativeResiduesOperation():
	TaskOperation(),
	reference_pose_( /* NULL */ ),
	verbose_( false ),
	prevent_repacking_( 0 )
{
}

/// @brief destructor
RestrictNativeResiduesOperation::~RestrictNativeResiduesOperation() {}

/// @brief clone
core::pack::task::operation::TaskOperationOP
RestrictNativeResiduesOperation::clone() const {
	return core::pack::task::operation::TaskOperationOP( new RestrictNativeResiduesOperation( *this ) );
}

core::pose::PoseCOP
RestrictNativeResiduesOperation::reference_pose() const
{
	return reference_pose_;
}

void
RestrictNativeResiduesOperation::reference_pose( core::pose::PoseCOP pose )
{
	reference_pose_ = pose;
}

void
RestrictNativeResiduesOperation::reference_pose( core::pose::Pose const & pose )
{
	reference_pose_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( pose ) ) );
}

bool
RestrictNativeResiduesOperation::verbose() const
{
	return( verbose_ );
}

void
RestrictNativeResiduesOperation::verbose( bool const verb )
{
	verbose_ = verb;
}

bool
RestrictNativeResiduesOperation::prevent_repacking() const
{
	return( prevent_repacking_ );
}

void
RestrictNativeResiduesOperation::prevent_repacking( bool const prev )
{
	prevent_repacking_ = prev;
}

/// @brief Loop over the residues in the current pose and restrict those that match the
/// reference pose (ie, native residues) to repacking.
void
RestrictNativeResiduesOperation::apply( Pose const & pose, PackerTask & task ) const
{
	runtime_assert( reference_pose() != 0 );
	core::Size total_residue_ref;
	core::pose::Pose asym_ref_pose;
	if(core::pose::symmetry::is_symmetric( *reference_pose() )) {
		core::pose::symmetry::extract_asymmetric_unit( *reference_pose(), asym_ref_pose);
  	for (core::Size i = 1; i <= asym_ref_pose.total_residue(); ++i) {
    	if (asym_ref_pose.residue_type(i).name() == "VRT") {
				asym_ref_pose.conformation().delete_residue_slow(asym_ref_pose.total_residue());
			}
		}
		total_residue_ref = asym_ref_pose.total_residue();
	} else {
		total_residue_ref = reference_pose()->total_residue();
		asym_ref_pose = *reference_pose();
	}
	core::Size total_residue;
	core::pose::Pose asym_pose;
	if (core::pose::symmetry::is_symmetric( pose )) {
		core::pose::symmetry::extract_asymmetric_unit(pose, asym_pose);
  	for (core::Size i = 1; i <= asym_pose.total_residue(); ++i) {
    	if (asym_pose.residue_type(i).name() == "VRT") {
				asym_pose.conformation().delete_residue_slow(asym_pose.total_residue());
			}
		}
		total_residue = asym_pose.total_residue();
	} else {
		total_residue = pose.total_residue();
		asym_pose = pose;
	}
	if( total_residue_ref != total_residue )
		utility_exit_with_message( "Reference pose and current pose have a different number of residues" );
	std::string select_non_native_resis("select non_native_resis, resi ");
	core::Size designable_resis = 0;
	for( core::Size resi=1; resi<=total_residue; ++resi ) {
		if ( asym_ref_pose.residue(resi).name3() == asym_pose.residue(resi).name3() ) {
			if ( prevent_repacking_ ) task.nonconst_residue_task(resi).prevent_repacking();
			else task.nonconst_residue_task(resi).restrict_to_repacking();
		} else {
			if ( verbose_ ) {
				core::Size output_resi = resi;
				if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
					output_resi = pose.pdb_info()->number( resi );
				}
				select_non_native_resis.append(ObjexxFCL::string_of(output_resi) + "+");
			}
			designable_resis++;
		}
	}
	// turn off everything in symmetric copies
	for( core::Size ir = total_residue+1; ir <= pose.total_residue(); ++ir){
		 task.nonconst_residue_task(ir).prevent_repacking();
	}

	if( designable_resis == 0 ) {
		TR<<"Warning: No designable residues identified in pose."<<std::endl;
	} else {
		TR<<designable_resis<<" non-native, designable residues found in pose"<<std::endl;
	}
	if ( verbose_ ) {
		TR<<select_non_native_resis<<std::endl;
	}
} // apply

void
RestrictNativeResiduesOperation::parse_tag( TagCOP tag , DataMap & )
{
	verbose( tag->getOption< bool >( "verbose", 0 ) );
	prevent_repacking( tag->getOption< bool >( "prevent_repacking", 0 ) );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string reference_pdb = "";

	if( option[ in::file::native ].user() ){
		reference_pdb = option[ in::file::native ]();
	} else if ( tag->hasOption("pdbname") ) {
		reference_pdb = tag->getOption< std::string >("pdbname");
	}
	if ( reference_pdb != "" ) {
		core::pose::PoseOP temp_pose( new core::pose::Pose );
		core::import_pose::pose_from_pdb( *temp_pose, reference_pdb );
		reference_pose( temp_pose );
		TR<<"Using pdb "<<reference_pdb<<" as reference."<<std::endl;
	}	else {
		throw utility::excn::EXCN_RosettaScriptsOption( "Native PDB not specified." );
	}
}

} // TaskOperations
} // toolbox
} // protocols

