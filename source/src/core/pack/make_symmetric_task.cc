// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pose/symmetry/util.hh
/// @brief utility functions for handling of symmetric conformations
/// @author Ingemar Andre

// Unit headers
#include <core/pack/make_symmetric_task.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/id/AtomID.hh>

// Package Headers
#include <core/pose/symmetry/util.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

#include <utility/vector1.hh>

#include <basic/Tracer.hh>


static basic::Tracer TR( "core.pack.make_symmetric_task" );

namespace core {
namespace pack {

void
make_symmetric_PackerTask_by_truncation(
	pose::Pose const & pose,
	pack::task::PackerTaskOP task
)
{
	using namespace conformation::symmetry;
	using namespace pose::symmetry;

	debug_assert( core::pose::symmetry::is_symmetric( pose ) );

	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( !symm_info->chi_is_independent(i) ) {
			task->nonconst_residue_task( i ).prevent_repacking();
		}
	}
}

task::PackerTaskOP
make_new_symmetric_PackerTask_by_truncation(
	pose::Pose const & pose,
	task::PackerTaskCOP non_symmetric_task
){
	using namespace core::pack::task;
	debug_assert( core::pose::symmetry::is_symmetric( pose ) );
	PackerTaskOP new_task = non_symmetric_task->clone();
	make_symmetric_PackerTask_by_truncation(pose,new_task);
	return new_task;
}

task::PackerTaskOP
make_new_symmetric_PackerTask_by_union(
	pose::Pose const & pose,
	task::PackerTaskCOP non_symmetric_task
){
	using namespace core::pack::task;
	debug_assert( core::pose::symmetry::is_symmetric( pose ) );

	PackerTaskOP new_task = TaskFactory::create_packer_task(pose);
	PackerTask_ const & o(dynamic_cast<PackerTask_ const &>(*non_symmetric_task));
	PackerTask_       & n(dynamic_cast<PackerTask_       &>(*new_task));

	if ( !o.symmetrize_by_union() ) utility_exit_with_message("incorrect PackerTask symmetrization request");

	n.update_commutative(o);

	conformation::symmetry::SymmetricConformation const & SymmConf( dynamic_cast<conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
	conformation::symmetry::SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	for ( Size i = 1; i <= symm_info->num_total_residues_without_pseudo(); ++i ) {
		Size const ifollow = symm_info->chi_follows(i);
		if ( ifollow != 0 && ifollow != i ) {
			n.update_residue_union(ifollow,o.residue_task(i));
		}
	}

	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( !symm_info->chi_is_independent(i) ) {
			n.nonconst_residue_task( i ).prevent_repacking();
		}
	}

	return new_task;
}

task::PackerTaskOP
make_new_symmetric_PackerTask_by_intersection(
	pose::Pose const & pose,
	task::PackerTaskCOP non_symmetric_task
){
	using namespace core::pack::task;
	debug_assert( core::pose::symmetry::is_symmetric( pose ) );

	PackerTaskOP new_task = TaskFactory::create_packer_task(pose);
	PackerTask_ const & o(dynamic_cast<PackerTask_ const &>(*non_symmetric_task));
	PackerTask_       & n(dynamic_cast<PackerTask_       &>(*new_task));

	if ( !o.symmetrize_by_intersection() ) utility_exit_with_message("incorrect PackerTask symmetrization request");

	n.update_commutative(o);

	conformation::symmetry::SymmetricConformation const & SymmConf( dynamic_cast<conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
	conformation::symmetry::SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	for ( Size i = 1; i <= symm_info->num_total_residues_without_pseudo(); ++i ) {
		Size const ifollow = symm_info->chi_follows(i);
		if ( ifollow != 0 && ifollow != i ) {
			n.update_residue_intersection(ifollow,o.residue_task(i));
		}
	}

	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( !symm_info->chi_is_independent(i) ) {
			n.nonconst_residue_task( i ).prevent_repacking();
		}
	}

	return new_task;
}

task::PackerTaskOP
make_new_symmetric_PackerTask_by_requested_method(
	pose::Pose const & pose,
	task::PackerTaskCOP non_symmetric_task
){
	using namespace core::pack::task;
	debug_assert( core::pose::symmetry::is_symmetric( pose ) );
	PackerTask_ const & o(dynamic_cast<PackerTask_ const &>(*non_symmetric_task));

	if ( o.symmetrize_by_union() ) {
		return make_new_symmetric_PackerTask_by_union(pose,non_symmetric_task);
	}
	if ( o.symmetrize_by_intersection() ) {
		return make_new_symmetric_PackerTask_by_intersection(pose,non_symmetric_task);
	}
	// TR << "YOU HAVE NOT SPECIFIED HOW YOUR PACKERTASK SHOULD BE SYMMETRIZED, TRUNCATING IT!" << std::endl;
	// return make_new_symmetric_PackerTask_by_truncation(pose,non_symmetric_task);
	return non_symmetric_task->clone();
}

} // pack
} // core
