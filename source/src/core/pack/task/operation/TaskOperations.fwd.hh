// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/TaskOperation.fwd.hh
/// @brief  Forward declaration of a class that performs an operation on a packer task,
///         usually, by a PackerTaskFactory right after the task's construction
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_task_operation_TaskOperations_fwd_hh
#define INCLUDED_core_pack_task_operation_TaskOperations_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

class RestrictToRepacking;
class RestrictResidueToRepacking;
class InitializeFromCommandline;
class RestrictAbsentCanonicalAAS;
class DisallowIfNonnative;
class ReadResfile;
class SetRotamerCouplings;
class AppendRotamer;
class AppendRotamerSet;
class PreserveCBeta;
class PreventRepacking;
class IncludeCurrent;
class RotamerExplosion;
class ExtraChiCutoff;

typedef utility::pointer::shared_ptr< RestrictToRepacking > RestrictToRepackingOP;
typedef utility::pointer::shared_ptr< RestrictResidueToRepacking > RestrictResidueToRepackingOP;
typedef utility::pointer::shared_ptr< RestrictAbsentCanonicalAAS > RestrictAbsentCanonicalAASOP;
typedef utility::pointer::shared_ptr< DisallowIfNonnative > DisallowIfNonnativeOP;
typedef utility::pointer::shared_ptr< InitializeFromCommandline > InitializeFromCommandlineOP;
typedef utility::pointer::shared_ptr< ReadResfile > ReadResfileOP;
typedef utility::pointer::shared_ptr< SetRotamerCouplings > SetRotamerCouplingsOP;
typedef utility::pointer::shared_ptr< AppendRotamer > AppendRotamerOP;
typedef utility::pointer::shared_ptr< AppendRotamerSet > AppendRotamerSetOP;
typedef utility::pointer::shared_ptr< PreserveCBeta > PreserveCBetaOP;
typedef utility::pointer::shared_ptr< PreventRepacking > PreventRepackingOP;
typedef utility::pointer::shared_ptr< IncludeCurrent > IncludeCurrentOP;
typedef utility::pointer::shared_ptr< RotamerExplosion > RotamerExplosionOP;
typedef utility::pointer::shared_ptr< ExtraChiCutoff > ExtraChiCutoffOP;

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
