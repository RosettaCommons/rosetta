// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/TaskOperationCreators.hh
/// @brief  Declaration for the class that connects TaskOperations with the TaskOperationFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_TaskOperationCreators_hh
#define INCLUDED_core_pack_task_operation_TaskOperationCreators_hh

#include <core/pack/task/operation/TaskOperationCreator.hh>

#include <core/pack/task/operation/TaskOperation.fwd.hh>

#include <string>


namespace core {
namespace pack {
namespace task {
namespace operation {

class RestrictYSDesignCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "RestrictYSDesign"; }
};

class PreventRepackingCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "PreventRepacking"; }
};

class PreserveCBetaCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "PreserveCBeta"; }
};

class AppendRotamerSetCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "AppendRotamerSet"; }
};

class AppendResidueRotamerSetCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "AppendResidueRotamerSet"; }
};

class AppendRotamerCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "AppendRotamer"; }
};

class SetRotamerCouplingsCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "SetRotamerCouplings"; }
};

class SetRotamerLinksCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "SetRotamerLinks"; }
};


class ReadResfileCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "ReadResfile"; }
};

class ReadResfileAndObeyLengthEventsCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "ReadResfileAndObeyLengthEvents"; }
};


class IncludeCurrentCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "IncludeCurrent"; }
};


class InitializeExtraRotsFromCommandlineCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "InitializeExtraRotsFromCommandline"; }
};

class InitializeFromCommandlineCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "InitializeFromCommandline"; }
};

class ExtraRotamersGenericCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "ExtraRotamersGeneric"; }
};

class RotamerExplosionCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "RotamerExplosionCreator"; }
};

class RestrictAbsentCanonicalAASCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "RestrictAbsentCanonicalAAS"; }
};

class RestrictResidueToRepackingCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "RestrictResidueToRepacking"; }
};

class RestrictToRepackingCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "RestrictToRepacking"; }
};

class DisallowIfNonnativeCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "DisallowIfNonnative"; }
};

class ExtraRotamersCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "ExtraRotamers"; }
};

class ExtraChiCutoffCreator : public TaskOperationCreator {
public:
	virtual TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "ExtraChiCutoff"; }
};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
