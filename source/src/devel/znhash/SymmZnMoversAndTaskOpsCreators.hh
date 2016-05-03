// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/znhash/SymmZnMoversAndTaskOpsCreators.hh
/// @brief  Declaration of creator classes for the mover and task operation factories
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Bryan Der (bder@email.unc.edu)

#ifndef INCLUDED_devel_znhash_SymmZnMoversAndTaskOpsCreators_HH
#define INCLUDED_devel_znhash_SymmZnMoversAndTaskOpsCreators_HH

#include <protocols/moves/MoverCreator.hh>
#include <core/pack/task/operation/TaskOperationCreator.hh>


namespace devel {
namespace znhash {

class InsertZincCoordinationRemarkLinesCreator : public protocols::moves::MoverCreator
{
public:
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
};

class DisableZnCoordinationResiduesTaskOpCreator : public core::pack::task::operation::TaskOperationCreator
{
public:
	virtual core::pack::task::operation::TaskOperationOP create_task_operation() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
	virtual std::string keyname() const;

};

class LoadZnCoordNumHbondCalculatorMoverCreator : public protocols::moves::MoverCreator
{
public:
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
};

}
}

#endif
