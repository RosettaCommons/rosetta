// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
	// XRW TEMP  protocols::moves::MoverOP create_mover() const override;
	// XRW TEMP  std::string keyname() const override;
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

class DisableZnCoordinationResiduesTaskOpCreator : public core::pack::task::operation::TaskOperationCreator
{
public:
	core::pack::task::operation::TaskOperationOP create_task_operation() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
	std::string keyname() const override;

};

class LoadZnCoordNumHbondCalculatorMoverCreator : public protocols::moves::MoverCreator
{
public:
	// XRW TEMP  protocols::moves::MoverOP create_mover() const override;
	// XRW TEMP  std::string keyname() const override;
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

}
}

#endif
