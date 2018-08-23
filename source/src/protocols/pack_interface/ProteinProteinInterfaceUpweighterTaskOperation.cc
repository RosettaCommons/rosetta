// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/pack_interface/ProteinProteinInterfaceUpweighterTaskOperation.cc
/// @brief upweight the interface energies
/// @author Longxing (longxing@uw.edu)

#include <protocols/pack_interface/ProteinProteinInterfaceUpweighterTaskOperationCreator.hh>
#include <protocols/pack_interface/ProteinProteinInterfaceUpweighterTaskOperation.hh>

#include <protocols/toolbox/IGEdgeReweighters.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/IGEdgeReweightContainer.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <utility/pointer/access_ptr.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace pack_interface {

using namespace core::pack::task::operation;
using namespace utility::tag;

static basic::Tracer tr( "protocols.pack_interface.ProteinProteinInterfaceUpweighterTaskOperation" );


core::pack::task::operation::TaskOperationOP
ProteinProteinInterfaceUpweighterTaskOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new ProteinProteinInterfaceUpweighter );
}

void ProteinProteinInterfaceUpweighterTaskOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ProteinProteinInterfaceUpweighter::provide_xml_schema( xsd );
}

std::string ProteinProteinInterfaceUpweighterTaskOperationCreator::keyname() const
{
	return ProteinProteinInterfaceUpweighter::keyname();
}


ProteinProteinInterfaceUpweighter::ProteinProteinInterfaceUpweighter()
{
	ppi_packer_weight_=1.0;
}

ProteinProteinInterfaceUpweighter::ProteinProteinInterfaceUpweighter(core::Real weight_in)
{
	ppi_packer_weight_= weight_in;
}

ProteinProteinInterfaceUpweighter::~ProteinProteinInterfaceUpweighter() = default;

core::pack::task::operation::TaskOperationOP
ProteinProteinInterfaceUpweighter::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new ProteinProteinInterfaceUpweighter(*this) );
}

void
ProteinProteinInterfaceUpweighter::parse_tag( TagCOP tag , DataMap & )
{
	if ( tag->hasOption("interface_weight") ) ppi_packer_weight_ = tag->getOption< core::Real >( "interface_weight", 1.0 );
}

void ProteinProteinInterfaceUpweighter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes.push_back( XMLSchemaAttribute::attribute_w_default(  "interface_weight", xsct_real, "scaling factor for the interaction energy of the residue pairs across the interface",  "1.0"  ) );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "scaling factor for the interaction energy of the residue pairs across the interface." );
}

/// @brief Change a packer task in some way.  The input pose is the one to which the input
/// task will be later applied.
void ProteinProteinInterfaceUpweighter::apply(
	Pose const &,
	PackerTask & task) const
{
	if ( ppi_packer_weight_ != 1.0 ) {
		core::pack::task::IGEdgeReweighterOP ppi_up( new protocols::toolbox::IGInterfaceEdgeUpweighter( ppi_packer_weight_ ) );
		core::pack::task::IGEdgeReweightContainerOP IGreweight = task.set_IGEdgeReweights();
		IGreweight->add_reweighter( ppi_up );

		tr.Info << "Packer Energies between the interface residues are upweighted by factor " << ppi_packer_weight_ << "." << std::endl;

	}

} //apply

}//namespace pack_interface
}//namespace protocols
