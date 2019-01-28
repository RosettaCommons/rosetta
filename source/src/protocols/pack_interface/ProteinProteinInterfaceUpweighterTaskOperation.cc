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
#include <core/scoring/dssp/Dssp.hh>
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
	return utility::pointer::make_shared< ProteinProteinInterfaceUpweighter >();
}

void ProteinProteinInterfaceUpweighterTaskOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ProteinProteinInterfaceUpweighter::provide_xml_schema( xsd );
}

std::string ProteinProteinInterfaceUpweighterTaskOperationCreator::keyname() const
{
	return ProteinProteinInterfaceUpweighter::keyname();
}


ProteinProteinInterfaceUpweighter::ProteinProteinInterfaceUpweighter() :
	ppi_packer_weight_(1.0),
	chains_to_ignore_loops_("")
{
}

ProteinProteinInterfaceUpweighter::ProteinProteinInterfaceUpweighter(core::Real weight_in) :
	ppi_packer_weight_(weight_in),
	chains_to_ignore_loops_("")
{
}

ProteinProteinInterfaceUpweighter::ProteinProteinInterfaceUpweighter(core::Real weight_in, std::string const & skip_loop_in_chain) :
	ppi_packer_weight_(weight_in),
	chains_to_ignore_loops_(skip_loop_in_chain)
{
}

ProteinProteinInterfaceUpweighter::~ProteinProteinInterfaceUpweighter() = default;

core::pack::task::operation::TaskOperationOP
ProteinProteinInterfaceUpweighter::clone() const
{
	return utility::pointer::make_shared< ProteinProteinInterfaceUpweighter >(*this);
}

void
ProteinProteinInterfaceUpweighter::parse_tag( TagCOP tag , DataMap & )
{
	if ( tag->hasOption("interface_weight") ) ppi_packer_weight_ = tag->getOption< core::Real >( "interface_weight", 1.0 );
	if ( tag->hasOption("skip_loop_in_chain") ) chains_to_ignore_loops_ = tag->getOption< std::string >( "skip_loop_in_chain", "");
}

void ProteinProteinInterfaceUpweighter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes.push_back( XMLSchemaAttribute::attribute_w_default(  "interface_weight", xsct_real, "scaling factor for the interaction energy of the residue pairs across the interface",  "1.0"  ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "skip_loop_in_chain", xs_string, "comma seperated chain idenfitiers for chains that the weight related to its loops are retained to 1.0", "" ) );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "scaling factor for the interaction energy of the residue pairs across the interface." );
}

/// @brief Change a packer task in some way.  The input pose is the one to which the input
/// task will be later applied.
void ProteinProteinInterfaceUpweighter::apply(
	Pose const & pose,
	PackerTask & task) const
{
	if ( ppi_packer_weight_ != 1.0 ) {
		std::string sec_str = "";
		if ( chains_to_ignore_loops_.size() != 0 ) {
			core::scoring::dssp::Dssp dssp( pose );
			sec_str = dssp.get_dssp_secstruct();
		}
		core::pack::task::IGEdgeReweighterOP ppi_up( utility::pointer::make_shared<protocols::toolbox::IGInterfaceEdgeUpweighter>( ppi_packer_weight_, chains_to_ignore_loops_, sec_str ) );
		core::pack::task::IGEdgeReweightContainerOP IGreweight = task.set_IGEdgeReweights();
		IGreweight->add_reweighter( ppi_up );

		tr.Info << "Packer Energies between the interface residues are upweighted by factor " << ppi_packer_weight_ << "." << std::endl;

	}

} //apply

}//namespace pack_interface
}//namespace protocols
