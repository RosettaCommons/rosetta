// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/pack_interface/ProteinProteinInterfaceUpweighterTaskOperation.hh
/// @brief upweight protein interface energies
/// @author longxing (longxing@uw.edu)


#ifndef INCLUDED_protocols_pack_interface_ProteinProteinInterfaceUpweighterTaskOperations_hh
#define INCLUDED_protocols_pack_interface_ProteinProteinInterfaceUpweighterTaskOperations_hh

#include <protocols/pack_interface/ProteinProteinInterfaceUpweighterTaskOperation.fwd.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>
#include <set>
#include <string>


namespace protocols {
namespace pack_interface {

/// @brief Class to alter a packer task to speficially upweight the protein-protein interaction energies
class ProteinProteinInterfaceUpweighter: public core::pack::task::operation::TaskOperation
{

public:
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;

public:

	ProteinProteinInterfaceUpweighter();
	ProteinProteinInterfaceUpweighter(core::Real weight_in);
	~ProteinProteinInterfaceUpweighter() override;

	TaskOperationOP clone() const override;

	/// @brief Change a packer task in some way.  The input pose is the one to which the input
	/// task will be later applied.
	void apply( Pose const & pose, PackerTask & task) const override;

	void parse_tag( TagCOP, DataMap & ) override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "ProteinProteinInterfaceUpweighter"; }

	core::Real get_weight() const { return ppi_packer_weight_;}

	/// @brief set the weight of the interaction across the interface
	void set_weight(core::Real const weight_in) {ppi_packer_weight_ = weight_in;}

private:
	/// Reweight protein-protein interaction by a factor of ppi_packer_weight_.
	core::Real ppi_packer_weight_;
};


}//namespace pack_interface
}//namespace protocols

#endif
