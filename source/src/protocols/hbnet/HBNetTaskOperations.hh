// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/hbnet/HBNetTaskOperations.cc
/// @brief task operations for HBNet
/// @author Scott Boyken (sboyken@gmail.com)

#ifndef INCLUDED_protocols_hbnet_HBNetTaskOperations_hh
#define INCLUDED_protocols_hbnet_HBNetTaskOperations_hh

#include <protocols/hbnet/HBNetTaskOperations.fwd.hh>
//#include <protocols/motifs/LigandMotifSearch.fwd.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>
#include <set>
#include <string>


namespace protocols {
namespace hbnet {


/// @brief sets allowed residue types to constrain HBNet residues in downstream design;
///         add this Taskop to any design movers downstream of HBNet
class ConstrainHBondNetwork : public core::pack::task::operation::TaskOperation
{

public:
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;

public:

	ConstrainHBondNetwork();
	~ConstrainHBondNetwork();

	virtual TaskOperationOP clone() const;

	virtual void apply( Pose const & pose, PackerTask & task ) const;

	virtual void parse_tag( TagCOP tag , DataMap & );

	static std::string keyname();

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};


}//namespace hbnet
}//namespace protocols

#endif
