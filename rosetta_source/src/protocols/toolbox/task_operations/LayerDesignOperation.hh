// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/task_operations/LayerDesignOperation.hh
/// @brief Design residues with selected amino acids depending on the enviroment: layer.
/// The layer of each residue is assigned as core, boundary, or surface, which are defined by
/// accessible surface of mainchain + CB. If resfile is read before calling this operation,
/// this operation is not applied for the residues defined by PIKAA.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

//  The following are using amino acids for each layer
/// @CORE
//    Loop: AFILPVWY
//  Strand:  FIL VWY
//   Helix: AFIL VWY ( P only at the beginning of helix )
//   HelixCapping: DNST
//
/// @BOUDNARY
//    Loop: ADEFGIKLNPQRSTVWY
//  Strand:  DEF IKLN QRSTVWY
//   Helix: ADE  IKLN QRSTV Y ( P only at the beginning of helix )
//   HelixCapping: DNST
//
/// @SURFACE
//    Loop: DEGHKNPQRST
//  Strand: DE HKN QRST
//   Helix: DE HKN QRST  ( P only at the beginning of helix )
//   HelixCapping: DNST



#ifndef INCLUDED_protocols_toolbox_task_operations_LayerDesignOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_LayerDesignOperation_hh

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/toolbox/task_operations/LayerDesignOperation.fwd.hh>
#include <protocols/flxbb/SelectResiduesByLayer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <iostream>
#include <utility/vector1.hh>

// Utility Headers

using namespace core::pack::task;

namespace protocols{
namespace toolbox{
namespace task_operations{

class LayerDesignOperation : public core::pack::task::operation::TaskOperation {
public:


	typedef std::string String;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
	typedef utility::tag::TagPtr TagPtr;

	typedef protocols::flxbb::SelectResiduesByLayer SelectResiduesByLayer;
	typedef protocols::flxbb::SelectResiduesByLayerOP SelectResiduesByLayerOP;

public:


	/// @brief default constructor
	LayerDesignOperation();

	/// @brief value constructor
	LayerDesignOperation( bool dsgn_core, bool dsgn_boundary, bool dsgn_surface );

	/// @brief destructor
	virtual ~LayerDesignOperation();

	/// @brief make clone
	virtual TaskOperationOP clone() const;


public:


	/// @brief layer to be designed
	void design_layer( bool const dsgn_core, bool const dsgn_boundary, bool const dsgn_surface );

	/// @brief accessible surface for evaluating residues are in surface or not
	void sasa_surface( Real const r, String const ss="" );

	/// @brief accessible surface for evaluating residues are in core or not
	void sasa_core( Real const r, String const ss="" );

	/// @brief set pore radius for colculating asa
	void pore_radius( Real ps );

	/// @brief set verbose
	void set_verbose( bool const  b ) { verbose_ = b; }

	/// @brief use original sequence for not designed layer
	void use_original_seq()
	{
		use_original_ = true;
	}


public:


	void parse_tag( TagPtr tag );


public:


	/// @brief apply
	virtual void apply( Pose const & pose, PackerTask & task ) const;


private:


	/// @brief add helix capping ?
	bool add_helix_capping_;

	/// @brief use original sequence for not designed layer ? otherwise the residues will be changed to ala
	bool use_original_;

	/// @brief
	bool verbose_;

	// define the layer each residue belong to
	SelectResiduesByLayerOP srbl_;


};


} // TaskOperations
} // toolbox
} // protocols
#endif
