// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/flxbb/LayerDesignOperation.hh
/// @brief Design residues with selected amino acids depending on the enviroment: layer.
/// The layer of each residue is assigned as core, boundary, or surface, which are defined by
/// accessible surface of mainchain + CB. If resfile is read before calling this operation,
/// this operation is not applied for the residues defined by PIKAA.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )
/// @modified Javier Castellanos (javiercv@uw.edu )

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



#ifndef INCLUDED_protocols_flxbb_LayerDesignOperation_hh
#define INCLUDED_protocols_flxbb_LayerDesignOperation_hh

#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/AA.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/flxbb/LayerDesignOperation.fwd.hh>
#include <protocols/toolbox/SelectResiduesByLayer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>
// AUTO-REMOVED #include <iostream>
#include <utility/vector1.hh>

#include <string>
#include <map>


// Utility Headers

using namespace core::pack::task;

namespace protocols {
namespace flxbb {

class LayerDesignOperation : public core::pack::task::operation::TaskOperation {
public:

	enum LayerOperationType {
		DESIGN,
		NO_DESIGN,
		OMIT
	};

	enum LayerSpecificationType {
		DESIGNABLE,
		PACKABLE,
		FIXED
	};

	typedef std::string String;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
	typedef utility::tag::TagCOP TagCOP;
	typedef std::map< std::string, TaskOperationOP > TaskLayers;
	typedef std::pair< std::string, TaskOperationOP > TaskLayer;

	// Layer Residues is a map of maps, the first key is the layer(core, boundary, intermediate)
	// and the second key is the secondary structure (L,E,H). The values are string of one letter
	// code aminoacids to be used in each layer.
	typedef std::map<  std::string, std::string > LayerDefinitions;
	typedef std::pair< std::string, std::string > LayerDefinition;
	typedef std::map< std::string,  LayerDefinitions > LayerResidues;
	typedef std::pair< std::string, LayerDefinitions > Layer;

public:


	/// @brief default constructor
	LayerDesignOperation();

	/// @brief value constructor
	LayerDesignOperation( bool dsgn_core, bool dsgn_boundary, bool dsgn_surface );

	LayerDesignOperation( LayerDesignOperation const & rval );

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

	/// @brief set restrict_restypes
	void set_restrict_restypes( bool const  b ) { restrict_restypes_ = b; }

	/// @brief use original sequence for not designed layer
	void use_original_seq()
	{
		use_original_ = true;
	}

	///@brief make pymol scripts showing the different layers
	void  make_pymol_script(bool value) { make_pymol_script_ = value; }

public:


	void parse_tag( TagCOP tag , DataMap & );


public:


	/// @brief apply
	virtual void apply( Pose const & pose, PackerTask & task ) const;

private:
	utility::vector1<bool> get_restrictions(std::string const & layer, std::string const & default_layer, std::string const &  ss_type) const;
	void set_default_layer_residues();


private:

	///@brief utility function to transform a vector of position into a pymol selection command
	std::string pos2select( utility::vector1< Size > const & pos) const;
	///@brief write a pymol command with the different layers as selections
	void write_pymol_script( Pose const & pos, toolbox::SelectResiduesByLayerOP srbl, std::map< std::string, utility::vector1<bool> > const & layer_specification,bool las_ligand, std::string const & filename ) const;

	/// @brief add helix capping ?
	bool add_helix_capping_;

	/// @brief use original sequence for not designed layer ?
	bool use_original_;

	///@brief
	bool repack_non_designed_residues_;

	bool verbose_;

	bool restrict_restypes_;

	bool make_pymol_script_;

	LayerResidues layer_residues_;

	std::map< std::string, bool > design_layer_;
	std::map< std::string, LayerSpecificationType > layer_specification_;
	std::map< std::string, LayerOperationType > layer_operation_;

	TaskLayers task_layers_;

	// define the layer each residue belong to
	toolbox::SelectResiduesByLayerOP srbl_;

	// for defining secondary structure
	protocols::jd2::parser::BluePrintOP blueprint_;
};

// utility class for chaining together task operations
class CombinedTaskOperation : public core::pack::task::operation::TaskOperation {
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef utility::vector1< TaskOperationOP > VecTaskOP;
public:
	CombinedTaskOperation(VecTaskOP ops);
	/// @brief apply
	virtual void apply( Pose const & pose, PackerTask & task ) const;
	/// @brief make clone
	virtual TaskOperationOP clone() const { return TaskOperationOP( new CombinedTaskOperation( *this ) ); }

private:
	VecTaskOP task_operations_;
};

typedef utility::pointer::shared_ptr< CombinedTaskOperation > CombinedTaskOperationOP;

} // flxbb
} // protocols

#endif
