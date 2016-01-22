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

//unit headers
#include <protocols/flxbb/LayerDesignOperation.fwd.hh>

//protocols headers
#include <protocols/jd2/parser/BluePrint.fwd.hh>
#include <core/select/util/SelectResiduesByLayer.fwd.hh>

//core headers
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/types.hh>

//utility headers
#include <utility/tag/Tag.fwd.hh>
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
	typedef std::map< std::string, std::string > LayerDefinitions;
	typedef std::pair< std::string, std::string > LayerDefinition;
	typedef std::map< std::string,  LayerDefinitions > LayerResidues;
	typedef std::pair< std::string, LayerDefinitions > Layer;

	/// @brief List of noncanonicals that can be used in each secondary structure type.
	/// @details This maps secondary structure string ("L","E","H") to allowed NCAAs (vector of strings of three-letter codes).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	typedef std::map < std::string, utility::vector1 <std::string> > LayerNCDefinitions;
	/// @brief List of LayerNCDefinitions that can be used in each layer.
	/// @details This maps layer string ("core","boundary","surface") to LayerNCDefinitions (which in turn map secondary structure
	/// strings to vectors of noncanonical residue three-letter codes).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	typedef std::map < std::string, LayerNCDefinitions > LayerNCResidues;

	static utility::vector1< std::string > const SS_TYPES;

public:

	/// @brief default constructor
	///
	LayerDesignOperation();

	/// @brief value constructor
	///
	LayerDesignOperation( bool dsgn_core, bool dsgn_boundary, bool dsgn_surface );

	/// @brief Copy constructor.
	///
	LayerDesignOperation( LayerDesignOperation const & rval );

	/// @brief destructor
	virtual ~LayerDesignOperation();

	/// @brief make clone
	virtual TaskOperationOP clone() const;

public:

	/// @brief add a new layer which will be created based on a task operation
	void add_layer(
		std::string const & layer_name,
		core::pack::task::operation::TaskOperationOP task,
		LayerOperationType const operation,
		LayerSpecificationType const specification );

	/// @brief gets the residues allowed for a given secondary structure in a given layer
	std::string const & layer_residues(
		std::string const & layer_name,
		std::string const & ss_name ) const;

	/// @brief gets the residues allowed for a given secondary structure in a given layer
	utility::vector1< std::string > const & layer_nc_residues(
		std::string const & layer_name,
		std::string const & ss_name ) const;

	/// @brief sets residues allowed for a particular secondary structure in the given layer
	void set_layer_residues(
		std::string const & layer_name,
		std::string const & ss_name,
		std::string const & residues );

	/// @brief copies residues allowed from src to dest
	void copy_layer_residues(
		std::string const & src_layer,
		std::string const & dest_layer );

	/// @brief sets nc residues allowed for a particular secondary structure in the given layer
	void set_nc_layer_residues(
		std::string const & layer_name,
		std::string const & ss_name,
		std::string const & residues );
	void set_nc_layer_residues(
		std::string const & layer_name,
		std::string const & ss_name,
		utility::vector1< std::string > const & residues );

	/// @brief sets layer operation
	void set_layer_operation( std::string const & layer_name, LayerOperationType const operation );

	/// @brief sets layer specification
	void set_layer_specification( std::string const & layer_name, LayerSpecificationType const specification );

	/// @brief layer to be designed in internal SRBL object
	void design_layer( bool const dsgn_core, bool const dsgn_boundary, bool const dsgn_surface );

	/// @brief accessible surface for evaluating residues are in surface or not
	void sasa_surface( Real const r, String const ss="" );

	/// @brief accessible surface for evaluating residues are in core or not
	void sasa_core( Real const r, String const ss="" );

	/// @brief set pore radius for colculating asa
	void pore_radius( Real ps );

	/// @brief set verbose
	void set_verbose( bool const b ) { verbose_ = b; }

	/// @brief set restrict_restypes
	void set_restrict_restypes( bool const b ) { restrict_restypes_ = b; }

	/// @brief Set whether to use sidechain neighbors to determine core/boundary/surface (default=false)
	void set_use_sidechain_neighbors( bool const value );

	/// @brief Set the midpoint of the sigmoidal distance falloff for the sidechain neighbors method.
	///
	void set_sc_neighbor_dist_midpoint( core::Real const &value );

	/// @brief Set the factor by which neighbor counts are divided when using the sidechain neighbors method.
	///
	void set_sc_neighbor_denominator( core::Real const &value );

	/// @brief Set a parameter in the calculation that the sidechain neighbors algorithm uses.
	/// @details See core::select::util::SelectResiduesByLayer class for details.
	void set_sc_neighbor_angle_shift_factor( core::Real const &value );

	/// @brief Set another parameter (the angle exponent) in the calculation that the sidechain neighbors algorithm uses.
	/// @details See core::select::util::SelectResiduesByLayer class for details.
	void set_sc_neighbor_angle_exponent( core::Real const &value );

	/// @brief Set another parameter (the distance exponent) in the calculation that the sidechain neighbors algorithm uses.
	/// @details See core::select::util::SelectResiduesByLayer class for details.
	void set_sc_neighbor_dist_exponent( core::Real const &value );

	/// @brief sets names of layers to design
	void set_design_layers( utility::vector1< std::string > const & layers );

	/// @brief use original sequence for not designed layer
	void use_original_seq() { use_original_ = true; }

	/// @brief make pymol scripts showing the different layers
	void make_pymol_script( bool value ) { make_pymol_script_ = value; }

public:

	void parse_tag( TagCOP tag, DataMap & );

	void parse_layer_tag( TagCOP layer_tag, DataMap & datamap );

	void parse_layer_secstruct_tag(
		TagCOP secstruct_tag,
		DataMap & datamap,
		std::string const & layer_name );

	/// @brief Sets whether residues set with PIKAA and NATRO in previously-applied resfiles are ignored.
	/// @brief Default behaviour is now FALSE to preserve commutativity.
	void set_ignore_pikaa_natro( bool const val ) { ignore_pikaa_natro_ = val; return; }

	/// @brief Returns whether residues set with PIKAA and NATRO in previously-applied resfiles are ignored.
	/// @brief Default behaviour is now FALSE to preserve commutativity.
	bool ignore_pikaa_natro() const { return ignore_pikaa_natro_; }

public:

	/// @brief apply
	virtual void apply( Pose const & pose, PackerTask & task ) const;

private:
	/// @brief gets a list of allowed amino acids, and includes the current amino acid in addition to those allowed by layer
	utility::vector1< bool >
	get_restrictions(
		std::string const & layer,
		std::string const & default_layer,
		std::string const & ss_type,
		char const current_aa ) const;

	utility::vector1< bool >
	get_restrictions(
		std::string const & layer,
		std::string const & /*default_layer*/,
		std::string const & ss_type) const;

	void set_default_layer_residues();

	/// @brief Take a string consisting of comma-separated three-letter codes and parse it, storing separate three-letter codes in a given
	/// utility::vector1 of strings.
	void parse_ncaa_list(
		std::string const & str,
		utility::vector1< std::string > & storage_vect );

	/// @brief Remove a list of residue types from another list.
	///
	void exclude_ncaas(
		utility::vector1< std::string > const & aas_to_exclude,
		utility::vector1< std::string > & storage_vect );

	void set_layer_residues(
		std::string const & layer_name,
		std::string const & ss_name,
		std::string const & residues,
		LayerResidues & layer_residues );

	void init_nc_layerdefinitions( std::string const & layer_name );

	/// @brief Set whether this TaskOperation uses a whole, symmetric pose to define layers (true) or
	/// just the asymmetric unit (false).  Default is true.
	void set_use_symmetry( bool const val ) { use_symmetry_=val; return; }

	/// @brief Get whether this TaskOperation uses a whole, symmetric pose to define layers (true) or
	/// just the asymmetric unit (false).
	inline bool use_symmetry() const { return use_symmetry_; }

private:

	/// @brief utility function to transform a vector of position into a pymol selection command
	std::string pos2select( utility::vector1< Size > const & pos ) const;

	/// @brief write a pymol command with the different layers as selections
	void write_pymol_script(
		Pose const & pos,
		core::select::util::SelectResiduesByLayerOP srbl,
		std::map< std::string, utility::vector1< bool > > const & layer_specification,
		bool las_ligand,
		std::string const & filename ) const;

	/// @brief add helix capping ?
	bool add_helix_capping_;

	/// @brief use original sequence for not designed layer ?
	bool use_original_;

	/// @brief
	bool repack_non_designed_residues_;

	bool verbose_;

	bool restrict_restypes_;

	bool make_pymol_script_;

	LayerResidues layer_residues_;

	/// @brief Map of string for layer ("core", "boundary", "surface") to LayerNCDefinitions map (which in turn maps secondary structure to list of allowed noncanonicals).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	LayerNCResidues layer_nc_residues_;

	/// @brief Should residues for which PIKAA or NATRO information was defined in a previously-applied resfile
	/// be ignored by LayerDesign?  Default is now FALSE.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	bool ignore_pikaa_natro_;

	std::map< std::string, bool > design_layer_;
	std::map< std::string, LayerSpecificationType > layer_specification_;
	std::map< std::string, LayerOperationType > layer_operation_;

	TaskLayers task_layers_;

	// define the layer each residue belong to
	core::select::util::SelectResiduesByLayerOP srbl_;

	// for defining secondary structure
	protocols::jd2::parser::BluePrintOP blueprint_;

	/// @brief Should LayerDesign work with symmetry?  Default false.
	/// @details If false, the asymmetric unit is extracted and used for layer setup.
	/// If true (default), the whole, symmetric pose is used for layer setup.
	bool use_symmetry_;
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
