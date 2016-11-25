// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/flxbb/LayerDesignOperation.cc
/// @brief Design residues with selected amino acids depending on the enviroment: layer.
/// The layer of each residue is assigned as core, boundary, or surface, which are defined by
/// accessible surface of mainchain + CB
/// @details If ignore_pikaa_natro(true) is invoked, then if a resfile is read before calling this
/// operation, this operation is not applied for the residues defined by PIKAA.  Note that this
/// breaks TaskOperation commutativity, so it is NOT the default behaviour.  In RosettaScripts,
/// the option ignore_pikaa_natro=true can be added to the <LayerDesign> tag to turn this on.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )
/// @modified Javier Castellanos (javiercv@uw.edu )
/// @modified Vikram K. Mulligan (vmullig@uw.edu) -- support for noncanonicals.

//  The following are using amino acids for each layer
/// @CORE
//    Loop: AFILPVWY
//  Strand:  FIL VWY
//   Helix: AFIL VWY ( P only at the beginning of helix )
//   HelixCapping: DNST
//
/// @BOUNDARY
//    Loop: ADEFGIKLNPQRSTVWY
//  Strand:  DEF IKLN QRSTVWY
//   Helix: ADE  IKLN QRSTVWY ( P only at the beginning of helix )
//   HelixCapping: DNST
//
/// @SURFACE
//    Loop: DEGHKNPQRST
//  Strand: DE HKN QRST
//   Helix: DE HKN QRST  ( P only at the beginning of helix )
//   HelixCapping: DNST


// Unit Headers
#include <protocols/flxbb/LayerDesignOperation.hh>
#include <protocols/flxbb/LayerDesignOperationCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <basic/Tracer.hh>
#include <core/select/util/SelectResiduesByLayer.hh>
#include <core/pose/symmetry/util.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/jd2/parser/BluePrint.hh>

// Utility Headers
#include <protocols/flxbb/utility.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <utility/vector1.hh>
#include <string>
#include <algorithm>

using basic::T;
using basic::Error;
using basic::Warning;


#include <utility/vector0.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <ObjexxFCL/format.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <boost/assign/list_inserter.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;
using namespace core::pack::task::operation;
using namespace utility::tag;

static THREAD_LOCAL basic::Tracer TR( "protocols.flxbb.LayerDesignOperation" );

// AMW: Note--there's no Creator for CombinedTaskOperation
// So I'm not doing anything with it yet.

namespace protocols {
namespace flxbb {

utility::vector1< std::string > const LayerDesignOperation::SS_TYPES =
boost::assign::list_of ("all")("Helix")("HelixCapping")("HelixStart")("Loop")("Strand");

//bool
//valid_ss_type( std::string const & query ) {
// for ( auto ss : SS_TYPES ) {
//  if ( ss == query ) return true;
// }
// return false;
//}

TaskOperationOP
LayerDesignOperationCreator::create_task_operation() const
{
	return TaskOperationOP( new LayerDesignOperation );
}

void LayerDesignOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LayerDesignOperation::provide_xml_schema( xsd );
}

std::string LayerDesignOperationCreator::keyname() const
{
	return LayerDesignOperation::keyname();
}

CombinedTaskOperation::CombinedTaskOperation( VecTaskOP  ops ):
	task_operations_( ops )
{ }

void
CombinedTaskOperation::apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const {
	using core::pack::task::operation::TaskOperationOP;
	BOOST_FOREACH ( TaskOperationOP t_op, task_operations_ ) {
		t_op->apply( pose, task );
	}
}


/// @brief default constructor
///
LayerDesignOperation::LayerDesignOperation():
	TaskOperation(),
	add_helix_capping_( true ),
	use_original_( true ),
	repack_non_designed_residues_( true ),
	verbose_( false ),
	restrict_restypes_( true ),
	make_pymol_script_( false ),
	ignore_pikaa_natro_( false ),
	srbl_( core::select::util::SelectResiduesByLayerOP( new core::select::util::SelectResiduesByLayer ) ),
	blueprint_( /* NULL */ ),
	use_symmetry_(true)
{
	set_default_layer_residues();
}

/// @brief value constructor
///
LayerDesignOperation::LayerDesignOperation( bool dsgn_core, bool dsgn_boundary, bool dsgn_surface ):
	TaskOperation(),
	add_helix_capping_( true ),
	use_original_( true ),
	repack_non_designed_residues_( true ),
	verbose_( false ),
	restrict_restypes_( true ),
	make_pymol_script_( false ),
	ignore_pikaa_natro_( false ),
	srbl_( core::select::util::SelectResiduesByLayerOP( new core::select::util::SelectResiduesByLayer ) ),
	blueprint_( /* NULL */ ),
	use_symmetry_(true)
{
	design_layer( dsgn_core, dsgn_boundary, dsgn_surface );
	/*
	design_layer_[std::string("core")]=dsgn_core;
	design_layer_[std::string("boundary")]=dsgn_boundary;
	design_layer_[std::string("surface")]=dsgn_surface;
	*/
	set_default_layer_residues();
}

/// @brief value constructor
LayerDesignOperation::LayerDesignOperation( bool dsgn_core, bool dsgn_boundary, bool dsgn_surface, bool Nt, bool Ct, Real core_score, Real surf_score, std::string coreAA, std::string boundaryAA, std::string surfaceAA ):
	TaskOperation(),
	add_helix_capping_( true ),
	use_original_( true ),
	repack_non_designed_residues_( true ),
	verbose_( false ),
	restrict_restypes_( true ),
	make_pymol_script_( false ),
	srbl_( core::select::util::SelectResiduesByLayerOP( new core::select::util::SelectResiduesByLayer ) ),
	blueprint_( nullptr )
{
	design_layer( dsgn_core, dsgn_boundary, dsgn_surface );
	design_layer_[std::string("core")]=dsgn_core;
	design_layer_[std::string("boundary")]=dsgn_boundary;
	design_layer_[std::string("surface")]=dsgn_surface;
	design_layer_[std::string("Nterm")]=Nt;
	design_layer_[std::string("Cterm")]=Ct;
	srbl_->sasa_surface( surf_score, String(""));
	srbl_->sasa_core( core_score, String(""));

	set_default_layer_residues();

	//if these flags are set, swap out the default residies
	if ( coreAA != "" ) {
		layer_residues_[String("core")][String("all")] = coreAA;
		layer_residues_[String("core")][String("Loop")] = coreAA;
		layer_residues_[String("core")][String("Strand")] = coreAA;
		layer_residues_[String("core")][String("Helix")] = coreAA;
		//  layer_residues_[String("core")][String("HelixStart")] = String("P");
		//making sure Pro is available
		if ( coreAA.find("P") == std::string::npos ) coreAA = coreAA + "P";
		layer_residues_[String("core")][String("HelixStart")] = coreAA;
	}
	if ( boundaryAA != "" ) {
		layer_residues_[String("boundary")][String("all")] = boundaryAA;
		layer_residues_[String("boundary")][String("Loop")] = boundaryAA;
		layer_residues_[String("boundary")][String("Strand")] = boundaryAA;
		layer_residues_[String("boundary")][String("Helix")] = boundaryAA;
		// layer_residues_[String("boundary")][String("HelixStart")] = String("P");
		//making sure Pro is available
		if ( boundaryAA.find("P") == std::string::npos ) boundaryAA = boundaryAA + "P";
		layer_residues_[String("boundary")][String("HelixStart")] = boundaryAA;
	}
	if ( surfaceAA != "" ) {
		layer_residues_[String("surface")][String("all")] = surfaceAA;
		layer_residues_[String("surface")][String("Loop")] = surfaceAA;
		layer_residues_[String("surface")][String("Strand")] = surfaceAA;
		layer_residues_[String("surface")][String("Helix")] = surfaceAA;
		//  layer_residues_[String("surface")][String("HelixStart")] = String("P");
		//making sure Pro is available
		if ( surfaceAA.find("P") == std::string::npos ) surfaceAA = surfaceAA + "P";
		layer_residues_[String("surface")][String("HelixStart")] = surfaceAA;
	}



}

/// @brief Copy constructor.
///
LayerDesignOperation::LayerDesignOperation( LayerDesignOperation const & )= default;


/// @brief destructor
LayerDesignOperation::~LayerDesignOperation() = default;

/// @brief clone
core::pack::task::operation::TaskOperationOP
LayerDesignOperation::clone() const {
	return core::pack::task::operation::TaskOperationOP( new LayerDesignOperation( *this ) );
}

/// @brief layer to be designed
void
LayerDesignOperation::design_layer( bool const dsgn_core, bool const dsgn_boundary, bool const dsgn_surface )
{
	srbl_->set_design_layer( dsgn_core, dsgn_boundary, dsgn_surface );
}

/// @brief accessible surface for evaluating residues are in surface or not
void
LayerDesignOperation::sasa_surface( Real const r, String const & ss )
{
	srbl_->sasa_surface( r, ss );
}

/// @brief accessible surface for evaluating residues are in core or not
void
LayerDesignOperation::sasa_core( Real const r, String const & ss )
{
	srbl_->sasa_core( r, ss );
}

/// @brief set pore radius for colculating asa
void
LayerDesignOperation::pore_radius( Real ps )
{
	srbl_->pore_radius( ps );
}

std::string
LayerDesignOperation::pos2select( utility::vector1< core::Size > const & pos ) const
{
	std::string str;
	if ( pos.empty() ) {
		return str;
	} else {
		str = boost::lexical_cast< std::string >( pos[1] );
	}
	for ( Size i = 2; i <= pos.size(); i++ ) {
		str += "+" + boost::lexical_cast< std::string >( pos[i] );
	}
	return str;
}

void
LayerDesignOperation::write_pymol_script( core::pose::Pose const & pose, core::select::util::SelectResiduesByLayerOP srbl, std::map< std::string, utility::vector1<bool> > const & layer_specification, bool has_ligand, std::string const & filename ) const
{
	using utility::io::ozstream;
	typedef utility::vector1<Size> VecSize;
	TR << "Writing pymol script with the layer information to "<< filename << std::endl;
	ozstream pymol( filename );
	TR << "Dumping pose for the script as layer_design_input.pdb" << std::endl;
	pose.dump_pdb("layer_design_input.pdb");

	// importing necessary pymol modules
	pymol << "from pymol import cmd" << std::endl;
	//load strucuture into pymol
	pymol << "cmd.load('layer_design_input.pdb')" << std::endl;
	pymol << "cmd.hide('all') " << std::endl;
	pymol << "cmd.show('spheres')" << std::endl;
	// make the selections for the basic layers
	VecSize core_residues = srbl->selected_core_residues();
	pymol << "cmd.select('core', 'resi  " << pos2select( core_residues )<< "')" << std::endl;
	pymol << "cmd.color('red', 'core')" << std::endl;
	VecSize boundary_residues = srbl->selected_boundary_residues();
	pymol << "cmd.select('boundary', 'resi  " << pos2select( boundary_residues )<< "')" << std::endl;
	pymol << "cmd.color('orange', 'boundary')" << std::endl;
	VecSize surface_residues = srbl->selected_surface_residues();
	pymol << "cmd.select('surface', 'resi  " << pos2select( surface_residues ) << "')" <<  std::endl;
	pymol << "cmd.color('yellow', 'surface')" << std::endl;
	// make the selections for the task layers
	utility::vector0< std::string > colors;
	colors.push_back( "green" );
	colors.push_back( "magenta" );
	colors.push_back( "cyan" );
	colors.push_back( "purpleblue" );
	colors.push_back( "hotpink" );
	colors.push_back( "olive" );
	Size layer = 0;
	for ( auto const & it : layer_specification ) {
		utility::vector1< Size > pos;
		for ( Size i = 1; i <= pose.size(); i++ ) {
			if ( it.second[ i ] ) {
				pos.push_back( i );
			}
		}
		pymol << "cmd.select('"<< it.first  <<"', 'resi  " << pos2select( pos )<< "')" << std::endl;
		pymol << "cmd.color('" << colors[ layer % 6  ] << "','" << it.first << "')" << std::endl;
		layer += 1;
	}
	if ( has_ligand ) {
		pymol << "cmd.select('ligand', 'organic')" << std::endl;
		pymol << "cmd.show('sticks','ligand')" << std::endl;
		pymol << "cmd.color('gray','ligand')" << std::endl;
	}
}

/// @details Indirection to get around c++11 flexibility in initializing
/// maps.  Used in set_default_layer_residues below.  Idea stolen from
/// this thread on github:
/// https://github.com/ethz-asl/libpointmatcher/issues/13
template< class T >
LayerDesignOperation::LayerDefinitions
makeMap( T const & map_initializer )
{
	LayerDesignOperation::LayerDefinitions m = map_initializer;
	return m;
}

void
LayerDesignOperation::init_nc_layerdefinitions( std::string const & layer_name )
{
	LayerNCDefinitions & def = layer_nc_residues_[ layer_name ];
	if ( ( layer_name == "Nterm" ) || ( layer_name == "Cterm" ) ) {
		def[ "all" ];
	} else {
		for ( auto const & ss : SS_TYPES ) {
			def[ ss ];
		}
	}
}

void
LayerDesignOperation::set_default_layer_residues() {
	TR << "Initializing the layers with the default residues" << std::endl;
	boost::assign::insert(layer_residues_) //If the defaults are modified here, be sure to set up defaults below for the layer_nc_residues_ object that defines noncanonicals.
		(std::string("core"), makeMap( boost::assign::map_list_of
		(std::string("all"),            std::string("AFILPVWYDNST"))
		(std::string("Loop"),            std::string("AFILPVWY"))
		(std::string("Strand"),           std::string("FILVWY"))
		(std::string("Helix"),            std::string("AFILVWY"))
		(std::string("HelixStart"),       std::string("AFILVWYP"))
		(std::string("HelixCapping"),     std::string("DNST"))
		))
		(std::string("boundary"), makeMap( boost::assign::map_list_of
		(std::string("all"),              std::string("ADEFGIKLNPQRSTVWY"))
		(std::string("Loop"),             std::string("ADEFGIKLNPQRSTVWY"))
		(std::string("Strand"),           std::string("DEFIKLNQRSTVWY"))
		(std::string("Helix"),            std::string("ADEIKLNQRSTVWY"))
		(std::string("HelixStart"),       std::string("ADEIKLNQRSTVWYP"))
		(std::string("HelixCapping"),     std::string("DNST"))
		))
		(std::string("surface"), makeMap( boost::assign::map_list_of
		(std::string("all"),              std::string("DEGHKNPQRST"))
		(std::string("Loop"),             std::string("DEGHKNPQRST"))
		(std::string("Strand"),           std::string("DEHKNQRST"))
		(std::string("Helix"),            std::string("DEHKNQRST"))
		(std::string("HelixStart"),       std::string("DEHKNQRSTP"))
		(std::string("HelixCapping"),     std::string("DNST"))
		))
		(std::string("Nterm"), makeMap( boost::assign::map_list_of
		(std::string("all"),              std::string("ACDEFGHIKLMNPQRSTVWY"))
		//(std::string("all"),     std::string("DNST"))
		))
		(std::string("Cterm"), makeMap( boost::assign::map_list_of
		(std::string("all"),              std::string("ACDEFGHIKLMNPQRSTVWY"))
		//(std::string("all"),              std::string("DEGNPQRK"))
		));

	// Note: by default, the layer_nc_residues_ object contains empty NCAA lists, since we're not designing with noncanonicals by default. --VKM
	init_nc_layerdefinitions( "core" );
	init_nc_layerdefinitions( "boundary" );
	init_nc_layerdefinitions( "surface" );
	init_nc_layerdefinitions( "Nterm" );
	init_nc_layerdefinitions( "Cterm" );

	boost::assign::insert(design_layer_)
		(std::string("core"), false )
		(std::string("boundary"), false )
		(std::string("surface"), false )
		(std::string("Nterm"), false )
		(std::string("Cterm"), false );

	boost::assign::insert(layer_specification_)
		(std::string("core"), DESIGNABLE)
		(std::string("boundary"), DESIGNABLE)
		(std::string("surface"), DESIGNABLE)
		(std::string("Nterm"), DESIGNABLE)
		(std::string("Cterm"), DESIGNABLE);

	boost::assign::insert(layer_operation_)
		(std::string("core"), DESIGN)
		(std::string("boundary"), DESIGN)
		(std::string("surface"), DESIGN)
		(std::string("Nterm"), DESIGN)
		(std::string("Cterm"), DESIGN);
}

utility::vector1<bool>
LayerDesignOperation::get_restrictions( std::string const & layer, std::string const & /*default_layer*/, std::string const & ss_type) const {
	// if the layer doesn't specify the required ss used the default layer one.
	utility::vector1< bool >  restrict_to_aa( chemical::num_canonical_aas, false );
	BOOST_FOREACH ( char restype, layer_residues_.find( layer )->second.find( ss_type )->second ) {
		restrict_to_aa[ chemical::aa_from_oneletter_code( restype ) ] = true;
	}
	return restrict_to_aa;
}

/// @brief Take a string consisting of comma-separated three-letter codes and parse it, storing separate three-letter codes in a given
/// utility::vector1 of strings.
void LayerDesignOperation::parse_ncaa_list( std::string const &str, utility::vector1<std::string> &storage_vect ) {
	for ( core::Size i=0; i<str.length()-2; ++i ) {
		if ( str[i]!=',' && str[i]!=' ' && str[i]!='\n' && str[i]!='\t' ) {
			storage_vect.push_back( str.substr(i,3) );
			i+=2;
		}
	}
	return;
}

/// @brief Remove a list of residue types from another list.
///
void LayerDesignOperation::exclude_ncaas( utility::vector1<std::string> const &aas_to_exclude, utility::vector1<std::string> &storage_vect ) {
	for ( core::Size i=1,imax=aas_to_exclude.size(); i<=imax; ++i ) { //Loop through all aas to exclude.
		for ( auto j = storage_vect.begin(); j<storage_vect.end(); ++j ) {
			if ( (*j)==aas_to_exclude[i] ) {
				storage_vect.erase(j);
			}
		}
	}
	return;
}

void
LayerDesignOperation::set_design_layers( utility::vector1< std::string > const & layers )
{
	BOOST_FOREACH ( std::string const & layer, layers ) {
		design_layer_[ layer ] = true;
	}
}

/// @brief Utility function to convert a string vector into a comma-separated list of strings.
///
std::string print_string_vector( utility::vector1< std::string > const & vect )
{
	std::string outstring("");
	for ( core::Size i=1, imax=vect.size(); i<=imax; ++i ) {
		if ( i > 1 ) outstring += ",";
		outstring += vect[i];
	}
	return outstring;
}

void
LayerDesignOperation::set_use_sidechain_neighbors( bool const value )
{
	srbl_->use_sidechain_neighbors( value );
	if ( srbl_->use_sidechain_neighbors() ) {
		srbl_->sasa_core( 5.2 );
		srbl_->sasa_surface( 2.0 );
	}
}

/// @brief Set the midpoint of the sigmoidal distance falloff for the sidechain neighbors method.
///
void
LayerDesignOperation::set_sc_neighbor_dist_midpoint(
	core::Real const &value
) {
	srbl_->set_dist_midpoint(value);
	return;
}

/// @brief Set the factor by which neighbor counts are divided when using the sidechain neighbors method.
///
void
LayerDesignOperation::set_sc_neighbor_denominator(
	core::Real const &value
) {
	srbl_->set_rsd_neighbor_denominator(value);
	return;
}

/// @brief Set a parameter in the calculation that the sidechain neighbors algorithm uses.
/// @details See core::select::util::SelectResiduesByLayer class for details.
void
LayerDesignOperation::set_sc_neighbor_angle_shift_factor( core::Real const &value ) {
	srbl_->set_angle_shift_factor(value);
	return;
}

/// @brief Set another parameter (the angle exponent) in the calculation that the sidechain neighbors algorithm uses.
/// @details See core::select::util::SelectResiduesByLayer class for details.
void
LayerDesignOperation::set_sc_neighbor_angle_exponent( core::Real const &value ) {
	srbl_->set_angle_exponent(value);
	return;
}

/// @brief Set another parameter (the distance exponent) in the calculation that the sidechain neighbors algorithm uses.
/// @details See core::select::util::SelectResiduesByLayer class for details.
void
LayerDesignOperation::set_sc_neighbor_dist_exponent( core::Real const &value ) {
	srbl_->set_dist_exponent(value);
	return;
}


/// @brief sets layer operation
void
LayerDesignOperation::set_layer_operation( std::string const & layer_name, LayerOperationType const operation )
{
	layer_operation_[ layer_name ] = operation;
}

/// @brief sets layer specification
void
LayerDesignOperation::set_layer_specification( std::string const & layer_name, LayerSpecificationType const specification )
{
	layer_specification_[ layer_name ] = specification;
}

/// @brief layer to be designed
void
LayerDesignOperation::add_layer(
	std::string const & layer_name,
	core::pack::task::operation::TaskOperationOP task,
	LayerOperationType const operation,
	LayerSpecificationType const specification )
{
	task_layers_[ layer_name ] = task;
	design_layer_[ layer_name ] = true;

	layer_nc_residues_[ layer_name ] ;
	layer_operation_[ layer_name ] = operation;
	layer_specification_[ layer_name ] = specification;
	LayerDefinitions newdef;
	layer_residues_[ layer_name ] = std::map< std::string, std::string >();
}

/// @brief gets the residues allowed for a given secondary structure in a given layer
std::string const &
LayerDesignOperation::layer_residues(
	std::string const & layer_name,
	std::string const & ss_name ) const
{
	auto it = layer_residues_.find( layer_name );
	if ( it == layer_residues_.end() ) {
		throw utility::excn::EXCN_Msg_Exception( "Layer " + layer_name + " was not found when trying to determine valid residue types." );
	}
	debug_assert( it != layer_residues_.end() );

	auto it2 = it->second.find( ss_name );
	if ( it2 == it->second.end() ) {
		throw utility::excn::EXCN_Msg_Exception( "SS type " + ss_name + " in layer " + layer_name + " was not found when trying to determine valid residue types." );
	}
	debug_assert( it2 != it->second.end() );

	return it2->second;
}

/// @brief gets the residues allowed for a given secondary structure in a given layer
utility::vector1< std::string > const &
LayerDesignOperation::layer_nc_residues(
	std::string const & layer_name,
	std::string const & ss_name ) const
{
	auto it = layer_nc_residues_.find( layer_name );
	if ( it == layer_nc_residues_.end() ) {
		throw utility::excn::EXCN_Msg_Exception( "Layer " + layer_name + " was not found when trying to determine valid residue types." );
	}
	debug_assert( it != layer_nc_residues_.end() );

	auto it2 = it->second.find( ss_name );
	if ( it2 == it->second.end() ) {
		throw utility::excn::EXCN_Msg_Exception( "SS type " + ss_name + " in layer " + layer_name + " was not found when trying to determine valid ncaa residue types." );
	}
	debug_assert( it2 != it->second.end() );

	return it2->second;
}

/// @brief copies residues allowed from src to dest
void
LayerDesignOperation::copy_layer_residues(
	std::string const & src_layer,
	std::string const & dest_layer )
{
	LayerResidues::const_iterator src = layer_residues_.find( src_layer );
	if ( src == layer_residues_.end() ) {
		throw utility::excn::EXCN_Msg_Exception( "LayerDesign: Layer named " + src_layer + " does not exist when trying to copy residue information to layer " + dest_layer );
	}

	auto lr = layer_residues_.find( dest_layer );
	if ( lr == layer_residues_.end() ) {
		layer_residues_.insert( std::make_pair( dest_layer, src->second ) );
	} else {
		lr->second = src->second;
	}
}

void
LayerDesignOperation::set_layer_residues(
	std::string const & layer_name,
	std::string const & ss_name,
	std::string const & residues )
{
	set_layer_residues( layer_name, ss_name, residues, layer_residues_ );
}

void
LayerDesignOperation::set_nc_layer_residues(
	std::string const & layer_name,
	std::string const & ss_name,
	std::string const & residues )
{
	layer_nc_residues_[ layer_name ][ ss_name ].clear(); //Clear the list of NCAAS
	parse_ncaa_list( residues, layer_nc_residues_[ layer_name ][ ss_name ] ); //Add these NCAAs.
}

void
LayerDesignOperation::set_nc_layer_residues(
	std::string const & layer_name,
	std::string const & ss_name,
	utility::vector1< std::string > const & residues )
{
	layer_nc_residues_[ layer_name ][ ss_name ] = residues;
}

void
LayerDesignOperation::set_layer_residues(
	std::string const & layer_name,
	std::string const & ss_name,
	std::string const & residues,
	LayerResidues & layer_residues )
{
	auto it = layer_residues.find( layer_name );
	if ( it == layer_residues.end() ) {
		it = layer_residues_.insert( std::make_pair( layer_name, LayerDefinitions() ) ).first;
	}
	debug_assert( it != layer_residues_.end() );

	auto it2 = it->second.find( ss_name );
	if ( it2 == it->second.end() ) {
		it2 = it->second.insert( std::make_pair( ss_name, residues ) ).first;
		debug_assert( it2->second == residues );
	} else {
		it2->second = residues;
	}
}

/// @brief apply
void
LayerDesignOperation::apply( Pose const & input_pose, PackerTask & task ) const
{
	using core::pack::task::PackerTask_;
	typedef std::map< std::string, utility::vector1<bool> > LayerSpecification;

	Pose pose;

	// find the designable residues for the different task layers
	LayerSpecification layer_specification;
	BOOST_FOREACH ( const TaskLayers::value_type& task_pair, task_layers_ ) {
		if ( TR.visible() ) TR << "Residues  for task layer " << task_pair.first << ": " <<std::endl;
		PackerTask_ layer_task(input_pose);
		task_pair.second->apply(input_pose, layer_task);
		utility::vector1< bool > designable_residues( layer_task.designing_residues() );
		utility::vector1< bool > packable_residues( layer_task.repacking_residues() );
		utility::vector1< bool > fixed_residues;
		runtime_assert( designable_residues.size() == packable_residues.size() );
		for ( Size i = 1; i <= designable_residues.size(); i++ ) {
			bool fixed( false );
			if ( designable_residues[i] ) {
				if ( TR.visible() ) TR << "\t- residue " << i << " is designable" << std::endl;
			} else if ( ! designable_residues[i] && packable_residues[i] ) {
				if ( TR.visible() ) TR << "\t- residue " << i << " is repackable" << std::endl;
			} else {
				if ( TR.visible() ) TR << "\t- residue " << i << " is fixed" << std::endl;
				fixed = true;
			}
			fixed_residues.push_back( fixed );
		}
		auto type_it = layer_specification_.find( task_pair.first );
		runtime_assert( type_it != layer_specification_.end() );
		LayerSpecificationType const type = (*type_it).second;
		if ( type == DESIGNABLE ) {
			layer_specification[ task_pair.first ] = designable_residues;
		} else if ( type == PACKABLE ) {
			layer_specification[ task_pair.first ] = packable_residues;
		} else if ( type == FIXED ) {
			layer_specification[ task_pair.first ] = fixed_residues;
		}
	}

	// symmetry check
	if ( core::pose::symmetry::is_symmetric( input_pose ) ) {
		if ( use_symmetry() ) {
			TR << "Symmetry detected; will be used in defining layers.  (To disable this, set \"use_symmetry=false\" in RosettaScripts.)" << std::endl;
			pose=input_pose;
		} else {
			TR << "Symmetry detected, extracting asymmetric unit." << std::endl;
			core::pose::symmetry::extract_asymmetric_unit( input_pose, pose , false );
		}
	} else {
		pose = input_pose;
	}

	std::string secstruct;
	if ( blueprint_ ) {
		secstruct = blueprint_->secstruct();
	} else {
		core::scoring::dssp::Dssp dssp( pose );
		dssp.dssp_reduced();
		secstruct = dssp.get_dssp_secstruct();
	}

	// we need to add a SS identifier for the ligand if there is one
	utility::vector1<Size> ligands = protocols::flxbb::find_ligands( pose );
	bool has_ligand = false;
	std::replace( secstruct.begin(), secstruct.end(), ' ', 'L'); // replace all ' ' to 'L'
	if ( TR.visible() ) TR << "secstruct is:" << secstruct << std::endl;
	for ( Size i=1; i <= ligands.size(); i++ ) {
		has_ligand = true;
	}


	srbl_->compute( pose, secstruct );

	// make a pymol script for visualizing the layers
	if ( make_pymol_script_ && !utility::file::file_exists( "layers.py" ) ) {
		if ( TR.visible() ) TR << "writing pymol script with the layer specification and saving it as layers.py" << std::endl;
		write_pymol_script(input_pose, srbl_, layer_specification, has_ligand, "layers.py");
	}

	// find the position of residues of helix capping and intial residue of helix
	bool flag( false );
	utility::vector1< bool > helix_capping( pose.size(), false );
	utility::vector1< bool > initial_helix( pose.size(), false );
	for ( Size i=1; i<=pose.size(); ++i ) {
		// if this residue is not a protein residue, we shouldn't process it further
		if ( ! pose.residue( i ).is_protein() ) continue;
		char ss( secstruct[i-1] );
		if ( ss == 'H' && flag == false && i != 1 ) {
			initial_helix[ i ] = true;
			helix_capping[ i-1 ] = true;
			flag = true;
		}
		if ( ss != 'H' && flag == true ) {
			flag = false;
		}
	}

	// terminal residues set to be all aa
	//utility::vector1<bool> restrict_to_aa( 20, true );
	//task.nonconst_residue_task( 1 ).restrict_absent_canonical_aas( restrict_to_aa );
	//task.nonconst_residue_task( pose.size() ).restrict_absent_canonical_aas( restrict_to_aa );

	TR << "---------------------------------------" << std::endl;
	for ( Size i=1; i<=pose.size(); ++i ) {
		// if the residue is not a protein, we should continue on to the next one, but make the non-protein residue repackable only
		// Removed by VKM -- the user can decide whether to apply design to non-protein residues or not.
		if ( ! pose.residue( i ).is_protein() ) { /*task.nonconst_residue_task( i ).restrict_to_repacking();*/ continue; }

		std::string srbl_layer = srbl_->layer(i);

		// treat peptide edges diferently because dssp always assigns them as loop
		if ( pose.residue( i ).is_upper_terminus() ) {
			srbl_layer = "Cterm";
		}
		if ( pose.residue( i ).is_lower_terminus() ) {
			srbl_layer = "Nterm";
		}

		// check if the residue is specified in any of the task defined layer and
		// append that layer into the active layers
		utility::vector1< std::string > active_layers;
		BOOST_FOREACH ( const LayerSpecification::value_type& layer_pair, layer_specification ) {
			if ( (layer_pair.second)[i] == true ) {
				active_layers.push_back( layer_pair.first );
			}
		}

		char ss( secstruct[i-1] );
		if ( TR.visible() ) {
			TR << "Residue " << i << std::endl;
			TR << "    ss=" << ss << " ";
		}

		if ( srbl_->use_sidechain_neighbors() ) {
			if ( TR.visible() ) TR << "Neighbors=" << ObjexxFCL::format::F( 5, 2, srbl_->rsd_sasa( i ) ) << std::endl;
		} else {
			if ( TR.visible() ) TR << "    Sasa=" << ObjexxFCL::format::F( 6, 2, srbl_->rsd_sasa( i ) ) << std::endl;
		}

		if ( TR.visible() ) TR << "    basic layer = " << srbl_layer << std::endl;

		// If there are no active layers and the working layer is designable
		// append the working layer
		if ( design_layer_.find(srbl_layer)->second || !active_layers.empty() ) { // srbl_layer is set to be designed
			if ( active_layers.empty() ) {
				active_layers.push_back( srbl_layer );
			}
		} else { // srbl_layer is not set to be designed, the sidechain will be reapcked or left untouched
			if ( repack_non_designed_residues_ ) {
				task.nonconst_residue_task( i ).restrict_to_repacking();
				if ( TR.visible() ) TR << "    restricting aminoacid to repacking" << std::endl;
			} else {
				task.nonconst_residue_task( i ).prevent_repacking();
				if ( TR.visible() ) TR << "    preventing aminoacid from  repacking" << std::endl;
			}
		}
		if ( TR.visible() ) {
			TR << "    Active layers: ";
			BOOST_FOREACH ( std::string s, active_layers ) {
				TR << s << "    " ;
			}
			if ( active_layers.empty() ) {
				TR << "none" << std::endl;
			} else {
				TR << std::endl;
			}
			TR << "---------------------------------------" << std::endl;
		}

		// skip the residue if this position is defined as PIKAA, NATRO or NATAA in the resfile
		// V. Mulligan 28 Jan 2015: NOTE that this violates one of the rules for TaskOperations: commutativity.  The behaviour changes depending
		// on whether a ReadResfile TaskOperation is applied before or after the LayerDesign TaskOperation.  I'm changing this to make it the
		// non-default behaviour, but keeping it as an option.
		if ( ignore_pikaa_natro() ) {
			const std::string resfile_cmd =  task.residue_task( i ).command_string();
			if ( resfile_cmd.find( "PIKAA" ) != std::string::npos || resfile_cmd.find( "NATRO" ) != std::string::npos ) {
				if ( TR.visible() && verbose_ ) {
					TR << " ,Resfile info is used. (" << resfile_cmd  << ")"<< std::endl;
				}
				continue;
			}
		}

		// operations precedence: 1) OMIT, 2) NO_DESIGN, 3) DESIGN
		// skip the residue if this position is included in a layer with operation specified as "OMIT"
		bool omit( false );
		bool design( true );
		BOOST_FOREACH ( std::string& layer, active_layers ) {
			auto type_it = layer_operation_.find( layer );
			runtime_assert( type_it != layer_operation_.end() );
			LayerOperationType const operation = (*type_it).second;
			if ( operation == OMIT ) {
				omit = true;
			} else if ( operation == NO_DESIGN ) {
				design = false;
			}
		}
		if ( omit ) {
			if ( TR.visible() ) TR << "Omitting residue " << i << std::endl;
			continue;
		}
		if ( !design ) {
			if ( TR.visible() ) TR << "Not designing residue " << i << " as specified by no_design in the XML tag." << std::endl;
			if ( repack_non_designed_residues_ ) {
				task.nonconst_residue_task( i ).restrict_to_repacking();
			} else {
				task.nonconst_residue_task( i ).prevent_repacking();
			}
			continue;
		}
		BOOST_FOREACH ( std::string& layer, active_layers ) {
			std::string sstype("");
			if ( layer == "Nterm" || layer == "Cterm" ) {
				sstype="all";
			} else if ( helix_capping[ i ] == true && add_helix_capping_ ) {
				sstype="HelixCapping";
			} else if ( initial_helix[i] ) {
				sstype="HelixStart";
			} else if ( ss == 'E' ) {
				sstype="Strand";
			} else if ( ss == 'L' ) {
				sstype="Loop";
			} else if ( ss == 'H' ) {
				sstype="Helix";
			}
			if ( sstype!="" ) {
				//Restrict allowed canonical residues:
				task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( layer, srbl_layer, sstype) );
				//Add allowed noncanonical residues:
				auto cur_nc_layer = layer_nc_residues_.find(layer);
				if ( cur_nc_layer!=layer_nc_residues_.end() ) {
					auto cur_nc_defs = cur_nc_layer->second.find(sstype);
					if ( cur_nc_defs!=cur_nc_layer->second.end() ) {
						if ( !cur_nc_defs->second.empty() ) { //If there are NCAAs defined for this layer:
							for ( auto it=cur_nc_defs->second.begin(); it < cur_nc_defs->second.end(); ++it ) { //Loop through all NCAAs defined for this layer.
								task.nonconst_residue_task(i).allow_noncanonical_aa( *it ); //Add the current NCAA to the list allowed for this residue.
							}
						}
					}
				}
			}
			if ( TR.visible() ) TR << i << " done " << std::endl << std::flush;
		} // for each layer
	} // for( i )
} // apply

void
LayerDesignOperation::parse_tag( TagCOP tag , DataMap & datamap )
{
	//typedef std::pair< std::string, bool > DesignLayerPair;

	use_original_ = tag->getOption< bool >( "use_original_non_designed_layer", 1 );

	String design_layers = tag->getOption< String >( "layer", "core_boundary_surface_Nterm_Cterm" );
	if ( design_layers == "all" ) {
		design_layers = "core_boundary_surface_Nterm_Cterm";
	}
	if ( design_layers == "other" || design_layers == "user" ) {
		design_layers = "";
	}
	utility::vector1< String > layers( utility::string_split( design_layers, '_' ) );
	set_design_layers( layers );

	repack_non_designed_residues_ = tag->getOption< bool >("repack_non_design", 1);

	//Should residues that have been set to design or not design in a resfile with PIKAA or NATRO,
	//where the resfile is applied BEFORE the LayerDesign operation, be ignored?  Default false to
	//preserve commutativity.
	set_ignore_pikaa_natro( tag->getOption<bool>( "ignore_pikaa_natro", false ) );

	srbl_->set_design_layer( true, true, true );

	if ( tag->hasOption( "use_sidechain_neighbors" ) ) {
		set_use_sidechain_neighbors( tag->getOption< bool >( "use_sidechain_neighbors" ) );
	}

	if ( tag->hasOption( "sc_neighbor_dist_midpoint" ) ) {
		set_sc_neighbor_dist_midpoint( tag->getOption<core::Real>("sc_neighbor_dist_midpoint") );
	}

	if ( tag->hasOption( "sc_neighbor_denominator" ) ) {
		set_sc_neighbor_denominator( tag->getOption<core::Real>("sc_neighbor_denominator") );
	}

	if ( tag->hasOption( "sc_neighbor_angle_shift_factor") ) {
		set_sc_neighbor_angle_shift_factor( tag->getOption<core::Real>("sc_neighbor_angle_shift_factor") );
	}

	if ( tag->hasOption( "sc_neighbor_angle_exponent") ) {
		set_sc_neighbor_angle_exponent( tag->getOption<core::Real>("sc_neighbor_angle_exponent") );
	}

	if ( tag->hasOption( "sc_neighbor_dist_exponent") ) {
		set_sc_neighbor_dist_exponent( tag->getOption<core::Real>("sc_neighbor_dist_exponent") );
	}

	if ( tag->hasOption( "pore_radius" ) ) {
		srbl_->pore_radius( tag->getOption< Real >( "pore_radius" ) );
	}

	if ( tag->hasOption( "core" ) ) {
		srbl_->sasa_core( tag->getOption< Real >( "core" ) );
	}
	if ( tag->hasOption( "surface" ) ) {
		srbl_->sasa_surface( tag->getOption< Real >( "surface" ) );
	}

	if ( tag->hasOption( "core_E" ) ) {
		srbl_->sasa_core( tag->getOption< Real >( "core_E" ), "E" );
	}
	if ( tag->hasOption( "core_L" ) ) {
		srbl_->sasa_core( tag->getOption< Real >( "core_L" ), "L" );
	}
	if ( tag->hasOption( "core_H" ) ) {
		srbl_->sasa_core( tag->getOption< Real >( "core_H" ), "H" );
	}

	if ( tag->hasOption( "surface_E" ) ) {
		srbl_->sasa_surface( tag->getOption< Real >( "surface_E" ), "E" );
	}
	if ( tag->hasOption( "surface_L" ) ) {
		srbl_->sasa_surface( tag->getOption< Real >( "surface_L" ), "L" );
	}
	if ( tag->hasOption( "surface_H" ) ) {
		srbl_->sasa_surface( tag->getOption< Real >( "surface_H" ), "H" );
	}

	if ( tag->hasOption( "make_rasmol_script" ) ) {
		srbl_->make_rasmol_format_file( tag->getOption< bool >( "make_rasmol_script" ) );
	}

	if ( tag->hasOption( "blueprint" ) ) {
		blueprint_ = protocols::jd2::parser::BluePrintOP( new protocols::jd2::parser::BluePrint( tag->getOption< std::string >("blueprint") ) );
	}

	set_use_symmetry( tag->getOption<bool>("use_symmetry", true) );

	set_verbose( tag->getOption< bool >( "verbose", false ) );
	set_restrict_restypes( tag->getOption< bool >( "restrict_restypes", true ) );
	make_pymol_script( tag->getOption< bool >("make_pymol_script", false) );

	for ( auto layer_tag : tag->getTags() ) {
		parse_layer_tag( layer_tag, datamap );
	}

	if ( TR.visible() ) {
		TR << "Layers to be designed:";
		for ( auto l_p : design_layer_ )  {
			if ( l_p.second ) {
				TR << "\t" << l_p.first;
			}
		}
		TR << std::endl;
		// print the layer definitions
		for ( auto layer : layer_residues_ ) {
			TR << "Layer " << layer.first << std::endl;
			for ( auto layer_def : layer.second ) {
				TR << "\t" << ObjexxFCL::format::LJ(15,layer_def.first) << "aa = " << layer_def.second << "\t ncaa = " << print_string_vector( layer_nc_residues_[ layer.first ][ layer_def.first ] ) << std::endl;
			}
			TR << std::endl;
		}
	}
}

std::string
layer_design_ss_layer_naming_func( std::string const & element_name )
{
	return "layer_design_ss_layer_" + element_name + "_type";
}

std::string
layer_design_ss_layer_group() {
	return "layer_design_ss_layer";
}

std::string
layer_design_ss_layer_or_taskop_group() {
	return "layer_design_ss_layer_or_taskop";
}

// AMW: Add a common_simple_type for an underscore separated string list?
// TO DO: Fix this after we've changed the format of LayerDesign
void LayerDesignOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using core::pack::task::operation::TaskOperationFactory;
	AttributeList attributes;

	attributes
		+ optional_name_attribute()
		+ XMLSchemaAttribute::attribute_w_default(  "use_original_non_designed_layer", xsct_rosetta_bool, "Restrict to repacking the non design layers",  "true"  )
		+ XMLSchemaAttribute::attribute_w_default(  "layer", xs_string, "Layer to be designed, other ex. core_surface means only design core and surface layer, other refers to the additional layers defined with packertasks",  "core_boundary_surface_Nterm_Cterm"  )
		+ XMLSchemaAttribute::attribute_w_default(  "repack_non_design", xsct_rosetta_bool, "If true, side chains will be repacked, left untouched if otherwise.",  "true"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ignore_pikaa_natro", xsct_rosetta_bool, "If true, and if a resfile is read before applying this TaskOperation, ignore any residues that have been set in the resfile with the PIKAA, NATRO, or NATAA commands.",  "false"  )
		+ XMLSchemaAttribute( "use_sidechain_neighbors", xsct_rosetta_bool , "If true, assign a residue's layers based on counting the number CA atoms from other residues within a cone in front of the residue's ca-cb vector. Because this option is no longer SASA based, the layer assignments will always be identical regardless of the protein sequence; i.e. layers could be assigned based on a polyalanine backbone and it would make no difference. This option changes the defaults for core and surface to neighbors less than 2 --surface-- and neighbors greater than 5.2 --core. HOWEVER, these defaults will be overwritten if core and surface are manually specified in declaring the taskoperation! So make sure you do not specify new core and surface settings appropriate for SASA when you are actually counting neighboring residues." )
		+ XMLSchemaAttribute( "sc_neighbor_dist_midpoint", xsct_real , "These values fine-tune the behavior of the sidechain neighbors residue-counting logic. Typically, a user need not change these from default values. For details on these, see the LayerSelector ResidueSelector's documentation." )
		+ XMLSchemaAttribute( "sc_neighbor_denominator", xsct_real , "These values fine-tune the behavior of the sidechain neighbors residue-counting logic. Typically, a user need not change these from default values. For details on these, see the LayerSelector ResidueSelector's documentation." )
		+ XMLSchemaAttribute( "sc_neighbor_angle_shift_factor", xsct_real , "These values fine-tune the behavior of the sidechain neighbors residue-counting logic. Typically, a user need not change these from default values. For details on these, see the LayerSelector ResidueSelector's documentation." )
		+ XMLSchemaAttribute( "sc_neighbor_angle_exponent", xsct_real , "These values fine-tune the behavior of the sidechain neighbors residue-counting logic. Typically, a user need not change these from default values. For details on these, see the LayerSelector ResidueSelector's documentation." )
		+ XMLSchemaAttribute( "sc_neighbor_dist_exponent", xsct_real , "These values fine-tune the behavior of the sidechain neighbors residue-counting logic. Typically, a user need not change these from default values. For details on these, see the LayerSelector ResidueSelector's documentation." )
		+ XMLSchemaAttribute( "pore_radius", xsct_real , "pore radius for calculating accessible surface area" )
		+ XMLSchemaAttribute( "core", xsct_real , "residues of which asa is less than core are defined as core" )
		+ XMLSchemaAttribute( "surface", xsct_real , "residues of which asa is greater than surface are defined as surface" )
		+ XMLSchemaAttribute( "core_E", xsct_real , "XRW TO DO" )
		+ XMLSchemaAttribute( "core_L", xsct_real , "XRW TO DO" )
		+ XMLSchemaAttribute( "core_H", xsct_real , "XRW TO DO" )
		+ XMLSchemaAttribute( "surface_E", xsct_real , "XRW TO DO" )
		+ XMLSchemaAttribute( "surface_L", xsct_real , "XRW TO DO" )
		+ XMLSchemaAttribute( "surface_H", xsct_real , "XRW TO DO" )
		+ XMLSchemaAttribute( "make_rasmol_script", xsct_rosetta_bool , "if true, write a rasmol script coloring the residues by the three basic layers, core, boundary and surface." )
		+ XMLSchemaAttribute( "blueprint", xs_string , "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default(  "use_symmetry", xsct_rosetta_bool, "More recently, the use_symmetry option has been added to permit LayerDesign to be symmetry-aware. If use_symmetry is set to true --the default, layers are defined for symmetric poses using the full, symmetric pose.",  "true"  )
		+ XMLSchemaAttribute::attribute_w_default(  "verbose", xsct_rosetta_bool, "print to tracer",  "false"  )
		+ XMLSchemaAttribute::attribute_w_default(  "restrict_restypes", xsct_rosetta_bool, "XRW TO DO",  "true"  )
		+ XMLSchemaAttribute::attribute_w_default(  "make_pymol_script", xsct_rosetta_bool, "if true, write a pymol script coloring the residues by the three basic layer and the aditional taskoperation defined layers..",  "false"  );

	XMLSchemaRepeatableCTNodeOP root_node(new XMLSchemaRepeatableCTNode);

	root_node->set_element_w_attributes(
		keyname(), attributes,
		"Note: The LayerDesign TaskOperation will likely be deprecated at some point in the future in favour of the"
		" LayerSelector ResidueSelector. It is strongly recommended that users start to switch over to the LayerSelector"
		" ResidueSelector, which permits greater flexibility in selecting residues. Design residues with selected amino"
		" acids depending on the enviroment--accessible surface area. The layer of each residue is assigned to one of the"
		" three basic layers--core, boundary or surface-- depending on the accessible surface area of mainchain + CB, or"
		" depending on the number of neighbours in a cone extending along the CA-CB vector --if the use_sidechain_neighbors"
		" option is used." );

	root_node->set_root_node_naming_func(
		&core::pack::task::operation::complex_type_name_for_task_op);

	root_node->set_kids_naming_func( & layer_design_ss_layer_naming_func );

	// child layer beneath LayerDesign:
	// reads "core", "boundary", "surface", "Nterm", "Cterm"
	// grandchild layer beneath this layer reads
	// SS_TYPES?

	// grandchild layer beneath CombinedTasks:
	// reads "all", "Helix", "Strand", "Loop"

	// grandchild layer beneath TaskLayer:
	// reads SS_TYPES??

	// ------------------- SS layer ----------------------- //
	XMLSchemaRestriction layer_design_operation_behavior;
	layer_design_operation_behavior.name( "layer_design_operation_behavior" );
	layer_design_operation_behavior.base_type( xs_string );
	layer_design_operation_behavior.add_restriction( xsr_enumeration, "design" );
	layer_design_operation_behavior.add_restriction( xsr_enumeration, "no_design" );
	layer_design_operation_behavior.add_restriction( xsr_enumeration, "omit" );
	xsd.add_top_level_element( layer_design_operation_behavior );

	XMLSchemaRestriction layer_design_specification_behavior;
	layer_design_specification_behavior.name( "layer_design_specification_behavior" );
	layer_design_specification_behavior.base_type( xs_string );
	layer_design_specification_behavior.add_restriction( xsr_enumeration, "designable" );
	layer_design_specification_behavior.add_restriction( xsr_enumeration, "repackable" );
	layer_design_specification_behavior.add_restriction( xsr_enumeration, "fixed" );
	xsd.add_top_level_element( layer_design_specification_behavior );


	typedef XMLSchemaAttribute Attr;
	AttributeList ss_layer_attributes;
	ss_layer_attributes
		+ Attr( "copy_layer", xs_string, "XRW TO DO" )
		+ Attr( "aa", xs_string, "XRW TO DO" )
		+ Attr( "ncaa", xs_string, "XRW TO DO" )
		+ Attr( "append", xs_string, "XRW TO DO" )
		+ Attr( "ncaa_append", xs_string, "XRW TO DO" )
		+ Attr( "exclude", xs_string, "XRW TO DO" )
		+ Attr( "ncaa_exclude", xs_string, "XRW TO DO" )
		+ Attr( "operation", "layer_design_operation_behavior", "XRW TO DO" )
		+ Attr( "specification", "layer_design_specification_behavior", "XRW TO DO" );

	// Define the group
	XMLSchemaModelGroup ss_group;
	ss_group.group_name( layer_design_ss_layer_group() );
	XMLSchemaModelGroupOP ss_choice( new XMLSchemaModelGroup );
	ss_choice->type( xsmgt_choice );
	XMLSchemaModelGroup ss_or_taskop_group;
	ss_or_taskop_group.group_name( layer_design_ss_layer_or_taskop_group() );
	XMLSchemaModelGroupOP ss_or_taskop_choice( new XMLSchemaModelGroup );
	ss_or_taskop_choice->type( xsmgt_choice );

	for ( auto ss_type : SS_TYPES ) {
		XMLSchemaComplexTypeGenerator ss_ct_gen;
		ss_ct_gen.element_name( ss_type )
			.complex_type_naming_func( & layer_design_ss_layer_naming_func )
			.description( "XRW TO DO" )
			.add_attributes( ss_layer_attributes )
			.write_complex_type_to_schema( xsd );

		XMLSchemaElementOP ss_element( new XMLSchemaElement );
		ss_element->name( ss_type ).type_name( layer_design_ss_layer_naming_func( ss_type ));
		ss_choice->append_particle( ss_element );
		ss_or_taskop_choice->append_particle( ss_element );
	}
	ss_group.append_particle( ss_choice );

	ss_or_taskop_group.append_particle( ss_or_taskop_choice );
	XMLSchemaModelGroupOP taskop_group_ref( new XMLSchemaModelGroup );
	taskop_group_ref->group_name( TaskOperationFactory::task_operation_xml_schema_group_name() );
	ss_or_taskop_choice->append_particle( taskop_group_ref );

	xsd.add_top_level_element( ss_group );
	xsd.add_top_level_element( ss_or_taskop_group );

	// ---------------------- TaskLayer complex type ------------------- //
	// The first subelement must be a task operation;
	// all subsequent subelements must be one of the SS subelements
	XMLSchemaSimpleSubelementList task_layer_task_op_element;
	task_layer_task_op_element.add_group_subelement( & TaskOperationFactory::task_operation_xml_schema_group_name );
	XMLSchemaSimpleSubelementList task_layer_ss_elements;
	task_layer_ss_elements.add_group_subelement( & layer_design_ss_layer_group );
	XMLSchemaComplexTypeGenerator task_layer_ctgen;
	task_layer_ctgen.element_name( "TaskLayer" )
		.complex_type_naming_func( & layer_design_ss_layer_naming_func )
		.description( "XRW TO DO" )
		.add_ordered_subelement_set_as_required( task_layer_task_op_element )
		.add_ordered_subelement_set_as_repeatable( task_layer_ss_elements )
		.write_complex_type_to_schema( xsd );

	// Everything below follows the repeatable subelement idiom;

	XMLSchemaRepeatableCTNodeOP core_node(new XMLSchemaRepeatableCTNode);
	XMLSchemaRepeatableCTNodeOP boundary_node(new XMLSchemaRepeatableCTNode);
	XMLSchemaRepeatableCTNodeOP surface_node(new XMLSchemaRepeatableCTNode);
	XMLSchemaRepeatableCTNodeOP nterm_node(new XMLSchemaRepeatableCTNode);
	XMLSchemaRepeatableCTNodeOP cterm_node(new XMLSchemaRepeatableCTNode);
	XMLSchemaRepeatableCTNodeOP combined_tasks_node(new XMLSchemaRepeatableCTNode);
	XMLSchemaRepeatableCTNodeOP task_layer_node(new XMLSchemaRepeatableCTNode);

	root_node->add_child(core_node)
		.add_child(boundary_node)
		.add_child(surface_node)
		.add_child(nterm_node)
		.add_child(cterm_node)
		.add_child(combined_tasks_node)
		.add_child(task_layer_node);

	XMLSchemaRepeatableCTNodeOP ss_layer_node( new XMLSchemaRepeatableCTNode );
	ss_layer_node->set_group_subelement( & layer_design_ss_layer_group );

	AttributeList empty_attribs; // this list is empty!
	core_node->set_element_w_attributes( "core", empty_attribs, "XRW TO DO" ).add_child( ss_layer_node );
	boundary_node->set_element_w_attributes( "boundary", empty_attribs, "XRW TO DO" ).add_child( ss_layer_node );
	surface_node->set_element_w_attributes( "surface", empty_attribs, "XRW TO DO" ).add_child( ss_layer_node );
	nterm_node->set_element_w_attributes( "Nterm", empty_attribs, "XRW TO DO" ).add_child( ss_layer_node );
	cterm_node->set_element_w_attributes( "Cterm", empty_attribs, "XRW TO DO" ).add_child( ss_layer_node );

	XMLSchemaRepeatableCTNodeOP ss_layer_or_taskop_node( new XMLSchemaRepeatableCTNode );
	ss_layer_or_taskop_node->set_group_subelement( & layer_design_ss_layer_or_taskop_group );

	AttributeList combined_task_attribs; combined_task_attribs + required_name_attribute();
	combined_tasks_node->set_element_w_attributes( "CombinedTasks", combined_task_attribs, "XRW TO DO");

	//combined_tasks_node->set_kids_naming_func( & layer_design_ss_layer_or_taskop_group );
	combined_tasks_node->add_child( ss_layer_or_taskop_node );

	task_layer_node->set_already_defined_element( "TaskLayer", & layer_design_ss_layer_naming_func );

	root_node->recursively_write_ct_to_schema(xsd);
}

void
LayerDesignOperation::parse_layer_tag( utility::tag::TagCOP layer_tag, DataMap & datamap )
{
	using core::pack::task::operation::TaskOperationFactory;
	std::string layer = layer_tag->getName(); // core, residue, boundary or taskoperation
	if ( layer == "core" || layer =="boundary" || layer == "surface" ||  layer == "Nterm" ||  layer == "Cterm" ) {
		TR << "redefining default layer " << layer << std::endl;
		for ( utility::tag::TagCOP const secstruct_tag : layer_tag->getTags() ) {
			parse_layer_secstruct_tag( secstruct_tag, datamap, layer );
		}
	} else if ( layer == "CombinedTasks" ) {
		std::string const comb_name = layer_tag->getOption< std::string >( "name" );
		layer = comb_name;
		TR << "Making a combined task named " << comb_name << std::endl;

		// 1st pass over subtags: read all the task operations.
		utility::vector1< TaskOperationOP > task_ops;
		for ( auto task_tag : layer_tag->getTags() ) {
			std::string task_op_type = task_tag->getName();
			if ( task_op_type == "all" || task_op_type == "Helix" || task_op_type == "Strand" || task_op_type == "Loop" ) {
				continue;
			} else if ( TaskOperationFactory::get_instance()->has_type( task_op_type ) ) {
				TaskOperationOP task = TaskOperationFactory::get_instance()->newTaskOperation( task_op_type, datamap, task_tag );
				task_ops.push_back( task );
			} else {
				utility::excn::EXCN_Msg_Exception( "Illegal subtag \"" + task_op_type + "\" of CombinedTasks subtag of CombinedTasks. Expected \"all\", \"Helix\", \"Strand\", \"Loop\" or a valid TaskOperation name" );
			}
		}
		CombinedTaskOperationOP comb( new CombinedTaskOperation( task_ops ) );
		add_layer( comb_name, comb, DESIGN, DESIGNABLE );

		// 2nd pass over subtags: read all the other acceptible subtags.
		for ( utility::tag::TagCOP const secstruct_tag : layer_tag->getTags() ) {
			std::string secstruct_type = secstruct_tag->getName();
			if ( secstruct_type == "all" || secstruct_type == "Helix" || secstruct_type == "Strand" || secstruct_type == "Loop" ) {
				parse_layer_secstruct_tag( secstruct_tag, datamap, layer );
			}
		}
	} else if ( layer == "TaskLayer" ) {
		bool first_tag = true;
		std::string const task_op_type = layer;
		std::string task_name;
		for ( utility::tag::TagCOP const task_layer_subtag : layer_tag->getTags() ) {
			std::string const task_op_type = task_layer_subtag->getName();
			if ( first_tag ) {
				first_tag = false;
				if ( TaskOperationFactory::get_instance()->has_type( task_op_type ) ) {
					task_name = task_layer_subtag->getOption< std::string >( "name" );
					TR << "Defining new layer from task type " << layer << " named " << task_name << std::endl;
					TaskOperationOP task = TaskOperationFactory::get_instance()->newTaskOperation( task_op_type, datamap, task_layer_subtag );
					add_layer( task_name, task, DESIGN, DESIGNABLE );
				} else {
					throw utility::excn::EXCN_Msg_Exception( "Invalid TaskOperation name \"" + task_op_type + "\" in TaskLayer subtag of LayerDesign" );
				}
			} else {
				parse_layer_secstruct_tag( task_layer_subtag, datamap, task_name );
			}
		}

	} else {
		utility_exit_with_message( "Invalid layer " + layer + ", valid layers are core, boundary, surface, TaskOperations, CombinedTasks, or TaskLayer" );
	}

	// layer_tag is not accessed from here forward.

	// check if layer ss is defined and if not fill it up
	BOOST_FOREACH ( Layer const & layer, layer_residues_ ) {
		// fill up empty canonical amino acid lists with all amino acids
		std::string all_layers_residues = "";
		LayerDefinitions::const_iterator d = layer_residues_[ layer.first ].find( "all" );
		if ( d == layer_residues_[ layer.first ].end() ) {
			all_layers_residues = "ARNDCEQGHILKMFPSTWYV";
		} else {
			all_layers_residues = d->second;
		}

		// for all ss types except all, fill in residues
		BOOST_FOREACH ( std::string const & ss_def_name, SS_TYPES ) {
			if ( ss_def_name == "all" ) {
				continue;
			}

			auto ld_it = layer_residues_[ layer.first ].find( ss_def_name );
			if ( ld_it == layer_residues_[ layer.first ].end() ) {
				TR << "layer " << layer.first << " has no specification for residues in " << ss_def_name << ".  The layer will be filled with the residues defined for all secondary structure types." << std::endl;
				set_layer_residues( layer.first, ss_def_name, all_layers_residues );
			}
		}

		// fill up empty non-canonical amino acid lists with "all" list if it exists
		utility::vector1< std::string > all_layers_ncaa_residues;
		LayerNCDefinitions::const_iterator nc = layer_nc_residues_[ layer.first ].find( "all" );
		if ( nc != layer_nc_residues_[ layer.first ].end() ) {
			all_layers_ncaa_residues = layer_nc_residues( layer.first, "all" );
		}

		// skip if "all" layer is empty
		if ( ! all_layers_ncaa_residues.size() ) {
			continue;
		}

		BOOST_FOREACH ( std::string const & ss_def_name, SS_TYPES ) {
			auto ldncaa_it = layer_nc_residues_[ layer.first ].find( ss_def_name );
			if ( ( ldncaa_it == layer_nc_residues_[ layer.first ].end() ) || ( !ldncaa_it->second.size() ) ) {
				TR << "layer " << layer.first << " has no specification for noncanonical residues in " << ss_def_name << ".  The layer will be filled with the noncanonical residues defined for all secondary structure types." << std::endl;
				set_nc_layer_residues( layer.first, ss_def_name, all_layers_ncaa_residues );
			}
		}
	}
}

/// @brief returns a string containing sorted, unique residue one-letter codes
std::string unique_chars( std::string const orig )
{
	std::set< char > temp_def_res_set( orig.begin(), orig.end() );
	return std::string( temp_def_res_set.begin(), temp_def_res_set.end() );
}

utility::vector1< std::string > unique_strs(
	utility::vector1< std::string > const & orig,
	utility::vector1< std::string > const & orig2 )
{
	std::set< std::string > strset( orig.begin(), orig.end() );
	for ( auto const & s : orig2 ) {
		strset.insert( s );
	}
	return utility::vector1< std::string >( strset.begin(), strset.end() );
}

void
LayerDesignOperation::parse_layer_secstruct_tag(
	TagCOP secstruct_tag,
	DataMap &,
	std::string const & layer_name )
{
	std::string secstruct = secstruct_tag->getName(); // Strand, Helix, Loop, HelixCapping

	if ( secstruct_tag->hasOption( "copy_layer" ) ) {
		std::string const layer_to_copy = secstruct_tag->getOption< std::string >( "copy_layer" );
		TR << "Copying definitions from layer " << layer_to_copy << " to layer " << layer_name << std::endl;
		copy_layer_residues( layer_to_copy, layer_name );
	}
	if ( secstruct_tag->hasOption( "aa" ) ) {
		std::string const aas = unique_chars( secstruct_tag->getOption< std::string >( "aa" ) );
		TR << "Assigning residues " << aas << " to layer " << layer_name << " for ss " << secstruct << std::endl;
		set_layer_residues( layer_name, secstruct, unique_chars( aas ) );
		if ( secstruct == "all" ) {
			LayerDefinitions & ld = layer_residues_[ layer_name ];
			for ( LayerDefinitions::const_iterator l = ld.begin(); l != ld.end(); ++l ) {
				if ( l->first != "all" ) {
					set_layer_residues( layer_name, l->first, aas );
				}
			}
		}
	}
	if ( secstruct_tag->hasOption( "ncaa" ) ) {
		std::string const aas = secstruct_tag->getOption< std::string >( "ncaa" );
		set_nc_layer_residues( layer_name, secstruct, aas );
		TR << "Assigning noncanonical residues " << aas << " to layer " << layer_name << " for ss " << secstruct << std::endl;
		if ( secstruct == "all" ) {
			LayerNCDefinitions & ld = layer_nc_residues_[ layer_name ];
			for ( LayerNCDefinitions::const_iterator l = ld.begin(); l != ld.end(); ++l ) {
				if ( l->first != "all" ) {
					set_nc_layer_residues( layer_name, l->first, aas );
				}
			}
		}
	}
	if ( secstruct_tag->hasOption( "append" ) ) {
		std::string const aas = secstruct_tag->getOption< std::string >( "append" );
		LayerDefinitions::const_iterator prevaas = layer_residues_[ layer_name ].find( secstruct );
		if ( prevaas == layer_residues_[ layer_name ].end() ) {
			set_layer_residues( layer_name, secstruct, unique_chars( aas ) );
		} else {
			set_layer_residues( layer_name, secstruct, unique_chars( aas + prevaas->second ) );
		}
		if ( secstruct == "all" ) {
			TR << "Appending residues " << aas << " to layer " << layer_name << std::endl;
			auto lrs =  layer_residues_.find( layer_name );
			for ( auto & ld : lrs->second ) {
				std::string const & cur_aas = ld.second;
				// prevent appending to all twice
				if ( ld.first != "all" ) {
					set_layer_residues( layer_name, ld.first, unique_chars( cur_aas + aas ) );
				}
			}
		}
	}
	if ( secstruct_tag->hasOption( "ncaa_append" ) ) {
		LayerNCDefinitions & defs = layer_nc_residues_[ layer_name ];
		std::string const aas = secstruct_tag->getOption< std::string >( "ncaa_append" );
		utility::vector1< std::string > aas_vec;
		parse_ncaa_list( aas, aas_vec );
		defs[ secstruct ] = unique_strs( aas_vec, defs[ secstruct ] );

		TR << "Appending noncanonical residues " << aas << " to layer " << layer_name << std::endl;
		for ( auto & def : defs ) {
			def.second = unique_strs( aas_vec, def.second );
		}
	}
	if ( secstruct_tag->hasOption( "exclude" ) ) {
		std::string const aas = secstruct_tag->getOption< std::string >( "exclude" );
		LayerDefinitions::const_iterator cur_aas = layer_residues_[ layer_name ].find( secstruct );
		if ( cur_aas != layer_residues_[ layer_name ].end() ) {
			std::set< char > temp_def_res_set( cur_aas->second.begin(), cur_aas->second.end() );
			BOOST_FOREACH ( char aa, aas ) {
				temp_def_res_set.erase( aa );
			}
			set_layer_residues( layer_name, secstruct, std::string( temp_def_res_set.begin(), temp_def_res_set.end() ) );
		} else {
			if ( secstruct != "all" ) {
				throw utility::excn::EXCN_Msg_Exception("LayerDesignOperation: Trying to exclude aas from a secstruct that doesn't exist. Layer = " + layer_name + " secstruct = " + secstruct );
			}
		}

		TR << "Excluding residues " << aas << " from layer " << layer_name << " secstruct " << secstruct << std::endl;
		if ( secstruct == "all" ) {
			auto lrs =  layer_residues_.find( layer_name );
			for ( auto ld = lrs->second.begin(); ld != lrs->second.end(); ld++ ) {
				if ( ld->first == "all" ) {
					continue;
				}
				std::set< char > temp_def_res_set( ld->second.begin(), ld->second.end() );
				BOOST_FOREACH ( char aa, aas ) {
					temp_def_res_set.erase( aa );
				}
				set_layer_residues( lrs->first, ld->first, std::string( temp_def_res_set.begin(), temp_def_res_set.end() ) );
			}
		}
	}
	if ( secstruct_tag->hasOption( "ncaa_exclude" ) ) {
		LayerNCDefinitions & defs = layer_nc_residues_[ layer_name ];
		std::string const aas = secstruct_tag->getOption< std::string >( "ncaa_exclude" );
		TR << "Excluding noncanonical residues " << aas << " from layer " << layer_name << std::endl;
		utility::vector1< std::string > res_to_exclude;
		parse_ncaa_list( aas, res_to_exclude );
		exclude_ncaas( res_to_exclude, defs[ secstruct ] );
		for ( auto & def : defs ) {
			exclude_ncaas( res_to_exclude, def.second );
		}
	}

	// Set layer behaivour
	if ( secstruct_tag->hasOption( "operation" ) ) {
		std::string const operation_str = secstruct_tag->getOption< std::string >( "operation" );
		TR << "Apply operation  " << operation_str << " to layer " << layer_name << std::endl;
		if ( operation_str == "design" ) {
			set_layer_operation( layer_name, DESIGN );
		} else if ( operation_str == "no_design" ) {
			set_layer_operation( layer_name, NO_DESIGN );
		} else if ( operation_str == "omit" ) {
			set_layer_operation( layer_name, OMIT );
		} else {
			utility_exit_with_message( "Invalid operation " + operation_str + " specified for layer " + layer_name + ". Valid options are \"design\", \"no_design\", and \"omit\"." );
		}
	}

	if ( secstruct_tag->hasOption( "specification" ) ) {
		std::string const specification = secstruct_tag->getOption< std::string >( "specification" );
		TR << "Apply specification  " << specification << " to layer " << layer_name << std::endl;
		if ( specification == "designable" ) {
			set_layer_specification( layer_name, DESIGNABLE );
		} else if ( specification == "repackable" ) {
			set_layer_specification( layer_name, PACKABLE );
		} else if ( specification == "fixed" ) {
			set_layer_specification( layer_name, FIXED );
		} else {
			utility_exit_with_message( "Invalid specification " + specification + " for layer " + layer_name + ". Valid options are \"designable\", \"packable\", and \"fixed\"." );
		}
	}
}

} // flxbb
} // protocols
