// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/flxbb/LayerDesignOperation.cc
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
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <basic/Tracer.hh>
#include <protocols/toolbox/SelectResiduesByLayer.hh>
#include <core/pose/symmetry/util.hh>

// Utility Headers
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <string>

using basic::T;
using basic::Error;
using basic::Warning;


#include <utility/vector0.hh>
#include <ObjexxFCL/format.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <boost/assign/list_inserter.hpp> 
#include <boost/assign/list_of.hpp> 
#include <boost/foreach.hpp>



using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;

static basic::Tracer TR("protocols.flxbb.LayerDesignOperation");

namespace protocols {
namespace flxbb {

core::pack::task::operation::TaskOperationOP
LayerDesignOperationCreator::create_task_operation() const
{
	return new LayerDesignOperation;
}


/// @brief default constructor
LayerDesignOperation::LayerDesignOperation():
	TaskOperation(),
	add_helix_capping_( true ),
	use_original_( true ),
	verbose_( false ),
	restrict_restypes_( true ),
	srbl_( new toolbox::SelectResiduesByLayer )
{
	set_default_layer_residues();
}

/// @brief value constructor
LayerDesignOperation::LayerDesignOperation( bool dsgn_core, bool dsgn_boundary, bool dsgn_surface ):
	TaskOperation(),
	add_helix_capping_( true ),
	use_original_( true ),
	verbose_( false ),
	restrict_restypes_( true ),
	srbl_( new toolbox::SelectResiduesByLayer )
{
	design_layer( dsgn_core, dsgn_boundary, dsgn_surface );
	set_default_layer_residues();
}

/// @brief destructor
LayerDesignOperation::~LayerDesignOperation() {
}

/// @brief clone
core::pack::task::operation::TaskOperationOP
LayerDesignOperation::clone() const {
	return new LayerDesignOperation( *this );
}

/// @brief layer to be designed
void
LayerDesignOperation::design_layer( bool const dsgn_core, bool const dsgn_boundary, bool const dsgn_surface )
{
	srbl_->set_design_layer( dsgn_core, dsgn_boundary, dsgn_surface );
}

/// @brief accessible surface for evaluating residues are in surface or not
void
LayerDesignOperation::sasa_surface( Real const r, String const ss )
{
	srbl_->sasa_surface( r, ss );
}

/// @brief accessible surface for evaluating residues are in core or not
void
LayerDesignOperation::sasa_core( Real const r, String const ss )
{
	srbl_->sasa_core( r, ss );
}

/// @brief set pore radius for colculating asa
void
LayerDesignOperation::pore_radius( Real ps )
{
	srbl_->pore_radius( ps );
}

void
LayerDesignOperation::set_default_layer_residues() {
	TR << "initializing the layer with the default residues" << std::endl;
	boost::assign::insert(layer_residues_)
					("core", boost::assign::map_list_of
							("Loop", 						"AFILPVWY") 
							("Strand",	 				"FILVWY") 
							("Helix", 					"AFILVWY") 
							("HelixStart", 			"AFILVWYP") 
							("HelixCapping", 		"DNST") 
					)
					("boundary", boost::assign::map_list_of
							("Loop", 						"ADEFGIKLNPQRSTVWY") 
							("Strand",	 				"DEFIKLNQRSTVWY") 
							("Helix", 					"ADEIKLNQRSTVWY") 
							("HelixStart", 			"ADEIKLNQRSTVWYP") 
							("HelixCapping", 		"DNST") 
					)
					("surface", boost::assign::map_list_of
							("Loop", 						"DEGHKNPQRST") 
							("Strand",	 				"DEHKNQRST") 
							("Helix", 					"DEHKNQRST") 
							("HelixStart",			"DEHKNQRSTP") 
							("HelixCapping", 		"DNST") 
					)
					("interface", boost::assign::map_list_of
							("Loop", 						"DEGHKNPQRST") 
							("Strand",	 				"DEHKNQRST") 
							("Helix", 					"DEHKNQRST") 
							("HelixStart",			"DEHKNQRSTP") 
							("HelixCapping", 		"DNST") 
					);
	
	boost::assign::insert(design_layer_)
					("core", false)
					("boundary", false)
					("surface", false);
}

utility::vector1<bool> 
LayerDesignOperation::get_restrictions( std::string const & layer, std::string const & default_layer, std::string const & ss_type) const {	
	// if the layer doesn't specify the required ss used the default layer one.
  std::string used_layer = ( layer_residues_.find(layer)->second.count(ss_type) != 0 ) ? layer : default_layer;
	utility::vector1<bool>  restrict_to_aa( chemical::num_canonical_aas, false );
	BOOST_FOREACH(char restype, layer_residues_.find(layer)->second.find(ss_type)->second){
		restrict_to_aa[chemical::aa_from_oneletter_code( restype )] = true;
	}
	return restrict_to_aa;
}

/// @brief apply
void
LayerDesignOperation::apply( Pose const & input_pose, PackerTask & task ) const
{
	using core::pack::task::PackerTask_;
	typedef std::map< std::string, utility::vector1<bool> > LayerSpecification;

	Pose pose;
	// symmetry check
	if(core::pose::symmetry::is_symmetric( input_pose ) ) {
		TR << "Symmetry detected, extracting asymmetric unit." << std::endl;
		core::pose::symmetry::extract_asymmetric_unit( input_pose, pose /*, false*/ );
	 } else {
		pose = input_pose;
  }
	// calc SelectResiduesByLayer
	srbl_->compute( pose, "" );

	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	// find the position of residues of helix capping and intial residue of helix
	bool flag( false );
	utility::vector1< bool > helix_capping( pose.total_residue(), false );
	utility::vector1< bool > initial_helix( pose.total_residue(), false );
	for( Size i=1; i<=pose.total_residue(); ++i ) {
		char ss( dssp.get_dssp_secstruct( i ) );
		if( ss == 'H' && flag == false && i != 1 ) {
			initial_helix[ i ] = true;
			helix_capping[ i-1 ] = true;
			flag = true;
		}
		if( ss != 'H' && flag == true ) {
			flag = false;
		}
	}
	// find the designable residues for the different task layers
	LayerSpecification layer_specification;
	BOOST_FOREACH(const TaskLayers::value_type& task_pair, task_layers_) { 
		TR << "Residues for task layer " << task_pair.first << ": ";
		PackerTask_ layer_task(pose);
  	task_pair.second->apply(pose, layer_task);
		utility::vector1< bool > designable_residues( layer_task.designing_residues() );
		layer_specification[ task_pair.first ] = designable_residues;
	}

	// terminal residues set to be all aa
	utility::vector1<bool> restrict_to_aa( 20, true );
	task.nonconst_residue_task( 1 ).restrict_absent_canonical_aas( restrict_to_aa );
	task.nonconst_residue_task( pose.total_residue() ).restrict_absent_canonical_aas( restrict_to_aa );

	for( Size i=2; i<=pose.total_residue()-1; ++i ) {

		const std::string srbl_layer = srbl_->layer(i);
		// check if the residue is specified in any of the task defined layer
		utility::vector1< std::string > active_layers;
		BOOST_FOREACH(const LayerSpecification::value_type& layer_pair, layer_specification) {
			if( (layer_pair.second)[i] == true) 
				active_layers.push_back( layer_pair.first );
		}
		// If there are no active layers and the working layer is designable
		// append the working layer
		if( active_layers.empty() && design_layer_.find(srbl_layer)->second ) {
			active_layers.push_back(srbl_layer);
		} else {
			if(use_original_)
  			task.nonconst_residue_task( i ).restrict_to_repacking();

		}

		char ss( dssp.get_dssp_secstruct( i ) );
		if( verbose_ ) {
			TR << " Resnum=" << i << " ,SS=" << ss << " "
				 << " Sasa=" << ObjexxFCL::fmt::F( 6, 2, srbl_->rsd_sasa( i ) );
		}
		TR << std::endl;

		// skip the residue if this position is defined as PIKAA by resfile
		if( task.residue_task( i ).command_string().find( "PIKAA" ) != std::string::npos ){
			if( verbose_ ) {
				TR << " ,Resfile info is used." << std::endl;
			}
			continue;
		}
		
		BOOST_FOREACH(std::string& layer, active_layers) {
  	  if( helix_capping[ i ] == true && add_helix_capping_ ) {
  			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( layer, srbl_layer, "HelixCapping") );
  			if( verbose_ ) 
  				TR << " ,Helix Capping " << std::endl;
  
  		} else if( initial_helix[i] ) {
  			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( layer, srbl_layer, "HelixStart") );
  
  		} else if( ss == 'E') {
  			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( layer, srbl_layer, "Strand") );
  
  		} else if( ss == 'L') {
  			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( layer, srbl_layer, "Loop") );
  
  		} else if( ss == 'H') {
  			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( layer, srbl_layer, "Helix") );
  
  		} 
		}
	} // for( i )
} // apply

void
LayerDesignOperation::parse_tag( TagPtr tag )
{
	using core::pack::task::operation::TaskOperationFactory;
	use_original_ = tag->getOption< bool >( "use_original_non_designed_layer", 1 );

	String design_layers = tag->getOption< String >( "layer", "core_boundary_surface" );
	utility::vector1< String > layers( utility::string_split( design_layers, '_' ) );
	BOOST_FOREACH(std::string &  layer, layers) {
		design_layer_[ layer ] = true;
	}
	

	srbl_->set_design_layer( true, true, true);

	if( tag->hasOption( "pore_radius" ) ) {
		srbl_->pore_radius( tag->getOption< Real >( "pore_radius" ) );
	}

	if( tag->hasOption( "core" ) ) {
		srbl_->sasa_core( tag->getOption< Real >( "core" ) );
	}
	if( tag->hasOption( "surface" ) ) {
		srbl_->sasa_surface( tag->getOption< Real >( "surface" ) );
	}

	if( tag->hasOption( "core_E" ) ) {
		srbl_->sasa_core( tag->getOption< Real >( "core_E" ), "E" );
	}
	if( tag->hasOption( "core_L" ) ) {
		srbl_->sasa_core( tag->getOption< Real >( "core_L" ), "L" );
	}
	if( tag->hasOption( "core_H" ) ) {
		srbl_->sasa_core( tag->getOption< Real >( "core_H" ), "H" );
	}

	if( tag->hasOption( "surface_E" ) ) {
		srbl_->sasa_surface( tag->getOption< Real >( "surface_E" ), "E" );
	}
	if( tag->hasOption( "surface_L" ) ) {
		srbl_->sasa_surface( tag->getOption< Real >( "surface_L" ), "L" );
	}
	if( tag->hasOption( "surface_H" ) ) {
		srbl_->sasa_surface( tag->getOption< Real >( "surface_H" ), "H" );
	}

	if( tag->hasOption( "make_rasmol_script" ) ) {
		srbl_->make_rasmol_format_file( tag->getOption< bool >( "make_rasmol_script" ) );
	}

	set_verbose( tag->getOption< bool >( "verbose", false ) );
	set_restrict_restypes( tag->getOption< bool >( "restrict_restypes", true ) );

	BOOST_FOREACH( utility::tag::TagPtr const layer_tag, tag->getTags() ){
		std::string layer = layer_tag->getName(); // core, residue, boundary
		if( TaskOperationFactory::get_instance()->has_type(layer) ) {
			std::string task_op_type = layer;
			std::string task_name = layer_tag->getOption< std::string >("name");
			TR << "Defining new layer from task type "<< layer << " named " << task_name << std::endl;
			layer = task_name;
			TaskOperationOP task = TaskOperationFactory::get_instance()->newTaskOperation(task_op_type, layer_tag);
			// store the task to use it in apply to find the designable residues for the layer
			// and add the extra layer to layer_residues_.
			task_layers_[ task_name ] = task;
			layer_residues_[ layer ] = std::map< std::string, std::string >();
		} 

		BOOST_FOREACH( utility::tag::TagPtr const secstruct_tag, layer_tag->getTags() ){
			std::string secstruct = secstruct_tag->getName(); // Strand, Helix, Loop, HelixCapping
			std::string aas = secstruct_tag->getOption< std::string >("aa");
			TR << "Setting layer residues for " << layer << " "<< secstruct <<" to " << aas << std::endl;
			layer_residues_[ layer ][ secstruct ] = aas;
		}
	}
}

} // flxbb
} // protocols
