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
/// @modified Javier Castellanos (javiercv@uw.edu )

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
#include <protocols/toolbox/SelectResiduesByLayer.hh>
#include <core/pose/symmetry/util.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/jd2/parser/BluePrint.hh>

// Utility Headers
#include <protocols/flxbb/utility.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <string>

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

static basic::Tracer TR("protocols.flxbb.LayerDesignOperation");

namespace protocols {
namespace flxbb {

core::pack::task::operation::TaskOperationOP
LayerDesignOperationCreator::create_task_operation() const
{
	return new LayerDesignOperation;
}

CombinedTaskOperation::CombinedTaskOperation( VecTaskOP  ops ):
	task_operations_( ops )
{ }

void
CombinedTaskOperation::apply(core::pose::Pose const & pose, PackerTask & task) const {
	using core::pack::task::operation::TaskOperationOP;
	BOOST_FOREACH( TaskOperationOP t_op, task_operations_ )
		t_op->apply( pose, task );
}


/// @brief default constructor
LayerDesignOperation::LayerDesignOperation():
	TaskOperation(),
	add_helix_capping_( true ),
	use_original_( true ),
	repack_non_designed_residues_( true ),
	verbose_( false ),
	restrict_restypes_( true ),
	make_pymol_script_( false ),
	srbl_( new toolbox::SelectResiduesByLayer ),
	blueprint_( NULL )
{
	set_default_layer_residues();
}

/// @brief value constructor
LayerDesignOperation::LayerDesignOperation( bool dsgn_core, bool dsgn_boundary, bool dsgn_surface ):
	TaskOperation(),
	add_helix_capping_( true ),
	use_original_( true ),
	repack_non_designed_residues_( true ),
	verbose_( false ),
	restrict_restypes_( true ),
	make_pymol_script_( false ),
	srbl_( new toolbox::SelectResiduesByLayer ),
	blueprint_( NULL )
{
	design_layer( dsgn_core, dsgn_boundary, dsgn_surface );
	set_default_layer_residues();
}

LayerDesignOperation::LayerDesignOperation( LayerDesignOperation const & rval ):
	 add_helix_capping_( rval.add_helix_capping_ ),
	 use_original_( rval.use_original_ ),
	 repack_non_designed_residues_( rval.repack_non_designed_residues_ ),
	 verbose_( rval.verbose_ ),
	 restrict_restypes_( rval.restrict_restypes_ ),
	 make_pymol_script_( rval.make_pymol_script_ ),
	 layer_residues_( rval.layer_residues_ ),
	 design_layer_( rval.design_layer_ ),
	 layer_specification_( rval.layer_specification_ ),
	 layer_operation_( rval.layer_operation_ ),
	 task_layers_( rval.task_layers_ ),
	 srbl_( rval.srbl_ ),
	 blueprint_( rval.blueprint_ )
{
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

std::string
LayerDesignOperation::pos2select( utility::vector1< core::Size > const & pos ) const
{
	std::string str;
	if( pos.empty() )
		return str;
	else
		str = boost::lexical_cast< std::string	>( pos[1] );
		for(Size i = 2; i <= pos.size(); i++)
			str += "+" + boost::lexical_cast< std::string >( pos[i] );
		return str;
}

void
LayerDesignOperation::write_pymol_script( core::pose::Pose const & pose, toolbox::SelectResiduesByLayerOP srbl, std::map< std::string, utility::vector1<bool> > const & layer_specification, bool has_ligand, std::string const & filename ) const
{
	using utility::io::ozstream;
	typedef utility::vector1<Size> VecSize;
	typedef std::map< std::string, utility::vector1<bool> > LayerSpecification;
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
	for(LayerSpecification::const_iterator it = layer_specification.begin(); it != layer_specification.end(); it++){
		utility::vector1< Size > pos;
		for(Size i = 1; i <= pose.total_residue(); i++)
			if( it->second[ i ])
				pos.push_back( i );
		pymol << "cmd.select('"<< it->first  <<"', 'resi  " << pos2select( pos )<< "')" << std::endl;
		pymol << "cmd.color('" << colors[ layer % 6  ] << "','" << it->first << "')" << std::endl;
		layer += 1;
	}
	if(has_ligand) {
		pymol << "cmd.select('ligand', 'organic')" << std::endl;
		pymol << "cmd.show('sticks','ligand')" << std::endl;
		pymol << "cmd.color('gray','ligand')" << std::endl;
	}
}


void
LayerDesignOperation::set_default_layer_residues() {
	TR << "Initializing the layers with the default residues" << std::endl;
	boost::assign::insert(layer_residues_)
		      ("core", boost::assign::map_list_of
					 		("all",							"AFILPVWYDNST")
							("Loop", 						"AFILPVWY")
							("Strand",	 				"FILVWY")
							("Helix", 					"AFILVWY")
							("HelixStart", 			"AFILVWYP")
							("HelixCapping", 		"DNST")
					)
					("boundary", boost::assign::map_list_of
							("all", 						"ADEFGIKLNPQRSTVWY")
							("Loop", 						"ADEFGIKLNPQRSTVWY")
							("Strand",	 				"DEFIKLNQRSTVWY")
							("Helix", 					"ADEIKLNQRSTVWY")
							("HelixStart", 			"ADEIKLNQRSTVWYP")
							("HelixCapping", 		"DNST")
					)
					("surface", boost::assign::map_list_of
							("all", 						"DEGHKNPQRST")
							("Loop", 						"DEGHKNPQRST")
							("Strand",	 				"DEHKNQRST")
							("Helix", 					"DEHKNQRST")
							("HelixStart",			"DEHKNQRSTP")
							("HelixCapping", 		"DNST")
					)
					("Nterm", boost::assign::map_list_of
					 ("all", "ACDEFGHIKLMNPQRSTVWY")
					)
					("Cterm", boost::assign::map_list_of
					 ("all", "ACDEFGHIKLMNPQRSTVWY")
					);

	boost::assign::insert(design_layer_)
		("core", false )
		("boundary", false )
		("surface", false )
		("Nterm", false )
		("Cterm", false );

	boost::assign::insert(layer_specification_)
		("core", DESIGNABLE)
		("boundary", DESIGNABLE)
		("surface", DESIGNABLE)
		("Nterm", DESIGNABLE)
		("Cterm", DESIGNABLE);

	boost::assign::insert(layer_operation_)
		("core", DESIGN)
		("boundary", DESIGN)
		("surface", DESIGN)
		("Nterm", DESIGN)
		("Cterm", DESIGN);
}

utility::vector1<bool>
LayerDesignOperation::get_restrictions( std::string const & layer, std::string const & default_layer, std::string const & ss_type) const {
	// if the layer doesn't specify the required ss used the default layer one.
  std::string used_layer = ( layer_residues_.find(layer)->second.count(ss_type) != 0 ) ? layer : default_layer;
	utility::vector1<bool>  restrict_to_aa( chemical::num_canonical_aas, false );
	BOOST_FOREACH(char restype, layer_residues_.find(layer)->second.find(ss_type)->second) {
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

	// find the designable residues for the different task layers
	LayerSpecification layer_specification;
	BOOST_FOREACH(const TaskLayers::value_type& task_pair, task_layers_) {
		TR << "Residues  for task layer " << task_pair.first << ": " <<std::endl;
		PackerTask_ layer_task(input_pose);
  	task_pair.second->apply(input_pose, layer_task);
		utility::vector1< bool > designable_residues( layer_task.designing_residues() );
		utility::vector1< bool > packable_residues( layer_task.repacking_residues() );
		utility::vector1< bool > fixed_residues;
		runtime_assert( designable_residues.size() == packable_residues.size() );
		for(Size i = 1; i <= designable_residues.size(); i++) {
			bool fixed( false );
			if( designable_residues[i] ) {
				TR << "\t- residue " << i << " is designable" << std::endl;
			} else if ( ! designable_residues[i] && packable_residues[i] ) {
				TR << "\t- residue " << i << " is repackable" << std::endl;
			} else {
				TR << "\t- residue " << i << " is fixed" << std::endl;
				fixed = true;
			}
			fixed_residues.push_back( fixed );
		}
		std::map< std::string, LayerSpecificationType >::const_iterator type_it = layer_specification_.find( task_pair.first );
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
	if(core::pose::symmetry::is_symmetric( input_pose ) ) {
		TR << "Symmetry detected, extracting asymmetric unit." << std::endl;
		core::pose::symmetry::extract_asymmetric_unit( input_pose, pose , false );
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
  TR << "secstruct is:" << secstruct << std::endl;
	for( Size i=1; i <= ligands.size(); i++) {
		TR << "adding an L to the secstruct string due to unknown AA" << std::endl;
		secstruct += 'L';
		has_ligand = true;
	}
	srbl_->compute( pose, secstruct );

	// make a pymol script for visualizing the layers
	if( make_pymol_script_ && !utility::file::file_exists( "layers.py" ) ) {
		TR << "writing pymol script with the layer specification and saving it as layers.py" << std::endl;
		write_pymol_script(input_pose, srbl_, layer_specification, has_ligand, "layers.py");
	}

	// find the position of residues of helix capping and intial residue of helix
	bool flag( false );
	utility::vector1< bool > helix_capping( pose.total_residue(), false );
	utility::vector1< bool > initial_helix( pose.total_residue(), false );
	for( Size i=1; i<=pose.total_residue(); ++i ) {
		// if this residue is not a protein residue, we shouldn't process it further
		if ( ! pose.residue( i ).is_protein() ) continue;
		char ss( secstruct[i] );
		if( ss == 'H' && flag == false && i != 1 ) {
			initial_helix[ i ] = true;
			helix_capping[ i-1 ] = true;
			flag = true;
		}
		if( ss != 'H' && flag == true ) {
			flag = false;
		}
	}

	// terminal residues set to be all aa
	//utility::vector1<bool> restrict_to_aa( 20, true );
	//task.nonconst_residue_task( 1 ).restrict_absent_canonical_aas( restrict_to_aa );
	//task.nonconst_residue_task( pose.total_residue() ).restrict_absent_canonical_aas( restrict_to_aa );

	TR << "---------------------------------------" << std::endl;
	for( Size i=1; i<=pose.total_residue(); ++i ) {
		// if the residue is not a protein, we should continue on to the next one, but make the non-protein residue repackable only
		if ( ! pose.residue( i ).is_protein() ){ task.nonconst_residue_task( i ).restrict_to_repacking(); continue; }


		std::string srbl_layer = srbl_->layer(i);

		// treat peptide edges diferently because dssp always assigns them as loop
		if( pose.residue( i ).is_upper_terminus() )
			srbl_layer = "Cterm";
		if( pose.residue( i ).is_lower_terminus() )
			srbl_layer = "Nterm";

		// check if the residue is specified in any of the task defined layer and
		// append that layer into the active layers
		utility::vector1< std::string > active_layers;
		BOOST_FOREACH(const LayerSpecification::value_type& layer_pair, layer_specification) {
			if( (layer_pair.second)[i] == true)
				active_layers.push_back( layer_pair.first );
		}

		char ss( secstruct[i] );
		TR << "Residue " << i << std::endl;
		TR << "    ss=" << ss << " "
			 << "    Sasa=" << ObjexxFCL::format::F( 6, 2, srbl_->rsd_sasa( i ) ) << std::endl;
		TR << "    basic layer = " << srbl_layer << std::endl;

		// If there are no active layers and the working layer is designable
		// append the working layer
		if( design_layer_.find(srbl_layer)->second || !active_layers.empty() ) { // srbl_layer is set to be designed
			if( active_layers.empty() )
				active_layers.push_back( srbl_layer );
		} else { // srbl_layer is not set to be designed, the sidechain will be reapcked or left untouched
			if( repack_non_designed_residues_ ) {
  			task.nonconst_residue_task( i ).restrict_to_repacking();
				TR << "    restricting aminoacid to repacking" << std::endl;
			}
			else {
  			task.nonconst_residue_task( i ).prevent_repacking();
				TR << "    prenventing aminoacid from  repacking" << std::endl;
			}
		}
		TR << "    Active layers: ";
		BOOST_FOREACH(std::string s, active_layers){
			TR << s << "    " ;
		}
		if(active_layers.empty())
			TR << "none" << std::endl;
		else
			TR << std::endl;
		TR << "---------------------------------------" << std::endl;


		// skip the residue if this position is defined as PIKAA, NATRO or NATAA in the resfile
		const std::string resfile_cmd =  task.residue_task( i ).command_string();
		if( resfile_cmd.find( "PIKAA" ) != std::string::npos ||	resfile_cmd.find( "NATRO" ) != std::string::npos ){
			if( verbose_ ) {
				TR << " ,Resfile info is used. (" << resfile_cmd  << ")"<< std::endl;
			}
			continue;
		}

		// operations precedence: 1) OMIT, 2) NO_DESIGN, 3) DESIGN
		// skip the residue if this position is included in a layer with operation specified as "OMIT"
 		bool omit( false );
		bool design( true );
		BOOST_FOREACH(std::string& layer, active_layers) {
			std::map< std::string, LayerOperationType >::const_iterator type_it = layer_operation_.find( layer );
			runtime_assert( type_it != layer_operation_.end() );
			LayerOperationType const operation = (*type_it).second;
			if ( operation == OMIT ) {
				omit = true;
			} else if ( operation == NO_DESIGN ) {
				design = false;
			}
		}
		if ( omit ) {
			TR << "Omitting residue " << i << std::endl;
			continue;
		}
		if ( !design ) {
			TR << "Not designing residue " << i << " as specified by no_design in the XML tag." << std::endl;
			if ( repack_non_designed_residues_ ) {
				task.nonconst_residue_task( i ).restrict_to_repacking();
			} else {
				task.nonconst_residue_task( i ).prevent_repacking();
			}
			continue;
		}
		BOOST_FOREACH(std::string& layer, active_layers) {
			if( layer == "Nterm" || layer == "Cterm") {
  			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( layer, srbl_layer, "all") );
			} else if( helix_capping[ i ] == true && add_helix_capping_ ) {
  			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( layer, srbl_layer, "HelixCapping") );

  		} else if( initial_helix[i] ) {
  			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( layer, srbl_layer, "HelixStart") );
  		} else if( ss == 'E') {
  			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( layer, srbl_layer, "Strand") );

  		} else if( ss == 'L') {
  			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( layer, srbl_layer, "Loop") );

  		} else if( ss == 'H') {
  			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( get_restrictions( layer, srbl_layer, "Helix") );

  		}
			TR << i << " done " << std::endl << std::flush;
		} // for each layer
	} // for( i )
} // apply

void
LayerDesignOperation::parse_tag( TagCOP tag , DataMap & datamap )
{
	using core::pack::task::operation::TaskOperationFactory;
	typedef std::pair< std::string, bool > DesignLayerPair;

	use_original_ = tag->getOption< bool >( "use_original_non_designed_layer", 1 );

	String design_layers = tag->getOption< String >( "layer", "core_boundary_surface_Nterm_Cterm" );
	if (design_layers == "all")
		design_layers = "core_boundary_surface_Nterm_Cterm";
	if (design_layers == "other" || design_layers == "user")
		design_layers = "";
	utility::vector1< String > layers( utility::string_split( design_layers, '_' ) );
	BOOST_FOREACH(std::string &  layer, layers) {
		design_layer_[ layer ] = true;
	}

	repack_non_designed_residues_ = tag->getOption< bool >("repack_non_design", 1);

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

	if( tag->hasOption( "blueprint" ) ) {
		blueprint_ = new protocols::jd2::parser::BluePrint( tag->getOption< std::string >("blueprint") );
	}

	set_verbose( tag->getOption< bool >( "verbose", false ) );
	set_restrict_restypes( tag->getOption< bool >( "restrict_restypes", true ) );
	make_pymol_script( tag->getOption< bool >("make_pymol_script", false) );

	BOOST_FOREACH( utility::tag::TagCOP const layer_tag, tag->getTags() ){

		std::string layer = layer_tag->getName(); // core, residue, boundary or taskoperation
		if( layer == "core" || layer =="boundary" || layer == "surface" ||  layer == "Nterm" ||  layer == "Cterm" )
        {
            TR << "redefining default layer " << layer << std::endl;
		} else if(layer == "CombinedTasks" ) {
			std::string comb_name = layer_tag->getOption< std::string >("name");
			TR << "Making a combined task named "<< comb_name << std::endl;
			utility::vector1< TaskOperationOP > task_ops;
			BOOST_FOREACH( utility::tag::TagCOP const task_tag, layer_tag->getTags() ) {
				std::string task_op_type = task_tag->getName();
				TaskOperationOP task = TaskOperationFactory::get_instance()->newTaskOperation(task_op_type, datamap, task_tag);
				task_ops.push_back( task );
			}
			CombinedTaskOperationOP comb = new CombinedTaskOperation( task_ops );
			task_layers_[ comb_name ] = comb;
			design_layer_[ comb_name ] = true;
			layer_residues_[ comb_name ] = std::map< std::string, std::string >();
		} else if( TaskOperationFactory::get_instance()->has_type(layer) ) {
			std::string task_op_type = layer;
			std::string task_name = layer_tag->getOption< std::string >("name");
			TR << "Defining new layer from task type "<< layer << " named " << task_name << std::endl;
			layer = task_name;
			TaskOperationOP task = TaskOperationFactory::get_instance()->newTaskOperation(task_op_type, datamap, layer_tag);
			// store the task to use it in apply to find the designable residues for the layer
			// and add the extra layer to layer_residues_.
			task_layers_[ task_name ] = task;
			design_layer_[ task_name ] = true;
			layer_residues_[ task_name ] = std::map< std::string, std::string >();
            
		} else {
			utility_exit_with_message( "Invalid layer " + layer + ", valid layers are core, boundary, surface, TaskOperations or CombinedTasks" );
		}

		BOOST_FOREACH( utility::tag::TagCOP const secstruct_tag, layer_tag->getTags() ) {
			std::string secstruct = secstruct_tag->getName(); // Strand, Helix, Loop, HelixCapping
			if( secstruct == "all" &&  secstruct_tag->hasOption("copy_layer") ) {
				const std::string layer_to_copy = secstruct_tag->getOption< std::string >("copy_layer");
				TR << "Copying definitions from layer " << layer_to_copy << " to layer " << layer << std::endl;
				layer_residues_[ layer ] = layer_residues_[ layer_to_copy ];
			}
			if( secstruct == "all" &&  secstruct_tag->hasOption("aa") ) {
				const std::string aas = secstruct_tag->getOption< std::string >("aa");
				LayerResidues::iterator lrs =  layer_residues_.find( layer );
				TR << "Assigning residues " << aas << " to layer " << lrs->first << std::endl;
				for(LayerDefinitions::iterator ld = lrs->second.begin(); ld != lrs->second.end(); ld++) {
					std::set<char> temp_def_res_set;
					temp_def_res_set.insert( aas.begin(), aas.end() );
					layer_residues_[ lrs->first ][ ld->first ] = std::string(temp_def_res_set.begin(), temp_def_res_set.end() );
				}
			}
			if( secstruct == "all" &&  secstruct_tag->hasOption("append") ) {
				const std::string aas = secstruct_tag->getOption< std::string >("append");
				LayerResidues::iterator lrs =  layer_residues_.find( layer );
				TR << "Appending residues " << aas << " to layer " << lrs->first << std::endl;
				for(LayerDefinitions::iterator ld = lrs->second.begin(); ld != lrs->second.end(); ld++) {
					std::set<char> temp_def_res_set( ld->second.begin(), ld->second.end());
					temp_def_res_set.insert( aas.begin(), aas.end() );
					layer_residues_[ lrs->first ][ ld->first ] = std::string(temp_def_res_set.begin(), temp_def_res_set.end() );
				}
			}

			if( secstruct == "all" &&  secstruct_tag->hasOption("exclude") ) {
				const std::string aas = secstruct_tag->getOption< std::string >("exclude");
				LayerResidues::iterator lrs =  layer_residues_.find( layer );
				TR << "Excluding residues " << aas << " from layer " << lrs->first << std::endl;
				for(LayerDefinitions::iterator ld = lrs->second.begin(); ld != lrs->second.end(); ld++) {
					std::set<char> temp_def_res_set( ld->second.begin(), ld->second.end());
					BOOST_FOREACH(char aa, aas)
						temp_def_res_set.erase(aa);
					layer_residues_[ lrs->first ][ ld->first ] = std::string(temp_def_res_set.begin(), temp_def_res_set.end() );
				}
			}

			if( secstruct_tag->hasOption("aa") ) {
			std::string aas = secstruct_tag->getOption< std::string >("aa");
			TR << "Setting layer residues for " << layer << " "<< secstruct <<" to " << aas << std::endl;
			layer_residues_[ layer ][ secstruct ] = aas;
			}

			if( secstruct_tag->hasOption("append") ) {
			std::string aas = secstruct_tag->getOption< std::string >("append");
			TR << "Appending residues "<< aas << " to layer " << layer << " "<< secstruct << std::endl;
			const std::string layer_res = layer_residues_[ layer ][ secstruct ];
			std::set<char> temp_def_res_set( layer_res.begin(), layer_res.end());
			temp_def_res_set.insert( aas.begin(), aas.end() );
			layer_residues_[ layer ][ secstruct ] = std::string(temp_def_res_set.begin(), temp_def_res_set.end() );
			}

			if( secstruct_tag->hasOption("exclude") ) {
			std::string aas = secstruct_tag->getOption< std::string >("exclude");
			TR << "Excluding residues "<< aas << " to layer " << layer << " "<< secstruct << std::endl;
			const std::string layer_res = layer_residues_[ layer ][ secstruct ];
			std::set<char> temp_def_res_set( layer_res.begin(), layer_res.end());
			BOOST_FOREACH(char aa, aas)
				temp_def_res_set.erase(aa);
			layer_residues_[ layer ][ secstruct ] = std::string(temp_def_res_set.begin(), temp_def_res_set.end() );
			}

		}
        
        std::string const operation_str = layer_tag->getOption< std::string >("operation", "design" );
        TR << "operation=" << operation_str << std::endl;
        if ( operation_str == "design" ) {
            layer_operation_[ layer ] = DESIGN;
        } else if ( operation_str == "no_design" ) {
            layer_operation_[ layer ] = NO_DESIGN;
        } else if ( operation_str == "omit" ) {
            layer_operation_[ layer ] = OMIT;
        } else {
            utility_exit_with_message( "Invalid operation " + operation_str + " specified for layer " + layer + ". Valid options are \"design\", \"no_design\", and \"omit\"." );
        }

        
        std::string const specification = layer_tag->getOption< std::string >("specification", "designable");
        if ( specification == "designable" ) {
            layer_specification_[ layer ] = DESIGNABLE;
        } else if ( specification == "repackable" ) {
            layer_specification_[ layer ] = PACKABLE;
        } else if ( specification == "fixed" ) {
            layer_specification_[ layer ] = FIXED;
        } else {
            utility_exit_with_message( "Invalid specification " + specification + " for layer " + layer + ". Valid options are \"designable\", \"packable\", and \"fixed\"." );
        }
	}

	// fill empty the empty layers of the task layers with the residues at the 'all' layer
	std::set< std::string > ss_def_names; // pick the ss names from the core layer
	LayerResidues::const_iterator default_layer_it = layer_residues_.find("core");
	BOOST_FOREACH(LayerDefinition const & layer_def, default_layer_it->second )
		ss_def_names.insert( layer_def.first );
	ss_def_names.erase("all");

	// check if layer ss is defined and if not fill it up
	BOOST_FOREACH(TaskLayer  const & task_layer, task_layers_) {
		const std::string all_layers_residues = (layer_residues_.count(task_layer.first)) ? layer_residues_[ task_layer.first ][ "all" ] : "ARNDCEQGHILKMFPSTWYV";
		BOOST_FOREACH( std::string const & ss_def_name, ss_def_names){
			LayerDefinitions::iterator ld_it = layer_residues_[ task_layer.first ].find( ss_def_name );
			if( ld_it == layer_residues_[ task_layer.first ].end() ) {
				TR << "layer " << task_layer.first << " has no specification for residues in " << ss_def_name << ", the layer will be filled with the residues defined for all secondary structure types." << std::endl;
				layer_residues_[ task_layer.first ][ ss_def_name ] = all_layers_residues;
			}
		}
	}


	TR << "Layers to be designed:";
	BOOST_FOREACH( DesignLayerPair l_p, design_layer_)  {
		if( l_p.second )
			TR << "\t" << l_p.first;
	}
	TR << std::endl;
	// print the layer definitions
	TR << std::endl;
	BOOST_FOREACH( Layer const &  layer, layer_residues_) {
		TR << "Layer " << layer.first << std::endl;
		BOOST_FOREACH( LayerDefinition const & layer_def, layer.second ) {
			TR << "\t" << ObjexxFCL::format::LJ(15,layer_def.first) << "aa = " << layer_def.second << std::endl;
		}
		TR << std::endl;
	}
}

} // flxbb
} // protocols
