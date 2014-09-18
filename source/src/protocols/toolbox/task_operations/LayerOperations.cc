// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/flxbb/LayerOperations.cc
/// @brief operation defines either core, surface or periphery. It solely determines which residue to pack/design, 
///	   does not define identities.
/// @author Eva-Maria Strauch (evas01@uw.edu)


// Unit Headers
#include <protocols/toolbox/task_operations/LayerOperations.hh>
#include <protocols/toolbox/task_operations/LayerOperationsCreator.hh>

// Project Headers
// AUTO-REMOVED #include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rigid/RigidBodyMover.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask_.hh>
#include <core/scoring/dssp/Dssp.hh>
// AUTO-REMOVED #include <core/scoring/sasa.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
#include <basic/Tracer.hh>
//#include <protocols/toolbox/LayerOperations.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/conformation/Conformation.hh>

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

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;


static thread_local basic::Tracer TR( "protocols.toolbox.LayerOperations" );

namespace protocols {
namespace toolbox {
namespace task_operations  {
    
    
/// @brief destructor
LayerOperations::~LayerOperations() {}

/// @brief clone
core::pack::task::operation::TaskOperationOP
LayerOperations::clone() const {
    return new LayerOperations( *this );
}

    
core::pack::task::operation::TaskOperationOP
LayerOperations::create_task_operation() const
{
	return new LayerOperations;
}

/// @brief default constructor
LayerOperations::LayerOperations():
pick_core_( false ),
pick_boundary_( false ),
pick_surface_( false ),
pore_radius_( 2.0 ),
make_rasmol_format_file_( false )
{
    initialize( 20.0, 40.0 );
}


/// @brief value constructor
LayerOperations::LayerOperations( bool const pick_core, bool const pick_boundary, bool const pick_surface ):
pick_core_( pick_core ),
pick_boundary_( pick_boundary ),
pick_surface_( pick_surface ),
pore_radius_( 2.0 ),
make_rasmol_format_file_( false )
{
    initialize( 20.0, 40.0 );
}

/// @brief value constructor
LayerOperations::LayerOperations( String const pick ):
pick_core_( false ),
pick_boundary_( false ),
pick_surface_( false ),
pore_radius_( 2.0 ),
make_rasmol_format_file_( false )
{
    initialize( 20.0, 40.0 );
    utility::vector1< String > layers( utility::string_split( pick, '_' ) );
    for( utility::vector1< String >::const_iterator iter = layers.begin(); iter != layers.end() ; ++iter) {
        String layer(*iter);
        if ( layer == "core" ) {
            pick_core_ = true;
        } else if ( layer == "surface" ) {
            pick_surface_ = true;
        } else if ( layer == "boundary" ) {
            pick_boundary_ = true;
        } else {
            TR << "Error!, wrong specification of layer_mode " << layer << std::endl;
            runtime_assert( false );
        }
    } // utility::vector1
}
/// @brief copy constructor
LayerOperations::LayerOperations( LayerOperations const & rval ):
pick_core_( rval.pick_core_ ),
pick_boundary_( rval.pick_boundary_ ),
pick_surface_( rval.pick_surface_ ),
pore_radius_( rval.pore_radius_ ),
burial_( rval.burial_ ),
surface_( rval.surface_ ),
selected_core_residues_( rval.selected_core_residues_ ),
selected_boundary_residues_( rval.selected_boundary_residues_ ),
selected_surface_residues_( rval.selected_surface_residues_ ),
make_rasmol_format_file_( rval.make_rasmol_format_file_ ),
rsd_sasa_( rval.rsd_sasa_ ),
rsd_layer_( rval.rsd_layer_ )
{}

/// @brief accessible surface for evaluating residues are in surface or not
void
LayerOperations::sasa_surface( Real const r, String const ss )
{
    if ( ss == "" ) {
        surface_[ 'H' ] = r;
        surface_[ 'L' ] = r;
        surface_[ 'E' ] = r;
    } else {
        runtime_assert( ss.length() == 1 );
        runtime_assert( ss.at( 0 ) == 'L' || ss.at( 0 ) == 'E' || ss.at( 0 ) == 'H' );
        surface_[ ss.at( 0 ) ] = r;
    }
}

/// @brief accessible surface for evaluating residues are in core or not
void
LayerOperations::sasa_core( Real const r, String const ss )
{
    if ( ss == "" ) {
        burial_[ 'H' ] = r;
        burial_[ 'L' ] = r;
        burial_[ 'E' ] = r;
    } else {
        runtime_assert( ss.length() == 1 );
        runtime_assert( ss.at( 0 ) == 'L' || ss.at( 0 ) == 'E' || ss.at( 0 ) == 'H' );
        burial_[ ss.at( 0 ) ] = r;
    }
}

/// @brief
void
LayerOperations::initialize( Real burial, Real surface )
{
    surface_.insert( std::map< char, Real >::value_type( 'H', surface ) );
    surface_.insert( std::map< char, Real >::value_type( 'L', surface ) );
    surface_.insert( std::map< char, Real >::value_type( 'E', surface ) );
    burial_.insert( std::map< char, Real >::value_type( 'H', burial ) );
    burial_.insert( std::map< char, Real >::value_type( 'L', burial ) );
    burial_.insert( std::map< char, Real >::value_type( 'E', burial ) );
}


/// @brief amino acid types excluded for selection
//void
//LayerOperations::exclude_aatypes_for_selection( utility::vector1< AA > const & aas ){
//    excluded_aatypes_for_selection_ = aas;
//}


/// @brief accessbile surface are of each residue
core::Real
LayerOperations::rsd_sasa( Size const i ) const{
    return rsd_sasa_[ i ];
}

/// @brief return defined layer for each residue
std::string
LayerOperations::layer( Size const i ) const{
    return rsd_layer_[ i ];
}

/// @brief selected residues on boundary
utility::vector1< core::Size > const &
LayerOperations::selected_boundary_residues() const{
    return selected_boundary_residues_;
}


/// @brief selected residues in core
utility::vector1< core::Size > const &
LayerOperations::selected_core_residues() const{
    return selected_core_residues_;
}


/// @brief selected residues on surface
utility::vector1< core::Size > const &
LayerOperations::selected_surface_residues() const{
    return selected_surface_residues_;
}


/// @brief return accessible surface area for each residue
utility::vector1< core::Real > const
LayerOperations::calc_rsd_sasa( Pose const & pose ) const {
    
    // define atom_map for main-chain and CB
    core::id::AtomID_Map< bool > atom_map;
    core::pose::initialize_atomid_map( atom_map, pose, false );
    for ( Size ir = 1; ir <= pose.total_residue(); ++ir ) {
        for ( Size j = 1; j<=5; ++j ) {
            core::id::AtomID atom( j, ir );
            atom_map.set( atom, true );
        }
    }
    
    // calc sasa
    core::id::AtomID_Map< Real > atom_sasa;
    utility::vector1< Real > rsd_sasa;
    core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, pore_radius_, false, atom_map );
    
    return rsd_sasa;
} // calc_residue_sasa
    
    
/// @brief
utility::vector1< core::Size > const
LayerOperations::compute( Pose const & pose )
{

    String secstruct = pose.secstruct();
    runtime_assert( pose.total_residue() == secstruct.length() );
    
    // clear
    rsd_sasa_.clear();
    rsd_layer_.clear();
    selected_core_residues_.clear();
    selected_surface_residues_.clear();
    selected_boundary_residues_.clear();
    
    std::ofstream output;
    if( make_rasmol_format_file_ ) {
        output.open( "srb.ras" ,std::ios::out );
    }
    
    TR << " pore_radius : " <<  pore_radius_ << std::endl;
    TR << " core ( E, L, H ): " << burial_[ 'E' ] << ' ' << burial_[ 'L' ] << ' ' << burial_[ 'H' ] << std::endl;
    TR << " surface (E, L, H ): " << surface_[ 'E' ] << ' ' << surface_[ 'L' ] << ' ' << surface_[ 'H' ] << std::endl;
    
    rsd_sasa_ = calc_rsd_sasa( pose );
    rsd_layer_.resize( pose.total_residue() );
    
    Size start_res = 1;
    Size end_res = pose.total_residue();
    
    if( chain_ > 0 ){
        start_res = pose.conformation().chain_begin( chain_ );
		end_res   =	pose.conformation().chain_end( chain_ );
    }
    
    utility::vector1< Size > selected_residues;
    for( Size iaa=start_res; iaa<= end_res ; iaa++ ) {
        
        char ss = secstruct.at( iaa-1 );
        runtime_assert( ss == 'L' || ss =='E' || ss=='H' );
        rsd_layer_[ iaa ] = "";
        
        if ( rsd_sasa_[ iaa ] <= burial_[ ss ] && pick_core_ ) {
            selected_core_residues_.push_back( iaa );
            selected_residues.push_back( iaa );
            
            if( make_rasmol_format_file_ ) {
                output << "select " << iaa << std::endl;
                output << "color blue " << std::endl;
            }
            rsd_layer_[ iaa ] = "core";
            
        } else if ( rsd_sasa_[ iaa ] >= surface_[ ss ] && pick_surface_ ) {
            selected_surface_residues_.push_back( iaa );
            selected_residues.push_back( iaa );
            
            if( make_rasmol_format_file_ ) {
                output << "select " << iaa << std::endl;
                output << "color red " << std::endl;
            }
            rsd_layer_[ iaa ] = "surface";
            
        } else if ( rsd_sasa_[ iaa ] < surface_[ ss ] && rsd_sasa_[ iaa ] > burial_[ ss ] && pick_boundary_ ) {
            selected_boundary_residues_.push_back( iaa );
            selected_residues.push_back( iaa );
            
            if( make_rasmol_format_file_ ) {
                output << "select " << iaa << std::endl;
                output << "color green " << std::endl;
            }
            rsd_layer_[ iaa ] = "boundary";
        }
    }
    
    return selected_residues;
    
} // compute

bool 
is_part ( utility::vector1< core::Size > vec, core::Size pos ){
    bool is_part_of_this = false;
    for( core::Size i = 1; i <= vec.size(); i++ )
        if( i == pos ) is_part_of_this = true;
    return is_part_of_this;
}
    
void
LayerOperations::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const{
    
    using namespace core::pack::task::operation;
    
    // separate bound partners if specified
	protocols::rigid::RigidBodyTransMoverOP translate;
	if( (jump_ > 0) && (pose.conformation().num_chains() > 1) ) {
		TR<<"translating along jump #"<<jump_<<std::endl;
		translate = new protocols::rigid::RigidBodyTransMover( pose, jump_ ) ;
		translate->step_size( 1000.0 );
		translate->apply( pose );
	}
    
	// calc layers based on residue asa
	utility::vector1< core::Size > selected_residues( compute( pose )) ;
        
    ///for some unfathomable reason OperateOnCertainResidues defaults to applying to all residues if none are defined, so you have to be careful here...
	OperateOnCertainResidues oocr_repacking, oocr_prevent_repacking;
    
    utility::vector1 <core::Size > repacking_residues;
    utility::vector1 <core::Size > norepack_residues;
    
    // since task operations are all about restricting and once a freedom has been restricted it cannot be regained. hence all negative here...
    for (Size res = 1; res <= pose.total_residue() ; res++ ){
        if ( !repack_rest_ ){
        	if (!is_part(selected_residues, res)) norepack_residues.push_back(res);
        }
        if ( design_ ){
        	if(!is_part(selected_residues, res)) repacking_residues.push_back(res);
        }
    }
	
    if( repacking_residues.size() ){
		oocr_repacking.op( new RestrictToRepackingRLT );
		oocr_repacking.residue_indices( repacking_residues );
		oocr_repacking.apply( pose, task );
	}
	if( norepack_residues.size() ){
		oocr_prevent_repacking.op( new PreventRepackingRLT );
		oocr_prevent_repacking.residue_indices( norepack_residues );
		oocr_prevent_repacking.apply( pose, task );
	}
    
    // move back together
	if( (jump_ > 0) && (pose.conformation().num_chains() > 1) ) {
		translate->trans_axis().negate();
		translate->apply( pose );
	}
}


void
LayerOperations::parse_tag( TagCOP tag , DataMap & )
{
	//use_original_ = tag->getOption< bool >( "use_original_non_designed_layer", 1 );

	String design_layers = tag->getOption< String >( "layer" );
	utility::vector1< String > layers( utility::string_split( design_layers, '_' ) );

    design_ = tag->getOption< bool > ("design" , 1 ); 
    
    repack_rest_ = tag->getOption< bool >("repack_rest" , 1 );
    
    chain_ = tag->getOption< Size > ("chain", 0 );
    
    jump_ = tag->getOption< Size > ("jump", 0 );

    //if ( jump_ > 0 && pose.conformation().num_chains() == 1 )
    //	utility_exit_with_message("need at least two chains for jump > 0 ");
        
	for ( utility::vector1< String >::const_iterator iter = layers.begin(); iter != layers.end() ; ++iter) {
		String layer(*iter);

		if ( layer == "core" ) {
			dsgn_core = true;
		} else if ( layer == "surface" ) {
			dsgn_surface = true;
		} else if ( layer == "boundary" ) {
			dsgn_boundary = true;
		} else {
			TR << "Error! wrong specification of layer_mode " << layer << std::endl;
			TR << "All layers are designed. " << std::endl;
			dsgn_core = true;
			dsgn_surface = true;
			dsgn_boundary = true;
		}
	}

	srbl_->set_design_layer( dsgn_core, dsgn_boundary, dsgn_surface );

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

	//set_verbose( tag->getOption< bool >( "verbose", false ) );


}

}
}
}
