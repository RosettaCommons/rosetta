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
// AUTO-REMOVED #include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask_.hh>
#include <core/scoring/dssp/Dssp.hh>
// AUTO-REMOVED #include <core/scoring/sasa.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <protocols/flxbb/SelectResiduesByLayer.hh>


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
	srbl_( new SelectResiduesByLayer )
{}


/// @brief value constructor
LayerDesignOperation::LayerDesignOperation( bool dsgn_core, bool dsgn_boundary, bool dsgn_surface ):
	TaskOperation(),
	add_helix_capping_( true ),
	use_original_( true ),
	verbose_( false ),
	srbl_( new SelectResiduesByLayer )
{
	design_layer( dsgn_core, dsgn_boundary, dsgn_surface );
}


/// @brief destructor
LayerDesignOperation::~LayerDesignOperation() {}

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


/// @brief apply
void
LayerDesignOperation::apply( Pose const & pose, PackerTask & task ) const
{
	// calc dssp
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();

	// calc SelectResiduesByLayer
	srbl_->compute( pose, dssp.get_dssp_secstruct() );

	// find the position of residues of helix capping and intial residue of helix
	bool flag( false );
	utility::vector1< bool > helix_capping( pose.total_residue(), false );
	utility::vector1< bool > initial_helix( pose.total_residue(), false );
	for( Size i=1; i<=pose.total_residue(); ++i ){
		char ss( dssp.get_dssp_secstruct( i ) );
		if( ss == 'H' && flag == false && i != 1 ){
			initial_helix[ i ] = true;
			helix_capping[ i-1 ] = true;
			flag = true;
		}
		if( ss != 'H' && flag == true ){
			flag = false;
		}
	}

	// terminal residues set to be allaa
	utility::vector1<bool> restrict_to_aa( 20, true );
	task.nonconst_residue_task( 1 ).restrict_absent_canonical_aas( restrict_to_aa );
	task.nonconst_residue_task( pose.total_residue() ).restrict_absent_canonical_aas( restrict_to_aa );

	for( Size i=2; i<=pose.total_residue()-1; ++i ) {

		char ss( dssp.get_dssp_secstruct( i ) );
		if( verbose_ ) {
			TR << "Resnum=" << i << " ,SS=" << ss << " "
				 << "Sasa=" << ObjexxFCL::fmt::F( 6, 2, srbl_->rsd_sasa( i ) );
		}

		// skip the residue if this position is defined as PIKAA by resfile
		if( task.residue_task( i ).command_string().find( "PIKAA" ) != std::string::npos ){
			if( verbose_ ) {
				TR << " ,Resfile info is used." << std::endl;
			}
			continue;
		}

		if( srbl_->layer( i ) == "core" ) {

			if( helix_capping[ i ] == true && add_helix_capping_ ) {
				utility::vector1<bool> restrict_to_aa( chemical::num_canonical_aas, false );
				restrict_to_aa[chemical::aa_from_name( "ASP" )] = true;
				restrict_to_aa[chemical::aa_from_name( "ASN" )] = true;
				restrict_to_aa[chemical::aa_from_name( "THR" )] = true;
				restrict_to_aa[chemical::aa_from_name( "SER" )] = true;
				task.nonconst_residue_task( i ).restrict_absent_canonical_aas( restrict_to_aa );
				if( verbose_ ) {
					TR << " ,Helix Capping " << std::endl;
				}
			} else {
				utility::vector1<bool> restrict_to_aa( chemical::num_canonical_aas, false );
				restrict_to_aa[chemical::aa_from_name( "TRP" )] = true;
				restrict_to_aa[chemical::aa_from_name( "PHE" )] = true;
				restrict_to_aa[chemical::aa_from_name( "ILE" )] = true;
				restrict_to_aa[chemical::aa_from_name( "LEU" )] = true;
				restrict_to_aa[chemical::aa_from_name( "VAL" )] = true;
				restrict_to_aa[chemical::aa_from_name( "TYR" )] = true;
				restrict_to_aa[chemical::aa_from_name( "MET" )] = true;

				if( ss == 'E' ) {

				} else if( ss == 'H' ) {
					restrict_to_aa[chemical::aa_from_name( "ALA" )] = true;
					if( initial_helix[ i ] ) {
						restrict_to_aa[chemical::aa_from_name( "PRO" )] = true;
					}
				} else {
					restrict_to_aa[chemical::aa_from_name( "ALA" )] = true;
					restrict_to_aa[chemical::aa_from_name( "PRO" )] = true;
				}
				task.nonconst_residue_task( i ).restrict_absent_canonical_aas( restrict_to_aa );
				if( verbose_ ) {
					TR << " " << srbl_->layer( i ) << std::endl;
				}
			}

		} else if ( srbl_->layer( i ) == "boundary" ) {

			if( helix_capping[ i ] == true && add_helix_capping_ ) {
				utility::vector1<bool> restrict_to_aa( chemical::num_canonical_aas, false );
				restrict_to_aa[chemical::aa_from_name( "ASP" )] = true;
				restrict_to_aa[chemical::aa_from_name( "ASN" )] = true;
				restrict_to_aa[chemical::aa_from_name( "THR" )] = true;
				restrict_to_aa[chemical::aa_from_name( "SER" )] = true;
				task.nonconst_residue_task( i ).restrict_absent_canonical_aas( restrict_to_aa );
				if( verbose_ ) {
					TR << " ,Helix Capping " << std::endl;
				}
			} else {

				utility::vector1<bool> restrict_to_aa( chemical::num_canonical_aas, true );
				if( ss == 'E' || ss == 'L' ) {

					restrict_to_aa[chemical::aa_from_name( "CYS" )] = false;
					restrict_to_aa[chemical::aa_from_name( "MET" )] = false;
					restrict_to_aa[chemical::aa_from_name( "HIS" )] = false;

				} else {

					restrict_to_aa[chemical::aa_from_name( "CYS" )] = false;
					restrict_to_aa[chemical::aa_from_name( "PHE" )] = false;
					restrict_to_aa[chemical::aa_from_name( "MET" )] = false;
					restrict_to_aa[chemical::aa_from_name( "HIS" )] = false;

				}

				if( ss == 'E' ){
					restrict_to_aa[chemical::aa_from_name( "GLY" )] = false;
					restrict_to_aa[chemical::aa_from_name( "ALA" )] = false;
					restrict_to_aa[chemical::aa_from_name( "PRO" )] = false;
				} else if( ss == 'H' ) {
					if( ! initial_helix[ i ] ) {
						restrict_to_aa[chemical::aa_from_name( "PRO" )] = false;
					}
					restrict_to_aa[chemical::aa_from_name( "GLY" )] = false;
				}
				task.nonconst_residue_task( i ).restrict_absent_canonical_aas( restrict_to_aa );
				if( verbose_ ) {
					TR << " " << srbl_->layer( i ) << std::endl;
				}
			}

		} else if ( srbl_->layer( i ) == "surface" ) {

			if( helix_capping[ i ] == true && add_helix_capping_ ) {
				utility::vector1<bool> restrict_to_aa( chemical::num_canonical_aas, false );
				restrict_to_aa[chemical::aa_from_name( "ASP" )] = true;
				restrict_to_aa[chemical::aa_from_name( "ASN" )] = true;
				restrict_to_aa[chemical::aa_from_name( "THR" )] = true;
				restrict_to_aa[chemical::aa_from_name( "SER" )] = true;
				task.nonconst_residue_task( i ).restrict_absent_canonical_aas( restrict_to_aa );
				if( verbose_ ) {
					TR << " ,Helix Capping " << std::endl;
				}
			} else {

				utility::vector1<bool> restrict_to_aa( chemical::num_canonical_aas, false );
				restrict_to_aa[chemical::aa_from_name( "GLU" )] = true;
				restrict_to_aa[chemical::aa_from_name( "ARG" )] = true;
				restrict_to_aa[chemical::aa_from_name( "ASP" )] = true;
				restrict_to_aa[chemical::aa_from_name( "LYS" )] = true;
				restrict_to_aa[chemical::aa_from_name( "HIS" )] = true;
				restrict_to_aa[chemical::aa_from_name( "ASN" )] = true;
				restrict_to_aa[chemical::aa_from_name( "GLN" )] = true;
				restrict_to_aa[chemical::aa_from_name( "SER" )] = true;
				restrict_to_aa[chemical::aa_from_name( "THR" )] = true;

				if( ss == 'E' ) {
				} else if ( ss == 'H' ) {
					if( initial_helix[ i ] == true ){
						restrict_to_aa[chemical::aa_from_name( "PRO" )] = true;
					}
				} else {
					restrict_to_aa[chemical::aa_from_name( "GLY" )] = true;
				}

				task.nonconst_residue_task( i ).restrict_absent_canonical_aas( restrict_to_aa );
				if( verbose_ ) {
					TR << " " << srbl_->layer( i ) << std::endl;
				}
			}

		} else {

			if( use_original_ ){
				task.nonconst_residue_task( i ).restrict_to_repacking();
				if( verbose_ ) {
					TR << " ,Original sequence used" << std::endl;
				}
			}else{
				utility::vector1<bool> restrict_to_aa( 20, false );
				restrict_to_aa[chemical::aa_from_name( "ALA" )] = true;
				task.nonconst_residue_task( i ).restrict_absent_canonical_aas( restrict_to_aa );
				if( verbose_ ) {
					TR << " ,No design" << std::endl;
				}
			}

		}

	} // for( i )
} // apply

void
LayerDesignOperation::parse_tag( TagPtr tag )
{
	use_original_ = tag->getOption< bool >( "use_original_non_designed_layer", 1 );

	String design_layers = tag->getOption< String >( "layer" );
	utility::vector1< String > layers( utility::string_split( design_layers, '_' ) );

	bool dsgn_core( false );
	bool dsgn_surface( false );
	bool dsgn_boundary( false );
	for ( utility::vector1< String >::const_iterator iter = layers.begin(); iter != layers.end() ; ++iter) {
		String layer(*iter);

		if ( layer == "core" ) {
			dsgn_core = true;
		} else if ( layer == "surface" ) {
			dsgn_surface = true;
		} else if ( layer == "boundary" ) {
			dsgn_boundary = true;
		} else {
			TR << "Error!, wrong specification of layer_mode " << layer << std::endl;
			TR << "Every layers are designed. " << std::endl;
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

	set_verbose( tag->getOption< bool >( "verbose", false ) );


}

} // flxbb
} // protocols

