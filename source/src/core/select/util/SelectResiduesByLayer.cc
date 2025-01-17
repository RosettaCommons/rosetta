// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/util/SelectResiduesByLayer.cc
/// @brief select residues depending on layer: core, boundary, and surface
/// The layer of residues are defined by accessible sufrace area of CB + mainchain,
/// or by the number of sidechains within a cone along the CA-CB vector.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )
/// @author Gabe Rocklin (sidechain neighbour selection)
/// @author Vikram K. Mulligan (vmullig@uw.edu -- moving this class to core and refactoring for noncanonicals)

// Unit Headers
#include <core/select/util/SelectResiduesByLayer.hh>

// Project Headers
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <basic/Tracer.hh>
#include <core/select/util/burial_utilities.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <fstream>

//Auto Headers
#include <core/pose/init_id_map.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.select.util.SelectResiduesByLayer" );

namespace core {
namespace select {
namespace util {

/// @brief default constructor
SelectResiduesByLayer::SelectResiduesByLayer() : utility::VirtualBase(),
	pick_core_( false ),
	pick_boundary_( false ),
	pick_surface_( false ),
	pore_radius_( 2.0 ),
	make_rasmol_format_file_( false ),
	use_sidechain_neighbors_( false ),
	is_selection_initialized_( false ),
	dist_midpoint_(9.0),
	rsd_neighbor_denominator_(1.0),
	angle_shift_factor_(0.5),
	angle_exponent_(2.0),
	dist_exponent_(1.0)
{
	initialize( 20.0, 40.0 );
}

/// @brief Copy constructor
///
SelectResiduesByLayer::SelectResiduesByLayer( SelectResiduesByLayer const &src ) :
	VirtualBase( src ),
	pick_core_(src.pick_core_),
	pick_boundary_(src.pick_boundary_),
	pick_surface_(src.pick_surface_),
	pore_radius_(src.pore_radius_),
	burial_(src.burial_),
	surface_(src.surface_),
	excluded_aatypes_for_selection_(src.excluded_aatypes_for_selection_),
	restricted_aatypes_for_selection_(src.restricted_aatypes_for_selection_),
	selected_core_residues_(src.selected_core_residues_),
	selected_boundary_residues_(src.selected_boundary_residues_),
	selected_surface_residues_(src.selected_surface_residues_),
	make_rasmol_format_file_(src.make_rasmol_format_file_),
	use_sidechain_neighbors_(src.use_sidechain_neighbors_),
	is_selection_initialized_(src.is_selection_initialized_),
	rsd_sasa_(src.rsd_sasa_),
	rsd_layer_(src.rsd_layer_),
	dist_midpoint_(src.dist_midpoint_),
	rsd_neighbor_denominator_(src.rsd_neighbor_denominator_),
	angle_shift_factor_(src.angle_shift_factor_),
	angle_exponent_(src.angle_exponent_),
	dist_exponent_(src.dist_exponent_)
{}

/// @brief value constructor
SelectResiduesByLayer::SelectResiduesByLayer(
	bool const pick_core,
	bool const pick_boundary,
	bool const pick_surface ) :
	utility::VirtualBase(),
	pick_core_( pick_core ),
	pick_boundary_( pick_boundary ),
	pick_surface_( pick_surface ),
	pore_radius_( 2.0 ),
	make_rasmol_format_file_( false ),
	use_sidechain_neighbors_( false ),
	is_selection_initialized_( false ),
	dist_midpoint_(9.0),
	rsd_neighbor_denominator_(1.0),
	angle_shift_factor_(0.5),
	angle_exponent_(2.0),
	dist_exponent_(1.0)
{
	initialize( 20.0, 40.0 );
}

/// @brief value constructor
SelectResiduesByLayer::SelectResiduesByLayer( String const &pick ) : utility::VirtualBase(),
	pick_core_( false ),
	pick_boundary_( false ),
	pick_surface_( false ),
	pore_radius_( 2.0 ),
	make_rasmol_format_file_( false ),
	use_sidechain_neighbors_( false ),
	is_selection_initialized_( false ),
	dist_midpoint_(9.0),
	rsd_neighbor_denominator_(1.0),
	angle_shift_factor_(0.5),
	angle_exponent_(2.0),
	dist_exponent_(1.0)
{
	initialize( 20.0, 40.0 );
	utility::vector1< String > layers( utility::string_split( pick, '_' ) );
	for ( utility::vector1< String >::const_iterator iter = layers.begin(); iter != layers.end() ; ++iter ) {
		String layer(*iter);
		if ( layer == "core" ) {
			pick_core_ = true;
		} else if ( layer == "surface" ) {
			pick_surface_ = true;
		} else if ( layer == "boundary" ) {
			pick_boundary_ = true;
		} else {
			TR.Fatal << "Error!, wrong specification of layer_mode " << layer << std::endl;
			utility_exit_with_message("Layer mode specification not understood.");
		}
	} // utility::vector1
}

/// @brief destructor
SelectResiduesByLayer::~SelectResiduesByLayer() = default;

/// @brief Clone operator.
/// @details Construct a copy of this object and return an owning pointer to the copy.
SelectResiduesByLayerOP
SelectResiduesByLayer::clone() const {
	return utility::pointer::make_shared< SelectResiduesByLayer >(*this);
}

/// @brief accessible surface for evaluating residues are in surface or not
void
SelectResiduesByLayer::sasa_surface( Real const &r, String const &ss )
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
SelectResiduesByLayer::sasa_core( Real const &r, String const &ss )
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
SelectResiduesByLayer::initialize( Real const &burial, Real const &surface )
{
	surface_.insert( std::map< char, Real >::value_type( 'H', surface ) );
	surface_.insert( std::map< char, Real >::value_type( 'L', surface ) );
	surface_.insert( std::map< char, Real >::value_type( 'E', surface ) );
	burial_.insert( std::map< char, Real >::value_type( 'H', burial ) );
	burial_.insert( std::map< char, Real >::value_type( 'L', burial ) );
	burial_.insert( std::map< char, Real >::value_type( 'E', burial ) );
}


/// @brief amino acid types excluded for selection
void
SelectResiduesByLayer::exclude_aatypes_for_selection( utility::vector1< AA > const & aas )
{
	excluded_aatypes_for_selection_ = aas;
}


/// @brief amino acid types restricted for selection
void
SelectResiduesByLayer::restrict_aatypes_for_selection( utility::vector1< AA > const & aas )
{
	restricted_aatypes_for_selection_ = aas;
}

bool
SelectResiduesByLayer::use_sidechain_neighbors() const
{
	return use_sidechain_neighbors_;
}

/// @brief accessbile surface are of each residue
core::Real
SelectResiduesByLayer::rsd_sasa( Size const i ) const
{
	return rsd_sasa_[ i ];
}

/// @brief return defined layer for each residue
std::string
SelectResiduesByLayer::layer( Size const i ) const
{
	return rsd_layer_[ i ];
}

/// @brief selected residues on boundary
utility::vector1< core::Size > const &
SelectResiduesByLayer::selected_boundary_residues() const
{
	return selected_boundary_residues_;
}


/// @brief selected residues in core
utility::vector1< core::Size > const &
SelectResiduesByLayer::selected_core_residues() const
{
	return selected_core_residues_;
}


/// @brief selected residues on surface
utility::vector1< core::Size > const &
SelectResiduesByLayer::selected_surface_residues() const
{
	return selected_surface_residues_;
}


/// @brief return accessible surface area for each residue
utility::vector1< core::Real > const
SelectResiduesByLayer::calc_rsd_sasa( Pose const & pose ) const {

	// define atom_map for main-chain and CB
	core::id::AtomID_Map< bool > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, false );
	for ( Size ir = 1; ir <= pose.size(); ++ir ) {
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
SelectResiduesByLayer::compute( Pose const & pose, String const &secstruct, bool const skip_dssp, bool const asu_only )
{
	String secstruct2(secstruct);

	// define secstruct if it is empty
	if ( secstruct2 == "" ) {
		if ( skip_dssp ) {
			for ( core::Size i=1, imax=pose.size(); i<=imax; ++i ) {
				secstruct2 = secstruct2 + "L";
			}
		} else {
			// calc dssp
			core::scoring::dssp::Dssp dssp( pose );
			dssp.dssp_reduced();
			secstruct2 = dssp.get_dssp_secstruct();
		}
	}
	runtime_assert( pose.size() == secstruct2.length() );

	// clear
	rsd_sasa_.clear();
	rsd_layer_.clear();
	selected_core_residues_.clear();
	selected_surface_residues_.clear();
	selected_boundary_residues_.clear();

	std::ofstream output;
	if ( make_rasmol_format_file_ ) {
		output.open( "srb.ras" ,std::ios::out );
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << " pore_radius : " <<  pore_radius_ << std::endl;
		TR.Debug << " core ( E, L, H ): " << burial_[ 'E' ] << ' ' << burial_[ 'L' ] << ' ' << burial_[ 'H' ] << std::endl;
		TR.Debug << " surface (E, L, H ): " << surface_[ 'E' ] << ' ' << surface_[ 'L' ] << ' ' << surface_[ 'H' ] << std::endl;
	}

	//calculate layer identity
	if ( use_sidechain_neighbors_ ) {
		rsd_sasa_ = calc_sc_neighbors( pose, angle_exponent(), angle_shift_factor(), dist_exponent(), dist_midpoint(), rsd_neighbor_denominator(), asu_only );
	} else {
		rsd_sasa_ = calc_rsd_sasa( pose );
	}

	core::Size target_pose_size;

	//check for symmetry
	bool sym_check = core::pose::symmetry::is_symmetric( pose );
	if ( sym_check && asu_only ) {
		TR.Info << "Pose is symmetric and asu_only=True. Returning residues on asymmetric unit ONLY." << std::endl;
		target_pose_size = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		TR.Debug << "Pose.size() = " << pose.size() << std::endl;
		TR.Debug << "ASU size = " << target_pose_size << std::endl;
	} else {
		target_pose_size = pose.size();
	}

	rsd_layer_.resize( target_pose_size );

	//apply layer identity to selection
	utility::vector1< Size > selected_residues;

	for ( Size iaa=1; iaa<=target_pose_size; iaa++ ) {

		//TR.Debug << "resi = " << iaa << std::endl;

		char ss = secstruct2.at( iaa-1 );
		runtime_assert( ss == 'L' || ss =='E' || ss=='H' );

		rsd_layer_[ iaa ] = "";

		bool flag( false );
		for ( Size i=1; i<=excluded_aatypes_for_selection_.size(); i++ ) {
			if ( pose.aa( iaa ) == excluded_aatypes_for_selection_[ i ] ) {
				flag = true;
				break;
			}
		}
		if ( flag ) continue;

		if ( restricted_aatypes_for_selection_.size() > 0 ) flag = true;
		for ( Size i=1; i<=restricted_aatypes_for_selection_.size(); i++ ) {
			if ( pose.aa( iaa ) == restricted_aatypes_for_selection_[ i ] ) {
				flag = false;
				break;
			}
		}
		if ( flag ) continue;

		if ( pick_core_ && ( (rsd_sasa_[ iaa ] <= burial_[ ss ] && !use_sidechain_neighbors_ ) ||
				(rsd_sasa_[ iaa ] >= burial_[ ss ] &&  use_sidechain_neighbors_ ) ) ) {

			selected_core_residues_.push_back( iaa );
			selected_residues.push_back( iaa );

			if ( make_rasmol_format_file_ ) {
				output << "select " << iaa << std::endl;
				output << "color blue " << std::endl;
			}

			rsd_layer_[ iaa ] = "core";

		} else if ( pick_surface_ && ( (rsd_sasa_[ iaa ] >= surface_[ ss ] && !use_sidechain_neighbors_ ) ||
				(rsd_sasa_[ iaa ] <= surface_[ ss ] &&  use_sidechain_neighbors_ ) ) ) {

			selected_surface_residues_.push_back( iaa );
			selected_residues.push_back( iaa );

			if ( make_rasmol_format_file_ ) {
				output << "select " << iaa << std::endl;
				output << "color red " << std::endl;
			}

			rsd_layer_[ iaa ] = "surface";

		} else if ( pick_boundary_ && ( ( rsd_sasa_[ iaa ] < surface_[ ss ] && rsd_sasa_[ iaa ] > burial_[ ss ] && !use_sidechain_neighbors_ ) ||
				( rsd_sasa_[ iaa ] > surface_[ ss ] && rsd_sasa_[ iaa ] < burial_[ ss ] &&  use_sidechain_neighbors_ ) ) ) {

			selected_boundary_residues_.push_back( iaa );
			selected_residues.push_back( iaa );

			if ( make_rasmol_format_file_ ) {
				output << "select " << iaa << std::endl;
				output << "color green " << std::endl;
			}

			rsd_layer_[ iaa ] = "boundary";

		}

	}

	is_selection_initialized_ = true;

	return selected_residues;

} // compute

} // util
} // select
} // core

#ifdef SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::util::SelectResiduesByLayer::save( Archive & arc ) const {
	arc( CEREAL_NVP( pick_core_ ) ); // bool
	arc( CEREAL_NVP( pick_boundary_ ) ); // bool
	arc( CEREAL_NVP( pick_surface_ ) ); // bool
	arc( CEREAL_NVP( cache_selection_ ) ); // bool
	arc( CEREAL_NVP( pore_radius_ ) ); // core::Real
	arc( CEREAL_NVP( burial_ ) ); // std::map< char, core::Real >
	arc( CEREAL_NVP( surface_ ) ); // std::map< char, core::Real >
	arc( CEREAL_NVP( excluded_aatypes_for_selection_ ) ); // utility::vector1< AA >
	arc( CEREAL_NVP( restricted_aatypes_for_selection_ ) ); // utility::vector1< AA >
	arc( CEREAL_NVP( selected_core_residues_ ) ); // utility::vector1< Size >
	arc( CEREAL_NVP( selected_boundary_residues_ ) ); // utility::vector1< Size >
	arc( CEREAL_NVP( selected_surface_residues_ ) ); // utility::vector1< Size >
	arc( CEREAL_NVP( make_rasmol_format_file_ ) ); // bool
	arc( CEREAL_NVP( use_sidechain_neighbors_ ) ); // bool
	arc( CEREAL_NVP( is_selection_initialized_ ) ); // bool
	arc( CEREAL_NVP( rsd_sasa_ ) ); // utility::vector1< Real >
	arc( CEREAL_NVP( rsd_layer_ ) ); // utility::vector1< String >
	arc( CEREAL_NVP( dist_midpoint_ ) ); // core::Real
	arc( CEREAL_NVP( rsd_neighbor_denominator_ ) ); // core::Real
	arc( CEREAL_NVP( angle_shift_factor_ ) ); // core::Real
	arc( CEREAL_NVP( angle_exponent_ ) ); // core::Real
	arc( CEREAL_NVP( dist_exponent_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::util::SelectResiduesByLayer::load( Archive & arc ) {
	arc( pick_core_ ); // bool
	arc( pick_boundary_ ); // bool
	arc( pick_surface_ ); // bool
	arc( cache_selection_ ); // bool
	arc( pore_radius_ ); // core::Real
	arc( burial_ ); // std::map< char, core::Real >
	arc( surface_ ); // std::map< char, core::Real >
	arc( excluded_aatypes_for_selection_ ); // utility::vector1< AA >
	arc( restricted_aatypes_for_selection_ ); // utility::vector1< AA >
	arc( selected_core_residues_ ); // utility::vector1< Size >
	arc( selected_boundary_residues_ ); // utility::vector1< Size >
	arc( selected_surface_residues_ ); // utility::vector1< Size >
	arc( make_rasmol_format_file_ ); // bool
	arc( use_sidechain_neighbors_ ); // bool
	arc( is_selection_initialized_ ); // bool
	arc( rsd_sasa_ ); // utility::vector1< Real >
	arc( rsd_layer_ ); // utility::vector1< String >
	arc( dist_midpoint_ ); // core::Real
	arc( rsd_neighbor_denominator_ ); // core::Real
	arc( angle_shift_factor_ ); // core::Real
	arc( angle_exponent_ ); // core::Real
	arc( dist_exponent_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::util::SelectResiduesByLayer );
CEREAL_REGISTER_TYPE( core::select::util::SelectResiduesByLayer )

CEREAL_REGISTER_DYNAMIC_INIT( core_select_util_SelectResiduesByLayer )
#endif // SERIALIZATION
