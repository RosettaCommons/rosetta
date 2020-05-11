// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/EPRSpinLabel.cc
/// @brief  This is a container class specific to the use of double electron-electron resonance data
/// @details This container manages the simulated distance distributions for the deer_decay and
///      deer_distance energy method. It also manages individual electron coordinate ensembles for
///      a given protein, although a reference to that pose is not stored in this object.
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

// Unit headers
#include <core/scoring/epr_deer/EPRSpinLabel.hh>

// Package headers
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/select/util.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/types.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>


// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/epr_deer.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>
#include <numeric/HomogeneousTransform.hh>

// C++ headers
#include <stdlib.h>
#include <iostream>
#include <set>

namespace core {
namespace scoring {
namespace epr_deer {

static basic::Tracer TR( "core.scoring.epr_deer.EPRSpinLabel" );

//////////////////////
// DEFAULT SPIN LABELS

// These coordinates describe the electron position on a dynamic MTSSL side chain in a coarse-grained manner.
// The XYZ coordinates are each assigned a weight (second member of the pair) to characterize a given electron's
// propensity to occupy that given area. This works as well as existing methods without relying on more expensive
// fullatom representations. Clashes are evaluated using a dummy atom (see below). Custom positions may also be
// provided. These are introduced using a homogeneous transform object - see corresponding cc file.
// The larger vector is more accurate but about 4-fold slower. For de novo folding, users are encouraged to use
// the spin label DEFAULT_FAST (the first one), whereas for things like comparative modeling DEFAULT (the second)
// is more appropriate.
static const
utility::vector1< PseudoElectron > mtssl13_ = {
std::make_pair( numeric::xyzVector< Real >( -5.767,  2.275, -0.295 ), 0.5447 ),
std::make_pair( numeric::xyzVector< Real >(  3.613, -6.337,  0.800 ), 0.7091 ),
std::make_pair( numeric::xyzVector< Real >( -1.705, -1.234,  5.859 ), 1.0000 ),
std::make_pair( numeric::xyzVector< Real >( -5.623, -4.733, -2.694 ), 0.2072 ),
std::make_pair( numeric::xyzVector< Real >( -2.224, -6.089,  5.277 ), 0.3287 ),
std::make_pair( numeric::xyzVector< Real >(  5.572, -2.864,  3.288 ), 0.4498 ),
std::make_pair( numeric::xyzVector< Real >(  1.108, -6.893, -2.574 ), 0.9798 ),
std::make_pair( numeric::xyzVector< Real >( -3.732,  1.880,  7.013 ), 0.4694 ),
std::make_pair( numeric::xyzVector< Real >(  2.653,  0.153,  6.667 ), 0.6293 ),
std::make_pair( numeric::xyzVector< Real >( -3.585, -7.244,  2.537 ), 0.1551 ),
std::make_pair( numeric::xyzVector< Real >( -1.877, -5.513, -1.528 ), 0.2841 ),
std::make_pair( numeric::xyzVector< Real >(  0.386, -6.957,  4.280 ), 0.3151 ),
std::make_pair( numeric::xyzVector< Real >( -0.407, -6.350, -5.107 ), 0.6001 )
};

static const
utility::vector1< PseudoElectron > mtssl50_unweighted_ = {
std::make_pair( numeric::xyzVector< core::Real > ( -3.02, -6.385, 4.598), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -3.171, -5.602, -1.129), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -5.965, 1.611, -1.727), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 1.721, -7.535, 3.201), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -7.095, -5.243, 2.081), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 5.401, -2.96, 4.29), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 1.535, -5.371, -4.155), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 2.451, -0.607, 6.988), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -8.324, 0.475, 3.255), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -6.721, -4.965, -2.153), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 1.768, -7.644, -1.366), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -5.965, -3.659, 5.684), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -3.628, 2.644, 6.662), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -7.029, -1.377, -5.071), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -8.731, -3.654, 0.004), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -3.698, -7.508, 1.932), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 4.807, -5.037, 1.085), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -5.315, 2.65, -0.775), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 0.822, -6.631, -3.274), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 0.941, -5.792, 5.57), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 2.066, -7.341, 1.191), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -2.384, -5.545, 6.003), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -0.505, -5.657, -5.723), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -8.803, -1.771, 0.792), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -6.585, -5.547, 4.04), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -5.742, -4.1, -4.064), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -1.205, -7.119, 3.747), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -6.679, 0.592, -4.703), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 3.965, -6.634, 0.124), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -0.491, -5.153, -2.001), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -7.825, 0.17, 5.212), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -0.309, -7.044, -4.492), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -2.704, -0.365, 5.652), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 5.743, -2.769, 2.286), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -6.208, 1.811, 0.086), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 2.854, 0.912, 6.346), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -4.644, -6.581, 2.998), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -5.152, 3.53, 0.742), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 0.087, -7.381, 4.601), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -4.32, -4.43, -1.535), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( 0.307, -7.925, -1.5), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -3.837, 1.117, 7.365), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -7.49, -4.188, -0.854), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -2.413, -7.644, 2.679), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -5.711, -5.438, -3.024), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -1.671, -6.489, 4.795), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -1.969, -5.785, -1.454), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -5.895, -5.548, 4.9), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -0.513, -1.394, 5.729), 0.02 ),
std::make_pair( numeric::xyzVector< core::Real >( -2.159, -1.769, 5.96), 0.02 )
};

/// @brief Constructor
EPRSpinLabel::EPRSpinLabel() {}

/// @brief Destructor
EPRSpinLabel::~EPRSpinLabel() {}

/// @brief Operator to return electrons from specific residue
utility::vector1< PseudoElectron > &
EPRSpinLabel::operator[](  std::pair< Size, std::string > const & res ) {
	return mapped_coords_[ res ];
}

/// @brief Allows const spin label data to be accessed
utility::vector1< PseudoElectron > const &
EPRSpinLabel::at( std::pair< Size, std::string > const & res ) const {
	return mapped_coords_.at( res );
}

/// @brief Returns a histogram between all coordinates for all residues - assumes complete labeling
std::map< Size, Real >
EPRSpinLabel::histogram(
	utility::vector1< std::pair< Size, std::string > > const & residues,
	Size const & bins_per_a,
	Real const & modifier // = 0.0
) {
	std::map< Size, Real > output;
	// This iterates through all pairs of residues. There should be ( n * ( n - 1 ) ) / 2 such pairs
	// In the typical case where only two residues are spin labeled, it will simply go through these loops once
	for ( Size i = 1; i < residues.size(); ++i ) {
		for ( Size j = i + 1; j <= residues.size(); ++j ) {
			auto distribution = normalize_distribution( histogram( residues[ i ], residues[ j ], bins_per_a, modifier ) );
			for ( auto const & dist : distribution ) {
				if ( output.find( dist.first ) == output.end() ) {
					output[ dist.first ] = 0.0;
				}
				output[ dist.first ] += dist.second;
			}
		}
	}
	return normalize_distribution( output );
}

/// @brief Returns a histogram between all coordinates between two residue/SL combinations
std::map< Size, Real >
EPRSpinLabel::histogram(
	std::pair< Size, std::string > const & res1,
	std::pair< Size, std::string > const & res2,
	Size const & bins_per_a,
	Real const & modifier, // = 0.0
	Real const & stdev // = 1.0
) {
	Size STDEV_RANGE = 4;
	if ( mapped_coords_.find( res1 ) == mapped_coords_.end() || mapped_coords_.find( res2 ) == mapped_coords_.end() ) {
		TR.Error << "Either " << res1.first << " or " << res2.first << " hasn't been labeled!" << std::endl;
	}
	std::map< Size, Real > output;
	auto const & res1_coords = mapped_coords_[ res1 ];
	auto const & res2_coords = mapped_coords_[ res2 ];
	for ( auto const & iter1 : res1_coords ) {
		for ( auto const & iter2 : res2_coords ) {
			Real dist = iter1.first.distance( iter2.first ) + ( modifier * bins_per_a );
			Real comb_weight = iter1.second * iter2.second;
			Size lowest_bin = std::max( round( dist * bins_per_a ) - round( STDEV_RANGE * bins_per_a ), 1.0 );
			Size highest_bin = round( dist * bins_per_a ) + int( STDEV_RANGE * bins_per_a );
			for ( Size bin = lowest_bin; bin <= highest_bin; ++bin ) {
				output[ bin ] += comb_weight * gauss( bin / Real( bins_per_a ), dist, stdev );
			}
		}
	}
	if ( output.empty() ) {
		TR.Error << "Distance distribution with no distances measured: " << res1.first << "\t" << res2.first << std::endl;
	}
	return output;
}

/// @brief Return a value for a gaussian distribution at a particular value, given a average and standard deviation
Real
EPRSpinLabel::gauss(
	Real const & dist,
	Real const & avg,
	Real const & stdev
) const {
	Real numerator = std::exp( -0.5 * ( pow( ( dist - avg ) / stdev, 2 ) ) );
	Real denominator = stdev * pow( 2.0 * numeric::constants::d::pi, 0.5 );
	return numerator / denominator;
}

/// @brief Label a residue with a certain spin label
void
EPRSpinLabel::label(
	Size const & res,
	std::string const & label,
	pose::Pose const & pose,
	bool const & skip_clash_eval // = false
) {
	auto pair = std::make_pair( res, label );
	if ( label == "MTSSL13" || label == "DEFAULT_FAST" ) {
		mapped_coords_[ pair ] = get_coords( res, pose, mtssl13_, skip_clash_eval );
	} else if ( label == "MTSSL50_UNWEIGHTED" || label == "DEFAULT" ) {
		mapped_coords_[ pair ] = get_coords( res, pose, mtssl50_unweighted_, skip_clash_eval );
	} else if ( label == "CUSTOM" ) {
		mapped_coords_[ pair ] = get_coords( res, pose, custom_coords_[ res ], true );
	} else {
		TR.Error << "Must specify DEFAULT, DEFAULT_FAST, or CUSTOM when declaring rotamers to use! Res : " << res << std::endl;
	}
}

/// @brief Normalize distribution so that the sum is equal to one
std::map< Size, Real >
EPRSpinLabel::normalize_distribution(
	std::map< Size, Real > sim_map
) const {
	if ( sim_map.empty() ) {
		TR.Error << "Found distance distribution with zero area under the curve!" << std::endl;
		return std::map< Size, Real >();
	}
	Size start_bin = std::max( sim_map.begin()->first, Size( 1 ) );
	Size end_bin = sim_map.rbegin()->first + 1;
	for ( Size i = start_bin; i <= end_bin; ++i ) {
		if ( sim_map.find( i ) == sim_map.end() ) {
			sim_map[ i ] = 0.0;
		}
	}
	Real baseline( 0.0 );
	for ( auto const & dist_weight : sim_map ) {
		baseline += dist_weight.second;
	}
	if ( baseline == 0 ) {
		TR.Error << "Found distance distribution with zero area under the curve!" << std::endl;
		return std::map< Size, Real >();
	}
	for ( auto& dist_weight : sim_map ) {
		dist_weight.second /= baseline;
	}
	return sim_map;
}

/// @brief Given a set of electrons, a pose, and a residue of interest, find viable coords
utility::vector1< PseudoElectron >
EPRSpinLabel::get_coords(
	Size const & res,
	pose::Pose const & pose,
	utility::vector1< PseudoElectron > const & electrons,
	bool const & skip_clash_eval // = false
) {
	assert( pose.residue( res ).is_protein() );

	// This is what allows the coordinates to go from the local coordinate frame to the global frame
	numeric::HomogeneousTransform< Real > res_to_global(
		numeric::HomogeneousTransform< Real >(
		pose.residue( res ).xyz( "N"  ),
		pose.residue( res ).xyz( "C"  ),
		pose.residue( res ).xyz( "CA" )
		).inverse() );
	numeric::xyzVector< Real > cb = res_to_global.to_local_coordinate( cb_coord_ );

	// Bring into the global coordinate frame
	utility::vector1< PseudoElectron > v_e;
	std::transform( electrons.begin(), electrons.end(), std::back_inserter( v_e ),
		[&]( PseudoElectron e ){ return std::make_pair( res_to_global.to_local_coordinate( e.first ), e.second ); } );
	if ( skip_clash_eval ) {
		return v_e;
	}

	for ( Real forgive_factor = 0.80; forgive_factor >= 0.00; forgive_factor -= 0.05 ) {
		auto v_e_temp( v_e );
		for ( auto & e : v_e_temp ) {
			e.second *= get_weight( res, center_of_mass( e.first, cb ), pose, forgive_factor );
		}
		Size n_viable = 0;
		for ( auto const & e_temp : v_e_temp ) {
			if ( e_temp.second > cutoff_ ) {
				n_viable += 1;
			}
		}
		if ( n_viable > 0 ) {
			return v_e_temp;
		}
	}
	// If at this point output has not returned, just return the original set of coordinates
	return v_e;
}

/// @brief retrieve weight for given coordinate
Real
EPRSpinLabel::get_weight(
	core::Size const & source_res,
	numeric::xyzVector< core::Real > const & center_of_mass,
	pose::Pose const & pose,
	core::Real const & forgive_factor
) {
	// Using the neighborgraph since it is a more efficient way of checking for clashes
	Real weight( 1.0 );
	auto const & neighborgraph = pose.energies().twelveA_neighbor_graph();
	auto const & atom_vdw = ScoringManager::get_instance()->get_AtomVDW( chemical::CENTROID );
	for (
			auto edge_iter = neighborgraph.get_node( source_res )->const_edge_list_begin();
			edge_iter != neighborgraph.get_node( source_res )->const_edge_list_end();
			++edge_iter
			) {
		Size const & res = ( *edge_iter )->get_other_node( source_res )->get_node_index();
		if ( pose.residue( res ).is_virtual_residue() ) {
			continue;
		}
		for ( Size ii = 1; ii <= pose.residue( res ).natoms(); ++ii ) {
			Real clash_d_sq = ( !pose.is_fullatom() && ii == pose.residue( res ).atom_index( "CEN" ) )
				? pow( 2.4 + forgive_factor * atom_vdw.approximate_vdw_radius( pose.residue( res ).atom_type_index( ii ) ), 2 )
				: pow( 2.4 + forgive_factor * ( pose.residue( res ).atom_type( ii ).lj_radius() ), 2 );
			Real distance_sq = center_of_mass.distance_squared( pose.residue( res ).atom( ii ).xyz() );
			weight *= ( pose.residue( res ).atom_is_backbone( ii ) ) ?
				Real( bool( distance_sq > clash_d_sq ) ) : std::min( distance_sq / clash_d_sq, 1.0 );
			if ( weight <= cutoff_ ) {
				return 0.0;
			}
		}
	}
	return weight;
}

/// @brief Retrieve cuttof for weights (if the weight is less than this, it is set to zero)
Real
EPRSpinLabel::cutoff() const {
	return cutoff_;
}

/// @brief Allows a custom set of electrons to be read without superimposition
void
EPRSpinLabel::load_custom_electrons(
	std::map< Size, utility::vector1< PseudoElectron > > const & all_coords
) {
	for ( auto const & res : all_coords ) {
		custom_coords_[ res.first ] = res.second;
	}
}

/// @brief   Returns the center of mass between a given electron coordinate and its CB
numeric::xyzVector< Real >
EPRSpinLabel::center_of_mass(
	numeric::xyzVector< Real > const & electron,
	numeric::xyzVector< Real > const & cb
) const {
	return ( electron - cb ) * 0.875 + cb;
}

} // namespace epr_deer
} // namespace scoring
} // namespace core
