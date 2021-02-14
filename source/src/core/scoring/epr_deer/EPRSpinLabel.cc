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
#include <core/scoring/epr_deer/util.hh>

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
#include <core/scoring/TenANeighborGraph.hh>
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
#include <utility/excn/Exceptions.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>
#include <numeric/HomogeneousTransform.hh>

#include <boost/algorithm/string.hpp>

// C++ headers
#include <stdlib.h>
#include <iostream>
#include <set>

namespace core {
namespace scoring {
namespace epr_deer {

/// @brief Tracer used for error messages
/// @details Global to avoid re-instantiating tracer with every new object
static basic::Tracer TR( "core.scoring.epr_deer.EPRSpinLabel" );

/// @brief Constructor
EPRSpinLabel::EPRSpinLabel() {
	init_vdw();
}

/// @brief Destructor
EPRSpinLabel::~EPRSpinLabel() {}

/// @brief Initialize object used to calculate centroid clashes
void
EPRSpinLabel::init_vdw() const {
	if ( !atom_vdw_ ) {
		atom_vdw_ = AtomVDWOP( new AtomVDW(
			ScoringManager::get_instance()->get_AtomVDW(
			chemical::CENTROID ) ) );
	}
}

/// @brief Read DB file for a given type of spin label
/// @param Name of file/SL type
/// @return PseudoSLs used for simulation of DEER distributions
utility::vector1< PseudoSL >
EPRSpinLabel::read_db_file(
	std::string const & name
) {

	utility::vector1< PseudoSL > output;

	// Pull the file
	std::string const fullname
		= path_ + ObjexxFCL::uppercased( name ) + ".txt";
	utility::io::izstream file_contents_db;
	basic::database::open( file_contents_db, fullname );

	// Iterate across each line
	std::string line;
	while ( !file_contents_db.eof() ) {
		getline( file_contents_db, line );
		utility::trim( line, " \t\n" );

		// This is in case the line contains problematic empty spaces
		utility::vector1< std::string > slraw, sl;
		boost::algorithm::split( slraw, line, boost::is_any_of( "\t " ) );
		for ( auto const & val : slraw ) {
			if ( val.size() > 0 ) {
				sl.push_back( val );
			}
		}

		// Pull XYZ and weight of each coord and convert to PseudoSL
		Real const x = std::stod( sl[ 1 ] );
		Real const y = std::stod( sl[ 2 ] );
		Real const z = std::stod( sl[ 3 ] );
		Real const w = std::stod( sl[ 4 ] );
		output.push_back( std::make_pair(
			numeric::xyzVector< Real >( x, y, z ), w ) );
	}
	return output;
}

/// @brief Operator to return nnn-const PseudoSL from specific residue
/// @param res: Residue info (number and spin label type)
/// @return Coords corresponding to residue
utility::vector1< PseudoSL > &
EPRSpinLabel::operator[](
	PairSizeString const & res
) {

	// Error in case these have not been labeled first
	if ( coords_.find( res ) == coords_.end() ) {
		throw CREATE_EXCEPTION( utility::excn::KeyError,
			"EPRSpinLabel::operator[]: index is out of range!" );
	} else {
		return coords_[ res ];
	}
}

/// @brief Return const PseudoSL from specific residue
/// @param res: Residue info (number and spin label type)
/// @return Coords corresponding to residue
utility::vector1< PseudoSL > const &
EPRSpinLabel::at( PairSizeString const & res ) const {

	// Error in case these have not been labeled first
	if ( coords_.find( res ) == coords_.end() ) {
		throw CREATE_EXCEPTION( utility::excn::KeyError,
			"EPRSpinLabel::operator[]: index is out of range!" );
	} else {
		return coords_.at( res );
	}
}

/// @brief Returns histogram between coordinate sets for residues
/// @param residues: Residues contributing to histogram
/// @param bins_per_a: Granularity of histogram (bins per angstrom)
/// @param mod: What to add to the X-axis
/// @param dist_ids: If a custom X-axis is used (default: empty)
/// @return Histogram with X- and Y-values being keys and values
/// @detail Note: the X-axis = bins_per_a * distance, rounded to an int
/// @detail Note: Equal labeling assumed for all sites
/// @detail This matters if residues.size() > 2
std::map< Size, Real >
EPRSpinLabel::histogram(
	utility::vector1< PairSizeString > const & residues,
	Size const & bins_per_a,
	int const & mod, // = 0
	Real const & stdev, // 1.0
	std::map< Size, Real > const & dist_ids // = {}
) {

	std::map< Size, Real > output;

	// Iterate through residue pairs
	for ( Size i = 1; i < residues.size(); ++i ) {
		for ( Size j = i + 1; j <= residues.size(); ++j ) {

			// Generate a histogram for this pair and add to output
			auto const distr = normalize( histogram( residues[ i ],
				residues[ j ], bins_per_a, mod, stdev, dist_ids ) );
			for ( auto const & dist : distr ) {
				add_to_map( output, dist.first, dist.second );
			}
		}
	}

	// Normalize so values add up to 1.0

	output = normalize( output );

	return output;
}

/// @brief Return histogram for pair of coordinate sets
/// @param res1_coords: Pair of PseudoSL coords for res1
/// @param res2_coords: Pair of PseudoSL coords for res2
/// @param bins_per_a: Bins per angstrom for distribution
/// @param mod: How much to shift X-axis
/// @param stdev: St deviation of gauss used to convolute pairwise dists
/// @param dist_ids: If a custom X-axis is used (default: empty)
/// @return Histogram with X- and Y-values being keys and values
/// @detail Note: the X-axis = bins_per_a * distance, rounded to an int
std::map< Size, Real >
EPRSpinLabel::histogram(
	utility::vector1< PseudoSL > const & res1_coords,
	utility::vector1< PseudoSL > const & res2_coords,
	Size const & bins_per_a,
	int const & mod, // = 0
	Real const & stdev, // = 1.0
	std::map< Size, Real > const & dist_ids // = {}
) {

	std::map< Size, Real > output;

	// How many st deviations away from the mean we go
	Size const STDEV_RANGE = 3;
	Real const width = STDEV_RANGE * stdev * bins_per_a;


	// Iterate over all pairs of coordinates
	for ( auto const & iter1 : res1_coords ) {
		for ( auto const & iter2 : res2_coords ) {

			// Get the dist (d) between and weight (w) of the two coords
			Real const w = iter1.second * iter2.second;
			Real const d = iter1.first.distance( iter2.first ) + mod;

			// If custom X-axis isn't used, define bin range and iterate
			if ( bins_per_a != 0 ) {
				Size const lo = std::max( 1.0, round( d * bins_per_a - width ) );
				Size const hi = round( d * bins_per_a + width );
				for ( Size bin = lo; bin <= hi; ++bin ) {
					Real const amp = gauss( bin / Real( bins_per_a ), d, stdev );
					add_to_map( output, bin, w * amp );
				}

				// Use custom X-axis if it is defined via dist_ids
			} else {

				for ( auto const & dist_id : dist_ids ) {
					Real const amp = gauss( dist_id.second, d, stdev );
					add_to_map( output, dist_id.first, w * amp );
				}
			}
		}
	}
	return normalize( output );
}

/// @brief Returns histogram between coordinate sets for residues
/// @param res1: Residue 1
/// @param res2: Residue 2
/// @param bins_per_a: Granularity of histogram (bins per angstrom)
/// @param mod: What to add to the X-axis
/// @param dist_ids: If a custom X-axis is used (default: empty)
/// @return Histogram with X- and Y-values being keys and values
/// @detail Note: the X-axis = bins_per_a * distance, rounded to an int
std::map< Size, Real >
EPRSpinLabel::histogram(
	PairSizeString const & res1,
	PairSizeString const & res2,
	Size const & bins_per_a,
	int const & mod, // = 0.0
	Real const & stdev, // = 1.0
	std::map< Size, Real > const & dist_ids // = {}
) {

	// Spot check to make sure we have coordinates for both residues
	for ( auto const & res : { res1, res2 } ) {
		if ( coords_.find( res ) == coords_.end() ) {
			TR.Error << "Residue " << res.first << " of SL type "
				<< res.second << " hasn't been labeled!" << std::endl;
		}
	}

	// Calculate histogram and do another spot check
	auto const output = histogram( coords_[ res1 ], coords_[ res2 ],
		bins_per_a, mod, stdev, dist_ids );
	if ( output.empty() ) {
		TR.Error << "Distance distribution with no distances measured: "
			<< res1.first << "\t" << res2.first << std::endl;
	}
	return output;
}

/// @brief Label a residue with a certain spin label
/// @param res: Residue index
/// @param label: SL type
/// @param pose: Pose used for superimposition
/// @param skip_clash_eval: Whether clash evaluation is skipped
/// @detail Clashes are skipped for custom coordinates for reasons
/// @detail  relating to the way they have been calculated
void
EPRSpinLabel::label(
	PairSizeString const & res_label,
	pose::Pose const & pose,
	bool const & skip_clash_eval // = false
) {

	// Aliases
	auto const & res = res_label.first;
	auto const & label = res_label.second;

	// Make sure computation isn't wasted by calling this multiple times
	if ( coords_.find( res_label ) == coords_.end() ) {

		// In the event that custom coords need to be used. This is because
		//  the rotamers were precomputed and do not need re-checked
		if ( res_label.second == "CUSTOM" ) {
			coords_[ res_label ] = calc_sl_for_res( res, pose,
				custom_coords_[ res ], true );

			// Otherwise proceed normally
		} else {

			// Check if the spin label type has even been read yet
			if ( deflt_coords_.find( label ) == deflt_coords_.end() ) {
				deflt_coords_[ label ] = read_db_file( label );
			}

			// Then go ahead and label
			coords_[ res_label ] = calc_sl_for_res( res, pose,
				deflt_coords_[ label ], skip_clash_eval );
		}
	}
}

/// @brief Normalize distribution so that the sum is equal to 1.0
/// @param sim_map: Simulated DEER distribution
/// @result Identical std::map except values add up to 1.0
std::map< Size, Real >
EPRSpinLabel::normalize(
	std::map< Size, Real > sim_map
) const {

	// Spot check that the contents aren't empty
	if ( sim_map.empty() ) {
		throw CREATE_EXCEPTION( utility::excn::KeyError,
			"Distribution has no elements!" );
	}

	// Add everything up
	Real total = 0.0;
	for ( auto const & x_y : sim_map ) {
		total += x_y.second;
	}

	// Another spot check to make sure we aren't dividing by zero...
	if ( total == 0.0 ) {
		throw CREATE_EXCEPTION( utility::excn::RangeError,
			"Total distribution area is zero!" );

		// ... or an infinitely large number...
	} else if ( std::isinf( std::abs( total ) ) ) {
		throw CREATE_EXCEPTION( utility::excn::RangeError,
			"Total distribution area is infinite!" );

	} else if ( std::isnan( std::abs( total ) ) ) {
		throw CREATE_EXCEPTION( utility::excn::RangeError,
			"Total distribution area is NaN!" );

		// Otherwise, proceed

		// Otherwise, proceed
	} else {

		// Divide by total and return
		for ( auto & x_y : sim_map ) {
			x_y.second /= total;
		}
		return sim_map;
	}
}

/// @brief Get positions of unpaired electrons at a residue
/// @param res: Residue number
/// @param pose: Pose for clash eval
/// @param sl_vec: Vector of PseudoSLs, which have positions of unpaired e
/// @param skip_clash_eval: Exactly what it suggested by the title
/// @param min_rad: Lowest radius for clash eval to check
/// @return Vector of PseudoSLs in local coordinate frame of residue
utility::vector1< PseudoSL >
EPRSpinLabel::calc_sl_for_res(
	Size const & res,
	pose::Pose const & pose,
	utility::vector1< PseudoSL > const & sl_vec,
	bool const & skip_clash_eval, // = false
	Real const & min_rad // = 0.0
) {

	// Spot check to make sure we are labeling an amino acid
	if ( !pose.residue( res ).is_protein() ) {
		throw CREATE_EXCEPTION( utility::excn::RangeError,
			"Attempting to spin label an non-amino acid residue!" );
	}

	// Create a local frame
	numeric::HomogeneousTransform< Real > const frame(
		numeric::HomogeneousTransform< Real >(
		pose.residue( res ).xyz( "N"  ),
		pose.residue( res ).xyz( "C"  ),
		pose.residue( res ).xyz( "CA" )
		).inverse() );
	numeric::xyzVector< Real > const cb = frame.to_local_coordinate( vrt_cb_ );

	// Bring into the global coordinate frame
	utility::vector1< PseudoSL > local_vec;
	std::transform( sl_vec.begin(), sl_vec.end(),
		std::back_inserter( local_vec ), [&]( PseudoSL e ){
			return std::make_pair( frame.to_local_coordinate( e.first ),
			e.second );
		} );

	// We can return this if skip_clash_eval is set to true
	if ( skip_clash_eval ) {
		return local_vec;
	}

	// Clash evaluation. This is the most time consuming step
	// We will start with a wide clash radius, and steadily decrease it
	//  if clashes prevent any PseudoSL from being placed.
	utility::vector1< PseudoSL > output;
	Real rad = 0.80;
	while ( output.size() == 0 && rad > min_rad ) {
		for ( auto const & e : local_vec ) {
			Real const w = weight( res, bulk( e.first, cb ), e.second, pose, rad );
			if ( w > cutoff_ ) {
				output.push_back( std::make_pair( e.first, w ) );
			}
		}
		rad -= 0.05;
	}

	// Return the SLs if it isn't empty
	if ( output.size() > 0 ) {
		return output;

		// If placing the SLs is impossible, just return the whole vector
	} else {
		return local_vec;
	}
}

/// @brief retrieve weight for given coordinate
/// @param res1: Residue over which the coordinate is being superimposed
/// @param clash_xyz: Coordinate used for clash calculation
/// @param w: Weight, passed by value since we need a new obj to modify
/// @param pose: Pose with all the residues we check for clash evaluation
/// @param vdw_rad: Radius of the clash_xyz atom to consider
/// @return Weight of PseudoSL at position given local environment of pose
Real
EPRSpinLabel::weight(
	Size const & res1,
	numeric::xyzVector< core::Real > const & clash_xyz,
	Real w,
	pose::Pose const & pose,
	Real const & rad
) {

	// Set initial check
	bool const & fa = pose.is_fullatom();

	// Neighborgraph provides most efficient way to check for clashes
	auto const & nbrs = pose.energies().tenA_neighbor_graph();

	// Iterate across all residues adjacent to res1 in neighbor graph
	auto it = nbrs.get_node( res1 )->const_edge_list_begin();
	for ( ; it != nbrs.get_node( res1 )->const_edge_list_end(); ++it ) {

		// The other residue
		Size const & res2 =
			( *it )->get_other_node( res1 )->get_node_index();
		if ( pose.residue( res2 ).is_virtual_residue() ) {
			continue;
		}

		// Iterate through atoms of this residue. Get dist & if backbone
		// Note that using distance squared halves compute time
		for ( Size ii = 1; ii <= pose.residue( res2 ).natoms(); ++ii ) {
			Real const d_sq = clash_xyz.distance_squared(
				pose.residue( res2 ).atom( ii ).xyz() );
			bool const & bb = pose.residue( res2 ).atom_is_backbone( ii );

			// If this is a centroid atom, get clash dist sq from AtomVDW
			Real vdw = 0.0;

			// NOTE: We need to have this funny-looking double loop
			// because Rosetta will crash if you ask for the "CEN" atom
			//  while in fullatom mode (fa)
			if ( fa ) {
				vdw = pose.residue( res2 ).atom_type( ii ).lj_radius();
			} else {
				if ( ii == pose.residue( res2 ).atom_index( "CEN" ) ) {
					vdw = atom_vdw_->approximate_vdw_radius(
						pose.residue( res2 ).atom_type_index( ii ) );
				} else {
					vdw = pose.residue( res2 ).atom_type( ii ).lj_radius();
				}
			}

			// Clash distance for comparison and modify weight
			Real const clash_d_sq = pow( 2.4 + rad * vdw, 2 );
			if ( bb && d_sq < clash_d_sq ) {
				return 0.0;
			} else {
				w *= std::min( 1.0, d_sq / clash_d_sq );
			}
		}
	}
	return w;
}

/// @brief Getter for cutoff for weights
/// @return Weight cutoff
Real
EPRSpinLabel::cutoff() const {
	return cutoff_;
}

/// @brief Allows a custom set of electrons to be read without superimposition
/// @param all_coords: Custom residue-specific PseudoSLs
void
EPRSpinLabel::load_custom_electrons(
	std::map< Size, utility::vector1< PseudoSL > > const & all_coords
) {
	for ( auto const & res : all_coords ) {
		custom_coords_[ res.first ] = res.second;
	}
}

/// @brief Goes through every residue in provided list and calculates
/// @param pose: Pose to label
/// @param residues: Residues that need to be labeled
/// @param skip_clash_eval: Whether clash evaluation should be skipped
void
EPRSpinLabel::label_everything(
	pose::Pose & pose,
	utility::vector1< PairSizeString > const & residues,
	bool const & skip_clash_eval
) {
	for ( auto const & res : residues ) {
		label( res, pose, skip_clash_eval );
	}
}

/// @brief Returns the center of mass between coordinate and CB
/// @param coord: Coordinate xyz
/// @param cb: CB coordinate
/// @return XYZ of bulk / center of mass of nitroxide ring
/// @detail See del Alamo et al 2020 Biophysical Journal for details on
///   why 0.875 was computed/chosen
numeric::xyzVector< Real >
EPRSpinLabel::bulk(
	numeric::xyzVector< Real > const & coord,
	numeric::xyzVector< Real > const & cb
) const {
	return ( coord - cb ) * 0.875 + cb;
}

} // namespace epr_deer
} // namespace scoring
} // namespace core
