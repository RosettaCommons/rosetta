// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/carbon_hbonds/AtomVDW.hh
/// @brief
/// @author Rhiju Das


// Unit Headers
#include <core/scoring/carbon_hbonds/CarbonHBondPotential.hh>
#include <core/scoring/rna/RNA_ScoringUtil.hh> // for get_fade_correction

// Package headers

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

// Utility headers
#include <utility/io/izstream.hh>

#include <numeric/conversions.hh>

#include <core/scoring/interpolation_util.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


namespace core {
namespace scoring {
namespace carbon_hbonds {

/// @details ctor, reads data file. Need to configure to allow alternate tables/atom_sets
CarbonHBondPotential::CarbonHBondPotential():
	bin_width_( 0 ),
	max_dist_( 0 ),
	aroC_scale_factor_( 2.0 ),
	num_carbon_donor_atoms_( 5 ),
	num_bins_( 45 ),
	rna_cos_theta_cutoff_( cos( numeric::conversions::radians( 60.0 ) ) /* 0.50 */),
	rna_cos_theta_fade_zone_( 0.2 ),
	rna_ch_o_bond_distance_( 2.35 )
{
	read_potential();
}

////////////////////////////////////////////////////////////////////////////////////////
//Why doesn't this helper function already exist in vector class?
Size
get_position_in_vector( utility::vector1< std::string> & vec, std::string const & element )
{
	Size count( 1 );
	for ( utility::vector1<std::string>::iterator iter = vec.begin(); iter < vec.end(); ++iter )	{
		if (*iter == element) return count;
		count++;
	}
	vec.push_back( element );
	return count;
}

////////////////////////////////////////////////////////////////////////////////////////
void
CarbonHBondPotential::read_potential()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::corrections::score;

	//Read in data file, and fill in private data.
	utility::io::izstream stream;

	// search in the local directory first
	stream.open( option[ ch_o_bond_potential ] );
	// then try database
	if ( !stream.good() ) {
		stream.close();
		basic::database::open( stream, option[ ch_o_bond_potential ] );
	}

	//This doesn't actually need to be hardwired -- if
	// I wasn't so lazy, I'd make it figure out the
	// dimensions based on the input file.
	for (Size i = 1; i <= num_carbon_donor_atoms_; i++ ){
		ObjexxFCL::FArray1D<Real> new_array; // =   new ObjexxFCL::FArray1D< Real> ;
		carbon_hbond_parameter_.push_back( new_array );
		carbon_hbond_parameter_[ i ].dimension( num_bins_ );
		carbon_hbond_parameter_[ i ] = 0.0;
	}

	if ( !stream.good() ) utility_exit_with_message( "Unable to open ch_o_bond_potential!" );

	std::string line, atom_name;
	Real bin_min, bin_max, value;
	bin_width_ =  0.0;

	while ( getline( stream, line ) ) {
		std::istringstream l( line );
		l >> atom_name >> bin_min >> bin_max >> value;
		//		atom_name = line.substr(0,4);

		Size const pos = get_position_in_vector( carbon_donors_, atom_name );

		if ( bin_width_ < 0.00001 ) bin_width_ = bin_max - bin_min;
		if (bin_max > max_dist_ ) max_dist_ = bin_max;

		Size const bin = static_cast<Size>( bin_max/bin_width_ + 0.001 );
		if (bin <= num_bins_ )		{
			carbon_hbond_parameter_[pos]( bin ) = value;
		}
	}

	// Force potential to be zero at boundary?
	//	for (Size i = 1; i <= num_carbon_donor_atoms_; i++ ) carbon_hbond_parameter_[i]( num_bins_ ) = 0.0;

	//Now read in atom type set, and lookup strings.
	core::chemical::AtomTypeSetCOP atom_type_set_cop =
		core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
	chemical::AtomTypeSet const & atom_type_set = *atom_type_set_cop;

	standard_atomtype_to_carbon_donor_index_.dimension( atom_type_set.n_atomtypes() );
	standard_atomtype_to_carbon_donor_index_ = 0;
	for (Size i = 1; i <= num_carbon_donor_atoms_; i++ ){
		Size const standard_atomtype_index = atom_type_set.atom_type_index( carbon_donors_[i] );
		//		std::cout << "HEY! " << i << " " << carbon_donors_[i] << " " << standard_atomtype_index << std::endl;
		standard_atomtype_to_carbon_donor_index_( standard_atomtype_index ) = i;
	}

	// Fudge for aromatic carbons.
	// Note this is totally based on intuition that DNA/RNA bases should make stronger carbon H-bonds than
	// what is inferred from proteins, due to activation by exocyclic atoms.
	// It would perhaps make sense to only apply the fudge for nucleic acids, and maybe then only
	// for appropriate positions.
	//
	Size const aro_index = standard_atomtype_to_carbon_donor_index_(	atom_type_set.atom_type_index( "aroC" ) );
	for (Size j = 1; j <= num_bins_; j++ ) {
		carbon_hbond_parameter_[ aro_index ]( j )  *= aroC_scale_factor_;
	}


	// Precalculate derivative. Would it make sense to do this and the interpolation via splines?
	// And is the potential smooth enough?
	//
	//	carbon_hbond_deriv_.dimension( carbon_hbond_parameter_.size1(), carbon_hbond_parameter_.size2() );
	//	carbon_hbond_deriv_ = 0.0;
	//	for (Size i = 1; i <= carbon_hbond_parameter_.size1(); i++ ) {
	//		for (Size j = 1; j < carbon_hbond_parameter_.size2(); j++ ) {
	//			carbon_hbond_deriv_(i,j) = ( carbon_hbond_parameter_(i,j+1) - carbon_hbond_parameter_(i,j) )/bin_width_;
	//		}
	//	}
}


//////////////////////////////////////////////////////////////////
Real
CarbonHBondPotential::get_potential_RNA(
	Vector const & r_vec,
	Vector const & z_i /*unit vector pointing from donor to its connected hydrogen*/,
	bool const & update_deriv,
	Vector & deriv_vector
) const
{

	Real const r = r_vec.length();
	Real const z = dot( z_i, r_vec );
	Real const cos_kappa = z / r;

  //Orientation dependence
	Real angle_fade_value( 1.0 ), angle_fade_deriv( 0.0 );
	core::scoring::rna::get_fade_correction( cos_kappa, rna_cos_theta_cutoff_, 1.5, rna_cos_theta_fade_zone_, angle_fade_value, angle_fade_deriv );
  Real score  = angle_fade_value;

  //Distance dependence
	Real const d_H_A = r_vec.length();

	Real const scaled_dist_H_A = d_H_A/ rna_ch_o_bond_distance_ ;
	Real const lj_style_distance_dependence = ( 1.0/std::pow(scaled_dist_H_A, 8) - 2.0/std::pow(scaled_dist_H_A, 4) );
	score *= lj_style_distance_dependence;

	Real distance_fade_value( 1.0 ), distance_fade_deriv( 0.0 );
	core::scoring::rna::get_fade_correction( d_H_A, 0, max_dist_, 0.5, distance_fade_value, distance_fade_deriv );
	score *= distance_fade_value;

	if ( update_deriv ) {
		Real const dEdcoskappa( angle_fade_deriv * lj_style_distance_dependence);

		Real dEdr = angle_fade_value * distance_fade_value * ( -8.0/rna_ch_o_bond_distance_ ) * ( 1.0/std::pow( scaled_dist_H_A, 9) - 1.0/std::pow( scaled_dist_H_A, 5) );
		dEdr += angle_fade_value * distance_fade_deriv * lj_style_distance_dependence;

		//Coordinate system ... z points along donor-hydrogen vector. x is orthogonal, and points towards acceptor.
		Vector x_i = r_vec - z_i * dot(r_vec,z_i);
		Real const x = x_i.length();
		x_i = x_i.normalize();

		Real const dEdx = dEdr * (x/r) + ( - x * z ) * dEdcoskappa / (r * r * r);
		Real const dEdz = dEdr * (z/r) + ( r*r - z*z) * dEdcoskappa / (r * r * r);

		deriv_vector = dEdx * x_i + dEdz * z_i;

	}

  return score;

}

/////////////////////////////////////////////////////////////////////////////////////////
Real
CarbonHBondPotential::get_potential(
	Size const & atom_type_index,
	Vector const & H_A_vector,
	Vector const & D_H_vector,
	Vector const & B_A_vector,
	bool calculate_deriv,
	Vector & deriv_vector
) const
{
	deriv_vector = 0.;

	Real dist_H_A = H_A_vector.length();
	Vector  H_A_norm_vector ( H_A_vector ); H_A_norm_vector.normalize();
	Vector  D_H_norm_vector ( D_H_vector ); D_H_norm_vector.normalize();
	Vector  B_A_norm_vector ( B_A_vector ); B_A_norm_vector.normalize();

	Real neg_cos_theta =  dot ( H_A_norm_vector, D_H_norm_vector );
	Real neg_cos_psi   = -dot ( H_A_norm_vector, B_A_norm_vector );

	Size const carbon_donor_index( standard_atomtype_to_carbon_donor_index_( atom_type_index ) );
	if ( carbon_donor_index == 0 ) return 0.0;

	Real dist_value( 0.0 );
	Real dEdHA( 0.0 );
	// We need to interpolate
	interpolate_value_and_deriv(
		carbon_hbond_parameter_[ carbon_donor_index ], bin_width_,
		dist_H_A, dist_value, dEdHA );

	Real theta_fade_value( 1.0 ), theta_fade_deriv( 0.0 );
	Real psi_fade_value( 1.0 ), psi_fade_deriv( 0.0 );
	core::scoring::rna::get_fade_correction( neg_cos_theta, 0.5, 1.5, 0.2, theta_fade_value, theta_fade_deriv );
	core::scoring::rna::get_fade_correction( neg_cos_psi,   0.4, 1.5, 0.2, psi_fade_value,   psi_fade_deriv );

	Real const value( dist_value * theta_fade_value * psi_fade_value );

	if ( ! calculate_deriv ) return value;

	deriv_vector = -dEdHA * theta_fade_value * psi_fade_value * H_A_norm_vector;

	Real neg_sin_theta = -(cross(D_H_norm_vector, H_A_norm_vector)).length();
	Vector theta_i = -cross( cross(D_H_norm_vector, H_A_norm_vector), H_A_norm_vector);
	theta_i.normalize();
	Vector theta_fade_gradient = (1./ dist_H_A) * theta_fade_deriv * neg_sin_theta * theta_i;
	deriv_vector += dist_value * psi_fade_value * theta_fade_gradient;

	Real neg_sin_psi = -(cross(B_A_norm_vector, H_A_norm_vector)).length();
	Vector psi_i = cross( cross(B_A_norm_vector, H_A_norm_vector), H_A_norm_vector);
	psi_i.normalize();
	Vector psi_fade_gradient = (1. / dist_H_A) * psi_fade_deriv * neg_sin_psi * psi_i;
	deriv_vector += dist_value * theta_fade_value * psi_fade_gradient;

	//Vector H_A_vector_incr = H_A_vector + 0.01 * psi_i;
	//Vector H_A_norm_vector_incr(H_A_vector_incr); H_A_norm_vector_incr.normalize();
	//Real neg_cos_psi_incr =  -dot ( H_A_norm_vector_incr, B_A_norm_vector );
	//Real psi_fade_value_incr( 1.0 ), psi_fade_deriv_incr( 0.0 );
	//core::scoring::rna::get_fade_correction( neg_cos_psi_incr, 0.3, 1.5, 0.3, psi_fade_value_incr, psi_fade_deriv_incr );
	//std::cout
	//<< " incr = "            << std::setw(8) << std::fixed << std::setprecision(5) << (neg_cos_psi_incr - neg_cos_psi)
	//<< " ∆E(psi) = "         << std::setw(8) << std::fixed << std::setprecision(5) << (psi_fade_value_incr - psi_fade_value)
	//<< " psi = "       << std::setw(8) << std::fixed << std::setprecision(5) << psi_i
	//<< " -cos(psi) = "       << std::setw(8) << std::fixed << std::setprecision(5) << neg_cos_psi
	//<< " -sin(psi) = "       << std::setw(8) << std::fixed << std::setprecision(5) << neg_sin_psi
	//<< " psi_fade_deriv1 = " << std::setw(8) << std::fixed << std::setprecision(5) << -(0.01 / dist_H_A) * psi_fade_deriv * neg_sin_psi
	//<< " psi_fade_deriv2 = " << std::setw(8) << std::fixed << std::setprecision(5) << -(0.01) * psi_fade_deriv * neg_sin_psi
	//<< " psi_fade_deriv3 = " << std::setw(8) << std::fixed << std::setprecision(5) << -(0.01) * psi_fade_deriv
	//<< std::endl;

	//∇{f(r)*g(theta)} = g(theta)∇f(r) + f(r)∇g(theta)
	//∇g(theta) = 1/r * ∂g(theta)/∂theta * theta_i
	//∂g(theta)/∂theta = -sin(theta) * theta_fade_deriv
	//theta_i = normalize(D_H_i X H_A_i) X H_A_i

	//	Real comparison_energy( 0.0 );
	//	Size const bin = static_cast<Size> ( dist_H_A / bin_width_ ) + 1;
	//	if ( bin <= num_bins_ )	comparison_energy = carbon_hbond_parameter_[carbon_donor_index]( bin );
	//	if ( value != 0.0 ){
	//		std::cout << "COMPARE: " << value << " " << comparison_energy << std::endl;
	//	}
	//	value = comparison_energy;

	return value;
}

/// @details Locally allocated deriv_vector is passed by reference into the more general implementation
/// of the get_potential function; moreover, the calculate-derivative flag is set to false in that call.
/// This is the "evaluate only the score" interface to the chbond potential.
Real
CarbonHBondPotential::get_potential(
	Size const & atom_type_index,
	Vector const & H_A_vector, //vector of hydrogen to acceptor
	Vector const & D_H_vector, //vector of donor hv atom to hydrogen
	Vector const & B_A_vector //vector of acceptor's base atom to acceptor
) const
{
	Vector deriv_vector;
	return get_potential(atom_type_index, H_A_vector, D_H_vector, B_A_vector, false, deriv_vector);
}


} // namespace carbon_hbonds
} // namespace scoring
} // namespace core

