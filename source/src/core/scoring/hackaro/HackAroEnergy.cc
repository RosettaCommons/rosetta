// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/elec/HackAroEnergy.cc
/// @brief  Electrostatics energy method for aromatic side chain (stand-in for pi/pi interactions).
/// @author Rhiju Das


// Unit headers
#include <core/scoring/hackaro/HackAroEnergy.hh>
#include <core/scoring/hackaro/HackAroEnergyCreator.hh>

// Package headers
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/Stub.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <numeric/xyzMatrix.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>

// ObjexxFCL headers

// C++


/////////////////////////////////////////////////////////////////////////////////////////
//  Trying to get more "T-shaped" edge-to-face interactions of aromatic side chains.
//  This is very similar to Kira's fa_plane in rosetta++, though I haven't bothered
//   to parameterize as carefully. This is just a hacky functional form whose derivatives
//   are easy to calculate and smooth.
//  Perhaps the best approach would be to calculate, e.g., benzene-benzene interactions
//   with high level quantum calculations -- but note that you need really careful treatment
//   of electron correlations to get dispersion right -- still not quite feasible for all rigid body
//   orientations.
// (I also tried using elec_aro_aro, i.e., coulomb's law between positively charged aromatic
//  hydrogens and the negative partial charges on aromatic carbons -- it didn't quite give me
//  what I wanted).  Rhiju, Dec. 2009.
/////////////////////////////////////////////////////////////////////////////////////////

using core::Real;
typedef  numeric::xyzMatrix< Real > Matrix;

namespace core {
namespace scoring {
namespace hackaro {


/// @details This must return a fresh instance of the HackAroEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
HackAroEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new HackAroEnergy );
}

ScoreTypes
HackAroEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( hack_aro );
	return sts;
}

/// c-tor
HackAroEnergy::HackAroEnergy() :
	parent( methods::EnergyMethodCreatorOP( new HackAroEnergyCreator ) )
{}

//clone
methods::EnergyMethodOP
HackAroEnergy::clone() const
{
	return methods::EnergyMethodOP( new HackAroEnergy );
}


void
HackAroEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


void
HackAroEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
HackAroEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	// Aromatic means: PHE, TRP, TYR (not HIS, currently).
	if ( rsd1.is_aromatic() && rsd2.is_aromatic() ) {
		residue_pair_energy_aro_aro( rsd1, rsd2, emap );
	}
}


///////////////////////////////////////////////////////////////////////////////
Vector
HackAroEnergy::get_centroid( conformation::Residue const & rsd ) const
{
	Vector centroid( 0.0 );
	Size numatoms = 0;
	for ( Size i=rsd.first_sidechain_atom(); i<= rsd.nheavyatoms(); ++i ) {
		centroid += rsd.xyz(i);
		numatoms++;
	}
	if ( numatoms > 0 ) {
		centroid /= static_cast< Real >( numatoms );
	} else { //Yo, is this a glycine?
		debug_assert( rsd.aa() == chemical::aa_gly );
		centroid = rsd.xyz( "CA" );
	}

	return centroid;
}

/////////////////////////////////////////////////////////////////////////////////
kinematics::Stub
HackAroEnergy::get_base_coordinate_system( conformation::Residue const & rsd, Vector const & centroid ) const
{
	using namespace chemical;
	Size res_type = rsd.aa();

	Vector x,y,z;

	// Make an axis pointing from base centroid to Watson-Crick edge.
	std::string WC_atom;
	if ( res_type == aa_phe ) WC_atom = " CZ ";
	if ( res_type == aa_tyr ) WC_atom = " CZ ";
	if ( res_type == aa_trp ) WC_atom = " CZ2";

	Vector const WC_coord (rsd.xyz( WC_atom ) );
	x = WC_coord - centroid;
	x.normalize();

	// Make a perpendicular axis pointing from centroid towards
	// Hoogstein edge (e.g., major groove in a double helix).
	std::string H_atom;
	if ( res_type == aa_phe ) H_atom = " CD1";
	if ( res_type == aa_tyr ) H_atom = " CD1";
	if ( res_type == aa_trp ) H_atom = " CD1";

	Vector const H_coord (rsd.xyz( H_atom ) );
	y = H_coord - centroid; //not orthonormal yet...
	z = cross(x, y);
	z.normalize(); // Should poSize roughly 5' to 3' if in a double helix.

	y = cross(z, x);
	y.normalize(); //not necessary but doesn't hurt.

	//  std::cout << "WC : " << WC_coord << "   H : " << H_coord << "    centroid: " << centroid << std::endl;

	return kinematics::Stub( Matrix::cols( x, y, z ), centroid );
}

///////////////////////////////////////////////////////////////
Real
HackAroEnergy::get_aro_axis_score_ANGLE(
	Real const cos_theta,
	Real & deriv ) const {

	// Favor "T-shaped" arrangements with theta = 90 degrees.
	// - sin^4 ( theta )
	Real const value = -1 * ( 1 - cos_theta * cos_theta) * ( 1 - cos_theta * cos_theta);
	deriv = 4 * ( 1 - cos_theta * cos_theta ) * cos_theta;

	return value;
}

///////////////////////////////////////////////////////////////
Real
HackAroEnergy::get_aro_axis_score_DIST(
	Real const dist,
	Real & deriv ) const{

	using namespace core::scoring::constraints;

	static Real const lower_bound_( -2.0 );
	static Real const upper_bound_(  7.0 );
	static Real const bound_zone_ (  2.0 ); // bonus is flat below 5.0 Angstroms.
	static Real const well_depth_ (  1.0 ); // angular dependence goes negative.
	static core::scoring::func::FuncOP dist_func_( new core::scoring::func::FadeFunc( lower_bound_, upper_bound_, bound_zone_, well_depth_) );

	Real const value = dist_func_->func( dist );
	deriv = dist_func_->dfunc( dist );

	return value;
}

//////////////////////////////////////////////////////////////////////////////////
void
HackAroEnergy::residue_pair_energy_aro_aro(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	EnergyMap & emap
) const
{
	debug_assert( rsd1.is_aromatic() );
	debug_assert( rsd2.is_aromatic() );

	Vector centroid1 = get_centroid( rsd1 );
	kinematics::Stub stub1 = get_base_coordinate_system( rsd1, centroid1 );

	Vector centroid2 = get_centroid( rsd2 );
	kinematics::Stub stub2 = get_base_coordinate_system( rsd2, centroid2 );

	Real const cos_theta = dot( stub1.M.col_z(), stub2.M.col_z() );
	Real const cen_dist =  ( centroid1 - centroid2 ).length();

	Real dummy( 0.0 );
	Real const angle_score = get_aro_axis_score_ANGLE( cos_theta, dummy );
	Real const dist_score  = get_aro_axis_score_DIST( cen_dist, dummy );

	Real const total_score = angle_score * dist_score;

	emap[ hack_aro ] += total_score;
}

/////////////////////////////////////////////////////////////
Size
chi2_torsion_atom_index( conformation::Residue const & rsd )
{
	chemical::AtomIndices const & atom_indices = rsd.chi_atoms( 2 /*chi # 2 must be torsion that controls aromatic base*/ );
	return atom_indices[ 4 ]; /* C2' ... C1' ... first base atom ... chi1 torsion atom*/
}

/////////////////////////////////////////////
void
HackAroEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	using namespace chemical;

	Size const pos1( atom_id.rsd() );
	Size const i   ( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( pos1 ) );

	if ( !rsd1.is_aromatic() ) return;

	// Only works for proteins!!!!!!!!!!!!
	if ( i != chi2_torsion_atom_index( rsd1 ) ) return;

	int const pos1_map( domain_map( pos1 ) );
	bool const pos1_fixed( pos1_map != 0 );

	// cached energies object
	Energies const & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

	// loop over *all* nbrs of rsd1 (not just upper or lower)
	for ( utility::graph::Graph::EdgeListConstIter
			iru  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
			irue = energy_graph.get_node( pos1 )->const_edge_list_end();
			iru != irue; ++iru ) {
		Size const pos2( (*iru)->get_other_ind( pos1 ) );

		if ( pos1_fixed && pos1_map == domain_map( pos2 ) ) continue; // fixed wrt one another
		conformation::Residue const & rsd2( pose.residue( pos2 ) );
		if ( !rsd2.is_aromatic() ) continue;

		debug_assert( pos2 != pos1 );
		eval_atom_derivative_aro_aro( rsd1, rsd2, weights, F1, F2 );
	} // loop over nbrs of rsd1
}


//////////////////////////////////////////////////
void
HackAroEnergy::eval_atom_derivative_aro_aro(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{

	debug_assert( rsd1.is_aromatic() );
	debug_assert( rsd2.is_aromatic() );

	Vector centroid1 = get_centroid( rsd1 );
	kinematics::Stub stub1 = get_base_coordinate_system( rsd1, centroid1 );

	Vector centroid2 = get_centroid( rsd2 );
	kinematics::Stub stub2 = get_base_coordinate_system( rsd2, centroid2 );

	Real const cos_theta = dot( stub1.M.col_z(), stub2.M.col_z() );
	Vector const d_i_j =  centroid1 - centroid2;
	Real const cen_dist =  d_i_j.length();
	Vector const r_i_j = d_i_j / cen_dist;

	Real angle_deriv( 0.0 ), dist_deriv( 0.0 );
	Real const angle_score = get_aro_axis_score_ANGLE( cos_theta, angle_deriv );
	Real const dist_score  = get_aro_axis_score_DIST( cen_dist, dist_deriv );

	// checked my previous wrk for rna_stack_axis, et.c., in RNA_LowResolutionPotential.cc
	Matrix const & M_i( stub1.M );
	Vector const & z_i = M_i.col_z();
	Matrix const & M_j( stub2.M );
	Vector const & z_j = M_j.col_z();

	Vector const f2 =
		-1.0 * angle_score * dist_deriv * r_i_j;

	Vector const f1 = cross( f2, centroid2 ) + dist_score * angle_deriv * cross( z_i, z_j );

	F1 -= weights[ hack_aro ] * f1;
	F2 -= weights[ hack_aro ] * f2;
}


/// @brief HackAroEnergy is context independent; no context graphs required
void
HackAroEnergy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}

Distance
HackAroEnergy::atomic_interaction_cutoff() const{
	return 4.0; // pad by 4.0 because these are long sidechains?
}
core::Size
HackAroEnergy::version() const
{
	return 1; // Initial versioning
}


}
}
}
