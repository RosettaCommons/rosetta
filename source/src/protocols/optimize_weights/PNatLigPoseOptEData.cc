// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/optimize_weights/PNatLigPoseOptEData.cc
///
/// @brief
/// @author Ian W. Davis


#ifdef USEMPI
#include <mpi.h>
#endif

#include <protocols/optimize_weights/PNatLigPoseOptEData.hh>

#include <core/scoring/ScoreType.hh>

#include <ObjexxFCL/format.hh>

#include <utility/vector1.functions.hh>

#include <ostream>
#include <sstream>
#include <string>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.optimize_weights.OptEData" );

using namespace core;
using namespace scoring;
using namespace ObjexxFCL::format;

namespace protocols {
namespace optimize_weights {


PNatLigPoseOptEData::PNatLigPoseOptEData():
	PNatStructureOptEData(),
	// Arbitrary values for now:
	kT_(1.0),
	multiplier_(1.0) // don't use this value; use the value passed in to the component weights file (ronj)
{
}


PNatLigPoseOptEData::~PNatLigPoseOptEData() {}


Real
PNatLigPoseOptEData::do_score(
	std::ostream & ostr,
	Multivec const & component_weights,
	Multivec const & vars,
	Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const ,//num_ref_dofs,
	int const ,//num_total_dofs,
	EnergyMap const & fixed_terms,
	ScoreTypes const & ,//score_list,
	ScoreTypes const & fixed_score_list,
	bool const print
) const
{
	using namespace core::optimization;
	using namespace utility;
	//std::cout << "In get_score() ... " << natives_.size() << " natives, " << decoys_.size() << " decoys" << std::endl;
	//std::cout << "Weights:";
	//for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) std::cout << " " << vars[ ii ];
	//std::cout << std::endl;

	if ( decoys_.size() == 0 || natives_.size() == 0 ) return 0.0; // wtf?

	utility::vector1< Real > decoy_energies( decoys_.size(), 0.0 );
	utility::vector1< Real > native_energies( natives_.size(), 0.0 );
	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
		for ( Size jj = 1; jj <= natives_.size(); ++jj ) {
			native_energies[ jj ] += vars[ ii ]  * natives_[ jj ]->free_data()[ ii ];
		}
		for ( Size jj = 1; jj <= decoys_.size(); ++jj ) {
			decoy_energies[ jj ] += vars[ ii ] * decoys_[ jj ]->free_data()[ ii ];
		}
	}
	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		for ( Size jj = 1; jj <= natives_.size(); ++jj ) {
			native_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * natives_[ jj ]->fixed_data()[ ii ];
		}
		for ( Size jj = 1; jj <= decoys_.size(); ++jj ) {
			decoy_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * decoys_[ jj ]->fixed_data()[ ii ];
		}
	}

	Real const best_native_energy = min( native_energies );
	Real const best_decoy_energy = min( decoy_energies );
	//std::cout << "Best native E = " << best_native_energy << " , best decoy E = " << best_decoy_energy << std::endl;
	Real const best_energy =  best_native_energy < best_decoy_energy ? best_native_energy : best_decoy_energy;
	for ( Size ii = 1; ii <= natives_.size(); ++ii ) {
		native_energies[ ii ] -= best_energy;
	}
	for ( Size ii = 1; ii <= decoys_.size(); ++ii ) {
		decoy_energies[ ii ] -= best_energy;
	}

	Real numerator(0.0), partition(0.0);
	Multivec dpartition( vars.size(), 0.0 ), dnumerator( vars.size(), 0.0 );

	Real const neginv_kT = (-1.0 / kT_);
	for ( Size ii(1); ii <= natives_.size(); ++ii ) {

		// Limit the improbability of each native to 1 in a million.
		// This prevents numerator ~ 0, which causes NANs and INFs in the derivatives.
		// It also limits the "force" that any one structure can exert on the minimization.
		Real const exp_term = std::max( 1e-6, std::exp( neginv_kT * native_energies[ ii ] ) );
		//if( exp_term > 1 || exp_term < 0 || std::isinf(exp_term) || std::isnan(exp_term) ) std::cout << "[" << tag() << "] native exp_term = " << exp_term << std::endl;
		numerator += exp_term;
		partition += exp_term;

		for ( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
			// note for derivatives: d/dw( e^-(E*w+...) ) = -E * e^-(E*w+...)
			Real e_dof_deriv( neginv_kT * natives_[ ii ]->free_data()[ e_dof ] * exp_term );
			//if( std::isinf(e_dof_deriv) || std::isnan(e_dof_deriv) ) std::cout << "[" << tag() << "," << e_dof << "] native e_dof_deriv = " << e_dof_deriv << "; Eterm = " << natives_[ ii ]->free_data()[ e_dof ]<< std::endl;
			dnumerator[ e_dof ] += e_dof_deriv;
			dpartition[ e_dof ] += e_dof_deriv;
		}
	}
	for ( Size ii(1); ii <= decoys_.size(); ++ii ) {

		// Because partition >= numerator in all cases, there is no minimum value for this term:
		Real const exp_term( std::exp( neginv_kT * decoy_energies[ ii ] ) );
		//if( exp_term > 1 || exp_term < 0 || std::isinf(exp_term) || std::isnan(exp_term) ) std::cout << "[" << tag() << "] decoy exp_term = " << exp_term << std::endl;
		partition += exp_term;

		// partitions for energy derivatives
		for ( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
			// note for derivatives: d/dw( e^-(E*w+...) ) = -E * e^-(E*w+...)
			Real e_dof_deriv( neginv_kT * decoys_[ ii ]->free_data()[ e_dof ] * exp_term );
			//if( std::isinf(e_dof_deriv) || std::isnan(e_dof_deriv) ) std::cout << "[" << tag() << "," << e_dof << "] decoy e_dof_deriv = " << e_dof_deriv << "; Eterm = " << decoys_[ ii ]->free_data()[ e_dof ]<< std::endl;
			dpartition[ e_dof ] += e_dof_deriv;
		}
	}

	// -1 is here just to make lower scores better -- don't need kT (I think)
	// std::abs() is here in case log(N/P) == 0, so that we get +0 instead of -0 (which seems to break the minimizer?)
	Real const total_score = std::abs( -1.0 * multiplier_ * std::log( numerator / partition ) );
	//std::cout << " total score: " << total_score << std::endl;

	// If score is small enough, don't compute derivatives -- just say they're zero.
	// This may help protect us from weird minimizer run-away when "perfection" is attainable.
	if ( total_score >= 1e-2 ) {
		// accumulate to passed-in derivative sums -- excludes reference energies
		//std::cout << "vars (dvars): ";
		for ( Size dof(1); dof <= num_energy_dofs; ++dof ) {
			Real const dP_P = dpartition[ dof ] / partition;
			Real const dN_N = dnumerator[ dof ] / numerator;
			Real const dE_dvar = multiplier_ * (dP_P - dN_N);
			// This error *should* never occur, thanks to the minimum value for numerator (above).
#ifndef WIN32
			// std::isinf and std::isnan seem to give Visual Studio problems for some reason.
			if ( std::isinf(dE_dvar) || std::isnan(dE_dvar) ) std::cout << "[" << tag() << "," << dof << "] final deriv = " << dE_dvar << "; " << dpartition[ dof ] << "/" << partition << " - " << dnumerator[ dof ] << "/" << numerator << std::endl;
			else dE_dvars[ dof ] += component_weights[ type() ] * dE_dvar;
			//std::cout << " " << vars[ dof ] << "(" << dE_dvar << ")";
#endif
		}
	}

	if ( print ) {
		ostr << "PNatLigPose " << tag() << X(1)
			<< " num: " << F(7,3,numerator) << " part: " << F(7,3,partition)
			<< " p: " << F(7,5,numerator / partition)
			<< " -lnp: " << F(6,4,-1.0 * std::log( numerator / partition ))
			<< " -compwt_lnp: " << F(6, 4, component_weights[ type() ] * (-1.0 * std::log( numerator / partition )) ) << std::endl;
	}

	return component_weights[ type() ] * total_score;
}


OptEPositionDataType
PNatLigPoseOptEData::type() const
{
	return prob_native_ligand_pose;
}


} // namespace optimize_weights
} // namespace protocols
