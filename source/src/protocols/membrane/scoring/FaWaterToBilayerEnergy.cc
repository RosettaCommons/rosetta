// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/membrane/scoring/FaWaterToBilayerEnergy.cc
/// @brief Implicit Lipid Membrane Model Water-to-Bilayer transfer energy (one-body)
/// @author  Rebecca Alford (ralford3@jhu.edu)

#ifndef INCLUDED_protocols_membrane_scoring_FaWaterToBilayerEnergy_cc
#define INCLUDED_protocols_membrane_scoring_FaWaterToBilayerEnergy_cc

// Unit headers
#include <protocols/membrane/scoring/FaWaterToBilayerEnergy.hh>
#include <protocols/membrane/scoring/FaWaterToBilayerEnergyCreator.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>

// Project Headers
#include <protocols/membrane/scoring/MEnvAtomParams.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/MembraneGeometry.hh>
#include <core/conformation/Atom.hh>


// Package headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/scoring/EnergyMap.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

// Utility Headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// C++ Headers

static basic::Tracer TR( "protocols.membrane.scoring.FaWaterToBilayerEnergy" );

namespace protocols {
namespace membrane {
namespace scoring {

/// @brief Return a fresh instance of the energy method
core::scoring::methods::EnergyMethodOP
FaWaterToBilayerEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< FaWaterToBilayerEnergy >( options );
}

/// @brief Return Score Types Required for Method
core::scoring::ScoreTypes
FaWaterToBilayerEnergyCreator::score_types_for_method() const {
	core::scoring::ScoreTypes sts;
	sts.push_back( core::scoring::fa_water_to_bilayer );
	return sts;
}

/// @brief Construct Energy Method from Etable
FaWaterToBilayerEnergy::FaWaterToBilayerEnergy( core::scoring::methods::EnergyMethodOptions const & options ):
	parent( utility::pointer::make_shared< FaWaterToBilayerEnergyCreator >() ),
	memb_lk_dgrefce_(),
	water_lk_dgrefce_(),
	atypes_list_(),
	use_fleming_de_( options.use_fleming_de() )
{
	std::string dbfile( "membrane/memb_fa_params_2019.txt" );
	using namespace basic;
	using namespace core;

	// #ifdef USEMPI
	//  utility_exit_with_message("fa_water_to_bilayer is not yet threadsafe for MPI mode!");
	// #endif

	// #ifdef MULTI_THREADED
	//  utility_exit_with_message("fa_water_to_bilayer is not yet threadsafe for MULTI_THREADED mode!");
	// #endif

	utility::io::izstream infile;
	TR << "Reading fa_water_to_bilayer parameters from the database" << std::endl;
	database::open( infile, dbfile );
	if ( !infile.good() ) {
		utility_exit_with_message( "Unable to open database file containing reference energies" );
	}

	std::string line;
	getline( infile, line );
	while ( getline( infile, line ) ) {
		utility::trim(line, "\t\n" );
		if ( line.empty() > 0 ) continue;
		if ( line.find("#",0) == 0 ) continue;

		std::istringstream l( line );
		std::string atype, memb_lk_dgrefce, water_lk_dgrefce;
		l >> atype >> water_lk_dgrefce >> memb_lk_dgrefce;

		if ( l.fail() ) {
			utility_exit_with_message( "bad line: " + line );
		}

		// Add data to maps
		atypes_list_.push_back( atype );
		memb_lk_dgrefce_.push_back( utility::from_string( memb_lk_dgrefce, core::Real(0.0) ) );
		water_lk_dgrefce_.push_back( utility::from_string( water_lk_dgrefce, core::Real(0.0) ) );
	}
}

/// @brief Clone Energy Method
core::scoring::methods::EnergyMethodOP
FaWaterToBilayerEnergy::clone() const {
	return core::scoring::methods::EnergyMethodOP( new FaWaterToBilayerEnergy( *this ) );
}


/// @brief Compute Per-Residue Energies
void
FaWaterToBilayerEnergy::residue_energy(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	core::scoring::EnergyMap & emap
) const {

	using namespace core;
	using namespace core::scoring;
	//geting rid of fixed values for each residue.
	if ( rsd.name3() == "MEM" || !rsd.is_protein() ) return;
	for ( core::Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i ) {
		MEnvAtomParamsCOP menv_params = get_menv_params_for_residue( pose, rsd, i );
		emap[ core::scoring::fa_water_to_bilayer ] += eval_fa_wtbe( *menv_params );
	}

}

/// @brief Fianlzie Total Per-Residue Energies
void
FaWaterToBilayerEnergy::finalize_total_energy(
	core::pose::Pose &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const {
	emap[ core::scoring::fa_water_to_bilayer ] += 0.0;
}

/// @brief Setup for Computing Derivatives
void
FaWaterToBilayerEnergy::setup_for_derivatives(
	core::pose::Pose &,
	core::scoring::ScoreFunction const & scfxn
) const {
	fa_wtbe_weight_ = scfxn.weights()[ core::scoring::fa_water_to_bilayer ];
}

/// @brief Evaluate Per-Atom Derivatives
void
FaWaterToBilayerEnergy::eval_atom_derivative(
	core::id::AtomID const & atom_id,
	core::pose::Pose const & pose,
	core::kinematics::DomainMap const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const &,
	core::Vector & F1,
	core::Vector & F2
) const {

	using namespace core;
	using namespace core::scoring;

	// Grab residue and atom numbers
	core::Size const i( atom_id.rsd() );
	core::Size const m( atom_id.atomno() );

	// Get the actual residue and atom
	conformation::Residue const & rsd1( pose.residue( i ) );
	if ( m > rsd1.nheavyatoms() ) return;
	Vector const heavy_atom_i( rsd1.xyz( m ) );

	// Get the relevant atom type parameters
	MEnvAtomParamsCOP p( get_menv_params_for_residue( pose, rsd1, m ) );

	// Set count pair weight
	Real cp_weight = 1.0;

	Real const deriv = p->dGfreeB() - p->dGfreeW();


	F1 += fa_wtbe_weight_ * cp_weight * deriv * p->f1();
	F2 += fa_wtbe_weight_ * cp_weight * deriv * p->f2();

}

/// @brief Fa_MbenvEnergy is context independent
void
FaWaterToBilayerEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}

/// @brief Setup Method for initial scoring
void
FaWaterToBilayerEnergy::setup_for_scoring(
	core::pose::Pose &,
	core::scoring::ScoreFunction const &
) const {}

/// @brief Evaluate Per-Atom Env term
core::Real
FaWaterToBilayerEnergy::eval_fa_wtbe(
	MEnvAtomParams const & p
) const {

	using namespace core;

	Real score( 0.0 );
	score = ( 1 - p.hydration() ) * ( p.dGfreeB() - p.dGfreeW() );
	return score;
}

MEnvAtomParamsCOP
FaWaterToBilayerEnergy::get_menv_params_for_residue(
	core::pose::Pose const & pose,
	core::conformation::Residue const & rsd,
	core::Size atomno
) const {

	using namespace core::pose;
	using namespace core::conformation::membrane;

	// Get the water-to-bilayer free energy components
	core::Size index( get_atype_index( rsd.atom_type( atomno ).name() ) );
	core::Real dGfreeW = water_lk_dgrefce_[ index ];
	core::Real dGfreeB = memb_lk_dgrefce_[ index ];

	if ( use_fleming_de_ && ( rsd.atom_type( atomno ).name()=="COO" || rsd.atom_type( atomno ).name()=="OOC" ) ) {
		if ( rsd.atom_type( atomno ).name()=="COO" ) {
			dGfreeB = -0.5225;
		} else {
			dGfreeB = -9.2649;
		}
		//TR << rsd.atom_type( atomno ).name() << dGfreeB << std::endl;
	}

	core::conformation::Conformation const & conf( pose.conformation() );
	core::conformation::membrane::MembraneGeometryCOP mp_geometry( conf.membrane_info()->membrane_geometry() );

	if ( !conf.membrane_info()->use_franklin() ) {
		TR.Warning << "Score term optimized with franklin transition function, but using IMM1 membrane transiton funtion!" << std::endl;
		TR.Warning << "Recommended to run with -mp:restore_lazaridis_imm_behavior false if using franklin2019.wts" << std::endl;
	}

	core::Real hydration = mp_geometry->f_transition( conf, rsd.seqpos(), atomno );
	core::Vector f1 = mp_geometry->f_transition_f1( conf, rsd.seqpos(), atomno );
	core::Vector f2 = mp_geometry->f_transition_f2( conf, rsd.seqpos(), atomno );

	MEnvAtomParamsCOP menv_parameter = MEnvAtomParamsCOP( new MEnvAtomParams( rsd.atom_type(atomno).name(), dGfreeW, dGfreeB, hydration, f1, f2 ) );

	return menv_parameter;

}

core::Size
FaWaterToBilayerEnergy::get_atype_index( std::string atype_name ) const {
	core::Size index(0);
	for ( core::Size ii = 1; ii <= atypes_list_.size(); ++ii ) {
		if ( atypes_list_[ii] == atype_name ) index = ii;
	}
	debug_assert( index != 0 );
	return index;
}

/// @brief Versioning
core::Size
FaWaterToBilayerEnergy::version() const {
	return 2;
}

} // membrane
} // scoring
} // core

#endif // INCLUDED_core_scoring_methods_FaWaterToBilayerEnergy_cc
