// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/Fa_MbenvEnergy.cc
/// @author Patrick Barth


// Unit headers
#include <core/energy_methods/Fa_MbenvEnergy.hh>
#include <core/energy_methods/Fa_MbenvEnergyCreator.hh>

// Package headers
#include <core/scoring/Membrane_FAPotential.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh> //pba
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/memb_etable/MembEtable.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh> //pba
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

//#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/FArray1.fwd.hh>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.energy_methods.Fa_MbenvEnergy" );

namespace core {
namespace energy_methods {



/// that created object reads in the three database files: p_aa, p_aa_pp, and p_aa_n.  That object is returned and then stored as
/// a private member variable here.
/// @details This must return a fresh instance of the Fa_MbsolvEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
Fa_MbenvEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< Fa_MbenvEnergy >(*( core::scoring::ScoringManager::get_instance()->memb_etable( options.etable_type() ).lock() ) );
}

core::scoring::ScoreTypes
Fa_MbenvEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( fa_mbenv );
	return sts;
}


Fa_MbenvEnergy::Fa_MbenvEnergy( core::scoring::etable::MembEtable const & memb_etable_in ):
	parent( utility::pointer::make_shared< Fa_MbenvEnergyCreator >() ),
	//memb_etable_(memb_etable_in),
	lk_dgrefce_(memb_etable_in.lk_dgrefce()),
	memb_lk_dgrefce_(memb_etable_in.memb_lk_dgrefce()),
	potential_( core::scoring::ScoringManager::get_instance()->get_Membrane_FAPotential() )
{}


core::scoring::methods::EnergyMethodOP
Fa_MbenvEnergy::clone() const {
	return utility::pointer::make_shared< Fa_MbenvEnergy >( *this );
}


/////////////////////////////////////////////////////////////////////////////
// methods for ContextDependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////


void
Fa_MbenvEnergy::setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	potential_.compute_fa_projection( pose );
}

void
Fa_MbenvEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	core::scoring::EnergyMap & emap ) const {

	for ( Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i ) {
		// TP3 waters should also not have fa_mbenv in hydrate/SPaDES protocol
		if ( basic::options::option[ basic::options::OptionKeys::score::water_hybrid_sf ] ) {
			if ( rsd.name() != "TP3" ) {
				emap[ core::scoring::fa_mbenv ] += eval_fa_mbenv( rsd.atom(i), core::scoring::Membrane_FAEmbed_from_pose( pose ).fa_proj(rsd.seqpos(),i));
			}
		} else { // default behavior
			emap[ core::scoring::fa_mbenv ] += eval_fa_mbenv( rsd.atom(i), core::scoring::Membrane_FAEmbed_from_pose( pose ).fa_proj(rsd.seqpos(),i));
		}

	}
}

////////////////////////////////////////////////
Real
Fa_MbenvEnergy::eval_fa_mbenv(
	conformation::Atom const & atom1,
	Real const & f1 ) const
{

	Real temp_score( 0.0 );
	temp_score = (1 - f1) * (memb_lk_dgrefce_(atom1.type()) - lk_dgrefce_(atom1.type()));
	return temp_score;
}

/////////////////////////////////////////////////////////////////////////////
// derivatives
/////////////////////////////////////////////////////////////////////////////
void
Fa_MbenvEnergy::setup_for_derivatives(
	pose::Pose & pose,
	core::scoring::ScoreFunction const & scfxn
) const
{
	potential_.compute_fa_projection( pose );
	fa_mbenv_weight_ = scfxn.weights()[ core::scoring::fa_mbenv ];
}

//////////////////////////////////////////////////////////////////////////////////////
void
Fa_MbenvEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & /*domain_map*/,
	core::scoring::ScoreFunction const &,// sfxn,
	core::scoring::EnergyMap const & /*weights*/,
	Vector & F1,
	Vector & F2
) const
{

	Size const i( atom_id.rsd() );
	Size const m( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( i ) );

	if ( m > rsd1.nheavyatoms() ) return;

	Vector const heavy_atom_i( rsd1.xyz( m ) );

	Real cp_weight = 1.0;

	Vector const center(MembraneEmbed_from_pose( pose ).center());

	Vector f1( 0.0 ), f2( 0.0 );

	Real const deriv = memb_lk_dgrefce_(rsd1.atom(m).type()) - lk_dgrefce_(rsd1.atom(m).type());
	Real dE_dZ_over_r = fa_mbenv_weight_ * deriv * core::scoring::Membrane_FAEmbed_from_pose( pose ).fa_proj_deriv(rsd1.seqpos(),m);

	Vector const d_ij = core::scoring::Membrane_FAEmbed_from_pose( pose ).fa_proj_coord(rsd1.seqpos(),m) - heavy_atom_i;
	Real const d_ij_norm = d_ij.length();
	if ( d_ij_norm == Real(0.0) ) return;

	Real const invd = 1.0 / d_ij_norm;
	f2 = d_ij * invd;
	f1 = core::scoring::Membrane_FAEmbed_from_pose( pose ).fa_proj_coord(rsd1.seqpos(),m).cross(heavy_atom_i);
	f1 *= invd;

	if ( dE_dZ_over_r != 0.0 ) {
		F1 += dE_dZ_over_r * cp_weight * f1;
		F2 += dE_dZ_over_r * cp_weight * f2;
	}
}

////////////////////////////////////////////////
void
Fa_MbenvEnergy::finalize_total_energy(
	pose::Pose & /*pose*/,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	//std::cout << "BEFORE emap[ core::scoring::fa_mbenv ] " << emap[ core::scoring::fa_mbenv ] << std::endl;
	emap[ core::scoring::fa_mbenv ] += 0.0;
	//std::cout << "AFTER emap[ core::scoring::fa_mbenv ] " << emap[ core::scoring::fa_mbenv ] << std::endl;
}

/// @brief Fa_MbenvEnergy is context independent; indicates that no context graphs are required
void
Fa_MbenvEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}

/// @details Pose must already contain a cenlist object or this method will fail.
core::scoring::Membrane_FAEmbed const &
Fa_MbenvEnergy::Membrane_FAEmbed_from_pose( pose::Pose const & pose ) const
{
	return *( utility::pointer::static_pointer_cast< core::scoring::Membrane_FAEmbed const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED ) ));
}

core::scoring::MembraneEmbed const &
Fa_MbenvEnergy::MembraneEmbed_from_pose( pose::Pose const & pose ) const
{
	debug_assert( pose.data().has( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED ) );
	return *( utility::pointer::static_pointer_cast< core::scoring::MembraneEmbed const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED ) ));
}

core::scoring::MembraneTopology const &
Fa_MbenvEnergy::MembraneTopology_from_pose( pose::Pose const & pose ) const
{
	return *( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
}
core::Size
Fa_MbenvEnergy::version() const
{
	return 1; // Initial versioning
}


} // scoring
} // core

