// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/Fa_MbenvEnergy.cc
/// @author Patrick Barth


// Unit headers
#include <core/scoring/methods/Fa_MbenvEnergy.hh>
#include <core/scoring/methods/Fa_MbenvEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/Membrane_FAPotential.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh> //pba
#include <core/scoring/methods/EnergyMethodOptions.hh> //pba
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/memb_etable/MembEtable.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh> //pba

//#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/FArray1.fwd.hh>

#include <utility/vector1.hh>
#include <basic/Tracer.hh> 

static thread_local basic::Tracer TR( "core.scoring.methods.Fa_MbEnvEnergy" );

namespace core {
namespace scoring {
namespace methods {


/// that created object reads in the three database files: p_aa, p_aa_pp, and p_aa_n.  That object is returned and then stored as
/// a private member variable here.
/// @details This must return a fresh instance of the Fa_MbsolvEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
Fa_MbenvEnergyCreator::create_energy_method(
  methods::EnergyMethodOptions const & options
) const {
  return methods::EnergyMethodOP( new Fa_MbenvEnergy(*( ScoringManager::get_instance()->memb_etable( options.etable_type() ).lock() ) ) );
}

ScoreTypes
Fa_MbenvEnergyCreator::score_types_for_method() const {
  ScoreTypes sts;
  sts.push_back( fa_mbenv );
  return sts;
}


Fa_MbenvEnergy::Fa_MbenvEnergy( etable::MembEtable const & memb_etable_in):
  parent( methods::EnergyMethodCreatorOP( new Fa_MbenvEnergyCreator ) ),
  memb_etable_(memb_etable_in),
  lk_dgrefce_(memb_etable_in.lk_dgrefce()),
  memb_lk_dgrefce_(memb_etable_in.memb_lk_dgrefce()),
	potential_( ScoringManager::get_instance()->get_Membrane_FAPotential() )
{}


EnergyMethodOP
Fa_MbenvEnergy::clone() const {
	return EnergyMethodOP( new Fa_MbenvEnergy( *this ) );
}


/////////////////////////////////////////////////////////////////////////////
// methods for ContextDependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

///
void
Fa_MbenvEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
  potential_.compute_fa_projection( pose );
}

void
Fa_MbenvEnergy::residue_energy(
  conformation::Residue const & rsd,
  pose::Pose const & pose,
  EnergyMap & emap ) const {

  for ( Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i ) {
      emap[ fa_mbenv ] += eval_fa_mbenv( rsd.atom(i), Membrane_FAEmbed_from_pose( pose ).fa_proj(rsd.seqpos(),i));
	  
  }
}

////////////////////////////////////////////////
Real
Fa_MbenvEnergy::eval_fa_mbenv(
  conformation::Atom const & atom1,
  Real const & f1 ) const
{

  Real temp_score( 0.0 );
  //Make this an input option for efficiency
  //bool const eval_deriv( true );

  // l1 and l2 are FArray LINEAR INDICES for fast lookup:
  // [ l1 ] == (disbin  ,attype2,attype1)
  // [ l2 ] == (disbin+1,attype2,attype1)

  temp_score = (1 - f1) * (memb_lk_dgrefce_(atom1.type()) - lk_dgrefce_(atom1.type()));
  return temp_score;

}

/////////////////////////////////////////////////////////////////////////////
// derivatives
/////////////////////////////////////////////////////////////////////////////
void
Fa_MbenvEnergy::setup_for_derivatives(
                                         pose::Pose & pose,
                                         ScoreFunction const & scfxn
) const
{
  potential_.compute_fa_projection( pose );
  fa_mbenv_weight_ = scfxn.weights()[ fa_mbenv ];
}

//////////////////////////////////////////////////////////////////////////////////////
void
Fa_MbenvEnergy::eval_atom_derivative(
  id::AtomID const & atom_id,
  pose::Pose const & pose,
  kinematics::DomainMap const & /*domain_map*/,
  ScoreFunction const &,// sfxn,
  EnergyMap const & /*weights*/,
  Vector & F1,
  Vector & F2
) const
{

  Size const i( atom_id.rsd() );
  Size const m( atom_id.atomno() );
  conformation::Residue const & rsd1( pose.residue( i ) );

  if ( m > rsd1.nheavyatoms() ) return;

  Vector const heavy_atom_i( rsd1.xyz( m ) );
	//  bool const pos1_fixed( domain_map( i ) != 0 );

        Real cp_weight = 1.0;

        Vector const center(MembraneEmbed_from_pose( pose ).center());

        Vector f1( 0.0 ), f2( 0.0 );

        Real const deriv = memb_lk_dgrefce_(rsd1.atom(m).type()) - lk_dgrefce_(rsd1.atom(m).type());
        Real dE_dZ_over_r = fa_mbenv_weight_ * deriv * Membrane_FAEmbed_from_pose( pose ).fa_proj_deriv(rsd1.seqpos(),m);

        Vector const d_ij = Membrane_FAEmbed_from_pose( pose ).fa_proj_coord(rsd1.seqpos(),m) - heavy_atom_i;
        Real const d_ij_norm = d_ij.length();
        if ( d_ij_norm == Real(0.0) ) return;

        Real const invd = 1.0 / d_ij_norm;
        f2 = d_ij * invd;
        f1 = Membrane_FAEmbed_from_pose( pose ).fa_proj_coord(rsd1.seqpos(),m).cross(heavy_atom_i);
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
  ScoreFunction const &,
  EnergyMap & emap
) const
{
  //std::cout << "BEFORE emap[ fa_mbenv ] " << emap[ fa_mbenv ] << std::endl;
  //emap[ fa_mbenv ] += Membrane_FAEmbed_from_pose( pose ).fa_penalty();
  emap[ fa_mbenv ] += 0.0;
  //std::cout << "AFTER emap[ fa_mbenv ] " << emap[ fa_mbenv ] << std::endl;
}

/// @brief Fa_MbenvEnergy is context independent; indicates that no context graphs are required
void
Fa_MbenvEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}

/// @details Pose must already contain a cenlist object or this method will fail.
Membrane_FAEmbed const &
Fa_MbenvEnergy::Membrane_FAEmbed_from_pose( pose::Pose const & pose ) const
{
  //using core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED;
  return *( utility::pointer::static_pointer_cast< core::scoring::Membrane_FAEmbed const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_FAEMBED ) ));
}

MembraneEmbed const &
Fa_MbenvEnergy::MembraneEmbed_from_pose( pose::Pose const & pose ) const
{
  //using core::pose::datacache::CacheableDataType::MEMBRANE_EMBED;
  assert( pose.data().has( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED ) );
  return *( utility::pointer::static_pointer_cast< core::scoring::MembraneEmbed const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_EMBED ) ));
}

MembraneTopology const &
Fa_MbenvEnergy::MembraneTopology_from_pose( pose::Pose const & pose ) const
{
  //using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
  return *( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
}
core::Size
Fa_MbenvEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

