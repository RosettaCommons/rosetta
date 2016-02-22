// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/methods/DNA_DihedralEnergy.cc
/// @brief  dna scoring
/// @author Phil Bradley
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/methods/DNA_DihedralEnergy.hh>
#include <core/scoring/methods/DNA_DihedralEnergyCreator.hh>

// Package Headers
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/DNA_DihedralPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/XYZ_Func.hh>


#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>

#include <core/id/TorsionID.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/types.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
//#include <core/options/util.hh>


// Utility headers
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <utility/tools/make_vector1.hh>


// C++


namespace core {
namespace scoring {
namespace methods {

static basic::Tracer TR("core.scoring.methods.DNA_DihedralEnergy" );

/// @details This must return a fresh instance of the LK_hack class,
/// never an instance already in use
methods::EnergyMethodOP
DNA_DihedralEnergyCreator::create_energy_method(
	EnergyMethodOptions const &// options
) const {
	return EnergyMethodOP( new DNA_DihedralEnergy() );
}

ScoreTypes
DNA_DihedralEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	//sts.push_back( dna_dihedral );
	sts.push_back( dna_dihedral_bb );
	sts.push_back( dna_dihedral_chi );
	sts.push_back( dna_dihedral_sugar );
	return sts;
}



/// ctor
DNA_DihedralEnergy::DNA_DihedralEnergy() :
	parent( EnergyMethodCreatorOP( new DNA_DihedralEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_DNA_DihedralPotential() )
{
	configure_from_options_system();
	//  add_score_type( dna_dihedral );
	//  add_score_type( dna_dihedral_bb );
	//  add_score_type( dna_dihedral_chi );
	//  add_score_type( dna_dihedral_sugar );
}
DNA_DihedralEnergy::DNA_DihedralEnergy( DNA_DihedralEnergy const & src ) :
	parent( EnergyMethodCreatorOP( new DNA_DihedralEnergyCreator ) ),
	potential_( src.potential_ )
{
	configure_from_options_system();
	//  add_score_type( dna_dihedral );
	//  add_score_type( dna_dihedral_bb );
	//  add_score_type( dna_dihedral_chi );
	//  add_score_type( dna_dihedral_sugar );
}

/// clone
EnergyMethodOP
DNA_DihedralEnergy::clone() const
{
	return EnergyMethodOP( new DNA_DihedralEnergy( *this ) );
}


void
DNA_DihedralEnergy::configure_from_options_system()
{
	score_delta_ = true;
	score_chi_ = true;
	//  score_delta_ = options::option[ options::OptionKeys::dna::specificity::score_delta_dna_dihedral_energy ];
	//  score_chi_   = options::option[ options::OptionKeys::dna::specificity::score_chi_dna_dihedral_energy ];
	//  use_jjh_dna_dihedral_potential_ =
	//   options::option[ options::OptionKeys::dna::specificity::use_jjh_dna_dihedral_potential ];

	/// HACK: if we're chi scoring and not using jjh, then use jjh just for chi
	//  use_jjh_dna_dihedral_potential_chi_only_ = ( score_chi_ && !use_jjh_dna_dihedral_potential_ );
}


bool
DNA_DihedralEnergy::defines_score_for_residue(
	conformation::Residue const & rsd
) const
{
	return rsd.is_DNA();
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

///
void
DNA_DihedralEnergy::residue_energy(
	conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	EnergyMap & emap
) const
{
	if ( rsd.is_DNA() ) {

		//   TR.Trace << "defines_score_for_residue: " << this->defines_score_for_residue( rsd ) <<
		//    " use_extended_residue_energy_interface: "<< this->use_extended_residue_energy_interface() <<
		//    " defines_dof_derivatives: " << this->defines_dof_derivatives( pose ) << std::endl;

		assert( this->defines_score_for_residue( rsd ) );
		assert( !this->use_extended_residue_energy_interface() );
		assert( this->defines_dof_derivatives( pose ) );

		Real bb_score( 0.0 );
		for ( Size tor=1; tor<= 6; ++tor ) {
			if ( tor == 1 && rsd.is_lower_terminus() ) continue;
			if ( tor >= 5 && rsd.is_upper_terminus() ) continue;
			if ( tor == 4 ) continue; // delta handled in sugar score
			Real score, dscore_dtor;
			potential_.eval_harmonic_backbone_torsion_score_and_deriv( tor, rsd, pose, score, dscore_dtor );
			bb_score += score;
		}

		/// need to know the sugar pucker for chi scoring and for delta + three sugar torsions
		std::pair< std::string, int > pucker;
		scoring::dna::get_sugar_pucker( rsd, pucker );
		Size const pucker_index( pucker.second );


		Real chi_score, dscore_dchi;
		potential_.eval_harmonic_sugar_pucker_dependent_chi_torsion_score_and_deriv( rsd, pose, pucker_index, chi_score,
			dscore_dchi );

		Real sugar_score( 0.0 );
		{ /// score the sugar deviation, this only goes into emap[ dna_dihedral_sugar ]
			utility::vector1< Real > torsions;
			scoring::dna::get_sugar_torsions( rsd, torsions );
			assert( torsions[1] == rsd.mainchain_torsion(4) ); // first one is delta
			for ( Size tor=1; tor<= 4; ++tor ) {
				Real score, dscore_dtor;
				potential_.eval_sugar_torsion_score_and_deriv( torsions[ tor ], tor, rsd, pucker_index, score, dscore_dtor );
				sugar_score += score;
			}
		}

		emap[ dna_dihedral_bb    ] += bb_score;
		emap[ dna_dihedral_chi   ] += chi_score;
		emap[ dna_dihedral_sugar ] += sugar_score;

		//   TR.Trace << "residue_energy: " << rsd.seqpos() << ' ' << rsd.name() << " pucker: "<< pucker.first << ' ' <<
		//    pucker.second << " bb_score: " << bb_score << " chi_score: " << chi_score << " sugar_score: " <<
		//    sugar_score << std::endl;

		if ( false ) { // sanity check on atom naming conventions, not sure this is still correct now that we renamed
			Real const h2dihedral( numeric::dihedral_degrees(rsd.xyz("H2''"),rsd.xyz("C2'"),rsd.xyz("C1'"),rsd.xyz("H2'")));
			if ( h2dihedral > 0.0 ) utility_exit_with_message("Bad H2'' -- H2' dihedral!");
			if ( !rsd.is_lower_terminus() && rsd.has_variant_type(chemical::CUTPOINT_UPPER) ) {
				Real const opdihedral( numeric::dihedral_degrees(rsd.xyz("OP1"),rsd.xyz("OVU1"),rsd.xyz("O5'"),rsd.xyz( "OP2")));
				if ( opdihedral < 0.0 ) {
					std::cerr << "Bad OP1 -- OP2 dihedral! " << opdihedral << ' ' << rsd.seqpos() << std::endl;
					utility_exit_with_message("Bad OP1 -- OP2 dihedral!");
				}
			}
			//std::cout << "standard dihedrals -120? " << h2dihedral << " +120? " << opdihedral << std::endl;
		}

	}
}


// ///
// void
// DNA_DihedralEnergy::residue_energy_old(
//  conformation::Residue const & rsd,
//  core::pose::Pose const & pose,
//  EnergyMap & emap
// ) const
// {
//  if ( rsd.is_DNA() ) {
//   if ( use_jjh_dna_dihedral_potential_ ) {
//    Real total_score(0.0);
//    for ( Size tor=1; tor<= 6; ++tor ) {
//     if ( !score_delta_ && tor == 4 ) continue;
//     Real score(0.0), dscore_dchi;
//     potential_.jjh_potential().eval_dna_bform_bb_torsion_score_residue( rsd, score, dscore_dchi, tor );
//     total_score += score;
//     emap[ dna_dihedral_bb ] = score;
//    }
//    if ( score_chi_ ) {
//     Real score(0.0), dscore_dchi;
//     potential_.jjh_potential().eval_dna_bform_chi_torsion_score_residue( rsd, score, dscore_dchi );
//     total_score += score;
//     emap[ dna_dihedral_chi ] = score;
//    }
//    emap[ dna_dihedral ] = total_score;


//   } else {
//    Real logprob_bb(0.0), logprob_chi(0.0);
//    for ( Size tor=1; tor<= 6; ++tor ) {
//     if ( !score_delta_ && tor == 4 ) continue;
//     Real const prob( potential_.torsion_probability( rsd, tor ) );
//     logprob_bb += std::log( prob );
//    }
//    if ( score_chi_ ) {
//     if ( use_jjh_dna_dihedral_potential_chi_only_ ) {
//      Real score(0.0), dscore_dchi;
//      potential_.jjh_potential().eval_dna_bform_chi_torsion_score_residue( rsd, score, dscore_dchi );
//      logprob_chi += ( -1 * score ); // -1 since we multiply by -1 below
//     } else {
//      Real const prob( potential_.torsion_probability( rsd, 7 ) );
//      logprob_chi += std::log( prob );
//     }
//    }

//    emap[ dna_dihedral     ] = -1.0 * ( logprob_bb + logprob_chi );
//    emap[ dna_dihedral_bb  ] = -1.0 * logprob_bb;
//    emap[ dna_dihedral_chi ] = -1.0 * logprob_chi;


//    { /// score the sugar deviation, this only goes into emap[ dna_dihedral_sugar ]
//     Real sugar_score( 0.0 );
//     std::pair< std::string, int > pucker;
//     scoring::dna::get_sugar_pucker( rsd, pucker );
//     Size const pucker_index( pucker.second );
//     utility::vector1< Real > torsions;
//     scoring::dna::get_sugar_torsions( rsd, torsions );
//     for ( Size tor=1; tor<= 4; ++tor ) {
//      Real score, dscore_dtor;
//      potential_.eval_sugar_torsion_score_and_deriv( torsions[ tor ], tor, rsd, pucker_index, score, dscore_dtor );
//      sugar_score += score;
//     }
//     emap[ dna_dihedral_sugar ] = sugar_score;
//    }
//   }

//   { // sanity check on atom naming conventions
//    Real const h2dihedral( numeric::dihedral_degrees(rsd.xyz("H2''"),rsd.xyz("C2'"),rsd.xyz("C1'"),rsd.xyz("H2'")));
//    if ( h2dihedral > 0.0 ) utility_exit_with_message("Bad H2'' -- H2' dihedral!");
//    if ( !rsd.is_lower_terminus() && rsd.has_variant_type(chemical::CUTPOINT_UPPER) ) {
//     Real const opdihedral( numeric::dihedral_degrees(rsd.xyz("OP1"),rsd.xyz("OVU1"),rsd.xyz("O5'"),rsd.xyz( "OP2")));
//     if ( opdihedral < 0.0 ) {
//      std::cerr << "Bad OP1 -- OP2 dihedral! " << opdihedral << ' ' << rsd.seqpos() << std::endl;
//      utility_exit_with_message("Bad OP1 -- OP2 dihedral!");
//     }
//    }
//    //std::cout << "standard dihedrals -120? " << h2dihedral << " +120? " << opdihedral << std::endl;
//   }

//  } else {
//   emap[ dna_dihedral ] = 0.0;
//   emap[ dna_dihedral_bb ] = 0.0;
//   emap[ dna_dihedral_chi ] = 0.0;
//   emap[ dna_dihedral_sugar ] = 0.0;
//  }
// }


///
Real
DNA_DihedralEnergy::eval_dof_derivative(
	id::DOF_ID const &, // dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	ScoreFunction const &,// sfxn,
	EnergyMap const & weights
) const
{
	Real deriv(0.0), score, dscore_dtor;

	bool const is_bb( tor_id.type() == id::BB ), is_chi( tor_id.type() == id::CHI );

	Size const tor( tor_id.torsion() );

	//  TR.Trace << "eval_dof_derivative: dof_id= " << dof_id << " tor_id= " << tor_id << std::endl;

	if ( tor_id.valid() && ( is_bb || is_chi ) && pose.residue( tor_id.rsd() ).is_DNA() ) {
		conformation::Residue const & rsd( pose.residue( tor_id.rsd() ) );
		if ( is_bb && tor != 4 ) {
			potential_.eval_harmonic_backbone_torsion_score_and_deriv( tor, rsd, pose, score, dscore_dtor );
			deriv = ( /*weights[ dna_dihedral ] + */ weights[ dna_dihedral_bb ] ) * dscore_dtor;
		} else {
			/// need to know the pucker for chi, delta, "chi2-4"
			std::pair< std::string, int > pucker;
			scoring::dna::get_sugar_pucker( rsd, pucker );
			Size const pucker_index( pucker.second );
			if ( is_chi && tor == 1 ) {
				// chi
				potential_.eval_harmonic_sugar_pucker_dependent_chi_torsion_score_and_deriv( rsd, pose, pucker_index, score,
					dscore_dtor );
				deriv = ( /*weights[ dna_dihedral ] + */ weights[ dna_dihedral_chi ] ) * dscore_dtor;
			} else {
				assert( ( is_bb && tor == 4 ) || ( is_chi && ( tor >= 2 && tor <= 4 ) ) );
				Size const sugar_tor( is_bb ? 1 : tor );
				utility::vector1< Real > torsions;
				scoring::dna::get_sugar_torsions( rsd, torsions );
				potential_.eval_sugar_torsion_score_and_deriv( torsions[ sugar_tor ], sugar_tor, rsd, pucker_index, score,
					dscore_dtor );

				deriv = weights[ dna_dihedral_sugar ] * dscore_dtor;
			}
		}
	}

	// note that the atomtree PHI dofs are in radians
	// use degrees since dE/dangle has angle in denominator
	return numeric::conversions::degrees( deriv );
}


Real
DNA_DihedralEnergy::eval_residue_dof_derivative(
	conformation::Residue const &,// rsd,
	ResSingleMinimizationData const &, // min_data,
	id::DOF_ID const & dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights
) const
{
	return eval_dof_derivative( dof_id, tor_id, pose, sfxn, weights );
}



void
DNA_DihedralEnergy::eval_residue_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,// min_data,
	pose::Pose const &,// pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	using utility::vector1;
	using utility::tools::make_vector1;
	using std::string;
	using id::AtomID;
	using namespace scoring::constraints;
	using namespace scoring::func;
	using numeric::conversions::radians;

	static vector1< string > atom_names
		( make_vector1( string("C5'"), string("C4'"), string("O4'"), string("C1'"), string("C2'"), string("H2''")));

	vector1< AtomID > atom_ids;
	Size const seqpos( rsd.seqpos() );
	for ( Size i=1; i<= atom_names.size(); ++i ) atom_ids.push_back( AtomID( rsd.atom_index( atom_names[i] ), seqpos ) );

	/// need to know the sugar pucker
	std::pair< string, int > pucker;
	scoring::dna::get_sugar_pucker( rsd, pucker );
	Size const pucker_index( pucker.second );

	//  TR.Trace << "eval_residue_derivatives: " << rsd.seqpos() << ' ' << rsd.name() << " pucker: "<< pucker.first << ' ' <<
	//   pucker.second << std::endl;

	ResidueXYZ const rsd_xyz( rsd );

	for ( Size tor=2; tor<= 4; ++tor ) { // delta is sugar torsion #1, in numbering used by the potential
		// this is silly and slow, create dihedral constraints to calc derivs
		Real mean,sdev;
		potential_.get_sugar_torsion_mean_and_sdev( tor, rsd, pucker_index, mean, sdev );
		DihedralConstraint const cst( atom_ids[tor-1], atom_ids[tor], atom_ids[tor+1], atom_ids[tor+2],
			func::FuncOP( new CircularHarmonicFunc( radians(mean), radians(sdev) )),
			dna_dihedral_sugar );

		if ( false ) { // debugging
			utility::vector1< Real > torsions;
			scoring::dna::get_sugar_torsions( rsd, torsions );
			Real const recompute_torsion( numeric::dihedral_degrees( rsd.xyz( atom_ids[tor-1].atomno() ),
				rsd.xyz( atom_ids[tor  ].atomno() ),
				rsd.xyz( atom_ids[tor+1].atomno() ),
				rsd.xyz( atom_ids[tor+2].atomno() ) ) );
			Real const dev( basic::subtract_degree_angles( recompute_torsion, torsions[tor] ) );
			TR.Trace << "recompute_torsion: " << dev << ' ' << recompute_torsion << ' ' << torsions[tor] << std::endl;
			runtime_assert( std::fabs( dev ) <1e-3 );
		}

		for ( Size i=tor-1; i<= tor+2; ++i ) {
			cst.fill_f1_f2( atom_ids[i], rsd_xyz,
				atom_derivs[ atom_ids[i].atomno() ].f1(), atom_derivs[ atom_ids[i].atomno() ].f2(), weights );
		}
	}



}

///
// Real
// DNA_DihedralEnergy::eval_dof_derivative_old(
//  id::DOF_ID const &,// dof_id,
//  id::TorsionID const & tor_id,
//  pose::Pose const & pose,
//  ScoreFunction const &,// sfxn,
//  EnergyMap const & weights
// ) const
// {
//  Real deriv(0.0);
//  bool const is_bb( tor_id.type() == id::BB ), is_chi( tor_id.type() == id::CHI );
//  if ( tor_id.valid() && ( is_bb || is_chi ) && pose.residue( tor_id.rsd() ).is_DNA() ) {
//   conformation::Residue const & rsd( pose.residue( tor_id.rsd() ) );
//   if ( ( is_bb  && ( score_delta_ || tor_id.torsion() != 4 ) ) ||
//      ( is_chi && score_chi_ && tor_id.torsion() == 1 ) ) {
//    // standard bb or chi torsion
//    Real scorefxn_weight( 0.0 ); //weights[ dna_dihedral ] );
//    if ( is_bb  ) scorefxn_weight += weights[ dna_dihedral_bb ];
//    if ( is_chi ) scorefxn_weight += weights[ dna_dihedral_chi ];
//    if ( use_jjh_dna_dihedral_potential_ ) {
//     Real score(0.0);
//     if ( is_bb ) {
//      potential_.jjh_potential().eval_dna_bform_bb_torsion_score_residue( rsd, score, deriv, tor_id.torsion());
//      deriv = scorefxn_weight * numeric::conversions::radians( deriv ); // already internally converted
//     } else {
//      potential_.jjh_potential().eval_dna_bform_chi_torsion_score_residue( rsd, score, deriv );
//      deriv = scorefxn_weight * -2.0 * deriv; // not already internally converted
//     }
//     /// not sure if we should be converting via degrees as below, have to check with numerical derivs...
//    } else {
//     if ( is_chi && use_jjh_dna_dihedral_potential_chi_only_ ) {
//      Real score(0.0);
//      potential_.jjh_potential().eval_dna_bform_chi_torsion_score_residue( rsd, score, deriv );
//      deriv = scorefxn_weight * -2.0 * deriv; // not already internally converted
//     } else {
//      Size const tor_index( is_bb ? tor_id.torsion() : 7 );
//      Real const prob ( potential_.torsion_probability      ( rsd, tor_index ) );
//      Real const dprob( potential_.torsion_probability_deriv( rsd, tor_index ) );
//      deriv = scorefxn_weight * -1.0 * dprob / prob;
//     }
//    }
//   }


//   if ( ( is_chi && tor_id.torsion() > 1 ) || ( is_bb && tor_id.torsion() == 4 ) ) {
//    if ( is_bb ) assert( !score_delta_ ); // otherwise we get the deriv wrong, should sum in from above...
//    /// a sugar dihedral
//    assert( tor_id.torsion() <= 4 );
//    Size const tor( is_bb ? 1 : tor_id.torsion() );

//    std::pair< std::string, int > pucker;
//    scoring::dna::get_sugar_pucker( rsd, pucker );
//    Size const pucker_index( pucker.second );
//    utility::vector1< Real > torsions;
//    scoring::dna::get_sugar_torsions( rsd, torsions );
//    Real score;
//    potential_.eval_sugar_torsion_score_and_deriv( torsions[ tor ], tor, rsd, pucker_index, score, deriv );

//    deriv *= weights[ dna_dihedral_sugar ];
//   }

//  }
//  // note that the atomtree PHI dofs are in radians
//  // use degrees since dE/dangle has angle in denominator
//  return numeric::conversions::degrees( deriv );
// }

/// @brief DNA_Dihedral Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
DNA_DihedralEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
)
const
{}


} // methods
} // scoring
} // core

