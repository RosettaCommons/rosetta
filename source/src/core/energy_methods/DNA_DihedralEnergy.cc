// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/energy_methods/DNA_DihedralEnergy.cc
/// @brief  dna scoring
/// @author Phil Bradley
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/energy_methods/DNA_DihedralEnergy.hh>
#include <core/energy_methods/DNA_DihedralEnergyCreator.hh>

// Package Headers
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/DNA_DihedralPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/PolymerBondedEnergyContainer.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/symmetry/util.hh>


#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>

#include <core/id/DOF_ID.fwd.hh>
#include <core/id/PartialAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

#include <numeric/deriv/dihedral_deriv.hh>


// Utility headers
#include <numeric/conversions.hh>


// C++


namespace core {
namespace energy_methods {


static basic::Tracer TR( "core.energy_methods.DNA_DihedralEnergy" );

/// @details This must return a fresh instance of the LK_hack class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
DNA_DihedralEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &// options
) const {
	return utility::pointer::make_shared< DNA_DihedralEnergy >();
}

core::scoring::ScoreTypes
DNA_DihedralEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( dna_dihedral_bb );
	sts.push_back( dna_dihedral_chi );
	sts.push_back( dna_dihedral_sugar );
	return sts;
}



/// ctor
DNA_DihedralEnergy::DNA_DihedralEnergy() :
	parent( utility::pointer::make_shared< DNA_DihedralEnergyCreator >() ),
	potential_( core::scoring::ScoringManager::get_instance()->get_DNA_DihedralPotential() )
{
	configure_from_options_system();
}
DNA_DihedralEnergy::DNA_DihedralEnergy( DNA_DihedralEnergy const & src ) :
	parent( utility::pointer::make_shared< DNA_DihedralEnergyCreator >() ),
	potential_( src.potential_ )
{
	configure_from_options_system();
}

/// clone
core::scoring::methods::EnergyMethodOP
DNA_DihedralEnergy::clone() const
{
	return utility::pointer::make_shared< DNA_DihedralEnergy >( *this );
}


void
DNA_DihedralEnergy::configure_from_options_system()
{
	score_delta_ = true;
	score_chi_ = true;
}


bool
DNA_DihedralEnergy::defines_intrares_energy_for_residue(
	conformation::Residue const & rsd
) const
{
	return rsd.is_DNA();
}

bool
DNA_DihedralEnergy::defines_residue_pair_energy(
	pose::Pose const &pose,
	Size i,
	Size j
) const {
	return ( pose.residue(i).is_DNA() && pose.residue(j).is_DNA() );
}

void
DNA_DihedralEnergy::setup_for_scoring(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &
) const
{
	using namespace core::scoring;
	using namespace core::scoring::methods;

	core::scoring::methods::LongRangeEnergyType const & lr_type( long_range_type() );
	core::scoring::Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == nullptr ) {
		create_new_lre_container = true;
	} else {
		core::scoring::LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		core::scoring::PolymerBondedEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::PolymerBondedEnergyContainer > ( lrc ) );
		if ( !dec || !dec->is_valid( pose ) ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		Size nres = pose.size();
		if ( core::pose::symmetry::is_symmetric(pose) ) {
			nres = core::pose::symmetry::symmetry_info(pose)->last_independent_residue();
		}

		TR << "Creating new peptide-bonded energy container (" << nres << ")" << std::endl;
		utility::vector1< ScoreType > s_types;
		s_types.push_back( dna_dihedral_bb );
		s_types.push_back( dna_dihedral_chi );
		s_types.push_back( dna_dihedral_sugar );
		core::scoring::LREnergyContainerOP new_dec( utility::pointer::make_shared< core::scoring::PolymerBondedEnergyContainer >( pose, s_types ) );
		energies.set_long_range_container( lr_type, new_dec );
	}
}

///
void
DNA_DihedralEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	if ( !rsd.is_DNA() ) return;

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
	utility::vector1<core::Real> puckerProbs;
	scoring::dna::get_sugar_pucker_distr( rsd, puckerProbs );

	Real chi_score=0;

	for ( Size pucker=1; pucker<=10; ++pucker ) {
		if ( puckerProbs[pucker]>1e-8 ) {
			Real chi_i_score, dscore_dchi_i;
			potential_.eval_harmonic_sugar_pucker_dependent_chi_torsion_score_and_deriv(
				rsd, pose, pucker-1, chi_i_score, dscore_dchi_i );
			chi_score += puckerProbs[pucker]*chi_i_score;
		}
	}

	Real sugar_score( 0.0 );
	utility::vector1< Real > sugar_torsions;
	scoring::dna::get_sugar_torsions( rsd, sugar_torsions );
	debug_assert( sugar_torsions[1] == rsd.mainchain_torsion(4) ); // first one is delta

	for ( Size pucker=1; pucker<=10; ++pucker ) {
		if ( puckerProbs[pucker]>1e-8 ) {
			for ( Size tor=1; tor<= 4; ++tor ) {
				Real score, dscore_dtor;
				potential_.eval_sugar_torsion_score_and_deriv(
					sugar_torsions[ tor ], tor, rsd, pucker-1, score, dscore_dtor );
				sugar_score += puckerProbs[pucker]*score;
			}
		}
	}

	emap[ core::scoring::dna_dihedral_bb    ] += bb_score;
	emap[ core::scoring::dna_dihedral_chi   ] += chi_score;
	emap[ core::scoring::dna_dihedral_sugar ] += sugar_score;

	//std::pair< std::string, int > puckerOld;
	//scoring::dna::get_sugar_pucker( rsd, puckerOld );
	//TR << "residue_energy: " << rsd.seqpos() << ' ' << rsd.name() << " pucker: "<< puckerOld.first << ' ' <<
	//  puckerOld.second << " bb_score: " << bb_score << " chi_score: " << chi_score << " sugar_score: " <<
	//  sugar_score << std::endl;
}

void
DNA_DihedralEnergy::residue_pair_energy(
	conformation::Residue const & ,//rsd1,
	conformation::Residue const & ,//rsd2,
	pose::Pose const & ,//pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & //emap
) const
{
	//
}


utility::vector1< id::PartialAtomID >
DNA_DihedralEnergy::atoms_with_dof_derivatives(
	conformation::Residue const & rsd,
	pose::Pose const &
) const
{
	utility::vector1< id::PartialAtomID > retlist;
	if ( rsd.is_DNA() ) {
		std::set< id::PartialAtomID > atoms;
		for ( Size tor=1; tor<= 6; ++tor ) {
			if ( tor == 1 && rsd.is_lower_terminus() ) continue;
			if ( tor >= 5 && rsd.is_upper_terminus() ) continue;
			//  if ( tor == 4 ) continue; // delta handled in sugar score
			// tor 4 is handled by the sugar torsion component, even if not by the mainchain torsion
			// component
			conformation::insert_partial_atom_ids_for_mainchain_torsion( rsd, tor, atoms );
		}

		for ( Size at : rsd.chi_atoms()[1] ) {
			atoms.insert( id::PartialAtomID( at, rsd.seqpos() ));
		}

		utility::vector1< utility::vector1< std::string > > sugar_torsion_atoms = scoring::dna::sugar_torsion_atom_names();
		for ( Size ii = 1; ii <= 4; ++ii ) {
			for ( Size jj = 1; jj <= 4; ++jj ) {
				atoms.insert( id::PartialAtomID( rsd.atom_index(sugar_torsion_atoms[ii][jj]), rsd.seqpos() ));
			}
		}

		retlist.resize(atoms.size());
		std::copy(atoms.begin(), atoms.end(), retlist.begin());
	}
	return retlist;

}

Real
DNA_DihedralEnergy::eval_intraresidue_dof_derivative(
	conformation::Residue const &,// rsd,
	core::scoring::ResSingleMinimizationData const &, // min_data,
	id::DOF_ID const & , //dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &, // sfxn,
	core::scoring::EnergyMap const & weights
) const
{
	Real deriv(0.0), score, dscore_dtor;

	bool const is_bb( tor_id.type() == id::BB ), is_chi( tor_id.type() == id::CHI );
	Size const tor( tor_id.torsion() );

	if ( !tor_id.valid() || ( !is_bb && !is_chi ) || !pose.residue( tor_id.rsd() ).is_DNA() ) return 0.0;

	conformation::Residue const & rsd( pose.residue( tor_id.rsd() ) );
	if ( is_bb && tor != 4 ) {
		if ( tor == 2 || tor == 3 ) {
			potential_.eval_harmonic_backbone_torsion_score_and_deriv( tor, rsd, pose, score, dscore_dtor );
			deriv = weights[ core::scoring::dna_dihedral_bb ] * dscore_dtor;
		}
		// other torsion derivatives in res-pair
	} else {
		utility::vector1<core::Real> puckerProbs;
		scoring::dna::get_sugar_pucker_distr( rsd, puckerProbs );
		if ( is_chi && tor == 1 ) {
			// chi
			for ( Size pucker=1; pucker<=10; ++pucker ) {
				if ( puckerProbs[pucker]>1e-8 ) {
					Real chi_i_score, dscore_dchi_i;
					potential_.eval_harmonic_sugar_pucker_dependent_chi_torsion_score_and_deriv(
						rsd, pose, pucker-1, chi_i_score, dscore_dchi_i );
					deriv += weights[ core::scoring::dna_dihedral_chi ] * puckerProbs[pucker] * dscore_dchi_i;
				}
			}
		} else {
			debug_assert( ( is_bb && tor == 4 ) || ( is_chi && ( tor >= 2 && tor <= 4 ) ) );
			Size const sugar_tor( is_bb ? 1 : tor );
			utility::vector1< Real > sugar_torsions;
			scoring::dna::get_sugar_torsions( rsd, sugar_torsions );
			for ( Size pucker=1; pucker<=10; ++pucker ) {
				if ( puckerProbs[pucker]>1e-8 ) {
					Real score, dscore_dtor;
					potential_.eval_sugar_torsion_score_and_deriv(
						sugar_torsions[ sugar_tor ], sugar_tor, rsd, pucker-1, score, dscore_dtor );
					deriv += weights[ core::scoring::dna_dihedral_sugar ] * puckerProbs[pucker]* dscore_dtor;
				}
			}
		}
	}

	// note that the atomtree PHI dofs are in radians
	// use degrees since dE/dangle has angle in denominator
	return numeric::conversions::degrees( deriv );
}


void
DNA_DihedralEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	core::scoring::ResSingleMinimizationData const &,// min_data,
	pose::Pose const & pose,
	core::scoring::EnergyMap const & weights,
	utility::vector1< core::scoring::DerivVectorPair > &atom_derivs
) const
{
	using namespace utility;
	using namespace utility::tools;
	using std::string;

	// we need to compute the dPuckerProb/dx part of the potential here
	vector1< string > names = {
		string( "C1'" ), string( "C2'" ), string( "C3'" ), string( "C4'" ), string( "O4'" ) };

	utility::vector1<core::Real> Es(10,0.0);
	utility::vector1< Real > sugar_torsions;
	scoring::dna::get_sugar_torsions( rsd, sugar_torsions );
	for ( Size pucker=1; pucker<=10; ++pucker ) {
		// chi
		Real chi_i_score, dscore_dchi_i;
		potential_.eval_harmonic_sugar_pucker_dependent_chi_torsion_score_and_deriv(
			rsd, pose, pucker-1, chi_i_score, dscore_dchi_i );
		Es[pucker] += weights[ core::scoring::dna_dihedral_chi ] * chi_i_score;

		// sugar
		for ( Size tor=1; tor<= 4; ++tor ) {
			Real score, dscore_dtor;
			potential_.eval_sugar_torsion_score_and_deriv(
				sugar_torsions[ tor ], tor, rsd, pucker-1, score, dscore_dtor );
			Es[pucker] += weights[ core::scoring::dna_dihedral_sugar ] * score;
		}
	}

	//fd: if this is too slow, we could replace with analytic derivatives...
	utility::vector1<core::Real> puckerProbsP, puckerProbsM;
	conformation::Residue rsdCopy = rsd;
	for ( core::Size ii=1; ii<=5; ++ii ) {
		Size atm_i = rsd.atom_index(names[ii]);
		Vector x_i = rsdCopy.xyz(atm_i);

		for ( core::Size jj=0; jj<3; ++jj ) {
			x_i[jj] += 0.0001;
			rsdCopy.set_xyz(atm_i,x_i);
			scoring::dna::get_sugar_pucker_distr( rsdCopy, puckerProbsP );
			x_i[jj] -= 0.0002;
			rsdCopy.set_xyz(atm_i,x_i);
			scoring::dna::get_sugar_pucker_distr( rsdCopy, puckerProbsM );
			x_i[jj] += 0.0001;
			rsdCopy.set_xyz(atm_i,x_i);

			for ( Size pucker=1; pucker<=10; ++pucker ) {
				Real deriv_ij = (puckerProbsP[pucker]-puckerProbsM[pucker])/0.0002;
				//TR << rsd.seqpos()
				atom_derivs[ atm_i ].f2()[jj] += deriv_ij * Es[pucker];
			}
		}
		Vector deriv_ii = atom_derivs[ atm_i ].f2();
		atom_derivs[ atm_i ].f1() = x_i.cross( x_i-deriv_ii );
	}

	static vector1< string > atom_names = {
		string("C5'"), string("C4'"), string("O4'"), string("C1'"), string("C2'"), string("H2''") };
	vector1< Size > atom_indices(atom_names.size());
	for ( Size ii=1; ii<=atom_names.size(); ++ii ) {
		atom_indices[ii] = rsd.atom_index( atom_names[ii] );
	}

	utility::vector1<core::Real> puckerProbs;
	scoring::dna::get_sugar_pucker_distr( rsd, puckerProbs );
	for ( Size tor=2; tor<= 4; ++tor ) { // delta is sugar torsion #1, in numbering used by the potential
		Real dscore_dtor = 0;
		for ( Size pucker=1; pucker<=10; ++pucker ) {
			if ( puckerProbs[pucker]>1e-8 ) {
				Real score_i, dscore_i_dtor;
				potential_.eval_sugar_torsion_score_and_deriv(
					sugar_torsions[ tor ], tor, rsd, pucker-1, score_i, dscore_i_dtor );
				dscore_dtor +=  weights[ core::scoring::dna_dihedral_sugar ] * puckerProbs[pucker]* dscore_i_dtor;
			}
		}

		// rads -> degs
		dscore_dtor = numeric::conversions::degrees( dscore_dtor );

		Size a1=atom_indices[tor-1], a2=atom_indices[tor], a3=atom_indices[tor+1], a4=atom_indices[tor+2];

		Vector f1(0.0), f2(0.0);
		Real phi;
		numeric::deriv::dihedral_p1_cosine_deriv(
			rsd.xyz( a1 ), rsd.xyz( a2 ), rsd.xyz( a3 ), rsd.xyz( a4 ), phi, f1, f2 );
		atom_derivs[ a1 ].f1() += dscore_dtor * f1;
		atom_derivs[ a1 ].f2() += dscore_dtor * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv(
			rsd.xyz( a1 ), rsd.xyz( a2 ), rsd.xyz( a3 ), rsd.xyz( a4 ), phi, f1, f2 );
		atom_derivs[ a2 ].f1() += dscore_dtor * f1;
		atom_derivs[ a2 ].f2() += dscore_dtor * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv(
			rsd.xyz( a4 ), rsd.xyz( a3 ), rsd.xyz( a2 ), rsd.xyz( a1 ), phi, f1, f2 );
		atom_derivs[ a3 ].f1() += dscore_dtor * f1;
		atom_derivs[ a3 ].f2() += dscore_dtor * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p1_cosine_deriv(
			rsd.xyz( a4 ), rsd.xyz( a3 ), rsd.xyz( a2 ), rsd.xyz( a1 ), phi, f1, f2 );
		atom_derivs[ a4 ].f1() += dscore_dtor * f1;
		atom_derivs[ a4 ].f2() += dscore_dtor * f2;
	}
}


void
DNA_DihedralEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsdA,
	conformation::Residue const & rsdB,
	core::scoring::ResSingleMinimizationData const &,
	core::scoring::ResSingleMinimizationData const &,
	core::scoring::ResPairMinimizationData const &,
	pose::Pose const & pose,
	core::scoring::EnergyMap const & weights,
	utility::vector1< core::scoring::DerivVectorPair > & rA_atom_derivs,
	utility::vector1< core::scoring::DerivVectorPair > & rB_atom_derivs
) const
{
	// we compute bb tors 5 & 6 on rsd1, and bb tor 1 on rsd2 here
	id::AtomID id1, id2, id3, id4;
	Real deriv, score, dscore_dtor1, phi;
	Real dscore_dtor2, dscore_deps2, dscore_dzeta2;
	Real dscore_dtor5, dscore_deps5, dscore_dzeta5;
	Real dscore_dtor6, dscore_deps6, dscore_dzeta6;
	Vector f1(0.0), f2(0.0);

	conformation::Residue const &rsd1 = rsdA.seqpos() < rsdB.seqpos() ? rsdA : rsdB;
	conformation::Residue const &rsd2 = rsdA.seqpos() < rsdB.seqpos() ? rsdB : rsdA;

	utility::vector1< core::scoring::DerivVectorPair > &r1_atom_derivs = rsdA.seqpos() < rsdB.seqpos() ? rA_atom_derivs : rB_atom_derivs;
	utility::vector1< core::scoring::DerivVectorPair > &r2_atom_derivs = rsdA.seqpos() < rsdB.seqpos() ? rB_atom_derivs : rA_atom_derivs;

	// depsilon and dzeta parts of beta
	potential_.eval_harmonic_backbone_torsion_score_and_deriv( 1, rsd2, pose, score, dscore_dtor1 );
	potential_.eval_harmonic_backbone_torsion_score_and_deriv( 2, rsd2, pose, score, dscore_dtor2, dscore_deps2, dscore_dzeta2 ); // NOTE RSD2!
	potential_.eval_harmonic_backbone_torsion_score_and_deriv( 5, rsd1, pose, score, dscore_dtor5, dscore_deps5, dscore_dzeta5 );
	potential_.eval_harmonic_backbone_torsion_score_and_deriv( 6, rsd1, pose, score, dscore_dtor6, dscore_deps6, dscore_dzeta6 );

	core::id::TorsionID tor5( rsd1.seqpos(), id::BB, 5 );
	bool tor5_invalid = pose.conformation().get_torsion_angle_atom_ids( tor5, id1, id2, id3, id4 );
	if ( !tor5_invalid ) {
		deriv = weights[ core::scoring::dna_dihedral_bb ] * numeric::conversions::degrees( dscore_dtor5+dscore_deps2+dscore_deps5+dscore_deps6 );
		numeric::deriv::dihedral_p1_cosine_deriv(rsd1.xyz(id1.atomno()), rsd1.xyz(id2.atomno()), rsd1.xyz(id3.atomno()), rsd2.xyz(id4.atomno()), phi, f1, f2);
		r1_atom_derivs[ id1.atomno() ].f1() += deriv * f1;
		r1_atom_derivs[ id1.atomno() ].f2() += deriv * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv(rsd1.xyz(id1.atomno()), rsd1.xyz(id2.atomno()), rsd1.xyz(id3.atomno()), rsd2.xyz(id4.atomno()), phi, f1, f2);
		r1_atom_derivs[ id2.atomno() ].f1() += deriv * f1;
		r1_atom_derivs[ id2.atomno() ].f2() += deriv * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv(rsd2.xyz(id4.atomno()), rsd1.xyz(id3.atomno()), rsd1.xyz(id2.atomno()), rsd1.xyz(id1.atomno()), phi, f1, f2);
		r1_atom_derivs[ id3.atomno() ].f1() += deriv * f1;
		r1_atom_derivs[ id3.atomno() ].f2() += deriv * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p1_cosine_deriv(rsd2.xyz(id4.atomno()), rsd1.xyz(id3.atomno()), rsd1.xyz(id2.atomno()), rsd1.xyz(id1.atomno()), phi, f1, f2);
		r2_atom_derivs[ id4.atomno() ].f1() += deriv * f1;
		r2_atom_derivs[ id4.atomno() ].f2() += deriv * f2;
	}

	core::id::TorsionID tor6( rsd1.seqpos(), id::BB, 6 );
	bool tor6_invalid = pose.conformation().get_torsion_angle_atom_ids( tor6, id1, id2, id3, id4 );
	if ( !tor6_invalid ) {
		deriv = weights[ core::scoring::dna_dihedral_bb ] * numeric::conversions::degrees( dscore_dtor6+dscore_dzeta2+dscore_dzeta5+dscore_dzeta6  );
		numeric::deriv::dihedral_p1_cosine_deriv(rsd1.xyz(id1.atomno()), rsd1.xyz(id2.atomno()), rsd2.xyz(id3.atomno()), rsd2.xyz(id4.atomno()), phi, f1, f2);
		r1_atom_derivs[ id1.atomno() ].f1() += deriv * f1;
		r1_atom_derivs[ id1.atomno() ].f2() += deriv * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv(rsd1.xyz(id1.atomno()), rsd1.xyz(id2.atomno()), rsd2.xyz(id3.atomno()), rsd2.xyz(id4.atomno()), phi, f1, f2);
		r1_atom_derivs[ id2.atomno() ].f1() += deriv * f1;
		r1_atom_derivs[ id2.atomno() ].f2() += deriv * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv(rsd2.xyz(id4.atomno()), rsd2.xyz(id3.atomno()), rsd1.xyz(id2.atomno()), rsd1.xyz(id1.atomno()), phi, f1, f2);
		r2_atom_derivs[ id3.atomno() ].f1() += deriv * f1;
		r2_atom_derivs[ id3.atomno() ].f2() += deriv * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p1_cosine_deriv(rsd2.xyz(id4.atomno()), rsd2.xyz(id3.atomno()), rsd1.xyz(id2.atomno()), rsd1.xyz(id1.atomno()), phi, f1, f2);
		r2_atom_derivs[ id4.atomno() ].f1() += deriv * f1;
		r2_atom_derivs[ id4.atomno() ].f2() += deriv * f2;
	}

	core::id::TorsionID tor1( rsd2.seqpos(), id::BB, 1 );
	bool tor1_invalid = pose.conformation().get_torsion_angle_atom_ids( tor1, id1, id2, id3, id4 );
	if ( !tor1_invalid ) {
		// rads -> degs
		deriv = weights[ core::scoring::dna_dihedral_bb ] * numeric::conversions::degrees( dscore_dtor1 );
		numeric::deriv::dihedral_p1_cosine_deriv(rsd1.xyz(id1.atomno()), rsd2.xyz(id2.atomno()), rsd2.xyz(id3.atomno()), rsd2.xyz(id4.atomno()), phi, f1, f2);
		r1_atom_derivs[ id1.atomno() ].f1() += deriv * f1;
		r1_atom_derivs[ id1.atomno() ].f2() += deriv * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv(rsd1.xyz(id1.atomno()), rsd2.xyz(id2.atomno()), rsd2.xyz(id3.atomno()), rsd2.xyz(id4.atomno()), phi, f1, f2);
		r2_atom_derivs[ id2.atomno() ].f1() += deriv * f1;
		r2_atom_derivs[ id2.atomno() ].f2() += deriv * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv(rsd2.xyz(id4.atomno()), rsd2.xyz(id3.atomno()), rsd2.xyz(id2.atomno()), rsd1.xyz(id1.atomno()), phi, f1, f2);
		r2_atom_derivs[ id3.atomno() ].f1() += deriv * f1;
		r2_atom_derivs[ id3.atomno() ].f2() += deriv * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p1_cosine_deriv(rsd2.xyz(id4.atomno()), rsd2.xyz(id3.atomno()), rsd2.xyz(id2.atomno()), rsd1.xyz(id1.atomno()), phi, f1, f2);
		r2_atom_derivs[ id4.atomno() ].f1() += deriv * f1;
		r2_atom_derivs[ id4.atomno() ].f2() += deriv * f2;
	}
}

core::scoring::methods::LongRangeEnergyType
DNA_DihedralEnergy::long_range_type() const { return core::scoring::methods::dna_dihedral_lr; }

void
DNA_DihedralEnergy::indicate_required_context_graphs(utility::vector1< bool > & ) const {}


} // scoring
} // core

