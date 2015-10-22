// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CenRotDunEnergy.cc
/// @brief  CenRot version of centroid dunbrack energy
/// @author Yuan Liu


// Unit headers
#include <core/pack/dunbrack/cenrot/CenRotDunEnergy.hh>
#include <core/pack/dunbrack/cenrot/CenRotDunEnergyCreator.hh>

// Package headers
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>

#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/angle_deriv.hh>
#include <numeric/deriv/dihedral_deriv.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <numeric/conversions.hh>

//#include <core/kinematics/Jump.hh>
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

// Utility headers

// C++

namespace core {
namespace pack {
namespace dunbrack {
namespace cenrot {

using namespace scoring;
using namespace scoring::methods;

/// @details This must return a fresh instance of the CenRotDunEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
CenRotDunEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new CenRotDunEnergy );
}

ScoreTypes
CenRotDunEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( cen_rot_dun );
	return sts;
}


/// c-tor
CenRotDunEnergy::CenRotDunEnergy() :
	parent( EnergyMethodCreatorOP( new CenRotDunEnergyCreator ) )
{}


/// clone
EnergyMethodOP
CenRotDunEnergy::clone() const {
	return EnergyMethodOP( new CenRotDunEnergy );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
void CenRotDunEnergy::setup_for_scoring( pose::Pose &, ScoreFunction const & ) const {
	// compute interpolated number of neighbors at various distance cutoffs
	//pose.update_residue_neighbors();
	//potential_.compute_centroid_environment( pose );
}

using namespace core::pack::dunbrack;
using namespace core::pack::dunbrack::cenrot;
using namespace core::pack::rotamers;

void CenRotDunEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,  //pose,
	EnergyMap & emap
) const {

	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) return;
	if ( rsd.aa() > core::chemical::num_canonical_aas ) return;
	debug_assert(rsd.residue_type_set().name()==chemical::CENTROID_ROT);

	/// accumulate total energies
	Real dun_score( 0.0 );

	//cal single residue cenrot lib
	SingleResidueRotamerLibraryFactory const & rotlibfact = *SingleResidueRotamerLibraryFactory::get_instance();
	SingleResidueRotamerLibraryCOP residue_rotamer_library(rotlibfact.get( rsd.type() ));

	if ( residue_rotamer_library == 0 ) return;

	SingleResidueCenrotLibraryCOP residue_cenrot_library(
		utility::pointer::dynamic_pointer_cast< SingleResidueCenrotLibrary const >( residue_rotamer_library )
	);

	RotamerLibraryScratchSpace scratch;
	dun_score = residue_cenrot_library->rotamer_energy( rsd, scratch );

	emap[ cen_rot_dun ] += dun_score;
} // residue_energy

// dunbrack term not only contains dU/dr but also
// has dU/dphi and dU/dpsi term, cz it's bb dependent

bool CenRotDunEnergy::defines_dof_derivatives( pose::Pose const & ) const { return true; }

/// @brief Evaluate the phi/psi and chi dihedral derivatives
/// for the input residue.
Real CenRotDunEnergy::eval_residue_dof_derivative(
	conformation::Residue const & rsd,
	scoring::ResSingleMinimizationData const &,
	id::DOF_ID const &,
	id::TorsionID const & tor_id,
	pose::Pose const &,
	scoring::ScoreFunction const &,
	scoring::EnergyMap const & weights
) const {
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	debug_assert(rsd.residue_type_set().name()==chemical::CENTROID_ROT);
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) )  return 0.0;
	if ( rsd.is_virtual_residue() ) return 0.0;
	if ( ! tor_id.valid() )  return 0.0;
	
	debug_assert( rsd.seqpos() == tor_id.rsd() );
	
	Real deriv(0.0);
	
	SingleResidueRotamerLibraryFactory const & rotlibfact = *SingleResidueRotamerLibraryFactory::get_instance();
	SingleResidueRotamerLibraryCOP residue_rotamer_library(rotlibfact.get( rsd.type() ));
	
	if ( residue_rotamer_library==0 ) return 0.0;
	
	SingleResidueCenrotLibraryCOP residue_cenrot_library(
			utility::pointer::dynamic_pointer_cast< SingleResidueCenrotLibrary const >(residue_rotamer_library) );
	
	if ( residue_cenrot_library != 0 && rsd.is_protein() && tor_id.type() == id::BB ) {
		RotamerLibraryScratchSpace scratch;
		residue_cenrot_library->eval_rotameric_energy_bb_dof_deriv( rsd, scratch );
		
		if ( tor_id.torsion() <= DUNBRACK_MAX_BBTOR ) {
			//for backbone torsion angles: phi, psi, omega?
			deriv = scratch.dE_dbb()[tor_id.torsion()];
		}
	}

	return numeric::conversions::degrees( weights[ cen_rot_dun ] * deriv );
}

/// @brief Deprecated.
Real CenRotDunEnergy::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	scoring::ScoreFunction const &,
	scoring::EnergyMap const & weights
) const {
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	debug_assert(pose.residue( tor_id.rsd() ).residue_type_set().name()==chemical::CENTROID_ROT);
	if ( pose.residue( tor_id.rsd() ).has_variant_type( core::chemical::REPLONLY ) ) {
		return 0.0;
	}

	if ( ! tor_id.valid() )  return 0.0;
	if ( pose.residue( tor_id.rsd() ).is_virtual_residue() ) return 0.0;
	
	Real deriv(0.0);
	
	SingleResidueRotamerLibraryFactory const & rotlibfact = *SingleResidueRotamerLibraryFactory::get_instance();
	SingleResidueRotamerLibraryCOP residue_rotamer_library(rotlibfact.get( pose.residue( tor_id.rsd() ).type() ));
	
	if ( residue_rotamer_library==0 ) return 0.0;
	
	SingleResidueCenrotLibraryCOP residue_cenrot_library(
			utility::pointer::dynamic_pointer_cast< SingleResidueCenrotLibrary const >(residue_rotamer_library) );
	
	if ( residue_cenrot_library && pose.residue_type( tor_id.rsd() ).is_protein() && tor_id.type() == id::BB ) {
		
		RotamerLibraryScratchSpace scratch;
		residue_cenrot_library->eval_rotameric_energy_bb_dof_deriv(pose.residue( tor_id.rsd() ), scratch);
		
		if ( tor_id.torsion() <= DUNBRACK_MAX_BBTOR ) {
			//for backbone torsion angles: phi, psi, omega?
			deriv = scratch.dE_dbb()[tor_id.torsion()];
		}
	}
	
	return numeric::conversions::degrees( weights[ cen_rot_dun ] * deriv );
}

void CenRotDunEnergy::eval_residue_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const &, //pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const {
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) return;
	if ( rsd.aa() > core::chemical::num_canonical_aas ) return;
	debug_assert(rsd.residue_type_set().name()==chemical::CENTROID_ROT);

	Real weight = weights[ cen_rot_dun ];

	//cal single residue cenrot lib
	SingleResidueRotamerLibraryFactory const & rotlibfact = *SingleResidueRotamerLibraryFactory::get_instance();
	SingleResidueRotamerLibraryCOP residue_rotamer_library(rotlibfact.get( rsd.type() ));

	if ( residue_rotamer_library==0 ) return;

	SingleResidueCenrotLibraryCOP residue_cenrot_library(
		utility::pointer::dynamic_pointer_cast< SingleResidueCenrotLibrary const >(residue_rotamer_library)
	);

	//get xyz of all the 4 atoms
	//   D -- C -- B -- A
	//   N -- CA - CB - CEN
	Size nA = rsd.atom_index("CEN");
	Size nB = rsd.atom_index("CB");
	Size nC = rsd.atom_index("CA");
	Size nD = rsd.atom_index("N");
	Vector const rA (rsd.atom(nA).xyz());
	Vector const rB (rsd.atom(nB).xyz());
	Vector const rC (rsd.atom(nC).xyz());
	Vector const rD (rsd.atom(nD).xyz());

	RotamerLibraryScratchSpace scratch;
	residue_cenrot_library->rotamer_energy_deriv(rsd, scratch);
	Real4 const & dE_dchi(scratch.dE_dchi());
	Real dE_ddis(dE_dchi[1]);
	Real dE_dang(dE_dchi[2]);
	Real dE_ddih(dE_dchi[3]);

	Vector f1_A, f1_B, f1_C, f1_D;
	Vector f2_A, f2_B, f2_C, f2_D;

	using namespace numeric::deriv;

	//dis contrib
	Real dis_cen_cb(0.0);
	distance_f1_f2_deriv(rA, rB, dis_cen_cb, f1_A, f2_A);
	Real w_dis = weight * dE_ddis;
	atom_derivs[ nA ].f1() += w_dis * f1_A;
	atom_derivs[ nA ].f2() += w_dis * f2_A;
	atom_derivs[ nB ].f1() -= w_dis * f1_A;
	atom_derivs[ nB ].f2() -= w_dis * f2_A;

	//ang contrib
	Real theta_cen_cb_ca(0.0);
	angle_p1_deriv(rA, rB, rC, theta_cen_cb_ca, f1_A, f2_A);
	angle_p2_deriv(rA, rB, rC, theta_cen_cb_ca, f1_B, f2_B);
	angle_p1_deriv(rC, rB, rA, theta_cen_cb_ca, f1_C, f2_C);
	Real w_ang = weight * dE_dang;
	atom_derivs[ nA ].f1() += w_ang * f1_A;
	atom_derivs[ nA ].f2() += w_ang * f2_A;
	atom_derivs[ nB ].f1() += w_ang * f1_B;
	atom_derivs[ nB ].f2() += w_ang * f2_B;
	atom_derivs[ nC ].f1() += w_ang * f1_C;
	atom_derivs[ nC ].f2() += w_ang * f2_C;


	//dih contrib
	Real cosdih(0.0);
	dihedral_p1_cosine_deriv(rA, rB, rC, rD, cosdih, f1_A, f2_A);
	dihedral_p2_cosine_deriv(rA, rB, rC, rD, cosdih, f1_B, f2_B);
	dihedral_p2_cosine_deriv(rD, rC, rB, rA, cosdih, f1_C, f2_C);
	dihedral_p1_cosine_deriv(rD, rC, rB, rA, cosdih, f1_D, f2_D);
	Real w_dih = weight * dE_ddih;
	atom_derivs[ nA ].f1() += w_dih * f1_A;
	atom_derivs[ nA ].f2() += w_dih * f2_A;
	atom_derivs[ nB ].f1() += w_dih * f1_B;
	atom_derivs[ nB ].f2() += w_dih * f2_B;
	atom_derivs[ nC ].f1() += w_dih * f1_C;
	atom_derivs[ nC ].f2() += w_dih * f2_C;
	atom_derivs[ nD ].f1() += w_dih * f1_D;
	atom_derivs[ nD ].f2() += w_dih * f2_D;
}

void CenRotDunEnergy::finalize_total_energy(
	pose::Pose &, //pose,
	ScoreFunction const &,
	EnergyMap &
) const {
}

core::Size CenRotDunEnergy::version() const {
	return 1; // Initial versioning
}

}
}
}
}
