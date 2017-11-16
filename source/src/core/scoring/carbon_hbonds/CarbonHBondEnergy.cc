// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/carbon_hbonds/CarbonHBondEnergy.fwd.hh
/// @brief  Hydrogen bond energy method forward declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay
/// @author Rhiju Das
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)

// Unit Headers
#include <core/scoring/carbon_hbonds/CarbonHBondEnergy.hh>
#include <core/scoring/carbon_hbonds/CarbonHBondEnergyCreator.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/carbon_hbonds/CarbonHBondPotential.hh>

// Project headers
#include <ObjexxFCL/format.hh>

#include <core/conformation/Residue.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/chemical/AtomType.hh>

#include <basic/Tracer.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

#include <core/scoring/DerivVectorPair.hh>
#include <utility/vector1.hh>
#include <boost/bind.hpp>

//Auto Headers
#include <core/scoring/EnergyMap.hh>


static basic::Tracer tr( "core.scoring.carbon_hbonds.CarbonHBondEnergy" );


namespace core {
namespace scoring {
namespace carbon_hbonds {

using namespace ObjexxFCL::format;

/// @details This must return a fresh instance of the CarbonHBondEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
CarbonHBondEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new CarbonHBondEnergy );
}

ScoreTypes
CarbonHBondEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( ch_bond );
	sts.push_back( ch_bond_sc_sc );
	sts.push_back( ch_bond_bb_sc );
	sts.push_back( ch_bond_bb_bb );
	return sts;
}


/// @brief copy c-tor
CarbonHBondEnergy::CarbonHBondEnergy() :
	parent( methods::EnergyMethodCreatorOP( new CarbonHBondEnergyCreator ) ),
	carbon_hbond_potential_( ScoringManager::get_instance()->get_CarbonHBondPotential() ),
	max_dis_( carbon_hbond_potential_.max_dis() ),
	max_dis2_( max_dis_*max_dis_ ),
	path_dist_cutoff_( 4 ),
	orientation_dep_rna_ch_o_bonds_( ! basic::options::option[ basic::options::OptionKeys::score::rna::disable_orientation_dependent_rna_ch_o_bonds ]),
	verbose_( false )
{}

/// copy ctor
CarbonHBondEnergy::CarbonHBondEnergy( CarbonHBondEnergy const & src ):
	parent( src ),
	carbon_hbond_potential_( ScoringManager::get_instance()->get_CarbonHBondPotential() ),
	max_dis_( carbon_hbond_potential_.max_dis() ),
	max_dis2_( max_dis_*max_dis_ ),
	path_dist_cutoff_( src.path_dist_cutoff_ ),
	orientation_dep_rna_ch_o_bonds_(  src.orientation_dep_rna_ch_o_bonds_ ),
	verbose_( src.verbose_ )
{}

/// clone
methods::EnergyMethodOP
CarbonHBondEnergy::clone() const
{
	return methods::EnergyMethodOP( new CarbonHBondEnergy( *this ) );
}


void
CarbonHBondEnergy::setup_for_scoring( pose::Pose & /*pose*/, ScoreFunction const & ) const
{
	// nothing for now.
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// Everything in here.
void
CarbonHBondEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	Real bb_bb(0.0);
	Real bb_sc(0.0);
	Real sc_sc(0.0);
	Real ch_bond_E =
		res_res_carbon_hbond_one_way( rsd1, rsd2, bb_bb, bb_sc, sc_sc) +
		res_res_carbon_hbond_one_way( rsd2, rsd1, bb_bb, bb_sc, sc_sc ) ;
	emap[ ch_bond_bb_bb ] += bb_bb;
	emap[ ch_bond_bb_sc ] += bb_sc;
	emap[ ch_bond_sc_sc ] += sc_sc;
	// store the energies
	emap[ ch_bond ] += ch_bond_E;
	// std::cout << "RUNNING SUM: " << rsd1.seqpos() <<  " " << rsd2.seqpos() << " " << ch_bond_E << " " << emap[ ch_bond ] << std::endl;
}

bool
CarbonHBondEnergy::minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

void
CarbonHBondEnergy::backbone_backbone_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	Real ch_bond_E =
		bb_bb_carbon_hbond_one_way( rsd1, rsd2 ) +
		bb_bb_carbon_hbond_one_way( rsd2, rsd1 ) ;

	// store the energies
	emap[ ch_bond ] += ch_bond_E;
	emap[ ch_bond_bb_bb ] += ch_bond_E;
	// std::cout << "RUNNING SUM: " << rsd1.seqpos() <<  " " << rsd2.seqpos() << " " << ch_bond_E << " " << emap[ ch_bond ] << std::endl;
}


void
CarbonHBondEnergy::backbone_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	Real ch_bond_E =
		bb_sc_carbon_hbond_one_way( rsd1, rsd2 ) +
		sc_bb_carbon_hbond_one_way( rsd2, rsd1 ) ;

	// store the energies
	emap[ ch_bond ] += ch_bond_E;
	emap[ ch_bond_bb_sc ] += ch_bond_E;

	// std::cout << "RUNNING SUM: " << rsd1.seqpos() <<  " " << rsd2.seqpos() << " " << ch_bond_E << " " << emap[ ch_bond ] << std::endl;
}


void
CarbonHBondEnergy::sidechain_sidechain_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap & emap
) const {
	Real ch_bond_E =
		sc_sc_carbon_hbond_one_way( rsd1, rsd2 ) +
		sc_sc_carbon_hbond_one_way( rsd2, rsd1 ) ;

	// store the energies
	emap[ ch_bond ] += ch_bond_E;
	emap[ ch_bond_sc_sc ] += ch_bond_E;
	// std::cout << "RUNNING SUM: " << rsd1.seqpos() <<  " " << rsd2.seqpos() << " " << ch_bond_E << " " << emap[ ch_bond ] << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
// Look more than four atoms away.
bool
CarbonHBondEnergy::path_distance_OK(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Size const ii,
	Size const jj
) const {
	Size const i( rsd1.seqpos() );
	Size const j( rsd2.seqpos() );
	// std::cout << i << " " << j << " " << ii << " " << jj << std::endl;
	if ( i == j && Size( rsd1.path_distance(ii,jj) ) <= path_dist_cutoff_ ) return false;
	else if ( rsd1.is_bonded( rsd2 ) ) {
		Size const path_size =
			rsd1.path_distance( ii, rsd1.connect_atom( rsd2 ) ) +
			rsd2.path_distance( jj, rsd2.connect_atom( rsd1 ) ) + 1;
		if ( path_size <= path_dist_cutoff_ ) return false;
	}
	return true;
}


/////////////////////////////////////////////////////////////////////
Real
CarbonHBondEnergy::res_res_carbon_hbond_one_way(
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Real & bb_bb,
	Real & bb_sc,
	Real & sc_sc
) const {
	Real res_res_energy( 0.0 ), energy( 0.0 );

	// Here we go -- cycle through non-polar hydrogens in don_aa, and all acceptors.
	for ( Size const don_h_atm : don_rsd.Hpos_apolar() ) {
		//  std::cout << "Apolar hydrogen: " << don_rsd.atom_name( don_h_atm ) << " in " <<  don_rsd.name1() << don_rsd.seqpos() << std::endl;

		for ( Size const acc_atm : acc_rsd.accpt_pos() ) {
			//check here whether is backbone

			/// square-distance check outside the call to get_atom_atom_carbon_hbond_energy
			if ( don_rsd.xyz( don_h_atm ).distance_squared( acc_rsd.xyz( acc_atm ) ) > max_dis2_ ) continue;

			if ( get_atom_atom_carbon_hbond_energy( don_h_atm, don_rsd,
					acc_atm, acc_rsd, energy ) ) {
				if ( don_rsd.atom_is_backbone(don_h_atm) && acc_rsd.atom_is_backbone(acc_atm) ) {
					//emap[ch_bond_bb_bb]+=energy;
					bb_bb +=energy;
				} else if ( !don_rsd.atom_is_backbone(don_h_atm) && !acc_rsd.atom_is_backbone(acc_atm) ) {
					//emap[ch_bond_sc_sc]+=energy;
					sc_sc +=energy;
				} else {
					//emap[ch_bond_bb_sc]+=energy;
					bb_sc +=energy;
				}
				res_res_energy += energy;
			}
		}
	}

	return res_res_energy;
}

Real
CarbonHBondEnergy::bb_bb_carbon_hbond_one_way(
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd
) const {

	Real res_res_energy( 0.0 ), energy( 0.0 );

	// Here we go -- cycle through non-polar hydrogens in don_aa, and all acceptors.
	for ( Size const don_h_atm : don_rsd.Hpos_apolar() ) {
		if ( don_h_atm >= don_rsd.first_sidechain_hydrogen() ) continue;
		//  std::cout << "Apolar hydrogen: " << don_rsd.atom_name( don_h_atm ) << " in " <<  don_rsd.name1() << don_rsd.seqpos() << std::endl;

		for ( Size const acc_atm : acc_rsd.accpt_pos() ) {
			if ( acc_atm > acc_rsd.last_backbone_atom() ) continue;

			/// square-distance check outside the call to get_atom_atom_carbon_hbond_energy
			if ( don_rsd.xyz( don_h_atm ).distance_squared( acc_rsd.xyz( acc_atm ) ) > max_dis2_ ) continue;

			get_atom_atom_carbon_hbond_energy( don_h_atm, don_rsd,
				acc_atm,   acc_rsd,
				energy );
			res_res_energy += energy;
		}
	}

	return res_res_energy;
}

Real
CarbonHBondEnergy::sc_bb_carbon_hbond_one_way(
	conformation::Residue const & don_rsd, // sidechain atoms on donor
	conformation::Residue const & acc_rsd  // backbone atoms on acceptor
) const
{

	Real res_res_energy( 0.0 ), energy( 0.0 );

	// Here we go -- cycle through non-polar hydrogens in don_aa, and all acceptors.
	for ( Size const don_h_atm : don_rsd.Hpos_apolar() ) {
		if ( don_h_atm < don_rsd.first_sidechain_hydrogen() ) continue;
		//  std::cout << "Apolar hydrogen: " << don_rsd.atom_name( don_h_atm ) << " in " <<  don_rsd.name1() << don_rsd.seqpos() << std::endl;

		for ( Size const acc_atm : acc_rsd.accpt_pos() ) {
			if ( acc_atm > acc_rsd.last_backbone_atom() ) continue;

			/// square-distance check outside the call to get_atom_atom_carbon_hbond_energy
			if ( don_rsd.xyz( don_h_atm ).distance_squared( acc_rsd.xyz( acc_atm ) ) > max_dis2_ ) continue;

			get_atom_atom_carbon_hbond_energy( don_h_atm, don_rsd,
				acc_atm,   acc_rsd,
				energy );
			res_res_energy += energy;
		}
	}

	return res_res_energy;
}

Real
CarbonHBondEnergy::bb_sc_carbon_hbond_one_way(
	conformation::Residue const & don_rsd, // backbone atoms on donor
	conformation::Residue const & acc_rsd  // sidechain atoms on acceptor
) const {
	Real res_res_energy( 0.0 ), energy( 0.0 );

	// Here we go -- cycle through non-polar hydrogens in don_aa, and all acceptors.
	for ( chemical::AtomIndices::const_iterator
			hnum  = don_rsd.Hpos_apolar().begin(),
			hnume = don_rsd.Hpos_apolar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atm( *hnum );
		if ( don_h_atm >= don_rsd.first_sidechain_hydrogen() ) continue;
		//  std::cout << "Apolar hydrogen: " << don_rsd.atom_name( don_h_atm ) << " in " <<  don_rsd.name1() << don_rsd.seqpos() << std::endl;

		for ( chemical::AtomIndices::const_iterator
				anum  = acc_rsd.accpt_pos().begin(),
				anume = acc_rsd.accpt_pos().end(); anum != anume; ++anum ) {

			Size const acc_atm( *anum );
			if ( acc_atm <= acc_rsd.last_backbone_atom() ) continue;

			/// square-distance check outside the call to get_atom_atom_carbon_hbond_energy
			if ( don_rsd.xyz( don_h_atm ).distance_squared( acc_rsd.xyz( acc_atm ) ) > max_dis2_ ) continue;

			get_atom_atom_carbon_hbond_energy( don_h_atm, don_rsd,
				acc_atm,   acc_rsd,
				energy );
			res_res_energy += energy;
		}
	}

	return res_res_energy;
}

Real
CarbonHBondEnergy::sc_sc_carbon_hbond_one_way(
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd
) const {
	Real res_res_energy( 0.0 ), energy( 0.0 );

	// Here we go -- cycle through non-polar hydrogens in don_aa, and all acceptors.
	for ( chemical::AtomIndices::const_iterator
			hnum  = don_rsd.Hpos_apolar().begin(),
			hnume = don_rsd.Hpos_apolar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atm( *hnum );
		if ( don_h_atm < don_rsd.first_sidechain_hydrogen() ) continue;
		//  std::cout << "Apolar hydrogen: " << don_rsd.atom_name( don_h_atm ) << " in " <<  don_rsd.name1() << don_rsd.seqpos() << std::endl;

		for ( chemical::AtomIndices::const_iterator
				anum  = acc_rsd.accpt_pos().begin(),
				anume = acc_rsd.accpt_pos().end(); anum != anume; ++anum ) {

			Size const acc_atm( *anum );
			if ( acc_atm <= acc_rsd.last_backbone_atom() ) continue;

			/// square-distance check outside the call to get_atom_atom_carbon_hbond_energy
			if ( don_rsd.xyz( don_h_atm ).distance_squared( acc_rsd.xyz( acc_atm ) ) > max_dis2_ ) continue;

			get_atom_atom_carbon_hbond_energy( don_h_atm, don_rsd, acc_atm, acc_rsd, energy );
			res_res_energy += energy;
		}
	}

	return res_res_energy;
}


void
CarbonHBondEnergy::res_res_carbon_hbond_derivs_one_way(
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & don_atom_derivs,
	utility::vector1< DerivVectorPair > & acc_atom_derivs
) const {
	bool const eval_bbbb( weights[ ch_bond ] != 0.0 || weights[ ch_bond_bb_bb ] != 0.0 );
	bool const eval_bbsc( weights[ ch_bond ] != 0.0 || weights[ ch_bond_bb_sc ] != 0.0 );
	bool const eval_scsc( weights[ ch_bond ] != 0.0 || weights[ ch_bond_bb_sc ] != 0.0 ); //Mistake on this line? Parin S. (sripakpa@stanford.edu) Jan 11, 2012
	bool const eval_bb( eval_bbsc || eval_bbbb );
	bool const eval_sc( eval_bbsc || eval_scsc );

	// Here we go -- cycle through non-polar hydrogens in don_aa, and all acceptors.
	for ( Size const don_h_atm : don_rsd.Hpos_apolar() ) {
		bool const don_h_bb( don_rsd.atom_is_backbone( don_h_atm ) );
		if ( don_h_bb ) {
			if ( ! eval_bb ) continue;
		} else {
			if ( ! eval_sc ) continue;
		}

		//  std::cout << "Apolar hydrogen: " << don_rsd.atom_name( don_h_atm ) << " in " <<  don_rsd.name1() << don_rsd.seqpos() << std::endl;

		for ( Size const acc_atm : acc_rsd.accpt_pos() ) {
			bool const acc_bb( acc_rsd.atom_is_backbone( acc_atm ) );

			/// skip the derivative evaluation if we have a zero weight for the interaction we're examining
			if ( acc_bb ) {
				if ( ! eval_bb ) continue;
				if ( don_h_bb ) {
					if ( ! eval_bbbb ) continue;
				} else {
					if ( ! eval_bbsc ) continue;
				}
			} else { // acceptor is sidechain
				if ( ! eval_sc ) continue;
				if ( don_h_bb ) {
					if ( ! eval_bbsc ) continue;
				} else {
					if ( ! eval_scsc ) continue;
				}
			}

			/// square-distance check outside the call to get_atom_atom_carbon_hbond_energy
			if ( don_rsd.xyz( don_h_atm ).distance_squared( acc_rsd.xyz( acc_atm ) ) > max_dis2_ ) continue;

			Vector f2; Real energy;
			if ( get_atom_atom_carbon_hbond_energy( don_h_atm, don_rsd,
					acc_atm, acc_rsd, energy, true, f2 ) ) {
				//// 1. scale the f2 vector by the weight that applies to this interaction
				if ( don_h_bb && acc_bb ) {
					f2 *= weights[ ch_bond ] + weights[ ch_bond_bb_bb ];
				} else if ( ! don_h_bb && ! acc_bb ) {
					f2 *= weights[ ch_bond ] + weights[ ch_bond_sc_sc ];
				} else {
					f2 *= weights[ ch_bond ] + weights[ ch_bond_bb_sc ];
				}

				if ( use_orientation_dep_rna_ch_o_bonds(don_rsd, acc_rsd) ) {
					//The standard code doesn't appear to work properly for the RNA case (i.e. fail the numerical_derivative_check() test)
					//I am including a special version for RNA. Parin S. (sripakpa@stanford.edu).  Jan 11, 2012

					Vector const f1 = cross( f2, acc_rsd.xyz( acc_atm ) );

					don_atom_derivs[ don_h_atm ].f2() += f2;
					don_atom_derivs[ don_h_atm ].f1() += f1 ;

					acc_atom_derivs[ acc_atm ].f2() -= f2;
					acc_atom_derivs[ acc_atm ].f1() -= f1;

				} else {

					/// 2. f2 is the force vector on the hydrogen; compute f1 by taking the cross product
					/// with the coordinate of the acceptor atom
					don_atom_derivs[ don_h_atm ].f2() += f2;
					Vector f1_H = cross( f2, acc_rsd.xyz( acc_atm )  );
					don_atom_derivs[ don_h_atm ].f1() += f1_H;

					/// 3. Since f2 is the force vector on the hydrogen, negate it to get the force
					/// vector on the acceptor; compute f1 by taking the cross product
					/// with the coordinate of the hydrogen atom
					f2 *= -1;
					acc_atom_derivs[ acc_atm ].f2() += f2;
					Vector f1_Acc = cross( f2, don_rsd.xyz( don_h_atm ) );
					acc_atom_derivs[ acc_atm ].f1() += f1_Acc;
				}

				/*std::cout << "  chbond deriv: " << energy << " " << don_rsd.seqpos() << " " << don_h_atm << " " << don_rsd.atom_name(don_h_atm) << " "
				<< acc_rsd.seqpos() << " " << acc_atm <<  " " << acc_rsd.atom_name( acc_atm ) << " " << don_h_bb << " " << acc_bb << std::endl;
				std::cout << "   don f1: " << don_atom_derivs[ don_h_atm ].f1().x() << " " << don_atom_derivs[ don_h_atm ].f1().y() << " " << don_atom_derivs[ don_h_atm ].f1().z() << std::endl;
				std::cout << "   don f2: " << don_atom_derivs[ don_h_atm ].f2().x() << " " << don_atom_derivs[ don_h_atm ].f2().y() << " " << don_atom_derivs[ don_h_atm ].f2().z() << std::endl;
				std::cout << "   acc f1: " << acc_atom_derivs[ acc_atm ].f1().x() << " " << acc_atom_derivs[ acc_atm ].f1().y() << " " << acc_atom_derivs[ acc_atm ].f1().z() << std::endl;
				std::cout << "   acc f2: " << acc_atom_derivs[ acc_atm ].f2().x() << " " << acc_atom_derivs[ acc_atm ].f2().y() << " " << acc_atom_derivs[ acc_atm ].f2().z() << std::endl;*/

			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
CarbonHBondEnergy::get_atom_atom_carbon_hbond_energy(
	Size const don_h_atm,
	conformation::Residue const & don_rsd,
	Size const acc_atm,
	conformation::Residue const & acc_rsd,
	Real & energy,
	bool const update_deriv /*= false*/,
	Vector & f2 /*=ZERO_VECTOR*/ // Only computes the force vector between don and acceptor -- calling code must compute f1
) const {
	energy = 0.0;
	f2 = 0.0;

	Size const don_atm( don_rsd.atom_base( don_h_atm ) );
	Size const base_atm( acc_rsd.atom_base( acc_atm ) );

	if ( !path_distance_OK( don_rsd, acc_rsd, don_atm, acc_atm ) ) return false; // Look more than four atoms away.

	////////////////////////////////////////////////////////////////////////////
	//No virtual atom seem to reach to point, probably becuase virtual_atom is neither Hpos_apolar nor accpt_atom.
	//But should still have these checks here just to be safe. [Parin S. (sripakpa@stanford.edu) Jan 11, 2012]
	if ( acc_rsd.is_virtual( acc_atm ) ) return false;
	if ( don_rsd.is_virtual( don_atm ) ) return false;
	if ( don_rsd.is_virtual( don_h_atm) ) return false;

	////////////////////////////////////////////////////////////////////////////

	Vector const & don_h_atm_xyz( don_rsd.atom( don_h_atm ).xyz() );
	Vector const & don_atm_xyz( don_rsd.atom( don_atm ).xyz() );
	Vector const & acc_atm_xyz( acc_rsd.atom( acc_atm ).xyz() );
	Vector const & base_atm_xyz( acc_rsd.atom( base_atm ).xyz() );


	Vector H_A_vector = acc_atm_xyz - don_h_atm_xyz;
	Vector D_H_vector = don_h_atm_xyz - don_atm_xyz;
	Vector B_A_vector = acc_atm_xyz - base_atm_xyz;


	if ( use_orientation_dep_rna_ch_o_bonds(don_rsd, acc_rsd) ) {

		Vector const r_H_A( acc_atm_xyz  - don_h_atm_xyz );
		Vector const z_D_H( ( don_h_atm_xyz - don_atm_xyz ).normalize() );

		energy = carbon_hbond_potential_.get_potential_RNA( r_H_A, z_D_H, update_deriv, f2 );
		// for some reason, get_potential_RNA() returns the derivative for the acceptor,
		// whereas get_potential() returns the derivative for the hydrogen.  Since this function
		// should return the derivative for the hydrogen, multiply f2 by -1.
		f2 *= -1;

		if ( verbose_ && energy < -0.05 ) {
			Real const angle_DH_A = numeric::conversions::degrees( angle_radians( don_atm_xyz, don_h_atm_xyz, acc_atm_xyz ) );
			tr <<"CHbond [RNA]: "<< don_rsd.name1() << I(3,don_rsd.seqpos())<<
				" atom "<< don_rsd.atom_name( don_h_atm )<< " [ " <<
				don_rsd.atom_name( don_atm) <<
				" ] bonded to acc_res " <<
				acc_rsd.name1()<< I(3, acc_rsd.seqpos()) <<
				" atom "<< acc_rsd.atom_name( acc_atm ) <<
				" with energy "<< F(8,3,energy) << " [" << F(8,3,H_A_vector.length()) << " Angstroms; "
				<< angle_DH_A << " degrees ]" << std::endl;
		}

	} else {

		energy = carbon_hbond_potential_.get_potential(
			don_rsd.atom_type_index( don_atm ), H_A_vector, D_H_vector, B_A_vector, update_deriv, f2);

		if ( verbose_ && energy < -0.05 ) {
			tr <<"CHbond: "<< don_rsd.name1() << I(3,don_rsd.seqpos())<<
				" atom "<< don_rsd.atom_name( don_h_atm )<< " [ " <<
				don_rsd.atom_name( don_atm) <<
				" ] bonded to acc_res " <<
				acc_rsd.name1()<< I(3, acc_rsd.seqpos()) <<
				" atom "<< acc_rsd.atom_name( acc_atm ) <<
				" with energy "<< F(8,3,energy) << " [" << F(8,3,H_A_vector.length()) << "]" << std::endl;
		}
	}

	return true;
}

void
CarbonHBondEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const &,
	pose::Pose const &, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {
	res_res_carbon_hbond_derivs_one_way( rsd1, rsd2, weights, r1_atom_derivs, r2_atom_derivs );
	res_res_carbon_hbond_derivs_one_way( rsd2, rsd1, weights, r2_atom_derivs, r1_atom_derivs );
}

void
CarbonHBondEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const &,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const {
	/// intra-residue carbon hbond derivatives; one call to
	res_res_carbon_hbond_derivs_one_way( rsd, rsd, weights, atom_derivs, atom_derivs );
}


////////////////////////////////////////////////////////////////////////////////////////////
Distance
CarbonHBondEnergy::atomic_interaction_cutoff() const
{
	//return hbonds::MAX_R + 1.35; // MAGIC NUMBER
	return 4.35; //MAX_R is no longer a constant
}

////////////////////////////////////////////////////////////////////////////////////////////
bool
CarbonHBondEnergy::defines_intrares_energy( EnergyMap const & /*weights*/ ) const
{
	return true;
}

void
CarbonHBondEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap & emap
) const {
	Real bb_bb(0.0);
	Real bb_sc(0.0);
	Real sc_sc(0.0);
	Real const res_energy = res_res_carbon_hbond_one_way( rsd, rsd, bb_bb, bb_sc, sc_sc );
	emap[ ch_bond ] += res_energy;
	emap[ ch_bond_bb_bb ] += bb_bb;
	emap[ ch_bond_bb_sc ] += bb_sc;
	emap[ ch_bond_sc_sc ] += sc_sc;

	// std::cout << "INTRARES" << rsd.seqpos() << " " << res_energy << " " << emap[ch_bond] << std::endl;
}

/// @brief CarbonHBondEnergy is not context sensitive
void
CarbonHBondEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const {
	/*nothing*/
}

core::Size
CarbonHBondEnergy::version() const
{
	return 1; // Initial versioning
}

////////////////////////////////////////////////////////////////////////////////////////////
bool
CarbonHBondEnergy::use_orientation_dep_rna_ch_o_bonds(conformation::Residue const & don_rsd, conformation::Residue const & acc_rsd) const
{

	return ( orientation_dep_rna_ch_o_bonds_ && don_rsd.is_RNA() && acc_rsd.is_RNA() );

}
////////////////////////////////////////////////////////////////////////////////////////////

} // carbon_hbonds
} // scoring
} // core
