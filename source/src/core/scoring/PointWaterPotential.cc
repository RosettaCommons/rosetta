// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/PointWaterPotential.hh
/// @brief  Statistical point water potential
/// @author Frank DiMaio

// Unit Headers
#include <core/scoring/PointWaterPotential.hh>

// Package Headers
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/AtomType.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <basic/basic.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/deriv/angle_deriv.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/angle.functions.hh>
#include <numeric/interpolation/spline/BicubicSpline.hh>

#include <core/scoring/DerivVectorPair.hh>
#include <core/chemical/AA.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/database/open.hh>

namespace core {
namespace scoring {

PointWaterPotential::PointWaterPotential() {
	read_pointwater_tables( );
}

//
//

core::Real
PointWaterPotential::eval_pointwater_score(
	core::chemical::AA const res_aa,
	core::conformation::Residue const &res,
	core::Vector O_pos
) const {
	using basic::subtract_degree_angles;

	core::Real score = 0;

	if ( res.is_protein() ) {
		for ( int i=1; i<=(int) potential_ids_.size(); ++i ) {
			core::chemical::AA tgt_aa = potential_ids_[i].aa;
			if ( tgt_aa != core::chemical::aa_unk && tgt_aa != res_aa ) continue;

			if ( !res.has(potential_ids_[i].name1) || !res.has(potential_ids_[i].name2) ) continue;

			core::Size atm1_idx = res.atom_index(potential_ids_[i].name1);
			core::Size atm2_idx = res.atom_index(potential_ids_[i].name2);
			core::Vector atm1 = res.atom( atm1_idx ).xyz();
			core::Vector atm2 = res.atom( atm2_idx ).xyz();
			core::Real d = (atm1 - O_pos).length();
			if ( d>6 ) continue;
			core::Real theta = numeric::angle_degrees( atm2,atm1,O_pos );
			core::Real score_i = aa_potentials_[i].F(d,theta);
			score += score_i;
		}
	} else {
		for ( core::chemical::AtomIndices::const_iterator
				anum  = res.accpt_pos().begin(),
				anume = res.accpt_pos().end(); anum != anume; ++anum ) {
			core::Size const aatm( *anum );
			core::Size const abase1( res.atom_base( aatm ) );
			core::Vector atm1 = res.atom( aatm ).xyz();
			core::Vector atm2 = res.atom( abase1 ).xyz();
			core::Real d = (atm1 - O_pos).length();
			if ( d>6 ) continue;
			core::Real theta = numeric::angle_degrees( atm2,atm1,O_pos );
			core::Real score_i = 0;
			core::chemical::Hybridization hybrid = res.atom_type( aatm ).hybridization();
			if ( hybrid == core::chemical::SP2_HYBRID ) {
				score_i = lig_potential_acc_sp2_.F(d,theta);
			} else if ( hybrid == core::chemical::SP3_HYBRID ) {
				score_i = lig_potential_acc_sp3_.F(d,theta);
			} // TODO: add ring
			score += score_i;
		}

		for ( core::chemical::AtomIndices::const_iterator
				hnum  = res.Hpos_polar().begin(),
				hnume = res.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			core::Size const hatm( *hnum );
			core::Size const abase1( res.atom_base( hatm ) );
			core::Vector atm1 = res.atom( hatm ).xyz();
			core::Vector atm2 = res.atom( abase1 ).xyz();
			core::Real d = (atm1 - O_pos).length();
			if ( d>6 ) continue;
			core::Real theta = numeric::angle_degrees( atm2,atm1,O_pos );
			core::Real score_i = lig_potential_don_.F(d,theta);
			score += score_i;

		}
	}

	return score;
}


Real
PointWaterPotential::eval_pointwater_derivs(
	core::chemical::AA const res_aa1,
	core::conformation::Residue const &res,
	core::Vector O_pos,
	utility::vector1< DerivVectorPair > & r_atom_derivs,
	utility::vector1< DerivVectorPair > & O_derivs,
	core::Real wt
) const {
	using basic::subtract_degree_angles;

	core::Real score = 0;
	if ( res.is_protein() ) {
		for ( int i=1; i<=(int) potential_ids_.size(); ++i ) {
			AA tgt_aa = potential_ids_[i].aa;
			if ( tgt_aa != core::chemical::aa_unk && tgt_aa != res_aa1 ) continue;

			if ( !res.has(potential_ids_[i].name1) || !res.has(potential_ids_[i].name2) ) continue;

			core::Size atm1_idx = res.atom_index(potential_ids_[i].name1);
			core::Size atm2_idx = res.atom_index(potential_ids_[i].name2);

			if ( atm1_idx == 0 || atm2_idx == 0 ) continue;

			core::Vector atm1 = res.atom( atm1_idx ).xyz();
			core::Vector atm2 = res.atom( atm2_idx ).xyz();

			core::Real d = (atm1 - O_pos).length();
			core::Real theta = numeric::angle_degrees( atm2,atm1,O_pos );

			if ( d>6 ) continue;

			core::Real score_i = aa_potentials_[i].F(d,theta);
			core::Real dE_dd = wt * aa_potentials_[i].dFdx(d,theta);
			core::Real dE_dtheta = wt * numeric::conversions::degrees( aa_potentials_[i].dFdy(d,theta) );

			// assign to atoms [ASSUME O is ATOM 1]
			// LENGTH
			Vector f1(0.0), f2(0.0);
			numeric::deriv::distance_f1_f2_deriv( atm1, O_pos, d, f1, f2 );
			r_atom_derivs[ atm1_idx ].f1() += dE_dd * f1;
			r_atom_derivs[ atm1_idx ].f2() += dE_dd * f2;
			O_derivs[ 1 ].f1() -= dE_dd * f1;
			O_derivs[ 1 ].f2() -= dE_dd * f2;

			// ANGLE
			numeric::deriv::angle_p1_deriv( atm2, atm1, O_pos, theta, f1, f2 );
			r_atom_derivs[ atm2_idx ].f1() += dE_dtheta * f1;
			r_atom_derivs[ atm2_idx ].f2() += dE_dtheta * f2;

			numeric::deriv::angle_p2_deriv( atm2, atm1, O_pos, theta, f1, f2 );
			r_atom_derivs[ atm1_idx ].f1() += dE_dtheta * f1;
			r_atom_derivs[ atm1_idx ].f2() += dE_dtheta * f2;

			numeric::deriv::angle_p1_deriv(O_pos, atm1, atm2, theta, f1, f2 );
			O_derivs[ 1 ].f1() += dE_dtheta * f1;
			O_derivs[ 1 ].f2() += dE_dtheta * f2;

			score += score_i;
		}
	} else {
		for ( core::chemical::AtomIndices::const_iterator
				anum  = res.accpt_pos().begin(),
				anume = res.accpt_pos().end(); anum != anume; ++anum ) {
			core::Size const aatm( *anum );
			core::Size const abase1( res.atom_base( aatm ) );
			core::Vector atm1 = res.atom( aatm ).xyz();
			core::Vector atm2 = res.atom( abase1 ).xyz();
			core::Real d = (atm1 - O_pos).length();
			if ( d>6 ) continue;
			core::Real theta = numeric::angle_degrees( atm2,atm1,O_pos );
			core::Real score_i=0, dE_dd=0, dE_dtheta=0;
			core::chemical::Hybridization hybrid = res.atom_type( aatm ).hybridization();
			if ( hybrid == core::chemical::SP2_HYBRID ) {
				score_i = lig_potential_acc_sp2_.F(d,theta);
				dE_dd = wt * lig_potential_acc_sp2_.dFdx(d,theta);
				dE_dtheta = wt * numeric::conversions::degrees( lig_potential_acc_sp2_.dFdy(d,theta) );
			} else if ( hybrid == core::chemical::SP3_HYBRID ) {
				score_i = lig_potential_acc_sp3_.F(d,theta);
				dE_dd = wt * lig_potential_acc_sp3_.dFdx(d,theta);
				dE_dtheta = wt * numeric::conversions::degrees( lig_potential_acc_sp3_.dFdy(d,theta) );
			} // TODO: add ring


			// assign to atoms [ASSUME O is ATOM 1]
			// LENGTH
			Vector f1(0.0), f2(0.0);
			numeric::deriv::distance_f1_f2_deriv( atm1, O_pos, d, f1, f2 );
			r_atom_derivs[ aatm ].f1() += dE_dd * f1;
			r_atom_derivs[ aatm ].f2() += dE_dd * f2;
			O_derivs[ 1 ].f1() -= dE_dd * f1;
			O_derivs[ 1 ].f2() -= dE_dd * f2;

			// ANGLE
			numeric::deriv::angle_p1_deriv( atm2, atm1, O_pos, theta, f1, f2 );
			r_atom_derivs[ abase1 ].f1() += dE_dtheta * f1;
			r_atom_derivs[ abase1 ].f2() += dE_dtheta * f2;

			numeric::deriv::angle_p2_deriv( atm2, atm1, O_pos, theta, f1, f2 );
			r_atom_derivs[ aatm ].f1() += dE_dtheta * f1;
			r_atom_derivs[ aatm ].f2() += dE_dtheta * f2;

			numeric::deriv::angle_p1_deriv(O_pos, atm1, atm2, theta, f1, f2 );
			O_derivs[ 1 ].f1() += dE_dtheta * f1;
			O_derivs[ 1 ].f2() += dE_dtheta * f2;

			score += score_i;
		}

		for ( core::chemical::AtomIndices::const_iterator
				hnum  = res.Hpos_polar().begin(),
				hnume = res.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			core::Size const hatm( *hnum );
			core::Size const abase1( res.atom_base( hatm ) );
			core::Vector atm1 = res.atom( hatm ).xyz();
			core::Vector atm2 = res.atom( abase1 ).xyz();
			core::Real d = (atm1 - O_pos).length();
			if ( d>6 ) continue;
			core::Real theta = numeric::angle_degrees( atm2,atm1,O_pos );
			core::Real score_i = lig_potential_don_.F(d,theta);
			core::Real dE_dd = wt * lig_potential_acc_sp2_.dFdx(d,theta);
			core::Real dE_dtheta = wt * numeric::conversions::degrees( lig_potential_acc_sp2_.dFdy(d,theta) );

			// assign to atoms [ASSUME O is ATOM 1]
			// LENGTH
			Vector f1(0.0), f2(0.0);
			numeric::deriv::distance_f1_f2_deriv( atm1, O_pos, d, f1, f2 );
			r_atom_derivs[ hatm ].f1() += dE_dd * f1;
			r_atom_derivs[ hatm ].f2() += dE_dd * f2;
			O_derivs[ 1 ].f1() -= dE_dd * f1;
			O_derivs[ 1 ].f2() -= dE_dd * f2;

			// ANGLE
			numeric::deriv::angle_p1_deriv( atm2, atm1, O_pos, theta, f1, f2 );
			r_atom_derivs[ abase1 ].f1() += dE_dtheta * f1;
			r_atom_derivs[ abase1 ].f2() += dE_dtheta * f2;

			numeric::deriv::angle_p2_deriv( atm2, atm1, O_pos, theta, f1, f2 );
			r_atom_derivs[ hatm ].f1() += dE_dtheta * f1;
			r_atom_derivs[ hatm ].f2() += dE_dtheta * f2;

			numeric::deriv::angle_p1_deriv(O_pos, atm1, atm2, theta, f1, f2 );
			O_derivs[ 1 ].f1() += dE_dtheta * f1;
			O_derivs[ 1 ].f2() += dE_dtheta * f2;

			score += score_i;
		}
	}

	return score;
}


///
/// load tables
void
PointWaterPotential::read_pointwater_tables( ) {
	using namespace core::chemical;
	using namespace basic::options;

	std::string pwat_database = option[OptionKeys::corrections::water::pointwater_params];
	//  std::cout << "POINTWATER database location = " << pwat_database << std::endl;

	///////
	// aa-specific
	potential_ids_.push_back( tableID( aa_unk, "H", "N")  ); // use aa_unk to mean all AAs
	potential_ids_.push_back( tableID( aa_unk, "O", "C")  ); // use aa_unk to mean all AAs
	potential_ids_.push_back( tableID( aa_arg, "HE", "NE")  );
	potential_ids_.push_back( tableID( aa_arg, "1HH1", "NH1")  );
	potential_ids_.push_back( tableID( aa_arg, "2HH1", "NH1")  );
	potential_ids_.push_back( tableID( aa_arg, "1HH2", "NH2")  );
	potential_ids_.push_back( tableID( aa_arg, "2HH2", "NH2")  );
	potential_ids_.push_back( tableID( aa_asn, "1HD2", "ND2")  );
	potential_ids_.push_back( tableID( aa_asn, "2HD2", "ND2")  );
	potential_ids_.push_back( tableID( aa_asn, "OD1", "CG")  );
	potential_ids_.push_back( tableID( aa_asp, "OD1", "CG")  );
	potential_ids_.push_back( tableID( aa_asp, "OD2", "CG")  );
	potential_ids_.push_back( tableID( aa_gln, "1HE2", "NE2")  );
	potential_ids_.push_back( tableID( aa_gln, "2HE2", "NE2")  );
	potential_ids_.push_back( tableID( aa_gln, "OE1", "CG")  );
	potential_ids_.push_back( tableID( aa_glu, "OE1", "CD")  );
	potential_ids_.push_back( tableID( aa_glu, "OE2", "CD")  );
	potential_ids_.push_back( tableID( aa_his, "HD1", "ND1")  );
	potential_ids_.push_back( tableID( aa_his, "HE2", "NE2")  );
	potential_ids_.push_back( tableID( aa_lys, "1HZ", "NZ")  );
	potential_ids_.push_back( tableID( aa_lys, "2HZ", "NZ")  );
	potential_ids_.push_back( tableID( aa_lys, "3HZ", "NZ")  );
	potential_ids_.push_back( tableID( aa_ser, "OG", "CB")  );
	potential_ids_.push_back( tableID( aa_thr, "OG1", "CB")  );
	potential_ids_.push_back( tableID( aa_trp, "HE1", "NE1")  );
	potential_ids_.push_back( tableID( aa_tyr, "OH", "CZ")  );

	aa_potentials_.resize( potential_ids_.size() );

	for ( int i=1; i<=(int)potential_ids_.size(); ++i ) {
		std::string aaname = (potential_ids_[i].aa != core::chemical::aa_unk) ?
			core::chemical::name_from_aa( potential_ids_[i].aa ) : "ALL";
		std::string name = "E_"+aaname+"_"+potential_ids_[i].name1+"_"+potential_ids_[i].name2+".dat";

		utility::io::izstream instream;
		//    basic::database::open( instream, "scoring/score_functions/pointwater/"+name);
		basic::database::open( instream, "scoring/score_functions/pointwater/"+pwat_database+"/"+name);

		ObjexxFCL::FArray2D<Real> pot_i;
		read_table_from_stream( instream, pot_i );
		dampen_longrange_energies( pot_i );
		setup_interpolation( pot_i, aa_potentials_[i] );
	}

	///////
	// ligand (general)
	{
		utility::io::izstream instream;
		basic::database::open( instream, "scoring/score_functions/pointwater/"+pwat_database+"/E_LIG_DON.dat");
		ObjexxFCL::FArray2D<Real> pot_i;
		read_table_from_stream( instream, pot_i );
		dampen_longrange_energies( pot_i );
		setup_interpolation( pot_i, lig_potential_don_ );
	}
	{
		utility::io::izstream instream;
		basic::database::open( instream, "scoring/score_functions/pointwater/"+pwat_database+"/E_LIG_ACC_SP2.dat");
		ObjexxFCL::FArray2D<Real> pot_i;
		read_table_from_stream( instream, pot_i );
		dampen_longrange_energies( pot_i );
		setup_interpolation( pot_i, lig_potential_acc_sp2_ );
	}
	{
		utility::io::izstream instream;
		basic::database::open( instream, "scoring/score_functions/pointwater/"+pwat_database+"/E_LIG_ACC_SP3.dat");
		ObjexxFCL::FArray2D<Real> pot_i;
		read_table_from_stream( instream, pot_i );
		dampen_longrange_energies( pot_i );
		setup_interpolation( pot_i, lig_potential_acc_sp3_ );
	}
}

void
PointWaterPotential::setup_interpolation(
	ObjexxFCL::FArray2D< Real > & x,
	numeric::interpolation::spline::BicubicSpline  &sx
) {
	using namespace numeric;
	using namespace numeric::interpolation::spline;

	// 24 x 24 with one cell padding all around ==> 26 x 26
	MathMatrix< Real > energy_vals( 26, 26 );
	for ( Size jj = 0; jj < 26; ++jj ) {
		for ( Size kk = 0; kk < 26; ++kk ) {
			energy_vals( jj, kk ) = x(jj+1,kk+1);
		}
	}

	// i index --> d
	// j index --> theta
	BorderFlag periodic_boundary[2] = { e_FirstDer, e_Natural };
	Real start_vals[2] = { -0.125, -3.75 };
	Real deltas[2] = {0.25, 7.5};         // 7.5 deg, 0.125A steps
	bool lincont[2] = {false,false};      //meaningless argument for a bicubic spline with periodic boundary conditions
	std::pair< Real, Real > unused[2];
	unused[0] = std::make_pair( 0.0, 0.0 );
	unused[1] = std::make_pair( 0.0, 0.0 );
	sx.train( periodic_boundary, start_vals, deltas, energy_vals, lincont, unused );
}

void
PointWaterPotential::read_table_from_stream(
	utility::io::izstream &stream,
	ObjexxFCL::FArray2D< Real > &table
) {
	table.dimension(26,26);

	int rowctr=2;
	core::Real buff;
	std::string line;
	while ( getline( stream, line ) ) {
		std::istringstream l(line);
		for ( int j=2; j<=25; ++j ) {
			l >> buff;
			table( rowctr, j) = buff;
		}
		table( rowctr, 1)  = table( rowctr, 2);
		table( rowctr, 26) = table( rowctr, 25);
		rowctr++;
	}
	for ( int j=1; j<=26; ++j ) {
		table(1,j) = table(2,j);
		table(26,j) = table(25,j);
	}
}

void
PointWaterPotential::dampen_longrange_energies(
	ObjexxFCL::FArray2D< Real > &table
) {
	// cap
	//for (int i=1; i<=26; ++i) {
	// for (int j=1; j<=26; ++j) {
	//  if (table(i,j)>0.0) table(i,j) = 0.0;
	// }
	//}

	// fade from full weight at 4A to 0 at 5A
	for ( int i=1; i<=26; ++i ) {
		core::Real scale=1;

		if ( i==17 ) scale = 0.75;
		if ( i==18 ) scale = 0.5;
		if ( i==19 ) scale = 0.25;
		if ( i>=20 ) scale = 0.0;

		for ( int j=1; j<=26; ++j ) {
			table(i,j) *= scale;
		}
	}
}

}
}
