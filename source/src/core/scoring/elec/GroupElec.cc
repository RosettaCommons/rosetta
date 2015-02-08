// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/GroupElec.cc
/// @brief  Electrostatic energy with a distance-dependant dielectric
/// @author Hahnbeom Park


// Unit headers
#include <core/scoring/elec/GroupElec.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>

// Package headers
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/ResidueNeighborList.hh>
//#include <core/scoring/hbonds/HBondSet.hh>

#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// Project headers
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/conformation/RotamerSetBase.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// C++ headers
#include <iostream>

static basic::Tracer TR("core.scoring.elec.GroupElec");

namespace core {
namespace scoring {
namespace elec {

using namespace core;

GroupElec::GroupElec( methods::EnergyMethodOptions const & options ):
	fade_type_( options.grpelec_fade_type() ),
	fade_param1_( options.grpelec_fade_param1() ),
	fade_param2_( options.grpelec_fade_param2() ),
	fade_hbond_( options.grpelec_fade_hbond() ),
	group_file_( options.elec_group_file() ),
	grp_cpfxn_( options.grp_cpfxn() )
{}

GroupElec::GroupElec( GroupElec const & src ): ReferenceCount(),
	fade_type_( src.fade_type_ ),
	fade_param1_( src.fade_param1_ ),
	fade_param2_( src.fade_param2_ ),
	fade_hbond_( src.fade_hbond_ ), //option[ score::elec_fade_hbond ]();
	group_file_( src.group_file_ ), //option[ score::elec_group_file ]();
	grp_cpfxn_( src.grp_cpfxn_ )
{}

GroupElec::~GroupElec(){}

void
GroupElec::initialize( etable::coulomb::Coulomb const &coulomb )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
	coulomb_ = coulomb.clone();

	build_groupinfo( group_file_ );
	if( option[ score::elec_group_extrafile ].user() ){
		std::string extra_file = option[ score::elec_group_extrafile ]();
		build_groupinfo( extra_file, true );
	}

	cpfxn_weight_.resize( 3 ); cpfxn_weight_[1] = 0.0; cpfxn_weight_[2] = 0.2; cpfxn_weight_[3] = 1.0;
	if( option[ score::grpelec_cpfxn_weight ].user() ) cpfxn_weight_ = option[ score::grpelec_cpfxn_weight ]();

}


ResElecGroup const &
GroupElec::get_group( core::chemical::ResidueType const &rsdtype ) const 
{
  std::map< std::string const, ResElecGroup >::const_iterator it 
    = rsdgrps_.find( rsdtype.name() );

  if( it == rsdgrps_.end() ){
    TR.Debug << "Building extra group on " << rsdtype.name() << std::endl;

    // otherwise assign new group
    ResElecGroup resgrp;
    resgrp.resize( 0 );

    // if one wants to use whole residue as group
    /*
      core::Size nheavy( 0.0 );
      utility::vector1< Size > resgrp;
      for( core::Size iatm = 1; iatm <= rsdtype->natoms(); ++iatm ){
      resgrp.push_back( iatm );
      if( !rsdtype->atom_is_hydrogen( iatm ) && !rsdtype->is_virtual( iatm ) ) nheavy++;
      }
      grp.grps.push_back( resgrp );
      grp.nheavy.push_back( nheavy );
    */

    // or, each atom as group
    for( core::Size iatm = 1; iatm <= rsdtype.natoms(); ++iatm ){
      if( rsdtype.is_virtual( iatm ) ) continue;

      ElecGroup grp;
      grp.atms.push_back( iatm );
      grp.comatms.push_back( iatm );
      grp.n_acceptor = 0.0;
      grp.n_donor = 0.0;
      grp.qeps = 0.0;
      resgrp.push_back( grp );
    }

    rsdgrps_[ rsdtype.name() ] = resgrp;
		it = rsdgrps_.find( rsdtype.name() );
  }
  //return rsdgrps_.at( rsdtype.name() );
	return it->second;
}

void
GroupElec::build_groupinfo( std::string const group_file,
														bool const extra )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// I will eventually move this to energy_method_options; for now here for convenience
	utility::vector1< Real > qeps = option[ score::grpelec_max_qeps ]();

	utility::io::izstream instream;
	std::string line;

	chemical::ResidueTypeSetCOP rsdtypeset =
		chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	TR << "Reading group information from: " << group_file << std::endl;
	if( extra ){
		instream.open( group_file );
		if( !instream )
			utility_exit_with_message( "No extra group_def file found: "+group_file );

	} else {
		basic::database::open( instream, group_file );
	}

	std::string resname;
	bool skip_residue( false );

	while( instream ){
		getline( instream, line );
		std::string s1("");
		std::istringstream linestream( line );
		linestream >> s1 ;

		if( s1.compare( "RESI" ) == 0 ){
			skip_residue = false;
			std::string s2("");
			linestream >> s2;
			resname = s2;

			if( !rsdtypeset->has_name( s2 ) ){
				skip_residue = true;
				continue;
			}
			//rsdtype = &rsdtypeset->name_map( resname );

			TR.Debug << "Adding group info on residue " << resname << std::endl;

		} else if( s1.compare( "GROUP" ) == 0 && !skip_residue ){
			TR.Debug << "Group: ";

			core::chemical::ResidueType const rsdtype = rsdtypeset->name_map( resname );

			ElecGroup grp;
			std::string atmname;
		  Size n_donor( 0 ), n_acceptor( 0 ), polar_type( 0 );

			linestream >> polar_type;
			if( polar_type >= 0 && polar_type < 4 ){
				grp.qeps = qeps[ polar_type+1 ];
			} else {
				TR << "PolarType index exceeds boundary, skip! line: " << line << std::endl;
			}

			while( linestream >> atmname ){
				if( !rsdtype.has( atmname ) )
					utility_exit_with_message( "No atom found in rsd/atm "+resname+"/"+atmname );
				Size const iatm = rsdtype.atom_index( atmname );
				grp.atms.push_back( iatm );

				// assign if any atom among group is donor/acceptor
				if( rsdtype.atom_type( iatm ).is_donor() ) n_donor ++;
				if( rsdtype.atom_type( iatm ).is_acceptor() ) n_acceptor ++;

				TR.Debug << " " << atmname;
			}

			// COM atms - take only Hpolar or acceptor for COM atoms
			if( n_donor > 0 || n_acceptor > 0 ){
				for( core::Size iatm = 1; iatm <= grp.atms.size(); ++iatm ){
					core::Size atmno = grp.atms[iatm];

					if( rsdtype.atom_is_polar_hydrogen( atmno ) ||
							rsdtype.heavyatom_is_an_acceptor( atmno ) )
						grp.comatms.push_back( atmno );
				} 
			}

			// otherwise fill in heavy atoms
			// why is Npro donor anyway? that assignment was causing weird behavior, though fixed now
			if( grp.comatms.size() == 0 ){
				for( core::Size iatm = 1; iatm <= grp.atms.size(); ++iatm ){
					core::Size atmno = grp.atms[iatm];
					if( atmno <= rsdtype.nheavyatoms() ) grp.comatms.push_back( atmno );
				}
			}

			TR.Debug << " / comatms: " << n_donor << " " << n_acceptor << ": ";
			for( core::Size i = 1; i <= grp.comatms.size(); ++i ) TR.Debug << " " << grp.comatms[i];
			TR.Debug << std::endl;

			grp.n_donor = n_donor;
			grp.n_acceptor =  n_acceptor;
			rsdgrps_[resname].push_back( grp );
		}

	}
}

Vector
GroupElec::get_grpdis2( 
               conformation::Residue const & rsd1,
							 conformation::Residue const & rsd2,
							 utility::vector1< Size > const &com1atms,
							 utility::vector1< Size > const &com2atms,
							 core::Vector &com1,
							 core::Vector &com2
			         ) const
{
  com1 = core::Vector( 0.0 );
  com2 = core::Vector( 0.0 );

  //core::Size n1( 0 ), n2( 0 );
  for( core::Size i = 1; i <= com1atms.size(); ++i ) com1 += rsd1.xyz( com1atms[i] ); 
  for( core::Size i = 1; i <= com2atms.size(); ++i ) com2 += rsd2.xyz( com2atms[i] ); 

  com1 /= (core::Real)( com1atms.size() );
  com2 /= (core::Real)( com2atms.size() );

  Vector dcom = com1 - com2;
  return dcom;
}

Real
GroupElec::eval_respair_group_coulomb( 
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
	//hbonds::HBondSet const & hbond_set
) const
{
  using namespace etable::count_pair;

	Real score( 0.0 );
	Real d2;
	ResElecGroup const &resgrp1 = get_group( rsd1.type() );
	ResElecGroup const &resgrp2 = get_group( rsd2.type() );

	// default "subtract"
	bool use_subtract( true ), use_shift( false );

	if( fade_type().compare( "shift" ) == 0 ){
		use_subtract = false; use_shift = true;
	} else if( fade_type().compare( "grpsubtract" ) == 0 ){
		use_subtract = true; use_shift = true;
	}

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	bool const is_bonded = ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) );
	bool intrares( rsd1.seqpos() == rsd2.seqpos() );

	/*
	TR << "entering score: " << rsd1.seqpos() << " " << rsd2.seqpos() << " " 
		 << rsd1.name() << " " << rsd2.name() << " " 
		 << resgrp1.size() << " " << resgrp2.size() << std::endl;
	*/

	for ( Size ii = 1; ii <= resgrp1.size(); ++ii ) {
		ElecGroup const &grp1 = resgrp1[ii];
		for ( Size jj = 1; jj <= resgrp2.size(); ++jj ) {
			if( intrares && ii >= jj ) continue;

			ElecGroup const &grp2 = resgrp2[jj];

			core::Vector com1, com2, dcom;
			dcom = get_grpdis2( rsd1, rsd2, grp1.comatms, grp2.comatms,
													com1, com2 );

			core::Real const grpdis2 = dcom.length_squared();

			//if( !use_subtract && (grpdis2 > coulomb().max_dis2()) ) continue;
			if( use_shift && (grpdis2 > coulomb().max_dis2()) ) continue;
			//if( grpdis2 > coulomb().max_dis2() ) continue;

			Real grp_cpweight( 1.0 );
			Size path_dist( 0 );
			if( is_bonded && grp_cpfxn_ )
				grp_cpweight = get_grp_countpair( grp1.atms, grp2.atms, cpfxn, 
																					path_dist );

			core::Real dE_dr( 0.0 );
			core::Real group_score( 0.0 );

			Real dsw_dr( 1.0 ), sw( 1.0 );
			//if( !use_subtract )
			if( use_shift )
				sw = eval_grp_trunc( false, grpdis2, false, dsw_dr );

			Size const n1( grp1.atms.size() );
			Size const n2( grp2.atms.size() );
			for ( Size kk = 1; kk <= n1; ++kk ){
				core::Size const &atm1( grp1.atms[kk] );
				Real const &q1( rsd1.atomic_charge( atm1 ) );

				for ( Size ll = 1; ll <= n2; ++ll ){
					core::Size const &atm2( grp2.atms[ll] );
					Real const &q2( rsd2.atomic_charge( atm2 ) );
					d2 = rsd1.xyz(atm1).distance_squared( rsd2.xyz(atm2) );

					if ( use_subtract && d2 > coulomb().max_dis2() ) continue;

					bool is_count( true );
					Real atom_cpweight( 1.0 );
					if( is_bonded && !grp_cpfxn_ ){
						path_dist = 0;
						is_count = cpfxn->count( atm1, atm2, atom_cpweight, path_dist );
					}

					if( !is_count ) continue;
					//atom_cpweight = 1.0;

					Real atompair_score( 0.0 );
					if( use_subtract ){
						atompair_score = coulomb().eval_atom_atom_fa_elecE( rsd1.xyz( atm1 ), q1, 
																																rsd2.xyz( atm2 ), q2 );
					} else {
						atompair_score = eval_standard_coulomb( q1, q2, d2, false, dE_dr );
					}

					/*
					TR << "kk/ll "  << ii << " " << jj << " " << kk << " " << ll << " " << std::sqrt(d2) 
						 << " " << atompair_score << " " << path_dist << " " << atom_cpweight
						 <<std::endl;
					*/

					/*
					TR << "Residue " << rsd1.seqpos() << " atom " << rsd1.atom_name(atm1) << " to Residue " << rsd2.seqpos() << " atom " << rsd2.atom_name(atm2)
						 << " q1 " << rsd1.atomic_charge(atm1) << " q2 " << rsd2.atomic_charge(atm2)
						 << " dist " << std::sqrt(d2) << " energy " << atompair_score << " cpwt " << atom_cpweight << std::endl;
					*/

					group_score += atom_cpweight*atompair_score;
				}
			}

			group_score *= grp_cpweight;

			// converge to constant for hbonding groups
			bool is_hbond_pair = (grp1.n_donor*grp2.n_acceptor + grp2.n_donor*grp1.n_acceptor > 0);
			//bool is_hbond_pair_from_HBscore = check_hbond_pairing_from_HBscore( grp1, grp2, hbond_set );
			//Real score_exp( group_score );

			Real dw_dE( 0.0 );
			//bool do_fade( false );
			if( is_hbond_pair && fade_hbond_ ) 
				fade_hbonding_group_score( grp1, grp2, group_score, dw_dE );

			/*
			if( std::abs(group_score) > 1.0e-5 && do_fade )
				printf("Score,Grp: %3d %3d %3d %3d %8.3f %2d %2d %8.5f %8.5f %8.5f\n", int(rsd1.seqpos()), int(rsd2.seqpos()),
							 int(ii), int(jj), std::sqrt(grpdis2), do_fade,
							 is_hbond_pair, grp_cpweight, score_exp, group_score);
			*/

			/*
			if( group_score_exp < -1.0 )
				TR << "do fade? " << rsd1.seqpos() << " " << rsd2.seqpos()
					 << " " << ii << " " << jj << ": " << do_fade << " "  << group_score_exp << " " << group_score << std::endl;
			*/

			group_score *= sw;

			score += group_score;
		} // grp2
	} // grp1

  return score;
}

bool
GroupElec::fade_hbonding_group_score( ElecGroup const &grp1,
																			ElecGroup const &grp2,
																			Real &group_score,
																			Real &dw_dE ) const
{
	if( grp1.qeps < 1e-3 || grp2.qeps < 1e-3 ) return false;

	Real Emin = -332.0637*grp1.qeps*grp2.qeps/4.0/coulomb().die();
	Real Efade = 0.9*Emin;
	dw_dE = 0.0;
	bool do_fade( false );

	//TR << "Fade?" << group_score << " " << Emin << " " << Efade << std::endl;
  if( group_score < Emin ){
		group_score = Emin;
		do_fade = true;

	} else if( group_score < Efade ){ // continuous
		// E' = 0.1Emin*(-x^3+x^2+x), x = 10.0*(E-0.9Emin)/Emin
		// E'(0) = 0, E'(1) = 0.1Emin
		// dE'dE = 0.1Emin*(-3x^2+2x+1)*dxdE, dxdE = 10/Emin => dE'dE = -3x^2+2x+1
		// dE'dE (0) = 0.1Emin*(10.0/Emin) = 1.0, dE'dx(1) = 0

		Real x = 10.0*(group_score-Efade)/Emin;
		dw_dE = -3.0*x*x + 2.0*x + 1.0;
		group_score = 0.9*Emin + 0.1*Emin*( -x*x*x + x*x + x );
		do_fade = false;
	}

	return do_fade;
}

void
GroupElec::eval_respair_group_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	//hbonds::HBondSet const & hbond_set,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs,
	core::Real const elec_weight,
	core::Real & Erespair
	) const 
{
  using namespace etable::count_pair;

	// default "subtract"
	bool use_subtract( true ), use_shift( false );
	if( fade_type().compare( "shift" ) == 0 ){
		use_subtract = false; use_shift = true;
	} else if( fade_type().compare( "grpsubtract" ) == 0 ){
		use_subtract = true; use_shift = true;
	}

	Erespair = 0.0;
	ResElecGroup const &resgrp1 = get_group( rsd1.type() );
	ResElecGroup const &resgrp2 = get_group( rsd2.type() );

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
	bool const is_bonded = ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) );
	bool intrares( rsd1.seqpos() == rsd2.seqpos() );

	// dummy array
	Vector Iv( 0.0 );
	utility::vector1< Vector > f2r1( r1_atom_derivs.size(), Iv );
	utility::vector1< Vector > f2r2( r2_atom_derivs.size(), Iv );

	for ( Size ii = 1; ii <= resgrp1.size(); ++ii ) {
		ElecGroup const &grp1 = resgrp1[ii];
		core::Size const &ncom1 = grp1.comatms.size();

		for ( Size jj = 1; jj <= resgrp2.size(); ++jj ) {
			if( intrares && (ii >= jj ) ) continue;

			ElecGroup const &grp2 = resgrp2[jj];
			core::Size const &ncom2 = grp2.comatms.size();

			core::Vector com1, com2, dcom;
			dcom = get_grpdis2( rsd1, rsd2, grp1.comatms, grp2.comatms,
													com1, com2 );
			core::Real const grpdis2 = dcom.length_squared();

			if ( use_shift && (grpdis2 > coulomb().max_dis2()) ) continue;
			//if ( !use_subtract && (grpdis2 > coulomb().max_dis2()) ) continue;
			//if ( grpdis2 > coulomb().max_dis2() ) continue;

			Real grp_cpweight( 1.0 );
			Size path_dist( 0 );

			if( is_bonded && grp_cpfxn_ )
				grp_cpweight = get_grp_countpair( grp1.atms, grp2.atms, cpfxn, 
																					path_dist );

			//core::Real grpdis = std::sqrt( grpdis2 );

			core::Size const &n1 = grp1.atms.size();
			core::Size const &n2 = grp2.atms.size();

			// get group-truncation info first
			Real dsw_dr( 1.0 ), sw( 1.0 );
			//if( !use_subtract )
			if( use_shift )
				sw = eval_grp_trunc( false, grpdis2, true, dsw_dr );

			// 1. derivative on Coulomb part: dE*sw
			core::Real group_score( 0.0 );

			utility::vector1< Vector > v1( n1, Iv );
			utility::vector1< Vector > v2( n2, Iv );

			for ( Size kk = 1; kk <= n1; ++kk ){
				core::Size const &atm1( grp1.atms[kk] );
				Vector const & atom1xyz( rsd1.xyz( atm1 ) );
				Real const &q1( rsd1.atomic_charge( atm1 ) );

				for ( Size ll = 1; ll <= n2; ++ll ){
					core::Size const &atm2( grp2.atms[ll] );
					Vector const & atom2xyz( rsd2.xyz( atm2 ) );
					Real const &q2( rsd2.atomic_charge( atm2 ) );

					Vector f2 = ( atom1xyz - atom2xyz );
					Real const &dis2( f2.length_squared() );
					Real dE_dr( 0.0 );

					if ( use_subtract && dis2 > coulomb().max_dis2() ) continue;

					bool is_count( true );
					Real atom_cpweight( 1.0 );
					if( is_bonded && !grp_cpfxn_ ){
						path_dist = 0;
						is_count = cpfxn->count( atm1, atm2, atom_cpweight, path_dist );
					}

					if( !is_count ) continue;
					atom_cpweight = 1.0;

					// untruncated energy
					Real atompair_score( 0.0 );
					if( use_subtract ){
						atompair_score = coulomb().eval_atom_atom_fa_elecE( rsd1.xyz( atm1 ), q1, 
																																rsd2.xyz( atm2 ), q2 );
						dE_dr = coulomb().eval_dfa_elecE_dr_over_r( dis2, q1, q2 );
					} else {
						atompair_score = eval_standard_coulomb( q1, q2, dis2, true, dE_dr );
					}

					Real sfxn_weight = atom_cpweight*elec_weight;

					group_score += atom_cpweight*atompair_score;

					f2 *= dE_dr*sw*sfxn_weight;
					v1[kk] += f2;
					v2[ll] -= f2;
				}
			}

			group_score *= grp_cpweight;
			Erespair += group_score;

			// converge to constant for hbonding groups
			bool is_hbond_pair = (grp1.n_donor*grp2.n_acceptor + grp2.n_donor*grp1.n_acceptor > 0);
			Real dw_dE( 1.0 );
			bool do_fade( false );

			if( is_hbond_pair && fade_hbond_ )
				do_fade = fade_hbonding_group_score( grp1, grp2, group_score, dw_dE );

			/*
			TR << "Deriv,groupscore: "  << rsd1.seqpos() << " " << rsd2.seqpos() 
				 << " " << ii << " " << jj << " " << std::sqrt(grpdis2) 
				 << " " << grp_cpfxn_ << " " << grp_cpweight << " " << path_dist
				 << " " << do_fade << " " << dw_dE << " " << sw << " " << dsw_dr
				 << " " << group_score << std::endl;
			*/

			if( dw_dE > 0.0 ){
				for ( Size kk = 1; kk <= n1; ++kk ){ v1[kk] *= dw_dE; }
 				for ( Size kk = 1; kk <= n2; ++kk ){ v2[kk] *= dw_dE; }
			}

			if( !do_fade ){ // long enough
				// E = E
				for ( Size kk = 1; kk <= n1; ++kk ){
					core::Size const &atm1( grp1.atms[kk] );
					f2r1[ atm1 ] += v1[kk]*grp_cpweight;
				}
				for ( Size ll = 1; ll <= n2; ++ll ){
					core::Size const &atm2( grp2.atms[ll] );
					f2r2[ atm2 ] += v2[ll]*grp_cpweight;
				}
			}
			// 2. derivative on truncation part: E*dsw
			// below will matter if applying different weight on bb/sc...

			//if( !use_subtract ){
			if( use_shift ){
				if( ncom1 > 0 ){
					core::Real const c_grp1_heavy = 1.0/((core::Real)(ncom1));

					for ( Size kk = 1; kk <= ncom1; ++kk ){
						core::Size const &atm1( grp1.comatms[kk] );
						Vector f2 = c_grp1_heavy*group_score*dsw_dr*dcom*elec_weight;
						f2r1[ atm1 ] += f2;
					}
				}

				if( ncom2 > 0 ){
					core::Real const c_grp2_heavy = 1.0/((core::Real)(ncom2));

					for ( Size kk = 1; kk <= ncom2; ++kk ){
						core::Size const &atm2( grp2.comatms[kk] );
						Vector f2 = -c_grp2_heavy*group_score*dsw_dr*dcom*elec_weight;
						f2r2[ atm2 ] += f2;
					}
				}
			}

		} // grp2
	} // grp2
  //TR << "end deriv" << std::endl;

	// finally get f1
	for( Size ii = 1; ii <= rsd1.natoms(); ++ii ){
		Vector const & atom1xyz( rsd1.xyz( ii ) );
		Vector const & f2 =	f2r1[ii];
		r1_atom_derivs[ ii ].f2() += f2;
		r1_atom_derivs[ ii ].f1() += atom1xyz.cross( -f2 );
	}
	for( Size ii = 1; ii <= rsd2.natoms(); ++ii ){
		Vector const & atom1xyz( rsd2.xyz( ii ) );
		Vector const & f2 =	f2r2[ii];
		r2_atom_derivs[ ii ].f2() += f2;
		r2_atom_derivs[ ii ].f1() += atom1xyz.cross( -f2 );
	}

	Erespair *= elec_weight;

}

Real
GroupElec::get_grp_countpair( 
							 utility::vector1< Size > const &grp1atms,
							 utility::vector1< Size > const &grp2atms,
							 etable::count_pair::CountPairFunctionCOP cpfxn,
							 Size &path_dist
							 ) const
{
	// hard code for now
	path_dist = 0;
	Size path_dist_min( 99 ), path_dist_max( 0 );
	//Real path_dist_avrg( 0.0 );
	for ( Size kk = 1; kk <= grp1atms.size(); ++kk ){
		core::Size const atm1( grp1atms[kk] );
		for ( Size ll = 1; ll <= grp2atms.size(); ++ll ){
			core::Size const atm2( grp2atms[ll] );

			Real weight;
			cpfxn->count( atm1, atm2, weight, path_dist );
			if( path_dist < path_dist_min ) path_dist_min = path_dist;
			if( path_dist > path_dist_max ) path_dist_max = path_dist;
		}
	}

	//if( grp_cpfxn_mode_.compare( "max" ) == 0 ){
	//path_dist = path_dist_max;
	//} else if( grp_cpfxn_mode_.compare( "average" ) == 0 ){
	//path_dist = Size(path_dist_avrg);
	//} else {
	path_dist = path_dist_min;
	//}

	if( path_dist > 5 ) return 1.0; // >1-6
	else if( path_dist < 3 ) return 0.0; // <1-4
	else return cpfxn_weight_[path_dist-2];
}

} // namespace elec
} // namespace scoring
} // namespace core

