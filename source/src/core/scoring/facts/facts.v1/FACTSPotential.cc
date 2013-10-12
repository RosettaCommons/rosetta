// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// @file:   core/scoring/facts/FACTSPotential.cc
// @brief:  The definitions of 3 classes of the FACTS algorithm resides here (see devel/khorvash/FACTSPotential.hh
// @author: Massih Khoravash
// @author: Hahnbeom Park

// Unit headers
#include <core/scoring/facts/FACTSPotential.hh>

// Project headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AA.hh>

//#include <core/chemical/MMAtomTypeSet.hh>
//#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>

#include <basic/prof.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <cassert>
#include <utility/assert.hh>

static basic::Tracer TR("core.scoring.FACTSPotential");

using namespace std;

# define Math_PI 3.14159265358979323846

namespace core {
namespace scoring {

/**************************************************************************************************/
/*                                                                                                */
/*    @brief: The FACTSRsdTypeInfo class provides all the constants and parameters for given      */
/*            residue type                                                                        */
/*                                                                                                */
/**************************************************************************************************/
void FACTSRsdTypeInfo::create_info( chemical::ResidueType const & rsd )
{
	initialize_parameters( rsd );
	initialize_selfpair( rsd );
}

// This function initializes native parameters that are used for empirical function calculations
void FACTSRsdTypeInfo::initialize_parameters( chemical::ResidueType const & rsd ){

  using namespace basic::options;
  using namespace basic::options::OptionKeys;
	
	// First define natoms
	natoms_ = rsd.natoms();

	//grab data from "database/chemical/atom_type_sets/fa_standard/extras/facts_params.txt
	Size const FACTS_RADIUS_INDEX( rsd.atom_type_set().extra_parameter_index( "FACTS_RADIUS" ) );
	Size const FACTS_CUT_INDEX   ( rsd.atom_type_set().extra_parameter_index( "FACTS_CUT"    ) );
	Size const FACTS_B1_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_B1"     ) );
	Size const FACTS_B2_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_B2"     ) );
	Size const FACTS_D1_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_D1"     ) );
	Size const FACTS_D2_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_D2"     ) );
	Size const FACTS_A0_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_A0"     ) );
	Size const FACTS_A1_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_A1"     ) );
	Size const FACTS_A2_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_A2"     ) );
	Size const FACTS_A3_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_A3"     ) );
	Size const FACTS_C0_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_C0"     ) );
	Size const FACTS_C1_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_C1"     ) );
	Size const FACTS_C2_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_C2"     ) );
	Size const FACTS_C3_INDEX    ( rsd.atom_type_set().extra_parameter_index( "FACTS_C3"     ) );

	Size FACTS_ALPHA_INDEX;
	Size const asp_patch( option[ score::facts_asp_patch ]() );
	if ( asp_patch == 2 ){
		FACTS_ALPHA_INDEX = rsd.atom_type_set().extra_parameter_index( "FACTS_ALPHA2"  );
	} else if ( asp_patch == 3 ){
		FACTS_ALPHA_INDEX = rsd.atom_type_set().extra_parameter_index( "FACTS_ALPHA3"  );
	} else {
		FACTS_ALPHA_INDEX = rsd.atom_type_set().extra_parameter_index( "FACTS_ALPHA"  );
	}

	// Initialize array sizes & put in parameters
	q_.resize( natoms() );
	COradius2_.resize( natoms() );
	volume_.resize( natoms() );
	alpha_.resize( natoms() );
	a0_.resize( natoms() );
	a1_.resize( natoms() );
	a2_.resize( natoms() );
	a3_.resize( natoms() );
	b1_.resize( natoms() );
	b2_.resize( natoms() );
	c0_.resize( natoms() );
	c1_.resize( natoms() );
	c2_.resize( natoms() );
	c3_.resize( natoms() );
	d1_.resize( natoms() );
	d2_.resize( natoms() );
	not_using_.resize( natoms() );
	is_chargedH_.resize( natoms() );

	for(Size i = 1; i <= natoms(); ++i){
		core::chemical::AtomType const &type = rsd.atom_type(i);

		// Partial charge
		q_[i] = rsd.atomic_charge( i );

		// Residue polarity
		string atmname( type.atom_type_name() );
		bool is_chargedH( atmname.compare( 0, 4, "Hpol" ) == 0 && 
											(	rsd.aa() == core::chemical::aa_arg ||
												rsd.aa() == core::chemical::aa_lys ||
												rsd.aa() == core::chemical::aa_his ) );

		is_chargedH_[i] = is_chargedH;
		Real vdw_radius;

		//Corrections for atomic parameters
		// 1. HIS aromatic carbons to be consistent with CHARMM definition
		if( rsd.aa() == core::chemical::aa_his && atmname == "aroC" ){
			vdw_radius = 1.8;
			alpha_[i] = 0.02;
			COradius2_[i] = 8.37057*8.37057;
			b1_[i] =  0.853399e+1; b2_[i] = -0.109161e+1;
			d1_[i] = -0.227929e+5; d2_[i] = -0.825630e+1;
			a0_[i] = -0.123827e+3; a1_[i] =  0.123827e+3; a2_[i] = 0.185613e-2; a3_[i] =  0.347902e+3;
			c0_[i] =  0.315402e+2; c1_[i] = -0.313502e+2; c2_[i] = 0.567384e-3; c3_[i] = -0.142830e+5;

  	// Otherwise just use default
		} else {
			vdw_radius = type.extra_parameter( FACTS_RADIUS_INDEX );
			if ( vdw_radius <= 1.0e-6 ){
				not_using_[i] = true;
			} else {
				not_using_[i] = false;
			}

			alpha_[i] = type.extra_parameter( FACTS_ALPHA_INDEX );
			COradius2_[i] = type.extra_parameter( FACTS_CUT_INDEX )*type.extra_parameter( FACTS_CUT_INDEX );
			if( COradius2_[i] <= 1.0e-3 ){
				TR << "Unrealistic cutoff radii: set to new cut " << 100.0 << std::endl;
				COradius2_[i] = 10.0*10.0;
			}

			b1_[i] = type.extra_parameter( FACTS_B1_INDEX );
			b2_[i] = type.extra_parameter( FACTS_B2_INDEX );
			d1_[i] = type.extra_parameter( FACTS_D1_INDEX );
			d2_[i] = type.extra_parameter( FACTS_D2_INDEX );
			a0_[i] = type.extra_parameter( FACTS_A0_INDEX );
			a1_[i] = type.extra_parameter( FACTS_A1_INDEX );
			a2_[i] = type.extra_parameter( FACTS_A2_INDEX );
			a3_[i] = type.extra_parameter( FACTS_A3_INDEX );
			c0_[i] = type.extra_parameter( FACTS_C0_INDEX );
			c1_[i] = type.extra_parameter( FACTS_C1_INDEX );
			c2_[i] = type.extra_parameter( FACTS_C2_INDEX );
			c3_[i] = type.extra_parameter( FACTS_C3_INDEX );
		}
		volume_[i] = (4.0/3.0) * Math_PI * vdw_radius * vdw_radius * vdw_radius;
	}

	// Reduce Born radii for charged polarH: sidechain Hpol of ARG/LYS/HIS
	/*
	for(Size i = 1; i <= natoms(); ++i){
		if( option[ score::facts_saltbridge_correction ].user() &&is_chargedH ){
			a0_[i] *= option[ score::facts_saltbridge_correction ]();
			a1_[i] *= option[ score::facts_saltbridge_correction ]();
			TR.Debug << rsd.seqpos() << " " << atmname << " " << a0_[i] << std::endl;
		}
	}
	*/

}

void FACTSRsdTypeInfo::initialize_selfpair( chemical::ResidueType const & rsd ){

	selfpair_.resize( rsd.natoms() );

	bool const plane_to_self( basic::options::option[ basic::options::OptionKeys::score::facts_plane_to_self ]() );

	for( Size atm1 = 1; atm1 <= rsd.natoms(); ++atm1 ){
		selfpair_[atm1].resize( rsd.natoms(), false );

		for( Size atm2 = 1; atm2 <= rsd.natoms(); ++atm2 ){
			if( rsd.path_distance(atm1,atm2) <= 2 ){
				selfpair_[atm1][atm2] = true;

			} else if( plane_to_self && rsd.aa() == core::chemical::aa_arg ) {
				// Add xHHx - NHx / xHHx - xHHx / xHHx - NE / xHHx - HE / NHx - HE
				// or in other words, add all pairs in H or E position
				if( rsd.atom_name(atm1).compare( 2, 1, "H" ) == 0 && 
						(	rsd.atom_name(atm2).compare( 2, 1, "H" ) == 0 || rsd.atom_name(atm2).compare( 2, 1, "E" ) == 0) )
					selfpair_[atm1][atm2] = true;

			} else if( plane_to_self && rsd.aa() == core::chemical::aa_gln ) {
				// Add all pairs in E position
				if( rsd.atom_name(atm1).compare( 2, 1, "E" ) == 0 && rsd.atom_name(atm2).compare( 2, 1, "E" ) == 0 )
					selfpair_[atm1][atm2] = true;

 			} else if( plane_to_self && rsd.aa() == core::chemical::aa_asn ) {
				// Add all pairs in D position
				if( rsd.atom_name(atm1).compare( 2, 1, "D" ) == 0 && rsd.atom_name(atm2).compare( 2, 1, "D" ) == 0 )
					selfpair_[atm1][atm2] = true;

 			} else if( plane_to_self && rsd.aa() == core::chemical::aa_his ) {
				// Add all pairs in D or E position
				if( rsd.atom_name(atm1).compare( 2, 1, "D" ) == 0 && 
						( rsd.atom_name(atm2).compare( 2, 1, "D" ) == 0 || rsd.atom_name(atm2).compare( 2, 1, "E" ) == 0 ) )
					selfpair_[atm1][atm2] = true;

			}
		}
	}
} // END void initialize_selfpair
/// FACTSRsdTypeInfo

/**************************************************************************************************/
/*                                                                                                */
/*    @brief: The FACTSResidueInfo class provides all the functions, constants and parameters     */
/*            for different atoms, which are required to calculate the solvation free energy of   */
/*                      of a molecule embedded in water using FACTS method                        */
/*                                                                                                */
/**************************************************************************************************/

//This function initializes all the values for FACTS original parameters, atomic volume, Ai, Bi, esolvE, sasa...
void FACTSResidueInfo::initialize( conformation::Residue const & rsd,  
																	 bool const is_rotamer )
{
	natoms_ = rsd.natoms();

	// This variable is used in res_res_burial and evaluate_polar(nonpolar)_energy
	flag_for_calculation_ = utility::vector1< bool >( natoms(), false );

	// Initialize Arrays
	Vector i( 0.0 );
	nmtr_ = utility::vector1< Vector >( natoms(), i ); // dnmtr of Ai (equation 4 on page 704 of FACTS paper)
	dnmtr_ = utility::vector1< Real >( natoms(), 1.0 ); // dnmtr of Bi (equation 4 on page 704 of FACTS paper)	
	Ai_ = utility::vector1< Real >( natoms(), 0.0 );// Ai (equation 3 on page 704 of FACTS paper)						
	Bi_ = utility::vector1< Real >( natoms(), 0.0 );// Bi (equation 4 on page 704 of FACTS paper)						
	Ci_ = utility::vector1< Real >( natoms(), 0.0 );// Ci (equation 6 on page 704 of FACTS paper)						
	Di_ = utility::vector1< Real >( natoms(), 0.0 );// Di (equation 10 on page 704 of FACTS paper)						

	esolvE_ = utility::vector1< Real >( natoms(), 0.0 ); // DeltaGi (equation 7 on page 704 of FACTS paper)			
	sasa_ = utility::vector1< Real >( natoms(), 0.0 ); // atomic SASA (equation 11 on page 706 of FACTS paper)	
	BR_ = utility::vector1< Real >( natoms(), 0.0 ); // BornRadius

	// auxiliary arrays for calculating derivatives
	if( !is_rotamer ){
		dG_dCi_ = utility::vector1< Real >( natoms() ); 
		dSA_dDi_ = utility::vector1< Real >( natoms() );
		dE_dBR_ = utility::vector1< Real >( natoms() );
		dE_drij2_ = utility::vector1< Vector >( natoms(), i );
		polarF2d_ = utility::vector1< Vector >( natoms(), i );
		polarF2BR_ = utility::vector1< Vector >( natoms(), i );
		nonpolarF2_ = utility::vector1< Vector >( natoms(), i );
	}

} // END void FACTSResidueInfo::initialize( conformation::Residue const & rsd){

void FACTSResidueInfo::allocate_selfneigh(){
	selfneigh_ = utility::vector1< SelfNeighInfo >( natoms() ); // Array for containing Self-neighbors
	for( Size i_atm = 1; i_atm <= natoms(); ++i_atm ){
		selfneigh_[i_atm].nneigh = 0;
		selfneigh_[i_atm].resID.resize( MAXNEIGH );
		selfneigh_[i_atm].atmID.resize( MAXNEIGH );
	}
}

void FACTSResidueInfo::store_xyz( Residue const &rsd ){
	xyz_.resize( 0 );
	for( Size iatm = 1; iatm <= rsd.natoms(); ++iatm ){
		xyz_.push_back( rsd.xyz(iatm) );
	}
}

/**************************************************************************************************/
/*                                                                                                */
/*    @breif: The class FACTSPoseInfo                                                             */
/*                                                                                                */
/**************************************************************************************************/

// Constructor
FACTSPoseInfo::FACTSPoseInfo( FACTSPoseInfo const & src ):CacheableData() 
{
	Size const src_size( src.size() );

	residue_info_.resize( src_size );
	placeholder_residue_.resize( src_size );
	placeholder_info_.resize( src_size );

	for ( Size i=1; i<= src_size; ++i ) {
		residue_info_[i] = src.residue_info_[i]->clone();
		if ( src.placeholder_residue_[i] ) {
			placeholder_residue_[i] = src.placeholder_residue_[i]->clone();
			placeholder_info_[i] = src.placeholder_info_[i]->clone();
		} else {
			placeholder_residue_[i] = 0;
			placeholder_info_[i] = 0;
		}
	}
	being_packed_ = src.being_packed_;
}

void FACTSPoseInfo::initialize( pose::Pose const & pose )
{
	Size const nres( pose.total_residue() );

	residue_info_.resize( nres, 0 );
	placeholder_residue_.resize( nres, 0 );
	placeholder_info_.resize( nres, 0 );

	for ( Size i=1; i<= nres; ++i ) {
		if ( !residue_info_[i] ){
			residue_info_[i] = new FACTSResidueInfo( pose.residue(i) );
		}else	{
			residue_info_[i]->initialize( pose.residue(i) );
		}
	}
	being_packed_.clear();
}

void FACTSPoseInfo::set_placeholder( Size const i, ResidueOP rsd, FACTSResidueInfoOP info )
{
	placeholder_residue_[ i ] = rsd;
	placeholder_info_[ i ] = info;
}

void FACTSPoseInfo::set_repack_list( utility::vector1< bool > const & repacking_residues )
{
	being_packed_.resize( size(), false );
	for ( Size i=1; i<= size(); ++i ) {
		being_packed_[i] = repacking_residues[ i ];
	}
}

bool FACTSPoseInfo::is_changed( pose::Pose const &pose ){
	for( Size ires = 1; ires <= pose.total_residue(); ++ ires ){

		Size const natom( pose.residue(ires).natoms() );
		utility::vector1<Vector> const facts_xyz = residue_info( ires ).xyz();

		if( natom != facts_xyz.size() )
			return true;

		for( Size iatm = 1; iatm <= natom; ++ iatm ){
			Vector const dxyz = facts_xyz[iatm] - pose.residue(ires).xyz(iatm);
			Real const d2 = dxyz.dot(dxyz);
			if( d2 > 1.0e-6 ) return true;
		}

	}
	return false;
}

/**************************************************************************************************/
/*                                                                                                */
/*    @brief: The  class    FACTSRotamerSetInfo                                                   */
/*                                                                                                */
/**************************************************************************************************/

void FACTSRotamerSetInfo::initialize( RotamerSet const & rotamer_set )
{
	Size const nrot( rotamer_set.num_rotamers() );
	residue_info_.resize( nrot );
	for ( Size i=1; i<= nrot; ++i ) {
		residue_info_[i] = new FACTSResidueInfo( *rotamer_set.rotamer(i), true );
	}
}

/**************************************************************************************************/
/*                                                                                                */
/*    @breif: The FACTSPotential class provides all the functions, constants, and parameters      */
/*            common to all atoms required to calculate the free energy of solvation of a         */
/*               (macro)molecule embedded in a continuum solvent using FACTS method               */
/*                                                                                                */
/**************************************************************************************************/
FACTSPotential::FACTSPotential ():
	MultiplicitiveFactor_(332.07156)
{
	set_default();
}

void
FACTSPotential::set_default(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Kappa_ = option[ score::facts_kappa ]();
	GBpair_cut_ = option[ score::facts_GBpair_cut ]();
 	//Tau_ = (1.0/option[ score::hackelec_die ]() ) - (1.0/78.5);
 	Tau_ = (1.0/option[ score::hackelec_die ]() ) - (1.0/20.0);
	do_apprx = option[ score::facts_apprx ]();
	min_dis_ = option[ score::hackelec_min_dis ]();
	selfenergy_scale_ = option[ score::facts_selfenergy_scale ]();
	intrares_scale_ = option[ score::facts_intrares_scale ]();
	saltbridge_correction_ = option[ score::facts_saltbridge_correction ]();
	dshift2_ = option[ score::facts_dshift ]()*option[ score::facts_dshift ]();

	if( do_apprx ){
		TR << "FACTS Approximation turned on." << endl;
	}
	TR << "FACTS ASP set using: " << option[ score::facts_asp_patch ]()<< endl;

	runtime_assert( option[score::hackelec_die]() == 1.0 );
	runtime_assert( option[score::hackelec_r_option]() );
}

// Dummy thing... I don't know why this is required
void FACTSPotential::build_placeholders(
																				pose::Pose const & pose,
																				FACTSPoseInfo & facts_info
																				) const
{
	Size const nres( pose.total_residue() );

	chemical::ResidueTypeSet const & residue_set( pose.residue(1).residue_type_set() );

	for ( Size i=1; i<= nres; ++i ) {
		if ( facts_info.being_packed(i) ) {
			Residue const & existing_rsd( pose.residue(i) );
			// build a placeholder at this position
			if ( existing_rsd.is_protein() ) {
				chemical::ResidueTypeCAP protein_placeholder_residue_type( &( residue_set.name_map("GB_AA_PLACEHOLDER") ) );
				// use appropriate termini variants if necessary:
				if ( existing_rsd.is_lower_terminus() ) {
					protein_placeholder_residue_type =
						&(residue_set.get_residue_type_with_variant_added( *protein_placeholder_residue_type,
																															 chemical::LOWER_TERMINUS ) );
				}
				if ( existing_rsd.is_upper_terminus() ) {
					protein_placeholder_residue_type =
						&(residue_set.get_residue_type_with_variant_added( *protein_placeholder_residue_type,
																															 chemical::UPPER_TERMINUS ) );
				}

				conformation::ResidueOP rsd
					( conformation::ResidueFactory::create_residue( *protein_placeholder_residue_type, existing_rsd,
																													pose.conformation() ) );
				FACTSResidueInfoOP rsd_info( new FACTSResidueInfo( *rsd ) );

				Size const dummy_index( rsd->atom_index("DUMM") );
				//rsd_info->set_atomic_radius( dummy_index,  dummy_radius_);
				assert( std::fabs( rsd->xyz("CA").distance( rsd->xyz( dummy_index )) - 2.44 ) < 1e-2 );
				facts_info.set_placeholder( i, rsd, rsd_info );
			}
		}
	}
}

//This function "pre-calculates" all the atomic contents (Born radius, sasa,...) and energy values
void FACTSPotential::setup_for_scoring(pose::Pose & pose, bool const & packing) const
{
	Size res1;
	Size res2;

	//timeval t1, t2, t3, t4;
	//gettimeofday(&t1, NULL );

	PROF_START( basic::FACTS_GET_ALL_BORN_RADII );
	Size const nres( pose.total_residue() );
	FACTSPoseInfoOP facts_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ) ) {
		facts_info = static_cast< FACTSPoseInfo* > 
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )() );
	} else {
		facts_info = new FACTSPoseInfo();
	}

	facts_info->initialize( pose );

	// Call this everytime in setup_for_scoring
	create_rsdtypeinfo( pose );

	if( !packing ){
		for ( res1 = 1; res1 <= nres; ++res1 ) {
			FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );
			facts1.allocate_selfneigh();
			facts1.store_xyz( pose.residue(res1) );
		}
	}

	Energies const & energies( pose.energies() );
	EnergyGraph const & energy_graph( energies.energy_graph() );

	// 1. First get Born radius, Solvation energy, SASA for all atoms
	for ( res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd1( pose.residue( res1 ) );
		Size natoms1 = rsd1.natoms();
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );

		if( facts1.natoms() == 0) continue;

		res_res_burial_for_scoring( rsd1, facts1, rsd1, facts1 );

		for ( graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node( res1 )->const_edge_list_begin(),
						irue = energy_graph.get_node( res1 )->const_edge_list_end();
					iru != irue; ++iru ) {
			Size const res2( (*iru)->get_other_ind( res1 ) );
			
			// Iter over upper only 
			if (res2 < res1) continue;
			
			conformation::Residue const & rsd2( pose.residue( res2 ) );
			res_res_burial_for_scoring( rsd1, facts1, pose.residue( res2 ), facts_info->residue_info( res2 ) );
		}
	}

	//gettimeofday(&t2, NULL );
	//Real elapsedTime1 = (t2.tv_sec - t1.tv_sec) * 1000.0;
	//elapsedTime1 += (t2.tv_usec - t1.tv_usec) / 1000.0;

	// 2. Refresh Born radii / SASA
	for ( res1 = 1; res1 <= nres; ++res1 ){
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );

		FACTSRsdTypeInfoCOP factstype1 = get_factstype( pose.residue( res1 ).type() );
		get_self_terms( factstype1, facts1, packing );
	}

	//gettimeofday(&t3, NULL );
	//Real elapsedTime2 = (t3.tv_sec - t2.tv_sec) * 1000.0;
	//elapsedTime2 += (t3.tv_usec - t2.tv_usec) / 1000.0;

	// 3. Then Pre-calculate all GB-pair-related parts and store them in FACTSINFO
	for ( res1 = 1; res1 <= nres; ++ res1 ){
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );
		Residue const & rsd1( pose.residue( res1 ) );
		facts1.GBpair_ = utility::vector1< Real >( nres, 0.0 );

		calculate_GBpair_exact( rsd1, rsd1, facts1, facts1 );

		for ( graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node( res1 )->const_edge_list_begin(),
						irue = energy_graph.get_node( res1 )->const_edge_list_end();
					iru != irue; ++iru ) {

			Size const res2( (*iru)->get_other_ind( res1 ) );
			// Iter over upper only 
			if (res2 < res1) continue;

			FACTSResidueInfo & facts2( facts_info->residue_info( res2 ) );
			Residue const & rsd2( pose.residue( res2 ) );
			calculate_GBpair_exact( rsd1, rsd2, facts1, facts2 );

			//if ( do_apprx ){
			//	calculate_GBpair_apprx( rsd1, rsd2, facts1, facts2 );
			//} else {
			//  calculate_GBpair_exact( rsd1, rsd2, facts1, facts2, packing );
			//}
		}
	}

	// 3. Finally store everything into pose
	pose.data().set( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO, facts_info );
	PROF_STOP( basic::FACTS_GET_ALL_BORN_RADII );

	/*
	// Log for Born radii and SASA
	for( res1 = 1; res1 <= nres; ++ res1 ){
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );
		Residue const & rsd1( pose.residue( res1 ) );

		for( Size atm1 = 1; atm1 <= pose.residue(res1).natoms(); ++atm1 ){

			Real norm1 = sqrt(facts1.polarF2d(atm1)[0]*facts1.polarF2d(atm1)[0] +
												facts1.polarF2d(atm1)[1]+facts1.polarF2d(atm1)[1]+
												facts1.polarF2d(atm1)[2]+facts1.polarF2d(atm1)[2]);
			Real norm2 = sqrt(facts1.polarF2BR(atm1)[0]*facts1.polarF2BR(atm1)[0] +
												facts1.polarF2BR(atm1)[1]+facts1.polarF2BR(atm1)[1]+
												facts1.polarF2BR(atm1)[2]+facts1.polarF2BR(atm1)[2]);

			string atmname = pose.residue(res1).atom_type(atm1).atom_type_name();

			if( atmname.compare("HNbb") == 0 || atmname.compare("Hpol") == 0 || 
					atmname.compare("OOC") == 0 || atmname.compare("OCbb") == 0){
				cout << "Res/Atm/Type/" << setw(4) << res1 << setw(4) << atm1 << " " << setw(4) << atmname;
				cout << "| " << setw(10) << facts1.BR(atm1);
				cout << "| " << setw(10) << facts1.sasa(atm1);
				cout << "| " << setw(10) << norm1 << " " << setw(10) << norm2 << std::endl;
				//cout << "|" << setw(10) << facts1.b1(atm1) << "|" << setw(10) << facts1.b2(atm1);
				//cout << "|" << setw(10) << facts1.d1(atm1) << "|" << setw(10) << facts1.d2(atm1);
				//cout << "|" << setw(10) << facts1.Ai(atm1) << "|" << setw(10) << facts1.Bi(atm1);
				//cout << "|" << setw(10) << facts1.Ci(atm1) << "|" << setw(10) << facts1.Di(atm1) << endl;
			}
		}
	}
	*/


	//gettimeofday(&t4, NULL );
	//Real elapsedTime3 = (t4.tv_sec - t3.tv_sec) * 1000.0;
	//elapsedTime3 += (t4.tv_usec - t3.tv_usec) / 1000.0;
	//cout << "Setup for scoring 1/2/3: " << elapsedTime1 << " " << elapsedTime2 << " " << elapsedTime3 << " ms." << endl;

}// END setup_for_scoring

void FACTSPotential::create_rsdtypeinfo( pose::Pose & pose ) const
{
	for( Size i_res = 1; i_res <= pose.total_residue(); ++i_res ){
		core::chemical::ResidueType const &rsdtype = pose.residue(i_res).type();

		std::map< chemical::ResidueType const *, FACTSRsdTypeInfoCOP >::const_iterator iter = 
			FACTSrsdtypemap_.find( &rsdtype );

		// Add new if not found
		if ( iter == FACTSrsdtypemap_.end() ) {
			TR << "Adding new FACTS residue type info: " << rsdtype.name() << std::endl;

			FACTSRsdTypeInfoOP rsdtypeinfo = new FACTSRsdTypeInfo;
			rsdtypeinfo->create_info( rsdtype );

			FACTSrsdtypemap_[ &rsdtype ] = rsdtypeinfo;

			// Assert if it is added
			std::map< chemical::ResidueType const *, FACTSRsdTypeInfoCOP >::const_iterator iter2 = 
				FACTSrsdtypemap_.find( & rsdtype );
			assert ( iter2 != FACTSrsdtypemap_.end() );
		}
	}
}

FACTSRsdTypeInfoCOP 
FACTSPotential::get_factstype( core::chemical::ResidueType const &rsdtype ) const
{

	std::map< chemical::ResidueType const *, FACTSRsdTypeInfoCOP >::const_iterator iter = FACTSrsdtypemap_.find( &rsdtype );
	if( iter == FACTSrsdtypemap_.end() ){
		TR << "Adding new FACTS residue type info: " << rsdtype.name() << std::endl;
		
		FACTSRsdTypeInfoOP rsdtypeinfo = new FACTSRsdTypeInfo;
		rsdtypeinfo->create_info( rsdtype );
		FACTSrsdtypemap_[ &rsdtype ] = rsdtypeinfo;
		
		// Assert if it is added
		std::map< chemical::ResidueType const *, FACTSRsdTypeInfoCOP >::const_iterator iter2 = 
			FACTSrsdtypemap_.find( & rsdtype );
		assert ( iter2 != FACTSrsdtypemap_.end() );
		return iter2->second;

	} else {
		return iter->second;
	}
}

// This function evaluates Ai (volume of each atom) and components of Bi (symmetry of each atom),
// which are converted to Born Radius & SASA in the next step
void FACTSPotential::res_res_burial(
																		conformation::Residue const & rsd1,
																		FACTSResidueInfo & facts1,
																		conformation::Residue const & rsd2,
																		FACTSResidueInfo const & facts2
																		) const
{
	bool const same_res( rsd1.seqpos() == rsd2.seqpos() );
	Size natoms1 = rsd1.natoms();
	Size natoms2 = rsd2.natoms();
	Vector xyz2, xyz1, dyxz;

	FACTSRsdTypeInfoCOP factstype1 = get_factstype( rsd1.type() );
	FACTSRsdTypeInfoCOP factstype2 = get_factstype( rsd2.type() );

	// START: This part calculates the volume and symmetry of the solute around atom atm
	// and maintains it in Ai_ and nmtr_ & dnmtr_
	for ( Size atm1 = 1; atm1 <= natoms1; ++atm1 ) {
		xyz1 = rsd1.xyz(atm1);

		if ( factstype1->not_using(atm1) ) continue;

		for ( Size atm2 = 1; atm2 <= natoms2; ++atm2 ) {
			if (same_res && (atm1 == atm2 )) continue;

			xyz2 = (rsd2.xyz( atm2 ));
			Vector const dxyz( xyz1 - xyz2 );

			if ( factstype2->not_using(atm2) || factstype2->volume(atm2) < 1e-3) continue;

			Real dis2 = xyz1.distance_squared( xyz2 );

			// take special care for distance to prevent from exploding
			if( dis2 != dis2 ) continue;
			if( dis2 < 0.01 ) dis2 = 0.01;

			Real dis = std::sqrt(dis2);

			Real CutOff_sqr = factstype1->COradius2(atm1);

			//this is a redundant check as there is a stricter check on the next if statment ...
			if( dis2 >= CutOff_sqr ) continue;
			facts1.flag_for_calculation_[atm1] = true;

			// Equation 5 on page 704 of FACTS paper
			Real theta_sqrt = 1.0 - (dis2 / CutOff_sqr);
			Real thetaij = theta_sqrt*theta_sqrt;

			// The term within the sigma of equation 3 on page 704 of FACTS'
			Real Vi = factstype2->volume(atm2) * thetaij;

			// 1. Ai in equation 3 of page 704 of FACTS paper
			facts1.Ai_[atm1] += Vi; 

			// 2. Bi: the xyz coordinate of the nmtr of equation 4 on page 704 of FACTS paper
			facts1.nmtr_[atm1]  += (Vi/dis2)*dxyz;
			facts1.dnmtr_[atm1] +=  Vi/dis;
		}
	}
}//end res_res_burial

// This function has same logic with res_res_burial, but modified for efficient scoring 
// and for derivative evaluation
void FACTSPotential::res_res_burial_for_scoring(
																		conformation::Residue const & rsd1,
																		FACTSResidueInfo & facts1,
																		conformation::Residue const & rsd2,
																		FACTSResidueInfo & facts2
																		) const
{
	bool const same_res( rsd1.seqpos() == rsd2.seqpos() );
	//Size natoms1 = rsd1.natoms();
	//Size natoms2 = rsd2.natoms();
	Real thetaij, theta_sqrt;
	Real MAX_SELFDCUT2 = min(100.0, GBPair_cut()*GBPair_cut() );

	FACTSRsdTypeInfoCOP factstype1 = get_factstype( rsd1.type() );
	FACTSRsdTypeInfoCOP factstype2 = get_factstype( rsd2.type() );
	Size natoms1 = factstype1->natoms();
	Size natoms2 = factstype2->natoms();

	runtime_assert( natoms1 == rsd1.natoms() );
	runtime_assert( natoms2 == rsd2.natoms() );

	// START: This part calculates the volume and symmetry of the solute around atom atm
	// and maintains it in Ai_ and nmtr_ & dnmtr_
	for ( Size atm1 = 1; atm1 <= natoms1; ++atm1 ) {
		Vector const &xyz1 = rsd1.xyz(atm1);
		if ( factstype1->not_using(atm1) ) continue;

		Real CutOff_sqr1 = min( factstype1->COradius2(atm1), MAX_SELFDCUT2 );
		SelfNeighInfo &neigh1 = facts1.selfneigh_[atm1];
		Real const &V1 = factstype1->volume(atm1);

		// Iterate only for upper diagonal :)
		for ( Size atm2 = 1; atm2 <= natoms2; ++atm2 ) {
			if (same_res && (atm1 >= atm2)) continue;
			if ( factstype2->not_using(atm2) ) continue;

			Vector const &xyz2 = rsd2.xyz( atm2 );
			Vector const dxyz( xyz1 - xyz2 );
			Real dis2 = dxyz[0]*dxyz[0] + dxyz[1]*dxyz[1] + dxyz[2]*dxyz[2];

			// take special care for distance to prevent from exploding
			if( dis2 != dis2 ) continue;
			if( dis2 < 0.01 ) dis2 = 0.01;

			if ( dis2 > MAX_SELFDCUT2 ) continue;

			Real dis = std::sqrt(dis2);
			Real CutOff_sqr2 = min( factstype2->COradius2(atm2),	MAX_SELFDCUT2 );

			// Consideration for the first atom
			if( dis2 <= CutOff_sqr1 && factstype2->volume(atm2) > 1e-3 ){
				if( neigh1.nneigh >= facts1.MAXNEIGH ) continue;
				facts1.flag_for_calculation_[atm1] = true;
				neigh1.nneigh++;
				neigh1.resID[ neigh1.nneigh ] = rsd2.seqpos();
				neigh1.atmID[ neigh1.nneigh ] = atm2;

				// Equation 5 on page 704 of FACTS paper
				theta_sqrt = 1.0 - (dis2 / CutOff_sqr1);
				thetaij = theta_sqrt*theta_sqrt;
				// The term within the sigma of equation 3 on page 704 of FACTS'
				Real Vi = factstype2->volume(atm2) * thetaij;

				// 1. Ai in equation 3 of page 704 of FACTS paper
				facts1.Ai_[atm1] += Vi; 

				// 2. Bi: the xyz coordinate of the nmtr of equation 4 on page 704 of FACTS paper
				facts1.nmtr_[atm1]  += (Vi/dis2)*dxyz;
				facts1.dnmtr_[atm1] +=  Vi/dis;
			}

			// Consideration for the second atom
			if( dis2 <= CutOff_sqr2 && factstype1->volume(atm1) > 1e-3 ){

				SelfNeighInfo &neigh2 = facts2.selfneigh_[atm2];

				if (neigh2.nneigh >= facts2.MAXNEIGH) continue;

				facts2.flag_for_calculation_[atm2] = true;
				neigh2.nneigh++;
				neigh2.resID[ neigh2.nneigh ] = rsd1.seqpos();
				neigh2.atmID[ neigh2.nneigh ] = atm1;

				// Equation 5 on page 704 of FACTS paper
				theta_sqrt = 1.0 - (dis2 / CutOff_sqr2);
				thetaij = theta_sqrt*theta_sqrt;

				// The term within the sigma of equation 3 on page 704 of FACTS'
				Real const Vi = V1*thetaij;

				// 1. Ai in equation 3 of page 704 of FACTS paper
				facts2.Ai_[atm2] += Vi; 

				// 2. Bi: the xyz coordinate of the nmtr of equation 4 on page 704 of FACTS paper
				facts2.nmtr_[atm2]  -= (Vi/dis2)*dxyz;
				facts2.dnmtr_[atm2] +=  Vi/dis; // the dnmtr of equation 4 on page 704 of FACTS paper
			}
		}
	}
}//end res_res_burial_for_scoring

// Converts Ai & Bi (actually its nmtr & dnmtr) into Born Radius & SASA
void FACTSPotential::get_self_terms(
																		FACTSRsdTypeInfoCOP factstype1,
																		FACTSResidueInfo & facts1,
																		bool const packing
																		) const
{
	for( Size atm1 = 1; atm1<=factstype1->natoms(); atm1++ ){
		// Confirm B denominator isn't any strange value
		if( !(facts1.dnmtr(atm1) >= 1.0 ) ) facts1.dnmtr_[atm1] = 1.0;

		if (facts1.flag_for_calculation(atm1) ){
			// Bi needs to be calculated after all vectors are collected in res_res_burial
			// Bi in equation 4 of page 704 of FACTS paper
			facts1.Bi_[atm1] = (facts1.nmtr(atm1)).norm() / facts1.dnmtr(atm1); 

			// 1. Born Radius & Self Polar energy!
			facts1.Ci_[atm1] = facts1.Ai(atm1) + (factstype1->b1( atm1 ) * facts1.Bi(atm1)) +
				(factstype1->b2(atm1)  * facts1.Ai(atm1) * facts1.Bi(atm1) );

			Real expterm = exp( -factstype1->a2(atm1) *	(facts1.Ci(atm1) - factstype1->a3(atm1)) );
			Real tmp = factstype1->a1(atm1)/(1.0 + expterm);

			// Equation 7 of page 704 of FACTS paper
			facts1.esolvE_[atm1] = factstype1->a0(atm1) + tmp;

			facts1.BR_[atm1] = -0.5*Tau()*MultiplicitiveFactor()/facts1.esolvE_[atm1];

			// Save partial derivative for nonpolar term here
			if( !packing )
				facts1.dG_dCi_[atm1] = factstype1->a2(atm1)*expterm*tmp/( 1.0 + expterm );

			// 2. SASA!
			// Equation 10 of page 706 of FACTS paper
			facts1.Di_[atm1] = facts1.Ai(atm1) + (factstype1->d1( atm1 ) * facts1.Bi(atm1)) +
				(factstype1->d2(atm1)  * facts1.Ai(atm1) * facts1.Bi(atm1) );

			// Equation 11 of page 706 of FACTS paper
			expterm = exp( -factstype1->c2(atm1) * (facts1.Di(atm1) - factstype1->c3(atm1)) );
			tmp = factstype1->c1(atm1)/(1.0 + expterm);

			facts1.sasa_[atm1] = factstype1->c0(atm1) + tmp;

			// Save partial derivative for nonpolar term here
			if( !packing )
				facts1.dSA_dDi_[atm1] = factstype1->c2(atm1)*expterm*tmp/( 1.0 + expterm );

			// Assert before proceed
			runtime_assert( facts1.BR(atm1) < 1.0e10 );
			runtime_assert( facts1.BR(atm1) > 1.0e-1 );
			runtime_assert( std::abs( facts1.esolvE(atm1) ) < 1.0e4 );
			runtime_assert( facts1.sasa(atm1) < 1.0e10 );
			runtime_assert( facts1.sasa(atm1) > 0.0 );
		}
	}
}

// Calculate polar interaction between rsd1 & rsd2 using Born Radius information
void FACTSPotential::calculate_GBpair_exact( 
												 conformation::Residue const & rsd1,
												 conformation::Residue const & rsd2,
												 FACTSResidueInfo & facts1, 
												 FACTSResidueInfo & facts2
												 ) const 
{
	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );
	bool adjacent = ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) );

	Real cut_off_square = GBpair_cut_ * GBpair_cut_;
	Real GBpair = 0.0;

	using namespace core::scoring::etable::count_pair;

	CountPairFunctionOP cpfxn13 =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
	CountPairFunctionOP cpfxn14 =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	FACTSRsdTypeInfoCOP factstype1 = get_factstype( rsd1.type() );
	FACTSRsdTypeInfoCOP factstype2 = get_factstype( rsd2.type() );

	for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++atm1 ) {
		Vector const &xyz1 = rsd1.xyz(atm1);
		Real const &q1 = factstype1->q(atm1);

		if( !facts1.flag_for_calculation(atm1) ||
				factstype1->not_using(atm1) || std::fabs( q1 ) < 1.0e-6 ) continue;

		Real dE_dBR = 0.0;
		Size path_dist = 0;
		Real selfscore = 0.0;
		Real pairscore = 0.0;

		for ( Size atm2 = 1; atm2 <= rsd2.natoms(); ++atm2 ) {
			Real const &q2 = factstype2->q( atm2 );

			if( !facts2.flag_for_calculation(atm2) ||
					factstype2->not_using(atm2) || std::fabs( q2 ) < 1.0e-6 ) continue;

			/// Hacking part
			// Note: self_pair precedes everything, then same_res and adjacent

			// 1. Selfpair definition: up to 1-3 (connected by angle) whatever respair relation is
			Real dumm = 1.0;
			bool self_pair = 
				(adjacent && !(cpfxn13->count( atm1, atm2, dumm, path_dist ))) ||
				(same_res && factstype1->selfpair(atm1,atm2));

			// Adjacent respair 1-4, 1-5: special care to avoid overlap with Rama term
			// 1-4 weight is 0.0, 1-5 weight is 0.2, and else 1.0
			// This is just to be consistent with the "hacked way" in hack_elec
			Real cpweight = 1.0;
			if( adjacent ){
				bool is_cp = cpfxn14->count( atm1, atm2, cpweight, path_dist );
				if( !self_pair && !is_cp ) cpweight = 0.0; // 1-4
			}

			// Collect the relationship information to correct the scales
			Real scale( 1.0 );
			if( self_pair ){
				scale = selfenergy_scale_;
			} else if ( same_res ){
				scale = intrares_scale_;
			} else if ( adjacent ){
				scale = cpweight; // This would vary from 0.0 to 1.0
			}

			// Start evaluation
			Vector const &xyz2 = rsd2.xyz( atm2 );
			Real dis2 = xyz1.distance_squared( xyz2 );

			if ( !(dis2 < cut_off_square) ) continue;

			Real dis = std::sqrt(dis2);
			Vector dxyz;
			Real BRi = facts1.BR(atm1);
			Real BRj = facts2.BR(atm2);
			Real dshift2( 0.0 );

			// Fading short-distance solvation effect, to be in harmonied with hack_elec
			// Consider self-energy (and pseudo self-energy b/w 1-2, 1-3) to be free from hacking.
			if( self_pair ){
				dshift2 = 0.0;
				dxyz = xyz1 - xyz2;

			} else {
				if ( dis < min_dis_ ){
					dis = min_dis_;
					dis2 = dis*dis;
					dxyz[0] = 0.0;
					dxyz[1] = 0.0;
					dxyz[2] = 0.0;
				} else {
					dxyz = xyz1 - xyz2;
				}

				dshift2 = std::min( dshift2_, dis2 );
				if( factstype1->is_chargedH(atm1) ) BRi *= saltbridge_correction_;
				if( factstype2->is_chargedH(atm2) ) BRj *= saltbridge_correction_;
			}

			// Main run
			Real BRij = BRi*BRj;
			Real tmp1 = dis2/Kappa();
			Real tmp2 = exp(-tmp1/BRij);
			Real tmp3 = dis2 + BRij*tmp2 - dshift2;

			// To avoid dividing by zero 
			if ( !(tmp3 > 1.0e-3) ) continue;

			Real fpair = scale*MultiplicitiveFactor()*Tau()*(q1*q2/sqrt(tmp3));

			// Shift function (required for truncation at cut_off)
			Real sf1 = 1.0 - dis2/cut_off_square;
			Real sf2 = sf1*sf1;

			// Derivative stuffs
			Real g1 = 0.5*fpair/tmp3;
			Real g2 = 2.0 - 2.0*tmp2/Kappa();
			Real dE_drij = g1*g2;
			Real dsf2_drij = 4.0*sf1/cut_off_square;
			Real dEsf_drij = dE_drij*sf2 + fpair*dsf2_drij;

			/* // just for debugging
			if( self_pair ){
				if( same_res ){
					selfscore -= 0.5*fpair*sf2;
				} else {
					selfscore -= fpair*sf2;
				}
			} else {
				pairscore -= fpair*sf2;
			}
			*/

			//string a1 = rsd1.atom_name(atm1);
			//string a2 = rsd2.atom_name(atm2);

			// Make sure that fpair isn't any weird value 
			if( fpair != fpair || std::abs(fpair) > 1.0e6 ){
				TR << "Bad pair interaction score! " << fpair << std::endl;
				continue;
			}

			// Store here and reuse at res_pair scoring & derivative call
			if ( same_res ){
				GBpair -= 0.5*fpair*sf2;

				/*
				if( std::abs( fpair*sf2 ) > 2.0 && !self_pair ){
					TR.debug << "Intra: " << std::setw(4) << rsd1.seqpos() << " " << a1;
					TR.debug << " " << std::setw(4) << rsd2.seqpos() << " " << a2;
					TR.debug << " " << -0.5*fpair*sf2 << std::endl;
				}
				*/

				facts1.dE_drij2_[atm1] += dEsf_drij*dxyz;
				facts1.dE_dBR_[atm1] -= 0.5*sf2*g1*tmp2*(BRj + tmp1/BRi);
				facts2.dE_dBR_[atm2] -= 0.5*sf2*g1*tmp2*(BRi + tmp1/BRj);

			} else {
				GBpair -= fpair*sf2;

				/*
				if( adjacent ){
					TR.Debug << "Facts Pair.dres1: " << std::setw(4) << rsd1.seqpos() << " " << a1;
					TR.Debug << " " << std::setw(4) << rsd2.seqpos() << " " << a2;
					TR.Debug << " " << -0.5*fpair*sf2;
					TR.Debug << " " << std::setw(4) << rsd1.atom_name( atm1 );
					TR.Debug << " " << std::setw(4) << rsd2.atom_name( atm2 );
					TR.Debug << " " << cpweight << " " << dumm;
					TR.Debug << " " << std::endl;
				}
				*/

				facts1.dE_drij2_[atm1] += dEsf_drij*dxyz;
				facts2.dE_drij2_[atm2] -= dEsf_drij*dxyz;
				facts1.dE_dBR_[atm1] -= sf2*g1*tmp2*(BRj + tmp1/BRi);
				facts2.dE_dBR_[atm2] -= sf2*g1*tmp2*(BRi + tmp1/BRj);
			}
		}//atm2

		if( same_res ){
			string atypename = rsd1.atom_type(atm1).atom_type_name();
			string atmname = rsd1.atom_name(atm1);
			TR.Debug << "Facts Self/BR/Ai/Bi/Ci/Di: " << std::setw(4) << rsd1.seqpos() << " " << rsd1.name();
			TR.Debug << " " << atmname << " " << atypename;
			TR.Debug << " " << std::setw(10) << selfscore;
			TR.Debug << " " << std::setw(10) << facts1.BR(atm1);
			TR.Debug << " " << std::setw(10) << facts1.Ai(atm1);
			TR.Debug << " " << std::setw(10) << facts1.Bi(atm1);
			TR.Debug << " " << std::setw(10) << facts1.Ci(atm1);
			TR.Debug << " " << std::setw(10) << facts1.Di(atm1);
			TR.Debug << std::endl;
		}
	}//atm1

	facts1.GBpair_[rsd2.seqpos()] = GBpair;

	if( !same_res ){
		TR.Debug << "Facts Pair:" << std::setw(4) << rsd1.seqpos() << " " << rsd1.name();
		TR.Debug << " " << rsd2.seqpos() << " " << rsd2.name();
		TR.Debug << " " << std::setw(10) << GBpair;
		TR.Debug << std::endl;
	}

} //END FACTSPotential::calculate_GBpair

// Calculate polar interaction between rsd1 & rsd2 using Born Radius information
// Note that this function is applied when rsd1=rsd2; use exact instead to reduce error
void FACTSPotential::calculate_GBpair_apprx(
												 conformation::Residue const & rsd1,
												 conformation::Residue const & rsd2,
												 FACTSResidueInfo & facts1,
												 FACTSResidueInfo & facts2 ) const
{
	Real cut_off_square = GBpair_cut_ * GBpair_cut_;
	Real GBpair = 0.0;
	Real Kappa_sqrt = sqrt(Kappa());
	Real a12 = 0.865924478668;
	Real b12 = 0.4399086973;

	FACTSRsdTypeInfoCOP factstype1 = get_factstype( rsd1.type() );
	FACTSRsdTypeInfoCOP factstype2 = get_factstype( rsd2.type() );

	for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++atm1 ) {
		Vector const &xyz1 = rsd1.xyz(atm1);
		Real const &q1 = factstype1->q(atm1);

		if( factstype1->not_using(atm1) || std::fabs( q1 ) < 1.0e-6 ) continue;
		Real dE_dBR = 0.0;

		for ( Size atm2 = 1; atm2 <= rsd2.natoms(); ++atm2 ) {
			Real const &q2 = factstype2->q( atm2 );

			if( factstype2->not_using(atm2) || std::fabs( q2 ) < 1.0e-6 ) continue;

			Vector const &xyz2 = rsd2.xyz( atm2 );
			Real dis2 = xyz1.distance_squared( xyz2 );

			if ( !(dis2 < cut_off_square ) ) continue;

			Real const &BRi = facts1.BR(atm1);
			Real const &BRj = facts2.BR(atm2);

			Real BRij = BRi*BRj;
			Real m = sqrt(BRi*BRj);
			Real dis = std::sqrt(dis2);
			Vector dxyz = xyz1 - xyz2;

			Real fpair_exact = -MultiplicitiveFactor()*Tau()*q1*q2/
				sqrt(dis2+BRij*exp(-dis2/(Kappa()*BRij)))*(1.0-dis2/cut_off_square)*(1.0-dis2/cut_off_square);

			if ( dis < m ){ // Inside approximator
				Real arg = 1.0/((a12+b12)*m);
				Real sfm = (1.0 - BRij/cut_off_square);
				Real fm = arg*sfm*sfm;
				Real gm = -arg*(sfm*sfm*a12*arg + sfm*4.0*m/cut_off_square);

				Real c = (gm - 2.0*fm/m + 2.0/BRij)/BRij;
				Real d = (-gm + 3.0*fm/m - 3.0/BRij)/m;

				Real fpair = -MultiplicitiveFactor()*Tau()*q1*q2*(c*dis2*dis + d*dis2 + 1.0/m);
				Real gpair = -MultiplicitiveFactor()*Tau()*q1*q2*(c*dis + d);

				// Store here and reuse at res_pair scoring & derivative call
				GBpair += fpair;
				facts1.dE_drij2_[atm1] += gpair*dxyz;
				facts2.dE_drij2_[atm2] -= gpair*dxyz;

				// Forget about dE_dBR for now... too messy
				//facts1.dE_dBR_[atm1] += ;
				//facts2.dE_dBR_[atm2] += ;

			}	else { // Outside approximater
				Real sf1 = 1.0 - dis2/cut_off_square;
				Real dsf2_drij = 4.0*sf1/cut_off_square;
				Real farg = a12*dis + b12*m;
				Real fpair = -MultiplicitiveFactor()*Tau()*(q1*q2/farg);
				Real farg2 = fpair*sf1*sf1/farg;

				// Shift function (required for truncation at cut_off)
				Real dEsf_drij = -farg2/dis - fpair*dsf2_drij;

				// Store here and reuse at res_pair scoring & derivative call
				GBpair += fpair*sf1*sf1;

				facts1.dE_drij2_[atm1] += dEsf_drij*dxyz;
				facts2.dE_drij2_[atm2] -= dEsf_drij*dxyz;
				facts1.dE_dBR_[atm1] += 0.5*farg2*b12*BRj / m;
				facts2.dE_dBR_[atm2] += 0.5*farg2*b12*BRi / m;

			}
		}//atm2

	}//atm1

	facts1.GBpair_[rsd2.seqpos()] = GBpair;
} //END FACTSPotential::calculate_GBpair

// Calculate derivatives for both polar & nonpolar interactions
// "res_res_burial, get_self_terms, calculate_GBpair" should precede this function,
// otherwise derivative will be inconsistent to current structure
void FACTSPotential::setup_for_derivatives( pose::Pose & pose ) const
{
	FACTSPoseInfoOP facts_info;
	Vector virtualcrd;
	Vector cross_v;

	//timeval t1, t2;
	//gettimeofday(&t1, NULL );

	if ( pose.data().has( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ) ) {
		facts_info = static_cast< FACTSPoseInfo* > 
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )() );
	} else {
		facts_info = new FACTSPoseInfo();
	}

	// Check whether information is changed from the given pose - if changed, call setup_for_scoring again
	if(	facts_info->is_changed( pose ) ){
		TR.Debug << "Pose changed since last scoring, call setup_for_scoring..." << std::endl;
		setup_for_scoring( pose, false );
		facts_info = static_cast< FACTSPoseInfo* > 
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )() );
	}

	//
	for ( Size res1 = 1; res1 <= facts_info->size(); ++res1){
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );
		FACTSRsdTypeInfoCOP factstype1 = get_factstype( pose.residue(res1).type() );

		for ( Size atm1 = 1; atm1 <= facts1.natoms(); ++atm1){
			Vector const &crd1 = pose.residue(res1).xyz(atm1);

			// Warning: Below should be problematic for SASA... must be checked
			if ( (fabs(facts1.esolvE(atm1)) <= 1e-6) || (factstype1->not_using(atm1)) ) continue;

			// 1. Derivative for distance dependent term
			facts1.polarF2d_[atm1] += facts1.dE_drij2_[atm1];

			// 2. Derivative for Born Radius & SASA dependent term
			Real dB_dBdnmtr = -facts1.Bi(atm1)/facts1.dnmtr(atm1);
			Real dB_dBnmtr  = 1.0/(facts1.nmtr(atm1).norm() * facts1.dnmtr(atm1));
			Real dBR_dG = facts1.BR(atm1)/facts1.esolvE(atm1);

			Real CutOff_sqr = factstype1->COradius2(atm1);
			
			SelfNeighInfo const &neigh = facts1.selfneigh(atm1);
			Size const &n_neigh = neigh.nneigh;

			// Self-Neighbor needs to be defined (this is defined at res_res_burial, prior to this function call)
			for (Size i_j = 1; i_j <= n_neigh; ++i_j){
				Size const &res2 = neigh.resID[i_j];
				Size const &atm2 = neigh.atmID[i_j];

				Vector const &crd2 = pose.residue(res2).xyz(atm2);

				FACTSResidueInfo & facts2( facts_info->residue_info( res2 ) );
				FACTSRsdTypeInfoCOP factstype2 = get_factstype( pose.residue(res2).type() );

				Vector const dxyz( crd1 - crd2 );
				Real const dis2 ( dxyz.length_squared() );
				Real const dis ( std::sqrt(dis2) );
				Real const i_dis2 ( 1.0/dis2 );

				Real theta_sqrt = 1.0 - (dis2 / CutOff_sqr);
				Real dtheta_tmp1 = factstype2->volume(atm2) * theta_sqrt*theta_sqrt;
				Real dtheta_tmp2 = 4.0*factstype2->volume(atm2) *theta_sqrt/CutOff_sqr;
				Real tmp1 = dtheta_tmp1*i_dis2; // = theta*Vj*xij/r ;
				Real tmp2 = (tmp1 + dtheta_tmp2)/dis;

				// Derivative for Ai
				Real dAi_drij =  -dtheta_tmp2;

				// Derivative for Bi
				Real arg12 = -2.0*dtheta_tmp1*i_dis2 - dtheta_tmp2;
				Vector argv1 = i_dis2*arg12*dxyz;
				Real dBn_drij = argv1[0]*facts1.nmtr(atm1)[0] + argv1[1]*facts1.nmtr(atm1)[1] + argv1[2]*facts1.nmtr(atm1)[2];
				Vector dBn_drij2 = facts1.nmtr(atm1)*dtheta_tmp1*i_dis2;
				Vector dBi_drij = (dB_dBnmtr*dBn_drij - dB_dBdnmtr*tmp2)*dxyz + dB_dBnmtr*dBn_drij2;

				// 2-1. Derivatives for Polar Interaction
				Vector dCi_drij_for_F2 = dAi_drij*dxyz + factstype1->b1(atm1)*dBi_drij + 
					factstype1->b2(atm1)*(facts1.Ai(atm1)*dBi_drij + facts1.Bi(atm1)*dAi_drij*dxyz);

				Vector dE_drij_for_F2 = facts1.dE_dBR(atm1)*dBR_dG*facts1.dG_dCi(atm1)*dCi_drij_for_F2;

				// Assert that derivative isn't any weird value
				Real drv_dot1( dE_drij_for_F2.dot( dE_drij_for_F2 ) );
				if( drv_dot1 != drv_dot1 || drv_dot1 > 1.0e6 ){
					TR << "Bad polar derivatives! " << drv_dot1 << std::endl;
					continue;
				} else {
					facts1.polarF2BR_[atm1] += dE_drij_for_F2;
					facts2.polarF2BR_[atm2] -= dE_drij_for_F2;
				}

				// 2-2. Derivatives for NonPolar Interaction
				Vector dDi_drij_for_F2 = dAi_drij*dxyz + factstype1->d1(atm1)*dBi_drij + 
					factstype1->d2(atm1)*(facts1.Ai(atm1)*dBi_drij + facts1.Bi(atm1)*dAi_drij*dxyz);

				Vector dSA_drij_for_F2 = factstype1->alpha(atm1)*facts1.dSA_dDi(atm1)*dDi_drij_for_F2;

				// Assert that derivative isn't any weird value
				Real drv_dot2( dSA_drij_for_F2.dot( dSA_drij_for_F2 ) );
				if( drv_dot2 != drv_dot2 || drv_dot2 > 1.0e6 ){
					TR << "Bad nonpolar derivatives! " << drv_dot2 << std::endl;
					continue;
				} else {
					facts1.nonpolarF2_[atm1] += dSA_drij_for_F2;
					facts2.nonpolarF2_[atm2] -= dSA_drij_for_F2;
				}
			}
		}
	}

	pose.data().set( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO, facts_info );
	PROF_STOP( basic::FACTS_GET_ALL_BORN_RADII );

	//gettimeofday(&t2, NULL );
	//Real elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
	//elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;
	//cout << "Setup for derivative 1/2: " << elapsedTime << " ms." << endl;

}

// Called at scoring step - for polar energy
// Just reuse scores calculated at setup_for_scoring
Real FACTSPotential::evaluate_polar_energy(Residue const & rsd1,
																					 FACTSResidueInfo const & facts1,
																					 Residue const & rsd2
																					 ) const {

	return facts1.GBpair( rsd2.seqpos() );
}

// Called at scoring step - for nonpolar energy
// Just reuse scores calculated at setup_for_scoring
Real FACTSPotential::evaluate_nonpolar_energy(Residue const & rsd1,
																							FACTSResidueInfo const & facts1,
																							Residue const & rsd2
																							) const {
	Real E_SA = 0.0;
	//FACTSRsdTypeInfo const factstype1 = get_factstype( rsd1.type() );
	FACTSRsdTypeInfoCOP factstype1 = get_factstype( rsd1.type() );

	if ( rsd1.seqpos() == rsd2.seqpos() ){
		for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++ atm1 ){
			E_SA += factstype1->alpha(atm1)*facts1.sasa(atm1);

			string atmname = rsd1.atom_type(atm1).atom_type_name();
			TR.Debug << "Facts SA: " << std::setw(4) << rsd1.seqpos()  << " " << rsd1.name();
			TR.Debug << " " << atmname;
			TR.Debug << " " << std::setw(10) << facts1.sasa(atm1) << " " << std::setw(10) << facts1.sasa(atm1)*factstype1->alpha(atm1);
			TR.Debug << std::endl;
		}
	}
	return E_SA;
}

void FACTSPotential::eval_atom_polar_derivative(
																					id::AtomID const & id,
																					Real const weight,
																					pose::Pose const & pose,
																					kinematics::DomainMap const & domain_map,
																					bool const exclude_DNA_DNA,
																					Vector & F1,
																					Vector & F2
																					) const
{
	FACTSPoseInfo const & facts_info( static_cast< FACTSPoseInfo const & >
																		( pose.data().get( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )));
	Size const atm1( id.atomno() );
	Size const res1( id.rsd() );

	// Pose stuff - this is not supported currently
	//int const i_map( domain_map( res2 ) ); 
	//bool const i_fixed( i_map != 0 );

	FACTSResidueInfo const & facts1( facts_info.residue_info( res1 ) );

	Vector tmpv = weight * (facts1.polarF2d(atm1) + facts1.polarF2BR(atm1));
	F2 += tmpv;

	Vector const &crd1 = pose.residue(res1).xyz(atm1);
	Vector virtualcrd = -tmpv + crd1;
	F1 += crd1.cross( virtualcrd );
}

// Called during minimization; Just call derivatives calculated at setup_for_derivative
void FACTSPotential::eval_atom_nonpolar_derivative(
																					id::AtomID const & id,
																					Real const weight,
																					pose::Pose const & pose,
																					kinematics::DomainMap const & domain_map,
																					bool const exclude_DNA_DNA,
																					Vector & F1,
																					Vector & F2
																					) const
{
	FACTSPoseInfo const & facts_info( static_cast< FACTSPoseInfo const & >
																		( pose.data().get( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )));
	Size const atm1( id.atomno() );
	Size const res1( id.rsd() );

	// Pose stuff - this is not supported currently
	//int const i_map( domain_map( res2 ) ); 
	//bool const i_fixed( i_map != 0 );

	FACTSResidueInfo const & facts1( facts_info.residue_info( res1 ) );

	Vector tmpv = weight * facts1.nonpolarF2(atm1);
	F2 += tmpv;

	Vector const &crd1 = pose.residue(res1).xyz(atm1);
	Vector virtualcrd = -tmpv + crd1;
	F1 += crd1.cross( virtualcrd );
}

/// Note: when called at the beginning of rotamer_trials, task.being_packed(i) will be false for all i
/// this ensures that we use all the information we have to compute the current set of radii
void FACTSPotential::setup_for_packing(
  pose::Pose & pose,
	utility::vector1< bool > const & repacking_residues ) const
{
	PROF_START( basic::FACTS_SETUP_FOR_PACKING );

	FACTSPoseInfoOP facts_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ) ) {
		facts_info = static_cast< FACTSPoseInfo* >
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )() );
	} else {
		facts_info = new FACTSPoseInfo();
		setup_for_scoring( pose, true );
	}

	/// store info about which positions are moving
	facts_info->set_repack_list( repacking_residues );

	/*
	build_placeholders( pose, *facts_info );
	// For non-packed residues
	get_template_born_radii( pose, *facts_info );
	*/

	pose.data().set( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO, facts_info );
	PROF_STOP( basic::FACTS_SETUP_FOR_PACKING );
}

void FACTSPotential::get_template_born_radii(pose::Pose const & pose, FACTSPoseInfo & facts_info) const{
	Size const nres( pose.total_residue() );
	runtime_assert( facts_info.size() == nres );

	for ( Size i=1; i<= nres; ++i ) {
		if ( facts_info.being_packed( i ) ) continue;

		Residue const & rsd1( pose.residue( i ) );
		FACTSResidueInfo & facts1( facts_info.residue_info( i ) );
		//FACTSRsdTypeInfo const factstype1 = get_factstype( rsd1.type() );
		FACTSRsdTypeInfoCOP factstype1 = get_factstype( rsd1.type() );
		runtime_assert( rsd1.natoms()<1 || std::fabs(facts1.Ai(1)) < 1e-3 );

		for ( Size j=1; j<= nres; ++j ) {
			// we are not using placeholder in FACTS - this can be changed if we want to do "Design"
			//if ( facts_info.being_packed(j) ) {
			//	res_res_burial( rsd1, facts1, facts_info.placeholder_residue(j), facts_info.placeholder_info(j) );
			//} else {
			res_res_burial( rsd1, facts1, pose.residue(j), facts_info.residue_info(j) );
			//}
		}
		get_self_terms( factstype1, facts1, true );
	}
}

/// called eg after a rotamer substitution is accepted during rotamer trials
void FACTSPotential::update_residue_for_packing(pose::Pose & pose,Size const seqpos) const
{
	FACTSPoseInfo & facts_info( static_cast< FACTSPoseInfo & >
															( pose.data().get( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ) ) );
	FACTSResidueInfo & facts( facts_info.residue_info( seqpos ) );

	Residue const & rsd( pose.residue( seqpos ) );
	facts.initialize( rsd );

	get_single_rotamer_born_radii( rsd, pose, facts_info, facts );
}

void FACTSPotential::get_rotamers_born_radii(pose::Pose const & pose, conformation::RotamerSetBase & rotamer_set) const {

	FACTSPoseInfo const & facts_info_pose( static_cast< FACTSPoseInfo const & >
																				 (pose.data().get( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )));

	// this will get cached in the rotamer set
	// this call should initialize the residue_info objects with the appropriate Residue info
	FACTSRotamerSetInfoOP facts_info_rotamers( new FACTSRotamerSetInfo( rotamer_set ) );

	for ( Size n=1; n<= rotamer_set.num_rotamers(); ++n ) {
		get_single_rotamer_born_radii( *rotamer_set.rotamer(n), 
																	 pose, facts_info_pose, 
																	 facts_info_rotamers->residue_info( n ) );
	}

	rotamer_set.data().set( core::conformation::RotamerSetCacheableDataType::FACTS_ROTAMER_SET_INFO, facts_info_rotamers );
}

void FACTSPotential::get_single_rotamer_born_radii(Residue const & rsd1,
																									 pose::Pose const & pose,
																									 FACTSPoseInfo const & facts_info,
																									 FACTSResidueInfo  & facts1) 	const
{
	Size natoms1, natoms2;
	//FACTSRsdTypeInfo const factstype1 = get_factstype( rsd1.type() );
	FACTSRsdTypeInfoCOP factstype1 = get_factstype( rsd1.type() );

	assert( rsd1.natoms()<1 || std::fabs(facts1.Ai(1)) < 1e-3 );
	for (Size res2=1; res2<= pose.total_residue(); ++res2 ) {
		// we are not using placeholder in FACTS - this can be changed if we want to do "Design"
		//if ( facts_info.being_packed( res2 ) ) {
		//	res_res_burial( rsd1, facts1, facts_info.placeholder_residue( res2 ), facts_info.placeholder_info( res2));
		//} else {
		res_res_burial( rsd1, facts1, pose.residue( res2 ), facts_info.residue_info( res2 ) );
		//}
		
		get_self_terms( factstype1, facts1, true );

	}//end for loop 1
}

// Given precalculated born radius, called at packing - for polar energy
Real FACTSPotential::evaluate_polar_otf_energy(Residue const & rsd1,
																							 FACTSResidueInfo const & facts1,
																							 Residue const & rsd2,
																							 FACTSResidueInfo const & facts2,
																							 utility::vector1< Real > const & dBRi1,
																							 utility::vector1< Real > const & dBRi2,
																							 bool do_correction
																							 ) const {

	Real score = 0.0;
	Real cut_off_square = GBpair_cut_ * GBpair_cut_;
	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );
	bool adjacent = ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) );

	Real BRi, BRj;

	using namespace core::scoring::etable::count_pair;

	CountPairFunctionOP cpfxn13 =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
	CountPairFunctionOP cpfxn14 =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	FACTSRsdTypeInfoCOP factstype1 = get_factstype( rsd1.type() );
	FACTSRsdTypeInfoCOP factstype2 = get_factstype( rsd2.type() );

	for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++atm1 ) {
		Vector xyz1 = rsd1.xyz(atm1);
		Real const &q1 = factstype1->q(atm1);

		if( factstype1->not_using(atm1) || std::fabs( q1 ) < 1.0e-6 ) continue;

		Size path_dist = 0;

		for ( Size atm2 = 1; atm2 <= rsd2.natoms(); ++atm2 ) {
			Real const &q2 = factstype2->q( atm2 );

			if( factstype2->not_using(atm2) || std::fabs( q2 ) < 1.0e-6 ) continue;

			/// Hacking part
			// Note: self_pair precedes everything, then same_res and adjacent

			// 1. Selfpair definition: up to 1-3 (connected by angle) whatever respair relation is
			Real dumm = 1.0;
			bool self_pair = 
				(adjacent && !(cpfxn13->count( atm1, atm2, dumm, path_dist ))) ||
				(same_res && factstype1->selfpair(atm1,atm2));

			// Adjacent respair 1-4, 1-5: special care to avoid overlap with Rama term
			// 1-4 weight is 0.0, 1-5 weight is 0.2, and else 1.0
			// This is just to be consistent with the "hacked way" in hack_elec
			Real cpweight = 1.0;
			if( adjacent ){
				bool is_cp = cpfxn14->count( atm1, atm2, cpweight, path_dist );
				if( !self_pair && !is_cp ) cpweight = 0.0; // 1-4
			}

			// Collect the relationship information to correct the scales
			Real scale( 1.0 );
			if( self_pair ){
				scale = selfenergy_scale_;
			} else if ( same_res ){
				scale = intrares_scale_;
			} else if ( adjacent ){
				scale = cpweight; // This would vary from 0.0 to 1.0
			}

			Vector const &xyz2 = rsd2.xyz( atm2 );
			Real dis2 = xyz1.distance_squared( xyz2 );

			if ( !(dis2 < cut_off_square) ) continue;

			Real dis = std::sqrt(dis2);
			Real BRi = facts1.BR(atm1);
			Real BRj = facts2.BR(atm2);
			Real dshift2( 0.0 );

			// Fading short-distance solvation effect, to be in harmonied with hack_elec
			// Consider self-energy (and pseudo self-energy b/w 1-2, 1-3) to be free from hacking.
			if( self_pair ){
				dshift2 = 0.0;

			}	else {
				if ( dis < min_dis_ ){
					dis = min_dis_;
					dis2 = dis*dis;
				}

				dshift2 = std::min( dshift2_, dis2 );
				if( factstype1->is_chargedH(atm1) ) BRi *= saltbridge_correction_;
				if( factstype2->is_chargedH(atm2) ) BRj *= saltbridge_correction_;
			}

			Real BRij = BRi*BRj;
			Real tmp1 = dis2/Kappa();
			Real tmp2 = exp(-tmp1/BRij);
			Real tmp3 = dis2 + BRij*tmp2 - dshift2;

			// To avoid dividing by zero 
			if ( !(tmp3 > 1.0e-3) ) continue;

			Real fpair = scale*MultiplicitiveFactor()*Tau()*(q1*q2/sqrt(tmp3));

			// Shift function (required for truncation at cut_off)
			Real sf1 = 1.0 - dis2/cut_off_square;
			Real sf2 = sf1*sf1;
			
			if ( same_res ){
				score -= 0.5*sf2*fpair;
			} else {
				score -= sf2*fpair;
			}
		}
	}
	return score;
}

/*
// Given rotamer pair, calculate BR & SA change induced by each other
// This is not being used currently
void FACTSPotential::evaluate_context_change_for_packing(
    Residue const & rsd1_ref,
		Residue const & rsd1,
		FACTSResidueInfo const & facts1,
		Residue const & rsd2_ref,
		Residue const & rsd2,
		FACTSResidueInfo const & facts2,
		utility::vector1< Real > & dBRi1,
		utility::vector1< Real > & dBRi2,
		utility::vector1< Real > & dSAi1,
		utility::vector1< Real > & dSAi2
		) const {

	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );
	Real fct = 1.0;
	if ( same_res ) {
		fct = 0.5;
	}

	Real theta_sqrt, thetaij, thetaij_ref;

	utility::vector1< Real > Ai1 = facts1.Ai();
	utility::vector1< Vector > nmtr1 = facts1.nmtr();
	utility::vector1< Real > dnmtr1 = facts1.dnmtr();
	utility::vector1< Real > Ai2 = facts2.Ai();
	utility::vector1< Vector > nmtr2 = facts2.nmtr();
	utility::vector1< Real > dnmtr2 = facts2.dnmtr();

	for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++atm1 ) {
		Vector xyz1 = rsd1.xyz(atm1);
		Real CutOff_sqr1 = facts1.COradius(atm1)*facts1.COradius(atm1);

		for ( Size atm2 = 1; atm2 <= rsd2.natoms(); ++atm2 ) {
			if ( same_res && atm1 == atm2 ) continue;

			Real CutOff_sqr2 = facts2.COradius(atm2)*facts2.COradius(atm2);
			Vector xyz2 = rsd2.xyz( atm2 );
			Vector dxyz = xyz1 - xyz2;

			Real dis2 = dxyz[0]*dxyz[0] + dxyz[1]*dxyz[1] + dxyz[2]*dxyz[2];
			Real dis = sqrt(dis2);

			// rotamer1
			Vector dxyz_ref1 = xyz1 - rsd2_ref.xyz(atm2);
			Real dis2_ref1 = dxyz_ref1[0]*dxyz_ref1[0] + dxyz_ref1[1]*dxyz_ref1[1] 
				+ dxyz_ref1[2]*dxyz_ref1[2];

			if ( abs(facts1.esolvE(atm1)) > 1e-6 && abs(dis2 - dis2_ref1) > 1e-3 
					 && ( dis2 < CutOff_sqr1 || dis2_ref1 < CutOff_sqr1 )) {
				Real dis_ref = sqrt(dis2_ref1);

				if ( dis2 < CutOff_sqr1 ){
					theta_sqrt = 1.0 - (dis2 / CutOff_sqr1);
					thetaij = theta_sqrt*theta_sqrt;
				}	else {
					thetaij = 0.0;
				}

				if ( dis2_ref1 < CutOff_sqr1 ){
					theta_sqrt = 1.0 - (dis2_ref1 / CutOff_sqr1);
					thetaij_ref = theta_sqrt*theta_sqrt;
				}	else {
					thetaij_ref = 0.0;
				}

				Ai1[atm1] += fct*facts2.volume(atm2)*(thetaij - thetaij_ref);

				nmtr1[atm1] += fct*facts2.volume(atm2)*(thetaij*dxyz/dis2 
																					 - thetaij_ref*dxyz_ref1/dis2_ref1);
				dnmtr1[atm1] += fct*facts2.volume(atm2)*(thetaij/dis - thetaij_ref/dis_ref);

			}

			// rotamer2
			Vector dxyz_ref2 = xyz2 - rsd1_ref.xyz(atm1);
			Real dis2_ref2 = dxyz_ref2[0]*dxyz_ref2[0] + dxyz_ref2[1]*dxyz_ref2[1] 
				+ dxyz_ref2[2]*dxyz_ref2[2];

			if ( abs(facts2.esolvE(atm2)) > 1e-6 && abs(dis2_ref2 - dis2) > 1e-3 
					 && (dis2 < CutOff_sqr2 || dis2_ref2 < CutOff_sqr2) ) {
				Real dis_ref = sqrt(dis2_ref2);

				if ( dis2 < CutOff_sqr2 ){
					theta_sqrt = 1.0 - (dis2 / CutOff_sqr2);
					thetaij = theta_sqrt*theta_sqrt;
				}	else {
					thetaij = 0.0;
				}
				
				if ( dis2_ref2 < CutOff_sqr2 ){
					theta_sqrt = 1.0 - (dis2_ref2 / CutOff_sqr2);
					thetaij_ref = theta_sqrt*theta_sqrt;
				}	else {
					thetaij_ref = 0.0;
				}

				Ai2[atm2] += fct*facts1.volume(atm1)*(thetaij - thetaij_ref);
				nmtr2[atm2] += fct*facts1.volume(atm1)*(thetaij*(-dxyz)/dis2 
																					 - thetaij_ref*dxyz_ref2/dis2_ref2);
				dnmtr2[atm2] += fct*facts1.volume(atm1)*(thetaij/dis - thetaij_ref/dis_ref);
			}

		}
	}

	// Re-evalute Ci and Di
	for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++atm1 ) {
		Real const &Ai = Ai1[atm1];
		Real const &Bi = nmtr1[atm1].norm()/dnmtr1[atm1];

		if ( abs(facts1.esolvE(atm1)) > 1e-6 ){
 			Real Ci = Ai + facts1.b1(atm1)*Bi + facts1.b2(atm1)*Ai*Bi;
			Real Gi = facts1.a0(atm1) + facts1.a1(atm1)/
				(1.0 + exp( -facts1.a2(atm1) *	(Ci - facts1.a3(atm1)) ) );
			dBRi1[atm1] = -0.5*Tau()*MultiplicitiveFactor()/Gi;
		}
		Real dDi = Ai + facts1.d1(atm1)*Bi + facts1.d2(atm1)*Ai*Bi - facts1.Di(atm1);
		dSAi1[atm1] = facts1.alpha(atm1)*facts1.dSA_dDi(atm1)*dDi;

	}

	for ( Size atm2 = 1; atm2 <= rsd2.natoms(); ++atm2 ) {
		Real const &Ai = Ai2[atm2];
		Real const &Bi = nmtr2[atm2].norm()/dnmtr2[atm2];

		if ( abs(facts2.esolvE(atm2)) > 1e-6 ){
			Real Ci = Ai + facts2.b1(atm2)*Bi + facts2.b2(atm2)*Ai*Bi;
			Real Gi = facts2.a0(atm2) + facts2.a1(atm2)/
			(1.0 + exp( -facts2.a2(atm2) *	(Ci - facts2.a3(atm2)) ) );
			dBRi2[atm2] = -0.5*Tau()*MultiplicitiveFactor()/Gi;
		}

		Real dDi = Ai + facts2.d1(atm2)*Bi + facts2.d2(atm2)*Ai*Bi - facts2.Di(atm2);
		dSAi2[atm2] = facts2.alpha(atm2)*facts2.dSA_dDi(atm2)*dDi;

	}
}

Real FACTSPotential::polar_energy_pack_corrector( Residue const & ref_rsd,
																									Residue const & rsd,
																									FACTSResidueInfo const & facts_info
																									) const
{
	Real score_correction = 0.0;

	for( Size iatm = 1; iatm <= rsd.natoms(); ++iatm ){
		Vector const &dxyz = rsd.xyz(iatm) - ref_rsd.xyz(iatm);
		Vector const &drv = facts_info.polarF2BR(iatm); //+ facts_info.polarF2d(iatm);

		Real dabs = dxyz[0]*dxyz[0] + dxyz[1]*dxyz[1] + dxyz[2]*dxyz[2];
		
		if( dabs > 1.0e-3 && dabs < 0.25 ){
			score_correction += drv[0]*dxyz[0] + drv[1]*dxyz[1] + drv[2]*dxyz[2];
		}
	}
	return score_correction;
}
*/

} // namespace scoring
} // namespace core
