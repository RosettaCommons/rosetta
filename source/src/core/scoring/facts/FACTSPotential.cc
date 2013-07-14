// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
 // (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file:   core/scoring/facts/FACTSPotential.cc
/// @brief:  The definitions of 3 classes of the FACTS algorithm resides here (see devel/khorvash/FACTSPotential.hh
/// @author: Hahnbeom Park

// Unit headers
#include <core/scoring/facts/FACTSPotential.hh>

// Project headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
//#include <core/pack/task/PackerTask.hh>
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


// fast math -- from https://code.google.com/p/fastapprox/downloads/detail?name=fastapprox-0.3.2.tar.gz
inline float
fastpow2 (float p) {
	float offset = (p < 0.0) ? 1.0f : 0.0f;
	float clipp = (p < -126.0) ? -126.0f : p;
	int w = (int)(clipp);
	float z = clipp - (float)(w) + offset;
	union { uint32_t i; float f; } v = { (uint32_t) ( (1 << 23) * (clipp + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z) ) };

	return v.f;
}

inline float
fastlog2 (float x) {
	union { float f; uint32_t i; } vx = { x };
	union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
	float y = vx.i;
	y *= 1.1920928955078125e-7f;

	return y - 124.22551499f
			- 1.498030302f * mx.f
			- 1.72587999f / (0.3520887068f + mx.f);
}

inline float fastexp (float p) { return fastpow2 (1.442695040f * p); }
inline float fastpow (float x, float p) { return fastpow2 (p * fastlog2 (x)); }


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
		q_[i] = rsd.atom( i ).charge();

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
			if( COradius2_[i] <= 1.0e-3 && !not_using_[i] ){ // Don't do this for Virtual atoms
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

//This function initializes all the values for FACTS original parameters, atomic volume, Ai, Bi, esolvE, sasa...
void FACTSResidueInfo::initialize(
		conformation::Residue const & rsd,
		FACTSRsdTypeInfoCOP restypeinfo,
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
		dsolv_dBR_ = utility::vector1< Real >( natoms() );
		elecF2_ = utility::vector1< Vector >( natoms(), i );
		solvF2d_ = utility::vector1< Vector >( natoms(), i );
		solvF2BR_ = utility::vector1< Vector >( natoms(), i );
		sasaF2_ = utility::vector1< Vector >( natoms(), i );
	}

	restypeinfo_ = restypeinfo;
} // END void FACTSResidueInfo::initialize( conformation::Residue const & rsd){

void FACTSResidueInfo::refresh_energy_cache( Size const nres ){
	E_elec_ = utility::vector1< Real >( nres, 0.0 );
	E_solv_ = utility::vector1< Real >( nres, 0.0 );
}

void FACTSResidueInfo::store_xyz( Residue const &rsd ){
	xyz_.resize( 0 );
	for( Size iatm = 1; iatm <= rsd.natoms(); ++iatm ){
		xyz_.push_back( rsd.xyz(iatm) );
	}
}


// Constructor
FACTSPoseInfo::FACTSPoseInfo() :
	context_derivative_empty_(true)
{
}

// Constructor
FACTSPoseInfo::FACTSPoseInfo( FACTSPoseInfo const & src ) : CacheableData()
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
	context_derivative_empty_ = src.context_derivative_empty_;
}

void FACTSPoseInfo::initialize( pose::Pose const & pose, FACTSRsdTypeMap &rsdtypemap )
{
	Size const nres( pose.total_residue() ); // maybe faster for symm if it is made symm-aware

	residue_info_.resize( nres, 0 );
	placeholder_residue_.resize( nres, 0 );
	placeholder_info_.resize( nres, 0 );

	for ( Size i=1; i<= nres; ++i ) {
		if ( !residue_info_[i] ) {
			residue_info_[i] = new FACTSResidueInfo();
			residue_info_[i]->set_enumeration_shell( true );
		}

		// Initialize only if the residue is in enumeration_shell
		// otherwise keep information stored previously
		if( residue_info_[i]->enumeration_shell() ){
			// initialize residuetypeinfo if it has not been
			core::chemical::ResidueType const &rsdtype = pose.residue(i).type();
			FACTSRsdTypeMap::const_iterator it = rsdtypemap.find( &rsdtype );

			if ( it == rsdtypemap.end() ) {
				TR << "Adding new FACTS residue type info: " << rsdtype.name() << std::endl;
				FACTSRsdTypeInfoOP rsdtypeinfo = new FACTSRsdTypeInfo;
				rsdtypeinfo->create_info( rsdtype );
				rsdtypemap[ &rsdtype ] = rsdtypeinfo;
				it = rsdtypemap.find( &rsdtype );
			}
			residue_info_[i]->initialize( pose.residue(i), it->second );
		}
	}
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

// Store change in residue level
void
FACTSPoseInfo::update_enumeration_shell( pose::Pose const &pose,
																				 bool const enumerate_second_shell ){

	Energies const & energies( pose.energies() );
	EnergyGraph const & energy_graph( energies.energy_graph() );

	// First check change in coordinate in residue level
	for( Size ires = 1; ires <= pose.total_residue(); ++ ires ){

		FACTSResidueInfo & facts1( residue_info( ires ) );
		facts1.set_changed( false );
		facts1.set_enumeration_shell( false );

		Size const natom( pose.residue(ires).natoms() );
		utility::vector1<Vector> const facts_xyz = residue_info( ires ).xyz();

		// Check residue conformation change by looking at coordinate
		// Is there better way than this?
		if( natom != facts_xyz.size() ){
			facts1.set_changed( true );
			facts1.set_enumeration_shell( true );

		} else {
			for( Size iatm = 1; iatm <= natom; ++ iatm ){
				Vector const dxyz = facts_xyz[iatm] - pose.residue(ires).xyz(iatm);
				Real const d2 = dxyz.dot(dxyz);
				if( d2 > 1.0e-6 ){
					facts1.set_changed( true );
					facts1.set_enumeration_shell( true );
					break;
				}
			}
		}
	}

	// Decide whether to expand enumeration to the second shell:
	// say A-B-C where A,B,C are residues and A-B < 10 A, B-C < 10 A, A-C > 10 A
	// if A's conformation changes, it will affect B's context, which will again change B-C interaction energy
	// In order to calculate energy induced by A's change, 
	// one should enumerate over second shell also otherwise B-C interaction would have error

	if( !enumerate_second_shell ) return;

	// Then iter over neighbors to expand shell where its structure change can
	// affect context difference
	for( Size res1 = 1; res1 <= pose.total_residue(); ++res1 ) {
		FACTSResidueInfo & facts1( residue_info( res1 ) );

		// Propagate change info into its first neighbor shell
		if( facts1.changed() ){
			facts1.set_enumeration_shell( true );
			for ( graph::Graph::EdgeListConstIter
							iru  = energy_graph.get_node( res1 )->const_edge_list_begin(),
							irue = energy_graph.get_node( res1 )->const_edge_list_end();
						iru != irue; ++iru ) {
				Size const res2( (*iru)->get_other_ind( res1 ) );
				FACTSResidueInfo & facts2( residue_info( res2 ) );
				facts2.set_enumeration_shell( true );
			}
		}
	}
}


void FACTSRotamerSetInfo::initialize( RotamerSet const & rotamer_set, FACTSRsdTypeMap &rsdtypemap )
{
	Size const nrot( rotamer_set.num_rotamers() );
	residue_info_.resize( nrot );
	for ( Size i=1; i<= nrot; ++i ) {
		core::chemical::ResidueType const &rsdtype = rotamer_set.rotamer(i)->type();
		FACTSRsdTypeMap::const_iterator it = rsdtypemap.find( &rsdtype );
		if ( it == rsdtypemap.end() ) {
			TR << "Adding new FACTS residue type info: " << rsdtype.name() << std::endl;
			FACTSRsdTypeInfoOP rsdtypeinfo = new FACTSRsdTypeInfo;
			rsdtypeinfo->create_info( rsdtype );
			rsdtypemap[ &rsdtype ] = rsdtypeinfo;
			it = rsdtypemap.find( &rsdtype );
		}

		residue_info_[i] = new FACTSResidueInfo( *rotamer_set.rotamer(i), it->second, true );
	}
}


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
	min_dis_ = option[ score::facts_min_dis ]();
 	Tau_ = (1.0/option[ score::facts_die ]() ) - (1.0/20.0);
	inv_die_ = 1.0/option[ score::facts_die ]();
	do_apprx = option[ score::facts_apprx ]();
	elec_sh_exponent_ = option[ score::facts_elec_sh_exponent ]();
	selfenergy_scale_ = option[ score::facts_selfenergy_scale ]();
	intrares_scale_ = option[ score::facts_intrares_scale ]();
	saltbridge_correction_ = option[ score::facts_saltbridge_correction ]();
	dshift2_ = option[ score::facts_dshift ]()*option[ score::facts_dshift ]();

	if( do_apprx ) TR << "FACTS Approximation turned on." << endl;

	TR << "FACTS ASP set using: " << option[ score::facts_asp_patch ]()<< endl;
}


//This function "pre-calculates" all the atomic contents (Born radius, sasa,...) and energy values
void FACTSPotential::setup_for_scoring(pose::Pose & pose, bool const & packing) const
{
	Size res1;
	Size res2;

	timeval t1, t2, t3, t4;
	gettimeofday(&t1, NULL );

	PROF_START( basic::FACTS_GET_ALL_BORN_RADII );
	Size const nres( pose.total_residue() );
	FACTSPoseInfoOP facts_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO ) ) {
		facts_info = static_cast< FACTSPoseInfo* >
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )() );
		facts_info->update_enumeration_shell( pose, true );
	} else {
		facts_info = new FACTSPoseInfo();
	}

	facts_info->initialize( pose, FACTSrsdtypemap_ );

	// Set enumeration shell / store current coordinate
	TR.Debug << "Enumeration shell: ";
	Size nshell( 0 );
	for ( res1 = 1; res1 <= nres; ++res1 ) {
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );

		if( packing ){  // Nothing to take care of; just do full enumeration
			facts1.set_enumeration_shell( true );

		} else { // Scoring, minimizing; refresh energy cache for enumeration shell
			if( facts1.enumeration_shell() ){
				facts1.refresh_energy_cache( nres );
				TR.Debug << " " << res1;
				nshell++;
			}
			facts1.store_xyz( pose.residue(res1) );
		}
	}
	TR.Debug << std::endl;
	TR.Debug << "nres_shell/nres_total " << nshell << "/" << nres << std::endl;

	Energies const & energies( pose.energies() );
	EnergyGraph const & energy_graph( energies.energy_graph() );

	// 1. First get Born radius, Solvation energy, SASA for all atoms
	for ( res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd1( pose.residue( res1 ) );
		Size natoms1 = rsd1.natoms();
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );

		// Iter res1 over enumeration_shell only
		if( facts1.natoms() == 0 ) continue;

		if( facts1.enumeration_shell() )
			res_res_burial_for_scoring( rsd1, facts1, rsd1, facts1 );

		// changed again from upperedge to edge to take care of enumeration_shell stuffs
		for ( graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node( res1 )->const_edge_list_begin(),
						irue = energy_graph.get_node( res1 )->const_edge_list_end();
					iru != irue; ++iru ) {
			Size const res2( (*iru)->get_other_ind( res1 ) );
			FACTSResidueInfo & facts2( facts_info->residue_info( res2 ) );

			// pass if both are not in enumeration shell
			if ( !facts1.enumeration_shell() && !facts2.enumeration_shell() ) continue;

			// pass double counting cases
			if ( facts1.enumeration_shell() && facts2.enumeration_shell() && res2 <= res1 ) continue;

			res_res_burial_for_scoring( rsd1, facts1, pose.residue( res2 ), facts2 );
		}
	}

	//gettimeofday(&t2, NULL );
	//Real elapsedTime1 = (t2.tv_sec - t1.tv_sec) * 1000.0;
	//elapsedTime1 += (t2.tv_usec - t1.tv_usec) / 1000.0;

	// 2. Refresh Born radii / SASA
	for ( res1 = 1; res1 <= nres; ++res1 ){
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );

		if( facts1.enumeration_shell() ){
			FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
			get_self_terms( factstype1, facts1, packing );

		}
	}

	//gettimeofday(&t3, NULL );
	//Real elapsedTime2 = (t3.tv_sec - t2.tv_sec) * 1000.0;
	//elapsedTime2 += (t3.tv_usec - t2.tv_usec) / 1000.0;

	// 3. Then Pre-calculate all GB-pair-related parts and store them in FACTSINFO
	Size nenum( 0 );
	Size nenum_full( 0 );
	for ( res1 = 1; res1 <= nres; ++ res1 ){
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );
		Residue const & rsd1( pose.residue( res1 ) );

		if( facts1.enumeration_shell() )
			calculate_GBpair_exact( rsd1, rsd1, facts1, facts1 );

		// changed again from upperedge to edge to take care of enumeration_shell stuffs
		for ( graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node( res1 )->const_edge_list_begin(),
						irue = energy_graph.get_node( res1 )->const_edge_list_end();
					iru != irue; ++iru ) {

			Size const res2( (*iru)->get_other_ind( res1 ) );
			Residue const & rsd2( pose.residue( res2 ) );
			FACTSResidueInfo & facts2( facts_info->residue_info( res2 ) );

			nenum_full ++;
			// option 1(slow). Enumerate all the shell-nonshell pairs 
			if ( !facts1.enumeration_shell() && !facts2.enumeration_shell() ) continue;

			// option 2(fast). Enumerate only shell-shell pairs: this seems like having error on derivatives...
			//if ( !facts1.enumeration_shell() || !facts2.enumeration_shell() ) continue;

			// pass double counting cases
			if ( facts1.enumeration_shell() && facts2.enumeration_shell() && res2 <= res1 ) continue;

			nenum ++;
			calculate_GBpair_exact( rsd1, rsd2, facts1, facts2 );
		}
	}
	TR.Debug << "nrespair for GBpair: " << nenum << "/" << nenum_full << std::endl;

	// 3. Finally store everything into pose
	pose.data().set( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO, facts_info );
	PROF_STOP( basic::FACTS_GET_ALL_BORN_RADII );

	gettimeofday(&t4, NULL );
	Real elapsedTime3 = (t4.tv_sec - t3.tv_sec) * 1000.0;
	elapsedTime3 += (t4.tv_usec - t3.tv_usec) / 1000.0;
	//cout << "Setup for scoring 1/2/3: " << elapsedTime1 << " " << elapsedTime2 << " " << elapsedTime3 << " ms." << endl;

}// END setup_for_scoring

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

	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();

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

	Real thetaij, theta_sqrt;
	Real MAX_SELFDCUT2 = min(100.0, GBPair_cut()*GBPair_cut() );

	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();
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
		//SelfNeighInfo &neigh1 = facts1.selfneigh_[atm1];
		Real const &V1 = factstype1->volume(atm1);

		// Iterate only for upper diagonal :)
		for ( Size atm2 = 1; atm2 <= natoms2; ++atm2 ) {
			if (same_res && (atm1 >= atm2)) continue;
			if ( factstype2->not_using(atm2) ) continue;

			Vector const &xyz2 = rsd2.xyz( atm2 );
			Vector const dxyz( xyz1 - xyz2 );
			Real dis2 = dxyz[0]*dxyz[0] + dxyz[1]*dxyz[1] + dxyz[2]*dxyz[2];

			// take special care for distance to prevent from exploding
			if( dis2 != dis2 ){
				continue;
			} else if( dis2 < 0.01 ) {
				dis2 = 0.01;
			} else if ( dis2 > MAX_SELFDCUT2 ) {
				continue;
			}

			Real dis = std::sqrt(dis2);
			Real CutOff_sqr2 = min( factstype2->COradius2(atm2),	MAX_SELFDCUT2 );

			// Consideration for the first atom
			if( facts1.enumeration_shell() && 
					dis2 <= CutOff_sqr1 && factstype2->volume(atm2) > 1e-3 ){
				//if( neigh1.nneigh >= facts1.MAXNEIGH ) continue;
				facts1.flag_for_calculation_[atm1] = true;
				//neigh1.nneigh++;
				//neigh1.resID[ neigh1.nneigh ] = rsd2.seqpos();
				//neigh1.atmID[ neigh1.nneigh ] = atm2;

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
			if( facts2.enumeration_shell() &&
					dis2 <= CutOff_sqr2 && factstype1->volume(atm1) > 1e-3 ){

				//SelfNeighInfo &neigh2 = facts2.selfneigh_[atm2];
				//if (neigh2.nneigh >= facts2.MAXNEIGH) continue;

				facts2.flag_for_calculation_[atm2] = true;
				//neigh2.nneigh++;
				//neigh2.resID[ neigh2.nneigh ] = rsd1.seqpos();
				//neigh2.atmID[ neigh2.nneigh ] = atm1;

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
			runtime_assert( facts1.sasa(atm1) >= 0.0 );
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
	using namespace core::scoring::etable::count_pair;

	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );
	bool adjacent = ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) );

	Real cut_off_square = GBpair_cut_ * GBpair_cut_;
	Real Esolv( 0.0 ), Eelec( 0.0 );

	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();

	CountPairFunctionOP cpfxn13 = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
	CountPairFunctionOP cpfxn14 = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++atm1 ) {
		Vector const &xyz1 = rsd1.xyz(atm1);
		Real const &q1 = factstype1->q(atm1);
		if( !facts1.flag_for_calculation(atm1) || factstype1->not_using(atm1) || std::fabs( q1 ) < 1.0e-6 ) continue;
		Size path_dist( 0 );

		for ( Size atm2 = 1; atm2 <= rsd2.natoms(); ++atm2 ) {
			Real const &q2 = factstype2->q( atm2 );
			if( !facts2.flag_for_calculation(atm2) || factstype2->not_using(atm2) || std::fabs( q2 ) < 1.0e-6 ) continue;

			// 1. Selfpair definition: up to 1-3 (connected by angle) whatever respair relation is
			Real dumm = 1.0;
			bool self_pair =
				(adjacent && !(cpfxn13->count( atm1, atm2, dumm, path_dist ))) ||
				(same_res && factstype1->selfpair(atm1,atm2));

			// Adjacent respair 1-4, 1-5: special care to avoid overlap with Rama term
			// 1-4 weight is 0.0, 1-5 weight is 0.2, and else 1.0
			// This is just to be consistent with the "hacked way" in fa_elec
			Real cpweight = 1.0;
			if( adjacent ){
				bool is_cp = cpfxn14->count( atm1, atm2, cpweight, path_dist );
				if( !self_pair && !is_cp ) cpweight = 0.0; // 1-4
			}

			// Collect the relationship information to correct the scales
			Real scale_solv( 1.0 );
			Real scale_elec( 1.0 );
			if( self_pair ){
				scale_solv = selfenergy_scale_; scale_elec = 0.0;
			} else if ( same_res ){
				scale_solv = intrares_scale_; scale_elec = 0.0;
			} else if ( adjacent ){
				scale_solv = scale_elec = cpweight;
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

			// Fading short-distance solvation effect, to be in harmonied with fa_elec
			// Consider self-energy (and pseudo self-energy b/w 1-2, 1-3) to be free from hacking.
			if( self_pair ){
				dshift2 = 0.0;
				dxyz = xyz1 - xyz2;
			} else {
				if ( dis < min_dis_ ){
					dis = min_dis_;
					dis2 = dis*dis;
					dxyz[0] = dxyz[1] = dxyz[2] = 0.0;
				} else {
					dxyz = xyz1 - xyz2;
				}

				dshift2 = std::min( dshift2_, dis2 );
				if( factstype1->is_chargedH(atm1) ) BRi *= saltbridge_correction_;
				if( factstype2->is_chargedH(atm2) ) BRj *= saltbridge_correction_;
			}

			Real BRij = BRi*BRj;
			Real tmp1 = dis2/Kappa();
			//Real tmp2 = exp(-tmp1/BRij);
			Real tmp2 = (Real)fastexp((float)(-tmp1/BRij));
			Real tmp3 = dis2 + BRij*tmp2 - dshift2;
			if ( !(tmp3 > 1.0e-3) ) continue;

			Real const arg = MultiplicitiveFactor()*q1*q2;
			Real const fsolv = scale_solv*arg*Tau()/std::sqrt(tmp3);
			Real const felec = scale_elec*arg*inv_die()/dis;

			// Shift function (required for truncation at cut_off)
			Real sf1 = 1.0 - dis2/cut_off_square;
			Real sf2 = sf1*sf1;
			//Real sf_elec = std::pow( sf1, elec_sh_exponent_ );
			Real sf_elec = (Real)fastpow( (float)sf1, (float)elec_sh_exponent_ );

			// Derivative stuffs
			Real g1 = 0.5*fsolv/tmp3;
			Real g2 = 2.0 - 2.0*tmp2/Kappa();
			Real dsolv_drij = g1*g2;
			Real delec_drij = -felec/dis2;

			Real dsf2_drij = 4.0*sf1/cut_off_square;
			Real dsf_elec_drij = -2.0*elec_sh_exponent_*sf_elec/(sf1*cut_off_square);

			Real dsolvsf_drij = dsolv_drij*sf2 + fsolv*dsf2_drij;
			Real delecsf_drij = delec_drij*sf_elec + felec*dsf_elec_drij;

			// Make sure that fpair isn't any weird value
			if( fsolv != fsolv || std::abs(fsolv) > 1.0e6 || felec != felec || std::abs(felec) > 1.0e6 ){
				TR << "Bad pair interaction score(fsolv/felec)! " << fsolv << " " << felec << std::endl;
				continue;
			}

			// Store here and reuse at res_pair scoring & derivative call
			if ( same_res ){ // No elec part
				Esolv -= 0.5*fsolv*sf2;

				facts1.solvF2d_[atm1] += dsolvsf_drij*dxyz;
				facts1.dsolv_dBR_[atm1] -= 0.5*sf2*g1*tmp2*(BRj + tmp1/BRi);
				facts2.dsolv_dBR_[atm2] -= 0.5*sf2*g1*tmp2*(BRi + tmp1/BRj);

			} else {
				Esolv -= fsolv*sf2;
				Eelec += felec*sf_elec;

				facts1.elecF2_[atm1] += delecsf_drij*dxyz;
				facts2.elecF2_[atm2] -= delecsf_drij*dxyz;

				facts1.solvF2d_[atm1] += dsolvsf_drij*dxyz;
				facts2.solvF2d_[atm2] -= dsolvsf_drij*dxyz;
				facts1.dsolv_dBR_[atm1] -= sf2*g1*tmp2*(BRj + tmp1/BRi);
				facts2.dsolv_dBR_[atm2] -= sf2*g1*tmp2*(BRi + tmp1/BRj);
			}
		}//atm2
	}//atm1

	// Finally store into facts residue info
	facts1.E_solv_[rsd2.seqpos()] = Esolv;
	facts1.E_elec_[rsd2.seqpos()] = Eelec;
} //END FACTSPotential::calculate_GBpair


// Calculate derivatives for both polar & nonpolar interactions
// Note:
// "res_res_burial, get_self_terms, calculate_GBpair" should precede this function,
// otherwise derivative will be inconsistent to current structure
//
// Then derivative for distance dependent part will be took into account during setup_for_scoring
// This is only for Context-dependent part: Born Radius & SASA dependent terms
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

	// This is to make sure that context-dependent derivative arrays are initialized at least once
	bool full_update( facts_info->context_derivative_empty() );

	// Iter
	// This will iterate full N x N, rather than the upper edge, because cutoff depends on atomic type
	Real MAX_SELFDCUT2 = min(100.0, GBPair_cut()*GBPair_cut() );
	Energies const & energies( pose.energies() );
	EnergyGraph const & energy_graph( energies.energy_graph() );

	for ( Size res1 = 1; res1 <= facts_info->size(); ++res1){
		FACTSResidueInfo & facts1( facts_info->residue_info( res1 ) );
		FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
		core::conformation::Residue const &rsd1 = pose.residue(res1);

		// Pass if res1 is not in enumeration_shell and not full_update
		if( !(full_update || facts1.enumeration_shell() ) ) continue;

		for ( Size atm1 = 1; atm1 <= facts1.natoms(); ++atm1){
			Vector const &crd1 = rsd1.xyz(atm1);

			// Warning: Below can be dangerous for SASA... 
			if ( (fabs(facts1.esolvE(atm1)) <= 1e-6) || (factstype1->not_using(atm1)) ) continue;

			Real dB_dBdnmtr = -facts1.Bi(atm1)/facts1.dnmtr(atm1);
			Real dB_dBnmtr  = 1.0/(facts1.nmtr(atm1).norm() * facts1.dnmtr(atm1));
			Real dBR_dG = facts1.BR(atm1)/facts1.esolvE(atm1);

			Real CutOff_sqr = factstype1->COradius2(atm1);
			Real CutOff_sqr1 = min( factstype1->COradius2(atm1), MAX_SELFDCUT2 );

			// intra-res
			for ( Size atm2 = 1; atm2 <= facts1.natoms(); ++atm2){
				if( atm1 == atm2 ) continue; // itself has nothing to do with its context

				Vector const &crd2 = rsd1.xyz(atm2);
				Vector const dxyz( crd1 - crd2 );

				atom_atom_context_derivative( facts1, facts1, atm1, atm2, dxyz,
						dB_dBdnmtr, dB_dBnmtr, dBR_dG,
						CutOff_sqr1, CutOff_sqr,
						full_update );
			}

			// inter-res
			for ( graph::Graph::EdgeListConstIter
							iru  = energy_graph.get_node( res1 )->const_edge_list_begin(),
							irue = energy_graph.get_node( res1 )->const_edge_list_end();
						iru != irue; ++iru ) {
				Size res2( (*iru)->get_other_ind( res1 ) );

				FACTSResidueInfo & facts2( facts_info->residue_info( res2 ) );
				core::conformation::Residue const &rsd2 = pose.residue(res2);

				for ( Size atm2 = 1; atm2 <= facts2.natoms(); ++atm2){
					Vector const &crd2 = rsd2.xyz(atm2);
					Vector const dxyz( crd1 - crd2 );

					atom_atom_context_derivative( facts1, facts2, atm1, atm2, dxyz,
							dB_dBdnmtr, dB_dBnmtr, dBR_dG,
							CutOff_sqr1, CutOff_sqr,
							full_update );

				}

			} // end inter-res

		} // atm1
	} // res1

	FACTSResidueInfo & facts1( facts_info->residue_info( 1 ) );

	// Save status before storing
	facts_info->context_derivative_empty_ = false;

	pose.data().set( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO, facts_info );
	PROF_STOP( basic::FACTS_GET_ALL_BORN_RADII );
}

// Calculate atom_atom context derivatives for "atm1" brought by "atm2"
// All the information is being stored into facts1 & facts2
void
FACTSPotential::atom_atom_context_derivative( FACTSResidueInfo & facts1,
		FACTSResidueInfo & facts2,
		Size const & atm1,
		Size const & atm2,
		Vector const & dxyz,
		Real const & dB_dBdnmtr,
		Real const & dB_dBnmtr,
		Real const & dBR_dG,
		Real const & CutOff_sqr1,
		Real const & CutOff_sqr,
		bool const full_update
		) const
{
	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();
	Real const dis2 ( dxyz.length_squared() );

	if( dis2 > CutOff_sqr1 || factstype2->volume(atm2) <= 1e-3 ) return;

	Real const dis ( std::sqrt(dis2) );
	Real const i_dis2 ( 1.0/dis2 );

	Real const theta_sqrt = 1.0 - (dis2 / CutOff_sqr);
	Real const dtheta_tmp1 = factstype2->volume(atm2) * theta_sqrt*theta_sqrt;
	Real const dtheta_tmp2 = 4.0*factstype2->volume(atm2) *theta_sqrt/CutOff_sqr;
	Real const tmp1 = dtheta_tmp1*i_dis2; // = theta*Vj*xij/r ;
	Real const tmp2 = (tmp1 + dtheta_tmp2)/dis;
	
	// Derivative for Ai
	Real const dAi_drij = -dtheta_tmp2;

	// Derivative for Bi
	Real const arg12 = -2.0*dtheta_tmp1*i_dis2 - dtheta_tmp2;
	Vector const argv1 = i_dis2*arg12*dxyz;
	Real const dBn_drij = argv1[0]*facts1.nmtr(atm1)[0] + argv1[1]*facts1.nmtr(atm1)[1] + argv1[2]*facts1.nmtr(atm1)[2];
	Vector const dBn_drij2 = facts1.nmtr(atm1)*dtheta_tmp1*i_dis2;
	Vector const dBi_drij = (dB_dBnmtr*dBn_drij - dB_dBdnmtr*tmp2)*dxyz + dB_dBnmtr*dBn_drij2;
		
	// 2-1. Derivatives for Polar Interaction
	Vector const dCi_drij_for_F2 = dAi_drij*dxyz + factstype1->b1(atm1)*dBi_drij +
		factstype1->b2(atm1)*(facts1.Ai(atm1)*dBi_drij + facts1.Bi(atm1)*dAi_drij*dxyz);
	
	Vector const dsolv_drij_for_F2 = facts1.dsolv_dBR(atm1)*dBR_dG*facts1.dG_dCi(atm1)*dCi_drij_for_F2;
	
	// Assert that derivative isn't any weird value
	Real const drv_dot1( dsolv_drij_for_F2.dot( dsolv_drij_for_F2 ) );
	if( drv_dot1 != drv_dot1 || drv_dot1 > 1.0e10 ){
		TR << "Bad solvation derivatives! " << drv_dot1 << std::endl;
		return;
	} else {
		facts1.solvF2BR_[atm1] += dsolv_drij_for_F2;
		if( full_update || facts2.enumeration_shell() ) 
			facts2.solvF2BR_[atm2] -= dsolv_drij_for_F2;
	}
	
	// 2-2. Derivatives for NonPolar Interaction
	Vector dDi_drij_for_F2 = dAi_drij*dxyz + factstype1->d1(atm1)*dBi_drij +
		factstype1->d2(atm1)*(facts1.Ai(atm1)*dBi_drij + facts1.Bi(atm1)*dAi_drij*dxyz);
	
	Vector dSA_drij_for_F2 = factstype1->alpha(atm1)*facts1.dSA_dDi(atm1)*dDi_drij_for_F2;
	
	// Assert that derivative isn't any weird value
	Real drv_dot2( dSA_drij_for_F2.dot( dSA_drij_for_F2 ) );
	if( drv_dot2 != drv_dot2 || drv_dot2 > 1.0e10 ){
		TR << "Bad nonpolar derivatives! " << drv_dot2 << std::endl;
		return;
	} else {
		facts1.sasaF2_[atm1] += dSA_drij_for_F2;
		if( full_update || facts2.enumeration_shell() ) 
			facts2.sasaF2_[atm2] -= dSA_drij_for_F2;
	}
}

// Called at scoring step - for polar energy
// Just reuse scores calculated at setup_for_scoring
void FACTSPotential::evaluate_polar_energy(Residue const & /*rsd1*/,
		 FACTSResidueInfo const & facts1,
		 Residue const & rsd2,
		 Real & E_elec,
		 Real & E_solv
		 ) const {

	E_elec = facts1.E_elec( rsd2.seqpos() );
	E_solv = facts1.E_solv( rsd2.seqpos() );
}

// Called at scoring step - for nonpolar energy
// Just reuse scores calculated at setup_for_scoring
Real FACTSPotential::evaluate_nonpolar_energy(Residue const & rsd1,
		FACTSResidueInfo const & facts1,
		Residue const & rsd2
		) const {
	Real E_SA = 0.0;
	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();

	if ( rsd1.seqpos() == rsd2.seqpos() ){
		for ( Size atm1 = 1; atm1 <= rsd1.natoms(); ++ atm1 ){
			E_SA += factstype1->alpha(atm1)*facts1.sasa(atm1);

			/*
			string atmname = rsd1.atom_type(atm1).atom_type_name();
			TR.Debug << "Facts SA: " << std::setw(4) << rsd1.seqpos()  << " " << rsd1.name();
			TR.Debug << " " << atmname;
			TR.Debug << " " << std::setw(10) << facts1.sasa(atm1) << " " << std::setw(10) << facts1.sasa(atm1)*factstype1->alpha(atm1);
			TR.Debug << std::endl;
			*/
		}
	}
	return E_SA;
}

void FACTSPotential::eval_atom_polar_derivative(
		id::AtomID const & id,
		Real const weight_elec,
		Real const weight_solv,
		pose::Pose const & pose,
		kinematics::DomainMap const & /*domain_map*/,
		bool const /*exclude_DNA_DNA*/,
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

	// Just copy from the scratch saved in facts1 pre-calculated at setup_for_derivates 
	Vector tmpv =
		weight_elec * facts1.elecF2(atm1)
		+ weight_solv * (facts1.solvF2d(atm1) + facts1.solvF2BR(atm1));
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
		kinematics::DomainMap const & /*domain_map*/,
		bool const /*exclude_DNA_DNA*/,
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

	// Just copy from the scratch saved in facts1 pre-calculated at setup_for_derivates
	Vector tmpv = weight * facts1.sasaF2(atm1);
	F2 += tmpv;

	Vector const &crd1 = pose.residue(res1).xyz(atm1);
	Vector virtualcrd = -tmpv + crd1;
	F1 += crd1.cross( virtualcrd );
}

// Note: when called at the beginning of rotamer_trials, task.being_packed(i) will be false for all i
// this ensures that we use all the information we have to compute the current set of radii
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
		FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
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
	core::chemical::ResidueType const &rsdtype = rsd.type();
	FACTSRsdTypeMap::const_iterator it = FACTSrsdtypemap_.find( &rsdtype );

	if ( it == FACTSrsdtypemap_.end() ) {
		TR << "Adding new FACTS residue type info: " << rsdtype.name() << std::endl;
		FACTSRsdTypeInfoOP rsdtypeinfo = new FACTSRsdTypeInfo;
		rsdtypeinfo->create_info( rsdtype );
		FACTSrsdtypemap_[ &rsdtype ] = rsdtypeinfo;
		it = FACTSrsdtypemap_.find( &rsdtype );
	}

	facts.initialize( rsd, it->second );
	get_single_rotamer_born_radii( rsd, pose, facts_info, facts );
}

void FACTSPotential::get_rotamers_born_radii(pose::Pose const & pose, conformation::RotamerSetBase & rotamer_set) const {

	FACTSPoseInfo const & facts_info_pose( static_cast< FACTSPoseInfo const & >
			(pose.data().get( core::pose::datacache::CacheableDataType::FACTS_POSE_INFO )));

	// this will get cached in the rotamer set
	// this call should initialize the residue_info objects with the appropriate Residue info
	FACTSRotamerSetInfoOP facts_info_rotamers( new FACTSRotamerSetInfo( rotamer_set, FACTSrsdtypemap_ ) );

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
	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();

	assert( rsd1.natoms()<1 || std::fabs(facts1.Ai(1)) < 1e-3 );
	for (Size res2=1; res2<= pose.total_residue(); ++res2 ) {
		// we are not using placeholder in FACTS - this can be changed if we want to do "Design"
		//if ( facts_info.being_packed( res2 ) ) {
		//	res_res_burial( rsd1, facts1, facts_info.placeholder_residue( res2 ), facts_info.placeholder_info( res2));
		//} else {
		res_res_burial( rsd1, facts1, pose.residue( res2 ), facts_info.residue_info( res2 ) );
		//}
	}//end for loop 1
	get_self_terms( factstype1, facts1, true );
}

// Given precalculated born radius, called at packing - for polar energy
void FACTSPotential::evaluate_polar_otf_energy(Residue const & rsd1,
		FACTSResidueInfo const & facts1,
		Residue const & rsd2,
		FACTSResidueInfo const & facts2,
		utility::vector1< Real > const & /*dBRi1*/,
		utility::vector1< Real > const & /*dBRi2*/,
		Real & E_elec,
		Real & E_solv,
		bool /*do_correction*/
		) const {

	// Initialize
	E_elec = 0.0;
	E_solv = 0.0;

	Real cut_off_square = GBpair_cut_ * GBpair_cut_;
	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );
	bool adjacent = ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) );

	Real BRi, BRj;

	using namespace core::scoring::etable::count_pair;

	CountPairFunctionOP cpfxn13 =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
	CountPairFunctionOP cpfxn14 =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	FACTSRsdTypeInfoCOP factstype1 = facts1.restypeinfo();
	FACTSRsdTypeInfoCOP factstype2 = facts2.restypeinfo();

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
			// This is just to be consistent with the "hacked way" in fa_elec
			Real cpweight = 1.0;
			if( adjacent ){
				bool is_cp = cpfxn14->count( atm1, atm2, cpweight, path_dist );
				if( !self_pair && !is_cp ) cpweight = 0.0; // 1-4
			}

			// Collect the relationship information to correct the scales
			Real scale_solv( 1.0 );
			Real scale_elec( 1.0 );
			if( self_pair ){
				scale_solv = selfenergy_scale_;
				scale_elec = 0.0;
			} else if ( same_res ){
				scale_solv = intrares_scale_;
				scale_elec = 0.0;
			} else if ( adjacent ){
				scale_solv = cpweight; // This would vary from 0.0 to 1.0
				scale_elec = cpweight;
			}

			Vector const &xyz2 = rsd2.xyz( atm2 );
			Real dis2 = xyz1.distance_squared( xyz2 );

			if ( !(dis2 < cut_off_square) ) continue;

			Real dis = std::sqrt(dis2);
			Real BRi = facts1.BR(atm1);
			Real BRj = facts2.BR(atm2);
			Real dshift2( 0.0 );

			// Fading short-distance solvation effect, to be in harmonied with fa_elec
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
			//Real tmp2 = exp(-tmp1/BRij);
			Real tmp2 = (Real)fastexp((float)(-tmp1/BRij));
			Real tmp3 = dis2 + BRij*tmp2 - dshift2;

			// To avoid dividing by zero
			if ( !(tmp3 > 1.0e-3) ) continue;

			Real const arg = MultiplicitiveFactor()*q1*q2;

			// Shift function (required for truncation at cut_off)
			Real sf1 = 1.0 - dis2/cut_off_square;
			Real sf2 = sf1*sf1;
			//Real sf_elec = std::pow( sf1, elec_sh_exponent_ );
			Real sf_elec = (Real)fastpow( (float)sf1, (float)elec_sh_exponent_ );

			if ( same_res ){
				//E_elec += 0.5*sf2*felec; no elec in same residue
				Real fsolv = scale_solv*arg*Tau()/std::sqrt(tmp3);
				E_solv -= 0.5*sf2*fsolv;

			} else {
				Real fsolv = scale_solv*arg*Tau()/std::sqrt(tmp3);
				Real felec = scale_elec*arg*inv_die()/dis;

				E_elec += sf_elec*felec;
				E_solv -= sf2*fsolv;
			}
		}
	}
}

} // namespace scoring
} // namespace core
