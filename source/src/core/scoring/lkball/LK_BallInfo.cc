// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/LK_BallInfo.cc
/// @brief  Orientation dependent variant of the LK Solvation using
/// @author Phil Bradley

// Unit headers
#include <core/scoring/lkball/LK_BallInfo.hh>

// // Package headers
//#include <core/pack/rotamer_set/WaterPackingInfo.hh>

// // Project headers
#include <core/chemical/RestypeDestructionEvent.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <basic/options/option.hh>


// HACK
#include <core/io/pdb/pdb_writer.hh>

#include <fstream> // HACK
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// #include <core/scoring/constraints/AngleConstraint.hh>

#include <core/chemical/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <core/kinematics/Stub.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>

#include <numeric/conversions.hh>



// #include <utility/vector1.functions.hh> // HACK

//#ifdef WIN32
//#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
//#endif

#ifdef    SERIALIZATION
// Project serialization headers
#include <core/chemical/ResidueType.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/fixedsizearray1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace lkball {

//bool const sidechain_only_hack( false ); // make this configurable

/// @details Auto-generated virtual destructor
LKB_ResidueInfo::~LKB_ResidueInfo() = default;


static basic::Tracer TR("core.scoring.methods.LK_BallInfo" );


/// LAZY
using namespace chemical;
using namespace pose;
using namespace conformation;

static Real const optimal_water_distance( 2.65 ); /// note that this number is re-defined in hbonds.cc !!

#ifdef MULTI_THREADED

utility::thread::ReadWriteMutex LKBallDatabase::lkball_db_mutex_;

#endif

/////////////////////////////////////////////////////////////////////////////

WaterBuilder::WaterBuilder(
	Vector const & water,
	conformation::Residue const & rsd,
	Size const atom1,
	Size const atom2,
	Size const atom3
):
	atom1_( atom1 ),
	atom2_( atom2 ),
	atom3_( atom3 ),
	xyz_local_( kinematics::Stub( rsd.xyz( atom1 ), rsd.xyz( atom2 ), rsd.xyz( atom3 ) ).global2local( water ) )
{
	runtime_assert( atom1 != atom2 && atom1 != atom3 && atom2 != atom3 );
}

bool WaterBuilder::operator==(WaterBuilder const&o) const
{
	return atom1_ == o.atom1_ && atom2_ == o.atom2_ && atom3_ == o.atom3_ && xyz_local_ == o.xyz_local_;
}


Vector
WaterBuilder::build( conformation::Residue const & rsd ) const
{
	return kinematics::Stub( rsd.xyz( atom1_ ), rsd.xyz( atom2_ ), rsd.xyz( atom3_ ) ).local2global( xyz_local_ );
}

// fpd get the derivative of water movement w.r.t. base atom movement
// fpd too slow?
void
WaterBuilder::derivatives(
	conformation::Residue const & rsd,
	numeric::xyzMatrix <Real> &dw_da1,
	numeric::xyzMatrix <Real> &dw_da2,
	numeric::xyzMatrix <Real> &dw_da3 ) const
{
	Real NUM_H = 1e-8; //?
	Vector w000=kinematics::Stub( rsd.xyz( atom1_ ), rsd.xyz( atom2_ ), rsd.xyz( atom3_ ) ).local2global( xyz_local_ );
	Vector w100=kinematics::Stub( rsd.xyz( atom1_ )+Vector(NUM_H,0.0,0.0), rsd.xyz( atom2_ ), rsd.xyz( atom3_ ) ).local2global( xyz_local_ );
	Vector w010=kinematics::Stub( rsd.xyz( atom1_ )+Vector(0.0,NUM_H,0.0), rsd.xyz( atom2_ ), rsd.xyz( atom3_ ) ).local2global( xyz_local_ );
	Vector w001=kinematics::Stub( rsd.xyz( atom1_ )+Vector(0.0,0.0,NUM_H), rsd.xyz( atom2_ ), rsd.xyz( atom3_ ) ).local2global( xyz_local_ );

	Vector w200=kinematics::Stub( rsd.xyz( atom1_ ), rsd.xyz( atom2_ )+Vector(NUM_H,0.0,0.0), rsd.xyz( atom3_ ) ).local2global( xyz_local_ );
	Vector w020=kinematics::Stub( rsd.xyz( atom1_ ), rsd.xyz( atom2_ )+Vector(0.0,NUM_H,0.0), rsd.xyz( atom3_ ) ).local2global( xyz_local_ );
	Vector w002=kinematics::Stub( rsd.xyz( atom1_ ), rsd.xyz( atom2_ )+Vector(0.0,0.0,NUM_H), rsd.xyz( atom3_ ) ).local2global( xyz_local_ );

	Vector w300=kinematics::Stub( rsd.xyz( atom1_ ), rsd.xyz( atom2_ ), rsd.xyz( atom3_ )+Vector(NUM_H,0.0,0.0) ).local2global( xyz_local_ );
	Vector w030=kinematics::Stub( rsd.xyz( atom1_ ), rsd.xyz( atom2_ ), rsd.xyz( atom3_ )+Vector(0.0,NUM_H,0.0) ).local2global( xyz_local_ );
	Vector w003=kinematics::Stub( rsd.xyz( atom1_ ), rsd.xyz( atom2_ ), rsd.xyz( atom3_ )+Vector(0.0,0.0,NUM_H) ).local2global( xyz_local_ );

	dw_da1.col_x( (w100-w000) / NUM_H ); dw_da1.col_y( (w010-w000) / NUM_H ); dw_da1.col_z( (w001-w000) / NUM_H );
	dw_da2.col_x( (w200-w000) / NUM_H ); dw_da2.col_y( (w020-w000) / NUM_H ); dw_da2.col_z( (w002-w000) / NUM_H );
	dw_da3.col_x( (w300-w000) / NUM_H ); dw_da3.col_y( (w030-w000) / NUM_H ); dw_da3.col_z( (w003-w000) / NUM_H );
}

/// The next two functions were taken from protocols/water/rotamer_building_functions.cc with slight modifications
///
Vector
build_optimal_water_O_on_donor(
	Vector const & hxyz,
	Vector const & dxyz
)
{
	core::Real donor_length = optimal_water_distance;
	if ( basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_donor ].user() ) {
		donor_length = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_donor ]();
	}
	return ( dxyz + donor_length * ( hxyz - dxyz ).normalized() );
}


utility::vector1< Vector >
build_optimal_water_Os_on_acceptor(
	Size const acc_atm,
	conformation::Residue const & acc_rsd,
	Vector const & a_xyz,
	Vector b1_xyz, // local copy
	Vector b2_xyz,// local copy
	chemical::Hybridization const & hybrid
) {
	using numeric::conversions::radians;
	using namespace chemical;

	Real const distance( optimal_water_distance ); // acceptor--O distance

	//Real theta(0);
	utility::vector1< Real > params_sp2, params_sp3, params_ring;

	if ( !basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp2 ].user() ) {
		params_sp2.push_back(distance);params_sp2.push_back(120.0);params_sp2.push_back(0.0);
		params_sp2.push_back(distance);params_sp2.push_back(120.0);params_sp2.push_back(180.0);
	} else {
		params_sp2 = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp2 ]();
	}
	if ( !basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp3 ].user() ) {
		params_sp3.push_back(distance);params_sp3.push_back(109.0);params_sp3.push_back(120.0);
		params_sp3.push_back(distance);params_sp3.push_back(109.0);params_sp3.push_back(240.0);
	} else {
		params_sp3 = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_sp3 ]();
	}
	if ( !basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_ring ].user() ) {
		params_ring.push_back(distance);params_ring.push_back(180.0);params_ring.push_back(0.0);
	} else {
		params_ring = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_waters_ring ]();
	}

	// detect special case of DNA phosphate hydration:
	utility::vector1< Vector > waters;

	std::string const & acc_atm_name( acc_rsd.atom_name( acc_atm ) );
	if ( acc_rsd.is_DNA() && acc_rsd.atom_is_backbone( acc_atm ) &&
			( acc_atm_name == " OP2" || acc_atm_name == " OP1" ) ) {

		// special case: hydration of the DNA phosphate group
		b1_xyz = acc_rsd.xyz( "P" );
		// One DNA variant (methylated phosphate) lacks OP2.
		if ( acc_rsd.has( "OP1") && acc_atm_name == " OP2" ) {
			b2_xyz = ( ( acc_atm_name == " OP2" ) ? acc_rsd.xyz( "OP1" ) : acc_rsd.xyz( "OP2" ) );
			kinematics::Stub stub( a_xyz, b1_xyz, b2_xyz );

			// these numbers are taken from "Hydration of the Phosphate Group in Double-Helical DNA",
			// Schneider, Patel, and Berman, Biophysical Journal Vol. 75 2422-2434
			waters.push_back( stub.spherical( radians( 45.0 ), radians( 180.0 - 125.0 ), distance ) );
			waters.push_back( stub.spherical( radians( 160.0 ), radians( 180.0 - 125.0 ), distance ) );
			waters.push_back( stub.spherical( radians( 280.0 ), radians( 180.0 - 125.0 ), distance ) );
		} else if ( acc_rsd.has( "OP2") && acc_atm_name == " OP1" ) {
			b2_xyz = ( ( acc_atm_name == " OP2" ) ? acc_rsd.xyz( "OP1" ) : acc_rsd.xyz( "OP2" ) );
			kinematics::Stub stub( a_xyz, b1_xyz, b2_xyz );

			// these numbers are taken from "Hydration of the Phosphate Group in Double-Helical DNA",
			// Schneider, Patel, and Berman, Biophysical Journal Vol. 75 2422-2434
			waters.push_back( stub.spherical( radians( 45.0 ), radians( 180.0 - 125.0 ), distance ) );
			waters.push_back( stub.spherical( radians( 160.0 ), radians( 180.0 - 125.0 ), distance ) );
			waters.push_back( stub.spherical( radians( 280.0 ), radians( 180.0 - 125.0 ), distance ) );
		}


	} else {

		if ( hybrid == SP2_HYBRID ) {
			if ( params_sp2.size() > 3*MAX_N_WATERS_PER_ATOM ) {
				std::ostringstream oss;
				oss << "More than " << MAX_N_WATERS_PER_ATOM << " waters were requested for atom "
					<< acc_rsd.atom_name( acc_atm ) << " on Residue " << acc_rsd.name() << ". The "
					<< "MAX_N_WATERS_PER_ATOM constant in src/core/scoring/lkball/LK_BallInfo.hh must "
					<< "be increased and the code recompiled to accommodate this many waters.";
				utility_exit_with_message( oss.str() );
			}
			kinematics::Stub stub( a_xyz, b1_xyz, b2_xyz );
			for ( core::Size i=0; i<params_sp2.size()/3; ++i ) {
				waters.push_back( stub.spherical( radians( params_sp2[3*i+3] ), radians( 180.0-params_sp2[3*i+2] ), params_sp2[3*i+1] ) );
			}
		} else if ( hybrid == SP3_HYBRID ) {
			if ( params_sp3.size() > 3*MAX_N_WATERS_PER_ATOM ) {
				std::ostringstream oss;
				oss << "More than " << MAX_N_WATERS_PER_ATOM << " waters were requested for atom "
					<< acc_rsd.atom_name( acc_atm ) << " on Residue " << acc_rsd.name() << ". The "
					<< "MAX_N_WATERS_PER_ATOM constant in src/core/scoring/lkball/LK_BallInfo.hh must "
					<< "be increased and the code recompiled to accommodate this many waters.";
				utility_exit_with_message( oss.str() );
			}
			kinematics::Stub stub( a_xyz, b1_xyz, b2_xyz );
			for ( core::Size i=0; i<params_sp3.size()/3; ++i ) {
				waters.push_back( stub.spherical( radians( params_sp3[3*i+3] ), radians( 180.0-params_sp3[3*i+2] ), params_sp3[3*i+1] ) );
			}
		} else if ( hybrid == RING_HYBRID ) {
			if ( params_ring.size() > 3*MAX_N_WATERS_PER_ATOM ) {
				std::ostringstream oss;
				oss << "More than " << MAX_N_WATERS_PER_ATOM << " waters were requested for atom "
					<< acc_rsd.atom_name( acc_atm ) << " on Residue " << acc_rsd.name() << ". The "
					<< "MAX_N_WATERS_PER_ATOM constant in src/core/scoring/lkball/LK_BallInfo.hh must "
					<< "be increased and the code recompiled to accommodate this many waters.";
				utility_exit_with_message( oss.str() );
			}

			b1_xyz = 0.5 * ( b1_xyz + b2_xyz );
			kinematics::Stub stub( a_xyz, b1_xyz, b2_xyz );
			for ( core::Size i=0; i<params_ring.size()/3; ++i ) {
				waters.push_back( stub.spherical( radians( params_ring[3*i+3] ), radians( 180.0-params_ring[3*i+2] ), params_ring[3*i+1] ) );
			}
		} else {
			TR.Error << "Bad hybridization type for acceptor " << hybrid << std::endl;
			utility_exit();
		}

	}

	return waters;
}


void
setup_water_builders_for_residue_type(
	ResidueType const & rsd_type,
	bool const sidechain_only,
	utility::vector1< WaterBuilders > & rsd_water_builders
)
{
	using namespace conformation;
	using namespace chemical;

	bool const no_SP3_acceptor_waters
		( false );//basic::options::option[ basic::options::OptionKeys::dna::specificity::no_SP3_acceptor_waters ] );

	bool no_acceptor_waters
		( false );//basic::options::option[ basic::options::OptionKeys::dna::specificity::no_acceptor_waters ] );

	bool const no_acceptor_waters_except_OH
		( false );//basic::options::option[ basic::options::OptionKeys::dna::specificity::no_acceptor_waters_except_OH ] );

	bool const lk_ball_subset1
		( false );//basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_subset1 ] );

	bool const dump_waters_pdb( false ); // for debugging purposes

	ResidueOP rsd( ResidueFactory::create_residue( rsd_type ) );

	rsd_water_builders.clear();
	rsd_water_builders.resize( rsd_type.natoms() ); // initialize to empty vector1s

	if ( lk_ball_subset1 ) { // only ARG ASN GLN sidechain donors
		if ( rsd->aa() != chemical::aa_arg &&
				rsd->aa() != chemical::aa_asn &&
				rsd->aa() != chemical::aa_gln ) return;
		no_acceptor_waters = true;
		runtime_assert( sidechain_only ); // sanity
	}

	utility::vector1< Vector > all_waters; // for debugging

	// donors
	for ( auto
			hnum  = rsd->Hpos_polar().begin(),
			hnume = rsd->Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const hatm( *hnum );
		if ( sidechain_only && rsd_type.atom_is_backbone( hatm ) ) continue;

		Size const datm( rsd->atom_base( hatm ) );
		Size datm_base( rsd->atom_base( datm ) );
		if ( datm_base == hatm ) { // could happen with water
			datm_base = 0;
			chemical::AtomIndices const & datm_nbrs( rsd_type.nbrs( datm ) );
			for ( Size ii=1; ii<= datm_nbrs.size(); ++ii ) {
				if ( datm_nbrs[ii] == hatm ) continue;
				else datm_base = datm_nbrs[ii];
			}
			runtime_assert( datm_base );
		}

		Vector const water( build_optimal_water_O_on_donor( rsd->xyz( hatm ), rsd->xyz( datm ) ) );
		rsd_water_builders[ datm ].push_back( WaterBuilder( water, *rsd, hatm, datm, datm_base ) );
		TR.Trace << "adding ideal water to rsd_type= " << rsd_type.name() << " anchored to atom " <<
			rsd_type.atom_name( datm ) << " datm_base= " << rsd_type.atom_name( datm_base ) << std::endl;
		all_waters.push_back( water );
	}

	// acceptors
	if ( !no_acceptor_waters ) {
		for ( auto
				anum  = rsd->accpt_pos().begin(),
				anume = rsd->accpt_pos().end(); anum != anume; ++anum ) {
			Size const aatm( *anum );
			if ( sidechain_only && rsd_type.atom_is_backbone( aatm ) ) continue;
			if ( rsd_type.atom_type( aatm ).hybridization() == chemical::SP3_HYBRID &&
					no_SP3_acceptor_waters ) {
				TR.Trace << "not adding acceptor waters on atom: " << rsd->atom_name( aatm ) << ' ' <<
					rsd->name() << std::endl;
				continue;
			}
			if ( no_acceptor_waters_except_OH && rsd_type.atom_type( aatm ).name() != "OH" ) continue;
			Size const abase1( rsd->atom_base( aatm ) ), abase2( rsd->abase2( aatm ) );
			utility::vector1< Vector > const waters
				( build_optimal_water_Os_on_acceptor( aatm, *rsd,
				rsd->xyz( aatm ), rsd->xyz( abase1 ), rsd->xyz( abase2 ),
				rsd->atom_type( aatm ).hybridization() ) );
			for ( Size i=1; i<= waters.size(); ++i ) {
				rsd_water_builders[ aatm ].push_back( WaterBuilder( waters[i], *rsd, aatm, abase1, abase2 ) );
				TR.Trace << "adding ideal water to rsd_type= " << rsd_type.name() << " anchored to atom " <<
					rsd_type.atom_name( aatm ) << " abase1= " << rsd_type.atom_name( abase1 ) <<
					" abase2= " << rsd_type.atom_name( abase2 ) << std::endl;
				all_waters.push_back( waters[i] );
			}
		}
	}

	/// let's do a sanity check on something we assume down below
	/*
	if ( rsd_type.aa() != aa_h2o ) {
	for ( Size i=1; i<= rsd_type.nheavyatoms(); ++i ) {
	//want to confirm that all waters are matched to a hydrogen
	if ( ! rsd_type.heavyatom_has_polar_hydrogens(i) ) continue;

	WaterBuilders const & water_builders( rsd_water_builders[i] );
	for ( Size j=1; j<= water_builders.size(); ++j ) {
	// look for a polar hydrogen that matches to this guy
	Size hpos(0);
	if ( rsd_type.atom_is_polar_hydrogen( water_builders[j].atom1() ) ) {
	hpos = water_builders[j].atom1();
	}
	if ( rsd_type.atom_is_polar_hydrogen( water_builders[j].atom2() ) ) {
	runtime_assert( !hpos );
	hpos = water_builders[j].atom2();
	}
	if ( rsd_type.atom_is_polar_hydrogen( water_builders[j].atom3() ) ) {
	runtime_assert( !hpos );
	hpos = water_builders[j].atom3();
	}
	runtime_assert( hpos );
	TR.Trace << "donor-water-anchor: " << j << ' ' << rsd_type.name() << ' ' << rsd_type.atom_name(hpos ) <<
	std::endl;
	}
	}
	}
	*/

	if ( dump_waters_pdb && !all_waters.empty() ) { // HACKING -- dump a pdb containing just this residue:
		std::ofstream out( std::string( rsd->name() + "_waters.pdb" ).c_str() );
		Size atom_number( 0 );
		rsd->chain( 1 );
		rsd->seqpos( 1 );
		io::pdb::dump_pdb_residue( *rsd, atom_number, out );

		/// now add all the waters
		char const chain = 'A';
		for ( Size i=1; i<= all_waters.size(); ++i ) {
			using namespace ObjexxFCL::format;
			++atom_number;
			std::string const atom_name( "W" + ObjexxFCL::lead_zero_string_of( i, 3 ) );
			out << "ATOM  " << I(5,atom_number) << ' ' << atom_name << ' ' <<
				rsd->name3() << ' ' << chain << I(4,rsd->seqpos() ) << "    " <<
				F(8,3,all_waters[i](1)) <<
				F(8,3,all_waters[i](2)) <<
				F(8,3,all_waters[i](3)) <<
				F(6,2,1.0) << F(6,2,0.0) << '\n';
		}
		out.close();
	}
}

WaterBuilderForRestype::WaterBuilderForRestype(
	WaterBuildersList const & builders,
	utility::vector1< AtomWeights > const & atom_weights
) :
	n_waters_( 0 ),
	n_waters_for_atom_( builders.size() ),
	water_offset_for_atom_( builders.size(), 0 ),
	builders_( builders ),
	atom_weights_( atom_weights )
{
	Size natoms = builders.size();
	for ( Size ii = 1; ii <= natoms; ++ii ) {
		Size ii_nwaters = builders[ ii ].size();
		water_offset_for_atom_[ ii ] = n_waters_;
		n_waters_for_atom_[ ii ] = ii_nwaters;
		n_waters_ += ii_nwaters;
	}
}



LKBallDatabase::LKBallDatabase() = default; // Empty initialization of the list

LKBallDatabase::~LKBallDatabase() {
	// We need to de-register the destruction observers from each remaining ResidueType
	// Otherwise when they get destroyed they'll attempt to notify this (no longer existent) database
	// (Even though this is a Singleton object, there's destruction order issues at the end of the run.)
	for ( auto entry: water_builders_map_ ) {
		entry.first->detach_destruction_obs( &LKBallDatabase::restype_destruction_observer, this );
	}
	// We assume that for every entry the water_builders_map_ has, the atom_weights_map_ also has one
}

bool
LKBallDatabase::has( chemical::ResidueType const & rsd_type ) const {
	ResidueType const * const address( &rsd_type );
#if defined MULTI_THREADED
	utility::thread::ReadLockGuard lock( lkball_db_mutex_ );
#endif
	return water_builders_map_.count( address );
}

/// called the first time we encounter a given ResidueType
void
LKBallDatabase::initialize_residue_type( ResidueType const & rsd_type )
{
	using namespace conformation;
	using namespace chemical;

	bool const sidechain_only = !(basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_for_bb ]());

	ResidueType const * const address( &rsd_type );

#if defined MULTI_THREADED
	utility::thread::WriteLockGuard lock( lkball_db_mutex_ );
#endif

	if ( water_builders_map_.find( address ) == water_builders_map_.end() ) {
		TR.Trace << "initialize_residue_type: " << rsd_type.name() << std::endl;

		// Attach a destruction observer so we can clean up the database if the residue gets destroyed
		rsd_type.attach_destruction_obs( &LKBallDatabase::restype_destruction_observer, this );

		WaterBuildersList builders;
		setup_water_builders_for_residue_type( rsd_type, sidechain_only, builders );

		utility::vector1< AtomWeights > atom_wts;
		setup_atom_weights( rsd_type, builders, atom_wts );

		WaterBuilderForRestypeOP restype_builder = utility::pointer::make_shared< WaterBuilderForRestype >( builders, atom_wts );
		water_builders_map_[ address ] = restype_builder;
	}
}

void
LKBallDatabase::setup_atom_weights(
	chemical::ResidueType const & rsd_type,
	WaterBuildersList const & rsd_water_builders, // for sanity
	utility::vector1< AtomWeights > & atom_wts
)
{
	using utility::vector1;
	using ObjexxFCL::stripped;

	chemical::AtomTypeSet const & atom_set( rsd_type.atom_type_set() );
	std::string atom_wts_tag( "_RATIO23.0_DEFAULT" );
	if ( basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_wtd_tag ].user() ) {
		atom_wts_tag = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_wtd_tag ]();
	}

	utility::vector1< Real > lk_ball_wtd_prefactors( 6, 1.0 ); //
	if ( basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_wtd_prefactors ].user() ) {
		lk_ball_wtd_prefactors = basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_wtd_prefactors ]();
		if ( lk_ball_wtd_prefactors.size() != 6 ) {
			utility_exit_with_message("Option -lk_ball_wtd_prefactors should provide 6 numbers: <don-iso> <don-ball> <acc-iso> <acc-ball> <don+acc-iso> <don+acc-ball>");
		}
	}
	// wts of 1.0 for non-donor, non-acceptor atoms; they will presumably not have waters attached...
	lk_ball_wtd_prefactors.insert( lk_ball_wtd_prefactors.begin(), 1.0 );
	lk_ball_wtd_prefactors.insert( lk_ball_wtd_prefactors.begin(), 1.0 );

	Size const lkbi_atom_wt_index( atom_set.extra_parameter_index( "LK_BALL_ISO_ATOM_WEIGHT"+atom_wts_tag ) );
	Size const  lkb_atom_wt_index( atom_set.extra_parameter_index( "LK_BALL_ATOM_WEIGHT"+atom_wts_tag ) );

	atom_wts.clear(); atom_wts.resize( rsd_type.natoms() );
	for ( Size i=1; i<= rsd_type.natoms(); ++i ) {
		atom_wts[i][1] = 0;
		atom_wts[i][2] = 0;
	}

	runtime_assert( lk_ball_wtd_prefactors.size() == 8 );
	for ( Size i=1; i<= rsd_type.nheavyatoms(); ++i ) {
		bool const is_donor( rsd_type.atom_type(i).is_donor() ), is_acceptor( rsd_type.atom_type(i).is_acceptor() );
		Size const prefactor_offset( 4 * is_acceptor + 2 * is_donor );
		Real const lkbi_prefactor( lk_ball_wtd_prefactors[ prefactor_offset + 1 ] ), lkb_prefactor( lk_ball_wtd_prefactors[ prefactor_offset + 2 ] );

		atom_wts[i][1] = lkbi_prefactor * rsd_type.atom_type(i).extra_parameter( lkbi_atom_wt_index );
		atom_wts[i][2] =  lkb_prefactor * rsd_type.atom_type(i).extra_parameter(  lkb_atom_wt_index );
		if ( !rsd_water_builders[i].empty() ) {
			using ObjexxFCL::format::F;
			TR.Trace << "lk_ball_wtd atom_wts: " <<
				" iso_prefactor: " << F(9,3,lkbi_prefactor) <<
				" iso_wt: " << F(9,3,atom_wts[i][1]) <<
				" ball_prefactor: " << F(9,3,lkb_prefactor) <<
				" ball_wt: " << F(9,3,atom_wts[i][2]) <<
				' ' << rsd_type.atom_name(i) << ' ' << rsd_type.atom_type(i).name() << ' ' << rsd_type.name() << std::endl;
		}
	}
}

WaterBuilderForRestypeCOP
LKBallDatabase::get_water_builder_for_restype( chemical::ResidueType const & rsd_type ) const {

	ResidueType const * const address( &rsd_type );

	WaterBuildersForRestypeMap::const_iterator it;
	{
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard lock( lkball_db_mutex_ );
#endif
		it = ( water_builders_map_.find( address ) );
	}

	if ( it == water_builders_map_.end() ) {
		utility_exit_with_message("LKB_ResidueInfo::initialize has not been called");
	}

	return it->second;
}


void
LKBallDatabase::reset_arrays_danger_expert_only()
{
#if defined MULTI_THREADED
	utility::thread::WriteLockGuard lock( lkball_db_mutex_ );
#endif

	// Disconnect the destruction observers first
	for ( auto entry: water_builders_map_ ) {
		entry.first->detach_destruction_obs( &LKBallDatabase::restype_destruction_observer, this );
	}
	water_builders_map_.clear();
}

void
LKBallDatabase::restype_destruction_observer( core::chemical::RestypeDestructionEvent const & event ) {
#if defined MULTI_THREADED
	utility::thread::WriteLockGuard lock( lkball_db_mutex_ );
#endif

	// Remove the soon-to-be-destroyed object from the database caches
	// std::map::erase() is robust if the key is not in the map
	water_builders_map_.erase( event.restype );
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

LKB_ResidueInfo::LKB_ResidueInfo(
	conformation::Residue const & rsd,
	bool compute_derivs
)
{
	initialize( rsd.type() );
	build_waters( rsd, compute_derivs );
}

LKB_ResidueInfo::LKB_ResidueInfo( LKB_ResidueInfo const & /*src*/ ) = default;

LKB_ResidueInfo::LKB_ResidueInfo() = default;

LKB_ResidueInfo &
LKB_ResidueInfo::operator=( LKB_ResidueInfo const & ) = default;

void
LKB_ResidueInfo::build_waters( Residue const & rsd, bool compute_derivs )
{
	if ( !this->matches_residue_type( rsd.type() ) ) {
		utility_exit_with_message("LKB_ResidueInfo::build_waters: mismatch: "+rsd_type_->name()+" "+rsd.type().name() );
	}

	if ( !has_waters_ ) return;

	if ( compute_derivs ) {
		dwater_datom_ready_ = true;
		dwater_datom1_.resize( water_builders_->n_waters() );
		dwater_datom2_.resize( water_builders_->n_waters() );
		dwater_datom3_.resize( water_builders_->n_waters() );
	}

	for ( Size ii=1; ii <= rsd.nheavyatoms(); ++ii ) {
		Size ii_offset = water_builders_->water_offset_for_atom()[ ii ];
		WaterBuilders const & ii_water_builders( water_builders_->builders()[ ii ] );
		for ( Size jj=1, jj_end = ii_water_builders.size(); jj <= jj_end; ++jj ) {
			Size jj_water_ind = jj + ii_offset;
			waters_[ jj_water_ind ] = ii_water_builders[jj].build( rsd );
			if ( compute_derivs ) {
				ii_water_builders[ jj ].derivatives( rsd,
					dwater_datom1_[ jj_water_ind ],
					dwater_datom2_[ jj_water_ind ],
					dwater_datom3_[ jj_water_ind ]);
			}
		}
	}
}

WaterBuilders const &
LKB_ResidueInfo::get_water_builder( conformation::Residue const & rsd, Size heavyatom ) const
{
	if ( !matches_residue_type( rsd.type() ) ) {
		utility_exit_with_message("LKB_ResidueInfo::get_water_builder: mismatch: "+rsd_type_->name()+" "+rsd.type().name() );
	}

	if ( !has_waters_ ) {
		utility_exit_with_message("LKB_ResidueInfo::get_water_builder: no water builders!");
	}

	if ( heavyatom > water_builders_->builders().size() ) {
		utility_exit_with_message("LKB_ResidueInfo::get_water_builder called on unrecognized atom");
	}

	return water_builders_->builders()[ heavyatom ];
}

/// resize the waters_ array
/// set has_waters_
/// setup atom_weights_
///fpd setup dwater_datom*_
void
LKB_ResidueInfo::initialize( ResidueType const & rsd )
{
	rsd_type_ = rsd.get_self_ptr();
	dwater_datom_ready_ = false;

	LKBallDatabase & lk_db( * LKBallDatabase::get_instance() );
	if ( ! lk_db.has( rsd ) ) {
		lk_db.initialize_residue_type( rsd );
	}

	water_builders_ = LKBallDatabase::get_instance()->get_water_builder_for_restype( rsd );
	waters_.resize( water_builders_->n_waters(), Vector( 0.0 ) );
	has_waters_ = waters_.size() > 0;

}

basic::datacache::CacheableDataOP
LKB_ResidueInfo::clone() const
{
	return utility::pointer::make_shared< LKB_ResidueInfo >( *this );
}


LKB_ResiduesInfo::LKB_ResiduesInfo( LKB_ResiduesInfo const & src ):
	CacheableData(src)
{
	residues_info_.clear();
	for ( Size i=1; i<= src.size(); ++i ) {
		residues_info_.push_back( utility::pointer::make_shared< LKB_ResidueInfo >( src[i] ) );
	}
}

basic::datacache::CacheableDataOP
LKB_ResiduesInfo::clone() const
{
	return utility::pointer::make_shared< LKB_ResiduesInfo >( *this );
}

bool
LKB_ResidueInfo::matches_residue_type( chemical::ResidueType const & rsd_type ) const {
	return ( rsd_type_.get() == &(rsd_type) );
}

chemical::ResidueType const &
LKB_ResidueInfo::residue_type() const { return *rsd_type_; }


}
}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::lkball::LKB_ResidueInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( rsd_type_ ) );
	// EXEMPT water_builders_
	arc( CEREAL_NVP( waters_ ) ); // utility::vector1< WaterCoords >
	arc( CEREAL_NVP( dwater_datom_ready_ ) ); // bool
	arc( CEREAL_NVP( dwater_datom1_ ) );
	arc( CEREAL_NVP( dwater_datom2_ ) );
	arc( CEREAL_NVP( dwater_datom3_ ) );
	arc( CEREAL_NVP( has_waters_ ) ); // _Bool
}

/// @brief Manually generated deserialization method:
/// First deserialize the residue type, then query the LKBallDatabase for
/// the parameters for this residue type.
template< class Archive >
void
core::scoring::lkball::LKB_ResidueInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( rsd_type_ );

	LKBallDatabase & lk_db( * LKBallDatabase::get_instance() );
	if ( ! lk_db.has( *rsd_type_ ) ) {
		lk_db.initialize_residue_type( *rsd_type_ );
	}
	water_builders_ = LKBallDatabase::get_instance()->get_water_builder_for_restype( *rsd_type_ );

	arc( waters_ ); // utility::vector1<WaterCoords>
	arc( dwater_datom_ready_ );
	arc( dwater_datom1_ );
	arc( dwater_datom2_ );
	arc( dwater_datom3_ );
	arc( has_waters_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::lkball::LKB_ResidueInfo );
CEREAL_REGISTER_TYPE( core::scoring::lkball::LKB_ResidueInfo )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::lkball::LKB_ResiduesInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( residues_info_ ) ); // utility::vector1<LKB_ResidueInfoOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::lkball::LKB_ResiduesInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( residues_info_ ); // utility::vector1<LKB_ResidueInfoOP>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::lkball::LKB_ResiduesInfo );
CEREAL_REGISTER_TYPE( core::scoring::lkball::LKB_ResiduesInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_lkball_LK_BallInfo )
#endif // SERIALIZATION
