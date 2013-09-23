// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   core/scoring/methods/LK_BallEnergy.hh
/// @brief  LK Solvation using hemisphere culling class declaration
/// @author David Baker
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/LK_BallInfo.hh>

// // Package headers
// #include <core/scoring/methods/LK_BallEnergy.hh>
// #include <core/scoring/ScoringManager.hh>
// #include <core/scoring/NeighborList.hh>
// #include <core/scoring/EnergyGraph.hh>
// #include <core/scoring/etable/Etable.hh>
// #include <core/scoring/etable/count_pair/CountPairFunction.hh>
// #include <core/scoring/etable/count_pair/CountPairFactory.hh>
// #include <core/scoring/etable/count_pair/types.hh>

// // Project headers
#include <core/pose/Pose.hh>
// #include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
//#include <core/pack/rotamer_set/WaterPackingInfo.hh>
#include <basic/options/util.hh>
#include <basic/database/open.hh>


#include <core/io/pdb/pose_io.hh> // HACK
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
#include <numeric/xyz.functions.hh>

#include <utility/io/izstream.hh>
#include <utility/tools/make_vector1.hh>

#include <sstream>

// #include <utility/vector1.functions.hh> // HACK

//#ifdef WIN32
//#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
//#endif

namespace core {
namespace scoring {
namespace methods {

	//bool const sidechain_only_hack( false ); // make this configurable

/// @details Auto-generated virtual destructor
LKB_ResidueInfo::~LKB_ResidueInfo() {}


static basic::Tracer TR("core.scoring.methods.LK_BallInfo" );


/// LAZY
using namespace chemical;
using namespace pose;
using namespace conformation;

	static Real const optimal_water_distance( 2.65 ); /// note that this number is re-defined in hbonds.cc !!

/// Not doing backbone waters on protein or DNA
// inline
// bool
// residue_type_has_waters( ResidueType const & rsd_type ) {
// 	return ( !sidechain_only_hack || ( rsd_type.is_polar() || rsd_type.is_charged() ||
// 																		 ( rsd_type.is_aromatic() &&
// 																			 ( rsd_type.aa() == aa_trp || rsd_type.aa() == aa_tyr ) ) ) );
// }

LKB_ResidueInfo::WaterBuilderMap LKB_ResidueInfo::water_builder_map_;
LKB_ResidueInfo::AtomWeightsMap LKB_ResidueInfo::atom_weights_map_;
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
{}

Vector
WaterBuilder::build( conformation::Residue const & rsd ) const
{
	return kinematics::Stub( rsd.xyz( atom1_ ), rsd.xyz( atom2_ ), rsd.xyz( atom3_ ) ).local2global( xyz_local_ );
}

/// The next two functions were taken from protocols/water/rotamer_building_functions.cc with slight modifications
///
Vector
build_optimal_water_O_on_donor(
	Vector const & hxyz,
	Vector const & dxyz
)
{
	return ( dxyz + optimal_water_distance * ( hxyz - dxyz ).normalized() );
}


utility::vector1< Vector >
build_optimal_water_Os_on_acceptor(
	Size const acc_atm,
	conformation::Residue const & acc_rsd,
	Vector const & a_xyz,
	Vector b1_xyz, // local copy
	Vector b2_xyz,// local copy
	chemical::Hybridization const & hybrid
)
{
	using numeric::conversions::radians;
	using namespace chemical;

	Real const distance( optimal_water_distance ); // acceptor--O distance

	Real theta(0);
	utility::vector1< Real > phi_list;

	// detect special case of DNA phosphate hydration:
	std::string const & acc_atm_name( acc_rsd.atom_name( acc_atm ) );

	if ( acc_rsd.is_DNA() && acc_rsd.atom_is_backbone( acc_atm ) &&
			 ( acc_atm_name == " OP2" || acc_atm_name == " OP1" ) ) {
		// special case: hydration of the DNA phosphate group
		b1_xyz = acc_rsd.xyz( "P" );
		b2_xyz = ( ( acc_atm_name == " OP2" ) ? acc_rsd.xyz( "OP1" ) : acc_rsd.xyz( "OP2" ) );
		// these numbers are taken from "Hydration of the Phosphate Group in Double-Helical DNA",
		// Schneider, Patel, and Berman, Biophysical Journal Vol. 75 2422-2434
		theta = 180.0 - 125.0;
		phi_list.push_back(  45.0 );
		phi_list.push_back( 160.0 );
		phi_list.push_back( 280.0 );
	} else {

		switch( hybrid ) {
		case SP2_HYBRID:
			theta = 180.0 - 120.0;
			phi_list.push_back(   0.0 );
			phi_list.push_back( 180.0 );
			break;
		case SP3_HYBRID:
			theta = 180.0 - 109.0;
			phi_list.push_back( 120.0 );
			phi_list.push_back( 240.0 );
			break;
		case RING_HYBRID:
			b1_xyz = 0.5 * ( b1_xyz + b2_xyz );
			theta = 0.0;
			phi_list.push_back( 0.0 ); // doesnt matter
			break;
		default:
			TR.Error << "Bad hybridization type for acceptor " << hybrid << std::endl;
			utility_exit();
		}
	}

	utility::vector1< Vector > waters;
	kinematics::Stub stub( a_xyz, b1_xyz, b2_xyz );
	for ( Size i=1; i<= phi_list.size(); ++i ) {
		waters.push_back( stub.spherical( radians( phi_list[i] ), radians( theta ), distance ) );
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

	//bool const no_lk_ball_for_SP2( options::option[ options::OptionKeys::dna::specificity::no_lk_ball_for_SP2 ] );
	bool const no_SP3_acceptor_waters( false );
	//( basic::options::option[ basic::options::OptionKeys::dna::specificity::no_SP3_acceptor_waters ] );

	bool const dump_waters_pdb( false ); // for debugging purposes

	ResidueOP rsd( ResidueFactory::create_residue( rsd_type ) );

	rsd_water_builders.clear();
	rsd_water_builders.resize( rsd_type.natoms() ); // initialize to empty vector1s

	utility::vector1< Vector > all_waters; // for debugging

	// donors
	for ( chemical::AtomIndices::const_iterator
					hnum  = rsd->Hpos_polar().begin(),
					hnume = rsd->Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const hatm( *hnum );
		if ( sidechain_only && rsd_type.atom_is_backbone( hatm ) ) continue;

		Size const datm( rsd->atom_base( hatm ) );
		Size const datm_base( rsd->atom_base( datm ) );

		Vector const water( build_optimal_water_O_on_donor( rsd->xyz( hatm ), rsd->xyz( datm ) ) );
		rsd_water_builders[ datm ].push_back( WaterBuilder( water, *rsd, hatm, datm, datm_base ) );
		TR.Trace << "adding ideal water to rsd_type= " << rsd_type.name() << " anchored to atom " <<
			rsd_type.atom_name( datm ) << " datm_base= " << rsd_type.atom_name( datm_base ) << std::endl;
		all_waters.push_back( water );
	}

	// acceptors
	for ( chemical::AtomIndices::const_iterator
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

LKB_ResidueInfo::LKB_ResidueInfo(
																 pose::Pose const &,
																 conformation::Residue const & rsd
																 )
{
	initialize( rsd.type() ); // sets atom wts
	build_waters( rsd );
}

LKB_ResidueInfo::LKB_ResidueInfo()
{
}



/// called the first time we encounter a given ResidueType
void
LKB_ResidueInfo::initialize_residue_type( ResidueType const & rsd_type ) const
{
	using namespace conformation;
	using namespace chemical;

	//assert( residue_type_has_waters( rsd_type ) );

	TR.Trace << "initialize_residue_type: " << rsd_type.name() << std::endl;

	bool const sidechain_only( true ); // consider making this option-configurable
	// bool const sidechain_only( ! basic::options::option[ basic::options::OptionKeys::dna::specificity::lk_ball_for_bb ] );

	ResidueTypeCOP const address( &rsd_type );
	assert( ! water_builder_map_.count( address ) );

	water_builder_map_[ address ]; // create entry in map
	utility::vector1< WaterBuilders > & rsd_water_builders( water_builder_map_.find( address )->second );

	setup_water_builders_for_residue_type( rsd_type, sidechain_only, rsd_water_builders );

	atom_weights_map_[ address ];
	utility::vector1< utility::vector1< Real > > & atom_wts( atom_weights_map_.find( address )->second );

	setup_atom_weights( rsd_type, rsd_water_builders, atom_wts );
}

void
LKB_ResidueInfo::setup_atom_weights(
																		chemical::ResidueType const & rsd_type,
																		utility::vector1< WaterBuilders > const & rsd_water_builders, // for sanity
																		utility::vector1< utility::vector1< Real > > & atom_wts
																		) const
{
	using utility::vector1;
	using ObjexxFCL::stripped;

	if ( true ) { // new way

		chemical::AtomTypeSet const & atom_set( rsd_type.atom_type_set() );
		Size const atom_wt_index( atom_set.extra_parameter_index( "LK_BALL_ATOM_WEIGHT" ) );

		atom_wts.clear(); atom_wts.resize( rsd_type.natoms() );
		for ( Size i=1; i<= rsd_type.natoms(); ++i ) atom_wts[i] = utility::tools::make_vector1( Real(0.0), Real(1.0) );

		for ( Size i=1; i<= rsd_type.nheavyatoms(); ++i ) {
			atom_wts[i][1] = 0.0; // the iso weight
			atom_wts[i][2] = rsd_type.atom_type(i).extra_parameter( atom_wt_index );
			if ( !rsd_water_builders[i].empty() ) {
				TR.Trace << "lk_ball_wtd atom_wts: " <<
					ObjexxFCL::format::F(9,3,atom_wts[i][1]) << ObjexxFCL::format::F(9,3,atom_wts[i][2]) << ' ' <<
					rsd_type.atom_name(i) << ' ' << rsd_type.name() << std::endl;
			}
		}

	} else { // the OLD WAY
		static bool init( false );
		static vector1< vector1< Real > > residue_sidechain_weights;
		if ( !init ) {
			init = true;
			residue_sidechain_weights.clear();
			residue_sidechain_weights.resize( chemical::num_aa_types );

			for ( Size i=1; i<= chemical::num_aa_types; ++i ) residue_sidechain_weights[i].clear();

			utility::io::izstream data;
			basic::database::open( data, "scoring/lk_ball_wtd_weights.txt" );
			std::string line;
			while ( getline( data, line ) ) {
				std::istringstream l( line );
				std::string tmp, tag1,tag2,tag3,tag4,name1s;
				Real polwt, donwt, accwt;
				l >> tag1 >> name1s >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >>
					tag2 >> polwt >>
					tag3 >> donwt >>
					tag4 >> accwt;
				if ( l.fail() || tag1 != "FITSOL" || tag2 != "polwt:" || tag3 != "donwt:" || tag4 != "accwt:" ) {
					utility_exit_with_message("parse error: "+line);
				}
				chemical::AA const aa( aa_from_oneletter_code( name1s[0] ) );
				residue_sidechain_weights[ aa ] = utility::tools::make_vector1( polwt, donwt, accwt );
			}
		} // initialize aa wts from file

		// initialize weights to 0.0 !!!!!!!!!!!!!!!!
		atom_wts.clear(); atom_wts.resize( rsd_type.natoms() );
		for ( Size i=1; i<= rsd_type.natoms(); ++i ) atom_wts[i] = utility::tools::make_vector1( Real(1.0), Real(0.0) );

		vector1< Real > sc_wts( residue_sidechain_weights[rsd_type.aa()] );

		if ( sc_wts.empty() ) {
			sc_wts = utility::tools::make_vector1( 1.0, 0.0, 0.0 );
			for ( Size i=1; i<= rsd_type.nheavyatoms(); ++i ) {
			if ( !rsd_water_builders[i].empty() && !rsd_type.atom_is_backbone(i) ) {
					TR.Trace << "[WARNING] Using default lk_ball_wtd atom_weights of iso=1, aniso=0 for " <<
						rsd_type.atom_name(i) << ' ' << rsd_type.atom_type(i).name() << ' ' << rsd_type.name() << std::endl;
				}
			}
		}

		Real const polwt( sc_wts[1] ), donwt( sc_wts[2] ), accwt( sc_wts[3] );

		for ( Size i=1; i<= rsd_type.nheavyatoms(); ++i ) {
			Size const nwaters( rsd_water_builders[i].size() );
			if ( !nwaters ) continue;
			std::string const atomname( stripped( rsd_type.atom_name(i) ) ),
				atomtypename( stripped( rsd_type.atom_type(i).name()));

			if ( rsd_type.atom_is_backbone(i) ) {
				if        ( rsd_type.is_lower_terminus() && atomname == "N" && rsd_type.aa() != aa_pro ) {
					runtime_assert( nwaters == 3 );
					runtime_assert( atomtypename == "Nlys" );
					atom_wts[i][1] = 1.053; atom_wts[i][2] = 1.955; // taken from Lysine sidechain
				} else if ( rsd_type.is_lower_terminus() && atomname == "N" && rsd_type.aa() == aa_pro ) {
					runtime_assert( nwaters == 2 );
					//runtime_assert( atomtypename == "Nlys" );
					atom_wts[i][1] = 1.0; atom_wts[i][2] = 0.0; // these were just made up, not sure if I understand Nterm-pro
				} else if ( rsd_type.is_upper_terminus() && ( atomname == "O" || atomname == "OXT" ) ) {
					runtime_assert( nwaters == 2 );
					runtime_assert( atomtypename == "OOC" );
					atom_wts[i][1] = 1.888; atom_wts[i][2] = 1.584; // taken from Asp sidechain
				} else if ( atomname == "N" && rsd_type.aa() != aa_pro ) { // these next 2 are fitted by apbs/ajob1c.bb.py:
					// python ajob1c.bb.py tmp.ros_pairs.hiq_flipv2.cp3_md9.log
					// using APBS backbone solvation data computed in apbs/job12
					// need to double-check that we get the same results from job14, which had correct Nterm H locations
					runtime_assert( nwaters == 1 );
					runtime_assert( atomtypename == "Nbb" );
					atom_wts[i][1] = 0.26; atom_wts[i][2] = 0.53;
				} else if ( atomname == "O") {
					runtime_assert( nwaters == 2 );
					runtime_assert( atomtypename == "OCbb" );
					atom_wts[i][1] = 1.48; atom_wts[i][2] = 0.0;
				} else {
					runtime_assert( nwaters == 0 );
					continue;
				}
			} else {
				chemical::AtomType const & at( rsd_type.atom_type(i) );
				atom_wts[i][1] = polwt; // iso wt
				atom_wts[i][2] = 0.0; // aniso wt
				if ( at.is_acceptor() ) atom_wts[i][2] = accwt; // aniso wt
				if ( at.is_donor() ) atom_wts[i][2] = donwt; // aniso wt, DONOR GOES LAST ==> SP3 uses donor weight
			}
			TR.Trace << "lk_ball_wtd atom_weight: " << rsd_type.name1() << ' ' <<
				ObjexxFCL::format::A(6,atomname) <<
				ObjexxFCL::format::A(6,atomtypename) <<
				" nwaters: " << nwaters <<
				" iso_wt: " << ObjexxFCL::format::F(9,3,atom_wts[i][1] ) <<
				" aniso_wt: " << ObjexxFCL::format::F(9,3,atom_wts[i][2] ) <<
				' ' << rsd_type.name() << std::endl;
		}
	}
}


void
LKB_ResidueInfo::build_waters( Residue const & rsd )
{
	if ( !has_waters_ ) return;

	// waters_ array has already been dimensioned properly
	ResidueTypeCOP const address( &( rsd.type() ) );

	WaterBuilderMap::const_iterator it( water_builder_map_.find( address ) );
	if ( it == water_builder_map_.end() ) {
		utility_exit_with_message("LKB_ResidueInfo::initialize has not been called");
	}

	for ( Size i=1; i<= rsd.nheavyatoms(); ++i ) {
		WaterBuilders const & water_builders( it->second[ i ] );
		for ( Size j=1, j_end = water_builders.size(); j<= j_end; ++j ) {
			waters_[i][j] = water_builders[j].build( rsd );
		}
	}
}

/// resize the waters_ array
/// set has_waters_
/// setup atom_weights_

void
LKB_ResidueInfo::initialize( ResidueType const & rsd )
{
	waters_.clear(); waters_.resize( rsd.nheavyatoms() );
	has_waters_ = false;

	ResidueTypeCOP const address( &rsd );

	WaterBuilderMap::const_iterator it( water_builder_map_.find( address ) );
	if ( it == water_builder_map_.end() ) {
		initialize_residue_type( rsd );
		it = water_builder_map_.find( address );
	}

	for ( Size i=1; i<= rsd.nheavyatoms(); ++i ) {
		WaterBuilders const & water_builders( it->second[ i ] );
		if ( water_builders.empty() ) {
			waters_[i].clear();
		} else {
			has_waters_ = true;
			waters_[i].resize( water_builders.size() );
		}
	}

	// get the atom_weights_
	atom_weights_ = atom_weights_map_.find( address )->second;
}

LKB_ResidueInfo::LKB_ResidueInfo( LKB_ResidueInfo const & src ):
	CacheableData(),
	waters_( src.waters_ ),
	atom_weights_( src.atom_weights_ ), // added 5/20/13
	has_waters_( src.has_waters_ )
{}

basic::datacache::CacheableDataOP
LKB_ResidueInfo::clone() const
{
	return new LKB_ResidueInfo( *this );
}


LKB_ResiduesInfo::LKB_ResiduesInfo( LKB_ResiduesInfo const & src ):
	CacheableData()
{
	residues_info_.clear();
	for ( Size i=1; i<= src.size(); ++i ) {
		residues_info_.push_back( new LKB_ResidueInfo( src[i] ) );
	}
}

basic::datacache::CacheableDataOP
LKB_ResiduesInfo::clone() const
{
	return new LKB_ResiduesInfo( *this );
}



}
}
}
