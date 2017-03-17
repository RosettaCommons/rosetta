// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/CartesianBondedEnergy.cc
/// @brief  Harmonic bondangle/bondlength/torsion constraints
/// @author Frank DiMaio
/// @author modified by Vikram K. Mulligan to allow D-amino acid minimization and cyclic geometry.

// Unit headers
#include <core/scoring/methods/CartBondedParameters.hh>
#include <core/scoring/methods/CartesianBondedEnergy.hh>
#include <core/scoring/methods/CartesianBondedEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AA.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/mm/MMBondAngleLibrary.hh>
#include <core/scoring/mm/MMBondLengthLibrary.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/scoring/PolymerBondedEnergyContainer.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/deriv/angle_deriv.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/dihedral_deriv.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/angle.functions.hh>
#include <core/scoring/DerivVectorPair.hh>

// options
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/optE.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

// C++ headers
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/string.hpp>

#include <utility/vector1.hh>

#include <core/pose/PDBInfo.hh>

// cutoff for reporting outliers
#define CUTOFF 1

namespace boost {
namespace tuples {

std::size_t hash_value(atm_name_quad const& e) {
	std::size_t seed = 0;
	boost::hash_combine(seed, e.get<0>());
	boost::hash_combine(seed, e.get<1>());
	boost::hash_combine(seed, e.get<2>());
	boost::hash_combine(seed, e.get<3>());
	boost::hash_combine(seed, e.get<4>());
	return seed;
}
std::size_t hash_value(atm_name_triple const& e) {
	std::size_t seed = 0;
	boost::hash_combine(seed, e.get<0>());
	boost::hash_combine(seed, e.get<1>());
	boost::hash_combine(seed, e.get<2>());
	boost::hash_combine(seed, e.get<3>());
	return seed;
}
std::size_t hash_value(atm_name_pair const& e) {
	std::size_t seed = 0;
	boost::hash_combine(seed, e.get<0>());
	boost::hash_combine(seed, e.get<1>());
	boost::hash_combine(seed, e.get<2>());
	return seed;
}
std::size_t hash_value(atm_name_single const& e) {
	std::size_t seed = 0;
	boost::hash_combine(seed, e.get<0>());
	boost::hash_combine(seed, e.get<1>());
	return seed;
}

bool operator==(atm_name_quad const& a,atm_name_quad const& b) {
	return a.get<0>() == b.get<0>() && a.get<1>() == b.get<1>() &&
		a.get<2>() == b.get<2>() && a.get<3>() == b.get<3>() &&
		a.get<4>() == b.get<4>();
}
bool operator==(atm_name_triple const& a,atm_name_triple const& b) {
	return a.get<0>() == b.get<0>() && a.get<1>() == b.get<1>() &&
		a.get<2>() == b.get<2>() && a.get<3>() == b.get<3>();
}
bool operator==(atm_name_pair const& a,atm_name_pair const& b) {
	return a.get<0>() == b.get<0>() && a.get<1>() == b.get<1>() &&
		a.get<2>() == b.get<2>();
}
bool operator==(atm_name_single const& a,atm_name_single const& b) {
	return a.get<0>() == b.get<0>() && a.get<1>() == b.get<1>();
}

}
}


namespace core {
namespace scoring {
namespace methods {


static THREAD_LOCAL basic::Tracer TR( "core.scoring.CartesianBondedEnergy" );
static THREAD_LOCAL basic::Tracer GEOMETRIES( "core.scoring.CartesianBondedEnergy.GEOMETRIES" );

// default spring constants
//  --> only if not specified in config files
static const Real K_LENGTH=300.0;
static const Real K_ANGLE=80.0;
static const Real K_TORSION=20.0;
static const Real K_TORSION_PROTON=10.0;  // proton chi
static const Real K_TORSION_IMPROPER=40.0;

IdealParametersDatabaseOP CartesianBondedEnergy::db_=NULL;

#if defined MULTI_THREADED
utility::thread::ReadWriteMutex CartesianBondedEnergy::params_db_mutex_;
utility::thread::ReadWriteMutex CartesianBondedEnergy::restype_db_mutex_;
#endif


//////////////////////
/// EnergyMethod Creator
methods::EnergyMethodOP
CartesianBondedEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new CartesianBondedEnergy( options ) );
}

ScoreTypes
CartesianBondedEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( cart_bonded );
	sts.push_back( cart_bonded_angle );
	sts.push_back( cart_bonded_length );
	sts.push_back( cart_bonded_torsion );
	sts.push_back( cart_bonded_ring );
	sts.push_back( cart_bonded_proper );
	sts.push_back( cart_bonded_improper );
	return sts;
}


////////////////////////
/// helper function
std::string get_restag( core::chemical::ResidueType const & restype ) {
	using namespace core::chemical;

	if ( core::chemical::is_canonical_D_aa(restype.aa()) ) {
		std::string rsdname = restype.name();
		if ( rsdname.substr(0, rsdname.find(chemical::PATCH_LINKER)) == "DHIS_D" ) rsdname="HIS_D"; //If this is a DHIS_D, return HIS_D.
		else rsdname=core::chemical::name_from_aa( core::chemical::get_L_equivalent( restype.aa() ) ); //Otherwise, for D-amino acids, return the L-equivalent.
		return rsdname;
	} else if ( !restype.is_protein() && !restype.is_NA() ) {
		return restype.name3();
	} else {
		std::string rsdname = restype.name();
		rsdname = rsdname.substr( 0, rsdname.find(chemical::PATCH_LINKER) );
		return rsdname;
	}
}

ResidueCartBondedParameters::ResidueCartBondedParameters() :
	bb_N_index_( 0 ),
	bb_CA_index_( 0 ),
	bb_C_index_( 0 ),
	bb_O_index_( 0 ),
	bb_H_index_( 0 ),
	pro_CD_index_( 0 )
{}

ResidueCartBondedParameters::~ResidueCartBondedParameters() {}

void ResidueCartBondedParameters::add_length_parameter(  Size2 atom_inds, CartBondedParametersCOP params )
{
	length_params_.push_back( std::make_pair( atom_inds, params ));
}

void ResidueCartBondedParameters::add_angle_parameter(   Size3 atom_inds, CartBondedParametersCOP params )
{
	angle_params_.push_back( std::make_pair( atom_inds, params ));
}

void ResidueCartBondedParameters::add_torsion_parameter( Size4 atom_inds, CartBondedParametersCOP params )
{
	torsion_params_.push_back( std::make_pair( atom_inds, params ));
}

void ResidueCartBondedParameters::add_improper_parameter( Size4 atom_inds, CartBondedParametersCOP params )
{
	improper_params_.push_back( std::make_pair( atom_inds, params ));
}

void ResidueCartBondedParameters::add_bbdep_length_parameter(  Size2 atom_inds, CartBondedParametersCOP params )
{
	bbdep_length_params_.push_back( std::make_pair( atom_inds, params ));
}

void ResidueCartBondedParameters::add_bbdep_angle_parameter(   Size3 atom_inds, CartBondedParametersCOP params )
{
	bbdep_angle_params_.push_back( std::make_pair( atom_inds, params ));
}

void ResidueCartBondedParameters::add_lower_connect_angle_params( Size3 atom_inds, CartBondedParametersCOP params )
{
	lower_connect_angle_params_.push_back( std::make_pair( atom_inds, params ));
}

void ResidueCartBondedParameters::add_upper_connect_angle_params( Size3 atom_inds, CartBondedParametersCOP params )
{
	upper_connect_angle_params_.push_back( std::make_pair( atom_inds, params ));
}

void ResidueCartBondedParameters::bb_N_index( Size index )  { bb_N_index_  = index; }
void ResidueCartBondedParameters::bb_CA_index( Size index ) { bb_CA_index_ = index; }
void ResidueCartBondedParameters::bb_C_index( Size index )  { bb_C_index_  = index; }
void ResidueCartBondedParameters::bb_O_index( Size index )  { bb_O_index_  = index; }
void ResidueCartBondedParameters::bb_H_index( Size index )  { bb_H_index_  = index; }
void ResidueCartBondedParameters::pro_CD_index( Size index )  { pro_CD_index_  = index; }

void ResidueCartBondedParameters::ca_cprev_n_h_interres_improper_params(
	CartBondedParametersCOP params
)
{
	ca_cprev_n_h_interres_improper_params_ = params;
}

void ResidueCartBondedParameters::oprev_cprev_n_h_interres_improper_params(
	CartBondedParametersCOP params
)
{
	oprev_cprev_n_h_interres_improper_params_ = params;
}


void ResidueCartBondedParameters::ca_nnext_c_o_interres_improper_params(
	CartBondedParametersCOP params
)
{
	ca_nnext_c_o_interres_improper_params_ = params;
}

void ResidueCartBondedParameters::pro_cd_cprev_n_ca_interres_improper_params(
	CartBondedParametersCOP params
)
{
	pro_cd_cprev_n_ca_interres_improper_params_ = params;
}

void ResidueCartBondedParameters::cprev_n_bond_length_params(
	CartBondedParametersCOP params
)
{
	cprev_n_bond_length_params_ = params;
}


////////////////////////
// constructors
IdealParametersDatabase::IdealParametersDatabase(Real k_len_in, Real k_ang_in, Real k_tors_in, Real k_tors_prot_in, Real k_tors_improper_in) {
	// figure out the defaults given priorities:
	//   1: energy methods options (constructor args)
	//   2: command line/options system
	//   3: fallback defaults
	Real k_len=K_LENGTH, k_ang=K_ANGLE, k_tors=K_TORSION, k_tors_prot=K_TORSION_PROTON, k_tors_improper=K_TORSION_IMPROPER;
	utility::vector1<core::Real> params;
	if ( basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user() ) {
		params = basic::options::option[ basic::options::OptionKeys::score::bonded_params ]();
	}
	if ( k_len_in >= 0 ) {
		k_len = k_len_in;
	} else if ( basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user() ) {
		if ( params.size() >= 1 ) k_len = params[1];
	}

	if ( k_ang_in >= 0 ) {
		k_ang = k_ang_in;
	} else if ( basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user() ) {
		if ( params.size() >= 2 ) k_ang = params[2];
	}

	if ( basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user() ) {
		if ( params.size() >= 3 ) k_tors = params[3];
		if ( params.size() >= 4 ) k_tors_prot = params[4];
		if ( params.size() >= 5 ) k_tors_improper = params[5];
	}
	if ( k_tors_in >= 0 ) {
		k_tors = k_tors_in;
	}
	if ( k_tors_prot_in >= 0 ) {
		k_tors_prot = k_tors_prot_in;
	}
	if ( k_tors_improper_in >= 0 ) {
		k_tors_improper = k_tors_improper_in;
	}

	// now that we resolved the defaults, initialize
	init(k_len, k_ang, k_tors, k_tors_prot, k_tors_improper);
}

// read bb-dep and bb-indep parameter files from the rosetta database
void
IdealParametersDatabase::init(
	Real k_len_in,
	Real k_ang_in,
	Real k_tors_in,
	Real k_tors_prot_in,
	Real k_tors_improper_in
) {
	TR << "Initializing IdealParametersDatabase with default Ks="
		<< k_len_in << " , " << k_ang_in << " , " << k_tors_in << " , " << k_tors_prot_in << " , " << k_tors_improper_in << std::endl;
	k_length_ = k_len_in;
	k_angle_ = k_ang_in;
	k_torsion_ = k_tors_in;
	k_torsion_proton_ = k_tors_prot_in;
	k_torsion_improper_ = k_tors_improper_in;

	// read bb-independent parameters
	std::string libpath = basic::options::option[ basic::options::OptionKeys::score::bonded_params_dir ]();
	bbdep_bond_params_ = basic::options::option[ basic::options::OptionKeys::corrections::score::bbdep_bond_params ]();
	bbdep_bond_devs_ = basic::options::option[ basic::options::OptionKeys::corrections::score::bbdep_bond_devs ]();

	bool const symm_gly( basic::options::option[ basic::options::OptionKeys::score::symmetric_gly_tables ]() );

	read_length_database( libpath+"/default-lengths.txt", symm_gly );
	read_angle_database( libpath+"/default-angles.txt", symm_gly );
	read_torsion_database( libpath+"/default-torsions.txt", symm_gly );
	read_improper_database( libpath+"/default-improper.txt", symm_gly ); // used to be called "torsion" until July 2016

	if ( !bbdep_bond_params_ ) return;

	// read bb-dep parameters
	//NOTE this must be read after the bbindep params are read
	//   as deviations for low count rama bins are smoothed using the bb-indep deviations
	read_bbdep_table( libpath+"/bbdep-all-but-gly-pro-xpro-ile-val.txt", bondlengths_bbdep_def_, bondangles_bbdep_def_, "ALA", false );
	read_bbdep_table( libpath+"/bbdep-graphdata-gly.txt", bondlengths_bbdep_gly_, bondangles_bbdep_gly_, "GLY", symm_gly );
	read_bbdep_table( libpath+"/bbdep-graphdata-ile-val.txt", bondlengths_bbdep_valile_, bondangles_bbdep_valile_, "VAL", false );
	read_bbdep_table( libpath+"/bbdep-graphdata-pro.txt", bondlengths_bbdep_pro_, bondangles_bbdep_pro_, "PRO", false );
	read_bbdep_table( libpath+"/bbdep-graphdata-xpro.txt", bondlengths_bbdep_prepro_, bondangles_bbdep_prepro_, "ALA", false );
}


void
IdealParametersDatabase::read_length_database(
	std::string const & infile,
	bool const symmetrize_gly
) {
	std::string line;
	std::string name3, atom1, atom2;
	Real mu_d, K_d;

	utility::io::izstream instream;
	atm_name_pair tuple;
	basic::database::open( instream, infile);
	while ( instream ) {
		getline( instream, line );
		if ( line[0] == '#' ) continue;

		std::istringstream linestream(line);

		linestream >> name3 >> atom1 >> atom2 >> mu_d >> K_d;
		tuple = boost::make_tuple( name3, atom1, atom2 );
		CartBondedParametersOP params_i( new BBIndepCartBondedParameters(mu_d, K_d) );
		bondlengths_indep_.insert( std::make_pair( tuple, params_i ) );
		tuple = boost::make_tuple( name3, atom2, atom1 );
		params_i = CartBondedParametersOP( new BBIndepCartBondedParameters(mu_d, K_d) );
		bondlengths_indep_.insert( std::make_pair( tuple, params_i ) );
	}
	TR << "Read " << bondlengths_indep_.size() << " bb-independent lengths." << std::endl;

	if ( symmetrize_gly ) {
		TR << "Symmetrizing glycine bond lengths table." << std::endl;
		atm_name_pair const tuple1( boost::make_tuple( "GLY", "CA", "1HA" ) );
		atm_name_pair const tuple2( boost::make_tuple( "GLY", "CA", "2HA" ) );
		CartBondedParametersOP param1, param2;
		boost::unordered_map<atm_name_pair,CartBondedParametersOP>::iterator it = bondlengths_indep_.find( tuple1 );
		boost::unordered_map<atm_name_pair,CartBondedParametersOP>::iterator it2 = bondlengths_indep_.find( tuple2 );
		debug_assert ( it != bondlengths_indep_.end() );
		debug_assert ( it2 != bondlengths_indep_.end() );
		param1 = it->second;
		param2 = it2->second;
		core::Real const avgmu( (param1->mu(0.0, 0.0) + param2->mu(0.0, 0.0)) / 2.0 );
		core::Real const avgK( (param1->K(0.0, 0.0) + param2->K(0.0, 0.0)) / 2.0 );
		bondlengths_indep_.erase(tuple1); bondlengths_indep_.erase(tuple2);
		bondlengths_indep_.insert( std::make_pair(tuple1, CartBondedParametersOP( new BBIndepCartBondedParameters(avgmu, avgK) ) ) );
		bondlengths_indep_.insert( std::make_pair(tuple2, CartBondedParametersOP( new BBIndepCartBondedParameters(avgmu, avgK) ) ) );
	}
}

void
IdealParametersDatabase::read_angle_database(
	std::string const & infile,
	bool const symmetrize_gly
) {
	std::string line;
	std::string name3, atom1, atom2, atom3;
	Real mu_d, K_d;

	utility::io::izstream instream;
	atm_name_triple tuple;
	basic::database::open( instream, infile);
	while ( instream ) {
		getline( instream, line );
		if ( line[0] == '#' ) continue;

		std::istringstream linestream(line);

		linestream >> name3 >> atom1 >> atom2 >> atom3 >> mu_d >> K_d;
		tuple = boost::make_tuple( name3, atom1, atom2, atom3 );
		CartBondedParametersOP params_i( new BBIndepCartBondedParameters(mu_d, K_d) );
		bondangles_indep_.insert( std::make_pair( tuple, params_i) );
		tuple = boost::make_tuple( name3, atom3, atom2, atom1 );
		bondangles_indep_.insert( std::make_pair( tuple, params_i) );
	}
	TR << "Read " << bondangles_indep_.size() << " bb-independent angles." << std::endl;

	if ( symmetrize_gly ) {
		TR << "Symmetrizing glycine bond angles table." << std::endl;
		atm_name_triple const tuple1( boost::make_tuple( "GLY", "N", "CA", "1HA" ) );
		atm_name_triple const tuple2( boost::make_tuple( "GLY", "N", "CA", "2HA" ) );
		core::Real const avgmu( (bondangles_indep_[ tuple1 ]->mu(0.0, 0.0) + bondangles_indep_[ tuple2 ]->mu(0.0, 0.0)) / 2.0 );
		core::Real const avgK( (bondangles_indep_[ tuple1 ]->K(0.0, 0.0) + bondangles_indep_[ tuple2 ]->K(0.0, 0.0)) / 2.0 );
		bondangles_indep_[ tuple1 ] = CartBondedParametersOP( new BBIndepCartBondedParameters(avgmu, avgK) );
		bondangles_indep_[ tuple2 ] = CartBondedParametersOP( new BBIndepCartBondedParameters(avgmu, avgK) );
	}
}

void
IdealParametersDatabase::read_torsion_database(
	std::string const & infile,
	bool const //symmetrize_gly
) {
	std::string line;
	std::string name3, atom1, atom2, atom3, atom4;
	Real mu_d, K_d;
	Size period;

	utility::io::izstream instream;
	atm_name_quad tuple;
	basic::database::open( instream, infile);
	while ( instream ) {
		getline( instream, line );
		if ( line[0] == '#' ) continue;

		std::istringstream linestream(line);

		linestream >> name3 >> atom1 >> atom2 >> atom3 >> atom4 >> mu_d >> K_d >> period;
		tuple = boost::make_tuple( name3, atom1, atom2, atom3, atom4 );
		CartBondedParametersOP params_i( new BBIndepCartBondedParameters(mu_d, K_d, period) );
		torsions_indep_.insert( std::make_pair( tuple, params_i) );
	}
	TR << "Read " << torsions_indep_.size() << " bb-independent torsions." << std::endl;
}

void
IdealParametersDatabase::read_improper_database(
	std::string const & infile,
	bool const symmetrize_gly
) {
	std::string line;
	std::string name3, atom1, atom2, atom3, atom4;
	Real mu_d, K_d;
	Size period;

	utility::io::izstream instream;
	atm_name_quad tuple;
	basic::database::open( instream, infile);
	while ( instream ) {
		getline( instream, line );
		if ( line[0] == '#' ) continue;

		std::istringstream linestream(line);

		linestream >> name3 >> atom1 >> atom2 >> atom3 >> atom4 >> mu_d >> K_d >> period;
		tuple = boost::make_tuple( name3, atom1, atom2, atom3, atom4 );
		CartBondedParametersOP params_i( new BBIndepCartBondedParameters(mu_d, K_d, period) );
		impropers_indep_.insert( std::make_pair( tuple, params_i) );
	}
	TR << "Read " << impropers_indep_.size() << " bb-independent improper tors." << std::endl;

	//hpark  extra impropers from the command line
	if ( basic::options::option[ basic::options::OptionKeys::score::extra_improper_file ].user() ) {
		std::string extra_file( basic::options::option[ basic::options::OptionKeys::score::extra_improper_file ]().c_str() );
		atm_name_quad tuple;
		std::string line;
		Size const size_before = impropers_indep_.size();

		utility::io::izstream instream( extra_file );
		if ( !instream.good() ) utility_exit_with_message( "Unable to open file: " + extra_file + '\n' );

		while ( instream.getline( line ) ) {
			if ( line[0] == '#' ) continue;
			std::istringstream l( line );
			l >> name3 >> atom1 >> atom2 >> atom3 >> atom4 >> mu_d >> K_d >> period;
			tuple = boost::make_tuple( name3, atom1, atom2, atom3, atom4 );
			CartBondedParametersOP params_i( new BBIndepCartBondedParameters(mu_d, K_d, period) );
			impropers_indep_.insert( std::make_pair( tuple, params_i) );
		}
		TR << "Read " << impropers_indep_.size() - size_before << " extra improper tors.";
		TR << extra_file << std::endl;
	}

	if ( symmetrize_gly ) {
		TR << "Symmetrizing glycine improper dihedrals table." << std::endl;
		//We'll ADD parameters for glycine for this.
		atm_name_quad const tuple1( boost::make_tuple( "GLY", "CA", "C", "N", "H" ) );
		atm_name_quad const tuple2( boost::make_tuple( "GLY", "O", "C", "N", "H" ) );
		atm_name_quad const tuple3( boost::make_tuple( "*", "CA", "C", "N", "H" ) );
		atm_name_quad const tuple4( boost::make_tuple( "*", "O", "C", "N", "H" ) );
		core::Real const mu( numeric::constants::d::pi );
		core::Real const K( 40 );
		impropers_indep_[ tuple1 ] = CartBondedParametersOP( new BBIndepCartBondedParameters(mu, K) );
		impropers_indep_[ tuple2 ] = CartBondedParametersOP( new BBIndepCartBondedParameters(mu, K) );
		impropers_indep_[ tuple3 ] = CartBondedParametersOP( new BBIndepCartBondedParameters(mu, K) );
		impropers_indep_[ tuple4 ] = CartBondedParametersOP( new BBIndepCartBondedParameters(mu, K) );
	}
}

/// @brief Read bb-dependent tables
/// @details Smooth using bbdep data corresponding to residue 'resbase'
void
IdealParametersDatabase::read_bbdep_table(
	std::string const &filename,
	boost::unordered_map< atm_name_single, CartBondedParametersOP > &bondlengths,
	boost::unordered_map< atm_name_pair, CartBondedParametersOP > &bondangles,
	std::string const &resbase,
	bool const symmetrize_table
) {
	using numeric::constants::d::pi;
	using namespace ObjexxFCL;

	//std::cout << "READ_BBDEP_TABLE() CALLED." << std::endl; //DELETE ME

	Real DEV_SCALE = 0.1; // scaling factor between bbindep and bbdep spring constants
	core::Size M = 10; // m estimate to smooth low count bins

	utility::io::izstream instream;
	basic::database::open( instream, filename );

	//PhiStart PhiStop PsiStart PsiStop Observations PhiAvg(i) PhiDev(i) PsiAvg(i) PsiDev(i) OmeAvg(i) OmeDev(i) C(-1)-NAvg(i) C(-1)-NDev(i) N-CAAvg(i)
	// N-CADev(i) CA-CBAvg(i) CA-CBDev(i) CA-CAvg(i) CA-CDev(i) C-OAvg(i) C-ODev(i) C(-1)-N-CAAvg(i) C(-1)-N-CADev(i) N-CA-CBAvg(i) N-CA-CBDev(i)
	// N-CA-CAvg(i) N-CA-CDev(i) CB-CA-CAvg(i) CB-CA-CDev(i) CA-C-OAvg(i) CA-C-ODev(i) CA-C-N(+1)Avg(i) CA-C-N(+1)Dev(i) O-C-N(+1)Avg(i) O-C-N(+1)Dev(i)
	// Chi1Avg(i) Chi1Dev(i) Chi2Avg(i) Chi2Dev(i) Chi3Avg(i) Chi3Dev(i) Chi4Avg(i) Chi4Dev(i) ZetaAvg(i) ZetaDev(i)
	core::Real phiL,phiH,psiL,psiH,phiAvg,phiDev,psiAvg,psiDev,OmegaAvg,OmegaDev;
	core::Real CNavg,CNdev, NCAavg,NCAdev, CACBavg,CACBdev, CACavg,CACdev, COavg,COdev;
	core::Real CNCAavg,CNCAdev, NCACBavg,NCACBdev, NCACavg,NCACdev, CBCACavg,CBCACdev;
	core::Real CACOavg,CACOdev, CACNavg,CACNdev, OCNavg,OCNdev;
	core::Real c1avg,c1dev,c2avg,c2dev,c3avg,c3dev,c4avg,c4dev,Zavg,Zdev;
	core::Size nobs;

	// the bond angles
	// their bb-indep angles
	CartBondedParametersOP temp;
	temp = bondlengths_indep_[boost::make_tuple(resbase,"C", "N")];
	if ( !temp ) temp = bondlengths_indep_[boost::make_tuple("*","C", "N")];
	if ( !temp ) temp = CartBondedParametersOP( new BBIndepCartBondedParameters(0,0) );
	core::Real CN_indep_avg=temp->mu(0,0),CN_indep_K=temp->K(0,0);

	temp = bondlengths_indep_[boost::make_tuple(resbase,"N", "CA")];
	if ( !temp ) temp = bondlengths_indep_[boost::make_tuple("*","N", "CA")];
	if ( !temp ) temp = CartBondedParametersOP( new BBIndepCartBondedParameters(0,0) );
	core::Real NCA_indep_avg=temp->mu(0,0),NCA_indep_K=temp->K(0,0);

	temp = bondlengths_indep_[boost::make_tuple(resbase,"CA", "CB")];
	if ( !temp ) temp = bondlengths_indep_[boost::make_tuple("*","CA", "CB")];
	if ( !temp ) temp = CartBondedParametersOP( new BBIndepCartBondedParameters(0,0) );
	core::Real CACB_indep_avg=temp->mu(0,0),CACB_indep_K=temp->K(0,0);

	temp = bondlengths_indep_[boost::make_tuple(resbase,"CA", "C")];
	if ( !temp ) temp = bondlengths_indep_[boost::make_tuple("*","CA", "C")];
	if ( !temp ) temp = CartBondedParametersOP( new BBIndepCartBondedParameters(0,0) );
	core::Real CAC_indep_avg=temp->mu(0,0),CAC_indep_K=temp->K(0,0);

	temp = bondlengths_indep_[boost::make_tuple(resbase,"C", "O")];
	if ( !temp ) temp = bondlengths_indep_[boost::make_tuple("*","C", "O")];
	if ( !temp ) temp = CartBondedParametersOP( new BBIndepCartBondedParameters(0,0) );
	core::Real CO_indep_avg=temp->mu(0,0),CO_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"C", "N", "CA")];
	if ( !temp ) temp = bondangles_indep_[boost::make_tuple("*","C", "N", "CA")];
	if ( !temp ) temp = CartBondedParametersOP( new BBIndepCartBondedParameters(0,0) );
	core::Real CNCA_indep_avg=temp->mu(0,0),CNCA_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"N", "CA", "CB")];
	if ( !temp ) temp = bondangles_indep_[boost::make_tuple("*","N", "CA", "CB")];
	if ( !temp ) temp = CartBondedParametersOP( new BBIndepCartBondedParameters(0,0) );
	core::Real NCACB_indep_avg=temp->mu(0,0),NCACB_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"N", "CA", "C")];
	if ( !temp ) temp = bondangles_indep_[boost::make_tuple("*","N", "CA", "C")];
	if ( !temp ) temp = CartBondedParametersOP( new BBIndepCartBondedParameters(0,0) );
	core::Real NCAC_indep_avg=temp->mu(0,0),NCAC_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"CB", "CA", "C")];
	if ( !temp ) temp = bondangles_indep_[boost::make_tuple("*","CB", "CA", "C")];
	if ( !temp ) temp = CartBondedParametersOP( new BBIndepCartBondedParameters(0,0) );
	core::Real CBCAC_indep_avg=temp->mu(0,0),CBCAC_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"CA", "C", "O")];
	if ( !temp ) temp = bondangles_indep_[boost::make_tuple("*","CA", "C", "O")];
	if ( !temp ) temp = CartBondedParametersOP( new BBIndepCartBondedParameters(0,0) );
	core::Real CACO_indep_avg=temp->mu(0,0),CACO_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"CA", "C", "N")];
	if ( !temp ) temp = bondangles_indep_[boost::make_tuple("*","CA", "C", "N")];
	if ( !temp ) temp = CartBondedParametersOP( new BBIndepCartBondedParameters(0,0) );
	core::Real CACN_indep_avg=temp->mu(0,0),CACN_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"O", "C", "N")];
	if ( !temp ) temp = bondangles_indep_[boost::make_tuple("*","O", "C", "N")];
	if ( !temp ) temp = CartBondedParametersOP( new BBIndepCartBondedParameters(0,0) );
	core::Real OCN_indep_avg=temp->mu(0,0),OCN_indep_K=temp->K(0,0);

	// all the relevant tables
	FArray2D<core::Real> CNavg_tbl(36,36,0.0),CNdev_tbl(36,36,0.0);
	FArray2D<core::Real> NCAavg_tbl(36,36,0.0),NCAdev_tbl(36,36,0.0);
	FArray2D<core::Real> CACBavg_tbl(36,36,0.0),CACBdev_tbl(36,36,0.0);
	FArray2D<core::Real> CACavg_tbl(36,36,0.0),CACdev_tbl(36,36,0.0);
	FArray2D<core::Real> COavg_tbl(36,36,0.0),COdev_tbl(36,36,0.0);
	FArray2D<core::Real> CNCAavg_tbl(36,36,0.0),CNCAdev_tbl(36,36,0.0);
	FArray2D<core::Real> NCACBavg_tbl(36,36,0.0),NCACBdev_tbl(36,36,0.0);
	FArray2D<core::Real> NCACavg_tbl(36,36,0.0),NCACdev_tbl(36,36,0.0);
	FArray2D<core::Real> CBCACavg_tbl(36,36,0.0),CBCACdev_tbl(36,36,0.0);
	FArray2D<core::Real> CACOavg_tbl(36,36,0.0),CACOdev_tbl(36,36,0.0);
	FArray2D<core::Real> CACNavg_tbl(36,36,0.0),CACNdev_tbl(36,36,0.0);
	FArray2D<core::Real> OCNavg_tbl(36,36,0.0),OCNdev_tbl(36,36,0.0);
	FArray2D<core::Size> Ns(DynamicIndexRange(36),DynamicIndexRange(36),(core::Size)0);

	std::string fileline;
	while ( getline(instream, fileline) ) {
		if ( fileline[0] == '#' ) continue;
		std::stringstream(fileline) >> phiL >> phiH >> psiL >> psiH >> nobs >> phiAvg >> phiDev >> psiAvg >> psiDev >> OmegaAvg >> OmegaDev >>  CNavg >> CNdev
			>> NCAavg >> NCAdev >>  CACBavg >> CACBdev >>  CACavg >> CACdev >>  COavg >> COdev
			>> CNCAavg >> CNCAdev >>  NCACBavg >> NCACBdev >>  NCACavg >> NCACdev >>  CBCACavg >> CBCACdev
			>> CACOavg >> CACOdev >>  CACNavg >> CACNdev >>  OCNavg >> OCNdev
			>> c1avg >> c1dev >> c2avg >> c2dev >> c3avg >> c3dev >> c4avg >> c4dev >> Zavg >> Zdev;

		Size phibin = (Size)std::floor(phiL/10.0 + 0.5); //Ranges from 0 to 35 for phibins starting from 0.0 to 35.0, respectively.
		Size psibin = (Size)std::floor(psiL/10.0 + 0.5);

		CNavg_tbl(phibin+1,psibin+1) = CNavg; // Plus 1 because FArray2D is 1-based, not 0-based.
		NCAavg_tbl(phibin+1,psibin+1) = NCAavg;
		CACBavg_tbl(phibin+1,psibin+1) = CACBavg;
		CACavg_tbl(phibin+1,psibin+1) = CACavg;
		COavg_tbl(phibin+1,psibin+1) = COavg;
		CNCAavg_tbl(phibin+1,psibin+1) = CNCAavg * pi/180;
		NCACBavg_tbl(phibin+1,psibin+1) = NCACBavg * pi/180;
		NCACavg_tbl(phibin+1,psibin+1) = NCACavg * pi/180;
		CBCACavg_tbl(phibin+1,psibin+1) = CBCACavg * pi/180;
		CACOavg_tbl(phibin+1,psibin+1) = CACOavg * pi/180;
		CACNavg_tbl(phibin+1,psibin+1) = CACNavg * pi/180;
		OCNavg_tbl(phibin+1,psibin+1) = OCNavg * pi/180;
		Ns(phibin+1,psibin+1) = nobs;

		if ( !bbdep_bond_devs_ || nobs<=3 || CN_indep_K==0 ) CNdev_tbl(phibin+1,psibin+1) = CN_indep_K;
		else CNdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CN_indep_K) + (nobs-1)*( (2*(CNdev)*(CNdev))/DEV_SCALE) )/ (M+nobs-1) );

		if ( !bbdep_bond_devs_ || nobs<=3 || NCA_indep_K==0 ) NCAdev_tbl(phibin+1,psibin+1) = NCA_indep_K;
		else NCAdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/NCA_indep_K) + (nobs-1)*( (2*(NCAdev)*(NCAdev))/DEV_SCALE) )/ (M+nobs-1) );

		if ( !bbdep_bond_devs_ || nobs<=3 || CACB_indep_K==0 ) CACBdev_tbl(phibin+1,psibin+1) = CACB_indep_K;
		else CACBdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CACB_indep_K) + (nobs-1)*( (2*(CACBdev)*(CACBdev))/DEV_SCALE) )/ (M+nobs-1) );

		if ( !bbdep_bond_devs_ || nobs<=3 || CAC_indep_K==0 ) CACdev_tbl(phibin+1,psibin+1) = CAC_indep_K;
		else CACdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CAC_indep_K) + (nobs-1)*( (2*(CACdev)*(CACdev))/DEV_SCALE) )/ (M+nobs-1) );

		if ( !bbdep_bond_devs_ || nobs<=3 || CO_indep_K==0 ) COdev_tbl(phibin+1,psibin+1) = CO_indep_K;
		else COdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CO_indep_K) + (nobs-1)*( (2*(COdev)*(COdev))/DEV_SCALE) )/ (M+nobs-1) );

		if ( !bbdep_bond_devs_ || nobs<=3 || CNCA_indep_K==0 ) CNCAdev_tbl(phibin+1,psibin+1) = CNCA_indep_K;
		else CNCAdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CNCA_indep_K) + (nobs-1)*( (2*(CNCAdev*pi/180)*(CNCAdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );

		if ( !bbdep_bond_devs_ || nobs<=3 || NCACB_indep_K==0 ) NCACBdev_tbl(phibin+1,psibin+1) = NCACB_indep_K;
		else NCACBdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/NCACB_indep_K) + (nobs-1)*( (2*(NCACBdev*pi/180)*(NCACBdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );

		if ( !bbdep_bond_devs_ || nobs<=3 || NCAC_indep_K==0 ) NCACdev_tbl(phibin+1,psibin+1) = NCAC_indep_K;
		else NCACdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/NCAC_indep_K) + (nobs-1)*( (2*(NCACdev*pi/180)*(NCACdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );

		if ( !bbdep_bond_devs_ || nobs<=3 || CBCAC_indep_K==0 ) CBCACdev_tbl(phibin+1,psibin+1) = CBCAC_indep_K;
		else CBCACdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CBCAC_indep_K) + (nobs-1)*( (2*(CBCACdev*pi/180)*(CBCACdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );

		if ( !bbdep_bond_devs_ || nobs<=3 || CACO_indep_K==0 ) CACOdev_tbl(phibin+1,psibin+1) = CACO_indep_K;
		else CACOdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CACO_indep_K) + (nobs-1)*( (2*(CACOdev*pi/180)*(CACOdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );

		if ( !bbdep_bond_devs_ || nobs<=3 || CACN_indep_K==0 ) CACNdev_tbl(phibin+1,psibin+1) = CACN_indep_K;
		else CACNdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CACN_indep_K) + (nobs-1)*( (2*(CACNdev*pi/180)*(CACNdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );

		if ( !bbdep_bond_devs_ || nobs<=3 || OCN_indep_K==0 ) OCNdev_tbl(phibin+1,psibin+1) = OCN_indep_K;
		else OCNdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/OCN_indep_K) + (nobs-1)*( (2*(OCNdev*pi/180)*(OCNdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );
	}

	// fill in missing values
	for ( int i=1; i<=36; ++i ) {
		for ( int j=1; j<=36; ++j ) {
			if ( Ns(i,j) != 0 ) continue;

			CNavg_tbl(i,j)    = CN_indep_avg;    CNdev_tbl(i,j)    = CN_indep_K;
			NCAavg_tbl(i,j)   = NCA_indep_avg;   NCAdev_tbl(i,j)   = NCA_indep_K;
			CACBavg_tbl(i,j)  = CACB_indep_avg;  CACBdev_tbl(i,j)  = CACB_indep_K;
			CACavg_tbl(i,j)   = CAC_indep_avg;   CACdev_tbl(i,j)   = CAC_indep_K;
			COavg_tbl(i,j)    = CO_indep_avg;    COdev_tbl(i,j)    = CO_indep_K;
			CNCAavg_tbl(i,j)  = CNCA_indep_avg;  CNCAdev_tbl(i,j)  = CNCA_indep_K;
			NCACBavg_tbl(i,j) = NCACB_indep_avg; NCACBdev_tbl(i,j) = NCACB_indep_K;
			NCACavg_tbl(i,j)  = NCAC_indep_avg;  NCACdev_tbl(i,j)  = NCAC_indep_K;
			CBCACavg_tbl(i,j) = CBCAC_indep_avg; CBCACdev_tbl(i,j) = CBCAC_indep_K;
			CACOavg_tbl(i,j)  = CACO_indep_avg;  CACOdev_tbl(i,j)  = CACO_indep_K;
			CACNavg_tbl(i,j)  = CACN_indep_avg;  CACNdev_tbl(i,j)  = CACN_indep_K;
			OCNavg_tbl(i,j)   = OCN_indep_avg;   OCNdev_tbl(i,j)   = OCN_indep_K;
		}
	}

	// Symmetrize the tables, if we're doing that (added by VKM, 18 Sept. 2016):
	if ( symmetrize_table ) {
		symmetrize_tables( CNavg_tbl );
		symmetrize_tables( NCAavg_tbl );
		symmetrize_tables( CACBavg_tbl );
		symmetrize_tables( CACavg_tbl );
		symmetrize_tables( COavg_tbl );
		symmetrize_tables( CNCAavg_tbl );
		symmetrize_tables( NCACBavg_tbl );
		symmetrize_tables( NCACavg_tbl );
		symmetrize_tables( CBCACavg_tbl );
		symmetrize_tables( CACOavg_tbl );
		symmetrize_tables( CACNavg_tbl );
		symmetrize_tables( OCNavg_tbl );
		symmetrize_tables( CNdev_tbl );
		symmetrize_tables( NCAdev_tbl );
		symmetrize_tables( CACBdev_tbl );
		symmetrize_tables( CACdev_tbl );
		symmetrize_tables( COdev_tbl );
		symmetrize_tables( CNCAdev_tbl );
		symmetrize_tables( NCACBdev_tbl );
		symmetrize_tables( NCACdev_tbl );
		symmetrize_tables( CBCACdev_tbl );
		symmetrize_tables( CACOdev_tbl );
		symmetrize_tables( CACNdev_tbl );
		symmetrize_tables( OCNdev_tbl );
	}

	// make the database entries
	bondlengths[boost::make_tuple("N","C")] = bondlengths[boost::make_tuple("C","N")] = CartBondedParametersOP( new BBDepCartBondedParameters(CNavg_tbl, CNdev_tbl ,"C-N") );
	bondlengths[boost::make_tuple("CA","N")] = bondlengths[boost::make_tuple("N","CA")] = CartBondedParametersOP( new BBDepCartBondedParameters(NCAavg_tbl, NCAdev_tbl,"CA-N") );
	bondlengths[boost::make_tuple("CB","CA")] = bondlengths[boost::make_tuple("CA","CB")] = CartBondedParametersOP( new BBDepCartBondedParameters(CACBavg_tbl, CACBdev_tbl,"CB-CA") );
	bondlengths[boost::make_tuple("C","CA")] = bondlengths[boost::make_tuple("CA","C")] = CartBondedParametersOP( new BBDepCartBondedParameters(CACavg_tbl, CACdev_tbl,"C-CA") );
	bondlengths[boost::make_tuple("O","C")] = bondlengths[boost::make_tuple("C","O")] = CartBondedParametersOP( new BBDepCartBondedParameters(COavg_tbl, COdev_tbl,"C-O") );
	bondangles[boost::make_tuple("CA","N","C")] = bondangles[boost::make_tuple("C","N","CA")] = CartBondedParametersOP( new BBDepCartBondedParameters(CNCAavg_tbl, CNCAdev_tbl, "C-N-CA") );
	bondangles[boost::make_tuple("CB","CA","N")] = bondangles[boost::make_tuple("N","CA","CB")] = CartBondedParametersOP( new BBDepCartBondedParameters(NCACBavg_tbl, NCACBdev_tbl, "N-CA-CB") );
	bondangles[boost::make_tuple("C","CA","N")] = bondangles[boost::make_tuple("N","CA","C")] = CartBondedParametersOP( new BBDepCartBondedParameters(NCACavg_tbl, NCACdev_tbl, "N-CA-C") );
	bondangles[boost::make_tuple("C","CA","CB")] = bondangles[boost::make_tuple("CB","CA","C")] = CartBondedParametersOP( new BBDepCartBondedParameters(CBCACavg_tbl, CBCACdev_tbl, "CB-CA-C") );
	bondangles[boost::make_tuple("O","C","CA")] = bondangles[boost::make_tuple("CA","C","O")] = CartBondedParametersOP( new BBDepCartBondedParameters(CACOavg_tbl, CACOdev_tbl, "CA-C-O") );
	bondangles[boost::make_tuple("N","C","CA")] = bondangles[boost::make_tuple("CA","C","N")] = CartBondedParametersOP( new BBDepCartBondedParameters(CACNavg_tbl, CACNdev_tbl, "CA-C-N") );
	bondangles[boost::make_tuple("N","C","O")] = bondangles[boost::make_tuple("O","C","N")] = CartBondedParametersOP( new BBDepCartBondedParameters(OCNavg_tbl, OCNdev_tbl, "O-C-N") );
}

/// @brief Symmetrize the glycine backbone-dependent table.
/// @details Only called if the score::symmetric_gly_tables option is used.  Intended for design
/// with glyceine in a mixed D/L context (in which there should be no preference for a left-handed
/// conformation over a right).
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
IdealParametersDatabase::symmetrize_tables(
	ObjexxFCL::FArray2D<core::Real> &table
) {

	//std::cout << "SYMMETRIZE_TABLES CALLED." << std::endl; //DELETE ME

	if ( TR.Debug.visible() ) {
		TR.Debug << "\nSymmetrizing gly table.  Pre-symm:";
		for ( core::Size i( table.l1() ), imax( table.u1() ); i<=imax; ++i ) {
			for ( core::Size j( table.l2() ), jmax( table.u2() ); j<=jmax; ++j ) {
				TR.Debug << table(i, j);
				if ( j<jmax ) TR.Debug << "\t";
			}
			TR.Debug << "\n";
		}
		TR.Debug << std::endl;
	}

	ObjexxFCL::FArray2D< core::Real > newtable( table ); //Copy the old table.  Slightly inefficient, but easier to make sure that I'm doing what I think I'm doing.

	for ( core::Size i( table.l1() ), imax( table.u1() ); i<=imax; ++i ) {
		for ( core::Size j( table.l2() ), jmax( table.u2() ); j<=jmax; ++j ) {
			core::Size const iother( imax - i + 1 );
			core::Size const jother( jmax - j + 1 );
			newtable( i, j ) = (table(i, j) + table(iother, jother)) / 2.0; //Slightly inefficient, since we do the same calculation twice for equivalent symmetric points, but not a big deal.
		}
	}

	table = newtable; //Copy the new table back

	if ( TR.Debug.visible() ) {
		TR.Debug << "\nSymmetrizing gly table.  Post-symm:";
		for ( core::Size i( table.l1() ), imax( table.u1() ); i<=imax; ++i ) {
			for ( core::Size j( table.l2() ), jmax( table.u2() ); j<=jmax; ++j ) {
				TR.Debug << table(i, j);
				if ( j<jmax ) TR.Debug << "\t";
			}
			TR.Debug << "\n";
		}
		TR.Debug << std::endl;
	}

}

void
IdealParametersDatabase::lookup_bondangle_buildideal(
	core::chemical::ResidueType const & restype,
	int atm1,
	int atm2,
	int atm3,
	Real & Ktheta,
	Real & theta0
) {
	Ktheta = k_angle_;

	// Create mini-conformation for idealized residue
	conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );

	// get angle
	numeric::xyzVector<core::Real> x,y,z; // atom coords
	if ( atm1 > 0 ) {
		x = newres->atom( atm1 ).xyz();
	} else {
		x = newres->residue_connection(-atm1).icoor().build( restype );
	}

	y = newres->atom( atm2 ).xyz();

	if ( atm3 > 0 ) {
		z = newres->atom( atm3 ).xyz();
	} else {
		z = newres->residue_connection(-atm3).icoor().build( restype );
	}

	theta0 = numeric::angle_radians ( x,y,z );

	// EXCEPTIONS
	// ignore proline (L- or D-proline) C-ND which is handled by pro_close
	if ( restype.aa() == core::chemical::aa_pro || restype.aa() == core::chemical::aa_dpr ) {
		bool hasCD = (atm1>0 && restype.atom_name(atm1)==" CD ")
			|| (atm2>0 && restype.atom_name(atm2)==" CD ")
			|| (atm3>0 && restype.atom_name(atm3)==" CD ");
		bool hasN = (atm1>0 && restype.atom_name(atm1)==" N  ")
			|| (atm2>0 && restype.atom_name(atm2)==" N  ")
			|| (atm3>0 && restype.atom_name(atm3)==" N  ");
		if ( hasCD && hasN ) {
			Ktheta = theta0 = 0.0;
		}
	}

	// fpd ignore centroid angle in ALA (L- or D-alanine) and GLY
	if ( (restype.aa() == core::chemical::aa_ala || restype.aa() == core::chemical::aa_dal || restype.aa() == core::chemical::aa_gly) &&
			( (atm1>0 && restype.atom_name(atm1) == " CEN") || (atm3>0 && restype.atom_name(atm3) == " CEN") ) ) {
		Ktheta = theta0 = 0.0;
	}

	//fpd  cutpoint variants
	if ( restype.has_variant_type(chemical::CUTPOINT_UPPER) &&
			( (atm1>0 && restype.atom_name(atm1) == "OVU1") || (atm3>0 && restype.atom_name(atm3) == "OVU1") ) ) {
		Ktheta = theta0 = 0.0;
	}
	if ( restype.has_variant_type(chemical::CUTPOINT_LOWER) &&
			( (atm1>0 && restype.atom_name(atm1) == "OVL1") || (atm3>0 && restype.atom_name(atm3) == "OVL1") ) ) {
		Ktheta = theta0 = 0.0;
	}
}


void
IdealParametersDatabase::lookup_bondlength_buildideal(
	core::chemical::ResidueType const & restype,
	int atm1,
	int atm2,
	Real &Kd,
	Real &d0
)
{
	Kd = k_length_;

	// Create mini-conformation for idealized residue
	core::conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );

	// get length
	numeric::xyzVector<core::Real> x,y,z; // atom coords
	if ( atm1 > 0 ) {
		x = newres->atom( atm1 ).xyz();
	} else {
		x = newres->residue_connection(-atm1).icoor().build( restype );
	}

	if ( atm2 > 0 ) {
		y = newres->atom( atm2 ).xyz();
	} else {
		y = newres->residue_connection(-atm2).icoor().build( restype );
	}

	d0 = (x-y).length();

	//  ignore proline (L- or D-proline) C-ND which is handled by pro_close
	if ( (restype.aa() == core::chemical::aa_pro || restype.aa() == core::chemical::aa_dpr) &&
			(atm1>0 && atm2>0) &&
			( ( restype.atom_name(atm1) == " CD " && restype.atom_name(atm2) == " N  ") ||
			( restype.atom_name(atm2) == " CD " && restype.atom_name(atm1) == " N  ") ) ) {
		Kd = d0 = 0.0;
	}
}


CartBondedParametersCOP
IdealParametersDatabase::lookup_improper(
	core::chemical::ResidueType const & restype,
	std::string const & atm1_name,
	std::string const & atm2_name,
	std::string const & atm3_name,
	std::string const & atm4_name
)
{
	using namespace core::chemical;

	// use 'annotated sequence' to id this restype
	std::string restag = get_restag( restype );

	// there is no bb-dep torsion database
	// lookup in bb-indep table
	// this can probably be made way faster
	atm_name_quad tuple1( restag, atm1_name,atm2_name,atm3_name,atm4_name );
	boost::unordered_map<atm_name_quad,CartBondedParametersOP>::iterator b_it = impropers_indep_.find( tuple1 );
	if ( b_it != impropers_indep_.end() ) {
		return b_it->second;
	}

	atm_name_quad tuple2( "*", atm1_name,atm2_name,atm3_name,atm4_name );
	b_it = impropers_indep_.find( tuple2 );
	if ( b_it != impropers_indep_.end() ) {
		return b_it->second;
	}

	// if we don't find this improper torsion in the table, it's unconstrained
	return NULL;
}

/// Angle Database
///   build from ideal if not found
CartBondedParametersCOP
IdealParametersDatabase::lookup_angle(
	core::chemical::ResidueType const & restype,
	bool pre_proline,
	std::string const & atm1_name,
	std::string const & atm2_name,
	std::string const & atm3_name,
	int atm1idx,
	int atm2idx,
	int atm3idx
) {
	using namespace core::chemical;

	std::string restag = get_restag( restype );
	atm_name_triple tuple( restag, atm1_name,atm2_name,atm3_name );

	// 1 lookup in bb-dep table
	// figure out what table to look in
	boost::unordered_map< atm_name_pair, CartBondedParametersOP > *angle_table;
	if ( restype.aa() == core::chemical::aa_gly ) {
		angle_table = &bondangles_bbdep_gly_;
	} else if ( restype.aa() == core::chemical::aa_pro || restype.aa() == core::chemical::aa_dpr ) { //L- or D-proline
		angle_table = &bondangles_bbdep_pro_;
	} else if ( pre_proline ) {
		angle_table = &bondangles_bbdep_prepro_;
	} else if ( restype.aa() == core::chemical::aa_ile || restype.aa() == core::chemical::aa_val || restype.aa() == core::chemical::aa_dil || restype.aa() == core::chemical::aa_dva ) { //L- or D-isoleucine or valine
		angle_table = &bondangles_bbdep_valile_;
	} else {
		angle_table = &bondangles_bbdep_def_;
	}
	atm_name_pair alt_tuple( atm1_name,atm2_name,atm3_name );
	boost::unordered_map<atm_name_pair,CartBondedParametersOP>::iterator a_it = angle_table->find( alt_tuple );
	if ( a_it != angle_table-> end() ) {
		return a_it->second;
	}

	// 2 lookup in bb-indep table
	boost::unordered_map<atm_name_triple,CartBondedParametersOP>::iterator b_it = bondangles_indep_.find( tuple );
	if ( b_it != bondangles_indep_.end() ) {
		return b_it->second;
	}

	atm_name_triple tuple_alt( "*", atm1_name,atm2_name,atm3_name );
	b_it = bondangles_indep_.find( tuple_alt );
	if ( b_it != bondangles_indep_.end() ) {
		return b_it->second;
	}

	// 3 build from ideal, add to bb-indep table
	Real Ktheta, theta0;
	lookup_bondangle_buildideal( restype, atm1idx, atm2idx, atm3idx, Ktheta, theta0 );
	bondangles_indep_[ tuple ] = CartBondedParametersOP( new BBIndepCartBondedParameters( theta0, Ktheta ) );
	TR << "Adding undefined angle "
		<< restag << ": " << tuple.get<1>() << "," << tuple.get<2>() << "," << tuple.get<3>() << " to DB with"
		<< " theta0 = " << theta0 << " , Ktheta = " << Ktheta << std::endl;

	return bondangles_indep_[ tuple ];
}

/// BondLength Database
///   build from ideal if not found
CartBondedParametersCOP
IdealParametersDatabase::lookup_length(
	chemical::ResidueType const & restype,
	bool pre_proline,
	std::string const & atm1_name,
	std::string const & atm2_name,
	int atm1idx,
	int atm2idx
) {
	// use 'annotated sequence' to id this restype
	std::string restag = get_restag( restype );
	atm_name_pair tuple( restag, atm1_name,atm2_name);

	// 1 lookup in bb-dep table
	// figure out what table to look in
	boost::unordered_map< atm_name_single, CartBondedParametersOP > *length_table;
	if ( restype.aa() == core::chemical::aa_gly ) {
		length_table = &bondlengths_bbdep_gly_;
	} else if ( restype.aa() == core::chemical::aa_pro || restype.aa() == core::chemical::aa_dpr ) { //D- or L-proline
		length_table = &bondlengths_bbdep_pro_;
	} else if ( restype.aa() == core::chemical::aa_ile || restype.aa() == core::chemical::aa_val || restype.aa() == core::chemical::aa_dil || restype.aa() == core::chemical::aa_dva ) { //D- or L-isoleucine or valine
		length_table = &bondlengths_bbdep_valile_;
	} else if ( pre_proline ) {
		length_table = &bondlengths_bbdep_prepro_;
	} else {
		length_table = &bondlengths_bbdep_def_;
	}

	atm_name_single alt_tuple( atm1_name,atm2_name );
	boost::unordered_map<atm_name_single,CartBondedParametersOP>::iterator a_it = length_table->find( alt_tuple );
	if ( a_it != length_table->end() ) {
		return a_it->second;
	}

	// 2 lookup in bb-indep table
	boost::unordered_map<atm_name_pair,CartBondedParametersOP>::iterator b_it = bondlengths_indep_.find( tuple );
	if ( b_it != bondlengths_indep_.end() ) {
		return b_it->second;
	}

	atm_name_pair tuple_alt( "*", atm1_name,atm2_name );
	b_it = bondlengths_indep_.find( tuple_alt );
	if ( b_it != bondlengths_indep_.end() ) {
		return b_it->second;
	}

	// 3 build from ideal, add to bb-indep table
	Real Kd, d0;
	lookup_bondlength_buildideal( restype, atm1idx, atm2idx, Kd, d0 );
	bondlengths_indep_[ tuple ] = CartBondedParametersOP( new BBIndepCartBondedParameters( d0, Kd ) );
	TR << "Adding undefined length "
		<< restag << ": " << tuple.get<1>() << "," << tuple.get<2>() << "," << " to DB with"
		<< " d0 = " << d0 << " , Kd = " << Kd << std::endl;

	return bondlengths_indep_[ tuple ];
}


// old-style interface to database
// somewhat slow as it will always build the ideal residue on the fly
void
IdealParametersDatabase::lookup_torsion_legacy(
	core::chemical::ResidueType const & restype,
	int atm1,
	int atm2,
	int atm3,
	int atm4,
	Real & Kphi,
	Real & phi0,
	Real &phi_step
) {
	using namespace core::chemical;

	Kphi=k_torsion_; // default
	phi_step=0;

	if ( !restype.is_protein() ) {
		Kphi = phi0 = 0;
		return;
	}

	// use 'annotated sequence' to id this restype
	//std::string restag = get_restag( restype );

	// Create mini-conformation for idealized residue
	core::conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );
	core::pose::Pose newpose;
	newpose.append_residue_by_bond(*newres);

	// figure out if we need to constrain this torsion
	bool need_to_constrain=true, proton_chi=false;

	// backbone -- we never need to constrain intrares backbones (for now)
	if ( atm2 <= (int)newres->last_backbone_atom() && atm3 <= (int)newres->last_backbone_atom() ) {
		need_to_constrain = false;
	}

	// chi
	for ( Size j=1, j_end = restype.nchi(); j<= j_end; ++j ) {
		id::AtomID n1,n2,n3,n4;
		newpose.conformation().get_torsion_angle_atom_ids( id::TorsionID(1,id::CHI,j) ,n1,n2,n3,n4);
		if ( ((int)n2.atomno() == atm2 && (int)n3.atomno() == atm3) || ((int)n2.atomno() == atm3 && (int)n3.atomno() == atm2) ) {
			if ( restype.is_proton_chi( j ) ) {
				proton_chi = true;
			} else {
				need_to_constrain = false;
			}
		}
	}

	// check if this corresponds to a DOF_ID
	if ( !need_to_constrain ) {
		Kphi=0.0;
		phi0=0.0;
	} else {
		// get angle
		numeric::xyzVector<core::Real> x,y,z, w; // atom coords
		x = newres->atom( atm1 ).xyz();
		y = newres->atom( atm2 ).xyz();
		z = newres->atom( atm3 ).xyz();
		w = newres->atom( atm4 ).xyz();

		phi0 = numeric::dihedral_radians ( x,y,z,w );

		//////// exceptions!
		// fpd ignore proline C-ND which is handled by pro_close
		if ( restype.aa() == core::chemical::aa_pro || restype.aa() == core::chemical::aa_dpr ) { //D- or L-proline
			bool hasCD = (restype.atom_name(atm1)==" CD ") || (restype.atom_name(atm2)==" CD ")
				|| (restype.atom_name(atm3)==" CD ") || (restype.atom_name(atm4)==" CD ");
			bool hasN = (restype.atom_name(atm1)==" N  ") || (restype.atom_name(atm2)==" N  ")
				|| (restype.atom_name(atm3)==" N  ") || (restype.atom_name(atm4)==" N  ");
			if ( hasCD && hasN ) {
				Kphi=0.0; phi0=0.0;
			}
		}
		if ( restype.has_variant_type(chemical::CUTPOINT_UPPER) &&
				( ( restype.atom_name(atm2) == " N  " && restype.atom_name(atm3) == " CA ") ||
				( restype.atom_name(atm3) == " N  " && restype.atom_name(atm2) == " CA ")  ) ) {
			Kphi=0.0; phi0=0.0;
		}
		if ( restype.atom_name(atm1) == " CEN" || restype.atom_name(atm4) == " CEN" ) {
			Kphi=0.0; phi0=0.0;
		}
		if ( (restype.aa() == aa_cys || restype.aa() == aa_dcs) && restype.has_variant_type( chemical::DISULFIDE )  && //D- or L-cystine
				(restype.atom_name(atm2) == " SG " || restype.atom_name(atm3) == " SG ") ) {
			Kphi=0.0; phi0=0.0;
		}
		if ( (restype.aa() == aa_arg || restype.aa() == aa_dar)  && //D- or L-arginine
				(   restype.atom_name(atm1) == " NH1" || restype.atom_name(atm4) == " NH1"
				|| restype.atom_name(atm1) == " NH2" || restype.atom_name(atm4) == " NH2") ) {
			phi_step = numeric::constants::d::pi;
		}
		if ( restype.has_variant_type(chemical::LOWER_TERMINUS_VARIANT) &&
				( ( restype.atom_name(atm2) == " N  " && restype.atom_name(atm3) == " CA ") ||
				( restype.atom_name(atm3) == " N  " && restype.atom_name(atm2) == " CA ")  ) ) {
			Kphi=k_torsion_proton_;
			phi_step = 2.0*numeric::constants::d::pi/3.0;
		}

		if ( proton_chi ) {
			Kphi=k_torsion_proton_;
			utility::vector1< Real > const & v = restype.proton_chi_samples( 1 );
			phi0 = numeric::constants::d::deg2rad*v[1];
			phi_step = std::fabs( numeric::constants::d::deg2rad*(v[2]-v[1]) );  // i guess this assumes we have two adjacent samples
		}
	}
}

// old-style interface to database
// somewhat slow as it will always build the ideal residue on the fly
void
IdealParametersDatabase::lookup_angle_legacy(
	core::pose::Pose const & /*pose*/,
	core::conformation::Residue const & res,
	int atm1,
	int atm2,
	int atm3,
	Real & Ktheta,
	Real & theta0
) {
	using namespace core::chemical;

	Ktheta=k_angle_;
	ResidueType const & restype = res.type();
	//std::string restag = get_restag( restype );

	// Create mini-conformation for idealized residue
	conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );

	// use these for mm Ktheta lookup
	//Size mmatm1, mmatm2, mmatm3;

	// get angle
	numeric::xyzVector<core::Real> x,y,z; // atom coords
	if ( atm1 > 0 ) {
		x = newres->atom( atm1 ).xyz();
	} else {
		x = newres->residue_connection( -atm1 ).icoor().build( restype );
	}
	if ( atm2 > 0 ) {
		y = newres->atom( atm2 ).xyz();
	} else {
		y = newres->residue_connection( -atm2 ).icoor().build( restype );
	}
	if ( atm3 > 0 ) {
		z = newres->atom( atm3 ).xyz();
	} else {
		z = newres->residue_connection( -atm3 ).icoor().build( restype );
	}

	theta0 = numeric::angle_radians ( x,y,z );

	//fpd ignore proline C-ND which is handled by pro_close
	if ( restype.aa() == core::chemical::aa_pro || restype.aa() == core::chemical::aa_dpr ) { //D- or L-proline
		bool hasCD = (atm1>0 && restype.atom_name(atm1)==" CD ")
			|| (atm2>0 && restype.atom_name(atm2)==" CD ")
			|| (atm3>0 && restype.atom_name(atm3)==" CD ");
		bool hasN = (atm1>0 && restype.atom_name(atm1)==" N  ")
			|| (atm2>0 && restype.atom_name(atm2)==" N  ")
			|| (atm3>0 && restype.atom_name(atm3)==" N  ");
		if ( hasCD && hasN ) {
			Ktheta = theta0 = 0.0;
		}
	}

	if ( (restype.aa() == core::chemical::aa_ala || restype.aa() == core::chemical::aa_dal || restype.aa() == core::chemical::aa_gly) && //D- or L-alanine, or glycine
			( (atm1>0 && restype.atom_name(atm1) == " CEN") || (atm3>0 && restype.atom_name(atm3) == " CEN") ) ) {
		Ktheta = theta0 = 0.0;
	}
	if ( restype.has_variant_type(chemical::CUTPOINT_UPPER) &&
			( (atm1>0 && restype.atom_name(atm1) == "OVU1") || (atm3>0 && restype.atom_name(atm3) == "OVU1") ) ) {
		Ktheta = theta0 = 0.0;
	}
	if ( restype.has_variant_type(chemical::CUTPOINT_LOWER) &&
			( (atm1>0 && restype.atom_name(atm1) == "OVL1") || (atm3>0 && restype.atom_name(atm3) == "OVL1") ) ) {
		Ktheta = theta0 = 0.0;
	}

	return;
}

// old-style interface to database
// somewhat slow as it will always build the ideal residue on the fly
void
IdealParametersDatabase::lookup_length_legacy(
	core::pose::Pose const & /*pose*/,
	core::conformation::Residue const & res,
	int atm1,
	int atm2,
	Real &Kd,
	Real &d0
)
{
	using namespace core::chemical;

	Kd=k_length_;
	ResidueType const & restype = res.type();
	//std::string restag = get_restag( restype );

	// Create mini-conformation for idealized residue
	core::conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );

	// get length
	numeric::xyzVector<core::Real> x,y,z; // atom coords
	if ( atm1 > 0 ) {
		x = newres->atom( atm1 ).xyz();
	} else {
		x = newres->residue_connection( -atm1 ).icoor().build( restype );
	}
	if ( atm2 > 0 ) {
		y = newres->atom( atm2 ).xyz();
	} else {
		y = newres->residue_connection( -atm2 ).icoor().build( restype );
	}

	d0 = (x-y).length();

	if ( (restype.aa() == core::chemical::aa_pro || restype.aa() == core::chemical::aa_dpr) && //D- or L-proline
			(atm1>0 && atm2>0) &&
			( ( restype.atom_name(atm1) == " CD " && restype.atom_name(atm2) == " N  ") ||
			( restype.atom_name(atm2) == " CD " && restype.atom_name(atm1) == " N  ") ) ) {
		Kd = d0 = 0.0;
	}

	return;
}

ResidueCartBondedParameters const &
IdealParametersDatabase::parameters_for_restype(
	core::chemical::ResidueType const & restype,
	bool prepro
)
{
	std::map< chemical::ResidueType const *, ResidueCartBondedParametersOP > const * resdatamap;
	std::map< chemical::ResidueType const *, ResidueCartBondedParametersOP >::const_iterator iter;

	{
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( CartesianBondedEnergy::restype_db_mutex_ );
#endif

		resdatamap = prepro ? &prepro_restype_data_ : &nonprepro_restype_data_;
		iter = resdatamap->find( & restype );
	}

	if ( iter == resdatamap->end() ) {
#if defined MULTI_THREADED
		utility::thread::WriteLockGuard lock( CartesianBondedEnergy::restype_db_mutex_ );
#endif
		if ( resdatamap->find( & restype ) == resdatamap->end() ) {
			create_parameters_for_restype( restype, prepro );
		}
		iter = resdatamap->find( & restype );
	}
	return *iter->second;
}


void
IdealParametersDatabase::create_parameters_for_restype(
	core::chemical::ResidueType const & rsd_type,
	bool prepro
)
{
	ResidueCartBondedParametersOP restype_params( new ResidueCartBondedParameters );

	std::string const rsdname = get_restag(rsd_type);  //fpd don't use seperate logic here

	// Iter over parameters - this would be fast enough as far as parameter size is small enough
	for ( boost::unordered_map<atm_name_quad,CartBondedParametersOP>::iterator b_it =
			impropers_indep_.begin(); b_it != impropers_indep_.end(); ++b_it ) {

		atm_name_quad const &tuple( b_it->first );

		if ( rsdname != tuple.get<0>() ) continue;

		CartBondedParametersCOP tor_params = b_it->second;

		// Also skip if any atom does not exist
		if ( !rsd_type.has( tuple.get<1>() ) || !rsd_type.has( tuple.get<2>() ) ||
				!rsd_type.has( tuple.get<3>() ) || !rsd_type.has( tuple.get<4>() ) ) continue;

		ResidueCartBondedParameters::Size4 ids;
		ids[1] = rsd_type.atom_index( tuple.get<1>() );
		ids[2] = rsd_type.atom_index( tuple.get<2>() );
		ids[3] = rsd_type.atom_index( tuple.get<3>() );
		ids[4] = rsd_type.atom_index( tuple.get<4>() );

		restype_params->add_improper_parameter( ids, tor_params );
	}

	typedef boost::unordered_multimap< atm_name_quad, CartBondedParametersOP >::const_iterator tors_iterator;
	for ( tors_iterator it = torsions_indep_.begin(); it != torsions_indep_.end();
			it = torsions_indep_.equal_range(it->first).second ) {

		atm_name_quad const &tuple( it-> first );

		if ( rsdname != tuple.get<0>() ) continue;

		// Also skip if any atom does not exist
		if ( !rsd_type.has( tuple.get<1>() ) || !rsd_type.has( tuple.get<2>() ) ||
				!rsd_type.has( tuple.get<3>() ) || !rsd_type.has( tuple.get<4>() ) ) continue;

		std::pair< tors_iterator, tors_iterator > range = torsions_indep_.equal_range( tuple );
		tors_iterator it2;

		for ( it2 = range.first; it2 != range.second; ++it2 ) {
			CartBondedParametersCOP tor_params = it2->second;

			ResidueCartBondedParameters::Size4 ids;
			ids[1] = rsd_type.atom_index( tuple.get<1>() );
			ids[2] = rsd_type.atom_index( tuple.get<2>() );
			ids[3] = rsd_type.atom_index( tuple.get<3>() );
			ids[4] = rsd_type.atom_index( tuple.get<4>() );
			restype_params->add_torsion_parameter( ids, tor_params );
		}
	}

	// for each angle in the residue
	for ( Size bondang = 1; bondang <= rsd_type.num_bondangles(); ++bondang ) {
		// get ResidueType ints
		Size rt1 = ( rsd_type.bondangle( bondang ) ).key1();
		Size rt2 = ( rsd_type.bondangle( bondang ) ).key2();
		Size rt3 = ( rsd_type.bondangle( bondang ) ).key3();
		ResidueCartBondedParameters::Size3 ids;
		ids[ 1 ] = rt1; ids[ 2 ] = rt2; ids[ 3 ] = rt3;

		// lookup Ktheta and theta0
		std::string atm1name=rsd_type.atom_name(rt1); boost::trim(atm1name);
		std::string atm2name=rsd_type.atom_name(rt2); boost::trim(atm2name);
		std::string atm3name=rsd_type.atom_name(rt3); boost::trim(atm3name);
		CartBondedParametersCOP ang_params = lookup_angle( rsd_type, prepro,
			atm1name,atm2name,atm3name, rt1,rt2,rt3 );

		if ( ang_params->is_null() ) continue;

		restype_params->add_angle_parameter( ids, ang_params );
	}

	// for each bond
	for ( Size atm_i=1; atm_i<=rsd_type.natoms(); ++atm_i ) {
		chemical::AtomIndices atm_nbrs = rsd_type.nbrs( atm_i );
		for ( Size j=1; j<=atm_nbrs.size(); ++j ) {
			Size atm_j = atm_nbrs[j];
			if ( atm_i >= atm_j ) continue; // only score each bond once -- use restype index to define ordering

			ResidueCartBondedParameters::Size2 ids;
			ids[1] = atm_i; ids[2] = atm_j;

			// lookup Kd and d0
			std::string atm1name = rsd_type.atom_name(atm_i); boost::trim(atm1name);
			std::string atm2name = rsd_type.atom_name(atm_j); boost::trim(atm2name);
			CartBondedParametersCOP len_params = lookup_length(rsd_type, prepro,
				atm1name,atm2name, atm_i, atm_j );

			if ( len_params->is_null() ) continue;

			restype_params->add_length_parameter( ids, len_params );

		}
	}

	// Virtual shadows: for the sake of working with virtual atoms that are
	// used to keep intra-residue rings closed, apply a loose constraint between
	// the virtual shadow on top of the atom it's shadowing.
	for ( Size ii = 1; ii <= rsd_type.natoms(); ++ii ) {
		Size ii_shadowee = rsd_type.atom_being_shadowed(ii);
		if ( ii_shadowee == 0 ) continue;
		ResidueCartBondedParameters::Size2 ids;
		ids[1] = ii; ids[2] = ii_shadowee;

		// weak constraint; distance of 0, spring constant of 10
		CartBondedParametersOP params_ii( new BBIndepCartBondedParameters( 0, 10 ) );
		restype_params->add_length_parameter( ids, params_ii );
	}


	/// Keep track of the bond angles and lengths that are used in the backbone dependent calculations

	// loop over all bb-dep angles and bonds;
	// hardcode the bbdep angles and lengths to save a bit of time
	// connection IDs are hardcoded ... is this a problem? APL: Yes, c- or n- terminal disulfides would give you
	// trouble, but the connection IDs don't have to be hard coded.

	//fpd protein only
	if ( rsd_type.is_protein() ) {
		/// backbone dependent bond lengths
		bool is_nterm = ( (rsd_type.aa() <= chemical::num_canonical_aas || core::chemical::is_canonical_D_aa(rsd_type.aa())) && rsd_type.is_lower_terminus()); //Modified by VKM to check for D-amino acids
		bool is_cterm = ( (rsd_type.aa() <= chemical::num_canonical_aas || core::chemical::is_canonical_D_aa(rsd_type.aa())) && rsd_type.is_upper_terminus()); //Modified by VKM to check for D-amino acids
		for ( int i=1; i<=5; ++i ) {
			if ( i==1 && is_nterm ) continue;
			if ( i==3 && ( rsd_type.aa() == core::chemical::aa_gly || rsd_type.aa() == core::chemical::aa_b3g /*"beta-glycine"*/ ) ) continue;

			std::string atm1,atm2;
			int rt1 = 0;
			int rt2 = 0;
			if ( i==1 ) { atm1="C";  atm2="N";  rt1 = -rsd_type.lower_connect_id(); rt2 = rsd_type.atom_index(" N  ");}
			if ( i==2 ) { atm1="N";  atm2="CA"; rt1 = rsd_type.atom_index(" N  "); rt2 = rsd_type.atom_index(" CA "); }
			if ( i==3 ) { atm1="CA"; atm2="CB"; rt1 = rsd_type.atom_index(" CA "); rt2 = rsd_type.atom_index(" CB "); }
			if ( i==4 ) { atm1="CA"; atm2="C";  rt1 = rsd_type.atom_index(" CA "); rt2 = rsd_type.atom_index(" C  "); }
			if ( i==5 ) { atm1="C";  atm2="O";  rt1 = rsd_type.atom_index(" C  "); rt2 = rsd_type.atom_index(" O  "); }

			ResidueCartBondedParameters::Size2 ids;
			ids[1] = (core::Size)std::max(rt1,0);
			ids[2] = (core::Size)std::max(rt2,0);

			CartBondedParametersCOP len_params = lookup_length( rsd_type, is_nterm ? false : prepro, atm1, atm2, rt1, rt2 );
			restype_params->add_bbdep_length_parameter( ids, len_params );
		}

		// backbone dependent bond angles
		for ( int i=1; i<=7; ++i ) {
			if ( (i==2 || i==4) && ( rsd_type.aa() == core::chemical::aa_gly || rsd_type.aa() == core::chemical::aa_b3g /*"beta-glycine"*/ ) ) continue;
			if ( i==1 && is_nterm ) continue;
			if ( (i==6 || i==7) && is_cterm ) continue;

			std::string atm1,atm2,atm3;
			int rt1 = 0;
			int rt2 = 0;
			int rt3 = 0;
			if ( i==1 ) { atm1="C";  atm2="N";  atm3="CA"; rt1=-rsd_type.lower_connect_id();  rt2=rsd_type.atom_index(" N  "); rt3=rsd_type.atom_index(" CA ");}
			if ( i==2 ) { atm1="N";  atm2="CA"; atm3="CB"; rt1=rsd_type.atom_index(" N  "); rt2=rsd_type.atom_index(" CA "); rt3=rsd_type.atom_index(" CB "); }
			if ( i==3 ) { atm1="N";  atm2="CA"; atm3="C";  rt1=rsd_type.atom_index(" N  "); rt2=rsd_type.atom_index(" CA "); rt3=rsd_type.atom_index(" C  "); }
			if ( i==4 ) { atm1="CB"; atm2="CA"; atm3="C";  rt1=rsd_type.atom_index(" CB "); rt2=rsd_type.atom_index(" CA "); rt3=rsd_type.atom_index(" C  "); }
			if ( i==5 ) { atm1="CA"; atm2="C";  atm3="O";  rt1=rsd_type.atom_index(" CA "); rt2=rsd_type.atom_index(" C  "); rt3=rsd_type.atom_index(" O  "); }
			if ( i==6 ) { atm1="CA"; atm2="C";  atm3="N";  rt1=rsd_type.atom_index(" CA "); rt2=rsd_type.atom_index(" C  "); rt3=-rsd_type.upper_connect_id(); }
			if ( i==7 ) { atm1="O";  atm2="C";  atm3="N";  rt1=rsd_type.atom_index(" O  "); rt2=rsd_type.atom_index(" C  "); rt3=-rsd_type.upper_connect_id(); }

			ResidueCartBondedParameters::Size3 ids;
			ids[1] = (core::Size)std::max(rt1,0); ids[2] = (core::Size)std::max(rt2,0); ids[3] = (core::Size)std::max(rt3,0);

			CartBondedParametersCOP ang_params = lookup_angle(rsd_type, prepro, atm1,atm2,atm3, rt1,rt2,rt3 );
			restype_params->add_bbdep_angle_parameter( ids, ang_params );
		}
	}

	// atom pairs near the upper connect:
	{
		if ( rsd_type.upper_connect_id() != 0 ) {
			utility::vector1< chemical::two_atom_set > const & upconn_ats(
				rsd_type.atoms_within_one_bond_of_a_residue_connection( rsd_type.upper_connect_id() ));
			for ( Size ii = 1; ii <= upconn_ats.size(); ++ii ) {

				std::string atm1name=rsd_type.atom_name( upconn_ats[ii].key2() ); boost::trim(atm1name);
				std::string atm2name=rsd_type.atom_name( upconn_ats[ii].key1() ); boost::trim(atm2name);
				std::string atm3name="N";
				if ( rsd_type.is_RNA() ) atm3name="P";

				// lookup Ktheta and theta0
				CartBondedParametersCOP ang_params = lookup_angle( rsd_type, prepro,
					atm1name, atm2name, atm3name, upconn_ats[ii].key2(), upconn_ats[ii].key1(),
					-1 * ( (int) rsd_type.upper_connect_id() ));
				ResidueCartBondedParameters::Size3 inds;
				inds[1] = upconn_ats[ii].key2();
				inds[2] = upconn_ats[ii].key1();
				inds[3] = 0;
				restype_params->add_upper_connect_angle_params( inds, ang_params );
			}
		}
	}

	// atom pairs near the lower connect:
	{
		if ( rsd_type.lower_connect_id() != 0 ) {
			utility::vector1< chemical::two_atom_set > const & loconn_ats(
				rsd_type.atoms_within_one_bond_of_a_residue_connection( rsd_type.lower_connect_id() ));
			for ( Size ii = 1; ii <= loconn_ats.size(); ++ii ) {

				std::string atm1name=rsd_type.atom_name( loconn_ats[ii].key2() ); boost::trim(atm1name);
				std::string atm2name=rsd_type.atom_name( loconn_ats[ii].key1() ); boost::trim(atm2name);
				std::string atm3name="C";
				if ( rsd_type.is_RNA() ) atm3name="O3'";

				// lookup Ktheta and theta0
				CartBondedParametersCOP ang_params = lookup_angle( rsd_type, prepro,
					atm1name, atm2name, atm3name, loconn_ats[ii].key2(), loconn_ats[ii].key1(),
					-1 * ( (int) rsd_type.lower_connect_id() ));
				ResidueCartBondedParameters::Size3 inds;
				inds[1] = loconn_ats[ii].key2();
				inds[2] = loconn_ats[ii].key1();
				inds[3] = 0;
				restype_params->add_lower_connect_angle_params( inds, ang_params );
			}
		}
	}

	/// Atom IDs
	if ( rsd_type.has( "N" ) ) {
		restype_params->bb_N_index( rsd_type.atom_index( "N" ) );
	}
	if ( rsd_type.has( "CA" ) ) {
		restype_params->bb_CA_index( rsd_type.atom_index( "CA" ) );
	}
	if ( rsd_type.has( "C" ) ) {
		restype_params->bb_C_index( rsd_type.atom_index( "C" ) );
	}
	if ( rsd_type.has( "O" ) ) {
		restype_params->bb_O_index( rsd_type.atom_index( "O" ) );
	}
	if ( rsd_type.has( "H" ) ) {
		restype_params->bb_H_index( rsd_type.atom_index( "H" ) );
	}
	if ( (rsd_type.aa() == core::chemical::aa_pro || rsd_type.aa() == core::chemical::aa_dpr) && rsd_type.has( "CD" ) ) {
		restype_params->pro_CD_index( rsd_type.atom_index( "CD" ) );
	}

	/// oprev_cprev_n_h improper torsion
	if ( restype_params->bb_N_index() != 0 && restype_params->bb_H_index() != 0 ) {
		CartBondedParametersCOP tor_params = lookup_improper( rsd_type, "O", "C", "N", "H" );
		restype_params->oprev_cprev_n_h_interres_improper_params( tor_params );
	}

	/// ca_cprev_n_h improper torsion
	if ( restype_params->bb_N_index() != 0 && restype_params->bb_H_index() != 0 && restype_params->bb_CA_index() != 0 ) {
		CartBondedParametersCOP tor_params = lookup_improper( rsd_type, "CA", "C", "N", "H" );
		restype_params->ca_cprev_n_h_interres_improper_params( tor_params );
	}

	// ca_nnext_c_o improper torsion
	if ( restype_params->bb_C_index() != 0 && restype_params->bb_O_index() != 0 && restype_params->bb_CA_index() != 0 ) {
		CartBondedParametersCOP tor_params = lookup_improper( rsd_type, "CA", "N", "C", "O" );
		restype_params->ca_nnext_c_o_interres_improper_params( tor_params );
	}

	// ca_nnext_c_o improper torsion
	if ( restype_params->bb_N_index() != 0 && restype_params->pro_CD_index() != 0 && restype_params->bb_CA_index() != 0 ) {
		CartBondedParametersCOP tor_params = lookup_improper( rsd_type, "CD", "C", "N", "CA" );
		restype_params->pro_cd_cprev_n_ca_interres_improper_params( tor_params );
	}

	// caprev_n bond length
	if ( restype_params->bb_N_index() != 0 && rsd_type.lower_connect_id() != 0 ) {
		CartBondedParametersCOP len_params = lookup_length( rsd_type, prepro, "N", "C", restype_params->bb_N_index(), -1 * ((int) rsd_type.lower_connect_id() ) );
		debug_assert( len_params );
		restype_params->cprev_n_bond_length_params( len_params );
	}

	if ( prepro ) {
		prepro_restype_data_[ & rsd_type ] = restype_params;
	} else {
		nonprepro_restype_data_[ & rsd_type ] = restype_params;
	}

}


//////////////////////
/// EnergyMethod
CartesianBondedEnergy::CartesianBondedEnergy( methods::EnergyMethodOptions const & options ) :
	parent( methods::EnergyMethodCreatorOP( new CartesianBondedEnergyCreator ) ),
	pro_nv_("NV")
{
	// if flag _or_ energy method wants a linear potential, make the potential linear
	linear_bonded_potential_ =
		basic::options::option[ basic::options::OptionKeys::score::linear_bonded_potential ]() ||
		options.get_cartesian_bonded_linear();

	// initialize databases
	core::Real cartbonded_len, cartbonded_ang, cartbonded_tors, cartbonded_proton , cartbonded_improper;
	options.get_cartesian_bonded_parameters( cartbonded_len, cartbonded_ang, cartbonded_tors, cartbonded_proton , cartbonded_improper );

	if ( !db_ ) {
#if defined MULTI_THREADED
		utility::thread::WriteLockGuard lock( params_db_mutex_ );
#endif
		if ( !db_ ) {
			db_ = IdealParametersDatabaseOP( new IdealParametersDatabase(cartbonded_len, cartbonded_ang, cartbonded_tors, cartbonded_proton , cartbonded_improper) );
		}
	}
}

CartesianBondedEnergy::CartesianBondedEnergy( CartesianBondedEnergy const & src ) : parent( src ) {
	linear_bonded_potential_ = src.linear_bonded_potential_;
	pro_nv_ = src.pro_nv_;
}

CartesianBondedEnergy::~CartesianBondedEnergy() {}

EnergyMethodOP
CartesianBondedEnergy::clone() const {
	return EnergyMethodOP( new CartesianBondedEnergy( *this ) );
}

void
CartesianBondedEnergy::setup_for_derivatives(
	pose::Pose &,
	ScoreFunction const &
) const {
	//idealize_proline_nvs(pose);
}

void
CartesianBondedEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const {
	using namespace methods;

	// create LR energy container
	LongRangeEnergyType const & lr_type( long_range_type() );
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		PolymerBondedEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::PolymerBondedEnergyContainer > ( lrc ) );
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
		s_types.push_back( cart_bonded );
		s_types.push_back( cart_bonded_angle );
		s_types.push_back( cart_bonded_length );
		s_types.push_back( cart_bonded_torsion );
		s_types.push_back( cart_bonded_ring );
		s_types.push_back( cart_bonded_proper );
		s_types.push_back( cart_bonded_improper );
		LREnergyContainerOP new_dec( new PolymerBondedEnergyContainer( pose, s_types ) );
		energies.set_long_range_container( lr_type, new_dec );
	}
}

void
CartesianBondedEnergy::idealize_proline_nvs(
	pose::Pose & pose
) const {

	// Idealize the NV atom of the prolines in the pose. This allows
	// for safe switching between cartesian and non-cartesian poses
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( pose.is_fullatom() && (pose.residue(i).aa() == core::chemical::aa_pro || pose.residue(i).aa() == core::chemical::aa_dpr) ) {
			core::Size nv_index = pose.residue(i).atom_index(pro_nv_);
			core::id::AtomID nv_id(nv_index, i);
			Vector nv_coords = pose.residue(i).build_atom_ideal((int)nv_index, pose.conformation());
			pose.set_xyz(nv_id, nv_coords);
		}
	}
}


bool
CartesianBondedEnergy::defines_residue_pair_energy(
	pose::Pose const &,
	Size ,
	Size
) const {
	// is this fn. called?
	// VKM -- 10 Sept 2016: No, no it doesn't seem to be.
	// FD in that case, since logic below is no longer correct, just return true
	return ( true );
}


void
CartesianBondedEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sf,
	EnergyMap & emap
) const
{
	bool res1first = rsd2.seqpos() > rsd1.seqpos();

	// override ordering if one residue has an upper connection to another (cyclic)
	if ( rsd1.has_upper_connect() && rsd1.connected_residue_at_resconn( rsd1.type().upper_connect_id() ) == rsd2.seqpos() ) {
		res1first = true;
	} else if ( rsd2.has_upper_connect() && rsd2.connected_residue_at_resconn( rsd2.type().upper_connect_id() ) == rsd1.seqpos() ) {
		res1first = false;
	}

	if ( res1first ) {
		residue_pair_energy_sorted( rsd1, rsd2, pose, sf, emap );
	} else { //Assumes that residue 1 is connected to residue 2's C-terminus.
		residue_pair_energy_sorted( rsd2, rsd1, pose, sf, emap );
	}
}


void
CartesianBondedEnergy::eval_intrares_energy(
	conformation::Residue const &rsd,
	pose::Pose const &pose,
	ScoreFunction const &,
	EnergyMap &emap
) const
{
	using namespace numeric;

	//(fpd) NOTE: this must agree with logic in PolymerBondedEnergyContainer::other_res_index
	//    that is, this function must be evaluated for all residues that do not take part in a two-body energy
	bool dont_score_this = ( rsd.type().is_polymer() && rsd.has_upper_connect() );
	if ( dont_score_this ) {
		core::Size const other_res( rsd.residue_connection_partner( rsd.upper_connect().index() ) );
		dont_score_this = (
			other_res != 0
			&& pose.residue(other_res).type().is_polymer()
			&& pose.residue(other_res).has_lower_connect()
			&& pose.residue(other_res).residue_connection_partner( pose.residue(other_res).lower_connect().index() ) == rsd.seqpos()
		);
	}
	if ( dont_score_this ) return;

	const core::Real d_multiplier = core::chemical::is_canonical_D_aa(rsd.aa()) ? -1.0 : 1.0 ;
	Real phi=0,psi=0;
	if ( rsd.is_protein() ) {
		phi = nonnegative_principal_angle_degrees( d_multiplier * rsd.mainchain_torsion(1));
		psi = nonnegative_principal_angle_degrees( d_multiplier * rsd.mainchain_torsion(2));
	}

	ResidueCartBondedParameters const * rsdparams;
	{
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
		rsdparams = &(db_->parameters_for_restype( rsd.type(), false ));
	}
	if ( rsd.aa() != core::chemical::aa_vrt ) {
		eval_singleres_energy(rsd, *rsdparams, phi, psi, pose, emap ); // calls singleres improper
	}
}


void
CartesianBondedEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & /*res_data_cache*/,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const {
	using namespace numeric;

	//(fpd) NOTE: this must agree with logic in PolymerBondedEnergyContainer::other_res_index
	//    that is, this function must be evaluated for all residues that do not take part in a two-body energy
	bool dont_score_this = ( rsd.type().is_polymer() && rsd.has_upper_connect() );
	if ( dont_score_this ) {
		core::Size const other_res( rsd.residue_connection_partner( rsd.upper_connect().index() ) );
		dont_score_this = (
			other_res != 0
			&& pose.residue(other_res).type().is_polymer()
			&& pose.residue(other_res).has_lower_connect()
			&& pose.residue(other_res).residue_connection_partner( pose.residue(other_res).lower_connect().index() ) == rsd.seqpos()
		);
	}
	if ( dont_score_this ) return;

	const core::Real d_multiplier = core::chemical::is_canonical_D_aa(rsd.aa()) ? -1.0 : 1.0 ;
	Real phi=0,psi=0;
	if ( rsd.is_protein() ) {
		phi = nonnegative_principal_angle_degrees( d_multiplier * rsd.mainchain_torsion(1));
		psi = nonnegative_principal_angle_degrees( d_multiplier * rsd.mainchain_torsion(2));
	}

	ResidueCartBondedParameters const * rsdparams;
	{
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
		rsdparams = &(db_->parameters_for_restype( rsd.type(), false ));
	}

	if ( rsd.aa() != core::chemical::aa_vrt ) {
		eval_singleres_derivatives( rsd, *rsdparams, phi, psi, weights, atom_derivs );
	}
}

///////////////////////////
///
void
CartesianBondedEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const & min1,
	ResSingleMinimizationData const & min2,
	ResPairMinimizationData const & min12,
	pose::Pose const & pose,
	EnergyMap const & wts,
	utility::vector1< DerivVectorPair > & r1_derivs,
	utility::vector1< DerivVectorPair > & r2_derivs
) const {
	bool res1first = rsd2.seqpos() > rsd1.seqpos();

	// override ordering if one residue has an upper connection to another (cyclic)
	if ( rsd1.has_upper_connect() && rsd1.connected_residue_at_resconn( rsd1.type().upper_connect_id() ) == rsd2.seqpos() ) {
		res1first = true;
	} else if ( rsd2.has_upper_connect() && rsd2.connected_residue_at_resconn( rsd2.type().upper_connect_id() ) == rsd1.seqpos() ) {
		res1first = false;
	}

	if ( res1first ) {
		eval_residue_pair_derivatives_sorted( rsd1, rsd2, min1, min2, min12, pose, wts, r1_derivs, r2_derivs );
	} else { //Assumes that residue 1 is connected to residue 2's C-terminus.
		eval_residue_pair_derivatives_sorted( rsd1, rsd2, min1, min2, min12, pose, wts, r1_derivs, r2_derivs );
	}
}

void
CartesianBondedEnergy::eval_residue_pair_derivatives_sorted(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & /*min_data*/,
	pose::Pose const & /*pose*/,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	using namespace numeric;

	bool preproline = (rsd2.aa()==core::chemical::aa_pro || rsd2.aa()==core::chemical::aa_dpr); //Is rsd2 either D-proline or L-proline?

	//Multipliers for D-amino acids:
	const core::Real d_multiplier1 = core::chemical::is_canonical_D_aa(rsd1.aa()) ? -1.0 : 1.0 ;
	const core::Real d_multiplier2 = core::chemical::is_canonical_D_aa(rsd2.aa()) ? -1.0 : 1.0 ;


	ResidueCartBondedParameters const *res1params, *res2params;
	{
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
		res1params = &(db_->parameters_for_restype( rsd1.type(), preproline ) );
		res2params = &(db_->parameters_for_restype( rsd2.type(), false ) );
	}

	// get phi,psis for bb-dep angles/lengths
	Real phi1=0,psi1=0,phi2=0,psi2=0;
	if ( rsd1.is_protein() ) {
		phi1 = nonnegative_principal_angle_degrees( d_multiplier1*rsd1.mainchain_torsion(1)); //Inverted if D-amino acid
		psi1 = nonnegative_principal_angle_degrees( d_multiplier1*rsd1.mainchain_torsion(2)); //Inverted if D-amino acid
	}
	if ( rsd2.is_protein() ) {
		phi2 = nonnegative_principal_angle_degrees( d_multiplier2*rsd2.mainchain_torsion(1)); //Inverted if D-amino acid
		psi2 = nonnegative_principal_angle_degrees( d_multiplier2*rsd2.mainchain_torsion(2)); //Inverted if D-amino acid
	}

	// get subcomponents
	if ( rsd1.aa() != core::chemical::aa_vrt ) {
		eval_singleres_derivatives( rsd1, *res1params, phi1, psi1, weights, r1_atom_derivs );
	}

	// If residue1 and 2 are the same (signal used by the PolymerBondedEnergyContainer for residues that aren't polymer-bonded), stop here to avoid double-counting.
	if ( rsd1.seqpos() == rsd2.seqpos() ) return;

	if ( rsd1.aa() == core::chemical::aa_vrt ) return;

	// bail out if the residues aren't bonded or we cross a cutpoint
	if ( !rsd1.is_bonded(rsd2) ) { return; }
	if ( rsd1.has_variant_type(core::chemical::CUTPOINT_LOWER) ) { return; }
	if ( rsd2.has_variant_type(core::chemical::CUTPOINT_UPPER) ) { return; }

	eval_interresidue_improper_derivatives( rsd1, rsd2, *res1params, *res2params, weights, r1_atom_derivs, r2_atom_derivs );

	eval_interresidue_angle_derivs_two_from_rsd1(
		rsd1, rsd2, *res1params, *res2params, phi1, psi1,
		weights, r1_atom_derivs, r2_atom_derivs );

	eval_interresidue_angle_derivs_two_from_rsd2(
		rsd1, rsd2, *res1params, *res2params, phi2, psi2,
		weights, r1_atom_derivs, r2_atom_derivs );

	eval_interresidue_bond_length_derivs(
		rsd1, rsd2, *res1params, *res2params, phi2, psi2,
		weights, r1_atom_derivs, r2_atom_derivs );
}

/// @details Evaluate dE/dphi and dE/dpsi for the backbone-dependent
/// terms, if they are active.
Real
CartesianBondedEnergy::eval_intraresidue_dof_derivative(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & /*min_data*/,
	id::DOF_ID const & /*dof_id*/,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	ScoreFunction const & /*sfxn*/,
	EnergyMap const & weights
) const {
	using namespace numeric;
	using numeric::constants::d::pi;

	// save some time if we're only doing bb-indep
	{
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
		if ( !db_->bbdep_bond_params() ) return 0.0;
	}

	if ( !tor_id.valid() || tor_id.type()!=id::BB || tor_id.torsion() > 2  || !rsd.is_protein() ) {
		return 0.0;
	}

	core::Size const iplus1_resid( rsd.has_upper_connect() ? rsd.connected_residue_at_resconn( rsd.type().upper_connect_id() ) : 0 ); //Index of residue connected at C-terminal connection; 0 if no connection there.

	//i+1 residue is either D-proline or L-proline
	bool const preproline ( iplus1_resid ? pose.residue( iplus1_resid ).aa() == core::chemical::aa_pro || pose.residue( iplus1_resid ).aa() == core::chemical::aa_dpr : false ); //Is this a residue preceding a proline?

	// phi/psi
	Real phi=0,psi=0;
	if ( rsd.is_protein() ) {
		phi = nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1));
		psi = nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2));
	}

	ResidueCartBondedParameters const *resparams;
	bool bbdepbonds;
	{
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
		resparams = &(db_->parameters_for_restype( rsd.type(), preproline ) );
		bbdepbonds = db_->bbdep_bond_devs();
	}

	core::Real deriv=0.0;

	/// Backbone dependent bond lengths
	utility::vector1< ResidueCartBondedParameters::length_parameter > const & lps = resparams->bbdep_length_parameters();
	for ( Size ii = 1; ii <= lps.size(); ++ii ) {
		ResidueCartBondedParameters::Size2 const & atids( lps[ ii ].first );
		CartBondedParameters const & len_params( *lps[ ii ].second );

		Real const K = len_params.K(phi,psi);
		Real const mu = len_params.mu(phi,psi);
		Vector at1xyz, at2xyz( rsd.xyz( atids[2] ));
		if ( atids[1] == 0 ) {
			// connection atom from previous residue
			chemical::ResidueType const & prevrestype = pose.residue_type( rsd.seqpos() - 1 );
			Size prev_c_atom = prevrestype.upper_connect().atomno();
			at1xyz = pose.xyz( id::AtomID( prev_c_atom, rsd.seqpos() - 1 ));
		} else {
			at1xyz = rsd.xyz( atids[1] );
		}

		Real const d = at1xyz.distance( at2xyz );

		core::Real dmu_dtor, dK_dtor=0.0;
		core::Real dscore_dK=0.0, dscore_dmu;

		if ( tor_id.torsion()==1 ) {
			dmu_dtor = len_params.dmu_dphi(phi,psi);
		} else {
			dmu_dtor = len_params.dmu_dpsi(phi,psi);
		}
		if ( linear_bonded_potential_ && std::fabs(d - mu)>1 ) {
			dscore_dmu = -K * ((d - mu)>0 ? 0.5 : -0.5);
		} else {
			dscore_dmu = -K * (d - mu);
		}
		if ( bbdepbonds ) {
			if ( tor_id.torsion()==1 ) {
				dK_dtor = len_params.dK_dphi(phi,psi);
			} else {
				dK_dtor = len_params.dK_dpsi(phi,psi);
			}
			if ( linear_bonded_potential_ && std::fabs(d - mu)>1 ) {
				dscore_dK = 0.5*std::fabs(d - mu);
			} else {
				dscore_dK = 0.5*(d - mu)*(d - mu);   // currently we have no -log(K) term in our score
			}
		}

		// derivatives w.r.t. phi/psi
		deriv += (weights[ cart_bonded_length ] + weights[ cart_bonded ]) * (dscore_dmu*dmu_dtor + dscore_dK*dK_dtor);
	}

	/// Backbone-dependent bond angles
	utility::vector1< ResidueCartBondedParameters::angle_parameter > const & aps = resparams->bbdep_angle_parameters();
	for ( Size ii = 1; ii <= aps.size(); ++ii ) {
		ResidueCartBondedParameters::Size3 const & atids( aps[ ii ].first );
		CartBondedParameters const & ang_params( *aps[ ii ].second );

		Real const K = ang_params.K(phi,psi);
		Real const mu = ang_params.mu(phi,psi);
		Vector at1xyz, at2xyz( rsd.xyz( atids[2] )), at3xyz;
		if ( atids[1] == 0 ) {
			// connection atom from previous residue
			chemical::ResidueType const & prevrestype = pose.residue_type( rsd.seqpos() - 1 );
			Size prev_c_atom = prevrestype.upper_connect().atomno();
			at1xyz = pose.xyz( id::AtomID( prev_c_atom, rsd.seqpos() - 1 ));
		} else {
			at1xyz = rsd.xyz( atids[1] );
		}
		if ( atids[3] == 0 ) {
			// connection atom from the next residue
			chemical::ResidueType const & nextrestype = pose.residue_type( rsd.seqpos() + 1 );
			Size next_n_atom = nextrestype.lower_connect().atomno();
			at3xyz = pose.xyz( id::AtomID( next_n_atom, rsd.seqpos() + 1 ));
		} else {
			at3xyz = rsd.xyz( atids[3] );
		}

		Real const angle = numeric::angle_radians(at1xyz,at2xyz,at3xyz);

		core::Real dmu_dtor, dK_dtor=0.0;
		core::Real dscore_dK=0.0, dscore_dmu;

		if ( tor_id.torsion()==1 ) {
			dmu_dtor = ang_params.dmu_dphi(phi,psi);
		} else {
			dmu_dtor = ang_params.dmu_dpsi(phi,psi);
		}
		if ( linear_bonded_potential_ && std::fabs(angle - mu)>1 ) {
			dscore_dmu = -K * ((angle - mu)>0 ? 0.5 : -0.5);
		} else {
			dscore_dmu = -K * (angle - mu);
		}

		if ( bbdepbonds ) {
			if ( tor_id.torsion()==1 ) {
				dK_dtor = ang_params.dK_dphi(phi,psi);
			} else {
				dK_dtor = ang_params.dK_dpsi(phi,psi);
			}
			if ( linear_bonded_potential_ && std::fabs(angle - mu)>1 ) {
				dscore_dK = 0.5*std::fabs(angle - mu);
			} else {
				dscore_dK = 0.5*(angle - mu)*(angle - mu);   // currently we have no -log(K) term in our score
			}
		}

		// derivatives w.r.t. phi/psi
		deriv += (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * (dscore_dmu*dmu_dtor + dscore_dK*dK_dtor);
	}

	return 180/pi * deriv;
}

/// @brief CartesianBondedEnergy does not have an atomic interation threshold
Distance
CartesianBondedEnergy::atomic_interaction_cutoff() const {
	return 0.0;
}


/// @brief CartesianBondedEnergy is context independent; indicates that no context graphs are required
void
CartesianBondedEnergy::indicate_required_context_graphs(utility::vector1< bool > & ) const {}


methods::LongRangeEnergyType
CartesianBondedEnergy::long_range_type() const { return methods::cart_bonded_lr; }

void
CartesianBondedEnergy::residue_pair_energy_sorted(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {

	using namespace numeric;

	bool const preproline(
		rsd1.seqpos() != rsd2.seqpos() && //These are two different residues?
		( rsd1.has_upper_connect() && rsd1.connected_residue_at_resconn( rsd1.type().upper_connect_id() ) == rsd2.seqpos() ) && //The second residue is connected to the first residue's C-terminus?
		(pose.residue( rsd2.seqpos() ).aa() == core::chemical::aa_pro || pose.residue( rsd2.seqpos() ).aa() == core::chemical::aa_dpr) //The second residue is D/L proline?
	);


	//Multipliers for D-amino acids:
	const core::Real d_multiplier1 = core::chemical::is_canonical_D_aa(rsd1.aa()) ? -1.0 : 1.0 ;
	const core::Real d_multiplier2 = core::chemical::is_canonical_D_aa(rsd2.aa()) ? -1.0 : 1.0 ;

	ResidueCartBondedParameters const *res1params, *res2params;
	{
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
		res1params = &(db_->parameters_for_restype( rsd1.type(), preproline ) );
		res2params = &(db_->parameters_for_restype( rsd2.type(), false ) );
	}

	Real phi1=0,psi1=0,phi2=0,psi2=0;
	if ( rsd1.is_protein() ) {
		phi1 = nonnegative_principal_angle_degrees( d_multiplier1 * rsd1.mainchain_torsion(1));
		psi1 = nonnegative_principal_angle_degrees( d_multiplier1 * rsd1.mainchain_torsion(2));
	}
	if ( rsd2.is_protein() ) {
		phi2 = nonnegative_principal_angle_degrees( d_multiplier2 * rsd2.mainchain_torsion(1));
		psi2 = nonnegative_principal_angle_degrees( d_multiplier2 * rsd2.mainchain_torsion(2));
	}

	// get one body component (but which has two-body influence based on whether or not rsd2 is a proline)
	if ( rsd1.aa() != core::chemical::aa_vrt ) {
		eval_singleres_energy(rsd1, *res1params, phi1, psi1, pose, emap ); // calls singleres improper
	}

	// If residue1 and 2 are the same, stop here to avoid double-counting.
	if ( !rsd1.is_bonded(rsd2) ) return;
	if ( rsd1.has_variant_type(core::chemical::CUTPOINT_LOWER) ) return;
	if ( rsd2.has_variant_type(core::chemical::CUTPOINT_UPPER) ) return;

	/// evaluate all the inter-residue energy components
	eval_residue_pair_energies( rsd1, rsd2, *res1params, *res2params, phi1, psi1, phi2, psi2, pose, emap );
}

void
CartesianBondedEnergy::eval_singleres_energy(
	conformation::Residue const & rsd,
	ResidueCartBondedParameters const & resparams,
	Real phi, //Must be inverted for D-amino acids
	Real psi, //Must be inverted for D-amino acids
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	using numeric::constants::d::pi;
	using namespace numeric;

	debug_assert ( rsd.aa() != core::chemical::aa_vrt );

	eval_singleres_torsion_energies( rsd, resparams, pose, emap );
	eval_singleres_improper_energies( rsd, resparams, phi, psi, pose, emap );
	eval_singleres_angle_energies(   rsd, resparams, phi, psi, pose, emap );
	eval_singleres_length_energies(  rsd, resparams, phi, psi, pose, emap );
	try {
		eval_singleres_ring_energies( rsd, pose, emap );//removed resparams
	} catch ( utility::excn::EXCN_Base const & e ) {
		pose.dump_pdb("badstructure.pdb");
		std::cout << "caught exception " << e.msg() << std::endl;
		exit(0);
	}
//eval_singleres_ring_energies( rsd, pose, emap );//removed resparams
}

void
CartesianBondedEnergy::eval_singleres_ring_energies(
	conformation::Residue const & rsd,
	//ResidueCartBondedParameters const & resparams,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	using namespace core::chemical;
	using numeric::constants::d::pi;

	Size const n_rings( rsd.type().n_rings() );
	if ( n_rings==0 ) return;

	core::Real Ktheta, Kphi;
	{
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
		Ktheta = db_->k_torsion();  // for now use default torsion (TO DO: add spring constants to DB!)
		Kphi = db_->k_angle();  // for now use default torsion (TO DO: add spring constants to DB!)
	}

	for ( core::uint jj( 1 ); jj <= n_rings; ++jj ) {
		// get the conformer of the ring
		core::chemical::rings::RingConformer rc = rsd.ring_conformer( jj, 180.0 );

		// now constrain each element of the ring
		utility::vector1< core::Size > atms = rsd.type().ring_atoms( jj );

		for ( core::Size ii=1; ii<=atms.size(); ++ii ) {
			core::Size atm1 = atms[(ii+atms.size()-2)%atms.size()+1];
			core::Size atm2 = atms[ii];
			core::Size atm3 = atms[ii%atms.size()+1];
			core::Size atm4 = atms[(ii+1)%atms.size()+1];

			// 1 constrain angle
			Real theta0 = rc.tau_angles[ii]*pi/180.0;
			Real angle = numeric::angle_radians( rsd.xyz(atm1), rsd.xyz(atm2), rsd.xyz(atm3) );
			Real const energy_angle = eval_score( angle, Ktheta, theta0 );

			// Send a message to the user about a bad angle, if necessary.
			// Make sure not to send output to a tracer in the middle of
			// scoring unless that tracer is visible, since that can be very expensive
			if ( energy_angle > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " <<
					pose.pdb_info()->number(rsd.seqpos()) << " RING angle: " <<
					get_restag(rsd.type()) << " : " <<
					rsd.atom_name( atm1 ) << " , " << rsd.atom_name( atm2 ) << " , " <<
					rsd.atom_name( atm3 ) << "   (" <<
					Ktheta << ") " << 180/pi * angle << " " << 180/pi * theta0 << "    sc=" <<
					energy_angle << std::endl;
			}

			// accumulate the energy
			emap[ cart_bonded_ring ] += energy_angle;
			emap[ cart_bonded ] += energy_angle; // potential double counting*/

			// 2 constrain torsion
			Real phi0 = rc.nu_angles[ii]*pi/180.0;
			angle = numeric::dihedral_radians(
				rsd.xyz( atm1 ), rsd.xyz( atm2 ), rsd.xyz( atm3 ), rsd.xyz( atm4 ) );
			Real del_phi = basic::subtract_radian_angles(angle, phi0);

			core::Real const energy_torsion = eval_score( del_phi, Kphi, 0 );

			if ( energy_torsion > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " <<
					pose.pdb_info()->number(rsd.seqpos()) << " RING torsion: " <<
					get_restag(rsd.type()) << " : " <<
					rsd.atom_name( atm1 ) << " , " << rsd.atom_name( atm2 ) << " , " <<
					rsd.atom_name( atm3 ) << " , " << rsd.atom_name( atm4 ) << "   (" <<
					Kphi << ") " << 180/pi * angle << " " << 180/pi * phi0 << "    sc="  << energy_torsion << std::endl;
			}

			emap[ cart_bonded_ring ] += energy_torsion;
			emap[ cart_bonded ] += energy_torsion; // potential double counting*/

		}
	}
}


void
CartesianBondedEnergy::eval_singleres_improper_energies(
	conformation::Residue const & rsd,
	ResidueCartBondedParameters const & resparams,
	Real const phi,
	Real const psi,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	using namespace core::chemical;
	using numeric::constants::d::pi;

	const core::Real d_multiplier = core::chemical::is_canonical_D_aa(rsd.aa()) ? -1.0 : 1.0 ; //Multiplier for D-amino acid derivatives

	utility::vector1< ResidueCartBondedParameters::torsion_parameter > const & itps(
		resparams.improper_parameters() );

	for ( Size ii = 1, iiend = itps.size(); ii <= iiend; ++ii ) {
		ResidueCartBondedParameters::Size4 const & atids( itps[ ii ].first );
		CartBondedParameters const & imp_params( *itps[ ii ].second );
		Real Kphi = imp_params.K(phi,psi);
		Real phi0 = d_multiplier * imp_params.mu(phi,psi);
		Real phi_step=2 * pi / imp_params.period();
		Real angle = numeric::dihedral_radians(
			rsd.xyz( atids[1] ), rsd.xyz( atids[2] ), rsd.xyz( atids[3] ), rsd.xyz( atids[4] ) );
		Real del_phi = basic::subtract_radian_angles( angle, phi0 );
		del_phi = basic::periodic_range( del_phi, phi_step );

		core::Real const energy_improper = eval_score( del_phi, Kphi, 0 );

		// Send a message to the user about a bad angle, if necessary.
		// Make sure not to send output to a tracer in the middle of
		// scoring unless that tracer is visible, since that can be very expensive
		if ( energy_improper > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
			TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " <<
				pose.pdb_info()->number(rsd.seqpos()) << " improper torsion: " <<
				get_restag(rsd.type()) << " : " <<
				rsd.atom_name( atids[1] ) << " , " << rsd.atom_name( atids[2] ) << " , " <<
				rsd.atom_name( atids[3] ) << " , " << rsd.atom_name( atids[4] ) << "   (" <<
				Kphi << ") " << angle << " " << phi0 << "    sc="  << energy_improper << std::endl;
		}

		emap[ cart_bonded_improper ] += energy_improper;
		emap[ cart_bonded_torsion ] += energy_improper;
		emap[ cart_bonded ] += energy_improper; // potential double counting*/
	}
}

/// @brief helper function to handle intrares bond torsions
void
CartesianBondedEnergy::eval_singleres_torsion_energies(
	conformation::Residue const & rsd,
	ResidueCartBondedParameters const & resparams,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	using namespace core::chemical;
	using numeric::constants::d::pi;

	const core::Real d_multiplier = core::chemical::is_canonical_D_aa(rsd.aa()) ? -1.0 : 1.0 ; //Multiplier for D-amino acid derivatives

	utility::vector1< ResidueCartBondedParameters::torsion_parameter > const & tps( resparams.torsion_parameters() );

	for ( Size ii = 1, iiend = tps.size(); ii <= iiend; ++ii ) {
		ResidueCartBondedParameters::Size4 const & atids( tps[ ii ].first );
		CartBondedParameters const & tor_params( *tps[ ii ].second );


		Real Kphi = tor_params.K(0,0);
		Real phi0 = d_multiplier * tor_params.mu(0,0);
		Real angle = numeric::dihedral_radians(
			rsd.xyz( atids[1] ), rsd.xyz( atids[2] ), rsd.xyz( atids[3] ), rsd.xyz( atids[4] ) );

		core::Real const energy_torsion = Kphi*(cos(tor_params.period()*angle - phi0) + 1.0);

		// Send a message to the user about a bad angle, if necessary.
		// Make sure not to send output to a tracer in the middle of
		// scoring unless that tracer is visible, since that can be very expensive
		if ( TR.Debug.visible() && pose.pdb_info() ) {
			TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " <<
				pose.pdb_info()->number(rsd.seqpos()) << " intrares torsion: " <<
				get_restag(rsd.type()) << " : " <<
				rsd.atom_name( atids[1] ) << " , " << rsd.atom_name( atids[2] ) << " , " <<
				rsd.atom_name( atids[3] ) << " , " << rsd.atom_name( atids[4] ) << "   (" <<
				Kphi << ") " << angle << " " << phi0 << "    sc="  << energy_torsion << std::endl;
		}

		emap[ cart_bonded_proper ] += energy_torsion;
		emap[ cart_bonded_torsion ] += energy_torsion;
		emap[ cart_bonded ] += energy_torsion; // potential double counting*/
	}
}

/// @brief helper function to handle intrares bond angles
void
CartesianBondedEnergy::eval_singleres_angle_energies(
	conformation::Residue const & rsd,
	ResidueCartBondedParameters const & resparams,
	Real const phi,
	Real const psi,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	using namespace core::chemical;

	utility::vector1< ResidueCartBondedParameters::angle_parameter > const & aps( resparams.angle_parameters() );

	for ( Size ii = 1, iiend = aps.size(); ii <= iiend; ++ii ) {
		ResidueCartBondedParameters::Size3 const & atids( aps[ ii ].first );
		CartBondedParameters const & ang_params( *aps[ ii ].second );

		// ring angle?  Let cart_bonded_ring handle it
		Size const n_rings( rsd.type().n_rings() );
		bool skip_this_angle=false;
		for ( core::uint i( 1 ); i <= n_rings && !skip_this_angle; ++i ) {
			if ( rsd.type().is_ring_atom( i, atids[1] ) && rsd.type().is_ring_atom( i, atids[2] )
					&& rsd.type().is_ring_atom( i, atids[2] ) ) {
				skip_this_angle=true;
			}
		}
		if ( skip_this_angle ) continue;

		Real Ktheta = ang_params.K(phi,psi);
		Real theta0 = ang_params.mu(phi,psi);
		Real angle = numeric::angle_radians( rsd.xyz(atids[1]), rsd.xyz(atids[2]), rsd.xyz(atids[3]) );

		Real const energy_angle = eval_score( angle, Ktheta, theta0 );

		// Send a message to the user about a bad angle, if necessary.
		// Make sure not to send output to a tracer in the middle of
		// scoring unless that tracer is visible, since that can be very expensive
		if ( energy_angle > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
			TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " <<
				pose.pdb_info()->number(rsd.seqpos()) << " intrares angle: " <<
				get_restag(rsd.type()) << " : " <<
				rsd.atom_name( atids[1] ) << " , " << rsd.atom_name( atids[2] ) << " , " <<
				rsd.atom_name( atids[3] ) << "   (" <<
				Ktheta << ") " << angle << " " << theta0 << "    sc=" <<
				energy_angle << std::endl;
		}

		// accumulate the energy
		emap[ cart_bonded_angle ] += energy_angle;
		emap[ cart_bonded ] += energy_angle; // potential double counting*/
	}
}

/// @brief helper function to handle intrares bond lengths
void
CartesianBondedEnergy::eval_singleres_length_energies(
	conformation::Residue const & rsd,
	ResidueCartBondedParameters const & resparams,
	Real const phi,
	Real const psi,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	using namespace core::chemical;

	utility::vector1< ResidueCartBondedParameters::length_parameter > const & lps( resparams.length_parameters() );

	for ( Size ii = 1, iiend = lps.size(); ii <= iiend; ++ii ) {
		ResidueCartBondedParameters::Size2 const & atids( lps[ ii ].first );
		CartBondedParameters const & len_params( *lps[ ii ].second );
		Real const Kd = len_params.K(phi,psi);
		Real const d0 = len_params.mu(phi,psi);
		Real const d = rsd.xyz(atids[1]).distance( rsd.xyz( atids[2] ));

		Real const energy_length = eval_score( d, Kd, d0 );
		// Send a message to the user about a bad angle, if necessary.
		// Make sure not to send output to a tracer in the middle of
		// scoring unless that tracer is visible, since that can be very expensive
		if ( energy_length > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
			TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " << pose.pdb_info()->number(rsd.seqpos()) << " intrares angle: " <<
				get_restag( rsd.type() ) << " : " <<
				rsd.atom_name( atids[1] ) << " , " << rsd.atom_name( atids[2] ) << "   ("
				<< Kd << ")   " << d << "  " << d0
				<< "   sc=" << energy_length << std::endl;
		}

		// accumulate the energy
		emap[ cart_bonded ] += energy_length;
		emap[ cart_bonded_length ] += energy_length;
	}
}

/// @details rsd1.seqpos < rsd2.seqpos, and rsd1&rsd2 should be chemically
/// bonded.
void
CartesianBondedEnergy::eval_residue_pair_energies(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResidueCartBondedParameters const & rsd1params,
	ResidueCartBondedParameters const & rsd2params,
	Real phi1,
	Real psi1,
	Real phi2,
	Real psi2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	eval_interresidue_improper_energy( rsd1, rsd2, rsd1params, rsd2params, pose, emap );
	eval_interresidue_angle_energies_two_from_rsd1( rsd1, rsd2, rsd1params, rsd2params, phi1, psi1, pose, emap );
	eval_interresidue_angle_energies_two_from_rsd2( rsd1, rsd2, rsd1params, rsd2params, phi2, psi2, pose, emap );
	eval_interresidue_bond_energy( rsd1, rsd2, rsd1params, rsd2params, phi2, psi2, pose, emap );
}

void
CartesianBondedEnergy::eval_interresidue_angle_energies_two_from_rsd1(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResidueCartBondedParameters const & rsd1params,
	ResidueCartBondedParameters const & rsd2params,
	Real phi1,
	Real psi1,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	using namespace core::chemical;
	bool rsd1_is_protein = (rsd1.aa() <= num_canonical_aas || core::chemical::is_canonical_D_aa(rsd1.aa()));
	bool rsd2_is_protein = (rsd2.aa() <= num_canonical_aas || core::chemical::is_canonical_D_aa(rsd2.aa()));
	bool rsd12_peptide_bonded = !rsd1.is_upper_terminus() && rsd1.residue_connection_partner( rsd1.upper_connect().index() ) == rsd2.seqpos();
	if ( rsd1_is_protein && rsd2_is_protein && rsd12_peptide_bonded ) {

		utility::vector1< ResidueCartBondedParameters::angle_parameter > const & aps( rsd1params.upper_connect_angle_params() );
		for ( Size ii = 1, iiend = aps.size(); ii <= iiend; ++ii ) {
			ResidueCartBondedParameters::Size3 const & atids( aps[ ii ].first );
			CartBondedParameters const & ang_params( *aps[ ii ].second );
			Real Ktheta = ang_params.K(phi1,psi1);
			Real theta0 = ang_params.mu(phi1,psi1);
			Real angle = numeric::angle_radians(
				rsd1.xyz(atids[1]), rsd1.xyz(atids[2]),
				rsd2.xyz( rsd2params.bb_N_index() ) );

			Real const energy_angle = eval_score( angle, Ktheta, theta0 );
			// Send a message to the user about a bad angle, if necessary.
			// Make sure not to send output to a tracer in the middle of
			// scoring unless that tracer is visible, since that can be very expensive
			if ( energy_angle > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd1.seqpos()
					<< " pdbpos: " << pose.pdb_info()->number(rsd1.seqpos()) << " angle rsd1: " << get_restag(rsd1.type())
					<< ":" << rsd1.atom_name( atids[1] ) << " , "
					<< rsd1.atom_name( atids[2] ) << " , " << rsd2.atom_name( rsd2params.bb_N_index() )<< "   ("
					<< Ktheta << ")   " << angle << "  " << theta0
					<< "   sc=" << energy_angle << std::endl;
			}

			// accumulate the energy
			emap[ cart_bonded_angle ] += energy_angle;
			emap[ cart_bonded ] += energy_angle; // potential double counting*/
		}

	} else {
		// At the time of the writing of this comment, the CartesianBondEnergy is used only
		// to evaluate protein energetics, and the only inter-residue bonds it is evaluated
		// on is peptide bonds.  If that persists, then this "else" clause will never be executed.
		// If it does get executed, then this else clause ought to be optimized.
		// Before the above "if" and the below "else" clauses got separated, the only code
		// executed was the "else" clause and it was rediculously slow.  I mean much much
		// slower than evaluating the etable energies.  REEAAALLLLY SLOW.

		// Update -- 11 July 2013 (V.K. Mulligan):
		// At the time of the writing of THIS comment, the CartesianBondEnergy is used only
		// to evaluate protein energetics OR energetics of D-, L-, or D/L-peptides.  If THAT
		// persists, then this "else" clause will never be executed.  Since it will never be
		// executed, I'm not going to try to figure out whether the comment above is true or
		// not.  It may well be the case that what follows is REEAAALLLLY SLOW.

		utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );
		for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
			Size const resconn_id1( r1_resconn_ids[ii] );
			Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

			Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
			Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

			// compute bond-angle energies with two atoms on rsd1
			utility::vector1< chemical::two_atom_set > const & rsd1_atoms_wi1_bond_of_ii(
				rsd1.type().atoms_within_one_bond_of_a_residue_connection( resconn_id1 ));
			for ( Size jj = 1; jj <= rsd1_atoms_wi1_bond_of_ii.size(); ++jj ) {
				debug_assert( rsd1_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno1 );
				Size const res1_lower_atomno = rsd1_atoms_wi1_bond_of_ii[ jj ].key2();

				Real const angle = numeric::angle_radians(
					rsd1.atom( res1_lower_atomno ).xyz(),
					rsd1.atom( resconn_atomno1 ).xyz(),
					rsd2.atom( resconn_atomno2 ).xyz() );
				std::string atm1name=rsd1.atom_name(res1_lower_atomno); boost::trim(atm1name);
				std::string atm2name=rsd1.atom_name(resconn_atomno1); boost::trim(atm2name);
				std::string atm3name=rsd2.atom_name(resconn_atomno2); boost::trim(atm3name);

				// lookup Ktheta and theta0
				CartBondedParametersCOP ang_params;
				{
#if defined MULTI_THREADED
					utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
					ang_params =
						db_->lookup_angle(rsd1.type(), (rsd2.aa() == chemical::aa_pro || rsd2.aa() == chemical::aa_dpr /*D- or L-proline*/),
						atm1name, atm2name, atm3name, res1_lower_atomno, resconn_atomno1, -resconn_id1);
				}

				if ( ang_params->is_null() ) continue;
				Real Ktheta=ang_params->K(phi1,psi1), theta0=ang_params->mu(phi1,psi1);

				Real const energy_angle = eval_score( angle, Ktheta, theta0 );

				// Send a message to the user about a bad angle, if necessary.
				// Make sure not to send output to a tracer in the middle of
				// scoring unless that tracer is visible, since that can be very expensive
				if ( energy_angle > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
					TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd1.seqpos()
						<< " pdbpos: " << pose.pdb_info()->number(rsd1.seqpos()) << " angle rsd1: " << get_restag(rsd1.type())
						<< ":" << rsd1.atom_name( res1_lower_atomno ) << " , "
						<< rsd1.atom_name( resconn_atomno1 ) << " , " << rsd2.atom_name( resconn_atomno2 )<< "   ("
						<< Ktheta << ")   " << angle << "  " << theta0
						<< "   sc=" << energy_angle << std::endl;
				}

				// accumulate the energy
				emap[ cart_bonded_angle ] += energy_angle;
				emap[ cart_bonded ] += energy_angle; // potential double counting*/
			}
		}
	}
}

void
CartesianBondedEnergy::eval_interresidue_angle_energies_two_from_rsd2(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResidueCartBondedParameters const & rsd1params,
	ResidueCartBondedParameters const & rsd2params,
	Real phi2,
	Real psi2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	using namespace core::chemical;
	bool rsd1_is_protein = (rsd1.aa() <= num_canonical_aas || core::chemical::is_canonical_D_aa(rsd1.aa()));
	bool rsd2_is_protein = (rsd2.aa() <= num_canonical_aas || core::chemical::is_canonical_D_aa(rsd2.aa()));
	bool rsd12_peptide_bonded = !rsd1.is_upper_terminus() && rsd1.residue_connection_partner( rsd1.upper_connect().index() ) == rsd2.seqpos();
	if ( rsd1_is_protein && rsd2_is_protein && rsd12_peptide_bonded ) {

		utility::vector1< ResidueCartBondedParameters::angle_parameter > const & aps( rsd2params.lower_connect_angle_params() );
		for ( Size ii = 1, iiend = aps.size(); ii <= iiend; ++ii ) {
			ResidueCartBondedParameters::Size3 const & atids( aps[ ii ].first );
			CartBondedParameters const & ang_params( *aps[ ii ].second );
			Real Ktheta = ang_params.K(phi2,psi2);
			Real theta0 = ang_params.mu(phi2,psi2);
			Real angle = numeric::angle_radians(
				rsd2.xyz(atids[1]), rsd2.xyz(atids[2]),
				rsd1.xyz( rsd1params.bb_C_index() ) );

			Real const energy_angle = eval_score( angle, Ktheta, theta0 );

			// Send a message to the user about a bad angle, if necessary.
			// Make sure not to send output to a tracer in the middle of
			// scoring unless that tracer is visible, since that can be very expensive
			if ( energy_angle > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd2.seqpos()
					<< " pdbpos: " << pose.pdb_info()->number(rsd2.seqpos()) << " angle rsd2: " << get_restag(rsd2.type())
					<< ":" << rsd2.atom_name( atids[1] ) << " , "
					<< rsd2.atom_name( atids[2] ) << " , " << rsd1.atom_name( rsd1params.bb_C_index()  )
					<< "   (" << Ktheta << ")   " << angle << "  " << theta0
					<< "   sc=" << energy_angle << std::endl;
			}

			// accumulate the energy
			emap[ cart_bonded_angle ] += energy_angle;
			emap[ cart_bonded ] += energy_angle; // potential double counting*/
		}

	} else {
		// At the time of the writing of this comment, the CartesianBondEnergy is used only
		// to evaluate protein energetics, and the only inter-residue bonds it is evaluated
		// on is peptide bonds.  If that persists, then this "else" clause will never be executed.
		// If it does get executed, then this else clause ought to be optimized.
		// Before the above "if" and the below "else" clauses got separated, the only code
		// executed was the "else" clause and it was rediculously slow.  I mean much much
		// slower than evaluating the etable energies.  REEAAALLLLY SLOW.

		// Update -- 11 July 2013 (V.K. Mulligan):
		// At the time of the writing of THIS comment, the CartesianBondEnergy is used only
		// to evaluate protein energetics OR energetics of D-, L-, or D/L-peptides.  If THAT
		// persists, then this "else" clause will never be executed.  Since it will never be
		// executed, I'm not going to try to figure out whether the comment above is true or
		// not.  It may well be the case that what follows is REEAAALLLLY SLOW.

		utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );
		for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
			Size const resconn_id1( r1_resconn_ids[ii] );
			Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

			Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
			Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

			/// compute bond-angle energies with two atoms on rsd2
			utility::vector1< chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
				rsd2.type().atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));
			for ( Size jj = 1; jj <= rsd2_atoms_wi1_bond_of_ii.size(); ++jj ) {
				debug_assert( rsd2_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno2 );
				Size const res2_lower_atomno = rsd2_atoms_wi1_bond_of_ii[ jj ].key2();

				// lookup Ktheta and theta0
				// set pre-proline to false here -- it refers to rsd2+1 which would make this three-body
				std::string atm1name=rsd2.atom_name(res2_lower_atomno); boost::trim(atm1name);
				std::string atm2name=rsd2.atom_name(resconn_atomno2); boost::trim(atm2name);
				std::string atm3name=rsd1.atom_name(resconn_atomno1); boost::trim(atm3name);

				CartBondedParametersCOP ang_params;
				{
#if defined MULTI_THREADED
					utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
					ang_params =
						db_->lookup_angle(rsd2.type(), false,
						atm1name, atm2name, atm3name, res2_lower_atomno, resconn_atomno2, -resconn_id2);
				}

				if ( ang_params->is_null() ) continue;
				Real Ktheta=ang_params->K(phi2,psi2), theta0=ang_params->mu(phi2,psi2);

				if ( Ktheta == 0.0 ) continue;
				Real const angle = numeric::angle_radians(
					rsd2.atom( res2_lower_atomno ).xyz(),
					rsd2.atom( resconn_atomno2 ).xyz(),
					rsd1.atom( resconn_atomno1 ).xyz() );

				Real const energy_angle = eval_score( angle, Ktheta, theta0 );

				// Send a message to the user about a bad angle, if necessary.
				// Make sure not to send output to a tracer in the middle of
				// scoring unless that tracer is visible, since that can be very expensive
				if ( energy_angle > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
					TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd2.seqpos()
						<< " pdbpos: " << pose.pdb_info()->number(rsd2.seqpos()) << " angle rsd2: " << get_restag(rsd2.type())
						<< ":" << rsd2.atom_name( res2_lower_atomno ) << " , "
						<< rsd2.atom_name( resconn_atomno2 ) << " , " << rsd1.atom_name( resconn_atomno1 )
						<< "   (" << Ktheta << ")   " << angle << "  " << theta0
						<< "   sc=" << energy_angle << std::endl;
				}

				// accumulate the energy
				emap[ cart_bonded_angle ] += energy_angle;
				emap[ cart_bonded ] += energy_angle; // potential double counting*/

			}
		}
	}
}

void
CartesianBondedEnergy::eval_interresidue_bond_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResidueCartBondedParameters const & rsd1params,
	ResidueCartBondedParameters const & rsd2params,
	Real phi2,
	Real psi2,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	using namespace core::chemical;
	bool rsd1_is_protein = (rsd1.aa() <= num_canonical_aas || core::chemical::is_canonical_D_aa(rsd1.aa()));
	bool rsd2_is_protein = (rsd2.aa() <= num_canonical_aas || core::chemical::is_canonical_D_aa(rsd2.aa()));
	bool rsd12_peptide_bonded = !rsd1.is_upper_terminus() && rsd1.residue_connection_partner( rsd1.upper_connect().index() ) == rsd2.seqpos();
	if ( rsd1_is_protein && rsd2_is_protein && rsd12_peptide_bonded ) {

		/////////////
		/// finally, compute the bondlength across the interface
		Real length = rsd1.xyz( rsd1params.bb_C_index() ).distance( rsd2.xyz( rsd2params.bb_N_index() ) );

		// lookup Kd and d0
		CartBondedParametersCOP len_params = rsd2params.cprev_n_bond_length_params();

		if ( len_params->is_null() ) return;
		Real Kd=len_params->K(phi2,psi2), d0=len_params->mu(phi2,psi2);

		Real const energy_length = eval_score( length, Kd, d0 );

		// Send a message to the user about a bad distance, if necessary.
		// Make sure not to send output to a tracer in the middle of
		// scoring unless that tracer is visible, since that can be very expensive
		if ( energy_length > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
			TR.Debug << pose.pdb_info()->name()
				<< " pdbpos rsd1: " << pose.pdb_info()->number(rsd1.seqpos()) << " length rsd1 rsd2: " << rsd1.seqpos() << " -- " << rsd2.seqpos() << "  "
				<< get_restag(rsd1.type()) << ":" << rsd1.atom_name( rsd1params.bb_C_index() ) << " , " << rsd2.atom_name( rsd2params.bb_N_index() )<< "   ("
				<< Kd << ")   " << length << "  " << d0
				<< "   sc=" << energy_length << std::endl;
		}

		// accumulate the energy
		emap[ cart_bonded ] += energy_length;
		emap[ cart_bonded_length ] += energy_length;

	} else {
		// At the time of the writing of this comment, the CartesianBondEnergy is used only
		// to evaluate protein energetics, and the only inter-residue bonds it is evaluated
		// on is peptide bonds.  If that persists, then this "else" clause will never be executed.
		// If it does get executed, then this else clause ought to be optimized.
		// Before the above "if" and the below "else" clauses got separated, the only code
		// executed was the "else" clause and it was rediculously slow.  I mean much much
		// slower than evaluating the etable energies.  REEAAALLLLY SLOW.

		// Update -- 11 July 2013 (V.K. Mulligan):
		// At the time of the writing of THIS comment, the CartesianBondEnergy is used only
		// to evaluate protein energetics OR energetics of D-, L-, or D/L-peptides.  If THAT
		// persists, then this "else" clause will never be executed.  Since it will never be
		// executed, I'm not going to try to figure out whether the comment above is true or
		// not.  It may well be the case that what follows is REEAAALLLLY SLOW.

		utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );
		for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
			Size const resconn_id1( r1_resconn_ids[ii] );
			Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

			Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
			Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

			/////////////
			/// finally, compute the bondlength across the interface
			Real length =
				( rsd2.atom( resconn_atomno2 ).xyz() - rsd1.atom( resconn_atomno1 ).xyz() ).length();

			// lookup Kd and d0
			// again, pre-pro == false, definitions are based on pro as rsd2+1
			std::string atm1name=rsd2.atom_name(resconn_atomno2); boost::trim(atm1name);
			std::string atm2name=rsd1.atom_name(resconn_atomno1); boost::trim(atm2name);

			CartBondedParametersCOP len_params;
			{
#if defined MULTI_THREADED
				utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
				len_params = db_->lookup_length(rsd2.type(), false, atm1name, atm2name, resconn_atomno2, -resconn_id2);
			}

			if ( len_params->is_null() ) continue;
			Real Kd=len_params->K(phi2,psi2), d0=len_params->mu(phi2,psi2);

			Real const energy_length = eval_score( length, Kd, d0 );

			// Send a message to the user about a bad distance, if necessary.
			// Make sure not to send output to a tracer in the middle of
			// scoring unless that tracer is visible, since that can be very expensive
			if ( energy_length > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
				TR.Debug << pose.pdb_info()->name()
					<< " pdbpos rsd1: " << pose.pdb_info()->number(rsd1.seqpos()) << " length rsd1 rsd2: " << rsd1.seqpos() << " -- " << rsd2.seqpos() << "  "
					<< get_restag(rsd1.type()) << ":" << rsd1.atom_name( resconn_atomno1 ) << " , " << rsd2.atom_name( resconn_atomno2 )<< "   ("
					<< Kd << ")   " << length << "  " << d0
					<< "   sc=" << energy_length << std::endl;
			}

			// accumulate the energy
			emap[ cart_bonded ] += energy_length;
			emap[ cart_bonded_length ] += energy_length;
		}
	}

}

void
CartesianBondedEnergy::eval_interresidue_improper_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResidueCartBondedParameters const & rsd1params,
	ResidueCartBondedParameters const & rsd2params,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	using namespace core::chemical;
	using numeric::constants::d::pi;

	if ( !rsd1.is_protein() || !rsd2.is_protein() || !rsd1.is_polymer_bonded( rsd2 ) ) return;

	// exit if this is not a backbone connection
	if ( rsd1.is_upper_terminus() || rsd1.residue_connection_partner( rsd1.upper_connect().index() ) != rsd2.seqpos() ) return;

	const core::Real d_multiplier2 = core::chemical::is_canonical_D_aa(rsd2.aa()) ? -1.0 : 1.0 ; //Multiplier for D-amino acid derivatives

	// note: this function is hardcoded for the particular planar groups crossing residue connections
	//       a more general scheme would be quite nice

	// backbone CA-Cprev-N-H
	if ( (rsd2.aa() != aa_pro && rsd2.aa() != aa_dpr /*Not D- or L-proline*/) && rsd2params.bb_H_index() != 0 ) {
		CartBondedParametersCOP tor_params = rsd2params.ca_cprev_n_h_interres_improper_params();
		if ( !tor_params->is_null() ) {
			Real const Kphi = tor_params->K(0,0);
			Real const phi0 = d_multiplier2 * tor_params->mu(0,0);
			Real const phi_step = 2*pi/tor_params->period();
			Real angle = numeric::dihedral_radians(
				rsd2.xyz( rsd2params.bb_CA_index() ),
				rsd1.xyz( rsd1params.bb_C_index() ),
				rsd2.xyz( rsd2params.bb_N_index() ),
				rsd2.xyz( rsd2params.bb_H_index() )
			);
			Real del_phi = basic::subtract_radian_angles(angle, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );

			Real energy_improper = eval_score( del_phi, Kphi, 0 );

			// Send a message to the user about a bad angle, if necessary.
			// Make sure not to send output to a tracer in the middle of
			// scoring unless that tracer is visible, since that can be very expensive
			if ( energy_improper > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd1.seqpos() << " pdbpos: " <<
					pose.pdb_info()->number(rsd1.seqpos()) << " improper torsion: " <<
					get_restag(rsd1.type()) << " : " <<
					rsd2.atom_name( rsd2params.bb_CA_index() ) << " , " << rsd1.atom_name( rsd1params.bb_C_index() ) << " , " <<
					rsd2.atom_name( rsd2params.bb_N_index() )  << " , " << rsd2.atom_name( rsd2params.bb_H_index() ) << "   (" <<
					Kphi << ") " << angle << " " << phi0 << "    sc=" << energy_improper << std::endl;
			}

			emap[ cart_bonded ] += energy_improper;
			emap[ cart_bonded_torsion ] += energy_improper;
			emap[ cart_bonded_improper ] += energy_improper;
		}

	}

	// backbone Oprev-Cprev-N-H
	if ( (rsd2.aa() != aa_pro && rsd2.aa() != aa_dpr /*Not D- or L-proline*/) && rsd2params.bb_H_index() != 0
			&& rsd1params.bb_O_index() != 0 ) {
		CartBondedParametersCOP tor_params = rsd2params.oprev_cprev_n_h_interres_improper_params();
		if ( !tor_params->is_null() ) {
			Real const Kphi = tor_params->K(0,0);
			Real const phi0 = d_multiplier2 * tor_params->mu(0,0);
			Real const phi_step = 2*pi/tor_params->period();
			Real angle = numeric::dihedral_radians(
				rsd1.xyz( rsd1params.bb_O_index() ),
				rsd1.xyz( rsd1params.bb_C_index() ),
				rsd2.xyz( rsd2params.bb_N_index() ),
				rsd2.xyz( rsd2params.bb_H_index() )
			);
			Real del_phi = basic::subtract_radian_angles(angle, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );

			Real energy_improper = eval_score( del_phi, Kphi, 0 );

			// Send a message to the user about a bad angle, if necessary.
			// Make sure not to send output to a tracer in the middle of
			// scoring unless that tracer is visible, since that can be very expensive
			if ( energy_improper > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd1.seqpos() << " pdbpos: " <<
					pose.pdb_info()->number(rsd1.seqpos()) << " improper torsion: " <<
					get_restag(rsd1.type()) << " : " <<
					rsd1.atom_name( rsd1params.bb_O_index() ) << " , " << rsd1.atom_name( rsd1params.bb_C_index() ) << " , " <<
					rsd2.atom_name( rsd2params.bb_N_index() )  << " , " << rsd2.atom_name( rsd2params.bb_H_index() ) << "   (" <<
					Kphi << ") " << angle << " " << phi0 << "    sc=" << energy_improper << std::endl;
			}

			emap[ cart_bonded ] += energy_improper;
			emap[ cart_bonded_torsion ] += energy_improper;
			emap[ cart_bonded_improper ] += energy_improper;
		}

	}


	// backbone CA-Nnext-C-O
	{
		CartBondedParametersCOP tor_params = rsd1params.ca_nnext_c_o_interres_improper_params();
		if ( !tor_params->is_null() ) {
			Real const Kphi = tor_params->K(0,0);
			Real const phi0 = d_multiplier2 * tor_params->mu(0,0);
			Real const phi_step = 2*pi/tor_params->period();
			Real angle = numeric::dihedral_radians(
				rsd1.xyz( rsd1params.bb_O_index() ),
				rsd1.xyz( rsd1params.bb_C_index() ),
				rsd2.xyz( rsd2params.bb_N_index() ),
				rsd1.xyz( rsd1params.bb_CA_index() )
			);
			Real del_phi = basic::subtract_radian_angles(angle, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );

			Real energy_improper = eval_score( del_phi, Kphi, 0 );

			// Send a message to the user about a bad angle, if necessary.
			// Make sure not to send output to a tracer in the middle of
			// scoring unless that tracer is visible, since that can be very expensive
			if ( energy_improper > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd1.seqpos() << " pdbpos: " <<
					pose.pdb_info()->number(rsd1.seqpos()) << " improper torsion: " <<
					get_restag(rsd1.type()) << " : " <<
					rsd1.atom_name( rsd1params.bb_O_index() ) << " , " << rsd1.atom_name( rsd1params.bb_C_index() ) << " , " <<
					rsd2.atom_name( rsd2params.bb_N_index() ) << " , " << rsd1.atom_name( rsd1params.bb_CA_index() ) << "   (" <<
					Kphi << ") " << angle << " " << phi0 << "    sc=" << energy_improper << std::endl;
			}

			emap[ cart_bonded ] += energy_improper;
			emap[ cart_bonded_torsion ] += energy_improper;
			emap[ cart_bonded_improper ] += energy_improper;
		}
	}

	// proline N planarity
	if ( rsd2.aa() == core::chemical::aa_pro || rsd2.aa() == core::chemical::aa_dpr ) { //D- or L-proline
		CartBondedParametersCOP tor_params = rsd2params.pro_cd_cprev_n_ca_interres_improper_params();
		if ( tor_params && !tor_params->is_null() ) {
			Real const Kphi = tor_params->K(0,0);
			Real const phi0 = d_multiplier2 * tor_params->mu(0,0);
			Real const phi_step=2*pi/tor_params->period();
			Real angle = numeric::dihedral_radians(
				rsd2.xyz( rsd2params.pro_CD_index() ),
				rsd1.xyz( rsd1params.bb_C_index() ),
				rsd2.xyz( rsd2params.bb_N_index() ),
				rsd2.xyz( rsd2params.bb_CA_index() )
			);
			Real del_phi = basic::subtract_radian_angles(angle, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );

			Real energy_improper = eval_score( del_phi, Kphi, 0 );

			// Send a message to the user about a bad angle, if necessary.
			// Make sure not to send output to a tracer in the middle of
			// scoring unless that tracer is visible, since that can be very expensive
			if ( energy_improper > CUTOFF && TR.Debug.visible() && pose.pdb_info() ) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd1.seqpos() << " pdbpos: " <<
					pose.pdb_info()->number(rsd1.seqpos()) << " improper torsion: " <<
					get_restag(rsd1.type()) << " : " <<
					rsd1.atom_name( rsd1params.bb_O_index() ) << " , " << rsd1.atom_name( rsd1params.bb_C_index() ) << " , " <<
					rsd2.atom_name( rsd2params.bb_N_index() ) << " , " << rsd1.atom_name( rsd1params.bb_CA_index() ) << "   (" <<
					Kphi << ") " << angle << " " << phi0 << "    sc=" << energy_improper << std::endl;
			}

			emap[ cart_bonded ] += energy_improper;
			emap[ cart_bonded_torsion ] += energy_improper;
			emap[ cart_bonded_improper ] += energy_improper;
		}
	}
}


void
CartesianBondedEnergy::eval_singleres_derivatives(
	conformation::Residue const & rsd,
	ResidueCartBondedParameters const & resparams,
	Real phi,
	Real psi,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r_atom_derivs
) const {
	using numeric::constants::d::pi;
	using namespace numeric;

	chemical::ResidueType const & rsd_type = rsd.type();
	if ( rsd_type.aa() == core::chemical::aa_vrt ) return;

	eval_singleres_torsion_derivatives( rsd, resparams, weights, r_atom_derivs );
	eval_singleres_improper_derivatives( rsd, resparams, phi, psi, weights, r_atom_derivs );
	eval_singleres_angle_derivatives(   rsd, resparams, phi, psi, weights, r_atom_derivs );
	eval_singleres_length_derivatives(  rsd, resparams, phi, psi, weights, r_atom_derivs );
	eval_singleres_ring_derivatives(  rsd, weights, r_atom_derivs );//BF. note removed resparams as unused variable.

}

/// @brief helper function to handle intrares bond torsions
void
CartesianBondedEnergy::eval_singleres_ring_derivatives(
	conformation::Residue const & rsd,
	//ResidueCartBondedParameters const & resparams,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r_atom_derivs
) const
{
	using namespace core::chemical;
	using numeric::constants::d::pi;

	Size const n_rings( rsd.type().n_rings() );
	if ( n_rings==0 ) return;

	core::Real Ktheta, Kphi;
	{
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
		Ktheta = db_->k_torsion();  // for now use default torsion (TO DO: add spring constants to DB!)
		Kphi = db_->k_angle();  // for now use default torsion (TO DO: add spring constants to DB!)
	}

	Real const weight = weights[ cart_bonded_ring ] + weights[ cart_bonded ];

	for ( core::uint jj( 1 ); jj <= n_rings; ++jj ) {
		// get the conformer of the ring
		core::chemical::rings::RingConformer rc = rsd.ring_conformer( jj, 180.0 );

		// now constrain each element of the ring
		utility::vector1< core::Size > atms = rsd.type().ring_atoms( jj );

		for ( core::Size ii=1; ii<=atms.size(); ++ii ) {
			core::Size rt1 = atms[(ii+atms.size()-2)%atms.size()+1];
			core::Size rt2 = atms[ii];
			core::Size rt3 = atms[ii%atms.size()+1];
			core::Size rt4 = atms[(ii+1)%atms.size()+1];

			// 1  angle
			Real theta0 = rc.tau_angles[ii]*pi/180.0;

			Vector f1(0.0), f2(0.0);
			Real theta(0.0), dE_dtheta;

			numeric::deriv::angle_p1_deriv( rsd.xyz(rt1), rsd.xyz(rt2), rsd.xyz(rt3), theta, f1, f2 );
			if ( linear_bonded_potential_ && std::fabs(theta - theta0)>1 ) {
				dE_dtheta = weight * Ktheta * ((theta > theta0)>0? 0.5 : -0.5);
			} else {
				dE_dtheta = weight * Ktheta * (theta - theta0);
			}

			r_atom_derivs[ rt1 ].f1() += dE_dtheta * f1;
			r_atom_derivs[ rt1 ].f2() += dE_dtheta * f2;

			numeric::deriv::angle_p2_deriv( rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), theta, f1, f2 );
			r_atom_derivs[ rt2 ].f1() += dE_dtheta * f1;
			r_atom_derivs[ rt2 ].f2() += dE_dtheta * f2;

			numeric::deriv::angle_p1_deriv( rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), theta, f1, f2 );
			r_atom_derivs[ rt3 ].f1() += dE_dtheta * f1;
			r_atom_derivs[ rt3 ].f2() += dE_dtheta * f2;

			// 2  torsion
			Real phi0 = rc.nu_angles[ii]*pi/180.0;

			f1 = f2 = Vector(0.0);
			Real phi=0, dE_dphi;

			numeric::deriv::dihedral_p1_cosine_deriv(
				rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), rsd.xyz( rt4 ), phi, f1, f2 );
			Real del_phi = basic::subtract_radian_angles(phi, phi0);
			if ( linear_bonded_potential_ && std::fabs(del_phi)>1 ) {
				dE_dphi = weight * Kphi * (del_phi>0? 0.5 : -0.5);
			} else {
				dE_dphi = weight * Kphi * del_phi;
			}

			r_atom_derivs[ rt1 ].f1() += dE_dphi * f1;
			r_atom_derivs[ rt1 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), rsd.xyz( rt4 ), phi, f1, f2 );
			r_atom_derivs[ rt2 ].f1() += dE_dphi * f1;
			r_atom_derivs[ rt2 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				rsd.xyz( rt4 ), rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), phi, f1, f2 );

			r_atom_derivs[ rt3 ].f1() += dE_dphi * f1;
			r_atom_derivs[ rt3 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv(
				rsd.xyz( rt4 ), rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), phi, f1, f2 );
			r_atom_derivs[ rt4 ].f1() += dE_dphi * f1;
			r_atom_derivs[ rt4 ].f2() += dE_dphi * f2;
		}
	}
}


/// @brief helper function to handle intrares bond angles
void
CartesianBondedEnergy::eval_singleres_angle_derivatives(
	conformation::Residue const & rsd,
	ResidueCartBondedParameters const & resparams,
	Real const phi,
	Real const psi,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r_atom_derivs
) const
{
	using namespace core::chemical;

	utility::vector1< ResidueCartBondedParameters::angle_parameter > const & aps( resparams.angle_parameters() );
	Real const weight = weights[ cart_bonded_angle ] + weights[ cart_bonded ];

	for ( Size ii = 1, iiend = aps.size(); ii <= iiend; ++ii ) {
		ResidueCartBondedParameters::Size3 const & atids( aps[ ii ].first );
		Size const rt1( atids[1] ), rt2( atids[2] ), rt3( atids[3] );
		CartBondedParameters const & ang_params( *aps[ ii ].second );
		if ( ang_params.is_null() ) continue;

		// ring angle?  Let cart_bonded_ring handle it
		Size const n_rings( rsd.type().n_rings() );
		bool skip_this_angle=false;
		for ( core::uint i( 1 ); i <= n_rings && !skip_this_angle; ++i ) {
			if ( rsd.type().is_ring_atom( i, atids[1] ) && rsd.type().is_ring_atom( i, atids[2] )
					&& rsd.type().is_ring_atom( i, atids[2] ) ) {
				skip_this_angle=true;
			}
		}
		if ( skip_this_angle ) continue;

		Real Ktheta = ang_params.K(phi,psi);
		Real theta0 = ang_params.mu(phi,psi);

		Vector f1(0.0), f2(0.0);
		Real theta(0.0), dE_dtheta;

		numeric::deriv::angle_p1_deriv( rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), theta, f1, f2 );
		if ( linear_bonded_potential_ && std::fabs(theta - theta0)>1 ) {
			dE_dtheta = weight * Ktheta * ((theta > theta0)>0? 0.5 : -0.5);
		} else {
			dE_dtheta = weight * Ktheta * (theta - theta0);
		}

		r_atom_derivs[ rt1 ].f1() += dE_dtheta * f1;
		r_atom_derivs[ rt1 ].f2() += dE_dtheta * f2;

		numeric::deriv::angle_p2_deriv( rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), theta, f1, f2 );
		r_atom_derivs[ rt2 ].f1() += dE_dtheta * f1;
		r_atom_derivs[ rt2 ].f2() += dE_dtheta * f2;

		numeric::deriv::angle_p1_deriv( rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), theta, f1, f2 );
		r_atom_derivs[ rt3 ].f1() += dE_dtheta * f1;
		r_atom_derivs[ rt3 ].f2() += dE_dtheta * f2;
	}
}

/// @brief helper function to handle intrares bond lengths
void
CartesianBondedEnergy::eval_singleres_length_derivatives(
	conformation::Residue const & rsd,
	ResidueCartBondedParameters const & resparams,
	Real const phi,
	Real const psi,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r_atom_derivs
) const
{
	using namespace core::chemical;

	utility::vector1< ResidueCartBondedParameters::length_parameter > const & lps( resparams.length_parameters() );
	Real const weight = weights[ cart_bonded_length ] + weights[ cart_bonded ];

	for ( Size ii = 1, iiend = lps.size(); ii <= iiend; ++ii ) {
		ResidueCartBondedParameters::Size2 const & atids( lps[ ii ].first );
		Size const rt1( atids[1] ), rt2( atids[2] );
		CartBondedParameters const & len_params( *lps[ ii ].second );
		if ( len_params.is_null() ) continue;

		Real const Kd = len_params.K(phi,psi);
		Real const d0 = len_params.mu(phi,psi);

		Vector f1(0.0), f2(0.0);
		Real d=0, dE_dd;

		numeric::deriv::distance_f1_f2_deriv( rsd.xyz( rt1 ), rsd.xyz( rt2 ), d, f1, f2 );
		if ( linear_bonded_potential_ && std::fabs(d - d0)>1 ) {
			dE_dd = weight * Kd * ((d - d0)>0? 0.5 : -0.5);
		} else {
			dE_dd = weight * Kd * (d - d0);
		}

		r_atom_derivs[ rt1 ].f1() += dE_dd * f1;
		r_atom_derivs[ rt1 ].f2() += dE_dd * f2;

		r_atom_derivs[ rt2 ].f1() -= dE_dd * f1;
		r_atom_derivs[ rt2 ].f2() -= dE_dd * f2;

	}
}

/// @brief helper function to handle intrares bond torsions
void
CartesianBondedEnergy::eval_singleres_improper_derivatives(
	conformation::Residue const & rsd,
	ResidueCartBondedParameters const & resparams,
	Real const phi,
	Real const psi,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r_atom_derivs
) const
{
	using namespace core::chemical;
	using numeric::constants::d::pi;

	const core::Real d_multiplier = core::chemical::is_canonical_D_aa(rsd.aa()) ? -1.0 : 1.0 ; //Multiplier for D-amino acid derivatives

	utility::vector1< ResidueCartBondedParameters::torsion_parameter > const & tps( resparams.improper_parameters() );

	Real const weight = weights[ cart_bonded_improper ] + weights[ cart_bonded_torsion ] + weights[ cart_bonded ];

	for ( Size ii = 1, iiend = tps.size(); ii <= iiend; ++ii ) {
		ResidueCartBondedParameters::Size4 const & atids( tps[ ii ].first );
		Size const rt1( atids[1] ), rt2( atids[2] ), rt3( atids[3] ), rt4( atids[4] );
		CartBondedParameters const & imp_params( *tps[ ii ].second );

		Real Kphi = imp_params.K(phi,psi);
		Real phi0 = d_multiplier * imp_params.mu(phi,psi);
		Real phi_step = 2 * pi / imp_params.period();

		Vector f1(0.0), f2(0.0);
		Real angle(0.0), dE_dphi;

		numeric::deriv::dihedral_p1_cosine_deriv( rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), rsd.xyz( rt4 ), angle, f1, f2 );
		Real del_phi = basic::subtract_radian_angles(angle, phi0);
		if ( phi_step>0 ) del_phi = basic::periodic_range( del_phi, phi_step );
		if ( linear_bonded_potential_ && std::fabs(del_phi)>1 ) {
			dE_dphi = weight * Kphi * (del_phi>0? 0.5 : -0.5);
		} else {
			dE_dphi = weight * Kphi * del_phi;
		}
		r_atom_derivs[ rt1 ].f1() += dE_dphi * f1;
		r_atom_derivs[ rt1 ].f2() += dE_dphi * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv( rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), rsd.xyz( rt4 ), angle, f1, f2 );
		r_atom_derivs[ rt2 ].f1() += dE_dphi * f1;
		r_atom_derivs[ rt2 ].f2() += dE_dphi * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv( rsd.xyz( rt4 ), rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), angle, f1, f2 );

		r_atom_derivs[ rt3 ].f1() += dE_dphi * f1;
		r_atom_derivs[ rt3 ].f2() += dE_dphi * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p1_cosine_deriv( rsd.xyz( rt4 ), rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), angle, f1, f2 );
		r_atom_derivs[ rt4 ].f1() += dE_dphi * f1;
		r_atom_derivs[ rt4 ].f2() += dE_dphi * f2;

	}
}

// deriv impropers
void
CartesianBondedEnergy::eval_singleres_torsion_derivatives(
	conformation::Residue const & rsd,
	ResidueCartBondedParameters const & resparams,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r_atom_derivs
) const {
	using namespace core::chemical;
	using numeric::constants::d::pi;

	const core::Real d_multiplier = core::chemical::is_canonical_D_aa(rsd.aa()) ? -1.0 : 1.0 ; //Multiplier for D-amino acid derivatives

	utility::vector1< ResidueCartBondedParameters::torsion_parameter > const & itps(
		resparams.torsion_parameters() );
	Real const weight = weights[ cart_bonded_proper ] + weights[ cart_bonded_torsion ] + weights[ cart_bonded ];

	for ( Size ii = 1, iiend = itps.size(); ii <= iiend; ++ii ) {
		ResidueCartBondedParameters::Size4 const & atids( itps[ ii ].first );
		Size const rt1( atids[1] ), rt2( atids[2] ), rt3( atids[3] ), rt4( atids[4] );
		CartBondedParameters const & tor_params( *itps[ ii ].second );

		// ring tors?  Let cart_bonded_ring handle it
		Size const n_rings( rsd.type().n_rings() );
		bool skip_this_torsion=false;
		for ( core::uint i( 1 ); i <= n_rings && !skip_this_torsion; ++i ) {
			if ( rsd.type().is_ring_atom( i, atids[1] ) && rsd.type().is_ring_atom( i, atids[2] )
					&& rsd.type().is_ring_atom( i, atids[2] ) && rsd.type().is_ring_atom( i, atids[3] ) ) {
				skip_this_torsion=true;
			}
		}
		if ( skip_this_torsion ) continue;

		Real Kphi = tor_params.K(0,0);
		Real phi0 = d_multiplier * tor_params.mu(0,0);

		Vector f1(0.0), f2(0.0);
		Real phi = numeric::dihedral_radians(
			rsd.xyz( atids[1] ), rsd.xyz( atids[2] ), rsd.xyz( atids[3] ), rsd.xyz( atids[4] ) );

		core::Real const dE_dphi = -Kphi*weight*tor_params.period()*sin(tor_params.period()*phi - phi0);

		numeric::deriv::dihedral_p1_cosine_deriv(
			rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), rsd.xyz( rt4 ), phi, f1, f2 );
		r_atom_derivs[ rt1 ].f1() += dE_dphi * f1;
		r_atom_derivs[ rt1 ].f2() += dE_dphi * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv(
			rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), rsd.xyz( rt4 ), phi, f1, f2 );
		r_atom_derivs[ rt2 ].f1() += dE_dphi * f1;
		r_atom_derivs[ rt2 ].f2() += dE_dphi * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv(
			rsd.xyz( rt4 ), rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), phi, f1, f2 );

		r_atom_derivs[ rt3 ].f1() += dE_dphi * f1;
		r_atom_derivs[ rt3 ].f2() += dE_dphi * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p1_cosine_deriv(
			rsd.xyz( rt4 ), rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), phi, f1, f2 );
		r_atom_derivs[ rt4 ].f1() += dE_dphi * f1;
		r_atom_derivs[ rt4 ].f2() += dE_dphi * f2;
	}
}


void
CartesianBondedEnergy::eval_interresidue_angle_derivs_two_from_rsd1(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResidueCartBondedParameters const & rsd1params,
	ResidueCartBondedParameters const & rsd2params,
	Real phi1,
	Real psi1,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	using namespace core::chemical;
	bool rsd1_is_protein = (rsd1.aa() <= num_canonical_aas || core::chemical::is_canonical_D_aa(rsd1.aa()));
	bool rsd2_is_protein = (rsd2.aa() <= num_canonical_aas || core::chemical::is_canonical_D_aa(rsd2.aa()));
	bool rsd12_peptide_bonded = !rsd1.is_upper_terminus() && rsd1.residue_connection_partner( rsd1.upper_connect().index() ) == rsd2.seqpos();
	if ( rsd1_is_protein && rsd2_is_protein && rsd12_peptide_bonded ) {

		/// Assumption: rsd1 and rsd2 share a peptide bond and only a peptide bond.
		// amw: can we simply remove this bad assumption (e.g. for ncbb/ ncccs?
		//debug_assert( rsd2.residue_connection_partner( rsd2.lower_connect().index() ) == rsd1.seqpos() );
		// amw: for triazolamers let's try just that
		//debug_assert( rsd1.connections_to_residue( rsd2 ).size() == 1 );

		utility::vector1< ResidueCartBondedParameters::angle_parameter > const & aps( rsd1params.upper_connect_angle_params() );
		for ( Size ii = 1, iiend = aps.size(); ii <= iiend; ++ii ) {
			ResidueCartBondedParameters::Size3 const & atids( aps[ ii ].first );
			CartBondedParameters const & ang_params( *aps[ ii ].second );
			Real Ktheta = ang_params.K(phi1,psi1);
			Real theta0 = ang_params.mu(phi1,psi1);

			// temp -- rename these variables, ok?
			Size const res1_lower_atomno = atids[1], resconn_atomno1 = atids[2], resconn_atomno2 = rsd2params.bb_N_index();

			Vector f1(0.0), f2(0.0);
			Real theta(0.0), dE_dtheta;
			numeric::deriv::angle_p1_deriv(
				rsd1.xyz( res1_lower_atomno ), rsd1.xyz( resconn_atomno1 ), rsd2.xyz( resconn_atomno2 ), theta, f1, f2 );
			if ( linear_bonded_potential_ && std::fabs(theta - theta0)>1 ) {
				dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * ((theta - theta0)>0? 0.5 : -0.5);
			} else {
				dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * (theta - theta0);
			}
			r1_atom_derivs[ res1_lower_atomno ].f1() += dE_dtheta * f1;
			r1_atom_derivs[ res1_lower_atomno ].f2() += dE_dtheta * f2;

			numeric::deriv::angle_p2_deriv(
				rsd1.xyz( res1_lower_atomno ), rsd1.xyz( resconn_atomno1 ), rsd2.xyz( resconn_atomno2 ), theta, f1, f2 );
			r1_atom_derivs[ resconn_atomno1 ].f1() += dE_dtheta * f1;
			r1_atom_derivs[ resconn_atomno1 ].f2() += dE_dtheta * f2;

			numeric::deriv::angle_p1_deriv(
				rsd2.xyz( resconn_atomno2 ), rsd1.xyz( resconn_atomno1 ), rsd1.xyz( res1_lower_atomno ), theta, f1, f2 );
			r2_atom_derivs[ resconn_atomno2 ].f1() += dE_dtheta * f1;
			r2_atom_derivs[ resconn_atomno2 ].f2() += dE_dtheta * f2;
		}

	} else {
		// At the time of the writing of this comment, the CartesianBondEnergy is used only
		// to evaluate protein energetics, and the only inter-residue bonds it is evaluated
		// on is peptide bonds.  If that persists, then this "else" clause will never be executed.
		// If it does get executed, then this else clause ought to be optimized.
		// Before the above "if" and the below "else" clauses got separated, the only code
		// executed was the "else" clause and it was rediculously slow.  I mean much much
		// slower than evaluating the etable energies.  REEAAALLLLY SLOW.

		// Update (11 July 2013 by V.K. Mulligan): The above now executes for both D- and L-amino
		// acid peptides and polypeptides.  What follows is for anything else.

		utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );
		for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
			Size const resconn_id1( r1_resconn_ids[ii] );
			Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

			Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
			Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

			// compute bond-angle energies with two atoms on rsd1
			utility::vector1< chemical::two_atom_set > const & rsd1_atoms_wi1_bond_of_ii(
				rsd1.type().atoms_within_one_bond_of_a_residue_connection( resconn_id1 ));
			for ( Size jj = 1; jj <= rsd1_atoms_wi1_bond_of_ii.size(); ++jj ) {
				debug_assert( rsd1_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno1 );
				Size const res1_lower_atomno = rsd1_atoms_wi1_bond_of_ii[ jj ].key2();

				//Real const angle = numeric::angle_radians( // why compute the angle if you're not going to use it?
				// rsd1.atom( res1_lower_atomno ).xyz(),
				// rsd1.atom( resconn_atomno1 ).xyz(),
				// rsd2.atom( resconn_atomno2 ).xyz() );
				std::string atm1name=rsd1.atom_name(res1_lower_atomno); boost::trim(atm1name);
				std::string atm2name=rsd1.atom_name(resconn_atomno1); boost::trim(atm2name);
				std::string atm3name=rsd2.atom_name(resconn_atomno2); boost::trim(atm3name);

				CartBondedParametersCOP ang_params;
				{
#if defined MULTI_THREADED
					utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
					ang_params =
						db_->lookup_angle(rsd1.type(), (rsd2.aa() == chemical::aa_pro || rsd2.aa() == chemical::aa_dpr /*D- or L-proline*/),
						atm1name, atm2name, atm3name, res1_lower_atomno, resconn_atomno1, -resconn_id1);
				}


				if ( ang_params->is_null() ) continue;
				Real Ktheta=ang_params->K(phi1,psi1), theta0=ang_params->mu(phi1,psi1);

				Vector f1(0.0), f2(0.0);
				Real theta(0.0), dE_dtheta;
				numeric::deriv::angle_p1_deriv(
					rsd1.xyz( res1_lower_atomno ), rsd1.xyz( resconn_atomno1 ), rsd2.xyz( resconn_atomno2 ), theta, f1, f2 );
				if ( linear_bonded_potential_ && std::fabs(theta - theta0)>1 ) {
					dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * ((theta - theta0)>0? 0.5 : -0.5);
				} else {
					dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * (theta - theta0);
				}
				r1_atom_derivs[ res1_lower_atomno ].f1() += dE_dtheta * f1;
				r1_atom_derivs[ res1_lower_atomno ].f2() += dE_dtheta * f2;

				numeric::deriv::angle_p2_deriv(
					rsd1.xyz( res1_lower_atomno ), rsd1.xyz( resconn_atomno1 ), rsd2.xyz( resconn_atomno2 ), theta, f1, f2 );
				r1_atom_derivs[ resconn_atomno1 ].f1() += dE_dtheta * f1;
				r1_atom_derivs[ resconn_atomno1 ].f2() += dE_dtheta * f2;

				numeric::deriv::angle_p1_deriv(
					rsd2.xyz( resconn_atomno2 ), rsd1.xyz( resconn_atomno1 ), rsd1.xyz( res1_lower_atomno ), theta, f1, f2 );
				r2_atom_derivs[ resconn_atomno2 ].f1() += dE_dtheta * f1;
				r2_atom_derivs[ resconn_atomno2 ].f2() += dE_dtheta * f2;
			}
		}
	}
}

void
CartesianBondedEnergy::eval_interresidue_angle_derivs_two_from_rsd2(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResidueCartBondedParameters const & rsd1params,
	ResidueCartBondedParameters const & rsd2params,
	Real phi2,
	Real psi2,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	using namespace core::chemical;
	bool rsd1_is_protein = (rsd1.aa() <= num_canonical_aas || core::chemical::is_canonical_D_aa(rsd1.aa()));
	bool rsd2_is_protein = (rsd2.aa() <= num_canonical_aas || core::chemical::is_canonical_D_aa(rsd2.aa()));
	bool rsd12_peptide_bonded = !rsd1.is_upper_terminus() && rsd1.residue_connection_partner( rsd1.upper_connect().index() ) == rsd2.seqpos();
	if ( rsd1_is_protein && rsd2_is_protein && rsd12_peptide_bonded ) {

		/// Assumption: rsd1 and rsd2 share a peptide bond and only a peptide bond.
		// amw: can we simply remove this bad assumption (e.g. for ncbb/ ncccs?
		//debug_assert( rsd2.residue_connection_partner( rsd2.lower_connect().index() ) == rsd1.seqpos() );
		// amw: for triazolamers let's try just that
		//debug_assert( rsd1.connections_to_residue( rsd2 ).size() == 1 );

		utility::vector1< ResidueCartBondedParameters::angle_parameter > const & aps( rsd2params.lower_connect_angle_params() );
		for ( Size ii = 1, iiend = aps.size(); ii <= iiend; ++ii ) {
			ResidueCartBondedParameters::Size3 const & atids( aps[ ii ].first );
			CartBondedParameters const & ang_params( *aps[ ii ].second );
			Real Ktheta = ang_params.K(phi2,psi2);
			Real theta0 = ang_params.mu(phi2,psi2);

			// temp -- rename these variables, ok?
			Size const res2_lower_atomno = atids[1], resconn_atomno2 = atids[2], resconn_atomno1 = rsd1params.bb_C_index();

			Vector f1(0.0), f2(0.0);
			Real theta(0.0), dE_dtheta;
			numeric::deriv::angle_p1_deriv(
				rsd2.xyz( res2_lower_atomno ), rsd2.xyz( resconn_atomno2 ), rsd1.xyz( resconn_atomno1 ), theta, f1, f2 );
			if ( linear_bonded_potential_ && std::fabs(theta - theta0)>1 ) {
				dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * ((theta - theta0)>0? 0.5 : -0.5);
			} else {
				dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * (theta - theta0);
			}
			r2_atom_derivs[ res2_lower_atomno ].f1() += dE_dtheta * f1;
			r2_atom_derivs[ res2_lower_atomno ].f2() += dE_dtheta * f2;

			numeric::deriv::angle_p2_deriv(
				rsd2.xyz( res2_lower_atomno ), rsd2.xyz( resconn_atomno2 ), rsd1.xyz( resconn_atomno1 ), theta, f1, f2 );
			r2_atom_derivs[ resconn_atomno2 ].f1() += dE_dtheta * f1;
			r2_atom_derivs[ resconn_atomno2 ].f2() += dE_dtheta * f2;

			numeric::deriv::angle_p1_deriv(
				rsd1.xyz( resconn_atomno1 ), rsd2.xyz( resconn_atomno2 ), rsd2.xyz( res2_lower_atomno ), theta, f1, f2 );
			r1_atom_derivs[ resconn_atomno1 ].f1() += dE_dtheta * f1;
			r1_atom_derivs[ resconn_atomno1 ].f2() += dE_dtheta * f2;
		}

	} else {
		// At the time of the writing of this comment, the CartesianBondEnergy is used only
		// to evaluate protein energetics, and the only inter-residue bonds it is evaluated
		// on is peptide bonds.  If that persists, then this "else" clause will never be executed.
		// If it does get executed, then this else clause ought to be optimized.
		// Before the above "if" and the below "else" clauses got separated, the only code
		// executed was the "else" clause and it was rediculously slow.  I mean much much
		// slower than evaluating the etable energies.  REEAAALLLLY SLOW.

		// Update (11 July 2013 by V.K. Mulligan): The above executes for both D- and L-amino acid
		// peptides and polypeptides, now.  What follows is for everything else.

		utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );
		for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
			Size const resconn_id1( r1_resconn_ids[ii] );
			Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

			Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
			Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

			utility::vector1< chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
				rsd2.type().atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));
			for ( Size jj = 1; jj <= rsd2_atoms_wi1_bond_of_ii.size(); ++jj ) {
				debug_assert( rsd2_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno2 );
				Size const res2_lower_atomno = rsd2_atoms_wi1_bond_of_ii[ jj ].key2();

				// lookup Ktheta and theta0
				// set pre-proline to false here -- it refers to rsd2+1 which would make this three-body
				std::string atm1name=rsd2.atom_name(res2_lower_atomno); boost::trim(atm1name);
				std::string atm2name=rsd2.atom_name(resconn_atomno2); boost::trim(atm2name);
				std::string atm3name=rsd1.atom_name(resconn_atomno1); boost::trim(atm3name);

				CartBondedParametersCOP ang_params;
				{
#if defined MULTI_THREADED
					utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
					ang_params =
						db_->lookup_angle(rsd2.type(), false, atm1name, atm2name, atm3name, res2_lower_atomno, resconn_atomno2, -resconn_id2);
				}

				if ( ang_params->is_null() ) continue;
				Real Ktheta=ang_params->K(phi2,psi2), theta0=ang_params->mu(phi2,psi2);

				Vector f1(0.0), f2(0.0);
				Real theta(0.0), dE_dtheta;
				numeric::deriv::angle_p1_deriv(
					rsd2.xyz( res2_lower_atomno ), rsd2.xyz( resconn_atomno2 ), rsd1.xyz( resconn_atomno1 ), theta, f1, f2 );
				if ( linear_bonded_potential_ && std::fabs(theta - theta0)>1 ) {
					dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * ((theta - theta0)>0? 0.5 : -0.5);
				} else {
					dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * (theta - theta0);
				}
				r2_atom_derivs[ res2_lower_atomno ].f1() += dE_dtheta * f1;
				r2_atom_derivs[ res2_lower_atomno ].f2() += dE_dtheta * f2;

				numeric::deriv::angle_p2_deriv(
					rsd2.xyz( res2_lower_atomno ), rsd2.xyz( resconn_atomno2 ), rsd1.xyz( resconn_atomno1 ), theta, f1, f2 );
				r2_atom_derivs[ resconn_atomno2 ].f1() += dE_dtheta * f1;
				r2_atom_derivs[ resconn_atomno2 ].f2() += dE_dtheta * f2;

				numeric::deriv::angle_p1_deriv(
					rsd1.xyz( resconn_atomno1 ), rsd2.xyz( resconn_atomno2 ), rsd2.xyz( res2_lower_atomno ), theta, f1, f2 );
				r1_atom_derivs[ resconn_atomno1 ].f1() += dE_dtheta * f1;
				r1_atom_derivs[ resconn_atomno1 ].f2() += dE_dtheta * f2;
			}
		}
	}
}

void
CartesianBondedEnergy::eval_interresidue_bond_length_derivs(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResidueCartBondedParameters const & rsd1params,
	ResidueCartBondedParameters const & rsd2params,
	Real phi2,
	Real psi2,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	using namespace core::chemical;
	bool rsd1_is_protein = (rsd1.aa() <= num_canonical_aas || core::chemical::is_canonical_D_aa(rsd1.aa()));
	bool rsd2_is_protein = (rsd2.aa() <= num_canonical_aas || core::chemical::is_canonical_D_aa(rsd2.aa()));
	bool rsd12_peptide_bonded = !rsd1.is_upper_terminus() && rsd1.residue_connection_partner( rsd1.upper_connect().index() ) == rsd2.seqpos();
	if ( rsd1_is_protein && rsd2_is_protein && rsd12_peptide_bonded ) {

		// lookup Kd and d0
		CartBondedParametersCOP len_params = rsd2params.cprev_n_bond_length_params();

		if ( len_params->is_null() ) return;
		Real const Kd = len_params->K(phi2,psi2);
		Real const d0 = len_params->mu(phi2,psi2);

		Size const r1at = rsd1params.bb_C_index();
		Size const r2at = rsd2params.bb_N_index();

		Vector f1(0.0), f2(0.0);
		Real d=0, dE_dd;

		numeric::deriv::distance_f1_f2_deriv( rsd2.xyz( r2at ), rsd1.xyz( r1at ), d, f1, f2 );
		if ( linear_bonded_potential_ && std::fabs(d - d0)>1 ) {
			dE_dd = (weights[ cart_bonded_length ] + weights[ cart_bonded ]) * Kd * ((d - d0)>0? 0.5 : -0.5);
		} else {
			dE_dd = (weights[ cart_bonded_length ] + weights[ cart_bonded ]) * Kd * (d - d0);
		}
		r2_atom_derivs[ r2at ].f1() += dE_dd * f1;
		r2_atom_derivs[ r2at ].f2() += dE_dd * f2;

		r1_atom_derivs[ r1at ].f1() -= dE_dd * f1;
		r1_atom_derivs[ r1at ].f2() -= dE_dd * f2;


	} else {
		// At the time of the writing of this comment, the CartesianBondEnergy is used only
		// to evaluate protein energetics, and the only inter-residue bonds it is evaluated
		// on is peptide bonds.  If that persists, then this "else" clause will never be executed.
		// If it does get executed, then this else clause ought to be optimized.
		// Before the above "if" and the below "else" clauses got separated, the only code
		// executed was the "else" clause and it was rediculously slow.  I mean much much
		// slower than evaluating the etable energies.  REEAAALLLLY SLOW.

		// Updated 11 July 2013 by V.K. Mulligan: the above code is for both D- and L- amino acid
		// peptides, now, while what follows is for everything else.

		utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );
		for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
			Size const resconn_id1( r1_resconn_ids[ii] );
			Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

			Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
			Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

			// [3] Bonds across connection
			std::string atm1name=rsd2.atom_name(resconn_atomno2); boost::trim(atm1name);
			std::string atm2name=rsd1.atom_name(resconn_atomno1); boost::trim(atm2name);

			CartBondedParametersCOP len_params;
			{
#if defined MULTI_THREADED
				utility::thread::ReadLockGuard lock( params_db_mutex_ );
#endif
				len_params = db_->lookup_length(rsd2.type(), false, atm1name, atm2name, resconn_atomno2, -resconn_id2);
			}

			if ( len_params->is_null() ) continue;
			Real Kd=len_params->K(phi2,psi2), d0=len_params->mu(phi2,psi2);

			Vector f1(0.0), f2(0.0);
			Real d=0, dE_dd;

			numeric::deriv::distance_f1_f2_deriv( rsd2.xyz( resconn_atomno2 ), rsd1.xyz( resconn_atomno1 ), d, f1, f2 );
			if ( linear_bonded_potential_ && std::fabs(d - d0)>1 ) {
				dE_dd = (weights[ cart_bonded_length ] + weights[ cart_bonded ]) * Kd * ((d - d0)>0? 0.5 : -0.5);
			} else {
				dE_dd = (weights[ cart_bonded_length ] + weights[ cart_bonded ]) * Kd * (d - d0);
			}

			r2_atom_derivs[ resconn_atomno2 ].f1() += dE_dd * f1;
			r2_atom_derivs[ resconn_atomno2 ].f2() += dE_dd * f2;

			numeric::deriv::distance_f1_f2_deriv( rsd1.xyz( resconn_atomno1 ), rsd2.xyz( resconn_atomno2 ), d, f1, f2 );
			r1_atom_derivs[ resconn_atomno1 ].f1() += dE_dd * f1;
			r1_atom_derivs[ resconn_atomno1 ].f2() += dE_dd * f2;
		}
	}
}

// deriv impropers
void
CartesianBondedEnergy::eval_interresidue_improper_derivatives(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	ResidueCartBondedParameters const & rsd1params,
	ResidueCartBondedParameters const & rsd2params,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {
	using namespace core::chemical;
	using numeric::constants::d::pi;

	//Multipliers for D-amino acids:
	//const core::Real d_multiplier1 = core::chemical::is_canonical_D_aa(res1.aa()) ? -1.0 : 1.0 ;
	const core::Real d_multiplier2 = core::chemical::is_canonical_D_aa(res2.aa()) ? -1.0 : 1.0 ;

	// backbone C-N-CA-H
	if ( !res1.is_protein() || !res2.is_protein() || !res1.is_polymer_bonded( res2 ) ) return;
	Real weight = weights[ cart_bonded_improper ] + weights[ cart_bonded_torsion ] + weights[ cart_bonded ];

	// backbone Oprev-Cprev-N-H
	if ( (res2.aa() != aa_pro && res2.aa() != aa_dpr /*NOT D- or L-proline*/) && rsd2params.bb_H_index() != 0 ) {
		CartBondedParametersCOP tor_params = rsd2params.oprev_cprev_n_h_interres_improper_params();
		if ( !tor_params->is_null() ) {
			Real const Kphi = tor_params->K(0,0);
			Real const phi0 = d_multiplier2 * tor_params->mu(0,0);
			Real const phi_step = 2*pi/tor_params->period();

			Vector f1(0.0), f2(0.0);
			Real phi=0, dE_dphi;

			Size const atm1( rsd1params.bb_O_index() );
			Size const atm2( rsd1params.bb_C_index() );
			Size const atm3( rsd2params.bb_N_index() );
			Size const atm4( rsd2params.bb_H_index() );

			numeric::deriv::dihedral_p1_cosine_deriv(
				res1.xyz( atm1 ), res1.xyz( atm2 ), res2.xyz( atm3 ), res2.xyz( atm4 ), phi, f1, f2 );
			Real del_phi = basic::subtract_radian_angles(phi, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );
			if ( linear_bonded_potential_ && std::fabs(del_phi)>1 ) {
				dE_dphi = weight * Kphi * (del_phi>0? 0.5 : -0.5);
			} else {
				dE_dphi = weight * Kphi * del_phi;
			}
			r1_atom_derivs[ atm1 ].f1() += dE_dphi * f1;
			r1_atom_derivs[ atm1 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res1.xyz( atm1 ), res1.xyz( atm2 ), res2.xyz( atm3 ), res2.xyz( atm4 ), phi, f1, f2 );
			r1_atom_derivs[ atm2 ].f1() += dE_dphi * f1;
			r1_atom_derivs[ atm2 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res2.xyz( atm4 ), res2.xyz( atm3 ), res1.xyz( atm2 ), res1.xyz( atm1 ), phi, f1, f2 );
			r2_atom_derivs[ atm3 ].f1() += dE_dphi * f1;
			r2_atom_derivs[ atm3 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv(
				res2.xyz( atm4 ), res2.xyz( atm3 ), res1.xyz( atm2 ), res1.xyz( atm1 ), phi, f1, f2 );
			r2_atom_derivs[ atm4 ].f1() += dE_dphi * f1;
			r2_atom_derivs[ atm4 ].f2() += dE_dphi * f2;
		}
	}


	// backbone CA-Cprev-N-H
	if ( (res2.aa() != aa_pro && res2.aa() != aa_dpr /*NOT D- or L-proline*/) && rsd2params.bb_H_index() != 0 ) {
		CartBondedParametersCOP tor_params = rsd2params.ca_cprev_n_h_interres_improper_params();
		if ( !tor_params->is_null() ) {
			Real const Kphi = tor_params->K(0,0);
			Real const phi0 = d_multiplier2 * tor_params->mu(0,0);
			Real const phi_step = 2*pi/tor_params->period();

			Vector f1(0.0), f2(0.0);
			Real phi=0, dE_dphi;

			Size const atm1( rsd2params.bb_CA_index() );
			Size const atm2( rsd1params.bb_C_index() );
			Size const atm3( rsd2params.bb_N_index() );
			Size const atm4( rsd2params.bb_H_index() );

			numeric::deriv::dihedral_p1_cosine_deriv(
				res2.xyz( atm1 ), res1.xyz( atm2 ), res2.xyz( atm3 ), res2.xyz( atm4 ), phi, f1, f2 );
			Real del_phi = basic::subtract_radian_angles(phi, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );
			if ( linear_bonded_potential_ && std::fabs(del_phi)>1 ) {
				dE_dphi = weight * Kphi * (del_phi>0? 0.5 : -0.5);
			} else {
				dE_dphi = weight * Kphi * del_phi;
			}
			r2_atom_derivs[ atm1 ].f1() += dE_dphi * f1;
			r2_atom_derivs[ atm1 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res2.xyz( atm1 ), res1.xyz( atm2 ), res2.xyz( atm3 ), res2.xyz( atm4 ), phi, f1, f2 );
			r1_atom_derivs[ atm2 ].f1() += dE_dphi * f1;
			r1_atom_derivs[ atm2 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res2.xyz( atm4 ), res2.xyz( atm3 ), res1.xyz( atm2 ), res2.xyz( atm1 ), phi, f1, f2 );
			r2_atom_derivs[ atm3 ].f1() += dE_dphi * f1;
			r2_atom_derivs[ atm3 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv(
				res2.xyz( atm4 ), res2.xyz( atm3 ), res1.xyz( atm2 ), res2.xyz( atm1 ), phi, f1, f2 );
			r2_atom_derivs[ atm4 ].f1() += dE_dphi * f1;
			r2_atom_derivs[ atm4 ].f2() += dE_dphi * f2;
		}
	}

	// backbone CA-Nnext-C-O
	{
		CartBondedParametersCOP tor_params = rsd1params.ca_nnext_c_o_interres_improper_params();
		if ( !tor_params->is_null() ) {  // however, if it is in the db but with 0 weight, that is OK
			core::Size const atm1 = rsd1params.bb_CA_index();
			core::Size const atm2 = rsd2params.bb_N_index();
			core::Size const atm3 = rsd1params.bb_C_index();
			core::Size const atm4 = rsd1params.bb_O_index();

			Real const Kphi = tor_params->K(0,0);
			Real const phi0 = d_multiplier2 * tor_params->mu(0,0);
			Real const phi_step=2*pi/tor_params->period();

			Vector f1(0.0), f2(0.0);
			Real phi=0, dE_dphi;

			numeric::deriv::dihedral_p1_cosine_deriv(
				res1.xyz( atm1 ), res2.xyz( atm2 ), res1.xyz( atm3 ), res1.xyz( atm4 ), phi, f1, f2 );
			Real del_phi = basic::subtract_radian_angles(phi, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );
			if ( linear_bonded_potential_ && std::fabs(del_phi)>1 ) {
				dE_dphi = (weights[ cart_bonded_improper ] + weights[ cart_bonded_torsion ] + weights[ cart_bonded ])
					* Kphi * (del_phi>0? 0.5 : -0.5);
			} else {
				dE_dphi = (weights[ cart_bonded_improper ] + weights[ cart_bonded_torsion ] + weights[ cart_bonded ])
					* Kphi * del_phi;
			}
			r1_atom_derivs[ atm1 ].f1() += dE_dphi * f1;
			r1_atom_derivs[ atm1 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res1.xyz( atm1 ), res2.xyz( atm2 ), res1.xyz( atm3 ), res1.xyz( atm4 ), phi, f1, f2 );
			r2_atom_derivs[ atm2 ].f1() += dE_dphi * f1;
			r2_atom_derivs[ atm2 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res1.xyz( atm4 ), res1.xyz( atm3 ), res2.xyz( atm2 ), res1.xyz( atm1 ), phi, f1, f2 );
			r1_atom_derivs[ atm3 ].f1() += dE_dphi * f1;
			r1_atom_derivs[ atm3 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv(
				res1.xyz( atm4 ), res1.xyz( atm3 ), res2.xyz( atm2 ), res1.xyz( atm1 ), phi, f1, f2 );
			r1_atom_derivs[ atm4 ].f1() += dE_dphi * f1;
			r1_atom_derivs[ atm4 ].f2() += dE_dphi * f2;
		}
	}

	// proline N planarity
	if ( res2.aa() == core::chemical::aa_pro || res2.aa() == core::chemical::aa_dpr /*D- or L-proline*/ ) {
		CartBondedParametersCOP tor_params = rsd2params.pro_cd_cprev_n_ca_interres_improper_params();
		if ( tor_params && !tor_params->is_null() ) {
			core::Size const atm1 = rsd2params.pro_CD_index();
			core::Size const atm2 = rsd1params.bb_C_index();
			core::Size const atm3 = rsd2params.bb_N_index();
			core::Size const atm4 = rsd2params.bb_CA_index();

			Real const Kphi = tor_params->K(0,0);
			Real const phi0 = d_multiplier2 * tor_params->mu(0,0);
			Real const phi_step=2*pi/tor_params->period();

			Vector f1(0.0), f2(0.0);
			Real phi=0, dE_dphi;

			numeric::deriv::dihedral_p1_cosine_deriv(
				res2.xyz( atm1 ), res1.xyz( atm2 ), res2.xyz( atm3 ), res2.xyz( atm4 ), phi, f1, f2 );
			Real del_phi = basic::subtract_radian_angles(phi, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );
			if ( linear_bonded_potential_ && std::fabs(del_phi)>1 ) {
				dE_dphi = weight * Kphi * (del_phi>0? 0.5 : -0.5);
			} else {
				dE_dphi = weight * Kphi * del_phi;
			}
			r2_atom_derivs[ atm1 ].f1() += dE_dphi * f1;
			r2_atom_derivs[ atm1 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res2.xyz( atm1 ), res1.xyz( atm2 ), res2.xyz( atm3 ), res2.xyz( atm4 ), phi, f1, f2 );
			r1_atom_derivs[ atm2 ].f1() += dE_dphi * f1;
			r1_atom_derivs[ atm2 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res2.xyz( atm4 ), res2.xyz( atm3 ), res1.xyz( atm2 ), res2.xyz( atm1 ), phi, f1, f2 );
			r2_atom_derivs[ atm3 ].f1() += dE_dphi * f1;
			r2_atom_derivs[ atm3 ].f2() += dE_dphi * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv(
				res2.xyz( atm4 ), res2.xyz( atm3 ), res1.xyz( atm2 ), res2.xyz( atm1 ), phi, f1, f2 );
			r2_atom_derivs[ atm4 ].f1() += dE_dphi * f1;
			r2_atom_derivs[ atm4 ].f2() += dE_dphi * f2;
		}
	}
}


Real
CartesianBondedEnergy::eval_score(
	Real val,  // actual value
	Real K,    // spring constant
	Real val0  // ideal value
) const
{
	Real const absdiff = std::fabs(val-val0);
	if ( linear_bonded_potential_ && absdiff > 1 ) {
		return 0.5 * K * absdiff;
	} else {
		return 0.5 * K * absdiff * absdiff;
	}
}

core::Size
CartesianBondedEnergy::version() const {
	return 2; //fd  Fix a typo in params file, update HIS_D params
}

} // namespace methods
} // namespace scoring
} // namespace core
