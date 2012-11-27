// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CartesianBondedEnergy.cc
/// @brief  Harmonic bondangle/bondlength/torsion constraints
/// @author

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
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
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

#include <core/scoring/PeptideBondedEnergyContainer.hh>
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


static basic::Tracer TR("core.scoring.CartesianBondedEnergy");
static basic::Tracer GEOMETRIES("core.scoring.CartesianBondedEnergy.GEOMETRIES");

// default spring constants
//  --> only if not specified in config files
static const Real K_LENGTH=300.0;
static const Real K_ANGLE=80.0;
static const Real K_TORSION=20.0;
static const Real K_TORSION_PROTON=10.0;  // proton chi
static const Real K_TORSION_IMPROPER=40.0;

IdealParametersDatabaseOP CartesianBondedEnergy::db_=NULL;


//////////////////////
/// EnergyMethod Creator
methods::EnergyMethodOP
CartesianBondedEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new CartesianBondedEnergy( options );
}

ScoreTypes
CartesianBondedEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( cart_bonded );
	sts.push_back( cart_bonded_angle );
	sts.push_back( cart_bonded_length );
	sts.push_back( cart_bonded_torsion );
	return sts;
}


////////////////////////
/// helper function
std::string get_restag( core::chemical::ResidueType const & restype ) {
	using namespace core::chemical;

	return restype.name3();
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
	if (basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user())
		params = basic::options::option[ basic::options::OptionKeys::score::bonded_params ]();
	if ( k_len_in >= 0 ) {
		k_len = k_len_in;
	} else if (basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user()) {
		if (params.size() >= 1) k_len = params[1];
	}

	if (k_ang_in >= 0) {
		k_ang = k_ang_in;
	} else if (basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user()) {
		if (params.size() >= 2) k_ang = params[2];
	}

	if (basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user()) {
		if (params.size() >= 3) k_tors = params[3];
		if (params.size() >= 4) k_tors_prot = params[4];
		if (params.size() >= 5) k_tors_improper = params[5];
	}
	if (k_tors_in >= 0)
		k_tors = k_tors_in;
	if (k_tors_prot_in >= 0)
		k_tors_prot = k_tors_prot_in;
	if (k_tors_improper_in >= 0)
		k_tors_improper = k_tors_improper_in;

	// now that we resolved the defaults, initialize
	init(k_len, k_ang, k_tors, k_tors_prot, k_tors_improper);
}

// read bb-dep and bb-indep parameter files from the rosetta database
void 
IdealParametersDatabase::init(Real k_len_in, Real k_ang_in, Real k_tors_in, Real k_tors_prot_in, Real k_tors_improper_in) {
	TR << "Initializing IdealParametersDatabase with default Ks="
		<< k_len_in << " , " << k_ang_in << " , " << k_tors_in << " , " << k_tors_prot_in << " , " << k_tors_improper_in << std::endl;
	k_length_ = k_len_in;
	k_angle_ = k_ang_in;
	k_torsion_ = k_tors_in;
	k_torsion_proton_ = k_tors_prot_in;
	k_torsion_improper_ = k_tors_improper_in;

	// read bb-independent parameters
	std::string name3, atom1, atom2, atom3, atom4;
	Real mu_d, K_d;
	Size period;
	std::string libpath = basic::options::option[ basic::options::OptionKeys::score::bonded_params_dir ]();
	bbdep_bond_params_ = basic::options::option[ basic::options::OptionKeys::corrections::score::bbdep_bond_params ]();
	bbdep_bond_devs_ = basic::options::option[ basic::options::OptionKeys::corrections::score::bbdep_bond_devs ]();

	{
		utility::io::izstream instream;
		atm_name_pair tuple;
		basic::database::open( instream, libpath+"/default-lengths.txt");
		while (instream) {
			instream >> name3 >> atom1 >> atom2 >> mu_d >> K_d;
			tuple = boost::make_tuple( name3, atom1, atom2 );
			CartBondedParametersOP params_i = new BBIndepCartBondedParameters(mu_d, K_d);
			bondlengths_indep_.insert( std::make_pair( tuple, params_i ) );
			tuple = boost::make_tuple( name3, atom2, atom1 );
			params_i = new BBIndepCartBondedParameters(mu_d, K_d);
			bondlengths_indep_.insert( std::make_pair( tuple, params_i ) );
		}
		TR << "Read " << bondlengths_indep_.size() << " bb-independent lengths." << std::endl;
	}

	{
		utility::io::izstream instream;
		atm_name_triple tuple;
		basic::database::open( instream, libpath+"/default-angles.txt");
		while (instream) {
			instream >> name3 >> atom1 >> atom2 >> atom3 >> mu_d >> K_d;
			tuple = boost::make_tuple( name3, atom1, atom2, atom3 );
			CartBondedParametersOP params_i = new BBIndepCartBondedParameters(mu_d, K_d);
			bondangles_indep_.insert( std::make_pair( tuple, params_i) );
			tuple = boost::make_tuple( name3, atom3, atom2, atom1 );
			bondangles_indep_.insert( std::make_pair( tuple, params_i) );
		}
		TR << "Read " << bondangles_indep_.size() << " bb-independent lengths." << std::endl;
	}

	{
		utility::io::izstream instream;
		atm_name_quad tuple;
		basic::database::open( instream, libpath+"/default-torsions.txt");
		while (instream) {
			instream >> name3 >> atom1 >> atom2 >> atom3 >> atom4 >> mu_d >> K_d >> period;
			tuple = boost::make_tuple( name3, atom1, atom2, atom3, atom4 );
			CartBondedParametersOP params_i = new BBIndepCartBondedParameters(mu_d, K_d, period);
			torsions_indep_.insert( std::make_pair( tuple, params_i) );
			tuple = boost::make_tuple( name3, atom4, atom3, atom2, atom1 );
			torsions_indep_.insert( std::make_pair( tuple, params_i) );
		}
		TR << "Read " << torsions_indep_.size() << " bb-independent lengths." << std::endl;
	}

	if (!bbdep_bond_params_) return;

	// read bb-dep parameters
	//NOTE this must be read after the bbindep params are read 
	//   as deviations for low count rama bins are smoothed using the bb-indep deviations
	read_bbdep_table( libpath+"/bbdep-all-but-gly-pro-xpro-ile-val.txt", bondlengths_bbdep_def_, bondangles_bbdep_def_, "ALA" );
	read_bbdep_table( libpath+"/bbdep-graphdata-gly.txt", bondlengths_bbdep_gly_, bondangles_bbdep_gly_, "GLY" );
	read_bbdep_table( libpath+"/bbdep-graphdata-ile-val.txt", bondlengths_bbdep_valile_, bondangles_bbdep_valile_, "VAL" );
	read_bbdep_table( libpath+"/bbdep-graphdata-pro.txt", bondlengths_bbdep_pro_, bondangles_bbdep_pro_, "PRO" );
	read_bbdep_table( libpath+"/bbdep-graphdata-xpro.txt", bondlengths_bbdep_prepro_, bondangles_bbdep_prepro_, "ALA" );
}


// Read bb independent tables
// smooth using bbdep data corresponding to residue 'resbase'
void
IdealParametersDatabase::read_bbdep_table(
			std::string filename,
			boost::unordered_map< atm_name_single, CartBondedParametersOP > &bondlengths,
			boost::unordered_map< atm_name_pair, CartBondedParametersOP > &bondangles,
			std::string resbase )
{
	using numeric::constants::f::pi;
	using namespace ObjexxFCL;

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
	if (!temp) temp = bondlengths_indep_[boost::make_tuple("*","C", "N")];
	if (!temp) temp = new BBIndepCartBondedParameters(0,0);
	core::Real CN_indep_avg=temp->mu(0,0),CN_indep_K=temp->K(0,0);

	temp = bondlengths_indep_[boost::make_tuple(resbase,"N", "CA")];
	if (!temp) temp = bondlengths_indep_[boost::make_tuple("*","N", "CA")];
	if (!temp) temp = new BBIndepCartBondedParameters(0,0);
	core::Real NCA_indep_avg=temp->mu(0,0),NCA_indep_K=temp->K(0,0);

	temp = bondlengths_indep_[boost::make_tuple(resbase,"CA", "CB")];
	if (!temp) temp = bondlengths_indep_[boost::make_tuple("*","CA", "CB")];
	if (!temp) temp = new BBIndepCartBondedParameters(0,0);
	core::Real CACB_indep_avg=temp->mu(0,0),CACB_indep_K=temp->K(0,0);

	temp = bondlengths_indep_[boost::make_tuple(resbase,"CA", "C")];
	if (!temp) temp = bondlengths_indep_[boost::make_tuple("*","CA", "C")];
	if (!temp) temp = new BBIndepCartBondedParameters(0,0);
	core::Real CAC_indep_avg=temp->mu(0,0),CAC_indep_K=temp->K(0,0);

	temp = bondlengths_indep_[boost::make_tuple(resbase,"C", "O")];
	if (!temp) temp = bondlengths_indep_[boost::make_tuple("*","C", "O")];
	if (!temp) temp = new BBIndepCartBondedParameters(0,0);
	core::Real CO_indep_avg=temp->mu(0,0),CO_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"C", "N", "CA")];
	if (!temp) temp = bondangles_indep_[boost::make_tuple("*","C", "N", "CA")];
	if (!temp) temp = new BBIndepCartBondedParameters(0,0);
	core::Real CNCA_indep_avg=temp->mu(0,0),CNCA_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"N", "CA", "CB")];
	if (!temp) temp = bondangles_indep_[boost::make_tuple("*","N", "CA", "CB")];
	if (!temp) temp = new BBIndepCartBondedParameters(0,0);
	core::Real NCACB_indep_avg=temp->mu(0,0),NCACB_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"N", "CA", "C")];
	if (!temp) temp = bondangles_indep_[boost::make_tuple("*","N", "CA", "C")];
	if (!temp) temp = new BBIndepCartBondedParameters(0,0);
	core::Real NCAC_indep_avg=temp->mu(0,0),NCAC_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"CB", "CA", "C")];
	if (!temp) temp = bondangles_indep_[boost::make_tuple("*","CB", "CA", "C")];
	if (!temp) temp = new BBIndepCartBondedParameters(0,0);
	core::Real CBCAC_indep_avg=temp->mu(0,0),CBCAC_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"CA", "C", "O")];
	if (!temp) temp = bondangles_indep_[boost::make_tuple("*","CA", "C", "O")];
	if (!temp) temp = new BBIndepCartBondedParameters(0,0);
	core::Real CACO_indep_avg=temp->mu(0,0),CACO_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"CA", "C", "N")];
	if (!temp) temp = bondangles_indep_[boost::make_tuple("*","CA", "C", "N")];
	if (!temp) temp = new BBIndepCartBondedParameters(0,0);
	core::Real CACN_indep_avg=temp->mu(0,0),CACN_indep_K=temp->K(0,0);

	temp = bondangles_indep_[boost::make_tuple(resbase,"O", "C", "N")];
	if (!temp) temp = bondangles_indep_[boost::make_tuple("*","O", "C", "N")];
	if (!temp) temp = new BBIndepCartBondedParameters(0,0);
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
	while(getline(instream, fileline)) {
		if(fileline[0] == '#') continue;
		std::stringstream(fileline) >> phiL >> phiH >> psiL >> psiH >> nobs >> phiAvg >> phiDev >> psiAvg >> psiDev >> OmegaAvg >> OmegaDev >>  CNavg >> CNdev
			>> NCAavg >> NCAdev >>  CACBavg >> CACBdev >>  CACavg >> CACdev >>  COavg >> COdev
			>> CNCAavg >> CNCAdev >>  NCACBavg >> NCACBdev >>  NCACavg >> NCACdev >>  CBCACavg >> CBCACdev
			>> CACOavg >> CACOdev >>  CACNavg >> CACNdev >>  OCNavg >> OCNdev
			>> c1avg >> c1dev >> c2avg >> c2dev >> c3avg >> c3dev >> c4avg >> c4dev >> Zavg >> Zdev;

		Size phibin = (Size)std::floor(phiL/10.0 + 0.5);
		Size psibin = (Size)std::floor(psiL/10.0 + 0.5);

		CNavg_tbl(phibin+1,psibin+1) = CNavg;
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

		if (!bbdep_bond_devs_ || nobs<=3 || CN_indep_K==0) CNdev_tbl(phibin+1,psibin+1) = CN_indep_K;
		else CNdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CN_indep_K) + (nobs-1)*( (2*(CNdev)*(CNdev))/DEV_SCALE) )/ (M+nobs-1) );

		if (!bbdep_bond_devs_ || nobs<=3 || NCA_indep_K==0 ) NCAdev_tbl(phibin+1,psibin+1) = NCA_indep_K;
		else NCAdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/NCA_indep_K) + (nobs-1)*( (2*(NCAdev)*(NCAdev))/DEV_SCALE) )/ (M+nobs-1) );

		if (!bbdep_bond_devs_ || nobs<=3 || CACB_indep_K==0 ) CACBdev_tbl(phibin+1,psibin+1) = CACB_indep_K;
		else CACBdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CACB_indep_K) + (nobs-1)*( (2*(CACBdev)*(CACBdev))/DEV_SCALE) )/ (M+nobs-1) );

		if (!bbdep_bond_devs_ || nobs<=3 || CAC_indep_K==0 ) CACdev_tbl(phibin+1,psibin+1) = CAC_indep_K;
		else CACdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CAC_indep_K) + (nobs-1)*( (2*(CACdev)*(CACdev))/DEV_SCALE) )/ (M+nobs-1) );

		if (!bbdep_bond_devs_ || nobs<=3 || CO_indep_K==0 ) COdev_tbl(phibin+1,psibin+1) = CO_indep_K;
		else COdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CO_indep_K) + (nobs-1)*( (2*(COdev)*(COdev))/DEV_SCALE) )/ (M+nobs-1) );

		if (!bbdep_bond_devs_ || nobs<=3 || CNCA_indep_K==0 ) CNCAdev_tbl(phibin+1,psibin+1) = CNCA_indep_K;
		else CNCAdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CNCA_indep_K) + (nobs-1)*( (2*(CNCAdev*pi/180)*(CNCAdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );

		if (!bbdep_bond_devs_ || nobs<=3 || NCACB_indep_K==0 ) NCACBdev_tbl(phibin+1,psibin+1) = NCACB_indep_K;
		else NCACBdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/NCACB_indep_K) + (nobs-1)*( (2*(NCACBdev*pi/180)*(NCACBdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );

		if (!bbdep_bond_devs_ || nobs<=3 || NCAC_indep_K==0 ) NCACdev_tbl(phibin+1,psibin+1) = NCAC_indep_K;
		else NCACdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/NCAC_indep_K) + (nobs-1)*( (2*(NCACdev*pi/180)*(NCACdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );

		if (!bbdep_bond_devs_ || nobs<=3 || CBCAC_indep_K==0 ) CBCACdev_tbl(phibin+1,psibin+1) = CBCAC_indep_K;
		else CBCACdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CBCAC_indep_K) + (nobs-1)*( (2*(CBCACdev*pi/180)*(CBCACdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );

		if (!bbdep_bond_devs_ || nobs<=3 || CACO_indep_K==0 ) CACOdev_tbl(phibin+1,psibin+1) = CACO_indep_K;
		else CACOdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CACO_indep_K) + (nobs-1)*( (2*(CACOdev*pi/180)*(CACOdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );

		if (!bbdep_bond_devs_ || nobs<=3 || CACN_indep_K==0 ) CACNdev_tbl(phibin+1,psibin+1) = CACN_indep_K;
		else CACNdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/CACN_indep_K) + (nobs-1)*( (2*(CACNdev*pi/180)*(CACNdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );

		if (!bbdep_bond_devs_ || nobs<=3 || OCN_indep_K==0 ) OCNdev_tbl(phibin+1,psibin+1) = OCN_indep_K;
		else OCNdev_tbl(phibin+1,psibin+1) = 1.0 / ( (M*(1/OCN_indep_K) + (nobs-1)*( (2*(OCNdev*pi/180)*(OCNdev*pi/180))/DEV_SCALE) )/ (M+nobs-1) );
	}

	// fill in missing values
	for (int i=1; i<=36; ++i) {
		for (int j=1; j<=36; ++j) {
			if (Ns(i,j) == 0) {
				CNavg_tbl(i,j) = CN_indep_avg; CNdev_tbl(i,j) = CN_indep_K;
				NCAavg_tbl(i,j) = NCA_indep_avg; NCAdev_tbl(i,j) = NCA_indep_K;
				CACBavg_tbl(i,j) = CACB_indep_avg; CACBdev_tbl(i,j) = CACB_indep_K;
				CACavg_tbl(i,j) = CAC_indep_avg; CACdev_tbl(i,j) = CAC_indep_K;
				COavg_tbl(i,j) = CO_indep_avg; COdev_tbl(i,j) = CO_indep_K;
				CNCAavg_tbl(i,j) = CNCA_indep_avg; CNCAdev_tbl(i,j) = CNCA_indep_K;
				NCACBavg_tbl(i,j) = NCACB_indep_avg; NCACBdev_tbl(i,j) = NCACB_indep_K;
				NCACavg_tbl(i,j) = NCAC_indep_avg; NCACdev_tbl(i,j) = NCAC_indep_K;
				CBCACavg_tbl(i,j) = CBCAC_indep_avg; CBCACdev_tbl(i,j) = CBCAC_indep_K;
				CACOavg_tbl(i,j) = CACO_indep_avg; CACOdev_tbl(i,j) = CACO_indep_K;
				CACNavg_tbl(i,j) = CACN_indep_avg; CACNdev_tbl(i,j) = CACN_indep_K;
				OCNavg_tbl(i,j) = OCN_indep_avg; OCNdev_tbl(i,j) = OCN_indep_K;
			}
		}
	}

	// make the database entries
	bondlengths[boost::make_tuple("N","C")] = bondlengths[boost::make_tuple("C","N")] = new BBDepCartBondedParameters(CNavg_tbl, CNdev_tbl ,"C-N");
	bondlengths[boost::make_tuple("CA","N")] = bondlengths[boost::make_tuple("N","CA")] = new BBDepCartBondedParameters(NCAavg_tbl, NCAdev_tbl,"CA-N");
	bondlengths[boost::make_tuple("CB","CA")] = bondlengths[boost::make_tuple("CA","CB")] = new BBDepCartBondedParameters(CACBavg_tbl, CACBdev_tbl,"CB-CA");
	bondlengths[boost::make_tuple("C","CA")] = bondlengths[boost::make_tuple("CA","C")] = new BBDepCartBondedParameters(CACavg_tbl, CACdev_tbl,"C-CA");
	bondlengths[boost::make_tuple("O","C")] = bondlengths[boost::make_tuple("C","O")] = new BBDepCartBondedParameters(COavg_tbl, COdev_tbl,"C-O");
	bondangles[boost::make_tuple("CA","N","C")] = bondangles[boost::make_tuple("C","N","CA")] = new BBDepCartBondedParameters(CNCAavg_tbl, CNCAdev_tbl, "C-N-CA");
	bondangles[boost::make_tuple("CB","CA","N")] = bondangles[boost::make_tuple("N","CA","CB")] = new BBDepCartBondedParameters(NCACBavg_tbl, NCACBdev_tbl, "N-CA-CB");
	bondangles[boost::make_tuple("C","CA","N")] = bondangles[boost::make_tuple("N","CA","C")] = new BBDepCartBondedParameters(NCACavg_tbl, NCACdev_tbl, "N-CA-C");
	bondangles[boost::make_tuple("C","CA","CB")] = bondangles[boost::make_tuple("CB","CA","C")] = new BBDepCartBondedParameters(CBCACavg_tbl, CBCACdev_tbl, "CB-CA-C");
	bondangles[boost::make_tuple("O","C","CA")] = bondangles[boost::make_tuple("CA","C","O")] = new BBDepCartBondedParameters(CACOavg_tbl, CACOdev_tbl, "CA-C-O");
	bondangles[boost::make_tuple("N","C","CA")] = bondangles[boost::make_tuple("CA","C","N")] = new BBDepCartBondedParameters(CACNavg_tbl, CACNdev_tbl, "CA-C-N");
	bondangles[boost::make_tuple("N","C","O")] = bondangles[boost::make_tuple("O","C","N")] = new BBDepCartBondedParameters(OCNavg_tbl, OCNdev_tbl, "O-C-N");
}



void
IdealParametersDatabase::lookup_bondangle_buildideal(
		core::conformation::Residue const & res,
		int atm1, int atm2, int atm3, Real &Ktheta, Real &theta0) {
	Ktheta = k_angle_;

	// Create mini-conformation for idealized residue
	chemical::ResidueType const & restype = res.type();
	conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );
	
	// get angle
	numeric::xyzVector<core::Real> x,y,z; // atom coords
	if (atm1 > 0) {
		x = newres->atom( atm1 ).xyz();
	} else {
		x = newres->residue_connection( -atm1 ).icoor().build( restype );
	}

	if (atm2 > 0) {
		y = newres->atom( atm2 ).xyz();
	} else {
		y = newres->residue_connection( -atm2 ).icoor().build( restype );
	}

	if (atm3 > 0) {
		z = newres->atom( atm3 ).xyz();
	} else {
		z = newres->residue_connection( -atm3 ).icoor().build( restype );
	}

	theta0 = numeric::angle_radians ( x,y,z );
	
	// EXCEPTIONS
	// ignore proline C-ND which is handled by pro_close
	if (restype.aa() == core::chemical::aa_pro) {
		bool hasCD = (atm1>0 && restype.atom_name(atm1)==" CD ")
	    || (atm2>0 && restype.atom_name(atm2)==" CD ")
		  || (atm3>0 && restype.atom_name(atm3)==" CD ");
		bool hasN = (atm1>0 && restype.atom_name(atm1)==" N  ")
	    || (atm2>0 && restype.atom_name(atm2)==" N  ")
		  || (atm3>0 && restype.atom_name(atm3)==" N  ");
		if (hasCD && hasN)
			Ktheta = theta0 = 0.0;
	}

	// fpd ignore centroid angle in ALA and GLY
	if ( (restype.aa() == core::chemical::aa_ala || restype.aa() == core::chemical::aa_gly) &&
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
			core::conformation::Residue const & res,
			int atm1, int atm2, Real &Kd, Real &d0 ) {
	Kd = k_length_;

	// Create mini-conformation for idealized residue
	chemical::ResidueType const & restype = res.type();
	core::conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );

	// get length
	numeric::xyzVector<core::Real> x,y,z; // atom coords
	if (atm1 > 0) {
		x = newres->atom( atm1 ).xyz();
	}
	else {
		x = newres->residue_connection( -atm1 ).icoor().build( restype );
	}

	if (atm2 > 0) {
		y = newres->atom( atm2 ).xyz();
	}
	else {
		y = newres->residue_connection( -atm2 ).icoor().build( restype );
	}

	d0 = (x-y).length();
	 
	//  ignore proline C-ND which is handled by pro_close
	if (restype.aa() == core::chemical::aa_pro &&
	    (atm1>0 && atm2>0) &&
	    ( ( restype.atom_name(atm1) == " CD " && restype.atom_name(atm2) == " N  ") ||
	      ( restype.atom_name(atm2) == " CD " && restype.atom_name(atm1) == " N  ") ) ) {
		Kd = d0 = 0.0;
	}
}
	


//////////////////////
/// Torsion Database
// lookup ideal intra-res torsion
// torsion DB is unique is that it never tries to build from ideal
CartBondedParametersCAP
IdealParametersDatabase::lookup_torsion(
			core::conformation::Residue const & res,
			std::string atm1_name, std::string atm2_name, std::string atm3_name, std::string atm4_name ) {
	using namespace core::chemical;

	// use 'annotated sequence' to id this restype
	chemical::ResidueType const & restype = res.type();
	std::string restag = get_restag( restype );

	// there is no bb-dep torsion database
	// lookup in bb-indep table
	// this can probably be made way faster
	atm_name_quad tuple1( restag, atm1_name,atm2_name,atm3_name,atm4_name );
	boost::unordered_map<atm_name_quad,CartBondedParametersOP>::iterator b_it = torsions_indep_.find( tuple1 );
	if ( b_it != torsions_indep_.end() ) {
		return CartBondedParametersCAP(*b_it->second);
	}

	atm_name_quad tuple2( "*", atm1_name,atm2_name,atm3_name,atm4_name );
	b_it = torsions_indep_.find( tuple2 );
	if ( b_it != torsions_indep_.end() ) {
		return CartBondedParametersCAP(*b_it->second);
	}

	// if we don't find this torsion in the table, it's unconstrained
	return NULL;
}


/// Angle Database
///   build from ideal if not found
CartBondedParametersCAP												
IdealParametersDatabase::lookup_angle(
			core::conformation::Residue const & res,
			bool pre_proline,
			std::string atm1_name, std::string atm2_name, std::string atm3_name,
			int atm1idx, int atm2idx, int atm3idx ) {
	using namespace core::chemical;

	// use 'annotated sequence' to id this restype
	chemical::ResidueType const &restype = res.type();
	std::string restag = get_restag( restype );
	atm_name_triple tuple( restag, atm1_name,atm2_name,atm3_name );

	// 1 lookup in bb-dep table
	// figure out what table to look in
	boost::unordered_map< atm_name_pair, CartBondedParametersOP > *angle_table;
	if (restype.aa() == core::chemical::aa_gly) {
		angle_table = &bondangles_bbdep_gly_;
	} else if (restype.aa() == core::chemical::aa_pro) {
		angle_table = &bondangles_bbdep_pro_;
	} else if ( pre_proline ) {
		angle_table = &bondangles_bbdep_prepro_;
	} else if (restype.aa() == core::chemical::aa_ile || restype.aa() == core::chemical::aa_val) {
		angle_table = &bondangles_bbdep_valile_;
	} else {
		angle_table = &bondangles_bbdep_def_;
	}
	atm_name_pair alt_tuple( atm1_name,atm2_name,atm3_name );
	boost::unordered_map<atm_name_pair,CartBondedParametersOP>::iterator a_it = angle_table->find( alt_tuple );
	if ( a_it != angle_table-> end() ) {
		 return CartBondedParametersCAP(*a_it->second);
	}

	// 2 lookup in bb-indep table
	boost::unordered_map<atm_name_triple,CartBondedParametersOP>::iterator b_it = bondangles_indep_.find( tuple );
	if ( b_it != bondangles_indep_.end() ) {
		 return CartBondedParametersCAP(*b_it->second);
	}

	atm_name_triple tuple_alt( "*", atm1_name,atm2_name,atm3_name );
	b_it = bondangles_indep_.find( tuple_alt );
	if ( b_it != bondangles_indep_.end() ) {
		return CartBondedParametersCAP(*b_it->second);
	}

	// 3 build from ideal, add to bb-indep table
	Real Ktheta, theta0;
	lookup_bondangle_buildideal( res, atm1idx, atm2idx, atm3idx, Ktheta, theta0 );
	bondangles_indep_[ tuple ] = new BBIndepCartBondedParameters( theta0, Ktheta );
	TR << "Adding undefined angle "
	   << tuple.get<0>() << ": " << tuple.get<1>() << "," << tuple.get<2>() << "," << tuple.get<3>() << " to DB with"
	   << " theta0 = " << theta0 << " , Ktheta = " << Ktheta << std::endl; 

	return CartBondedParametersCAP(*bondangles_indep_[ tuple ]);
}

/// BondLength Database
///   build from ideal if not found
CartBondedParametersCAP
IdealParametersDatabase::lookup_length(
			conformation::Residue const & res,
			bool pre_proline,
			std::string atm1_name, std::string atm2_name,
			int atm1idx, int atm2idx ) {
	chemical::ResidueType const &restype = res.type();

	// use 'annotated sequence' to id this restype
	std::string restag = get_restag( restype );
	atm_name_pair tuple( restag, atm1_name,atm2_name);

	// 1 lookup in bb-dep table
	// figure out what table to look in
	boost::unordered_map< atm_name_single, CartBondedParametersOP > *length_table;
	if (restype.aa() == core::chemical::aa_gly) {
		length_table = &bondlengths_bbdep_gly_;
	} else if (restype.aa() == core::chemical::aa_pro) {
		length_table = &bondlengths_bbdep_pro_;
	} else if (restype.aa() == core::chemical::aa_ile || restype.aa() == core::chemical::aa_val) {
		length_table = &bondlengths_bbdep_valile_;
	} else if ( pre_proline ) {
		length_table = &bondlengths_bbdep_prepro_;
	} else {
		length_table = &bondlengths_bbdep_def_;
	}

	atm_name_single alt_tuple( atm1_name,atm2_name );
	boost::unordered_map<atm_name_single,CartBondedParametersOP>::iterator a_it = length_table->find( alt_tuple );
	if ( a_it != length_table->end() ) {
		 return CartBondedParametersCAP(*a_it->second);
	}

	// 2 lookup in bb-indep table
	boost::unordered_map<atm_name_pair,CartBondedParametersOP>::iterator b_it = bondlengths_indep_.find( tuple );
	if ( b_it != bondlengths_indep_.end() ) {
		 return CartBondedParametersCAP(*b_it->second);
	}

	atm_name_pair tuple_alt( "*", atm1_name,atm2_name );
	b_it = bondlengths_indep_.find( tuple_alt );
	if ( b_it != bondlengths_indep_.end() ) {
		return CartBondedParametersCAP(*b_it->second);
	}

	// 3 build from ideal, add to bb-indep table
	Real Kd, d0;
	lookup_bondlength_buildideal( res, atm1idx, atm2idx, Kd, d0 );
	bondlengths_indep_[ tuple ] = new BBIndepCartBondedParameters( d0, Kd );
	TR << "Adding undefined length "
	   << tuple.get<0>() << ": " << tuple.get<1>() << "," << tuple.get<2>() << "," << " to DB with"
	   << " d0 = " << d0 << " , Kd = " << Kd << std::endl; 

	return CartBondedParametersCAP(*bondlengths_indep_[ tuple ]);
}



// old-style interface to database
// somewhat slow as it will always build the ideal residue on the fly
void													
IdealParametersDatabase::lookup_torsion_legacy( core::chemical::ResidueType const & restype,
    int atm1, int atm2, int atm3, int atm4, Real &Kphi, Real &phi0, Real &phi_step ) {
	using namespace core::chemical;

	Kphi=k_torsion_; // default
	phi_step=0;

	if (!restype.is_protein()) {
		Kphi = phi0 = 0;
		return;
	}

	// use 'annotated sequence' to id this restype
	std::string restag = get_restag( restype );

	// Create mini-conformation for idealized residue
	core::conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );
	core::pose::Pose newpose;
	newpose.append_residue_by_bond(*newres);

	// figure out if we need to constrain this torsion
	bool need_to_constrain=true, proton_chi=false;

	// backbone -- we never need to constrain intrares backbones (for now)
	if ( atm2 <= (int)newres->last_backbone_atom() && atm3 <= (int)newres->last_backbone_atom()) {
			need_to_constrain = false;
	}

	// chi
	for ( Size j=1, j_end = restype.nchi(); j<= j_end; ++j ) {
		id::AtomID n1,n2,n3,n4;
		newpose.conformation().get_torsion_angle_atom_ids( id::TorsionID(1,id::CHI,j) ,n1,n2,n3,n4);
		if ( ((int)n2.atomno() == atm2 && (int)n3.atomno() == atm3) || ((int)n2.atomno() == atm3 && (int)n3.atomno() == atm2) ) {
			if (restype.is_proton_chi( j ) ) {
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
		if (restype.aa() == core::chemical::aa_pro) {
			bool hasCD = (restype.atom_name(atm1)==" CD ") || (restype.atom_name(atm2)==" CD ")
								|| (restype.atom_name(atm3)==" CD ") || (restype.atom_name(atm4)==" CD ");
			bool hasN = (restype.atom_name(atm1)==" N  ") || (restype.atom_name(atm2)==" N  ")
							 || (restype.atom_name(atm3)==" N  ") || (restype.atom_name(atm4)==" N  ");
			if (hasCD && hasN) {
				Kphi=0.0; phi0=0.0;
			}
		}
		if ( restype.has_variant_type(chemical::CUTPOINT_UPPER) &&
				 ( ( restype.atom_name(atm2) == " N  " && restype.atom_name(atm3) == " CA ") ||
					 ( restype.atom_name(atm3) == " N  " && restype.atom_name(atm2) == " CA ")  ) ) {
			Kphi=0.0; phi0=0.0;
		}
		if ( restype.atom_name(atm1) == " CEN" || restype.atom_name(atm4) == " CEN") {
			Kphi=0.0; phi0=0.0;
		}
		if ( restype.aa() == aa_cys && restype.has_variant_type( chemical::DISULFIDE )  &&
				(restype.atom_name(atm2) == " SG " || restype.atom_name(atm3) == " SG ") ) {
			Kphi=0.0; phi0=0.0;
		}
		if ( restype.aa() == aa_arg  &&
				(   restype.atom_name(atm1) == " NH1" || restype.atom_name(atm4) == " NH1" 
				 || restype.atom_name(atm1) == " NH2" || restype.atom_name(atm4) == " NH2") ) {
			phi_step = numeric::constants::f::pi;
		}
		if ( restype.has_variant_type(chemical::LOWER_TERMINUS) && 
				 ( ( restype.atom_name(atm2) == " N  " && restype.atom_name(atm3) == " CA ") ||
					 ( restype.atom_name(atm3) == " N  " && restype.atom_name(atm2) == " CA ")  ) ) {
			Kphi=k_torsion_proton_;
			phi_step = 2.0*numeric::constants::f::pi/3.0;
		}

		if (proton_chi) {
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
IdealParametersDatabase::lookup_angle_legacy( core::pose::Pose const & pose, core::conformation::Residue const & res,
	  int atm1, int atm2, int atm3, Real &Ktheta, Real &theta0) {
	using namespace core::chemical;

	Ktheta=k_angle_;
	ResidueType const & restype = res.type();
	std::string restag = get_restag( restype );

	// Create mini-conformation for idealized residue
	conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );
	
	// use these for mm Ktheta lookup
	Size mmatm1, mmatm2, mmatm3;						
	
	// get angle
	numeric::xyzVector<core::Real> x,y,z; // atom coords
	if (atm1 > 0) {
		x = newres->atom( atm1 ).xyz();
	}
	else {
		x = newres->residue_connection( -atm1 ).icoor().build( restype );
	}
	if (atm2 > 0) {
		y = newres->atom( atm2 ).xyz();
	}
	else {
		y = newres->residue_connection( -atm2 ).icoor().build( restype );
	}
	if (atm3 > 0) {
		z = newres->atom( atm3 ).xyz();
	}
	else {
		z = newres->residue_connection( -atm3 ).icoor().build( restype );
	}

	theta0 = numeric::angle_radians ( x,y,z );
	
	//fpd ignore proline C-ND which is handled by pro_close
	if (restype.aa() == core::chemical::aa_pro) {
		bool hasCD = (atm1>0 && restype.atom_name(atm1)==" CD ")
	    || (atm2>0 && restype.atom_name(atm2)==" CD ")
		  || (atm3>0 && restype.atom_name(atm3)==" CD ");
		bool hasN = (atm1>0 && restype.atom_name(atm1)==" N  ")
	    || (atm2>0 && restype.atom_name(atm2)==" N  ")
		  || (atm3>0 && restype.atom_name(atm3)==" N  ");
		if (hasCD && hasN)
			Ktheta = theta0 = 0.0;
	}

	if ( (restype.aa() == core::chemical::aa_ala || restype.aa() == core::chemical::aa_gly) &&
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
IdealParametersDatabase::lookup_length_legacy( core::pose::Pose const & pose, core::conformation::Residue const & res,
		int atm1, int atm2, Real &Kd, Real &d0 )
{
	using namespace core::chemical;

	Kd=k_length_;
	ResidueType const & restype = res.type();
	std::string restag = get_restag( restype );

	// Create mini-conformation for idealized residue
	core::conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );

	// get length
	numeric::xyzVector<core::Real> x,y,z; // atom coords
	if (atm1 > 0) {
		x = newres->atom( atm1 ).xyz();
	}
	else {
		x = newres->residue_connection( -atm1 ).icoor().build( restype );
	}
	if (atm2 > 0) {
		y = newres->atom( atm2 ).xyz();
	}
	else {
		y = newres->residue_connection( -atm2 ).icoor().build( restype );
	}

	d0 = (x-y).length();
	  
	if (restype.aa() == core::chemical::aa_pro &&
	    (atm1>0 && atm2>0) &&
	    ( ( restype.atom_name(atm1) == " CD " && restype.atom_name(atm2) == " N  ") ||
	      ( restype.atom_name(atm2) == " CD " && restype.atom_name(atm1) == " N  ") ) ) {
		Kd = d0 = 0.0;
	}

	return;
}



//////////////////////
/// EnergyMethod
CartesianBondedEnergy::CartesianBondedEnergy( methods::EnergyMethodOptions const & options ) :
	parent( new CartesianBondedEnergyCreator ) {
	// if flag _or_ energy method wants a linear potential, make the potential linear
	linear_bonded_potential_ = 
		basic::options::option[ basic::options::OptionKeys::score::linear_bonded_potential ]() ||
		options.get_cartesian_bonded_linear();

	// initialize databases
	core::Real cartbonded_len, cartbonded_ang, cartbonded_tors, cartbonded_proton , cartbonded_improper;
	options.get_cartesian_bonded_parameters( cartbonded_len, cartbonded_ang, cartbonded_tors, cartbonded_proton , cartbonded_improper );
	if (!db_)
		db_ = new IdealParametersDatabase(cartbonded_len, cartbonded_ang, cartbonded_tors, cartbonded_proton , cartbonded_improper);
}

CartesianBondedEnergy::CartesianBondedEnergy( CartesianBondedEnergy const & src ) : parent( src ) {
	linear_bonded_potential_ = src.linear_bonded_potential_;
}

CartesianBondedEnergy::~CartesianBondedEnergy() {}

EnergyMethodOP
CartesianBondedEnergy::clone() const {
	return new CartesianBondedEnergy( *this );
}


methods::LongRangeEnergyType
CartesianBondedEnergy::long_range_type() const { return methods::cart_bonded_lr; }

void
CartesianBondedEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const {
	using namespace methods;

	// create LR energy container
	LongRangeEnergyType const & lr_type( long_range_type() );
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		PeptideBondedEnergyContainerOP dec( static_cast< PeptideBondedEnergyContainer * > ( lrc.get() ) );
		Size nres = pose.total_residue();
		if( core::pose::symmetry::is_symmetric(pose) )
			nres = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		if ( dec->size() != nres ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		Size nres = pose.total_residue();
		if( core::pose::symmetry::is_symmetric(pose) )
			nres = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		TR << "Creating new peptide-bonded energy container (" << nres << ")" << std::endl;
		utility::vector1< ScoreType > s_types;
		s_types.push_back( cart_bonded );
		s_types.push_back( cart_bonded_angle );
		s_types.push_back( cart_bonded_length );
		s_types.push_back( cart_bonded_torsion );
		LREnergyContainerOP new_dec = new PeptideBondedEnergyContainer( nres, s_types );
		energies.set_long_range_container( lr_type, new_dec );
	}
}

bool
CartesianBondedEnergy::defines_residue_pair_energy(
	pose::Pose const &,
	Size res1,
	Size res2
) const {
	// is this fn. called?
	return ( res1 == (res2+1) || res1 == (res2-1) );
}


void
CartesianBondedEnergy::eval_improper_torsions(
   conformation::Residue const & rsd1,
	 conformation::Residue const & rsd2,
	 pose::Pose const & pose,
	 ScoreFunction const & ,
	 EnergyMap & emap
) const {
	using namespace core::chemical;
	using numeric::constants::f::pi;

	// backbone C-N-CA-H
	assert( rsd1.seqpos() < rsd2.seqpos() );
	conformation::Residue const &res_low = rsd1;
	conformation::Residue const &res_high = rsd2;

	core::Real energy_torsion = 0.0;

	if (!res_low.is_protein() || !res_high.is_protein()) return;

	if (res_high.aa() != aa_pro && res_high.has(" H  ") ) {  // check for terminal variants
		core::Size atm1 = res_high.atom_index(" CA ");
		core::Size atm2 = res_low.atom_index(" C  ");
		core::Size atm3 = res_high.atom_index(" N  ");
		core::Size atm4 = res_high.atom_index(" H  ");

		CartBondedParametersCAP tor_params = db_->lookup_torsion(res_low, "CA","C","N","H");
		runtime_assert( tor_params ); //fpd  fail if this is not in DB
		if (!tor_params->is_null()) {  // however, if it is in the db but with 0 weight, that is OK
			Real Kphi = tor_params->K(0,0), phi0=tor_params->mu(0,0), phi_step=2*pi/tor_params->period();
			Real angle = numeric::dihedral_radians
				( res_high.atom( atm1 ).xyz(), res_low.atom( atm2 ).xyz(), res_high.atom( atm3 ).xyz(), res_high.atom( atm4 ).xyz() );

			Real del_phi = basic::subtract_radian_angles(angle, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );
	
			if (pose.pdb_info() && 0.5*Kphi*del_phi*del_phi > CUTOFF) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << res_low.seqpos() << " pdbpos: " << pose.pdb_info()->number(res_low.seqpos()) << " improper torsion: " <<
					res_low.name() << " : " <<
					res_high.atom_name( atm1 ) << " , " << res_low.atom_name( atm2 ) << " , " <<
					res_high.atom_name( atm3 ) << " , " << res_high.atom_name( atm4 ) << "   (" <<
					Kphi << ") " << angle << " " << phi0 << "    sc=" << 0.5*Kphi*del_phi*del_phi << std::endl;
			}
	
			if (linear_bonded_potential_ && std::fabs(del_phi)>1) 
				energy_torsion += 0.5*Kphi*std::fabs(del_phi);
			else 
				energy_torsion += 0.5*Kphi*del_phi*del_phi;
		}
	}

	{
		core::Size atm1 = res_low.atom_index(" CA ");
		core::Size atm2 = res_high.atom_index(" N  ");
		core::Size atm3 = res_low.atom_index(" C  ");
		core::Size atm4 = res_low.atom_index(" O  ");

		CartBondedParametersCAP tor_params = db_->lookup_torsion(res_low, "CA", "N", "C", "O");
		runtime_assert( tor_params ); //fpd  fail if this is not in DB
		if (!tor_params->is_null()) {  // however, if it is in the db but with 0 weight, that is OK
			Real Kphi = tor_params->K(0,0), phi0=tor_params->mu(0,0), phi_step=2*pi/tor_params->period();
			Real angle = numeric::dihedral_radians
				( res_low.atom( atm1 ).xyz(), res_high.atom( atm2 ).xyz(), res_low.atom( atm3 ).xyz(), res_low.atom( atm4 ).xyz() );
	
			Real del_phi = basic::subtract_radian_angles(angle, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );
	
			if (pose.pdb_info() && 0.5*Kphi*del_phi*del_phi > CUTOFF) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << res_low.seqpos() << " pdbpos: " << pose.pdb_info()->number(res_low.seqpos()) << " improper torsion: " <<
					res_low.name() << " : " <<
					res_low.atom_name( atm1 ) << " , " << res_high.atom_name( atm2 ) << " , " <<
					res_low.atom_name( atm3 ) << " , " << res_low.atom_name( atm4 ) << "   (" <<
					Kphi << ") " << angle << " " << phi0 << "    sc="  << 0.5*Kphi*del_phi*del_phi << std::endl;
			}
	
			if (linear_bonded_potential_ && std::fabs(del_phi)>1) 
				energy_torsion += 0.5*Kphi*std::fabs(del_phi);
			else 
				energy_torsion += 0.5*Kphi*del_phi*del_phi;
		}
	}

	emap[ cart_bonded_torsion ] += energy_torsion;
	emap[ cart_bonded ] += energy_torsion; // potential double counting
}

void
CartesianBondedEnergy::eval_improper_torsions(
   conformation::Residue const & rsd,
	 pose::Pose const & pose,
	 ScoreFunction const & ,
	 EnergyMap & emap
) const {
	using namespace core::chemical;
	using numeric::constants::f::pi;

	if (!pose.is_fullatom()) return;

	std::string atm1, atm2, atm3, atm4;
	core::Size rt1, rt2, rt3, rt4;
	if (rsd.aa() == aa_asp) {
		// asp CB-CG-OD1-OD2
		atm1 = "CB"; atm2 = "CG"; atm3 = "OD1"; atm4 = "OD2";
	} else if (rsd.aa() == aa_glu) {
		// glu CG-CD-OE1-OE2
		atm1 = "CG"; atm2 = "CD"; atm3 = "OE1"; atm4 = "OE2";
	} else if (rsd.aa() == aa_asn) {
		// asn CB-CG-OD-ND
		atm1 = "CB"; atm2 = "CG"; atm3 = "OD1"; atm4 = "ND2";
	} else if (rsd.aa() == aa_gln) {
		// gln CG-CD-OE-NE
		atm1 = "CG"; atm2 = "CD"; atm3 = "OE1"; atm4 = "NE2";
	} else {
		return;
	}
	rt1 = rsd.atom_index(atm1);
	rt2 = rsd.atom_index(atm2);
	rt3 = rsd.atom_index(atm3);
	rt4 = rsd.atom_index(atm4);

	CartBondedParametersCAP tor_params = db_->lookup_torsion(rsd, atm1, atm2, atm3, atm4);
	runtime_assert( tor_params ); //fpd  fail if this is not in DB
	if (tor_params->is_null())
		return;  // however, if it is in the db but with 0 weight, that is OK

	Real Kphi = tor_params->K(0,0), phi0=tor_params->mu(0,0), phi_step=2*pi/tor_params->period();
	Real angle = numeric::dihedral_radians
		( rsd.atom( rt1 ).xyz(), rsd.atom( rt2 ).xyz(), rsd.atom( rt3 ).xyz(), rsd.atom( rt4 ).xyz() );
	Real del_phi = basic::subtract_radian_angles(angle, phi0);
	del_phi = basic::periodic_range( del_phi, phi_step );

	if (pose.pdb_info() && 0.5*Kphi*del_phi*del_phi > CUTOFF) {
		TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " << pose.pdb_info()->number(rsd.seqpos()) << " improper torsion: " <<
			rsd.name() << " : " <<
			atm1 << " , " << atm2 << " , " <<
			atm3 << " , " << atm4 << "   (" <<
			Kphi << ") " << angle << " " << phi0 << "    sc="  << 0.5*Kphi*del_phi*del_phi << std::endl;
	}

	core::Real energy_torsion;
	if (linear_bonded_potential_ && std::fabs(del_phi)>1) 
		energy_torsion = 0.5*Kphi*std::fabs(del_phi);
	else 
		energy_torsion = 0.5*Kphi*del_phi*del_phi;

	emap[ cart_bonded_torsion ] += energy_torsion;
	emap[ cart_bonded ] += energy_torsion; // potential double counting
}


// because of preproline potential, everything is done in res_pair
void
CartesianBondedEnergy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const
{
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
	if ( rsd1.seqpos() < rsd2.seqpos() ) {
		residue_pair_energy_sorted( rsd1, rsd2, pose, sf, emap );
	} else {
		residue_pair_energy_sorted( rsd2, rsd1, pose, sf, emap );
	}
}

void
CartesianBondedEnergy::residue_pair_energy_sorted(
   conformation::Residue const & rsd1,
	 conformation::Residue const & rsd2,
	 pose::Pose const & pose,
	 ScoreFunction const & sf,
	 EnergyMap & emap
) const {

	using namespace numeric;

	assert( rsd2.seqpos() > rsd1.seqpos() );

	core::Size resid = rsd1.seqpos();
	bool is_cterm = (pose.fold_tree().is_cutpoint( resid ));
	bool preproline = (!is_cterm && rsd2.aa() == core::chemical::aa_pro);

	if ( rsd1.aa() == core::chemical::aa_vrt) return;

	// get subcomponents
	eval_singleres_energy(rsd1, pose, sf, emap, preproline); // calls singleres improper

	// last residue won't ever be rsd1, so we need to explicitly call eval_singleres
	Size nres = pose.total_residue();
	if( core::pose::symmetry::is_symmetric(pose) ) nres = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	if (rsd2.seqpos() == nres) {
		eval_singleres_energy(rsd2, pose, sf, emap, false);
	}

	// bail out if the residues aren't bonded or we cross a cutpoint
	if (!rsd1.is_bonded(rsd2)) { return; }
	if ( is_cterm ) { return; }

	eval_improper_torsions( rsd1, rsd2, pose, sf, emap );

	Real energy_angle=0, energy_length=0;
	chemical::ResidueType const & rsd1_type = rsd1.type();
	chemical::ResidueType const & rsd2_type = rsd2.type();

	// get phi,psis for bb-dep angles/lengths
	Real phi1=0,psi1=0,phi2=0,psi2=0;
	if (rsd1.is_protein()) {
		phi1 = nonnegative_principal_angle_degrees( rsd1.mainchain_torsion(1));
		psi1 = nonnegative_principal_angle_degrees( rsd1.mainchain_torsion(2));
	}
	if (rsd2.is_protein()) {
		phi2 = nonnegative_principal_angle_degrees( rsd2.mainchain_torsion(1));
		psi2 = nonnegative_principal_angle_degrees( rsd2.mainchain_torsion(2));
	}

	utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );
	for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
		Size const resconn_id1( r1_resconn_ids[ii] );
		Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

		Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
		Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

		// compute bond-angle energies with two atoms on rsd1
		utility::vector1< chemical::two_atom_set > const & rsd1_atoms_wi1_bond_of_ii(
			rsd1_type.atoms_within_one_bond_of_a_residue_connection( resconn_id1 ));
		for ( Size jj = 1; jj <= rsd1_atoms_wi1_bond_of_ii.size(); ++jj ) {
			assert( rsd1_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno1 );
			Size const res1_lower_atomno = rsd1_atoms_wi1_bond_of_ii[ jj ].key2();

			Real const angle = numeric::angle_radians(
				rsd1.atom( res1_lower_atomno ).xyz(),
				rsd1.atom( resconn_atomno1 ).xyz(),
				rsd2.atom( resconn_atomno2 ).xyz() );
			std::string atm1name=rsd1.atom_name(res1_lower_atomno); boost::trim(atm1name);
			std::string atm2name=rsd1.atom_name(resconn_atomno1); boost::trim(atm2name);
			std::string atm3name=rsd2.atom_name(resconn_atomno2); boost::trim(atm3name);

			// lookup Ktheta and theta0
			CartBondedParametersCAP ang_params = db_->lookup_angle(rsd1, preproline,
				atm1name, atm2name, atm3name, res1_lower_atomno, resconn_atomno1, -resconn_id1);
			if (ang_params->is_null()) continue;
			Real Ktheta=ang_params->K(phi1,psi1), theta0=ang_params->mu(phi1,psi1);

			if (pose.pdb_info() && 0.5*Ktheta*(angle-theta0)*(angle-theta0) > CUTOFF) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd1.seqpos() 
				   << " pdbpos: " << pose.pdb_info()->number(rsd1.seqpos()) << " angle rsd1: " << rsd1.name() << ":" << rsd1_type.atom_name( res1_lower_atomno ) << " , "
				   << rsd1_type.atom_name( resconn_atomno1 ) << " , " << rsd2_type.atom_name( resconn_atomno2 )<< "   ("
				   << Ktheta << ")   " << angle << "  " << theta0
				   << "   sc=" << 0.5*Ktheta*(angle-theta0)*(angle-theta0) << std::endl;
			}

			// accumulate the energy
			if (linear_bonded_potential_ && std::fabs(angle-theta0)>1) {
				energy_angle += 0.5*Ktheta*std::fabs(angle-theta0);
			} else {
				energy_angle += 0.5*Ktheta*(angle-theta0) * (angle-theta0);
			}
		}

		/// compute bond-angle energies with two atoms on rsd2
		utility::vector1< chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
			rsd2_type.atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));
		for ( Size jj = 1; jj <= rsd2_atoms_wi1_bond_of_ii.size(); ++jj ) {
			assert( rsd2_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno2 );
			Size const res2_lower_atomno = rsd2_atoms_wi1_bond_of_ii[ jj ].key2();

			// lookup Ktheta and theta0
			// set pre-proline to false here -- it refers to rsd2+1 which would make this three-body
			std::string atm1name=rsd2.atom_name(res2_lower_atomno); boost::trim(atm1name);
			std::string atm2name=rsd2.atom_name(resconn_atomno2); boost::trim(atm2name);
			std::string atm3name=rsd1.atom_name(resconn_atomno1); boost::trim(atm3name);
			CartBondedParametersCAP ang_params = db_->lookup_angle(rsd2, false,
				atm1name, atm2name, atm3name, res2_lower_atomno, resconn_atomno2, -resconn_id2);
			if (ang_params->is_null()) continue;
			Real Ktheta=ang_params->K(phi2,psi2), theta0=ang_params->mu(phi2,psi2);

			if (Ktheta == 0.0) continue;
			Real const angle = numeric::angle_radians(
				rsd2.atom( res2_lower_atomno ).xyz(),
				rsd2.atom( resconn_atomno2 ).xyz(),
				rsd1.atom( resconn_atomno1 ).xyz() );

			if (pose.pdb_info() && 0.5*Ktheta*(angle-theta0)*(angle-theta0) > CUTOFF) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd2.seqpos() 
				   << " pdbpos: " << pose.pdb_info()->number(rsd2.seqpos()) << " angle rsd2: " << rsd2.name() << ":" << rsd2_type.atom_name( res2_lower_atomno ) << " , "
			     << rsd2_type.atom_name( resconn_atomno2 ) << " , " << rsd1_type.atom_name( resconn_atomno1 )<< "   ("
				   << Ktheta << ")   " << angle << "  " << theta0
				   << "   sc=" << 0.5*Ktheta*(angle-theta0)*(angle-theta0) << std::endl;
			}
			
			// accumulate the energy
			if (linear_bonded_potential_ && std::fabs(angle-theta0)>1) {
				energy_angle += 0.5*Ktheta*std::fabs(angle-theta0);
			} else {
				energy_angle += 0.5*Ktheta*(angle-theta0) * (angle-theta0);
			}
		}

		/////////////
		/// finally, compute the bondlength across the interface
		Real length =
			( rsd2.atom( resconn_atomno2 ).xyz() - rsd1.atom( resconn_atomno1 ).xyz() ).length();

		// lookup Kd and d0
		// again, pre-pro == false, definitions are based on pro as rsd2+1
		std::string atm1name=rsd2.atom_name(resconn_atomno2); boost::trim(atm1name);
		std::string atm2name=rsd1.atom_name(resconn_atomno1); boost::trim(atm2name);
		CartBondedParametersCAP len_params = db_->lookup_length(rsd2, false,
				atm1name, atm2name, resconn_atomno2, -resconn_id2);

		if (len_params->is_null()) continue;
		Real Kd=len_params->K(phi2,psi2), d0=len_params->mu(phi2,psi2);

		if (pose.pdb_info() && 0.5*Kd*(length-d0) * (length-d0) > CUTOFF) {
			TR.Debug << pose.pdb_info()->name() 
			   << " pdbpos rsd1: " << pose.pdb_info()->number(rsd1.seqpos()) << " length rsd1 rsd2: " << rsd1.seqpos() << " -- " << rsd2.seqpos() << "  "
			   << rsd1.name() << ":" << rsd1_type.atom_name( resconn_atomno1 ) << " , " << rsd2_type.atom_name( resconn_atomno2 )<< "   ("
			   << Kd << ")   " << length << "  " << d0
			   << "   sc=" << 0.5*Kd*(length-d0)*(length-d0) << std::endl;
		}

		// accumulate the energy
		if (linear_bonded_potential_ && std::fabs(length-d0)>1) {
			energy_length += 0.5*Kd*std::fabs(length-d0);
		} else {
			energy_length += 0.5*Kd*(length-d0)*(length-d0);
		}
	}

	emap[ cart_bonded_angle ] += energy_angle;
	emap[ cart_bonded_length ] += energy_length;

	// fpd note this is double counting if both cart_bonded and cart_bonded_* are set
	emap[ cart_bonded ] += energy_angle;
	emap[ cart_bonded ] += energy_length;
}

void
CartesianBondedEnergy::eval_singleres_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sf,
	EnergyMap & emap,
	bool preproline
) const {
	using numeric::constants::f::pi;
	using namespace numeric;

	Real energy_angle=0, energy_length=0, energy_torsion=0;	
	chemical::ResidueType const & rsd_type = rsd.type();

	if ( rsd_type.aa() == core::chemical::aa_vrt) return;

	eval_improper_torsions( rsd, pose, sf, emap );

	// phi/psi
	Real phi=0,psi=0;
	if (rsd.is_protein()) {
		phi = nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1));
		psi = nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2));
	}

	for ( Size dihe = 1; dihe <= rsd_type.ndihe(); ++dihe ){
		// get ResidueType ints
		int rt1 = ( rsd_type.dihedral( dihe ) ).key1();
		int rt2 = ( rsd_type.dihedral( dihe ) ).key2();
		int rt3 = ( rsd_type.dihedral( dihe ) ).key3();
		int rt4 = ( rsd_type.dihedral( dihe ) ).key4();

		// lookup Kphi and phi0
		std::string atm1name=rsd.atom_name(rt1); boost::trim(atm1name);
		std::string atm2name=rsd.atom_name(rt2); boost::trim(atm2name);
		std::string atm3name=rsd.atom_name(rt3); boost::trim(atm3name);
		std::string atm4name=rsd.atom_name(rt4); boost::trim(atm4name);
		CartBondedParametersCAP tor_params = db_->lookup_torsion(rsd,
			atm1name,atm2name,atm3name,atm4name );

		if (!tor_params || tor_params->is_null()) continue;

		Real Kphi=tor_params->K(phi,psi), phi0=tor_params->mu(phi,psi), phi_step=2*pi/tor_params->period();

		// get angle
		Real angle = numeric::dihedral_radians
			( rsd.atom( rt1 ).xyz(), rsd.atom( rt2 ).xyz(),
				rsd.atom( rt3 ).xyz(), rsd.atom( rt4 ).xyz() );

		// accumulate the energy
		Real del_phi = basic::subtract_radian_angles(angle, phi0);
		if (phi_step>0) del_phi = basic::periodic_range( del_phi, phi_step );

 		if (pose.pdb_info() && 0.5*Kphi*del_phi*del_phi > CUTOFF) {
 			TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " << pose.pdb_info()->number(rsd.seqpos()) << " intrares torsion: " <<
 				rsd_type.name() << " : " <<
 			  rsd.atom_name( rt1 ) << " , " << rsd.atom_name( rt2 ) << " , " <<
 			  rsd.atom_name( rt3 ) << " , " << rsd.atom_name( rt4 )<< "   ("
			  << Kphi << ")   " << angle << "  " << phi0
			  << "   sc=" << 0.5*Kphi*del_phi*del_phi << std::endl;
 		}

		if (linear_bonded_potential_ && std::fabs(del_phi)>1) 
			energy_torsion += 0.5*Kphi*std::fabs(del_phi);
		else 
			energy_torsion += 0.5*Kphi*del_phi*del_phi;
	}

	// for each angle in the residue
	for ( Size bondang = 1; bondang <= rsd_type.num_bondangles(); ++bondang ) {
		// get ResidueType ints
		Size rt1 = ( rsd_type.bondangle( bondang ) ).key1();
		Size rt2 = ( rsd_type.bondangle( bondang ) ).key2();
		Size rt3 = ( rsd_type.bondangle( bondang ) ).key3();

		// lookup Ktheta and theta0
		std::string atm1name=rsd.atom_name(rt1); boost::trim(atm1name);
		std::string atm2name=rsd.atom_name(rt2); boost::trim(atm2name);
		std::string atm3name=rsd.atom_name(rt3); boost::trim(atm3name);
		CartBondedParametersCAP ang_params = db_->lookup_angle(rsd, preproline,
			atm1name,atm2name,atm3name, rt1,rt2,rt3 );

		if (ang_params->is_null()) continue;
		Real Ktheta=ang_params->K(phi,psi), theta0=ang_params->mu(phi,psi);

		// get angle
		Real const angle = numeric::angle_radians(
			rsd.atom( rt1 ).xyz(),
			rsd.atom( rt2 ).xyz(),
			rsd.atom( rt3 ).xyz() );

 		if ( pose.pdb_info() && 0.5*Ktheta*(angle-theta0) * (angle-theta0) > CUTOFF) {
 			TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " << pose.pdb_info()->number(rsd.seqpos()) << " intrares angle: " <<
 				rsd_type.name() << " : " <<
 			  rsd.atom_name( rt1 ) << " , " << rsd.atom_name( rt2 ) << " , " <<
 			  rsd.atom_name( rt3 ) << "   ("
			  << Ktheta << ")   " << angle << "  " << theta0
			  << "   sc=" << 0.5*Ktheta*(angle-theta0) * (angle-theta0) << std::endl;
 		}

		// accumulate the energy
		if (linear_bonded_potential_ && std::fabs(angle - theta0)>1) {
			energy_angle += 0.5*Ktheta*std::fabs(angle-theta0);
		} else {
			energy_angle += 0.5*Ktheta*(angle-theta0) * (angle-theta0);
		}
	}

	// for each bond
	for (Size atm_i=1; atm_i<=rsd_type.natoms(); ++atm_i) {
		chemical::AtomIndices atm_nbrs = rsd_type.nbrs( atm_i );
		for (Size j=1; j<=atm_nbrs.size(); ++j) {
			Size atm_j = atm_nbrs[j];
			if ( atm_i<atm_j ) { // only score each bond once -- use restype index to define ordering

				// lookup Kd and d0
				std::string atm1name=rsd.atom_name(atm_i); boost::trim(atm1name);
				std::string atm2name=rsd.atom_name(atm_j); boost::trim(atm2name);
				CartBondedParametersCAP len_params = db_->lookup_length(rsd, preproline,
					atm1name,atm2name, atm_i,atm_j );
				if (len_params->is_null()) continue;

				Real Kd=len_params->K(phi,psi), d0=len_params->mu(phi,psi);
				Real const d = ( rsd.xyz( atm_i )-rsd.xyz( atm_j ) ).length();

				if ( pose.pdb_info() && 0.5*Kd*(d-d0)*(d-d0) > CUTOFF) {
					TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " << pose.pdb_info()->number(rsd.seqpos()) << " intrares angle: " <<
						rsd_type.name() << " : " <<
						rsd.atom_name( atm_i ) << " , " << rsd.atom_name( atm_j ) << "   ("
						<< Kd << ")   " << d << "  " << d0
						<< "   sc=" << 0.5*Kd*(d-d0)*(d-d0) << std::endl;
				}

				// accumulate the energy
				if (linear_bonded_potential_ && std::fabs(d - d0)>1) {
					energy_length += 0.5*Kd*std::fabs(d-d0);
				} else {
					energy_length += 0.5*Kd*(d-d0)*(d-d0);
				}
			}
		}
	}

	// add energy to emap
	emap[ cart_bonded_angle ] += energy_angle;
	emap[ cart_bonded_length ] += energy_length;
	emap[ cart_bonded_torsion ] += energy_torsion;

	// fpd note this is dounble counting if both cart_bonded and cart_bonded_* are set
	emap[ cart_bonded ] += energy_angle;
	emap[ cart_bonded ] += energy_length;
	emap[ cart_bonded ] += energy_torsion; 
}


// dof (bbdep) derivatives
Real
CartesianBondedEnergy::eval_intraresidue_dof_derivative(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		id::DOF_ID const & dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights
) const {
	using namespace numeric;
	using numeric::constants::f::pi;

	// save some time if we're only doing bb-indep
	if (!db_->bbdep_bond_params()) return 0.0;

	if ( !tor_id.valid() || tor_id.type()!=id::BB || tor_id.torsion() > 2  || !rsd.is_protein() ) 
		return 0.0;

	core::Size resid = rsd.seqpos();
	bool is_nterm = ((resid==1) || pose.fold_tree().is_cutpoint( resid-1 ));
	bool is_cterm = ((resid==pose.total_residue()) || pose.fold_tree().is_cutpoint( resid ));
	bool preproline = (!is_cterm && pose.residue( resid+1 ).aa() == core::chemical::aa_pro);

	// phi/psi
	Real phi=0,psi=0;
	if (rsd.is_protein()) {
		phi = nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1));
		psi = nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2));
	}

	core::Real deriv=0.0;

	// loop over all bb-dep angles and bonds, summing dmu_dphi and dK_dphi
	// hardcode the bbdep angles and lengths to save a bit of time
	// connection IDs are hardcoded ... is this a problem?
	for (int i=1; i<=5; ++i) {
		if (i==1 && is_nterm) continue;
		if (i==3 && rsd.aa() == core::chemical::aa_gly) continue;

		std::string atm1,atm2;
		core::Size rt1,rt2;
		if (i==1) { atm1="C"; atm2="N"; rt1=-1; rt2=rsd.atom_index(" N  ");}
		if (i==2) { atm1="N"; atm2="CA"; rt1=rsd.atom_index(" N  "); rt2=rsd.atom_index(" CA "); }
		if (i==3) { atm1="CA"; atm2="CB"; rt1=rsd.atom_index(" CA "); rt2=rsd.atom_index(" CB "); }
		if (i==4) { atm1="CA"; atm2="C"; rt1=rsd.atom_index(" CA "); rt2=rsd.atom_index(" C  "); }
		if (i==5) { atm1="C"; atm2="O"; rt1=rsd.atom_index(" C  "); rt2=rsd.atom_index(" O  "); }

		Vector atom1xyz, atom2xyz = rsd.atom( rt2 ).xyz();
		if (i==1)
			atom1xyz = pose.residue( resid-1 ).atom(" C  ").xyz();
		else 
			atom1xyz = rsd.atom( rt1 ).xyz();
		CartBondedParametersCAP len_params = db_->lookup_length(rsd, (i==1)?false:preproline, atm1,atm2, rt1,rt2 );

		Real const d = ( atom2xyz-atom1xyz ).length();
		core::Real mu,K, dmu_dtor, dK_dtor=0.0;
		core::Real dscore_dK=0.0, dscore_dmu;

		K = len_params->K(phi,psi);
		mu = len_params->mu(phi,psi);

		if (tor_id.torsion()==1) {
			dmu_dtor = len_params->dmu_dphi(phi,psi);
		} else {
			dmu_dtor = len_params->dmu_dpsi(phi,psi);
		}
		if (linear_bonded_potential_ && std::fabs(d - mu)>1) {
			dscore_dmu = -K * ((d - mu)>0 ? 1:-1);
		} else {
			dscore_dmu = -K * (d - mu);
		}
		if (db_->bbdep_bond_devs()) {
			if (tor_id.torsion()==1) {
				dK_dtor = len_params->dK_dphi(phi,psi);
			} else {
				dK_dtor = len_params->dK_dpsi(phi,psi);
			}
			if (linear_bonded_potential_ && std::fabs(d - mu)>1) {
				dscore_dK = 0.5*std::fabs(d - mu);
			} else {
				dscore_dK = 0.5*(d - mu)*(d - mu);   // currently we have no -log(K) term in our score
			}
		}

		// derivatives w.r.t. phi/psi
		deriv += (weights[ cart_bonded_length ] + weights[ cart_bonded ]) * (dscore_dmu*dmu_dtor + dscore_dK*dK_dtor);
	}

	for (int i=1; i<=7; ++i) {
		if ((i==2 || i==4) && rsd.aa() == core::chemical::aa_gly) continue;
		if (i==1 && is_nterm) continue;
		if ((i==6 || i==7) && is_cterm) continue;

		std::string atm1,atm2,atm3;
		core::Size rt1,rt2,rt3;
		if (i==1) { atm1="C"; atm2="N";atm3="CA"; rt1=-1; rt2=rsd.atom_index(" N  "); rt3=rsd.atom_index(" CA ");}
		if (i==2) { atm1="N"; atm2="CA";atm3="CB"; rt1=rsd.atom_index(" N  "); rt2=rsd.atom_index(" CA "); rt3=rsd.atom_index(" CB "); }
		if (i==3) { atm1="N"; atm2="CA";atm3="C"; rt1=rsd.atom_index(" N  "); rt2=rsd.atom_index(" CA "); rt3=rsd.atom_index(" C  "); }
		if (i==4) { atm1="CB"; atm2="CA";atm3="C"; rt1=rsd.atom_index(" CB "); rt2=rsd.atom_index(" CA "); rt3=rsd.atom_index(" C  "); }
		if (i==5) { atm1="CA"; atm2="C";atm3="O"; rt1=rsd.atom_index(" CA "); rt2=rsd.atom_index(" C  "); rt3=rsd.atom_index(" O  "); }
		if (i==6) { atm1="CA"; atm2="C";atm3="N"; rt1=rsd.atom_index(" CA "); rt2=rsd.atom_index(" C  "); rt3=-2; }
		if (i==7) { atm1="O"; atm2="C";atm3="N"; rt1=rsd.atom_index(" O  "); rt2=rsd.atom_index(" C  "); rt3=-2; }

		Vector atom1xyz, atom2xyz = rsd.atom( rt2 ).xyz(), atom3xyz;
		if (i==1)
			atom1xyz = pose.residue( resid-1 ).atom(" C  ").xyz();
		else 
			atom1xyz = rsd.atom( rt1 ).xyz();
		if (i==6 || i==7)
			atom3xyz = pose.residue( resid+1 ).atom(" N  ").xyz();
		else 
			atom3xyz = rsd.atom( rt3 ).xyz();

		CartBondedParametersCAP ang_params = db_->lookup_angle(rsd, preproline, atm1,atm2,atm3, rt1,rt2,rt3 );

		Real const angle = numeric::angle_radians(atom1xyz,atom2xyz,atom3xyz);

		core::Real mu,K, dmu_dtor, dK_dtor=0.0;
		core::Real dscore_dK=0.0, dscore_dmu;

		K = ang_params->K(phi,psi);
		mu = ang_params->mu(phi,psi);

		if (tor_id.torsion()==1) {
			dmu_dtor = ang_params->dmu_dphi(phi,psi);
		} else {
			dmu_dtor = ang_params->dmu_dpsi(phi,psi);
		}
		if (linear_bonded_potential_ && std::fabs(angle - mu)>1) {
			dscore_dmu = -K * ((angle - mu)>0 ? 1:-1);
		} else {
			dscore_dmu = -K * (angle - mu);
		}

		if (db_->bbdep_bond_devs()) {
			if (tor_id.torsion()==1) {
				dK_dtor = ang_params->dK_dphi(phi,psi);
			} else {
				dK_dtor = ang_params->dK_dpsi(phi,psi);
			}
			if (linear_bonded_potential_ && std::fabs(angle - mu)>1) {
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


///
void
CartesianBondedEnergy::eval_singleres_derivatives(
			conformation::Residue const & rsd,
			EnergyMap const & weights,
			utility::vector1< DerivVectorPair > & r_atom_derivs,
			bool preproline
) const {
	using numeric::constants::f::pi;
	using namespace numeric;

	chemical::ResidueType const & rsd_type = rsd.type();
	if ( rsd_type.aa() == core::chemical::aa_vrt) return;

	// impropers
	eval_improper_torsions_derivative( rsd, weights, r_atom_derivs );

	// phi/psi
	Real phi=0,psi=0;
	if (rsd.is_protein()) {
		phi = nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1));
		psi = nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2));
	}

	// [1] torsions
	for ( Size dihe = 1; dihe <= rsd_type.ndihe(); ++dihe ){
		// get ResidueType ints
		int rt1 = ( rsd_type.dihedral( dihe ) ).key1();
		int rt2 = ( rsd_type.dihedral( dihe ) ).key2();
		int rt3 = ( rsd_type.dihedral( dihe ) ).key3();
		int rt4 = ( rsd_type.dihedral( dihe ) ).key4();

		// lookup Kphi and phi0
		std::string atm1name=rsd.atom_name(rt1); boost::trim(atm1name);
		std::string atm2name=rsd.atom_name(rt2); boost::trim(atm2name);
		std::string atm3name=rsd.atom_name(rt3); boost::trim(atm3name);
		std::string atm4name=rsd.atom_name(rt4); boost::trim(atm4name);
		CartBondedParametersCAP tor_params = db_->lookup_torsion(rsd,
			atm1name,atm2name,atm3name,atm4name );

		if (!tor_params || tor_params->is_null()) continue;

		Real Kphi=tor_params->K(phi,psi), phi0=tor_params->mu(phi,psi), phi_step=2*pi/tor_params->period();

		Vector f1(0.0), f2(0.0);
		Real angle(0.0), dE_dphi;

		numeric::deriv::dihedral_p1_cosine_deriv( rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), rsd.xyz( rt4 ), angle, f1, f2 );
		Real del_phi = basic::subtract_radian_angles(angle, phi0);
		if (phi_step>0) del_phi = basic::periodic_range( del_phi, phi_step );
		if (linear_bonded_potential_ && std::fabs(del_phi)>1)
			dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * (del_phi>0? 1 : -1);
		else
			dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * del_phi;
		r_atom_derivs[ rt1 ].f1() += dE_dphi * f1;
		r_atom_derivs[ rt1 ].f2() += dE_dphi * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv( rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), rsd.xyz( rt4 ), angle, f1, f2 );
		r_atom_derivs[ rt2 ].f1() += dE_dphi * f1;
		r_atom_derivs[ rt2 ].f2() += dE_dphi * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p2_cosine_deriv( rsd.xyz( rt4 ), rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), angle, f1, f2 );
		del_phi = basic::subtract_radian_angles(angle, phi0);
		if (phi_step>0) del_phi = basic::periodic_range( del_phi, phi_step );
		if (linear_bonded_potential_ && std::fabs(del_phi)>1)
			dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * (del_phi>0? 1 : -1);
		else
			dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * del_phi;
		r_atom_derivs[ rt3 ].f1() += dE_dphi * f1;
		r_atom_derivs[ rt3 ].f2() += dE_dphi * f2;

		f1 = f2 = Vector(0.0);
		numeric::deriv::dihedral_p1_cosine_deriv( rsd.xyz( rt4 ), rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), angle, f1, f2 );
		r_atom_derivs[ rt4 ].f1() += dE_dphi * f1;
		r_atom_derivs[ rt4 ].f2() += dE_dphi * f2;
	}


	// [2] angles
	for ( Size bondang = 1; bondang <= rsd_type.num_bondangles(); ++bondang ) {
		// get ResidueType ints
		Size rt1 = ( rsd_type.bondangle( bondang ) ).key1();
		Size rt2 = ( rsd_type.bondangle( bondang ) ).key2();
		Size rt3 = ( rsd_type.bondangle( bondang ) ).key3();

		// lookup Ktheta and theta0
		std::string atm1name=rsd.atom_name(rt1); boost::trim(atm1name);
		std::string atm2name=rsd.atom_name(rt2); boost::trim(atm2name);
		std::string atm3name=rsd.atom_name(rt3); boost::trim(atm3name);
		CartBondedParametersCAP ang_params = db_->lookup_angle(rsd, preproline,
			atm1name,atm2name,atm3name, rt1,rt2,rt3 );

		if (ang_params->is_null()) continue;

		Real Ktheta=ang_params->K(phi,psi), theta0=ang_params->mu(phi,psi);

		Vector f1(0.0), f2(0.0);
		Real theta(0.0), dE_dtheta;

		numeric::deriv::angle_p1_deriv( rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), theta, f1, f2 );
		if (linear_bonded_potential_ && std::fabs(theta - theta0)>1)
			dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * ((theta - theta0)>0? 1 : -1);
		else
			dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * (theta - theta0);
		r_atom_derivs[ rt1 ].f1() += dE_dtheta * f1;
		r_atom_derivs[ rt1 ].f2() += dE_dtheta * f2;

		numeric::deriv::angle_p2_deriv( rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), theta, f1, f2 );
		r_atom_derivs[ rt2 ].f1() += dE_dtheta * f1;
		r_atom_derivs[ rt2 ].f2() += dE_dtheta * f2;

		numeric::deriv::angle_p1_deriv( rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), theta, f1, f2 );
		r_atom_derivs[ rt3 ].f1() += dE_dtheta * f1;
		r_atom_derivs[ rt3 ].f2() += dE_dtheta * f2;
	}

	// [3] bonds
	for (Size atm_i=1; atm_i<=rsd_type.natoms(); ++atm_i) {
		chemical::AtomIndices atm_nbrs = rsd_type.nbrs( atm_i );
		for (Size j=1; j<=atm_nbrs.size(); ++j) {
			Size atm_j = atm_nbrs[j];
			if ( atm_i<atm_j ) { // only score each bond once -- use restype index to define ordering
				// lookup Kd and d0
				std::string atm1name=rsd.atom_name(atm_i); boost::trim(atm1name);
				std::string atm2name=rsd.atom_name(atm_j); boost::trim(atm2name);
				CartBondedParametersCAP len_params = db_->lookup_length(rsd, preproline,
					atm1name,atm2name, atm_i,atm_j );
				if (len_params->is_null()) continue;

				Real Kd=len_params->K(phi,psi), d0=len_params->mu(phi,psi);

				Vector f1(0.0), f2(0.0);
				Real d=0, dE_dd;

				numeric::deriv::distance_f1_f2_deriv( rsd.xyz( atm_i ), rsd.xyz( atm_j ), d, f1, f2 );
				if (linear_bonded_potential_ && std::fabs(d - d0)>1)
					dE_dd = (weights[ cart_bonded_length ] + weights[ cart_bonded ]) * Kd * ((d - d0)>0? 1 : -1);
				else
					dE_dd = (weights[ cart_bonded_length ] + weights[ cart_bonded ]) * Kd * (d - d0);
				r_atom_derivs[ atm_i ].f1() += dE_dd * f1;
				r_atom_derivs[ atm_i ].f2() += dE_dd * f2;

				numeric::deriv::distance_f1_f2_deriv( rsd.xyz( atm_j ), rsd.xyz( atm_i ), d, f1, f2 );
				r_atom_derivs[ atm_j ].f1() += dE_dd * f1;
				r_atom_derivs[ atm_j ].f2() += dE_dd * f2;
			}
		}
	}
}


///////////////////////////
///
void
CartesianBondedEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	using namespace numeric;

	assert( rsd2.seqpos() > rsd1.seqpos() );
	bool preproline = (rsd2.aa()==core::chemical::aa_pro);

	chemical::ResidueType const & rsd1_type = rsd1.type();
	chemical::ResidueType const & rsd2_type = rsd2.type();

	if ( rsd1.aa() == core::chemical::aa_vrt) return;

	// get subcomponents
	eval_singleres_derivatives( rsd1, weights, r1_atom_derivs, preproline );

	// cterm special case
	Size nres = pose.total_residue();
	if( core::pose::symmetry::is_symmetric(pose) ) nres = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	if (rsd2.seqpos() == nres) {
		eval_singleres_derivatives(rsd2, weights, r2_atom_derivs, false );
	}

	// bail out if the residues aren't bonded or we cross a cutpoint
	if (!rsd1.is_bonded(rsd2)) { return; }
	if ( pose.fold_tree().is_cutpoint( rsd1.seqpos() ) ) { return; }

	eval_improper_torsions_derivative( rsd1, rsd2, weights, r1_atom_derivs, r2_atom_derivs );

	// get phi,psis for bb-dep angles/lengths
	Real phi1=0,psi1=0,phi2=0,psi2=0;
	if (rsd1.is_protein()) {
		phi1 = nonnegative_principal_angle_degrees( rsd1.mainchain_torsion(1));
		psi1 = nonnegative_principal_angle_degrees( rsd1.mainchain_torsion(2));
	}
	if (rsd2.is_protein()) {
		phi2 = nonnegative_principal_angle_degrees( rsd2.mainchain_torsion(1));
		psi2 = nonnegative_principal_angle_degrees( rsd2.mainchain_torsion(2));
	}

	utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );
	for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
		Size const resconn_id1( r1_resconn_ids[ii] );
		Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

		Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
		Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

		/// [1] angles involving two atoms on this rsd1
		utility::vector1< chemical::two_atom_set > const & rsd1_atoms_wi1_bond_of_ii(
			rsd1_type.atoms_within_one_bond_of_a_residue_connection( resconn_id1 ));
		for ( Size jj = 1; jj <= rsd1_atoms_wi1_bond_of_ii.size(); ++jj ) {
			assert( rsd1_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno1 );
			Size const res1_lower_atomno = rsd1_atoms_wi1_bond_of_ii[ jj ].key2();

			Real const angle = numeric::angle_radians(
				rsd1.atom( res1_lower_atomno ).xyz(),
				rsd1.atom( resconn_atomno1 ).xyz(),
				rsd2.atom( resconn_atomno2 ).xyz() );
			std::string atm1name=rsd1.atom_name(res1_lower_atomno); boost::trim(atm1name);
			std::string atm2name=rsd1.atom_name(resconn_atomno1); boost::trim(atm2name);
			std::string atm3name=rsd2.atom_name(resconn_atomno2); boost::trim(atm3name);

			// lookup Ktheta and theta0
			CartBondedParametersCAP ang_params = db_->lookup_angle(rsd1, preproline,
				atm1name, atm2name, atm3name, res1_lower_atomno, resconn_atomno1, -resconn_id1);
			if (ang_params->is_null()) continue;
			Real Ktheta=ang_params->K(phi1,psi1), theta0=ang_params->mu(phi1,psi1);

			Vector f1(0.0), f2(0.0);
			Real theta(0.0), dE_dtheta;
			numeric::deriv::angle_p1_deriv( 
			   rsd1.xyz( res1_lower_atomno ), rsd1.xyz( resconn_atomno1 ), rsd2.xyz( resconn_atomno2 ), theta, f1, f2 );
			if (linear_bonded_potential_ && std::fabs(theta - theta0)>1)
				dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * ((theta - theta0)>0? 1 : -1);
			else
				dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * (theta - theta0);
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


		/// [2] angles involving two atoms on this rsd2
		utility::vector1< chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
			rsd2_type.atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));
		for ( Size jj = 1; jj <= rsd2_atoms_wi1_bond_of_ii.size(); ++jj ) {
			assert( rsd2_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno2 );
			Size const res2_lower_atomno = rsd2_atoms_wi1_bond_of_ii[ jj ].key2();

			// lookup Ktheta and theta0
			// set pre-proline to false here -- it refers to rsd2+1 which would make this three-body
			std::string atm1name=rsd2.atom_name(res2_lower_atomno); boost::trim(atm1name);
			std::string atm2name=rsd2.atom_name(resconn_atomno2); boost::trim(atm2name);
			std::string atm3name=rsd1.atom_name(resconn_atomno1); boost::trim(atm3name);
			CartBondedParametersCAP ang_params = db_->lookup_angle(rsd2, false,
				atm1name, atm2name, atm3name, res2_lower_atomno, resconn_atomno2, -resconn_id2);
			if (ang_params->is_null()) continue;
			Real Ktheta=ang_params->K(phi2,psi2), theta0=ang_params->mu(phi2,psi2);

			Vector f1(0.0), f2(0.0);
			Real theta(0.0), dE_dtheta;
			numeric::deriv::angle_p1_deriv( 
			   rsd2.xyz( res2_lower_atomno ), rsd2.xyz( resconn_atomno2 ), rsd1.xyz( resconn_atomno1 ), theta, f1, f2 );
			if (linear_bonded_potential_ && std::fabs(theta - theta0)>1)
				dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * ((theta - theta0)>0? 1 : -1);
			else
				dE_dtheta = (weights[ cart_bonded_angle ] + weights[ cart_bonded ]) * Ktheta * (theta - theta0);
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

		// [3] Bonds across connection
		std::string atm1name=rsd2.atom_name(resconn_atomno2); boost::trim(atm1name);
		std::string atm2name=rsd1.atom_name(resconn_atomno1); boost::trim(atm2name);
		CartBondedParametersCAP len_params = db_->lookup_length(rsd2, false,
				atm1name, atm2name, resconn_atomno2, -resconn_id2);

		if (len_params->is_null()) continue;
		Real Kd=len_params->K(phi2,psi2), d0=len_params->mu(phi2,psi2);

		Vector f1(0.0), f2(0.0);
		Real d=0, dE_dd;

		numeric::deriv::distance_f1_f2_deriv( rsd2.xyz( resconn_atomno2 ), rsd1.xyz( resconn_atomno1 ), d, f1, f2 );
		if (linear_bonded_potential_ && std::fabs(d - d0)>1)
			dE_dd = (weights[ cart_bonded_length ] + weights[ cart_bonded ]) * Kd * ((d - d0)>0? 1 : -1);
		else
			dE_dd = (weights[ cart_bonded_length ] + weights[ cart_bonded ]) * Kd * (d - d0);
		r2_atom_derivs[ resconn_atomno2 ].f1() += dE_dd * f1;
		r2_atom_derivs[ resconn_atomno2 ].f2() += dE_dd * f2;

		numeric::deriv::distance_f1_f2_deriv( rsd1.xyz( resconn_atomno1 ), rsd2.xyz( resconn_atomno2 ), d, f1, f2 );
		r1_atom_derivs[ resconn_atomno1 ].f1() += dE_dd * f1;
		r1_atom_derivs[ resconn_atomno1 ].f2() += dE_dd * f2;
	}
}


// deriv impropers
void
CartesianBondedEnergy::eval_improper_torsions_derivative(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {
	using namespace core::chemical;
	using numeric::constants::f::pi;

	assert ( res1.seqpos() < res2.seqpos() );

	// backbone C-N-CA-H
	if (!res1.is_protein() || !res2.is_protein()) return;

	if (res2.aa() != aa_pro && res2.has(" H  ") ) {  // check for terminal variants
		core::Size atm1 = res2.atom_index(" CA ");
		core::Size atm2 = res1.atom_index(" C  ");
		core::Size atm3 = res2.atom_index(" N  ");
		core::Size atm4 = res2.atom_index(" H  ");

		CartBondedParametersCAP tor_params = db_->lookup_torsion(res1, "CA","C","N","H");
		runtime_assert( tor_params ); //fpd  fail if this is not in DB
		if (!tor_params->is_null()) {  // however, if it is in the db but with 0 weight, that is OK
			Real Kphi = tor_params->K(0,0), phi0=tor_params->mu(0,0), phi_step=2*pi/tor_params->period();
	
			Vector f1(0.0), f2(0.0);
			Real phi=0, dE_dphi;

			numeric::deriv::dihedral_p1_cosine_deriv(
				res2.xyz( atm1 ), res1.xyz( atm2 ), res2.xyz( atm3 ), res2.xyz( atm4 ), phi, f1, f2 );
			Real del_phi = basic::subtract_radian_angles(phi, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );
			if (linear_bonded_potential_ && std::fabs(del_phi)>1)
				dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * (del_phi>0? 1 : -1);
			else
				dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * del_phi;
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
			del_phi = basic::subtract_radian_angles(phi, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );
			if (linear_bonded_potential_ && std::fabs(del_phi)>1)
				dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * (del_phi>0? 1 : -1);
			else
				dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * del_phi;
			r2_atom_derivs[ atm3 ].f1() += dE_dphi * f1;
			r2_atom_derivs[ atm3 ].f2() += dE_dphi * f2;
	
			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv( 
				res2.xyz( atm4 ), res2.xyz( atm3 ), res1.xyz( atm2 ), res2.xyz( atm1 ), phi, f1, f2 );
			r2_atom_derivs[ atm4 ].f1() += dE_dphi * f1;
			r2_atom_derivs[ atm4 ].f2() += dE_dphi * f2;
		}
	}

	{
		core::Size atm1 = res1.atom_index(" CA ");
		core::Size atm2 = res2.atom_index(" N  ");
		core::Size atm3 = res1.atom_index(" C  ");
		core::Size atm4 = res1.atom_index(" O  ");

		CartBondedParametersCAP tor_params = db_->lookup_torsion(res1, "CA", "N", "C", "O");
		runtime_assert( tor_params ); //fpd  fail if this is not in DB
		if (!tor_params->is_null()) {  // however, if it is in the db but with 0 weight, that is OK
			Real Kphi = tor_params->K(0,0), phi0=tor_params->mu(0,0), phi_step=2*pi/tor_params->period();
	
			Vector f1(0.0), f2(0.0);
			Real phi=0, dE_dphi;

			numeric::deriv::dihedral_p1_cosine_deriv(
				res1.xyz( atm1 ), res2.xyz( atm2 ), res1.xyz( atm3 ), res1.xyz( atm4 ), phi, f1, f2 );
			Real del_phi = basic::subtract_radian_angles(phi, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );
			if (linear_bonded_potential_ && std::fabs(del_phi)>1)
				dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * (del_phi>0? 1 : -1);
			else
				dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * del_phi;
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
			del_phi = basic::subtract_radian_angles(phi, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );
			if (linear_bonded_potential_ && std::fabs(del_phi)>1)
				dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * (del_phi>0? 1 : -1);
			else
				dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * del_phi;
			r1_atom_derivs[ atm3 ].f1() += dE_dphi * f1;
			r1_atom_derivs[ atm3 ].f2() += dE_dphi * f2;
	
			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv( 
				res1.xyz( atm4 ), res1.xyz( atm3 ), res2.xyz( atm2 ), res1.xyz( atm1 ), phi, f1, f2 );
			r1_atom_derivs[ atm4 ].f1() += dE_dphi * f1;
			r1_atom_derivs[ atm4 ].f2() += dE_dphi * f2;
		}
	}
}


// deriv impropers
void
CartesianBondedEnergy::eval_improper_torsions_derivative(
	conformation::Residue const & rsd,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r_atom_derivs
) const {
	using namespace core::chemical;
	using numeric::constants::f::pi;

	if (rsd.type().residue_type_set().name() != core::chemical::FA_STANDARD) return;

	std::string atm1, atm2, atm3, atm4;
	core::Size rt1, rt2, rt3, rt4;
	if (rsd.aa() == aa_asp) {
		// asp CB-CG-OD1-OD2
		atm1 = "CB"; atm2 = "CG"; atm3 = "OD1"; atm4 = "OD2";
	} else if (rsd.aa() == aa_glu) {
		// glu CG-CD-OE1-OE2
		atm1 = "CG"; atm2 = "CD"; atm3 = "OE1"; atm4 = "OE2";
	} else if (rsd.aa() == aa_asn) {
		// asn CB-CG-OD-ND
		atm1 = "CB"; atm2 = "CG"; atm3 = "OD1"; atm4 = "ND2";
	} else if (rsd.aa() == aa_gln) {
		// gln CG-CD-OE-NE
		atm1 = "CG"; atm2 = "CD"; atm3 = "OE1"; atm4 = "NE2";
	} else {
		return;
	}
	rt1 = rsd.atom_index(atm1);
	rt2 = rsd.atom_index(atm2);
	rt3 = rsd.atom_index(atm3);
	rt4 = rsd.atom_index(atm4);

	Vector f1(0.0), f2(0.0);
	Real phi=0, dE_dphi;

	CartBondedParametersCAP tor_params = db_->lookup_torsion(rsd, atm1, atm2, atm3, atm4);
	runtime_assert( tor_params ); //fpd  fail if this is not in DB
	if (tor_params->is_null())
		return;  // however, if it is in the db but with 0 weight, that is OK

	Real Kphi = tor_params->K(0,0), phi0=tor_params->mu(0,0), phi_step=2*pi/tor_params->period();

	numeric::deriv::dihedral_p1_cosine_deriv(
		rsd.xyz( rt1 ), rsd.xyz( rt2 ), rsd.xyz( rt3 ), rsd.xyz( rt4 ), phi, f1, f2 );
	Real del_phi = basic::subtract_radian_angles(phi, phi0);
	del_phi = basic::periodic_range( del_phi, phi_step );
	if (linear_bonded_potential_ && std::fabs(del_phi)>1)
		dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * (del_phi>0? 1 : -1);
	else
		dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * del_phi;
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
	del_phi = basic::subtract_radian_angles(phi, phi0);
	del_phi = basic::periodic_range( del_phi, phi_step );
	if (linear_bonded_potential_ && std::fabs(del_phi)>1)
		dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * (del_phi>0? 1 : -1);
	else
		dE_dphi = (weights[ cart_bonded_torsion ] + weights[ cart_bonded ]) * Kphi * del_phi;
	r_atom_derivs[ rt3 ].f1() += dE_dphi * f1;
	r_atom_derivs[ rt3 ].f2() += dE_dphi * f2;

	f1 = f2 = Vector(0.0);
	numeric::deriv::dihedral_p1_cosine_deriv( 
		rsd.xyz( rt4 ), rsd.xyz( rt3 ), rsd.xyz( rt2 ), rsd.xyz( rt1 ), phi, f1, f2 );
	r_atom_derivs[ rt4 ].f1() += dE_dphi * f1;
	r_atom_derivs[ rt4 ].f2() += dE_dphi * f2;
}


/// @brief CartesianBondedEnergy does not have an atomic interation threshold
Distance
CartesianBondedEnergy::atomic_interaction_cutoff() const {
	return 0.0;
}

/// @brief CartesianBondedEnergy is context independent; indicates that no context graphs are required
void
CartesianBondedEnergy::indicate_required_context_graphs(utility::vector1< bool > & ) const {}

core::Size
CartesianBondedEnergy::version() const {
	return 1; // Initial versioning
}

} // namespace methods
} // namespace scoring
} // namespace core
