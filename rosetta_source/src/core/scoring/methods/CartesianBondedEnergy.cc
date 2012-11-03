// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CartesianBondedEnergy.cc
/// @brief  Harmonic constraints on rosetta's ideal parameters; structure of this code borrows heavily from MMBondAngle & MMTorsionAngle code
/// @author

// Unit headers
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

#include <core/scoring/PeptideBondedEnergyContainer.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>

// Utility headers
#include <utility/string_util.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/deriv/angle_deriv.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/dihedral_deriv.hh>

// options
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/optE.OptionKeys.gen.hh>

// C++ headers
#include <iostream>

#include <utility/vector1.hh>

#include <core/pose/PDBInfo.hh>
#define CUTOFF 10 

namespace boost {
namespace tuples {

// hash functions
std::size_t hash_value(residx_atm_quad const& e) {
	std::size_t seed = 0;
	boost::hash_combine(seed, e.get<0>());
	boost::hash_combine(seed, e.get<1>());
	boost::hash_combine(seed, e.get<2>());
	boost::hash_combine(seed, e.get<3>());
	boost::hash_combine(seed, e.get<4>());
	return seed;
}
std::size_t hash_value(residx_atm_triple const& e) {
	std::size_t seed = 0;
	boost::hash_combine(seed, e.get<0>());
	boost::hash_combine(seed, e.get<1>());
	boost::hash_combine(seed, e.get<2>());
	boost::hash_combine(seed, e.get<3>());
	return seed;
}
std::size_t hash_value(residx_atm_pair const& e) {
	std::size_t seed = 0;
	boost::hash_combine(seed, e.get<0>());
	boost::hash_combine(seed, e.get<1>());
	boost::hash_combine(seed, e.get<2>());
	return seed;
}

bool operator==(residx_atm_quad const& a,residx_atm_quad const& b) {
	return a.get<0>() == b.get<0>() && a.get<1>() == b.get<1>() &&
	       a.get<2>() == b.get<2>() && a.get<3>() == b.get<3>() &&
	       a.get<4>() == b.get<4>();
}
bool operator==(residx_atm_triple const& a,residx_atm_triple const& b) {
	return a.get<0>() == b.get<0>() && a.get<1>() == b.get<1>() &&
	       a.get<2>() == b.get<2>() && a.get<3>() == b.get<3>();
}
bool operator==(residx_atm_pair const& a,residx_atm_pair const& b) {
	return a.get<0>() == b.get<0>() && a.get<1>() == b.get<1>() &&
	       a.get<2>() == b.get<2>();
}

}
}



namespace core {
namespace scoring {
namespace methods {

/// @details Auto-generated virtual destructor
BondLengthDatabase::~BondLengthDatabase() {}

/// @details Auto-generated virtual destructor
BondAngleDatabase::~BondAngleDatabase() {}

/// @details Auto-generated virtual destructor
TorsionDatabase::~TorsionDatabase() {}


static basic::Tracer TR("core.scoring.CartesianBondedEnergy");

// default spring constants
static const Real K_LENGTH=300.0;
static const Real K_ANGLE=150.0;
static const Real K_TORSION=40.0;
static const Real K_TORSION_PROTON=40.0;  // proton chi

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

	// use 'annotated sequence' to id this restype
	char c = restype.name1();
	std::string restag(1,c);
	if (  ( !oneletter_code_specifies_aa(c) || name_from_aa( aa_from_oneletter_code(c) ) != restype.name() ) ) {
		restag = restag + '[' + restype.name() + ']';
	}
	// don't pool centroid and FA
	if (restype.residue_type_set().name() != "fa_standard")
		restag += "_"+restype.residue_type_set().name();
	return restag;
}


////////////////////////
// constructors
BondLengthDatabase::BondLengthDatabase(Real k_length_in) :
	mm_bondlength_library_( ScoringManager::get_instance()->get_MMBondLengthLibrary() ) {
	// check energy method options first, then command flag, then go to default
	if (k_length_in >= 0) {
		k_length = k_length_in;
	} else if (basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user()) {
		utility::vector1<core::Real> params = basic::options::option[ basic::options::OptionKeys::score::bonded_params ]();
		k_length = params[1];
	} else {
		k_length = K_LENGTH;
	}
	TR.Debug << "Set Klength to " << k_length << std::endl;
}

BondAngleDatabase::BondAngleDatabase(Real k_angle_in) :
	mm_bondangle_library_( ScoringManager::get_instance()->get_MMBondAngleLibrary() ) {
	// check energy method options first, then command flag, then go to default
	if (k_angle_in >= 0) {
		k_angle = k_angle_in;
	} else if (basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user()) {
		utility::vector1<core::Real> params = basic::options::option[ basic::options::OptionKeys::score::bonded_params ]();
		k_angle = params[2];
	} else {
		k_angle = K_ANGLE;
	}
	TR.Debug << "Set Kangle to " << k_angle << std::endl;
}

TorsionDatabase::TorsionDatabase(Real k_tors_in, Real k_prot_tors_in) {
	// check energy method options first, then command flag, then go to default
	k_torsion = K_TORSION;
	k_torsion_proton = K_TORSION_PROTON;
	if (basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user()) {
		utility::vector1<core::Real> params = basic::options::option[ basic::options::OptionKeys::score::bonded_params ]();
		k_torsion = params[3];
		k_torsion_proton = params[4];
	}
	if (k_tors_in >= 0)
		k_torsion = k_tors_in;
	if (k_prot_tors_in >= 0)
		k_torsion_proton = k_prot_tors_in;

	TR.Debug << "Set Ktors to " << k_torsion << std::endl;
	TR.Debug << "Set Ktorsprot to " << k_torsion_proton << std::endl;
}


//////////////////////
/// Torsion Database
// lookup ideal intra-res torsion; insert in DB if not there
void												
TorsionDatabase::lookup(
		core::chemical::ResidueType const & restype,
		int atm1, int atm2, int atm3, int atm4, Real &Kphi, Real &phi0, Real &phi_step ) {
	using namespace core::chemical;

	Kphi=k_torsion; // default
	phi_step=0;

	//fpd we should fix this to work with DNA/ligands rather than ignoring them
	if (!restype.is_protein()) {
		Kphi = phi0 = 0;
		return;
	}

	// use 'annotated sequence' to id this restype
	std::string restag = get_restag( restype );

	// lookup in table
	residx_atm_quad tuple( restag, atm1,atm2,atm3,atm4 );
	boost::unordered_map<residx_atm_quad,core::Real>::iterator b_it = torsions_.find( tuple );
	if ( b_it != torsions_.end() ) {
		phi0 = b_it->second;

		boost::unordered_map<residx_atm_quad,core::Real>::iterator k_it = Kphis_.find( tuple );
		if ( k_it != torsions_.end() )
			Kphi = k_it->second;

		boost::unordered_map<residx_atm_quad,core::Real>::iterator s_it = torsion_steps_.find( tuple );
		if ( s_it != torsions_.end() )
			phi_step = s_it->second;

	} else {
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
				if (hasCD && hasN)
				Kphi=0.0;
				phi0=0.0;
			}

			//fpd  cutpoint variants
			if ( restype.has_variant_type(chemical::CUTPOINT_UPPER) &&
					 ( ( restype.atom_name(atm2) == " N  " && restype.atom_name(atm3) == " CA ") ||
						 ( restype.atom_name(atm3) == " N  " && restype.atom_name(atm2) == " CA ")  ) ) {
				Kphi=0.0;
				phi0=0.0;
			}
		
			//fpd ignore centroid torsion
			if ( restype.atom_name(atm1) == " CEN" || restype.atom_name(atm4) == " CEN") {
				Kphi=0.0;
				phi0=0.0;
			}

			//fpd cyd
			if ( restype.aa() == aa_cys && restype.has_variant_type( chemical::DISULFIDE )  &&
					(restype.atom_name(atm2) == " SG " || restype.atom_name(atm3) == " SG ") ) {
				Kphi=0.0;
				phi0=0.0;
			}

			//fpd arg NH1-NH2 flip
			if ( restype.aa() == aa_arg  &&
					(   restype.atom_name(atm1) == " NH1" || restype.atom_name(atm4) == " NH1" 
					 || restype.atom_name(atm1) == " NH2" || restype.atom_name(atm4) == " NH2") ) {
				torsion_steps_[ tuple ] = phi_step = numeric::constants::f::pi;
			}

			// N terminal hydrogens
			if ( restype.has_variant_type(chemical::LOWER_TERMINUS) && 
					 ( ( restype.atom_name(atm2) == " N  " && restype.atom_name(atm3) == " CA ") ||
						 ( restype.atom_name(atm3) == " N  " && restype.atom_name(atm2) == " CA ")  ) ) {
				Kphi=k_torsion_proton;
				torsion_steps_[ tuple ] = phi_step = 2.0*numeric::constants::f::pi/3.0;
			}

			if (proton_chi) {
				Kphi=k_torsion_proton;
				//if (restype.aa()==core::chemical::aa_tyr) {
				//	torsion_steps_[ tuple ] = phi_step = numeric::constants::f::pi;
				//} else {
				//	torsion_steps_[ tuple ] = phi_step = 2.0*numeric::constants::f::pi/3.0;
				//}

				//fpd it turns out that the .params file of SER and THR put the 
				//      HG 60 degrees off from where it should be
				// query proton_chi to get the real value
				utility::vector1< Real > const & v = restype.proton_chi_samples( 1 );
				phi0 = numeric::constants::d::deg2rad*v[1];
				phi_step = std::fabs( numeric::constants::d::deg2rad*(v[2]-v[1]) );  // i guess this assumes we have two adjacent samples
				torsion_steps_[ tuple ] = phi_step;
			}
		}
		Kphis_[ tuple ] = Kphi;
		torsions_[ tuple ] = phi0;

		//TR << "torsion " << restype.name() << "  " 
		//   << restype.atom_name(atm1) << " " <<  restype.atom_name(atm2) << " "
		//   << restype.atom_name(atm3) << " " <<  restype.atom_name(atm4) << " --- " << phi0 << " , " << Kphi << std::endl;
	}
}


//////////////////////
/// BondAngle Database
// lookup ideal bondangle; insert in DB if not there
// atm ids > 0 ==> atom index
// atm ids < 0 ==> residue connection id
void
BondAngleDatabase::lookup( 
	pose::Pose const & pose, conformation::Residue const & res, int atm1, int atm2, int atm3, Real &Ktheta, Real &theta0 ) {
	using namespace core::chemical;

	Ktheta=k_angle;
	ResidueType const & restype = res.type();
	std::string restag = get_restag( restype );
 
	// lookup in table
	residx_atm_triple tuple( restag, atm1,atm2,atm3 );
	boost::unordered_map<residx_atm_triple,core::Real>::iterator b_it = bondangles_.find( tuple );

	if ( b_it != bondangles_.end() ) {
		theta0 = b_it->second;
		boost::unordered_map<residx_atm_triple,core::Real>::iterator k_it = Kangles_.find( tuple );
		if ( k_it != Kangles_.end() )
			Ktheta = k_it->second;
		return;
	}

	// Create mini-conformation for idealized residue
	conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );
	
	// if interres angle, get neighbor info 
	Size neighb_mmatm;
	int const resconn = std::min( std::min( atm1, atm2 ), atm3 );	
	if ( resconn < 0) {		
		Size const neighb = res.residue_connection_partner( -resconn );
		conformation::Residue const neighb_res = pose.residue( neighb );
		Size const neighb_resconn = res.residue_connection_conn_id( -resconn );
		ResidueType const & neighb_restype = neighb_res.type();
		Size const neighb_atm = neighb_restype.residue_connection( neighb_resconn ).atomno();
		neighb_mmatm = neighb_restype.atom( neighb_atm ).mm_atom_type_index();
	}

	// use these for mm Ktheta lookup
	Size mmatm1, mmatm2, mmatm3;						
	
	// get angle
	numeric::xyzVector<core::Real> x,y,z; // atom coords
	if (atm1 > 0) {
		mmatm1 = restype.atom( atm1 ).mm_atom_type_index();
		x = newres->atom( atm1 ).xyz();
	}
	else {
		mmatm1 = neighb_mmatm;
		x = newres->residue_connection( -atm1 ).icoor().build( restype );
	}
	if (atm2 > 0) {
		mmatm2 = restype.atom( atm2 ).mm_atom_type_index();
		y = newres->atom( atm2 ).xyz();
	}
	else {
		mmatm2 = neighb_mmatm;
		y = newres->residue_connection( -atm2 ).icoor().build( restype );
	}
	if (atm3 > 0) {
		mmatm3 = restype.atom( atm3 ).mm_atom_type_index();
		z = newres->atom( atm3 ).xyz();
	}
	else {
		mmatm3 = neighb_mmatm;
		z = newres->residue_connection( -atm3 ).icoor().build( restype );
	}

	theta0 = numeric::angle_radians ( x,y,z );
	
	// get Ktheta from CHARMM
	if ( basic::options::option[ basic::options::OptionKeys::score::bonded_params_CHARMM_angle ].user() ) {

		mm::mm_bondangle_atom_tri mm_atomtype_set( mmatm1, mmatm2, mmatm3);
		mm::mm_bondangle_library_citer_pair pair = mm_bondangle_library_.lookup(
			mm_atomtype_set.key1(),
			mm_atomtype_set.key2(),
			mm_atomtype_set.key3() 
		);
		Ktheta = ((pair.first)->second).key1();			// look in core/mm/MMBondAngleLibrary.hh to figure out this one
		//TR << "Set Ktheta to " << Ktheta << std::endl;
	}
		
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

	// fpd ignore centroid angle in ALA and GLY
	if ( (restype.aa() == core::chemical::aa_ala || restype.aa() == core::chemical::aa_gly) &&
	     ( (atm1>0 && restype.atom_name(atm1) == " CEN") || (atm3>0 && restype.atom_name(atm3) == " CEN") ) ) {
		Ktheta = theta0 = 0.0;
	}

	//fpd  cutpoint variants (?????)
	if ( restype.has_variant_type(chemical::CUTPOINT_UPPER) &&
			 ( (atm1>0 && restype.atom_name(atm1) == "OVU1") || (atm3>0 && restype.atom_name(atm3) == "OVU1") ) ) {
		Ktheta = theta0 = 0.0;
	}
	if ( restype.has_variant_type(chemical::CUTPOINT_LOWER) &&
			 ( (atm1>0 && restype.atom_name(atm1) == "OVL1") || (atm3>0 && restype.atom_name(atm3) == "OVL1") ) ) {
		Ktheta = theta0 = 0.0;
	}

	bondangles_[ tuple ] = theta0;
	Kangles_[ tuple ] = Ktheta;

	return;
}

//////////////////////
/// BondLength Database
// lookup ideal bondangle; insert in DB if not there
// atm ids > 0 ==> atom index
// atm ids < 0 ==> residue connection id
void
BondLengthDatabase::lookup
	( pose::Pose const & pose, conformation::Residue const & res, int atm1, int atm2, Real &Kd, Real &d0 ) {
	using namespace core::chemical;

	Kd=k_length;
	ResidueType const & restype = res.type();
	std::string restag = get_restag( restype );

	// lookup in table
	residx_atm_pair tuple( restag, atm1,atm2 );
	boost::unordered_map<residx_atm_pair,core::Real>::iterator b_it = bondlengths_.find( tuple );

	if ( b_it != bondlengths_.end() ) {
		d0 = b_it->second;
		boost::unordered_map<residx_atm_pair,core::Real>::iterator k_it = Klengths_.find( tuple );
		if ( k_it != Klengths_.end() )
			Kd = k_it->second;
		return;
	}

	// Create mini-conformation for idealized residue
	core::conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );

	// if interres angle, get neighbor info 
	Size neighb_mmatm;
	int const resconn = std::min( atm1, atm2 );	
	if ( resconn < 0) {		
		Size const neighb = res.residue_connection_partner( -resconn );
		conformation::Residue const neighb_res = pose.residue( neighb );
		Size const neighb_resconn = res.residue_connection_conn_id( -resconn );
		ResidueType const & neighb_restype = neighb_res.type();
		Size const neighb_atm = neighb_restype.residue_connection( neighb_resconn ).atomno();
		neighb_mmatm = neighb_restype.atom( neighb_atm ).mm_atom_type_index();
		//TR.Debug << "resconn atom name: " << neighb_res.atom_name( neighb_atm ) << "rsd atom name: " << res.atom_name( std::max(atm1, atm2)) << std::endl;
		//TR.Debug << "rsd atom no: " << std::max(atm1, atm2) << " resconn atom no: " << neighb_atm << std::endl;
		//TR.Debug << "rsd mm_atom no: " << restype.mm_atom_type_index( atm1 ) << " resconn mm_atom no: " << neighb_mmatm << std::endl;
	}

	// use these for mm Ktheta lookup
	Size mmatm1, mmatm2;

	// get length
	numeric::xyzVector<core::Real> x,y,z; // atom coords
	if (atm1 > 0) {
		mmatm1 = restype.atom( atm1 ).mm_atom_type_index();
		x = newres->atom( atm1 ).xyz();
	}
	else {
		mmatm1 = neighb_mmatm;
		x = newres->residue_connection( -atm1 ).icoor().build( restype );
	}
	if (atm2 > 0) {
		mmatm2 = restype.atom( atm2 ).mm_atom_type_index();
		y = newres->atom( atm2 ).xyz();
	}
	else {
		mmatm2 = neighb_mmatm;
		y = newres->residue_connection( -atm2 ).icoor().build( restype );
	}

	d0 = (x-y).length();
	  
	// get Kd from CHARMM
	if ( basic::options::option[ basic::options::OptionKeys::score::bonded_params_CHARMM_length ].user() ) {
		mm::mm_bondlength_atom_pair mm_atomtype_set( mmatm1, mmatm2 );
		mm::mm_bondlength_library_citer_pair pair = mm_bondlength_library_.lookup(
			mm_atomtype_set.key1(),
			mm_atomtype_set.key2()
		);
		Kd = ((pair.first)->second).key1();			// look in core/mm/MMBondLengthLibrary.hh to figure out this one

		/*// lookup is supposed to return 0 for virtual and exit with error if param not found. fatal error is disabled for now.
		if ( Kd == 0 and !( restype.atom_type(atm1).is_virtual() || restype.atom_type(atm2).is_virtual() ) ) {
			Kd = 300;								//arbitrary number
			TR << "Unknown Kd set to " << Kd << std::endl;
		}
		else
			TR.Debug << "Set Kd to " << Kd << std::endl;*/
	}

	// fpd ignore proline C-ND which is handled by pro_close
	if (restype.aa() == core::chemical::aa_pro &&
	    (atm1>0 && atm2>0) &&
	    ( ( restype.atom_name(atm1) == " CD " && restype.atom_name(atm2) == " N  ") ||
	      ( restype.atom_name(atm2) == " CD " && restype.atom_name(atm1) == " N  ") ) ) {
		Kd = d0 = 0.0;
	}

	bondlengths_[ tuple ] = d0;
	Klengths_[ tuple ] = Kd;

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
	options.get_cartesian_bonded_parameters( cartbonded_len_, cartbonded_ang_, cartbonded_tors_, cartbonded_proton_ );
	db_angle_ = new BondAngleDatabase(cartbonded_len_);
	db_length_ = new BondLengthDatabase(cartbonded_ang_);
	db_torsion_ = new TorsionDatabase(cartbonded_tors_, cartbonded_proton_);
}

CartesianBondedEnergy::CartesianBondedEnergy( CartesianBondedEnergy const & src ) : parent( src ) {
	linear_bonded_potential_ = src.linear_bonded_potential_;
	db_angle_ = src.db_angle_;
	db_length_ = src.db_length_;
	db_torsion_ = src.db_torsion_;
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
		if ( dec->size() != pose.total_residue() ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		TR << "Creating new peptide-bonded energy container (" << pose.total_residue() << ")" << std::endl;
		utility::vector1< ScoreType > s_types;
		s_types.push_back( cart_bonded );
		s_types.push_back( cart_bonded_angle );
		s_types.push_back( cart_bonded_length );
		s_types.push_back( cart_bonded_torsion );
		LREnergyContainerOP new_dec = new PeptideBondedEnergyContainer( pose.total_residue(), s_types );
		energies.set_long_range_container( lr_type, new_dec );
	}
}

///
bool
CartesianBondedEnergy::defines_intrares_energy( EnergyMap const & ) const {
	return true;
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

	// backbone C-N-CA-H
	conformation::Residue const &res_low = (rsd1.seqpos() < rsd2.seqpos())? rsd1 : rsd2;
	conformation::Residue const &res_high = (rsd1.seqpos() < rsd2.seqpos())? rsd2 : rsd1;

	core::Real energy_torsion = 0.0;

	if (!res_low.is_protein() || !res_high.is_protein()) return;

	if (res_high.aa() != aa_pro && res_high.has(" H  ") ) {  // check for terminal variants
		core::Size atm1 = res_high.atom_index(" CA ");
		core::Size atm2 = res_low.atom_index(" C  ");
		core::Size atm3 = res_high.atom_index(" N  ");
		core::Size atm4 = res_high.atom_index(" H  ");

		Real Kphi=db_torsion_->ktors(), phi0=numeric::constants::f::pi, phi_step=numeric::constants::f::pi;

		Real angle = numeric::dihedral_radians
			( res_high.atom( atm1 ).xyz(), res_low.atom( atm2 ).xyz(), res_high.atom( atm3 ).xyz(), res_high.atom( atm4 ).xyz() );
		Real del_phi = basic::subtract_radian_angles(angle, phi0);
		del_phi = basic::periodic_range( del_phi, phi_step );

		if (pose.pdb_info() && 0.5*Kphi*del_phi*del_phi > CUTOFF) {
			TR.Debug << pose.pdb_info()->name() << " seqpos: " << res_low.seqpos() << " pdbpos: " << pose.pdb_info()->number(res_low.seqpos()) << " improper torsion: " <<
				res_low.name() << " : " <<
				res_high.atom_name( atm1 ) << " , " << res_low.atom_name( atm2 ) << " , " <<
				res_high.atom_name( atm3 ) << " , " << res_high.atom_name( atm4 ) << "   " <<
				Kphi << " " << 0.5*Kphi*del_phi*del_phi << std::endl;
		}

		if (linear_bonded_potential_ && std::fabs(del_phi)>1) 
			energy_torsion += 0.5*Kphi*std::fabs(del_phi);
		else 
			energy_torsion += 0.5*Kphi*del_phi*del_phi;
	}

	{
		core::Size atm1 = res_low.atom_index(" CA ");
		core::Size atm2 = res_high.atom_index(" N  ");
		core::Size atm3 = res_low.atom_index(" C  ");
		core::Size atm4 = res_low.atom_index(" O  ");

		Real Kphi=db_torsion_->ktors(), phi0=numeric::constants::f::pi, phi_step=numeric::constants::f::pi;

		Real angle = numeric::dihedral_radians
			( res_low.atom( atm1 ).xyz(), res_high.atom( atm2 ).xyz(), res_low.atom( atm3 ).xyz(), res_low.atom( atm4 ).xyz() );
		Real del_phi = basic::subtract_radian_angles(angle, phi0);
		del_phi = basic::periodic_range( del_phi, phi_step );

		if (pose.pdb_info() && 0.5*Kphi*del_phi*del_phi > CUTOFF) {
			TR.Debug << pose.pdb_info()->name() << " seqpos: " << res_low.seqpos() << " pdbpos: " << pose.pdb_info()->number(res_low.seqpos()) << " improper torsion: " <<
				res_low.name() << " : " <<
				res_low.atom_name( atm1 ) << " , " << res_high.atom_name( atm2 ) << " , " <<
				res_low.atom_name( atm3 ) << " , " << res_low.atom_name( atm4 ) << "   " <<
				Kphi << " " << 0.5*Kphi*del_phi*del_phi << std::endl;
		}

		if (linear_bonded_potential_ && std::fabs(del_phi)>1) 
			energy_torsion += 0.5*Kphi*std::fabs(del_phi);
		else 
			energy_torsion += 0.5*Kphi*del_phi*del_phi;
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

	if (!pose.is_fullatom()) return;

	core::Size atm1, atm2, atm3, atm4;
	// asp CB-CG-OD1-OD2
	if (rsd.aa() == aa_asp) {
		atm1 = rsd.atom_index(" CB ");
		atm2 = rsd.atom_index(" CG ");
		atm3 = rsd.atom_index(" OD1");
		atm4 = rsd.atom_index(" OD2");
	} else if (rsd.aa() == aa_glu) {
		// glu CG-CD-OE1-OE2
		atm1 = rsd.atom_index(" CG ");
		atm2 = rsd.atom_index(" CD ");
		atm3 = rsd.atom_index(" OE1");
		atm4 = rsd.atom_index(" OE2");
	} else if (rsd.aa() == aa_asn) {
		// asn CB-CG-OD-ND
		atm1 = rsd.atom_index(" CB ");
		atm2 = rsd.atom_index(" CG ");
		atm3 = rsd.atom_index(" OD1");
		atm4 = rsd.atom_index(" ND2");
	} else if (rsd.aa() == aa_gln) {
		// gln CG-CD-OE-NE
		atm1 = rsd.atom_index(" CG ");
		atm2 = rsd.atom_index(" CD ");
		atm3 = rsd.atom_index(" OE1");
		atm4 = rsd.atom_index(" NE2");
	} else {
		return;
	}

	Real Kphi=db_torsion_->ktors(), phi0=numeric::constants::f::pi, phi_step=numeric::constants::f::pi;

	Real angle = numeric::dihedral_radians
		( rsd.atom( atm1 ).xyz(), rsd.atom( atm2 ).xyz(), rsd.atom( atm3 ).xyz(), rsd.atom( atm4 ).xyz() );
	Real del_phi = basic::subtract_radian_angles(angle, phi0);
	del_phi = basic::periodic_range( del_phi, phi_step );

	if (pose.pdb_info() && 0.5*Kphi*del_phi*del_phi > CUTOFF) {
		TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " << pose.pdb_info()->number(rsd.seqpos()) << " improper torsion: " <<
			rsd.name() << " : " <<
			rsd.atom_name( atm1 ) << " , " << rsd.atom_name( atm2 ) << " , " <<
			rsd.atom_name( atm3 ) << " , " << rsd.atom_name( atm4 ) << "   " <<
			Kphi << " " << 0.5*Kphi*del_phi*del_phi << std::endl;
	}

	core::Real energy_torsion;
	if (linear_bonded_potential_ && std::fabs(del_phi)>1) 
		energy_torsion = 0.5*Kphi*std::fabs(del_phi);
	else 
		energy_torsion = 0.5*Kphi*del_phi*del_phi;

	emap[ cart_bonded_torsion ] += energy_torsion;
	emap[ cart_bonded ] += energy_torsion; // potential double counting
}


// deriv impropers
void
CartesianBondedEnergy::eval_improper_torsions_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const & ,
	ScoreFunction const & ,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	using namespace core::chemical;

	core::chemical::ResidueType const & restype( pose.residue_type( id.rsd() ) );
	core::conformation::Residue const & res( pose.residue( id.rsd() ) );
	Size const atomno( id.atomno());
	std::string atm_name = restype.atom_name(atomno);
	if ( !restype.is_protein() ) return;

	core::Size atm1, atm2, atm3, atm4;
	numeric::xyzVector< core::Real > atm1xyz, atm2xyz, atm3xyz, atm4xyz;

	Vector f1(0.0), f2(0.0);
	Real phi(0.0);
	Real dE_dphi;

	// INTER-RES
	{
		if ( atm_name == " C  " || atm_name == " N  " || atm_name == " CA " || atm_name == " H  ") {
			bool applyThis = true;
			if (atm_name == " C  " && (id.rsd() == pose.total_residue() || !pose.residue( id.rsd()+1 ).is_protein()) )
				applyThis=false;
			if (atm_name != " C  " && (id.rsd() == 1 || !pose.residue( id.rsd()-1 ).is_protein()) )
				applyThis=false;

			if (applyThis) {
				core::conformation::Residue const & res_low = (atm_name == " C  ")? res : pose.residue( id.rsd()-1 );
				core::conformation::Residue const & res_high = (atm_name == " C  ")? pose.residue( id.rsd()+1 ) : res;
				if (res_high.aa() != aa_pro && res_low.is_bonded(res_high) && res_high.has(" H  ")) {
					atm1 = res_high.atom_index(" CA ");
					atm2 = res_low.atom_index(" C  "); 
					atm3 = res_high.atom_index(" N  ");
					atm4 = res_high.atom_index(" H  ");
					atm1xyz = res_high.xyz( atm1 ); atm2xyz = res_low.xyz( atm2 ); atm3xyz = res_high.xyz( atm3 ); atm4xyz = res_high.xyz( atm4 );

					if ( atm1 == atomno ) {
						numeric::deriv::dihedral_p1_cosine_deriv( atm1xyz, atm2xyz, atm3xyz, atm4xyz, phi, f1, f2 );
					} else if ( atm2 == atomno ) {
						numeric::deriv::dihedral_p2_cosine_deriv( atm1xyz, atm2xyz, atm3xyz, atm4xyz, phi, f1, f2 );
					} else if ( atm3 == atomno ) {
						numeric::deriv::dihedral_p2_cosine_deriv( atm4xyz, atm3xyz, atm2xyz, atm1xyz, phi, f1, f2 );
					} else if ( atm4 == atomno ) {
						numeric::deriv::dihedral_p1_cosine_deriv( atm4xyz, atm3xyz, atm2xyz, atm1xyz, phi, f1, f2 );
					}

					Real Kphi=db_torsion_->ktors(), phi0=numeric::constants::f::pi, phi_step=numeric::constants::f::pi;
					Real del_phi = basic::subtract_radian_angles(phi, phi0);
					del_phi = basic::periodic_range( del_phi, phi_step );
					if (linear_bonded_potential_ && std::fabs(del_phi)>1) {
						dE_dphi = weights[ cart_bonded_torsion ] * Kphi * (del_phi>0? 1 : -1);
						dE_dphi += weights[ cart_bonded ] * Kphi * (del_phi>0? 1 : -1);
					} else {
						dE_dphi = weights[ cart_bonded_torsion ] * Kphi * del_phi;
						dE_dphi += weights[ cart_bonded ] * Kphi * del_phi;
					}
					F1 += dE_dphi * f1;
					F2 += dE_dphi * f2;
				} 
			}
		}

		if ( atm_name == " C  " || atm_name == " N  " || atm_name == " CA " || atm_name == " O  ") {
			bool applyThis = true;
			if (atm_name != " N  " && (id.rsd() == pose.total_residue() || !pose.residue( id.rsd()+1 ).is_protein()) )
				applyThis=false;
			if (atm_name == " N  " && (id.rsd() == 1 || !pose.residue( id.rsd()-1 ).is_protein()) )
				applyThis=false;

			if (applyThis) {
				core::conformation::Residue const & res_low = (atm_name != " N  ")? res : pose.residue( id.rsd()-1 );
				core::conformation::Residue const & res_high = (atm_name != " N  ")? pose.residue( id.rsd()+1 ) : res;
				if (res_low.is_bonded(res_high)) {
					atm1 = res_low.atom_index(" CA ");
					atm2 = res_high.atom_index(" N  ");
					atm3 = res_low.atom_index(" C  ");
					atm4 = res_low.atom_index(" O  ");
					atm1xyz = res_low.xyz( atm1 ); atm2xyz = res_high.xyz( atm2 ); atm3xyz = res_low.xyz( atm3 ); atm4xyz = res_low.xyz( atm4 );

					if ( atm1 == atomno ) {
						numeric::deriv::dihedral_p1_cosine_deriv( atm1xyz, atm2xyz, atm3xyz, atm4xyz, phi, f1, f2 );
					} else if ( atm2 == atomno ) {
						numeric::deriv::dihedral_p2_cosine_deriv( atm1xyz, atm2xyz, atm3xyz, atm4xyz, phi, f1, f2 );
					} else if ( atm3 == atomno ) {
						numeric::deriv::dihedral_p2_cosine_deriv( atm4xyz, atm3xyz, atm2xyz, atm1xyz, phi, f1, f2 );
					} else if ( atm4 == atomno ) {
						numeric::deriv::dihedral_p1_cosine_deriv( atm4xyz, atm3xyz, atm2xyz, atm1xyz, phi, f1, f2 );
					}

					Real Kphi=db_torsion_->ktors(), phi0=numeric::constants::f::pi, phi_step=numeric::constants::f::pi;
					Real del_phi = basic::subtract_radian_angles(phi, phi0);
					del_phi = basic::periodic_range( del_phi, phi_step );
					if (linear_bonded_potential_ && std::fabs(del_phi)>1) {
						dE_dphi = weights[ cart_bonded_torsion ] * Kphi * (del_phi>0? 1 : -1);
						dE_dphi += weights[ cart_bonded ] * Kphi * (del_phi>0? 1 : -1);
					} else {
						dE_dphi = weights[ cart_bonded_torsion ] * Kphi * del_phi;
						dE_dphi += weights[ cart_bonded ] * Kphi * del_phi;
					}
					F1 += dE_dphi * f1;
					F2 += dE_dphi * f2;
				}
			} 
		}
	}

	// INTRA-RES
	{
		bool applyThis = true;
		if (!pose.is_fullatom()) return;

		if (res.aa() == aa_asp) {
			// asp CB-CG-OD1-OD2
			atm1 = res.atom_index(" CB ");
			atm2 = res.atom_index(" CG ");
			atm3 = res.atom_index(" OD1");
			atm4 = res.atom_index(" OD2");
		} else if (res.aa() == aa_glu) {
			// glu CG-CD-OE1-OE2
			atm1 = res.atom_index(" CG ");
			atm2 = res.atom_index(" CD ");
			atm3 = res.atom_index(" OE1");
			atm4 = res.atom_index(" OE2");
		} else if (res.aa() == aa_asn) {
			// asn CB-CG-OD-ND
			atm1 = res.atom_index(" CB ");
			atm2 = res.atom_index(" CG ");
			atm3 = res.atom_index(" OD1");
			atm4 = res.atom_index(" ND2");
		} else if (res.aa() == aa_gln) {
			// gln CG-CD-OE-NE
			atm1 = res.atom_index(" CG ");
			atm2 = res.atom_index(" CD ");
			atm3 = res.atom_index(" OE1");
			atm4 = res.atom_index(" NE2");
		} else {
			applyThis = false;
		}

		if (applyThis) {
			atm1xyz = res.xyz( atm1 );
			atm2xyz = res.xyz( atm2 );
			atm3xyz = res.xyz( atm3 );
			atm4xyz = res.xyz( atm4 );
		
			if ( atm1 == atomno ) {
				numeric::deriv::dihedral_p1_cosine_deriv( atm1xyz, atm2xyz, atm3xyz, atm4xyz, phi, f1, f2 );
			} else if ( atm2 == atomno ) {
				numeric::deriv::dihedral_p2_cosine_deriv( atm1xyz, atm2xyz, atm3xyz, atm4xyz, phi, f1, f2 );
			} else if ( atm3 == atomno ) {
				numeric::deriv::dihedral_p2_cosine_deriv( atm4xyz, atm3xyz, atm2xyz, atm1xyz, phi, f1, f2 );
			} else if ( atm4 == atomno ) {
				numeric::deriv::dihedral_p1_cosine_deriv( atm4xyz, atm3xyz, atm2xyz, atm1xyz, phi, f1, f2 );
			} else {
				applyThis = false;
			}
		}

		if (applyThis) {
			Real Kphi=db_torsion_->ktors(), phi0=numeric::constants::f::pi, phi_step=numeric::constants::f::pi;
			Real del_phi = basic::subtract_radian_angles(phi, phi0);
			del_phi = basic::periodic_range( del_phi, phi_step );
			if (linear_bonded_potential_ && std::fabs(del_phi)>1) {
				dE_dphi = weights[ cart_bonded_torsion ] * Kphi * (del_phi>0? 1 : -1);
				dE_dphi += weights[ cart_bonded ] * Kphi * (del_phi>0? 1 : -1);
			} else {
				dE_dphi = weights[ cart_bonded_torsion ] * Kphi * del_phi;
				dE_dphi += weights[ cart_bonded ] * Kphi * del_phi;
			}

			F1 += dE_dphi * f1;
			F2 += dE_dphi * f2;
		}
	}
}


void
CartesianBondedEnergy::residue_pair_energy(
   conformation::Residue const & rsd1,
	 conformation::Residue const & rsd2,
	 pose::Pose const & pose,
	 ScoreFunction const & sf,
	 EnergyMap & emap
) const {
	// bail out if the residues aren't bonded
	if (!rsd1.is_bonded(rsd2)) return;

	//fpd is we are not connected in the fold tree then don't constrain across the adjacent residues
	if ( pose.fold_tree().is_cutpoint( std::min( rsd1.seqpos(), rsd2.seqpos() ) ) ) return;

	Real energy_angle = 0;
	Real energy_length = 0;

	// get residue types
	chemical::ResidueType const & rsd1_type = rsd1.type();
	chemical::ResidueType const & rsd2_type = rsd2.type();

	utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );

	for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {
		Size const resconn_id1( r1_resconn_ids[ii] );
		Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

		// only look at rsd1->rsd2 connections
		//if (rsd1.connected_residue_at_resconn( resconn_id1 ) != rsd2.seqpos() ) continue;

		Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
		Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

		/////////////
		/// compute the bond-angle energies from pairs of atoms within-1 bond on rsd1 with
		/// the the connection atom on rsd2.
		utility::vector1< chemical::two_atom_set > const & rsd1_atoms_wi1_bond_of_ii(
			rsd1_type.atoms_within_one_bond_of_a_residue_connection( resconn_id1 ));
		for ( Size jj = 1; jj <= rsd1_atoms_wi1_bond_of_ii.size(); ++jj ) {
			assert( rsd1_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno1 );
			Size const res1_lower_atomno = rsd1_atoms_wi1_bond_of_ii[ jj ].key2();

			Real const angle = numeric::angle_radians(
				rsd1.atom( res1_lower_atomno ).xyz(),
				rsd1.atom( resconn_atomno1 ).xyz(),
				rsd2.atom( resconn_atomno2 ).xyz() );

			// lookup Ktheta and theta0
			Real Ktheta, theta0;
			db_angle_->lookup( pose, rsd1, res1_lower_atomno, resconn_atomno1, -resconn_id1, Ktheta, theta0 );

			if (Ktheta == 0.0) continue;

			if (pose.pdb_info() && 0.5*Ktheta*(angle-theta0)*(angle-theta0) > CUTOFF) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd1.seqpos() << " pdbpos: " << pose.pdb_info()->number(rsd1.seqpos()) << " angle rsd1: " << rsd1.name() << ":" << rsd1_type.atom_name( res1_lower_atomno ) << " , "
				   << rsd1_type.atom_name( resconn_atomno1 ) << " , " << rsd2_type.atom_name( resconn_atomno2 )
				   << "    " << angle << "  " << theta0 << "   "
				   << Ktheta << "    " << 0.5*Ktheta*(angle-theta0)*(angle-theta0) << std::endl;
			}

			// accumulate the energy
			if (linear_bonded_potential_ && std::fabs(angle-theta0)>1) {
				energy_angle += 0.5*Ktheta*std::fabs(angle-theta0);
			} else {
				energy_angle += 0.5*Ktheta*(angle-theta0) * (angle-theta0);
			}
		}

		/////////////
		/// compute the bond-angle energies from pairs of atoms within-1 bond on rsd2 with
		/// the the connection atom on rsd1.
		utility::vector1< chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
			rsd2_type.atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));
		for ( Size jj = 1; jj <= rsd2_atoms_wi1_bond_of_ii.size(); ++jj ) {
			assert( rsd2_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno2 );
			Size const res2_lower_atomno = rsd2_atoms_wi1_bond_of_ii[ jj ].key2();

			// lookup Ktheta and theta0
			Real Ktheta, theta0;
			db_angle_->lookup( pose, rsd2, res2_lower_atomno, resconn_atomno2, -resconn_id2, Ktheta, theta0 );

			if (Ktheta == 0.0) continue;
			Real const angle = numeric::angle_radians(
				rsd2.atom( res2_lower_atomno ).xyz(),
				rsd2.atom( resconn_atomno2 ).xyz(),
				rsd1.atom( resconn_atomno1 ).xyz() );

			if (pose.pdb_info() && 0.5*Ktheta*(angle-theta0)*(angle-theta0) > CUTOFF) {
				TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd2.seqpos() << " pdbpos: " << pose.pdb_info()->number(rsd2.seqpos()) << " angle rsd2: " << rsd2.name() << ":" << rsd2_type.atom_name( res2_lower_atomno ) << " , "
			   << rsd2_type.atom_name( resconn_atomno2 ) << " , " << rsd1_type.atom_name( resconn_atomno1 )
				 << "    " << angle << "  " << theta0 << "   "
				 << Ktheta << "    " << 0.5*Ktheta*(angle-theta0)*(angle-theta0) << std::endl;
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

		//TR.Debug << "interres resconn atom name: " << rsd1.atom_name( resconn_atomno1 ) << "rsd atom name: " << rsd2.atom_name( resconn_atomno2 ) << std::endl;
		//TR.Debug << "interres rsd1 atom no: " << resconn_atomno1 << " rsd2 (resconn) atom no: " << resconn_atomno2 << std::endl;
		// lookup Ktheta and theta0
		Real Kd, d0;
		db_length_->lookup( pose, rsd1, resconn_atomno1, -resconn_id1, Kd, d0 );

		if (pose.pdb_info() && 0.5*Kd*(length-d0) * (length-d0) > CUTOFF) {
			TR.Debug << pose.pdb_info()->name() << " pdbpos rsd1: " << pose.pdb_info()->number(rsd1.seqpos()) << " length rsd1 rsd2: " << rsd1.seqpos() << " -- " << rsd2.seqpos() << "  "
			   << rsd1.name() << ":" << rsd1_type.atom_name( resconn_atomno1 ) << " , " << rsd2_type.atom_name( resconn_atomno2 )
			   << "    " << length << " [" << d0 << "]" 
			   << Kd << "    " << 0.5*Kd*std::fabs(length-d0) << std::endl;
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

	eval_improper_torsions( rsd1, rsd2, pose, sf, emap );
}

void
CartesianBondedEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sf,
	EnergyMap & emap
) const
{
	Real energy_angle = 0;
	Real energy_length = 0;
	Real energy_torsion = 0;	

	// get residue type
	chemical::ResidueType const & rsd_type = rsd.type();

	// for each torsion _that doesn't correspond to a DOF_ID in the pose_
	for ( Size dihe = 1; dihe <= rsd_type.ndihe(); ++dihe ){
		// get ResidueType ints
		int rt1 = ( rsd_type.dihedral( dihe ) ).key1();
		int rt2 = ( rsd_type.dihedral( dihe ) ).key2();
		int rt3 = ( rsd_type.dihedral( dihe ) ).key3();
		int rt4 = ( rsd_type.dihedral( dihe ) ).key4();

		// lookup Ktheta and theta0
		Real Kphi, phi0, phi_step;
		db_torsion_->lookup( rsd.type(), rt1, rt2, rt3, rt4, Kphi, phi0, phi_step );
		if (Kphi == 0.0) continue;

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
 			  rsd.atom_name( rt3 ) << " , " << rsd.atom_name( rt4 ) << "   " <<
			  Kphi << " " << 0.5*Kphi*del_phi*del_phi << std::endl;
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

		// check for vrt
		//if ( rsd_type.atom_type(rt1).is_virtual()
		//       || rsd_type.atom_type(rt2).is_virtual()
		//       || rsd_type.atom_type(rt3).is_virtual() )
		if ( rsd_type.aa() == core::chemical::aa_vrt)
			continue;

		// lookup Ktheta and theta0
		Real Ktheta, theta0;
		db_angle_->lookup( pose, rsd, rt1, rt2, rt3, Ktheta, theta0 );
		if (Ktheta == 0.0) continue;

		// get angle
		Real const angle = numeric::angle_radians(
			rsd.atom( rt1 ).xyz(),
			rsd.atom( rt2 ).xyz(),
			rsd.atom( rt3 ).xyz() );

 		if ( pose.pdb_info() && 0.5*Ktheta*(angle-theta0) * (angle-theta0) > CUTOFF) {
 			TR.Debug << pose.pdb_info()->name() << " seqpos: " << rsd.seqpos() << " pdbpos: " << pose.pdb_info()->number(rsd.seqpos()) << " intrares angle: " <<
 				rsd_type.name() << " : " <<
 			  rsd.atom_name( rt1 ) << " , " << rsd.atom_name( rt2 ) << " , " <<
 			  rsd.atom_name( rt3 ) << "   " << angle << "  " << theta0 << "     " <<
 			  Ktheta << " " << 0.5*Ktheta*(angle-theta0) * (angle-theta0) << std::endl;
 		}
		// accumulate the energy
		if (linear_bonded_potential_ && std::fabs(angle - theta0)>1) {
			energy_angle += 0.5*Ktheta*std::fabs(angle-theta0);
		} else {
			energy_angle += 0.5*Ktheta*(angle-theta0) * (angle-theta0);
		}
	}

	// for each bond in the residue
	// for each bonded atom
	for (Size atm_i=1; atm_i<=rsd_type.natoms(); ++atm_i) {
		chemical::AtomIndices atm_nbrs = rsd_type.nbrs( atm_i );
		for (Size j=1; j<=atm_nbrs.size(); ++j) {
			Size atm_j = atm_nbrs[j];
			if ( atm_i<atm_j ) { // only score each bond once -- use restype index to define ordering
				// check for vrt
				//if ( rsd_type.atom_type(atm_i).is_virtual() || rsd_type.atom_type(atm_j).is_virtual() )
				if ( rsd_type.aa() == core::chemical::aa_vrt)
					continue;

				//TR << "intrares resconn atom name: " << rsd.atom_name( atm_i ) << "rsd atom name: " << rsd.atom_name( atm_j ) << std::endl;
				//TR << "intrares rsd atom no: " << atm_i << " rsd atom no: " << atm_j << std::endl;
				// lookup Ktheta and theta0
				Real Kd, d0;
				db_length_->lookup( pose, rsd, atm_i, atm_j, Kd, d0 );
				if (Kd == 0.0) continue;

				Real const d = ( rsd.atom( atm_i ).xyz()-rsd.atom( atm_j ).xyz() ).length();

				if ( 0.5*Kd*(d-d0)*(d-d0) > CUTOFF) {
					TR.Debug << "length-ideal" << 
						rsd_type.name() << " : " <<
						rsd.atom_name( atm_i ) << " , " << rsd.atom_name( atm_j ) << "   " <<
						Kd << " " << 0.5*Kd*(d-d0)*(d-d0) << std::endl;
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

	eval_improper_torsions( rsd, pose, sf, emap );
}


///////////////////////////
///
void
CartesianBondedEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const & sf,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	Vector LF1( 0.0 );
	Vector LF2( 0.0 );

	core::chemical::ResidueType const & restype( pose.residue_type( id.rsd() ) );
	core::conformation::Residue const & res( pose.residue( id.rsd() ) );
	Size const atomno( id.atomno());

	if ( restype.aa() == core::chemical::aa_vrt) return;

	// (A0) intra-res torsions
	utility::vector1< Size > const & diheds( restype.dihedrals_for_atom( id.atomno() ) );
	for ( Size ii = 1, ii_end = diheds.size(); ii <= ii_end; ++ii ) {
		chemical::dihedral_atom_set const & ii_dihed( restype.dihedral( diheds[ ii ] ) );
		assert( ii_dihed.key1() == atomno || ii_dihed.key2() == atomno || ii_dihed.key3() == atomno  || ii_dihed.key4() == atomno );

		// lookup Kphi and phi0
		Real Kphi, phi0, phi_step;
		db_torsion_->lookup( res.type(), ii_dihed.key1(), ii_dihed.key2(), ii_dihed.key3(), ii_dihed.key4(), Kphi, phi0, phi_step);
		if (Kphi == 0.0) continue;

		Vector f1(0.0), f2(0.0);
		Real phi(0.0);
		if ( ii_dihed.key1() == atomno ) {
			numeric::deriv::dihedral_p1_cosine_deriv(
				res.xyz( ii_dihed.key1() ),
				res.xyz( ii_dihed.key2() ),
				res.xyz( ii_dihed.key3() ),
				res.xyz( ii_dihed.key4() ),
				phi, f1, f2 );
		} else if ( ii_dihed.key2() == atomno ) {
			numeric::deriv::dihedral_p2_cosine_deriv(
				res.xyz( ii_dihed.key1() ),
				res.xyz( ii_dihed.key2() ),
				res.xyz( ii_dihed.key3() ),
				res.xyz( ii_dihed.key4() ),
				phi, f1, f2 );
		} else if ( ii_dihed.key3() == atomno ) {
			numeric::deriv::dihedral_p2_cosine_deriv(
				res.xyz( ii_dihed.key4() ),
				res.xyz( ii_dihed.key3() ),
				res.xyz( ii_dihed.key2() ),
				res.xyz( ii_dihed.key1() ),
				phi, f1, f2 );
		} else {
			numeric::deriv::dihedral_p1_cosine_deriv(
				res.xyz( ii_dihed.key4() ),
				res.xyz( ii_dihed.key3() ),
				res.xyz( ii_dihed.key2() ),
				res.xyz( ii_dihed.key1() ),
				phi, f1, f2 );
		}

		Real del_phi = basic::subtract_radian_angles(phi, phi0);
		if (phi_step>0) del_phi = basic::periodic_range( del_phi, phi_step );
		Real dE_dphi;

		if (linear_bonded_potential_ && std::fabs(del_phi)>1) {
			dE_dphi = weights[ cart_bonded_torsion ] * Kphi * (del_phi>0? 1 : -1);
			dE_dphi += weights[ cart_bonded ] * Kphi * (del_phi>0? 1 : -1);
		} else {
			dE_dphi = weights[ cart_bonded_torsion ] * Kphi * del_phi;
			dE_dphi += weights[ cart_bonded ] * Kphi * del_phi;
		}

		LF1 += dE_dphi * f1;
		LF2 += dE_dphi * f2;
	}

	/// (A1) intra-res angles
	utility::vector1< Size > const & angs( restype.bondangles_for_atom( id.atomno() ) );
	for ( Size ii = 1, ii_end = angs.size(); ii <= ii_end; ++ii ) {
		chemical::bondangle_atom_set const & ii_bangle( restype.bondangle( angs[ ii ] ) );
		assert( ii_bangle.key1() == atomno || ii_bangle.key2() == atomno || ii_bangle.key3() == atomno );

		// lookup Ktheta and theta0
		Real Ktheta, theta0;
		db_angle_->lookup( pose, res, ii_bangle.key1(), ii_bangle.key2(), ii_bangle.key3(), Ktheta, theta0);
		if (Ktheta == 0.0) continue;

		Vector f1(0.0), f2(0.0);
		Real theta(0.0);
		if ( ii_bangle.key1() == atomno ) {
			numeric::deriv::angle_p1_deriv(
				res.xyz( ii_bangle.key1() ), res.xyz( ii_bangle.key2() ), res.xyz( ii_bangle.key3() ),
				theta, f1, f2 );
		} else if ( ii_bangle.key2() == atomno ) {
			numeric::deriv::angle_p2_deriv(
				res.xyz( ii_bangle.key1() ), res.xyz( ii_bangle.key2() ), res.xyz( ii_bangle.key3() ),
				theta, f1, f2 );
		} else {
			numeric::deriv::angle_p1_deriv(
				res.xyz( ii_bangle.key3() ), res.xyz( ii_bangle.key2() ), res.xyz( ii_bangle.key1() ),
				theta, f1, f2 );
		}

		Real dE_dtheta;

		if (linear_bonded_potential_ && std::fabs(theta - theta0)>1) {
			dE_dtheta = weights[ cart_bonded_angle ] * Ktheta * ((theta - theta0)>0? 1 : -1);
			dE_dtheta += weights[ cart_bonded ] * Ktheta * ((theta - theta0)>0? 1 : -1);
		} else {
			dE_dtheta = weights[ cart_bonded_angle ] * Ktheta * (theta - theta0);
			dE_dtheta += weights[ cart_bonded ] * Ktheta * (theta - theta0);
		}

		LF1 += dE_dtheta * f1;
		LF2 += dE_dtheta * f2;
	}

	/// (A2) intra-res bonds
	chemical::AtomIndices atm_nbrs = restype.nbrs( atomno );
	for (Size j=1; j<=atm_nbrs.size(); ++j) {
		Size atm2 = atm_nbrs[j];

		// lookup Kd and d0
		Real Kd, d0;
		db_length_->lookup( pose, res, atomno, atm2, Kd, d0 );
		if (Kd == 0.0) continue;

		Vector f1(0.0), f2(0.0);
		Real d=0;
		numeric::deriv::distance_f1_f2_deriv( res.xyz( atomno ), res.xyz( atm2 ), d, f1, f2 );

		Real dE_dd;

		if (linear_bonded_potential_ && std::fabs(d - d0)>1) {
			dE_dd = weights[ cart_bonded_length ] * Kd * ((d - d0)>0? 1 : -1);
			dE_dd += weights[ cart_bonded ] * Kd * ((d - d0)>0? 1 : -1);
		} else {
			dE_dd = weights[ cart_bonded_length ] * Kd * (d - d0);
			dE_dd += weights[ cart_bonded ] * Kd * (d - d0);
		}

		LF1 += dE_dd * f1;
		LF2 += dE_dd * f2;
	}

	/// (B) angles involving two atoms on this residue.
	utility::vector1< std::pair< Size, Size > > const & interres_wi1_for_this_atom(
		restype.within1bonds_sets_for_atom( atomno ));
	for ( Size ii = 1; ii <= interres_wi1_for_this_atom.size(); ++ii ) {
		Size const ii_resconn   = interres_wi1_for_this_atom[ ii ].first;
		Size const ii_whichpair = interres_wi1_for_this_atom[ ii ].second;
		chemical::two_atom_set const & ii_pair = restype.
			atoms_within_one_bond_of_a_residue_connection( ii_resconn )[ ii_whichpair ];
		assert( ii_pair.key1() == atomno || ii_pair.key2() == atomno );

		/// Find the neighbor residue and atom
		Size const ii_neighb      = res.residue_connection_partner( ii_resconn );

		if (ii_neighb == 0) continue;  //fpd undefined connection id

		Size const neighb_resconn = res.residue_connection_conn_id( ii_resconn );
		conformation::Residue const & neighb_res( pose.residue( ii_neighb ));
		chemical::ResidueType const & neighb_restype( pose.residue_type( ii_neighb ) );
		Size const neighb_atom    = neighb_restype.residue_connection( neighb_resconn ).atomno();

		//fpd chainbreak variants also mess things up ... explicitly check for chainbreaks
		if ( pose.fold_tree().is_cutpoint( std::min( neighb_res.seqpos(), res.seqpos() ) ) ) continue;

		// lookup Ktheta and theta0
		Real Ktheta, theta0;
		db_angle_->lookup( pose, res, -ii_resconn, ii_pair.key1(), ii_pair.key2(), Ktheta, theta0 );
		if (Ktheta == 0.0) continue;

		Vector f1(0.0), f2(0.0);
		Real theta(0.0);
		if ( ii_pair.key1() == atomno ) {
			numeric::deriv::angle_p2_deriv(
				neighb_res.xyz( neighb_atom ),
				res.xyz( ii_pair.key1() ),
				res.xyz( ii_pair.key2() ),
				theta, f1, f2 );
		} else { // ( ii_pair.key2() == atomno )
			numeric::deriv::angle_p1_deriv(
				res.xyz( ii_pair.key2() ),
				res.xyz( ii_pair.key1() ),
				neighb_res.xyz( neighb_atom ),
				theta, f1, f2 );
		}

		Real dE_dtheta;
		if (linear_bonded_potential_  && std::fabs(theta - theta0)>1 ) {
			dE_dtheta = weights[ cart_bonded_angle ] *  Ktheta * ((theta - theta0)>0? 1 : -1);
			dE_dtheta += weights[ cart_bonded ] *  Ktheta * ((theta - theta0)>0? 1 : -1);
		} else {
			dE_dtheta  = weights[ cart_bonded_angle ] * Ktheta * (theta - theta0);
			dE_dtheta  += weights[ cart_bonded ] * Ktheta * (theta - theta0);
		}

		LF1 += dE_dtheta * f1;
		LF2 += dE_dtheta * f2;
	}

	/// (C) angles involving two atoms on the neighbor residues.
	utility::vector1< Size > const & connections( restype.residue_connections_for_atom( atomno ) );
	for ( Size ii = 1; ii <= connections.size(); ++ii ) {
		Size const ii_resconn   = connections[ ii ];

		/// Find the neighbor residue and connection atom
		Size const ii_neighb      = res.residue_connection_partner( ii_resconn );

		if (ii_neighb == 0) continue;  //fpd undefined connection id

		Size const neighb_resconn = res.residue_connection_conn_id( ii_resconn );
		conformation::Residue const & neighb_res( pose.residue( ii_neighb ));
		chemical::ResidueType const & neighb_restype( pose.residue_type( ii_neighb ) );
		Size const neighb_atom1    = neighb_restype.residue_connection( neighb_resconn ).atomno();

		//fpd chainbreak variants also mess things up ... explicitly check for chainbreaks
		if ( pose.fold_tree().is_cutpoint( std::min( neighb_res.seqpos(), res.seqpos() ) ) ) continue;

		utility::vector1< chemical::two_atom_set > const & neighb_atoms_wi1(
			neighb_restype.atoms_within_one_bond_of_a_residue_connection( neighb_resconn ));

		for ( Size jj = 1; jj <= neighb_atoms_wi1.size(); ++jj ) {

			chemical::two_atom_set const & neighb_pair = neighb_atoms_wi1[ jj ];
			assert( neighb_pair.key1() == neighb_atom1 );

			Size const neighb_atom2 = neighb_pair.key2();

			Vector f1(0.0), f2(0.0);
			Real theta(0.0);

			// lookup Ktheta and theta0
			Real Ktheta, theta0;
			db_angle_->lookup( pose, neighb_res, -neighb_resconn, neighb_atom1, neighb_atom2, Ktheta, theta0 );
			if (Ktheta == 0.0) continue;

			numeric::deriv::angle_p1_deriv(
				res.xyz( atomno ),
				neighb_res.xyz( neighb_atom1 ),
				neighb_res.xyz( neighb_atom2 ),
				theta, f1, f2 );

			Real dE_dtheta;
			if (linear_bonded_potential_ && std::fabs(theta - theta0)>1) {
				dE_dtheta = weights[ cart_bonded_angle ] *  Ktheta * ((theta - theta0)>0? 1 : -1);
				dE_dtheta += weights[ cart_bonded ] *  Ktheta * ((theta - theta0)>0? 1 : -1);
			} else {
				dE_dtheta  = weights[ cart_bonded_angle ] * Ktheta * (theta - theta0);
				dE_dtheta  += weights[ cart_bonded ] * Ktheta * (theta - theta0);
			}

			LF1 += dE_dtheta * f1;
			LF2 += dE_dtheta * f2;
		}
	}

	// (D) Bonds across connections
	for ( Size ii = 1; ii <= connections.size(); ++ii ) {
		Size const ii_resconn   = connections[ ii ];

		Size const ii_neighb      = res.residue_connection_partner( ii_resconn );

		if (ii_neighb == 0) continue;  //fpd undefined connection id

		Size const neighb_resconn = res.residue_connection_conn_id( ii_resconn );
		conformation::Residue const & neighb_res( pose.residue( ii_neighb ));
		chemical::ResidueType const & neighb_restype( pose.residue_type( ii_neighb ) );
		Size const neighb_atom1    = neighb_restype.residue_connection( neighb_resconn ).atomno();

		// check for vrt
		//if ( restype.atom_type(atomno).is_virtual() || neighb_restype.atom_type(neighb_atom1).is_virtual() )
		//	continue;

		//fpd chainbreak variants also mess things up ... explicitly check for chainbreaks
		if ( pose.fold_tree().is_cutpoint( std::min( neighb_res.seqpos(), res.seqpos() ) ) ) continue;

		// lookup Kd and d0
		Real Kd, d0;
		db_length_->lookup( pose, res, atomno, -ii_resconn, Kd, d0 );
		if (Kd == 0.0) continue;

		Vector f1(0.0), f2(0.0);
		Real d=0;
		numeric::deriv::distance_f1_f2_deriv( res.xyz( atomno ), neighb_res.xyz( neighb_atom1 ), d, f1, f2 );

		Real dE_dd;

		if (linear_bonded_potential_ && std::fabs(d - d0)>1) {
			dE_dd = weights[ cart_bonded_length ] * Kd * ((d - d0)>0? 1 : -1);
			dE_dd += weights[ cart_bonded ] * Kd * ((d - d0)>0? 1 : -1);
		} else {
			dE_dd = weights[ cart_bonded_length ] * Kd * (d - d0);
			dE_dd += weights[ cart_bonded ] * Kd * (d - d0);
		}

		LF1 += dE_dd * f1;
		LF2 += dE_dd * f2;
	}

	F1 += LF1;
	F2 += LF2;

	eval_improper_torsions_derivative( id, pose, domain_map, sf, weights, F1, F2);
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
