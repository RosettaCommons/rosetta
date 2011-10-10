// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CartesianBondedEnergy.cc
/// @brief
/// @author

// Unit headers
#include <core/scoring/methods/CartesianBondedEnergy.hh>
#include <core/scoring/methods/CartesianBondedEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/kinematics/AtomTree.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>  // for is_virtual()
#include <core/chemical/util.hh>  // for is_virtual()
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>

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

// C++ headers
#include <iostream>

//Auto Headers
#include <core/id/AtomID.hh>



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


static basic::Tracer TR("core.scoring.CartesianBondedEnergy");

// default spring constants
static const Real K_ANGLE=600.0;
static const Real K_BOND=600.0;
static const Real K_TORSION=300.0;
static const Real K_TORSION_PROTON=60.0;  // proton chi

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
BondLengthDatabase::BondLengthDatabase() {
	k_bond = K_BOND;
	if (basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user()) {
		utility::vector1<core::Real> params = basic::options::option[ basic::options::OptionKeys::score::bonded_params ]();
		k_bond = params[1];
	}
}

BondAngleDatabase::BondAngleDatabase() {
	k_angle = K_ANGLE;
	if (basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user()) {
		utility::vector1<core::Real> params = basic::options::option[ basic::options::OptionKeys::score::bonded_params ]();
		k_angle = params[2];
	}
}

TorsionDatabase::TorsionDatabase() {
	k_torsion = K_TORSION;
	k_torsion_proton = K_TORSION_PROTON;
	if (basic::options::option[ basic::options::OptionKeys::score::bonded_params ].user()) {
		utility::vector1<core::Real> params = basic::options::option[ basic::options::OptionKeys::score::bonded_params ]();
		k_torsion = params[3];
		k_torsion_proton = params[4];
	}
}


//////////////////////
/// Torsion Database
// lookup ideal torsion; insert in DB if not there
void
TorsionDatabase::lookup(
		core::chemical::ResidueType const & restype,
		int atm1, int atm2, int atm3, int atm4, Real &Kphi, Real &phi0, Real &phi_step ) {
	using namespace core::chemical;

	Kphi=k_torsion; // default
	phi_step=0;

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

		return;
	}

	// Create mini-conformation for idealized residue
	core::conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );
	core::pose::Pose newpose;
	newpose.append_residue_by_bond(*newres);

	// check if this corresponds to a DOF_ID
	if ( newpose.atom_tree().torsion_angle_dof_id(
				 id::AtomID(atm1,1),id::AtomID(atm2,1),id::AtomID(atm3,1),id::AtomID(atm4,1)
			 ) != id::BOGUS_DOF_ID) {
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

		////////
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
	
		//fpd  torsion_angle_dof_id should pick these up (?)
		if ( restype.aa()==core::chemical::aa_trp || restype.aa()==core::chemical::aa_phe
				 || restype.aa()==core::chemical::aa_tyr || restype.aa()==core::chemical::aa_his ) {
			if ( ( restype.atom_name(atm2) == " CB " && restype.atom_name(atm3) == " CG ") ||
					 ( restype.atom_name(atm3) == " CB " && restype.atom_name(atm2) == " CG ") ) {
				Kphi=0.0;
				phi0=0.0;
			}
		}

		//fpd  cutpoint variants
		if ( restype.has_variant_type(chemical::CUTPOINT_UPPER) &&
				 ( ( restype.atom_name(atm2) == " N  " && restype.atom_name(atm3) == " CA ") ||
					 ( restype.atom_name(atm3) == " N  " && restype.atom_name(atm2) == " CA ")  ) ) {
			Kphi=0.0;
			phi0=0.0;
		}
	
		// fpd ignore centroid torsion
		if ( restype.atom_name(atm1) == " CEN" || restype.atom_name(atm4) == " CEN") {
			Kphi=0.0;
			phi0=0.0;
		}

		//fpd  proton CHIs
		if ( restype.aa()==core::chemical::aa_cys && 
				 ( ( restype.atom_name(atm2) == " CB " && restype.atom_name(atm3) == " SG ") ||
					 ( restype.atom_name(atm3) == " CB " && restype.atom_name(atm2) == " SG ")  ) ) {
			phi_step = numeric::constants::f::pi_2_over_3;
			torsion_steps_[ tuple ] = phi_step;
			Kphi=k_torsion_proton;
		}
		if (restype.aa()==core::chemical::aa_thr &&
				 ( ( restype.atom_name(atm2) == " CB " && restype.atom_name(atm3) == " OG1") ||
					 ( restype.atom_name(atm3) == " CB " && restype.atom_name(atm2) == " OG1")  ) ) {
			phi_step = numeric::constants::f::pi_2_over_3;
			torsion_steps_[ tuple ] = phi_step;
			Kphi=k_torsion_proton;
		}
		if (restype.aa()==core::chemical::aa_ser &&
				 ( ( restype.atom_name(atm2) == " CB " && restype.atom_name(atm3) == " OG ") ||
					 ( restype.atom_name(atm3) == " CB " && restype.atom_name(atm2) == " OG ")  ) ) {
			phi_step = numeric::constants::f::pi_2_over_3;
			torsion_steps_[ tuple ] = phi_step;
			Kphi=k_torsion_proton;
		}
		if (restype.aa()==core::chemical::aa_tyr &&
				 ( ( restype.atom_name(atm2) == " CZ " && restype.atom_name(atm3) == " OH ") ||
					 ( restype.atom_name(atm3) == " OH " && restype.atom_name(atm2) == " CZ ")  ) ) {
			phi_step = numeric::constants::f::pi;
			torsion_steps_[ tuple ] = phi_step;
			Kphi=k_torsion_proton;
		}
	}
	torsions_[ tuple ] = phi0;
	Kphis_[ tuple ] = Kphi;

	return;
}
//////////////////////
/// BondAngle Database
// lookup ideal bondangle; insert in DB if not there
// atm ids > 0 ==> atom index
// atm ids < 0 ==> residue connection id
void
BondAngleDatabase::lookup(
		core::chemical::ResidueType const & restype, int atm1, int atm2, int atm3, Real &Ktheta, Real &theta0 ) {
	using namespace core::chemical;

	Ktheta=k_angle;

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
	core::conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );

	// get angle
	numeric::xyzVector<core::Real> x,y,z; // atom coords
	if (atm1 > 0)
		x = newres->atom( atm1 ).xyz();
	else
		x = newres->residue_connection( -atm1 ).icoor().build( restype );
	if (atm2 > 0)
		y = newres->atom( atm2 ).xyz();
	else
		y = newres->residue_connection( -atm2 ).icoor().build( restype );
	if (atm3 > 0)
		z = newres->atom( atm3 ).xyz();
	else
		z = newres->residue_connection( -atm3 ).icoor().build( restype );

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
  ( core::chemical::ResidueType const & restype, int atm1, int atm2, Real &Kd, Real &d0 ) {
	using namespace core::chemical;

	Kd=k_bond;

	std::string restag = get_restag( restype );

	// lookup in table
	residx_atm_pair tuple( restag, atm1,atm2 );
	boost::unordered_map<residx_atm_pair,core::Real>::iterator b_it = bondlengths_.find( tuple );

	if ( b_it != bondlengths_.end() ) {
		d0 = b_it->second;
		boost::unordered_map<residx_atm_pair,core::Real>::iterator k_it = Kbonds_.find( tuple );
		if ( k_it != Kbonds_.end() )
			Kd = k_it->second;
		return;
	}

	// Create mini-conformation for idealized residue
	core::conformation::ResidueOP newres = conformation::ResidueFactory::create_residue( restype );

	// get length
	numeric::xyzVector<core::Real> x,y,z; // atom coords
	if (atm1 > 0)
		x = newres->atom( atm1 ).xyz();
	else
		x = newres->residue_connection( -atm1 ).icoor().build( restype );
	if (atm2 > 0)
		y = newres->atom( atm2 ).xyz();
	else
		y = newres->residue_connection( -atm2 ).icoor().build( restype );

	d0 = (x-y).length();

	// fpd ignore proline C-ND which is handled by pro_close
	if (restype.aa() == core::chemical::aa_pro &&
	    (atm1>0 && atm2>0) &&
	    ( ( restype.atom_name(atm1) == " CD " && restype.atom_name(atm2) == " N  ") ||
	      ( restype.atom_name(atm2) == " CD " && restype.atom_name(atm1) == " N  ") ) ) {
		Kd = d0 = 0.0;
	}

	bondlengths_[ tuple ] = d0;
	Kbonds_[ tuple ] = Kd;

	return;
}


//////////////////////
/// EnergyMethod
CartesianBondedEnergy::CartesianBondedEnergy( methods::EnergyMethodOptions const & options ) :
	parent( new CartesianBondedEnergyCreator ) {
	linear_bonded_potential_ = basic::options::option[ basic::options::OptionKeys::score::linear_bonded_potential ]();
}

CartesianBondedEnergy::CartesianBondedEnergy( CartesianBondedEnergy const & src ) : parent( src ) {
	linear_bonded_potential_ = src.linear_bonded_potential_;
	db_angle_ = src.db_angle_;
	db_length_ = src.db_length_;
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
		LREnergyContainerOP new_dec = new PeptideBondedEnergyContainer( pose.total_residue(), cart_bonded );
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
	pose::Pose const & pose,
	Size res1,
	Size res2
) const {
	// is this fn. called?
	return ( res1 == (res2+1) || res1 == (res2-1) );
}


void
CartesianBondedEnergy::residue_pair_energy(
   conformation::Residue const & rsd1,
	 conformation::Residue const & rsd2,
	 pose::Pose const & pose,
	 ScoreFunction const & ,
	 EnergyMap & emap
) const {
	// bail out if the residues aren't bonded
	if (!rsd1.is_bonded(rsd2)) return;

	Real energy = 0;

	// get residue types
	chemical::ResidueType const & rsd1_type = rsd1.type();
	chemical::ResidueType const & rsd2_type = rsd2.type();

	utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );

	for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {

		Size const resconn_id1( r1_resconn_ids[ii] );
		Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

		Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
		Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

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
			db_angle_.lookup( rsd1.type(), res1_lower_atomno, resconn_atomno1, -resconn_id1, Ktheta, theta0 );

			if (Ktheta == 0.0) continue;

			//if (0.5*Ktheta*(angle-theta0)*(angle-theta0) > 10.0) {
			//	TR << rsd1.name() << ":" << rsd1_type.atom_name( res1_lower_atomno ) << " , "
			//	   << rsd1_type.atom_name( resconn_atomno1 ) << " , " << rsd2_type.atom_name( resconn_atomno2 )
			//	   << "    " << angle << "  " << theta0 << "   "
			//	   << "    " << 0.5*Ktheta*(angle-theta0)*(angle-theta0) << std::endl;
			//}

			// accumulate the energy
			if (linear_bonded_potential_ && std::fabs(angle-theta0)>1) {
				energy += 0.5*Ktheta*std::fabs(angle-theta0);
			} else {
				energy += 0.5*Ktheta*(angle-theta0) * (angle-theta0);
			}
		}

		/// compute the bond-angle energies from pairs of atoms within-1 bond on rsd2 with
		/// the the connection atom on rsd1.
		utility::vector1< chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
			rsd2_type.atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));
		for ( Size jj = 1; jj <= rsd2_atoms_wi1_bond_of_ii.size(); ++jj ) {
			assert( rsd2_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno2 );
			Size const res2_lower_atomno = rsd2_atoms_wi1_bond_of_ii[ jj ].key2();

			// lookup Ktheta and theta0
			Real Ktheta, theta0;
			db_angle_.lookup( rsd2.type(), res2_lower_atomno, resconn_atomno2, -resconn_id2, Ktheta, theta0 );

			if (Ktheta == 0.0) continue;
			Real const angle = numeric::angle_radians(
				rsd2.atom( res2_lower_atomno ).xyz(),
				rsd2.atom( resconn_atomno2 ).xyz(),
				rsd1.atom( resconn_atomno1 ).xyz() );

			//if (0.5*Ktheta*(angle-theta0)*(angle-theta0) > 10.0) {
			//	TR << rsd2.name() << ":" << rsd2_type.atom_name( res2_lower_atomno ) << " , "
			//	   << rsd2_type.atom_name( resconn_atomno2 ) << " , " << rsd1_type.atom_name( resconn_atomno1 )
			//	   << "    " << 0.5*Ktheta*(angle-theta0)*(angle-theta0) << std::endl;
			//}

			// accumulate the energy
			if (linear_bonded_potential_ && std::fabs(angle-theta0)>1) {
				energy += 0.5*Ktheta*std::fabs(angle-theta0);
			} else {
				energy += 0.5*Ktheta*(angle-theta0) * (angle-theta0);
			}
		}

		/// finally, compute the bondlength across the interface
		Real const length =
			( rsd2.atom( resconn_atomno2 ).xyz() - rsd1.atom( resconn_atomno1 ).xyz() ).length();

		// lookup Ktheta and theta0
		Real Kd, d0;
		db_length_.lookup( rsd1.type(), resconn_atomno1, -resconn_id1, Kd, d0 );

		//if (0.5*Kd*(length-d0) * (length-d0) > 10.0) {
		//	TR << rsd1.seqpos() << " -- " << rsd2.seqpos() << "  "
		//	   << rsd1.name() << ":" << rsd1_type.atom_name( resconn_atomno1 ) << " , " << rsd2_type.atom_name( resconn_atomno2 )
		//	   << "    " << length << " [" << d0 << "]" 
		//	   << "    " << 0.5*Kd*std::fabs(length-d0) << std::endl;
		//}

		// accumulate the energy
		if (linear_bonded_potential_ && std::fabs(length-d0)>1) {
			energy += 0.5*Kd*std::fabs(length-d0);
		} else {
			energy += 0.5*Kd*(length-d0)*(length-d0);
		}
	}

	emap[ cart_bonded ] += energy;
}

void
CartesianBondedEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & ,
	EnergyMap & emap
) const
{
	Real energy = 0;

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
		db_torsion_.lookup( rsd.type(), rt1, rt2, rt3, rt4, Kphi, phi0, phi_step );
		if (Kphi == 0.0) continue;

		// get angle
		Real angle = numeric::dihedral_radians
			( rsd.atom( rt1 ).xyz(), rsd.atom( rt2 ).xyz(),
				rsd.atom( rt3 ).xyz(), rsd.atom( rt4 ).xyz() );

		// accumulate the energy
		Real del_phi = basic::subtract_radian_angles(angle, phi0);
		if (phi_step>0) del_phi = basic::periodic_range( del_phi, phi_step );

 		//if (0.5*Kphi*del_phi*del_phi > 10) {
 		//	TR <<
 		//		rsd_type.name() << " : " <<
 		//	  rsd.atom_name( rt1 ) << " , " << rsd.atom_name( rt2 ) << " , " <<
 		//	  rsd.atom_name( rt3 ) << " , " << rsd.atom_name( rt4 ) << "   " <<
 		//	  0.5*Kphi*del_phi*del_phi << std::endl;
 		//}

		if (linear_bonded_potential_ && std::fabs(del_phi)>1) {
			energy += 0.5*Kphi*std::fabs(del_phi);
		} else {
			energy += 0.5*Kphi*del_phi*del_phi;
		}
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
		db_angle_.lookup( rsd.type(), rt1, rt2, rt3, Ktheta, theta0 );
		if (Ktheta == 0.0) continue;

		// get angle
		Real const angle = numeric::angle_radians(
			rsd.atom( rt1 ).xyz(),
			rsd.atom( rt2 ).xyz(),
			rsd.atom( rt3 ).xyz() );

 		//if ( 0.5*Ktheta*(angle-theta0) * (angle-theta0) > 10) {
 		//	TR <<
 		//		rsd_type.name() << " : " <<
 		//	  rsd.atom_name( rt1 ) << " , " << rsd.atom_name( rt2 ) << " , " <<
 		//	  rsd.atom_name( rt3 ) << "   " << angle << "  " << theta0 << "     " <<
 		//	  0.5*Ktheta*(angle-theta0) * (angle-theta0) << std::endl;
 		//}
		// accumulate the energy
		if (linear_bonded_potential_ && std::fabs(angle - theta0)>1) {
			energy += 0.5*Ktheta*std::fabs(angle-theta0);
		} else {
			energy += 0.5*Ktheta*(angle-theta0) * (angle-theta0);
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

				// lookup Ktheta and theta0
				Real Kd, d0;
				db_length_.lookup( rsd.type(), atm_i, atm_j, Kd, d0 );
				if (Kd == 0.0) continue;

				Real const d = ( rsd.atom( atm_i ).xyz()-rsd.atom( atm_j ).xyz() ).length();

				//if ( 0.5*Kd*(d-d0)*(d-d0) > 10) {
				//	TR <<
				//		rsd_type.name() << " : " <<
				//		rsd.atom_name( atm_i ) << " , " << rsd.atom_name( atm_j ) << "   " <<
				//		0.5*Kd*(d-d0)*(d-d0) << std::endl;
				//}

				// accumulate the energy
				if (linear_bonded_potential_ && std::fabs(d - d0)>1) {
					energy += 0.5*Kd*std::fabs(d-d0);
				} else {
					energy += 0.5*Kd*(d-d0)*(d-d0);
				}
			}
		}
	}

	// add energy to emap
	emap[ cart_bonded ] += energy;
}


///////////////////////////
///
void
CartesianBondedEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
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

		// lookup Kphi and phi0
		Real Kphi, phi0, phi_step;
		db_torsion_.lookup( res.type(), ii_dihed.key1(), ii_dihed.key2(), ii_dihed.key3(), ii_dihed.key4(), Kphi, phi0, phi_step);

		if (Kphi == 0.0) continue;

		Real del_phi = basic::subtract_radian_angles(phi, phi0);
		if (phi_step>0) del_phi = basic::periodic_range( del_phi, phi_step );
		Real dE_dphi;

		if (linear_bonded_potential_ && std::fabs(del_phi)>1) {
			dE_dphi = weights[ cart_bonded ] * Kphi * (del_phi>0? 1 : -1);
		} else {
			dE_dphi = weights[ cart_bonded ] * Kphi * del_phi;
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
		db_angle_.lookup( res.type(), ii_bangle.key1(), ii_bangle.key2(), ii_bangle.key3(), Ktheta, theta0);
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
			dE_dtheta = weights[ cart_bonded ] * Ktheta * ((theta - theta0)>0? 1 : -1);
		} else {
			dE_dtheta = weights[ cart_bonded ] * Ktheta * (theta - theta0);
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
		db_length_.lookup( res.type(), atomno, atm2, Kd, d0 );
		if (Kd == 0.0) continue;

		Vector f1(0.0), f2(0.0);
		Real d=0;
		numeric::deriv::distance_f1_f2_deriv( res.xyz( atomno ), res.xyz( atm2 ), d, f1, f2 );


		Real dE_dd;

		if (linear_bonded_potential_ && std::fabs(d - d0)>1) {
			dE_dd = weights[ cart_bonded ] * Kd * ((d - d0)>0? 1 : -1);
		} else {
			dE_dd = weights[ cart_bonded ] * Kd * (d - d0);
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

		// lookup Ktheta and theta0
		Real Ktheta, theta0;
		db_angle_.lookup( res.type(), -ii_resconn, ii_pair.key1(), ii_pair.key2(), Ktheta, theta0 );
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
			dE_dtheta = weights[ cart_bonded ] *  Ktheta * ((theta - theta0)>0? 1 : -1);
		} else {
			dE_dtheta  = weights[ cart_bonded ] * Ktheta * (theta - theta0);
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
			db_angle_.lookup( neighb_res.type(), -neighb_resconn, neighb_atom1, neighb_atom2, Ktheta, theta0 );
			if (Ktheta == 0.0) continue;

			numeric::deriv::angle_p1_deriv(
				res.xyz( atomno ),
				neighb_res.xyz( neighb_atom1 ),
				neighb_res.xyz( neighb_atom2 ),
				theta, f1, f2 );


			Real dE_dtheta;
			if (linear_bonded_potential_ && std::fabs(theta - theta0)>1) {
				dE_dtheta = weights[ cart_bonded ] *  Ktheta * ((theta - theta0)>0? 1 : -1);
			} else {
				dE_dtheta  = weights[ cart_bonded ] * Ktheta * (theta - theta0);
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

		// lookup Kd and d0
		Real Kd, d0;
		db_length_.lookup( res.type(), atomno, -ii_resconn, Kd, d0 );
		if (Kd == 0.0) continue;

		Vector f1(0.0), f2(0.0);
		Real d=0;
		numeric::deriv::distance_f1_f2_deriv( res.xyz( atomno ), neighb_res.xyz( neighb_atom1 ), d, f1, f2 );

		Real dE_dd;

		if (linear_bonded_potential_ && std::fabs(d - d0)>1) {
			dE_dd = weights[ cart_bonded ] * Kd * ((d - d0)>0? 1 : -1);
		} else {
			dE_dd = weights[ cart_bonded ] * Kd * (d - d0);
		}

		LF1 += dE_dd * f1;
		LF2 += dE_dd * f2;
	}

	F1 += LF1;
	F2 += LF2;
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
