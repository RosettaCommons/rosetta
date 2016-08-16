// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
/// @file AtomType.cc
///
/// @brief
/// A class for defining atom parameters, known as atom_types
///
/// @details
/// This class contains the "chemical" information for atoms. This does not contain the actual
/// xyz coordinates of the class (xyz found in core/conformation/Atom.hh. The atom_type properties
/// are assigned by the class AtomTypeSet which is initiated from the ChemicalManager. Atom type properties
/// are currently are read in from the file located chemical/atom_type_sets/fa_standard/atom_properties.txt.
/// These properties contain the the properties of LJ_RADIUS, LJ_WDEPTH, LK_DGRFREE, LK_LAMBDA, LK_VOLUME.
/// These properties are used in the scoring function fa_atr, fa_rep, fa_sol, which is located in the Etable
/// (core/scoring/etable/Etable.hh)
/// Additional parameters are acceptor/donor, hybridization, and orbital parameters.
///
/// @author Phil Bradley
/// @author Steven Combs - comments
/////////////////////////////////////////////////////////////////////////


// Rosetta headers
#include <core/chemical/AtomType.hh>

// Utility headers
#include <utility/exit.hh>

namespace core {
namespace chemical {

/// @details S-H bond length in CYS.
Real const MAX_CHEMICAL_BOND_TO_HYDROGEN_LENGTH = { 1.35 };

void
AtomType::print(
	std::ostream & out
) const {

	out
		<< "Atom Type: " << name() << std::endl
		<< "\telement: " << element() << std::endl
		<< "\tLennard Jones: radius=" << lj_radius() << " wdepth=" << lj_wdepth() << std::endl
		<< "\tLazaridis Karplus: lambda=" << lk_lambda() << " "
		<< "volume=" << lk_volume() << " "
		<< "dgfree=" << lk_dgfree() << std::endl
		<< "\tproperties: "
		<< (is_acceptor() ? "ACCEPTOR " : "")
		<< (is_donor() ? "DONOR " : "")
		<< (is_polar_hydrogen() ? "POLAR_HYDROGEN " : "")
		<< (is_h2o() ? "H2O " : "")
		<< (is_aromatic() ? "AROMATIC " : "")
		<< (atom_has_orbital() ? "ORBITALS " : "");
	switch(hybridization()){
	case SP2_HYBRID : out << "SP2_HYBRID "; break;
	case SP3_HYBRID : out << "SP3_HYBRID "; break;
	case RING_HYBRID : out << "RING_HYBRID "; break;
	case UNKNOWN_HYBRID : break;
	default :
		utility_exit_with_message("Attempting retrive hydrid for atom type '" + name() +
			"', however the hybridization type is not recognized.");
	}
	out << std::endl;
	out << "Extra Parameters:";
	for ( Size i = 1; i <= extra_parameters_.size(); ++i ) {
		out << " " << extra_parameters_[i];
	}
	out << std::endl;
}

std::ostream &
operator<< (std::ostream & out, const AtomType & atom_type ){
	atom_type.print( out );
	return out;
}


///////////////////////////////////////////////////////////////////////////////
/// @brief set LJ and LK solvation parameter for this atom type
///
/// @details currently parameters are "LJ_RADIUS","LJ_WDEPTH","LK_VOLUME",
/// "LK_DGFREE","LK_LAMBDA".It will abort if the parameter name is not
/// Recognized. Supplemented by membrane specific solvation parameters:
/// "MEMB_LK_DGFREE","MEMB_LK_DGREFCE","LK_DGREFCE". These are the header files
/// in atom_properties.txt
void
AtomType::set_parameter(
	std::string const & param,
	Real const setting
)
{
	if ( param == "LJ_RADIUS" ) {
		lj_radius_ = setting;
	} else if ( param == "LJ_WDEPTH" ) {
		lj_wdepth_ = setting;
	} else if ( param == "LK_VOLUME" ) {
		lk_volume_ = setting;
	} else if ( param == "LK_DGFREE" ) {
		lk_dgfree_ = setting;
	} else if ( param == "LK_LAMBDA" ) {
		lk_lambda_ = setting;
		/*} else if ( param == "MEMB_LK_DGFREE" ) {  //pba
		memb_lk_dgfree_ = setting;
		} else if ( param == "LK_DGREFCE" ) {      //pba
		lk_dgrefce_ = setting;
		} else if ( param == "MEMB_LK_DGREFCE" ) { //pba
		memb_lk_dgrefce_ = setting;*/
	} else {
		utility_exit_with_message( "unrecognized atomtype parameter "+param );
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @brief set relevant properties for this atom type
///
/// @details currently properties are "ACCEPTOR","DONOR","POLAR_HYDROGEN",
/// "H2O", and hybridization types including "SP2_HYBRID", "SP3_HYBRID" and
/// "RING_HYBRID". It will abort if the property name is not recognized. To add
/// properties, edit atom_properties.txt and add your property to the last column
/// then add code here that will read the property.
///
void
AtomType::set_property(
	std::string const & property,
	bool const setting
)
{
	if ( property == "ACCEPTOR" ) {
		is_acceptor_ = setting;
	} else if ( property == "DONOR" ) {
		is_donor_ = setting;
	} else if ( property == "POLAR_HYDROGEN" ) {
		is_polar_hydrogen_ = setting;
	} else if ( property == "AROMATIC" ) {
		is_aromatic_ = setting;
	} else if ( property == "H2O" ) {
		is_h2o_ = setting;
	} else if ( property == "ORBITALS" ) { //is the atom type orbital? defined in atom_properties.txt
		atom_has_orbitals_ = setting;
	} else if ( property == "VIRTUAL" ) { //is the atom type virtual? defined in atom_properties.txt
		atom_is_virtual_ = setting;
	} else if ( property == "REPULSIVE" ) {
		atom_is_repulsive_ = setting;
	} else if ( property == "SP2_HYBRID" ) {
		hybridization_ = SP2_HYBRID;
	} else if ( property == "SP3_HYBRID" ) {
		hybridization_ = SP3_HYBRID;
	} else if ( property == "RING_HYBRID" ) {
		hybridization_ = RING_HYBRID;
	} else {
		utility_exit_with_message( "unrecognized atomtype property "+property );
	}
}

void
AtomType::clear_properties() {
	is_acceptor_ = false;
	is_donor_ = false;
	is_polar_hydrogen_ = false;
	is_aromatic_ = false;
	is_h2o_ = false;
	atom_has_orbitals_ = false;
	atom_is_virtual_ = false;
	hybridization_ = UNKNOWN_HYBRID;
	extra_parameters_.clear();
}

void
AtomType::add_property(
	std::string const & property
) {
	if ( property == "ACCEPTOR" ) {
		is_acceptor_ = true;
	} else if ( property == "DONOR" ) {
		is_donor_ = true;
	} else if ( property == "POLAR_HYDROGEN" ) {
		is_polar_hydrogen_ = true;
	} else if ( property == "AROMATIC" ) {
		is_aromatic_ = true;
	} else if ( property == "H2O" ) {
		is_h2o_ = true;
	} else if ( property == "ORBITALS" ) {
		atom_has_orbitals_ = true;
	} else if ( property == "VIRTUAL" ) {
		atom_is_virtual_ = true;
	} else if ( property == "REPULSIVE" ) {
		atom_is_repulsive_ = true;
	} else if ( property == "SP2_HYBRID" ) {
		hybridization_ = SP2_HYBRID;
	} else if ( property == "SP3_HYBRID" ) {
		hybridization_ = SP3_HYBRID;
	} else if ( property == "RING_HYBRID" ) {
		hybridization_ = RING_HYBRID;
	} else {
		utility_exit_with_message("Attempting to set non-existant property '" + property +
			"' on atom type '" + name() + "'.");
	}
}


utility::vector1< std::string >
AtomType::get_all_properties() const {
	utility::vector1< std::string > properties;
	if ( is_acceptor() ) properties.push_back("ACCEPTOR");
	if ( is_donor() ) properties.push_back("DONOR");
	if ( is_polar_hydrogen() ) properties.push_back("POLAR_HYDROGEN");
	if ( is_aromatic() ) properties.push_back("AROMATIC");
	if ( is_h2o() ) properties.push_back("H2O");
	if ( atom_has_orbital() ) properties.push_back("ORBITALS");
	if ( is_virtual() ) properties.push_back("VIRTUAL");
	if ( is_repulsive() ) properties.push_back("REPULSIVE");

	switch(hybridization()){
	case SP2_HYBRID : properties.push_back("SP2_HYBRID"); break;
	case SP3_HYBRID : properties.push_back("SP3_HYBRID"); break;
	case RING_HYBRID : properties.push_back("RING_HYBRID"); break;
	case UNKNOWN_HYBRID : break;
	default :
		utility_exit_with_message("Attempting retrive hydrid for atom type '" + name() +
			"', however the hybridization type is not recognized.");
	}
	return properties;
}

void
AtomType::set_all_extra_parameters(
	utility::vector1< Real > const & extra_parameters
) {
	extra_parameters_ = extra_parameters;
}


} // pose
} // core
