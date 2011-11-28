// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin AtomType
///
/// @brief
/// A class for defining atom parameters, known as atom_types
///
/// @detailed
/// This class contains the "chemical" information for atoms. This does not contain the actual
/// xyz coordinates of the class (xyz found in core/conformation/Atom.hh. The atom_type properties
/// are assigned by the class AtomTypeSet which is initiated from the ChemicalManager. Atom type properties
/// are currently are read in from the file located chemical/atom_type_sets/fa_standard/atom_properties.txt.
/// These properties contain the the properties of LJ_RADIUS, LJ_WDEPTH, LK_DGRFREE, LK_LAMBDA, LK_VOLUME.
/// These properties are used in the scoring function fa_atr, fa_rep, fa_sol, which is located in the Etable (core/scoring/etable/Etable.hh)
/// Additional parameters are acceptor/donor, hybridzation, and orbital paramaters.
///
///
///
/// @authors
/// Phil Bradley
/// Steven Combs - comments
///
///
/// @last_modified December 6 2010
/////////////////////////////////////////////////////////////////////////






// Rosetta headers
#include <core/chemical/AtomType.hh>


// ObjexxFCL headers
//#include <ObjexxFCL/ObjexxFCL.hh>
//#include <ObjexxFCL/string.functions.hh>

// Numeric headers


// Utility headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>

// C++ headers


namespace core {
namespace chemical {

/// @details S-H bond length in CYS.
Real const MAX_CHEMICAL_BOND_TO_HYDROGEN_LENGTH = { 1.35 };

///@brief is atom type virtual?
bool AtomType::is_virtual() const
{
	assert( (name_ == "VIRT") ? atom_is_virtual_ : true ); // Raise an error if an atom type named VIRT is not virtual.
  return (atom_is_virtual_);
}

///////////////////////////////////////////////////////////////////////////////
/// @brief set LJ and LK solvation parameter for this atom type
///
/// @details currently parameters are "LJ_RADIUS","LJ_WDEPTH","LK_VOLUME",
/// "LK_DGFREE","LK_LAMBDA".It will abort if the parameter name is not
/// recoganized. Supplemented by membrane specific solvation parameters:
/// "MEMB_LK_DGFREE","MEMB_LK_DGREFCE","LK_DGREFCE". These are the header files
/// in atom_properties.txt
///
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
/// "RING_HYBRID". It will abort if the property name is not recoganized. To add
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
	} else if(property == "AROMATIC"){
		is_aromatic_ = setting;
	} else if ( property == "H2O" ) {
		is_h2o_ = setting;
	} else if (property == "ORBITALS"){ //is the atom type orbital? defined in atom_properties.txt
		atom_has_orbitals_ = setting;
	} else if(property == "VIRTUAL"){ //is the atom type virtual? defined in atom_properties.txt
		atom_is_virtual_ = setting;
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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


} // pose
} // core
