// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_Mg_KnowledgeBasedPotential.hh
/// @brief
/// @author Rhiju Das


// Unit Headers
#include <core/scoring/rna/RNA_Mg_KnowledgeBasedPotential.hh>

// Package headers
#include <core/scoring/rna/RNA_ScoringInfo.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>


static basic::Tracer tr( "core.scoring.rna.RNA_Mg_KnowledgeBasedPotential" );

namespace core {
namespace scoring {
namespace rna {

/// @details ctor
RNA_Mg_KnowledgeBasedPotential::RNA_Mg_KnowledgeBasedPotential():
	gaussian_parameter_phosphate_oxygen_( -5.93, 2.22, 0.55 ), // amplitude, center, width
	gaussian_parameter_imine_           ( -3.41, 2.40, 0.33 ),
	gaussian_parameter_exocyclic_oxygen_( -4.20, 2.33, 0.43 ), // Reduced this from 5.2 to 3.4 after seeing too many close interactions.
	gaussian_parameter_o2prime_          ( -3.88, 2.54, 0.34 ),
	gaussian_parameter_phosphate_p_     ( 5.00, 2.00, 0.25 ), //this is a penalty!
	gaussian_parameter_polar_H_         ( 5.00, 2.25, 0.50 ), //this is a penalty!
	gaussian_parameter_nonpolar_H_      ( 5.00, 2.00, 0.50 ), //this is a penalty!
	//gaussian_parameter_aromatic_H_      ( 5.00, 2.00, 0.50 ), //assume no penalty here. Just a guess.

	gaussian_parameter_phosphate_oxygen_indirect_( -1.96, 4.07, 0.51 ), // longer-range, water-mediated interactions
	gaussian_parameter_imine_indirect_           ( -1.24, 4.18, 0.28 ),
	gaussian_parameter_exocyclic_oxygen_indirect_( -1.93, 4.08, 0.45 ),
	gaussian_parameter_o2prime_indirect_          ( -1.28, 4.10, 0.37 ),

	gaussian_parameter_costheta_phosphate_oxygen_( 1.00, -0.91, 0.49 ), // amplitude, center, width of angular 'form factor'
	gaussian_parameter_costheta_imine_           ( 1.00, -0.97, 0.18 ), // note how sharp this is!
	gaussian_parameter_costheta_exocyclic_oxygen_( 1.00, -0.79, 0.34 ),
	gaussian_parameter_costheta_o2prime_          ( 1.00, -0.91, 0.49 ), // not enough stats, copy from phosphate_oxygen
	gaussian_parameter_costheta_polar_H_         ( 1.00, -1.00, 0.25 ), //guess
	gaussian_parameter_costheta_nonpolar_H_      ( 1.00, -1.00, 0.25 ), //guess

	gaussian_parameter_costheta_phosphate_oxygen_indirect_( 1.00, -0.74, 0.55 ), // amplitude, center, width of angular 'form factor'
	gaussian_parameter_costheta_imine_indirect_           ( 1.00, -0.94, 0.46 ),
	gaussian_parameter_costheta_exocyclic_oxygen_indirect_( 1.00, -0.56, 0.61 ),
	gaussian_parameter_costheta_o2prime_indirect_          ( 1.00, -0.74, 0.55 ) // not enough stats, copy from phosphate_oxygen

{
}

//////////////////////////////////////////////////////////////////
GaussianParameter
RNA_Mg_KnowledgeBasedPotential::get_mg_potential_gaussian_parameter( core::conformation::Residue const & rsd, Size const j ) const{
	bool is_phosphate_oxygen( false );
	return get_mg_potential_gaussian_parameter( rsd, j, is_phosphate_oxygen );
}


/////////////////////////////////
// This is actually pretty general, and would even work for non-RNA residues, but
// for now stick to RNA where I've tried to derive a reasonable low resolution potential
// from the available PDB statistics in the RNA09 set.
GaussianParameter
RNA_Mg_KnowledgeBasedPotential::get_mg_potential_gaussian_parameter( core::conformation::Residue const & rsd, Size const j, bool & is_phosphate_oxygen ) const{

	is_phosphate_oxygen = false;

	if ( rsd.is_RNA() ) {

		std::string const atom_type_name = rsd.atom_type( j ).name();

		//tr << "atom check: " << j << " " << atom_type_name << "." << std::endl;

		// This information should probably go into a special RNA_Mg_Potential.cc function or something.
		if ( atom_type_name == "OOC" ) { // This is a  OP2 or OP1 nonbridging phosphate oxygen
			is_phosphate_oxygen = true;
			return gaussian_parameter_phosphate_oxygen_;
		} else if ( atom_type_name == "Nhis" ) { // imine nitrogens
			return gaussian_parameter_imine_;
		} else if ( atom_type_name == "OCbb" ) { // exocyclic oxygens
			return gaussian_parameter_exocyclic_oxygen_;
		} else if ( rsd.atom_name( j ) == " O2'" ){
			return gaussian_parameter_o2prime_;
		} else if ( rsd.atom_name( j ) == " P  " ){
			return gaussian_parameter_phosphate_p_;
		} else if ( atom_type_name == "Hpol" ){
			return gaussian_parameter_polar_H_;
		} else if ( atom_type_name == "Hapo" ) { //|| atom_type_name == "Haro" ){
			return gaussian_parameter_nonpolar_H_;
		}
	}

	return GaussianParameter( 0.0, 0.0, 0.0 );
}


/////////////////////////////////
GaussianParameter
RNA_Mg_KnowledgeBasedPotential::get_mg_potential_indirect_gaussian_parameter( core::conformation::Residue const & rsd, Size const j ) const{

	if ( rsd.is_RNA() ) {

		std::string const atom_type_name = rsd.atom_type( j ).name();

		// This information should probably go into a special RNA_Mg_Potential.cc function or something.
		if ( atom_type_name == "OOC" ) { // This is a  OP2 or OP1 nonbridging phosphate oxygen
			return gaussian_parameter_phosphate_oxygen_indirect_;
		} else if ( atom_type_name == "Nhis" ) { // imine nitrogens
			return gaussian_parameter_imine_indirect_;
		} else if ( atom_type_name == "OCbb" ) { // exocyclic oxygens
			return gaussian_parameter_exocyclic_oxygen_indirect_;
		} else if ( rsd.atom_name( j ) == " O2'" ){
			return gaussian_parameter_o2prime_indirect_;
		}
	}

	return GaussianParameter( 0.0, 0.0, 0.0 );
}


/////////////////////////////////
GaussianParameter
RNA_Mg_KnowledgeBasedPotential::get_mg_potential_costheta_gaussian_parameter( core::conformation::Residue const & rsd, Size const j ) const{

	if ( rsd.is_RNA() ) {

		std::string const atom_type_name = rsd.atom_type( j ).name();

		// This information should probably go into a special RNA_Mg_Potential.cc function or something.
		if ( atom_type_name == "OOC" ) { // This is a  OP2 or OP1 nonbridging phosphate oxygen
			return gaussian_parameter_costheta_phosphate_oxygen_;
		} else if ( atom_type_name == "Nhis" ) { // imine nitrogens
			return gaussian_parameter_costheta_imine_;
		} else if ( atom_type_name == "OCbb" ) { // exocyclic oxygens
			return gaussian_parameter_costheta_exocyclic_oxygen_;
		} else if ( rsd.atom_name( j ) == " O2'" ){
			return gaussian_parameter_costheta_o2prime_;
		} else if ( atom_type_name == "Hpol" ){
			return gaussian_parameter_costheta_polar_H_;
		} else if ( atom_type_name == "Hapo" ) { //|| atom_type_name == "Haro" ){
			return gaussian_parameter_costheta_nonpolar_H_;
		}
	}

	return GaussianParameter( 0.0, 0.0, 0.0 );
}

/////////////////////////////////
GaussianParameter
RNA_Mg_KnowledgeBasedPotential::get_mg_potential_costheta_indirect_gaussian_parameter( core::conformation::Residue const & rsd, Size const j ) const{

	if ( rsd.is_RNA() ) {

		std::string const atom_type_name = rsd.atom_type( j ).name();

		// This information should probably go into a special RNA_Mg_Potential.cc function or something.
		if ( atom_type_name == "OOC" ) { // This is a  OP2 or OP1 nonbridging phosphate oxygen
			return gaussian_parameter_costheta_phosphate_oxygen_indirect_;
		} else if ( atom_type_name == "Nhis" ) { // imine nitrogens
			return gaussian_parameter_costheta_imine_indirect_;
		} else if ( atom_type_name == "OCbb" ) { // exocyclic oxygens
			return gaussian_parameter_costheta_exocyclic_oxygen_indirect_;
		} else if ( rsd.atom_name( j ) == " O2'" ){
			return gaussian_parameter_costheta_o2prime_indirect_;
		}
	}

	return GaussianParameter( 0.0, 0.0, 0.0 );
}


//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_Mg_KnowledgeBasedPotential::setup_info_for_mg_calculation( pose::Pose & pose ) const
{
	//We don't know a priori which atom numbers correspond to which
	// atom names (e.g., O2' on an adenosine could be different depending
	// on whether its at a chainbreak, terminus, etc.)
	//Better to do a quick setup every time to pinpoint atoms that require
	//  monitoring for Mg binding.

	rna::RNA_ScoringInfo & rna_scoring_info( rna::nonconst_rna_scoring_info_from_pose( pose ) );

	if ( rna_scoring_info.mg_calculation_annotated_sequence() == pose.annotated_sequence() ) {
		//		tr << "matching sequence -- early return from setup_info_for_mg_calculation" << std::endl;
		return; // should be up to date
	}
	rna_scoring_info.set_mg_calculation_annotated_sequence( pose.annotated_sequence() );

	utility::vector1< bool > & is_magnesium = rna_scoring_info.nonconst_is_magnesium();
	utility::vector1< utility::vector1< Size > > &
			atom_numbers( rna_scoring_info.nonconst_atom_numbers_for_mg_calculation() );

	Size const total_residue( pose.total_residue() );

	atom_numbers.resize( total_residue );
	is_magnesium.resize( total_residue );

	for ( Size i = 1; i <= total_residue; i++ ) {

		conformation::Residue const & rsd( pose.residue( i ) );
		is_magnesium[ i ] = false;
		atom_numbers[ i ].clear();

		if ( rsd.is_RNA() ) {


			// we go over all atoms, because we are putting in some repulsions (from polar hydrogens & phosphorus)
			for ( Size j = 1; j <= rsd.natoms(); j++ ){
				GaussianParameter gaussian_parameter = get_mg_potential_gaussian_parameter( rsd, j );
				//tr << j << rsd.atom_name(j) << " ==> gaussian parameter " << gaussian_parameter.center << std::endl;
				if ( gaussian_parameter.center > 0.0 ) atom_numbers[ i ].push_back( j );
			}

		} else if ( rsd.name3() == " MG" ){
			is_magnesium[ i ] = true;
		}

		//tr << rsd.name3() << ' ' << i << ' ' << is_magnesium[ i ] << std::endl;

	}
}


} // namespace rna
} // namespace scoring
} // namespace core

