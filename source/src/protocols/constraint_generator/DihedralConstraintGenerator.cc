// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/DihedralConstraintGenerator.cc
/// @brief A cst generator that creates Dihedral constraints for specified residues using a residue selector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/constraint_generator/DihedralConstraintGenerator.hh>
#include <protocols/constraint_generator/DihedralConstraintGeneratorCreator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>
#include <protocols/constraint_generator/util.hh>



#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>

#include <protocols/rosetta_scripts/util.hh>

// Basic/Utility headers
#include <numeric/conversions.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

static basic::Tracer TR( "protocols.constraint_generator.DihedralConstraintGenerator" );

namespace protocols {
namespace constraint_generator {
using namespace core::scoring::constraints;
using namespace core::scoring::func;
using namespace core::id;
using namespace core::select::residue_selector;

protocols::constraint_generator::ConstraintGeneratorOP
DihedralConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new DihedralConstraintGenerator );
}

std::string
DihedralConstraintGeneratorCreator::keyname() const
{
	return DihedralConstraintGenerator::class_name();
}

DihedralConstraintGenerator::DihedralConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( DihedralConstraintGenerator::class_name() )
{
}

DihedralConstraintGenerator::~DihedralConstraintGenerator()
{}

protocols::constraint_generator::ConstraintGeneratorOP
DihedralConstraintGenerator::clone() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new DihedralConstraintGenerator( *this ) );
}

std::string
DihedralConstraintGenerator::class_name()
{
	return "DihedralConstraintGenerator";
}

void
DihedralConstraintGenerator::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

void
DihedralConstraintGenerator::set_torsion_type(core::id::MainchainTorsionType torsion){
	torsion_ = torsion;
	custom_torsion_.clear();
}

void
DihedralConstraintGenerator::set_sd_degree(core::Real sd){
	sd_ = sd;
}

void
DihedralConstraintGenerator::set_custom_dihedral( utility::vector1< core::id::AtomID > const & custom_torsion ){
	assert( custom_torsion.size() == 4 );

	custom_torsion_ = custom_torsion;
}

core::scoring::constraints::ConstraintCOPs
DihedralConstraintGenerator::apply( core::pose::Pose const & pose) const
{

	utility::vector1< ConstraintCOP > constraints;

	core::Real sd_rad = numeric::conversions::radians( sd_ );
	core::Real current_torsion_angle = 0;

	utility::vector1< AtomID > local_custom_torsion = custom_torsion_;

	if ( parsed_custom_torsion_ ) {
		local_custom_torsion = parse_custom_torsion( parsed_atoms_, parsed_resnums_, pose);
	}

	if ( local_custom_torsion.size() != 0 ) {

		TR << "Using a Custom Torsion" << std::endl;
		numeric::xyzVector< core::Real>  atom1_xyz = pose.residue(local_custom_torsion[1].rsd()).xyz( local_custom_torsion[1].atomno());
		numeric::xyzVector< core::Real>  atom2_xyz = pose.residue(local_custom_torsion[2].rsd()).xyz( local_custom_torsion[2].atomno());
		numeric::xyzVector< core::Real>  atom3_xyz = pose.residue(local_custom_torsion[3].rsd()).xyz( local_custom_torsion[3].atomno());
		numeric::xyzVector< core::Real>  atom4_xyz = pose.residue(local_custom_torsion[4].rsd()).xyz( local_custom_torsion[4].atomno());

		current_torsion_angle = numeric::dihedral_degrees(atom1_xyz, atom2_xyz, atom3_xyz, atom4_xyz);

		core::Real current_torsion_angle_radians = numeric::conversions::radians( current_torsion_angle);
		CircularHarmonicFuncOP circ_func( new core::scoring::func::CircularHarmonicFunc(current_torsion_angle_radians, sd_rad) );
		DihedralConstraintOP cst( new DihedralConstraint(local_custom_torsion[1], local_custom_torsion[2], local_custom_torsion[3], local_custom_torsion[4], circ_func) );

		constraints.push_back(cst);
		return constraints;
	}
	if ( selector_ == nullptr ) {
		utility_exit_with_message("ResidueSelector required for the DihedralConstraintGenerator!");
	}
	utility::vector1< bool > subset = selector_->apply( pose );
	utility::vector1< core::Size > cst_residues = selection_positions( subset );
	for ( core::Size i : cst_residues ) {
		if ( pose.residue_type( i ).is_carbohydrate() ) {
			utility::vector1< AtomID > ref_atoms = core::conformation::carbohydrates::get_reference_atoms( torsion_, pose.conformation(), i);

			current_torsion_angle = core::conformation::carbohydrates::get_glycosidic_torsion( torsion_, pose.conformation(), i );

			core::Real current_torsion_angle_radians = numeric::conversions::radians( current_torsion_angle);
			CircularHarmonicFuncOP circ_func( new core::scoring::func::CircularHarmonicFunc(current_torsion_angle_radians, sd_rad) );
			DihedralConstraintOP cst( new DihedralConstraint(ref_atoms[1], ref_atoms[2], ref_atoms[3], ref_atoms[4], circ_func) );

			constraints.push_back(cst);
		} else {

			AtomID atm1;
			AtomID atm2;
			AtomID atm3;
			AtomID atm4;

			TorsionID torsion_id = TorsionID(i, BB, core::Size( torsion_ ));
			bool fail = pose.conformation().get_torsion_angle_atom_ids(torsion_id, atm1, atm2, atm3, atm4);
			if ( ! fail ) {
				TR.Debug << "Applying dihedral constraints to residue " << i << std::endl;

				current_torsion_angle = pose.conformation().torsion_angle(atm1, atm2, atm3, atm4);

				//core::Real current_torsion_angle_radians = numeric::conversions::radians( current_torsion_angle); Already in radians
				CircularHarmonicFuncOP circ_func( new core::scoring::func::CircularHarmonicFunc(current_torsion_angle, sd_rad) );
				DihedralConstraintOP cst( new DihedralConstraint(atm1, atm2, atm3, atm4, circ_func) );
				constraints.push_back( cst );
			}
		}
	}

	return constraints;
}

void
DihedralConstraintGenerator::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data){


	set_sd_degree( tag->getOption< core::Real >("sd", sd_ ));
	if ( tag->hasOption("dihedral") ) {
		std::string dih = tag->getOption< std::string >( "dihedral" );
		set_torsion_type( parse_torsion_type( dih ) );
		parsed_custom_torsion_ = false;
	}

	if ( tag->hasOption( "dihedral_atoms") ) {
		parsed_atoms_ = tag->getOption< std::string >( "dihedral_atoms ");
		parsed_custom_torsion_ = true;
	}
	if ( tag->hasOption(" dihedral_resnums") ) {
		parsed_resnums_ = tag->getOption< std::string >( "dihedral_resnums ");
		parsed_custom_torsion_ = true;
	}

	//Input error checking.
	if ( parsed_custom_torsion_ ) {
		if ( (parsed_atoms_ == "") || ( parsed_resnums_ == "") ) {
			utility_exit_with_message(" If setting a custom dihedral angle, both dihedral_atoms and dihedral_resnums must be set!");
		} else if ( tag->hasOption("dihedral") ) {
			utility_exit_with_message(" Custom dihedral angle with a set dihedral cannot be done.  Choose one or the other.");
		} else if ( tag->hasOption("residue_selector") ) {
			utility_exit_with_message("Custom dihedral angle not currently compatible with a residue selector. ");
		}
	}

	if ( tag->hasOption( "residue_selector") ) {

		ResidueSelectorCOP selector = parse_residue_selector( tag, data );
		if ( selector ) set_residue_selector( selector );
	}
}


void
DihedralConstraintGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	DihedralConstraintGenerator::provide_xml_schema( xsd );
}

void
DihedralConstraintGenerator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	AttributeList attlist;

	attlist + XMLSchemaAttribute("dihedral", xs_string, "The cannonical dihedral being set.  Currently, only Backbone dihedrals are accepted here.  Works for proteins and carbohydrates ");

	attlist + XMLSchemaAttribute("sd", xsct_real, "The standard deviation used for the CircularHarmonic Func.  Default is 16 degrees, THis which was found by taking the mean SD of all dihedral angles of either PHI or PSI for each North (Antibody) CDR Cluster (CDRs are the main antibody loops).  This is a fairly tight constraint and allows a bit of movement while not changing overall struture much.");


	attlist + XMLSchemaAttribute("dihedral_atoms", xs_string, "Comma-separated list of atom names.  FOUR atoms used for calculation of a custom dihedral.  Must also pass dihedral_residues.");

	attlist + XMLSchemaAttribute("dihedral_residues", xs_string, "Comma-separated list of resnums.  Must be FOUR residues to define a custom dihedral. See RosettaScripts resnum definition.  Rosetta or pose numbering (ex. 23 or 42A).  Ranges are accepted as well.   These are completely arbitrary and DO NOT have to be in order. Must also pass dihedral_atoms.");

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	std::string description = "A cst generator that creates Dihedral constraints for specified residues using a residue selector.\n"
		" Uses CircularHarmonic constraints, since CircularGaussian func does not exist. \n"
		" By default, works on Protein and carbohydrate BackBone dihedrals (see dihedral option), however CUSTOM ARBITRARY DIHEDRALS can be set. \n"
		"  See the dihedral_atoms and dihedral_residues tags to set a custom dihedral.\n\n"
		" Will only work on ONE type of dihedral angle to allow complete customization.";

	protocols::rosetta_scripts::attributes_for_parse_residue_selector( attlist );
	ConstraintGeneratorFactory::xsd_constraint_generator_type_definition_w_attributes(
		xsd,
		class_name(),
		description,
		attlist );
}

///@brief Turn two strings, comma-separated with atom names and resnums into a custom dihedral vector.
utility::vector1< core::id::AtomID >
parse_custom_torsion( std::string const & atom_names_list, std::string const & resnums, core::pose::Pose const & pose ){
	TR << "Parsing Custom Torsion" << std::endl;

	utility::vector1< AtomID > local_custom_torsion;
	utility::vector1< std::string > atom_names = utility::string_split_simple( atom_names_list, ',');
	TR << utility::to_string(atom_names) << std::endl;
	assert( atom_names.size() == 4 );

	utility::vector1< core::Size > residues = core::pose::get_resnum_list_ordered( resnums, pose );
	TR << utility::to_string(residues) << std::endl;

	assert( residues.size() == 4 );

	//Create the AtomIDs.

	for ( core::Size i = 1; i <= 4; ++i ) {

		core::Size res = residues[ i ];
		AtomID atmid = AtomID( pose.residue( res).atom_index( atom_names[ i ]), res );
		local_custom_torsion.push_back( atmid );
	}
	//TR << "Created AtomIDs" << std::endl;
	return local_custom_torsion;
}

///@brief Get the MainchainTorsionType from a string (phi/psi/omega/omega2/omega3).
core::id::MainchainTorsionType
parse_torsion_type( std::string const & dih){

	std::string local_dih = utility::lower( dih );
	if ( local_dih == "phi" ) {
		return phi_dihedral;
	} else if ( local_dih == "psi" ) {
		return psi_dihedral;
	} else if ( local_dih == "omega" ) {
		return omega_dihedral;
	} else if ( local_dih == "omega2" ) {
		return omega2_dihedral ;
	} else if ( local_dih == "omega3" ) {
		return omega3_dihedral;
	} else {
		utility_exit_with_message("DihedralConstraintGenerator cannot understand torsion type.  Please set a custom torsion! "+dih);
	}
}


} //protocols
} //constraint_generator
