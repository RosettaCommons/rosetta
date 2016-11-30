// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Adds all the required constraints and bonds to cyclize a pose.
/// @author Parisa Hosseinzadeh (parisah@uw.edu)

//Core Headers
#include <core/id/AtomID.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/Func.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/kinematics/FoldTree.hh>

// Protocols Headers
#include <protocols/cyclic_peptide/CreateAngleConstraint.fwd.hh>
#include <protocols/cyclic_peptide/CreateAngleConstraint.hh>
#include <protocols/cyclic_peptide/CreateDistanceConstraint.fwd.hh>
#include <protocols/cyclic_peptide/CreateDistanceConstraint.hh>
#include <protocols/cyclic_peptide/CreateTorsionConstraint.fwd.hh>
#include <protocols/cyclic_peptide/CreateTorsionConstraint.hh>
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/rosetta_scripts/util.hh>

//Numeric Headers
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>

//Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

//Mover Headers
#include <protocols/cyclic_peptide/PeptideCyclizeMover.hh>
#include <protocols/cyclic_peptide/PeptideCyclizeMoverCreator.hh>

//Stream headers+basic headers
#include <basic/Tracer.hh>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <sstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


//General construction
static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.PeptideCyclizeMover" );

namespace protocols {
namespace cyclic_peptide {

moves::MoverOP PeptideCyclizeMover::clone() const
{
	return moves::MoverOP( new PeptideCyclizeMover( *this ) );
}

moves::MoverOP PeptideCyclizeMover::fresh_instance() const
{
	return moves::MoverOP( new PeptideCyclizeMover );
}

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP PeptideCyclizeMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new PeptideCyclizeMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP PeptideCyclizeMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return PeptideCyclizeMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP PeptideCyclizeMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "PeptideCyclizeMover";
// XRW TEMP }

///@brief destructor of the mover
//
PeptideCyclizeMover::~PeptideCyclizeMover(){}

/// @brief Constructor for PeptideCyclizeMover.
///
PeptideCyclizeMover::PeptideCyclizeMover() :
	Mover("PeptideCyclizeMover"),
	selector_()
{
	distance_assigned_ = false;
	bond_assigned_ = false;
	angle_assigned_ = false;
	torsion_assigned_ = false;
	//set_rosetta_scripts_tag( utility::tag::TagOP( new utility::tag::Tag() ) ); //do I need this?
}

//these are setters. Set user-defined values.
//use to set user defined bonds
void PeptideCyclizeMover::set_bond(core::Size res1, std::string atom1, core::Size res2, std::string atom2,  bool add_termini, bool rebuild_fold_tree){

	res1_.push_back(res1);
	atom1_.push_back(atom1);
	res2_.push_back(res2);
	atom2_.push_back(atom2);
	add_termini_.push_back(add_termini);
	rebuild_fold_tree_.push_back(rebuild_fold_tree);

	counter_+=1;

	//printf ("I am setting bond for you\n");
	bond_assigned_=true;

	return;
}
//use to set user defined distance constraints
void PeptideCyclizeMover::set_distance(core::Size res1, std::string atom1, core::Size res2, std::string atom2, std::string cst_func){

	res1_dist_.push_back(res1);
	atom1_dist_.push_back(atom1);
	res2_dist_.push_back(res2);
	atom2_dist_.push_back(atom2);
	cst_func_dist_.push_back(cst_func);

	distance_assigned_=true;
	//printf ("Now I am setting distance for you\n");

	return;
}
//setting user defined angle constraints
void PeptideCyclizeMover::set_angle(core::Size res_center, std::string atom_center, core::Size res1, std::string atom1, core::Size res2, std::string atom2, std::string cst_func){

	res1_angle_.push_back(res1);
	atom1_angle_.push_back(atom1);
	res2_angle_.push_back(res2);
	atom2_angle_.push_back(atom2);
	res_center_.push_back(res_center);
	atom_center_.push_back(atom_center);
	cst_func_angle_.push_back(cst_func);

	angle_assigned_=true;
	//printf ("setting angle dude\n");

	return;
}
//setting user defined torsion constaints
void PeptideCyclizeMover::set_torsion(core::Size res1, std::string atom1, core::Size res2, std::string atom2, core::Size res3, std::string atom3, core::Size res4, std::string atom4,std::string cst_func){

	res1_torsion_.push_back(res1);
	atom1_torsion_.push_back(atom1);
	res2_torsion_.push_back(res2);
	atom2_torsion_.push_back(atom2);
	res3_torsion_.push_back(res3);
	atom3_torsion_.push_back(atom3);
	res4_torsion_.push_back(res4);
	atom4_torsion_.push_back(atom4);
	cst_func_torsion_.push_back(cst_func);

	torsion_assigned_=true;
	//printf ("time for setting torsion\n");

	return;
}

//setting selectors
void PeptideCyclizeMover::set_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	if ( selector_in ) {
		selector_ = selector_in;
	} else {
		utility_exit_with_message("Error in protocols::cyclic_peptide::PeptideCyclizeMover::set_selector(): Null pointer passed to function!");
	}
	return;
}

//using this function, you can call the default function, i.e. peptide bond between the two termini of the pose.
void PeptideCyclizeMover::set_default
(){

	bond_assigned_=false;
	distance_assigned_=false;
	angle_assigned_=false;
	torsion_assigned_=false;

	return;
}


///////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.

void PeptideCyclizeMover::apply( core::pose::Pose & pose )
{
	core::select::residue_selector::ResidueSubset subset( pose.total_residue(), true );
	if ( pose.total_residue()==0 ) {
		return;
	} else {
		if ( selector_ ) {
			//printf ("setting residue selector\n");
			subset = selector_->apply( pose );
		}
	}

	get_all(subset,pose); //get all the values of constaintrs


	protocols::cyclic_peptide::CreateDistanceConstraintOP distance(new protocols::cyclic_peptide::CreateDistanceConstraint);
	distance->set(res1_dist_,atom1_dist_,res2_dist_,atom2_dist_,cst_func_dist_);
	distance -> apply(pose);
	std::cout << "setting distance " << cst_func_dist_ << " between " << res1_dist_ << atom1_dist_ << " and " << res2_dist_ << atom2_dist_ << std::endl;

	protocols::cyclic_peptide::CreateAngleConstraintOP angle1(new protocols::cyclic_peptide::CreateAngleConstraint);
	protocols::cyclic_peptide::CreateAngleConstraintOP angle2(new protocols::cyclic_peptide::CreateAngleConstraint);

	if ( !angle_assigned_ ) {
		angle2->set(res_center2_,atom_center2_,res1_angle_,atom1_angle2_,res2_angle_,atom2_angle2_,cst_func_angle2_);
		angle2->apply(pose);
		angle1->set(res_center_,atom_center_,res1_angle_,atom1_angle_,res2_angle_,atom2_angle_,cst_func_angle_);
		angle1->apply(pose);
		//printf("this is where angle constraints being set, if not assigned\n");
		std::cout << "we are first setting the angle function " << cst_func_angle_ << " between " << res1_angle_ << atom1_angle_ << " and " << res_center_ << atom_center_ << " and " << res2_angle_ << atom2_angle_ << std::endl << "and then the function " << cst_func_angle2_ << " between " << res1_angle_ << atom1_angle2_ << " and " << res_center2_ << atom_center2_ << " and " << res2_angle_ << atom2_angle2_ << std::endl;
	} else {

		std::cout << "setting angles you assigned" << std::endl;
		std::cout <<res1_angle_ << atom1_angle_ << res_center_<< atom_center_ << res2_angle_ << atom2_angle_ << cst_func_angle_ << std::endl;
		angle1->set(res_center_,atom_center_,res1_angle_,atom1_angle_,res2_angle_,atom2_angle_,cst_func_angle_);
		angle1->apply(pose);
	}

	protocols::cyclic_peptide::CreateTorsionConstraintOP torsion(new protocols::cyclic_peptide::CreateTorsionConstraint);
	torsion->set(res1_torsion_,atom1_torsion_,res2_torsion_,atom2_torsion_,res3_torsion_,atom3_torsion_,res4_torsion_,atom4_torsion_,cst_func_torsion_);
	torsion->apply(pose);

	std::cout << "we are setting the torsion function " << cst_func_torsion_ << " between " << res1_torsion_ << atom1_torsion_ << " and " << res2_torsion_ << atom2_torsion_ << " and " << res3_torsion_ << atom3_torsion_ << " and " << res4_torsion_ << atom4_torsion_ << std::endl;

	if ( !bond_assigned_ ) {
		protocols::cyclic_peptide::DeclareBondOP bond_close(new protocols::cyclic_peptide::DeclareBond);
		bond_close->set(res1_[1],atom1_[1],res2_[1],atom2_[1],add_termini_[1],false,0,0,rebuild_fold_tree_[1]);
		bond_close->apply(pose);
	} else {
		for ( core::Size i=1; i<=counter_; i++ ) {
			protocols::cyclic_peptide::DeclareBondOP bond_close(new protocols::cyclic_peptide::DeclareBond);
			bond_close->set(res1_[i],atom1_[i],res2_[i],atom2_[i],add_termini_[i],false,0,0,rebuild_fold_tree_[i]);
			bond_close->apply(pose);
		}
	}


}

////////////////////////////////////////////////////////////////////////////////


/// @brief Returns the name of this mover ("PeptideCyclizeMover").

// XRW TEMP std::string
// XRW TEMP PeptideCyclizeMover::get_name() const {
// XRW TEMP  return "PeptideCyclizeMover";
// XRW TEMP }

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
PeptideCyclizeMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
)
{
	std::string const name( tag->getOption<std::string>( "name" ));

	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );//enabling addition of multiple user defined tags
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
	counter_=0;
	for ( tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() == "Bond" ) {
			//printf ("you have assigned a Bond to be formed\n");
			set_bond ((*tag_it)->getOption<core::Size> ("res1"),(*tag_it)->getOption<std::string> ("atom1"),(*tag_it)->getOption<core::Size> ("res2"),(*tag_it)->getOption<std::string> ("atom2"), (*tag_it)->getOption< bool >( "add_termini" ),(*tag_it)->getOption< bool >( "rebuild_fold_tree", false ));
		}
		if ( (*tag_it)->getName() == "Distance" ) {
			//printf ("getting distance constraints from user\n");
			set_distance ((*tag_it)->getOption<core::Size> ("res1"),(*tag_it)->getOption<std::string> ("atom1"),(*tag_it)->getOption<core::Size> ("res2"),(*tag_it)->getOption<std::string> ("atom2"),(*tag_it)->getOption<std::string> ("cst_func",""));
		}
		if ( (*tag_it)->getName() == "Angle" ) {
			//printf ("getting angle constraints from user \n");
			set_angle ((*tag_it)->getOption<core::Size> ("res_center"),(*tag_it)->getOption<std::string> ("atom_center"),(*tag_it)->getOption<core::Size> ("res1"),(*tag_it)->getOption<std::string> ("atom1"),(*tag_it)->getOption<core::Size> ("res2"),(*tag_it)->getOption<std::string> ("atom2"),(*tag_it)->getOption<std::string> ("cst_func",""));
		}
		if ( (*tag_it)->getName() == "Torsion" ) {
			std::cout << "getting torsion constraints from user" << std::endl;
			set_torsion ((*tag_it)->getOption<core::Size> ("res1"),(*tag_it)->getOption<std::string> ("atom1"),(*tag_it)->getOption<core::Size> ("res2"),(*tag_it)->getOption<std::string> ("atom2"),(*tag_it)->getOption<core::Size> ("res3"),(*tag_it)->getOption<std::string> ("atom3"),(*tag_it)->getOption<core::Size> ("res4"),(*tag_it)->getOption<std::string> ("atom4"),(*tag_it)->getOption<std::string> ("cst_func"));
		}
	}

	if ( tag->hasOption("residue_selector") ) {
		//printf ("residue selector being read in\n");
		set_selector( protocols::rosetta_scripts::parse_residue_selector( tag, data ) );
	}


}

//Private functions:

///@brief getters
//these are getters. They basically assign all the private variables.

void PeptideCyclizeMover::get_values ( ) {
	//generating a residue type and getting the values for ideal bond length and angles from it
	core::chemical::ResidueTypeSetCOP standard_residues = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	core::chemical ::ResidueType new_rsd =standard_residues->name_map("ALA");

	core::chemical::AtomICoor icoor_lower ( new_rsd.lower_connect().icoor());
	core::chemical::AtomICoor icoor_upper ( new_rsd.upper_connect().icoor());

	angle2_=numeric::constants::d::pi-icoor_lower.theta(); //the angle for the peptide is 180-angle from lower_connect. Basic Math :P
	angle1_=numeric::constants::d::pi-icoor_upper.theta();
	distance_=icoor_lower.d();

}
//This function sets the values if they are not already assigned.
void PeptideCyclizeMover::get_all(core::select::residue_selector::ResidueSubset subset, core::pose::Pose const & pose){

	core::Size const & nres=pose.total_residue();

	//calling to find angles if they are not assigned.
	get_values();

	std::cout << "the bond is " << distance_ << " and the angles are " << angle1_ << " and " << angle2_ << std::endl;

	if ( !bond_assigned_ ) {
		core::Size starter=0;
		core::Size counter=0;
		for ( Size i=1; i<= nres; ++i ) {
			if ( subset[i] ) {
				counter=i;
				starter+=1;
				if ( starter==1 ) {
					res2_.push_back(counter);
					atom2_.push_back("N");
				}
			}
		}
		res1_.push_back(counter);
		atom1_.push_back("C");
		add_termini_.push_back(false);
		rebuild_fold_tree_.push_back(false);
		std::cout << "The last and first residue to close are "<<  res1_[1] << atom1_[1] << " & " << res2_[1] <<  atom2_[1] << std::endl;
	}

	if ( !distance_assigned_ ) {
		core::Size starter=0;
		core::Size counter=0;
		for ( Size i=1; i<= nres; ++i ) {
			if ( subset[i] ) {
				counter=i;
				starter+=1;
				if ( starter==1 ) {
					res2_dist_.push_back(counter);
					atom2_dist_.push_back("N");
				}
			}
		}
		res1_dist_.push_back(counter);
		atom1_dist_.push_back("C");
		std::string dist_func_;
		std::ostringstream sstm;//to add back the integer into string and use it as a function
		sstm <<"HARMONIC " << distance_ << " 0.01";
		dist_func_=sstm.str();
		cst_func_dist_.push_back(dist_func_);
		std::cout << "this is what I am setting for you for distance with function of " << dist_func_ << std::endl;
	}

	if ( !angle_assigned_ ) {
		core::Size starter=0;
		core::Size counter=0;
		for ( Size i=1; i<= nres; ++i ) {
			if ( subset[i] ) {
				counter=i;
				starter+=1;
				if ( starter==1 ) {
					res2_angle_.push_back(counter);
					atom2_angle_.push_back("N");
					res_center2_.push_back(counter);
					atom_center2_.push_back("N");
					atom2_angle2_.push_back("CA");
				}
			}
		}
		res1_angle_.push_back(counter);
		atom1_angle_.push_back("CA");
		res_center_.push_back(counter);
		atom_center_.push_back("C");
		atom1_angle2_.push_back("C");
		std::string angle1_func_;
		std::ostringstream sstm1;
		sstm1 <<"HARMONIC " << angle1_ << " 0.01";
		angle1_func_=sstm1.str();
		std::string angle2_func_;
		std::ostringstream sstm2;
		sstm2 <<"HARMONIC " << angle2_ << " 0.01";
		angle2_func_=sstm2.str();
		cst_func_angle_.push_back(angle1_func_);
		cst_func_angle2_.push_back(angle2_func_);
	}

	if ( !torsion_assigned_ ) {
		core::Size starter=0;
		core::Size counter=0;
		for ( Size i=1; i<= nres; ++i ) {
			if ( subset[i] ) {
				counter=i;
				starter+=1;
				if ( starter==1 ) {
					res3_torsion_.push_back(counter);
					atom3_torsion_.push_back("N");
					res4_torsion_.push_back(counter);
					atom4_torsion_.push_back("CA");
				}
			}
		}
		res1_torsion_.push_back(counter);
		atom1_torsion_.push_back("CA");
		res2_torsion_.push_back(counter);
		atom2_torsion_.push_back("C");
		cst_func_torsion_.push_back("CIRCULARHARMONIC 3.141592654 0.005");
	}

}

std::string PeptideCyclizeMover::get_name() const {
	return mover_name();
}

std::string PeptideCyclizeMover::mover_name() {
	return "PeptideCyclizeMover";
}

std::string subtag_for_cyclize( std::string const & foo ) {
	return "subtag_for_cyclize_" + foo + "_complex_type";
}

void PeptideCyclizeMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist, bond_attributes, distance_attributes, angle_attributes, torsion_attributes;

	XMLSchemaSimpleSubelementList ssl;

	bond_attributes + XMLSchemaAttribute( "res1", xsct_non_negative_integer, "Residue one" )
		+ XMLSchemaAttribute( "atom1", xs_string, "Atom one" )
		+ XMLSchemaAttribute( "res2", xsct_non_negative_integer, "Residue two" )
		+ XMLSchemaAttribute( "atom2", xs_string, "Atom two" )
		+ XMLSchemaAttribute( "add_termini", xsct_rosetta_bool, "Add terminal types where necessary" )
		+ XMLSchemaAttribute::attribute_w_default( "rebuild_fold_tree", xsct_rosetta_bool, "Rebuild the fold tree around this bond", "false" );

	distance_attributes + XMLSchemaAttribute( "res1", xsct_non_negative_integer, "Residue one" )
		+ XMLSchemaAttribute( "atom1", xs_string, "Atom one" )
		+ XMLSchemaAttribute( "res2", xsct_non_negative_integer, "Residue two" )
		+ XMLSchemaAttribute( "atom2", xs_string, "Atom two" )
		+ XMLSchemaAttribute( "cst_func", xs_string, "Function to use as a constraint" );

	angle_attributes + XMLSchemaAttribute( "res_center", xsct_non_negative_integer, "Central residue" )
		+ XMLSchemaAttribute( "atom_center", xs_string, "Central atom" )
		+ XMLSchemaAttribute( "res1", xsct_non_negative_integer, "Residue one" )
		+ XMLSchemaAttribute( "atom1", xs_string, "Atom one" )
		+ XMLSchemaAttribute( "res2", xsct_non_negative_integer, "Residue two" )
		+ XMLSchemaAttribute( "atom2", xs_string, "Atom two" )
		+ XMLSchemaAttribute( "cst_func", xs_string, "Function to use as a constraint" );

	torsion_attributes + XMLSchemaAttribute( "res1", xsct_non_negative_integer, "Residue one" )
		+ XMLSchemaAttribute( "atom1", xs_string, "Atom one" )
		+ XMLSchemaAttribute( "res2", xsct_non_negative_integer, "Residue two" )
		+ XMLSchemaAttribute( "atom2", xs_string, "Atom two" )
		+ XMLSchemaAttribute( "res3", xsct_non_negative_integer, "Residue three" )
		+ XMLSchemaAttribute( "atom3", xs_string, "Atom three" )
		+ XMLSchemaAttribute( "res4", xsct_non_negative_integer, "Residue four" )
		+ XMLSchemaAttribute( "atom4", xs_string, "Atom four" )
		+ XMLSchemaAttribute( "cst_func", xs_string, "Function to use as a constraint" );

	ssl.add_simple_subelement( "Bond", bond_attributes, "Tags describing a bond from the macrocycle"/*, 0 minoccurs*/ )
		.add_simple_subelement( "Distance", distance_attributes, "Tags describing a distance from the macrocycle"/*, 0 minoccurs*/ )
		.add_simple_subelement( "Angle", angle_attributes, "Tags describing a angle from the macrocycle"/*, 0 minoccurs*/ )
		.add_simple_subelement( "Torsion", torsion_attributes, "Tags describing a torsion from the macrocycle"/*, 0 minoccurs*/ )
		.complex_type_naming_func( & subtag_for_cyclize );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, ssl );
}

std::string PeptideCyclizeMoverCreator::keyname() const {
	return PeptideCyclizeMover::mover_name();
}

protocols::moves::MoverOP
PeptideCyclizeMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PeptideCyclizeMover );
}

void PeptideCyclizeMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PeptideCyclizeMover::provide_xml_schema( xsd );
}




} // moves
} // protocols
