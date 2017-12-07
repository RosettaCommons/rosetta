// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/MetalContactsConstraintGenerator.cc
/// @brief This constraint generator takes residue selectors for a residue/residues containing metal ion(s) and for residue(s) for which to set up contacts. It allows users to specify which base atoms will be used to define angles/dihedrals to constrain; ideal values for angles/dihedrals/distances; and an option to constrain to native values.
/// @author guffysl (guffy@email.unc.edu)

// Unit headers
#include <protocols/constraint_generator/MetalContactsConstraintGenerator.hh>
#include <protocols/constraint_generator/MetalContactsConstraintGeneratorCreator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>


//Protocols


//Core
#include <core/util/metalloproteins_util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/id/AtomID.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

static basic::Tracer TR( "protocols.constraint_generator.MetalContactsConstraintGenerator" );

namespace protocols {
namespace constraint_generator {

protocols::constraint_generator::ConstraintGeneratorOP
MetalContactsConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new MetalContactsConstraintGenerator );
}

std::string
MetalContactsConstraintGeneratorCreator::keyname() const
{
	return MetalContactsConstraintGenerator::class_name();
}

MetalContactsConstraintGenerator::MetalContactsConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( MetalContactsConstraintGenerator::class_name() )
{
}

MetalContactsConstraintGenerator::~MetalContactsConstraintGenerator()
{}

MetalContactsConstraintGenerator::MetalContactsConstraintGenerator( MetalContactsConstraintGenerator const & src ):
	protocols::constraint_generator::ConstraintGenerator( src )
{
	dist_cutoff_multiplier_ = src.dist_cutoff_multiplier_;
	use_ligand_selector_ = src.use_ligand_selector_;
	ligand_atom_name_ = src.ligand_atom_name_;
	ligand_selector_ = src.ligand_selector_;
	ligand_resnum_string_ = src.ligand_resnum_string_;
	use_contact_selector_ = src.use_contact_selector_;
	contact_selector_ = src.contact_selector_;
	contact_resnum_string_ = src.contact_resnum_string_;
	base_atom_name_ = src.base_atom_name_;
	base_base_atom_name_ = src.base_base_atom_name_;
	ideal_distance_ = src.ideal_distance_;
	ideal_angle_about_contact_ = src.ideal_angle_about_contact_;
	ideal_dihedral_about_contact_ = src.ideal_dihedral_about_contact_;
	ideal_angle_about_metal_ = src.ideal_angle_about_metal_;
	ideal_dihedral_about_metal_ = src.ideal_dihedral_about_metal_;
	ideal_dihedral_3_ = src.ideal_dihedral_3_;
	score_against_internal_contacts_ = src.score_against_internal_contacts_;
	constrain_to_closest_= src.constrain_to_closest_;
}


protocols::constraint_generator::ConstraintGeneratorOP
MetalContactsConstraintGenerator::clone() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new MetalContactsConstraintGenerator( *this ) );
}

std::string
MetalContactsConstraintGenerator::class_name()
{
	return "MetalContactsConstraintGenerator";
}

void
MetalContactsConstraintGenerator::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap)
{
	dist_cutoff_multiplier_ = tag->getOption< core::Real >( "dist_cutoff_multiplier", 1.0 );
	ligand_atom_name_ = tag->getOption< std::string >( "ligand_atom_name" );
	if ( tag->hasOption( "base_atom_name" ) ) {
		base_atom_name_ = tag->getOption< std::string >( "base_atom_name" );
	} else {
		base_atom_name_ = "";
	}
	if ( tag->hasOption( "base_base_atom_name" ) ) {
		base_base_atom_name_ = tag->getOption< std::string >( "base_base_atom_name" );
	} else {
		base_base_atom_name_ = "";
	}
	ideal_distance_ = tag->getOption< core::Real >( "ideal_distance", -1.0 ); //Negative value indicates to use native distance

	if ( tag->hasOption( "ideal_angle_about_contact" ) ) {
		utility::vector1< std::string > angle_strings = utility::string_split( tag->getOption< std::string >( "ideal_angle_about_contact" ), ',' );
		for ( std::string val: angle_strings ) {
			core::Real n(0);
			std::istringstream ss( val );
			ss >> n;
			ideal_angle_about_contact_.push_back( n );
		}
	}
	if ( tag->hasOption( "ideal_dihedral_about_contact" ) ) {
		utility::vector1< std::string > dihedral_strings = utility::string_split( tag->getOption< std::string >( "ideal_dihedral_about_contact" ), ',' );
		for ( std::string val: dihedral_strings ) {
			core::Real n(0);
			std::istringstream ss( val );
			ss >> n;
			ideal_dihedral_about_contact_.push_back( n );
		}
	}
	if ( tag->hasOption( "ideal_angle_about_metal" ) ) {
		utility::vector1< std::string > angle_strings = utility::string_split( tag->getOption< std::string >( "ideal_angle_about_metal" ), ',' );
		for ( std::string val: angle_strings ) {
			core::Real n(0);
			std::istringstream ss( val );
			ss >> n;
			ideal_angle_about_metal_.push_back( n );
		}
	}

	if ( tag->hasOption( "ideal_dihedral_about_metal" ) ) {
		utility::vector1< std::string > dihedral_strings = utility::string_split( tag->getOption< std::string >( "ideal_dihedral_about_metal" ), ',' );
		for ( std::string val: dihedral_strings ) {
			core::Real n(0);
			std::istringstream ss( val );
			ss >> n;
			ideal_dihedral_about_metal_.push_back( n );
		}
	}

	if ( tag->hasOption( "ideal_dihedral_3" ) ) {
		utility::vector1< std::string > dihedral_strings = utility::string_split( tag->getOption< std::string >( "ideal_dihedral_3" ), ',' );
		for ( std::string val: dihedral_strings ) {
			core::Real n(0);
			std::istringstream ss( val );
			ss >> n;
			ideal_dihedral_3_.push_back( n );
		}
	}

	score_against_internal_contacts_ = tag->getOption< bool >( "score_against_internal_contacts", false );
	constrain_to_closest_ = tag->getOption< bool >( "constrain_to_closest", true );
	//Res selectors and resnums
	//Ligand
	if ( tag->hasOption( "ligand_selector" ) && tag->hasOption( "ligand_resnum" ) ) {
		utility_exit_with_message( "Ligand must be specified with either a selector or a residue number, not both!" );
	}
	if ( !tag->hasOption( "ligand_selector" ) && !tag->hasOption( "ligand_resnum" ) ) {
		utility_exit_with_message( "You must specify a ligand residue!" );
	}
	if ( tag->hasOption( "ligand_selector" ) ) {
		use_ligand_selector_ = true;
		std::string selector_string = tag->getOption< std::string >( "ligand_selector" );
		ligand_selector_= datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_string );
	} else {
		ligand_resnum_string_ = tag->getOption< std::string >( "ligand_resnum" );
		use_ligand_selector_ = false;
	}

	//Contacts
	if ( tag->hasOption( "contact_selector" ) && tag->hasOption( "contact_resnums" ) ) {
		utility_exit_with_message( "Contact residues may be specified with a selector or a residue number, but not both!" );
	}
	if ( tag->hasOption( "contact_selector" ) ) {
		use_contact_selector_ = true;
		std::string selector_string = tag->getOption< std::string >( "contact_selector" );
		contact_selector_= datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_string );
	} else if ( tag->hasOption( "contact_resnums" ) ) {
		contact_resnum_string_ = tag->getOption< std::string >( "contact_resnums" );
		use_contact_selector_ = false;
	} else {
		use_contact_selector_ = false;
	}
}

core::scoring::constraints::ConstraintCOPs
MetalContactsConstraintGenerator::apply( core::pose::Pose const & pose ) const
{
	core::chemical::AtomTypeSetCAP atom_types_;
	if ( pose.is_fullatom() ) {
		atom_types_ = core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
	} else {
		atom_types_ = core::chemical::ChemicalManager::get_instance()->atom_type_set("centroid");
	}
	core::scoring::constraints::ConstraintCOPs csts;
	//Get ligand resnum
	core::Size ligand_resnum = get_ligand_resnum( pose );
	//Get contact resnums
	std::set< core::Size > contact_resnums; //If this is empty, we'll constrain all contacts
	get_contact_resnums( pose, contact_resnums );



	//TODO--allow selection of multiple ligands?







	std::set< core::id::AtomID > already_added_contacts;
	//For the specified ligand atom:
	core::id::AtomID metal_atom_id = core::id::AtomID( pose.residue( ligand_resnum ).atom_index( ligand_atom_name_), ligand_resnum );

	//Define funcs

	//Distance will have a harmonic func; all others will be AmbiguousConstraints with as many CircularHarmonicFunc as there are distances
	core::scoring::func::FuncOP distance_func;
	if ( ideal_distance_ > 0 ) {
		distance_func = core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc( ideal_distance_, 0.05 ) );
	}
	utility::vector1< core::scoring::func::FuncOP > angle_about_contact_func;
	for ( core::Real angle: ideal_angle_about_contact_ ) {
		angle = numeric::conversions::radians< core::Real >( angle );
		TR.Debug << "Ideal angle: " << angle << std::endl;
		angle_about_contact_func.push_back( core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( angle, 0.1 ) ) );
	}
	utility::vector1< core::scoring::func::FuncOP > dihedral_about_contact_func;
	for ( core::Real dihedral: ideal_dihedral_about_contact_ ) {
		dihedral = numeric::conversions::radians< core::Real >( dihedral );
		TR.Debug << "Ideal dihedral: " << dihedral << std::endl;
		dihedral_about_contact_func.push_back( core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( dihedral, 0.2 ) ) );
	}
	utility::vector1< core::scoring::func::FuncOP > angle_about_metal_func;
	for ( core::Real angle: ideal_angle_about_metal_ ) {
		angle = numeric::conversions::radians< core::Real >( angle );
		TR.Debug << "Ideal angle2: " << angle << std::endl;
		angle_about_metal_func.push_back( core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( angle, 0.1 ) ) );
	}
	utility::vector1< core::scoring::func::FuncOP > dihedral_about_metal_func;
	for ( core::Real dihedral: ideal_dihedral_about_metal_ ) {
		dihedral = numeric::conversions::radians< core::Real >( dihedral );
		TR.Debug << "Ideal dihedral2: " << dihedral << std::endl;
		dihedral_about_metal_func.push_back( core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( dihedral, 0.2 ) ) );
	}

	utility::vector1< core::scoring::func::FuncOP > dihedral_3_func;
	for ( core::Real dihedral: ideal_dihedral_3_ ) {
		dihedral = numeric::conversions::radians< core::Real >( dihedral );
		TR.Debug << "Ideal dihedral3: " << dihedral << std::endl;
		dihedral_3_func.push_back( core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( dihedral, 0.2 ) ) );
	}
	already_added_contacts.clear();
	//Find out what internal contact atoms are and add to already_added_contacts
	if ( score_against_internal_contacts_ ) {
		//These are found differently, by looking at the internal connectivity of the ligand as specified in the params file
		utility::vector1< core::Size > bonded_atoms = pose.residue( ligand_resnum ).bonded_neighbor( pose.residue( ligand_resnum ).atom_index( ligand_atom_name_ ) );
		//Add all these to already_added_contacts so we'll have constraints against them
		for ( core::Size atno: bonded_atoms ) {
			if ( ( *atom_types_.lock() )[ pose.residue( ligand_resnum ).atom( atno ).type() ].atom_type_name() == "VIRT"  ) {
				continue;
			}
			if ( atno == pose.residue( ligand_resnum ).atom_index( ligand_atom_name_ ) ) {
				continue;
			}
			already_added_contacts.insert( core::id::AtomID( atno, ligand_resnum ) );
		}
	}
	utility::vector1< core::id::AtomID > metal_contact_ids = core::util::find_metalbinding_atoms_helper( pose, metal_atom_id, dist_cutoff_multiplier_ );



	TR.Debug << "Outputting initial measures for contacts" << std::endl;
	core::id::AtomID base_atom;
	core::id::AtomID base_base_atom;
	core::id::AtomID base_2;
	core::id::AtomID base_base_2;
	for ( core::id::AtomID contact: metal_contact_ids ) {
		TR.Debug << "CONTACT\tDISTANCE\tANGLE\tDIHEDRAL" << std::endl;
		if ( base_atom_name_ != "" ) {
			base_atom = core::id::AtomID( pose.residue( contact.rsd() ).atom_index( base_atom_name_ ), contact.rsd() );
		} else {
			base_atom = core::id::AtomID( pose.residue( contact.rsd() ).atom_base( contact.atomno() ), contact.rsd() );
		}
		if ( base_base_atom_name_ != "" ) {
			base_base_atom = core::id::AtomID( pose.residue( contact.rsd() ).atom_index( base_base_atom_name_ ), contact.rsd() );
		} else {
			base_base_atom = core::id::AtomID( pose.residue( contact.rsd() ).atom_base( base_atom.atomno() ), contact.rsd() );
		}
		numeric::xyzVector< core::Real > base_base_coords = pose.residue( base_base_atom.rsd() ).atom( base_base_atom.atomno() ).xyz();
		numeric::xyzVector< core::Real > base_coords = pose.residue( base_atom.rsd() ).atom( base_atom.atomno() ).xyz();
		numeric::xyzVector< core::Real > contact_coords = pose.residue( contact.rsd() ).atom( contact.atomno() ).xyz();
		numeric::xyzVector< core::Real > metal_coords = pose.residue( metal_atom_id.rsd() ).atom( metal_atom_id.atomno() ).xyz();

		core::Real distance = contact_coords.distance( metal_coords );
		core::Real angle = numeric::angle_radians( base_coords, contact_coords, metal_coords );
		core::Real dihedral = numeric::dihedral_radians( base_base_coords, base_coords, contact_coords, metal_coords );
		TR << contact.rsd() << "\t" << distance << "\t" << angle << "\t" << dihedral << std::endl;

		TR.Debug << "CONTACT 1\tCONTACT 2\tANGLE\tDIHEDRAL 1\tDIHEDRAL 2" << std::endl;
		for ( core::id::AtomID contact2: metal_contact_ids ) {
			if ( base_atom_name_ != "" ) {
				base_2 = core::id::AtomID( pose.residue( contact2.rsd() ).atom_index( base_atom_name_ ), contact2.rsd() );
			} else {
				base_2 = core::id::AtomID( pose.residue( contact2.rsd() ).atom_base( contact2.atomno() ), contact2.rsd() );
			}
			if ( base_base_atom_name_ != "" ) {
				base_base_2 = core::id::AtomID( pose.residue( contact2.rsd() ).atom_index( base_base_atom_name_ ), contact2.rsd() );
			} else {
				base_base_2 = core::id::AtomID( pose.residue( contact2.rsd() ).atom_base( base_2.atomno() ), contact2.rsd() );
			}
			numeric::xyzVector< core::Real > base_2_coords = pose.residue( base_2.rsd() ).atom( base_2.atomno() ).xyz();
			numeric::xyzVector< core::Real > contact_2_coords = pose.residue( contact2.rsd() ).atom( contact2.atomno() ).xyz();
			core::Real angle2 = numeric::angle_radians( contact_coords, metal_coords, contact_2_coords );
			core::Real dihedral2 = numeric::dihedral_radians( base_coords, contact_coords, metal_coords, contact_2_coords );
			core::Real dihedral3 = numeric::dihedral_radians( contact_coords, metal_coords, contact_2_coords, base_2_coords );
			TR << contact.rsd() << "\t" << contact2.rsd() << "\t" << angle2 << "\t" << dihedral2 << "\t" << dihedral3 << std::endl;
		}
	}



	//TEMP
	TR << "ATOMS IDENTIFIED FOR CONSTRAINTS:" << std::endl;
	TR << "Metal: Residue " << metal_atom_id.rsd() << " atom " << metal_atom_id.atomno() << std::endl;
	for ( core::id::AtomID id: metal_contact_ids ) {
		TR << "Residue " << id.rsd() << " atom " << id.atomno() << std::endl;
	}

	//For each contact:
	for ( core::id::AtomID contact: metal_contact_ids ) {

		if ( ( *atom_types_.lock() )[ pose.residue( contact.rsd() ).atom( contact.atomno() ).type() ].atom_type_name() == "VIRT"  ) {
			continue;
		}


		//Skip any contacts that aren't included in contact_resnums (unless it's empty)
		if ( contact_resnums.size() != 0 && contact_resnums.find( contact.rsd() ) == contact_resnums.end() ) {
			TR << "Residue " << contact.rsd() << " not included in specified contact residues. Skipping contact." << std::endl;
			continue;
		}
		TR << "Begin constraint for residue " << contact.rsd() << " atom " << contact.atomno() << std::endl;
		//Get the base and base_base atoms
		core::id::AtomID base_atom;
		core::id::AtomID base_base_atom;
		bool valid_base = true;
		bool valid_base_base = true;
		if ( base_atom_name_ != "" ) {
			base_atom = core::id::AtomID( pose.residue( contact.rsd() ).atom_index( base_atom_name_ ), contact.rsd() );
		} else {
			base_atom = core::id::AtomID( pose.residue( contact.rsd() ).atom_base( contact.atomno() ), contact.rsd() );
			if ( base_atom.atomno() == contact.atomno() ) {
				valid_base = false;
			}
		}
		if ( base_base_atom_name_ != "" ) {
			base_base_atom = core::id::AtomID( pose.residue( contact.rsd() ).atom_index( base_base_atom_name_ ), contact.rsd() );
		} else {
			base_base_atom = core::id::AtomID( pose.residue( contact.rsd() ).atom_base( base_atom.atomno() ), contact.rsd() );
			if ( base_base_atom.atomno() == contact.atomno() || base_base_atom.atomno() == base_atom.atomno() ) {
				valid_base_base = false;
			}
		}


		//Define distance func if using current distances
		if ( ideal_distance_ <= 0 ) {
			numeric::xyzVector< core::Real > dist_at1_coords = pose.residue( contact.rsd() ).atom( contact.atomno() ).xyz();
			numeric::xyzVector< core::Real > dist_at2_coords = pose.residue( metal_atom_id.rsd() ).atom( metal_atom_id.atomno() ).xyz();
			core::Real distance = dist_at1_coords.distance( dist_at2_coords );
			TR << "Constraining to distance " << distance << std::endl;
			distance_func = core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc( distance, 0.05 ) );
		}
		//Distance constraint
		core::scoring::constraints::ConstraintCOP distance_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AtomPairConstraint( contact, metal_atom_id, distance_func, core::scoring::metalbinding_constraint ) );
		csts.push_back( distance_cst );


		//Angle constraint about contact
		core::scoring::constraints::ConstraintCOPs contact_angle_csts;
		numeric::xyzVector< core::Real > angle_at1_coords = pose.residue( base_atom.rsd() ).atom( base_atom.atomno() ).xyz();
		numeric::xyzVector< core::Real > angle_at2_coords = pose.residue( contact.rsd() ).atom( contact.atomno() ).xyz();
		numeric::xyzVector< core::Real > angle_at3_coords = pose.residue( metal_atom_id.rsd() ).atom( metal_atom_id.atomno() ).xyz();
		core::Real current_angle = numeric::angle_radians( angle_at1_coords, angle_at2_coords, angle_at3_coords );
		//If using current values, define here
		if ( angle_about_contact_func.size() == 0 && valid_base ) {
			//Define func
			TR << "Constraining to angle (in radians ) " << current_angle << std::endl;
			core::scoring::func::FuncOP anglefunc = core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( current_angle, .1 ) ); //~5 degrees
			//Define constraint
			core::scoring::constraints::ConstraintCOP angle_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AngleConstraint( base_atom, contact, metal_atom_id, anglefunc, core::scoring::metalbinding_constraint ) );
			//Add to set
			contact_angle_csts.push_back( angle_cst );
		} else if ( constrain_to_closest_ ) {
			//We will use the constraint corresponding to the ideal value closest to the current value
			utility::vector1< core::Real > diff_vector;
			for ( core::scoring::func::FuncOP anglefunc: angle_about_contact_func ) {
				diff_vector.push_back( std::abs( std::dynamic_pointer_cast< core::scoring::func::CircularHarmonicFunc>( anglefunc )->x0() - current_angle ) );
			}
			core::scoring::func::FuncOP anglefunc = angle_about_contact_func.at( utility::arg_min( diff_vector ) );
			//Define constraint
			core::scoring::constraints::ConstraintCOP angle_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AngleConstraint( base_atom, contact, metal_atom_id, anglefunc, core::scoring::metalbinding_constraint ) );
			//Add to set
			contact_angle_csts.push_back( angle_cst );
		} else {
			for ( core::scoring::func::FuncOP anglefunc: angle_about_contact_func ) {
				if ( !valid_base ) break;
				core::scoring::constraints::ConstraintCOP angle_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AngleConstraint( base_atom, contact, metal_atom_id, anglefunc, core::scoring::metalbinding_constraint ) );
				contact_angle_csts.push_back( angle_cst );
			}
		}
		if ( valid_base ) {
			core::scoring::constraints::ConstraintCOP contact_angle_ambig_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AmbiguousConstraint( contact_angle_csts ) );
			csts.push_back( contact_angle_ambig_cst );
		}






		//Dihedral constraint about contact
		if ( valid_base && valid_base_base ) {
			core::scoring::constraints::ConstraintCOPs contact_dihedral_csts;
			numeric::xyzVector< core::Real > di_at1_coords = pose.residue( base_base_atom.rsd() ).atom( base_base_atom.atomno() ).xyz();
			numeric::xyzVector< core::Real > di_at2_coords = pose.residue( base_atom.rsd() ).atom( base_atom.atomno() ).xyz();
			numeric::xyzVector< core::Real > di_at3_coords = pose.residue( contact.rsd() ).atom( contact.atomno() ).xyz();
			numeric::xyzVector< core::Real > di_at4_coords = pose.residue( metal_atom_id.rsd() ).atom( metal_atom_id.atomno() ).xyz();
			core::Real current_dihedral = numeric::dihedral_radians( di_at1_coords, di_at2_coords, di_at3_coords, di_at4_coords );
			//If using current values
			if ( dihedral_about_contact_func.size() == 0 ) {
				//Define func
				TR << "Constraining to dihedral (in radians ) " << current_dihedral << std::endl;
				core::scoring::func::FuncOP dihedralfunc = core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( current_dihedral, .2 ) ); //~10 degrees
				//Define constraint
				core::scoring::constraints::ConstraintCOP dihedral_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::DihedralConstraint( base_base_atom, base_atom, contact, metal_atom_id, dihedralfunc, core::scoring::metalbinding_constraint ) );
				//Add to set
				contact_dihedral_csts.push_back( dihedral_cst );
			} else if ( constrain_to_closest_ ) {
				//We will use the constraint corresponding to the ideal value closest to the current value
				utility::vector1< core::Real > diff_vector;
				for ( core::scoring::func::FuncOP dihedralfunc: dihedral_about_contact_func ) {
					diff_vector.push_back( std::abs( std::dynamic_pointer_cast< core::scoring::func::CircularHarmonicFunc>( dihedralfunc )->x0() - current_dihedral ) );
				}
				core::scoring::func::FuncOP dihedralfunc = dihedral_about_contact_func.at( utility::arg_min( diff_vector ) );
				//Define constraint
				core::scoring::constraints::ConstraintCOP dihedral_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::DihedralConstraint( base_base_atom, base_atom, contact, metal_atom_id, dihedralfunc, core::scoring::metalbinding_constraint ) );
				//Add to set
				contact_dihedral_csts.push_back( dihedral_cst );
			} else {
				for ( core::scoring::func::FuncOP dihedralfunc: dihedral_about_contact_func ) {
					core::scoring::constraints::ConstraintCOP dihedral_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::DihedralConstraint( base_base_atom, base_atom, contact, metal_atom_id, dihedralfunc, core::scoring::metalbinding_constraint ) );
					contact_dihedral_csts.push_back( dihedral_cst );
				}
			}
			core::scoring::constraints::ConstraintCOP contact_dihedral_ambig_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AmbiguousConstraint( contact_dihedral_csts ) );
			csts.push_back( contact_dihedral_ambig_cst );
		}





		//Angle constraints about metal--make an AmbiguousConstraint & use one for each func
		for ( core::id::AtomID other_contact: already_added_contacts ) {
			TR << "Begin other constraint for residue " << other_contact.rsd() << " atom " << other_contact.atomno() << std::endl;
			core::scoring::constraints::ConstraintCOPs angle_csts;
			numeric::xyzVector< core::Real > an_at1_coords = pose.residue( contact.rsd() ).atom( contact.atomno() ).xyz();
			numeric::xyzVector< core::Real > an_at2_coords = pose.residue( metal_atom_id.rsd() ).atom( metal_atom_id.atomno() ).xyz();
			numeric::xyzVector< core::Real > an_at3_coords = pose.residue( other_contact.rsd() ).atom( other_contact.atomno() ).xyz();
			core::Real current_angle2 = numeric::angle_radians( an_at1_coords, an_at2_coords, an_at3_coords );

			//If using current values
			if ( angle_about_metal_func.size() == 0 ) {
				//Define func
				TR << "Constraining to angle (in radians ) " << current_angle2 << std::endl;
				core::scoring::func::FuncOP anglefunc = core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( current_angle2, .1 ) );
				//Define constraint
				core::scoring::constraints::ConstraintCOP angle_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AngleConstraint( contact, metal_atom_id, other_contact, anglefunc, core::scoring::metalbinding_constraint ) );
				//Add to set
				angle_csts.push_back( angle_cst );
			} else if ( constrain_to_closest_ ) {
				utility::vector1< core::Real > diff_vector;
				for ( core::scoring::func::FuncOP anglefunc: angle_about_metal_func ) {
					diff_vector.push_back( std::abs( std::dynamic_pointer_cast< core::scoring::func::CircularHarmonicFunc>( anglefunc )->x0() - current_angle2 ) );
				}
				core::scoring::func::FuncOP anglefunc = angle_about_metal_func.at( utility::arg_min( diff_vector ) );
				core::scoring::constraints::ConstraintCOP angle_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AngleConstraint( contact, metal_atom_id, other_contact, anglefunc, core::scoring::metalbinding_constraint ) );
				//Add to set
				angle_csts.push_back( angle_cst );
			} else {
				for ( core::scoring::func::FuncOP anglefunc: angle_about_metal_func ) {
					core::scoring::constraints::ConstraintCOP angle_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AngleConstraint( contact, metal_atom_id, other_contact, anglefunc, core::scoring::metalbinding_constraint ) );
					angle_csts.push_back( angle_cst );
				}
			}
			core::scoring::constraints::ConstraintCOP angle_ambig_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AmbiguousConstraint( angle_csts ) );
			csts.push_back( angle_ambig_cst );















			//Dihedral constraints about metal--sticking in same for loop to avoid looping twice
			if ( valid_base ) {
				core::scoring::constraints::ConstraintCOPs dihedral_csts;
				numeric::xyzVector< core::Real > d2_at1_coords = pose.residue( base_atom.rsd() ).atom( base_atom.atomno() ).xyz();
				numeric::xyzVector< core::Real > d2_at2_coords = pose.residue( contact.rsd() ).atom( contact.atomno() ).xyz();
				numeric::xyzVector< core::Real > d2_at3_coords = pose.residue( metal_atom_id.rsd() ).atom( metal_atom_id.atomno() ).xyz();
				numeric::xyzVector< core::Real > d2_at4_coords = pose.residue( other_contact.rsd() ).atom( other_contact.atomno() ).xyz();

				core::Real current_dihedral2 = numeric::dihedral_radians( d2_at1_coords, d2_at2_coords, d2_at3_coords, d2_at4_coords );

				//If using current values
				if ( dihedral_about_metal_func.size() == 0 ) {
					//Define func
					TR << "Constraining to dihedral (in radians ) " << current_dihedral2 << std::endl;
					core::scoring::func::FuncOP dihedralfunc = core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( current_dihedral2, 0.2 ) );
					//Define constraint
					core::scoring::constraints::ConstraintCOP dihedral_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::DihedralConstraint( base_atom, contact, metal_atom_id, other_contact, dihedralfunc, core::scoring::metalbinding_constraint ) );
					//Add to set
					dihedral_csts.push_back( dihedral_cst );
				} else if ( constrain_to_closest_ ) {
					utility::vector1< core::Real > diff_vector;
					for ( core::scoring::func::FuncOP dihedralfunc: dihedral_about_metal_func ) {
						diff_vector.push_back( std::abs( std::dynamic_pointer_cast< core::scoring::func::CircularHarmonicFunc>( dihedralfunc )->x0() - current_dihedral2 ) );
					}
					core::scoring::func::FuncOP dihedralfunc = dihedral_about_metal_func.at( utility::arg_min( diff_vector ) );
					core::scoring::constraints::ConstraintCOP dihedral_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::DihedralConstraint( base_atom, contact, metal_atom_id, other_contact, dihedralfunc, core::scoring::metalbinding_constraint ) );
					//Add to set
					dihedral_csts.push_back( dihedral_cst );
				} else {
					for ( core::scoring::func::FuncOP dihedralfunc: dihedral_about_metal_func ) {
						core::scoring::constraints::ConstraintCOP dihedral_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::DihedralConstraint( base_atom, contact, metal_atom_id, other_contact, dihedralfunc, core::scoring::metalbinding_constraint ) );
						dihedral_csts.push_back( dihedral_cst );
					}
				}
				core::scoring::constraints::ConstraintCOP dihedral_ambig_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AmbiguousConstraint( dihedral_csts ) );
				csts.push_back( dihedral_ambig_cst );
			}
			//Dihedral 3 constraints
			core::scoring::constraints::ConstraintCOPs dihedral_3_csts;
			//Get AtomID for base of other_contact
			core::id::AtomID other_base;
			if ( base_atom_name_ != "" ) {
				TR << "WARNING: Using Rosetta base for other contact base atom." << std::endl;
			}
			other_base = core::id::AtomID( pose.residue( other_contact.rsd() ).atom_base( other_contact.atomno() ), other_contact.rsd() );
			if ( other_base.atomno() == other_contact.atomno() ) {
				continue;
			}

			numeric::xyzVector< core::Real > d3_at1_coords = pose.residue( contact.rsd() ).atom( contact.atomno() ).xyz();
			numeric::xyzVector< core::Real > d3_at2_coords = pose.residue( metal_atom_id.rsd() ).atom( metal_atom_id.atomno() ).xyz();
			numeric::xyzVector< core::Real > d3_at3_coords = pose.residue( other_contact.rsd() ).atom( other_contact.atomno() ).xyz();
			numeric::xyzVector< core::Real > d3_at4_coords = pose.residue( other_base.rsd() ).atom( other_base.atomno() ).xyz();
			core::Real current_dihedral3 = numeric::dihedral_radians( d3_at1_coords, d3_at2_coords, d3_at3_coords, d3_at4_coords );

			//If using current values
			if ( dihedral_3_func.size() == 0 ) {
				//Define func
				TR << "Constraining to dihedral (in radians ) " << current_dihedral3 << std::endl;
				core::scoring::func::FuncOP dihedral3func = core::scoring::func::FuncOP( new core::scoring::func::CircularHarmonicFunc( current_dihedral3, 0.2 ) );
				//Define constraint
				core::scoring::constraints::ConstraintCOP dihedral_3_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::DihedralConstraint( contact, metal_atom_id, other_contact, other_base, dihedral3func, core::scoring::metalbinding_constraint ) );
				//Add to set
				dihedral_3_csts.push_back( dihedral_3_cst );
			} else if ( constrain_to_closest_ ) {
				utility::vector1< core::Real > diff_vector;
				for ( core::scoring::func::FuncOP dihedral3func: dihedral_3_func ) {
					diff_vector.push_back( std::abs( std::dynamic_pointer_cast< core::scoring::func::CircularHarmonicFunc>( dihedral3func )->x0() - current_dihedral3 ) );
				}
				core::scoring::func::FuncOP dihedral3func = dihedral_3_func.at( utility::arg_min( diff_vector ) );
				core::scoring::constraints::ConstraintCOP dihedral_3_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::DihedralConstraint( contact, metal_atom_id, other_contact, other_base, dihedral3func, core::scoring::metalbinding_constraint ) );
				//Add to set
				dihedral_3_csts.push_back( dihedral_3_cst );
			} else {
				for ( core::scoring::func::FuncOP dihedral3func: dihedral_3_func ) {
					core::scoring::constraints::ConstraintCOP dihedral_3_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::DihedralConstraint( contact, metal_atom_id, other_contact, other_base, dihedral3func, core::scoring::metalbinding_constraint ) );
					dihedral_3_csts.push_back( dihedral_3_cst );
				}
			}
			core::scoring::constraints::ConstraintCOP dihedral_3_ambig_cst = core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AmbiguousConstraint( dihedral_3_csts ) );
			csts.push_back( dihedral_3_ambig_cst );
		}

		//Add to ConstraintCOPs
		//Add to already_added_contacts
		TR << "Done with this atom" << std::endl;
		already_added_contacts.insert( contact );
	}
	return csts;
}

void
MetalContactsConstraintGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	MetalContactsConstraintGenerator::provide_xml_schema( xsd );
}

void
MetalContactsConstraintGenerator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "dist_cutoff_multiplier", xsct_real, "Multiply van der Waals radii of metal and contact atom by this value during contact detection.", "1.0" )
		+ XMLSchemaAttribute::required_attribute( "ligand_atom_name", xs_string, "Name of ligand metal atom you want to constrain." )
		+ XMLSchemaAttribute( "ligand_selector", xs_string, "Residue selector specifying which metal-containing ligand to constrain" )
		+ XMLSchemaAttribute( "contact_selector", xs_string, "Residue selector restricting which residues should be considered as potential ligand contacts." )
		+ XMLSchemaAttribute( "ligand_resnum", xsct_positive_integer, "Residue number for ligand to be constrained" )
		+ XMLSchemaAttribute( "contact_resnums", xsct_refpose_enabled_residue_number_cslist, "Residue numbers for residues that could be considered as contacts. If neither this option nor a residue selector is specified, then all residues are considered." )
		+ XMLSchemaAttribute( "base_atom_name", xs_string, "Name of atom to use as base of contact atoms for angles/dihedrals. Defaults to residue's base atom for contact atom." )
		+ XMLSchemaAttribute( "base_base_atom_name", xs_string, "Name of atom to use as base of base of contact atoms for angles/dihedrals. Defaults to residue's base atom." )
		+ XMLSchemaAttribute( "ideal_distance", xsct_real, "Ideal distance between constrained metal and contact atom. Defaults to current distance." )
		+ XMLSchemaAttribute( "ideal_angle_about_contact", xsct_real_cslist, "Comma-separated list of possible optimal angles, base-contact-metal. Defaults to current angle." )
		+ XMLSchemaAttribute( "ideal_dihedral_about_contact", xsct_real_cslist, "Comma-separated list of possible optimal dihedrals, base_base-base-contact-metal. Defaults to current dihedral." )
		+ XMLSchemaAttribute( "ideal_angle_about_metal", xsct_real_cslist, "Comma-separated list of possible optimal angles, contact-metal-other_contact. Defaults to current angle." )
		+ XMLSchemaAttribute( "ideal_dihedral_about_metal", xsct_real_cslist, "Comma-separated list of possible optimal dihedrals, base-contact-metal-other_contact. Defaults to current dihedral." )
		+ XMLSchemaAttribute( "ideal_dihedral_3", xsct_real_cslist, "Comma-separated list of possible optimal dihedrals, contact-metal-other_contact-other_base. Defaults to current dihedral." )
		+ XMLSchemaAttribute::attribute_w_default( "score_against_internal_contacts", xsct_rosetta_bool, "Should we score angles and dihedrals vs other atoms in the ligand?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "constrain_to_closest", xsct_rosetta_bool, "When multiple ideal values are provided for an angle/dihedral, constrain to the one that is closest to the current value", "true" );

	ConstraintGeneratorFactory::xsd_constraint_generator_type_definition_w_attributes(
		xsd,
		class_name(),
		"Generates distance, angle, and dihedral constraints for the specified metal atom in the selected ligand residue",
		attlist );
}


//GETTERS
core::Real
MetalContactsConstraintGenerator::get_dist_cutoff_multiplier() const{
	return dist_cutoff_multiplier_;
}

bool
MetalContactsConstraintGenerator::get_constrain_to_closest() const{
	return constrain_to_closest_;
}

std::string
MetalContactsConstraintGenerator::get_ligand_atom_name() const{
	return ligand_atom_name_;
}

bool
MetalContactsConstraintGenerator::get_use_ligand_selector() const{
	return use_ligand_selector_;
}

core::select::residue_selector::ResidueSelectorCOP
MetalContactsConstraintGenerator::get_ligand_selector() const{
	return ligand_selector_;
}

std::string
MetalContactsConstraintGenerator::get_ligand_resnum_string() const{
	return ligand_resnum_string_;
}

bool
MetalContactsConstraintGenerator::get_use_contact_selector() const{
	return use_contact_selector_;
}

core::select::residue_selector::ResidueSelectorCOP
MetalContactsConstraintGenerator::get_contact_selector() const{
	return contact_selector_;
}

std::string
MetalContactsConstraintGenerator::get_contact_resnum_string() const{
	return contact_resnum_string_;
}

std::string
MetalContactsConstraintGenerator::get_base_atom_name() const{
	return base_atom_name_;
}

std::string
MetalContactsConstraintGenerator::get_base_base_atom_name() const{
	return base_base_atom_name_;
}

core::Real
MetalContactsConstraintGenerator::get_ideal_distance() const{
	return ideal_distance_;
}

utility::vector1< core::Real >
MetalContactsConstraintGenerator::get_ideal_angle_about_contact() const{
	return ideal_angle_about_contact_;
}

utility::vector1< core::Real >
MetalContactsConstraintGenerator::get_ideal_dihedral_about_contact() const{
	return ideal_dihedral_about_contact_;
}

utility::vector1< core::Real >
MetalContactsConstraintGenerator::get_ideal_angle_about_metal() const{
	return ideal_angle_about_metal_;
}

utility::vector1< core::Real >
MetalContactsConstraintGenerator::get_ideal_dihedral_about_metal() const{
	return ideal_dihedral_about_metal_;
}

utility::vector1< core::Real >
MetalContactsConstraintGenerator::get_ideal_dihedral_3() const{
	return ideal_dihedral_3_;
}

bool
MetalContactsConstraintGenerator::get_score_against_internal_contacts() const{
	return score_against_internal_contacts_;
}

//SETTERS
void
MetalContactsConstraintGenerator::set_dist_cutoff_multiplier( core::Real setting){
	dist_cutoff_multiplier_ = setting;
}

void
MetalContactsConstraintGenerator::set_constrain_to_closest( bool setting ){
	constrain_to_closest_ = setting;
}

void
MetalContactsConstraintGenerator::set_ligand_atom_name( std::string setting){
	ligand_atom_name_ = setting;
}

void
MetalContactsConstraintGenerator::set_use_ligand_selector( bool setting){
	use_ligand_selector_ = setting;
}

void
MetalContactsConstraintGenerator::set_ligand_selector( core::select::residue_selector::ResidueSelectorCOP setting){
	ligand_selector_ = setting;
}

void
MetalContactsConstraintGenerator::set_ligand_resnum_string( std::string setting){
	ligand_resnum_string_ = setting;
}

void
MetalContactsConstraintGenerator::set_use_contact_selector( bool setting){
	use_contact_selector_ = setting;
}

void
MetalContactsConstraintGenerator::set_contact_selector(core::select::residue_selector::ResidueSelectorCOP setting){
	contact_selector_ = setting;
}

void
MetalContactsConstraintGenerator::set_contact_resnum_string( std::string setting){
	contact_resnum_string_ = setting;
}

void
MetalContactsConstraintGenerator::set_base_atom_name(std::string setting){
	base_atom_name_ = setting;
}

void
MetalContactsConstraintGenerator::set_base_base_atom_name(  std::string setting){
	base_base_atom_name_ = setting;
}

void
MetalContactsConstraintGenerator::set_ideal_distance( core::Real setting){
	ideal_distance_ = setting;
}

void
MetalContactsConstraintGenerator::set_ideal_angle_about_contact( utility::vector1< core::Real > setting ){
	ideal_angle_about_contact_ = setting;
}

void
MetalContactsConstraintGenerator::set_ideal_dihedral_about_contact( utility::vector1< core::Real > setting){
	ideal_dihedral_about_contact_ = setting;
}

void
MetalContactsConstraintGenerator::set_ideal_angle_about_metal( utility::vector1< core::Real > setting ){
	ideal_angle_about_metal_ = setting;
}

void
MetalContactsConstraintGenerator::set_ideal_dihedral_about_metal( utility::vector1< core::Real > setting){
	ideal_dihedral_about_metal_ = setting;
}

void
MetalContactsConstraintGenerator::set_ideal_dihedral_3( utility::vector1< core::Real > setting){
	ideal_dihedral_3_ = setting;
}

void
MetalContactsConstraintGenerator::set_score_against_internal_contacts(bool setting){
	score_against_internal_contacts_ = setting;
}


//PRIVATE METHODS

//Function to get ligand resnum
///@brief Uses private data to compute ligand resnums based on selector/resnum string and pose
core::Size
MetalContactsConstraintGenerator::get_ligand_resnum( core::pose::Pose const & pose ) const{
	bool found = false;
	core::Size return_val = 0;
	if ( ligand_selector_ && use_ligand_selector_ ) {
		core::select::residue_selector::ResidueSubset lig_subset( pose.total_residue(), false );
		lig_subset = ligand_selector_->apply( pose );
		for ( Size ii = 1; ii <= lig_subset.size(); ++ii ) {
			if ( lig_subset[ ii ] ) {
				if ( found ) {
					utility_exit_with_message( "ERROR: Ligand selector must specify only one residue!" );
				}
				return_val = ii;
				found = true;
			}
		}
		return return_val;
	}
	//Otherwise take from string
	return_val = core::pose::parse_resnum( ligand_resnum_string_, pose, true );
	return return_val;
}

//Function to get selector resnums
///@brief Uses private data to compute contact residue numbers based on selector/resnum string and pose; inserts them into the provided set
void
MetalContactsConstraintGenerator::get_contact_resnums( core::pose::Pose const & pose, std::set< core::Size > & output ) const{
	if ( contact_selector_ && use_contact_selector_ ) {
		core::select::residue_selector::ResidueSubset contact_subset( pose.total_residue(), false );
		contact_subset = contact_selector_->apply( pose );
		for ( Size ii = 1; ii <= contact_subset.size(); ++ii ) {
			if ( contact_subset[ ii ] ) {
				output.insert( ii );
				TR.Debug << "Found contact residue " << ii << std::endl;
			}
		}
	} else { // grab from string
		std::set< Size > const res_vec( get_resnum_list( contact_resnum_string_, pose ) );
		output.insert( res_vec.begin(), res_vec.end() );
		for ( core::Size num: output ) {
			if ( num > pose.total_residue() || num == 0 ) {
				std::stringstream err_msg;
				err_msg << "Residue " << num << " not found in pose!\n";
				throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg.str() );
			}
		}
	}
}






} //protocols
} //constraint_generator
