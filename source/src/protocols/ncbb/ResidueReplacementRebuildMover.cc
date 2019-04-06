// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ncbb/ResidueReplacementRebuildMover.cc
/// @brief A simple method to mutate a pose fairly destructively, and barely recover
/// @author Andy Watkins (andy.watkins2@gmail.com)

// Unit headers
#include <protocols/ncbb/ResidueReplacementRebuildMover.hh>
#include <protocols/ncbb/ResidueReplacementRebuildMoverCreator.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/variant_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <protocols/minimization_packing/MinMover.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.ncbb.ResidueReplacementRebuildMover" );

namespace protocols {
namespace ncbb {

using namespace core::pose;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ResidueReplacementRebuildMover::ResidueReplacementRebuildMover():
	protocols::moves::Mover( ResidueReplacementRebuildMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
ResidueReplacementRebuildMover::ResidueReplacementRebuildMover( ResidueReplacementRebuildMover const & src ):
	protocols::moves::Mover( src )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ResidueReplacementRebuildMover::~ResidueReplacementRebuildMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/*
ResidueReplacementRebuildMover::make_poses():
poses_for_later = []

pose = pose_from_pdb("1fok_0001.pdb")
for resi_chosen in range(43, 595):#pose.size()): # starting residue. actually, only do one for now.
print(pose.residue_type(resi_chosen-1).name())
print(pose.residue_type(resi_chosen).name())
print(pose.residue_type(resi_chosen+1).name())
#for resi_chosen in range(214, 215): # starting residue. actually, only do one for now.
for resi_len in range(4):
if resi_chosen + resi_len > 594: continue#pose.size()-1: continue
for resname in [ "ORTHO_POLYARAMID_ALA", "PRE_METHYLENE_ORTHO_POLYARAMID_ALA", "POST_METHYLENE_ORTHO_POLYARAMID_ALA",
"META_POLYARAMID_ALA", "PRE_METHYLENE_META_POLYARAMID_ALA", "POST_METHYLENE_META_POLYARAMID_ALA",
"PARA_POLYARAMID_ALA", "PRE_METHYLENE_PARA_POLYARAMID_ALA", "POST_METHYLENE_PARA_POLYARAMID_ALA",
"PRE_METHYLENE_POST_METHYLENE_ORTHO_POLYARAMID_ALA", "PRE_METHYLENE_POST_METHYLENE_META_POLYARAMID_ALA", "PRE_METHYLENE_POST_METHYLENE_PARA_POLYARAMID_ALA"]:
# copy pose residue by residue
print("Testing {}, substituting for {} and {} subsequent residues.".format(resname, resi_chosen, resi_len))

newpose = make_new_pose(pose, resi_chosen, resname, resi_len)

# ok, min
premin = sfxn(newpose)


movemap = rosetta.core.kinematics.MoveMap()
movemap.set_jump(True)
movemap.set_bb(False)
movemap.set_chi(True)
#adjust_movemap_backbone_for_aramid(movemap, newpose, resi_chosen)
for resi in range(max(2,resi_chosen-2), min(newpose.size(), resi_chosen+resi_len+2)):
if newpose.residue_type(resi).is_protein():
movemap.set_bb(resi, True)
elif newpose.residue_type(resi).is_aramid():
movemap.set_bb(resi, True)
movemap.set(DOF_ID(AtomID(newpose.residue_type(resi).atom_index('N'), resi), D), True)

#print(movemap)
minmover = rosetta.protocols.minimization_packing.MinMover(movemap, sfxn, "lbfgs_armijo_nonmonotone", 0.001, True)
minmover.cartesian(True)
print("  Built pose and about to minimize")
#newpose.dump_pdb('premin.pdb')
minmover.apply(newpose)
print("  Minimized from {} to {}.".format(premin, sfxn(newpose)))
#newpose.dump_pdb('min.pdb')
poses_for_later.append({'pose': Pose(newpose), 'score': sfxn(newpose), 'resi': resi_chosen, 'resname': resname, 'resi_len': resi_len})
return poses_for_later
*/

Pose
ResidueReplacementRebuildMover::make_new_pose(
	Pose const & pose
) const {
	Pose newpose;

	auto rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	newpose.append_residue_by_jump( pose.residue( 1 ), 1 );
	TR.Trace << "Appended first residue by jump." << std::endl;
	for ( Size resi = 2; resi < resi_chosen_; ++resi ) {
		if ( pose.residue_type(resi).is_protein() && pose.residue_type(resi-1).is_DNA() ) {
			newpose.append_residue_by_jump(pose.residue(resi), newpose.size());
		} else if ( pose.residue_type(resi).is_DNA() && pose.residue_type(resi-1).is_protein() ) {
			newpose.append_residue_by_jump(pose.residue(resi), newpose.size());
		} else if ( pose.residue_type(resi-1).has_variant_type(core::chemical::UPPER_TERMINUS_VARIANT) ) {
			newpose.append_residue_by_jump(pose.residue(resi), newpose.size());
		} else {
			newpose.append_residue_by_bond(pose.residue(resi), false);
		}
		TR.Trace << "Appended residue " << resi << " by bond or jump as appropriate. New pose size now " << newpose.size() << "." << std::endl;
	}

	// instead of 214, add an aramid
	newpose.append_residue_by_bond( core::conformation::Residue( rts->name_map( resname_ ), true), true);
	TR.Trace << "Appended residue " << resname_ << " by bond. New pose size now " << newpose.size() << "." << std::endl;

	//newpose.dump_pdb('prefoo.pdb')
	newpose.append_residue_by_jump( pose.residue( resi_chosen_+1+resi_len_ ), newpose.size() );
	TR.Trace << "Appended residue " << resi_chosen_+1+resi_len_  << " by jump. New pose size now " << newpose.size() << "." << std::endl;

	for ( Size resi = resi_chosen_ + 2 + resi_len_; resi <= pose.size(); ++resi ) {
		if ( pose.residue_type(resi).is_protein() and pose.residue_type(resi-1).is_DNA() ) {
			newpose.append_residue_by_jump(pose.residue(resi), newpose.size());
		} else if ( pose.residue_type(resi).is_DNA() and pose.residue_type(resi-1).is_protein() ) {
			newpose.append_residue_by_jump(pose.residue(resi), newpose.size());
		} else if ( pose.residue_type(resi-1).has_variant_type(core::chemical::UPPER_TERMINUS_VARIANT) ) {
			newpose.append_residue_by_jump(pose.residue(resi), newpose.size());
		} else {
			newpose.append_residue_by_bond(pose.residue(resi), false);
		}
		TR.Trace << "Appended residue " << resi << " by bond or jump as appropriate. New pose size now " << newpose.size() << "." << std::endl;
	}

	// OK. Now we have a setup.
	// Correctly add cutpoint variants to 214.
	TR.Debug << "About to add cutpoint variants to " << resi_chosen_ << " which is " << newpose.residue_type( resi_chosen_ ).name() << " and the next res" << std::endl;
	core::pose::correctly_add_cutpoint_variants( newpose, resi_chosen_ );

	return newpose;
}

/// @brief Apply the mover
void
ResidueReplacementRebuildMover::apply( core::pose::Pose & pose ) {

	using namespace protocols::minimization_packing;
	using namespace core::kinematics;

	TR << "About to replace residues " << resi_chosen_ << " through " << resi_chosen_ + resi_len_ << std::endl;
	TR << "with one " << resname_ << "." << std::endl;

	pose = make_new_pose( pose );

	if ( !scorefxn_ ) {
		TR.Error << "Scorefunction not specified?" << std::endl;
		utility_exit();
	}
	core::Real premin = ( *scorefxn_ )( pose );

	auto movemap = utility::pointer::make_shared< MoveMap >();
	movemap->set_jump( true );
	movemap->set_bb( false );
	movemap->set_chi( true );
	using namespace core::id;
	for ( Size resi = std::max( Size( 2 ), Size( resi_chosen_ - 2 ) );
			resi < std::min( Size( pose.size() ), Size( resi_chosen_ + resi_len_ + 2 ) );
			++resi ) {
		if ( pose.residue_type( resi ).is_protein() ) {
			movemap->set_bb( resi, true );
		} else if ( pose.residue_type( resi ).is_aramid() ) {
			movemap->set_bb( resi, true );
		}
		movemap->set( DOF_ID( AtomID( pose.residue_type( resi ).atom_index( "N" ), resi ), D ), true );
	}

	//print(movemap)
	auto minmover = utility::pointer::make_shared< MinMover >( movemap, scorefxn_, "lbfgs_armijo_nonmonotone", 0.001, true );
	minmover->cartesian( true );
	TR << "\tBuilt pose and about to minimize" << std::endl;
	//newpose.dump_pdb('premin.pdb')
	minmover->apply( pose );
	TR << "\tMinimized from " << premin << " to " << ( *scorefxn_ )( pose );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
ResidueReplacementRebuildMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
ResidueReplacementRebuildMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	if ( tag->hasOption("scorefxn") ) {
		score_function( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	}

	//resi_chosen_ =
	set_resi_chosen( tag->getOption< Size >( "seqpos" ) );
	//resi_len_ =
	set_resi_len( tag->getOption< Size >( "nres_additional", 0 ) );
	//resname_ =
	set_resname( tag->getOption< std::string >( "residue_name" ) );

}

void ResidueReplacementRebuildMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "seqpos", xsct_non_negative_integer, "The residue to delete" )
		+ XMLSchemaAttribute( "nres_additional", xsct_non_negative_integer, "Additional adjacent residues, beyond the seqpos specified, to delete" )
		+ XMLSchemaAttribute( "residue_name", xs_string, "The residue name for replacement" );

	rosetta_scripts::attributes_for_parse_score_function( attlist );


	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "DOCUMENTATION STRING", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ResidueReplacementRebuildMover::fresh_instance() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< ResidueReplacementRebuildMover >() );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ResidueReplacementRebuildMover::clone() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< ResidueReplacementRebuildMover >( *this ) );
}

std::string ResidueReplacementRebuildMover::get_name() const {
	return mover_name();
}

std::string ResidueReplacementRebuildMover::mover_name() {
	return "ResidueReplacementRebuildMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
ResidueReplacementRebuildMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< ResidueReplacementRebuildMover >() );
}

std::string
ResidueReplacementRebuildMoverCreator::keyname() const
{
	return ResidueReplacementRebuildMover::mover_name();
}

void ResidueReplacementRebuildMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueReplacementRebuildMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, ResidueReplacementRebuildMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //ncbb
