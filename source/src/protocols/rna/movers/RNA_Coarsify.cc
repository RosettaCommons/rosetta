// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/movers/RNA_Coarsify.cc
/// @brief Make a fullatom RNA pose into a coarse representation.
/// @author Andy Watkins (andy.watkins2@gmail.com)

// Unit headers
#include <protocols/rna/movers/RNA_Coarsify.hh>
#include <protocols/rna/movers/RNA_CoarsifyCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.rna.movers.RNA_Coarsify" );

namespace protocols {
namespace rna {
namespace movers {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
RNA_Coarsify::RNA_Coarsify():
	protocols::moves::Mover( RNA_Coarsify::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
RNA_Coarsify::RNA_Coarsify( RNA_Coarsify const & src ):
	protocols::moves::Mover( src )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
RNA_Coarsify::~RNA_Coarsify(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
RNA_Coarsify::apply( core::pose::Pose & pose ){
	// Take a pose and return a coarse version.

	core::pose::PoseOP coarse_pose = utility::pointer::make_shared< Pose >();
	core::import_pose::make_coarse_pose( pose, *coarse_pose );
	pose = *coarse_pose;

	// using namespace core;
	// using namespace core::scoring;
	// using namespace core::scoring::func;
	// using namespace core::scoring::constraints;

	// ScoreFunctionOP scorefxn = utility::pointer::make_shared< ScoreFunction >();
	// scorefxn->set_weight( core::scoring::coordinate_constraint, 1.0 );
	// scorefxn->set_weight( core::scoring::rna_vdw, 1.0 );


	// core::pose::PoseOP coarse_pose = utility::pointer::make_shared< Pose >();

	// auto coarse_rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( "coarse_rna" );
	// FuncOP constraint_func = utility::pointer::make_shared< FlatHarmonicFunc >( 0, 0.25, 1 );

	// for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
	//  auto const & coarse_rt = coarse_rts->name_map( pose.residue_type( ii ).name() );
	//  if ( ii == 1 || !pose.residue( ii ).is_bonded( ii - 1 ) ) {
	//   coarse_pose->append_residue_by_jump( core::conformation::Residue( coarse_rt, true ), 1 );
	//  } else {
	//   coarse_pose->append_residue_by_bond( core::conformation::Residue( coarse_rt, true ), true );
	//  }


	//  core::pose::addVirtualResAsRoot( *coarse_pose );
	//  // core::pose::addVirtualResAsRoot( copy_pose );
	//  core::Size const my_anchor( coarse_pose->fold_tree().root() );

	//  // ConstraintOP constraint = utility::pointer::make_shared< CoordinateConstraint >(
	//  //  core::id::AtomID( coarse_pose->residue_type( coarse_pose->size()-1 ).atom_index( "P" ), coarse_pose->size()-1 ),
	//  //  core::id::AtomID( 1, my_anchor ),
	//  //  pose.residue( coarse_pose->size()-1 ).xyz( pose.residue_type( coarse_pose->size()-1 ).atom_index( "P" ) ),
	//  //  constraint_func );
	//  // coarse_pose->add_constraint( constraint );

	//  // constraint = utility::pointer::make_shared< CoordinateConstraint >(
	//  //  core::id::AtomID( coarse_pose->residue_type( coarse_pose->size()-1 ).atom_index( "S" ), coarse_pose->size()-1 ),
	//  //  core::id::AtomID( 1, my_anchor ),
	//  //  pose.residue( coarse_pose->size()-1 ).xyz( pose.residue_type( coarse_pose->size()-1 ).atom_index( "C1'" ) ),
	//  //  constraint_func );
	//  // coarse_pose->add_constraint( constraint );

	//  // if ( coarse_rt.name3() == "  A" || coarse_rt.name3() == "  G" ) {
	//  //  constraint = utility::pointer::make_shared< CoordinateConstraint >(
	//  //   core::id::AtomID( coarse_pose->residue_type( coarse_pose->size()-1 ).atom_index( "CEN" ), coarse_pose->size()-1 ),
	//  //   core::id::AtomID( 1, my_anchor ),
	//  //   pose.residue( coarse_pose->size()-1 ).xyz( pose.residue_type( coarse_pose->size()-1 ).atom_index( "N9" ) ),
	//  //   constraint_func );
	//  //  coarse_pose->add_constraint( constraint );
	//  // } else {
	//  //  constraint = utility::pointer::make_shared< CoordinateConstraint >(
	//  //   core::id::AtomID( coarse_pose->residue_type( coarse_pose->size()-1 ).atom_index( "CEN" ), coarse_pose->size()-1 ),
	//  //   core::id::AtomID( 1, my_anchor ),
	//  //   pose.residue( coarse_pose->size()-1 ).xyz( pose.residue_type( coarse_pose->size()-1 ).atom_index( "N1" ) ),
	//  //   constraint_func );
	//  //  coarse_pose->add_constraint( constraint );
	//  // }

	//  for ( Size jj = 1; jj <= coarse_pose->size() - 1; ++jj ) {
	//   ConstraintOP constraint = utility::pointer::make_shared< CoordinateConstraint >(
	//    core::id::AtomID( coarse_pose->residue_type( jj ).atom_index( "P" ), jj ),
	//    core::id::AtomID( 1, my_anchor ),
	//    pose.residue( jj ).xyz( pose.residue_type( jj ).atom_index( "P" ) ),
	//    constraint_func );
	//   coarse_pose->add_constraint( constraint );

	//   constraint = utility::pointer::make_shared< CoordinateConstraint >(
	//    core::id::AtomID( coarse_pose->residue_type( jj ).atom_index( "S" ), jj ),
	//    core::id::AtomID( 1, my_anchor ),
	//    core::pose::rna::get_sugar_centroid( pose.residue( jj ) ),
	//    constraint_func );
	//   coarse_pose->add_constraint( constraint );
	//   constraint = utility::pointer::make_shared< CoordinateConstraint >(
	//    core::id::AtomID( coarse_pose->residue_type( jj ).atom_index( "CEN" ), jj ),
	//    core::id::AtomID( 1, my_anchor ),
	//    core::chemical::rna::get_rna_base_centroid( pose.residue( jj ) ),
	//    constraint_func );
	//   coarse_pose->add_constraint( constraint );
	//   // if ( coarse_rt.name3() == "  A" || coarse_rt.name3() == "  G" ) {
	//   //  constraint = utility::pointer::make_shared< CoordinateConstraint >(
	//   //   core::id::AtomID( coarse_pose->residue_type( jj ).atom_index( "CEN" ), jj ),
	//   //   core::id::AtomID( 1, my_anchor ),
	//   //   pose.residue( jj ).xyz( pose.residue_type( jj ).atom_index( "N9" ) ),
	//   //   constraint_func );
	//   //  coarse_pose->add_constraint( constraint );
	//   // } else {
	//   //  constraint = utility::pointer::make_shared< CoordinateConstraint >(
	//   //   core::id::AtomID( coarse_pose->residue_type( jj ).atom_index( "CEN" ), jj ),
	//   //   core::id::AtomID( 1, my_anchor ),
	//   //   pose.residue( jj ).xyz( pose.residue_type( jj ).atom_index( "N1" ) ),
	//   //   constraint_func );
	//   //  coarse_pose->add_constraint( constraint );
	//   // }
	//  }


	//  kinematics::MoveMapOP suite_mm = utility::pointer::make_shared< kinematics::MoveMap >();
	//  suite_mm->set_bb( true );
	//  suite_mm->set_chi( true );
	//  // for ( Size kk = 0; kk <= 4; ++kk ) {
	//  //  suite_mm->set_bb( coarse_pose->size() - kk, true );
	//  //  suite_mm->set_chi( coarse_pose->size() - kk, true );
	//  // }
	//  suite_mm->set_jump( true );

	//  protocols::minimization_packing::MinMoverOP minm = utility::pointer::make_shared< protocols::minimization_packing::MinMover >( suite_mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true );
	//  TR << "Minimized coarse pose " << ii << " from " << ( *scorefxn )( *coarse_pose ) << " to ";
	//  minm->apply( *coarse_pose );
	//  TR << ( *scorefxn )( *coarse_pose ) << "." << std::endl;

	//  core::pose::remove_virtual_residues( *coarse_pose );
	//  std::stringstream ss;
	//  ss << "min_" << ii << ".pdb";
	//  coarse_pose->dump_scored_pdb( ss.str(), *scorefxn );
	// }


	pose = *coarse_pose;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
RNA_Coarsify::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
RNA_Coarsify::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap&
) {

}
void RNA_Coarsify::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "DOCUMENTATION STRING", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
RNA_Coarsify::fresh_instance() const
{
	return utility::pointer::make_shared< RNA_Coarsify >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
RNA_Coarsify::clone() const
{
	return utility::pointer::make_shared< RNA_Coarsify >( *this );
}

std::string RNA_Coarsify::get_name() const {
	return mover_name();
}

std::string RNA_Coarsify::mover_name() {
	return "RNA_Coarsify";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
RNA_CoarsifyCreator::create_mover() const
{
	return utility::pointer::make_shared< RNA_Coarsify >();
}

std::string
RNA_CoarsifyCreator::keyname() const
{
	return RNA_Coarsify::mover_name();
}

void RNA_CoarsifyCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RNA_Coarsify::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, RNA_Coarsify const & mover )
{
	mover.show(os);
	return os;
}

} //movers
} //rna
} //protocols
