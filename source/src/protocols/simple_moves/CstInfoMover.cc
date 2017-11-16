// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/simple_moves/CstInfoMover.cc
/// @brief A Mover to output information about constraints
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/simple_moves/CstInfoMover.hh>
#include <protocols/simple_moves/CstInfoMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/MultiConstraint.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/izstream.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.simple_moves.CstInfoMover" );


namespace protocols {
namespace simple_moves {

CstInfoMover::CstInfoMover():
	protocols::moves::Mover( "CstInfoMover" ),
	cst_file_( std::string() ),
	dump_cst_file_( std::string() ),
	prefix_( "CST" ),
	recursive_( false )
{}

CstInfoMover::~CstInfoMover(){}

CstInfoMover::CstInfoMover( CstInfoMover const & src ):
	protocols::moves::Mover( src ),
	cst_file_( src.cst_file_ ),
	dump_cst_file_( src.dump_cst_file_ ),
	prefix_( src.prefix_ ),
	recursive_( src.recursive_ )
{}

void
CstInfoMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	cst_file( tag->getOption< std::string >( "cst_file", std::string() ) );
	dump_cst_file( tag->getOption< std::string >( "dump_cst_file", std::string() ) );
	prefix( tag->getOption< std::string >( "prefix", "CST" ) );
	recursive( tag->getOption<bool>( "recursive", false ) );
}

protocols::moves::MoverOP
CstInfoMover::clone() const{
	return protocols::moves::MoverOP( new CstInfoMover( *this ) );
}

/*
CstInfoMover & CstInfoMoveroperator=( CstInfoMover const & src){
return CstInfoMover( src );
}
*/


moves::MoverOP
CstInfoMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new CstInfoMover );
}

// XRW TEMP std::string
// XRW TEMP CstInfoMover::get_name() const {
// XRW TEMP  return "CstInfoMover";
// XRW TEMP }

void
CstInfoMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
	output << "CstInfoMover using file " << cst_file_ << " and recursive " << recursive_ << " with prefix " << prefix_ << std::endl;
}

std::ostream &operator<< (std::ostream &os, CstInfoMover const &mover)
{
	mover.show(os);
	return os;
}


void
CstInfoMover::apply( core::pose::Pose& pose ){
	using namespace core::scoring::constraints;

	ConstraintCOPs all_constraints;

	if ( cst_file_.empty() ) {
		TR << "Loading constraints from pose: Warning - constraint ordering may be inconsistent." << std::endl;
		all_constraints = get_constraints_from_pose( pose );
	} else {
		all_constraints = get_constraints_from_file( cst_file_, pose );
	}

	if ( !dump_cst_file_.empty() ) {
		ConstraintSetOP cst_set( new ConstraintSet );
		cst_set->add_constraints( all_constraints );
		ConstraintIO::get_instance()->write_constraints( dump_cst_file_, *cst_set, pose );
	}

	TR << "Read " << all_constraints.size() << " constraints. " << std::endl;
	if ( cst_file_.empty() ) {
		TR << "------------ Constraints read from pose -------------" << std::endl;
	} else {
		TR << "------------ Constraints read for file " << cst_file_ << "-------------" << std::endl;
	}
	add_info_for_csts( all_constraints, pose, prefix_ );
	TR << "----------------------------------------------------" << std::endl;

}

void
CstInfoMover::add_info_for_csts( core::scoring::constraints::ConstraintCOPs const & cstlist, core::pose::Pose & pose, std::string const & tag ) const {
	using namespace core::scoring::constraints;

	for ( core::Size ii( 1 ); ii <= cstlist.size(); ++ii ) {
		ConstraintCOP const & cst( cstlist[ii] );
		core::Real measurement( cst->dist( pose ) );
		core::Real score( cst->score( pose ) );

		// Add the values to the pose extra scores
		std::string newtag( tag + "_" + utility::to_string( ii ) );
		TR << "# Constraint for '" << newtag << "' is:" << std::endl;
		cst->show_def( TR, pose );
		TR << std::endl; // Space for multi-constraints
		core::pose::setPoseExtraScore( pose, newtag + "_measure", measurement );
		core::pose::setPoseExtraScore( pose, newtag + "_score", score );

		if ( recursive_ ) {
			MultiConstraintCOP multi_cst( utility::pointer::dynamic_pointer_cast< MultiConstraint const >( cst )  );
			if ( multi_cst ) {
				ConstraintCOPs const & subcsts( multi_cst->member_constraints() );
				add_info_for_csts( subcsts, pose, newtag );
			}
		}
	}

}

core::scoring::constraints::ConstraintCOPs
CstInfoMover::get_constraints_from_file( std::string const & filename, core::pose::Pose const & pose ) const {
	using namespace core::scoring::constraints;

	ConstraintCOPs all_constraints;

	ConstraintIO & cstio( *ConstraintIO::get_instance() );

	utility::io::izstream data( filename.c_str() );
	if ( !data ) {
		utility_exit_with_message( "ERROR: In CstInfoMover, cannot open constraint file " + filename );
	}

	while ( data ) {
		ConstraintCOP cst_op( cstio.read_individual_constraint_new( data, pose, cstio.get_func_factory() ) );
		if ( cst_op ) {
			all_constraints.push_back( cst_op );
		} else if ( ! data.eof() ) {
			TR.Error << "Unable to read constraint " << all_constraints.size() + 1 << " from file " << filename << std::endl;
			utility_exit_with_message( "In CstInfoMover, unable to read cst." );
		}
	}

	return all_constraints;
}


core::scoring::constraints::ConstraintCOPs
CstInfoMover::get_constraints_from_pose( core::pose::Pose const & pose ) const {
	using namespace core::scoring::constraints;

	if ( ! pose.constraint_set() ) {
		ConstraintCOPs all_constraints;
		TR << "Asked to get constraints from pose which apparently doesn't have any constraints loaded." << std::endl;
		return all_constraints; // Empty constraints
	}

	// From what I can tell about ConstraintSet::get_all_constraints(), it should return the constraints in a stable order,
	// and that two poses with identical constraints should get the same order of the constraints.
	// This is probably the best we can ask for for pose-stored constraints.

	pose.constraint_set()->show(TR.Debug);

	return pose.constraint_set()->get_all_constraints();
}

/////////////// Creator ///////////////

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP CstInfoMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new CstInfoMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CstInfoMoverCreator::keyname() const {
// XRW TEMP  return CstInfoMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CstInfoMover::mover_name(){
// XRW TEMP  return "CstInfoMover";
// XRW TEMP }

std::string CstInfoMover::get_name() const {
	return mover_name();
}

std::string CstInfoMover::mover_name() {
	return "CstInfoMover";
}

void CstInfoMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute(
		"cst_file", xs_string,
		"Which constraints to report. If not given (or given as an empty string), "
		"report on the constraints that are currently in the pose. "
		"For consistency between different structures, it's highly recommended to "
		"supply a constraints file. In-pose constraints are not stored in any "
		"defined order, there's no guarantee that the constraint numbering will "
		"match up structure-to-structure, even if every pose has identical constraints. "
		"(If a constraint file is used, the order of the constraints will be the "
		"same as the order of the constraints in the constraint file.)" )
		+ XMLSchemaAttribute::attribute_w_default(
		"dump_cst_file", xs_string,
		"File name to WRITE the constraints into. Constraints are only printed if a file "
		" name is provided",
		std::string() )
		+ XMLSchemaAttribute::attribute_w_default(
		"prefix", xs_string,
		"What prefix to give the values. Important if you have more than one CstInfoMover "
		"(or more than one application), as later values with the same prefix "
		"will overwrite earlier values.",
		"CST" )
		+ XMLSchemaAttribute::attribute_w_default(
		"recursive", xsct_rosetta_bool,
		"Should sub-constraints of a multi-constraint also be output",
		"false" );

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Report information about how well the pose reflects given constraint values",
		attlist );
}

std::string CstInfoMoverCreator::keyname() const {
	return CstInfoMover::mover_name();
}

protocols::moves::MoverOP
CstInfoMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new CstInfoMover );
}

void CstInfoMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CstInfoMover::provide_xml_schema( xsd );
}


} //protocols
} //simple_moves
