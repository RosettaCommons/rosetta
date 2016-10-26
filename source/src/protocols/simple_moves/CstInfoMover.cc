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

#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/MultiConstraint.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/izstream.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.CstInfoMover" );


namespace protocols {
namespace simple_moves {

CstInfoMover::CstInfoMover():
	protocols::moves::Mover( "CstInfoMover" ),
	prefix_( "CST" ),
	recursive_( false )
{}

CstInfoMover::~CstInfoMover(){}

CstInfoMover::CstInfoMover( CstInfoMover const & src ):
	protocols::moves::Mover( src ),
	cst_file_( src.cst_file_ ),
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
	cst_file( tag->getOption< std::string >( "cst_file", "" ) );
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

std::string
CstInfoMover::get_name() const {
	return "CstInfoMover";
}

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
			TR.Error << "ERROR Unable to read constraint " << all_constraints.size() + 1 << " from file " << filename << std::endl;
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

protocols::moves::MoverOP
CstInfoMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new CstInfoMover );
}

std::string
CstInfoMoverCreator::keyname() const {
	return CstInfoMoverCreator::mover_name();
}

std::string
CstInfoMoverCreator::mover_name(){
	return "CstInfoMover";
}

} //protocols
} //simple_moves


