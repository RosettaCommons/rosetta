// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Assigns a ConstraintSet to a pose. Reads and creats ConstraintSet from file via command line option -constraints::cst_file, unless a ConstraintSet is supplied via the constructor or the constraint_set() method.
/// @author ashworth

#include <protocols/constraint_movers/ConstraintSetMover.hh>
#include <protocols/constraint_movers/ConstraintSetMoverCreator.hh>


#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/datacache/DataMapObj.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

//Auto Headers

#include <utility/io/izstream.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.constraint_movers.ConstraintSetMover" );

namespace protocols {
namespace constraint_movers {

using namespace core;
using namespace basic::options;
using namespace scoring;
using namespace constraints;
using namespace utility::tag;




ConstraintSetMover::ConstraintSetMover()
: protocols::moves::Mover( ConstraintSetMover::mover_name() )
{
	read_options();
}

ConstraintSetMover::~ConstraintSetMover()= default;

ConstraintSetMover::ConstraintSetMover( std::string const & type )
: protocols::moves::Mover(type)
{
	read_options();
}

void
ConstraintSetMover::read_options()
{
	if ( option[ OptionKeys::constraints::cst_file ].user() ) {
		cst_file_ = option[ OptionKeys::constraints::cst_file ]().front();
	}

	if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		cst_fa_file_ = option[ OptionKeys::constraints::cst_fa_file ]().front();
	} else {
		cst_fa_file_ = cst_file_;
	}

	add_constraints_ = false; // This is never used, lol.
}

void
ConstraintSetMover::constraint_set( ConstraintSetCOP cst_set )
{
	constraint_set_low_res_ = utility::pointer::make_shared< ConstraintSet >( *cst_set );
	constraint_set_high_res_ = utility::pointer::make_shared< ConstraintSet >( *cst_set );
}

void
ConstraintSetMover::constraint_file( std::string const & cst_file )
{
	cst_file_ = cst_file;
	cst_fa_file_ = cst_file;
}


ConstraintSetOP ConstraintSetMover::constraint_set() { return constraint_set_low_res_; }
ConstraintSetCOP ConstraintSetMover::constraint_set() const { return constraint_set_low_res_; }

void
ConstraintSetMover::apply( Pose & pose )
{
	if ( cst_file_ == "print" ) {
		std::cout << "*****************showing constraint set************************"  <<std::endl;
		pose.constraint_set()->show_definition(std::cout,pose);
	} else {
		if ( !constraint_set_low_res_ && !pose.is_fullatom() ) {
			// uninitialized filename not tolerated, in order to avoid potential confusion
			if ( cst_file_.empty() ) {
				// We catch this error later if file_to_input_map_ is used
				if ( file_to_input_map_.empty() ) utility_exit_with_message("Can\'t read constraints from empty file!");
			} else if ( cst_file_ == "none" ) constraint_set_low_res_ = utility::pointer::make_shared< ConstraintSet >();
			// special case: set cst_file_ to "none" to effectively remove constraints from Pose
			else {
				constraint_set_low_res_ =
					ConstraintIO::get_instance()->read_constraints( cst_file_, utility::pointer::make_shared< ConstraintSet >(), pose );
				//ConstraintIO::get_instance()->read_constraints_new( cst_file_, new ConstraintSet, pose );
			}
		}
		if ( !constraint_set_high_res_ && pose.is_fullatom() ) {
			// uninitialized filename not tolerated, in order to avoid potential confusion
			if ( cst_file_.empty() ) {
				// We catch this error later if file_to_input_map_ is used
				if ( file_to_input_map_.empty() ) utility_exit_with_message("Can\'t read constraints from empty file!");
			} else if ( cst_fa_file_ == "none" ) constraint_set_high_res_ = utility::pointer::make_shared< ConstraintSet >();
			// special case: set cst_file_ to "none" to effectively remove constraints from Pose
			else {
				constraint_set_high_res_ =
					ConstraintIO::get_instance()->read_constraints( cst_fa_file_, utility::pointer::make_shared< ConstraintSet >(), pose );
				//ConstraintIO::get_instance()->read_constraints_new( cst_fa_file_, new ConstraintSet, pose );
			}
		}
		// Use the cst_set unless the input tag is in the file_to_input_map_
		ConstraintSetOP constraints_to_use = nullptr;
		if ( pose.is_fullatom() ) {
			constraints_to_use = constraint_set_high_res_;
		} else {
			constraints_to_use = constraint_set_low_res_;
		}
		if ( ! file_to_input_map_.empty() ) {
			if ( file_to_input_map_.count( current_input_name_ ) ) {
				std::string const & fname = file_to_input_map_.at( current_input_name_ );
				TR << "Loading constraints from " << fname << " for tag \"" << current_input_name_ << "\"" << std::endl;
				constraints_to_use = ConstraintIO::get_instance()->read_constraints( fname, ConstraintSetOP( new ConstraintSet ), pose );
			} else {
				// This serves a a way to figure out what the tag is and to not silently use the default if you give the wrong tag.
				TR << "Using default constraints for tag \"" << current_input_name_ << "\"" << std::endl;
				if ( ! constraints_to_use ) {
					utility_exit_with_message( "ConstraintSetMover: Default constraints were never set!" );
				}
			}
		}

		if ( pose.is_fullatom() ) {
			if ( add_constraints() ) {
				pose.add_constraints( constraints_to_use->get_all_constraints());
			} else {
				pose.constraint_set( constraints_to_use );
			}
			if ( TR.Debug.visible() ) {
				TR.Debug << "High-res constraints:" << std::endl;
				constraints_to_use->show_definition(TR.Debug, pose);
				TR.Debug.flush(); // Make sure the tracer is output if the definitions don't use a std::endl
			}
		} else {
			if ( add_constraints() ) {
				pose.add_constraints( constraints_to_use->get_all_constraints() );
			} else {
				pose.constraint_set( constraints_to_use );
			}
			if ( TR.Debug.visible() ) {
				TR.Debug << "Low-res constraints:" << std::endl;
				constraints_to_use->show_definition(TR.Debug, pose);
				TR.Debug.flush(); // Make sure the tracer is output if the definitions don't use a std::endl
			}
		}
	}
}


void ConstraintSetMover::set_cst_map_file( std::string const & file_name ) {
	std::string contents;
	utility::io::izstream in_file( file_name );
	if ( ! in_file ) {
		utility_exit_with_message( "ConstraintSetMover: couldn't open file: \"" + file_name + "\"" );
	}
	utility::slurp( in_file, contents );
	in_file.close();
	set_cst_map_file_contents( contents );
}

void ConstraintSetMover::set_cst_map_file_contents( std::string const & file_contents ) {

	std::stringstream ss( file_contents );

	std::string line;

	while ( getline( ss, line ) ) {

		std::string save_line = line;

		// delete comments and whitespace
		line = line.substr( 0, line.find("#") );
		utility::trim( line, " \n\r\t" );

		if ( line.length() == 0 ) continue;

		utility::vector1< std::string > split = utility::split_whitespace( line );

		if ( split.size() != 2 ) {
			utility_exit_with_message( "ConstraintSetMover: invalid line in cst_map_file: \"" + save_line + "\"" );
		}

		std::string const & file = split[1];
		std::string const & tag = split[2];

		if ( file_to_input_map_.count( tag ) ) {
			utility_exit_with_message( "ConstraintSetMover: tag specified twice in cst_map_file: \"" + tag + "\"" );
		}
		file_to_input_map_[ tag ] = file;
	}

}


protocols::moves::MoverOP ConstraintSetMover::clone() const { return utility::pointer::make_shared< protocols::constraint_movers::ConstraintSetMover >( *this ); }
protocols::moves::MoverOP ConstraintSetMover::fresh_instance() const { return utility::pointer::make_shared< ConstraintSetMover >(); }

void
ConstraintSetMover::register_options()
{
	option.add_relevant( OptionKeys::constraints::cst_file );
	option.add_relevant( OptionKeys::constraints::cst_fa_file );
}

void
ConstraintSetMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data
)
{
	if ( tag->hasOption("cst_file") ) cst_file_ = tag->getOption<std::string>("cst_file");
	if ( tag->hasOption("cst_fa_file") ) cst_fa_file_ = tag->getOption<std::string>("cst_fa_file");
	else cst_fa_file_=cst_file_;
	add_constraints( tag->getOption< bool >( "add_constraints", false ) );
	TR << "of type ConstraintSetMover with constraint file: " << cst_file_ << std::endl;
	if ( cst_fa_file_ != cst_file_ ) {
		TR << "of type ConstraintSetMover with fullatom constraint file: " << cst_fa_file_ << std::endl;
	}

	if ( tag->hasOption("cst_map_file") ) set_cst_map_file( tag->getOption<std::string>("cst_map_file") );
	current_input_name_ = "";
	if ( data.has( "strings", "current_input_name" ) ) {
		using StringWrapper = basic::datacache::DataMapObj< std::string >;
		auto ptr = data.get_ptr< StringWrapper >( "strings", "current_input_name" );
		runtime_assert( ptr );
		current_input_name_ = ptr->obj;
	}
}

std::string ConstraintSetMover::get_name() const {
	return mover_name();
}

std::string ConstraintSetMover::mover_name() {
	return "ConstraintSetMover";
}

void ConstraintSetMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute( "cst_file", xs_string, "(implicitly centroid compatible) constraint file; attempts to autodetect if pose is centroid and uses these")
		+ XMLSchemaAttribute( "cst_fa_file", xs_string, "(implicitly fullatom compatible) constraint file; attempts to autodetect if pose is fullatom and uses these.  If not supplied, uses the value in cst_file")
		+ XMLSchemaAttribute( "cst_map_file", xs_string, "A file to map constraint sets to input tags. With one entry per line, the format is: /path/to/cst/file.cst tag")
		+ XMLSchemaAttribute::attribute_w_default( "add_constraints", xsct_rosetta_bool, "if True, ADD these constraints to the existing ConstraintSet in the Pose; if False, REPLACE the constraints in the Pose with these.", "false");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Adds constraints from a constraint file to the pose.", attlist );
}

std::string ConstraintSetMoverCreator::keyname() const {
	return ConstraintSetMover::mover_name();
}

protocols::moves::MoverOP
ConstraintSetMoverCreator::create_mover() const {
	return utility::pointer::make_shared< ConstraintSetMover >();
}

void ConstraintSetMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ConstraintSetMover::provide_xml_schema( xsd );
}


} // moves
} // protocols
