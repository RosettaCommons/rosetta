// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/RenderGridsToKinemage.cc
/// @author Sam DeLuca

#include <protocols/qsar/RenderGridsToKinemage.hh>
#include <protocols/qsar/RenderGridsToKinemageCreator.hh>

#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <utility>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/tag/Tag.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/color_util.hh>

#include <utility/excn/Exceptions.hh>
#include <iostream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/qsar/scoring_grid/schema_util.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace qsar {

static basic::Tracer TR( "protocols.qsar.RenderGridsToKinemage" );





ColorGradient::ColorGradient(numeric::xyzVector<core::Real> const & value,
	core::Real const & lower,
	core::Real const & upper,
	std::string const & name)
{
	color_value = value;
	lower_bound = lower;
	upper_bound = upper;
	color_name = name;
}

RenderGridsToKinemage::RenderGridsToKinemage() :
	protocols::moves::Mover("RenderGridsToKinemage")
{}

RenderGridsToKinemage::RenderGridsToKinemage(
	scoring_grid::GridSetCOP grid_set_prototype,
	std::string const & grid_name,
	std::string const & filename):
	grid_set_prototype_(std::move( grid_set_prototype )),
	filename_( filename ),
	grid_name_( grid_name )
{}

RenderGridsToKinemage::RenderGridsToKinemage( RenderGridsToKinemage const & ) = default;

RenderGridsToKinemage::~RenderGridsToKinemage() = default;

protocols::moves::MoverOP RenderGridsToKinemage::clone() const
{
	return utility::pointer::make_shared< RenderGridsToKinemage >(*this);
}



void RenderGridsToKinemage::apply(core::pose::Pose & pose)
{
	debug_assert( grid_set_prototype_ );
	debug_assert( grid_name_ != "" );
	debug_assert( filename_ != "" );

	std::string chain( grid_set_prototype_->chain() );
	utility::vector1< core::Size > chain_residues( core::pose::get_resnums_for_chain( pose, chain ) );
	core::Vector center( core::pose::all_atom_center(pose, chain_residues) );

	scoring_grid::GridSetCOP grid_set = scoring_grid::GridManager::get_instance()->get_grids(*grid_set_prototype_, pose, center, chain);
	scoring_grid::SingleGridCOP grid = utility::pointer::dynamic_pointer_cast<scoring_grid::SingleGrid const>( grid_set->get_grid( grid_name_ ) );
	if ( ! grid ) {
		utility_exit_with_message( "RenderGridsToKinemage is currently unable to output the contents of a metagrid.  Sorry.");
	}

	setup_colors( *grid );

	//we can write multiple grids to one file, peek at the first line and see if the header has been written
	utility::io::izstream infile;
	infile.open(filename_);
	std::string first_word;
	infile >> first_word;
	infile.close();

	utility::io::ozstream outfile;
	if ( first_word == "@kinemage" ) {
		//header already written, append to existing file, don't write the header
		outfile.open(filename_,std::ios::app);
	} else {
		//header not written, open for overwrite, write the header
		outfile.open(filename_,std::ios::out);
		write_header(outfile);
	}
	write_colors(outfile);
	write_points(outfile, *grid);
	outfile.close();

}

void RenderGridsToKinemage::parse_my_tag(utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{

	grid_set_prototype_ = scoring_grid::parse_grid_set_from_tag( tag, data );

	if ( !tag->hasOption("grid_name") ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "RenderGridsToKinemage requires the'grid name' option");
	}

	grid_name_ = tag->getOption<std::string>("grid_name");
	filename_ = tag->getOption<std::string>("file_name");

	if ( !tag->hasOption("low_color") && !tag->hasOption("zero_color") && !tag->hasOption("high_color") ) {
		if ( tag->hasOption("color") ) {
			//one color
			color_mode_ = 1;
			color_ = numeric::comma_seperated_string_to_xyz<core::Real>(tag->getOption<std::string>("color"));

		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "To use RenderGridsToKinemage you must specify color for a 1 color grid plot, "
				"high_color and low_color for a two color gradient, or high_color, low_color and "
				"zero_color for a three color gradient centered at zero");
		}
	} else if ( tag->hasOption("low_color") && tag->hasOption("zero_color") && tag->hasOption("high_color") ) {
		//3 part gradient
		color_mode_ = 3;
		low_color_ = numeric::comma_seperated_string_to_xyz<core::Real>(tag->getOption<std::string>("low_color"));
		zero_color_ = numeric::comma_seperated_string_to_xyz<core::Real>(tag->getOption<std::string>("zero_color"));
		high_color_ = numeric::comma_seperated_string_to_xyz<core::Real>(tag->getOption<std::string>("high_color"));

	} else if ( tag->hasOption("low_color")&& tag->hasOption("high_color") ) {
		//2 part gradient
		color_mode_ = 2;
		low_color_ = numeric::comma_seperated_string_to_xyz<core::Real>(tag->getOption<std::string>("low_color"));
		high_color_ = numeric::comma_seperated_string_to_xyz<core::Real>(tag->getOption<std::string>("high_color"));
	} else {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "To use RenderGridsToKinemage you must specify color for a 1 color grid plot, "
			"high_color and low_color for a two color gradient, or high_color, low_color and "
			"zero_color for a three color gradient centered at zero");
	}

	gradient_bins_ = tag->getOption<core::Size>("gradient_bins",10);
	stride_ = tag->getOption<core::Size>("stride",2);
	std::cout << stride_ <<std::endl;

}


void RenderGridsToKinemage::setup_colors( scoring_grid::SingleGrid const & grid )
{
	if ( color_mode_ == 1 ) {
		setup_one_color_scheme(grid);
	} else if ( color_mode_ == 2 ) {
		setup_two_color_scheme(grid);
	} else if ( color_mode_ == 3 ) {
		setup_three_color_scheme(grid);
	} else {
		TR.Fatal << "Improper value of color mode found: " << color_mode_ << std::endl;
		utility_exit_with_message("Improper color mode in RenderGridsToKinemage."); //This should never happen
	}

}

void RenderGridsToKinemage::setup_one_color_scheme( scoring_grid::SingleGrid const & grid )
{
	core::Real grid_min = grid.get_min_value();
	core::Real grid_max = grid.get_max_value();

	ColorGradient new_color(color_,grid_min,grid_max,grid_name_+"_color");
	color_data_.push_back(new_color);
}

void RenderGridsToKinemage::setup_two_color_scheme( scoring_grid::SingleGrid const & grid )
{
	core::Real grid_min = grid.get_min_value();
	core::Real grid_max = grid.get_max_value();

	core::Real red_step = (high_color_.x() - low_color_.x())/static_cast<core::Real>(gradient_bins_);
	core::Real green_step = (high_color_.y() - low_color_.y())/static_cast<core::Real>(gradient_bins_);
	core::Real blue_step = (high_color_.z() - low_color_.z())/static_cast<core::Real>(gradient_bins_);

	core::Real value_step = (grid_max - grid_min)/gradient_bins_;


	core::Vector current_color = low_color_;
	core::Real current_low = grid_min;
	core::Real current_high = current_low+value_step;
	for ( core::Size bin = 1; bin <= gradient_bins_; ++bin ) {
		current_color.x(current_color.x()+red_step);
		current_color.y(current_color.y()+green_step);
		current_color.z(current_color.z()+blue_step);
		std::string color_bin_string(utility::to_string<core::Size>(bin));
		std::string color_name(grid_name_+"_"+color_bin_string+"_color");
		ColorGradient new_color(current_color,current_low,current_high,color_name);
		color_data_.push_back(new_color);
		current_low = current_high;
		current_high = current_low+value_step;
	}

}

void RenderGridsToKinemage::setup_three_color_scheme( scoring_grid::SingleGrid const & grid )
{
	core::Real grid_min = grid.get_min_value();
	core::Real grid_max = grid.get_max_value();
	if ( grid_min >= 0 ) {
		utility_exit_with_message("This grid does not have any scores greater than 0, a three color gradient makes no sense here");
	}

	core::Size half_bins = gradient_bins_/2;

	core::Real red_step_low = (zero_color_.x()-low_color_.x())/static_cast<core::Real>(half_bins);
	core::Real blue_step_low = (zero_color_.y()-low_color_.y())/static_cast<core::Real>(half_bins);
	core::Real green_step_low = (zero_color_.z()-low_color_.z())/static_cast<core::Real>(half_bins);

	core::Real red_step_high = (high_color_.x()-zero_color_.x())/static_cast<core::Real>(half_bins);
	core::Real blue_step_high = (high_color_.y()-zero_color_.y())/static_cast<core::Real>(half_bins);
	core::Real green_step_high = (high_color_.z()-zero_color_.z())/static_cast<core::Real>(half_bins);

	core::Real low_value_step = (-grid_min)/gradient_bins_;
	core::Real high_value_step = (grid_max)/gradient_bins_;

	core::Vector current_color = low_color_;
	core::Real current_low = grid_min;
	core::Real current_high = current_low+low_value_step;
	//set up colors between low_color and zero
	for ( core::Size bin = 1; bin <= half_bins; ++bin ) {
		current_color.x(current_color.x()+red_step_low);
		current_color.y(current_color.y()+green_step_low);
		current_color.z(current_color.z()+blue_step_low);
		std::string color_bin_string(utility::to_string<core::Size>(bin));
		std::string color_name(grid_name_+"_"+color_bin_string+"_color");
		ColorGradient new_color(current_color,current_low,current_high,color_name);
		color_data_.push_back(new_color);
		current_low = current_high;
		current_high = current_low+low_value_step;
	}


	//set up colors between zero and high color
	current_color = zero_color_;
	current_low = 0.0;
	current_high = current_low+high_value_step;

	for ( core::Size bin = 1; bin <= half_bins; ++bin ) {
		current_color.x(current_color.x()+red_step_high);
		current_color.y(current_color.y()+green_step_high);
		current_color.z(current_color.z()+blue_step_high);
		std::string color_bin_string(utility::to_string<core::Size>(bin+half_bins));
		ColorGradient new_color(current_color,current_low,current_high,grid_name_+"_"+color_bin_string+"_color");
		color_data_.push_back(new_color);
		current_low = current_high;
		current_high = current_low+high_value_step;
	}

}

void RenderGridsToKinemage::write_points(utility::io::ozstream & kin_file, scoring_grid::SingleGrid const & grid)
{
	//setup grid master
	std::string master_name = grid_name_;
	kin_file << "@master {" << master_name << "}" <<std::endl;
	//setup dotlist
	kin_file << "@dotlist {" << master_name << " dots} color=white master={" <<master_name << "}" <<std::endl;

	//loop through each color in color_data
	for ( core::Size color_index = 1; color_index <= color_data_.size(); ++color_index ) {
		//write dot for each point in color_range
		ColorGradient current_color(color_data_[color_index]);
		std::list<std::pair<core::Vector, core::Real> > point_list(
			grid.get_point_value_list_within_range(current_color.lower_bound,current_color.upper_bound,stride_));

		auto point_iterator = point_list.begin();
		for ( ; point_iterator != point_list.end(); ++point_iterator ) {
			core::Vector coords(point_iterator->first);
			core::Real value(point_iterator->second);

			//don't output zero points
			if ( value != 0.0 ) {
				kin_file <<"{" << master_name << "_points} " << current_color.color_name << " " <<
					coords.x() << " " << coords.y() << " " << coords.z() <<std::endl;
			}
		}
	}

}

void RenderGridsToKinemage::write_colors(utility::io::ozstream & kin_file)
{
	for ( core::Size color_index = 1; color_index <= color_data_.size(); ++color_index ) {
		ColorGradient current_color(color_data_[color_index]);
		//we asked the user to specify in RGB, but the kinemage needs to be in HSV
		core::Vector hsv_color(numeric::rgb_to_hsv(current_color.color_value));
		kin_file << "@hsvcolor {" << current_color.color_name << "} " <<
			hsv_color.x() << " "<< hsv_color.y()*100 << " " << hsv_color.z()*100 << std::endl;
	}
}

void RenderGridsToKinemage::write_header(utility::io::ozstream & kin_file)
{
	kin_file << "@kinemage 1" <<std::endl;
	kin_file << "@group dominant {dots}" <<std::endl;
}

std::string RenderGridsToKinemage::get_name() const {
	return mover_name();
}

std::string RenderGridsToKinemage::mover_name() {
	return "RenderGridsToKinemage";
}

std::string real_regex_pattern() {
	return "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?";
}

void RenderGridsToKinemage::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaRestriction three_dim_real_vector;
	three_dim_real_vector.name( "three_dim_real_vector" );
	three_dim_real_vector.base_type( xs_string );
	three_dim_real_vector.add_restriction( xsr_pattern, real_regex_pattern() + "," + real_regex_pattern() + "," + real_regex_pattern() );
	xsd.add_top_level_element( three_dim_real_vector );
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "grid_name", xs_string, "The name of the grid to dump to the file" )
		+ XMLSchemaAttribute::required_attribute( "file_name", xs_string, "The name of the file to which to dump the given grid" )
		+ XMLSchemaAttribute( "low_color", "three_dim_real_vector", "RGB values for color to use for low value" )
		+ XMLSchemaAttribute( "zero_color", "three_dim_real_vector", "RGB values for color to use for zero" )
		+ XMLSchemaAttribute( "high_color", "three_dim_real_vector", "RGB values for color to use for high value" )
		+ XMLSchemaAttribute( "color", "three_dim_real_vector", "RGB values for color to use" )
		+ XMLSchemaAttribute::attribute_w_default( "gradient_bins", xsct_non_negative_integer, "Size of bins to use", "10" )
		+ XMLSchemaAttribute::attribute_w_default( "stride", xsct_non_negative_integer, "Separation between possible colors", "2" );

	scoring_grid::attributes_for_parse_grid_set_from_tag( attlist, "The Grid Set from which to take the grid." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Visualization tool for grids", attlist );
}

std::string RenderGridsToKinemageCreator::keyname() const {
	return RenderGridsToKinemage::mover_name();
}

protocols::moves::MoverOP
RenderGridsToKinemageCreator::create_mover() const {
	return utility::pointer::make_shared< RenderGridsToKinemage >();
}

void RenderGridsToKinemageCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RenderGridsToKinemage::provide_xml_schema( xsd );
}


}
}


