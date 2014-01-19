// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/RenderGridsToKinemage.cc
/// @author Sam DeLuca


#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/RenderGridsToKinemage.hh>
#include <protocols/qsar/RenderGridsToKinemageCreator.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>

#include <core/pose/Pose.hh>
#include <protocols/qsar/scoring_grid/GridBase.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/tag/Tag.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/color_util.hh>

#include <utility/excn/Exceptions.hh>
#include <iostream>





namespace protocols {
namespace qsar {

std::string RenderGridsToKinemageCreator::keyname() const
{
	return RenderGridsToKinemageCreator::mover_name();
}


moves::MoverOP RenderGridsToKinemageCreator::create_mover() const
{
	return new RenderGridsToKinemage;
}

std::string RenderGridsToKinemageCreator::mover_name()
{
	return "RenderGridsToKinemage";
}

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
		protocols::moves::Mover("RenderGridsToKinemage"),
		color_mode_(0),
		gradient_bins_(10),
		stride_(1)
		//low_cut_(0.0),
		//high_cut_(0.0)
{

}

RenderGridsToKinemage::RenderGridsToKinemage(
		RenderGridsToKinemage const & mover
		) :
		protocols::moves::Mover(mover),
		filename_(mover.filename_),
		color_mode_(mover.color_mode_),
		gradient_bins_(mover.gradient_bins_),
		stride_(mover.stride_),
		grid_name_(mover.grid_name_),
		color_(mover.color_),
		low_color_(mover.low_color_),
		zero_color_(mover.zero_color_),
		high_color_(mover.high_color_),
		color_data_(mover.color_data_),
		grid_(mover.grid_)
{

}

RenderGridsToKinemage::~RenderGridsToKinemage()
{
	//
}

protocols::moves::MoverOP RenderGridsToKinemage::clone() const
{
	return new RenderGridsToKinemage(*this);
}

std::string RenderGridsToKinemage::get_name() const
{
	return "RenderGridsToKinemage";
}


void RenderGridsToKinemage::apply(core::pose::Pose & )
{

	setup_colors();

	//we can write multiple grids to one file, peek at the first line and see if the header has been written
	utility::io::izstream infile;
	infile.open(filename_);
	std::string first_word;
	infile >> first_word;
	infile.close();

	utility::io::ozstream outfile;
	if(first_word == "@kinemage")
	{
		//header already written, append to existing file, don't write the header
		outfile.open(filename_,std::ios::app);
	}else
	{
		//header not written, open for overwrite, write the header
		outfile.open(filename_,std::ios::out);
		write_header(outfile);
	}
	write_colors(outfile);
	write_points(outfile);
	outfile.close();

}

void RenderGridsToKinemage::parse_my_tag(utility::tag::TagCOP const tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
)
{

	if(!tag->hasOption("grid_name"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("RenderGridsToKinemage requires the'grid name' option");
	}

	grid_name_ = tag->getOption<std::string>("grid_name");
	grid_ = utility::pointer::dynamic_pointer_cast<scoring_grid::SingleGrid,scoring_grid::GridBase>(scoring_grid::GridManager::get_instance()->get_grid(grid_name_) );
	if( ! grid_ )
	{
		throw utility::excn::EXCN_RosettaScriptsOption("RenderGridsToKinemage is currently unable to output the contents of a metagrid.  Sorry.");
	}

	if(!tag->hasOption("file_name"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("RenderGridsToKinemage requires the 'file_name' option");
	}
	filename_ = tag->getOption<std::string>("file_name");

	if(!tag->hasOption("low_color") && !tag->hasOption("zero_color") && !tag->hasOption("high_color"))
	{
		if(tag->hasOption("color"))
		{
			//one color
			color_mode_ = 1;
			color_ = numeric::comma_seperated_string_to_xyz<core::Real>(tag->getOption<std::string>("color"));

		}else
		{
			throw utility::excn::EXCN_RosettaScriptsOption("To use RenderGridsToKinemage you must specify color for a 1 color grid plot, "
				"high_color and low_color for a two color gradient, or high_color, low_color and "
				"zero_color for a three color gradient centered at zero");
		}
	}else if(tag->hasOption("low_color") && tag->hasOption("zero_color") && tag->hasOption("high_color"))
	{
		//3 part gradient
		color_mode_ = 3;
		low_color_ = numeric::comma_seperated_string_to_xyz<core::Real>(tag->getOption<std::string>("low_color"));
		zero_color_ = numeric::comma_seperated_string_to_xyz<core::Real>(tag->getOption<std::string>("zero_color"));
		high_color_ = numeric::comma_seperated_string_to_xyz<core::Real>(tag->getOption<std::string>("high_color"));

	}else if(tag->hasOption("low_color")&& tag->hasOption("high_color"))
	{
		//2 part gradient
		color_mode_ = 2;
		low_color_ = numeric::comma_seperated_string_to_xyz<core::Real>(tag->getOption<std::string>("low_color"));
		high_color_ = numeric::comma_seperated_string_to_xyz<core::Real>(tag->getOption<std::string>("high_color"));
	}else
	{
		throw utility::excn::EXCN_RosettaScriptsOption("To use RenderGridsToKinemage you must specify color for a 1 color grid plot, "
			"high_color and low_color for a two color gradient, or high_color, low_color and "
			"zero_color for a three color gradient centered at zero");
	}

	gradient_bins_ = tag->getOption<core::Size>("gradient_bins",10);
	stride_ = tag->getOption<core::Size>("stride",2);
	std::cout << stride_ <<std::endl;

}


void RenderGridsToKinemage::setup_colors()
{
	if(color_mode_ == 1)
	{
		setup_one_color_scheme();
	}else if(color_mode_ == 2)
	{
		setup_two_color_scheme();
	}else if(color_mode_ == 3)
	{
		setup_three_color_scheme();
	}else
	{
		assert(false); //This should never happen
	}

}

void RenderGridsToKinemage::setup_one_color_scheme()
{
	core::Real grid_min =grid_->get_min_value();
	core::Real grid_max = grid_->get_max_value();

	ColorGradient new_color(color_,grid_min,grid_max,grid_name_+"_color");
	color_data_.push_back(new_color);
}

void RenderGridsToKinemage::setup_two_color_scheme()
{
	core::Real grid_min = grid_->get_min_value();
	core::Real grid_max = grid_->get_max_value();

	core::Real red_step = (high_color_.x() - low_color_.x())/static_cast<core::Real>(gradient_bins_);
	core::Real green_step = (high_color_.y() - low_color_.y())/static_cast<core::Real>(gradient_bins_);
	core::Real blue_step = (high_color_.z() - low_color_.z())/static_cast<core::Real>(gradient_bins_);

	core::Real value_step = (grid_max - grid_min)/gradient_bins_;


	core::Vector current_color = low_color_;
	core::Real current_low = grid_min;
	core::Real current_high = current_low+value_step;
	for(core::Size bin = 1; bin <= gradient_bins_; ++bin)
	{
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

void RenderGridsToKinemage::setup_three_color_scheme()
{
	core::Real grid_min = grid_->get_min_value();
	core::Real grid_max = grid_->get_max_value();
	if(grid_min >= 0)
	{
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
	for(core::Size bin = 1; bin <= half_bins; ++bin)
	{
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

	for(core::Size bin = 1; bin <= half_bins;++bin)
	{
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

void RenderGridsToKinemage::write_points(utility::io::ozstream & kin_file)
{
	//setup grid master
	std::string master_name = grid_name_;
	kin_file << "@master {" << master_name << "}" <<std::endl;
	//setup dotlist
	kin_file << "@dotlist {" << master_name << " dots} color=white master={" <<master_name << "}" <<std::endl;

	//loop through each color in color_data
	for(core::Size color_index = 1; color_index <= color_data_.size();++color_index)
	{
		//write dot for each point in color_range
		ColorGradient current_color(color_data_[color_index]);
		std::list<std::pair<core::Vector, core::Real> > point_list(
				grid_->get_point_value_list_within_range(current_color.lower_bound,current_color.upper_bound,stride_));

		std::list<std::pair<core::Vector, core::Real> >::iterator point_iterator = point_list.begin();
		for(; point_iterator != point_list.end(); ++point_iterator)
		{
			core::Vector coords(point_iterator->first);
			core::Real value(point_iterator->second);

			//don't output zero points
			if( value != 0.0)
			{
				kin_file <<"{" << master_name << "_points} " << current_color.color_name << " " <<
					coords.x() << " " << coords.y() << " " << coords.z() <<std::endl;
			}
		}
	}

}

void RenderGridsToKinemage::write_colors(utility::io::ozstream & kin_file)
{
	for(core::Size color_index = 1; color_index <= color_data_.size();++color_index)
	{
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

}
}


