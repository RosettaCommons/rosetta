// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/ElectrostaticMetaGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/ElectroStaticMetaGrid.hh>
#include <protocols/qsar/scoring_grid/ElectroStaticMetaGridCreator.hh>
#include <core/conformation/Residue.hh>

#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

#include <numeric/util.hh>


#include <map>

namespace protocols {
namespace qsar {
namespace scoring_grid {

static basic::Tracer TR("protocols.qsar.scoring_grid.ElectroStaticMetaGrid");

std::string ElectroStaticMetaGridCreator::keyname()const
{
	return ElectroStaticMetaGridCreator::grid_name();
}

GridBaseOP ElectroStaticMetaGridCreator::create_grid(utility::tag::TagPtr const tag) const
{
	GridBaseOP electrostatic_grid= new ElectroStaticMetaGrid();

	electrostatic_grid->parse_my_tag(tag);

	return electrostatic_grid;
}

std::string ElectroStaticMetaGridCreator::grid_name()
{
	return "ElectroStaticMetaGrid";
}


ElectroStaticMetaGrid::ElectroStaticMetaGrid() : type_("ElectroStaticMetaGrid"),chain_('A'),weight_(1.0)
{
	core::Real charges[] = {
		-0.760,-0.750,-0.660,-0.620,-0.610,-0.550,
		-0.530,-0.510,-0.470,-0.370,-0.270,-0.250,
		-0.180,-0.160,-0.130,-0.115,-0.100,-0.090,
		 0.000, 0.070, 0.095, 0.115, 0.310, 0.430,
		 0.510, 0.550, 0.620, 1.000, 1.500, 2.000, 3.000};

	//When we start supporting c++11 and I can use initializer lists for vectors im throwing
	//a party and you're all invited -sld
	charges_ = utility::vector1<core::Real>(charges,charges+sizeof(charges)/sizeof(core::Real));

}

void ElectroStaticMetaGrid::initialize(core::Vector const & center, core::Real width, core::Real resolution)
{
	//A list of 31 charges pulled from python/apps/public/molfile_to_params.py
	//Most charges that you find should be close to one of these.

	TR << "initializing electrostatic grid" << std::endl;

	//make grids for all the charges and initialize them
	for(int charge_index = 1; charge_index <= charges_.size(); ++charge_index)
	{
		core::Real current_charge = charges_[charge_index];
		add_new_charge_grid(current_charge);
		charge_grid_map_[current_charge].initialize(center,width,resolution);
	}
}


void ElectroStaticMetaGrid::refresh(core::pose::Pose const & pose, core::Vector const &  center )
{
	std::map<core::Real, ChargeGrid>::iterator charge_grid_it(charge_grid_map_.begin());
	TR << "refreshing electrostatic grid.  This could take a while"  <<std::endl;
	for(;charge_grid_it != charge_grid_map_.end();++charge_grid_it)
	{
		TR << "refreshing electrostatic grid with charge: " <<charge_grid_it->first <<std::endl;
		charge_grid_it->second.set_chain(chain_);
		charge_grid_it->second.refresh(pose,center);
	}

}

void ElectroStaticMetaGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void ElectroStaticMetaGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}

void ElectroStaticMetaGrid::parse_my_tag(utility::tag::TagPtr const tag)
{
	if (!tag->hasOption("weight")){
		utility_exit_with_message("Could not make AtrGrid: you must specify a weight when making a new grid");
	}
	weight_= tag->getOption<core::Real>("weight");
}


core::Real ElectroStaticMetaGrid::score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP qsar_map)
{
	core::Real score = 0.0;
	for(core::Size atom_index = 1; atom_index <= residue.natoms(); ++atom_index)
	{
		core::Real charge = residue.atomic_charge(atom_index);
		core::Real nearest_charge_in_map = numeric::find_nearest_value<core::Real>(charges_,charge);
		score += charge_grid_map_[nearest_charge_in_map].get_point(residue.xyz(atom_index));

	}

}

std::string ElectroStaticMetaGrid::get_type()
{
	return type_;
}

void ElectroStaticMetaGrid::set_chain(char chain)
{
	chain_ = chain;
}


void ElectroStaticMetaGrid::add_new_charge_grid(core::Real const & charge)
{
	charge_grid_map_.insert(std::make_pair(charge,ChargeGrid(charge,1.0)));
}

void ElectroStaticMetaGrid::dump_BRIX(std::string const & prefix)
{
	utility_exit_with_message("ElectroStaticMetaGrids are incapable of outputting BRIX grids right now sorry.");
}

}
}
}
