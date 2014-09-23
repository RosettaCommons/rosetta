// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/ShapeGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/ShapeGrid.hh>
#include <protocols/qsar/scoring_grid/ShapeGridCreator.hh>


#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>

#include <basic/database/open.hh>

#include <numeric/interpolation/util.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <numeric/kdtree/KDTree.hh>
#include <numeric/kdtree/nearest_neighbors.hh>

#include <utility/tag/Tag.hh>
#include <utility/io/izstream.hh>
#include <utility/tools/make_vector.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {
	
static thread_local basic::Tracer TR( "protocols.qsar.scoring_grid.ShapeGrid" );

std::string ShapeGridCreator::keyname() const
{
	return ShapeGridCreator::grid_name();
}

GridBaseOP ShapeGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP shape_grid( new ShapeGrid() );

	shape_grid->parse_my_tag(tag);

	return shape_grid;
}

GridBaseOP ShapeGridCreator::create_grid() const
{
	return GridBaseOP( new ShapeGrid() );
}


std::string ShapeGridCreator::grid_name()
{
	return "ShapeGrid";
}

ShapeGrid::ShapeGrid() : SingleGrid("ShapeGrid")
{
	distance_bin_width_ = 10.0/25.0;
	theta_bin_width_ = 2.0/25.0;
	phi_bin_width_ = 360.0/25.0;
	load_kbp_data();
}
	
ShapeGrid::~ShapeGrid()
{
	
}

	
void ShapeGrid::refresh(core::pose::Pose const & pose, core::Vector const & )
{
	this->fill_with_value(0.0);
	
	//First, loop through all the protein residues, create a KD tree node where coordinates are xyz coords, and data is pointer to residue.
	//This will be used to avoid having to traverse the entire protein for every point.
	core::Size chain_id = core::pose::get_chain_id_from_chain(get_chain(),pose);
	
	utility::vector1<utility::vector1<core::Real> > coord_list;
	utility::vector1<utility::pointer::ReferenceCountOP> data_list;
	for(core::Size residue_index = 1; residue_index <= pose.n_residue(); ++residue_index)
	{
		core::conformation::ResidueOP residue( new core::conformation::Residue(pose.residue(residue_index)) );
		if(residue->chain() == chain_id)
		{
			continue;
		}
		
		core::Vector coords(residue->nbr_atom_xyz());
		utility::vector1<core::Real> xyz_vector(3,0.0);
		xyz_vector[1] = coords.x();
		xyz_vector[2] = coords.y();
		xyz_vector[3] = coords.z();
		
		coord_list.push_back(xyz_vector);
		data_list.push_back(static_cast<utility::pointer::ReferenceCountOP>(residue));
	}
	
	numeric::kdtree::KDTree residue_tree(coord_list,data_list);
	
	//For every point in the grid, find the closest residues with centers within
	//10.0 A.  For each located residue, find score from KBP.
	//take the highest score from each residue and place in the grid.
    numeric::xyzVector<core::Size> dimensions = get_dimensions();
	for(core::Size x_index =0; x_index < dimensions.x(); ++x_index)
	{
		for(core::Size y_index = 0; y_index < dimensions.y(); ++y_index)
		{
			for(core::Size z_index = 0; z_index < dimensions.z(); ++z_index)
			{
				core::Vector query_coords (this->get_pdb_coords(x_index, y_index, z_index));
				utility::vector1<core::Real> kdcoords(3,0.0);
				kdcoords[1] = query_coords.x();
				kdcoords[2] = query_coords.y();
				kdcoords[3] = query_coords.z();
				numeric::kdtree::KDPointList nearest_residues(numeric::kdtree::nearest_neighbors(residue_tree, kdcoords, 10, 10.0));
				
				core::Real score_for_grid = get_point_score(nearest_residues, query_coords);
				this->set_point(query_coords, score_for_grid);
            }
        }
    }
}
	
void ShapeGrid::refresh(core::pose::Pose const & pose, core::Vector const & center,core::Size const & )
{
	refresh(pose, center);
}
	
void ShapeGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}
	
void ShapeGrid::parse_my_tag(utility::tag::TagCOP)
{
	
}
	
core::Real ShapeGrid::score(core::conformation::UltraLightResidue const & residue, core::Real const max_score, qsarMapOP)
{
	core::Real total_score = 0.0;
	for(core::Size atomno = 1; atomno <= residue.natoms();++atomno)
	{
		core::Vector coords(residue[atomno]);
		total_score += this->get_point(coords);
		if(total_score >= max_score)
		{
			break;
		}
	}
	return total_score;
}
	
core::Real ShapeGrid::atom_score(core::conformation::UltraLightResidue const & residue, core::Size atomno, qsarMapOP)
{
	return this->get_point(residue[atomno]);
}
	
core::Real ShapeGrid::score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP)
{
	core::Real total_score = 0.0;
	for(core::Size atomno = 1; atomno <= residue.natoms();++atomno)
	{
		core::Vector coords(residue.xyz(atomno));
		total_score += this->get_point(coords);
		if(total_score >= max_score)
		{
			break;
		}
	}
	return total_score;
}
	
core::Real ShapeGrid::atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapOP )
{
	return this->get_point(residue.xyz(atomno));
}

utility::json_spirit::Value ShapeGrid::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;
	
	//Serialize the base data for the grid.  we don't need to store the kbp because it is huge and doesnt ever change
	Pair base_data("base_data",SingleGrid::serialize());
	
	return Value(utility::tools::make_vector(base_data));
}

void ShapeGrid::deserialize(utility::json_spirit::mObject data)
{
	SingleGrid::deserialize(data["base_data"].get_obj());
}
	
core::Real ShapeGrid::get_point_score(numeric::kdtree::KDPointList const & nearest_residues,core::Vector const & query_coords)
{
	
	//For points that are out in empty space, set score to zero
	//We know the point is in space if it has no neighbors within 10 A.
	if(nearest_residues.size() == 0)
	{
		return 0.0;
	}
	
	core::Real total_score = 0.0;
	
	for(numeric::kdtree::KDPointList::const_iterator it = nearest_residues.begin();it != nearest_residues.end();++it)
	{
		numeric::kdtree::KDPointOP kd_point(*it);
		//This is a bit dangerous
		utility::pointer::ReferenceCountOP data_pointer(kd_point->data());
		core::conformation::ResidueOP residue(utility::pointer::dynamic_pointer_cast< core::conformation::Residue > ( data_pointer ));
		
		core::Vector ca(residue->xyz("CA"));
		core::Vector cb;
		if(residue->name3() == "GLY")
		{
			cb = residue->xyz("2HA");
		}else
		{
			cb = residue->xyz("CB");
		}
		core::Vector n(residue->xyz("N"));
		core::Vector ha;
		if(residue->name3() == "GLY")
		{
			ha = residue->xyz("1HA");
		}else
		{
			ha = residue->xyz("HA");
		}
		core::Real distance = cb.distance(query_coords);
		
		//if the distance > 8.0 A, the score will be 0.0
		if(distance > 8.0)
		{
			continue;
		}
		core::Real theta = cos(numeric::angle_radians(query_coords, cb, ca));
		core::Real phi = numeric::dihedral_degrees(ha, ca, cb, query_coords);
		std::string res_name3(residue->name3());
		core::Real score = get_score_from_angles(res_name3, distance, theta, phi);
		if(score < 0)
		{
			total_score += score;
		}
		//We don't want to penalize empty space, so a sigmoid smoothing function is applied based on the distance from CB.
		
		/*
		//if the distance < 4.0 A, the smoothing factor will be 1 and we can avoid computing an exp
		if(distance < 4.0)
		{
			total_score += score;
		}else
		{
			core::Real smooth_score = score*(-(1.0/(1.0+std::exp(-3.0*distance+18.0)))+1.0);
			total_score += smooth_score;
		}
		 */
	}
	return total_score/(core::Real)nearest_residues.size();
}

core::Real ShapeGrid::get_score_from_angles(std::string const & name3,core::Real distance, core::Real theta, core::Real phi)
{
	int distance_bin = static_cast<int>(std::floor(distance/distance_bin_width_));
	int theta_bin = static_cast<int>(std::floor((1.0+theta)/theta_bin_width_));
	int phi_bin = static_cast<int>(std::floor((180.0+phi)/phi_bin_width_));
	
	return kbp_data_[name3]->getValue(distance_bin, theta_bin, phi_bin);
	
}
	
void ShapeGrid::load_kbp_data()
{
	
	// The KBP data is stored as a list of lists in the form [name3,data_list].
	// Here, name3 is the 3letter abbreviation of the residue, data_list is a list in the form [distance,theta,phi]
	// the points are in array order.
	
	TR << "Loading Shape grid KBP data" << std::endl;
	std::string filename(basic::database::full_name("scoring/qsar/shape_histogram_data.js"));
	if(!utility::file::file_exists(filename))
    {
        utility_exit_with_message("cannot parse "+filename+" because it does not exist");
    }
    
	utility::io::izstream infile;
	infile.open(filename.c_str(),std::ifstream::in);
	utility::json_spirit::mValue kbp_data_list;
    if(!utility::json_spirit::read(infile,kbp_data_list))
    {
        infile.close();
        utility_exit_with_message("cannot parse JSON in file "+ filename);
    }
	
	infile.close();
	

	//The base and resolution don't really matter here because I'm just using this class as a convenient 3D array
	KBPGridOP empty_array( new core::grid::CartGrid<core::Real> );
	empty_array->setBase(0.0, 0.0, 0.0);
	empty_array->setDimensions(25,25,25,1.0,1.0,1.0);
	empty_array->setupZones();
	empty_array->zero();
	
	
	utility::json_spirit::mArray kbp_outer_array = kbp_data_list.get_array();
	for(utility::json_spirit::mArray::iterator res_data = kbp_outer_array.begin(); res_data != kbp_outer_array.end();++res_data)
	{
		utility::json_spirit::mArray current_res_data(res_data->get_array());
		std::string resname(current_res_data[0].get_str());
		
		KBPGridOP current_array( new core::grid::CartGrid<core::Real> );
		empty_array->clone(*current_array);
		
		
		utility::json_spirit::mArray current_res_list(current_res_data[1].get_array());
		utility::json_spirit::mArray::iterator res_list_it = current_res_list.begin();
		for(int distance_bin = 0; distance_bin < 25; ++distance_bin)
		{
			for(int theta_bin = 0; theta_bin < 25; ++theta_bin)
			{
				for(int phi_bin = 0; phi_bin < 25; ++phi_bin)
				{
					core::Real kbp_value = res_list_it->get_real();
					current_array->setValue(distance_bin, theta_bin, phi_bin, kbp_value);
					res_list_it++;
					
				}
			}
		}
		
		kbp_data_.insert(std::make_pair(resname,current_array));
	}
	TR << "Done loading Shape grid KBP data" << std::endl;
}
	
}
}
}
