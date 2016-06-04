// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/qsar/qsarMap.cc
/// @author Sam DeLuca

#include <protocols/qsar/qsarMap.hh>
//#include <protocols/qsar/qsarTypeManager.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace qsar {

/// @details Auto-generated virtual destructor
qsarMap::~qsarMap() {}

/// @details Auto-generated virtual destructor
qsarPoint::~qsarPoint() {}

static THREAD_LOCAL basic::Tracer qsarMapTracer( "protocols.qsar.qsarMap" );


qsarPoint::qsarPoint(std::string type, core::Real value, std::string name, core::conformation::ResidueOP residue ):
	type_(type), value_(value), atom_name_(name), residue_(residue)
{}

void qsarPoint::set_value(core::Real value)
{
	value_ = value;
}

core::Real qsarPoint::get_value()
{
	return value_;
}

std::string qsarPoint::get_type()
{
	return type_;
}

std::string qsarPoint::get_name()
{
	return atom_name_;
}

core::conformation::ResidueOP qsarPoint::get_residue()
{
	return residue_;
}

qsarMap::qsarMap(std::string map_name, core::conformation::ResidueOP residue) :
	map_name_(map_name), residue_(residue)
{ }

core::Size qsarMap::size()
{
	return qsar_map_.size();
}

void qsarMap::fill_with_value(core::Size value,utility::vector1<std::string> grids_to_use)
{
	this->clear();
	core::Size const heavy_atom_size(residue_->nheavyatoms());
	//std::cout << "atom size" << heavy_atom_size <<std::endl;

	for ( core::Size atom_index = 1; atom_index <= heavy_atom_size; ++atom_index ) {
		std::string atom_name(residue_->atom_name(atom_index));
		for ( core::Size grid_index = 1; grid_index <= grids_to_use.size(); ++grid_index ) {
			//qsarType current_type(static_cast<qsarType>(enum_index));
			qsarPointOP current_point( new qsarPoint(grids_to_use[grid_index],value,atom_name,residue_) );
			std::string point_name(grids_to_use[grid_index]+"_"+atom_name);
			this->add_point(point_name,current_point);
		}
	}

	//std::cout << "filled with value"<<std::endl;
}

// MolData no longer exists - if you need this, look into the string properties of the restype
//bool qsarMap::fill_from_mol_data(core::chemical::sdf::MolData mol_data)
//{
// //the mol data line has type qsar_map.  it is a collection of lines like this:
// //atomno\tqsar_type\tqsar_weight
// //atomno is a Size, qsar_type is a string, qsar_weight is a float
// this->clear();
// utility::vector1<std::string> data_lines(mol_data.get_mol_data_string_vector("qsar_map", '\n'));
// if(data_lines.size() == 0)
// {
//  return false;
// }
//
// for(core::Size index = 1; index <= data_lines.size(); ++index)
// {
//  std::string line(data_lines[index]);
//
//  utility::vector1<std::string> line_fields(utility::string_split(line, '\t'));
//  if(line_fields.size() != 3)
//  {
//   utility_exit_with_message("a qsar_map data line doesn't have 3 fields. something is wrong. Aborting.");
//  }
//
//  std::string atom_name(residue_->atom_name(utility::string2int(line_fields[1])));
//  std::string qsar_type(line_fields[2]);
//  //qsarType qsar_type(qsarTypeManager::qsar_type_from_name(line_fields[2]));
//  core::Real qsar_weight(utility::string2float(line_fields[3]));
//  qsarPointOP current_point(new qsarPoint(qsar_type,qsar_weight,atom_name,residue_));
//  std::string point_name(line_fields[3]+"_"+atom_name);
//  this->add_point(point_name,current_point);
// }
// return true;
//}

void qsarMap::add_point(std::string point_name, qsarPointOP new_point)
{
	std::pair<std::string, qsarPointOP> new_qsar_map_item(point_name,new_point);

	std::string new_point_type(new_point->get_type());
	//qsarType new_point_type(new_point->get_type());
	std::pair<std::string,qsarPointOP> new_type_map_item(new_point_type,new_point);

	core::conformation::ResidueOP point_residue(new_point->get_residue());
	std::string atom_name(new_point->get_name());
	core::Size atom_index(point_residue->atom_index(atom_name));

	std::pair<core::Size, qsarPointOP> new_atom_map_item(atom_index,new_point);

	qsar_map_.insert(new_qsar_map_item);
	type_map_.insert(new_type_map_item);
	atom_map_.insert(new_atom_map_item);
}

void qsarMap::clear()
{
	qsar_map_.clear();
	type_map_.clear();
	atom_map_.clear();
}

qsarPointOP qsarMap::get_point(std::string const & point_name)
{
	std::map<std::string,qsarPointOP>::iterator point(qsar_map_.find(point_name));
	if ( point == qsar_map_.end() ) {
		return 0;
	} else {
		return point->second;
	}
}

qsarPointOP qsarMap::get_point(core::Size const atom_id, std::string const & type)
{
	std::multimap<core::Size, qsarPointOP>::iterator lower_bound(atom_map_.lower_bound(atom_id));
	std::multimap<core::Size, qsarPointOP>::iterator upper_bound(atom_map_.upper_bound(atom_id));

	std::multimap<core::Size,qsarPointOP>::iterator current_point(lower_bound);
	for ( ; current_point !=upper_bound; ++current_point ) {
		qsarPointOP point = current_point->second;
		if ( point->get_type() == type ) {
			return point;
		}
	}
	qsarMapTracer << "couldn't find a point, returning 0" <<std::endl;
	return 0;

}

core::conformation::ResidueOP qsarMap::get_residue()
{
	return residue_;
}

utility::vector1<qsarPointOP> qsarMap::find_points_for_atom(core::Size const atom_id)
{
	utility::vector1<qsarPointOP> points;

	std::multimap<core::Size, qsarPointOP>::iterator lower_bound(atom_map_.lower_bound(atom_id));
	std::multimap<core::Size, qsarPointOP>::iterator upper_bound(atom_map_.upper_bound(atom_id));

	std::multimap<core::Size,qsarPointOP>::iterator current_point(lower_bound);
	for ( ; current_point != upper_bound; ++current_point ) {
		points.push_back(current_point->second);
	}

	return points;
}


utility::vector1<qsarPointOP> qsarMap::find_points_of_type(std::string const & type)
{
	utility::vector1<qsarPointOP> points;

	std::multimap<std::string, qsarPointOP>::iterator lower_bound(type_map_.lower_bound(type));
	std::multimap<std::string, qsarPointOP>::iterator upper_bound(type_map_.upper_bound(type));

	std::multimap<std::string,qsarPointOP>::iterator current_point(lower_bound);
	for ( ; current_point != upper_bound; ++current_point ) {
		points.push_back(current_point->second);
	}

	return points;
}

}
}
