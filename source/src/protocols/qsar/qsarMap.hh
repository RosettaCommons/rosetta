// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/qsar/qsarMap.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_qsarMap_hh
#define INCLUDED_protocols_qsar_qsarMap_hh

#include <map>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/qsar/qsarMap.fwd.hh>
//#include <protocols/qsar/qsarTypeManager.fwd.hh>

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <utility/vector1.hh>
#include <iostream>


namespace protocols {
namespace qsar {

class qsarPoint : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~qsarPoint() override;
	qsarPoint(std::string type, core::Real value, std::string name, core::conformation::ResidueOP residue );
	/// @brief set the value of the qsar point
	void set_value(core::Real value);


	/// @brief return the value of the qsar point
	core::Real get_value();

	/// @brief return the type of the qsar point
	std::string get_type();

	/// @brief return the name of the atom the qsar point is associated with
	std::string get_name();

	/// @brief return a pointer to the residue the qsar point is associated with
	core::conformation::ResidueOP get_residue();

private:
	std::string type_;
	core::Real value_;
	std::string atom_name_;
	core::conformation::ResidueOP residue_;
};


class qsarMap : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~qsarMap() override;
	qsarMap(std::string map_name, core::conformation::ResidueOP residue);

	/// @brief get size of map
	core::Size size();

	/// @brief initialize grid so that every point has a constant value for every atom and every qsarType. mostly for debugging
	void fill_with_value(core::Size value,utility::vector1<std::string> grids_to_use);

	// MolData no longer exists - if you need this, look into the string properties of the restype
	// brief initialize grid using data from the mol_data object attached to a residue (if availible)
	// bool fill_from_mol_data(core::chemical::sdf::MolData mol_data);

	/// @brief add a new qsar point to the map
	void add_point(std::string point_name, qsarPointOP new_point);

	/// @brief clear the qsar map
	void clear();

	/// @brief return a point in the map from the point name
	qsarPointOP get_point(std::string const & point_name);

	/// @brief return a point in the map from the atom_id and qsarType
	qsarPointOP get_point(core::Size const atom_id, std::string const & type);

	/// @brief return the residue associated with the qsarMap
	core::conformation::ResidueOP get_residue();

	/// @brief return a vector of points associated with a given atom_id
	utility::vector1<qsarPointOP> find_points_for_atom(core::Size const atom_id);

	/// @brief return a vector of all points associated with a given qsarType;
	utility::vector1<qsarPointOP> find_points_of_type(std::string const & type);


private:

	std::string map_name_;
	std::map<std::string,qsarPointOP> qsar_map_;
	std::multimap<std::string, qsarPointOP> type_map_;
	std::multimap<core::Size,qsarPointOP> atom_map_;

	core::conformation::ResidueOP residue_;

};

}
}


#endif
