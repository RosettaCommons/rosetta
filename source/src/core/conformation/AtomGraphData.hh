// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/conformation/AtomGraphData.hh
/// @author Sam DeLuca

#ifndef INCLUDED_core_conformation_AtomGraphData_hh
#define INCLUDED_core_conformation_AtomGraphData_hh

#include <core/conformation/PointGraphData.hh>

namespace core {
namespace conformation {

class AtomGraphVertexData : public PointGraphVertexData
{
public:
	AtomGraphVertexData() : PointGraphVertexData(), residue_id_(0),atom_radius_squared_(0.0) {}
	AtomGraphVertexData( numeric::xyzVector<core::Real> const & coors,std::string atom_name, core::Size residue_id)
	: PointGraphVertexData(coors),atom_name_(atom_name),residue_id_(residue_id) {}

	std::string & atom_name() {return atom_name_;}
	std::string const & atom_name() const {return atom_name_;}

	core::Size & residue_id() {return residue_id_;}
	core::Size const & residue_id() const {return residue_id_;}

	core::Real & atom_radius_squared() {return atom_radius_squared_;}
	core::Real const & atom_radius_squared() const {return atom_radius_squared_;}

private:
	std::string atom_name_;
	core::Size residue_id_;
	core::Real atom_radius_squared_;

};


class AtomGraphEdgeData : public PointGraphEdgeData
{
public:
	AtomGraphEdgeData() : PointGraphEdgeData() {}

	/// @brief inputs and outputs are distances squared
	AtomGraphEdgeData(platform::Real d2) : PointGraphEdgeData(d2) {}
};

}
}

#endif /* ATOMGRAPHDATA_HH_ */
