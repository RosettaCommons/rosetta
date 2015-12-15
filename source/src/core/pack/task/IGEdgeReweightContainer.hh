// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file edge reweighting for Interaction Graphs
/// @brief
/// @author Florian Richter, floric@u.washington.edu, june 08

#ifndef INCLUDED_core_pack_task_IGEdgeReweightContainer_hh
#define INCLUDED_core_pack_task_IGEdgeReweightContainer_hh

// Unit headers
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
// Package headers

// Project headers
#include <core/types.hh>
//#include <core/conformation/Residue.fwd.hh>
//#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
//#include <core/id/AtomID.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


// STL Headers

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace task {


/// @brief helper class for IGEdgeReweightContainer
class IGEdgeReweighter : public utility::pointer::ReferenceCount {

public:
	IGEdgeReweighter() : default_weight_(1.0) {}
	virtual ~IGEdgeReweighter();

	virtual
	Real
	get_edge_reweight(
		pose::Pose const & pose,
		PackerTask const & task,
		Size res1,
		Size res2
	) const = 0;

protected:

	//this is the default non-upweighted weight. all reweighters derived from this class
	//should return this variable for edges with no specific upweighting
	core::Real const default_weight_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief class that interfaces to containers holding IGedge weights between individual residues of the task
/// @brief note: this class only knows about the pose numbering, not about the numbers in the IG
class IGEdgeReweightContainer : public utility::pointer::ReferenceCount {

public:

	IGEdgeReweightContainer(Size nres );
	virtual ~IGEdgeReweightContainer();

	Real res_res_weight(
		pose::Pose const & pose,
		PackerTask const & task,
		Size res1id,
		Size res2id
	) const;

	void add_reweighter( IGEdgeReweighterOP reweighter );

	utility::vector1< IGEdgeReweighterOP >::const_iterator
	reweighters_begin() const {
		return edge_reweighters_.begin(); }

	utility::vector1< IGEdgeReweighterOP >::const_iterator
	reweighters_end() const {
		return edge_reweighters_.end(); }

private:

	utility::vector1< IGEdgeReweighterOP > edge_reweighters_;
	Size nres_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	IGEdgeReweightContainer();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //task
} //pack
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_IGEdgeReweightContainer )
#endif // SERIALIZATION


#endif
