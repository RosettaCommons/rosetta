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

#ifndef INCLUDED_protocols_toolbox_IGEdgeReweighters_hh
#define INCLUDED_protocols_toolbox_IGEdgeReweighters_hh

// Unit headers
#include <core/pack/task/IGEdgeReweightContainer.hh>
// Package headers
//#include <core/pose/Pose.fwd.hh>
//#include <core/pack/task/PackerTask.fwd.hh>
#include <set>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {

class IGLigandDesignEdgeUpweighter : public core::pack::task::IGEdgeReweighter{

public:
	IGLigandDesignEdgeUpweighter( core::Real weight_factor ){ weight_factor_ = weight_factor; }

	core::Real get_edge_reweight(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask const & task,
		Size res1,
		Size res2
	) const;

private:
	core::Real weight_factor_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	IGLigandDesignEdgeUpweighter();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


//template <class T>
class ResidueGroupIGEdgeUpweighter : public core::pack::task::IGEdgeReweighter{

	//typedef typename T::const_iterator t_it;

public:
	/*
	ResidueGroupIGEdgeUpweighter(
	core::Real weight_factor,
	T const & group1_in,
	T const & group2_in
	) : weight_factor_(weight_factor)
	{
	for( t_it g1it = group1_in.begin(); g1it != group1_in.end(); ++g1it ) group1_.insert( *g1it );
	for( t_it g2it = group2_in.begin(); g2it != group2_in.end(); ++g2it ) group2_.insert( *g2it );
	}
	*/
	// /*
	ResidueGroupIGEdgeUpweighter(
		core::Real weight_factor,
		utility::vector1< core::Size> const & group1_in,
		utility::vector1< core::Size> const & group2_in
	) : weight_factor_(weight_factor)
	{
		for ( utility::vector1< core::Size>::const_iterator g1it = group1_in.begin(); g1it != group1_in.end(); ++g1it ) group1_.insert( *g1it );
		for ( utility::vector1< core::Size>::const_iterator g2it = group2_in.begin(); g2it != group2_in.end(); ++g2it ) group2_.insert( *g2it );
	}
	// */

	core::Real get_edge_reweight(
		core::pose::Pose const &, // pose,
		core::pack::task::PackerTask const &, // task,
		Size res1,
		Size res2
	) const {
		if ( ( (group1_.find(res1) != group1_.end()) && (group2_.find(res2) != group2_.end() ) )
				||( (group2_.find(res1) != group2_.end()) && (group1_.find(res2) != group1_.end() ) ) ) {
			return weight_factor_;
		} else return default_weight_;
	}

private:

	core::Real weight_factor_;
	std::set< core::Size > group1_;
	std::set< core::Size > group2_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ResidueGroupIGEdgeUpweighter();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //namespace toolbox
} //namespace protocols


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_IGEdgeReweighters )
#endif // SERIALIZATION


#endif
