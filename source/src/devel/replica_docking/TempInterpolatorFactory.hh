// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Zhe Zhang

#ifndef INCLUDED_devel_replica_docking_TempInterpolatorFactory_hh
#define INCLUDED_devel_replica_docking_TempInterpolatorFactory_hh

// Unit Headers
#include <devel/replica_docking/TempInterpolatorFactory.fwd.hh>
#include <devel/replica_docking/TempInterpolator.fwd.hh>

// #include <protocols/moves/MoverCreator.hh>
// #include <protocols/moves/Mover.fwd.hh>
// #include <basic/datacache/DataMap.fwd.hh>
// #include <protocols/filters/Filter.fwd.hh> // Filters_map (typedef)

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
// #include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// c++ headers
#include <map>
#include <set>
#include <utility/vector1.hh>


namespace devel {
namespace replica_docking {


class TempInterpolatorFactory : public utility::pointer::ReferenceCount
{
public:
	//	typedef std::map< std::string, MoverCreatorOP > MoverMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;
	typedef core::pose::Pose Pose;

public:
	virtual ~TempInterpolatorFactory();

	static
	TempInterpolatorFactory * get_instance();

	TempInterpolatorBaseOP new_tempInterpolator( utility::tag::TagCOP const tag, core::Size n_level );

private:
	TempInterpolatorFactory();

	// Unimplemented -- uncopyable
	TempInterpolatorFactory( TempInterpolatorFactory const & );
	TempInterpolatorFactory const & operator = ( TempInterpolatorFactory const & );

private:
	static TempInterpolatorFactory * instance_;

};

} //namespace
} //namespace

#endif
