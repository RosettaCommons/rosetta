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

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// c++ headers
#include <map>
#include <set>
#include <utility/vector1.hh>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <thread>
#endif
#endif

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

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif
#endif

private:

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static TempInterpolatorFactory * create_singleton_instance();

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
