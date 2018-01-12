// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/PoseSelectorFactory.hh
/// @brief
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_moves_PoseSelectorFactory_hh
#define INCLUDED_protocols_moves_PoseSelectorFactory_hh

// Unit Headers
#include <protocols/rosetta_scripts/PoseSelectorFactory.fwd.hh>
#include <protocols/rosetta_scripts/PoseSelector.fwd.hh>
#include <protocols/rosetta_scripts/PoseSelector.hh>
#include <protocols/rosetta_scripts/PoseSelectorCreator.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
// c++ headers
#include <map>
#include <set>

namespace protocols {
namespace rosetta_scripts {

/// @brief This templated class will register an instance of an
/// PoseSelectorCreator (class T) with the PoseSelectorFactory.  It will ensure
/// that no PoseSelectorCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class PoseSelectorRegistrator : public utility::factory::WidgetRegistrator< PoseSelectorFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< PoseSelectorFactory, T > parent;
public:
	PoseSelectorRegistrator() : parent() {}
};


class PoseSelectorFactory : public utility::SingletonBase< PoseSelectorFactory >
{
public:
	friend class utility::SingletonBase< PoseSelectorFactory >;

	typedef std::map< std::string, PoseSelectorCreatorOP > PoseSelectorMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;
	typedef core::pose::Pose Pose;

public:
	virtual ~PoseSelectorFactory();

	void factory_register( PoseSelectorCreatorOP creator );

	/// @brief Create a PoseSelector given its identifying string
	PoseSelectorOP newPoseSelector( std::string const & );

	/// @brief return new PoseSelector by Tag parsing; the identifying string for the PoseSelector is in the Tag
	PoseSelectorOP
	newPoseSelector(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	void define_pose_selector_group( utility::tag::XMLSchemaDefinition & xsd ) const;
	static std::string pose_selector_group_name();
	static std::string
	complex_type_name_for_pose_selector( std::string const & selector_name );


private:
	PoseSelectorFactory();

	// Unimplemented -- uncopyable
	PoseSelectorFactory( PoseSelectorFactory const & );
	PoseSelectorFactory const & operator = ( PoseSelectorFactory const & );


private:

	PoseSelectorMap poseselector_creator_map_;
};

} //namespace rosetta_scripts
} //namespace protocols

#endif
