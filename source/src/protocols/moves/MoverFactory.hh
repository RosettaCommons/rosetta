// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/MoverFactory.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_moves_MoverFactory_hh
#define INCLUDED_protocols_moves_MoverFactory_hh

// Unit Headers
#include <protocols/moves/MoverFactory.fwd.hh>
#include <protocols/moves/MoverCreator.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh> // Filters_map (typedef)

#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// c++ headers
#include <map>
#include <set>

namespace protocols {
namespace moves {

/// @brief This templated class will register an instance of an
/// MoverCreator (class T) with the MoverFactory.  It will ensure
/// that no MoverCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class MoverRegistrator : public utility::factory::WidgetRegistrator< MoverFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< MoverFactory, T > parent;
public:
	MoverRegistrator() : parent() {}
};


class MoverFactory : public utility::SingletonBase< MoverFactory >
{
public:
	friend class utility::SingletonBase< MoverFactory >;

	typedef std::map< std::string, MoverCreatorOP > MoverMap;
	typedef utility::tag::Tag Tag;
	typedef utility::tag::TagCOP TagCOP;
	typedef core::pose::Pose Pose;

public:
	virtual ~MoverFactory();

	void factory_register( MoverCreatorOP creator );

	/// @brief Create a mover given its identifying string
	MoverOP newMover( std::string const & );

	/// @brief return new Mover by Tag parsing; the identifying string for the Mover is in the Tag
	MoverOP
	newMover(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		Movers_map const &,
		Pose const &
	);

	/// @brief Read access to the set of all MoverCreators; for unit testing purposes
	MoverMap const & mover_creator_map() const;

	void define_mover_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
	static std::string mover_xml_schema_group_name();

private:
	MoverFactory();

	// Unimplemented -- uncopyable
	MoverFactory( MoverFactory const & ) = delete;
	MoverFactory const & operator = ( MoverFactory const & ) = delete;

private:

	MoverMap mover_creator_map_;

	std::set< std::string > forbidden_names_; //certain names have historic meanings and shouldn't be assigned to new movers

};

} //namespace moves
} //namespace protocols

#endif
