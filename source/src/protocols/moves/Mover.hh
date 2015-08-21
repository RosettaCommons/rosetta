// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/Mover.hh
/// @brief: A mover is an object that can apply a conformational change to a pose.
///
/// @details: Each derived class should define its own apply() statement.
///   DO NOT JUST GO ADDING FUNCTIONS TO THIS CLASS.  YOU IMPLY EVERYONE WHO
///   DERIVES FROM IT IS OBLIGATED TO IMPLEMENT THAT FUNCTION, OR, IF THEY
///   DON'T, THAT THEY ARE IMPLICITLY AGREEING THAT THE BASE CLASS IMPLEMENTATION
///   IS CORRECT.
///
/// @author Monica Berrondo
/// @author Jeff Gray
/// @author Steven Lewis
/// @author Sarel Fleishman

#ifndef INCLUDED_protocols_moves_Mover_hh
#define INCLUDED_protocols_moves_Mover_hh

// Unit Headers
#include <protocols/moves/Mover.fwd.hh>

// Package headers
#include <protocols/moves/MoverStatus.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

//#include <protocols/jobdist/Jobs.hh>
#include <basic/datacache/DataMap.fwd.hh>
// ObjexxFCL Headers

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>
#include <map>
#include <list>
#include <iostream>

#include <protocols/jobdist/Jobs.fwd.hh>
#include <utility/vector1.hh>

#ifdef WIN32
#include <utility/tag/Tag.hh>
#endif

#ifdef USELUA
#include <lua.hpp>
#include <luabind/luabind.hpp>
#endif

#include <core/io/serialization/PipeMap.fwd.hh>
#include <utility/lua/LuaObject.hh>
#include <utility/lua/LuaIterator.hh>

namespace protocols {
namespace moves {

#ifdef USELUA
void lregister_Mover( lua_State * lstate );
void lregister_SerializableState( lua_State * lstate );
void SerializableState_set( SerializableStateSP state, std::string key, std::string val );
void SerializableState_set( SerializableStateSP state, std::string key, core::Real val );
std::string SerializableState_get( SerializableStateSP state, std::string key );
#endif


class Mover : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< Mover >
{
public:
	typedef utility::tag::TagCOP TagCOP;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseCOP PoseCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef std::list< std::string > Strings; // should match jd2::Job::Strings

public:
	Mover();
	virtual ~Mover();

	// Factory<Mover> functions
	// this really should be pure
	virtual MoverSP create();

	static std::string name() {
		return "UNDEFINED NAME";
	}

	// self pointers
	inline MoverCOP get_self_ptr() const { return shared_from_this(); }
	inline MoverOP  get_self_ptr() { return shared_from_this(); }
	inline MoverCAP get_self_weak_ptr() const { return MoverCAP( shared_from_this() ); }
	inline MoverAP  get_self_weak_ptr() { return MoverAP( shared_from_this() ); }

	// elscripts functions
	virtual void apply( core::io::serialization::PipeMap & pmap);

	// called right before mover is used , allowing mover to set settings based on state
	virtual void parse_state( SerializableState const & state );

	// called once, when mover is instantiated
	virtual void parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & score_fxns,
		utility::lua::LuaObject const & tasks,
		MoverCacheSP cache );

	virtual void save_state( SerializableState & state );

	/// @brief sets the type for a mover; name_ has been removed (2010/05/14)
	Mover( std::string const & type_name );

	Mover( Mover const & other );
	Mover& operator=( Mover const& other );


	virtual core::Real last_proposal_density_ratio() {return 1;} //ek added 2/25/10

	/// @brief Overload this static method if you access options within the mover.
	/// @details These options will end up in -help of your application if users of this mover call register_options.
	/// Do this recursively!
	/// If you use movers within your mover, call their register_options in your register_options() method.
	static void register_options() {}

	virtual void apply( Pose & ) = 0;

	std::string const & type() const { return type_; }

	void set_type( std::string const& setting ) {
		type_ = setting;
	}

	/// @brief A tag is a unique identifier used to identify structures produced
	/// by this Mover. get_current_tag() returns the tag, and set_current_tag( std::string tag )
	/// sets the tag.  This functionality is not intended for use with the 2008 job distributor.
	std::string get_current_tag() const {
		//  if ( jd2::jd2_used() ) return jd2::current_output_name();
		return current_tag_;
	}

	virtual void set_current_tag( std::string const & new_tag ) { current_tag_ = new_tag; }

	/// @brief setter for poses contained for rms
	virtual void set_input_pose( PoseCOP pose );

	/// @brief setter for native poses contained for rms ---- we should get rid of this method? it is widely used, but a bit unsafe
	virtual void set_native_pose( PoseCOP pose );

	PoseCOP get_input_pose() const;
	PoseCOP get_native_pose() const; // ---- we should get rid of this method? it is widely used, but a bit unsafe

	/// @brief: Unit test support function.  Apply one move to a given pose.
	///      Allows extra test specific functions to be called before applying
	virtual void test_move( Pose & pose ) {
		apply( pose );
	}

	//protected:
	//OL 4/23/08 made this public. it is not really a safety issue to have that
	//protected and it allows more detail in MC diagnosis
	void type( const std::string & type_in ) { type_ = type_in; }

	/// @brief Return a clone of the Mover object.
	virtual MoverOP clone() const /*= 0*/;

	//////////////////////////////////begin parser interface////////////////////////////

	/// @brief Called by MoverFactory when constructing new Movers. Takes care of the specific mover's parsing.
	virtual
	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		Movers_map const & movers,
		Pose const & pose
	);

	/// @brief Each derived class must specify its name.  The class name.
	virtual std::string get_name() const = 0;
	std::string get_type() const { return( type_ ); }

	///end parser interface, start Job Distributor interface/////////////
	//void output_intermediate_pose( Pose const & pose, std::string const & stage_tag );

	/// @brief returns status after an apply().  The job distributor (august 08 vintage) will check this function to see if your protocol wants to filter its results - if your protocol wants to say "that run was no good, skip it" then use the protected last_move_status(MoverStatus) to change the value that this function will return.
	MoverStatus get_last_move_status() const;

	/// @brief resets status to SUCCESS, meant to be used before an apply().  The job distributor (august 08 vintage) uses this to ensure non-accumulation of status across apply()s.
	void reset_status();

	///fpd
	/// @brief Mechanism by which a mover may return multiple output poses from a single input pose.
	virtual
	core::pose::PoseOP get_additional_output();

	/// @brief Strings container can be used to return miscellaneous info (as std::string) from a mover, such as notes about the results of apply(). The job distributor (Apr 09 vintage) will check this function to see if your protocol wants to add string info to the Job that ran this mover. One way this can be useful is that later, a JobOutputter may include/append this info to an output file.
	/// @brief clear_info is called by jd2 before calling apply
	virtual void clear_info() { info_.clear(); }

	/// @brief non-const accessor
	virtual Strings & info() { return info_; }

	/// @brief const accessor
	virtual Strings const & info() const { return info_; }

	/// @brief Inform the Job Distributor (August '08 vintage) whether this object needs to be freshly regenerated on
	/// each use.
	virtual bool reinitialize_for_each_job() const;

	/// @brief Inform the Job Distributor (August '08 vintage) whether this object needs to be regenerated when the input
	/// pose is about to change, (for example, if the Mover has special code on the first apply() that is only valid for
	/// that one input pose).
	virtual bool reinitialize_for_new_input() const;

	/// @brief Generates a new Mover object freshly created with the default ctor.
	virtual MoverOP fresh_instance() const /*= 0*/;

	///////////////////////////////end Job Distributor interface////////////////////////////////////////

	void set_current_job( protocols::jobdist::BasicJobCOP job );

	jobdist::BasicJobCOP get_current_job() const;

	/// @brief Outputs details about the Mover, including current settings.
	virtual void show(std::ostream & output=std::cout) const;

protected:
	/// @brief nonvirtual setter for MoverStatus last_status_.  Protected means that only the mover itself will be able to change its own status.  The job distributor (august 08 vintage) is aware of status set with this function and will do what the MoverStatus says.
	void set_last_move_status( MoverStatus status );

private:

	std::string type_;
	std::string current_tag_;

	PoseCOP input_pose_;
	PoseCOP native_pose_;

	/// @brief used to track if movers fail their filters - jobdist::BasicJob Distributor queries for this value
	MoverStatus last_status_;
	/// @brief miscellaneous info: optional notes about mover results. Fed to Job by JobDistributor (jd2), allowing JobOutputters to optionally report this information in some way.
	Strings info_;

	/// @brief this field can be NULL or it refers to the current JOB this mover is working on.
	jobdist::BasicJobCOP current_job_;

}; // end Mover base class

} // moves

/// @brief Insertion operator (overloaded so that the Mover can be "printed").
//std::ostream & operator<<(std::ostream & output, moves::Mover const & mover);

} // protocols

#endif //INCLUDED_protocols_moves_Mover_HH
