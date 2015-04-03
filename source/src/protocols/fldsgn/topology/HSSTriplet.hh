// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/topology/HSSTriplet.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_HSSTriplet_hh
#define INCLUDED_protocols_fldsgn_topology_HSSTriplet_hh

// Unit headers
#include <protocols/fldsgn/topology/HSSTriplet.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>

// Project headers
#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>
#include <string>

namespace protocols {
namespace fldsgn {
namespace topology {

class HSSTriplet : public utility::pointer::ReferenceCount {
public:


	typedef core::Size Size;
	typedef core::Real Real;
	typedef std::string String;
	typedef core::Vector Vector;

	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;
public:


	/// @brief default constructor
	HSSTriplet():
		helix_( 0 ),
		strand1_( 0 ),
		strand2_( 0 ),
		hsheet_dist_( 0.0 ),
		hs_angle_( 0.0 ),
		hs1_dist_( 0.0 ),
		hs2_dist_( 0.0 ),
		ss_dist_( 0.0 ),
		ss_orient_( "" ),
		hs1_orient_( "" ),
		hs2_orient_( "" ),
		left_handed_( false ),
		geometry_is_initialized_( false )
	{}

	/// @Brief value constructor
	HSSTriplet(
		Size const h,
		Size const s1,
		Size const s2
	):
		helix_( h ),
		strand1_( s1 ),
		strand2_( s2 ),
		hsheet_dist_( 0.0 ),
		hs_angle_( 0.0 ),
		hs1_dist_( 0.0 ),
		hs2_dist_( 0.0 ),
		ss_dist_( 0.0 ),
		ss_orient_( "" ),
		hs1_orient_( "" ),
		hs2_orient_( "" ),
		left_handed_( false ),
		geometry_is_initialized_( false )
	{}

	/// @brief value constructor
	HSSTriplet( String const & hss );

	/// @brief default destructor
	virtual ~HSSTriplet() ; // auto-removing definition from header{}

	/// @brief copy constructor
	HSSTriplet( HSSTriplet const & hss );


	/// @brief operator ==
	inline
	bool operator ==( HSSTriplet const & rval ) const {
		return ( helix_ == rval.helix_ && strand1_ == rval.strand1_ && strand2_ == rval.strand2_ );
	}


public:


	/// @brief IO Operator
	friend std::ostream & operator<<(std::ostream & out, const HSSTriplet & s );


public:


	inline
	Size helix() const
	{
		return helix_;
	}

	inline
	Size strand1() const
	{
		return strand1_;
	}

	inline
	Size strand2() const
	{
		return strand2_;
	}


public:


	/// @brief reutrn distance between sheet ( defined by the 2 strands ) and helix
	Real hsheet_dist() const;

	/// @brief return distance between sheet ( defined by the 2 strands ) and helix
	Real hs_angle() const;

	/// @brief distance between mid helix and midpoint of 1st strand
	Real hs1_dist() const;

	/// @brief distance between mid helix and midpoint of 2nd strand
	Real hs2_dist() const;

	/// @brief distance between midpoints of strands
	Real ss_dist() const;

	/// @brief orientation between strands
	String ss_orient() const;

	/// @brief orientation between helix and 1st strand
	String hs1_orient() const;

	/// @brief orientation between helix and 2nd strand
	String hs2_orient() const;

	/// @brief hsstriplet is left handed or not
	bool left_handed() const;

	/// @brief geometry is initialized or not
	inline
	bool geometry_is_initialized() const
	{
		return geometry_is_initialized_;
	}


public:


	void calc_geometry( SS_Info2_COP const ssinfo );


private:


	/// @brief helix of hsstriplet
	Size helix_;

	/// @brief 1st strand of hsstriplet
	Size strand1_;

	/// @brief 2nd strand of hsstriplet
	Size strand2_;

	/// @brief distance between sheet ( defined by the 2 strands ) and helix
	Real hsheet_dist_;

	/// @brief distance between sheet ( defined by the 2 strands ) and helix
	Real hs_angle_;

	/// @brief distance between mid helix and midpoint of strand1
	Real hs1_dist_;

	/// @brief distance between mid helix and midpoint of strand2
	Real hs2_dist_;

	/// @brief distance between midpoints of strands
	Real ss_dist_;

	/// @brief orientation between strands
	String ss_orient_;

	/// @brief orientation between helix and 1st strand
	String hs1_orient_;

	/// @brief orientation between helix and 2nd strand
	String hs2_orient_;

	/// @brief hsstriplet is left-handed or not
	bool left_handed_;

	/// @brief geometry is initialized
	bool geometry_is_initialized_;


}; // HSSTriplet


class HSSTripletSet : public utility::pointer::ReferenceCount {
public:


	typedef core::Size Size;
	typedef std::string String;


public:


	/// @brief default constructor
	HSSTripletSet();

	/// @brief value constructor
	HSSTripletSet( HSSTriplets const & s );

	/// @brief value constructor
	HSSTripletSet( String const & s );

	/// @brief copy constructor
	HSSTripletSet( HSSTripletSet const & s );

	/// @brief destructor
	virtual ~HSSTripletSet();


public:// operator


	/// @brief IO Operator
	friend std::ostream & operator<<(std::ostream & out, const HSSTripletSet & s );


public: //mutator


	/// @brief
	void add_hsstriplets( HSSTriplets const & s );

	/// @brief
	void push_back( HSSTripletOP const hssop );

	/// @brief
	void clear();


public: // accessor


	/// @brief return an interator that points to the first HSSTripletOP
	HSSIterator begin() {
		return hss_triplets_.begin();
	}

	/// @brief return an interator that points beyond the last HSSTripletOP
	HSSIterator end() {
		return hss_triplets_.end();
	}

	/// @brief return an const interator that points to the first HSSTripletOP
	HSSConstIterator begin() const {
		return hss_triplets_.begin();
	}

	/// @brief return an const interator that points beyond the last HSSTripletOP
	HSSConstIterator end() const {
		return hss_triplets_.end();
	}

	/// @brief return the size of vector of hss_triplets_
	Size size() const {
		return hss_triplets_.size();
	}


public: //accessor


	/// @brief
	HSSTripletOP hss_triplet( Size const helix );


	/// @brief
	HSSTriplets const & hss_triplets() const;


private: //data


	HSSTriplets hss_triplets_;


	std::map< Size, HSSTripletOP > helix2hss_;


}; // HSSTripletSet


} // namespace topology
} // namespace fldsgn
} // namespace protocols

#endif
