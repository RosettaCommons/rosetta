// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/topology/SS_Info2.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_SS_Info2_hh
#define INCLUDED_protocols_fldsgn_topology_SS_Info2_hh

/// Unit headers
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>

/// Package headers
#include <protocols/fldsgn/topology/BB_Pos.hh>

// Project headers
#include <core/types.hh>
#include <basic/datacache/CacheableData.hh>

/// Numeric headers
#include <numeric/xyzVector.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace fldsgn {
namespace topology {

class SS_Base : public utility::pointer::ReferenceCount {
public:


	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector Vector;
	typedef protocols::fldsgn::topology::BB_Pos BB_Pos;


public:// construct/destruct


	/// @brief default constructor
	SS_Base();

	/// @brief value constructor
	SS_Base( Size const & begin, Size const & end );

	/// @brief value constructor
	// SS_Base( Size const & begin, Size const & end, Vector const & v );

	/// @brief copy constructor
	SS_Base( SS_Base const & s );

	/// @brief destructor
	virtual ~SS_Base();


public:// accessors


	inline Size begin() const { return begin_; }

	inline Size end() const { return end_; }

	inline Size length() const { return end_ - begin_ + 1; }

	inline Vector orient() const { return orient_; }

	inline Vector Nend_orient() const { return Nend_orient_; }

	inline Vector Cend_orient() const { return Cend_orient_; }

	inline Vector Nend_pos() const { return Nend_pos_; }

	inline Vector Cend_pos() const { return Cend_pos_; }

	inline Vector mid_pos() const { return mid_pos_; }


public:// accessors


	inline bool is_geometry_initialized() const { return is_geometry_initialized_; }


protected: // setters


	/// @brief set vector from Ca of Nterminal to Ca of Cterminal
	void orient( Vector const & v ) { orient_ = v; }

	/// @brief set orient vector of N-terminal SS
	void Nend_orient( Vector const & v ) { Nend_orient_ = v; }

	/// @brief set orient vector of C-terminal SS
	void Cend_orient( Vector const & v ) { Cend_orient_ = v; }

	/// @brief set positional vector of N-terminal
	void Nend_pos( Vector const & v ) { Nend_pos_ = v; }

	/// @brief set positional vector of C-terminal
	void Cend_pos( Vector const & v ) { Cend_pos_ = v; }

	/// @brief set positional vector of mid point
	void mid_pos( Vector const & v ) { mid_pos_ = v; }

	/// @brief set geometry is initialized or not
	void is_geometry_initialized( bool const v ) { is_geometry_initialized_ = v; }


private: ///data


	/// @brief begintial residue of strand
	Size begin_;

	/// @brief end residue of strand
	Size end_;

	/// @brief
	bool is_geometry_initialized_;

	/// @brief orient vector from Ca of Nterminal to Ca of Cterminal
	Vector orient_;

	/// @brief orient vector of Nterminal end
	Vector Nend_orient_;

	/// @brief orient vector of Cterminal end
	Vector Cend_orient_;

	/// @brief positional vector of N- and C- teriminals
	Vector Nend_pos_, Cend_pos_;

	/// @brief positional vector of mid point
	Vector mid_pos_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

//////////////////////////////////////////////////////////////////////////////////////////////////////
class Strand : public SS_Base {


	typedef SS_Base Parent;


public:// construct/destruct


	/// @brief default constructor
	Strand();

	/// @brief value constructor
	Strand( Size const & begin, Size const & end );

	/// @brief copy constructor
	Strand( Strand const & s );

	/// @brief destructor
	~Strand();

	/// @brief
	friend std::ostream & operator<<(std::ostream & out, const Strand & st );

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


public:


	// @brief
	void calc_geometry( BB_Pos const & bbpos );


};

//////////////////////////////////////////////////////////////////////////////////////////////////////
class Helix : public SS_Base {


	typedef SS_Base Parent;
	typedef core::Vector Vector;


public: // constructor/destructor


	/// @brief default constructor
	Helix();

	/// @brief value constructor
	Helix( Size const & begin, Size const & end );

	/// @brief copy constructor
	Helix( Helix const & s );

	/// @brief destructor
	~Helix();

	/// @brief
	Real bend() const { return bend_; }

	/// @brief
	friend std::ostream & operator<<( std::ostream & out, const Helix & hx );


public: //


	/// @brief calc geometry of helix
	void calc_geometry( BB_Pos const & bbpos );


private: // data


	Real bend_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

//////////////////////////////////////////////////////////////////////////////////////////////////////
class Loop : public SS_Base {

	typedef std::string String;
	typedef SS_Base Parent;
	typedef core::Vector Vector;


public: // constructor/destructor


	/// @brief default constructor
	Loop();

	/// @brief value constructor
	Loop( Size const & begin, Size const & end, String const & type="" );

	/// @brief copy constructor
	Loop( Loop const & s );

	/// @brief destructor
	~Loop();

	/// @brief
	friend std::ostream & operator<<( std::ostream & out, const Loop & hx );


public:

	String type() const { return type_; };

public:

	/// @brief
	String type_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

//////////////////////////////////////////////////////////////////////////////////////////////////////
class SS_Info2 : public basic::datacache::CacheableData {


	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef std::string String;
	typedef protocols::fldsgn::topology::BB_Pos BB_Pos;


public: // constructor/destructor


	/// @brief default constructor
	SS_Info2();

	/// @brief value constructor
	SS_Info2( String const & secstruct );

	/// @brief value constructor
	SS_Info2( Pose const & pose, String const & secstruct = "" );

	/// @brief copy constructor
	SS_Info2( SS_Info2 const & s );

	/// @brief destructor
	~SS_Info2();


public:


	/// @brief make clone
	basic::datacache::CacheableDataOP clone() const;

	/// @brief output info of SS_Info2
	friend
	std::ostream & operator<<(std::ostream & out, const SS_Info2 & ssinfo );

	/// @brief initialize parameters of this class
	void initialize( String const & secstruct );

	/// @brief initialize parameters of this class
	void initialize( Pose const & pose, String const & secstruct = "" );


public:


	/// @brief get flag for telling whether bb_pos_ was initiliazed by pose or not
	inline
	bool
	bbpos_is_set() const
	{
		return bbpos_is_set_;
	}

	/// @brief string of secondary structure elements
	inline
	String
	secstruct() const
	{
		return secstruct_;
	}

	inline
	char
	secstruct( Size ii ) const
	{
		return secstruct_.at( ii-1 );
	}

	/// @brief get xyz-coordinates of backbone structure
	inline
	BB_Pos const &
	bb_pos() const
	{
		return bb_pos_;
	}

	/// @brief return strands
	inline
	Strands const &
	strands() const
	{
		return strands_;
	}

	/// @brief return helices
	inline
	Helices const &
	helices() const
	{
		return helices_;
	}

	/// @brief return loops
	inline
	Loops const &
	loops() const
	{
		return loops_;
	}

	/// @brief return owning pointer of strand given an index of strands
	StrandCOP const
	strand( Size is ) const
	{
		runtime_assert( is <= strands_.size() );
		return strands_[ is ];
	}

	/// @brief return owning pointer of helix given an index of helices
	HelixCOP const
	helix( Size ih ) const
	{
		runtime_assert( ih <= helices_.size() );
		return helices_[ ih ];
	}

	/// @brief return owning pointer of loop given an index of loops
	LoopCOP const
	loop( Size il ) const
	{
		runtime_assert( il <= loops_.size() );
		return loops_[ il ];
	}

	/// @brief return strand index in strands given a residue number
	inline
	Size
	strand_id( Size const nres ) const
	{
		return strand_id_[ nres ];
	}

	/// @brief return helix index in helices given a residue number
	inline
	Size
	helix_id( Size const nres ) const
	{
		return helix_id_[ nres ];
	}

	/// @brief return loop index in loops given a residue number
	inline
	Size
	loop_id( Size const nres ) const
	{
		return loop_id_[ nres ];
	}

	/// @brief return the index of secondary structure element given a residue number
	inline
	Size
	ss_element_id( Size const nres ) const
	{
		return ss_element_id_[ nres ];
	}


public:


	/// @brief set orientation vector of secondary structures given a pose
	void set_SSorient( Pose const & pose );

	/// @brief set orientation vector of secondary structures given a pose which is defined in the constructor
	void set_SSorient();

	/// @brief clear data
	void clear_data();


private:


	/// @brief
	void
	resize( Size const nres );

	/// @brief identify secondary structures
	void identify_ss( String const & secstruct );


private: //data


	/// @brief flag for telling whether bb_pos_ was initiliazed by pose or not
	bool bbpos_is_set_;

	/// @brief string of secondary structure elements
	String secstruct_;

	/// @brief xyz-coordinates of backbone
	BB_Pos bb_pos_;

	/// @brief vector of StrandOP
	Strands strands_;

	/// @brief vector for storing index of strand id for each residue position
	utility::vector1< Size > strand_id_;

	/// @brief vector of HelixOP
	Helices helices_;

	/// @brief vector for storing index of helix id for each residue position
	utility::vector1< Size > helix_id_;

	/// @brief vector of loops
	Loops loops_;
	utility::vector1< Size > loop_id_;

	/// @brief vector for storing index of secondary structure element for each residue position
	utility::vector1< Size > ss_element_id_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace topology
} // namespace fldsgn
} // naemspace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_fldsgn_topology_SS_Info2 )
#endif // SERIALIZATION


#endif
