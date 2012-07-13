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
#include <protocols/fldsgn/topology/StrandPairing.hh>

// Project headers
#include <core/types.hh>
#include <basic/datacache/CacheableData.hh>

/// Numeric headers
#include <numeric/xyzVector.hh>

// C++ headers
#include <string>
#include <map>
#include <utility/vector1.hh>


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

	/// @brief default constructor
	SS_Base( char const & name );
	
	/// @brief value constructor
	SS_Base( char const & name, Size const & begin, Size const & end );

	/// @brief value constructor
	// SS_Base( Size const & begin, Size const & end, Vector const & v );

	/// @brief copy constructor
	SS_Base(	SS_Base const & s );

	/// @brief destructor
	~SS_Base();


public:// accessors

	inline char name() const { return name_; }
	
	inline Size begin() const {	return begin_; }

	inline Size end() const {	return end_; }

	inline Size length() const { return end_ - begin_ + 1; }

	inline Vector orient() const { return orient_; }

	inline Vector Nend_orient() const { return Nend_orient_; }

	inline Vector Cend_orient() const { return Cend_orient_; }

	inline Vector Nend_pos() const { return Nend_pos_; }

	inline Vector Cend_pos() const { return Cend_pos_; }

	inline Vector mid_pos() const { return mid_pos_; }
	
	bool is_member( Size const s ) const;


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

	/// @brief character description of this secondary structure, H, E, L
	char name_;
	
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
	Strand(	Strand const & s );

	/// @brief destructor
	~Strand();

  /// @brief
	friend std::ostream & operator<<(std::ostream & out, const Strand & st );


public: // accessor
	
	
	// @brief strand has bulge ?
	bool has_bulge() const;
	
	// @brief return reidue numbers of bulges
	inline utility::vector1< Size > bulges() const { return bulges_; }
	
	// @brief
	bool set_bulge( Size const & s );

	
public:


	// @brief calculate geometry of strand
	void calc_geometry( BB_Pos const & bbpos );	
	

private:

	/// bulges_ are initialized by set_bulges of SS_Info2 based on StrandPairingSet
	/// @brief residue numbers of bulge
	utility::vector1< Size > bulges_;


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


};

//////////////////////////////////////////////////////////////////////////////////////////////////////
class SS_Info2 : public basic::datacache::CacheableData {


	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef std::string String;
	typedef protocols::fldsgn::topology::BB_Pos BB_Pos;
	typedef protocols::fldsgn::topology::StrandPairingSetOP StrandPairingSetOP;
	typedef protocols::fldsgn::topology::StrandPairingOP StrandPairingOP;
	

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

	/// @brief redifine this with abego information
	/// Currently, this extends strand definition to N-terminal side 
	/// if the abego of loop region is 2 consecutive beta value, B, Z, Y, P, S
	void redefine_with_abego();
	
	/// @brief set bulge in Strand
	void set_bulge( StrandPairingSetCOP const spairset );
	
	
public: //accessor


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
	
	/// @brief
	inline
	char
	secstruct( Size ii ) const
	{
		return secstruct_.at( ii-1 );		
	}
	
	/// @brief
	inline
	utility::vector1< String > 
	abego() const
	{
		return abego_;		
	}

	/// @brief get string of abego given a residue number
	inline
	String
	abego( Size const ii ) const
	{
		return abego_[ ii ];
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
	
	/// @brief return HE_elements
	inline
	utility::vector1<Size> const &
	HE_elements() const
	{
		return HE_elements_;				
	}

	/// @brief return owning pointer of strand given an index of strands
	StrandOP const
	strand( Size is ) const
	{
		runtime_assert( is <= strands_.size() );
		return strands_[ is ];
	}

	/// @brief return owning pointer of helix given an index of helices
	HelixOP const
	helix( Size ih ) const
	{
		runtime_assert( ih <= helices_.size() );
		return helices_[ ih ];
	}

	/// @brief return owning pointer of loop given an index of loops
	LoopOP const
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
	
	/// @brief return ss element OP given a ss element id
	SS_BaseCOP const
	ss_element_id2op( Size const id )
	{
		runtime_assert( ss_element_id_.back() >= id );
		return ss_element_id2op_[ id ];
	}
	
	/// @brief get bulges for given a strand id
	utility::vector1< Size > const bulges( Size const id ) const;

	/// @brief get SS_element id from helix, strand, or loop index
	Size get_ssid_by_hleid( char const & ss, Size const & index );
	
	/// @brief get helix, strand, or loop index from SS_element id
	String get_hleid_by_ssid( Size const & index );
	
	/// @brief get chirality for consecutive 3 secondary structur elements given a first ss element id
	char
	get_chiral( Size const id ) const;
	
	/// @brief get chirality for consecutive 3 secondary structur elements given a first ss element described H, E, or L w/ its index
	char
	get_chiral( char const & ss, Size const & index );

	
public: // mutator
	
	
	/// @brief set minimum lengtht of helix ( default 0 )
	void set_min_helix_length( Size const s ) { min_helix_length_ = s; }
	
	/// @brief set minimum lengtht of strand ( default 0 )
	void set_min_strand_length( Size const s ) { min_strand_length_ = s; }

	/// @brief set secondary structure
	void secstruct( Size const ii, char const & aa  );
	
	/// @brief set abego
	void abego( Size const ii, String const & aa );


public:


	/// @brief make string of HE elements
	String make_string_HE_elements() const;
	
	/// @brief make string of chiralities
	String make_string_chiralities() const;	
	
	/// @brief set orientation vector of secondary structures given a pose
	void set_SSorient( Pose const & pose );

	/// @brief set orientation vector of secondary structures given a pose which is defined in the constructor
	void set_SSorient();

	///	@brief set chirality of consecutive secondary structure elements
	void set_chiral_consecutive_sstriplets();
		
	/// @brief clear data
	void clear_data();
	
	
public: // utility
	
	/// @brief set secondary structure information into pose
	void set_ss_into_pose( Pose & pose );

private:


	/// @brief
	void
	resize( Size const nres );

	/// @brief identify secondary structures
	void identify_ss( String const & secstruct );
	
	/// @brief modify secstruct based on the length of ss elements, depending on min_helix_length_ and min_strand_length_
	void reduce_secstruct();
	
	/// @brief set orientation vector of helix
	void set_SSorient_helix();

	/// @brief set orientation vector of strand
	void set_SSorient_strand();

	
private: // helper function
	

	/// @brief relate helix, strand, loop index to sselement id
	void set_hleid_ssid( char const & ss, Size const & index, Size const & id );

												
private: //data


	/// @brief flag for telling whether bb_pos_ was initiliazed by pose or not
	bool bbpos_is_set_;

	/// @brief string of secondary structure elements
	String secstruct_;
	
	/// @brief minimum length of helix (default 0)
	Size min_helix_length_;

	/// @brief minimum length of strand (default 0)
	Size min_strand_length_;

	/// @brief string of abego
	utility::vector1< std::string > abego_;

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

	/// @brief vector for SS element id for each residue position
	utility::vector1< Size > ss_element_id_;
	
	/// @brief map from SS element id to SS element OPs
	mutable std::map< Size, SS_BaseOP > ss_element_id2op_;
	
	/// @brief map from helix, strand, or loop id to SS element id
	std::map< String, Size > hleid_to_ssid_;
	
	/// @brief map from helix, strand, or loop id to SS element id
	std::map< Size, String > ssid_to_hleid_;
 
	/// @brief vector for secondary structure elelement ids without loop elements
	utility::vector1< Size > HE_elements_;
	
	/// will be removed !!!! 
	/// @brief chirality of consecutive secondary structure elements
	utility::vector1< char > chiral_sstriplet_;
	
	/// @brief
	bool initialized_by_spairset_;

};

} // namespace topology
} // namespace fldsgn
} // naemspace protocols

#endif
