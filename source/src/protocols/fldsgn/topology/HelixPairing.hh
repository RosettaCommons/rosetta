// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/topology/HelixPairing.hh
/// @brief header file of HelixPairing.cc
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_HelixPairing_hh
#define INCLUDED_protocols_fldsgn_topology_HelixPairing_hh

// Unit headers
#include <protocols/fldsgn/topology/HelixPairing.fwd.hh>

// Project headers
#include <core/types.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <string>


namespace protocols {
namespace fldsgn {
namespace topology {

class HelixPairing : public utility::pointer::ReferenceCount {
public:

  typedef std::string String;
  typedef core::Size Size;
  typedef core::Real Real;
	typedef core::Vector Vector;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;

public:// construct/destruct


  /// @brief default constructor
  HelixPairing();

  /// @brief value constructor
  HelixPairing(
     Size const h1,
     Size const h2,
     char const o
								);


	/// @brief copy constructor
  HelixPairing( String const & hp );

	/// @brief copy constructor
  HelixPairing( HelixPairing const & hp );

  /// @brief default destructor
  virtual ~HelixPairing();

	/// @brief clone this object
	HelixPairingOP clone();

  /// @brief return strand pairing
	friend
	std::ostream & operator<<(std::ostream & out, const HelixPairing &hp);


public: //accessors


	/// @brief the strand number of the 1st strand in strand pairing
	inline Size h1() const	{	return h1_;	}

	/// @brief the strand number of the 2nd strand in strand pairing
	inline Size h2() const	{	return h2_; }

	/// @brief orientation, parallel or anti-parallel, of helix pairing
	inline char orient() const { return orient_; }

	/// @brief HelixPairing is descripbed as s1()-s2().orient()
	/// For example, 2-3.A means 2nd and 3rd helices make anti-parallel helix_pairing
	inline String name() const { return name_; }

	/// @brief
	inline Real dist() const { return dist_; }

	/// @brief
	inline Real cross_angle() const { return cross_angle_; }

	/// @brief helix angle projected on sheet helices belongs to
	inline Real align_angle() const { return align_angle_; }

	/// @brief is parallel
	bool is_parallel() const;


public:

	/// @brief
	void calc_geometry( SS_Info2_COP const ss_info );


	//private: // initialize


	void initialize();


private:  // data


  /// @brief Helix number of first strand in the strand pair
  Size h1_;

  /// @brief Helix number of second strand in the strand pair
  Size h2_;

  /// @brief two helices make a pair  by parallel, "P", anti parallel, "A", and if not defined, "N"
  char orient_;

	/// @brief helix_pairing as in the style: h1_-h2_.orient_
	String name_;

	/// @brief
	Real dist_;

	/// @brief
	Real cross_angle_;

	/// @brief
	Real align_angle_;

	/// @brief
	Size loop_length_;


};


////////////////////////////////////////////////////////////////////////////////////////////////////////////
class HelixPairingSet : public utility::pointer::ReferenceCount {
public: // typedef


	typedef std::string String;
	typedef core::Size Size;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;

public:// construct/destruct


  /// @brief default constructor
  HelixPairingSet();

  /// @brief value constructor
  HelixPairingSet( HelixPairings const & helix_pairings );

  /// @brief value constructor
  HelixPairingSet( String const & helix_pairings );

	/// @brief copy constructor
  HelixPairingSet( HelixPairingSet const & s );

  /// @brief default destructor
  virtual ~HelixPairingSet();

	/// @brief clone this object
	HelixPairingSetOP clone() const;

  /// @brief return strand pairing
	friend std::ostream & operator<<( std::ostream & out, const HelixPairingSet &s );


public: // mutators


	/// @brief add HelixPairingOP to StrandPairingSet
	void push_back( HelixPairingOP const hop );

	/// @brief clear data of this HelixPairingSet
	void clear();


public: // accessors


	/// @brief return one of the strand_pairings given a number
	HelixPairingOP helix_pairing( Size const s ) const;

	/// @brief return the pointer of the helix pairing, given the two helix numbers of h1, and h2
	/// if h1 and h2 does not make pairing, return 0
	HelixPairingOP helix_pairing( Size const h1, Size const h2 );

	/// @brief return all helix pairings
	HelixPairings const &	helix_pairings() const;

	/// @brief return the size of helix_pairings_
	Size size() const;


public:


	/// @brief calc geomtry of helix pairing
	void  calc_geometry( SS_Info2_COP ss_info );

	/// @brief the name of HelixPairingSet is expressed by the combination of helix pairings
	String name() const;



private:


	/// @brief create map _strand_pairings_
	void create_map_helix_pairings();


private:// data


	/// @brief vector1 including owning pointers of HelixPairing
	HelixPairings helix_pairings_;

	/// @brief the name of HelixPairingSet is expressed by the combination of helix pairings
	String hpairset_name_;

	/// @brief the total number of strands included in HelixPairingSet
	Size num_helices_;

	/// @brief whether the map_helix_pairings_ is initialized or not
	bool initialize_map_helix_pairings_;

	/// @brief 2D table of the pointer of helix pairing, which is sorted by the helix number
	utility::vector1< utility::vector1< HelixPairingOP > > map_helix_pairings_;


}; // HelixPairingSet

} // namespace topology
} // namespace fldsgn
} // namespace protocols

#endif
