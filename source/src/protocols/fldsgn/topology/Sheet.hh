// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/topology/Sheet.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_Sheet_hh
#define INCLUDED_protocols_fldsgn_topology_Sheet_hh

// unit headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/fldsgn/topology/Sheet.fwd.hh>

// project headers
#include <core/types.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


namespace protocols {
namespace fldsgn {
namespace topology {

class Sheet : public utility::pointer::ReferenceCount {
public:


	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef utility::vector1< Size > VecSize;
	typedef utility::vector1< int > VecInt;
	typedef utility::vector1< Real > VecReal;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;


public:// construct/destruct


	/// @brief default constructor
	Sheet();

	/// @brief value constructor
	Sheet( VecSize const & order_strands,	VecInt const & orient_strands, bool is_barrel );

	/// @brief copy constructor
	Sheet( Sheet const & s );

	/// @brief default destructor
	virtual ~Sheet();

	/// @brief clone this object
	SheetOP clone() {
		return SheetOP( new Sheet( *this ) );
	}

	/// @brief return strand pairing
	friend
	std::ostream & operator<<( std::ostream & out, const Sheet &s );


public:


	/// @brief intialize this class
	void initialize();


public: //accessors


	/// @brief the number strands inclued in
	Size num_strands() const {	return num_strands_;	}

	/// @brief is this barrel ?
	bool is_barrel() const { return is_barrel_; }

	VecSize order_strands() const {	return order_strands_; }

	Size order_strand( Size const s ) const { return order_strands_[ s ]; }

	VecInt orient_strands() const { return orient_strands_; }

	int orient_strand( Size const s ) const { return orient_strands_[ s ]; }

	Size strand_order( Size const s ) {	return strand_order_[ s ];	}

	VecInt ca_cb_orients() const { return ca_cb_orients_; }

	int ca_cb_orient( Size const s ) const { return ca_cb_orients_[ s ]; }


public:

	/// @brief
	int
	which_side( Vector const vec ) const;

	/// @brief calc surface areas only with beta-sheet
	utility::vector1< Real >
	calc_sasa_bothsides( Pose const & pose, SS_Info2_COP const ssinfo, Real pore_radius=1.5 );

	/// @brief calc geometry of sheet, sheet_plane_, sheet_center_, ca_cb_orients_
	void
	calc_geometry( SS_Info2_COP const ssinfo );


public: //

	bool is_member( Size const s );

private:  // data

	/// @brief
	Size num_strands_;

	/// @brief
	bool is_barrel_;

	/// @brief order of strand in sheet -> id of strand
	VecSize order_strands_;

	/// @brief order of strand in sheet -> id of strand
	VecInt orient_strands_;

	/// @brief id of strand -> order of strand in sheet
	std::map< Size, Size > strand_order_;

	/// @brief vector defining sheet plane
	Vector sheet_plane_;

	/// @brief "center" of sheet
	Vector sheet_center_;

	VecInt ca_cb_orients_;

	/// @brief geometries was calculated or not
	bool is_geometry_initialized_;

};


class SheetSet : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~SheetSet();


	typedef core::Size Size;
	typedef core::Real Real;
	typedef utility::vector1< int > VecInt;
	typedef utility::vector1< Size > VecSize;
	typedef utility::vector1< Real > VecReal;
	typedef protocols::fldsgn::topology::StrandPairingSet StrandPairingSet;
	typedef protocols::fldsgn::topology::StrandPairingSetCOP StrandPairingSetCOP;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;


public:// construct/destruct


	/// @brief default constructor
	SheetSet();

	/// @brief value constructor
	SheetSet( Sheets const & sheets );

	/// @brief value constructor
	SheetSet( SS_Info2_COP const ssinfo, StrandPairingSetCOP const spairset );

	/// @brief copy constructor
	SheetSet( SheetSet const & s );

	/// @brief return strand pairing
	friend
	std::ostream & operator<<( std::ostream & out, const SheetSet &s );


public:


	void initialize( SS_Info2_COP const ssinfo, StrandPairingSetCOP const spairset );


public: // mutators


	/// @brief
	void push_back( SheetOP const sop );

	/// @brief
	void clear();


public: // accessors


	/// @brief
	SheetOP sheet( Size const s ) const;

	/// @brief
	inline Sheets sheets() const { return sheets_; }

	/// @brief return number of sheets
	inline Size	num_sheets() const { return sheets_.size();	}

	/// @brief return number of sheets
	inline Size	size() const { return sheets_.size();	}

	/// @brief return the id of sheet that a given strand belongs to.
	Size which_sheet( Size const s ) const;

	/// @brief return strand pairing set
	StrandPairingSet spairset() const;

public:


	/// @brief
	void calc_geometry( SS_Info2_COP const ssinfo );


private://


	/// @brief
	void set_sheet_number() const;


private:  // data


	/// @brief sheet number given a strand
	mutable std::map< Size, Size > sheet_number_;

	/// @brief
	Sheets sheets_;

	/// @brief
	StrandPairingSet spairset_;


};

} // namespace topology
} // namespace fldsgn
} // namespace protocol

#endif
