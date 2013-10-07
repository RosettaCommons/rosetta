// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/vall/VallResidue.hh
/// @brief  class for managing a line of the Vall fragment library
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_VallResidue_hh
#define INCLUDED_core_fragment_picking_old_vall_VallResidue_hh


// unit headers
#include <core/fragment/picking_old/vall/VallResidue.fwd.hh>

// type headers
#include <core/types.hh>

// project headers
#include <core/chemical/AA.hh>
#include <core/fragment/BBTorsionSRFD.hh>

// utility headers
#include <utility/pointer/owning_ptr.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

// C++ headers
#include <iostream>

#include <utility/vector1.hh>




namespace core {
namespace fragment {
namespace picking_old {
namespace vall {


/// @brief class for managing a line of the Vall fragment library
class VallResidue {


public: // typedefs


	typedef core::Size Size;
	typedef core::Real Real;
	typedef std::string String;

	typedef core::fragment::BBTorsionSRFD BBTorsionSRFD;
	typedef core::fragment::BBTorsionSRFDOP BBTorsionSRFDOP;


public: // construct/destruct


	/// @brief default constructor
	VallResidue();


	/// @brief string constructor
	VallResidue( String const & line );


	/// @brief copy constructor
	VallResidue( VallResidue const & rval );


	/// @brief default destructor
	~VallResidue();


public: // copy assignment


	/// @brief copy assignment
	VallResidue & operator =( VallResidue const & rval );


public: // conversion


	/// @brief build a BBTorsionSRFD of the given type from this page
	/// @param[in] srfd_type BBTorsionSRFD::create() will be called from this object.
	/// @return A BBTorsionSRFD of the given type initialized with the backbone
	///  torsion information from this page.
	BBTorsionSRFDOP bbtorsion_srfd( BBTorsionSRFD const & srfd_type = BBTorsionSRFD() ) const;


public: // accessors


	/// @brief id of fragment source (e.g. pdb name)
	inline
	String const & id() const {
		return id_;
	}


	/// @brief one letter amino acid code
	inline
	char aa() const {
		return aa_;
	}


	/// @brief one letter secondary structure code
	inline
	char ss() const {
		return ss_;
	}


	/// @brief residue sequence number in source
	inline
	Size resi() const {
		return resi_;
	}


	/// @brief x-coordinate of C-alpha
	inline
	Real x() const {
		return x_;
	}


	/// @brief y-coordinate of C-alpha
	inline
	Real y() const {
		return y_;
	}


	/// @brief z-coordinate of C-alpha
	inline
	Real z() const {
		return z_;
	}


	/// @brief phi backbone torsion in degrees
	inline
	Real phi() const {
		return phi_;
	}


	/// @brief psi backbone torsion in degrees
	inline
	Real psi() const {
		return psi_;
	}


	/// @brief omega backbone torsion in degrees
	inline
	Real omega() const {
		return omega_;
	}


	/// @brief per amino acid profile data
	inline
	utility::vector1< Real > const & profile() const {
		return profile_;
	}


	/// @brief stores the 1-based indexing for accessing this residue
	///  via VallSection::operator []
	inline
	Size position_index() const {
		return position_index_;
	}


	/// @brief stores the 1-based indexing for accessing the VallSection
	///  this residue is located in via VallLibrary::operator []
	inline
	Size section_index() const {
		return section_index_;
	}


public: // mutators


	/// @brief id of fragment source (e.g. pdb name)
	inline
	void id( String const & s ) {
		id_ = s;
	}


	/// @brief one letter amino acid code
	inline
	void aa( char const c ) {
		aa_ = c;
	}


	/// @brief one letter secondary structure code
	inline
	void ss( char const c ) {
		ss_ = c;
	}


	/// @brief residue sequence number in source
	inline
	void resi( Size const i ) {
		resi_ = i;
	}


	/// @brief x-coordinate of C-alpha
	inline
	void x( Real const val ) {
		x_ = val;
	}


	/// @brief y-coordinate of C-alpha
	inline
	void y( Real const val ) {
		y_ = val;
	}


	/// @brief z-coordinate of C-alpha
	inline
	void z( Real const val ) {
		z_ = val;
	}


	/// @brief phi backbone torsion in degrees
	inline
	void phi( Real const val ) {
		phi_ = val;
	}


	/// @brief psi backbone torsion in degrees
	inline
	void psi( Real const val ) {
		psi_ = val;
	}


	/// @brief omega backbone torsion in degrees
	inline
	void omega( Real const val ) {
		omega_ = val;
	}


	/// @brief per amino acid profile data
	inline
	void profile( utility::vector1< Real > const & v ) {
		profile_ = v;
	}


	/// @brief stores the 1-based indexing for accessing this residue
	///  via VallSection::operator []
	inline
	void position_index( Size const idx ) {
		position_index_ = idx;
	}


	/// @brief stores the 1-based indexing for accessing the VallSection
	///  this residue is located in via VallLibrary::operator []
	inline
	void section_index( Size const idx ) {
		section_index_ = idx;
	}


public: // query


	/// @brief has profile info?
	inline
	bool has_profile() const {
		return profile_.size() > 0;
	}


public: // methods


	/// @brief fill internal data from string
	void fill_from_string( String const & line );


private: // static methods


	/// @brief return a vector specifying the order of profile data in Vall
	static utility::vector1< core::chemical::AA > order_vector();


	/// @brief return a formatting string for fill_from_string() dependent
	///  upon actual type of core::Real and core::Size
	/// @remarks This is necessary for sscanf; wrong type can give wrong
	///  input.
	static String format_string();
  static String format_string_2011();

private: // data


	/// @brief id of fragment source (e.g. pdb name)
	String id_;


	/// @brief one letter amino acid code
	char aa_;


	/// @brief one letter secondary structure code
	char ss_;


	/// @brief residue sequence number in source
	Size resi_;


	/// @brief x-coordinate of C-alpha
	Real x_;


	/// @brief y-coordinate of C-alpha
	Real y_;


	/// @brief z-coordinate of C-alpha
	Real z_;


	/// @brief phi backbone torsion in degrees
	Real phi_;


	/// @brief psi backbone torsion in degrees
	Real psi_;


	/// @brief omega backbone torsion in degrees
	Real omega_;


	/// @brief per amino acid profile data
	utility::vector1< Real > profile_;


	/// @brief stores the 1-based indexing for accessing this residue
	///  via VallSection::operator []
	Size position_index_;


	/// @brief stores the 1-based indexing for accessing the VallSection
	///  this residue is located in via VallLibrary::operator []
	Size section_index_;


private: // static data


	/// @brief order of amino acid profile data in Vall
	static utility::vector1< core::chemical::AA > order_;


	/// @brief formatting string for fill_from_string()
	static String format_;
	static String format_2011_;


};


/// @brief wrapper for a collection of VallResidue
class VallResidues : public utility::vector1< VallResidue > {


public: // typedefs


	/// @brief Page typedef required to satisfy concept
	typedef VallResidue Page;


};


} // vall
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_vall_VallResidue_HH */
