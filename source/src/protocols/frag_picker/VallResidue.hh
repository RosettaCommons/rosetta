// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/VallResidue.hh
/// @brief  class for managing a line of the Vall fragment library
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_frag_picker_VallResidue_hh
#define INCLUDED_protocols_frag_picker_VallResidue_hh

// unit headers
#include <protocols/frag_picker/VallResidue.fwd.hh>

// project headers
#include <protocols/frag_picker/ContactTypes.hh>

// type headers
#include <core/types.hh>

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// project headers
#include <core/chemical/AA.hh>

// C++ headers
#include <iostream>
#include <map>

//Auto Headers
#include <core/fragment/BBTorsionSRFD.fwd.hh>
#include <utility/vector1_bool.hh>


namespace protocols {
namespace frag_picker {

/// @brief class for managing a line of the Vall fragment library
class VallResidue: public utility::pointer::ReferenceCount {

public:
	// typedefs

	typedef std::string String;

	typedef core::fragment::BBTorsionSRFD BBTorsionSRFD;
	typedef core::fragment::BBTorsionSRFDOP BBTorsionSRFDOP;

public:
	// construct/destruct

	/// @brief default constructor
	VallResidue();

	/// @brief string constructor
	VallResidue(String const & line);

	/// @brief copy constructor
	VallResidue(VallResidue const & rval);

	/// @brief default destructor
	~VallResidue() override;

public:
	/// @brief copy assignment
	VallResidue & operator =(VallResidue const & rval);

public:
	/// @brief build a BBTorsionSRFD from this page
	BBTorsionSRFDOP bbtorsion_srfd() const;

public:
	/// @brief id of fragment source (e.g. pdb name)
	inline String const & id() const {
		return id_;
	}

	/// @brief sets the key id
	void key(core::Size key) { key_ = key; }

	/// @brief sets the key id
	core::Size key() { return key_; }

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

	/// @brief one letter secondary structure STR code
	inline
	char ss_str() const {
		return ss_str_;
	}

	/// @brief residue sequence number in source
	inline core::Size resi() const {
		return resi_;
	}

	/// @brief x-coordinate of C-alpha
	inline core::Real x() const {
		return x_;
	}

	/// @brief y-coordinate of C-alpha
	inline core::Real y() const {
		return y_;
	}

	/// @brief z-coordinate of C-alpha
	inline core::Real z() const {
		return z_;
	}

	/// @brief x-coordinate of C-beta
	inline core::Real cbx() const {
		return cbx_;
	}

	/// @brief y-coordinate of C-beta
	inline core::Real cby() const {
		return cby_;
	}

	/// @brief z-coordinate of C-beta
	inline core::Real cbz() const {
		return cbz_;
	}

	/// @brief x-coordinate of centroid
	inline core::Real cenx() const {
		return cenx_;
	}

	/// @brief y-coordinate of centroid
	inline core::Real ceny() const {
		return ceny_;
	}

	/// @brief z-coordinate of centroid
	inline core::Real cenz() const {
		return cenz_;
	}

	/// @brief solvent accessible area
	inline core::Real sa() const {
		return sa_;
	}

	/// @brief solvent accessible area normalized
	inline core::Real sa_norm() const {
		return sa_norm_;
	}

	/// @brief phi backbone torsion in degrees from DSSP
	inline core::Real dssp_phi() const {
		return dssp_phi_;
	}

	/// @brief psi backbone torsion in degrees from DSSP
	inline core::Real dssp_psi() const {
		return dssp_psi_;
	}

	/// @brief all-atom residue depth
	inline core::Real depth() const {
		return all_atom_residue_depth_;
	}

	/// @brief number of alignments
	inline core::Size nali() const {
		return nali_;
	}

	/// @brief phi backbone torsion in degrees
	inline core::Real phi() const {
		return phi_;
	}

	/// @brief psi backbone torsion in degrees
	inline core::Real psi() const {
		return psi_;
	}

	/// @brief omega backbone torsion in degrees
	inline core::Real omega() const {
		return omega_;
	}

	/// @brief b factor average for backbone
	inline core::Real bF() const {
		return bF_;
	}

	/// @brief per amino acid profile data
	inline utility::vector1<core::Real> const & profile() const {
		return profile_;
	}

	/// @brief per amino acid structure profile data
	inline utility::vector1<core::Real> const & profile_struct() const {
		return profile_struct_;
	}

	/// @brief secondary chemical shifts
	inline utility::vector1<core::Real> const & secondary_shifts() {
		return sec_shift_data_;
	}

	/// @brief per amino acid profile data
	inline void profile(utility::vector1<core::Real> const & v) {
		profile_ = v;
	}

	/// @brief per amino acid structure profile data
	inline void profile_struct(utility::vector1<core::Real> const & v) {
		profile_struct_ = v;
	}


public:
	// mutators

	/// @brief id of fragment source (e.g. pdb name)
	inline
	void id(String const & s) {
		id_ = s;
	}

	/// @brief one letter amino acid code
	inline
	void aa(char const c) {
		aa_ = c;
	}

	/// @brief one letter secondary structure code
	inline
	void ss(char const c) {
		ss_ = c;
	}

	/// @brief one letter secondary structure STR code
	inline
	void ss_str(char const c) {
		ss_str_ = c;
	}

	inline
	void set_ss_from_str() {
		switch(ss_str_) {
		case 'H' :
			ss_ = 'H';
			break;
		case 'G' :
			ss_ = 'H';
			break;
		case 'I' :
			ss_ = 'H';
			break;
		case 'A' :
			ss_ = 'E';
			break;
		case 'E' :
			ss_ = 'E';
			break;
		case 'M' :
			ss_ = 'E';
			break;
		case 'P' :
			ss_ = 'E';
			break;
		case 'Q' :
			ss_ = 'E';
			break;
		case 'Z' :
			ss_ = 'E';
			break;
		case 'B' :
			ss_ = 'E';
			break;
		default :
			ss_ = 'L';
		}
	}

	/// @brief residue sequence number in source
	inline
	void resi(core::Size const i) {
		resi_ = i;
	}

	/// @brief x-coordinate of C-alpha
	inline
	void x(core::Real const val) {
		x_ = val;
	}

	/// @brief y-coordinate of C-alpha
	inline
	void y(core::Real const val) {
		y_ = val;
	}

	/// @brief z-coordinate of C-alpha
	inline
	void z(core::Real const val) {
		z_ = val;
	}


	/// @brief x-coordinate of C-beta
	inline
	void cbx(core::Real const val) {
		cbx_ = val;
	}

	/// @brief y-coordinate of C-beta
	inline
	void cby(core::Real const val) {
		cby_ = val;
	}

	/// @brief z-coordinate of C-beta
	inline
	void cbz(core::Real const val) {
		cbz_ = val;
	}

	/// @brief x-coordinate of centroid
	inline
	void cenx(core::Real const val) {
		cenx_ = val;
	}

	/// @brief y-coordinate of centroid
	inline
	void ceny(core::Real const val) {
		ceny_ = val;
	}

	/// @brief z-coordinate of centroid
	inline
	void cenz(core::Real const val) {
		cenz_ = val;
	}


	/// @brief solvent accessible area
	inline
	void sa(core::Real const val) {
		sa_ = val;
	}

	/// @brief number of alignments
	inline
	void nali(core::Size const val) {
		nali_ = val;
	}

	/// @brief all-atom residue depth
	inline
	void depth(core::Real const depth) {
		all_atom_residue_depth_ = depth;
	}

	/// @brief phi backbone torsion in degrees
	inline
	void phi(core::Real const val) {
		phi_ = val;
	}

	/// @brief psi backbone torsion in degrees
	inline
	void psi(core::Real const val) {
		psi_ = val;
	}

	/// @brief omega backbone torsion in degrees
	inline
	void omega(core::Real const val) {
		omega_ = val;
	}

	/// @brief b factor average for backbone
	inline
	void bF(core::Real const val) {
		bF_ = val;
	}

public:
	// query

	/// @brief has profile info?
	inline
	bool has_profile() const {
		return profile_.size() > 0;
	}

	/// @brief has structure profile info?
	inline
	bool has_profile_struct() const {
		return profile_struct_.size() > 0;
	}

	inline
	bool has_chemical_shifts() const {
		return sec_shift_data_.size() > 0;
	}

public:
	// methods

	/// @brief fill internal data from string
	void fill_from_string(String const & line);

	/// @brief fill internal data from string
	void fill_from_string_cs(String const & line);

	/// @brief fill internal data from string
	void fill_from_string_version1(String const & line);

	void fill_from_string_residue_depth_version1(String const & line);

	core::Real distance_squared( VallResidueCOP r );

	core::Real distance_squared_cb( VallResidueCOP r );

	core::Real distance_squared_cen( VallResidueCOP r );

	core::Real distance_squared( VallResidueCOP r, ContactType const & type );

private:
	// static methods

	/// @brief return a vector specifying the order of profile data in Vall
	static utility::vector1<core::chemical::AA> order_vector();

	/// @brief return a formatting string for fill_from_string() dependent
	///  upon actual type of core::Real and core::Size
	/// @remarks This is necessary for sscanf; wrong type can give wrong
	///  input.
	static String format_string();
	/// @brief return a formatting string for fill_from_string() dependent
	///  upon actual type of core::Real and core::Size
	/// @remarks This is necessary for sscanf; wrong type can give wrong
	///  input.
	static String format_string_cs();

	/// @brief return a formatting string for fill_from_string() dependent
	///  upon actual type of core::Real and core::Size
	/// @remarks This is necessary for sscanf; wrong type can give wrong
	///  input.
	static String format_string_version1();

	static String format_string_residue_depth_version1();

private:
	// data

	/// @brief integer key for a residue is simply the line number in a vall file
	core::Size key_;

	/// @brief id of fragment source (e.g. pdb name)
	String id_;

	/// @brief one letter amino acid code
	char aa_;

	/// @brief one letter secondary structure code (3 letter alphabet)
	char ss_;

	/// @brief one letter secondary structure code (STR alphabet - Karplus et. al.)
	char ss_str_;

	/// @brief residue sequence number in source
	core::Size resi_;

	/// @brief b factor average for backbone
	core::Real bF_;

	/// @brief x-coordinate of C-alpha
	core::Real x_;

	/// @brief y-coordinate of C-alpha
	core::Real y_;

	/// @brief z-coordinate of C-alpha
	core::Real z_;

	/// @brief x-coordinate of C-beta
	core::Real cbx_;

	/// @brief y-coordinate of C-beta
	core::Real cby_;

	/// @brief z-coordinate of C-beta
	core::Real cbz_;

	/// @brief x-coordinate of centroid
	core::Real cenx_;

	/// @brief y-coordinate of centroid
	core::Real ceny_;

	/// @brief z-coordinate of centroid
	core::Real cenz_;

	/// @brief phi backbone torsion in degrees
	core::Real phi_;

	/// @brief psi backbone torsion in degrees
	core::Real psi_;

	/// @brief omega backbone torsion in degrees
	core::Real omega_;

	/// @brief solvent accessible area
	core::Real sa_;

	/// @brief normalized solvent accessible area
	core::Real sa_norm_;

	/// @brief phi backbone torsion in degrees from DSSP program
	core::Real dssp_phi_;

	/// @brief psi backbone torsion in degrees from DSSP program
	core::Real dssp_psi_;

	/// @brief number of alignments used for profile
	core::Size nali_;

	/// @brief per amino acid profile data
	utility::vector1<core::Real> profile_;

	/// @brief per amino acid structure profile data
	utility::vector1<core::Real> profile_struct_;

	utility::vector1<core::Real> sec_shift_data_;

	core::Real all_atom_residue_depth_;

private:
	// static data

	/// @brief order of amino acid profile data in Vall
	static utility::vector1<core::chemical::AA> order_;

};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_VallResidue_HH */
