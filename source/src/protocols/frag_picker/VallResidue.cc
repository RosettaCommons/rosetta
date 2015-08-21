// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/VallResidue.cc
/// @brief  class for managing a line of the Vall fragment library
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/frag_picker/VallResidue.hh>

// project headers
#include <core/chemical/AA.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <protocols/frag_picker/Faraggi_SA.hh>

// boost headers
#include <boost/type_traits.hpp>

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/exit.hh>

// C++ headers
#include <cstdio>
#include <iostream>
#include <sstream>

// Tracer
#include <basic/Tracer.hh>

namespace protocols {
namespace frag_picker {

static thread_local basic::Tracer TR( "protocols.frag_picker.VallResidue" );

// static initialization
utility::vector1<core::chemical::AA> VallResidue::order_ = order_vector();

/// @brief default constructor
VallResidue::VallResidue() :
	utility::pointer::ReferenceCount(), key_(0),id_(""),profile_(20),profile_struct_(20),sec_shift_data_(/* 0 */),
	position_index_(0), section_index_(0) {
}

/// @brief string constructor
VallResidue::VallResidue(String const & line) {
	if ( line.length() > 300 ) {
		fill_from_string_version1(line);
	} else if ( line.length() > 240 ) {
		fill_from_string_cs(line);
	} else if ( line.length() > 110 ) {
		fill_from_string(line);
	} else {
		fill_from_string_residue_depth_version1(line);
	}
}

/// @brief copy constructor
VallResidue::VallResidue(VallResidue const & rval) :
	utility::pointer::ReferenceCount(),
	key_(rval.key_),
	id_(rval.id_),
	aa_(rval.aa_),
	ss_(rval.ss_),
	ss_str_(rval.ss_str_),
	resi_(rval.resi_),
	bF_(rval.bF_),
	x_(rval.x_),
	y_(rval.y_),
	z_(rval.z_),
	cbx_(rval.cbx_),
	cby_(rval.cby_),
	cbz_(rval.cbz_),
	cenx_(rval.cenx_),
	ceny_(rval.ceny_),
	cenz_(rval.cenz_),
	phi_(rval.phi_),
	psi_(rval.psi_),
	omega_(rval.omega_),
	sa_(rval.sa_),
	sa_norm_(rval.sa_norm_),
	dssp_phi_(rval.dssp_phi_),
	dssp_psi_(rval.dssp_psi_),
	nali_(rval.nali_),
	profile_(rval.profile_),
	profile_struct_(rval.profile_struct_),
	sec_shift_data_(rval.sec_shift_data_),
	position_index_(rval.position_index_),
	section_index_(rval.section_index_),
	all_atom_residue_depth_(rval.all_atom_residue_depth_){
}

/// @brief default destructor
VallResidue::~VallResidue() {
}

/// @brief copy assignment
VallResidue & VallResidue::operator =(VallResidue const & rval) {
	if ( this != &rval ) {
		id_ = rval.id_;
		aa_ = rval.aa_;
		ss_ = rval.ss_;
		ss_str_ = rval.ss_str_;
		resi_ = rval.resi_;
		bF_ = rval.bF_;
		x_ = rval.x_;
		y_ = rval.y_;
		z_ = rval.z_;
		cbx_ = rval.cbx_;
		cby_ = rval.cby_;
		cbz_ = rval.cbz_;
		cenx_ = rval.cenx_;
		ceny_ = rval.ceny_;
		cenz_ = rval.cenz_;
		dssp_phi_ = rval.dssp_phi_;
		dssp_psi_ = rval.dssp_psi_;
		sa_ = rval.sa_;
		sa_norm_ = rval.sa_norm_;
		nali_ = rval.nali_;
		phi_ = rval.phi_;
		psi_ = rval.psi_;
		omega_ = rval.omega_;
		profile_ = rval.profile_;
		profile_struct_ = rval.profile_struct_;
		sec_shift_data_ = rval.sec_shift_data_;
		position_index_ = rval.position_index_;
		section_index_ = rval.section_index_;
	}
	return *this;
}

/// @brief build a BBTorsionSRFD from this page
VallResidue::BBTorsionSRFDOP VallResidue::bbtorsion_srfd() const {
	BBTorsionSRFDOP srfd( new BBTorsionSRFD() );

	srfd->set_sequence(aa_);
	srfd->set_secstruct(ss_);

	// currently, phi = 1, psi = 2, omega = 3
	// see BBTorsionSRFD and Pose implemenations for details
	srfd->set_torsion(1, phi_);
	srfd->set_torsion(2, psi_);
	srfd->set_torsion(3, omega_);
	srfd->set_coordinates( x_, y_, z_ );

	return srfd;
}

/// @brief fill internal data from string
/// @details Values are delimited by whitespace.  Ordering is:
///  <tt> id  aa  ss  resi  dummy  dummy  x  y  z  phi  psi  omega  dummy  dummy  dummy  dummy  (aa profile_info, 20 columns) BFactor </tt>
void VallResidue::fill_from_string(String const & line) {

	char id[] = { '\0', '\0', '\0', '\0', '\0', '\0' }; // 5 char + 1 null termination

	static String format_here = format_string();

	// Use sscanf here; Vall is huge and istringstream is way too slow.
	// *One* sscanf call per line, multiple calls decreases performance!
	std::sscanf(line.c_str(), format_here.c_str(), &id, &aa_, &ss_, &resi_, &x_,
		&y_, &z_, &phi_, &psi_, &omega_, &profile_[order_[1]],
		&profile_[order_[2]], &profile_[order_[3]], &profile_[order_[4]],
		&profile_[order_[5]], &profile_[order_[6]], &profile_[order_[7]],
		&profile_[order_[8]], &profile_[order_[9]], &profile_[order_[10]],
		&profile_[order_[11]], &profile_[order_[12]],
		&profile_[order_[13]], &profile_[order_[14]],
		&profile_[order_[15]], &profile_[order_[16]],
		&profile_[order_[17]], &profile_[order_[18]],
		&profile_[order_[19]], &profile_[order_[20]]);
	ss_str_ = ss_; // this vall format does not use an STR alphabet so just set as ss_
	id_ = id;
}

/// @brief fill internal data from string
/// @details Values are delimited by whitespace.  Ordering is:
///  <tt> id  aa  ss  resi  dummy  dummy  x  y  z  phi  psi  omega  dummy  dummy  dummy  dummy  (aa profile_info, 20 columns) (CS data, 12 columns)</tt>
void VallResidue::fill_from_string_cs(String const & line) {

	char id[] = { '\0', '\0', '\0', '\0', '\0', '\0' }; // 5 char + 1 null termination
	//std::cerr<<line<<std::endl;
	sec_shift_data_.resize(12);

	// Use sscanf here; Vall is huge and istringstream is way too slow.
	// *One* sscanf call per line, multiple calls decreases performance!

	static String format_here = format_string_cs();

	std::sscanf(line.c_str(), format_here.c_str(), &id, &aa_, &ss_, &resi_, &bF_,
		&x_, &y_, &z_, &phi_, &psi_, &omega_, &profile_[order_[1]],
		&profile_[order_[2]], &profile_[order_[3]], &profile_[order_[4]],
		&profile_[order_[5]], &profile_[order_[6]], &profile_[order_[7]],
		&profile_[order_[8]], &profile_[order_[9]], &profile_[order_[10]],
		&profile_[order_[11]], &profile_[order_[12]],
		&profile_[order_[13]], &profile_[order_[14]],
		&profile_[order_[15]], &profile_[order_[16]],
		&profile_[order_[17]], &profile_[order_[18]],
		&profile_[order_[19]], &profile_[order_[20]],
		&sec_shift_data_[1], &sec_shift_data_[2], &sec_shift_data_[3],
		&sec_shift_data_[4], &sec_shift_data_[5], &sec_shift_data_[6],
		&sec_shift_data_[7], &sec_shift_data_[8], &sec_shift_data_[9],
		&sec_shift_data_[10], &sec_shift_data_[11], &sec_shift_data_[12]);
	ss_str_ = ss_; // this vall format does not use an STR alphabet so just set as ss_
	id_ = id;
}

/// @brief fill internal data from string
/// @details Values are delimited by whitespace.  Ordering is:
///  <tt> id  aa  ss  resi  dummy  dummy  x  y  z cbz cby cbz cenx ceny cenz phi  psi  omega dssp_phi dssp_psi sa nali  (aa profile_info, 20 columns) BFactor </tt>
void VallResidue::fill_from_string_version1(String const & line) {

	char id[] = { '\0', '\0', '\0', '\0', '\0', '\0' }; // 5 char + 1 null termination

	static String format_here = format_string_version1();

	// Use sscanf here; Vall is huge and istringstream is way too slow.
	// *One* sscanf call per line, multiple calls decreases performance!
	std::sscanf(line.c_str(), format_here.c_str(), &id, &aa_, &ss_str_, &resi_, &bF_, &x_,
		&y_, &z_, &cbx_, &cby_, &cbz_, &cenx_, &ceny_, &cenz_, &phi_, &psi_, &omega_, &dssp_phi_, &dssp_psi_, &sa_, &nali_, &profile_[order_[1]],
		&profile_[order_[2]], &profile_[order_[3]], &profile_[order_[4]],
		&profile_[order_[5]], &profile_[order_[6]], &profile_[order_[7]],
		&profile_[order_[8]], &profile_[order_[9]], &profile_[order_[10]],
		&profile_[order_[11]], &profile_[order_[12]],
		&profile_[order_[13]], &profile_[order_[14]],
		&profile_[order_[15]], &profile_[order_[16]],
		&profile_[order_[17]], &profile_[order_[18]],
		&profile_[order_[19]], &profile_[order_[20]],
		&profile_struct_[order_[1]],&profile_struct_[order_[2]],
		&profile_struct_[order_[3]], &profile_struct_[order_[4]],
		&profile_struct_[order_[5]], &profile_struct_[order_[6]],
		&profile_struct_[order_[7]], &profile_struct_[order_[8]],
		&profile_struct_[order_[9]], &profile_struct_[order_[10]],
		&profile_struct_[order_[11]], &profile_struct_[order_[12]],
		&profile_struct_[order_[13]], &profile_struct_[order_[14]],
		&profile_struct_[order_[15]], &profile_struct_[order_[16]],
		&profile_struct_[order_[17]], &profile_struct_[order_[18]],
		&profile_struct_[order_[19]], &profile_struct_[order_[20]]);

	sa_norm_ = sa_/protocols::frag_picker::sa_faraggi_max(aa_);
	set_ss_from_str();
	id_ = id;
}

void VallResidue::fill_from_string_residue_depth_version1( String const & line ) {

	char id[] = { '\0', '\0', '\0', '\0', '\0', '\0' }; // 5 char + 1 null termination

	static String format_here = format_string_residue_depth_version1();

	std::sscanf(line.c_str(), format_here.c_str(), &id, &aa_, &resi_, &x_,
		&y_, &z_, &all_atom_residue_depth_);
	id_ = id;
	ss_ = 'L'; // placeholder
	ss_str_ = 'C'; // placeholder
}


/// @brief return a vector specifying the order of profile data in Vall
utility::vector1<core::chemical::AA> VallResidue::order_vector() {
	utility::vector1<core::chemical::AA> order(20);

	order[1] = core::chemical::aa_from_oneletter_code('A');
	order[2] = core::chemical::aa_from_oneletter_code('C');
	order[3] = core::chemical::aa_from_oneletter_code('D');
	order[4] = core::chemical::aa_from_oneletter_code('E');
	order[5] = core::chemical::aa_from_oneletter_code('F');
	order[6] = core::chemical::aa_from_oneletter_code('G');
	order[7] = core::chemical::aa_from_oneletter_code('H');
	order[8] = core::chemical::aa_from_oneletter_code('I');
	order[9] = core::chemical::aa_from_oneletter_code('K');
	order[10] = core::chemical::aa_from_oneletter_code('L');
	order[11] = core::chemical::aa_from_oneletter_code('M');
	order[12] = core::chemical::aa_from_oneletter_code('N');
	order[13] = core::chemical::aa_from_oneletter_code('P');
	order[14] = core::chemical::aa_from_oneletter_code('Q');
	order[15] = core::chemical::aa_from_oneletter_code('R');
	order[16] = core::chemical::aa_from_oneletter_code('S');
	order[17] = core::chemical::aa_from_oneletter_code('T');
	order[18] = core::chemical::aa_from_oneletter_code('V');
	order[19] = core::chemical::aa_from_oneletter_code('W');
	order[20] = core::chemical::aa_from_oneletter_code('Y');

	return order;
}

/// @brief return a formatting string for fill_from_string() dependent
///  upon actual type of core::Real and core::Size
/// @remarks This is necessary for sscanf; wrong type can give wrong
///  input.
VallResidue::String VallResidue::format_string() {
	using boost::is_same;

	// resolve types
	bool const is_ulong = is_same<Size, unsigned long>::value;
	bool const is_double = is_same<Real, double>::value;

	// format w/ the following order:
	std::ostringstream s;
	s << "%s" << " "; // (string) id
	s << "%c" << " "; // (char)   aa
	s << "%c" << " "; // (char)   ss
	s << (is_ulong ? "%lu" : "%u") << " "; // (int)    resi
	s << "%*s" << " "; // (string) dummy
	s << "%*s" << " ";// (string) dummy
	s << (is_double ? "%lf" : "%f") << " "; // (real)   x
	s << (is_double ? "%lf" : "%f") << " "; // (real)   y
	s << (is_double ? "%lf" : "%f") << " "; // (real)   z
	s << (is_double ? "%lf" : "%f") << " "; // (real)   phi
	s << (is_double ? "%lf" : "%f") << " "; // (real)   psi
	s << (is_double ? "%lf" : "%f") << " "; // (real)   omega
	s << "%*s" << " "; // (string) dummy
	s << "%*s" << " "; // (string) dummy
	s << "%*s" << " "; // (string) dummy
	s << "%*s"; // (string) dummy

	// (real)   (aa profile_info, 20 columns)
	for ( Size i = 0; i < 20; ++i ) {
		s << " " << (is_double ? "%lf" : "%f");
	}

	TR.Debug << "vall line format (no CS) is:\n"<<s.str()<<std::endl;

	return s.str();
}

VallResidue::String VallResidue::format_string_cs() {
	using boost::is_same;

	// resolve types
	bool const is_ulong = is_same<Size, unsigned long>::value;
	bool const is_double = is_same<Real, double>::value;

	// format w/ the following order:
	std::ostringstream s;
	s << "%s" << " "; // (string) id
	s << "%c" << " "; // (char)   aa
	s << "%c" << " "; // (char)   ss
	s << (is_ulong ? "%lu" : "%u") << " "; // (int)    resi
	s << "%*s" << " "; // (string) dummy
	s << "%*s" << " ";// (string) dummy
	s << (is_double ? "%lf" : "%f") << " "; // (real) b factor
	s << (is_double ? "%lf" : "%f") << " "; // (real)   x
	s << (is_double ? "%lf" : "%f") << " "; // (real)   y
	s << (is_double ? "%lf" : "%f") << " "; // (real)   z
	s << (is_double ? "%lf" : "%f") << " "; // (real)   phi
	s << (is_double ? "%lf" : "%f") << " "; // (real)   psi
	s << (is_double ? "%lf" : "%f") << " "; // (real)   omega

	// (real)   (aa profile_info, 20 columns)
	for ( Size i = 0; i < 20; ++i ) {
		s << " " << (is_double ? "%lf" : "%f");
	}

	// (real)   (CS data, 12 columns)
	for ( Size i = 0; i < 12; ++i ) {
		s << " " << (is_double ? "%lf" : "%f");
	}


	TR.Debug << "vall line format (WITH CS) is:\n"<<s.str()<<std::endl;
	return s.str();
}

/// @brief return a formatting string for fill_from_string_version1() dependent
///  upon actual type of core::Real and core::Size
/// @remarks This is necessary for sscanf; wrong type can give wrong
///  input.
VallResidue::String VallResidue::format_string_version1() {
	using boost::is_same;

	// resolve types
	bool const is_ulong = is_same<Size, unsigned long>::value;
	bool const is_double = is_same<Real, double>::value;

	// format w/ the following order:
	std::ostringstream s;
	s << "%s" << " "; // (string) id
	s << "%c" << " "; // (char)   aa
	s << "%c" << " "; // (char)   ss
	s << (is_ulong ? "%lu" : "%u") << " "; // (int)    resi
	s << (is_double ? "%lf" : "%f") << " "; // (real)   avg bfactor
	s << (is_double ? "%lf" : "%f") << " "; // (real)   x
	s << (is_double ? "%lf" : "%f") << " "; // (real)   y
	s << (is_double ? "%lf" : "%f") << " "; // (real)   z
	s << (is_double ? "%lf" : "%f") << " "; // (real)   cbx
	s << (is_double ? "%lf" : "%f") << " "; // (real)   cby
	s << (is_double ? "%lf" : "%f") << " "; // (real)   cbz
	s << (is_double ? "%lf" : "%f") << " "; // (real)   cenx
	s << (is_double ? "%lf" : "%f") << " "; // (real)   ceny
	s << (is_double ? "%lf" : "%f") << " "; // (real)   cenz
	s << (is_double ? "%lf" : "%f") << " "; // (real)   phi
	s << (is_double ? "%lf" : "%f") << " "; // (real)   psi
	s << (is_double ? "%lf" : "%f") << " "; // (real)   omega
	s << (is_double ? "%lf" : "%f") << " "; // (real)   dssp_phi
	s << (is_double ? "%lf" : "%f") << " "; // (real)   dssp_psi
	s << (is_double ? "%lf" : "%f") << " "; // (real) sa
	s << (is_ulong ? "%ld" : "%d") << " "; // (size) nali

	// (real)   (aa profile_info, 20 columns)
	for ( Size i = 0; i < 20; ++i ) {
		s << " " << (is_double ? "%lf" : "%f");
	}

	// (real)  (aa structure profile_info, 20 columns)
	for ( Size i = 0; i < 20; ++i ) {
		s << " " << (is_double ? "%lf" : "%f");
	}

	TR.Debug << "vall line format version1 is:\n"<<s.str()<<std::endl;

	return s.str();
}

VallResidue::String VallResidue::format_string_residue_depth_version1() {
	using boost::is_same;

	// resolve types
	bool const is_ulong = is_same<Size, unsigned long>::value;
	bool const is_double = is_same<Real, double>::value;

	// format w/ the following order:
	std::ostringstream s;
	s << "%s" << " "; // (string) id
	s << "%c" << " "; // (char)   aa
	s << (is_ulong ? "%lu" : "%u") << " "; // (int)    resi
	s << "%*s" << " "; // (string) dummy    actual pdb resi
	s << (is_double ? "%lf" : "%f") << " "; // (real)   x
	s << (is_double ? "%lf" : "%f") << " "; // (real)   y
	s << (is_double ? "%lf" : "%f") << " "; // (real)   z
	s << (is_double ? "%lf" : "%f") << " "; // (real)   all-atom residue depth

	TR.Debug << "vall line format residue depth  version1 is:\n"<<s.str()<<std::endl;

	return s.str();
}

core::Real VallResidue::distance_squared( VallResidueCOP r ) {
	core::Real x = x_-r->x();
	core::Real y = y_-r->y();
	core::Real z = z_-r->z();
	return x*x+y*y+z*z;
}

core::Real VallResidue::distance_squared_cb( VallResidueCOP r ) {
	// use CA for Gly
	if ( aa_ == 'G' ) {
		if ( r->aa() == 'G' ) {
			return distance_squared( r );
		} else {
			core::Real x = x_-r->cbx();
			core::Real y = y_-r->cby();
			core::Real z = z_-r->cbz();
			return x*x+y*y+z*z;
		}
	} else if ( r->aa() == 'G' ) {
		core::Real x = cbx_-r->x();
		core::Real y = cby_-r->y();
		core::Real z = cbz_-r->z();
		return x*x+y*y+z*z;
	} else {
		core::Real x = cbx_-r->cbx();
		core::Real y = cby_-r->cby();
		core::Real z = cbz_-r->cbz();
		return x*x+y*y+z*z;
	}
}

core::Real VallResidue::distance_squared_cen( VallResidueCOP r ) {
	core::Real x = cenx_-r->cenx();
	core::Real y = ceny_-r->ceny();
	core::Real z = cenz_-r->cenz();
	return x*x+y*y+z*z;
}

core::Real VallResidue::distance_squared( VallResidueCOP r, ContactType const & type ) {
	switch (type) {
	case CA :
		return distance_squared(r);
	case CB :
		return distance_squared_cb(r);
	case CEN :
		return distance_squared_cen(r);
	default :
		utility_exit_with_message("unkown contacts type for distance");
	}
}

} // frag_picker
} // protocols
