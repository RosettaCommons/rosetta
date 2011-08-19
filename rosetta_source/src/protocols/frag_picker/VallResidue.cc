// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/frag_picker/VallResidue.cc
/// @brief  class for managing a line of the Vall fragment library
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/frag_picker/VallResidue.hh>

// project headers
#include <core/chemical/AA.hh>
#include <core/fragment/BBTorsionSRFD.hh>

// boost headers
#include <boost/type_traits.hpp>

// utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <cstdio>
#include <iostream>
#include <sstream>

// Tracer
#include <basic/Tracer.hh>

namespace protocols {
namespace frag_picker {

static basic::Tracer TR("protocols.frag_picker.VallResidue");

// static initialization
utility::vector1<core::chemical::AA> VallResidue::order_ = order_vector();
VallResidue::String VallResidue::format_ = format_string_cs();

/// @brief default constructor
VallResidue::VallResidue() :
	utility::pointer::ReferenceCount(), key_(0),id_(""),profile_(20), sec_shift_data_(0),
			position_index_(0), section_index_(0) {
}

/// @brief string constructor
VallResidue::VallResidue(String const & line) {
	if (line.length() < 240)
		fill_from_string(line);
	else
		fill_from_string_cs(line);

}

/// @brief copy constructor
VallResidue::VallResidue(VallResidue const & rval) :
	utility::pointer::ReferenceCount(), key_(rval.key_),id_(rval.id_), aa_(rval.aa_), ss_(
			rval.ss_), resi_(rval.resi_), bF_(rval.bF_), x_(rval.x_), y_(rval.y_),
			z_(rval.z_), phi_(rval.phi_), psi_(rval.psi_), omega_(rval.omega_),
			profile_(rval.profile_), sec_shift_data_(rval.sec_shift_data_),
	position_index_(rval.position_index_), section_index_(rval.section_index_) {
}

/// @brief default destructor
VallResidue::~VallResidue() {
}

/// @brief copy assignment
VallResidue & VallResidue::operator =(VallResidue const & rval) {
	if (this != &rval) {
		id_ = rval.id_;
		aa_ = rval.aa_;
		ss_ = rval.ss_;
		resi_ = rval.resi_;
		bF_ = rval.bF_;
		x_ = rval.x_;
		y_ = rval.y_;
		z_ = rval.z_;
		phi_ = rval.phi_;
		psi_ = rval.psi_;
		omega_ = rval.omega_;
		profile_ = rval.profile_;
		sec_shift_data_ = rval.sec_shift_data_;
		position_index_ = rval.position_index_;
		section_index_ = rval.section_index_;
	}
	return *this;
}

/// @brief build a BBTorsionSRFD from this page
VallResidue::BBTorsionSRFDOP VallResidue::bbtorsion_srfd() const {
	BBTorsionSRFD * srfd = new BBTorsionSRFD();

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

	id_ = id;
}

/// @brief fill internal data from string
/// @details Values are delimited by whitespace.  Ordering is:
///  <tt> id  aa  ss  resi  dummy  dummy  x  y  z  phi  psi  omega  dummy  dummy  dummy  dummy  (aa profile_info, 20 columns) (CS data, 12 columns)</tt>
void VallResidue::fill_from_string_cs(String const & line) {

	char id[] = { '\0', '\0', '\0', '\0', '\0', '\0' }; // 5 char + 1 null termination
	//std::cerr<<line<<std::endl;
	sec_shift_data_.resize(12);

	//format_ = format_string_cs();

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

	id_ = id;
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
	for (Size i = 0; i < 20; ++i) {
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
	for (Size i = 0; i < 20; ++i) {
		s << " " << (is_double ? "%lf" : "%f");
	}

	// (real)   (CS data, 12 columns)
	for (Size i = 0; i < 12; ++i) {
		s << " " << (is_double ? "%lf" : "%f");
	}


        TR.Debug << "vall line format (WITH CS) is:\n"<<s.str()<<std::endl;
	return s.str();
}

core::Real VallResidue::distance_squared( VallResidueCOP r ) {
	core::Real x = x_-r->x();
	core::Real y = y_-r->y();
	core::Real z = z_-r->z();
	return x*x+y*y+z*z;
}

} // frag_picker
} // protocols
