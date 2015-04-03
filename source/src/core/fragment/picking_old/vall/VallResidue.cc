// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/vall/VallResidue.cc
/// @brief  class for managing a line of the Vall fragment library
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/fragment/picking_old/vall/VallResidue.hh>

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

#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {


// static initialization
utility::vector1< core::chemical::AA > VallResidue::order_ = order_vector();
VallResidue::String VallResidue::format_ = format_string();
VallResidue::String VallResidue::format_2011_ = format_string_2011();

/// @brief default constructor
VallResidue::VallResidue() :
	aa_( 0 ),
	ss_( 0 ),
	resi_( 0 ),
	x_( 0.0 ),
	y_( 0.0 ),
	z_( 0.0 ),
	phi_( 0.0 ),
	psi_( 0.0 ),
	omega_( 0.0 ),
	profile_( 20 ),
	position_index_( 0 ),
	section_index_( 0 )
{}


/// @brief string constructor
VallResidue::VallResidue( String const & line ) {
	fill_from_string( line );
}


/// @brief copy constructor
VallResidue::VallResidue( VallResidue const & rval ) :
	id_( rval.id_ ),
	aa_( rval.aa_ ),
	ss_( rval.ss_ ),
	resi_( rval.resi_ ),
	x_( rval.x_ ),
	y_( rval.y_ ),
	z_( rval.z_ ),
	phi_( rval.phi_ ),
	psi_( rval.psi_ ),
	omega_( rval.omega_ ),
	profile_( rval.profile_ ),
	position_index_( rval.position_index_ ),
	section_index_( rval.section_index_ )
{}


/// @brief default destructor
VallResidue::~VallResidue() {}


/// @brief copy assignment
VallResidue & VallResidue::operator =( VallResidue const & rval ) {
	if ( this != &rval ) {
		id_ = rval.id_;
		aa_ = rval.aa_;
		ss_ = rval.ss_;
		resi_ = rval.resi_;
		x_ = rval.x_;
		y_ = rval.y_;
		z_ = rval.z_;
		phi_ = rval.phi_;
		psi_ = rval.psi_;
		omega_ = rval.omega_;
		profile_ = rval.profile_;
		position_index_ = rval.position_index_;
		section_index_ = rval.section_index_;
	}
	return *this;
}


/// @brief build a BBTorsionSRFD of the given type from this page
/// @param[in] srfd_type BBTorsionSRFD::create() will be called from this object.
/// @return A BBTorsionSRFD of the given type initialized with the backbone
///  torsion information from this page.
VallResidue::BBTorsionSRFDOP VallResidue::bbtorsion_srfd( BBTorsionSRFD const & srfd_type ) const {
	BBTorsionSRFDOP srfd = utility::pointer::static_pointer_cast< BBTorsionSRFD > ( srfd_type.create() );

	srfd->set_sequence( aa_ );
	srfd->set_secstruct( ss_ );

	// currently, phi = 1, psi = 2, omega = 3
	// see BBTorsionSRFD and Pose implemenations for details
	srfd->set_torsion( 1, phi_ );
	srfd->set_torsion( 2, psi_ );
	srfd->set_torsion( 3, omega_ );

	return srfd;
}


/// @brief fill internal data from string
/// @details Values are delimited by whitespace.  Ordering is:
///  <tt> id  aa  ss  resi  dummy  dummy  x  y  z  phi  psi  omega  dummy  dummy  dummy  dummy  (aa profile_info, 20 columns) </tt>
void VallResidue::fill_from_string( String const & line ) {

	char id[] = { '\0', '\0', '\0', '\0', '\0', '\0' }; // 5 char + 1 null termination

	// Use sscanf here; Vall is huge and istringstream is way too slow.
	// *One* sscanf call per line, multiple calls decreases performance!
  //2006 Vall line is ~224 units char long
	//2011 Vall line is ~396 units long
	if(line.size() < 230){
		std::sscanf(
		line.c_str(),
		format_.c_str(),
		&id, &aa_, &ss_, &resi_,
		&x_, &y_, &z_,
		&phi_, &psi_, &omega_,
		&profile_[ order_[ 1 ] ],
		&profile_[ order_[ 2 ] ],
		&profile_[ order_[ 3 ] ],
		&profile_[ order_[ 4 ] ],
		&profile_[ order_[ 5 ] ],
		&profile_[ order_[ 6 ] ],
		&profile_[ order_[ 7 ] ],
		&profile_[ order_[ 8 ] ],
		&profile_[ order_[ 9 ] ],
		&profile_[ order_[ 10 ] ],
		&profile_[ order_[ 11 ] ],
		&profile_[ order_[ 12 ] ],
		&profile_[ order_[ 13 ] ],
		&profile_[ order_[ 14 ] ],
		&profile_[ order_[ 15 ] ],
		&profile_[ order_[ 16 ] ],
		&profile_[ order_[ 17 ] ],
		&profile_[ order_[ 18 ] ],
		&profile_[ order_[ 19 ] ],
		&profile_[ order_[ 20 ] ]
	);
	}
	if(line.size() > 300){
			Real bfactor;//all definded here are unused
		  Real cb_x;
			Real cb_y;
		  Real cb_z;
			Real cen_x;
		  Real cen_y;
			Real cen_z;
			Real dssp_phi;
			Real dssp_psi;
			Real dssp_sa; 
			Real n_ali; 
			utility::vector1< Real > profile_s(20);
			std::sscanf(
				line.c_str(),
				format_2011_.c_str(),
				&id,&aa_, &ss_, &resi_,&bfactor,
				&x_,&y_,&z_,
				&cb_x,&cb_y,&cb_z,
				&cen_x,&cen_y,&cen_z,
				&phi_, &psi_, &omega_,
				&dssp_phi, &dssp_psi, &dssp_sa, &n_ali,
				&profile_[ order_[ 1 ] ],
				&profile_[ order_[ 2 ] ],
				&profile_[ order_[ 3 ] ],
				&profile_[ order_[ 4 ] ],
				&profile_[ order_[ 5 ] ],
				&profile_[ order_[ 6 ] ],
				&profile_[ order_[ 7 ] ],
				&profile_[ order_[ 8 ] ],
				&profile_[ order_[ 9 ] ],
				&profile_[ order_[ 10 ] ],
				&profile_[ order_[ 11 ] ],
				&profile_[ order_[ 12 ] ],
				&profile_[ order_[ 13 ] ],
				&profile_[ order_[ 14 ] ],
				&profile_[ order_[ 15 ] ],
				&profile_[ order_[ 16 ] ],
				&profile_[ order_[ 17 ] ],
				&profile_[ order_[ 18 ] ],
				&profile_[ order_[ 19 ] ],
				&profile_[ order_[ 20 ] ],
				&profile_s[ order_[ 1 ] ],
				&profile_s[ order_[ 2 ] ],
				&profile_s[ order_[ 3 ] ],
				&profile_s[ order_[ 4 ] ],
				&profile_s[ order_[ 5 ] ],
				&profile_s[ order_[ 6 ] ],
				&profile_s[ order_[ 7 ] ],
				&profile_s[ order_[ 8 ] ],
				&profile_s[ order_[ 9 ] ],
				&profile_s[ order_[ 10 ] ],
				&profile_s[ order_[ 11 ] ],
				&profile_s[ order_[ 12 ] ],
				&profile_s[ order_[ 13 ] ],
				&profile_s[ order_[ 14 ] ],
				&profile_s[ order_[ 15 ] ],
				&profile_s[ order_[ 16 ] ],
				&profile_s[ order_[ 17 ] ],
				&profile_s[ order_[ 18 ] ],
				&profile_s[ order_[ 19 ] ],
				&profile_s[ order_[ 20 ] ]
		);
	}
	id_ = id;

}


/// @brief return a vector specifying the order of profile data in Vall
utility::vector1< core::chemical::AA > VallResidue::order_vector() {
	utility::vector1< core::chemical::AA > order( 20 );

	order[ 1] = core::chemical::aa_from_oneletter_code( 'A' );
	order[ 2] = core::chemical::aa_from_oneletter_code( 'C' );
	order[ 3] = core::chemical::aa_from_oneletter_code( 'D' );
	order[ 4] = core::chemical::aa_from_oneletter_code( 'E' );
	order[ 5] = core::chemical::aa_from_oneletter_code( 'F' );
	order[ 6] = core::chemical::aa_from_oneletter_code( 'G' );
	order[ 7] = core::chemical::aa_from_oneletter_code( 'H' );
	order[ 8] = core::chemical::aa_from_oneletter_code( 'I' );
	order[ 9] = core::chemical::aa_from_oneletter_code( 'K' );
	order[10] = core::chemical::aa_from_oneletter_code( 'L' );
	order[11] = core::chemical::aa_from_oneletter_code( 'M' );
	order[12] = core::chemical::aa_from_oneletter_code( 'N' );
	order[13] = core::chemical::aa_from_oneletter_code( 'P' );
	order[14] = core::chemical::aa_from_oneletter_code( 'Q' );
	order[15] = core::chemical::aa_from_oneletter_code( 'R' );
	order[16] = core::chemical::aa_from_oneletter_code( 'S' );
	order[17] = core::chemical::aa_from_oneletter_code( 'T' );
	order[18] = core::chemical::aa_from_oneletter_code( 'V' );
	order[19] = core::chemical::aa_from_oneletter_code( 'W' );
	order[20] = core::chemical::aa_from_oneletter_code( 'Y' );

	return order;
}


/// @brief return a formatting string for fill_from_string() dependent
///  upon actual type of core::Real and core::Size
/// @remarks This is necessary for sscanf; wrong type can give wrong
///  input.
VallResidue::String VallResidue::format_string() {
	using boost::is_same;

	// resolve types
	bool const is_ulong = is_same< Size, unsigned long >::value;
	bool const is_double = is_same< Real, double >::value;

	// format w/ the following order:
	std::ostringstream s;
	s << "%s" << " "; // (string) id
	s << "%c" << " "; // (char)   aa
	s << "%c" << " "; // (char)   ss
	s << ( is_ulong ? "%lu" : "%u" ) << " "; // (int)    resi
	s << "%*s" << " "; // (string) dummy
	s << "%*s" << " ";// (string) dummy
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   x
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   y
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   z
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   phi
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   psi
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   omega
	s << "%*s" << " "; // (string) dummy
	s << "%*s" << " "; // (string) dummy
	s << "%*s" << " "; // (string) dummy
	s << "%*s"; // (string) dummy

	// (real)   (aa profile_info, 20 columns)
	for ( Size i = 0; i < 20; ++i ) {
		s << " " << ( is_double ? "%lf" : "%f" );
	}

	return s.str();
}
/// @brief return a formatting string for fill_from_string() dependent
///  upon actual type of core::Real and core::Size 
//   This version is compatible with 2011 vall
/// @remarks This is necessary for sscanf; wrong type can give wrong
///  input.
VallResidue::String VallResidue::format_string_2011() {
	using boost::is_same;

	// resolve types
	bool const is_ulong = is_same< Size, unsigned long >::value;
	bool const is_double = is_same< Real, double >::value;
	// format w/ the following order:
	std::ostringstream s;
	s << "%s" << " "; // (string) id
	s << "%c" << " "; // (char)   aa
	s << "%c" << " "; // (char)   ss
	s << ( is_ulong ? "%lu" : "%u" ) << " "; // (int)    resi
  s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   bfactor
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   ca_x
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   ca_y
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   ca_z
  s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   cb_x
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   cb_y
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   cb_z
  s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   cen_x
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   cen_y
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   cen_z
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   phi
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   psi
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   omega
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   dssp_phi
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   dssp_psi
	s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   dssp_sa
  s << ( is_double ? "%lf" : "%f" ) << " "; // (real)   n_ali
	// (real)   (aa profile_info, 20 columns)
	for ( Size i = 0; i < 40; ++i ) {
		s << " " << ( is_double ? "%lf" : "%f" );
	}
	return s.str();
}

} // vall
} // picking_old
} // fragment
} // core
