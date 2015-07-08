// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @details
///
///
///
/// @author Bjorn Wallner

// Unit Headers
#include <protocols/jumping/DisulfPairingLibrary.hh>


// Package Headers
#include <protocols/jumping/util.hh>

// Project Headers
#include <core/types.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>

#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/JumpSRFD.hh>
#include <core/fragment/JumpingFrame.hh>

#include <basic/database/open.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2A.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/io/izstream.hh>
// #include <utility/io/ozstream.hh>
// #include <numeric/numeric.functions.hh>

#include <basic/Tracer.hh>
// #include <basic/options/option.hh>
// #include <basic/options/keys/OptionKeys.hh>

// numeric headers
#include <numeric/random/random.hh>


#include <numeric/xyzMatrix.hh>
#include <numeric/xyzMatrix.fwd.hh>
//// C++ headers
#include <cstdlib>
#include <string>

#include <protocols/jumping/DisulfPairingsList.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

static thread_local basic::Tracer tr( "protocols.jumping" );

// Singleton instance and mutex static data members
namespace utility {

using protocols::jumping::StandardDisulfPairingLibrary;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< StandardDisulfPairingLibrary >::singleton_mutex_{};
template <> std::atomic< StandardDisulfPairingLibrary * > utility::SingletonBase< StandardDisulfPairingLibrary >::instance_( 0 );
#else
template <> StandardDisulfPairingLibrary * utility::SingletonBase< StandardDisulfPairingLibrary >::instance_( 0 );
#endif

}

using core::Real;
using namespace core;
using namespace basic;
using namespace ObjexxFCL;
//using namespace basic::options;

namespace protocols {
namespace jumping {

/// @details Auto-generated virtual destructor
BaseDisulfPairingLibrary::~BaseDisulfPairingLibrary() {}

//------------------------------------------------------------------------------
// the x-axis of this coordinate system is along the p(*,2) -> p(*,1) bond vector
// the y-axis is in the p(*,*) plane, with positive dot product to p(*,2) -> p(*,3) vector

// m(*,1) is the x-axis unit vector of this new coord system expressed in
// the absolute coordinates defining the positions p(*,j)
// this corresponds to a column in the matrix m (non-CHarlie convention)
//
// thus multiplication by m can be interpreted as taking a vector in the
// local coordinate system defined by m and expressing it in the absolute
// coordinate system in which the coords p are defined
//
// by multiplication I mean non-charlie multiplication (row,col) indexing
//
// likewise, multiplication by m^t = m inverse can be thought of
// as expressing a vector given in absolute coords in terms of the
// local coordinate system defined by m


void
dis_get_coordinate_system(
	numeric::xyzMatrix_double const & p, //FArray2A_double p, // input
	numeric::xyzMatrix_double & m //FArray2A_double m // output
)
{
	using namespace numeric;

	xyzVector_double a1 = p.col_x() - p.col_y();
	xyzVector_double a2 = p.col_z() - p.col_y();
	a1.normalize();
	xyzVector_double a3 = cross( a1, a2 );
	a3.normalize();
	a2 = cross( a3, a1 );

	m = xyzMatrix_double::cols( a1, a2, a3 );
}

void
dis_get_ncac(
	FArray2A_float pos,
	numeric::xyzMatrix_double & p
)
{
	pos.dimension( 3, 5 );
	using namespace numeric;

	xyzVector_double n( &pos(1,1) );
	xyzVector_double ca( &pos(1,2) );
	xyzVector_double c( &pos(1,4) );

	p = xyzMatrix_double::cols( n, ca, c );

}

numeric::xyzMatrix_double
dis_get_ncac ( FArray2A_float pos )
{
	pos.dimension( 3, 5 );
	using namespace numeric;

	xyzVector_double n( &pos(1,1) );
	xyzVector_double ca( &pos(1,2) );
	xyzVector_double c( &pos(1,4) );

	return xyzMatrix_double::cols( n, ca, c );
}


//helper code to make an RT from two Epos
// does this live somewhere else in mini, haven't found it !
using core::kinematics::RT;
RT dis_RT_from_epos( FArray2A_float Epos1, FArray2A_float Epos2)
{
	/// rotation matrix, written in stub1 frame
	RT::Matrix rotation( 0.0 ); // 3x3
	/// tranlsation vector, written in stub1 frame
	RT::Vector translation( 0.0 ); // 3

	Size const MAX_POS( 5 ); // param::MAX_POS
	Epos1.dimension(3,MAX_POS);
	Epos2.dimension(3,MAX_POS);

	//bool const local_debug ( false );

	numeric::xyzMatrix_double p1, p2, m1, m2;

	// get coordinate systems from both residues
	dis_get_ncac(Epos1,p1);
	dis_get_ncac(Epos2,p2);
	dis_get_coordinate_system(p1,m1);
	dis_get_coordinate_system(p2,m2);

	// consider this:       |xx xy xz|
	// coordinate frame M = |yx yy yz|
	//                      |zx zy zz|
	// each column is a unit vector written in genuine frame.
	//
	// vector A in frame M can be rewritten as B in genuine frame
	// by the formula B = M x A, thus A = M^T x B
	// a simple example of this would be: the unit vector (1,0,0) in frame M
	// is actually (xx,yx,zx) in genuine frame. mathematically,
	// |xx|   |xx xy xz|   |1|
	// |yx| = |yx yy yz| x |0| ==> B = M x A
	// |zx|   |zx zy zz|   |0|
	//
	// the above formula has another layer of meaning: rotation
	// keeping the genuine frame fixed, a vector can be rotated by applying
	// matrix M onto it, e.g., (1,0,0) rotated to (xx,yx,zx)

	numeric::xyzVector_double e1( &Epos1(1,2) );
	numeric::xyzVector_double e2( &Epos2(1,2) );

	// ( e2 - e1 ) is the vector in genuine frame,
	// translation is the vector in m1 frame. so m1^T is multiplied.
	translation = m1.transposed() * ( e2 - e1 );

	// let's look at the rotation matrix
	// A, B, C are three vectors in genuine frame and are related by rotation
	// B = M1 x A; C = M2 x A;
	// ==> A = M1^T x B = M2^T x C
	// ==> C = M2 x M1^T x B
	// but note that C and B are both in genuine frame and we want a rotation
	// matrix to be applied onto a vector in M1 frame, so comes another step of
	// conversion -- left-multiply M1^T on both sides:
	// M1^T x C = M1^T x M2 x M1^T x B
	// C' = M1^T x M2 x B', as C' and B' are written in M1 frame.
	// so we get the rotation matrix as M1^T x M2.
	// but wait a minute, what Phil orginally got below is M2^T x M1 and it is
	// impossible for that to be wrong, then what happens?

	// It turns out when this rotation matrix is further applied to a vector,
	// it uses Charlies' (col,row) convention (see Dvect_multiply()
	// in RT::make_jump) which means there is one more transpose to do.
	// Now an agreement is reached:
	//  (M2^T x M1)^T = M1^T x (M2^T)^T = M1^T x M2
	// since xyzMatrix uses the normal (row,col) convention, we will switch to
	// rotation = M1^T x M2

	rotation = m1.transposed() * m2;
	/************************Phil's legacy code *********************/
	// rotation(j,*) is the j-th unit vector of 2nd coord sys written in 1st coord-sys
	// 	for ( int i=1; i<=3; ++i ) {
	// 			for ( int j=1; j<=3; ++j ) {
	// 				// DANGER: not sure about the order... ////////////////////////
	// 				// should sync with make_jump
	//  				rotation(j,i) = Ddotprod( m1(1,i), m2(1,j) ); // 2nd guess
	//  				//rotation(j,i) = Ddotprod( m1(1,j), m2(1,i) ); // 1st guess
	// 			}
	// 		}
	/************************Phil's legacy code ********************/
	RT rt;
	rt.set_translation( translation );
	rt.set_rotation( rotation );

	return rt;
}


DisulfTemplate::DisulfTemplate ( std::string const& s1, std::string const& s2, std::string const& s3 ) :
	phi  ( 2, 0.0 ),
	psi  ( 2, 0.0 ),
	omega( 2, 0.0 ),
	secstruct(2,'H')
{
	atoms_downstream_.reserve(3);
	atoms_downstream_.push_back( s1 );
	atoms_downstream_.push_back( s2 );
	atoms_downstream_.push_back( s3 );
	atoms_upstream_ = atoms_downstream_;
}

DisulfTemplate::DisulfTemplate ( std::string const& c, std::string const& s1, std::string const& s2, std::string const& s3 ) :
	phi  ( 2, 0.0 ),
	psi  ( 2, 0.0 ),
	omega( 2, 0.0 ),
	secstruct(2,'H')
{
	atoms_downstream_.reserve(4);
	atoms_downstream_.push_back( c );
	atoms_downstream_.push_back( s1 );
	atoms_downstream_.push_back( s2 );
	atoms_downstream_.push_back( s3 );
	atoms_upstream_ = atoms_downstream_;
}


///////////////////////////////////////////////////////////////////////////////
void
DisulfPairingLibrary::read_from_file( std::string const& fn)
{
	//const float MAX_NO_DIST ( 3.1 );
	std::string line,res1,res2;
	std::string tag,filename;
	int pos1,pos2,seq_sep;
	char ss1, ss2;
  //float o,p1,p2,mn_dist,mx_dist,
	//phi1,psi1,omega1,phi2,psi2,omega2;
	Size const MAX_POS( 5 ); // param::MAX_POS
	FArray2D_float Epos1(3,MAX_POS), Epos2(3,MAX_POS);
	utility::io::izstream data( fn ); //or from database file

	RT this_rt;

	while ( getline( data,line ) ) {
		std::istringstream is( line );


		is >> tag >> filename >> pos1 >> pos2 >> res1 >> res2 >>
			ss1 >> ss2 >> seq_sep >> this_rt;

		if ( is.fail() || tag != "DISULF" ) continue;

		// fill in a new pairing template
		protocols::jumping::DisulfTemplate t("CA","N","CA","C");
		t.rt_ = this_rt;

		//ANGLES NOT CURRENTLY DEFINED IN LIBRARY
		//TAKE THEM OUT OF TEMPLATE?

		//t.phi   (1) = phi1;
		//t.phi   (2) = phi2;
		//t.psi   (1) = psi1;
		//t.psi   (2) = psi2;
		//t.omega (1) = omega1;
		//t.omega (2) = omega2;

		all_pairings_.push_back( t );

		pairings_[ std::make_pair( 1, 1 ) ].push_back( t );

		++num_of_pairings_;
	}
}

/// @details puts all jump-geometries that fit the orientation and pleating into
/// list of FragData's. Try to reuse these FragData for different Frames that have same orientation and pleating
void DisulfPairingLibrary::create_jump_fragments(
	bool bWithTorsion,
	core::fragment::FragDataOPs& frags
) const {
	using namespace core::fragment;

	//	read_jump_templates(); // self-initializing
	runtime_assert( all_pairings_.size() > 0 );

	const DisulfTemplateList & templates
		( all_pairings_ );

	const int ntemplates ( templates.size() );
	const int iStart( 1 ); // in templates start residue is number 1
	const int iStop ( 2 ); // in templates stop residue is number 2
	frags.reserve( ntemplates );
	for ( DisulfTemplateList::const_iterator it=templates.begin(),	eit=templates.end();
				it!=eit; ++it ) {
		frags.push_back( core::fragment::FragDataOP( new FragData ) );
		if ( bWithTorsion ) {
			BBTorsionSRFDOP start( new BBTorsionSRFD( 3, 'E', 'X' ) );
			start->set_torsion( 1, it->phi( iStart ) );
			start->set_torsion( 2, it->psi( iStart ) );
			start->set_torsion( 3, it->omega( iStart ) );

			frags.back()->add_residue( start );
		}

		frags.back()->add_residue( SingleResidueFragDataOP( new UpJumpSRFD() ) );
		frags.back()->add_residue( SingleResidueFragDataOP( new DownJumpSRFD( it->rt_, it->atoms_downstream_, it->atoms_upstream_, 'X' ) ) );

		if ( bWithTorsion ) {
			BBTorsionSRFDOP stop( new BBTorsionSRFD( 3, 'E', 'X' ) );
			stop->set_torsion( 1, it->phi( iStop ) );
			stop->set_torsion( 2, it->psi( iStop ) );
			stop->set_torsion( 3, it->omega( iStop ) );

			frags.back()->add_residue( stop );
		}
		frags.back()->set_valid(); // yes there is data in this Fragment
	} // for-loop
} // create_jump_fragments

void
DisulfPairingLibrary::generate_jump_frags(
	DisulfPairingsList const& pairings,
	kinematics::MoveMap const& mm,
	bool bWithTorsion,
	core::fragment::FragSet& frags_accumulator
) const
{

	fragment::FragDataOPs frag_data;
	create_jump_fragments( bWithTorsion, frag_data );

	for ( Size jump_nr = 1; jump_nr <= pairings.size(); ++jump_nr) {

		//int const jump_nr ( jump_nr );
		int const startpos( pairings[ jump_nr ].pos1 );
		int const endpos( pairings[ jump_nr ].pos2 );

		if ( mm.get_bb( startpos ) && mm.get_bb( endpos ) ) {
			Size const length( bWithTorsion ? 4 : 2 );
			runtime_assert( length == frag_data.front()->size() );
			fragment::JumpingFrameOP frame = generate_empty_jump_frame( startpos, endpos, length );
			frame->add_fragment( frag_data );
			frags_accumulator.add( frame );
		} else {
			utility_exit_with_message("need to implement this: make ss-library fragments that only contain those torsions for residues\
						that can move according to movemap -- call this function with	\
						bWithTorsions = false ... and it works for now");
		}
	}
} // method


StandardDisulfPairingLibrary *
StandardDisulfPairingLibrary::create_singleton_instance()
{
	StandardDisulfPairingLibrary * instance = new StandardDisulfPairingLibrary;
	//instance = new StandardDisulfPairingLibrary();
	std::cout << "READING START" << std::endl;
	instance->read_from_file( basic::database::full_name("sampling/disulfide_jump_database_wip.dat") );
	std::cout << "READING END" << std::endl;
	return instance;
}

} // jumping
} // protocols
