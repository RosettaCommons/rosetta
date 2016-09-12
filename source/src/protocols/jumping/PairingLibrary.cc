// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @details
///
/// @author Bjorn Wallner
/// @author Christopher Miles (cmiles@uw.edu)

// Unit Headers
#include <protocols/jumping/PairingLibrary.hh>

// Package Headers
#include <protocols/jumping/util.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/id/NamedStubID.hh>

#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/JumpSRFD.hh>
#include <core/fragment/JumpingFrame.hh>
#include <core/scoring/dssp/PairingsList.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2A.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/thread/threadsafe_creation.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyz.functions.hh>

// C++ headers
#include <cstdlib>
#include <string>

// Singleton instance and mutex static data members
namespace utility {

using protocols::jumping::StandardPairingLibrary;

#if defined MULTI_THREADED
template <> std::mutex utility::SingletonBase< StandardPairingLibrary >::singleton_mutex_{};
template <> std::atomic< StandardPairingLibrary * > utility::SingletonBase< StandardPairingLibrary >::instance_( 0 );
#else
template <> StandardPairingLibrary * utility::SingletonBase< StandardPairingLibrary >::instance_( nullptr );
#endif

}

static THREAD_LOCAL basic::Tracer tr( "protocols.jumping" );

using core::Real;
using namespace core;
using namespace basic;
using namespace ObjexxFCL;

namespace protocols {
namespace jumping {

BasePairingLibrary::~BasePairingLibrary() = default;

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
get_coordinate_system(
	numeric::xyzMatrix_double const & p, // input
	numeric::xyzMatrix_double & m        // output
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
get_ncac(
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
get_ncac ( FArray2A_float pos )
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
RT RT_from_epos( FArray2A_float Epos1, FArray2A_float Epos2)
{
	/// rotation matrix, written in stub1 frame
	RT::Matrix rotation( 0.0 ); // 3x3
	/// tranlsation vector, written in stub1 frame
	RT::Vector translation( 0.0 ); // 3

	Size const MAX_POS( 5 ); // param::MAX_POS
	Epos1.dimension(3,MAX_POS);
	Epos2.dimension(3,MAX_POS);

	numeric::xyzMatrix_double p1, p2, m1, m2;

	// get coordinate systems from both residues
	get_ncac(Epos1,p1);
	get_ncac(Epos2,p2);
	get_coordinate_system(p1,m1);
	get_coordinate_system(p2,m2);

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
	RT rt;
	rt.set_translation( translation );
	rt.set_rotation( rotation );

	return rt;
}


PairingTemplate::PairingTemplate ( std::string const& s1, std::string const& s2, std::string const& s3 ) :
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

PairingTemplate::PairingTemplate ( std::string const& c, std::string const& s1, std::string const& s2, std::string const& s3 ) :
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
PairingLibrary::read_from_file( std::string const& fn)
{
	const float MAX_NO_DIST ( 3.1 );
	std::string line;
	std::string tag,filename;
	int pos1,pos2;
	float o,p1,p2,mn_dist,mx_dist,
		phi1,psi1,omega1,phi2,psi2,omega2;
	Size const MAX_POS( 5 ); // param::MAX_POS
	FArray2D_float Epos1(3,MAX_POS), Epos2(3,MAX_POS);
	utility::io::izstream data( fn ); //or from database file

	while ( getline( data,line ) ) {
		std::istringstream is( line );
		Vector n1, ca1, c1;
		Vector n2, ca2, c2;
		is >> tag >> filename >> pos1 >> pos2 >> mn_dist >> mx_dist >>
			o >> p1 >> p2 >>
			Epos1(1,1) >> Epos1(2,1) >> Epos1(3,1) >>
			Epos1(1,2) >> Epos1(2,2) >> Epos1(3,2) >>
			Epos1(1,4) >> Epos1(2,4) >> Epos1(3,4) >>
			Epos2(1,1) >> Epos2(2,1) >> Epos2(3,1) >>
			Epos2(1,2) >> Epos2(2,2) >> Epos2(3,2) >>
			Epos2(1,4) >> Epos2(2,4) >> Epos2(3,4) >>

			phi1 >> psi1 >> omega1 >>
			phi2 >> psi2 >> omega2;

		RT this_rt(RT_from_epos(Epos1,Epos2));

		if ( is.fail() || tag != "PAIR" ) continue;

		runtime_assert ( pos1 < pos2 && p1 * p2 > 0.0 &&
			std::abs(phi1) < 185 && std::abs(psi1) < 185 && std::abs(omega1) < 185 &&
			std::abs(phi2) < 185 && std::abs(psi2) < 185 && std::abs(omega2) < 185 );

		// filter
		// note that the filename contains info about the scop class so you could
		// in principle filter for b,c,or d class proteins individually
		if ( mx_dist > MAX_NO_DIST ||
				phi1 >  0.0 || phi2 >  0.0 ||
				psi1 < 50.0 || psi2 < 50.0 ||
				std::abs( omega1 ) < 90 ||
				std::abs( omega2 ) < 90 ) continue;

		// fill in a new beta-pairing template
		PairingTemplate t("CA","N","CA","C");
		t.rt_ = this_rt;
		t.phi   (1) = phi1;
		t.phi   (2) = phi2;
		t.psi   (1) = psi1;
		t.psi   (2) = psi2;
		t.omega (1) = omega1;
		t.omega (2) = omega2;

		this_rt.reverse();
		PairingTemplate t_reverse("CA","N","CA","C");
		t_reverse.rt_ = this_rt;
		t_reverse.phi   (1) = phi2; //reverse... also put torsions on the right side of jump
		t_reverse.phi   (2) = phi1;
		t_reverse.psi   (1) = psi2;
		t_reverse.psi   (2) = psi1;
		t_reverse.omega (1) = omega2;
		t_reverse.omega (2) = omega1;
		const int o_key( ( o  < 0.0 ) ? 1 : 2 ); // orientation
		const int p_key( ( p1 < 0.0 ) ? 1 : 2 ); // pleating
		Vector dNN = n1-ca1;
		Vector dNC = n1-c1;
		Vector dCCA = c1-ca1;
		pairings_[ std::make_pair( o_key, p_key ) ].push_back( t );
		pairings_[ std::make_pair( o_key, p_key ) ].push_back( t_reverse );

		//each "view" gives a different result, it cannot be decided a-priori which one will be more appropriate
		++num_of_pairings_;
	}
}

///////////////////////////////////////////////////////////////////////////////
void
PairingLibrary::read_from_file_no_filters( std::string const& fn)
{
	std::string line;
	std::string tag,filename;
	int pos1,pos2;
	float o,p1,p2,mn_dist,mx_dist,
		phi1,psi1,omega1,phi2,psi2,omega2;
	Size const MAX_POS( 5 ); // param::MAX_POS
	FArray2D_float Epos1(3,MAX_POS), Epos2(3,MAX_POS);
	utility::io::izstream data( fn ); //or from database file
	std::ofstream template_infofile("jump_TMH_templates.dat.info");
	while ( getline( data,line ) ) {
		std::istringstream is( line );
		Vector n1, ca1, c1;
		Vector n2, ca2, c2;
		is >> tag >> filename >> pos1 >> pos2 >> mn_dist >> mx_dist >>
			o >> p1 >> p2 >>

			n1.x() >> n1.y() >> n1.z() >>
			ca1.x() >> ca1.y() >> ca1.z() >>
			c1.x() >> c1.y() >> c1.z() >>

			n2.x() >> n2.y() >> n2.z() >>
			ca2.x() >> ca2.y() >> ca2.z() >>
			c2.x() >> c2.y() >> c2.z() >>

			phi1 >> psi1 >> omega1 >>
			phi2 >> psi2 >> omega2;

		RT this_rt( kinematics::Stub( ca1, n1, ca1, c1 ), kinematics::Stub( ca2, n2, ca2, c2) );

		// fill in a new beta-pairing template
		char ss1,ss2;
		if ( p1 == 'E' || p1 == 1 ) {
			ss1 = 'E';
		} else if ( p1 == 'H' || p1 == 2 ) {
			ss1 = 'H';
		} else if ( p1 == 'L' || p1 == 3 ) {
			ss1 = 'L';
		} else {
			std::cout << "bad secstruct: " << p1 << std::endl;
			continue;
		}
		if ( p2 == 'E' || p2 == 1 ) {
			ss2 = 'E';
		} else if ( p2 == 'H' || p2 == 2 ) {
			ss2 = 'H';
		} else if ( p2 == 'L' || p2 == 3 ) {
			ss2 = 'L';
		} else {
			std::cout << "bad secstruct: " << p2 << std::endl;
			continue;
		}

		template_infofile << fn << ' ' << this_rt << "\n";

		PairingTemplate t("CA","N","CA","C");
		t.rt_ = this_rt;
		t.phi   (1) = phi1;
		t.phi   (2) = phi2;
		t.psi   (1) = psi1;
		t.psi   (2) = psi2;
		t.omega (1) = omega1;
		t.omega (2) = omega2;
		t.secstruct(1)=ss1;
		t.secstruct(2)=ss2;

		// bw for TMH the jump library are specific to the positions. these are defined in the template file to 1 or 2.
		if ( pos1==0 && pos2 == 1 ) {
			const int o_key = (int)o;
			const int p_key = 0; //p1;

			// This is for generic use in case pos1,pos2 from pairings_file is not defined in the jump library.
			pairings_[ std::make_pair( o_key, p_key ) ].push_back( t );
			++num_of_pairings_;
		} else {
			const int o_key = (int)o;
			const int p_key = 0; //p1;
			std::cout << t.rt_ << "\n";

			// This is when we want to test a number of different jump for a particular position.
			pairings_[ std::make_pair( pos1, pos2 ) ].push_back( t );

			//Put all in generic library as well...
			pairings_[ std::make_pair( o_key, p_key ) ].push_back( t );
			++num_of_pairings_;
			++num_of_pairings_;
		}
	}
	template_infofile.close();
}

///////////////////////////////////////////////////////////////////////////////
kinematics::RT
PairingLibrary::get_random_beta_sheet_jump(
	int const orientation,
	int const pleating
) const
{
	runtime_assert( pairings_.size() > 0 );

	// key for looking up the template geometry:
	std::pair<int,int> key( orientation, pleating );

	// HACK to get it to compile -- fix this later
	runtime_assert ( pairings_.count( key ) == 1 );

	const PairingTemplateList & templates
		( pairings_.find( key )->second );

	const int ntemplates ( templates.size() );

	int const index( static_cast<int>( numeric::random::rg().uniform() * ntemplates ) );
	const PairingTemplate &t ( templates[ index ] );

	return t.rt_;
}

///////////////////////////////////////////////////////////////////////////////
kinematics::RT
PairingLibrary::get_random_tmh_jump(int const orientation,
	int const pos1,
	int const pos2
) const
{
	assert( pairings_.size() > 0 );

	// key for looking up the template geometry:
	std::pair<int,int> generic_key (orientation,0);
	std::pair<int,int> specific_key (pos1,pos2);
	// std::pair<int,int> key (specific_key); // Unused variable causes warning.

	const PairingTemplateList & templates
		( pairings_.find( specific_key )->second );

	const int ntemplates ( templates.size() );
	if ( ntemplates>0 ) {
		int const index( static_cast<int>( numeric::random::rg().uniform() * ntemplates ) );
		const PairingTemplate &t ( templates[ index ] );
		return t.rt_;
	} else { // use the generic key
		std::cout << "No key found for " << pos1 << ' ' << pos2 << " using the generic template\n";
		const PairingTemplateList & templates_generic( pairings_.find( generic_key )->second );
		const int ntemplates_generic ( templates_generic.size() );
		int const index( static_cast<int>( numeric::random::rg().uniform() * ntemplates_generic) );
		const PairingTemplate &t ( templates[ index ] );
		return t.rt_;
	}
}

///////////////////////////////////////////////////////////////////////////////
void
PairingLibrary::set_tmh_jump(core::pose::Pose pose,
	int const jump_number,
	int const orientation,
	int const pos1,
	int const pos2
) const
{
	assert( pairings_.size() > 0 );

	// key for looking up the template geometry:
	std::pair<int,int> generic_key (orientation,0);
	std::pair<int,int> specific_key (pos1,pos2);
	// std::pair<int,int> key (specific_key); // Unused variable causes warning.

	const PairingTemplateList & templates
		( pairings_.find( specific_key )->second );

	const int ntemplates ( templates.size() );
	if ( ntemplates>0 ) {
		int const index( static_cast<int>( numeric::random::rg().uniform() * ntemplates ) );
		const PairingTemplate &t ( templates[ index ] );
		std::cout << jump_number << "\n";
		std::cout << t.rt_;

		pose.set_phi(pos1,t.phi(1));
		pose.set_phi(pos2,t.phi(2));
		pose.set_psi(pos1,t.psi(1));
		pose.set_psi(pos2,t.psi(2));
		pose.set_omega(pos1,t.omega(1));
		pose.set_omega(pos2,t.omega(2));
		pose.set_secstruct(pos1,t.secstruct(1));
		pose.set_secstruct(pos2,t.secstruct(2));

		id::StubID up_stub(   core::pose::named_stub_id_to_stub_id( core::id::NamedStubID( "CA","N","CA","C", pos1 ), pose ) );
		id::StubID down_stub( core::pose::named_stub_id_to_stub_id( core::id::NamedStubID( "CA","N","CA","C", pos2 ), pose ) );
		pose.conformation().set_stub_transform( up_stub, down_stub, t.rt_ );
	} else { // use the generic key
		std::cout << "No key found for " << pos1 << ' ' << pos2 << " using the generic template\n";
		const PairingTemplateList & templates_generic( pairings_.find( generic_key )->second );
		const int ntemplates_generic ( templates_generic.size() );
		int const index( static_cast<int>( numeric::random::rg().uniform() * ntemplates_generic) );
		const PairingTemplate &t ( templates[ index ] );
		std::cout << jump_number << "\n";
		std::cout << t.rt_;

		pose.set_phi(pos1,t.phi(1));
		pose.set_phi(pos2,t.phi(2));
		pose.set_psi(pos1,t.psi(1));
		pose.set_psi(pos2,t.psi(2));
		pose.set_omega(pos1,t.omega(1));
		pose.set_omega(pos2,t.omega(2));
		pose.set_secstruct(pos1,t.secstruct(1));
		pose.set_secstruct(pos2,t.secstruct(2));

		id::StubID up_stub(   core::pose::named_stub_id_to_stub_id( core::id::NamedStubID( "CA","N","CA","C", pos1 ), pose ) );
		id::StubID down_stub( core::pose::named_stub_id_to_stub_id( core::id::NamedStubID( "CA","N","CA","C", pos2 ), pose ) );
		pose.conformation().set_stub_transform( up_stub, down_stub, t.rt_ );
	}
}

/// @details puts all jump-geometries that fit the orientation and pleating into
/// list of FragData's. Try to reuse these FragData for different Frames that have same orientation and pleating
void PairingLibrary::create_jump_fragments(
	int const orientation,
	int const pleating,
	bool bWithTorsion,
	core::fragment::FragDataOPs& frags
) const {
	using namespace core::fragment;

	// read_jump_templates(); // self-initializing
	runtime_assert( pairings_.size() > 0 );

	// key for looking up the template geometry:
	std::pair<int,int> key( orientation, pleating );

	// HACK to get it to compile -- fix this later
	runtime_assert ( pairings_.count( key ) == 1 );

	const PairingTemplateList & templates
		( pairings_.find( key )->second );

	const int ntemplates ( templates.size() );
	const int iStart( 1 ); // in templates start residue is number 1
	const int iStop ( 2 ); // in templates stop residue is number 2
	frags.reserve( ntemplates );
	for ( auto const & it : templates ) {
		frags.push_back( core::fragment::FragDataOP( new FragData ) );
		if ( bWithTorsion ) {
			BBTorsionSRFDOP start( new BBTorsionSRFD( 3, 'E', 'X' ) );
			start->set_torsion( 1, it.phi( iStart ) );
			start->set_torsion( 2, it.psi( iStart ) );
			start->set_torsion( 3, it.omega( iStart ) );

			frags.back()->add_residue( start );
		}

		frags.back()->add_residue( SingleResidueFragDataOP( new UpJumpSRFD() ) );
		frags.back()->add_residue( SingleResidueFragDataOP( new DownJumpSRFD( it.rt_, it.atoms_downstream_, it.atoms_upstream_, 'X' ) ) );

		if ( bWithTorsion ) {
			BBTorsionSRFDOP stop( new BBTorsionSRFD( 3, 'E', 'X' ) );
			stop->set_torsion( 1, it.phi( iStop ) );
			stop->set_torsion( 2, it.psi( iStop ) );
			stop->set_torsion( 3, it.omega( iStop ) );

			frags.back()->add_residue( stop );
		}
		frags.back()->set_valid(); // yes there is data in this Fragment
	} // for-loop
} // create_jump_fragments

void
PairingLibrary::generate_jump_frags(
	core::scoring::dssp::PairingsList const& pairings,
	kinematics::MoveMap const& mm,
	bool bWithTorsion,
	core::fragment::FragSet& frags_accumulator
) {
	runtime_assert( has_orientation_and_pleating( pairings ) );

	// find out how many different kind of fragments are we interested in:
	// max of four: A 1 , A 2, P 1, P 2
	typedef utility::vector1< Size > JumpList;
	typedef std::map< std::pair< Size, Size >, JumpList > JumpOrientations;
	JumpOrientations jump_kind;
	Size jump_nr ( 1 );
	for ( auto const & pairing : pairings ) {
		Size o_key ( pairing.Orientation() ); // < 0 ? 1 : 2 );
		Size p_key ( pairing.Pleating() ); // < 0 ? 1 : 2 );
		jump_kind[ std::make_pair( o_key, p_key ) ].push_back( jump_nr++ );
	}

	// now generate fragments for each of the maximum four JumpOrientations present
	for ( JumpOrientations::const_iterator it=jump_kind.begin(), eit=jump_kind.end();
			it!=eit;
			++it ) {
		Size o_key( it->first.first ); //orientation
		Size p_key( it->first.second ); //pleating ... believe me or not, it is in first.second
		fragment::FragDataOPs frag_data;
		create_jump_fragments( o_key, p_key, bWithTorsion, frag_data );
		for ( int jump_nr : it->second ) {
			int const startpos( pairings[ jump_nr ].Pos1() );
			int const endpos( pairings[ jump_nr ].Pos2() );

			if ( mm.get_bb( startpos ) && mm.get_bb( endpos ) ) {
				Size const length( bWithTorsion ? 4 : 2 );
				runtime_assert( length == frag_data.front()->size() );
				fragment::JumpingFrameOP frame = generate_empty_jump_frame( startpos, endpos, length );
				frame->add_fragment( frag_data );
				frags_accumulator.add( frame );
			} else {
				utility_exit_with_message("need to implement this: make ss-library fragments that only contain those torsions for residues "
					"that can move according to movemap -- call this function with\t"
					"bWithTorsions = false ... and it works for now");
			}
		} // for JumpList iteration
	} // loop over orientations and pleatings
} // method


StandardPairingLibrary *
StandardPairingLibrary::create_singleton_instance()
{
	auto * instance = new StandardPairingLibrary;
	instance->read_from_file( basic::database::full_name("scoring/score_functions/jump_templates_SSpairs_v2.dat") );
	return instance;
}

} // jumping
} // protocols
