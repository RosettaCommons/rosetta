// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_sic_dock_Rose_hh
#define INCLUDED_protocols_sic_dock_Rose_hh

#include <protocols/sic_dock/Rose.fwd.hh>

#include <platform/types.hh>
#include <utility/vector1.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <protocols/sic_dock/xyzStripeHashPose.hh>
#include <numeric/xyzTransform.hh>


namespace protocols {
namespace sic_dock {

class Rose {
public:
	typedef platform::Size Size;
	typedef platform::Real Real;
	typedef core::id::AtomID AID;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseCOP PoseCOP;
	typedef protocols::sic_dock::xyzStripeHashPose Hash;
	typedef protocols::sic_dock::xyzStripeHashPoseOP HashOP;
	typedef protocols::sic_dock::xyzStripeHashPoseCOP HashCOP;
	typedef numeric::xyzTransform<Real>  X;
	typedef numeric::xyzMatrix<Real>     M;
	typedef numeric::xyzVector<Real>     V;
	typedef X const   XC;
	typedef M const   MC;
	typedef V const   VC;
	typedef X const & XCR;
	typedef M const & MCR;
	typedef V const & VCR;
	typedef Rose const & RCR;

	PoseCOP p;
	HashCOP h;
	X       x;

	Rose(PoseCOP p);
	Rose(PoseCOP p, sic_dock::PoseCoordPickMode const & coord_picker );
	Rose(PoseCOP p, core::id::AtomID_Map<Real>  const & clash_atoms  );

	inline
	Rose(RCR r, XCR x_in) : p(r.p), h(r.h), x(x_in) {}

	virtual ~Rose() {}

	inline bool clashes (VCR point) const { return h->clash(~x*point); }
	bool clashes (RCR other) const;
	Size contacts(RCR other) const;

	friend inline bool operator==( RCR a, RCR b ){ return ( a.p==b.p && a.h==b.h && a.x.distance_squared(b.x) <= 0.000001 ); }
	friend inline bool operator!=( RCR a, RCR b ){ return ( a.p!=b.p || a.h!=b.h || a.x.distance_squared(b.x) >  0.000001 ); }

	// friend inline Rose operator +( RCR a, VCR b ){ return Rose( a, a.x+b   ); }
	inline Rose operator +(VCR b ){ return Rose( *this, x+b   ); }  // Refactoring to make it PyRosetta friendly

	// friend inline Rose operator +( VCR a, RCR b ){ return Rose( b, a  +b.x ); }

	//friend inline Rose operator -( RCR a, VCR b ){ return Rose( a, a.x-b   ); }
	inline Rose operator -( VCR b ){ return Rose( *this, x-b   ); }  // Refactoring to make it PyRosetta friendly

	//friend inline Rose operator *( RCR a, XCR b ){ return Rose( a, a.x*b   ); }
	inline Rose operator *( XCR b ){ return Rose( *this, x*b   ); }  // Refactoring to make it PyRosetta friendly

	//friend inline Rose operator *( RCR a, MCR b ){ return Rose( a, a.x*b   ); }
	inline Rose operator *( MCR b ){ return Rose( *this, x*b   ); }  // Refactoring to make it PyRosetta friendly

	// friend inline Rose operator -( VCR a, RCR b ){ return Rose( b, a  -b.x ); }
	friend inline Rose operator *( XCR a, RCR b ){ return Rose( b, a  *b.x ); }
	friend inline Rose operator *( MCR a, RCR b ){ return Rose( b, a  *b.x ); }

	// friend inline Rose operator /( RCR a, XCR b ){ return Rose( a, a.x/b   ); }
	// friend inline Rose operator /( XCR a, RCR b ){ return Rose( b, a  /b.x ); }

	inline XC res_anchor(Size const & ir) const { return x*X( p->xyz(AID(1,ir)), p->xyz(AID(2,ir)), p->xyz(AID(3,ir)) ); }
	inline XC   n_anchor() const { return res_anchor(      1       ); }
	inline XC   c_anchor() const { return res_anchor(p->n_residue()); }
	inline void align_n(XCR a){ x = a * ~X( p->xyz(AID(1,       1      )), p->xyz(AID(2,       1      )), p->xyz(AID(3,       1      )) ); }
	inline void align_c(XCR a){ x = a * ~X( p->xyz(AID(1,p->n_residue())), p->xyz(AID(2,p->n_residue())), p->xyz(AID(3,p->n_residue())) ); }

	void dump_pdb(std::ostream      & out  ) const;
	void dump_pdb(std::string const & fname) const;
	void dump_minimal_pdb(std::ostream & out, char chain='R');

	PoseCOP pose() const;

	//////////////// debug //////////////////
	bool clashes_naive(RCR other) const;
	Size contacts_naive(RCR other) const;
};



} // namespace sic_dock
} // namespace protocols

#endif
