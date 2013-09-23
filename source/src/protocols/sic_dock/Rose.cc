// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:

#include <protocols/sic_dock/Rose.hh>

#include <protocols/sic_dock/util.hh>

#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/sasa.hh>
#include <core/pose/util.hh>

namespace protocols {
namespace sic_dock {

using platform::Size;
using platform::Real;
using std::string;
using utility::vector1;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;
using ObjexxFCL::format::RJ;
using numeric::min;
using numeric::max;
using std::cout;
using std::cerr;
using std::endl;
typedef numeric::xyzVector<platform::Real> Vec;
typedef numeric::xyzMatrix<platform::Real> Mat;

Rose::Rose(PoseCOP pin                                                   ) : p(pin),h(new Hash(*p,         BB ,4.0)) {}
Rose::Rose(PoseCOP pin, sic_dock::PoseCoordPickMode const & coord_picker ) : p(pin),h(new Hash(*p,coord_picker,4.0)) {}
Rose::Rose(PoseCOP pin, core::id::AtomID_Map<Real>  const & clash_atoms  ) : p(pin),h(new Hash(*p,clash_atoms ,4.0)) {}



bool Rose::clashes(RCR o) const {
	// Real const & thresh_2 = h->grid_size2();
	RCR b( h->natom() >  o.h->natom() ? *this : o );
	RCR s( h->natom() <= o.h->natom() ? *this : o );
	// XC s2b( b.h->translation() + ( (s.x-s.h->translation()) / b.x ) );
	XC s2b( b.h->translation() + ( ~b.x * (s.x-s.h->translation()) ) );
	for(Hash::const_iterator i = s.h->begin(); i != s.h->end(); ++i){
		if( b.h->clash_raw( s2b * *i ) ) return true;
	}
	return false;
}

Size Rose::contacts(RCR o) const {
	Size count = 0;
	// Real const & thresh_2 = h->grid_size2();
	RCR b( h->natom() >  o.h->natom() ? *this : o );
	RCR s( h->natom() <= o.h->natom() ? *this : o );
	// (~d)*n
	// XC s2b( b.h->translation() + ( (s.x-s.h->translation()) / b.x ) );
	XC s2b( b.h->translation() + ( ~b.x * (s.x-s.h->translation()) ) );
	for(Hash::const_iterator i = s.h->begin(); i != s.h->end(); ++i){
		count += b.h->nbcount_raw( s2b * *i );
	}
	return count;
}

core::pose::PoseCOP Rose::pose() const {
	core::pose::PoseOP pose = new core::pose::Pose(*p);
	for(Size ir = 1; ir <= pose->n_residue(); ++ir) {
		for(Size ia = 1; ia <= pose->residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose->set_xyz( aid, x.xform(pose->xyz(aid)) );
		}
	}
	return pose;
}

void Rose::dump_pdb(std::ostream & out) const {
	pose()->dump_pdb(out);
}

void Rose::dump_pdb(std::string const & fname) const {
	pose()->dump_pdb(fname);
}

void Rose::dump_minimal_pdb(std::ostream & out, char chain){
	for(Size ir = 1; ir <= p->n_residue(); ++ir){
		V v;
		v = x * p->xyz(AID(1,ir)); out<<"ATOM  "<<I(5,3*ir-2)<<' '<<"  N "<<' '<<"GLY"<<' '<<chain<<I(4,ir)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
		v = x * p->xyz(AID(2,ir)); out<<"ATOM  "<<I(5,3*ir-1)<<' '<<" CA "<<' '<<"GLY"<<' '<<chain<<I(4,ir)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
		v = x * p->xyz(AID(3,ir)); out<<"ATOM  "<<I(5,3*ir-0)<<' '<<"  C "<<' '<<"GLY"<<' '<<chain<<I(4,ir)<<"    "<<F(8,3,v.x())<<F(8,3,v.y())<<F(8,3,v.z())<<F(6,2,1.0)<<F(6,2,1.0)<<endl;
	}
}


//////////////// debugging ///////////////////////

bool Rose::clashes_naive(RCR o) const {
	Real const thresh_2 = h->grid_size2();
	for(Hash::const_iterator i = o.h->begin(); i != o.h->end(); ++i){
		VC u( o.x.xform(*i-o.h->translation()) );
		for(Hash::const_iterator j = h->begin(); j != h->end(); ++j){
			VC v( x.xform(*j-h->translation()) );
			if( u.distance_squared(v) < thresh_2 ) return true;
		}
	}
	return false;
}

Size Rose::contacts_naive(RCR o) const {
	Size count = 0;
	Real const thresh_2 = h->grid_size2();
	for(Hash::const_iterator i = o.h->begin(); i != o.h->end(); ++i){
		VC u( o.x.xform(*i-o.h->translation()) );
		for(Hash::const_iterator j = h->begin(); j != h->end(); ++j){
			VC v( x.xform(*j-h->translation()) );
			if( u.distance_squared(v) <= thresh_2 ) count++;
		}
	}
	return count;
}

} // namespace sic_dock
} // namespace protocols
