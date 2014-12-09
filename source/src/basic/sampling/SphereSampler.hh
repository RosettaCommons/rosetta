// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_basic_sampling_SphereSampler_hh
#define INCLUDED_basic_sampling_SphereSampler_hh

#include <basic/database/open.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <numeric/xyzVector.hh>
#include <numeric/types.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>
#include <set>

namespace basic {
namespace sampling {

typedef numeric::Real Real;

class SphereNode {
private:
	int const level_;
	numeric::xyzVector<Real> const axis_;
	utility::vector1<SphereNode const *> parents_;
	utility::vector1<SphereNode const *> children_;
	utility::vector1<SphereNode const *> neighbors_;
public:
	SphereNode(int level, numeric::xyzVector<Real> axis)
		:	level_(level), axis_(axis)
	{
		neighbors_.reserve(6);
	}
	void reserve_child() { children_.reserve(7); }
	void add_child   (SphereNode const * c) { children_.push_back(c); }
	void add_neighbor(SphereNode const * n) { neighbors_.push_back(n); }
	void add_parent  (SphereNode const * p) { parents_.push_back(p); }
	numeric::Size num_children()    const { return children_.size(); }
	numeric::Size num_neighbor() const { return neighbors_.size(); }
	numeric::Size num_parents()  const { return parents_.size(); }
	SphereNode const & child   (numeric::Size i) const { return *(children_[i]); }
	SphereNode const & neighbor(numeric::Size i) const { return *(neighbors_[i]); }
	SphereNode const & parent  (numeric::Size i) const { return *(parents_[i]); }
	// SphereNode const * child_ptr   (numeric::Size i) const { return children_[i]; }
	// SphereNode const * neighbor_ptr(numeric::Size i) const { return neighbors_[i]; }
	int level() const { return level_; }
	numeric::xyzVector<Real> const & axis() const { return axis_; }
};

class SphereSampler {
private:
	utility::vector1<SphereNode*> allnodes_;
public:
	SphereSampler() : allnodes_(40974, (SphereNode*)(NULL) ){
		read_sphere_data();
		assert( sanity_check() );
	}
	~SphereSampler() {
		for(utility::vector1<SphereNode*>::iterator i = allnodes_.begin(); i != allnodes_.end(); ++i) {
			delete *i;
		}
	}

	void read_sphere_data() {
		utility::io::izstream in;
		if(!basic::database::open(in,"sampling/spheres/sphere_hierarchy.dat.gz") || !in.good()) {
			utility_exit_with_message("problem opening sampling/spheres/sphere_hierarchy.dat.gz");
		}
		utility::vector1<int> parent(40974);
		utility::vector1<utility::vector1<int> > children(40974,utility::vector1<int>(0)),neighbors(40974);
		for(int i = 1; i <= 40974; ++i) {
			int l,p;
			numeric::xyzVector<Real> a;
			if(!(in >> l >> p >> a.x() >> a.y() >> a.z())) utility_exit_with_message("problem with sampling/spheres/sphere_hierarchy.dat.gz");
			assert(1 <= l && l <= 7);
			assert(fabs(a.length()-1.0) < 0.00001);
			allnodes_[i] = new SphereNode(l,a);
			assert((1 <= p && p <= 40974) || p == 65535);  // I assume the parentheses should be this way? ~Labonte
			parent[i] = p;
			utility::vector1<int> c(7);
			if(!(in >> c[1])) utility_exit_with_message("problem with sampling/spheres/sphere_hierarchy.dat.gz");
			if(c[1] != 65535) {
				if(!(in >> c[2] >> c[3] >> c[4] >> c[5] >> c[6] >> c[7])) utility_exit_with_message("problem with sampling/spheres/sphere_hierarchy.dat.gz");
				if( c.back() == 65535 ) c.pop_back();
				children[i] = c;
			}
			utility::vector1<int> n(6);
			if(!(in >> n[1] >> n[2] >> n[3] >> n[4] >> n[5] >> n[6])) utility_exit_with_message("problem with sampling/spheres/sphere_hierarchy.dat.gz");
			if( n.back() == 65535 ) n.pop_back();
			neighbors[i] = n;
		}
		in.close();
		for(int i = 1; i <= 40974; ++i) {
			// neighbors
			for(utility::vector1<int>::iterator ic = neighbors[i].begin(); ic != neighbors[i].end(); ++ic) {
				allnodes_[i]->add_neighbor( allnodes_[*ic] );
			}
			// level 7 nodes are leaves
			if(allnodes_[i]->level() == 7 ) continue;
			// children
			allnodes_[i]->reserve_child();
			for(utility::vector1<int>::iterator ic = children[i].begin(); ic != children[i].end(); ++ic) {
				allnodes_[ i ]->add_child ( allnodes_[*ic] );
				allnodes_[*ic]->add_parent( allnodes_[ i ] );
			}
		}
	}

	inline Real covering_radius_degrees(int level) {
		switch(level) {
			case 1: return 37.378;
			case 2: return 21.0; // probably a little lower
			case 3: return 11.685637449452;
			case 4: return  5.814062961918;
			case 5: return  2.902260274250;
			case 6: return  1.449774552065;
			case 7: return  0.724517158473;
			default: utility_exit_with_message("bad level");
		}
	}

	void pdb_from_level(int level, std::string fname) {
		using namespace ObjexxFCL::format;
		utility::io::ozstream out(fname);
		int count = 0;
		for(numeric::Size i = 1; i <= num_sample(level); ++i) {
			SphereNode const & n( sample(level,i) );
			numeric::xyzVector<Real> const & p(n.axis());
			out<<"HETATM"<<I(5,++count)<<' '<<"NODE"<<' '<<"SPH"<<' '<<"A"<<I(4,1)<<"    ";
			out<<F(8,3,100.0*p.x())<<F(8,3,100.0*p.y())<<F(8,3,100.0*p.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
		}
		out.close();
	}

	bool sanity_check() {
		int n5=0,n6=0;
		for(utility::vector1<SphereNode*>::iterator i = allnodes_.begin(); i != allnodes_.end(); ++i) {
			SphereNode const & n(**i);
			if(n.level()==7) {
				assert(n.num_children()==0);
			} else {
				assert(n.num_children()==6 || n.num_children()==7);
				n6 += n.num_children()==6;
			}
			assert( n.num_neighbor() == 5 || n.num_neighbor() == 6);
			n5 += n.num_neighbor()==5;
			for(numeric::Size in = 1; in <= n.num_neighbor(); ++in) {
				SphereNode const & b(n.neighbor(in));
				assert(n.level()==b.level());
				bool is_neighbor = false;
				for(numeric::Size ib = 1; ib <= b.num_neighbor(); ++ib) {
					is_neighbor |= &b.neighbor(ib) == &n;
				}
				assert(is_neighbor);
			}
			for(numeric::Size ic = 1; ic <= n.num_children(); ++ic) {
				SphereNode const & c(n.child(ic));
				assert(n.level()+1==c.level());
				bool is_parent = false;
				for(numeric::Size ip = 1; ip <= c.num_parents(); ++ip) {
					is_parent |= &c.parent(ip) == &n;
				}
				assert(is_parent);
			}
			for(numeric::Size ip = 1; ip <= n.num_parents(); ++ip) {
				SphereNode const & p(n.parent(ip));
				assert(n.level()-1==p.level());
				bool is_child = false;
				for(numeric::Size ic = 1; ic <= p.num_children(); ++ic) {
					is_child |= &p.child(ic) == &n;
				}
				assert(is_child);
			}
		}
		for(numeric::Size l = 1; l <= 7; ++l) {
			for(numeric::Size i = 1; i <= num_sample(l); ++i) {
				assert( sample(l,i).level() == (int)l);
			}
		}
		assert( n5==7*12 );
		assert( n6==6*12 );
		return true;
	}

	inline numeric::Size num_sample(int level) const {
		switch(level) {
			case 1: return    12;
			case 2: return    32;
			case 3: return   122;
			case 4: return   482;
			case 5: return  1922;
			case 6: return  7682;
			case 7: return 30722;
			default: utility_exit_with_message("bad level");
		}
		return 0; // never
	}
	inline SphereNode const & sample(int l, int i) const { return *allnodes_[i+sample_offset(l)]; }
	numeric::Size sample_offset(int level) const {
		switch(level) {
			case 1: return     0;
			case 2: return    12;
			case 3: return    44;
			case 4: return   166;
			case 5: return   648;
			case 6: return  2570;
			case 7: return 10252;
			default: utility_exit_with_message("bad level");
		}
		return 0; // never
	}

	int num_sample1() const { return num_sample(1); }
	int num_sample2() const { return num_sample(2); }
	int num_sample3() const { return num_sample(3); }
	int num_sample4() const { return num_sample(4); }
	int num_sample5() const { return num_sample(5); }
	int num_sample6() const { return num_sample(6); }
	int num_sample7() const { return num_sample(7); }

	SphereNode const & sample1(int i) { return sample(1,i); }
	SphereNode const & sample2(int i) { return sample(2,i); }
	SphereNode const & sample3(int i) { return sample(3,i); }
	SphereNode const & sample4(int i) { return sample(4,i); }
	SphereNode const & sample5(int i) { return sample(5,i); }
	SphereNode const & sample6(int i) { return sample(6,i); }
	SphereNode const & sample7(int i) { return sample(7,i); }

};






}
}

#endif
