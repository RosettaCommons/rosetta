// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/sampling/orientations/QuaternionGrid.fwd.hh
/// @brief  provide an evenly spaced set of rotations
/// @author Will Sheffler <willsheffler@gmail.com>

// much of this comes from http://charles.karney.info/orientation/

#ifndef INCLUDED_basic_sampling_orientations_QuaternionGrid_hh
#define INCLUDED_basic_sampling_orientations_QuaternionGrid_hh

// unit headers
#include <basic/sampling/orientations/QuaternionGrid.fwd.hh>

// C++ headers
#include <iosfwd>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/types.hh>

// C++ headers
#include <map>

namespace basic {
namespace sampling {
namespace orientations {


// Minimal quaternion class -- isn't there a quaternion class in src/numeric?
class Quaternion {
public:
	numeric::Real w, x, y, z;
	Quaternion(numeric::Real ww = 1, numeric::Real xx = 0, numeric::Real yy = 0, numeric::Real zz = 0)
		: w(ww)    , x(xx)    , y(yy)    , z(zz) {}
	void Normalize();
	void Canonicalize() {
		Normalize();
		// Make first biggest element positive
		numeric::Real mag = w;
		if (std::abs(x) > std::abs(mag))      mag = x;
		if (std::abs(y) > std::abs(mag))      mag = y;
		if (std::abs(z) > std::abs(mag))      mag = z;
		if (mag < 0) {      w *= -1;      x *= -1;      y *= -1;      z *= -1;    }
		return;
	}
	// a.Times(b) returns a * b
	Quaternion Times(const Quaternion& q) const {
		numeric::Real
			mw = w*q.w - x*q.x - y*q.y - z*q.z,
			mx = w*q.x + x*q.w + y*q.z - z*q.y,
			my = w*q.y + y*q.w + z*q.x - x*q.z,
			mz = w*q.z + z*q.w + x*q.y - y*q.x;
		return Quaternion(mw, mx, my, mz);
	}
	void Print(std::ostream& s) const;
	numeric::xyzVector<numeric::Real> euler() const;
	void PrintEuler(std::ostream& s) const;
	numeric::xyzMatrix<numeric::Real> rotation_matrix() const ;
};

// Class to hold a set of orientations and weights
class QuatSet {
public:
	Quaternion Orientation(size_t i) const {
		return m_v[i];
	}
	numeric::Real Weight(size_t i) const {
		return m_w[i];
	}
	size_t Number() const {
		return m_v.size();
	}
	void Add(const Quaternion& q, numeric::Real w = 1) {
		Quaternion v(q);
		v.Canonicalize();
		m_v.push_back(v);
		m_w.push_back(w);
	}
	void Clear() {
		m_v.clear();
		m_w.clear();
	}
	void Print(std::ostream& s, bool euler = false, std::size_t prec = 6) const;
private:
	std::vector<Quaternion> m_v;
	std::vector<numeric::Real> m_w;
};

// The triple of grid indices
class Triple {
public:
	long a, b, c;
	Triple(long aa, long bb, long cc) :
		a(aa),
		b(bb),
		c(cc) {}
};

// Generate the permutations and sign changes for a Triple.
class Permute {
public:
	Permute(Triple x) {
		assert(x.a >= x.b && x.b >= x.c && x.c >= 0);
		m_arr.push_back(x);
		size_t n = 1;
		// Do the sign changes
		if (x.a != 0) {
			for (size_t i = 0; i < n; ++i)
				m_arr.push_back(Triple(-m_arr[i].a, m_arr[i].b, m_arr[i].c));
			n *= 2;
		}
		if (x.b != 0) {
			for (size_t i = 0; i < n; ++i)
				m_arr.push_back(Triple(m_arr[i].a, -m_arr[i].b, m_arr[i].c));
			n *= 2;
		}
		if (x.c != 0) {
			for (size_t i = 0; i < n; ++i)
				m_arr.push_back(Triple(m_arr[i].a, m_arr[i].b, -m_arr[i].c));
			n *= 2;
		}
		if (x.a == x.b && x.b == x.c)
			return;
		// With at least two distinct indices we can rotate the set thru 3
		// permuations.
		for (size_t i = 0; i < n; ++i) {
			m_arr.push_back(Triple(m_arr[i].b, m_arr[i].c, m_arr[i].a));
			m_arr.push_back(Triple(m_arr[i].c, m_arr[i].a, m_arr[i].b));
		}
		n *= 3;
		if (x.a == x.b || x.b == x.c)
			return;
		// With three distinct indices we can in addition interchange the
		// first two indices (to yield all 6 permutations of 3 indices).
		for (size_t i = 0; i < n; ++i) {
			m_arr.push_back(Triple(m_arr[i].b, m_arr[i].a, m_arr[i].c));
		}
		n *= 2;
	}
	size_t Number() const {
		return m_arr.size();
	}
	Triple Member(size_t i) const {
		return m_arr[i];
	}
private:
	std::vector<Triple> m_arr;
};


class QuaternionGrid : public utility::pointer::ReferenceCount {
public:

 	QuaternionGrid(std::string const & name, std::istream & in);
 	void print() const;
 	long num_samples() const { return ntot; }
 	Quaternion quaternion(long i) const { return s.Orientation(i-1); }
	numeric::Real maxrad() const {return maxrad_;}
 	numeric::Real weight(long i) const { return s.Weight(i-1); }
private:
 	numeric::xyzVector<numeric::Real> euler(long i) const { return s.Orientation(i-1).euler(); }
 	std::string name_;
	numeric::Real delta, sigma, maxrad_, coverage;
	size_t ncell, ntot, nent;
	QuatSet s;
	friend std::ostream & operator<<(std::ostream & out, QuaternionGrid const & q);
};
std::ostream & operator<<(std::ostream & out, QuaternionGrid const & q);

class QuatDBMetadata {
public:
	QuatDBMetadata(std::string _name, long _NN, numeric::Real _radius, numeric::Real _cover, numeric::Real _delta ) : name(_name),N(_NN),radius(_radius),cover(_cover),delta(_delta) {}
	std::string name;
	long N;
	numeric::Real radius,cover,delta;
};

class QuaternionGridManager : public utility::SingletonBase< QuaternionGridManager >
{
public:
	friend class utility::SingletonBase< QuaternionGridManager >;

	QuaternionGridCOP request_by_name(std::string const & name);
	QuaternionGridCOP request_by_size(long target_size);
	QuaternionGridCOP request_by_radius(numeric::Real target_radius);

private:
	QuaternionGridManager();
	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static QuaternionGridManager * create_singleton_instance();
	void fill_metadata();

private:

	utility::vector1<QuatDBMetadata> by_size_,by_radius_,by_cover_;
	std::map<std::string,QuaternionGridCOP> grids_;
};

}
}
}


#endif /*INCLUDED_basic_sampling_QuaternionGrid_HH*/
