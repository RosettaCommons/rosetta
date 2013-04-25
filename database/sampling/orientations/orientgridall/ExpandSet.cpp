// Expand an orientation set expressed as a set of grids into a set of
// explicit orientations.
//
// Written by Charles Karney
// Copyright (c) 2006 Sarnoff Corporation. All rights reserved.
//
// For more information, see
//
//    http://charles.karney.info/orientation/
//
// Compile with, e.g.,
//
//    g++ -O2 -o ExpandSet ExpandSet.cpp
//
// Run with
//
//   ./ExpandSet [-e] < grid-file > orientation-file
//
// If -e is specified, the orientations are written as Euler angles,
// otherwise they are written as quaternions.  Format of the grid file:
//
//     Any number of initial comment lines beginning with #
//     A line containing "format grid"
//     A line containing: delta sigma ntot ncell nent maxrad coverage
//     nent lines containing: k l m weight radius mult
//
// Here k >= l >= m >= 0.  delta and sigma are used to define the grid.
// ntot is the total number of orientations, ncell = ntot/24 is the
// number of orientations per cell of the 48-cell.  maxrad is the
// covering radius of the set and radius is the radius of the Voronoi
// cell.  Both are measured in degrees.  coverage is the coverage of the
// set, i.e., how much overlap there is when caps of radius maxrad are
// placed at each point; coverage = 1 means no overlap.
//
// For each triplet, [k l m], generate mult distinct permutations by
// changing the order and the signs of the elements.  Each [k l m] is
// converted to a point in a truncated cube [x y z] =
// [pind(k/2,delta,sigma) pind(l/2,delta,sigma) pind(m/2,delta,sigma)]
// Each [x y z] is converted to a unit quaternion via p = [1 x y z]; q =
// p/|p| to give ncell orientations.  Finally, the 24 rotational cube
// symmetries are applied to the results to yield ntot = orientations.
// The weights are normalized such that sum mult weight = sum mult =
// ncell

// Format of the orientation file
//
//     Any number of initial comment lines beginning with #
//     A line containing "format quaternion" or "format euler"
//     A line containing: ntot maxrad coverage
//     ntot lines containing: q0 q1 q2 q3 weight # for quaternions
//     ntot lines containing: alpha beta gamma weight # for euler.
//
// The weights are normalized such that sum weight = ntot.

#include <cassert>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

namespace {
  char rcsid[] = "$Id: ExpandSet.cpp 6126 2006-04-04 12:32:09Z ckarney $";
}

// Windows doesn't define M_PI in the standard header?
#if !defined(M_PI)
#define M_PI 3.1415926535897932384626433832795028841971694
#endif

using namespace std;

// Minimal quaternion class
class Quaternion {
public:
  double w, x, y, z;
  Quaternion(double ww = 1, double xx = 0, double yy = 0, double zz = 0)
    : w(ww)
    , x(xx)
    , y(yy)
    , z(zz) {}
  void Normalize() {
    double t = w*w + x*x + y*y + z*z;
    assert(t > 0);
    t = 1/sqrt(t);
    w *= t;
    x *= t;
    y *= t;
    z *= t;
    return;
  }
  void Canonicalize() {
    Normalize();
    // Make first biggest element positive
    double mag = w;
    if (abs(x) > abs(mag))
      mag = x;
    if (abs(y) > abs(mag))
      mag = y;
    if (abs(z) > abs(mag))
      mag = z;
    if (mag < 0) {
      w *= -1;
      x *= -1;
      y *= -1;
      z *= -1;
    }
    return;
  }
  // a.Times(b) returns a * b
  Quaternion Times(const Quaternion& q) const {
    double
      mw = w*q.w - x*q.x - y*q.y - z*q.z,
      mx = w*q.x + x*q.w + y*q.z - z*q.y,
      my = w*q.y + y*q.w + z*q.x - x*q.z,
      mz = w*q.z + z*q.w + x*q.y - y*q.x;
    return Quaternion(mw, mx, my, mz);
  }
  void Print(ostream& s) const;
  void PrintEuler(ostream& s) const;
};

// Class to hold a set of orientations and weights
class PackSet {
public:
  Quaternion Orientation(size_t i) const {
    return m_v[i];
  }
  double Weight(size_t i) const {
    return m_w[i];
  }
  size_t Number() const {
    return m_v.size();
  }
  void Add(const Quaternion& q, double w = 1) {
    Quaternion v(q);
    v.Canonicalize();
    m_v.push_back(v);
    m_w.push_back(w);
  }
  void Clear() {
    m_v.clear();
    m_w.clear();
  }
  void Print(ostream& s, bool euler = false, size_t prec = 6) const {
    for (size_t i = 0; i < Number(); ++i) {
      if (euler)
	m_v[i].PrintEuler(s);
      else
	m_v[i].Print(s);
      s << " " << fixed << setprecision(prec) << setw(prec + 2) << m_w[i] << endl;
    }
  }
private:
  vector<Quaternion> m_v;
  vector<double> m_w;
};

// The triple of grid indices
class Triple {
public:
  int a, b, c;
  Triple(int aa, int bb, int cc)
    : a(aa)
    , b(bb)
    , c(cc) {}
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
  vector<Triple> m_arr;
};

// The rotational symmetries of the cube.  (Not normalized, since
// PackSet.Add does this.)
static double CubeSyms[24][4] = {
  {1, 0, 0, 0},
  // 180 deg rotations about 3 axes
  {0, 1, 0, 0},
  {0, 0, 1, 0},
  {0, 0, 0, 1},
  // +/- 120 degree rotations about 4 leading diagonals
  {1, 1, 1, 1},
  {1, 1, 1,-1},
  {1, 1,-1, 1},
  {1, 1,-1,-1},
  {1,-1, 1, 1},
  {1,-1, 1,-1},
  {1,-1,-1, 1},
  {1,-1,-1,-1},
  // +/- 90 degree rotations about 3 axes
  {1, 1, 0, 0},
  {1,-1, 0, 0},
  {1, 0, 1, 0},
  {1, 0,-1, 0},
  {1, 0, 0, 1},
  {1, 0, 0,-1},
  // 180 degree rotations about 6 face diagonals
  {0, 1, 1, 0},
  {0, 1,-1, 0},
  {0, 1, 0, 1},
  {0, 1, 0,-1},
  {0, 0, 1, 1},
  {0, 0, 1,-1},
};

// Convert from index to position.  The sinh scaling tries to compensate
// for the bunching up that occurs when [1 x y z] is projected onto the
// unit sphere.
double pind(double ind, double delta, double sigma) {
  return (sigma == 0) ? ind * delta : sinh(sigma * ind * delta) / sigma;
}

int main(int argc, char* argv[], char*[]) {
  bool euler = false;
  if (argc > 1 && string(argv[1]) == "-e")
    euler = true;
  assert(cin.good());
  string line;
  while (cin.peek() == '#') {
    getline(cin, line);
    cout << line << endl;
  }
  assert(cin.good());
  getline(cin, line);
  assert(line == "format grid");
  cout << "format " << (euler ? "euler" : "quaternion") << endl;
  double delta, sigma, maxrad, coverage;
  size_t ncell, ntot, nent;
  cin >> delta >> sigma >> ntot >> ncell >> nent >> maxrad >> coverage;
  // Use extra digit of precision with weights and radii.  This also
  // triggers a memory minimizing expansion.
  const bool fine = delta < 0.05;
  cout << ntot << " " << fixed
       << setprecision(fine ? 3 : 2) << maxrad << " "
       << setprecision(5) << coverage << endl;
  PackSet s;
  size_t ncell1 = 0;
  for (size_t n = 0; n < nent; ++n) {
    int k, l, m;
    size_t mult;
    double r, w;
    assert(cin.good());
    cin >> k >> l >> m >> w >> r >> mult;
    Permute p(Triple(k, l, m));
    assert(mult == p.Number());
    for (size_t i = 0; i < mult; ++i) {
      Triple t = p.Member(i);
      s.Add(Quaternion(1.0,
		       pind(0.5 * t.a, delta, sigma),
		       pind(0.5 * t.b, delta, sigma),
		       pind(0.5 * t.c, delta, sigma)),
	    w);
    }
    ncell1 += mult;
    if (fine) {
      // Skip n = 0; that's already included.
      for (size_t n = 1; n < 24; ++n) {
	Quaternion q(CubeSyms[n][0], CubeSyms[n][1],
		     CubeSyms[n][2], CubeSyms[n][3]);
	for (size_t i = 0; i < mult; ++i)
	  s.Add(q.Times(s.Orientation(i)), s.Weight(i));
      }
      s.Print(cout, euler, fine ? 7 : 6);
      s.Clear();
    }
  }
  assert(cin.good());
  assert(ncell1 == ncell);
  if (!fine) {
    size_t nc = s.Number();
    assert(nc == ncell);
    for (size_t n = 1; n < 24; ++n) {
      Quaternion q(CubeSyms[n][0], CubeSyms[n][1],
		   CubeSyms[n][2], CubeSyms[n][3]);
      for (size_t i = 0; i < nc; ++i)
	s.Add(q.Times(s.Orientation(i)), s.Weight(i));
    }
    assert(s.Number() == ntot);
    s.Print(cout, euler, fine ? 7 : 6);
    s.Clear();
  }
  return 0;
}

void Quaternion::Print(ostream& s) const {
  s << fixed << setprecision(9) << setw(12) << w << " ";
  s << setw(12) << x << " ";
  s << setw(12) << y << " ";
  s << setw(12) << z;
}

void Quaternion::PrintEuler(ostream& s) const {
  // Print out orientation as a set of Euler angles, following the
  // convention given in
  //
  //    http://www.mhl.soton.ac.uk/research/help/Euler/index.html
  //
  // Rotation by Euler angles [a,b,c] is defined as rotation by c about
  // z axis, followed by rotation by b about y axis. followed by
  // rotation by a about z axis (again).
  //
  // Convert to rotation matrix (assume quaternion is already
  // normalized)
  double
    // m00 = 1 - 2*y*y - 2*z*z,
    m01 =     2*x*y - 2*z*w,
    m02 =     2*x*z + 2*y*w,
    // m10 =     2*x*y + 2*z*w,
    m11 = 1 - 2*x*x - 2*z*z,
    m12 =     2*y*z - 2*x*w,
    m20 =     2*x*z - 2*y*w,
    m21 =     2*y*z + 2*x*w,
    m22 = 1 - 2*x*x - 2*y*y;
  // Taken from Ken Shoemake, "Euler Angle Conversion", Graphics Gems
  // IV, Academic 1994.
  //
  //    http://vered.rose.utoronto.ca/people/david_dir/GEMS/GEMS.html
  double sy = sqrt(m02*m02 + m12*m12);
  //  double sy = sqrt(m10*m10 + m20*m20);
  double a,  b, c;
  b = atan2(sy, m22);
  if (sy > 16 * numeric_limits<double>::epsilon()) {
	a = atan2(m12, m02);
	c = atan2(m21, -m20);
  } else {
    a = atan2(-m01, m11);
    c = 0;
  }
  s << fixed << setprecision(9) << setw(12) << a << " "
    << setw(12) << b << " " << setw(12) << c;

#if !defined(NDEBUG)
  // Sanity check.  Convert from Euler angles back to a quaternion, q
  Quaternion q = Quaternion(cos(a/2), 0, 0, sin(a/2)). // a about z
    Times(Quaternion(cos(b/2), 0, sin(b/2), 0). // b about y
	  Times(Quaternion(cos(c/2), 0, 0, sin(c/2)))); // c about z
  // and check that q is parallel to *this.
  double t = abs(q.w * w + q.x * x + q.y * y + q.z * z);
  assert(t > 1 - 16 * numeric_limits<double>::epsilon());
#endif
}
