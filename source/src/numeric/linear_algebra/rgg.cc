#include <numeric/linear_algebra/rgg.hh>

namespace numeric {
namespace linear_algebra {

using namespace fem::major_types;
using fem::common;

// double epslon( {{{1

/// @brief Estimate unit roundoff in quantities of size x.
///
/// This program should function properly on all systems satisfying the 
/// following two assumptions:
///
/// 1.  The base used in representing floating point numbers is not a power of 
///     three.
/// 2.  The quantity a in statement 10 is represented to the accuracy used in 
///     floating point variables that are stored in memory.
///
/// The statement number 10 and the goto 10 are intended to force optimizing 
/// compilers to generate code satisfying assumption 2.  Under these 
/// assumptions, it should be true that:
///
/// -   `a` is not exactly equal to four-thirds.
/// -   `b` has a zero for its last bit or digit.
/// -   `c` is not exactly equal to one.
/// -   `eps` measures the separation of 1.0 from the next larger floating
///     point number.
///
/// The developers of EISPACK would appreciate being informed about any systems 
/// where these assumptions do not hold.  This version dated 4/6/83.  Fable was 
/// used to convert the original source to C++.

double epslon(
		double const& x) {

	double return_value = fem::double0;
	//double a = fem::double0;
	double b;
	double c;
	double eps;

	double a = 4.0e0 / 3.0e0;
	statement_10:
	b = a - 1.0e0;
	c = b + b + b;
	eps = fem::dabs(c - 1.0e0);
	if (eps == 0.0e0) {
		goto statement_10;
	}
	return_value = eps * fem::dabs(x);
	return return_value;
}

// void qzhes(  {{{1

/// @brief This subroutine is the first step of the QZ algorithm for solving 
/// generalized matrix eigenvalue problems.
///
/// This subroutine accepts a pair of real general matrices and reduces one of 
/// them to upper Hessenberg form and the other to upper triangular form using 
/// orthogonal transformations.  it is usually followed by `qzit`, `qzval` and, 
/// possibly, `qzvec`.
///
/// @param[in] nm Must be set to the row dimension of two-dimensional array 
/// parameters as declared in the calling program dimension statement.
///
/// @param[in] n The order of the matrices.
///
/// @param[in] a Contains a real general matrix.
///
/// @param[in] b Contains a real general matrix.
///
/// @param[in] matz Should be set to `true` if the right hand transformations
///      are to be accumulated for later use in computing eigenvectors, and to 
///      `false` otherwise.
///
/// @param[out] a Reduced to upper Hessenberg form.  The elements below the 
/// first subdiagonal have been set to zero.
///
/// @param[out] b Reduced to upper triangular form.  The elements below the 
/// main diagonal have been set to zero.
///
/// @param[out] z Contains the product of the right hand transformations if 
/// matz has been set to `true`.  Otherwise, `z` is not referenced.
///
/// Reference: Siam J. Numer. Anal. 10, 241-256(1973) by Moler and Stewart.
///
/// Questions and comments should be directed to Burton S. Garbow,
/// Mathematics and Computer Science Division, Argonne National Laboratory.
/// This version dated August 1983.  Fable was used to convert the original 
/// Fortran source into C++.

void qzhes(
		int const& nm,
		int const& n,
		arr_ref<double, 2> a,
		arr_ref<double, 2> b,
		bool const& matz,
		arr_ref<double, 2> z) {

  a(dimension(nm, n));
  b(dimension(nm, n));
  z(dimension(nm, n));
  int j = fem::int0;
  int i = fem::int0;
  int nm1 = fem::int0;
  int l = fem::int0;
  int l1 = fem::int0;
  double s = fem::double0;
  double r = fem::double0;
  double rho = fem::double0;
  double t = fem::double0;
  int nm2 = fem::int0;
  int k = fem::int0;
  int nk1 = fem::int0;
  int lb = fem::int0;
  double u1 = fem::double0;
  double u2 = fem::double0;
  double v1 = fem::double0;
  double v2 = fem::double0;

	//C     .......... initialize z ..........
  if (!matz) {
    goto statement_10;
  }

  FEM_DO_SAFE(j, 1, n) {
    FEM_DO_SAFE(i, 1, n) {
      z(i, j) = 0.0e0;
    }
    z(j, j) = 1.0e0;
  }
  //C     .......... reduce b to upper triangular form ..........
  statement_10:
  if (n <= 1) {
    goto statement_170;
  }
  nm1 = n - 1;

  FEM_DO_SAFE(l, 1, nm1) {
    l1 = l + 1;
    s = 0.0e0;

    FEM_DO_SAFE(i, l1, n) {
      s += fem::dabs(b(i, l));
    }

    if (s == 0.0e0) {
      goto statement_100;
    }
    s += fem::dabs(b(l, l));
    r = 0.0e0;

    FEM_DO_SAFE(i, l, n) {
      b(i, l) = b(i, l) / s;
      r += fem::pow2(b(i, l));
    }

    r = fem::dsign(fem::dsqrt(r), b(l, l));
    b(l, l) += r;
    rho = r * b(l, l);

    FEM_DO_SAFE(j, l1, n) {
      t = 0.0e0;

      FEM_DO_SAFE(i, l, n) {
        t += b(i, l) * b(i, j);
      }

      t = -t / rho;

      FEM_DO_SAFE(i, l, n) {
        b(i, j) += t * b(i, l);
      }

    }

    FEM_DO_SAFE(j, 1, n) {
      t = 0.0e0;

      FEM_DO_SAFE(i, l, n) {
        t += b(i, l) * a(i, j);
      }

      t = -t / rho;

      FEM_DO_SAFE(i, l, n) {
        a(i, j) += t * b(i, l);
      }

    }

    b(l, l) = -s * r;

    FEM_DO_SAFE(i, l1, n) {
      b(i, l) = 0.0e0;
    }

    statement_100:;
  }
  //C     .......... reduce a to upper hessenberg form, while
  //C                keeping b triangular ..........
  if (n == 2) {
    goto statement_170;
  }
  nm2 = n - 2;

  FEM_DO_SAFE(k, 1, nm2) {
    nk1 = nm1 - k;
    //C     .......... for l=n-1 step -1 until k+1 do -- ..........
    FEM_DO_SAFE(lb, 1, nk1) {
      l = n - lb;
      l1 = l + 1;
      //C     .......... zero a(l+1,k) ..........
      s = fem::dabs(a(l, k)) + fem::dabs(a(l1, k));
      if (s == 0.0e0) {
        goto statement_150;
      }
      u1 = a(l, k) / s;
      u2 = a(l1, k) / s;
      r = fem::dsign(fem::dsqrt(u1 * u1 + u2 * u2), u1);
      v1 = -(u1 + r) / r;
      v2 = -u2 / r;
      u2 = v2 / v1;

      FEM_DO_SAFE(j, k, n) {
        t = a(l, j) + u2 * a(l1, j);
        a(l, j) += t * v1;
        a(l1, j) += t * v2;
      }

      a(l1, k) = 0.0e0;

      FEM_DO_SAFE(j, l, n) {
        t = b(l, j) + u2 * b(l1, j);
        b(l, j) += t * v1;
        b(l1, j) += t * v2;
      }
      //C     .......... zero b(l+1,l) ..........
      s = fem::dabs(b(l1, l1)) + fem::dabs(b(l1, l));
      if (s == 0.0e0) {
        goto statement_150;
      }
      u1 = b(l1, l1) / s;
      u2 = b(l1, l) / s;
      r = fem::dsign(fem::dsqrt(u1 * u1 + u2 * u2), u1);
      v1 = -(u1 + r) / r;
      v2 = -u2 / r;
      u2 = v2 / v1;

      FEM_DO_SAFE(i, 1, l1) {
        t = b(i, l1) + u2 * b(i, l);
        b(i, l1) += t * v1;
        b(i, l) += t * v2;
      }

      b(l1, l) = 0.0e0;

      FEM_DO_SAFE(i, 1, n) {
        t = a(i, l1) + u2 * a(i, l);
        a(i, l1) += t * v1;
        a(i, l) += t * v2;
      }

      if (!matz) {
        goto statement_150;
      }

      FEM_DO_SAFE(i, 1, n) {
        t = z(i, l1) + u2 * z(i, l);
        z(i, l1) += t * v1;
        z(i, l) += t * v2;
      }

      statement_150:;
    }

  }

  statement_170:;
}

// void qzit(  {{{1

/// @brief This subroutine is the second step of the QZ algorithm
/// for solving generalized matrix eigenvalue problems.
///
/// This subroutine accepts a pair of real matrices, one of them
/// in upper Hessenberg form and the other in upper triangular form.
/// it reduces the Hessenberg matrix to quasi-triangular form using
/// orthogonal transformations while maintaining the triangular form
/// of the other matrix.  It is usually preceded by `qzhes` and
/// followed by `qzval` and, possibly, `qzvec`.
///
/// @param[in] nm Must be set to the row dimension of two-dimensional array 
/// parameters as declared in the calling program dimension statement.
///
/// @param[in] n The order of the matrices.
///
/// @param[in] a Contains a real upper Hessenberg matrix.
///
/// @param[in] b Contains a real upper triangular matrix.
///
/// @param[in] eps1 A tolerance used to determine negligible elements.  Zero or 
/// negative values may be given, in which case an element will be neglected 
/// only if it is less than roundoff error times the norm of its matrix.  If 
/// the input `eps1` is positive, then an element will be considered negligible 
/// if it is less than `eps1` times the norm of its matrix.  A positive value 
/// of `eps1` may result in faster execution, but less accurate results.
///
/// @param[in] matz Should be set to `true` if the right hand transformations 
/// are to be accumulated for later use in computing eigenvectors, and to 
/// `false` otherwise.
///
/// @param[in] z Contains, if matz has been set to `true`, the transformation 
/// matrix produced in the reduction by `qzhes`, if performed, or else the 
/// identity matrix.  if `matz` has been set to `false`, `z` is not referenced.
///
/// @param[out] a Reduced to quasi-triangular form.  The elements below the 
/// first subdiagonal are still zero and no two consecutive subdiagonal 
/// elements are nonzero.
///
/// @param[out] b Still in upper triangular form, although its elements have 
/// been altered.  The location `b(n,1)` is used to store `eps1` times the norm 
/// of `b` for later use by `qzval` and `qzvec`.
///
/// @param[out] z Contains the product of the right hand transformations (for 
/// both steps) if `matz` has been set to `true`.
///
/// @param[out] ierr Set to zero for normal return, `j` if the limit of 30*n 
/// iterations is exhausted while the j-th eigenvalue is being sought.
///
/// Reference: Siam J. Numer. Anal. 10, 241-256(1973) by Moler and Stewart,
/// as modified in technical note NASA tn D-7305(1973) by Ward.
///
/// Questions and comments should be directed to Burton S. Garbow,
/// Mathematics and Computer Science Division, Argonne National Laboratory.
/// This version dated August 1983.  Fable was used to convert the original 
/// Fortran source into C++.


void qzit(
		int const& nm,
		int const& n,
		arr_ref<double, 2> a,
		arr_ref<double, 2> b,
		double const& eps1,
		bool const& matz,
		arr_ref<double, 2> z,
		int& ierr) {

  a(dimension(nm, n));
  b(dimension(nm, n));
  z(dimension(nm, n));
  //double anorm = fem::double0;
  //double bnorm = fem::double0;
  int i = fem::int0;
  double ani = fem::double0;
  double bni = fem::double0;
  int j = fem::int0;
  double ep = fem::double0;
  double epsa = fem::double0;
  double epsb = fem::double0;
  int lor1 = fem::int0;
  int enorn = fem::int0;
  int en = fem::int0;
  int itn = fem::int0;
  int its = fem::int0;
  int na = fem::int0;
  int enm2 = fem::int0;
  int ish = fem::int0;
  int ll = fem::int0;
  int lm1 = fem::int0;
  int l = fem::int0;
  int ld = fem::int0;
  int l1 = fem::int0;
  double b11 = fem::double0;
  double s = fem::double0;
  double u1 = fem::double0;
  double u2 = fem::double0;
  double r = fem::double0;
  double v1 = fem::double0;
  double v2 = fem::double0;
  double t = fem::double0;
  double a11 = fem::double0;
  double a21 = fem::double0;
  double b22 = fem::double0;
  double b33 = fem::double0;
  double b44 = fem::double0;
  double a33 = fem::double0;
  double a34 = fem::double0;
  double a43 = fem::double0;
  double a44 = fem::double0;
  double b34 = fem::double0;
  double sh = fem::double0;
  double a1 = fem::double0;
  double a2 = fem::double0;
  double a12 = fem::double0;
  double a22 = fem::double0;
  double b12 = fem::double0;
  double a3 = fem::double0;
  int k = fem::int0;
  bool notlas = fem::bool0;
  int k1 = fem::int0;
  int k2 = fem::int0;
  int km1 = fem::int0;
  double u3 = fem::double0;
  double v3 = fem::double0;

  ierr = 0;
  //C     .......... compute epsa,epsb ..........
  double anorm = 0.0e0;
  double bnorm = 0.0e0;

  FEM_DO_SAFE(i, 1, n) {
    ani = 0.0e0;
    if (i != 1) {
      ani = fem::dabs(a(i, i - 1));
    }
    bni = 0.0e0;

    FEM_DO_SAFE(j, i, n) {
      ani += fem::dabs(a(i, j));
      bni += fem::dabs(b(i, j));
    }

    if (ani > anorm) {
      anorm = ani;
    }
    if (bni > bnorm) {
      bnorm = bni;
    }
  }

  if (anorm == 0.0e0) {
    anorm = 1.0e0;
  }
  if (bnorm == 0.0e0) {
    bnorm = 1.0e0;
  }
  ep = eps1;
  if (ep > 0.0e0) {
    goto statement_50;
  }
  //C     .......... use roundoff level if eps1 is zero ..........
  ep = epslon(1.0e0);
  statement_50:
  epsa = ep * anorm;
  epsb = ep * bnorm;
  //C     .......... reduce a to quasi-triangular form, while
  //C                keeping b triangular ..........
  lor1 = 1;
  enorn = n;
  en = n;
  itn = 30 * n;
  //C     .......... begin qz step ..........
  statement_60:
  if (en <= 2) {
    goto statement_1001;
  }
  if (!matz) {
    enorn = en;
  }
  its = 0;
  na = en - 1;
  enm2 = na - 1;
  statement_70:
  ish = 2;
  //C     .......... check for convergence or reducibility.
  //C                for l=en step -1 until 1 do -- ..........
  FEM_DO_SAFE(ll, 1, en) {
    lm1 = en - ll;
    l = lm1 + 1;
    if (l == 1) {
      goto statement_95;
    }
    if (fem::dabs(a(l, lm1)) <= epsa) {
      goto statement_90;
    }
  }

  statement_90:
  a(l, lm1) = 0.0e0;
  if (l < na) {
    goto statement_95;
  }
  //C     .......... 1-by-1 or 2-by-2 block isolated ..........
  en = lm1;
  goto statement_60;
  //C     .......... check for small top of b ..........
  statement_95:
  ld = l;
  statement_100:
  l1 = l + 1;
  b11 = b(l, l);
  if (fem::dabs(b11) > epsb) {
    goto statement_120;
  }
  b(l, l) = 0.0e0;
  s = fem::dabs(a(l, l)) + fem::dabs(a(l1, l));
  u1 = a(l, l) / s;
  u2 = a(l1, l) / s;
  r = fem::dsign(fem::dsqrt(u1 * u1 + u2 * u2), u1);
  v1 = -(u1 + r) / r;
  v2 = -u2 / r;
  u2 = v2 / v1;

  FEM_DO_SAFE(j, l, enorn) {
    t = a(l, j) + u2 * a(l1, j);
    a(l, j) += t * v1;
    a(l1, j) += t * v2;
    t = b(l, j) + u2 * b(l1, j);
    b(l, j) += t * v1;
    b(l1, j) += t * v2;
  }

  if (l != 1) {
    a(l, lm1) = -a(l, lm1);
  }
  lm1 = l;
  l = l1;
  goto statement_90;
  statement_120:
  a11 = a(l, l) / b11;
  a21 = a(l1, l) / b11;
  if (ish == 1) {
    goto statement_140;
  }
  //C     .......... iteration strategy ..........
  if (itn == 0) {
    goto statement_1000;
  }
  if (its == 10) {
    goto statement_155;
  }
  //C     .......... determine type of shift ..........
  b22 = b(l1, l1);
  if (fem::dabs(b22) < epsb) {
    b22 = epsb;
  }
  b33 = b(na, na);
  if (fem::dabs(b33) < epsb) {
    b33 = epsb;
  }
  b44 = b(en, en);
  if (fem::dabs(b44) < epsb) {
    b44 = epsb;
  }
  a33 = a(na, na) / b33;
  a34 = a(na, en) / b44;
  a43 = a(en, na) / b33;
  a44 = a(en, en) / b44;
  b34 = b(na, en) / b44;
  t = 0.5e0 * (a43 * b34 - a33 - a44);
  r = t * t + a34 * a43 - a33 * a44;
  if (r < 0.0e0) {
    goto statement_150;
  }
  //C     .......... determine single shift zeroth column of a ..........
  ish = 1;
  r = fem::dsqrt(r);
  sh = -t + r;
  s = -t - r;
  if (fem::dabs(s - a44) < fem::dabs(sh - a44)) {
    sh = s;
  }
  //C     .......... look for two consecutive small
  //C                sub-diagonal elements of a.
  //C                for l=en-2 step -1 until ld do -- ..........
  FEM_DO_SAFE(ll, ld, enm2) {
    l = enm2 + ld - ll;
    if (l == ld) {
      goto statement_140;
    }
    lm1 = l - 1;
    l1 = l + 1;
    t = a(l, l);
    if (fem::dabs(b(l, l)) > epsb) {
      t = t - sh * b(l, l);
    }
    if (fem::dabs(a(l, lm1)) <= fem::dabs(t / a(l1, l)) * epsa) {
      goto statement_100;
    }
  }

  statement_140:
  a1 = a11 - sh;
  a2 = a21;
  if (l != ld) {
    a(l, lm1) = -a(l, lm1);
  }
  goto statement_160;
  //C     .......... determine double shift zeroth column of a ..........
  statement_150:
  a12 = a(l, l1) / b22;
  a22 = a(l1, l1) / b22;
  b12 = b(l, l1) / b22;
  a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11) /
    a21 + a12 - a11 * b12;
  a2 = (a22 - a11) - a21 * b12 - (a33 - a11) - (a44 - a11) + a43 * b34;
  a3 = a(l1 + 1, l1) / b22;
  goto statement_160;
  //C     .......... ad hoc shift ..........
  statement_155:
  a1 = 0.0e0;
  a2 = 1.0e0;
  a3 = 1.1605e0;
  statement_160:
  its++;
  itn = itn - 1;
  if (!matz) {
    lor1 = ld;
  }
  //C     .......... main loop ..........
  FEM_DO_SAFE(k, l, na) {
    notlas = k != na && ish == 2;
    k1 = k + 1;
    k2 = k + 2;
    km1 = fem::max0(k - 1, l);
    ll = fem::min0(en, k1 + ish);
    if (notlas) {
      goto statement_190;
    }
    //C     .......... zero a(k+1,k-1) ..........
    if (k == l) {
      goto statement_170;
    }
    a1 = a(k, km1);
    a2 = a(k1, km1);
    statement_170:
    s = fem::dabs(a1) + fem::dabs(a2);
    if (s == 0.0e0) {
      goto statement_70;
    }
    u1 = a1 / s;
    u2 = a2 / s;
    r = fem::dsign(fem::dsqrt(u1 * u1 + u2 * u2), u1);
    v1 = -(u1 + r) / r;
    v2 = -u2 / r;
    u2 = v2 / v1;

    FEM_DO_SAFE(j, km1, enorn) {
      t = a(k, j) + u2 * a(k1, j);
      a(k, j) += t * v1;
      a(k1, j) += t * v2;
      t = b(k, j) + u2 * b(k1, j);
      b(k, j) += t * v1;
      b(k1, j) += t * v2;
    }

    if (k != l) {
      a(k1, km1) = 0.0e0;
    }
    goto statement_240;
    //C     .......... zero a(k+1,k-1) and a(k+2,k-1) ..........
    statement_190:
    if (k == l) {
      goto statement_200;
    }
    a1 = a(k, km1);
    a2 = a(k1, km1);
    a3 = a(k2, km1);
    statement_200:
    s = fem::dabs(a1) + fem::dabs(a2) + fem::dabs(a3);
    if (s == 0.0e0) {
      goto statement_260;
    }
    u1 = a1 / s;
    u2 = a2 / s;
    u3 = a3 / s;
    r = fem::dsign(fem::dsqrt(u1 * u1 + u2 * u2 + u3 * u3), u1);
    v1 = -(u1 + r) / r;
    v2 = -u2 / r;
    v3 = -u3 / r;
    u2 = v2 / v1;
    u3 = v3 / v1;

    FEM_DO_SAFE(j, km1, enorn) {
      t = a(k, j) + u2 * a(k1, j) + u3 * a(k2, j);
      a(k, j) += t * v1;
      a(k1, j) += t * v2;
      a(k2, j) += t * v3;
      t = b(k, j) + u2 * b(k1, j) + u3 * b(k2, j);
      b(k, j) += t * v1;
      b(k1, j) += t * v2;
      b(k2, j) += t * v3;
    }

    if (k == l) {
      goto statement_220;
    }
    a(k1, km1) = 0.0e0;
    a(k2, km1) = 0.0e0;
    //C     .......... zero b(k+2,k+1) and b(k+2,k) ..........
    statement_220:
    s = fem::dabs(b(k2, k2)) + fem::dabs(b(k2, k1)) + fem::dabs(b(k2, k));
    if (s == 0.0e0) {
      goto statement_240;
    }
    u1 = b(k2, k2) / s;
    u2 = b(k2, k1) / s;
    u3 = b(k2, k) / s;
    r = fem::dsign(fem::dsqrt(u1 * u1 + u2 * u2 + u3 * u3), u1);
    v1 = -(u1 + r) / r;
    v2 = -u2 / r;
    v3 = -u3 / r;
    u2 = v2 / v1;
    u3 = v3 / v1;

    FEM_DO_SAFE(i, lor1, ll) {
      t = a(i, k2) + u2 * a(i, k1) + u3 * a(i, k);
      a(i, k2) += t * v1;
      a(i, k1) += t * v2;
      a(i, k) += t * v3;
      t = b(i, k2) + u2 * b(i, k1) + u3 * b(i, k);
      b(i, k2) += t * v1;
      b(i, k1) += t * v2;
      b(i, k) += t * v3;
    }

    b(k2, k) = 0.0e0;
    b(k2, k1) = 0.0e0;
    if (!matz) {
      goto statement_240;
    }

    FEM_DO_SAFE(i, 1, n) {
      t = z(i, k2) + u2 * z(i, k1) + u3 * z(i, k);
      z(i, k2) += t * v1;
      z(i, k1) += t * v2;
      z(i, k) += t * v3;
    }
    //C     .......... zero b(k+1,k) ..........
    statement_240:
    s = fem::dabs(b(k1, k1)) + fem::dabs(b(k1, k));
    if (s == 0.0e0) {
      goto statement_260;
    }
    u1 = b(k1, k1) / s;
    u2 = b(k1, k) / s;
    r = fem::dsign(fem::dsqrt(u1 * u1 + u2 * u2), u1);
    v1 = -(u1 + r) / r;
    v2 = -u2 / r;
    u2 = v2 / v1;

    FEM_DO_SAFE(i, lor1, ll) {
      t = a(i, k1) + u2 * a(i, k);
      a(i, k1) += t * v1;
      a(i, k) += t * v2;
      t = b(i, k1) + u2 * b(i, k);
      b(i, k1) += t * v1;
      b(i, k) += t * v2;
    }

    b(k1, k) = 0.0e0;
    if (!matz) {
      goto statement_260;
    }

    FEM_DO_SAFE(i, 1, n) {
      t = z(i, k1) + u2 * z(i, k);
      z(i, k1) += t * v1;
      z(i, k) += t * v2;
    }

    statement_260:;
  }
  //C     .......... end qz step ..........
  goto statement_70;
  //C     .......... set error -- all eigenvalues have not
  //C                converged after 30*n iterations ..........
  statement_1000:
  ierr = en;
  //C     .......... save epsb for use by qzval and qzvec ..........
  statement_1001:
  if (n > 1) {
    b(n, 1) = epsb;
  }
}

// void qzval(  {{{1

/// @brief This subroutine is the third step of the QZ algorithm
/// for solving generalized matrix eigenvalue problems,
///
/// This subroutine accepts a pair of real matrices, one of them
/// in quasi-triangular form and the other in upper triangular form.
/// It reduces the quasi-triangular matrix further, so that any
/// remaining 2-by-2 blocks correspond to pairs of complex
/// eigenvalues, and returns quantities whose ratios give the
/// generalized eigenvalues.  It is usually preceded by `qzhes`
/// and `qzit` and may be followed by `qzvec`.
///
/// @param[in] nm Must be set to the row dimension of two-dimensional array 
/// parameters as declared in the calling program dimension statement.
///
/// @param[in] n The order of the matrices.
///
/// @param[in] a Contains a real upper quasi-triangular matrix.
///
/// @param[in] b Contains a real upper triangular matrix.  In addition, 
/// location `b(n,1)` contains the tolerance quantity `epsb` computed and saved 
/// in `qzit`.
///
/// @param[in] matz Should be set to `true` if the right hand transformations 
/// are to be accumulated for later use in computing eigenvectors, and to 
/// `false` otherwise.
///
/// @param[in] z Contains, if `matz` has been set to `true`, the transformation 
/// matrix produced in the reductions by `qzhes` and `qzit`, if performed, or 
/// else the identity matrix.  If `matz` has been set to `false`, `z` is not 
/// referenced.
///
/// @param[out] a Reduced further to a quasi-triangular matrix in which all 
/// nonzero subdiagonal elements correspond to pairs of complex eigenvalues.
///
/// @param[out] b Still in upper triangular form, although its elements have 
/// been altered.  `b(n,1)` is unaltered.
///
/// @param[out] alfr,alfi Contain the real and imaginary parts of the diagonal 
/// elements of the triangular matrix that would be obtained if a were reduced 
/// completely to triangular form by unitary transformations.  Non-zero values 
/// of `alfi` occur in pairs, the first member positive and the second 
/// negative.
///
/// @param[out] beta Contains the diagonal elements of the corresponding `b`, 
/// normalized to be real and non-negative.  The generalized eigenvalues are 
/// then the ratios \f$\frac{\mathrm{alfr} + i \times \mathrm{alfi}} 
/// {\mathrm{beta}}\f$.
///
/// @param[out] z Contains the product of the right hand transformations (for 
/// all three steps) if `matz` has been set to `true`.
///
/// Reference: Siam J. Numer. Anal. 10, 241-256(1973) by Moler and Stewart.
///
/// Questions and comments should be directed to Burton S. Garbow,
/// Mathematics and Computer Science Division, Argonne National Laboratory.
/// This version dated August 1983.  Fable was used to convert the original 
/// Fortran source into C++.

void qzval(
		int const& nm,
		int const& n,
		arr_ref<double, 2> a,
		arr_ref<double, 2> b,
		arr_ref<double> alfr,
		arr_ref<double> alfi,
		arr_ref<double> beta,
		bool const& matz,
		arr_ref<double, 2> z) {

  a(dimension(nm, n));
  b(dimension(nm, n));
  alfr(dimension(n));
  alfi(dimension(n));
  beta(dimension(n));
  z(dimension(nm, n));
  //double epsb = fem::double0;
  //int isw = fem::int0;
  int nn = fem::int0;
  int en = fem::int0;
  int na = fem::int0;
  double a1 = fem::double0;
  double a2 = fem::double0;
  double bn = fem::double0;
  double an = fem::double0;
  double a11 = fem::double0;
  double a12 = fem::double0;
  double a21 = fem::double0;
  double a22 = fem::double0;
  double b11 = fem::double0;
  double b12 = fem::double0;
  double b22 = fem::double0;
  double e = fem::double0;
  double ei = fem::double0;
  double s = fem::double0;
  double t = fem::double0;
  double c = fem::double0;
  double d = fem::double0;
  double u1 = fem::double0;
  double u2 = fem::double0;
  double r = fem::double0;
  double v1 = fem::double0;
  double v2 = fem::double0;
  int i = fem::int0;
  int j = fem::int0;
  double a11r = fem::double0;
  double a11i = fem::double0;
  double a12r = fem::double0;
  double a12i = fem::double0;
  double a22r = fem::double0;
  double a22i = fem::double0;
  double a1i = fem::double0;
  double a2i = fem::double0;
  double cz = fem::double0;
  double szr = fem::double0;
  double szi = fem::double0;
  double cq = fem::double0;
  double sqr = fem::double0;
  double sqi = fem::double0;
  double ssr = fem::double0;
  double ssi = fem::double0;
  double tr = fem::double0;
  double ti = fem::double0;
  double dr = fem::double0;
  double di = fem::double0;

  double epsb = b(n, 1);
  int isw = 1;
  //C     .......... find eigenvalues of quasi-triangular matrices.
  //C                for en=n step -1 until 1 do -- ..........
  FEM_DO_SAFE(nn, 1, n) {
    en = n + 1 - nn;
    na = en - 1;
    if (isw == 2) {
      goto statement_505;
    }
    if (en == 1) {
      goto statement_410;
    }
    if (a(en, na) != 0.0e0) {
      goto statement_420;
    }
    //C     .......... 1-by-1 block, one real root ..........
    statement_410:
    alfr(en) = a(en, en);
    if (b(en, en) < 0.0e0) {
      alfr(en) = -alfr(en);
    }
    beta(en) = fem::dabs(b(en, en));
    alfi(en) = 0.0e0;
    goto statement_510;
    //C     .......... 2-by-2 block ..........
    statement_420:
    if (fem::dabs(b(na, na)) <= epsb) {
      goto statement_455;
    }
    if (fem::dabs(b(en, en)) > epsb) {
      goto statement_430;
    }
    a1 = a(en, en);
    a2 = a(en, na);
    bn = 0.0e0;
    goto statement_435;
    statement_430:
    an = fem::dabs(a(na, na)) + fem::dabs(a(na, en)) + fem::dabs(a(en,
      na)) + fem::dabs(a(en, en));
    bn = fem::dabs(b(na, na)) + fem::dabs(b(na, en)) + fem::dabs(b(en, en));
    a11 = a(na, na) / an;
    a12 = a(na, en) / an;
    a21 = a(en, na) / an;
    a22 = a(en, en) / an;
    b11 = b(na, na) / bn;
    b12 = b(na, en) / bn;
    b22 = b(en, en) / bn;
    e = a11 / b11;
    ei = a22 / b22;
    s = a21 / (b11 * b22);
    t = (a22 - e * b22) / b22;
    if (fem::dabs(e) <= fem::dabs(ei)) {
      goto statement_431;
    }
    e = ei;
    t = (a11 - e * b11) / b11;
    statement_431:
    c = 0.5e0 * (t - s * b12);
    d = c * c + s * (a12 - e * b12);
    if (d < 0.0e0) {
      goto statement_480;
    }
    //C     .......... two real roots.
    //C                zero both a(en,na) and b(en,na) ..........
    e += (c + fem::dsign(fem::dsqrt(d), c));
    a11 = a11 - e * b11;
    a12 = a12 - e * b12;
    a22 = a22 - e * b22;
    if (fem::dabs(a11) + fem::dabs(a12) < fem::dabs(a21) + fem::dabs(a22)) {
      goto statement_432;
    }
    a1 = a12;
    a2 = a11;
    goto statement_435;
    statement_432:
    a1 = a22;
    a2 = a21;
    //C     .......... choose and apply real z ..........
    statement_435:
    s = fem::dabs(a1) + fem::dabs(a2);
    u1 = a1 / s;
    u2 = a2 / s;
    r = fem::dsign(fem::dsqrt(u1 * u1 + u2 * u2), u1);
    v1 = -(u1 + r) / r;
    v2 = -u2 / r;
    u2 = v2 / v1;

    FEM_DO_SAFE(i, 1, en) {
      t = a(i, en) + u2 * a(i, na);
      a(i, en) += t * v1;
      a(i, na) += t * v2;
      t = b(i, en) + u2 * b(i, na);
      b(i, en) += t * v1;
      b(i, na) += t * v2;
    }

    if (!matz) {
      goto statement_450;
    }

    FEM_DO_SAFE(i, 1, n) {
      t = z(i, en) + u2 * z(i, na);
      z(i, en) += t * v1;
      z(i, na) += t * v2;
    }

    statement_450:
    if (bn == 0.0e0) {
      goto statement_475;
    }
    if (an < fem::dabs(e) * bn) {
      goto statement_455;
    }
    a1 = b(na, na);
    a2 = b(en, na);
    goto statement_460;
    statement_455:
    a1 = a(na, na);
    a2 = a(en, na);
    //C     .......... choose and apply real q ..........
    statement_460:
    s = fem::dabs(a1) + fem::dabs(a2);
    if (s == 0.0e0) {
      goto statement_475;
    }
    u1 = a1 / s;
    u2 = a2 / s;
    r = fem::dsign(fem::dsqrt(u1 * u1 + u2 * u2), u1);
    v1 = -(u1 + r) / r;
    v2 = -u2 / r;
    u2 = v2 / v1;

    FEM_DO_SAFE(j, na, n) {
      t = a(na, j) + u2 * a(en, j);
      a(na, j) += t * v1;
      a(en, j) += t * v2;
      t = b(na, j) + u2 * b(en, j);
      b(na, j) += t * v1;
      b(en, j) += t * v2;
    }

    statement_475:
    a(en, na) = 0.0e0;
    b(en, na) = 0.0e0;
    alfr(na) = a(na, na);
    alfr(en) = a(en, en);
    if (b(na, na) < 0.0e0) {
      alfr(na) = -alfr(na);
    }
    if (b(en, en) < 0.0e0) {
      alfr(en) = -alfr(en);
    }
    beta(na) = fem::dabs(b(na, na));
    beta(en) = fem::dabs(b(en, en));
    alfi(en) = 0.0e0;
    alfi(na) = 0.0e0;
    goto statement_505;
    //C     .......... two complex roots ..........
    statement_480:
    e += c;
    ei = fem::dsqrt(-d);
    a11r = a11 - e * b11;
    a11i = ei * b11;
    a12r = a12 - e * b12;
    a12i = ei * b12;
    a22r = a22 - e * b22;
    a22i = ei * b22;
    if (fem::dabs(a11r) + fem::dabs(a11i) + fem::dabs(a12r) +
        fem::dabs(a12i) < fem::dabs(a21) + fem::dabs(a22r) +
        fem::dabs(a22i)) {
      goto statement_482;
    }
    a1 = a12r;
    a1i = a12i;
    a2 = -a11r;
    a2i = -a11i;
    goto statement_485;
    statement_482:
    a1 = a22r;
    a1i = a22i;
    a2 = -a21;
    a2i = 0.0e0;
    //C     .......... choose complex z ..........
    statement_485:
    cz = fem::dsqrt(a1 * a1 + a1i * a1i);
    if (cz == 0.0e0) {
      goto statement_487;
    }
    szr = (a1 * a2 + a1i * a2i) / cz;
    szi = (a1 * a2i - a1i * a2) / cz;
    r = fem::dsqrt(cz * cz + szr * szr + szi * szi);
    cz = cz / r;
    szr = szr / r;
    szi = szi / r;
    goto statement_490;
    statement_487:
    szr = 1.0e0;
    szi = 0.0e0;
    statement_490:
    if (an < (fem::dabs(e) + ei) * bn) {
      goto statement_492;
    }
    a1 = cz * b11 + szr * b12;
    a1i = szi * b12;
    a2 = szr * b22;
    a2i = szi * b22;
    goto statement_495;
    statement_492:
    a1 = cz * a11 + szr * a12;
    a1i = szi * a12;
    a2 = cz * a21 + szr * a22;
    a2i = szi * a22;
    //C     .......... choose complex q ..........
    statement_495:
    cq = fem::dsqrt(a1 * a1 + a1i * a1i);
    if (cq == 0.0e0) {
      goto statement_497;
    }
    sqr = (a1 * a2 + a1i * a2i) / cq;
    sqi = (a1 * a2i - a1i * a2) / cq;
    r = fem::dsqrt(cq * cq + sqr * sqr + sqi * sqi);
    cq = cq / r;
    sqr = sqr / r;
    sqi = sqi / r;
    goto statement_500;
    statement_497:
    sqr = 1.0e0;
    sqi = 0.0e0;
    //C     .......... compute diagonal elements that would result
    //C                if transformations were applied ..........
    statement_500:
    ssr = sqr * szr + sqi * szi;
    ssi = sqr * szi - sqi * szr;
    i = 1;
    tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22;
    ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22;
    dr = cq * cz * b11 + cq * szr * b12 + ssr * b22;
    di = cq * szi * b12 + ssi * b22;
    goto statement_503;
    statement_502:
    i = 2;
    tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22;
    ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21;
    dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22;
    di = -ssi * b11 - sqi * cz * b12;
    statement_503:
    t = ti * dr - tr * di;
    j = na;
    if (t < 0.0e0) {
      j = en;
    }
    r = fem::dsqrt(dr * dr + di * di);
    beta(j) = bn * r;
    alfr(j) = an * (tr * dr + ti * di) / r;
    alfi(j) = an * t / r;
    if (i == 1) {
      goto statement_502;
    }
    statement_505:
    isw = 3 - isw;
    statement_510:;
  }
  b(n, 1) = epsb;

}

// void qzvec(  {{{1

/// @brief This subroutine is the optional fourth step of the QZ algorithm for 
/// solving generalized matrix eigenvalue problems.
///
/// This subroutine accepts a pair of real matrices, one of them in 
/// quasi-triangular form (in which each 2-by-2 block corresponds to a pair of 
/// complex eigenvalues) and the other in upper triangular form.  It computes 
/// the eigenvectors of the triangular problem and transforms the results back 
/// to the original coordinate system.  It is usually preceded by `qzhes`, 
/// `qzit`, and `qzval`.
///
/// @param[in] nm Must be set to the row dimension of two-dimensional array 
/// parameters as declared in the calling program dimension statement.
///
/// @param[in] n The order of the matrices.
///
/// @param[in] a Contains a real upper quasi-triangular matrix.  Not altered by 
/// this subroutine.  Its subdiagonal elements provide information about the 
/// storage of the complex eigenvectors.
///
/// @param[in] b Contains a real upper triangular matrix.  In addition, 
/// location `b(n,1)` contains the tolerance quantity `epsb` computed and saved 
/// in `qzit`.
///
/// @param[in] alfr,alfi,beta  Vectors with components whose ratios 
/// \f$\frac{\mathrm{alfr} + i \times \mathrm{alfi} } {\mathrm{beta}}\f$ are 
/// the generalized eigenvalues.  They are usually obtained from `qzval`.
///
/// @param[in] z Contains the transformation matrix produced in the reductions 
/// by `qzhes`, `qzit`, and `qzval`, if performed.  If the eigenvectors of the 
/// triangular problem are desired, `z` must contain the identity matrix.
///
/// @param[out] b Destroyed.
///
/// @param[out] z Contains the real and imaginary parts of the eigenvectors.  
/// Each eigenvector is normalized so that the modulus of its largest component 
/// is 1.0.
///     - if `alfi(i) == 0.0`, the i-th eigenvalue is real and the i-th column
///       of z contains its eigenvector.
///     - if `alfi(i) != 0.0`, the i-th eigenvalue is complex.
///     - if `alfi(i) > 0.0`, the eigenvalue is the first of a complex pair and
///       the i-th and (i+1)-th columns of `z` contain its eigenvector.
///     - if `alfi(i) < 0.0`, the eigenvalue is the second of a complex pair
///       and the (i-1)-th and i-th columns of z contain the conjugate of its 
///       eigenvector.
///
/// Reference: Siam J. Numer. Anal. 10, 241-256(1973) by Moler and Stewart.
///
/// Questions and comments should be directed to Burton S. Garbow,
/// Mathematics and Computer Science Division, Argonne National Laboratory.
/// This version dated August 1983.  Fable was used to convert the original 
/// Fortran source into C++.


void qzvec(
		int const& nm,
		int const& n,
		arr_cref<double, 2> a,
		arr_ref<double, 2> b,
		arr_cref<double> alfr,
		arr_cref<double> alfi,
		arr_cref<double> beta,
		arr_ref<double, 2> z) {

  a(dimension(nm, n));
  b(dimension(nm, n));
  alfr(dimension(n));
  alfi(dimension(n));
  beta(dimension(n));
  z(dimension(nm, n));
  //double epsb = fem::double0;
  //int isw = fem::int0;
  int nn = fem::int0;
  int en = fem::int0;
  int na = fem::int0;
  int m = fem::int0;
  double alfm = fem::double0;
  double betm = fem::double0;
  int ii = fem::int0;
  int i = fem::int0;
  double w = fem::double0;
  double r = fem::double0;
  int j = fem::int0;
  double zz = fem::double0;
  double s = fem::double0;
  double t = fem::double0;
  double x = fem::double0;
  double y = fem::double0;
  double q = fem::double0;
  double almr = fem::double0;
  double almi = fem::double0;
  int enm2 = fem::int0;
  double w1 = fem::double0;
  double ra = fem::double0;
  double sa = fem::double0;
  double x1 = fem::double0;
  double z1 = fem::double0;
  double tr = fem::double0;
  double ti = fem::double0;
  double dr = fem::double0;
  double di = fem::double0;
  double rr = fem::double0;
  double d = fem::double0;
  double t1 = fem::double0;
  double t2 = fem::double0;
  int jj = fem::int0;
  int k = fem::int0;

  double epsb = b(n, 1);
  int isw = 1;
  //C     .......... for en=n step -1 until 1 do -- ..........
  FEM_DO_SAFE(nn, 1, n) {
    en = n + 1 - nn;
    na = en - 1;
    if (isw == 2) {
      goto statement_795;
    }
    if (alfi(en) != 0.0e0) {
      goto statement_710;
    }
    //C     .......... real vector ..........
    m = en;
    b(en, en) = 1.0e0;
    if (na == 0) {
      goto statement_800;
    }
    alfm = alfr(m);
    betm = beta(m);
    //C     .......... for i=en-1 step -1 until 1 do -- ..........
    FEM_DO_SAFE(ii, 1, na) {
      i = en - ii;
      w = betm * a(i, i) - alfm * b(i, i);
      r = 0.0e0;

      FEM_DO_SAFE(j, m, en) {
        r += (betm * a(i, j) - alfm * b(i, j)) * b(j, en);
      }

      if (i == 1 || isw == 2) {
        goto statement_630;
      }
      if (betm * a(i, i - 1) == 0.0e0) {
        goto statement_630;
      }
      zz = w;
      s = r;
      goto statement_690;
      statement_630:
      m = i;
      if (isw == 2) {
        goto statement_640;
      }
      //C     .......... real 1-by-1 block ..........
      t = w;
      if (w == 0.0e0) {
        t = epsb;
      }
      b(i, en) = -r / t;
      goto statement_700;
      //C     .......... real 2-by-2 block ..........
      statement_640:
      x = betm * a(i, i + 1) - alfm * b(i, i + 1);
      y = betm * a(i + 1, i);
      q = w * zz - x * y;
      t = (x * s - zz * r) / q;
      b(i, en) = t;
      if (fem::dabs(x) <= fem::dabs(zz)) {
        goto statement_650;
      }
      b(i + 1, en) = (-r - w * t) / x;
      goto statement_690;
      statement_650:
      b(i + 1, en) = (-s - y * t) / zz;
      statement_690:
      isw = 3 - isw;
      statement_700:;
    }
    //C     .......... end real vector ..........
    goto statement_800;
    //C     .......... complex vector ..........
    statement_710:
    m = na;
    almr = alfr(m);
    almi = alfi(m);
    betm = beta(m);
    //C     .......... last vector component chosen imaginary so that
    //C                eigenvector matrix is triangular ..........
    y = betm * a(en, na);
    b(na, na) = -almi * b(en, en) / y;
    b(na, en) = (almr * b(en, en) - betm * a(en, en)) / y;
    b(en, na) = 0.0e0;
    b(en, en) = 1.0e0;
    enm2 = na - 1;
    if (enm2 == 0) {
      goto statement_795;
    }
    //C     .......... for i=en-2 step -1 until 1 do -- ..........
    FEM_DO_SAFE(ii, 1, enm2) {
      i = na - ii;
      w = betm * a(i, i) - almr * b(i, i);
      w1 = -almi * b(i, i);
      ra = 0.0e0;
      sa = 0.0e0;

      FEM_DO_SAFE(j, m, en) {
        x = betm * a(i, j) - almr * b(i, j);
        x1 = -almi * b(i, j);
        ra += x * b(j, na) - x1 * b(j, en);
        sa += x * b(j, en) + x1 * b(j, na);
      }

      if (i == 1 || isw == 2) {
        goto statement_770;
      }
      if (betm * a(i, i - 1) == 0.0e0) {
        goto statement_770;
      }
      zz = w;
      z1 = w1;
      r = ra;
      s = sa;
      isw = 2;
      goto statement_790;
      statement_770:
      m = i;
      if (isw == 2) {
        goto statement_780;
      }
      //C     .......... complex 1-by-1 block ..........
      tr = -ra;
      ti = -sa;
      statement_773:
      dr = w;
      di = w1;
      //C     .......... complex divide (t1,t2) = (tr,ti) / (dr,di) ..........
      statement_775:
      if (fem::dabs(di) > fem::dabs(dr)) {
        goto statement_777;
      }
      rr = di / dr;
      d = dr + di * rr;
      t1 = (tr + ti * rr) / d;
      t2 = (ti - tr * rr) / d;
      switch (isw) {
        case 1: goto statement_787;
        case 2: goto statement_782;
        default: break;
      }
      statement_777:
      rr = dr / di;
      d = dr * rr + di;
      t1 = (tr * rr + ti) / d;
      t2 = (ti * rr - tr) / d;
      switch (isw) {
        case 1: goto statement_787;
        case 2: goto statement_782;
        default: break;
      }
      //C     .......... complex 2-by-2 block ..........
      statement_780:
      x = betm * a(i, i + 1) - almr * b(i, i + 1);
      x1 = -almi * b(i, i + 1);
      y = betm * a(i + 1, i);
      tr = y * ra - w * r + w1 * s;
      ti = y * sa - w * s - w1 * r;
      dr = w * zz - w1 * z1 - x * y;
      di = w * z1 + w1 * zz - x1 * y;
      if (dr == 0.0e0 && di == 0.0e0) {
        dr = epsb;
      }
      goto statement_775;
      statement_782:
      b(i + 1, na) = t1;
      b(i + 1, en) = t2;
      isw = 1;
      if (fem::dabs(y) > fem::dabs(w) + fem::dabs(w1)) {
        goto statement_785;
      }
      tr = -ra - x * b(i + 1, na) + x1 * b(i + 1, en);
      ti = -sa - x * b(i + 1, en) - x1 * b(i + 1, na);
      goto statement_773;
      statement_785:
      t1 = (-r - zz * b(i + 1, na) + z1 * b(i + 1, en)) / y;
      t2 = (-s - zz * b(i + 1, en) - z1 * b(i + 1, na)) / y;
      statement_787:
      b(i, na) = t1;
      b(i, en) = t2;
      statement_790:;
    }
    //C     .......... end complex vector ..........
    statement_795:
    isw = 3 - isw;
    statement_800:;
  }
  //C     .......... end back substitution.
  //C                transform to original coordinate system.
  //C                for j=n step -1 until 1 do -- ..........
  FEM_DO_SAFE(jj, 1, n) {
    j = n + 1 - jj;

    FEM_DO_SAFE(i, 1, n) {
      zz = 0.0e0;

      FEM_DO_SAFE(k, 1, j) {
        zz += z(i, k) * b(k, j);
      }

      z(i, j) = zz;
    }
  }
  //C     .......... normalize so that modulus of largest
  //C                component of each vector is 1.
  //C                (isw is 1 initially from before) ..........
  FEM_DO_SAFE(j, 1, n) {
    d = 0.0e0;
    if (isw == 2) {
      goto statement_920;
    }
    if (alfi(j) != 0.0e0) {
      goto statement_945;
    }

    FEM_DO_SAFE(i, 1, n) {
      if (fem::dabs(z(i, j)) > d) {
        d = fem::dabs(z(i, j));
      }
    }

    FEM_DO_SAFE(i, 1, n) {
      z(i, j) = z(i, j) / d;
    }

    goto statement_950;

    statement_920:
    FEM_DO_SAFE(i, 1, n) {
      r = fem::dabs(z(i, j - 1)) + fem::dabs(z(i, j));
      if (r != 0.0e0) {
        r = r * fem::dsqrt(fem::pow2((z(i, j - 1) / r)) + fem::pow2((z(i,
          j) / r)));
      }
      if (r > d) {
        d = r;
      }
    }

    FEM_DO_SAFE(i, 1, n) {
      z(i, j - 1) = z(i, j - 1) / d;
      z(i, j) = z(i, j) / d;
    }

    statement_945:
    isw = 3 - isw;
    statement_950:;
  }

}

// void rgg( {{{1

/// @brief This subroutine calls the recommended sequence of
/// subroutines from the eigensystem subroutine package (EISPACK)
/// to find the eigenvalues and eigenvectors (if desired)
/// for the real general generalized eigenproblem \f$ Ax = \lambda Bx \f$.
///
/// @param[in] nm  Must be set to the row dimension of the two-dimensional 
/// array parameters as declared in the calling program dimension statement.
///
/// @param[in] n The order of the matrices \f$A\f$ and \f$B\f$.
///
/// @param[in] a A real general matrix.
///
/// @param[in] b A real general matrix.
///
/// @param[in] matz  An integer variable set equal to zero if only eigenvalues 
/// are desired.  otherwise it is set to any non-zero integer for both 
/// eigenvalues and eigenvectors.
///
/// @param[out] alfr The real parts of the numerators of the eigenvalues.
///
/// @param[out] alfi The imaginary parts of the numerators of the eigenvalues.
///
/// @param[out] beta The denominators of the eigenvalues, which are thus given 
/// by the ratios \f$\frac{\mathrm{alfr} + i \times \mathrm{alfi}} 
/// {\mathrm{beta}}\f$.  Complex conjugate pairs of eigenvalues appear 
/// consecutively with the eigenvalue having the positive imaginary part first.
///
/// @param[out] z The real and imaginary parts of the eigenvectors if matz is 
/// not zero.  If the j-th eigenvalue is real, the j-th column of z contains 
/// its eigenvector.  If the j-th eigenvalue is complex with positive imaginary 
/// part, the j-th and (j+1)-th columns of z contain the real and imaginary 
/// parts of its eigenvector.  The conjugate of this vector is the eigenvector 
/// for the conjugate eigenvalue.
///
/// @param[out] ierr An integer output variable set equal to an error 
/// completion code described in the documentation for qzit.  The normal 
/// completion code is zero.
///
/// Questions and comments should be directed to Burton S. Garbow,
/// Mathematics and Computer Science Division, Argonne National Laboratory.
/// This version dated August 1983.  Fable was used to convert the original 
/// Fortran source into C++.

void rgg(
		int const& nm,
		int const& n,
		arr_ref<double, 2> a,
		arr_ref<double, 2> b,
		arr_ref<double> alfr,
		arr_ref<double> alfi,
		arr_ref<double> beta,
		int const& matz,
		arr_ref<double, 2> z,
		int& ierr) {

  a(dimension(nm, n));
  b(dimension(nm, n));
  alfr(dimension(n));
  alfi(dimension(n));
  beta(dimension(n));
  z(dimension(nm, n));
  bool tf = fem::bool0;

  if (n <= nm) {
    goto statement_10;
  }
  ierr = 10 * n;
  goto statement_50;

  statement_10:
  if (matz != 0) {
    goto statement_20;
  }
  //C     .......... find eigenvalues only ..........
  tf = false;
  qzhes(nm, n, a, b, tf, z);
  qzit(nm, n, a, b, 0.0e0, tf, z, ierr);
  qzval(nm, n, a, b, alfr, alfi, beta, tf, z);
  goto statement_50;
  //C     .......... find both eigenvalues and eigenvectors ..........
  statement_20:
  tf = true;
  qzhes(nm, n, a, b, tf, z);
  qzit(nm, n, a, b, 0.0e0, tf, z, ierr);
  qzval(nm, n, a, b, alfr, alfi, beta, tf, z);
  if (ierr != 0) {
    goto statement_50;
  }
  qzvec(nm, n, a, b, alfr, alfi, beta, z);
  statement_50:;
}
// }}}1

} // namespace linear_algebra
} // namespace numeric
