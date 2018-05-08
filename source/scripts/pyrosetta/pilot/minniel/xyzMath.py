"""
Easy 3D Linear Algebra, like xyz\* in rosetta
"""
from random    import gauss,uniform
from math      import pi,sqrt,sin,cos,acos,asin,atan2,degrees,radians,copysign
from itertools import chain,product,izip
import math
import operator as op

EPS = 0.00001
SQRTEPS = sqrt(EPS)
ATET = 54.735610317245360079 # asin(sr2/sr3)
AOCT = 35.264389682754668343 # asin(sr1/sr3)
AICS = 20.89774264557		 # asin(G/2/sr3)

def isint  (x): return type(x) is int
def isfloat(x): return type(x) is float
def isnum  (x): return isint(x) or isfloat(x)
def isiter (x): return hasattr(x,"__iter__")
def islist (x): return type(x) is list
def istuple(x): return type(x) is tuple

def isvec  (x): return type(x) is Vec
def ismat  (x): return type(x) is Mat
def isxform(x): return type(x) is Xform

# def isvec  (x): return str(type(x)).split('.')[-1] == "Vec'>"
# def ismat  (x): return str(type(x)).split('.')[-1] == "Mat'>"
# def isxform(x): return str(type(x)).split('.')[-1] == "Xform'>"

class Vec(object):
	"""a Vector like xyzVector<Real> in rosetta

	>>> v = Vec(1,2,3)
	>>> print v, 10*v
	(1.000000,2.000000,3.000000) (10.000000,20.000000,30.000000)

	multiplication is a dot prod at the moment
	>>> v*v
	14.0
	>>> assert Vec(1,0,-0) == Vec(1,-0,0)
	"""
	def __init__(self,x=0.0,y=None,z=None):
		if y is None:
			if isnum(x):
				self.x,self.y,self.z = (float(x),)*3
			elif isvec(x):
				self.x,self.y,self.z = x.x,x.y,x.z
			elif isiter(x):
				i = iter(x)
				self.x,self.y,self.z = i.next(),i.next(),i.next()
			else: raise NotImplementedError
		elif z is not None:
			assert isnum(x) and isnum(y) and isnum(z)
			self.x,self.y,self.z = float(x),float(y),float(z)
		else: raise NotImplementedError
		assert isfloat(self.x)
		assert isfloat(self.y)
		assert isfloat(self.z)
	def dot(u,v):     assert isvec(v); return u.x*v.x+u.y*v.y+u.z*v.z
	def normdot(u,v): assert isvec(v); return min(1.0,max(-1.0,u.dot(v)/u.length()/v.length()))
	def angle(u,v):
		assert isvec(v)
		d = u.normdot(v)
		if d > 1.0-EPS: return 0.0;
		if d < EPS-1.0: return pi
		return acos(d)
	def angle_degrees(u,v): return degrees(u.angle(v))
	def lineangle(u,v):
		assert isvec(v);
		if u.length() < SQRTEPS or v.length < SQRTEPS: return 0.0
		ang = abs(acos( u.normdot(v) ))
		return ang if ang < pi/2.0 else pi-ang
	def lineangle_degrees(u,v):
		if isvec(v): return degrees(u.lineangle(v))
		raise NotImplementedError
	def length(u):         return sqrt(u.dot(u))
	def length_squared(u): return      u.dot(u)
	def distance(u,v):         assert isvec(v); return (u-v).length()
	def distance_squared(u,v): assert isvec(v); return (u-v).length_squared()
	def cross(u,v): assert isvec(v); return Vec(u.y*v.z-u.z*v.y,u.z*v.x-u.x*v.z,u.x*v.y-u.y*v.x)
	def __mul__(u,a):
		if   isnum(a): return Vec(u.x*a,u.y*a,u.z*a)
		elif isvec(a): return u.dot(a)
		else: return a.__rmul__(u)
	def __rmul__(u,a): return u*a
	def __add__(u,v):
		if isvec(v): return Vec(u.x+v.x,u.y+v.y,u.z+v.z)
		return v.__radd__(u)
	def __sub__(u,v):
		if isvec(v): return Vec(u.x-v.x,u.y-v.y,u.z-v.z)
		return v.__rsub__(u)
	def __neg__(u): return Vec(-u.x,-u.y,-u.z)
	def __div__(u,a): return u*(1.0/a)
	def __str__(self): return "(%f,%f,%f)"%(self.x,self.y,self.z)
	def __repr__(self): return "Vec( %f, %f, %f )"%(self.x,self.y,self.z)
	def normalize(u):
		l = u.length()
		u.x /= l
		u.y /= l
		u.z /= l
	def normalized(u):
		v = Vec(u)
		v.normalize()
		return v
	def outer(u,v):
		assert isvec(v)
		return Mat( u.x*v.x, u.x*v.y, u.x*v.z,
						u.y*v.x, u.y*v.y, u.y*v.z,
						u.z*v.x, u.z*v.y, u.z*v.z      )
	def __eq__(self,other):
		assert isvec(other)
		return ( abs(self.x-other.x) < EPS and
					abs(self.y-other.y) < EPS and
					abs(self.z-other.z) < EPS )
	def rounded(self,sd):
		return Vec(round(self.x,sd), round(self.y,sd), round(self.z,sd) )
	def unit(v):
		if   abs(v.x) > SQRTEPS: return v/v.x
		elif abs(v.y) > SQRTEPS: return v/v.y
		elif abs(v.z) > SQRTEPS: return v/v.z
	def __len__(v):
		return 3
	def abs(v):
		return Vec(abs(v.x),abs(v.y),abs(v.z))
	def __getitem__(v,i):
		if i is 0: return v.x
		if i is 1: return v.x
		if i is 2: return v.x
		raise IndexError
	def tuple(v):
		return (v.x,v.y,v.z)
	def key(v):
		return v.abs().unit().rounded(6).tuple()

Ux = Vec(1,0,0)
Uy = Vec(0,1,0)
Uz = Vec(0,0,1)
V0 = Vec(0,0,0)

def randvec(n=None):
	if n is None: return Vec(gauss(0,1),gauss(0,1),gauss(0,1))
	return [Vec(gauss(0,1),gauss(0,1),gauss(0,1)) for i in range(n)]

def randveccube(r=1.0):
	return Vec(uniform(-1,1),uniform(-1,1),uniform(-1,1))*r

def randvecball(r=1.0):
	v =randveccube(r)
	while v.length_squared() > r*r: v = randveccube(r)
	return v

def randnorm(n=None):
	"""
	>>> assert abs(randnorm().length()-1.0) < 0.0000001
	"""
	if n is None: return randvec().normalized()
	return (randvec().normalized() for i in range(n))

def coplanar(x1,x2,x3,x4):
	"""
	>>> u,v,w = randvec(3)
	>>> a,b,c = (gauss(0,10) for i in range(3))
	>>> assert     coplanar(u, v, w, u + a*(u-v) + b*(v-w) + c*(w-u) )
	>>> assert not coplanar(u, v, w, u + a*(u-v) + b*(v-w) + c*(w-u) + randvec().cross(u-v) )
	"""
	return abs((x3-x1).dot((x2-x1).cross(x4-x3))) < SQRTEPS

def rmsd(l,m):
	"""
	>>> l,m = randvec(6),randvec(6)
	>>> rmsd(l,l)
	0.0
	"""
	rmsd = 0.0
	for u,v in izip(l,m):
		rmsd += u.distance_squared(v)
	return sqrt(rmsd)

class Mat(object):
	"""docstring for Mat

	>>> m = Mat(2,0,0,0,1,0,0,0,1)
	>>> print m
	Mat[ (2.000000,0.000000,0.000000), (0.000000,1.000000,0.000000), (0.000000,0.000000,1.000000) ]
	>>> print m*m
	Mat[ (4.000000,0.000000,0.000000), (0.000000,1.000000,0.000000), (0.000000,0.000000,1.000000) ]
	>>> print Mat(*range(1,10)) * Mat(*range(10,19))
	Mat[ (84.000000,90.000000,96.000000), (201.000000,216.000000,231.000000), (318.000000,342.000000,366.000000) ]
	>>> assert Mat(0.0,1.0,2.0,3,4,5,6,7,8) == Mat(-0,1,2,3,4,5.0,6.0,7.0,8.0)
	>>> print Mat(100,2,3,4,5,6,7,8,9).det()
	-297.0
	>>> m = Mat(100,2,3,4,5,6,7,8,9)
	>>> assert m * ~m == Imat
	"""
	def __init__(self, xx=None, xy=None, xz=None, yx=None, yy=None, yz=None, zx=None, zy=None, zz=None):
		super(Mat, self).__init__()
		if xx is None: # identity default
			self.xx, self.xy, self.xz = 1.0,0.0,0.0
			self.yx, self.yy, self.yz = 0.0,1.0,0.0
			self.zx, self.zy, self.zz = 0.0,0.0,1.0
		elif xy is None and ismat(xx):
			self.xx, self.xy, self.xz = xx.xx, xx.xy, xx.xz
			self.yx, self.yy, self.yz = xx.yx, xx.yy, xx.yz
			self.zx, self.zy, self.zz = xx.zx, xx.zy, xx.zz
		elif yx is None and isvec(xx) and isvec(xy) and isvec(xz):
			self.xx, self.xy, self.xz = xx.x, xy.x, xz.x
			self.yx, self.yy, self.yz = xx.y, xy.y, xz.y
			self.zx, self.zy, self.zz = xx.z, xy.z, xz.z
		elif isnum(xx):
			self.xx, self.xy, self.xz = float(xx), float(xy), float(xz)
			self.yx, self.yy, self.yz = float(yx), float(yy), float(yz)
			self.zx, self.zy, self.zz = float(zx), float(zy), float(zz)
		else:
			raise NotImplementedError
		assert isfloat(self.xx) and isfloat(self.xy) and isfloat(self.xz)
		assert isfloat(self.yx) and isfloat(self.yy) and isfloat(self.yz)
		assert isfloat(self.zx) and isfloat(self.zy) and isfloat(self.zz)
	def row(m,i):
		assert isint(i)
		if   i is 0: return Vec(m.xx,m.xy,m.xz)
		elif i is 1: return Vec(m.yx,m.yy,m.yz)
		elif i is 2: return Vec(m.zx,m.zy,m.zz)
		else: assert 0 <= i and i <= 2
	def col(m,i):
		assert isint(i)
		if   i is 0: return Vec(m.xx,m.yx,m.zx)
		elif i is 1: return Vec(m.xy,m.yy,m.zy)
		elif i is 2: return Vec(m.xz,m.yz,m.zz)
		else: assert 0 <= i and i <= 2
	def rowx(m): return m.row(0)
	def rowy(m): return m.row(1)
	def rowz(m): return m.row(2)
	def colx(m): return m.col(0)
	def coly(m): return m.col(1)
	def colz(m): return m.col(2)
	def __invert__(m):
		"""
		>>> from random import random
		>>> for i in range(10):
		...     m = random()*Mat(random(),random(),random(),random(),random(),random(),random(),random(),random())
		...     assert ~m*m == Imat
		...     assert m*~m == Imat
		"""
		return Mat(   m.zz*m.yy-m.zy*m.yz  , -(m.zz*m.xy-m.zy*m.xz) ,   m.yz*m.xy-m.yy*m.xz  ,
						-(m.zz*m.yx-m.zx*m.yz) ,   m.zz*m.xx-m.zx*m.xz  , -(m.yz*m.xx-m.yx*m.xz) ,
						  m.zy*m.yx-m.zx*m.yy  , -(m.zy*m.xx-m.zx*m.xy) ,   m.yy*m.xx-m.yx*m.xy  ) / m.det()
	def __mul__(m,rhs):
		if   isnum(rhs): return Mat( rhs*m.xx, rhs*m.xy, rhs*m.xz, rhs*m.yx, rhs*m.yy, rhs*m.yz, rhs*m.zx, rhs*m.zy, rhs*m.zz )
		elif isvec(rhs): return Vec( m.rowx()*rhs, m.rowy()*rhs, m.rowz()*rhs )
		elif ismat(rhs): return Mat( m.rowx()*rhs.colx(), m.rowx()*rhs.coly(), m.rowx()*rhs.colz(),
							m.rowy()*rhs.colx(), m.rowy()*rhs.coly(), m.rowy()*rhs.colz(),
							m.rowz()*rhs.colx(), m.rowz()*rhs.coly(), m.rowz()*rhs.colz() )
		else: return rhs.__rmul__(m)
	def __rmul__(m,v):
		if   isnum(v): return m*v
		elif isvec(v): return Vec( m.colx()*v, m.coly()*v, m.colz()*v )
	def __div__(m,v): return m*(1/v)
	def __add__(m,v):
		if   isnum(v): return Mat(v   +m.xx,v   +m.xy,v   +m.xz,v   +m.yx,v   +m.yy,v   +m.yz,v   +m.zx,v   +m.zy,v   +m.zz)
		elif ismat(v): return Mat(v.xx+m.xx,v.xy+m.xy,v.xz+m.xz,v.yx+m.yx,v.yy+m.yy,v.yz+m.yz,v.zx+m.zx,v.zy+m.zy,v.zz+m.zz)
		else: return v.__radd__(m)
	def __sub__(m,v): return m + -v
	def __neg__(m): return m * -1
	def __str__(m): return "Mat[ %s, %s, %s ]" % (str(m.rowx()),str(m.rowy()),str(m.rowz()))
	def __repr__(m): return "Mat( %s, %s, %s )" % (repr(m.colx()),repr(m.coly()),repr(m.colz()))
	def transpose(m):
		m = Mat( m.xx, m.yx, m.zx, m.xy, m.yy, m.zy, m.xz, m.yz, m.zz )
	def transposed(m): return Mat( m.xx, m.yx, m.zx, m.xy, m.yy, m.zy, m.xz, m.yz, m.zz )
	def det(m):
			 #   a11  (a33  a22- a32  a23)- a21 ( a33  a12- a32  a13)+ a31(  a23  a12- a22  a13)
		return m.xx*(m.zz*m.yy-m.zy*m.yz)-m.yx*(m.zz*m.xy-m.zy*m.xz)+m.zx*(m.yz*m.xy-m.yy*m.xz)
	def trace(m): return m.xx+m.yy+m.zz
	def add_diagonal(m,v): return Mat( v.x+m.xx, m.xy, m.xz, m.yx, v.y+m.yy, m.yz, m.zx, m.zy, v.z+m.zz )
	def is_rotation(m): return (m.colx().isnormal() and m.coly().isnormal() and m.colz().isnormal() and
				  m.rowx().isnormal() and m.rowy().isnormal() and m.rowz().isnormal()   )
	def __eq__(self,other): return ( abs(self.xx-other.xx) < EPS and
					abs(self.xy-other.xy) < EPS and
					abs(self.xz-other.xz) < EPS and
					abs(self.yx-other.yx) < EPS and
					abs(self.yy-other.yy) < EPS and
					abs(self.yz-other.yz) < EPS and
					abs(self.zx-other.zx) < EPS and
					abs(self.zy-other.zy) < EPS and
					abs(self.zz-other.zz) < EPS )
	def rotation_axis(R):
		"""
		>>> axis ,ang  = randnorm(),uniform(-pi,pi)
		>>> axis2,ang2 = rotation_matrix(axis,ang).rotation_axis()
		>>> assert abs( abs(ang) - abs(ang2) ) < EPS
		>>> assert axis == axis2 * copysign(1,ang*ang2)
		"""
		cos_theta = sin_cos_range((R.trace()-1.0)/2.0);
		if cos_theta > -1.0+EPS and cos_theta < 1.0-EPS:
			x = ( 1.0 if R.zy > R.yz else -1.0 ) * sqrt( max(0.0, ( R.xx - cos_theta ) / ( 1.0 - cos_theta ) ) )
			y = ( 1.0 if R.xz > R.zx else -1.0 ) * sqrt( max(0.0, ( R.yy - cos_theta ) / ( 1.0 - cos_theta ) ) )
			z = ( 1.0 if R.yx > R.xy else -1.0 ) * sqrt( max(0.0, ( R.zz - cos_theta ) / ( 1.0 - cos_theta ) ) )
			theta = acos( cos_theta );
			assert abs( x*x + y*y + z*z - 1 ) <= 0.01
			return Vec(x,y,z),theta
		elif cos_theta >= 1.0-EPS: return Vec(1.0,0.0,0.0),0.0
		else:
			nnT = (R+Imat)/2.0
			x,y,z = 0.0,0.0,0.0;
			if nnT.xx > EPS:
				x = sqrt( nnT.xx )
				y = nnT.yx / x
				z = nnT.zx / x
			elif nnT.yy > EPS:
				x = 0
				y = sqrt(nnT.yy)
				z = nnT.zy / y
			else:
				assert( nnT.zz > EPS );
				x = 0
				y = 0
				z = sqrt( nnT.zz )
			assert abs( x*x + y*y + z*z - 1.0 ) <= 0.01
			return Vec( x, y, z ),pi
	def euler_angles(self):
		FLOAT_PRECISION = 1e-5
		if self.zz >= 1 - FLOAT_PRECISION:
			e1 = math.acos( sin_cos_range( self.xx ) )
			e2 = 0.0
			e3 = 0.0
			return Vec(e1,e2,e3)
		if self.zz <= -1 + FLOAT_PRECISION:
			e1 = math.acos( sin_cos_range( self.xx ) )
			e2 = 0.0
			e3 = math.pi
			return Vec(e1,e2,e3)
		pos_sin_theta = math.sqrt( 1 - self.zz*self.zz ) # sin2theta = 1 - cos2theta.
		#  two values are possible here: my convention is to use positive theta only.
		#  corresponding theta between [0,pi/2] -> [0,90] since st > 0
		#  and asin returns value between [-pi/2, pi/2]
		e3 = math.asin( pos_sin_theta )
		#  decide whether the actual positive theta is between [pi/2, pi[ using the value of cos(theta)
		#  which happens to be the matrix element self.zz (and is thus signed).
		if self.zz < 0:
			e3 = math.pi - e3
		# e1 = math.atan2(  -UU(1,3), UU(2,3) )   #  between -Pi and Pi -> [-180,180]
		# e2 = math.atan2( UU(3,1), UU(3,2) )     #  between -Pi and Pi -> [-180, 180]
		#  this is atan( sin_phi * c, cos_phi * c  ) as opposed to Alex's atan( -sin_phi * c, -cos_phi * c ).
		e1 = math.atan2( self.zx, -self.zy )
		e2 = math.atan2( self.xz,  self.yz )
		return Vec(e1,e2,e3)

	def from_euler_angles(self,euler):
		phi   = euler.x()
		psi   = euler.y()
		theta = euler.z()
		cos_phi   = math.cos( phi )
		sin_phi   = math.sin( phi )
		cos_psi   = math.cos( psi )
		sin_psi   = math.sin( psi )
		cos_theta = math.cos( theta )
		sin_theta = math.sin( theta )
		self.xx =  cos_psi * cos_phi - cos_theta * sin_phi * sin_psi
		self.xy =  cos_psi * sin_phi + cos_theta * cos_phi * sin_psi
		self.xz =  sin_psi * sin_theta
		self.yx = -sin_psi * cos_phi - cos_theta * sin_phi * cos_psi
		self.yy = -sin_psi * sin_phi + cos_theta * cos_phi * cos_psi
		self.yz =  cos_psi * sin_theta
		self.zx =  sin_theta * sin_phi
		self.zy = -sin_theta * cos_phi
		self.zz =        cos_theta
		return self

Imat = Mat(1,0,0,0,1,0,0,0,1)

def projection_matrix(v):
	m = Mat( v.x * v.x, v.x * v.y, v.x * v.z, v.y * v.x, v.y * v.y, v.y * v.z, v.z * v.x, v.z * v.y, v.z * v.z )
	return m / v.dot(v)

def proj(u,v):
	"""
	>>> u = Vec(1,0,0); v = Vec(1,1,1)
	>>> proj(u,v)
	Vec( 1.000000, 0.000000, 0.000000 )
	"""
	return projection_matrix(u)*v

def projperp(u,v):
	"""
	>>> u = Vec(1,0,0); v = Vec(1,1,1)
	>>> projperp(u,v)
	Vec( 0.000000, 1.000000, 1.000000 )
	"""
	return v - proj(u,v)

def rotation_matrix(axis,angle):
	n = axis.normalized()
	sin_theta = sin( angle )
	cos_theta = cos( angle )
	R = projection_matrix(n)
	R *= 1.0 - cos_theta
	R.xx += cos_theta;       R.xy -= sin_theta * n.z; R.xz += sin_theta * n.y
	R.yx += sin_theta * n.z; R.yy += cos_theta;       R.yz -= sin_theta * n.x
	R.zx -= sin_theta * n.y; R.zy += sin_theta * n.x; R.zz += cos_theta
	return R;

def rotation_matrix_degrees(axis,angle):
	""" get a rotation matrix

	>>> rx180 = rotation_matrix_degrees(Vec(1,0,0),180.0)
	>>> rx90  = rotation_matrix_degrees(Vec(1,0,0),90.0)
	>>> print rx90*rx90 == rx180
	True
	>>> r = rotation_matrix_degrees(Vec(1,0,0),45.0)
	>>> print r
	Mat[ (1.000000,0.000000,0.000000), (0.000000,0.707107,-0.707107), (0.000000,0.707107,0.707107) ]
	>>> assert r*r == rx90
	>>> assert r*r*r*r == rx180
	>>> assert r*r*r*r*r*r*r*r == Imat
	>>> assert ~r == r.transposed()

	>>> ang = uniform(0,1)*360.0-180.0
	>>> v = randvec()
	>>> axs = randnorm()
	>>> while(abs(v.dot(axs))>0.9): axs = randnorm()
	>>> u = rotation_matrix_degrees(projperp(v,axs),ang)*v
	>>> assert abs(u.angle_degrees(v)-abs(ang)) < SQRTEPS
	>>> test_rotation_mat()
	test_rotation_mat PASS
	"""
	return rotation_matrix(axis,radians(angle))

def test_rotation_mat():
	import random
	for i in range(100):
		a0 = randnorm()
		t0 = uniform(-pi,pi)
		a,t = rotation_matrix(a0,t0).rotation_axis()
		if t0 < 0.01: continue
		if abs(t-pi) < EPS:
			if (abs(a.x-a0.x) < 0.001 and abs(a.y-a0.y) < 0.001 and abs(a.z-a0.z) < 0.001) or \
				(abs(a.x+a0.x) < 0.001 and abs(a.y+a0.y) < 0.001 and abs(a.z+a0.z) < 0.001):
				continue
			else:
				print a0
				print a
				return False
		if not abs(t-t0) < EPS or not (a.normalized()-a0.normalized()).length() < EPS:
			print a0.normalized(), t0
			print a.normalized() , t
			print "FAIL"
			return
	print "test_rotation_mat PASS"


def randrot(n=None):
 	if n is None: return rotation_matrix_degrees(randnorm(),uniform(0,1)*360)
 	return (rotation_matrix_degrees(randnorm(),uniform(0,1)*360) for i in range(n))

class Xform(object):
	"""Coordinate frame like rosetta Xform, behaves also as a rosetta Stub

	>>> x = Xform(R=Imat,t=Uz)
	>>> print x
	Xform( Mat[ (1.000000,0.000000,0.000000), (0.000000,1.000000,0.000000), (0.000000,0.000000,1.000000) ], (0.000000,0.000000,1.000000) )
	>>> assert (x*x) == Xform(R=Imat,t=2*Uz)
	>>> x = Xform(R=rotation_matrix_degrees(Vec(1,0,0),90.0),t=Vec(0,0,0))
	>>> print x
	Xform( Mat[ (1.000000,0.000000,0.000000), (0.000000,0.000000,-1.000000), (0.000000,1.000000,0.000000) ], (0.000000,0.000000,0.000000) )
	>>> assert x*x*x*x == Ixform
	>>> x.t = Ux
	>>> assert x*x*x*x == Xform(R=Imat,t=4*Ux)
	>>> x.t = Uz
	>>> print x
	Xform( Mat[ (1.000000,0.000000,0.000000), (0.000000,0.000000,-1.000000), (0.000000,1.000000,0.000000) ], (0.000000,0.000000,1.000000) )
	>>> assert x               == Xform(R=rotation_matrix_degrees(Ux, 90.0),t=Vec(0, 0,1))
	>>> assert x*x             == Xform(R=rotation_matrix_degrees(Ux,180.0),t=Vec(0,-1,1))
	>>> assert x*x*x           == Xform(R=rotation_matrix_degrees(Ux,270.0),t=Vec(0,-1,0))
	>>> assert x*x*x*x         == Xform(R=rotation_matrix_degrees(Ux,  0.0),t=Vec(0, 0,0))
	>>> assert x*x*x*x*x       == Xform(R=rotation_matrix_degrees(Ux, 90.0),t=Vec(0, 0,1))
	>>> assert x*x*x*x*x*x     == Xform(R=rotation_matrix_degrees(Ux,180.0),t=Vec(0,-1,1))
	>>> assert x*x*x*x*x*x*x   == Xform(R=rotation_matrix_degrees(Ux,270.0),t=Vec(0,-1,0))
	>>> assert x*x*x*x*x*x*x*x == Xform(R=rotation_matrix_degrees(Ux,  0.0),t=Vec(0, 0,0))
	>>> x = Xform(rotation_matrix_degrees(Vec(1,2,3),123),Vec(5,7,9))
	>>> assert ~x *  x == Ixform
	>>> assert  x * ~x == Ixform

	Frames / RTs are interchangable:

	>>> fr = Xform(rotation_matrix_degrees(Vec(1,2,3), 65.64),t=Vec(3,2,1))
	>>> to = Xform(rotation_matrix_degrees(Vec(7,5,3),105.44),t=Vec(10,9,8))
	>>> x = to/fr
	>>> assert to/Ixform ==  to
	>>> assert Ixform/fr == ~fr
	>>> assert (to * ~fr) * fr == to
	>>> assert x * fr == to

	>>> a1 = randnorm()
	>>> b1 = randnorm()
	>>> ang = uniform(0,1)*360.0-180.0
	>>> a2 = rotation_matrix_degrees(a1.cross(randnorm()),ang) * a1
	>>> b2 = rotation_matrix_degrees(b1.cross(randnorm()),ang) * b1
	>>> assert abs(angle(a1,a2) - angle(b1,b2)) < EPS
	>>> xa = Xform().from_two_vecs(a1,a2)
	>>> xb = Xform().from_two_vecs(b1,b2)
	>>> assert xa.tolocal(a1) == xb.tolocal(b1)
	>>> assert xa.tolocal(a2) == xb.tolocal(b2)
	>>> assert ~xa*a1 == ~xb*b1
	>>> assert ~xa*a2 == ~xb*b2
	>>> assert xb/xa*a1 == b1
	>>> assert xb/xa*a2 == b2

	add/sub with Vecs:

	>>> X = randxform()
	>>> u,v = randvec(2)
	>>> assert isxform(u+X) and isxform(X+u) and isxform(u-X) and isxform(X-u)
	>>> assert X*v+u == (u+X)*v
	>>> assert X*(v+u) == (X+u)*v
	>>> assert Xform(u)*X*v == (u+X)*v
	>>> assert X*Xform(u)*v == (X+u)*v
	>>> assert X*v-u == (u-X)*v
	>>> assert X*(v-u) == (X-u)*v

	mul,div with Mats:

	>>> R = randrot()
	>>> assert isxform(R*X) and isxform(X*R)
	>>> assert R*X*u == (R*X)*u == R*(X*u)
	>>> assert X*R*u == (X*R)*u == X*(R*u)
	>>> assert Xform(R)*X*u == Xform(R)*(X*u)
	>>> assert X*Xform(R)*u == X*(Xform(R,V0)*u)
	>>> assert X/X*v == v

	mul/div Xforms:

	>>> Y = randxform()
	>>> assert isxform(X/Y) and isxform(X*Y)
	>>> assert X/Y*v == X*~Y*v

	these don't work yet:

	>>> axis,ang,cen = randnorm(),uniform(-pi,pi),randvec()      #doctest: +SKIP
	>>> X = rotation_around(axis,ang,cen)                        #doctest: +SKIP
	>>> axis2,ang2,cen2 = X.rotation_center()                    #doctest: +SKIP
	>>> assert abs( abs(ang) - abs(ang2) ) < EPS                 #doctest: +SKIP
	>>> assert axis == axis2 * copysign(1,ang*ang2)              #doctest: +SKIP
	>>> print cen                                                #doctest: +SKIP
	>>> print cen2                                               #doctest: +SKIP

	>>> x =                Xform( Mat( Vec( 0.816587, -0.306018, 0.489427 ), Vec( 0.245040, 0.951487, 0.186086 ), Vec( -0.522629, -0.032026, 0.851959 ) ), Vec( 1.689794, 1.535762, -0.964428 ) )
	>>> assert repr(x) == "Xform( Mat( Vec( 0.816587, -0.306018, 0.489427 ), Vec( 0.245040, 0.951487, 0.186086 ), Vec( -0.522629, -0.032026, 0.851959 ) ), Vec( 1.689794, 1.535762, -0.964428 ) )"
	"""
	def __init__(self, R=None, t=None):
		super(Xform, self).__init__()
		if isvec(R) and t is None: R,t = Imat,R
		self.R = R if R else Imat
		self.t = t if t else V0
		assert ismat(self.R) and isvec(self.t)
	# def rotation_center(X):
	#    axis,ang = X.rotation_axis()
	#    cen = -(X.R-Imat).transposed()*X.t
	#    return axis,ang,cen
	def from_four_points(s,cen,a,b,c):
		s.t = cen
		e1 = (a-b).normalized()
		e3 = e1.cross(c-b).normalized()
		e2 = e1.cross(e3).normalized()
		# print "from_four_points"
		# print e1
		# print e2
		# print e3
		s.R = Mat(e1.x,e2.x,e3.x,e1.y,e2.y,e3.y,e1.z,e2.z,e3.z)
		return s
	def from_two_vecs(s,a,b):
		e1 = a.normalized()
		e2 = projperp(a,b).normalized()
		e3 = e1.cross(e2)
		return Xform( Mat(e1.x,e2.x,e3.x,e1.y,e2.y,e3.y,e1.z,e2.z,e3.z),V0)
	def tolocal(s,x): return s.R.transposed() * (x - s.t)
	def toglobal(s,x): return (s.R * x) + s.t
	def __invert__(self):
		R = ~self.R
		t = R * -self.t
		return Xform(R,t)
	def __mul__(X,o):
		if   isvec(o):   return X.R * o + X.t
		elif isxform(o): return Xform(X.R*o.R,X.R*(o.t) + X.t)
		elif ismat(o):   return Xform(X.R*o,X.t)
		elif islist(o):  return [X*x for x in o]
		elif istuple(o): return tuple([X*x for x in o])
		elif isiter(o):  return (X*x for x in o)
		else:            return o.__rmul__(X)
	def __rmul__(X,o):
		if ismat(o): return Xform(o*X.R,o*X.t)
		raise NotImplementedError
	def __div__(X,o):
		if isxform(o): return X*~o
		return o.__rdiv__(X)
	def __add__(X,v):
		if isvec(v): return Xform( X.R, X.t + X.R*v )
		return v.__radd__(X)
	def __radd__(X,v):
		if isvec(v): return Xform( X.R, X.t + v )
		raise NotImplementedError
	def __sub__(X,v):
		if isvec(v): return Xform( X.R, X.t - X.R*v )
		return v.__rsub__(X)
	def __rsub__(X,v):
		if isvec(v): return Xform( X.R, X.t - v )
		raise NotImplementedError
	def __eq__(self,other): return self.R==other.R and self.t==other.t
	def __str__ (self): return "Xform( %s, %s )" % (str(self.R),str(self.t))
	def __repr__(self): return "Xform( %s, %s )" % (repr(self.R),repr(self.t))
	def __eq__(X,Y):
		assert isxform(Y)
		return X.R == Y.R and X.t == Y.t
	def rotation_axis(X): return X.R.rotation_axis()
	# def rotation_axis_center(X):
	# 	"""
	# 	R(x-c)+c = Rx+c-Rc = Rx+t  ->  c-Rc=t  (I-R)c=t
	# 	singular! shit!
	# 	>>> x = rotation_around_degrees(Vec(1,1,1),90,Vec(0,2,5))
	# 	>>> print x.rotation_axis_center()

	# 	"""
	# 	axis, ang = X.R.rotation_axis()
	# 	p1before = Vec(0,0,0)
	# 	p2before = Vec(1,0,0)
	# 	p1after = X*p1before
	# 	p2after = X*p2before
	# 	vl1 = (p1after-p1before).cross(axis)
	# 	vl2 = (p2after-p2before).cross(axis)
	# 	cen = line_line_closest_points(vl1,p1after,vl2,p2after)
	# 	return axis,ang,cen
	def pretty(self):
		a,r = self.rotation_axis()
		if self.t.length() > EPS: return "Xform( axis=%s, ang=%f, dir=%s, dis=%f )"%(str(a),degrees(r),str(self.t.normalized()),self.t.length())
		else: return "Xform( axis=%s, ang=%f, dir=%s, dis=%f )"%(str(a),degrees(r),str(V0),0)

Ixform = Xform(Imat,V0)

def read_tokens(f):
	for line in f:
		 for token in line.split():
			  yield token

def read_xforms(fname,N=9e9,start=0):
	xforms = list()
	with open(fname) as fin:
		for i in range(start): fin.readline()
		while True:
			 x = list(chain(*(fin.readline().split() for i in range(4))))
			 if len(xforms) is N or not x: break
			 x = map(float,x)
			 xforms.append( Xform(Mat(*x[:9]),Vec(*x[9:])) )
	return xforms

def stub(cen=None, a=None, b=None, c=None):
	if cen is None: cen = a
	if c   is None: cen,a,b,c = cen,cen,a,b
	# print cen
	# print a
	# print b
	# print c
	return Xform().from_four_points(cen,a,b,c)

def randxform(n=None):
	if n is None: return Xform(randrot(),randvec())
	return (Xform(randrot(),randvec()) for i in range(n))

def rotation_around(axs,ang,cen):
	"""
	>>> x = rotation_around(Ux,1,Uy)
	>>> x * Uy
	Vec( 0.000000, 1.000000, 0.000000 )
	"""
	R = rotation_matrix(axs,ang)
	return Xform(R,R*-cen+cen)

def rotation_around_degrees(axs,ang,cen): return rotation_around(axs,radians(ang),cen)

def test():
	test_rotation_mat()


def dihedral(p1,p2,p3,p4):
	"""
	>>> dihedral_degrees(Ux,Uy,V0,Uz)
	90.0
	>>> dihedral_degrees(Ux,V0,Uy,Uz)
	-90.0
	"""
	a = ( p2 - p1 ).normalized()
	b = ( p3 - p2 ).normalized()
	c = ( p4 - p3 ).normalized()
	x = -a.dot(c) + a.dot(b) * b.dot(c)
	y =  a.dot( b.cross(c) );
	return atan2(y,x)

def dihedral_degrees(p1,p2,p3,p4): return degrees(dihedral(p1,p2,p3,p4))

def angle(p1,p2,p3=None):
	if p3 is None: return acos( p1.normalized().dot(p2.normalized()) )
	else:
			a = ( p2 - p1 ).normalized()
			b = ( p2 - p3 ).normalized()
			return acos( a.dot(b) )

def angle_degrees(p1,p2,p3=None):
	if p3 is None: return 180.0/math.pi*acos( p1.normalized().dot(p2.normalized()) )
	else:
			a = ( p2 - p1 ).normalized()
			b = ( p2 - p3 ).normalized()
			return 180.0/math.pi*acos( a.dot(b) )

def sin_cos_range(x):
	assert -1.001 < x < 1.001
	return min(1.0,max(-1.0,x))

def point_line_distance(p,a,c):
	"""
	>>> point_line_distance(V0,Uy,V0)
	0.0
	>>> round(point_line_distance(V0,Uy,Ux+Uz),8)
	1.41421356
	>>> round(point_line_distance(Ux,Ux,Vec(3,2,1)) , 8)
	2.23606798
	"""
	return projperp(a,p-c).length()

def line_line_angle(a1,a2):
	a = a1.angle(a2)
	return min(a,math.pi-a)

def line_line_angle_degrees(a1,a2):
	return line_line_angle(a1,a2)*180/math.pi

def line_line_distance(a1,c1,a2,c2):
	"""
	>>> line_line_distance(Ux,V0,Uy,V0)
	0.0
	>>> round(line_line_distance(Ux,Vec(0,1,2),Ux,Vec(3,2,1)) , 8)
	1.41421356
	>>> line_line_distance(Ux,10*Uy,Uz,99.0*Ux)
	10.0
	>>> X = randxform()
	>>> round(line_line_distance(X.R*Ux,X*Vec(0,1,2),X.R*Ux,X*Vec(3,2,1)) , 8)
	1.41421356
	"""
	a1 = a1.normalized()
	a2 = a2.normalized()
	if abs(a1.dot(a2)) > 0.9999: return projperp(a1,c2-c1).length()
	a = a1
	b = a2
	c = c2-c1
	n = abs(c.dot(a.cross(b)))
	d = a.cross(b).length()
	if abs(d) < EPS: return 0
	return n/d

def line_line_closest_points(A1,C1,A2,C2):
	"""
	>>> print line_line_closest_points(Ux,Ux,Uy,Uy)
	(Vec( 0.000000, 0.000000, 0.000000 ), Vec( 0.000000, 0.000000, 0.000000 ))

	>>> print line_line_closest_points(Ux,Uy,Uy,Uz)
	(Vec( 0.000000, 1.000000, 0.000000 ), Vec( 0.000000, 1.000000, 1.000000 ))
	"""
	# Line1 : C1 + t1*A1
	# Line2 : C2 + t2*A2
	# (Ucase: Vectors , Lcase: Scalars)
	# To find the points with minimun distance is equivalent to find Q1(t1)
	# & Q2(t2) /
	# Q2 - Q1 = C2+t2*A2 - C1-t1*A1 = k*(A2 x A1)
	# ( x mean Cross product)
	# Using some tricks and vector properties the solution is:
	# print type(A1)
	# print type(A2)
	# print A2.cross(A1)
	C21 = C2 - C1
	M = A2.cross(A1)
	m2 = M.dot(M)
	R = C21.cross(M/m2)
	t1 = R.dot(A2)
	t2 = R.dot(A1)
	Q1 = C1 + t1 * A1
	Q2 = C2 + t2 * A2
	return Q1,Q2

def skew_lines_center(A1,C1,A2,C2):
	"""
	>>> skew_lines_center(Ux,V0,Uy,V0)
	Vec( 0.000000, 0.000000, 0.000000 )

	>>> skew_lines_center(Ux,Uy,Uy,Ux)
	Vec( 1.000000, 1.000000, 0.000000 )

	>>> skew_lines_center(10*Ux,10*Uy,10*Uy,10*Uz)
	Vec( 0.000000, 10.000000, 5.000000 )

	# >>> skew_lines_center(Ux,Uy,Uz,Ux)
	"""
	p1,p2 = line_line_closest_points(A1,C1,A2,C2)
	return (p1+p2)/2.0

def skew_lines_relation(a1,c1,a2,c2):
	p1,p2 = line_line_closest_points(a1,c1,a2,c2)
	return dihedral( p1+a1, p1, p2, p2+a2 )

def skew_lines_relation_z(axis,cen):
	return skew_lines_relation(axis,cen,Vec(0,0,1),Vec(0,0,0))

def align_skew_line_pairs(aa1,ac1,aa2,ac2,ba1,bc1,ba2,bc2):
	"""
	>>> aa1 = Ux
	>>> ac1 = V0
	>>> aa2 = Uy
	>>> ac2 = V0
	>>> ba1 = Ux
	>>> bc1 = V0
	>>> ba2 = Uy
	>>> bc2 = V0
	>>> print align_skew_line_pairs(aa1,ac1,aa2,ac2,ba1,bc1,ba2,bc2)
	Xform( Mat[ (1.000000,-0.000000,-0.000000), (-0.000000,1.000000,0.000000), (0.000000,-0.000000,1.000000) ], (0.000000,0.000000,0.000000) )

	>>> aa1 = Uy
	>>> ac1 = V0
	>>> aa2 = Uz
	>>> ac2 = V0
	>>> ba1 = Ux
	>>> bc1 = Uz
	>>> ba2 = Uy
	>>> bc2 = Uz
	>>> print align_skew_line_pairs(aa1,ac1,aa2,ac2,ba1,bc1,ba2,bc2)
	Xform( Mat[ (0.000000,1.000000,0.000000), (0.000000,0.000000,1.000000), (1.000000,0.000000,0.000000) ], (0.000000,0.000000,1.000000) )
	"""
	acen = skew_lines_center(aa1,ac1,aa2,ac2)
	bcen = skew_lines_center(ba1,bc1,ba2,bc2)
	R = alignvectors_minangle(aa1,aa2,ba1,ba2).R
	return Xform(R,bcen-acen)

def  sindeg(x): return sin(math.pi/180.0*x)
def  cosdeg(x): return cos(math.pi/180.0*x)
def asindeg(x): return 180.0/math.pi*asin(x)
def acosdeg(x): return 180.0/math.pi*acos(x)

def rotation_around_dof_to_set_vec_vec_angle(dofaxis,tgt0,v1,v2):
	"""
	>>> from random import random
	>>> for i in range(10):
	... 	dof = randnorm()
	... 	tgt = random()*180
	... 	v1 = randnorm()
	... 	v2 = randnorm()
	... 	ANGs = rotation_around_dof_to_set_vec_vec_angle(dof,tgt,v1,v2)
	... 	for a in ANGs:
	... 		R = rotation_around_degrees(dof,a,Vec(0,0,0))
	... 		act = line_line_angle_degrees(v1,R*v2)
	... 		if 0.000001 < min(abs(act-tgt),abs(act-180+tgt)):
	... 			print a,tgt,act
	>>> print "if no other output, correct"
	if no other output, correct

	"""
	tgt0 = tgt0%360
	tgt0 = min(tgt0,360-tgt0)
	tgt0 = min(tgt0,180-tgt0)
	a1 = angle_degrees(dofaxis,v1)
	a2 = angle_degrees(dofaxis,v2)
	if a1 > 90.0: a1,v1 = 180.0-a1,-v1
	if a2 > 90.0: a2,v2 = 180.0-a2,-v2
	R = rotation_around_degrees((Vec(0,0,1)+dofaxis)/2.0,180.0,Vec(0,0,0))
	dofaxis = R*dofaxis
	v1 = R*v1
	v2 = R*v2
	ret0 = []
	for tgt in (tgt0,180-tgt0):
		if (-tgt <= a2+a1 <= +tgt) == (-tgt <= a2-a1 <= +tgt): continue
		sc1 = cosdeg(a1-tgt),sindeg(a1-tgt)
		sc2 = cosdeg(a1+tgt),sindeg(a1+tgt)
		# print a1-tgt
		# print a1+tgt
		# print "seg:",sc1
		# print "seg:",sc2
		m = (sc1[1]-sc2[1])/(sc1[0]-sc2[0])
		b = sc1[1]-m*sc1[0]
		# print "mx+b",m,b
		d = cos(a2/180*pi)
		h0 = m*d+b
		h = h0/sindeg(a2)
		# test = (-tgt <= a2+a1 <= +tgt) != (-tgt <= a2-a1 <= +tgt)
		# assert test == (-1 <= h <= 1)
		# if test: continue
		# print "h",h
		# if h < -1.0 or h > 1.0: continue
		xa = 90-asindeg(h)
		# print asindeg(h)
		di = -dihedral_degrees(v1,Vec(0,0,0),dofaxis,v2)
		# print -di
		ret0.append((di+xa)%360.0)
		ret0.append((di-xa)%360.0)
	if not ret0: return []
	ret0.sort()
	ret = [ret0[0]]
	for i in range(1,len(ret0)):
		if abs(ret0[i-1]-ret0[i]) > 1.0:
			ret.append(ret0[i])
	return ret


def ray_sphere_intersection(lin,l0in,cen,r):
	"""
	>>> v = randnorm()
	>>> v = Ux
	>>> assert v.distance( ray_sphere_intersection(v,V0,V0,1) ) < 0.00001
	>>> assert not ray_sphere_intersection(v,V0,-2*v,1.0)
	"""
	l = lin.normalized()
	l0 = l0in - cen
 	a = l.dot(l)
 	b = l.dot(l0)*2.0
 	c = l0.dot(l0) - r*r
 	disc = b * b - 4 * a * c;
 	if disc < 0: return None
 	distSqrt = sqrt(disc)
 	if b < 0: q = (-b - distSqrt)/2.0
 	else:     q = (-b + distSqrt)/2.0
 	t0 = q / a;
 	t1 = c / q;
 	if t0 > t1: t0,t1 = t1,t0
 	if t1 >= 0:
	 	if t0 >= 0: return l*t0+l0+cen
 		else:       return l*t1+l0+cen
 	return None

def line_plane_intersection(l,l0,n,p0):
	"""
	>>> l  = Ux
	>>> l0 = randvec()
	>>> n  = Ux
	>>> p0 = V0
	>>> assert line_plane_intersection(l,l0,n,p0)[1] == Vec(0,l0.y,l0.z)
	>>> n = randnorm()
	>>> p0 = randvec().cross(n)
	>>> l = randvec()
	>>> l0 = p0+l*gauss(0,10)
	>>> assert line_plane_intersection(l,l0,n,p0)[1] == p0
	"""
	n = n.normalized()
	d = (p0-l0).dot(n) / l.dot(n)
	return d,d*l+l0

def slide_to_make_lines_intersect(dof,l,l0,m,m0):
	"""
	>>> v = randvec()
	>>> assert abs(slide_to_make_lines_intersect(Ux,Uy,v,Uz,V0) + v.x ) < EPS
	>>> dof,l,l0,m,m0 = randvec(5)
	>>> d = slide_to_make_lines_intersect(dof,l,l0,m,m0)
	>>> l0 = l0 + d*dof
	>>> assert abs(line_line_distance(l,l0,m,m0)) < EPS
	"""
	n  = l.cross(m)
	p0 = m0
	d,i = line_plane_intersection(dof,l0,n,p0)
	assert ( (i-l0).normalized().dot(dof.normalized()) - 1.0 ) < EPS
	assert i-l0 == dof*d
	return d

def alignvector(a,b):
	"""
	>>> u = randvec()
	>>> v = randvec()
	>>> assert v.angle(alignvector(u,v)*u) < EPS
	"""
	return rotation_around(a.normalized()+b.normalized(),pi,V0)

def alignaroundaxis(axis,u,v):
	"""
	>>> axis = randnorm()
	>>> u = randvec()
	>>> angle = uniform(-pi,pi)
	>>> v = rotation_matrix(axis,angle)*u
	>>> uprime = alignaroundaxis(axis,u,v)*u
	>>> assert v.angle(uprime) < EPS
	>>> v = randvec()
	>>> uprime = alignaroundaxis(axis,u,v)*u
	>>> assert coplanar(V0,axis,v,uprime)
	"""
	return rotation_around(axis, -dihedral(u,axis,V0,v),V0 )

def alignvectors_minangle(a1,a2,b1,b2):
	"""
	exact alignment:

	>>> for i in range(10):
	... 	angdeg = uniform(-180,180)
	... 	a1 = randvec()
	... 	b1 = randnorm()*a1.length()
	... 	l2 = gauss(0,1)
	... 	a2 = rotation_matrix_degrees(a1.cross(randnorm()),angdeg) * a1 * l2
	... 	b2 = rotation_matrix_degrees(b1.cross(randnorm()),angdeg) * b1 * l2
	... 	assert abs(angle(a1,a2) - angle(b1,b2)) < EPS
	... 	Xa2b = alignvectors_minangle(a1,a2,b1,b2)
	... 	assert Xa2b.t.length() < EPS
	... 	assert (Xa2b*a1).distance(b1) < EPS
	... 	assert (Xa2b*a2).distance(b2) < EPS

	if angle(a1,a2) != angle(b1,2b), minimize deviation

	>>> a1,a2,b1,b2 = randvec(4)
	>>> Xa2b = alignvectors_minangle(a1,a2,b1,b2)
	>>> assert coplanar(b1,b2,Xa2b*a1,Xa2b*a2)
	>>> assert (b1.angle(a1)+b2.angle(a2)) > (b1.angle(Xa2b*a1)+b2.angle(Xa2b*a2))

	# >>> tgt1 = -Vec(0.816497,0.000000,0.577350)
	# >>> tgt2 = Vec(0.000000,0.000000,1.000000)
	# >>> orig1 = Vec(0.000000,0.000000,1.000000)
	# >>> orig2 = Vec(-0.723746,0.377967,-0.577350)
	# >>> print orig1.angle_degrees(orig2)
	# >>> print tgt1.angle_degrees(tgt2)
	# >>> x = alignvectors_minangle(orig1,orig2,tgt1,tgt2)
	# >>> print tgt1,x*orig1
	# >>> print tgt2,x*orig2
	"""
	aaxis = (a1.normalized()+a2.normalized())/2.0
	baxis = (b1.normalized()+b2.normalized())/2.0
	Xmiddle = alignvector(aaxis,baxis)
	assert (baxis).angle(Xmiddle*(aaxis)) < SQRTEPS
	Xaround = alignaroundaxis(baxis, Xmiddle*a1, b1 )#
	X = Xaround * Xmiddle
	assert (b1.angle(a1)+b2.angle(a2)) >= (b1.angle(X*a1)+b2.angle(X*a2))
	return X
	# not so good if angles don't match:
	# xa = Xform().from_two_vecs(a2,a1)
	# xb = Xform().from_two_vecs(b2,b1)
	# return xb/xa

def alignvectors(a1,a2,b1,b2):
	"same as alignvectors_minangle"
	return alignvectors_minangle(a1,a2,b1,b2)

# def alignvectors_kindamindis(a1,a2,b1,b2):
#    """
#    >>> ang = uniform(0,1)*360.0-180.0
#    >>> a1 = randvec()
#    >>> b1 = randnorm()*a1.length()
#    >>> l2 = gauss(0,1)
#    >>> a2 = rotation_matrix_degrees(a1.cross(randnorm()),ang) * a1 * l2
#    >>> b2 = rotation_matrix_degrees(b1.cross(randnorm()),ang) * b1 * l2
#    >>> assert abs(angle(a1,a2) - angle(b1,b2)) < EPS
#    >>> Xa2b = alignvectors(a1,a2,b1,b2)
#    >>> assert Xa2b.t.length() < EPS
#    >>> assert (Xa2b*a1).distance(b1) < EPS
#    >>> assert (Xa2b*a2).distance(b2) < EPS

#    >>> a1 = randvec()
#    >>> b1 = randvec()
#    >>> a2 = randvec()
#    >>> b2 = randvec()
#    >>> Xa2b = alignvectors(a1,a2,b1,b2)
#    >>> assert coplanar(b1,b2,Xa2b*a1,Xa2b*a2)
#    >>> if not (b1.distance(a1)+b2.distance(a2)) > (b1.distance(Xa2b*a1)+b2.distance(Xa2b*a2)):
#    ...   print b1
#    ...   print b2
#    ...   print a1
#    ...   print a2
#    ...   print Xa2b*a1
#    ...   print Xa2b*a2
#    """
#    Xmiddle = alignvector(a1+a2,b1+b2)
#    assert (b1+b2).angle(Xmiddle*(a1+a2)) < SQRTEPS
#    assert (b1+b2).angle(Xmiddle*a1+Xmiddle*a2) < SQRTEPS
#    Xaround = alignaroundaxis(b1+b2, Xmiddle*a1, b1 )#
#    return Xaround * Xmiddle
#    # xa = Xform().from_two_vecs(a2,a1)
#    # xb = Xform().from_two_vecs(b2,b1)
#    return xb/xa

def get_test_generators1():
	x1 = rotation_around_degrees(Vec(0,0,1),180,Vec(0,0,0))
	x2 = rotation_around_degrees(Vec(1,1,1),120,Vec(1,0,0))
	return x1,x2

def expand_xforms(G,N=3,c=Vec(1,3,10)):
	"""
	>>> G = get_test_generators1()
	>>> for x in expand_xforms(G): print x*Ux
	(-1.000000,0.000000,0.000000)
	(1.000000,0.000000,0.000000)
	(1.000000,-0.000000,0.000000)
	(-1.000000,0.000000,0.000000)
	(1.000000,-2.000000,0.000000)
	(1.000000,0.000000,0.000000)
	(-1.000000,2.000000,0.000000)
	(-1.000000,0.000000,0.000000)
	(1.000000,-2.000000,0.000000)
	(1.000000,-0.000000,-2.000000)
	"""
	seenit = set()
	for Xs in chain(G,*(product(G,repeat=n) for n in range(2,N+1))):
		X = Xs if isinstance(Xs,Xform) else reduce(Xform.__mul__,Xs)
		v = X*c
		key = (round(v.x,3),round(v.y,3),round(v.z,3))
		if key not in seenit:
			seenit.add(key)
			yield X

def find_identities(G,n=6,c=Vec(1,3,10)):
	"""
	>>> G = get_test_generators1()
	>>> for I in find_identities(G): print I.t
	(0.000000,0.000000,0.000000)
	(-2.000000,2.000000,2.000000)
	(2.000000,-2.000000,2.000000)
	"""
	for x in expand_xforms(G,n,c):
		if (abs(x.R.xx-1.0) < 0.0000001 and
			 abs(x.R.yy-1.0) < 0.0000001 and
			 abs(x.R.zz-1.0) < 0.0000001 ):
			yield x

def get_cell_bounds_orthogonal_only(G,n=6,c=Vec(1,3,10)):
	"""
	very slow... need to speed up

	>>> G = get_test_generators1()                  #doctest: +SKIP
	>>> get_cell_bounds_orthogonal_only(G[:2],12)   #doctest: +SKIP
	(4.0, 4.0, 4.0)
	"""
	mnx,mny,mnz = 9e9,9e9,9e9
	for i in (I.t for I in find_identities(G,n)):
		if abs(i.x) > SQRTEPS and abs(i.y) < SQRTEPS and abs(i.z) < SQRTEPS: mnx = min(mnx,abs(i.x))
		if abs(i.x) < SQRTEPS and abs(i.y) > SQRTEPS and abs(i.z) < SQRTEPS: mny = min(mny,abs(i.y))
		if abs(i.x) < SQRTEPS and abs(i.y) < SQRTEPS and abs(i.z) > SQRTEPS: mnz = min(mnz,abs(i.z))
	return round(mnx,3),round(mny,3),round(mnz,3)




AXIS_I2 = Vec(+1.00000000000000,+0.00000000000000,+0.00000000000000 ).normalized()
AXIS_I3 = Vec(+0.93417235896272,+0.00000000000000,+0.35682208977309 ).normalized()
AXIS_I5 = Vec(+0.85065080835204,+0.52573111211914,-0.00000000000000 ).normalized()
AXIS_O2 = Vec( 1.00000000000000, 1.00000000000000, 0.00000000000000 ).normalized()
AXIS_O3 = Vec( 1.00000000000000, 1.00000000000000, 1.00000000000000 ).normalized()
AXIS_O4 = Vec( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized()
AXIS_T2 = Vec( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized()
AXIS_T3 = Vec( 1.00000000000000, 1.00000000000000, 1.00000000000000 ).normalized()
ROT_I2 = rotation_around(AXIS_I2,math.pi*2/2,V0)
ROT_I3 = rotation_around(AXIS_I3,math.pi*2/3,V0)
ROT_I5 = rotation_around(AXIS_I5,math.pi*2/5,V0)
ROT_O2 = rotation_around(AXIS_O2,math.pi*2/2,V0)
ROT_O3 = rotation_around(AXIS_O3,math.pi*2/3,V0)
ROT_O4 = rotation_around(AXIS_O4,math.pi*2/4,V0)
ROT_T2 = rotation_around(AXIS_T2,math.pi*2/2,V0)
ROT_T3 = rotation_around(AXIS_T3,math.pi*2/3,V0)

SYMAXIS = dict()
SYMCEN  = dict()
SYMAXIS["ICS",5] = AXIS_I5
SYMCEN ["ICS",5] = Vec(0,0,0)
SYMAXIS["ICS",3] = AXIS_I3
SYMCEN ["ICS",3] = Vec(0,0,0)
SYMAXIS["ICS",2] = AXIS_I2
SYMCEN ["ICS",2] = Vec(0,0,0)
SYMAXIS["OCT",4] = AXIS_O4
SYMCEN ["OCT",4] = Vec(0,0,0)
SYMAXIS["OCT",3] = AXIS_O3
SYMCEN ["OCT",3] = Vec(0,0,0)
SYMAXIS["OCT",2] = AXIS_O2
SYMCEN ["OCT",2] = Vec(0,0,0)
SYMAXIS["TET",3] = AXIS_T3
SYMCEN ["TET",3] = Vec(0,0,0)
SYMAXIS["TET",2] = AXIS_T2
SYMCEN ["TET",2] = Vec(0,0,0)
SYMAXIS["P6",2] = Uz
SYMAXIS["P6",3] = Uz
SYMAXIS["P6",6] = Uz
SYMCEN ["P6",2] = Vec(0.50000000000,0.000000000,0.000000000)
SYMCEN ["P6",3] = Vec(0.50000000000,math.atan(math.pi/6.0)/2,0.000000000)
SYMCEN ["P6",6] = Vec(0.00000000000,0.000000000,0.000000000)

# for i in range(10): print "FIXME SYM FRAME GEN"
# SYMICS = tuple(expand_xforms((ROT_I2,ROT_I3),12))
# SYMOCT = tuple(expand_xforms((ROT_O2,ROT_O3),7))
# SYMTET = tuple(expand_xforms((ROT_T2,ROT_T3),7))

SYMTET = [
	Xform( Mat( Vec( 1.000000, 0.000000, 0.000000 ), Vec( 0.000000, -1.000000, 0.000000 ), Vec( 0.000000, -0.000000, -1.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ),
	Xform( Mat( Vec( 0.000000, 1.000000, -0.000000 ), Vec( -0.000000, 0.000000, 1.000000 ), Vec( 1.000000, -0.000000, 0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ),
	Xform( Mat( Vec( 1.000000, 0.000000, 0.000000 ), Vec( 0.000000, 1.000000, -0.000000 ), Vec( 0.000000, 0.000000, 1.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ),
	Xform( Mat( Vec( 0.000000, -1.000000, 0.000000 ), Vec( -0.000000, -0.000000, -1.000000 ), Vec( 1.000000, 0.000000, -0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ),
	Xform( Mat( Vec( 0.000000, 1.000000, -0.000000 ), Vec( 0.000000, -0.000000, -1.000000 ), Vec( -1.000000, 0.000000, -0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ),
	Xform( Mat( Vec( -0.000000, 0.000000, 1.000000 ), Vec( 1.000000, -0.000000, 0.000000 ), Vec( 0.000000, 1.000000, -0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ),
	Xform( Mat( Vec( 0.000000, -1.000000, 0.000000 ), Vec( 0.000000, 0.000000, 1.000000 ), Vec( -1.000000, -0.000000, 0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ),
	Xform( Mat( Vec( -0.000000, -0.000000, -1.000000 ), Vec( 1.000000, 0.000000, -0.000000 ), Vec( 0.000000, -1.000000, 0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ),
	Xform( Mat( Vec( 0.000000, -0.000000, -1.000000 ), Vec( -1.000000, -0.000000, -0.000000 ), Vec( -0.000000, 1.000000, -0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ),
	Xform( Mat( Vec( -0.000000, 0.000000, 1.000000 ), Vec( -1.000000, 0.000000, -0.000000 ), Vec( -0.000000, -1.000000, 0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ),
	Xform( Mat( Vec( -1.000000, -0.000000, -0.000000 ), Vec( -0.000000, 1.000000, 0.000000 ), Vec( 0.000000, 0.000000, -1.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ),
	Xform( Mat( Vec( -1.000000, 0.000000, -0.000000 ), Vec( -0.000000, -1.000000, 0.000000 ), Vec( -0.000000, 0.000000, 1.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) )
]
SYMOCT = [
	Xform( Mat( Vec(0.000000,1.000000,0.000000), Vec(1.000000,0.000000,-0.000000), Vec(-0.000000,0.000000,-1.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,-0.000000,1.000000), Vec(1.000000,0.000000,-0.000000), Vec(-0.000000,1.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(1.000000,0.000000,-0.000000), Vec(0.000000,1.000000,0.000000), Vec(0.000000,-0.000000,1.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(1.000000,0.000000,-0.000000), Vec(0.000000,-0.000000,1.000000), Vec(0.000000,-1.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.000000,0.000000,-1.000000), Vec(0.000000,1.000000,0.000000), Vec(1.000000,-0.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.000000,1.000000,0.000000), Vec(0.000000,-0.000000,1.000000), Vec(1.000000,0.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,1.000000,0.000000), Vec(-0.000000,0.000000,-1.000000), Vec(-1.000000,0.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,-0.000000,1.000000), Vec(-0.000000,1.000000,0.000000), Vec(-1.000000,-0.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,-1.000000,-0.000000), Vec(1.000000,0.000000,0.000000), Vec(0.000000,-0.000000,1.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(1.000000,-0.000000,-0.000000), Vec(-0.000000,0.000000,-1.000000), Vec(0.000000,1.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(1.000000,0.000000,0.000000), Vec(0.000000,-1.000000,-0.000000), Vec(-0.000000,0.000000,-1.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.000000,0.000000,-1.000000), Vec(1.000000,-0.000000,-0.000000), Vec(-0.000000,-1.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-1.000000,0.000000,0.000000), Vec(0.000000,1.000000,-0.000000), Vec(-0.000000,0.000000,-1.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-1.000000,-0.000000,0.000000), Vec(0.000000,0.000000,1.000000), Vec(-0.000000,1.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,-0.000000,1.000000), Vec(0.000000,-1.000000,-0.000000), Vec(1.000000,0.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,1.000000,-0.000000), Vec(-1.000000,0.000000,0.000000), Vec(0.000000,-0.000000,1.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,0.000000,1.000000), Vec(-1.000000,-0.000000,0.000000), Vec(0.000000,-1.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,-1.000000,-0.000000), Vec(-0.000000,-0.000000,1.000000), Vec(-1.000000,-0.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.000000,-1.000000,-0.000000), Vec(-0.000000,0.000000,-1.000000), Vec(1.000000,-0.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.000000,0.000000,-1.000000), Vec(-1.000000,0.000000,0.000000), Vec(0.000000,1.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,0.000000,-1.000000), Vec(-0.000000,-1.000000,-0.000000), Vec(-1.000000,0.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-1.000000,0.000000,0.000000), Vec(-0.000000,0.000000,-1.000000), Vec(-0.000000,-1.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-1.000000,-0.000000,-0.000000), Vec(0.000000,-1.000000,-0.000000), Vec(-0.000000,-0.000000,1.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,-1.000000,-0.000000), Vec(-1.000000,-0.000000,-0.000000), Vec(0.000000,0.000000,-1.000000) ), Vec(0.000000,0.000000,0.000000) ),
]
SYMICS = [
	Xform( Mat( Vec(1.000000,0.000000,0.000000), Vec(0.000000,-1.000000,-0.000000), Vec(0.000000,0.000000,-1.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.809017,-0.309017,0.500000), Vec(0.309017,-0.500000,-0.809017), Vec(0.500000,0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.809017,0.309017,0.500000), Vec(-0.309017,-0.500000,0.809017), Vec(0.500000,-0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(1.000000,0.000000,0.000000), Vec(0.000000,1.000000,0.000000), Vec(0.000000,-0.000000,1.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.809017,-0.309017,0.500000), Vec(-0.309017,0.500000,0.809017), Vec(-0.500000,-0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.809017,0.309017,-0.500000), Vec(0.309017,0.500000,0.809017), Vec(0.500000,-0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.809017,0.309017,-0.500000), Vec(-0.309017,-0.500000,-0.809017), Vec(-0.500000,0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.809017,0.309017,0.500000), Vec(0.309017,0.500000,-0.809017), Vec(-0.500000,0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.500000,-0.809017,0.309017), Vec(0.809017,0.309017,-0.500000), Vec(0.309017,0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.809017,-0.309017,-0.500000), Vec(-0.309017,0.500000,-0.809017), Vec(0.500000,0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.500000,-0.809017,0.309017), Vec(-0.809017,-0.309017,0.500000), Vec(-0.309017,-0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.809017,-0.309017,-0.500000), Vec(0.309017,-0.500000,0.809017), Vec(-0.500000,-0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.500000,0.809017,-0.309017), Vec(0.809017,-0.309017,0.500000), Vec(0.309017,-0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.309017,0.500000,0.809017), Vec(0.500000,-0.809017,0.309017), Vec(0.809017,0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.309017,-0.500000,0.809017), Vec(-0.500000,-0.809017,-0.309017), Vec(0.809017,-0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.500000,0.809017,-0.309017), Vec(-0.809017,0.309017,-0.500000), Vec(-0.309017,0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.309017,0.500000,0.809017), Vec(-0.500000,0.809017,-0.309017), Vec(-0.809017,-0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.309017,-0.500000,0.809017), Vec(0.500000,0.809017,0.309017), Vec(-0.809017,0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.500000,-0.809017,-0.309017), Vec(0.809017,0.309017,0.500000), Vec(-0.309017,-0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.309017,-0.500000,-0.809017), Vec(0.500000,0.809017,-0.309017), Vec(0.809017,-0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.309017,0.500000,-0.809017), Vec(-0.500000,0.809017,0.309017), Vec(0.809017,0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.500000,0.809017,0.309017), Vec(-0.809017,0.309017,0.500000), Vec(0.309017,-0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.500000,-0.809017,-0.309017), Vec(-0.809017,-0.309017,-0.500000), Vec(0.309017,0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.309017,-0.500000,-0.809017), Vec(-0.500000,-0.809017,0.309017), Vec(-0.809017,0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.309017,0.500000,-0.809017), Vec(0.500000,-0.809017,-0.309017), Vec(-0.809017,-0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.500000,0.809017,0.309017), Vec(0.809017,-0.309017,-0.500000), Vec(-0.309017,0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,0.000000,1.000000), Vec(1.000000,-0.000000,-0.000000), Vec(0.000000,1.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.309017,-0.500000,0.809017), Vec(0.500000,-0.809017,-0.309017), Vec(0.809017,0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,-1.000000,0.000000), Vec(-0.000000,-0.000000,-1.000000), Vec(1.000000,0.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,0.000000,1.000000), Vec(-1.000000,0.000000,0.000000), Vec(-0.000000,-1.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.309017,-0.500000,0.809017), Vec(-0.500000,0.809017,0.309017), Vec(-0.809017,-0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,-1.000000,0.000000), Vec(0.000000,0.000000,1.000000), Vec(-1.000000,-0.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,-0.000000,-1.000000), Vec(1.000000,0.000000,0.000000), Vec(0.000000,-1.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.309017,0.500000,-0.809017), Vec(0.500000,0.809017,0.309017), Vec(0.809017,-0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,1.000000,0.000000), Vec(0.000000,-0.000000,1.000000), Vec(1.000000,-0.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.309017,0.500000,0.809017), Vec(-0.500000,-0.809017,0.309017), Vec(0.809017,-0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,-0.000000,-1.000000), Vec(-1.000000,-0.000000,-0.000000), Vec(-0.000000,1.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.309017,0.500000,-0.809017), Vec(-0.500000,-0.809017,-0.309017), Vec(-0.809017,0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(0.000000,1.000000,0.000000), Vec(-0.000000,0.000000,-1.000000), Vec(-1.000000,0.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.309017,0.500000,0.809017), Vec(0.500000,0.809017,-0.309017), Vec(-0.809017,0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.500000,-0.809017,0.309017), Vec(0.809017,-0.309017,0.500000), Vec(-0.309017,0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.500000,-0.809017,-0.309017), Vec(0.809017,-0.309017,-0.500000), Vec(0.309017,-0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.309017,-0.500000,-0.809017), Vec(-0.500000,0.809017,-0.309017), Vec(0.809017,0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.500000,-0.809017,0.309017), Vec(-0.809017,0.309017,-0.500000), Vec(0.309017,-0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.500000,-0.809017,-0.309017), Vec(-0.809017,0.309017,0.500000), Vec(-0.309017,0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.309017,-0.500000,-0.809017), Vec(0.500000,-0.809017,0.309017), Vec(-0.809017,-0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.500000,0.809017,-0.309017), Vec(0.809017,0.309017,-0.500000), Vec(-0.309017,-0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.500000,0.809017,0.309017), Vec(0.809017,0.309017,0.500000), Vec(0.309017,0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.809017,0.309017,0.500000), Vec(0.309017,-0.500000,0.809017), Vec(0.500000,0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.809017,-0.309017,0.500000), Vec(-0.309017,-0.500000,-0.809017), Vec(0.500000,-0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.500000,0.809017,-0.309017), Vec(-0.809017,-0.309017,0.500000), Vec(0.309017,0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.500000,0.809017,0.309017), Vec(-0.809017,-0.309017,-0.500000), Vec(-0.309017,-0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.809017,0.309017,0.500000), Vec(-0.309017,0.500000,-0.809017), Vec(-0.500000,-0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.809017,-0.309017,0.500000), Vec(0.309017,0.500000,0.809017), Vec(-0.500000,0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.809017,-0.309017,-0.500000), Vec(0.309017,0.500000,-0.809017), Vec(0.500000,-0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.809017,0.309017,-0.500000), Vec(-0.309017,0.500000,0.809017), Vec(0.500000,0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.809017,-0.309017,-0.500000), Vec(-0.309017,-0.500000,0.809017), Vec(-0.500000,0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-0.809017,0.309017,-0.500000), Vec(0.309017,-0.500000,-0.809017), Vec(-0.500000,-0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-1.000000,-0.000000,0.000000), Vec(0.000000,-1.000000,-0.000000), Vec(0.000000,0.000000,1.000000) ), Vec(0.000000,0.000000,0.000000) ),
	Xform( Mat( Vec(-1.000000,-0.000000,0.000000), Vec(-0.000000,1.000000,0.000000), Vec(-0.000000,-0.000000,-1.000000) ), Vec(0.000000,0.000000,0.000000) ),
]

def cyclic_axis(coords):
	if max(len(x) for x in coords) != min(len(x) for x in coords):
		print "coord vectors not same length",
		return None
	cen = reduce(op.add,(x for c in coords for x in c),V0) / len(coords[0]*len(coords))
	guesses = []
	axis = Vec(0,0,0)
	diserr = 0
	rand = randvec()
	for xyzs in izip(*coords):
		a = reduce(op.add,(xyz-cen for xyz in xyzs))/len(xyzs)
		da = a if a.dot(axis) > 0 else -a
		axis += da
		guesses.append(da)
		dis = [a.distance(cen) for a in xyzs]
		avgdis = sum(dis)/len(dis)
		diserr += sum((x-avgdis)**2 for x in dis)
	axis.normalize()
	diserr = math.sqrt(diserr/len(coords)/len(coords[0]))
	angerr = math.sqrt(sum(projperp(axis,a).length_squared() for a in guesses)/len(guesses))
	return axis,cen,diserr,angerr

def symmetrize_xform(anchor,x,nf=None): #,debug=False):
	"""
	# >>> x = rotation_around_degrees(Uz,180,V0)
	# >>> assert symmetrize_xform(Ux,x,2) == x
	#>>> x = rotation_around_degrees(Vec(0,0,1),121,Vec(1,1,1))
	#>>> x,c = symmetrize_xform(Ux,x,3); print x.pretty(); print c
	# >>> print x.pretty()
	"""
	if not nf:
		axs,ang	= x.rotation_axis()
		nf    = int(round(2*pi/ang))
	sang = 360.0/nf
	axs,ang = x.rotation_axis()
	p1 = anchor
	p2 = x*p1
	axs = projperp(p2-p1,axs)
	cendir = axs.cross(p2-p1).normalized()
	cenlen = (p2-p1).length()/2/math.tan(sang/2/180*pi)
	# print "p1",p1
	# print "p2",p2
	# print "axs",axs
	# print "p2-p1",p2-p1
	# print "(p2-p1)/2",(p2-p1)/2
	# print "cendir",cendir,   cendir.angle_degrees(p2-p1)
	# print "cenlen",cenlen
	# print "math.tan(sang/2/180*pi)",math.tan(sang/2/180*pi)
	# print "(p2-p1).length()",(p2-p1).length()
	cen = (p1+p2)/2 + cendir*cenlen
	# print cen
	return axs,cen,nf
















if __name__ == '__main__':
	# xforms = read_xforms("/Users/sheffler/project/symgen_movie/fere.pdb_xforms.dat")
	# print len(xforms)
	# print xforms[-1]

	import doctest
	r = doctest.testmod()
	print r
