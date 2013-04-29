struct xyzmatrix
{
  float xx_, xy_, xz_;
  float yx_, yy_, yz_;
  float zx_, zy_, zz_;
};

struct xyzvector
{
  float x_, y_, z_;
};

/// Modify matrix 1, m1, by performing a right-multiply into it of m2
void
right_multiply_by_private(
  __private struct xyzmatrix * m1,
  __private struct xyzmatrix * m2
)
{
  float x, y, z; // Temporaries

  // First row
  x = ( m1->xx_ * m2->xx_ ) + ( m1->xy_ * m2->yx_ ) + ( m1->xz_ * m2->zx_ );
  y = ( m1->xx_ * m2->xy_ ) + ( m1->xy_ * m2->yy_ ) + ( m1->xz_ * m2->zy_ );
  z = ( m1->xx_ * m2->xz_ ) + ( m1->xy_ * m2->yz_ ) + ( m1->xz_ * m2->zz_ );
  m1->xx_ = x; m1->xy_ = y; m1->xz_ = z;

  // Second row
  x = ( m1->yx_ * m2->xx_ ) + ( m1->yy_ * m2->yx_ ) + ( m1->yz_ * m2->zx_ );
  y = ( m1->yx_ * m2->xy_ ) + ( m1->yy_ * m2->yy_ ) + ( m1->yz_ * m2->zy_ );
  z = ( m1->yx_ * m2->xz_ ) + ( m1->yy_ * m2->yz_ ) + ( m1->yz_ * m2->zz_ );
  m1->yx_ = x; m1->yy_ = y; m1->yz_ = z;

  // Third row
  x = ( m1->zx_ * m2->xx_ ) + ( m1->zy_ * m2->yx_ ) + ( m1->zz_ * m2->zx_ );
  y = ( m1->zx_ * m2->xy_ ) + ( m1->zy_ * m2->yy_ ) + ( m1->zz_ * m2->zy_ );
  z = ( m1->zx_ * m2->xz_ ) + ( m1->zy_ * m2->yz_ ) + ( m1->zz_ * m2->zz_ );
  m1->zx_ = x; m1->zy_ = y; m1->zz_ = z;
}

/// Modify matrix 1, m1, by performing a right-multiply into it of m2-transposed
void
right_multiply_by_transpose_private(
  __private struct xyzmatrix * m1,
  __private struct xyzmatrix * m2
)
{
  float x, y, z; // Temporaries

  // First row
  x = ( m1->xx_ * m2->xx_ ) + ( m1->xy_ * m2->xy_ ) + ( m1->xz_ * m2->xz_ );
  y = ( m1->xx_ * m2->yx_ ) + ( m1->xy_ * m2->yy_ ) + ( m1->xz_ * m2->yz_ );
  z = ( m1->xx_ * m2->zx_ ) + ( m1->xy_ * m2->zy_ ) + ( m1->xz_ * m2->zz_ );
  m1->xx_ = x; m1->xy_ = y; m1->xz_ = z;

  // Second row
  x = ( m1->yx_ * m2->xx_ ) + ( m1->yy_ * m2->xy_ ) + ( m1->yz_ * m2->xz_ );
  y = ( m1->yx_ * m2->yx_ ) + ( m1->yy_ * m2->yy_ ) + ( m1->yz_ * m2->yz_ );
  z = ( m1->yx_ * m2->zx_ ) + ( m1->yy_ * m2->zy_ ) + ( m1->yz_ * m2->zz_ );
  m1->yx_ = x; m1->yy_ = y; m1->yz_ = z;

  // Third row
  x = ( m1->zx_ * m2->xx_ ) + ( m1->zy_ * m2->xy_ ) + ( m1->zz_ * m2->xz_ );
  y = ( m1->zx_ * m2->yx_ ) + ( m1->zy_ * m2->yy_ ) + ( m1->zz_ * m2->yz_ );
  z = ( m1->zx_ * m2->zx_ ) + ( m1->zy_ * m2->zy_ ) + ( m1->zz_ * m2->zz_ );
  m1->zx_ = x; m1->zy_ = y; m1->zz_ = z;

}


/// Modify matrix 1, m1, by performing a left-multiply of matrix m2 into it.
void
left_multiply_by_private(
  __private struct xyzmatrix * m1,
  __private struct xyzmatrix * m2
)
{
  float x, y, z; // Temporaries

  // First column
  x = ( m2->xx_ * m1->xx_ ) + ( m2->xy_ * m1->yx_ ) + ( m2->xz_ * m1->zx_ );
  y = ( m2->yx_ * m1->xx_ ) + ( m2->yy_ * m1->yx_ ) + ( m2->yz_ * m1->zx_ );
  z = ( m2->zx_ * m1->xx_ ) + ( m2->zy_ * m1->yx_ ) + ( m2->zz_ * m1->zx_ );
  m1->xx_ = x; m1->yx_ = y; m1->zx_ = z;

  // Second column
  x = ( m2->xx_ * m1->xy_ ) + ( m2->xy_ * m1->yy_ ) + ( m2->xz_ * m1->zy_ );
  y = ( m2->yx_ * m1->xy_ ) + ( m2->yy_ * m1->yy_ ) + ( m2->yz_ * m1->zy_ );
  z = ( m2->zx_ * m1->xy_ ) + ( m2->zy_ * m1->yy_ ) + ( m2->zz_ * m1->zy_ );
  m1->xy_ = x; m1->yy_ = y; m1->zy_ = z;

  // Third column
  x = ( m2->xx_ * m1->xz_ ) + ( m2->xy_ * m1->yz_ ) + ( m2->xz_ * m1->zz_ );
  y = ( m2->yx_ * m1->xz_ ) + ( m2->yy_ * m1->yz_ ) + ( m2->yz_ * m1->zz_ );
  z = ( m2->zx_ * m1->xz_ ) + ( m2->zy_ * m1->yz_ ) + ( m2->zz_ * m1->zz_ );
  m1->xz_ = x; m1->yz_ = y; m1->zz_ = z;

}

/// Modify matrix 1, m1, by performing a left-multiply into it of m2-transposed
void
left_multiply_by_transpose_private(
  __private struct xyzmatrix * m1,
  __private struct xyzmatrix * m2
)
{
  float x, y, z; // Temporaries

  // First column
  x = ( m2->xx_ * m1->xx_ ) + ( m2->yx_ * m1->yx_ ) + ( m2->zx_ * m1->zx_ );
  y = ( m2->xy_ * m1->xx_ ) + ( m2->yy_ * m1->yx_ ) + ( m2->zy_ * m1->zx_ );
  z = ( m2->xz_ * m1->xx_ ) + ( m2->yz_ * m1->yx_ ) + ( m2->zz_ * m1->zx_ );
  m1->xx_ = x; m1->yx_ = y; m1->zx_ = z;

  // Second column
  x = ( m2->xx_ * m1->xy_ ) + ( m2->yx_ * m1->yy_ ) + ( m2->zx_ * m1->zy_ );
  y = ( m2->xy_ * m1->xy_ ) + ( m2->yy_ * m1->yy_ ) + ( m2->zy_ * m1->zy_ );
  z = ( m2->xz_ * m1->xy_ ) + ( m2->yz_ * m1->yy_ ) + ( m2->zz_ * m1->zy_ );
  m1->xy_ = x; m1->yy_ = y; m1->zy_ = z;

  // Third column
  x = ( m2->xx_ * m1->xz_ ) + ( m2->yx_ * m1->yz_ ) + ( m2->zx_ * m1->zz_ );
  y = ( m2->xy_ * m1->xz_ ) + ( m2->yy_ * m1->yz_ ) + ( m2->zy_ * m1->zz_ );
  z = ( m2->xz_ * m1->xz_ ) + ( m2->yz_ * m1->yz_ ) + ( m2->zz_ * m1->zz_ );
  m1->xz_ = x; m1->yz_ = y; m1->zz_ = z;

}


/// Set the value for one of the 9 elements in a 3x3 xyzmatrix given a row
/// index / column index pair, indexing from 0.
/// The first row is:     xx, xy, xz;
/// the second row is:    yx, yy, yz;
/// and the third row is: zx, zy, zz;
void set_xyzmatrix_value_private(
  __private struct xyzmatrix * m,
  int row,
  int column,
  float value
)
{
  int rowmajor = row*3 + column;
  switch ( rowmajor ) {
  case 0 :
    m->xx_ = value; break;
  case 1:
    m->xy_ = value; break;
  case 2:
    m->xz_ = value; break;
  case 3 :
    m->yx_ = value; break;
  case 4:
    m->yy_ = value; break;
  case 5:
    m->yz_ = value; break;
  case 6 :
    m->zx_ = value; break;
  case 7:
    m->zy_ = value; break;
  case 8:
    m->zz_ = value; break;
  default :
    // ERROR!
    // The GPU can't exactly throw an error, so
    // try to leave a message in the input xyzmatrix
    // saying something went wrong.
    m->xx_ = m->xy_ = m->xz_ = 0;
    m->yx_ = m->yy_ = m->yz_ = 0;
    m->zx_ = m->zy_ = m->zz_ = 0;
    break;
  }
}


/// Set the value for one of the 9 elements in a 3x3 xyzmatrix given a row
/// index / column index pair, indexing from 0.
/// The first row is:     xx, xy, xz;
/// the second row is:    yx, yy, yz;
/// and the third row is: zx, zy, zz;
float get_xyzmatrix_value_private(
  __private struct xyzmatrix * m,
  int row,
  int column
)
{
  int rowmajor = row*3 + column;
  switch ( rowmajor ) {
  case 0 :
    return m->xx_;
  case 1:
    return m->xy_;
  case 2:
    return m->xz_;
  case 3 :
    return m->yx_;
  case 4:
    return m->yy_;
  case 5:
    return m->yz_;
  case 6 :
    return m->zx_;
  case 7:
    return m->zy_;
  case 8:
    return m->zz_;
  }
  return 1/0; // error rowmajor > 8
}


/// Set the value for one of the 9 elements in a 3x3 xyzmatrix given a row
/// index / column index pair, indexing from 0.
/// The first row is:     xx, xy, xz;
/// the second row is:    yx, yy, yz;
/// and the third row is: zx, zy, zz;
void set_xyzmatrix_value_global(
  __global struct xyzmatrix * m,
  int row,
  int column,
  float value
)
{
  int rowmajor = row*3 + column;
  switch ( rowmajor ) {
  case 0 :
    m->xx_ = value; break;
  case 1:
    m->xy_ = value; break;
  case 2:
    m->xz_ = value; break;
  case 3 :
    m->yx_ = value; break;
  case 4:
    m->yy_ = value; break;
  case 5:
    m->yz_ = value; break;
  case 6 :
    m->zx_ = value; break;
  case 7:
    m->zy_ = value; break;
  case 8:
    m->zz_ = value; break;
  default :
    // ERROR!
    // The GPU can't exactly throw an error, so
    // try to leave a message in the input xyzmatrix
    // saying something went wrong.
    m->xx_ = m->xy_ = m->xz_ = 0;
    m->yx_ = m->yy_ = m->yz_ = 0;
    m->zx_ = m->zy_ = m->zz_ = 0;
    break;
  }
}


/// Set the value for one of the 9 elements in a 3x3 xyzmatrix given a row
/// index / column index pair, indexing from 0.
/// The first row is:     xx, xy, xz;
/// the second row is:    yx, yy, yz;
/// and the third row is: zx, zy, zz;
float get_xyzmatrix_value_global(
  __global struct xyzmatrix * m,
  int row,
  int column
)
{
  int rowmajor = row*3 + column;
  switch ( rowmajor ) {
  case 0 :
    return m->xx_;
  case 1:
    return m->xy_;
  case 2:
    return m->xz_;
  case 3 :
    return m->yx_;
  case 4:
    return m->yy_;
  case 5:
    return m->yz_;
  case 6 :
    return m->zx_;
  case 7:
    return m->zy_;
  case 8:
    return m->zz_;
  }
  return 1/0; // error rowmajor > 8
}

/// Set the value for one of the 9 elements in a 3x3 xyzmatrix given a row
/// index / column index pair, indexing from 0.
/// The first row is:     xx, xy, xz;
/// the second row is:    yx, yy, yz;
/// and the third row is: zx, zy, zz;
void set_xyzmatrix_value_local(
  __local struct xyzmatrix * m,
  int row,
  int column,
  float value
)
{
  int rowmajor = row*3 + column;
  switch ( rowmajor ) {
  case 0 :
    m->xx_ = value; break;
  case 1:
    m->xy_ = value; break;
  case 2:
    m->xz_ = value; break;
  case 3 :
    m->yx_ = value; break;
  case 4:
    m->yy_ = value; break;
  case 5:
    m->yz_ = value; break;
  case 6 :
    m->zx_ = value; break;
  case 7:
    m->zy_ = value; break;
  case 8:
    m->zz_ = value; break;
  default :
    // ERROR!
    // The GPU can't exactly throw an error, so
    // try to leave a message in the input xyzmatrix
    // saying something went wrong.
    m->xx_ = m->xy_ = m->xz_ = 0;
    m->yx_ = m->yy_ = m->yz_ = 0;
    m->zx_ = m->zy_ = m->zz_ = 0;
    break;
  }
}


/// Set the value for one of the 9 elements in a 3x3 xyzmatrix given a row
/// index / column index pair, indexing from 0.
/// The first row is:     xx, xy, xz;
/// the second row is:    yx, yy, yz;
/// and the third row is: zx, zy, zz;
float get_xyzmatrix_value_local(
  __local struct xyzmatrix * m,
  int row,
  int column
)
{
  int rowmajor = row*3 + column;
  switch ( rowmajor ) {
  case 0 :
    return m->xx_;
  case 1:
    return m->xy_;
  case 2:
    return m->xz_;
  case 3 :
    return m->yx_;
  case 4:
    return m->yy_;
  case 5:
    return m->yz_;
  case 6 :
    return m->zx_;
  case 7:
    return m->zy_;
  case 8:
    return m->zz_;
  }
  return 1/0; // error rowmajor > 8
}

/// set a value in an xyzvector, indexing from 0
void
set_xyzvector_value_private(
  __private struct xyzvector * v,
  int index,
  float value
)
{
  switch ( index ) {
  case 0 :
    v->x_ = value; break;
  case 1:
    v->y_ = value; break;
  case 2:
    v->z_ = value; break;
  } 
}

/// get a value out of an xyzvector, indexing from 0
float
get_xyzvector_value_private(
  __private struct xyzvector * v,
  int index
)
{
  switch ( index ) {
  case 0 :
    return v->x_;
  case 1 :
    return v->y_;
  case 2 :
    return v->z_;
  }
}


/// set a value in an xyzvector, indexing from 0
void
set_xyzvector_value_local(
  __local struct xyzvector * v,
  int index,
  float value
)
{
  switch ( index ) {
  case 0 :
    v->x_ = value; break;
  case 1:
    v->y_ = value; break;
  case 2:
    v->z_ = value; break;
  } 
}

/// get a value out of an xyzvector, indexing from 0
float
get_xyzvector_value_local(
  __local struct xyzvector * v,
  int index
)
{
  switch ( index ) {
  case 0 :
    return v->x_;
  case 1 :
    return v->y_;
  case 2 :
    return v->z_;
  }
}


/// set a value in an xyzvector, indexing from 0
void
set_xyzvector_value_global(
  __global struct xyzvector * v,
  int index,
  float value
)
{
  switch ( index ) {
  case 0 :
    v->x_ = value; break;
  case 1:
    v->y_ = value; break;
  case 2:
    v->z_ = value; break;
  } 
}

/// get a value out of an xyzvector, indexing from 0
float
get_xyzvector_value_global(
  __global struct xyzvector * v,
  int index
)
{
  switch ( index ) {
  case 0 :
    return v->x_;
  case 1 :
    return v->y_;
  case 2 :
    return v->z_;
  }
}

void
set_xyzmatrix_to_identity_private(
  __private struct xyzmatrix * m
)
{
  m->xx_ = 1.0f; m->xy_ = 0.0f; m->xz_ = 0.0f;
  m->yx_ = 0.0f; m->yy_ = 1.0f; m->yz_ = 0.0f;
  m->zx_ = 0.0f; m->zy_ = 0.0f; m->zz_ = 1.0f;
}

/// swap the off-diagonal entries
void
xyzmatrix_transpose_private(
  __private struct xyzmatrix * m
)
{
  float temp;
  temp = m->xy_;
  m->xy_ = m->yx_;
  m->yx_ = temp;

  temp = m->xz_;
  m->xz_ = m->zx_;
  m->zx_ = temp;

  temp = m->yz_;
  m->yz_ = m->zy_;
  m->zy_ = temp;
}


/// compute m*v
struct xyzvector
xyzmatrix_xyzvector_multiply_private(
  __private struct xyzmatrix * m,
  __private struct xyzvector * v
)
{
  struct xyzvector p;
  p.x_ = m->xx_ * v->x_ + m->xy_ * v->y_ + m->xz_ * v->z_;
  p.y_ = m->yx_ * v->x_ + m->yy_ * v->y_ + m->yz_ * v->z_;
  p.z_ = m->zx_ * v->x_ + m->zy_ * v->y_ + m->zz_ * v->z_;
  return p;
}

float
xyzvector_dot_product_private(
  __private struct xyzvector * v1,
  __private struct xyzvector * v2
)
{
  return v1->x_*v2->x_ + v1->y_*v2->y_ + v1->z_*v2->z_;
}

// compute the square magnitude of an xyzvector
float
xyzvector_square_magnitude_private(
  __private struct xyzvector * v
)
{
  return v->x_*v->x_ + v->y_*v->y_ + v->z_*v->z_;
}

/// compute the square distance between two points, held in xyzvectors
float
xyzvector_square_distance_private(
  __private struct xyzvector * v1,
  __private struct xyzvector * v2
)
{
  float sum2 = 0;
  float val;
  val = v1->x_ - v2->x_;
  sum2 += val*val;
  val = v1->y_ - v2->y_;
  sum2 += val*val;
  val = v1->z_ - v2->z_;
  sum2 += val*val;
  return sum2;
}

