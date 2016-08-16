// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#include "boost/python.hpp"

#include "numeric/xyzVector.hh"
#include "numeric/xyzTriple.hh"
#include "numeric/xyzMatrix.hh"
#include "numeric/types.hh"
#include "numeric/Quaternion.hh"
#include "numeric/EulerAngles.hh"
#include "numeric/BodyPosition.hh"
#include <numeric/io.hh>

#include <vector>

namespace bp = boost::python;
using namespace numeric;
using std::size_t;
using std::vector;

typedef bp::return_value_policy< bp::reference_existing_object > CP_REF;
typedef bp::return_value_policy< bp::copy_const_reference >      CP_CCR;
typedef bp::return_value_policy< bp::copy_non_const_reference >  CP_CNCR;

template<class T>
void instantiate_numeric_functions(std::string type_name)
{
  { xyzVector<T>   (*t)( xyzMatrix<T> const &, xyzVector<T> const &) = &(operator *); }
  { xyzVector<T>   (*t)( xyzMatrix<T> const &, xyzVector<T> const &) = &product; }
  { xyzVector<T> & (*t)( xyzMatrix<T> const &, xyzVector<T>       &) = &inplace_product; }
  { xyzVector<T>   (*t)( xyzMatrix<T> const &, xyzVector<T> const &) = &transpose_product; }
  { xyzVector<T> & (*t)( xyzMatrix<T> const &, xyzVector<T>       &) = &inplace_transpose_product; }
  { xyzMatrix<T>   (*t)( xyzVector<T> const &, xyzVector<T> const &) = &outer_product; }
  { xyzMatrix<T>   (*t)( xyzVector<T> const & ) = &projection_matrix; }

  { // numeric::rotation_matrix
    typedef xyzMatrix<T>   (*function_type)( xyzVector<T> const &, T const & );
    std::string docstring = "Rotation matrix for rotation about an axis by an angle in radians.";

    bp::def(
        ("rotation_matrix_" + type_name).c_str()
        , function_type(&rotation_matrix)
        , ( bp::arg("axis"), bp::arg("theta") )
        , bp::return_value_policy< bp::return_by_value >()
        , docstring);

    bp::def(
        "rotation_matrix"
        , function_type(&rotation_matrix)
        , ( bp::arg("axis"), bp::arg("theta") )
        , bp::return_value_policy< bp::return_by_value >()
        , docstring);
  }

  { // numeric::rotation_matrix_degrees
    typedef xyzMatrix<T>   (*function_type)( xyzVector<T> const &, T const & );
    std::string docstring = "Rotation matrix for rotation about an axis by an angle in degrees.";

    bp::def(
        ("rotation_matrix_degrees_" + type_name).c_str()
        , function_type( &rotation_matrix_degrees )
        , ( bp::arg("axis"), bp::arg("theta") )
        , bp::return_value_policy< bp::return_by_value >()
        , docstring);

    bp::def(
        "rotation_matrix_degrees"
        , function_type( &rotation_matrix_degrees )
        , ( bp::arg("axis"), bp::arg("theta") )
        , bp::return_value_policy< bp::return_by_value >()
        , docstring);
  }

  { // numeric::rotation_matrix_radians
    typedef xyzMatrix<T>   (*function_type)( xyzVector<T> const &, T const & );
    std::string docstring = "Rotation matrix for rotation about an axis by an angle in radians.";

    bp::def(
        ("rotation_matrix_radians_" + type_name).c_str()
        , function_type( &rotation_matrix_radians )
        , ( bp::arg("axis"), bp::arg("theta") )
        , bp::return_value_policy< bp::return_by_value >()
        , docstring);

    bp::def(
        "rotation_matrix_radians"
        , function_type( &rotation_matrix_radians )
        , ( bp::arg("axis"), bp::arg("theta") )
        , bp::return_value_policy< bp::return_by_value >()
        , docstring);
  }


  { //numeric::rotation_angle
    typedef T (*function_type)( xyzMatrix<T> const &);
    std::string docstring = \
         "Transformation from rotation matrix to magnitude of helical rotation, input matrix must be orthogonal.\nOrientation of axis chosen so that the angle of rotation is non-negative [0,pi].\nnumeric::rotation_axis returns both axis and angle of rotation.";

    bp::def(
        ("rotation_angle_" + type_name).c_str()
        , function_type( &rotation_angle )
        , ( bp::arg("rotation_matrix"))
        , bp::return_value_policy< bp::return_by_value >()
        , docstring);

    bp::def(
        "rotation_angle"
        , function_type( &rotation_angle )
        , ( bp::arg("rotation_matrix"))
        , bp::return_value_policy< bp::return_by_value >()
        , docstring);
  }
  { //numeric::rotation_axis_angle
    typedef xyzVector<T>   (*function_type)( xyzMatrix<T> const &);
    std::string docstring = \
          "Transformation from rotation matrix to compact axis-angle representation\nInput matrix must be orthogonal\nOrientation of axis chosen so that the angle of rotation is non-negative [0,pi]\nResulting vector will be oriented in axis of rotation with magnitude equal to magnitude of rotation.";

    bp::def(
        ("rotation_axis_angle_" + type_name).c_str()
        , function_type( &rotation_axis_angle )
        , ( bp::arg("rotation_matrix"))
        , bp::return_value_policy< bp::return_by_value >()
        , docstring);

    bp::def(
        "rotation_axis_angle"
        , function_type( &rotation_axis_angle )
        , ( bp::arg("rotation_matrix"))
        , bp::return_value_policy< bp::return_by_value >()
        , docstring);
  }

  { // numeric::rotation_matrix
    typedef xyzMatrix<T>   (*function_type)( xyzVector<T> const &);

    std::string docstring = "Rotation matrix for rotation from axis-angle representation.\nMagnitude of rotation (in radians) is taken as axis_angle.magnitude().";

    bp::def(
        ("rotation_matrix_" + type_name).c_str()
        , function_type( &rotation_matrix )
        , ( bp::arg("axis_angle") )
        , bp::return_value_policy< bp::return_by_value >()
        , docstring );

    bp::def(
        "rotation_matrix"
        , function_type( &rotation_matrix )
        , ( bp::arg("axis_angle") )
        , bp::return_value_policy< bp::return_by_value >()
        , docstring );
  }

  { xyzMatrix<T>   (*t)( T const & ) = &x_rotation_matrix; }
  { xyzMatrix<T>   (*t)( T const & ) = &x_rotation_matrix_radians; }
  { xyzMatrix<T>   (*t)( T const & ) = &x_rotation_matrix_degrees; }
  { xyzMatrix<T>   (*t)( T const & ) = &y_rotation_matrix; }
  { xyzMatrix<T>   (*t)( T const & ) = &y_rotation_matrix_radians; }
  { xyzMatrix<T>   (*t)( T const & ) = &y_rotation_matrix_degrees; }
  { xyzMatrix<T>   (*t)( T const & ) = &z_rotation_matrix; }
  { xyzMatrix<T>   (*t)( T const & ) = &z_rotation_matrix_radians; }
  { xyzMatrix<T>   (*t)( T const & ) = &z_rotation_matrix_degrees; }

  { xyzVector<T>   (*t)( xyzMatrix<T> const &, T &) = &rotation_axis; }

  { xyzVector<T>   (*t)( xyzMatrix<T> const &, T const & ) = &eigenvalue_jacobi; }
  { xyzVector<T>   (*t)( xyzMatrix<T> const &, T const &, xyzMatrix<T> & ) = &eigenvector_jacobi; }
  { void (*t)( xyzMatrix< T > const &, int const, int const, xyzMatrix<T> & ) = &jacobi_rotation; }

  { //::numeric::Quaternion< T >
      typedef bp::class_< numeric::Quaternion< T > > Quaternion_typename_exposer_t;
      Quaternion_typename_exposer_t Quaternion_typename_exposer = Quaternion_typename_exposer_t( std::string("Quaternion_" + type_name).c_str() );
      bp::scope Quaternion_typename_scope( Quaternion_typename_exposer );
      Quaternion_typename_exposer.def( bp::init< >() );
      Quaternion_typename_exposer.def( bp::init<  T  const &,  T  const &,  T  const &,  T  const &, bp::optional< bool > >(( bp::arg("w_a"), bp::arg("x_a"), bp::arg("y_a"), bp::arg("z_a"), bp::arg("precise")=(bool const)(true) )) );
      { //::numeric::Quaternion< T >::I

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::Quaternion< T > const & ( *I_function_type )(  );

          Quaternion_typename_exposer.def(
              "I"
              , I_function_type( &::numeric::Quaternion< T >::I )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::apply

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::Quaternion< T > & ( exported_class_t::*apply_function_type )( ::numeric::Quaternion< T > const &,bool const ) ;

          Quaternion_typename_exposer.def(
              "apply"
              , apply_function_type( &::numeric::Quaternion< T >::apply )
              , ( bp::arg("q"), bp::arg("precise")=(bool const)(true) )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::axis

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*axis_function_type )(  ) const;

          Quaternion_typename_exposer.def(
              "axis"
              , axis_function_type( &::numeric::Quaternion< T >::axis )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::axis

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*axis_function_type )( ::numeric::xyzVector< T > & ) const;

          Quaternion_typename_exposer.def(
              "axis"
              , axis_function_type( &::numeric::Quaternion< T >::axis )
              , ( bp::arg("u") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::conjugate

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::Quaternion< T > & ( exported_class_t::*conjugate_function_type )(  ) ;

          Quaternion_typename_exposer.def(
              "conjugate"
              , conjugate_function_type( &::numeric::Quaternion< T >::conjugate )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::dot

          typedef numeric::Quaternion< T > exported_class_t;
          typedef  T  ( exported_class_t::*dot_function_type )( ::numeric::Quaternion< T > const & ) const;

          Quaternion_typename_exposer.def(
              "dot"
              , dot_function_type( &::numeric::Quaternion< T >::dot )
              , ( bp::arg("q") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::dot_product

          typedef numeric::Quaternion< T > exported_class_t;
          typedef  T  ( exported_class_t::*dot_product_function_type )( ::numeric::Quaternion< T > const & ) const;

          Quaternion_typename_exposer.def(
              "dot_product"
              , dot_product_function_type( &::numeric::Quaternion< T >::dot_product )
              , ( bp::arg("q") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::identity

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::Quaternion< T > ( *identity_function_type )(  );

          Quaternion_typename_exposer.def(
              "identity"
              , identity_function_type( &::numeric::Quaternion< T >::identity )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::invert

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::Quaternion< T > & ( exported_class_t::*invert_function_type )(  ) ;

          Quaternion_typename_exposer.def(
              "invert"
              , invert_function_type( &::numeric::Quaternion< T >::invert )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::is_normalized

          typedef numeric::Quaternion< T > exported_class_t;
          typedef bool ( exported_class_t::*is_normalized_function_type )(  T  const & ) const;

          Quaternion_typename_exposer.def(
              "is_normalized"
              , is_normalized_function_type( &::numeric::Quaternion< T >::is_normalized )
              , ( bp::arg("tol")=numeric::NumericTraits< T >::quaternion_tolerance() )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::left_multiply_by

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::Quaternion< T > & ( exported_class_t::*left_multiply_by_function_type )( ::numeric::Quaternion< T > const &,bool const ) ;

          Quaternion_typename_exposer.def(
              "left_multiply_by"
              , left_multiply_by_function_type( &::numeric::Quaternion< T >::left_multiply_by )
              , ( bp::arg("q"), bp::arg("precise")=(bool const)(true) )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::left_multiply_by_inverse_of

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::Quaternion< T > & ( exported_class_t::*left_multiply_by_inverse_of_function_type )( ::numeric::Quaternion< T > const &,bool const ) ;

          Quaternion_typename_exposer.def(
              "left_multiply_by_inverse_of"
              , left_multiply_by_inverse_of_function_type( &::numeric::Quaternion< T >::left_multiply_by_inverse_of )
              , ( bp::arg("q"), bp::arg("precise")=(bool const)(true) )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::normalize

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::Quaternion< T > & ( exported_class_t::*normalize_function_type )(  ) ;

          Quaternion_typename_exposer.def(
              "normalize"
              , normalize_function_type( &::numeric::Quaternion< T >::normalize )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::normalize_if_needed

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::Quaternion< T > & ( exported_class_t::*normalize_if_needed_function_type )(  T  const & ) ;

          Quaternion_typename_exposer.def(
              "normalize_if_needed"
              , normalize_if_needed_function_type( &::numeric::Quaternion< T >::normalize_if_needed )
              , ( bp::arg("tol")=numeric::NumericTraits< T >::quaternion_tolerance() )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::not_normalized

          typedef numeric::Quaternion< T > exported_class_t;
          typedef bool ( exported_class_t::*not_normalized_function_type )(  T  const & ) const;

          Quaternion_typename_exposer.def(
              "not_normalized"
              , not_normalized_function_type( &::numeric::Quaternion< T >::not_normalized )
              , ( bp::arg("tol")=numeric::NumericTraits< T >::quaternion_tolerance() )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::right_multiply_by

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::Quaternion< T > & ( exported_class_t::*right_multiply_by_function_type )( ::numeric::Quaternion< T > const &,bool const ) ;

          Quaternion_typename_exposer.def(
              "right_multiply_by"
              , right_multiply_by_function_type( &::numeric::Quaternion< T >::right_multiply_by )
              , ( bp::arg("q"), bp::arg("precise")=(bool const)(true) )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::right_multiply_by_inverse_of

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::Quaternion< T > & ( exported_class_t::*right_multiply_by_inverse_of_function_type )( ::numeric::Quaternion< T > const &,bool const ) ;

          Quaternion_typename_exposer.def(
              "right_multiply_by_inverse_of"
              , right_multiply_by_inverse_of_function_type( &::numeric::Quaternion< T >::right_multiply_by_inverse_of )
              , ( bp::arg("q"), bp::arg("precise")=(bool const)(true) )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::swap

          typedef numeric::Quaternion< T > exported_class_t;
          typedef void ( exported_class_t::*swap_function_type )( ::numeric::Quaternion< T > & ) ;

          Quaternion_typename_exposer.def(
              "swap"
              , swap_function_type( &::numeric::Quaternion< T >::swap )
              , ( bp::arg("q") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::Quaternion< T >::to_identity

          typedef numeric::Quaternion< T > exported_class_t;
          typedef ::numeric::Quaternion< T > & ( exported_class_t::*to_identity_function_type )(  ) ;

          Quaternion_typename_exposer.def(
              "to_identity"
              , to_identity_function_type( &::numeric::Quaternion< T >::to_identity )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      Quaternion_typename_exposer.staticmethod( "I" );
      Quaternion_typename_exposer.staticmethod( "identity" );
      { //property "w"[fget=::numeric::Quaternion< T >::w]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  const & ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "w"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::w )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "x"[fget=::numeric::Quaternion< T >::x]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  const & ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "x"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::x )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "y"[fget=::numeric::Quaternion< T >::y]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  const & ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "y"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::y )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "z"[fget=::numeric::Quaternion< T >::z]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  const & ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "z"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::z )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "w_squared"[fget=::numeric::Quaternion< T >::w_squared]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "w_squared"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::w_squared )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "x_squared"[fget=::numeric::Quaternion< T >::x_squared]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "x_squared"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::x_squared )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "y_squared"[fget=::numeric::Quaternion< T >::y_squared]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "y_squared"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::y_squared )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "z_squared"[fget=::numeric::Quaternion< T >::z_squared]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "z_squared"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::z_squared )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "norm"[fget=::numeric::Quaternion< T >::norm]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "norm"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::norm )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "norm_squared"[fget=::numeric::Quaternion< T >::norm_squared]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "norm_squared"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::norm_squared )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "norm_error"[fget=::numeric::Quaternion< T >::norm_error]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "norm_error"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::norm_error )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "norm_squared_error"[fget=::numeric::Quaternion< T >::norm_squared_error]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "norm_squared_error"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::norm_squared_error )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "magnitude"[fget=::numeric::Quaternion< T >::magnitude]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "magnitude"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::magnitude )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "magnitude_squared"[fget=::numeric::Quaternion< T >::magnitude_squared]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "magnitude_squared"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::magnitude_squared )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "magnitude_error"[fget=::numeric::Quaternion< T >::magnitude_error]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "magnitude_error"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::magnitude_error )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "magnitude_squared_error"[fget=::numeric::Quaternion< T >::magnitude_squared_error]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "magnitude_squared_error"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::magnitude_squared_error )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "angle"[fget=::numeric::Quaternion< T >::angle]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "angle"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::angle )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "conjugated"[fget=::numeric::Quaternion< T >::conjugated]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef ::numeric::Quaternion< T > ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "conjugated"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::conjugated )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "inverse"[fget=::numeric::Quaternion< T >::inverse]

          typedef numeric::Quaternion< T > fget_class_t;

          typedef ::numeric::Quaternion< T > ( fget_class_t::*fget )(  ) const;

          Quaternion_typename_exposer.add_property(
              "inverse"
              , bp::make_function(
                    fget( &::numeric::Quaternion< T >::inverse )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      Quaternion_typename_exposer.def( bp::self != bp::self );
      Quaternion_typename_exposer.def( bp::self * bp::self );
      Quaternion_typename_exposer.def( bp::self == bp::self );

	  Quaternion_typename_exposer.def( bp::self_ns::str( bp::self ) );

  }

  { // ::numeric::EulerAngles<T>
    utility::py::wrap_access_pointer< ::numeric::EulerAngles<T> >(("EulerAngles_" + type_name).c_str());
    typedef boost::python::class_< ::numeric::EulerAngles<T>, boost::python::bases< ::numeric::xyzVector<T> > > EulerAngles_typename_exposer_type;
    EulerAngles_typename_exposer_type EulerAngles_typename_exposer(("EulerAngles_" + type_name).c_str(), "Euler angles 3-D orientation representation\n@remarks\nThe three euler angles (in radians) that describing a rotation operation\nof a Z axis rotation by the angle phi (position 1), followed by\nan X axis rotation by the angle theta (position 3), followed by another\nZ axis rotation by the angle psi (position 2).\nthis->code is a modified version of Alex Z's code from r++.\n@details\nThe range of phi is [ -pi, pi ];\nThe range of psi is [ -pi, pi ];\nThe range of theta is [ 0, pi ];\n", boost::python::init <  >() );
    EulerAngles_typename_exposer.def( boost::python::init< ::numeric::EulerAngles<T> const & > ( (boost::python::arg("")) , "Euler angles 3-D orientation representation\n@remarks\nThe three euler angles (in radians) that describing a rotation operation\nof a Z axis rotation by the angle phi (position 1), followed by\nan X axis rotation by the angle theta (position 3), followed by another\nZ axis rotation by the angle psi (position 2).\nthis->code is a modified version of Alex Z's code from r++.\n@details\nThe range of phi is [ -pi, pi ];\nThe range of psi is [ -pi, pi ];\nThe range of theta is [ 0, pi ];\n" ) );
    EulerAngles_typename_exposer.def( boost::python::init< ::numeric::xyzVector<T> const & > ( (boost::python::arg("v")) , "Copy constructor\n" ) );
    EulerAngles_typename_exposer.def( boost::python::init< T const &, T const &, T const & > ( (boost::python::arg("phi"), boost::python::arg("psi"), boost::python::arg("theta")) , "Triple value constructor\n" ) );
    EulerAngles_typename_exposer.def( boost::python::init< ::numeric::xyzMatrix<T> > ( (boost::python::arg("rotation_matrix")) , "Construct from rotation matrix.\n" ) );

    { // ::numeric::EulerAngles<T>::from_rotation_matrix
      typedef void ( ::numeric::EulerAngles<T>:: * from_rotation_matrix_function_type)(::numeric::xyzMatrix<T> __a0) ;

      EulerAngles_typename_exposer.def("from_rotation_matrix"
        , from_rotation_matrix_function_type( &::numeric::EulerAngles<T>::from_rotation_matrix )
        , ( boost::python::arg("matrix") )
        , "The equivalent rotation matrix representation of the euler angles would be:\nFIGURE 1:\nR = [\n      cos(psi)cos(phi)-cos(theta)sin(phi)sin(psi)        cos(psi)sin(phi)+cos(theta)cos(phi)sin(psi)      sin(psi)sin(theta)\n     -sin(psi)cos(phi)-cos(theta)sin(phi)cos(psi)       -sin(psi)sin(phi)+cos(theta)cos(phi)cos(psi)      cos(psi)sin(theta)\n                  sin(theta)sin(phi)                                 -sin(theta)cos(phi)                        cos(theta)\n]\nThe zz_ coordinate gives away theta.\nTheta may be computed as acos( zz_ ), or, as Alex does it, asin( sqrt( 1 - zz^2))\nSince there is redundancy in theta, this->function chooses a theta with a positive\nsin(theta): i.e. quadrants I and II.  Assuming we have a positive sin theta\npushes phi and psi into conforming angles.\nNOTE on theta: asin returns a value in the range [ -pi/2, pi/2 ], and we have artificially\ncreated a positive sin(theta), so we will get a asin( pos_sin_theta ), we have a value\nin the range [ 0, pi/2 ].  To convert this->into the actual angle theta, we examine the zz sign.\nIf zz is negative, we chose the quadrant II theta.\nThat is, asin( pos_sin_theta) returned an angle, call it theta'.  Now, if cos( theta ) is negative,\nthen we want to choose the positive x-axis rotation that's equivalent to -theta'.  To do so,\nwe reflect q through the y axis (See figure 2 below) to get p and then measure theta as pi - theta'.\nFIGURE 2:\n II        |         I\n           |\n   p.      |      .q (cos(-theta'), abs(sin(theta')))\n      .    |    .\ntheta'( .  |  .  )  theta' = asin( abs(sin(theta))\n-----------------------\n           |\n           |\n           |\n III       |        IV\n           |\n The angle between the positive x axis and p is pi - theta'.\nSince zx and zy contain only phi terms and a constant sin( theta ) term,\nphi is given by atan2( sin_phi, cos_phi ) = atan2( c*sin_phi, c*cos_phi ) = atan2( zx, -zy )\nfor c positive and non-zero.  If sin_theta is zero, or very close to zero, we're at gimbal lock.\nMoreover, since xz and yz contain only psi terms, psi may also be deduced using atan2.\nThere are 2 degenerate cases (gimbal lock)\n1. theta close to 0  (North Pole singularity), or\n2. theta close to pi (South Pole singularity)\nFor these, we take: phi=acos(xx), theta = 0 (resp. Pi/2), psi = 0\n" );
    }


    { // ::numeric::EulerAngles<T>::to_rotation_matrix
      typedef ::numeric::xyzMatrix<T> ( ::numeric::EulerAngles<T>:: * to_rotation_matrix_function_type)() ;

      EulerAngles_typename_exposer.def("to_rotation_matrix"
        , to_rotation_matrix_function_type( &::numeric::EulerAngles<T>::to_rotation_matrix )
        , "Construct rotation matrix from three euler angles that describe the frame.\nSee the description for from_rotation_matrix to understand\nthe Z-X-Z transformation convention.\n" );
    }


    { // ::numeric::EulerAngles<T>::phi
      typedef T & ( ::numeric::EulerAngles<T>:: * phi_function_type)() ;

      EulerAngles_typename_exposer.def("phi"
        , phi_function_type( &::numeric::EulerAngles<T>::phi )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "Value phi in radians\n" );
    }


    { // ::numeric::EulerAngles<T>::phi
      typedef void ( ::numeric::EulerAngles<T>:: * phi_function_type)(T const & __a0) ;

      EulerAngles_typename_exposer.def("phi"
        , phi_function_type( &::numeric::EulerAngles<T>::phi )
        , ( boost::python::arg("value") )
        , "Set value phi in radians\n" );
    }


    { // ::numeric::EulerAngles<T>::phi_radians
      typedef T & ( ::numeric::EulerAngles<T>:: * phi_radians_function_type)() ;

      EulerAngles_typename_exposer.def("phi_radians"
        , phi_radians_function_type( &::numeric::EulerAngles<T>::phi_radians )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "Value phi in radians\n" );
    }


    { // ::numeric::EulerAngles<T>::phi_radians
      typedef void ( ::numeric::EulerAngles<T>:: * phi_radians_function_type)(T const & __a0) ;

      EulerAngles_typename_exposer.def("phi_radians"
        , phi_radians_function_type( &::numeric::EulerAngles<T>::phi_radians )
        , ( boost::python::arg("value") )
        , "Set value phi in radians\n" );
    }


    { // ::numeric::EulerAngles<T>::phi_degrees
      typedef T ( ::numeric::EulerAngles<T>:: * phi_degrees_function_type)() ;

      EulerAngles_typename_exposer.def("phi_degrees"
        , phi_degrees_function_type( &::numeric::EulerAngles<T>::phi_degrees )
        , "Value phi in degrees\n" );
    }


    { // ::numeric::EulerAngles<T>::phi_degrees
      typedef void ( ::numeric::EulerAngles<T>:: * phi_degrees_function_type)(T const & __a0) ;

      EulerAngles_typename_exposer.def("phi_degrees"
        , phi_degrees_function_type( &::numeric::EulerAngles<T>::phi_degrees )
        , ( boost::python::arg("value") )
        , "Set value phi in degrees\n" );
    }


    { // ::numeric::EulerAngles<T>::psi
      typedef T & ( ::numeric::EulerAngles<T>:: * psi_function_type)() ;

      EulerAngles_typename_exposer.def("psi"
        , psi_function_type( &::numeric::EulerAngles<T>::psi )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "Value psi in radians\n" );
    }


    { // ::numeric::EulerAngles<T>::psi
      typedef void ( ::numeric::EulerAngles<T>:: * psi_function_type)(T const & __a0) ;

      EulerAngles_typename_exposer.def("psi"
        , psi_function_type( &::numeric::EulerAngles<T>::psi )
        , ( boost::python::arg("value") )
        , "Set value psi in radians\n" );
    }


    { // ::numeric::EulerAngles<T>::psi_radians
      typedef T & ( ::numeric::EulerAngles<T>:: * psi_radians_function_type)() ;

      EulerAngles_typename_exposer.def("psi_radians"
        , psi_radians_function_type( &::numeric::EulerAngles<T>::psi_radians )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "Value psi in radians\n" );
    }


    { // ::numeric::EulerAngles<T>::psi_radians
      typedef void ( ::numeric::EulerAngles<T>:: * psi_radians_function_type)(T const & __a0) ;

      EulerAngles_typename_exposer.def("psi_radians"
        , psi_radians_function_type( &::numeric::EulerAngles<T>::psi_radians )
        , ( boost::python::arg("value") )
        , "Set value psi in radians\n" );
    }


    { // ::numeric::EulerAngles<T>::psi_degrees
      typedef T ( ::numeric::EulerAngles<T>:: * psi_degrees_function_type)() ;

      EulerAngles_typename_exposer.def("psi_degrees"
        , psi_degrees_function_type( &::numeric::EulerAngles<T>::psi_degrees )
        , "Value psi in degrees\n" );
    }


    { // ::numeric::EulerAngles<T>::psi_degrees
      typedef void ( ::numeric::EulerAngles<T>:: * psi_degrees_function_type)(T const & __a0) ;

      EulerAngles_typename_exposer.def("psi_degrees"
        , psi_degrees_function_type( &::numeric::EulerAngles<T>::psi_degrees )
        , ( boost::python::arg("value") )
        , "Set value psi in degrees\n" );
    }


    { // ::numeric::EulerAngles<T>::theta
      typedef T & ( ::numeric::EulerAngles<T>:: * theta_function_type)() ;

      EulerAngles_typename_exposer.def("theta"
        , theta_function_type( &::numeric::EulerAngles<T>::theta )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "Value theta in radians\n" );
    }


    { // ::numeric::EulerAngles<T>::theta
      typedef void ( ::numeric::EulerAngles<T>:: * theta_function_type)(T const & __a0) ;

      EulerAngles_typename_exposer.def("theta"
        , theta_function_type( &::numeric::EulerAngles<T>::theta )
        , ( boost::python::arg("value") )
        , "Set value theta in radians\n" );
    }


    { // ::numeric::EulerAngles<T>::theta_radians
      typedef T & ( ::numeric::EulerAngles<T>:: * theta_radians_function_type)() ;

      EulerAngles_typename_exposer.def("theta_radians"
        , theta_radians_function_type( &::numeric::EulerAngles<T>::theta_radians )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "Value theta in radians\n" );
    }


    { // ::numeric::EulerAngles<T>::theta_radians
      typedef void ( ::numeric::EulerAngles<T>:: * theta_radians_function_type)(T const & __a0) ;

      EulerAngles_typename_exposer.def("theta_radians"
        , theta_radians_function_type( &::numeric::EulerAngles<T>::theta_radians )
        , ( boost::python::arg("value") )
        , "Set value theta in radians\n" );
    }


    { // ::numeric::EulerAngles<T>::theta_degrees
      typedef T ( ::numeric::EulerAngles<T>:: * theta_degrees_function_type)() ;

      EulerAngles_typename_exposer.def("theta_degrees"
        , theta_degrees_function_type( &::numeric::EulerAngles<T>::theta_degrees )
        , "Value theta in degrees\n" );
    }


    { // ::numeric::EulerAngles<T>::theta_degrees
      typedef void ( ::numeric::EulerAngles<T>:: * theta_degrees_function_type)(T const & __a0) ;

      EulerAngles_typename_exposer.def("theta_degrees"
        , theta_degrees_function_type( &::numeric::EulerAngles<T>::theta_degrees )
        , ( boost::python::arg("value") )
        , "Set value theta in degrees\n" );
    }


    { // ::numeric::EulerAngles<T>::from_degrees
      typedef ::numeric::EulerAngles<T> ( * from_degrees_function_type)(T __phi, T __psi, T __theta);

      EulerAngles_typename_exposer.def("from_degrees"
        , from_degrees_function_type( &::numeric::EulerAngles<T>::from_degrees )
        , ( boost::python::arg("phi"), boost::python::arg("psi"), boost::python::arg("theta") )
        , "Static constructor from degrees" );
    }

    { // ::numeric::EulerAngles<T>::from_degrees
      typedef ::numeric::EulerAngles<T> ( * from_degrees_function_type)(xyzVector<T> __vector);

      EulerAngles_typename_exposer.def("from_degrees"
        , from_degrees_function_type( &::numeric::EulerAngles<T>::from_degrees )
        , boost::python::arg("vector")
        , "Static constructor from degrees" );
    }

		EulerAngles_typename_exposer.staticmethod("from_degrees");

    { // ::numeric::EulerAngles<T>::from_radians
      typedef ::numeric::EulerAngles<T> ( * from_radians_function_type)(T __phi, T __psi, T __theta);

      EulerAngles_typename_exposer.def("from_radians"
        , from_radians_function_type( &::numeric::EulerAngles<T>::from_radians )
        , ( boost::python::arg("phi"), boost::python::arg("psi"), boost::python::arg("theta") )
        , "Static constructor from radians" );
    }

    { // ::numeric::EulerAngles<T>::from_radians
      typedef ::numeric::EulerAngles<T> ( * from_radians_function_type)(xyzVector<T> __vector);

      EulerAngles_typename_exposer.def("from_radians"
        , from_radians_function_type( &::numeric::EulerAngles<T>::from_radians )
        , boost::python::arg("vector")
        , "Static constructor from radians" );
    }

		EulerAngles_typename_exposer.staticmethod("from_radians");

    { // ::numeric::EulerAngles<T>::angular_distance_between
      typedef T ( * angular_distance_between_function_type)(::numeric::EulerAngles<T> __a0, ::numeric::EulerAngles<T> __a1);

      EulerAngles_typename_exposer.def("angular_distance_between"
        , angular_distance_between_function_type( &::numeric::EulerAngles<T>::angular_distance_between )
        , ( boost::python::arg("a1"), boost::python::arg("a2") )
        , "Get angular distance between two sets of Euler Angles." );
    }

		EulerAngles_typename_exposer.staticmethod("angular_distance_between");

  }
}

template<class T>
void instantiate_real_numeric_containers(std::string type_name)
{
  { // ::numeric::xyzTransform< T > 

    utility::py::wrap_access_pointer< ::numeric::xyzTransform< T > >( std::string("xyzTransform_" + type_name).c_str()  );

    typedef boost::python::class_< ::numeric::xyzTransform< T >, ::utility::pointer::shared_ptr< ::numeric::xyzTransform< T > > > xyzTransform_typename_exposer_type;

    xyzTransform_typename_exposer_type xyzTransform_typename_exposer( std::string("xyzTransform_" + type_name).c_str()  , "numeric/xyzTransform.fwd.hh:29", boost::python::init <  >() );
    xyzTransform_typename_exposer.def( boost::python::init< ::numeric::xyzTransform< T > const & > ( (boost::python::arg("")) , "numeric/xyzTransform.fwd.hh:29" ) );
    xyzTransform_typename_exposer.def( boost::python::init< ::numeric::xyzMatrix< T > const & > ( (boost::python::arg("rin")) , "numeric/xyzTransform.hh:55" ) );
    xyzTransform_typename_exposer.def( boost::python::init< ::numeric::xyzVector< T > const & > ( (boost::python::arg("tin")) , "numeric/xyzTransform.hh:56" ) );
    xyzTransform_typename_exposer.def( boost::python::init< ::numeric::xyzMatrix< T > const &, ::numeric::xyzVector< T > const & > ( (boost::python::arg("rin"), boost::python::arg("tin")) , "numeric/xyzTransform.hh:57" ) );
    xyzTransform_typename_exposer.def( boost::python::init< ::utility::fixedsizearray1<T,6ul> const & > ( (boost::python::arg("_rt6")) , "numeric/xyzTransform.hh:58" ) );
    xyzTransform_typename_exposer.def( boost::python::init< ::numeric::xyzVector< T > const &, ::numeric::xyzVector< T > const &, ::numeric::xyzVector< T > const & > ( (boost::python::arg("u"), boost::python::arg("v"), boost::python::arg("w")) , "numeric/xyzTransform.hh:71" ) );
    xyzTransform_typename_exposer.def( boost::python::init< ::numeric::xyzVector< T > const &, ::numeric::xyzVector< T > const &, ::numeric::xyzVector< T > const &, ::numeric::xyzVector< T > const & > ( (boost::python::arg("c"), boost::python::arg("u"), boost::python::arg("v"), boost::python::arg("w")) , "numeric/xyzTransform.hh:72" ) );
  
    { // ::numeric::xyzTransform< T >::from_four_points
      typedef ::numeric::xyzTransform< T > & ( ::numeric::xyzTransform< T >:: * from_four_points_function_type)(::numeric::xyzVector< T > const & __a0, ::numeric::xyzVector< T > const & __a1, ::numeric::xyzVector< T > const & __a2, ::numeric::xyzVector< T > const & __a3) ;
    
      xyzTransform_typename_exposer.def("from_four_points"
        , from_four_points_function_type( &::numeric::xyzTransform< T >::from_four_points )
        , ( boost::python::arg("c"), boost::python::arg("u"), boost::python::arg("v"), boost::python::arg("w") )
        , boost::python::return_value_policy< boost::python::reference_existing_object >()
        , "numeric/xyzTransform.hh:61" );
    }
    
  
    { // ::numeric::xyzTransform< T >::identity
      typedef ::numeric::xyzTransform< T > ( * identity_function_type)();
    
      xyzTransform_typename_exposer.def("identity"
        , identity_function_type( &::numeric::xyzTransform< T >::identity )
        , "numeric/xyzTransform.hh:74" );
    }
    
  
    { // ::numeric::xyzTransform< T >::BAD_XFORM
      typedef ::numeric::xyzTransform< T > ( * BAD_XFORM_function_type)();
    
      xyzTransform_typename_exposer.def("BAD_XFORM"
        , BAD_XFORM_function_type( &::numeric::xyzTransform< T >::BAD_XFORM )
        , "numeric/xyzTransform.hh:75" );
    }
    
  
    { // ::numeric::xyzTransform< T >::BAD_RT6
      typedef ::utility::fixedsizearray1<T,6ul> ( * BAD_RT6_function_type)();
    
      xyzTransform_typename_exposer.def("BAD_RT6"
        , BAD_RT6_function_type( &::numeric::xyzTransform< T >::BAD_RT6 )
        , "numeric/xyzTransform.hh:76" );
    }
    
  
    { // ::numeric::xyzTransform< T >::xx const
      typedef T const & ( ::numeric::xyzTransform< T >:: * xx_function_type)()  const;
    
      xyzTransform_typename_exposer.def("xx"
        , xx_function_type( &::numeric::xyzTransform< T >::xx )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:78" );
    }
    
  
    { // ::numeric::xyzTransform< T >::xy const
      typedef T const & ( ::numeric::xyzTransform< T >:: * xy_function_type)()  const;
    
      xyzTransform_typename_exposer.def("xy"
        , xy_function_type( &::numeric::xyzTransform< T >::xy )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:79" );
    }
    
  
    { // ::numeric::xyzTransform< T >::xz const
      typedef T const & ( ::numeric::xyzTransform< T >:: * xz_function_type)()  const;
    
      xyzTransform_typename_exposer.def("xz"
        , xz_function_type( &::numeric::xyzTransform< T >::xz )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:80" );
    }
    
  
    { // ::numeric::xyzTransform< T >::yx const
      typedef T const & ( ::numeric::xyzTransform< T >:: * yx_function_type)()  const;
    
      xyzTransform_typename_exposer.def("yx"
        , yx_function_type( &::numeric::xyzTransform< T >::yx )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:81" );
    }
    
  
    { // ::numeric::xyzTransform< T >::yy const
      typedef T const & ( ::numeric::xyzTransform< T >:: * yy_function_type)()  const;
    
      xyzTransform_typename_exposer.def("yy"
        , yy_function_type( &::numeric::xyzTransform< T >::yy )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:82" );
    }
    
  
    { // ::numeric::xyzTransform< T >::yz const
      typedef T const & ( ::numeric::xyzTransform< T >:: * yz_function_type)()  const;
    
      xyzTransform_typename_exposer.def("yz"
        , yz_function_type( &::numeric::xyzTransform< T >::yz )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:83" );
    }
    
  
    { // ::numeric::xyzTransform< T >::zx const
      typedef T const & ( ::numeric::xyzTransform< T >:: * zx_function_type)()  const;
    
      xyzTransform_typename_exposer.def("zx"
        , zx_function_type( &::numeric::xyzTransform< T >::zx )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:84" );
    }
    
  
    { // ::numeric::xyzTransform< T >::zy const
      typedef T const & ( ::numeric::xyzTransform< T >:: * zy_function_type)()  const;
    
      xyzTransform_typename_exposer.def("zy"
        , zy_function_type( &::numeric::xyzTransform< T >::zy )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:85" );
    }
    
  
    { // ::numeric::xyzTransform< T >::zz const
      typedef T const & ( ::numeric::xyzTransform< T >:: * zz_function_type)()  const;
    
      xyzTransform_typename_exposer.def("zz"
        , zz_function_type( &::numeric::xyzTransform< T >::zz )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:86" );
    }
    
  
    { // ::numeric::xyzTransform< T >::px const
      typedef T const & ( ::numeric::xyzTransform< T >:: * px_function_type)()  const;
    
      xyzTransform_typename_exposer.def("px"
        , px_function_type( &::numeric::xyzTransform< T >::px )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:87" );
    }
    
  
    { // ::numeric::xyzTransform< T >::py const
      typedef T const & ( ::numeric::xyzTransform< T >:: * py_function_type)()  const;
    
      xyzTransform_typename_exposer.def("py"
        , py_function_type( &::numeric::xyzTransform< T >::py )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:88" );
    }
    
  
    { // ::numeric::xyzTransform< T >::pz const
      typedef T const & ( ::numeric::xyzTransform< T >:: * pz_function_type)()  const;
    
      xyzTransform_typename_exposer.def("pz"
        , pz_function_type( &::numeric::xyzTransform< T >::pz )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:89" );
    }
    
  
    { // ::numeric::xyzTransform< T >::x const
      typedef T const & ( ::numeric::xyzTransform< T >:: * x_function_type)()  const;
    
      xyzTransform_typename_exposer.def("x"
        , x_function_type( &::numeric::xyzTransform< T >::x )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:90" );
    }
    
  
    { // ::numeric::xyzTransform< T >::y const
      typedef T const & ( ::numeric::xyzTransform< T >:: * y_function_type)()  const;
    
      xyzTransform_typename_exposer.def("y"
        , y_function_type( &::numeric::xyzTransform< T >::y )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:91" );
    }
    
  
    { // ::numeric::xyzTransform< T >::z const
      typedef T const & ( ::numeric::xyzTransform< T >:: * z_function_type)()  const;
    
      xyzTransform_typename_exposer.def("z"
        , z_function_type( &::numeric::xyzTransform< T >::z )
        , boost::python::return_value_policy< boost::python::copy_const_reference >()
        , "numeric/xyzTransform.hh:92" );
    }
    
  
    { // ::numeric::xyzTransform< T >::xx
      typedef T & ( ::numeric::xyzTransform< T >:: * xx_function_type)() ;
    
      xyzTransform_typename_exposer.def("xx"
        , xx_function_type( &::numeric::xyzTransform< T >::xx )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:93" );
    }
    
  
    { // ::numeric::xyzTransform< T >::xy
      typedef T & ( ::numeric::xyzTransform< T >:: * xy_function_type)() ;
    
      xyzTransform_typename_exposer.def("xy"
        , xy_function_type( &::numeric::xyzTransform< T >::xy )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:94" );
    }
    
  
    { // ::numeric::xyzTransform< T >::xz
      typedef T & ( ::numeric::xyzTransform< T >:: * xz_function_type)() ;
    
      xyzTransform_typename_exposer.def("xz"
        , xz_function_type( &::numeric::xyzTransform< T >::xz )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:95" );
    }
    
  
    { // ::numeric::xyzTransform< T >::yx
      typedef T & ( ::numeric::xyzTransform< T >:: * yx_function_type)() ;
    
      xyzTransform_typename_exposer.def("yx"
        , yx_function_type( &::numeric::xyzTransform< T >::yx )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:96" );
    }
    
  
    { // ::numeric::xyzTransform< T >::yy
      typedef T & ( ::numeric::xyzTransform< T >:: * yy_function_type)() ;
    
      xyzTransform_typename_exposer.def("yy"
        , yy_function_type( &::numeric::xyzTransform< T >::yy )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:97" );
    }
    
  
    { // ::numeric::xyzTransform< T >::yz
      typedef T & ( ::numeric::xyzTransform< T >:: * yz_function_type)() ;
    
      xyzTransform_typename_exposer.def("yz"
        , yz_function_type( &::numeric::xyzTransform< T >::yz )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:98" );
    }
    
  
    { // ::numeric::xyzTransform< T >::zx
      typedef T & ( ::numeric::xyzTransform< T >:: * zx_function_type)() ;
    
      xyzTransform_typename_exposer.def("zx"
        , zx_function_type( &::numeric::xyzTransform< T >::zx )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:99" );
    }
    
  
    { // ::numeric::xyzTransform< T >::zy
      typedef T & ( ::numeric::xyzTransform< T >:: * zy_function_type)() ;
    
      xyzTransform_typename_exposer.def("zy"
        , zy_function_type( &::numeric::xyzTransform< T >::zy )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:100" );
    }
    
  
    { // ::numeric::xyzTransform< T >::zz
      typedef T & ( ::numeric::xyzTransform< T >:: * zz_function_type)() ;
    
      xyzTransform_typename_exposer.def("zz"
        , zz_function_type( &::numeric::xyzTransform< T >::zz )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:101" );
    }
    
  
    { // ::numeric::xyzTransform< T >::px
      typedef T & ( ::numeric::xyzTransform< T >:: * px_function_type)() ;
    
      xyzTransform_typename_exposer.def("px"
        , px_function_type( &::numeric::xyzTransform< T >::px )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:102" );
    }
    
  
    { // ::numeric::xyzTransform< T >::py
      typedef T & ( ::numeric::xyzTransform< T >:: * py_function_type)() ;
    
      xyzTransform_typename_exposer.def("py"
        , py_function_type( &::numeric::xyzTransform< T >::py )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:103" );
    }
    
  
    { // ::numeric::xyzTransform< T >::pz
      typedef T & ( ::numeric::xyzTransform< T >:: * pz_function_type)() ;
    
      xyzTransform_typename_exposer.def("pz"
        , pz_function_type( &::numeric::xyzTransform< T >::pz )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:104" );
    }
    
  
    { // ::numeric::xyzTransform< T >::x
      typedef T & ( ::numeric::xyzTransform< T >:: * x_function_type)() ;
    
      xyzTransform_typename_exposer.def("x"
        , x_function_type( &::numeric::xyzTransform< T >::x )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:105" );
    }
    
  
    { // ::numeric::xyzTransform< T >::y
      typedef T & ( ::numeric::xyzTransform< T >:: * y_function_type)() ;
    
      xyzTransform_typename_exposer.def("y"
        , y_function_type( &::numeric::xyzTransform< T >::y )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:106" );
    }
    
  
    { // ::numeric::xyzTransform< T >::z
      typedef T & ( ::numeric::xyzTransform< T >:: * z_function_type)() ;
    
      xyzTransform_typename_exposer.def("z"
        , z_function_type( &::numeric::xyzTransform< T >::z )
        , boost::python::return_value_policy< boost::python::copy_non_const_reference>()
        , "numeric/xyzTransform.hh:107" );
    }
    
  
    { // ::numeric::xyzTransform< T >::inverse const
      typedef ::numeric::xyzTransform< T > ( ::numeric::xyzTransform< T >:: * inverse_function_type)()  const;
    
      xyzTransform_typename_exposer.def("inverse"
        , inverse_function_type( &::numeric::xyzTransform< T >::inverse )
        , "numeric/xyzTransform.hh:110" );
    }
    
  
    { // ::numeric::xyzTransform< T >::distance const
      typedef T ( ::numeric::xyzTransform< T >:: * distance_function_type)(::numeric::xyzTransform< T > const & __a0)  const;
    
      xyzTransform_typename_exposer.def("distance"
        , distance_function_type( &::numeric::xyzTransform< T >::distance )
        , ( boost::python::arg("b") )
        , "numeric/xyzTransform.hh:112" );
    }
    
  
    { // ::numeric::xyzTransform< T >::distance_squared const
      typedef T ( ::numeric::xyzTransform< T >:: * distance_squared_function_type)(::numeric::xyzTransform< T > const & __a0)  const;
    
      xyzTransform_typename_exposer.def("distance_squared"
        , distance_squared_function_type( &::numeric::xyzTransform< T >::distance_squared )
        , ( boost::python::arg("b") )
        , "numeric/xyzTransform.hh:113" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rot
      typedef ::numeric::xyzTransform< T > ( * rot_function_type)(::numeric::xyzMatrix< T > const & __a0, ::numeric::xyzVector< T > const & __a1, ::numeric::xyzVector< T > const & __a2);
    
      xyzTransform_typename_exposer.def("rot"
        , rot_function_type( &::numeric::xyzTransform< T >::rot )
        , ( boost::python::arg("rot"), boost::python::arg("o_cen"), boost::python::arg("cen") )
        , "numeric/xyzTransform.hh:131" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rot
      typedef ::numeric::xyzTransform< T > ( * rot_function_type)(::numeric::xyzMatrix< T > const & __a0, ::numeric::xyzVector< T > const & __a1);
    
      xyzTransform_typename_exposer.def("rot"
        , rot_function_type( &::numeric::xyzTransform< T >::rot )
        , ( boost::python::arg("rot"), boost::python::arg("cen")=(numeric::xyzVector<T>(0, 0, 0)) )
        , "numeric/xyzTransform.hh:132" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rot
      typedef ::numeric::xyzTransform< T > ( * rot_function_type)(::numeric::xyzVector< T > const & __a0, T const & __a1, ::numeric::xyzVector< T > const & __a2);
    
      xyzTransform_typename_exposer.def("rot"
        , rot_function_type( &::numeric::xyzTransform< T >::rot )
        , ( boost::python::arg("axs"), boost::python::arg("ang"), boost::python::arg("cen")=(numeric::xyzVector<T>(0, 0, 0)) )
        , "numeric/xyzTransform.hh:133" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rot_deg
      typedef ::numeric::xyzTransform< T > ( * rot_deg_function_type)(::numeric::xyzVector< T > const & __a0, T const & __a1, ::numeric::xyzVector< T > const & __a2);
    
      xyzTransform_typename_exposer.def("rot_deg"
        , rot_deg_function_type( &::numeric::xyzTransform< T >::rot_deg )
        , ( boost::python::arg("axs"), boost::python::arg("ang"), boost::python::arg("cen")=(numeric::xyzVector<T>(0, 0, 0)) )
        , "numeric/xyzTransform.hh:134" );
    }
    
  
    { // ::numeric::xyzTransform< T >::align
      typedef ::numeric::xyzTransform< T > ( * align_function_type)(::numeric::xyzVector< T > const & __a0, ::numeric::xyzVector< T > const & __a1);
    
      xyzTransform_typename_exposer.def("align"
        , align_function_type( &::numeric::xyzTransform< T >::align )
        , ( boost::python::arg("to"), boost::python::arg("from")=(numeric::xyzVector<T>(1, 0, 0)) )
        , "numeric/xyzTransform.hh:136" );
    }
    
  
    { // ::numeric::xyzTransform< T >::align_fast
      typedef ::numeric::xyzTransform< T > ( * align_fast_function_type)(::numeric::xyzVector< T > const & __a0, ::numeric::xyzVector< T > const & __a1);
    
      xyzTransform_typename_exposer.def("align_fast"
        , align_fast_function_type( &::numeric::xyzTransform< T >::align_fast )
        , ( boost::python::arg("to"), boost::python::arg("from")=(numeric::xyzVector<T>(1.0e+0, 0.0, 0.0)) )
        , "numeric/xyzTransform.hh:141" );
    }
    
  
    { // ::numeric::xyzTransform< T >::xform const
      typedef ::numeric::xyzVector< T > ( ::numeric::xyzTransform< T >:: * xform_function_type)(::numeric::xyzVector< T > const & __a0)  const;
    
      xyzTransform_typename_exposer.def("xform"
        , xform_function_type( &::numeric::xyzTransform< T >::xform )
        , ( boost::python::arg("v") )
        , "numeric/xyzTransform.hh:145" );
    }
    
  
    { // ::numeric::xyzTransform< T >::inv_xform const
      typedef ::numeric::xyzVector< T > ( ::numeric::xyzTransform< T >:: * inv_xform_function_type)(::numeric::xyzVector< T > const & __a0)  const;
    
      xyzTransform_typename_exposer.def("inv_xform"
        , inv_xform_function_type( &::numeric::xyzTransform< T >::inv_xform )
        , ( boost::python::arg("v") )
        , "numeric/xyzTransform.hh:146" );
    }
    
  
    { // ::numeric::xyzTransform< T >::to_quaternion const
      typedef void ( ::numeric::xyzTransform< T >:: * to_quaternion_function_type)(T & __a0, T & __a1, T & __a2, T & __a3)  const;
    
      xyzTransform_typename_exposer.def("to_quaternion"
        , to_quaternion_function_type( &::numeric::xyzTransform< T >::to_quaternion )
        , ( boost::python::arg("qw"), boost::python::arg("qx"), boost::python::arg("qy"), boost::python::arg("qz") )
        , "numeric/xyzTransform.hh:149" );
    }
    
  
    { // ::numeric::xyzTransform< T >::from_quaternion
      typedef void ( ::numeric::xyzTransform< T >:: * from_quaternion_function_type)(T const & __a0, T const & __a1, T const & __a2, T const & __a3) ;
    
      xyzTransform_typename_exposer.def("from_quaternion"
        , from_quaternion_function_type( &::numeric::xyzTransform< T >::from_quaternion )
        , ( boost::python::arg("qw"), boost::python::arg("qx"), boost::python::arg("qy"), boost::python::arg("qz") )
        , "numeric/xyzTransform.hh:168" );
    }
    
  
    { // ::numeric::xyzTransform< T >::euler_angles_rad const
      typedef ::numeric::xyzVector< T > ( ::numeric::xyzTransform< T >:: * euler_angles_rad_function_type)()  const;
    
      xyzTransform_typename_exposer.def("euler_angles_rad"
        , euler_angles_rad_function_type( &::numeric::xyzTransform< T >::euler_angles_rad )
        , "see numeric/HomogeneousTransform\n" );
    }
    
  
    { // ::numeric::xyzTransform< T >::euler_angles_deg const
      typedef ::numeric::xyzVector< T > ( ::numeric::xyzTransform< T >:: * euler_angles_deg_function_type)()  const;
    
      xyzTransform_typename_exposer.def("euler_angles_deg"
        , euler_angles_deg_function_type( &::numeric::xyzTransform< T >::euler_angles_deg )
        , "numeric/xyzTransform.hh:250" );
    }
    
  
    { // ::numeric::xyzTransform< T >::from_euler_angles_rad
      typedef ::numeric::xyzTransform< T > & ( ::numeric::xyzTransform< T >:: * from_euler_angles_rad_function_type)(T const & __a0, T const & __a1, T const & __a2) ;
    
      xyzTransform_typename_exposer.def("from_euler_angles_rad"
        , from_euler_angles_rad_function_type( &::numeric::xyzTransform< T >::from_euler_angles_rad )
        , ( boost::python::arg("phi"), boost::python::arg("psi"), boost::python::arg("theta") )
        , boost::python::return_value_policy< boost::python::reference_existing_object >()
        , "numeric/xyzTransform.hh:257" );
    }
    
  
    { // ::numeric::xyzTransform< T >::from_euler_angles_rad
      typedef ::numeric::xyzTransform< T > & ( ::numeric::xyzTransform< T >:: * from_euler_angles_rad_function_type)(::numeric::xyzVector< T > const & __a0) ;
    
      xyzTransform_typename_exposer.def("from_euler_angles_rad"
        , from_euler_angles_rad_function_type( &::numeric::xyzTransform< T >::from_euler_angles_rad )
        , ( boost::python::arg("euler") )
        , boost::python::return_value_policy< boost::python::reference_existing_object >()
        , "numeric/xyzTransform.hh:267" );
    }
    
  
    { // ::numeric::xyzTransform< T >::from_euler_angles_deg
      typedef ::numeric::xyzTransform< T > & ( ::numeric::xyzTransform< T >:: * from_euler_angles_deg_function_type)(T const & __a0, T const & __a1, T const & __a2) ;
    
      xyzTransform_typename_exposer.def("from_euler_angles_deg"
        , from_euler_angles_deg_function_type( &::numeric::xyzTransform< T >::from_euler_angles_deg )
        , ( boost::python::arg("phi"), boost::python::arg("psi"), boost::python::arg("theta") )
        , boost::python::return_value_policy< boost::python::reference_existing_object >()
        , "numeric/xyzTransform.hh:270" );
    }
    
  
    { // ::numeric::xyzTransform< T >::from_euler_angles_deg
      typedef ::numeric::xyzTransform< T > & ( ::numeric::xyzTransform< T >:: * from_euler_angles_deg_function_type)(::numeric::xyzVector< T > const & __a0) ;
    
      xyzTransform_typename_exposer.def("from_euler_angles_deg"
        , from_euler_angles_deg_function_type( &::numeric::xyzTransform< T >::from_euler_angles_deg )
        , ( boost::python::arg("euler") )
        , boost::python::return_value_policy< boost::python::reference_existing_object >()
        , "numeric/xyzTransform.hh:274" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rt6 const
      typedef ::utility::fixedsizearray1<T,6ul> ( ::numeric::xyzTransform< T >:: * rt6_function_type)()  const;
    
      xyzTransform_typename_exposer.def("rt6"
        , rt6_function_type( &::numeric::xyzTransform< T >::rt6 )
        , "numeric/xyzTransform.hh:280" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rt6
      typedef ::numeric::xyzTransform< T > & ( ::numeric::xyzTransform< T >:: * rt6_function_type)(::utility::fixedsizearray1<T,6ul> const & __a0) ;
    
      xyzTransform_typename_exposer.def("rt6"
        , rt6_function_type( &::numeric::xyzTransform< T >::rt6 )
        , ( boost::python::arg("rt6") )
        , boost::python::return_value_policy< boost::python::reference_existing_object >()
        , "numeric/xyzTransform.hh:306" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rt6
      typedef ::numeric::xyzTransform< T > & ( ::numeric::xyzTransform< T >:: * rt6_function_type)(T const & __a0, T const & __a1, T const & __a2, T const & __a3, T const & __a4, T const & __a5) ;
    
      xyzTransform_typename_exposer.def("rt6"
        , rt6_function_type( &::numeric::xyzTransform< T >::rt6 )
        , ( boost::python::arg("i"), boost::python::arg("j"), boost::python::arg("k"), boost::python::arg("l"), boost::python::arg("m"), boost::python::arg("n") )
        , boost::python::return_value_policy< boost::python::reference_existing_object >()
        , "numeric/xyzTransform.hh:315" );
    }
    
  
    { // ::numeric::xyzTransform< T >::hash64 const
      typedef ::uint64_t ( ::numeric::xyzTransform< T >:: * hash64_function_type)(T const & __a0, T const & __a1)  const;
    
      xyzTransform_typename_exposer.def("hash64"
        , hash64_function_type( &::numeric::xyzTransform< T >::hash64 )
        , ( boost::python::arg("cw")=(1.000000000000000055511151231257827021181583404541015625e-1), boost::python::arg("aw")=((3.6e+2 / 1.024e+3)) )
        , "numeric/xyzTransform.hh:325" );
    }
    
  
    { // ::numeric::xyzTransform< T >::symhash64 const
      typedef ::uint64_t ( ::numeric::xyzTransform< T >:: * symhash64_function_type)(T const & __a0, T const & __a1)  const;
    
      xyzTransform_typename_exposer.def("symhash64"
        , symhash64_function_type( &::numeric::xyzTransform< T >::symhash64 )
        , ( boost::python::arg("cw")=(1.000000000000000055511151231257827021181583404541015625e-1), boost::python::arg("aw")=((3.6e+2 / 1.024e+3)) )
        , "numeric/xyzTransform.hh:341" );
    }
    
  
    { // ::numeric::xyzTransform< T >::intersect3D_2Planes const
      typedef int ( ::numeric::xyzTransform< T >:: * intersect3D_2Planes_function_type)(typename ::numeric::xyzTransform< T >::Plane __a0, typename ::numeric::xyzTransform< T >::Plane __a1, typename ::numeric::xyzTransform< T >::Line * __a2)  const;
    
      xyzTransform_typename_exposer.def("intersect3D_2Planes"
        , intersect3D_2Planes_function_type( &::numeric::xyzTransform< T >::intersect3D_2Planes )
        , ( boost::python::arg("Pn1"), boost::python::arg("Pn2"), boost::python::arg("L") )
        , "numeric/xyzTransform.hh:364" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rotation_axis const
      typedef void ( ::numeric::xyzTransform< T >:: * rotation_axis_function_type)(::numeric::xyzVector< T > & __a0, ::numeric::xyzVector< T > & __a1, T & __a2)  const;
    
      xyzTransform_typename_exposer.def("rotation_axis"
        , rotation_axis_function_type( &::numeric::xyzTransform< T >::rotation_axis )
        , ( boost::python::arg("axis"), boost::python::arg("cen"), boost::python::arg("angle") )
        , "numeric/xyzTransform.hh:420" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rotation_axis const
      typedef ::numeric::xyzVector< T > ( ::numeric::xyzTransform< T >:: * rotation_axis_function_type)()  const;
    
      xyzTransform_typename_exposer.def("rotation_axis"
        , rotation_axis_function_type( &::numeric::xyzTransform< T >::rotation_axis )
        , "numeric/xyzTransform.hh:454" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rotation_angle_degrees const
      typedef T ( ::numeric::xyzTransform< T >:: * rotation_angle_degrees_function_type)()  const;
    
      xyzTransform_typename_exposer.def("rotation_angle_degrees"
        , rotation_angle_degrees_function_type( &::numeric::xyzTransform< T >::rotation_angle_degrees )
        , "numeric/xyzTransform.hh:459" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rotation_angle const
      typedef T ( ::numeric::xyzTransform< T >:: * rotation_angle_function_type)()  const;
    
      xyzTransform_typename_exposer.def("rotation_angle"
        , rotation_angle_function_type( &::numeric::xyzTransform< T >::rotation_angle )
        , "numeric/xyzTransform.hh:464" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rotation_cosine const
      typedef T ( ::numeric::xyzTransform< T >:: * rotation_cosine_function_type)()  const;
    
      xyzTransform_typename_exposer.def("rotation_cosine"
        , rotation_cosine_function_type( &::numeric::xyzTransform< T >::rotation_cosine )
        , "numeric/xyzTransform.hh:470" );
    }
    
  
    { // ::numeric::xyzTransform< T >::rotation_sine const
      typedef T ( ::numeric::xyzTransform< T >:: * rotation_sine_function_type)()  const;
    
      xyzTransform_typename_exposer.def("rotation_sine"
        , rotation_sine_function_type( &::numeric::xyzTransform< T >::rotation_sine )
        , "numeric/xyzTransform.hh:474" );
    }
    
  
    { // ::numeric::xyzTransform< T >::approx_lever_distance const
      typedef T ( ::numeric::xyzTransform< T >:: * approx_lever_distance_function_type)(::numeric::xyzTransform< T > const & __a0, T const & __a1)  const;
    
      xyzTransform_typename_exposer.def("approx_lever_distance"
        , approx_lever_distance_function_type( &::numeric::xyzTransform< T >::approx_lever_distance )
        , ( boost::python::arg("o"), boost::python::arg("lever")=(1.0e+0) )
        , "numeric/xyzTransform.hh:479" );
    }
    
  
    { // ::numeric::xyzTransform< T >::bad const
      typedef bool ( ::numeric::xyzTransform< T >:: * bad_function_type)()  const;
    
      xyzTransform_typename_exposer.def("bad"
        , bad_function_type( &::numeric::xyzTransform< T >::bad )
        , "numeric/xyzTransform.hh:486" );
    }
    
  
    { // ::numeric::xyzTransform< T >::badfast const
      typedef bool ( ::numeric::xyzTransform< T >:: * badfast_function_type)()  const;
    
      xyzTransform_typename_exposer.def("badfast"
        , badfast_function_type( &::numeric::xyzTransform< T >::badfast )
        , "numeric/xyzTransform.hh:492" );
    }
    
  
    xyzTransform_typename_exposer.staticmethod("rot_deg");
  
  
    xyzTransform_typename_exposer.staticmethod("align");
  
  
    xyzTransform_typename_exposer.staticmethod("BAD_XFORM");
  
  
    xyzTransform_typename_exposer.staticmethod("BAD_RT6");
  
  
    xyzTransform_typename_exposer.staticmethod("align_fast");
  
  
    xyzTransform_typename_exposer.staticmethod("rot");
  
  
    xyzTransform_typename_exposer.staticmethod("identity");
  
  }
}

template<class T>
void instantiate_numeric_containers(std::string type_name)
{
  { //::numeric::xyzMatrix< T >
      typedef bp::class_< numeric::xyzMatrix< T > > xyzMatrix_typename_exposer_t;
      xyzMatrix_typename_exposer_t xyzMatrix_typename_exposer = xyzMatrix_typename_exposer_t( std::string("xyzMatrix_" + type_name).c_str()  ); // "numeric___xyzMatrix_ T "
      bp::scope xyzMatrix_typename_scope( xyzMatrix_typename_exposer );
      xyzMatrix_typename_exposer.def( bp::init< >() );
      xyzMatrix_typename_exposer.def( bp::init<  T  const & >(( bp::arg("t") )) );
      { //::numeric::xyzMatrix< T >::I

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > const & ( *I_function_type )(  );

          xyzMatrix_typename_exposer.def(
              "I"
              , I_function_type( &::numeric::xyzMatrix< T >::I )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::add_diagonal

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*add_diagonal_function_type )(  T  const &, T  const &, T  const & ) ;

          xyzMatrix_typename_exposer.def(
              "add_diagonal"
              , add_diagonal_function_type( &::numeric::xyzMatrix< T >::add_diagonal )
              , ( bp::arg("xx_a"), bp::arg("yy_a"), bp::arg("zz_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::clear

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*clear_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "clear"
              , clear_function_type( &::numeric::xyzMatrix< T >::clear )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::col

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*col_function_type )( int const ) const;

          xyzMatrix_typename_exposer.def(
              "col"
              , col_function_type( &::numeric::xyzMatrix< T >::col )
              , ( bp::arg("i") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::col

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*col_function_type )( int const,::numeric::xyzVector< T > const & ) ;

          xyzMatrix_typename_exposer.def(
              "col"
              , col_function_type( &::numeric::xyzMatrix< T >::col )
              , ( bp::arg("i"), bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::col_x

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*col_x_function_type )(  ) const;

          xyzMatrix_typename_exposer.def(
              "col_x"
              , col_x_function_type( &::numeric::xyzMatrix< T >::col_x )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::col_x

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*col_x_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzMatrix_typename_exposer.def(
              "col_x"
              , col_x_function_type( &::numeric::xyzMatrix< T >::col_x )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::col_y

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*col_y_function_type )(  ) const;

          xyzMatrix_typename_exposer.def(
              "col_y"
              , col_y_function_type( &::numeric::xyzMatrix< T >::col_y )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::col_y

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*col_y_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzMatrix_typename_exposer.def(
              "col_y"
              , col_y_function_type( &::numeric::xyzMatrix< T >::col_y )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::col_z

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*col_z_function_type )(  ) const;

          xyzMatrix_typename_exposer.def(
              "col_z"
              , col_z_function_type( &::numeric::xyzMatrix< T >::col_z )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::col_z

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*col_z_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzMatrix_typename_exposer.def(
              "col_z"
              , col_z_function_type( &::numeric::xyzMatrix< T >::col_z )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::cols

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > ( *cols_function_type )(  T  const &, T  const &, T  const &, T  const &, T  const &, T  const &, T  const &, T  const &, T  const & );

          xyzMatrix_typename_exposer.def(
              "cols"
              , cols_function_type( &::numeric::xyzMatrix< T >::cols )
              , ( bp::arg("xx_a"), bp::arg("yx_a"), bp::arg("zx_a"), bp::arg("xy_a"), bp::arg("yy_a"), bp::arg("zy_a"), bp::arg("xz_a"), bp::arg("yz_a"), bp::arg("zz_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::diag

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > ( *diag_function_type )(  T  const &, T  const &, T  const & );

          xyzMatrix_typename_exposer.def(
              "diag"
              , diag_function_type( &::numeric::xyzMatrix< T >::diag )
              , ( bp::arg("xx_a"), bp::arg("yy_a"), bp::arg("zz_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::identity

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > ( *identity_function_type )(  );

          xyzMatrix_typename_exposer.def(
              "identity"
              , identity_function_type( &::numeric::xyzMatrix< T >::identity )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      xyzMatrix_typename_exposer.def( bp::self *= bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::self += bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::self -= bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::self /= bp::other<  T  >() );
      { //::numeric::xyzMatrix< T >::row

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*row_function_type )( int const ) const;

          xyzMatrix_typename_exposer.def(
              "row"
              , row_function_type( &::numeric::xyzMatrix< T >::row )
              , ( bp::arg("i") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::row

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*row_function_type )( int const,::numeric::xyzVector< T > const & ) ;

          xyzMatrix_typename_exposer.def(
              "row"
              , row_function_type( &::numeric::xyzMatrix< T >::row )
              , ( bp::arg("i"), bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::row_x

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*row_x_function_type )(  ) const;

          xyzMatrix_typename_exposer.def(
              "row_x"
              , row_x_function_type( &::numeric::xyzMatrix< T >::row_x )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::row_x

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*row_x_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzMatrix_typename_exposer.def(
              "row_x"
              , row_x_function_type( &::numeric::xyzMatrix< T >::row_x )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::row_y

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*row_y_function_type )(  ) const;

          xyzMatrix_typename_exposer.def(
              "row_y"
              , row_y_function_type( &::numeric::xyzMatrix< T >::row_y )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::row_y

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*row_y_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzMatrix_typename_exposer.def(
              "row_y"
              , row_y_function_type( &::numeric::xyzMatrix< T >::row_y )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::row_z

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*row_z_function_type )(  ) const;

          xyzMatrix_typename_exposer.def(
              "row_z"
              , row_z_function_type( &::numeric::xyzMatrix< T >::row_z )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::row_z

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*row_z_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzMatrix_typename_exposer.def(
              "row_z"
              , row_z_function_type( &::numeric::xyzMatrix< T >::row_z )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::rows

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > ( *rows_function_type )(  T  const &, T  const &, T  const &, T  const &, T  const &, T  const &, T  const &, T  const &, T  const & );

          xyzMatrix_typename_exposer.def(
              "rows"
              , rows_function_type( &::numeric::xyzMatrix< T >::rows )
              , ( bp::arg("xx_a"), bp::arg("xy_a"), bp::arg("xz_a"), bp::arg("yx_a"), bp::arg("yy_a"), bp::arg("yz_a"), bp::arg("zx_a"), bp::arg("zy_a"), bp::arg("zz_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::set_diagonal

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*set_diagonal_function_type )(  T  const &, T  const &, T  const & ) ;

          xyzMatrix_typename_exposer.def(
              "set_diagonal"
              , set_diagonal_function_type( &::numeric::xyzMatrix< T >::set_diagonal )
              , ( bp::arg("xx_a"), bp::arg("yy_a"), bp::arg("zz_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::subtract_diagonal

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*subtract_diagonal_function_type )(  T  const &, T  const &, T  const & ) ;

          xyzMatrix_typename_exposer.def(
              "subtract_diagonal"
              , subtract_diagonal_function_type( &::numeric::xyzMatrix< T >::subtract_diagonal )
              , ( bp::arg("xx_a"), bp::arg("yy_a"), bp::arg("zz_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::to_diag

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*to_diag_function_type )(  T  const &, T  const &, T  const & ) ;

          xyzMatrix_typename_exposer.def(
              "to_diag"
              , to_diag_function_type( &::numeric::xyzMatrix< T >::to_diag )
              , ( bp::arg("xx_a"), bp::arg("yy_a"), bp::arg("zz_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::to_identity

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*to_identity_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "to_identity"
              , to_identity_function_type( &::numeric::xyzMatrix< T >::to_identity )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::transpose

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*transpose_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "transpose"
              , transpose_function_type( &::numeric::xyzMatrix< T >::transpose )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::inverse

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*inverse_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "inverse"
              , inverse_function_type( &::numeric::xyzMatrix< T >::inverse )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::xx

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef  T  & ( exported_class_t::*xx_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "xx"
              , xx_function_type( &::numeric::xyzMatrix< T >::xx )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::xy

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef  T  & ( exported_class_t::*xy_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "xy"
              , xy_function_type( &::numeric::xyzMatrix< T >::xy )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::xz

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef  T  & ( exported_class_t::*xz_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "xz"
              , xz_function_type( &::numeric::xyzMatrix< T >::xz )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::yx

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef  T  & ( exported_class_t::*yx_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "yx"
              , yx_function_type( &::numeric::xyzMatrix< T >::yx )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::yy

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef  T  & ( exported_class_t::*yy_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "yy"
              , yy_function_type( &::numeric::xyzMatrix< T >::yy )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::yz

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef  T  & ( exported_class_t::*yz_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "yz"
              , yz_function_type( &::numeric::xyzMatrix< T >::yz )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::zero

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef ::numeric::xyzMatrix< T > & ( exported_class_t::*zero_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "zero"
              , zero_function_type( &::numeric::xyzMatrix< T >::zero )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::zx

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef  T  & ( exported_class_t::*zx_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "zx"
              , zx_function_type( &::numeric::xyzMatrix< T >::zx )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::zy

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef  T  & ( exported_class_t::*zy_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "zy"
              , zy_function_type( &::numeric::xyzMatrix< T >::zy )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzMatrix< T >::zz

          typedef numeric::xyzMatrix< T > exported_class_t;
          typedef  T  & ( exported_class_t::*zz_function_type )(  ) ;

          xyzMatrix_typename_exposer.def(
              "zz"
              , zz_function_type( &::numeric::xyzMatrix< T >::zz )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      xyzMatrix_typename_exposer.staticmethod( "I" );
      xyzMatrix_typename_exposer.staticmethod( "cols" );
      xyzMatrix_typename_exposer.staticmethod( "diag" );
      xyzMatrix_typename_exposer.staticmethod( "identity" );
      xyzMatrix_typename_exposer.staticmethod( "rows" );
      { //property "xx"[fget=::numeric::xyzMatrix< T >::xx, fset=::numeric::xyzMatrix< T >::xx]

          typedef numeric::xyzMatrix< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzMatrix_typename_exposer.add_property(
              "xx"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::xx )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzMatrix< T >::xx )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "xy"[fget=::numeric::xyzMatrix< T >::xy, fset=::numeric::xyzMatrix< T >::xy]

          typedef numeric::xyzMatrix< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzMatrix_typename_exposer.add_property(
              "xy"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::xy )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzMatrix< T >::xy )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "xz"[fget=::numeric::xyzMatrix< T >::xz, fset=::numeric::xyzMatrix< T >::xz]

          typedef numeric::xyzMatrix< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzMatrix_typename_exposer.add_property(
              "xz"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::xz )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzMatrix< T >::xz )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "yx"[fget=::numeric::xyzMatrix< T >::yx, fset=::numeric::xyzMatrix< T >::yx]

          typedef numeric::xyzMatrix< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzMatrix_typename_exposer.add_property(
              "yx"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::yx )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzMatrix< T >::yx )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "yy"[fget=::numeric::xyzMatrix< T >::yy, fset=::numeric::xyzMatrix< T >::yy]

          typedef numeric::xyzMatrix< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzMatrix_typename_exposer.add_property(
              "yy"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::yy )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzMatrix< T >::yy )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "yz"[fget=::numeric::xyzMatrix< T >::yz, fset=::numeric::xyzMatrix< T >::yz]

          typedef numeric::xyzMatrix< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzMatrix_typename_exposer.add_property(
              "yz"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::yz )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzMatrix< T >::yz )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "zx"[fget=::numeric::xyzMatrix< T >::zx, fset=::numeric::xyzMatrix< T >::zx]

          typedef numeric::xyzMatrix< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzMatrix_typename_exposer.add_property(
              "zx"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::zx )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzMatrix< T >::zx )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "zy"[fget=::numeric::xyzMatrix< T >::zy, fset=::numeric::xyzMatrix< T >::zy]

          typedef numeric::xyzMatrix< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzMatrix_typename_exposer.add_property(
              "zy"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::zy )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzMatrix< T >::zy )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "zz"[fget=::numeric::xyzMatrix< T >::zz, fset=::numeric::xyzMatrix< T >::zz]

          typedef numeric::xyzMatrix< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzMatrix_typename_exposer.add_property(
              "zz"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::zz )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzMatrix< T >::zz )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "is_zero"[fget=::numeric::xyzMatrix< T >::is_zero]

          typedef numeric::xyzMatrix< T > fget_class_t;

          typedef bool ( fget_class_t::*fget )(  ) const;

          xyzMatrix_typename_exposer.add_property(
              "is_zero"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::is_zero )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "is_identity"[fget=::numeric::xyzMatrix< T >::is_identity]

          typedef numeric::xyzMatrix< T > fget_class_t;

          typedef bool ( fget_class_t::*fget )(  ) const;

          xyzMatrix_typename_exposer.add_property(
              "is_identity"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::is_identity )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "det"[fget=::numeric::xyzMatrix< T >::det]

          typedef numeric::xyzMatrix< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzMatrix_typename_exposer.add_property(
              "det"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::det )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "trace"[fget=::numeric::xyzMatrix< T >::trace]

          typedef numeric::xyzMatrix< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzMatrix_typename_exposer.add_property(
              "trace"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::trace )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "transposed"[fget=::numeric::xyzMatrix< T >::transposed]

          typedef numeric::xyzMatrix< T > fget_class_t;

          typedef ::numeric::xyzMatrix< T > ( fget_class_t::*fget )(  ) const;

          xyzMatrix_typename_exposer.add_property(
              "transposed"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::transposed )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }

      { //property "inverted"[fget=::numeric::xyzMatrix< T >::inverted]

          typedef numeric::xyzMatrix< T > fget_class_t;

          typedef ::numeric::xyzMatrix< T > ( fget_class_t::*fget )(  ) const;

          xyzMatrix_typename_exposer.add_property(
              "inverted"
              , bp::make_function(
                    fget( &::numeric::xyzMatrix< T >::inverse )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }

      xyzMatrix_typename_exposer.def( bp::other<  T  >() != bp::self );
      xyzMatrix_typename_exposer.def( bp::self != bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::self != bp::self );
      xyzMatrix_typename_exposer.def( bp::other<  T  >() * bp::self );
      xyzMatrix_typename_exposer.def( bp::self * bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::self * bp::self );
      xyzMatrix_typename_exposer.def( bp::self * bp::other< numeric::xyzVector< T > >() );
      xyzMatrix_typename_exposer.def( bp::other<  T  >() + bp::self );
      xyzMatrix_typename_exposer.def( bp::self + bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::self + bp::self );
      xyzMatrix_typename_exposer.def( bp::other<  T  >() - bp::self );
      xyzMatrix_typename_exposer.def( bp::self - bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::self - bp::self );
      xyzMatrix_typename_exposer.def( bp::self / bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::other<  T  >() < bp::self );
      xyzMatrix_typename_exposer.def( bp::self < bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::self < bp::self );
      xyzMatrix_typename_exposer.def( bp::other<  T  >() <= bp::self );
      xyzMatrix_typename_exposer.def( bp::self <= bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::self <= bp::self );
      xyzMatrix_typename_exposer.def( bp::other<  T  >() == bp::self );
      xyzMatrix_typename_exposer.def( bp::self == bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::self == bp::self );
      xyzMatrix_typename_exposer.def( bp::other<  T  >() > bp::self );
      xyzMatrix_typename_exposer.def( bp::self > bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::self > bp::self );
      xyzMatrix_typename_exposer.def( bp::other<  T  >() >= bp::self );
      xyzMatrix_typename_exposer.def( bp::self >= bp::other<  T  >() );
      xyzMatrix_typename_exposer.def( bp::self >= bp::self );

	  xyzMatrix_typename_exposer.def( bp::self_ns::str( bp::self ) );
  }

  bp::implicitly_convertible<  T  const &, numeric::xyzMatrix< T > >();

  { //::numeric::xyzTriple< T >
      typedef bp::class_< numeric::xyzTriple< T > > xyzTriple_typename_exposer_t;
      xyzTriple_typename_exposer_t xyzTriple_typename_exposer = xyzTriple_typename_exposer_t( std::string("xyzTriple_" + type_name).c_str() ); //"numeric___xyzTriple_ T "
      bp::scope xyzTriple_typename_scope( xyzTriple_typename_exposer );
      xyzTriple_typename_exposer.def( bp::init< >() );
      xyzTriple_typename_exposer.def( bp::init<  T  const & >(( bp::arg("t") )) );
      xyzTriple_typename_exposer.def( bp::init<  T  const &,  T  const &,  T  const & >(( bp::arg("x_a"), bp::arg("y_a"), bp::arg("z_a") )) );
      { //::numeric::xyzTriple< T >::assign

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*assign_function_type )(  T  const &, T  const &, T  const & ) ;

          xyzTriple_typename_exposer.def(
              "assign"
              , assign_function_type( &::numeric::xyzTriple< T >::assign )
              , ( bp::arg("x_a"), bp::arg("y_a"), bp::arg("z_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::clear

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*clear_function_type )(  ) ;

          xyzTriple_typename_exposer.def(
              "clear"
              , clear_function_type( &::numeric::xyzTriple< T >::clear )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::cross

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > ( exported_class_t::*cross_function_type )( ::numeric::xyzTriple< T > const & ) const;

          xyzTriple_typename_exposer.def(
              "cross"
              , cross_function_type( &::numeric::xyzTriple< T >::cross )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::cross_product

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > ( exported_class_t::*cross_product_function_type )( ::numeric::xyzTriple< T > const & ) const;

          xyzTriple_typename_exposer.def(
              "cross_product"
              , cross_product_function_type( &::numeric::xyzTriple< T >::cross_product )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::distance

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef  T  ( exported_class_t::*distance_function_type )( ::numeric::xyzTriple< T > const & ) const;

          xyzTriple_typename_exposer.def(
              "distance"
              , distance_function_type( &::numeric::xyzTriple< T >::distance )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::distance_squared

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef  T  ( exported_class_t::*distance_squared_function_type )( ::numeric::xyzTriple< T > const & ) const;

          xyzTriple_typename_exposer.def(
              "distance_squared"
              , distance_squared_function_type( &::numeric::xyzTriple< T >::distance_squared )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::dot

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef  T  ( exported_class_t::*dot_function_type )( ::numeric::xyzTriple< T > const & ) const;

          xyzTriple_typename_exposer.def(
              "dot"
              , dot_function_type( &::numeric::xyzTriple< T >::dot )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::dot_product

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef  T  ( exported_class_t::*dot_product_function_type )( ::numeric::xyzTriple< T > const & ) const;

          xyzTriple_typename_exposer.def(
              "dot_product"
              , dot_product_function_type( &::numeric::xyzTriple< T >::dot_product )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::equal_length

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef bool ( exported_class_t::*equal_length_function_type )( ::numeric::xyzTriple< T > const & ) ;

          xyzTriple_typename_exposer.def(
              "equal_length"
              , equal_length_function_type( &::numeric::xyzTriple< T >::equal_length )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::inner_product

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef  T  ( exported_class_t::*inner_product_function_type )( ::numeric::xyzTriple< T > const & ) const;

          xyzTriple_typename_exposer.def(
              "inner_product"
              , inner_product_function_type( &::numeric::xyzTriple< T >::inner_product )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::is_normalized

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef bool ( exported_class_t::*is_normalized_function_type )(  ) const;

          xyzTriple_typename_exposer.def(
              "is_normalized"
              , is_normalized_function_type( &::numeric::xyzTriple< T >::is_normalized )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::is_normalized

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef bool ( exported_class_t::*is_normalized_function_type )(  T  const & ) const;

          xyzTriple_typename_exposer.def(
              "is_normalized"
              , is_normalized_function_type( &::numeric::xyzTriple< T >::is_normalized )
              , ( bp::arg("tol") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::is_unit

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef bool ( exported_class_t::*is_unit_function_type )(  ) const;

          xyzTriple_typename_exposer.def(
              "is_unit"
              , is_unit_function_type( &::numeric::xyzTriple< T >::is_unit )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::is_unit

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef bool ( exported_class_t::*is_unit_function_type )(  T  const & ) const;

          xyzTriple_typename_exposer.def(
              "is_unit"
              , is_unit_function_type( &::numeric::xyzTriple< T >::is_unit )
              , ( bp::arg("tol") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::longer

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef bool ( exported_class_t::*longer_function_type )( ::numeric::xyzTriple< T > const & ) ;

          xyzTriple_typename_exposer.def(
              "longer"
              , longer_function_type( &::numeric::xyzTriple< T >::longer )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::longer_or_equal

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef bool ( exported_class_t::*longer_or_equal_function_type )( ::numeric::xyzTriple< T > const & ) ;

          xyzTriple_typename_exposer.def(
              "longer_or_equal"
              , longer_or_equal_function_type( &::numeric::xyzTriple< T >::longer_or_equal )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::max

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*max_function_type )( ::numeric::xyzTriple< T > const & ) ;

          xyzTriple_typename_exposer.def(
              "max"
              , max_function_type( &::numeric::xyzTriple< T >::max )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::min

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*min_function_type )( ::numeric::xyzTriple< T > const & ) ;

          xyzTriple_typename_exposer.def(
              "min"
              , min_function_type( &::numeric::xyzTriple< T >::min )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::negate

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*negate_function_type )(  ) ;

          xyzTriple_typename_exposer.def(
              "negate"
              , negate_function_type( &::numeric::xyzTriple< T >::negate )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::negated

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > ( exported_class_t::*negated_function_type )(  ) const;

          xyzTriple_typename_exposer.def(
              "negated"
              , negated_function_type( &::numeric::xyzTriple< T >::negated )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::negated

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef void ( exported_class_t::*negated_function_type )( ::numeric::xyzTriple< T > & ) const;

          xyzTriple_typename_exposer.def(
              "negated"
              , negated_function_type( &::numeric::xyzTriple< T >::negated )
              , ( bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalize

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*normalize_function_type )(  ) ;

          xyzTriple_typename_exposer.def(
              "normalize"
              , normalize_function_type( &::numeric::xyzTriple< T >::normalize )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalize

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*normalize_function_type )(  T  const & ) ;

          xyzTriple_typename_exposer.def(
              "normalize"
              , normalize_function_type( &::numeric::xyzTriple< T >::normalize )
              , ( bp::arg("length_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalize_any

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*normalize_any_function_type )(  ) ;

          xyzTriple_typename_exposer.def(
              "normalize_any"
              , normalize_any_function_type( &::numeric::xyzTriple< T >::normalize_any )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalize_any

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*normalize_any_function_type )(  T  const & ) ;

          xyzTriple_typename_exposer.def(
              "normalize_any"
              , normalize_any_function_type( &::numeric::xyzTriple< T >::normalize_any )
              , ( bp::arg("length_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalize_or_zero

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*normalize_or_zero_function_type )(  ) ;

          xyzTriple_typename_exposer.def(
              "normalize_or_zero"
              , normalize_or_zero_function_type( &::numeric::xyzTriple< T >::normalize_or_zero )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalize_or_zero

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*normalize_or_zero_function_type )(  T  const & ) ;

          xyzTriple_typename_exposer.def(
              "normalize_or_zero"
              , normalize_or_zero_function_type( &::numeric::xyzTriple< T >::normalize_or_zero )
              , ( bp::arg("length_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalized

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef void ( exported_class_t::*normalized_function_type )( ::numeric::xyzTriple< T > & ) const;

          xyzTriple_typename_exposer.def(
              "normalized"
              , normalized_function_type( &::numeric::xyzTriple< T >::normalized )
              , ( bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalized

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef void ( exported_class_t::*normalized_function_type )(  T  const &,::numeric::xyzTriple< T > & ) const;

          xyzTriple_typename_exposer.def(
              "normalized"
              , normalized_function_type( &::numeric::xyzTriple< T >::normalized )
              , ( bp::arg("length_a"), bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalized

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > ( exported_class_t::*normalized_function_type )(  ) const;

          xyzTriple_typename_exposer.def(
              "normalized"
              , normalized_function_type( &::numeric::xyzTriple< T >::normalized )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalized

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > ( exported_class_t::*normalized_function_type )(  T  const & ) const;

          xyzTriple_typename_exposer.def(
              "normalized"
              , normalized_function_type( &::numeric::xyzTriple< T >::normalized )
              , ( bp::arg("length_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalized_any

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef void ( exported_class_t::*normalized_any_function_type )( ::numeric::xyzTriple< T > & ) const;

          xyzTriple_typename_exposer.def(
              "normalized_any"
              , normalized_any_function_type( &::numeric::xyzTriple< T >::normalized_any )
              , ( bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalized_any

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef void ( exported_class_t::*normalized_any_function_type )(  T  const &,::numeric::xyzTriple< T > & ) const;

          xyzTriple_typename_exposer.def(
              "normalized_any"
              , normalized_any_function_type( &::numeric::xyzTriple< T >::normalized_any )
              , ( bp::arg("length_a"), bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalized_any

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > ( exported_class_t::*normalized_any_function_type )(  ) const;

          xyzTriple_typename_exposer.def(
              "normalized_any"
              , normalized_any_function_type( &::numeric::xyzTriple< T >::normalized_any )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalized_any

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > ( exported_class_t::*normalized_any_function_type )(  T  const & ) const;

          xyzTriple_typename_exposer.def(
              "normalized_any"
              , normalized_any_function_type( &::numeric::xyzTriple< T >::normalized_any )
              , ( bp::arg("length_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalized_or_zero

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef void ( exported_class_t::*normalized_or_zero_function_type )( ::numeric::xyzTriple< T > & ) const;

          xyzTriple_typename_exposer.def(
              "normalized_or_zero"
              , normalized_or_zero_function_type( &::numeric::xyzTriple< T >::normalized_or_zero )
              , ( bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalized_or_zero

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef void ( exported_class_t::*normalized_or_zero_function_type )(  T  const &,::numeric::xyzTriple< T > & ) const;

          xyzTriple_typename_exposer.def(
              "normalized_or_zero"
              , normalized_or_zero_function_type( &::numeric::xyzTriple< T >::normalized_or_zero )
              , ( bp::arg("length_a"), bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalized_or_zero

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > ( exported_class_t::*normalized_or_zero_function_type )(  ) const;

          xyzTriple_typename_exposer.def(
              "normalized_or_zero"
              , normalized_or_zero_function_type( &::numeric::xyzTriple< T >::normalized_or_zero )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::normalized_or_zero

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > ( exported_class_t::*normalized_or_zero_function_type )(  T  const & ) const;

          xyzTriple_typename_exposer.def(
              "normalized_or_zero"
              , normalized_or_zero_function_type( &::numeric::xyzTriple< T >::normalized_or_zero )
              , ( bp::arg("length_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::not_equal_length

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef bool ( exported_class_t::*not_equal_length_function_type )( ::numeric::xyzTriple< T > const & ) ;

          xyzTriple_typename_exposer.def(
              "not_equal_length"
              , not_equal_length_function_type( &::numeric::xyzTriple< T >::not_equal_length )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      xyzTriple_typename_exposer.def( bp::self *= bp::other<  T  >() );
      xyzTriple_typename_exposer.def( bp::self += bp::other<  T  >() );
      xyzTriple_typename_exposer.def( -bp::self );
      xyzTriple_typename_exposer.def( bp::self -= bp::other<  T  >() );
      xyzTriple_typename_exposer.def( bp::self /= bp::other<  T  >() );
      { //::numeric::xyzTriple< T >::operator[]

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef  T  const & ( exported_class_t::*__getitem___function_type )( int const ) const;

          xyzTriple_typename_exposer.def(
              "__getitem__"
              , __getitem___function_type( &::numeric::xyzTriple< T >::at )
              , ( bp::arg("i") )
              , bp::return_value_policy< bp::copy_const_reference >() );

      }
      { //::numeric::xyzTriple< T >::operator[]

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef  T  & ( exported_class_t::*__getitem___function_type )( int const ) ;

          xyzTriple_typename_exposer.def(
              "__getitem__"
              , __getitem___function_type( &::numeric::xyzTriple< T >::at )
              , ( bp::arg("i") )
              , bp::return_value_policy< bp::copy_non_const_reference >() );

      }
      { //::numeric::xyzTriple< T >::project_normal

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*project_normal_function_type )( ::numeric::xyzTriple< T > const & ) ;

          xyzTriple_typename_exposer.def(
              "project_normal"
              , project_normal_function_type( &::numeric::xyzTriple< T >::project_normal )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::project_parallel

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*project_parallel_function_type )( ::numeric::xyzTriple< T > const & ) ;

          xyzTriple_typename_exposer.def(
              "project_parallel"
              , project_parallel_function_type( &::numeric::xyzTriple< T >::project_parallel )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::projected_normal

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > ( exported_class_t::*projected_normal_function_type )( ::numeric::xyzTriple< T > const & ) const;

          xyzTriple_typename_exposer.def(
              "projected_normal"
              , projected_normal_function_type( &::numeric::xyzTriple< T >::projected_normal )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::projected_normal

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef void ( exported_class_t::*projected_normal_function_type )( ::numeric::xyzTriple< T > const &,::numeric::xyzTriple< T > & ) const;

          xyzTriple_typename_exposer.def(
              "projected_normal"
              , projected_normal_function_type( &::numeric::xyzTriple< T >::projected_normal )
              , ( bp::arg("v"), bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::projected_parallel

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > ( exported_class_t::*projected_parallel_function_type )( ::numeric::xyzTriple< T > const & ) const;

          xyzTriple_typename_exposer.def(
              "projected_parallel"
              , projected_parallel_function_type( &::numeric::xyzTriple< T >::projected_parallel )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::projected_parallel

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef void ( exported_class_t::*projected_parallel_function_type )( ::numeric::xyzTriple< T > const &,::numeric::xyzTriple< T > & ) ;

          xyzTriple_typename_exposer.def(
              "projected_parallel"
              , projected_parallel_function_type( &::numeric::xyzTriple< T >::projected_parallel )
              , ( bp::arg("v"), bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::shorter

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef bool ( exported_class_t::*shorter_function_type )( ::numeric::xyzTriple< T > const & ) ;

          xyzTriple_typename_exposer.def(
              "shorter"
              , shorter_function_type( &::numeric::xyzTriple< T >::shorter )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::shorter_or_equal

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef bool ( exported_class_t::*shorter_or_equal_function_type )( ::numeric::xyzTriple< T > const & ) ;

          xyzTriple_typename_exposer.def(
              "shorter_or_equal"
              , shorter_or_equal_function_type( &::numeric::xyzTriple< T >::shorter_or_equal )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::x

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef  T  & ( exported_class_t::*x_function_type )(  ) ;

          xyzTriple_typename_exposer.def(
              "x"
              , x_function_type( &::numeric::xyzTriple< T >::x )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::y

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef  T  & ( exported_class_t::*y_function_type )(  ) ;

          xyzTriple_typename_exposer.def(
              "y"
              , y_function_type( &::numeric::xyzTriple< T >::y )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::z

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef  T  & ( exported_class_t::*z_function_type )(  ) ;

          xyzTriple_typename_exposer.def(
              "z"
              , z_function_type( &::numeric::xyzTriple< T >::z )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzTriple< T >::zero

          typedef numeric::xyzTriple< T > exported_class_t;
          typedef ::numeric::xyzTriple< T > & ( exported_class_t::*zero_function_type )(  ) ;

          xyzTriple_typename_exposer.def(
              "zero"
              , zero_function_type( &::numeric::xyzTriple< T >::zero )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //property "x"[fget=::numeric::xyzTriple< T >::x, fset=::numeric::xyzTriple< T >::x]

          typedef numeric::xyzTriple< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzTriple_typename_exposer.add_property(
              "x"
              , bp::make_function(
                    fget( &::numeric::xyzTriple< T >::x )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzTriple< T >::x )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "y"[fget=::numeric::xyzTriple< T >::y, fset=::numeric::xyzTriple< T >::y]

          typedef numeric::xyzTriple< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzTriple_typename_exposer.add_property(
              "y"
              , bp::make_function(
                    fget( &::numeric::xyzTriple< T >::y )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzTriple< T >::y )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "z"[fget=::numeric::xyzTriple< T >::z, fset=::numeric::xyzTriple< T >::z]

          typedef numeric::xyzTriple< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzTriple_typename_exposer.add_property(
              "z"
              , bp::make_function(
                    fget( &::numeric::xyzTriple< T >::z )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzTriple< T >::z )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "is_zero"[fget=::numeric::xyzTriple< T >::is_zero]

          typedef numeric::xyzTriple< T > fget_class_t;

          typedef bool ( fget_class_t::*fget )(  ) const;

          xyzTriple_typename_exposer.add_property(
              "is_zero"
              , bp::make_function(
                    fget( &::numeric::xyzTriple< T >::is_zero )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "length"[fget=::numeric::xyzTriple< T >::length]

          typedef numeric::xyzTriple< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzTriple_typename_exposer.add_property(
              "length"
              , bp::make_function(
                    fget( &::numeric::xyzTriple< T >::length )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "length_squared"[fget=::numeric::xyzTriple< T >::length_squared]

          typedef numeric::xyzTriple< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzTriple_typename_exposer.add_property(
              "length_squared"
              , bp::make_function(
                    fget( &::numeric::xyzTriple< T >::length_squared )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "norm"[fget=::numeric::xyzTriple< T >::norm]

          typedef numeric::xyzTriple< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzTriple_typename_exposer.add_property(
              "norm"
              , bp::make_function(
                    fget( &::numeric::xyzTriple< T >::norm )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "norm_squared"[fget=::numeric::xyzTriple< T >::norm_squared]

          typedef numeric::xyzTriple< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzTriple_typename_exposer.add_property(
              "norm_squared"
              , bp::make_function(
                    fget( &::numeric::xyzTriple< T >::norm_squared )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "magnitude"[fget=::numeric::xyzTriple< T >::magnitude]

          typedef numeric::xyzTriple< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzTriple_typename_exposer.add_property(
              "magnitude"
              , bp::make_function(
                    fget( &::numeric::xyzTriple< T >::magnitude )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "magnitude_squared"[fget=::numeric::xyzTriple< T >::magnitude_squared]

          typedef numeric::xyzTriple< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzTriple_typename_exposer.add_property(
              "magnitude_squared"
              , bp::make_function(
                    fget( &::numeric::xyzTriple< T >::magnitude_squared )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      xyzTriple_typename_exposer.def( bp::other<  T  >() != bp::self );
      xyzTriple_typename_exposer.def( bp::self != bp::other<  T  >() );
      xyzTriple_typename_exposer.def( bp::self != bp::self );
      xyzTriple_typename_exposer.def( bp::other<  T  >() * bp::self );
      xyzTriple_typename_exposer.def( bp::self * bp::other<  T  >() );
      xyzTriple_typename_exposer.def( bp::other<  T  >() + bp::self );
      xyzTriple_typename_exposer.def( bp::self + bp::other<  T  >() );
      xyzTriple_typename_exposer.def( bp::self + bp::self );
      xyzTriple_typename_exposer.def( bp::other<  T  >() - bp::self );
      xyzTriple_typename_exposer.def( bp::self - bp::other<  T  >() );
      xyzTriple_typename_exposer.def( bp::self - bp::self );
      xyzTriple_typename_exposer.def( bp::self / bp::other<  T  >() );
      xyzTriple_typename_exposer.def( bp::other<  T  >() < bp::self );
      xyzTriple_typename_exposer.def( bp::self < bp::other<  T  >() );
      xyzTriple_typename_exposer.def( bp::self < bp::self );
      xyzTriple_typename_exposer.def( bp::other<  T  >() <= bp::self );
      xyzTriple_typename_exposer.def( bp::self <= bp::other<  T  >() );
      xyzTriple_typename_exposer.def( bp::self <= bp::self );
      xyzTriple_typename_exposer.def( bp::other<  T  >() == bp::self );
      xyzTriple_typename_exposer.def( bp::self == bp::other<  T  >() );
      xyzTriple_typename_exposer.def( bp::self == bp::self );
      xyzTriple_typename_exposer.def( bp::other<  T  >() > bp::self );
      xyzTriple_typename_exposer.def( bp::self > bp::other<  T  >() );
      xyzTriple_typename_exposer.def( bp::self > bp::self );
      xyzTriple_typename_exposer.def( bp::other<  T  >() >= bp::self );
      xyzTriple_typename_exposer.def( bp::self >= bp::other<  T  >() );
      xyzTriple_typename_exposer.def( bp::self >= bp::self );

	  xyzTriple_typename_exposer.def( bp::self_ns::str( bp::self ) );
  }

  bp::implicitly_convertible<  T  const &, numeric::xyzTriple< T > >();

  { //::numeric::xyzVector< T >
      typedef bp::class_< numeric::xyzVector< T > > xyzVector_typename_exposer_t;
      xyzVector_typename_exposer_t xyzVector_typename_exposer = xyzVector_typename_exposer_t( std::string("xyzVector_" + type_name).c_str() ); // "numeric___xyzVector_ T "
      bp::scope xyzVector_typename_scope( xyzVector_typename_exposer );
      xyzVector_typename_exposer.def( bp::init< >() );
      xyzVector_typename_exposer.def( bp::init<  T  const & >(( bp::arg("t") )) );
      xyzVector_typename_exposer.def( bp::init<  T  const &,  T  const &,  T  const & >(( bp::arg("x_a"), bp::arg("y_a"), bp::arg("z_a") )) );
      { //::numeric::xyzVector< T >::assign

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*assign_function_type )(  T  const &, T  const &, T  const & ) ;

          xyzVector_typename_exposer.def(
              "assign"
              , assign_function_type( &::numeric::xyzVector< T >::assign )
              , ( bp::arg("x_a"), bp::arg("y_a"), bp::arg("z_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::clear

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*clear_function_type )(  ) ;

          xyzVector_typename_exposer.def(
              "clear"
              , clear_function_type( &::numeric::xyzVector< T >::clear )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::cross

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*cross_function_type )( ::numeric::xyzVector< T > const & ) const;

          xyzVector_typename_exposer.def(
              "cross"
              , cross_function_type( &::numeric::xyzVector< T >::cross )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::cross_product

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*cross_product_function_type )( ::numeric::xyzVector< T > const & ) const;

          xyzVector_typename_exposer.def(
              "cross_product"
              , cross_product_function_type( &::numeric::xyzVector< T >::cross_product )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::distance

          typedef numeric::xyzVector< T > exported_class_t;
          typedef  T  ( exported_class_t::*distance_function_type )( ::numeric::xyzVector< T > const & ) const;

          xyzVector_typename_exposer.def(
              "distance"
              , distance_function_type( &::numeric::xyzVector< T >::distance )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::distance_squared

          typedef numeric::xyzVector< T > exported_class_t;
          typedef  T  ( exported_class_t::*distance_squared_function_type )( ::numeric::xyzVector< T > const & ) const;

          xyzVector_typename_exposer.def(
              "distance_squared"
              , distance_squared_function_type( &::numeric::xyzVector< T >::distance_squared )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::dot

          typedef numeric::xyzVector< T > exported_class_t;
          typedef  T  ( exported_class_t::*dot_function_type )( ::numeric::xyzVector< T > const & ) const;

          xyzVector_typename_exposer.def(
              "dot"
              , dot_function_type( &::numeric::xyzVector< T >::dot )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::dot_product

          typedef numeric::xyzVector< T > exported_class_t;
          typedef  T  ( exported_class_t::*dot_product_function_type )( ::numeric::xyzVector< T > const & ) const;

          xyzVector_typename_exposer.def(
              "dot_product"
              , dot_product_function_type( &::numeric::xyzVector< T >::dot_product )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::equal_length

          typedef numeric::xyzVector< T > exported_class_t;
          typedef bool ( exported_class_t::*equal_length_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzVector_typename_exposer.def(
              "equal_length"
              , equal_length_function_type( &::numeric::xyzVector< T >::equal_length )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::inner_product

          typedef numeric::xyzVector< T > exported_class_t;
          typedef  T  ( exported_class_t::*inner_product_function_type )( ::numeric::xyzVector< T > const & ) const;

          xyzVector_typename_exposer.def(
              "inner_product"
              , inner_product_function_type( &::numeric::xyzVector< T >::inner_product )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::is_normalized

          typedef numeric::xyzVector< T > exported_class_t;
          typedef bool ( exported_class_t::*is_normalized_function_type )(  ) const;

          xyzVector_typename_exposer.def(
              "is_normalized"
              , is_normalized_function_type( &::numeric::xyzVector< T >::is_normalized )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::is_normalized

          typedef numeric::xyzVector< T > exported_class_t;
          typedef bool ( exported_class_t::*is_normalized_function_type )(  T  const & ) const;

          xyzVector_typename_exposer.def(
              "is_normalized"
              , is_normalized_function_type( &::numeric::xyzVector< T >::is_normalized )
              , ( bp::arg("tol") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::is_unit

          typedef numeric::xyzVector< T > exported_class_t;
          typedef bool ( exported_class_t::*is_unit_function_type )(  ) const;

          xyzVector_typename_exposer.def(
              "is_unit"
              , is_unit_function_type( &::numeric::xyzVector< T >::is_unit )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::is_unit

          typedef numeric::xyzVector< T > exported_class_t;
          typedef bool ( exported_class_t::*is_unit_function_type )(  T  const & ) const;

          xyzVector_typename_exposer.def(
              "is_unit"
              , is_unit_function_type( &::numeric::xyzVector< T >::is_unit )
              , ( bp::arg("tol") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::longer

          typedef numeric::xyzVector< T > exported_class_t;
          typedef bool ( exported_class_t::*longer_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzVector_typename_exposer.def(
              "longer"
              , longer_function_type( &::numeric::xyzVector< T >::longer )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::longer_or_equal

          typedef numeric::xyzVector< T > exported_class_t;
          typedef bool ( exported_class_t::*longer_or_equal_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzVector_typename_exposer.def(
              "longer_or_equal"
              , longer_or_equal_function_type( &::numeric::xyzVector< T >::longer_or_equal )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::max

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*max_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzVector_typename_exposer.def(
              "max"
              , max_function_type( &::numeric::xyzVector< T >::max )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::min

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*min_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzVector_typename_exposer.def(
              "min"
              , min_function_type( &::numeric::xyzVector< T >::min )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::negate

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*negate_function_type )(  ) ;

          xyzVector_typename_exposer.def(
              "negate"
              , negate_function_type( &::numeric::xyzVector< T >::negate )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::negated

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*negated_function_type )(  ) const;

          xyzVector_typename_exposer.def(
              "negated"
              , negated_function_type( &::numeric::xyzVector< T >::negated )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::negated

          typedef numeric::xyzVector< T > exported_class_t;
          typedef void ( exported_class_t::*negated_function_type )( ::numeric::xyzVector< T > & ) const;

          xyzVector_typename_exposer.def(
              "negated"
              , negated_function_type( &::numeric::xyzVector< T >::negated )
              , ( bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalize

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*normalize_function_type )(  ) ;

          xyzVector_typename_exposer.def(
              "normalize"
              , normalize_function_type( &::numeric::xyzVector< T >::normalize )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalize

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*normalize_function_type )(  T  const & ) ;

          xyzVector_typename_exposer.def(
              "normalize"
              , normalize_function_type( &::numeric::xyzVector< T >::normalize )
              , ( bp::arg("length_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalize_any

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*normalize_any_function_type )(  ) ;

          xyzVector_typename_exposer.def(
              "normalize_any"
              , normalize_any_function_type( &::numeric::xyzVector< T >::normalize_any )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalize_any

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*normalize_any_function_type )(  T  const & ) ;

          xyzVector_typename_exposer.def(
              "normalize_any"
              , normalize_any_function_type( &::numeric::xyzVector< T >::normalize_any )
              , ( bp::arg("length_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalize_or_zero

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*normalize_or_zero_function_type )(  ) ;

          xyzVector_typename_exposer.def(
              "normalize_or_zero"
              , normalize_or_zero_function_type( &::numeric::xyzVector< T >::normalize_or_zero )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalize_or_zero

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*normalize_or_zero_function_type )(  T  const & ) ;

          xyzVector_typename_exposer.def(
              "normalize_or_zero"
              , normalize_or_zero_function_type( &::numeric::xyzVector< T >::normalize_or_zero )
              , ( bp::arg("length_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalized

          typedef numeric::xyzVector< T > exported_class_t;
          typedef void ( exported_class_t::*normalized_function_type )( ::numeric::xyzVector< T > & ) const;

          xyzVector_typename_exposer.def(
              "normalized"
              , normalized_function_type( &::numeric::xyzVector< T >::normalized )
              , ( bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalized

          typedef numeric::xyzVector< T > exported_class_t;
          typedef void ( exported_class_t::*normalized_function_type )(  T  const &,::numeric::xyzVector< T > & ) const;

          xyzVector_typename_exposer.def(
              "normalized"
              , normalized_function_type( &::numeric::xyzVector< T >::normalized )
              , ( bp::arg("length_a"), bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalized

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*normalized_function_type )(  ) const;

          xyzVector_typename_exposer.def(
              "normalized"
              , normalized_function_type( &::numeric::xyzVector< T >::normalized )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalized

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*normalized_function_type )(  T  const & ) const;

          xyzVector_typename_exposer.def(
              "normalized"
              , normalized_function_type( &::numeric::xyzVector< T >::normalized )
              , ( bp::arg("length_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalized_any

          typedef numeric::xyzVector< T > exported_class_t;
          typedef void ( exported_class_t::*normalized_any_function_type )( ::numeric::xyzVector< T > & ) const;

          xyzVector_typename_exposer.def(
              "normalized_any"
              , normalized_any_function_type( &::numeric::xyzVector< T >::normalized_any )
              , ( bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalized_any

          typedef numeric::xyzVector< T > exported_class_t;
          typedef void ( exported_class_t::*normalized_any_function_type )(  T  const &,::numeric::xyzVector< T > & ) const;

          xyzVector_typename_exposer.def(
              "normalized_any"
              , normalized_any_function_type( &::numeric::xyzVector< T >::normalized_any )
              , ( bp::arg("length_a"), bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalized_any

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*normalized_any_function_type )(  ) const;

          xyzVector_typename_exposer.def(
              "normalized_any"
              , normalized_any_function_type( &::numeric::xyzVector< T >::normalized_any )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalized_any

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*normalized_any_function_type )(  T  const & ) const;

          xyzVector_typename_exposer.def(
              "normalized_any"
              , normalized_any_function_type( &::numeric::xyzVector< T >::normalized_any )
              , ( bp::arg("length_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalized_or_zero

          typedef numeric::xyzVector< T > exported_class_t;
          typedef void ( exported_class_t::*normalized_or_zero_function_type )( ::numeric::xyzVector< T > & ) const;

          xyzVector_typename_exposer.def(
              "normalized_or_zero"
              , normalized_or_zero_function_type( &::numeric::xyzVector< T >::normalized_or_zero )
              , ( bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalized_or_zero

          typedef numeric::xyzVector< T > exported_class_t;
          typedef void ( exported_class_t::*normalized_or_zero_function_type )(  T  const &,::numeric::xyzVector< T > & ) const;

          xyzVector_typename_exposer.def(
              "normalized_or_zero"
              , normalized_or_zero_function_type( &::numeric::xyzVector< T >::normalized_or_zero )
              , ( bp::arg("length_a"), bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalized_or_zero

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*normalized_or_zero_function_type )(  ) const;

          xyzVector_typename_exposer.def(
              "normalized_or_zero"
              , normalized_or_zero_function_type( &::numeric::xyzVector< T >::normalized_or_zero )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::normalized_or_zero

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*normalized_or_zero_function_type )(  T  const & ) const;

          xyzVector_typename_exposer.def(
              "normalized_or_zero"
              , normalized_or_zero_function_type( &::numeric::xyzVector< T >::normalized_or_zero )
              , ( bp::arg("length_a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::not_equal_length

          typedef numeric::xyzVector< T > exported_class_t;
          typedef bool ( exported_class_t::*not_equal_length_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzVector_typename_exposer.def(
              "not_equal_length"
              , not_equal_length_function_type( &::numeric::xyzVector< T >::not_equal_length )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      xyzVector_typename_exposer.def( bp::self *= bp::other<  T  >() );
      xyzVector_typename_exposer.def( bp::self += bp::other<  T  >() );
      xyzVector_typename_exposer.def( -bp::self );
      xyzVector_typename_exposer.def( bp::self -= bp::other<  T  >() );
      xyzVector_typename_exposer.def( bp::self /= bp::other<  T  >() );
      { //::numeric::xyzVector< T >::operator[]

          typedef numeric::xyzVector< T > exported_class_t;
          typedef  T  const & ( exported_class_t::*__getitem___function_type )( int const ) const;

          xyzVector_typename_exposer.def(
              "__getitem__"
              , __getitem___function_type( &::numeric::xyzVector< T >::at )
              , ( bp::arg("i") )
              , bp::return_value_policy< bp::copy_const_reference >() );

      }
      { //::numeric::xyzVector< T >::operator[]

          typedef numeric::xyzVector< T > exported_class_t;
          typedef  T  & ( exported_class_t::*__getitem___function_type )( int const ) ;

          xyzVector_typename_exposer.def(
              "__getitem__"
              , __getitem___function_type( &::numeric::xyzVector< T >::at )
              , ( bp::arg("i") )
              , bp::return_value_policy< bp::copy_non_const_reference >() );

      }
      { //::numeric::xyzVector< T >::project_normal

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*project_normal_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzVector_typename_exposer.def(
              "project_normal"
              , project_normal_function_type( &::numeric::xyzVector< T >::project_normal )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::project_parallel

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*project_parallel_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzVector_typename_exposer.def(
              "project_parallel"
              , project_parallel_function_type( &::numeric::xyzVector< T >::project_parallel )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::projected_normal

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*projected_normal_function_type )( ::numeric::xyzVector< T > const & ) const;

          xyzVector_typename_exposer.def(
              "projected_normal"
              , projected_normal_function_type( &::numeric::xyzVector< T >::projected_normal )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::projected_normal

          typedef numeric::xyzVector< T > exported_class_t;
          typedef void ( exported_class_t::*projected_normal_function_type )( ::numeric::xyzVector< T > const &,::numeric::xyzVector< T > & ) const;

          xyzVector_typename_exposer.def(
              "projected_normal"
              , projected_normal_function_type( &::numeric::xyzVector< T >::projected_normal )
              , ( bp::arg("v"), bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::projected_parallel

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > ( exported_class_t::*projected_parallel_function_type )( ::numeric::xyzVector< T > const & ) const;

          xyzVector_typename_exposer.def(
              "projected_parallel"
              , projected_parallel_function_type( &::numeric::xyzVector< T >::projected_parallel )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::projected_parallel

          typedef numeric::xyzVector< T > exported_class_t;
          typedef void ( exported_class_t::*projected_parallel_function_type )( ::numeric::xyzVector< T > const &,::numeric::xyzVector< T > & ) ;

          xyzVector_typename_exposer.def(
              "projected_parallel"
              , projected_parallel_function_type( &::numeric::xyzVector< T >::projected_parallel )
              , ( bp::arg("v"), bp::arg("a") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::shorter

          typedef numeric::xyzVector< T > exported_class_t;
          typedef bool ( exported_class_t::*shorter_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzVector_typename_exposer.def(
              "shorter"
              , shorter_function_type( &::numeric::xyzVector< T >::shorter )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::shorter_or_equal

          typedef numeric::xyzVector< T > exported_class_t;
          typedef bool ( exported_class_t::*shorter_or_equal_function_type )( ::numeric::xyzVector< T > const & ) ;

          xyzVector_typename_exposer.def(
              "shorter_or_equal"
              , shorter_or_equal_function_type( &::numeric::xyzVector< T >::shorter_or_equal )
              , ( bp::arg("v") )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::x

          typedef numeric::xyzVector< T > exported_class_t;
          typedef  T  & ( exported_class_t::*x_function_type )(  ) ;

          xyzVector_typename_exposer.def(
              "x"
              , x_function_type( &::numeric::xyzVector< T >::x )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::y

          typedef numeric::xyzVector< T > exported_class_t;
          typedef  T  & ( exported_class_t::*y_function_type )(  ) ;

          xyzVector_typename_exposer.def(
              "y"
              , y_function_type( &::numeric::xyzVector< T >::y )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::z

          typedef numeric::xyzVector< T > exported_class_t;
          typedef  T  & ( exported_class_t::*z_function_type )(  ) ;

          xyzVector_typename_exposer.def(
              "z"
              , z_function_type( &::numeric::xyzVector< T >::z )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //::numeric::xyzVector< T >::zero

          typedef numeric::xyzVector< T > exported_class_t;
          typedef ::numeric::xyzVector< T > & ( exported_class_t::*zero_function_type )(  ) ;

          xyzVector_typename_exposer.def(
              "zero"
              , zero_function_type( &::numeric::xyzVector< T >::zero )
              , bp::return_value_policy< bp::return_by_value >() );

      }
      { //property "x"[fget=::numeric::xyzVector< T >::x, fset=::numeric::xyzVector< T >::x]

          typedef numeric::xyzVector< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzVector_typename_exposer.add_property(
              "x"
              , bp::make_function(
                    fget( &::numeric::xyzVector< T >::x )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzVector< T >::x )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "y"[fget=::numeric::xyzVector< T >::y, fset=::numeric::xyzVector< T >::y]

          typedef numeric::xyzVector< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzVector_typename_exposer.add_property(
              "y"
              , bp::make_function(
                    fget( &::numeric::xyzVector< T >::y )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzVector< T >::y )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "z"[fget=::numeric::xyzVector< T >::z, fset=::numeric::xyzVector< T >::z]

          typedef numeric::xyzVector< T > exported_class_t;

          typedef  T  const & ( exported_class_t::*fget )(  ) const;
          typedef void ( exported_class_t::*fset )(  T  const & ) ;

          xyzVector_typename_exposer.add_property(
              "z"
              , bp::make_function(
                    fget( &::numeric::xyzVector< T >::z )
                  , bp::return_value_policy< bp::return_by_value >() )
              , bp::make_function(
                    fset( &::numeric::xyzVector< T >::z )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "is_zero"[fget=::numeric::xyzVector< T >::is_zero]

          typedef numeric::xyzVector< T > fget_class_t;

          typedef bool ( fget_class_t::*fget )(  ) const;

          xyzVector_typename_exposer.add_property(
              "is_zero"
              , bp::make_function(
                    fget( &::numeric::xyzVector< T >::is_zero )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "length"[fget=::numeric::xyzVector< T >::length]

          typedef numeric::xyzVector< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzVector_typename_exposer.add_property(
              "length"
              , bp::make_function(
                    fget( &::numeric::xyzVector< T >::length )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "length_squared"[fget=::numeric::xyzVector< T >::length_squared]

          typedef numeric::xyzVector< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzVector_typename_exposer.add_property(
              "length_squared"
              , bp::make_function(
                    fget( &::numeric::xyzVector< T >::length_squared )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "norm"[fget=::numeric::xyzVector< T >::norm]

          typedef numeric::xyzVector< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzVector_typename_exposer.add_property(
              "norm"
              , bp::make_function(
                    fget( &::numeric::xyzVector< T >::norm )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "norm_squared"[fget=::numeric::xyzVector< T >::norm_squared]

          typedef numeric::xyzVector< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzVector_typename_exposer.add_property(
              "norm_squared"
              , bp::make_function(
                    fget( &::numeric::xyzVector< T >::norm_squared )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "magnitude"[fget=::numeric::xyzVector< T >::magnitude]

          typedef numeric::xyzVector< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzVector_typename_exposer.add_property(
              "magnitude"
              , bp::make_function(
                    fget( &::numeric::xyzVector< T >::magnitude )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      { //property "magnitude_squared"[fget=::numeric::xyzVector< T >::magnitude_squared]

          typedef numeric::xyzVector< T > fget_class_t;

          typedef  T  ( fget_class_t::*fget )(  ) const;

          xyzVector_typename_exposer.add_property(
              "magnitude_squared"
              , bp::make_function(
                    fget( &::numeric::xyzVector< T >::magnitude_squared )
                  , bp::return_value_policy< bp::return_by_value >() )  );

      }
      xyzVector_typename_exposer.def( bp::other<  T  >() != bp::self );
      xyzVector_typename_exposer.def( bp::self != bp::other<  T  >() );
      xyzVector_typename_exposer.def( bp::self != bp::self );
      xyzVector_typename_exposer.def( bp::other<  T  >() * bp::self );
      xyzVector_typename_exposer.def( bp::self * bp::other<  T  >() );
      xyzVector_typename_exposer.def( bp::other<  T  >() + bp::self );
      xyzVector_typename_exposer.def( bp::self + bp::other<  T  >() );
      xyzVector_typename_exposer.def( bp::self + bp::self );
      xyzVector_typename_exposer.def( bp::other<  T  >() - bp::self );
      xyzVector_typename_exposer.def( bp::self - bp::other<  T  >() );
      xyzVector_typename_exposer.def( bp::self - bp::self );
      xyzVector_typename_exposer.def( bp::self / bp::other<  T  >() );
      xyzVector_typename_exposer.def( bp::other<  T  >() < bp::self );
      xyzVector_typename_exposer.def( bp::self < bp::other<  T  >() );
      xyzVector_typename_exposer.def( bp::self < bp::self );
      xyzVector_typename_exposer.def( bp::other<  T  >() <= bp::self );
      xyzVector_typename_exposer.def( bp::self <= bp::other<  T  >() );
      xyzVector_typename_exposer.def( bp::self <= bp::self );
      xyzVector_typename_exposer.def( bp::other<  T  >() == bp::self );
      xyzVector_typename_exposer.def( bp::self == bp::other<  T  >() );
      xyzVector_typename_exposer.def( bp::self == bp::self );
      xyzVector_typename_exposer.def( bp::other<  T  >() > bp::self );
      xyzVector_typename_exposer.def( bp::self > bp::other<  T  >() );
      xyzVector_typename_exposer.def( bp::self > bp::self );
      xyzVector_typename_exposer.def( bp::other<  T  >() >= bp::self );
      xyzVector_typename_exposer.def( bp::self >= bp::other<  T  >() );
      xyzVector_typename_exposer.def( bp::self >= bp::self );

	  xyzVector_typename_exposer.def( bp::self_ns::str( bp::self ) );
  }

  bp::implicitly_convertible<  T  const &, numeric::xyzVector< T > >();
}

void __numeric_by_hand_ending__()
{
	// instantiate template classes at the end of numeric bindings so they act as default classes for Python for template bindings (some of them might get binded automatically due to template specifications)
	instantiate_numeric_containers<float>("float");
	instantiate_numeric_containers<numeric::Real>("Real");
	instantiate_real_numeric_containers<float>("float");
	instantiate_real_numeric_containers<numeric::Real>("Real");

	instantiate_numeric_containers<numeric::Size>("Size");
	instantiate_numeric_containers<numeric::SSize>("SSize");


	instantiate_numeric_functions<float>("float");
	instantiate_numeric_functions<numeric::Real>("Real");
}
