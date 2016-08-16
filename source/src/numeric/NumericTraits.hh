// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/NumericTraits.hh
/// @brief  Numeric type traits
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_NumericTraits_hh
#define INCLUDED_numeric_NumericTraits_hh


namespace numeric {


/// @brief NumericTraits: Numeric type traits
template< typename T >
struct NumericTraits
{
	typedef T Type;


	/// @brief Zero
	inline static Type zero() { return Type( 0.0L ); }

	/// @brief One
	inline static Type one()  { return Type( 1.0L ); }

	/// @brief Two
	inline static Type two()  { return Type( 2.0L ); }

	/// @brief pi
	inline static Type pi()                 { return Type( 3.14159265358979323846264338327950288L ); }

	/// @brief 2*pi
	inline static Type pi_2()               { return Type( 6.28318530717958647692528676655900577L ); }

	/// @brief pi/2
	inline static Type pi_over_2()          { return Type( 1.57079632679489661923132169163975145L ); }

	/// @brief pi/3
	inline static Type pi_over_3()          { return Type( 1.04719755119659774615421446109316763L ); }

	/// @brief pi/4
	inline static Type pi_over_4()          { return Type( 0.785398163397448309615660845819875721L ); }

	/// @brief (2*pi)/3
	inline static Type pi_2_over_3()        { return Type( 2.09439510239319549230842892218633527L ); }

	/// @brief (3*pi)/4
	inline static Type pi_3_over_4()        { return Type( 2.35619449019234492884698253745962716L ); }

	/// @brief pi/180
	inline static Type pi_over_180()        { return Type( 0.0174532925199432957692369076848861271L ); }

	/// @brief pi/180
	inline static Type degrees_to_radians() { return Type( 0.0174532925199432957692369076848861271L ); }

	/// @brief pi/180
	inline static Type deg2rad()            { return Type( 0.0174532925199432957692369076848861271L ); }

	/// @brief 180/pi
	inline static Type radians_to_degrees() { return Type( 57.2957795130823208767981548141051703L ); }

	/// @brief 180/pi
	inline static Type rad2deg()            { return Type( 57.2957795130823208767981548141051703L ); }


}; // NumericTraits


/// @brief NumericTraits: Numeric type traits float specialization
template<>
struct NumericTraits< float >
{
	typedef float Type;


	/// @brief Zero
	inline static Type zero() { return Type( 0.0F ); }

	/// @brief One
	inline static Type one()  { return Type( 1.0F ); }

	/// @brief Two
	inline static Type two()  { return Type( 2.0F ); }

	/// @brief pi
	inline static Type pi()                 { return 3.14159265358979323846264338327950288F; }

	/// @brief 2*pi
	inline static Type pi_2()               { return 6.28318530717958647692528676655900577F; }

	/// @brief pi/2
	inline static Type pi_over_2()          { return 1.57079632679489661923132169163975145F; }

	/// @brief pi/3
	inline static Type pi_over_3()          { return 1.04719755119659774615421446109316763F; }

	/// @brief pi/4
	inline static Type pi_over_4()          { return 0.785398163397448309615660845819875721F; }

	/// @brief (2*pi)/3
	inline static Type pi_2_over_3()        { return 2.09439510239319549230842892218633527F; }

	/// @brief (3*pi)/4
	inline static Type pi_3_over_4()        { return 2.35619449019234492884698253745962716F; }

	/// @brief pi/180
	inline static Type pi_over_180()        { return 0.0174532925199432957692369076848861271F; }

	/// @brief pi/180
	inline static Type degrees_to_radians() { return 0.0174532925199432957692369076848861271F; }

	/// @brief pi/180
	inline static Type deg2rad()            { return 0.0174532925199432957692369076848861271F; }

	/// @brief 180/pi
	inline static Type radians_to_degrees() { return 57.2957795130823208767981548141051703F; }

	/// @brief 180/pi
	inline static Type rad2deg()            { return 57.2957795130823208767981548141051703F; }

	/// @brief Tolerance
	inline static Type tolerance() { return 1.0E-6F; }

	/// @brief Length tolerance
	inline static Type length_tolerance() { return 1.0E-6F; }

	/// @brief Angle tolerance (radians)
	inline static Type angle_tolerance() { return 1.0E-6F; }

	/// @brief Sine cosine range tolerance
	inline static Type sin_cos_tolerance() { return 1.0E-6F; }

	/// @brief Quaternion normalization tolerance
	inline static Type quaternion_tolerance() { return 1.0E-6F; }


}; // NumericTraits


/// @brief NumericTraits: Numeric type traits double specialization
template<>
struct NumericTraits< double >
{
	typedef double Type;


	/// @brief Zero
	inline static Type zero() { return Type( 0.0 ); }

	/// @brief One
	inline static Type one()  { return Type( 1.0 ); }

	/// @brief Two
	inline static Type two()  { return Type( 2.0 ); }

	/// @brief pi
	inline static Type pi()                 { return 3.14159265358979323846264338327950288; }

	/// @brief 2*pi
	inline static Type pi_2()               { return 6.28318530717958647692528676655900577; }

	/// @brief pi/2
	inline static Type pi_over_2()          { return 1.57079632679489661923132169163975145; }

	/// @brief pi/3
	inline static Type pi_over_3()          { return 1.04719755119659774615421446109316763; }

	/// @brief (2*pi)/3
	inline static Type pi_2_over_3()        { return 2.09439510239319549230842892218633527; }

	/// @brief (3*pi)/4
	inline static Type pi_3_over_4()        { return 2.35619449019234492884698253745962716; }

	/// @brief pi/4
	inline static Type pi_over_4()          { return 0.785398163397448309615660845819875721; }

	/// @brief pi/180
	inline static Type pi_over_180()        { return 0.0174532925199432957692369076848861271; }

	/// @brief pi/180
	inline static Type degrees_to_radians() { return 0.0174532925199432957692369076848861271; }

	/// @brief pi/180
	inline static Type deg2rad()            { return 0.0174532925199432957692369076848861271; }

	/// @brief 180/pi
	inline static Type radians_to_degrees() { return 57.2957795130823208767981548141051703; }

	/// @brief 180/pi
	inline static Type rad2deg()            { return 57.2957795130823208767981548141051703; }

	/// @brief Tolerance
	inline static Type tolerance() { return 1.0E-9; }

	/// @brief Length tolerance
	inline static Type length_tolerance() { return 1.0E-9; }

	/// @brief Angle tolerance (radians)
	inline static Type angle_tolerance() { return 1.0E-9; }

	/// @brief Sine cosine range tolerance
	inline static Type sin_cos_tolerance() { return 1.0E-9; }

	/// @brief Quaternion normalization tolerance
	inline static Type quaternion_tolerance() { return 1.0E-9; }


}; // NumericTraits


/// @brief NumericTraits: Numeric type traits long double specialization
template<>
struct NumericTraits< long double >
{
	typedef long double Type;


	/// @brief Zero
	inline static Type zero() { return Type( 0.0L ); }

	/// @brief One
	inline static Type one()  { return Type( 1.0L ); }

	/// @brief Two
	inline static Type two()  { return Type( 2.0L ); }

	/// @brief pi
	inline static Type pi()                 { return 3.14159265358979323846264338327950288L; }

	/// @brief 2*pi
	inline static Type pi_2()               { return 6.28318530717958647692528676655900577L; }

	/// @brief pi/2
	inline static Type pi_over_2()          { return 1.57079632679489661923132169163975145L; }

	/// @brief pi/3
	inline static Type pi_over_3()          { return 1.04719755119659774615421446109316763L; }

	/// @brief pi/4
	inline static Type pi_over_4()          { return 0.785398163397448309615660845819875721L; }

	/// @brief (2*pi)/3
	inline static Type pi_2_over_3()        { return 2.09439510239319549230842892218633527L; }

	/// @brief (3*pi)/4
	inline static Type pi_3_over_4()        { return 2.35619449019234492884698253745962716L; }

	/// @brief pi/180
	inline static Type pi_over_180()        { return 0.0174532925199432957692369076848861271L; }

	/// @brief pi/180
	inline static Type degrees_to_radians() { return 0.0174532925199432957692369076848861271L; }

	/// @brief pi/180
	inline static Type deg2rad()            { return 0.0174532925199432957692369076848861271L; }

	/// @brief 180/pi
	inline static Type radians_to_degrees() { return 57.2957795130823208767981548141051703L; }

	/// @brief 180/pi
	inline static Type rad2deg()            { return 57.2957795130823208767981548141051703L; }

	/// @brief Tolerance
	inline static Type tolerance() { return 1.0E-9L; }

	/// @brief Length tolerance
	inline static Type length_tolerance() { return 1.0E-9L; }

	/// @brief Angle tolerance (radians)
	inline static Type angle_tolerance() { return 1.0E-9L; }

	/// @brief Sine cosine range tolerance
	inline static Type sin_cos_tolerance() { return 1.0E-9L; }

	/// @brief Quaternion normalization tolerance
	inline static Type quaternion_tolerance() { return 1.0E-9L; }


}; // NumericTraits


} // namespace numeric


#endif // INCLUDED_numeric_NumericTraits_HH
