#ifndef INCLUDED_ObjexxFCL_FArray5A_hh
#define INCLUDED_ObjexxFCL_FArray5A_hh


// FArray5A: Fortran-Compatible 5D Argument Array
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/FArray5A.fwd.hh>
#include <ObjexxFCL/FArray5P.hh>
#include <ObjexxFCL/StaticIndexRange.hh>


namespace ObjexxFCL {


/// @brief FArray5A: Fortran-Compatible 5D Argument Array
template< typename T >
class FArray5A :
	public FArray5< T >
{


private: // Types


	typedef  FArray5< T >  Super;
	typedef  typename Super::real_FArray  real_FArray;
	typedef  typename Super::proxy_FArray  proxy_FArray;
	typedef  typename Super::arg_FArray  arg_FArray;
	typedef  internal::ProxySentinel  ProxySentinel;


public: // Types


	typedef  typename Super::Base  Base;
	typedef  typename Base::Section  Section;
	typedef  typename Super::IR  SIR;
	typedef  StaticIndexRange  IR;

	// STL Style
	typedef  typename Base::value_type  value_type;
	typedef  typename Base::reference  reference;
	typedef  typename Base::const_reference  const_reference;
	typedef  typename Base::pointer  pointer;
	typedef  typename Base::const_pointer  const_pointer;
	typedef  typename Base::size_type  size_type;
	typedef  typename Base::difference_type  difference_type;

	// C++ Style
	typedef  typename Base::Value  Value;
	typedef  typename Base::Reference  Reference;
	typedef  typename Base::ConstReference  ConstReference;
	typedef  typename Base::Pointer  Pointer;
	typedef  typename Base::ConstPointer  ConstPointer;
	typedef  typename Base::Size  Size;
	typedef  typename Base::Difference  Difference;

	using Super::array_;
	using Super::array_size_;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
	using Super::not_const_proxy;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	using Super::npos;
	using Super::sarray_;
	using Super::shift_;
	using Super::shift_set;
	using Super::size_set;
	using Super::s1_;
	using Super::s2_;
	using Super::s3_;
	using Super::s4_;


public: // Creation


	/// @brief Default Constructor
	inline
	FArray5A() :
		Super( ProxySentinel() )
	{}


	/// @brief Copy Constructor
	inline
	FArray5A( FArray5A const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ ),
		I4_( a.I4_ ),
		I5_( a.I5_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
	}


	/// @brief Non-Const Copy Constructor
	inline
	FArray5A( FArray5A & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ ),
		I4_( a.I4_ ),
		I5_( a.I5_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
	}


	/// @brief Real Constructor
	inline
	FArray5A( real_FArray const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ ),
		I4_( a.I4_ ),
		I5_( a.I5_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
	}


	/// @brief Non-Const Real Constructor
	inline
	FArray5A( real_FArray & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ ),
		I4_( a.I4_ ),
		I5_( a.I5_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
	}


	/// @brief Proxy Constructor
	inline
	FArray5A( proxy_FArray const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ ),
		I4_( a.I4_ ),
		I5_( a.I5_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
	}


	/// @brief Non-Const Proxy Constructor
	inline
	FArray5A( proxy_FArray & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ ),
		I4_( a.I4_ ),
		I5_( a.I5_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
	}


	/// @brief Super Constructor
	inline
	FArray5A( Super const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1() ),
		I2_( a.I2() ),
		I3_( a.I3() ),
		I4_( a.I4() ),
		I5_( a.I5() )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
	}


	/// @brief Non-Const Super Constructor
	inline
	FArray5A( Super & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1() ),
		I2_( a.I2() ),
		I3_( a.I3() ),
		I4_( a.I4() ),
		I5_( a.I5() )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
	}


	/// @brief Base Constructor
	inline
	FArray5A( Base const & a ) :
		Super( a, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( 1 ),
		I4_( 1 ),
		I5_( a.size() )
	{
		shift_set( 5 );
		s1_ = s2_ = s3_ = s4_ = 1;
	}


	/// @brief Non-Const Base Constructor
	inline
	FArray5A( Base & a ) :
		Super( a, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( 1 ),
		I4_( 1 ),
		I5_( a.size() )
	{
		shift_set( 5 );
		s1_ = s2_ = s3_ = s4_ = 1;
	}


	/// @brief Section Constructor
	inline
	FArray5A( Section const & s ) :
		Super( s, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( 1 ),
		I4_( 1 ),
		I5_( s.size() )
	{
		shift_set( 5 );
		s1_ = s2_ = s3_ = s4_ = 1;
	}


	/// @brief Non-Const Section Constructor
	inline
	FArray5A( Section & s ) :
		Super( s, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( 1 ),
		I4_( 1 ),
		I5_( s.size() )
	{
		shift_set( 5 );
		s1_ = s2_ = s3_ = s4_ = 1;
	}


	/// @brief Value Constructor
	inline
	FArray5A( T const & t ) :
		Super( t, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( 1 ),
		I4_( 1 ),
		I5_( star ) // Unbounded
	{
		shift_set( 5 );
		s1_ = s2_ = s3_ = s4_ = 1;
	}


	/// @brief Non-Const Value Constructor
	inline
	FArray5A( T & t ) :
		Super( t, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( 1 ),
		I4_( 1 ),
		I5_( star ) // Unbounded
	{
		shift_set( 5 );
		s1_ = s2_ = s3_ = s4_ = 1;
	}


	/// @brief Copy + IndexRange Constructor
	inline
	FArray5A( FArray5A const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Copy + IndexRange Constructor
	inline
	FArray5A( FArray5A & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Real + IndexRange Constructor
	inline
	FArray5A( real_FArray const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Real + IndexRange Constructor
	inline
	FArray5A( real_FArray & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Proxy + IndexRange Constructor
	inline
	FArray5A( proxy_FArray const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Proxy + IndexRange Constructor
	inline
	FArray5A( proxy_FArray & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Super + IndexRange Constructor
	inline
	FArray5A( Super const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Super + IndexRange Constructor
	inline
	FArray5A( Super & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Base + IndexRange Constructor
	inline
	FArray5A( Base const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Base + IndexRange Constructor
	inline
	FArray5A( Base & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Section + IndexRange Constructor
	inline
	FArray5A( Section const & s, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( s, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Section + IndexRange Constructor
	inline
	FArray5A( Section & s, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( s, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Value + IndexRange Constructor
	inline
	FArray5A( T const & t, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( t, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Value + IndexRange Constructor
	inline
	FArray5A( T & t, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) :
		Super( t, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a ),
		I5_( I5_a )
	{
		dimension_argument();
	}


	/// @brief Destructor
	inline
	virtual
	~FArray5A()
	{}


public: // Assignment


	/// @brief Copy Assignment
	inline
	FArray5A &
	operator =( FArray5A const & a )
	{
		if ( this != &a ) {
			if ( ! equal_dimension( a ) ) dimension( a );
			Base::operator =( a );
		}
		return *this;
	}


	/// @brief Super Assignment
	inline
	FArray5A &
	operator =( Super const & a )
	{
		if ( this != &a ) {
			if ( ! equal_dimension( a ) ) dimension( a );
			Base::operator =( a );
		}
		return *this;
	}


	/// @brief Super Assignment Template
	template< typename U >
	inline
	FArray5A &
	operator =( FArray5< U > const & a )
	{
		if ( ! equal_dimension( a ) ) dimension( a );
		Base::operator =( a );
		return *this;
	}


	/// @brief += Array Template
	template< typename U >
	inline
	FArray5A &
	operator +=( FArray5< U > const & a )
	{
		Super::operator +=( a );
		return *this;
	}


	/// @brief -= Array Template
	template< typename U >
	inline
	FArray5A &
	operator -=( FArray5< U > const & a )
	{
		Super::operator -=( a );
		return *this;
	}


	/// @brief = Value
	inline
	FArray5A &
	operator =( T const & t )
	{
		Super::operator =( t );
		return *this;
	}


	/// @brief += Value
	inline
	FArray5A &
	operator +=( T const & t )
	{
		Super::operator +=( t );
		return *this;
	}


	/// @brief -= Value
	inline
	FArray5A &
	operator -=( T const & t )
	{
		Super::operator -=( t );
		return *this;
	}


	/// @brief *= Value
	inline
	FArray5A &
	operator *=( T const & t )
	{
		Super::operator *=( t );
		return *this;
	}


	/// @brief /= Value
	inline
	FArray5A &
	operator /=( T const & t )
	{
		Super::operator /=( t );
		return *this;
	}


public: // Subscript


	/// @brief array( i1, i2, i3, i4, i5 ) const
	inline
	T const &
	operator ()( int const i1, int const i2, int const i3, int const i4, int const i5 ) const
	{
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) && ( I4_.contains( i4 ) ) && ( I5_.contains( i5 ) ) );
		return sarray_[ ( ( ( ( ( ( ( i5 * s4_ ) + i4 ) * s3_ ) + i3 ) * s2_ ) + i2 ) * s1_ ) + i1 ];
	}


	/// @brief array( i1, i2, i3, i4, i5 )
	inline
	T &
	operator ()( int const i1, int const i2, int const i3, int const i4, int const i5 )
	{
		proxy_const_assert( not_const_proxy() );
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) && ( I4_.contains( i4 ) ) && ( I5_.contains( i5 ) ) );
		return sarray_[ ( ( ( ( ( ( ( i5 * s4_ ) + i4 ) * s3_ ) + i3 ) * s2_ ) + i2 ) * s1_ ) + i1 ];
	}


	/// @brief Section Starting at array( i1, i2, i3, i4, i5 )
	inline
	Section const
	a( int const i1, int const i2, int const i3, int const i4, int const i5 ) const
	{
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) && ( I4_.contains( i4 ) ) && ( I5_.contains( i5 ) ) );
		size_type const offset( ( ( ( ( ( ( ( ( i5 * s4_ ) + i4 ) * s3_ ) + i3 ) * s2_ ) + i2 ) * s1_ ) + i1 ) - shift_ );
		return Section( static_cast< T const * >( array_ + offset ), ( array_size_ != npos ) ? array_size_ - offset : npos );
	}


	/// @brief Non-Const Section Starting at array( i1, i2, i3, i4, i5 )
	inline
	Section
	a( int const i1, int const i2, int const i3, int const i4, int const i5 )
	{
		proxy_const_assert( not_const_proxy() );
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) && ( I4_.contains( i4 ) ) && ( I5_.contains( i5 ) ) );
		size_type const offset( ( ( ( ( ( ( ( ( i5 * s4_ ) + i4 ) * s3_ ) + i3 ) * s2_ ) + i2 ) * s1_ ) + i1 ) - shift_ );
		return Section( array_ + offset, ( array_size_ != npos ) ? array_size_ - offset : npos );
	}


	/// @brief Linear Index
	inline
	size_type
	index( int const i1, int const i2, int const i3, int const i4, int const i5 ) const
	{
		return ( ( ( ( ( ( ( ( ( i5 * s4_ ) + i4 ) * s3_ ) + i3 ) * s2_ ) + i2 ) * s1_ ) + i1 ) - shift_ );
	}


public: // Predicate


	/// @brief Dimensions Initialized?
	inline
	bool
	dimensions_initialized() const
	{
		return true;
	}


	/// @brief Contains Indexed Element?
	inline
	bool
	contains( int const i1, int const i2, int const i3, int const i4, int const i5 ) const
	{
		return ( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) && ( I4_.contains( i4 ) ) && ( I5_.contains( i5 ) ) );
	}


	/// @brief Initializer Active?
	inline
	bool
	initializer_active() const
	{
		return false;
	}


public: // Inspector


	/// @brief IndexRange of Dimension 1
	inline
	IR const &
	I1() const
	{
		return I1_;
	}


	/// @brief Lower Index of Dimension 1
	inline
	int
	l1() const
	{
		return I1_.l();
	}


	/// @brief Upper Index of Dimension 1
	inline
	int
	u1() const
	{
		return I1_.u();
	}


	/// @brief IndexRange of Dimension 2
	inline
	IR const &
	I2() const
	{
		return I2_;
	}


	/// @brief Lower Index of Dimension 2
	inline
	int
	l2() const
	{
		return I2_.l();
	}


	/// @brief Upper Index of Dimension 2
	inline
	int
	u2() const
	{
		return I2_.u();
	}


	/// @brief IndexRange of Dimension 3
	inline
	IR const &
	I3() const
	{
		return I3_;
	}


	/// @brief Lower Index of Dimension 3
	inline
	int
	l3() const
	{
		return I3_.l();
	}


	/// @brief Upper Index of Dimension 3
	inline
	int
	u3() const
	{
		return I3_.u();
	}


	/// @brief IndexRange of Dimension 4
	inline
	IR const &
	I4() const
	{
		return I4_;
	}


	/// @brief Lower Index of Dimension 4
	inline
	int
	l4() const
	{
		return I4_.l();
	}


	/// @brief Upper Index of Dimension 4
	inline
	int
	u4() const
	{
		return I4_.u();
	}


	/// @brief IndexRange of Dimension 5
	inline
	IR const &
	I5() const
	{
		return I5_;
	}


	/// @brief Lower Index of Dimension 5
	inline
	int
	l5() const
	{
		return I5_.l();
	}


	/// @brief Upper Index of Dimension 5
	inline
	int
	u5() const
	{
		return I5_.u();
	}


	/// @brief Size of Dimension 5
	inline
	size_type
	size5() const
	{
		return I5_.size();
	}


public: // Modifier


	/// @brief Clear
	inline
	FArray5A &
	clear()
	{
		Super::clear();
		I1_.clear();
		I2_.clear();
		I3_.clear();
		I4_.clear();
		I5_.clear();
		return *this;
	}


	/// @brief Dimension by IndexRanges Even if Const
	inline
	FArray5A const &
	dim( IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a ) const
	{
		const_cast< FArray5A & >( *this ).dimension( I1_a, I2_a, I3_a, I4_a, I5_a );
		return *this;
	}


	/// @brief Dimension by Array Even if Const
	template< typename U >
	inline
	FArray5A const &
	dim( FArray5< U > const & a ) const
	{
		const_cast< FArray5A & >( *this ).dimension( a );
		return *this;
	}


	/// @brief Dimension by IndexRange
	inline
	FArray5A &
	dimension( IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a, IR const & I5_a )
	{
		I1_.assign_value_of( I1_a );
		I2_.assign_value_of( I2_a );
		I3_.assign_value_of( I3_a );
		I4_.assign_value_of( I4_a );
		I5_.assign_value_of( I5_a );
		dimension_argument();
		return *this;
	}


	/// @brief Dimension by Array
	template< typename U >
	inline
	FArray5A &
	dimension( FArray5< U > const & a )
	{
		I1_.assign_value_of( a.I1() );
		I2_.assign_value_of( a.I2() );
		I3_.assign_value_of( a.I3() );
		I4_.assign_value_of( a.I4() );
		I5_.assign_value_of( a.I5() );
		dimension_argument();
		return *this;
	}


	/// @brief Attach to Argument Array
	inline
	FArray5A &
	attach( FArray5A const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		I4_.assign_value_of( a.I4_ );
		I5_.assign_value_of( a.I5_ );
		return *this;
	}


	/// @brief Attach to Non-Const Argument Array
	inline
	FArray5A &
	attach( FArray5A & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		I4_.assign_value_of( a.I4_ );
		I5_.assign_value_of( a.I5_ );
		return *this;
	}


	/// @brief Attach to Real Array
	inline
	FArray5A &
	attach( real_FArray const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		I4_.assign_value_of( a.I4_ );
		I5_.assign_value_of( a.I5_ );
		return *this;
	}


	/// @brief Attach to Non-Const Real Array
	inline
	FArray5A &
	attach( real_FArray & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		I4_.assign_value_of( a.I4_ );
		I5_.assign_value_of( a.I5_ );
		return *this;
	}


	/// @brief Attach to Proxy Array
	inline
	FArray5A &
	attach( proxy_FArray const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		I4_.assign_value_of( a.I4_ );
		I5_.assign_value_of( a.I5_ );
		return *this;
	}


	/// @brief Attach to Non-Const Proxy Array
	inline
	FArray5A &
	attach( proxy_FArray & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		I4_.assign_value_of( a.I4_ );
		I5_.assign_value_of( a.I5_ );
		return *this;
	}


	/// @brief Attach to Super Array
	inline
	FArray5A &
	attach( Super const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
		I1_.assign_value_of( a.I1() );
		I2_.assign_value_of( a.I2() );
		I3_.assign_value_of( a.I3() );
		I4_.assign_value_of( a.I4() );
		I5_.assign_value_of( a.I5() );
		return *this;
	}


	/// @brief Attach to Non-Const Super Array
	inline
	FArray5A &
	attach( Super & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		s4_ = a.s4_;
		I1_.assign_value_of( a.I1() );
		I2_.assign_value_of( a.I2() );
		I3_.assign_value_of( a.I3() );
		I4_.assign_value_of( a.I4() );
		I5_.assign_value_of( a.I5() );
		return *this;
	}


	/// @brief Attach to Base Array
	inline
	FArray5A &
	attach( Base const & a )
	{
		Base::attach( a, 5 );
		s1_ = s2_ = s3_ = s4_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = 1;
		I4_ = 1;
		I5_ = a.size();
		return *this;
	}


	/// @brief Attach to Non-Const Base Array
	inline
	FArray5A &
	attach( Base & a )
	{
		Base::attach( a, 5 );
		s1_ = s2_ = s3_ = s4_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = 1;
		I4_ = 1;
		I5_ = a.size();
		return *this;
	}


	/// @brief Attach to Section
	inline
	FArray5A &
	attach( Section const & s )
	{
		Base::attach( s, 5 );
		s1_ = s2_ = s3_ = s4_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = 1;
		I4_ = 1;
		I5_ = s.size();
		return *this;
	}


	/// @brief Attach to Non-Const Section
	inline
	FArray5A &
	attach( Section & s )
	{
		Base::attach( s, 5 );
		s1_ = s2_ = s3_ = s4_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = 1;
		I4_ = 1;
		I5_ = s.size();
		return *this;
	}


	/// @brief Attach to Value
	inline
	FArray5A &
	attach( T const & t )
	{
		Base::attach( t, 5 );
		s1_ = s2_ = s3_ = s4_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = 1;
		I4_ = 1;
		I5_ = star; // Unbounded
		return *this;
	}


	/// @brief Attach to Non-Const Value
	inline
	FArray5A &
	attach( T & t )
	{
		Base::attach( t, 5 );
		s1_ = s2_ = s3_ = s4_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = 1;
		I4_ = 1;
		I5_ = star; // Unbounded
		return *this;
	}


	/// @brief Detach from Source Array
	inline
	FArray5A &
	detach()
	{
		Base::detach();
		s1_ = s2_ = s3_ = s4_ = 0;
		I1_.clear();
		I2_.clear();
		I3_.clear();
		I4_.clear();
		I5_.clear();
		return *this;
	}


protected: // Functions


	/// @brief Dimension by IndexRanges
	inline
	void
	dimension_assign( SIR const & I1_a, SIR const & I2_a, SIR const & I3_a, SIR const & I4_a, SIR const & I5_a )
	{
		I1_.assign_value_of( I1_a );
		I2_.assign_value_of( I2_a );
		I3_.assign_value_of( I3_a );
		I4_.assign_value_of( I4_a );
		I5_.assign_value_of( I5_a );
		dimension_argument();
	}


private: // Functions


	/// @brief Dimension by Current IndexRanges
	inline
	void
	dimension_argument()
	{
		assert( I1_.bounded_value() );
		assert( I2_.bounded_value() );
		assert( I3_.bounded_value() );
		assert( I4_.bounded_value() );
		s1_ = I1_.size();
		s2_ = I2_.size();
		s3_ = I3_.size();
		s4_ = I4_.size();
		if ( I5_.bounded_value() ) { // Bounded
			size_set( this->size_of( s1_, s2_, s3_, s4_, I5_.size() ) );
		} else if ( array_size_ != npos ) { // Unbounded with bounded data array
			size_type const slice_size( this->size_of( s1_, s2_, s3_, s4_ ) );
			if ( slice_size > 0 ) { // Infer upper index and size
				I5_.u( I5_.lz() + ( array_size_ / slice_size ) - 1 );
				size_set( this->size_of( slice_size, I5_.size() ) );
			} else {
				size_set( array_size_ );
			}
		} else { // Unbounded with unbounded data array
			size_set( npos );
		}
		shift_set( ( ( ( ( ( ( ( I5_.lz() * s4_ ) +  I4_.lz() ) * s3_ ) + I3_.lz() ) * s2_ ) + I2_.lz() ) * s1_ ) + I1_.lz() );
	}


private: // Data


	/// @brief Dimension 1 index range
	IR I1_;

	/// @brief Dimension 2 index range
	IR I2_;

	/// @brief Dimension 3 index range
	IR I3_;

	/// @brief Dimension 4 index range
	IR I4_;

	/// @brief Dimension 5 index range
	IR I5_;


}; // FArray5A


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray5A_HH
