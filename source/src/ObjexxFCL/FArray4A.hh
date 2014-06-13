#ifndef INCLUDED_ObjexxFCL_FArray4A_hh
#define INCLUDED_ObjexxFCL_FArray4A_hh


// FArray4A: Fortran-Compatible 4D Argument Array
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
#include <ObjexxFCL/FArray4A.fwd.hh>
#include <ObjexxFCL/FArray4P.hh>
#include <ObjexxFCL/StaticIndexRange.hh>


namespace ObjexxFCL {


/// @brief FArray4A: Fortran-Compatible 4D Argument Array
template< typename T >
class FArray4A :
	public FArray4< T >
{


private: // Types


	typedef  FArray4< T >  Super;
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


public: // Creation


	/// @brief Default Constructor
	inline
	FArray4A() :
		Super( ProxySentinel() )
	{}


	/// @brief Copy Constructor
	inline
	FArray4A( FArray4A const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ ),
		I4_( a.I4_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
	}


	/// @brief Non-Const Copy Constructor
	inline
	FArray4A( FArray4A & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ ),
		I4_( a.I4_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
	}


	/// @brief Real Constructor
	inline
	FArray4A( real_FArray const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ ),
		I4_( a.I4_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
	}


	/// @brief Non-Const Real Constructor
	inline
	FArray4A( real_FArray & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ ),
		I4_( a.I4_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
	}


	/// @brief Proxy Constructor
	inline
	FArray4A( proxy_FArray const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ ),
		I4_( a.I4_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
	}


	/// @brief Non-Const Proxy Constructor
	inline
	FArray4A( proxy_FArray & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ ),
		I4_( a.I4_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
	}


	/// @brief Super Constructor
	inline
	FArray4A( Super const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1() ),
		I2_( a.I2() ),
		I3_( a.I3() ),
		I4_( a.I4() )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
	}


	/// @brief Non-Const Super Constructor
	inline
	FArray4A( Super & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1() ),
		I2_( a.I2() ),
		I3_( a.I3() ),
		I4_( a.I4() )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
	}


	/// @brief Base Constructor
	inline
	FArray4A( Base const & a ) :
		Super( a, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( 1 ),
		I4_( a.size() )
	{
		shift_set( 4 );
		s1_ = s2_ = s3_ = 1;
	}


	/// @brief Non-Const Base Constructor
	inline
	FArray4A( Base & a ) :
		Super( a, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( 1 ),
		I4_( a.size() )
	{
		shift_set( 4 );
		s1_ = s2_ = s3_ = 1;
	}


	/// @brief Section Constructor
	inline
	FArray4A( Section const & s ) :
		Super( s, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( 1 ),
		I4_( s.size() )
	{
		shift_set( 4 );
		s1_ = s2_ = s3_ = 1;
	}


	/// @brief Non-Const Section Constructor
	inline
	FArray4A( Section & s ) :
		Super( s, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( 1 ),
		I4_( s.size() )
	{
		shift_set( 4 );
		s1_ = s2_ = s3_ = 1;
	}


	/// @brief Value Constructor
	inline
	FArray4A( T const & t ) :
		Super( t, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( 1 ),
		I4_( star ) // Unbounded
	{
		shift_set( 4 );
		s1_ = s2_ = s3_ = 1;
	}


	/// @brief Non-Const Value Constructor
	inline
	FArray4A( T & t ) :
		Super( t, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( 1 ),
		I4_( star ) // Unbounded
	{
		shift_set( 4 );
		s1_ = s2_ = s3_ = 1;
	}


	/// @brief Copy + IndexRange Constructor
	inline
	FArray4A( FArray4A const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Copy + IndexRange Constructor
	inline
	FArray4A( FArray4A & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Real + IndexRange Constructor
	inline
	FArray4A( real_FArray const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Real + IndexRange Constructor
	inline
	FArray4A( real_FArray & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Proxy + IndexRange Constructor
	inline
	FArray4A( proxy_FArray const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Proxy + IndexRange Constructor
	inline
	FArray4A( proxy_FArray & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Super + IndexRange Constructor
	inline
	FArray4A( Super const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Super + IndexRange Constructor
	inline
	FArray4A( Super & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Base + IndexRange Constructor
	inline
	FArray4A( Base const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Base + IndexRange Constructor
	inline
	FArray4A( Base & a, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Section + IndexRange Constructor
	inline
	FArray4A( Section const & s, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( s, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Section + IndexRange Constructor
	inline
	FArray4A( Section & s, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( s, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Value + IndexRange Constructor
	inline
	FArray4A( T const & t, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( t, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Value + IndexRange Constructor
	inline
	FArray4A( T & t, IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) :
		Super( t, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a ),
		I4_( I4_a )
	{
		dimension_argument();
	}


	/// @brief Destructor
	inline
	virtual
	~FArray4A()
	{}


public: // Assignment


	/// @brief Copy Assignment
	inline
	FArray4A &
	operator =( FArray4A const & a )
	{
		if ( this != &a ) {
			if ( ! equal_dimension( a ) ) dimension( a );
			Base::operator =( a );
		}
		return *this;
	}


	/// @brief Super Assignment
	inline
	FArray4A &
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
	FArray4A &
	operator =( FArray4< U > const & a )
	{
		if ( ! equal_dimension( a ) ) dimension( a );
		Base::operator =( a );
		return *this;
	}


	/// @brief += Array Template
	template< typename U >
	inline
	FArray4A &
	operator +=( FArray4< U > const & a )
	{
		Super::operator +=( a );
		return *this;
	}


	/// @brief -= Array Template
	template< typename U >
	inline
	FArray4A &
	operator -=( FArray4< U > const & a )
	{
		Super::operator -=( a );
		return *this;
	}


	/// @brief = Value
	inline
	FArray4A &
	operator =( T const & t )
	{
		Super::operator =( t );
		return *this;
	}


	/// @brief += Value
	inline
	FArray4A &
	operator +=( T const & t )
	{
		Super::operator +=( t );
		return *this;
	}


	/// @brief -= Value
	inline
	FArray4A &
	operator -=( T const & t )
	{
		Super::operator -=( t );
		return *this;
	}


	/// @brief *= Value
	inline
	FArray4A &
	operator *=( T const & t )
	{
		Super::operator *=( t );
		return *this;
	}


	/// @brief /= Value
	inline
	FArray4A &
	operator /=( T const & t )
	{
		Super::operator /=( t );
		return *this;
	}


public: // Subscript


	/// @brief array( i1, i2, i3, i4 ) const
	inline
	T const &
	operator ()( int const i1, int const i2, int const i3, int const i4 ) const
	{
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) && ( I4_.contains( i4 ) ) );
		return sarray_[ ( ( ( ( ( i4 * s3_ ) + i3 ) * s2_ ) + i2 ) * s1_ ) + i1 ];
	}


	/// @brief array( i1, i2, i3, i4 )
	inline
	T &
	operator ()( int const i1, int const i2, int const i3, int const i4 )
	{
		proxy_const_assert( not_const_proxy() );
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) && ( I4_.contains( i4 ) ) );
		return sarray_[ ( ( ( ( ( i4 * s3_ ) + i3 ) * s2_ ) + i2 ) * s1_ ) + i1 ];
	}


	/// @brief Section Starting at array( i1, i2, i3, i4 )
	inline
	Section const
	a( int const i1, int const i2, int const i3, int const i4 ) const
	{
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) && ( I4_.contains( i4 ) ) );
		size_type const offset( ( ( ( ( ( ( i4 * s3_ ) + i3 ) * s2_ ) + i2 ) * s1_ ) + i1 ) - shift_ );
		return Section( static_cast< T const * >( array_ + offset ), ( array_size_ != npos ) ? array_size_ - offset : npos );
	}


	/// @brief Non-Const Section Starting at array( i1, i2, i3, i4 )
	inline
	Section
	a( int const i1, int const i2, int const i3, int const i4 )
	{
		proxy_const_assert( not_const_proxy() );
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) && ( I4_.contains( i4 ) ) );
		size_type const offset( ( ( ( ( ( ( i4 * s3_ ) + i3 ) * s2_ ) + i2 ) * s1_ ) + i1 ) - shift_ );
		return Section( array_ + offset, ( array_size_ != npos ) ? array_size_ - offset : npos );
	}


	/// @brief Linear Index
	inline
	size_type
	index( int const i1, int const i2, int const i3, int const i4 ) const
	{
		return ( ( ( ( ( ( ( i4 * s3_ ) + i3 ) * s2_ ) + i2 ) * s1_ ) + i1 ) - shift_ );
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
	contains( int const i1, int const i2, int const i3, int const i4 ) const
	{
		return ( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) && ( I4_.contains( i4 ) ) );
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


	/// @brief Size of Dimension 4
	inline
	size_type
	size4() const
	{
		return I4_.size();
	}


public: // Modifier


	/// @brief Clear
	inline
	FArray4A &
	clear()
	{
		Super::clear();
		I1_.clear();
		I2_.clear();
		I3_.clear();
		I4_.clear();
		return *this;
	}


	/// @brief Dimension by IndexRanges Even if Const
	inline
	FArray4A const &
	dim( IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a ) const
	{
		const_cast< FArray4A & >( *this ).dimension( I1_a, I2_a, I3_a, I4_a );
		return *this;
	}


	/// @brief Dimension by Array Even if Const
	template< typename U >
	inline
	FArray4A const &
	dim( FArray4< U > const & a ) const
	{
		const_cast< FArray4A & >( *this ).dimension( a );
		return *this;
	}


	/// @brief Dimension by IndexRange
	inline
	FArray4A &
	dimension( IR const & I1_a, IR const & I2_a, IR const & I3_a, IR const & I4_a )
	{
		I1_.assign_value_of( I1_a );
		I2_.assign_value_of( I2_a );
		I3_.assign_value_of( I3_a );
		I4_.assign_value_of( I4_a );
		dimension_argument();
		return *this;
	}


	/// @brief Dimension by Array
	template< typename U >
	inline
	FArray4A &
	dimension( FArray4< U > const & a )
	{
		I1_.assign_value_of( a.I1() );
		I2_.assign_value_of( a.I2() );
		I3_.assign_value_of( a.I3() );
		I4_.assign_value_of( a.I4() );
		dimension_argument();
		return *this;
	}


	/// @brief Attach to Argument Array
	inline
	FArray4A &
	attach( FArray4A const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		I4_.assign_value_of( a.I4_ );
		return *this;
	}


	/// @brief Attach to Non-Const Argument Array
	inline
	FArray4A &
	attach( FArray4A & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		I4_.assign_value_of( a.I4_ );
		return *this;
	}


	/// @brief Attach to Real Array
	inline
	FArray4A &
	attach( real_FArray const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		I4_.assign_value_of( a.I4_ );
		return *this;
	}


	/// @brief Attach to Non-Const Real Array
	inline
	FArray4A &
	attach( real_FArray & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		I4_.assign_value_of( a.I4_ );
		return *this;
	}


	/// @brief Attach to Proxy Array
	inline
	FArray4A &
	attach( proxy_FArray const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		I4_.assign_value_of( a.I4_ );
		return *this;
	}


	/// @brief Attach to Non-Const Proxy Array
	inline
	FArray4A &
	attach( proxy_FArray & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		I4_.assign_value_of( a.I4_ );
		return *this;
	}


	/// @brief Attach to Super Array
	inline
	FArray4A &
	attach( Super const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		I1_.assign_value_of( a.I1() );
		I2_.assign_value_of( a.I2() );
		I3_.assign_value_of( a.I3() );
		I4_.assign_value_of( a.I4() );
		return *this;
	}


	/// @brief Attach to Non-Const Super Array
	inline
	FArray4A &
	attach( Super & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		s3_ = a.s3_;
		I1_.assign_value_of( a.I1() );
		I2_.assign_value_of( a.I2() );
		I3_.assign_value_of( a.I3() );
		I4_.assign_value_of( a.I4() );
		return *this;
	}


	/// @brief Attach to Base Array
	inline
	FArray4A &
	attach( Base const & a )
	{
		Base::attach( a, 4 );
		s1_ = s2_ = s3_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = 1;
		I4_ = a.size();
		return *this;
	}


	/// @brief Attach to Non-Const Base Array
	inline
	FArray4A &
	attach( Base & a )
	{
		Base::attach( a, 4 );
		s1_ = s2_ = s3_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = 1;
		I4_ = a.size();
		return *this;
	}


	/// @brief Attach to Section
	inline
	FArray4A &
	attach( Section const & s )
	{
		Base::attach( s, 4 );
		s1_ = s2_ = s3_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = 1;
		I4_ = s.size();
		return *this;
	}

	/// @brief Attach to Non-Const Section
	inline
	FArray4A &
	attach( Section & s )
	{
		Base::attach( s, 4 );
		s1_ = s2_ = s3_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = 1;
		I4_ = s.size();
		return *this;
	}


	/// @brief Attach to Value
	inline
	FArray4A &
	attach( T const & t )
	{
		Base::attach( t, 4 );
		s1_ = s2_ = s3_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = 1;
		I4_ = star; // Unbounded
		return *this;
	}


	/// @brief Attach to Non-Const Value
	inline
	FArray4A &
	attach( T & t )
	{
		Base::attach( t, 4 );
		s1_ = s2_ = s3_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = 1;
		I4_ = star; // Unbounded
		return *this;
	}


	/// @brief Detach from Source Array
	inline
	FArray4A &
	detach()
	{
		Base::detach();
		s1_ = s2_ = s3_ = 0;
		I1_.clear();
		I2_.clear();
		I3_.clear();
		I4_.clear();
		return *this;
	}


protected: // Functions


	/// @brief Dimension by IndexRanges
	inline
	void
	dimension_assign( SIR const & I1_a, SIR const & I2_a, SIR const & I3_a, SIR const & I4_a )
	{
		I1_.assign_value_of( I1_a );
		I2_.assign_value_of( I2_a );
		I3_.assign_value_of( I3_a );
		I4_.assign_value_of( I4_a );
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
		s1_ = I1_.size();
		s2_ = I2_.size();
		s3_ = I3_.size();
		if ( I4_.bounded_value() ) { // Bounded
			size_set( this->size_of( s1_, s2_, s3_, I4_.size() ) );
		} else if ( array_size_ != npos ) { // Unbounded with bounded data array
			size_type const slice_size( this->size_of( s1_, s2_, s3_ ) );
			if ( slice_size > 0 ) { // Infer upper index and size
				I4_.u( I4_.lz() + ( array_size_ / slice_size ) - 1 );
				size_set( this->size_of( slice_size, I4_.size() ) );
			} else {
				size_set( array_size_ );
			}
		} else { // Unbounded with unbounded data array
			size_set( npos );
		}
		shift_set( ( ( ( ( ( I4_.lz() * s3_ ) + I3_.lz() ) * s2_ ) + I2_.lz() ) * s1_ ) + I1_.lz() );
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


}; // FArray4A


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray4A_HH
