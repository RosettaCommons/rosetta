#ifndef INCLUDED_ObjexxFCL_FArray3A_hh
#define INCLUDED_ObjexxFCL_FArray3A_hh


// FArray3A: Fortran-Compatible 3D Argument Array
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
#include <ObjexxFCL/FArray3A.fwd.hh>
#include <ObjexxFCL/FArray3P.hh>
#include <ObjexxFCL/StaticIndexRange.hh>


namespace ObjexxFCL {


/// @brief FArray3A: Fortran-Compatible 3D Argument Array
template< typename T >
class FArray3A :
	public FArray3< T >
{


private: // Types


	typedef  FArray3< T >  Super;
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


public: // Creation


	/// @brief Default Constructor
	inline
	FArray3A() :
		Super( ProxySentinel() )
	{}


	/// @brief Copy Constructor
	inline
	FArray3A( FArray3A const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
	}


	/// @brief Non-Const Copy Constructor
	inline
	FArray3A( FArray3A & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
	}


	/// @brief Real Constructor
	inline
	FArray3A( real_FArray const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
	}


	/// @brief Non-Const Real Constructor
	inline
	FArray3A( real_FArray & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
	}


	/// @brief Proxy Constructor
	inline
	FArray3A( proxy_FArray const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
	}


	/// @brief Non-Const Proxy Constructor
	inline
	FArray3A( proxy_FArray & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		I3_( a.I3_ )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
	}


	/// @brief Super Constructor
	inline
	FArray3A( Super const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1() ),
		I2_( a.I2() ),
		I3_( a.I3() )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
	}


	/// @brief Non-Const Super Constructor
	inline
	FArray3A( Super & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1() ),
		I2_( a.I2() ),
		I3_( a.I3() )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		s2_ = a.s2_;
	}


	/// @brief Base Constructor
	inline
	FArray3A( Base const & a ) :
		Super( a, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( a.size() )
	{
		shift_set( 3 );
		s1_ = s2_ = 1;
	}


	/// @brief Non-Const Base Constructor
	inline
	FArray3A( Base & a ) :
		Super( a, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( a.size() )
	{
		shift_set( 3 );
		s1_ = s2_ = 1;
	}


	/// @brief Section Constructor
	inline
	FArray3A( Section const & s ) :
		Super( s, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( s.size() )
	{
		shift_set( 3 );
		s1_ = s2_ = 1;
	}


	/// @brief Non-Const Section Constructor
	inline
	FArray3A( Section & s ) :
		Super( s, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( s.size() )
	{
		shift_set( 3 );
		s1_ = s2_ = 1;
	}


	/// @brief Value Constructor
	inline
	FArray3A( T const & t ) :
		Super( t, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( star ) // Unbounded
	{
		shift_set( 3 );
		s1_ = s2_ = 1;
	}


	/// @brief Non-Const Value Constructor
	inline
	FArray3A( T & t ) :
		Super( t, ProxySentinel() ),
		I1_( 1 ),
		I2_( 1 ),
		I3_( star ) // Unbounded
	{
		shift_set( 3 );
		s1_ = s2_ = 1;
	}


	/// @brief Copy + IndexRange Constructor
	inline
	FArray3A( FArray3A const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Copy + IndexRange Constructor
	inline
	FArray3A( FArray3A & a, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Real + IndexRange Constructor
	inline
	FArray3A( real_FArray const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Real + IndexRange Constructor
	inline
	FArray3A( real_FArray & a, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Proxy + IndexRange Constructor
	inline
	FArray3A( proxy_FArray const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Proxy + IndexRange Constructor
	inline
	FArray3A( proxy_FArray & a, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Super + IndexRange Constructor
	inline
	FArray3A( Super const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Super + IndexRange Constructor
	inline
	FArray3A( Super & a, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Base + IndexRange Constructor
	inline
	FArray3A( Base const & a, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Base + IndexRange Constructor
	inline
	FArray3A( Base & a, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Section + IndexRange Constructor
	inline
	FArray3A( Section const & s, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( s, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Section + IndexRange Constructor
	inline
	FArray3A( Section & s, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( s, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Value + IndexRange Constructor
	inline
	FArray3A( T const & t, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( t, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Value + IndexRange Constructor
	inline
	FArray3A( T & t, IR const & I1_a, IR const & I2_a, IR const & I3_a ) :
		Super( t, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		I3_( I3_a )
	{
		dimension_argument();
	}


	/// @brief Destructor
	inline
	virtual
	~FArray3A()
	{}


public: // Assignment


	/// @brief Copy Assignment
	inline
	FArray3A &
	operator =( FArray3A const & a )
	{
		if ( this != &a ) {
			if ( ! equal_dimension( a ) ) dimension( a );
			Base::operator =( a );
		}
		return *this;
	}


	/// @brief Super Assignment
	inline
	FArray3A &
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
	FArray3A &
	operator =( FArray3< U > const & a )
	{
		if ( ! equal_dimension( a ) ) dimension( a );
		Base::operator =( a );
		return *this;
	}


	/// @brief += Array Template
	template< typename U >
	inline
	FArray3A &
	operator +=( FArray3< U > const & a )
	{
		Super::operator +=( a );
		return *this;
	}


	/// @brief -= Array Template
	template< typename U >
	inline
	FArray3A &
	operator -=( FArray3< U > const & a )
	{
		Super::operator -=( a );
		return *this;
	}


	/// @brief = Value
	inline
	FArray3A &
	operator =( T const & t )
	{
		Super::operator =( t );
		return *this;
	}


	/// @brief += Value
	inline
	FArray3A &
	operator +=( T const & t )
	{
		Super::operator +=( t );
		return *this;
	}


	/// @brief -= Value
	inline
	FArray3A &
	operator -=( T const & t )
	{
		Super::operator -=( t );
		return *this;
	}


	/// @brief *= Value
	inline
	FArray3A &
	operator *=( T const & t )
	{
		Super::operator *=( t );
		return *this;
	}


	/// @brief /= Value
	inline
	FArray3A &
	operator /=( T const & t )
	{
		Super::operator /=( t );
		return *this;
	}


public: // Subscript


	/// @brief array( i1, i2, i3 ) const
	inline
	T const &
	operator ()( int const i1, int const i2, int const i3 ) const
	{
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) );
		return sarray_[ ( ( ( i3 * s2_ ) + i2 ) * s1_ ) + i1 ];
	}


	/// @brief array( i1, i2, i3 )
	inline
	T &
	operator ()( int const i1, int const i2, int const i3 )
	{
		proxy_const_assert( not_const_proxy() );
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) );
		return sarray_[ ( ( ( i3 * s2_ ) + i2 ) * s1_ ) + i1 ];
	}


	/// @brief Section Starting at array( i1, i2, i3 )
	inline
	Section const
	a( int const i1, int const i2, int const i3 ) const
	{
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) );
		size_type const offset( ( ( ( ( i3 * s2_ ) + i2 ) * s1_ ) + i1 ) - shift_ );
		return Section( static_cast< T const * >( array_ + offset ), ( array_size_ != npos ) ? array_size_ - offset : npos );
	}


	/// @brief Non-Const Section Starting at array( i1, i2, i3 )
	inline
	Section
	a( int const i1, int const i2, int const i3 )
	{
		proxy_const_assert( not_const_proxy() );
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) );
		size_type const offset( ( ( ( ( i3 * s2_ ) + i2 ) * s1_ ) + i1 ) - shift_ );
		return Section( array_ + offset, ( array_size_ != npos ) ? array_size_ - offset : npos );
	}


	/// @brief Linear Index
	inline
	size_type
	index( int const i1, int const i2, int const i3 ) const
	{
		return ( ( ( ( ( i3 * s2_ ) + i2 ) * s1_ ) + i1 ) - shift_ );
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
	contains( int const i1, int const i2, int const i3 ) const
	{
		return ( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) && ( I3_.contains( i3 ) ) );
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


	/// @brief Size of Dimension 3
	inline
	size_type
	size3() const
	{
		return I3_.size();
	}


public: // Modifier


	/// @brief Clear
	inline
	FArray3A &
	clear()
	{
		Super::clear();
		I1_.clear();
		I2_.clear();
		I3_.clear();
		return *this;
	}


	/// @brief Dimension by IndexRanges Even if Const
	inline
	FArray3A const &
	dim( IR const & I1_a, IR const & I2_a, IR const & I3_a ) const
	{
		const_cast< FArray3A & >( *this ).dimension( I1_a, I2_a, I3_a );
		return *this;
	}


	/// @brief Dimension by Array Even if Const
	template< typename U >
	inline
	FArray3A const &
	dim( FArray3< U > const & a ) const
	{
		const_cast< FArray3A & >( *this ).dimension( a );
		return *this;
	}


	/// @brief Dimension by IndexRange
	inline
	FArray3A &
	dimension( IR const & I1_a, IR const & I2_a, IR const & I3_a )
	{
		I1_.assign_value_of( I1_a );
		I2_.assign_value_of( I2_a );
		I3_.assign_value_of( I3_a );
		dimension_argument();
		return *this;
	}


	/// @brief Dimension by Array
	template< typename U >
	inline
	FArray3A &
	dimension( FArray3< U > const & a )
	{
		I1_.assign_value_of( a.I1() );
		I2_.assign_value_of( a.I2() );
		I3_.assign_value_of( a.I3() );
		dimension_argument();
		return *this;
	}


	/// @brief Attach to Argument Array
	inline
	FArray3A &
	attach( FArray3A const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		return *this;
	}


	/// @brief Attach to Non-Const Argument Array
	inline
	FArray3A &
	attach( FArray3A & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		return *this;
	}


	/// @brief Attach to Real Array
	inline
	FArray3A &
	attach( real_FArray const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		return *this;
	}


	/// @brief Attach to Non-Const Real Array
	inline
	FArray3A &
	attach( real_FArray & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		return *this;
	}


	/// @brief Attach to Proxy Array
	inline
	FArray3A &
	attach( proxy_FArray const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		return *this;
	}


	/// @brief Attach to Non-Const Proxy Array
	inline
	FArray3A &
	attach( proxy_FArray & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		I1_.assign_value_of( a.I1_ );
		I2_.assign_value_of( a.I2_ );
		I3_.assign_value_of( a.I3_ );
		return *this;
	}


	/// @brief Attach to Super Array
	inline
	FArray3A &
	attach( Super const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		I1_.assign_value_of( a.I1() );
		I2_.assign_value_of( a.I2() );
		I3_.assign_value_of( a.I3() );
		return *this;
	}


	/// @brief Attach to Non-Const Super Array
	inline
	FArray3A &
	attach( Super & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		s2_ = a.s2_;
		I1_.assign_value_of( a.I1() );
		I2_.assign_value_of( a.I2() );
		I3_.assign_value_of( a.I3() );
		return *this;
	}


	/// @brief Attach to Base Array
	inline
	FArray3A &
	attach( Base const & a )
	{
		Base::attach( a, 3 );
		s1_ = s2_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = a.size();
		return *this;
	}


	/// @brief Attach to Non-Const Base Array
	inline
	FArray3A &
	attach( Base & a )
	{
		Base::attach( a, 3 );
		s1_ = s2_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = a.size();
		return *this;
	}


	/// @brief Attach to Section
	inline
	FArray3A &
	attach( Section const & s )
	{
		Base::attach( s, 3 );
		s1_ = s2_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = s.size();
		return *this;
	}


	/// @brief Attach to Non-Const Section
	inline
	FArray3A &
	attach( Section & s )
	{
		Base::attach( s, 3 );
		s1_ = s2_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = s.size();
		return *this;
	}


	/// @brief Attach to Value
	inline
	FArray3A &
	attach( T const & t )
	{
		Base::attach( t, 3 );
		s1_ = s2_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = star; // Unbounded
		return *this;
	}


	/// @brief Attach to Non-Const Value
	inline
	FArray3A &
	attach( T & t )
	{
		Base::attach( t, 3 );
		s1_ = s2_ = 1;
		I1_ = 1;
		I2_ = 1;
		I3_ = star; // Unbounded
		return *this;
	}


	/// @brief Detach from Source Array
	inline
	FArray3A &
	detach()
	{
		Base::detach();
		s1_ = s2_ = 0;
		I1_.clear();
		I2_.clear();
		I3_.clear();
		return *this;
	}


protected: // Functions


	/// @brief Dimension by IndexRanges
	inline
	void
	dimension_assign( SIR const & I1_a, SIR const & I2_a, SIR const & I3_a )
	{
		I1_.assign_value_of( I1_a );
		I2_.assign_value_of( I2_a );
		I3_.assign_value_of( I3_a );
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
		s1_ = I1_.size();
		s2_ = I2_.size();
		if ( I3_.bounded_value() ) { // Bounded
			size_set( this->size_of( s1_, s2_, I3_.size() ) );
		} else if ( array_size_ != npos ) { // Unbounded with bounded data array
			size_type const slice_size( this->size_of( s1_, s2_ ) );
			if ( slice_size > 0 ) { // Infer upper index and size
				I3_.u( I3_.lz() + ( array_size_ / slice_size ) - 1 );
				size_set( this->size_of( slice_size, I3_.size() ) );
			} else {
				size_set( array_size_ );
			}
		} else { // Unbounded with unbounded data array
			size_set( npos );
		}
		shift_set( ( ( ( I3_.lz() * s2_ ) + I2_.lz() ) * s1_ ) + I1_.lz() );
	}


private: // Data


	/// @brief Dimension 1 index range
	IR I1_;

	/// @brief Dimension 2 index range
	IR I2_;

	/// @brief Dimension 3 index range
	IR I3_;


}; // FArray3A


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray3A_HH
