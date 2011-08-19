#ifndef INCLUDED_ObjexxFCL_FArray1A_hh
#define INCLUDED_ObjexxFCL_FArray1A_hh


// FArray1A: Fortran-Compatible 1D Argument Array
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
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1P.hh>
#include <ObjexxFCL/StaticIndexRange.hh>


namespace ObjexxFCL {


/// @brief FArray1A: Fortran-Compatible 1D Argument Array
template< typename T >
class FArray1A :
	public FArray1< T >
{


private: // Types


	typedef  FArray1< T >  Super;
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


public: // Creation


	/// @brief Default Constructor
	inline
	FArray1A() :
		Super( ProxySentinel() )
	{}


	/// @brief Copy Constructor
	inline
	FArray1A( FArray1A const & a ) :
		Super( a, ProxySentinel() ),
		I_( a.I_ )
	{
		shift_set( a.shift_ );
	}


	/// @brief Non-Const Copy Constructor
	inline
	FArray1A( FArray1A & a ) :
		Super( a, ProxySentinel() ),
		I_( a.I_ )
	{
		shift_set( a.shift_ );
	}


	/// @brief Real Constructor
	inline
	FArray1A( real_FArray const & a ) :
		Super( a, ProxySentinel() ),
		I_( a.I_ )
	{
		shift_set( a.shift_ );
	}


	/// @brief Non-Const Real Constructor
	inline
	FArray1A( real_FArray & a ) :
		Super( a, ProxySentinel() ),
		I_( a.I_ )
	{
		shift_set( a.shift_ );
	}


	/// @brief Proxy Constructor
	inline
	FArray1A( proxy_FArray const & a ) :
		Super( a, ProxySentinel() ),
		I_( a.I_ )
	{
		shift_set( a.shift_ );
	}


	/// @brief Non-Const Proxy Constructor
	inline
	FArray1A( proxy_FArray & a ) :
		Super( a, ProxySentinel() ),
		I_( a.I_ )
	{
		shift_set( a.shift_ );
	}


	/// @brief Super Constructor
	inline
	FArray1A( Super const & a ) :
		Super( a, ProxySentinel() ),
		I_( a.I() )
	{
		shift_set( a.shift_ );
	}


	/// @brief Non-Const Super Constructor
	inline
	FArray1A( Super & a ) :
		Super( a, ProxySentinel() ),
		I_( a.I() )
	{
		shift_set( a.shift_ );
	}


	/// @brief Base Constructor
	inline
	FArray1A( Base const & a ) :
		Super( a, ProxySentinel() ),
		I_( a.size() )
	{
		shift_set( 1 );
	}


	/// @brief Non-Const Base Constructor
	inline
	FArray1A( Base & a ) :
		Super( a, ProxySentinel() ),
		I_( a.size() )
	{
		shift_set( 1 );
	}


	/// @brief Section Constructor
	inline
	FArray1A( Section const & s ) :
		Super( s, ProxySentinel() ),
		I_( s.size() )
	{
		shift_set( 1 );
	}


	/// @brief Non-Const Section Constructor
	inline
	FArray1A( Section & s ) :
		Super( s, ProxySentinel() ),
		I_( s.size() )
	{
		shift_set( 1 );
	}


	/// @brief Value Constructor
	inline
	FArray1A( T const & t ) :
		Super( t, ProxySentinel() ),
		I_( star ) // Unbounded
	{
		shift_set( 1 );
	}


	/// @brief Non-Const Value Constructor
	inline
	FArray1A( T & t ) :
		Super( t, ProxySentinel() ),
		I_( star ) // Unbounded
	{
		shift_set( 1 );
	}


	/// @brief Copy + IndexRange Constructor
	inline
	FArray1A( FArray1A const & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Copy + IndexRange Constructor
	inline
	FArray1A( FArray1A & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Real + IndexRange Constructor
	inline
	FArray1A( real_FArray const & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Real + IndexRange Constructor
	inline
	FArray1A( real_FArray & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Proxy + IndexRange Constructor
	inline
	FArray1A( proxy_FArray const & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Proxy + IndexRange Constructor
	inline
	FArray1A( proxy_FArray & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Super + IndexRange Constructor
	inline
	FArray1A( Super const & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Super + IndexRange Constructor
	inline
	FArray1A( Super & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Base + IndexRange Constructor
	inline
	FArray1A( Base const & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Base + IndexRange Constructor
	inline
	FArray1A( Base & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Section + IndexRange Constructor
	inline
	FArray1A( Section const & s, IR const & I_a ) :
		Super( s, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Section + IndexRange Constructor
	inline
	FArray1A( Section & s, IR const & I_a ) :
		Super( s, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Value + IndexRange Constructor
	inline
	FArray1A( T const & t, IR const & I_a ) :
		Super( t, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Non-Const Value + IndexRange Constructor
	inline
	FArray1A( T & t, IR const & I_a ) :
		Super( t, ProxySentinel() ),
		I_( I_a )
	{
		dimension_argument();
	}


	/// @brief Destructor
	inline
	virtual
	~FArray1A()
	{}


public: // Assignment


	/// @brief Copy Assignment
	inline
	FArray1A &
	operator =( FArray1A const & a )
	{
		if ( this != &a ) {
			if ( ! this->equal_dimension( a ) ) dimension( a );
			Base::operator =( a );
		}
		return *this;
	}


	/// @brief Super Assignment
	inline
	FArray1A &
	operator =( Super const & a )
	{
		if ( this != &a ) {
			if ( ! this->equal_dimension( a ) ) dimension( a );
			Base::operator =( a );
		}
		return *this;
	}


	/// @brief Super Assignment Template
	template< typename U >
	inline
	FArray1A &
	operator =( FArray1< U > const & a )
	{
		if ( ! equal_dimension( a ) ) dimension( a );
		Base::operator =( a );
		return *this;
	}


	/// @brief += Array Template
	template< typename U >
	inline
	FArray1A &
	operator +=( FArray1< U > const & a )
	{
		Super::operator +=( a );
		return *this;
	}


	/// @brief -= Array Template
	template< typename U >
	inline
	FArray1A &
	operator -=( FArray1< U > const & a )
	{
		Super::operator -=( a );
		return *this;
	}


	/// @brief = Value
	inline
	FArray1A &
	operator =( T const & t )
	{
		Super::operator =( t );
		return *this;
	}


	/// @brief += Value
	inline
	FArray1A &
	operator +=( T const & t )
	{
		Super::operator +=( t );
		return *this;
	}


	/// @brief -= Value
	inline
	FArray1A &
	operator -=( T const & t )
	{
		Super::operator -=( t );
		return *this;
	}


	/// @brief *= Value
	inline
	FArray1A &
	operator *=( T const & t )
	{
		Super::operator *=( t );
		return *this;
	}


	/// @brief /= Value
	inline
	FArray1A &
	operator /=( T const & t )
	{
		Super::operator /=( t );
		return *this;
	}


public: // Subscript


	/// @brief array( i ) const
	inline
	T const &
	operator ()( int const i ) const
	{
		assert( I_.contains( i ) );
		return sarray_[ i ];
	}


	/// @brief array( i )
	inline
	T &
	operator ()( int const i )
	{
		proxy_const_assert( not_const_proxy() );
		assert( I_.contains( i ) );
		return sarray_[ i ];
	}


	/// @brief Section Starting at array( i )
	inline
	Section const
	a( int const i ) const
	{
		assert( I_.contains( i ) );
		return Section( static_cast< T const * >( sarray_ + i ), ( array_size_ != npos ) ? array_size_ - ( i - shift_ ) : npos );
	}


	/// @brief Non-Const Section Starting at array( i )
	inline
	Section
	a( int const i )
	{
		proxy_const_assert( not_const_proxy() );
		assert( I_.contains( i ) );
		return Section( sarray_ + i, ( array_size_ != npos ) ? array_size_ - ( i - shift_ ) : npos );
	}


	/// @brief Linear Index
	inline
	size_type
	index( int const i ) const
	{
		return ( i - shift_ );
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
	contains( int const i ) const
	{
		return I_.contains( i );
	}


	/// @brief Initializer Active?
	inline
	bool
	initializer_active() const
	{
		return false;
	}


public: // Inspector


	/// @brief IndexRange
	inline
	IR const &
	I1() const
	{
		return I_;
	}


	/// @brief Lower Index
	inline
	int
	l1() const
	{
		return I_.l();
	}


	/// @brief Upper Index
	inline
	int
	u1() const
	{
		return I_.u();
	}


	/// @brief Size
	inline
	size_type
	size1() const
	{
		return I_.size();
	}


	/// @brief IndexRange
	inline
	IR const &
	I() const
	{
		return I_;
	}


	/// @brief Lower Index
	inline
	int
	l() const
	{
		return I_.l();
	}


	/// @brief Upper Index
	inline
	int
	u() const
	{
		return I_.u();
	}


public: // Modifier


	/// @brief Clear
	inline
	FArray1A &
	clear()
	{
		Super::clear();
		I_.clear();
		return *this;
	}


	/// @brief Dimension by IndexRanges Even if Const
	inline
	FArray1A const &
	dim( IR const & I_a ) const
	{
		const_cast< FArray1A & >( *this ).dimension( I_a );
		return *this;
	}


	/// @brief Dimension by Array Even if Const
	template< typename U >
	inline
	FArray1A const &
	dim( FArray1< U > const & a ) const
	{
		const_cast< FArray1A & >( *this ).dimension( a );
		return *this;
	}


	/// @brief Dimension by IndexRanges
	inline
	FArray1A &
	dimension( IR const & I_a )
	{
		I_.assign_value_of( I_a );
		dimension_argument();
		return *this;
	}


	/// @brief Dimension by Array
	template< typename U >
	inline
	FArray1A &
	dimension( FArray1< U > const & a )
	{
		I_.assign_value_of( a.I() );
		dimension_argument();
		return *this;
	}


	/// @brief Attach to Argument Array
	inline
	FArray1A &
	attach( FArray1A const & a )
	{
		Base::attach( a );
		I_.assign_value_of( a.I_ );
		return *this;
	}


	/// @brief Attach to Non-Const Argument Array
	inline
	FArray1A &
	attach( FArray1A & a )
	{
		Base::attach( a );
		I_.assign_value_of( a.I_ );
		return *this;
	}


	/// @brief Attach to Real Array
	inline
	FArray1A &
	attach( real_FArray const & a )
	{
		Base::attach( a );
		I_.assign_value_of( a.I_ );
		return *this;
	}


	/// @brief Attach to Non-Const Real Array
	inline
	FArray1A &
	attach( real_FArray & a )
	{
		Base::attach( a );
		I_.assign_value_of( a.I_ );
		return *this;
	}


	/// @brief Attach to Proxy Array
	inline
	FArray1A &
	attach( proxy_FArray const & a )
	{
		Base::attach( a );
		I_.assign_value_of( a.I_ );
		return *this;
	}


	/// @brief Attach to Non-Const Proxy Array
	inline
	FArray1A &
	attach( proxy_FArray & a )
	{
		Base::attach( a );
		I_.assign_value_of( a.I_ );
		return *this;
	}


	/// @brief Attach to Super Array
	inline
	FArray1A &
	attach( Super const & a )
	{
		Base::attach( a );
		I_.assign_value_of( a.I() );
		return *this;
	}


	/// @brief Attach to Non-Const Super Array
	inline
	FArray1A &
	attach( Super & a )
	{
		Base::attach( a );
		I_.assign_value_of( a.I() );
		return *this;
	}


	/// @brief Attach to Base Array
	inline
	FArray1A &
	attach( Base const & a )
	{
		Base::attach( a, 1 );
		I_ = a.size();
		return *this;
	}


	/// @brief Attach to Non-Const Base Array
	inline
	FArray1A &
	attach( Base & a )
	{
		Base::attach( a, 1 );
		I_ = a.size();
		return *this;
	}


	/// @brief Attach to Section
	inline
	FArray1A &
	attach( Section const & s )
	{
		Base::attach( s, 1 );
		I_ = s.size();
		return *this;
	}


	/// @brief Attach to Non-Const Section
	inline
	FArray1A &
	attach( Section & s )
	{
		Base::attach( s, 1 );
		I_ = s.size();
		return *this;
	}


	/// @brief Attach to Value
	inline
	FArray1A &
	attach( T const & t )
	{
		Base::attach( t, 1 );
		I_ = star; // Unbounded
		return *this;
	}


	/// @brief Attach to Non-Const Value
	inline
	FArray1A &
	attach( T & t )
	{
		Base::attach( t, 1 );
		I_ = star; // Unbounded
		return *this;
	}


	/// @brief Detach from Source Array
	inline
	FArray1A &
	detach()
	{
		Base::detach();
		I_.clear();
		return *this;
	}


protected: // Functions


	/// @brief Dimension by IndexRanges
	inline
	void
	dimension_assign( SIR const & I_a )
	{
		I_.assign_value_of( I_a );
		dimension_argument();
	}


private: // Functions


	/// @brief Dimension by Current IndexRanges
	inline
	void
	dimension_argument()
	{
		if ( I_.bounded_value() ) { // Bounded
			size_set( I_.size() );
		} else if ( array_size_ != npos ) { // Unbounded with bounded data array
			// Infer upper index and size
			I_.u( I_.lz() + array_size_ - 1 );
			size_set( I_.size() );
		} else { // Unbounded with unbounded data array
			size_set( npos );
		}
		shift_set( I_.lz() );
	}


private: // Data


	/// @brief Index range
	IR I_;


}; // FArray1A


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray1A_HH
