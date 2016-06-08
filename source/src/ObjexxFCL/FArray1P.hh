#ifndef INCLUDED_ObjexxFCL_FArray1P_hh
#define INCLUDED_ObjexxFCL_FArray1P_hh


// FArray1P: Fortran-Compatible 1D Proxy Array
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
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>


namespace ObjexxFCL {


/// @brief FArray1P: Fortran-Compatible 1D Proxy Array
template< typename T >
class FArray1P :
	public FArray1< T >,
	public ObserverMulti
{


private: // Types


	typedef  FArray1< T >  Super;
	typedef  typename Super::real_FArray  real_FArray;
	typedef  typename Super::proxy_FArray  proxy_FArray;
	typedef  internal::ProxySentinel  ProxySentinel;


private: // Friend


	friend class FArray1A< T >;


public: // Types


	typedef  typename Super::Base  Base;
	typedef  typename Base::Section  Section;
	typedef  typename Super::IR  SIR;
	typedef  DynamicIndexRange  IR;

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
	using Super::const_proxy;
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
	FArray1P() :
		Super( ProxySentinel() ),
		source_( 0 )
	{}


	/// @brief Copy Constructor
	inline
	FArray1P( FArray1P const & a ) :
		Super( a, ProxySentinel() ),
		ObserverMulti(),
		I_( a.I_ ),
		source_( &a )
	{
		shift_set( a.shift_ );
		insert_as_observer();
	}


	/// @brief Non-Const Copy Constructor
	inline
	FArray1P( FArray1P & a ) :
		Super( a, ProxySentinel() ),
		ObserverMulti(),
		I_( a.I_ ),
		source_( &a )
	{
		shift_set( a.shift_ );
		insert_as_observer();
	}


	/// @brief Real Constructor
	inline
	FArray1P( real_FArray const & a ) :
		Super( a, ProxySentinel() ),
		I_( a.I_ ),
		source_( &a )
	{
		shift_set( a.shift_ );
		insert_as_observer();
	}


	/// @brief Non-Const Real Constructor
	inline
	FArray1P( real_FArray & a ) :
		Super( a, ProxySentinel() ),
		I_( a.I_ ),
		source_( &a )
	{
		shift_set( a.shift_ );
		insert_as_observer();
	}


	/// @brief Super Constructor
	inline
	FArray1P( Super const & a ) :
		Super( a, ProxySentinel() ),
		I_( a.I() ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		shift_set( a.shift_ );
		insert_as_observer();
	}


	/// @brief Non-Const Super Constructor
	inline
	FArray1P( Super & a ) :
		Super( a, ProxySentinel() ),
		I_( a.I() ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		shift_set( a.shift_ );
		insert_as_observer();
	}


	/// @brief Base Constructor
	inline
	FArray1P( Base const & a ) :
		Super( a, ProxySentinel() ),
		I_( a.size() ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		shift_set( 1 );
		insert_as_observer();
	}


	/// @brief Non-Const Base Constructor
	inline
	FArray1P( Base & a ) :
		Super( a, ProxySentinel() ),
		I_( a.size() ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		shift_set( 1 );
		insert_as_observer();
	}


	/// @brief Section Constructor
	inline
	FArray1P( Section const & s ) :
		Super( s, ProxySentinel() ),
		I_( s.size() ),
		source_( 0 )
	{
		shift_set( 1 );
		insert_as_observer();
	}


	/// @brief Non-Const Section Constructor
	inline
	FArray1P( Section & s ) :
		Super( s, ProxySentinel() ),
		I_( s.size() ),
		source_( 0 )
	{
		shift_set( 1 );
		insert_as_observer();
	}


	/// @brief Value Constructor
	inline
	FArray1P( T const & t ) :
		Super( t, ProxySentinel() ),
		I_( star ), // Unbounded
		source_( 0 )
	{
		shift_set( 1 );
		insert_as_observer();
	}


	/// @brief Non-Const Value Constructor
	inline
	FArray1P( T & t ) :
		Super( t, ProxySentinel() ),
		I_( star ), // Unbounded
		source_( 0 )
	{
		shift_set( 1 );
		insert_as_observer();
	}


	/// @brief Copy + IndexRange Constructor
	inline
	FArray1P( FArray1P const & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a ),
		source_( &a )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Non-Const Copy + IndexRange Constructor
	inline
	FArray1P( FArray1P & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a ),
		source_( &a )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Real + IndexRange Constructor
	inline
	FArray1P( real_FArray const & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a ),
		source_( &a )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Non-Const Real + IndexRange Constructor
	inline
	FArray1P( real_FArray & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a ),
		source_( &a )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Super + IndexRange Constructor
	inline
	FArray1P( Super const & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Non-Const Super + IndexRange Constructor
	inline
	FArray1P( Super & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Base + IndexRange Constructor
	inline
	FArray1P( Base const & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Non-Const Base + IndexRange Constructor
	inline
	FArray1P( Base & a, IR const & I_a ) :
		Super( a, ProxySentinel() ),
		I_( I_a ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Section + IndexRange Constructor
	inline
	FArray1P( Section const & s, IR const & I_a ) :
		Super( s, ProxySentinel() ),
		I_( I_a ),
		source_( 0 )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Non-Const Section + IndexRange Constructor
	inline
	FArray1P( Section & s, IR const & I_a ) :
		Super( s, ProxySentinel() ),
		I_( I_a ),
		source_( 0 )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Value + IndexRange Constructor
	inline
	FArray1P( T const & t, IR const & I_a ) :
		Super( t, ProxySentinel() ),
		I_( I_a ),
		source_( 0 )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Non-Const Value + IndexRange Constructor
	inline
	FArray1P( T & t, IR const & I_a ) :
		Super( t, ProxySentinel() ),
		I_( I_a ),
		source_( 0 )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Destructor
	inline
	virtual
	~FArray1P()
	{
		if ( source_ ) source_->remove_observer( *this );
	}


public: // Assignment


	/// @brief Copy Assignment
	inline
	FArray1P &
	operator =( FArray1P const & a )
	{
		if ( this != &a ) {
			if ( ! equal_dimension( a ) ) dimension( a );
			Base::operator =( a );
		}
		return *this;
	}


	/// @brief Super Assignment
	inline
	FArray1P &
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
	FArray1P &
	operator =( FArray1< U > const & a )
	{
		if ( ! equal_dimension( a ) ) dimension( a );
		Base::operator =( a );
		return *this;
	}


	/// @brief += Array Template
	template< typename U >
	inline
	FArray1P &
	operator +=( FArray1< U > const & a )
	{
		Super::operator +=( a );
		return *this;
	}


	/// @brief -= Array Template
	template< typename U >
	inline
	FArray1P &
	operator -=( FArray1< U > const & a )
	{
		Super::operator -=( a );
		return *this;
	}


	/// @brief = Value
	inline
	FArray1P &
	operator =( T const & t )
	{
		Super::operator =( t );
		return *this;
	}


	/// @brief += Value
	inline
	FArray1P &
	operator +=( T const & t )
	{
		Super::operator +=( t );
		return *this;
	}


	/// @brief -= Value
	inline
	FArray1P &
	operator -=( T const & t )
	{
		Super::operator -=( t );
		return *this;
	}


	/// @brief *= Value
	inline
	FArray1P &
	operator *=( T const & t )
	{
		Super::operator *=( t );
		return *this;
	}


	/// @brief /= Value
	inline
	FArray1P &
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
		assert( I_.initialized() );
		return ( i - shift_ );
	}


public: // Predicate


	/// @brief Dimensions Initialized?
	inline
	bool
	dimensions_initialized() const
	{
		return I_.initialized();
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
	FArray1P &
	clear()
	{
		Super::clear();
		I_.clear_no_notify();
		source_ = 0;
		return *this;
	}


	/// @brief Dimension by IndexRanges Even if Const
	inline
	FArray1P const &
	dim( IR const & I_a ) const
	{
		const_cast< FArray1P & >( *this ).dimension( I_a );
		return *this;
	}


	/// @brief Dimension by Array Even if Const
	template< typename U >
	inline
	FArray1P const &
	dim( FArray1< U > const & a ) const
	{
		const_cast< FArray1P & >( *this ).dimension( a );
		return *this;
	}


	/// @brief Dimension by IndexRanges
	inline
	FArray1P &
	dimension( IR const & I_a )
	{
		I_.assign_no_notify( I_a );
		dimension_proxy();
		return *this;
	}


	/// @brief Dimension by Array
	template< typename U >
	inline
	FArray1P &
	dimension( FArray1< U > const & a )
	{
		I_.assign_no_notify( a.I() );
		dimension_proxy();
		return *this;
	}


	/// @brief Attach to Proxy Array
	inline
	FArray1P &
	attach( FArray1P const & a )
	{
		Base::attach( a );
		I_ = a.I_;
		if ( source_ ) source_->remove_observer( *this );
		source_ = &a;
		a.insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Non-Const Proxy Array
	inline
	FArray1P &
	attach( FArray1P & a )
	{
		Base::attach( a );
		I_ = a.I_;
		if ( source_ ) source_->remove_observer( *this );
		source_ = &a;
		a.insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Real Array
	inline
	FArray1P &
	attach( real_FArray const & a )
	{
		Base::attach( a );
		I_ = a.I_;
		if ( source_ ) source_->remove_observer( *this );
		source_ = &a;
		a.insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Non-Const Real Array
	inline
	FArray1P &
	attach( real_FArray & a )
	{
		Base::attach( a );
		I_ = a.I_;
		if ( source_ ) source_->remove_observer( *this );
		source_ = &a;
		a.insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Super Array
	inline
	FArray1P &
	attach( Super const & a )
	{
		Base::attach( a );
		I_ = a.I();
		if ( source_ ) source_->remove_observer( *this );
		source_ = dynamic_cast< SubjectMulti const * >( &a );
		if ( source_ ) source_->insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Non-Const Super Array
	inline
	FArray1P &
	attach( Super & a )
	{
		Base::attach( a );
		I_ = a.I();
		if ( source_ ) source_->remove_observer( *this );
		source_ = dynamic_cast< SubjectMulti const * >( &a );
		if ( source_ ) source_->insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Base Array
	inline
	FArray1P &
	attach( Base const & a )
	{
		Base::attach( a, 1 );
		I_ = a.size();
		if ( source_ ) source_->remove_observer( *this );
		source_ = dynamic_cast< SubjectMulti const * >( &a );
		if ( source_ ) source_->insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Non-Const Base Array
	inline
	FArray1P &
	attach( Base & a )
	{
		Base::attach( a, 1 );
		I_ = a.size();
		if ( source_ ) source_->remove_observer( *this );
		source_ = dynamic_cast< SubjectMulti const * >( &a );
		if ( source_ ) source_->insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Section
	inline
	FArray1P &
	attach( Section const & s )
	{
		Base::attach( s, 1 );
		I_ = s.size();
		if ( source_ ) source_->remove_observer( *this );
		source_ = 0;
		return *this;
	}


	/// @brief Attach to Non-Const Section
	inline
	FArray1P &
	attach( Section & s )
	{
		Base::attach( s, 1 );
		I_ = s.size();
		if ( source_ ) source_->remove_observer( *this );
		source_ = 0;
		return *this;
	}


	/// @brief Attach to Value
	inline
	FArray1P &
	attach( T const & t )
	{
		Base::attach( t, 1 );
		I_ = star; // Unbounded
		if ( source_ ) source_->remove_observer( *this );
		source_ = 0;
		return *this;
	}


	/// @brief Attach to Non-Const Value
	inline
	FArray1P &
	attach( T & t )
	{
		Base::attach( t, 1 );
		I_ = star; // Unbounded
		if ( source_ ) source_->remove_observer( *this );
		source_ = 0;
		return *this;
	}


	/// @brief Detach from Source Array
	inline
	FArray1P &
	detach()
	{
		Base::detach();
		I_.clear();
		if ( source_ ) source_->remove_observer( *this );
		source_ = 0;
		return *this;
	}


public: // Observer Modifier


	/// @brief Update
	inline
	void
	update()
	{
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		if ( source_ ) {
			if ( const_proxy() ) {
				this->update_to( *dynamic_cast< Base const * >( source_ ) );
			} else {
				this->update_to( *dynamic_cast< Base * >( const_cast< SubjectMulti * >( source_ ) ) );
			}
		}
#else
		if ( source_ ) this->update_to( *dynamic_cast< Base const * >( source_ ) );
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		dimension_proxy();
	}


	/// @brief Update for Destruction of a Subject
	inline
	void
	destructed( Subject const & subject )
	{
		if ( ( source_ ) && ( &subject == source_ ) ) { // Source array is being destructed
			Base::detach();
			I_.clear();
			source_ = 0;
		}
	}


protected: // Functions


	/// @brief Dimension by IndexRanges
	inline
	void
	dimension_assign( SIR const & I_a )
	{
		I_.assign_no_notify( I_a );
		dimension_proxy();
	}


private: // Functions


	/// @brief Dimension by Current IndexRanges
	inline
	void
	dimension_proxy()
	{
		if ( dimensions_initialized() ) {
			if ( I_.bounded_value() ) { // Bounded
				size_set( I_.size() );
			} else if ( array_size_ != npos ) { // Unbounded with bounded data array
				// Infer upper index and size
				I_.u_no_notify( I_.lz() + array_size_ - 1 );
				size_set( I_.size() );
			} else { // Unbounded with unbounded data array
				size_set( npos );
			}
			shift_set( I_.lz() );
		} else {
			size_set( 0 );
			shift_set( 0 );
		}
	}


	/// @brief Insert as Observer of the IndexRanges and Source Array
	inline
	void
	insert_as_observer()
	{
		I_.insert_observer( *this );
		if ( source_ ) source_->insert_observer( *this );
	}


	/* Unused private
	/// @brief Remove as Observer of the IndexRanges and Source Array
	inline
	void
	remove_as_observer()
	{
		I_.remove_observer( *this );
		if ( source_ ) source_->remove_observer( *this );
	}
	*/


private: // Data


	/// @brief Index range
	IR I_;

	/// @brief Pointer (non-owning) to source array (0 if unknown)
	SubjectMulti const * source_;


}; // FArray1P


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray1P_HH
