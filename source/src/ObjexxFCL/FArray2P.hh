#ifndef INCLUDED_ObjexxFCL_FArray2P_hh
#define INCLUDED_ObjexxFCL_FArray2P_hh


// FArray2P: Fortran-Compatible 2D Proxy Array
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
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>


namespace ObjexxFCL {


/// @brief FArray2P: Fortran-Compatible 2D Proxy Array
template< typename T >
class FArray2P :
	public FArray2< T >,
	public ObserverMulti
{


private: // Types


	typedef  FArray2< T >  Super;
	typedef  typename Super::real_FArray  real_FArray;
	typedef  typename Super::proxy_FArray  proxy_FArray;
	typedef  internal::ProxySentinel  ProxySentinel;


private: // Friend


	friend class FArray2A< T >;


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
	using Super::s1_;


public: // Creation


	/// @brief Default Constructor
	inline
	FArray2P() :
		Super( ProxySentinel() ),
		source_( 0 )
	{}


	/// @brief Copy Constructor
	inline
	FArray2P( FArray2P const & a ) :
		Super( a, ProxySentinel() ),
		ObserverMulti(),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		source_( &a )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		insert_as_observer();
	}


	/// @brief Non-Const Copy Constructor
	inline
	FArray2P( FArray2P & a ) :
		Super( a, ProxySentinel() ),
		ObserverMulti(),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		source_( &a )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		insert_as_observer();
	}


	/// @brief Real Constructor
	inline
	FArray2P( real_FArray const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		source_( &a )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		insert_as_observer();
	}


	/// @brief Non-Const Real Constructor
	inline
	FArray2P( real_FArray & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1_ ),
		I2_( a.I2_ ),
		source_( &a )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		insert_as_observer();
	}


	/// @brief Super Constructor
	inline
	FArray2P( Super const & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1() ),
		I2_( a.I2() ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		insert_as_observer();
	}


	/// @brief Non-Const Super Constructor
	inline
	FArray2P( Super & a ) :
		Super( a, ProxySentinel() ),
		I1_( a.I1() ),
		I2_( a.I2() ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		shift_set( a.shift_ );
		s1_ = a.s1_;
		insert_as_observer();
	}


	/// @brief Base Constructor
	inline
	FArray2P( Base const & a ) :
		Super( a, ProxySentinel() ),
		I1_( 1 ),
		I2_( a.size() ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		shift_set( 2 );
		s1_ = 1;
		insert_as_observer();
	}


	/// @brief Non-Const Base Constructor
	inline
	FArray2P( Base & a ) :
		Super( a, ProxySentinel() ),
		I1_( 1 ),
		I2_( a.size() ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		shift_set( 2 );
		s1_ = 1;
		insert_as_observer();
	}


	/// @brief Section Constructor
	inline
	FArray2P( Section const & s ) :
		Super( s, ProxySentinel() ),
		I1_( 1 ),
		I2_( s.size() ),
		source_( 0 )
	{
		shift_set( 2 );
		s1_ = 1;
		insert_as_observer();
	}


	/// @brief Non-Const Section Constructor
	inline
	FArray2P( Section & s ) :
		Super( s, ProxySentinel() ),
		I1_( 1 ),
		I2_( s.size() ),
		source_( 0 )
	{
		shift_set( 2 );
		s1_ = 1;
		insert_as_observer();
	}


	/// @brief Value Constructor
	inline
	FArray2P( T const & t ) :
		Super( t, ProxySentinel() ),
		I1_( 1 ),
		I2_( star ), // Unbounded
		source_( 0 )
	{
		shift_set( 2 );
		s1_ = 1;
		insert_as_observer();
	}


	/// @brief Non-Const Value Constructor
	inline
	FArray2P( T & t ) :
		Super( t, ProxySentinel() ),
		I1_( 1 ),
		I2_( star ), // Unbounded
		source_( 0 )
	{
		shift_set( 2 );
		s1_ = 1;
		insert_as_observer();
	}


	/// @brief Copy + IndexRange Constructor
	inline
	FArray2P( FArray2P const & a, IR const & I1_a, IR const & I2_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		source_( &a )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Non-Const Copy + IndexRange Constructor
	inline
	FArray2P( FArray2P & a, IR const & I1_a, IR const & I2_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		source_( &a )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Real + IndexRange Constructor
	inline
	FArray2P( real_FArray const & a, IR const & I1_a, IR const & I2_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		source_( &a )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Non-Const Real + IndexRange Constructor
	inline
	FArray2P( real_FArray & a, IR const & I1_a, IR const & I2_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		source_( &a )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Super + IndexRange Constructor
	inline
	FArray2P( Super const & a, IR const & I1_a, IR const & I2_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Non-Const Super + IndexRange Constructor
	inline
	FArray2P( Super & a, IR const & I1_a, IR const & I2_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Base + IndexRange Constructor
	inline
	FArray2P( Base const & a, IR const & I1_a, IR const & I2_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Non-Const Base + IndexRange Constructor
	inline
	FArray2P( Base & a, IR const & I1_a, IR const & I2_a ) :
		Super( a, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		source_( dynamic_cast< SubjectMulti const * >( &a ) )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Section + IndexRange Constructor
	inline
	FArray2P( Section const & s, IR const & I1_a, IR const & I2_a ) :
		Super( s, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		source_( 0 )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Non-Const Section + IndexRange Constructor
	inline
	FArray2P( Section & s, IR const & I1_a, IR const & I2_a ) :
		Super( s, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		source_( 0 )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Value + IndexRange Constructor
	inline
	FArray2P( T const & t, IR const & I1_a, IR const & I2_a ) :
		Super( t, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		source_( 0 )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Non-Const Value + IndexRange Constructor
	inline
	FArray2P( T & t, IR const & I1_a, IR const & I2_a ) :
		Super( t, ProxySentinel() ),
		I1_( I1_a ),
		I2_( I2_a ),
		source_( 0 )
	{
		dimension_proxy();
		insert_as_observer();
	}


	/// @brief Destructor
	inline
	virtual
	~FArray2P()
	{
		if ( source_ ) source_->remove_observer( *this );
	}


public: // Assignment


	/// @brief Copy Assignment
	inline
	FArray2P &
	operator =( FArray2P const & a )
	{
		if ( this != &a ) {
			if ( ! equal_dimension( a ) ) dimension( a );
			Base::operator =( a );
		}
		return *this;
	}


	/// @brief Super Assignment
	inline
	FArray2P &
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
	FArray2P &
	operator =( FArray2< U > const & a )
	{
		if ( ! equal_dimension( a ) ) dimension( a );
		Base::operator =( a );
		return *this;
	}


	/// @brief += Array Template
	template< typename U >
	inline
	FArray2P &
	operator +=( FArray2< U > const & a )
	{
		Super::operator +=( a );
		return *this;
	}


	/// @brief -= Array Template
	template< typename U >
	inline
	FArray2P &
	operator -=( FArray2< U > const & a )
	{
		Super::operator -=( a );
		return *this;
	}


	/// @brief *= Array Template
	template< typename U >
	inline
	FArray2P &
	operator *=( FArray2< U > const & a )
	{
		Super::operator *=( a );
		return *this;
	}


	/// @brief = Value
	inline
	FArray2P &
	operator =( T const & t )
	{
		Super::operator =( t );
		return *this;
	}


	/// @brief += Value
	inline
	FArray2P &
	operator +=( T const & t )
	{
		Super::operator +=( t );
		return *this;
	}


	/// @brief -= Value
	inline
	FArray2P &
	operator -=( T const & t )
	{
		Super::operator -=( t );
		return *this;
	}


	/// @brief *= Value
	inline
	FArray2P &
	operator *=( T const & t )
	{
		Super::operator *=( t );
		return *this;
	}


	/// @brief /= Value
	inline
	FArray2P &
	operator /=( T const & t )
	{
		Super::operator /=( t );
		return *this;
	}


public: // Subscript


	/// @brief array( i1, i2 ) const
	inline
	T const &
	operator ()( int const i1, int const i2 ) const
	{
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) );
		return sarray_[ ( i2 * s1_ ) + i1 ];
	}


	/// @brief array( i1, i2 )
	inline
	T &
	operator ()( int const i1, int const i2 )
	{
		proxy_const_assert( not_const_proxy() );
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) );
		return sarray_[ ( i2 * s1_ ) + i1 ];
	}


	/// @brief Section Starting at array( i1, i2 )
	inline
	Section const
	a( int const i1, int const i2 ) const
	{
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) );
		size_type const offset( ( ( i2 * s1_ ) + i1 ) - shift_ );
		return Section( static_cast< T const * >( array_ + offset ), ( array_size_ != npos ) ? array_size_ - offset : npos );
	}


	/// @brief Non-Const Section Starting at array( i1, i2 )
	inline
	Section
	a( int const i1, int const i2 )
	{
		proxy_const_assert( not_const_proxy() );
		assert( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) );
		size_type const offset( ( ( i2 * s1_ ) + i1 ) - shift_ );
		return Section( array_ + offset, ( array_size_ != npos ) ? array_size_ - offset : npos );
	}


	/// @brief Linear Index
	inline
	size_type
	index( int const i1, int const i2 ) const
	{
		assert( ( I1_.initialized() ) && ( I2_.initialized() ) );
		return ( ( ( i2 * s1_ ) + i1 ) - shift_ );
	}


public: // Predicate


	/// @brief Dimensions Initialized?
	inline
	bool
	dimensions_initialized() const
	{
		return ( ( I1_.initialized() ) && ( I2_.initialized() ) );
	}


	/// @brief Contains Indexed Element?
	inline
	bool
	contains( int const i1, int const i2 ) const
	{
		return ( ( I1_.contains( i1 ) ) && ( I2_.contains( i2 ) ) );
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


	/// @brief Size of Dimension 2
	inline
	size_type
	size2() const
	{
		return I2_.size();
	}


public: // Modifier


	/// @brief Clear
	inline
	FArray2P &
	clear()
	{
		Super::clear();
		I1_.clear_no_notify();
		I2_.clear_no_notify();
		source_ = 0;
		return *this;
	}


	/// @brief Dimension by IndexRanges Even if Const
	inline
	FArray2P const &
	dim( IR const & I1_a, IR const & I2_a ) const
	{
		const_cast< FArray2P & >( *this ).dimension( I1_a, I2_a );
		return *this;
	}


	/// @brief Dimension by Array Even if Const
	template< typename U >
	inline
	FArray2P const &
	dim( FArray2< U > const & a ) const
	{
		const_cast< FArray2P & >( *this ).dimension( a );
		return *this;
	}


	/// @brief Dimension by IndexRanges
	inline
	FArray2P &
	dimension( IR const & I1_a, IR const & I2_a )
	{
		I1_.assign_no_notify( I1_a );
		I2_.assign_no_notify( I2_a );
		dimension_proxy();
		return *this;
	}


	/// @brief Dimension by Array
	template< typename U >
	inline
	FArray2P &
	dimension( FArray2< U > const & a )
	{
		I1_.assign_no_notify( a.I1() );
		I2_.assign_no_notify( a.I2() );
		dimension_proxy();
		return *this;
	}


	/// @brief Attach to Proxy Array
	inline
	FArray2P &
	attach( FArray2P const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		I1_ = a.I1_;
		I2_ = a.I2_;
		if ( source_ ) source_->remove_observer( *this );
		source_ = &a;
		a.insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Non-Const Proxy Array
	inline
	FArray2P &
	attach( FArray2P & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		I1_ = a.I1_;
		I2_ = a.I2_;
		if ( source_ ) source_->remove_observer( *this );
		source_ = &a;
		a.insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Real Array
	inline
	FArray2P &
	attach( real_FArray const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		I1_ = a.I1_;
		I2_ = a.I2_;
		if ( source_ ) source_->remove_observer( *this );
		source_ = &a;
		a.insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Non-Const Real Array
	inline
	FArray2P &
	attach( real_FArray & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		I1_ = a.I1_;
		I2_ = a.I2_;
		if ( source_ ) source_->remove_observer( *this );
		source_ = &a;
		a.insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Super Array
	inline
	FArray2P &
	attach( Super const & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		I1_ = a.I1();
		I2_ = a.I2();
		if ( source_ ) source_->remove_observer( *this );
		source_ = dynamic_cast< SubjectMulti const * >( &a );
		if ( source_ ) source_->insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Non-Const Super Array
	inline
	FArray2P &
	attach( Super & a )
	{
		Base::attach( a );
		s1_ = a.s1_;
		I1_ = a.I1();
		I2_ = a.I2();
		if ( source_ ) source_->remove_observer( *this );
		source_ = dynamic_cast< SubjectMulti const * >( &a );
		if ( source_ ) source_->insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Base Array
	inline
	FArray2P &
	attach( Base const & a )
	{
		Base::attach( a, 2 );
		s1_ = 1;
		I1_ = 1;
		I2_ = a.size();
		if ( source_ ) source_->remove_observer( *this );
		source_ = dynamic_cast< SubjectMulti const * >( &a );
		if ( source_ ) source_->insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Non-Const Base Array
	inline
	FArray2P &
	attach( Base & a )
	{
		Base::attach( a, 2 );
		s1_ = 1;
		I1_ = 1;
		I2_ = a.size();
		if ( source_ ) source_->remove_observer( *this );
		source_ = dynamic_cast< SubjectMulti const * >( &a );
		if ( source_ ) source_->insert_observer( *this );
		return *this;
	}


	/// @brief Attach to Section
	inline
	FArray2P &
	attach( Section const & s )
	{
		Base::attach( s, 2 );
		s1_ = 1;
		I1_ = 1;
		I2_ = s.size();
		if ( source_ ) source_->remove_observer( *this );
		source_ = 0;
		return *this;
	}


	/// @brief Attach to Non-Const Section
	inline
	FArray2P &
	attach( Section & s )
	{
		Base::attach( s, 2 );
		s1_ = 1;
		I1_ = 1;
		I2_ = s.size();
		if ( source_ ) source_->remove_observer( *this );
		source_ = 0;
		return *this;
	}


	/// @brief Attach to Value
	inline
	FArray2P &
	attach( T const & t )
	{
		Base::attach( t, 2 );
		s1_ = 1;
		I1_ = 1;
		I2_ = star; // Unbounded
		if ( source_ ) source_->remove_observer( *this );
		source_ = 0;
		return *this;
	}


	/// @brief Attach to Non-Const Value
	inline
	FArray2P &
	attach( T & t )
	{
		Base::attach( t, 2 );
		s1_ = 1;
		I1_ = 1;
		I2_ = star; // Unbounded
		if ( source_ ) source_->remove_observer( *this );
		source_ = 0;
		return *this;
	}


	/// @brief Detach from Source Array
	inline
	FArray2P &
	detach()
	{
		Base::detach();
		s1_ = 0;
		I1_.clear();
		I2_.clear();
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
			s1_ = 0;
			I1_.clear();
			I2_.clear();
			source_ = 0;
		}
	}


protected: // Functions


	/// @brief Dimension by IndexRanges
	inline
	void
	dimension_assign( SIR const & I1_a, SIR const & I2_a )
	{
		I1_.assign_no_notify( I1_a );
		I2_.assign_no_notify( I2_a );
		dimension_proxy();
	}


private: // Functions


	/// @brief Dimension by Current IndexRanges
	inline
	void
	dimension_proxy()
	{
		assert( I1_.not_unbounded() );
		s1_ = I1_.size();
		if ( dimensions_initialized() ) {
			if ( I2_.bounded_value() ) { // Bounded
				size_set( this->size_of( s1_, I2_.size() ) );
			} else if ( array_size_ != npos ) { // Unbounded with bounded data array
				if ( s1_ > 0 ) { // Infer upper index and size
					I2_.u_no_notify( I2_.lz() + ( array_size_ / s1_ ) - 1 );
					size_set( this->size_of( s1_, I2_.size() ) );
				} else {
					size_set( array_size_ );
				}
			} else { // Unbounded with unbounded data array
				size_set( npos );
			}
			shift_set( ( I2_.lz() * s1_ ) + I1_.lz() );
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
		I1_.insert_observer( *this );
		I2_.insert_observer( *this );
		if ( source_ ) source_->insert_observer( *this );
	}

	/* // Unused private function
	/// @brief Remove as Observer of the IndexRanges and Source Array
	inline
	void
	remove_as_observer()
	{
		I1_.remove_observer( *this );
		I2_.remove_observer( *this );
		if ( source_ ) source_->remove_observer( *this );
	}
	*/

private: // Data


	/// @brief Dimension 1 index range
	IR I1_;

	/// @brief Dimension 2 index range
	IR I2_;

	/// @brief Pointer (non-owning) to source array (0 if unknown)
	SubjectMulti const * source_;


}; // FArray2P


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArray2P_HH
