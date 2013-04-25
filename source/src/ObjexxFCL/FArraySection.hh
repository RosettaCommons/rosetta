#ifndef INCLUDED_ObjexxFCL_FArraySection_hh
#define INCLUDED_ObjexxFCL_FArraySection_hh


// FArraySection: Fortran-Compatible Array Section Proxy
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
#include <ObjexxFCL/FArraySection.fwd.hh>


namespace ObjexxFCL {


// Forward Declarations
template< typename > class FArray;


/// @brief FArraySection: Fortran-Compatible Array Section Proxy
template< typename T >
class FArraySection
{


private: // Friend


	friend class FArray< T >;


public: // Types


	// STL style
	typedef  T  value_type;
	typedef  T &  reference;
	typedef  T const &  const_reference;
	typedef  T *  pointer;
	typedef  T const *  const_pointer;
	typedef  std::size_t  size_type;
	typedef  std::ptrdiff_t  difference_type;

	// C++ style
	typedef  T  Value;
	typedef  T &  Reference;
	typedef  T const &  ConstReference;
	typedef  T *  Pointer;
	typedef  T const *  ConstPointer;
	typedef  std::size_t  Size;
	typedef  std::ptrdiff_t  Difference;


public: // Creation


	/// @brief Copy Constructor
	inline
	FArraySection( FArraySection const & s ) :
		array_( s.array_ ),
		size_( s.size_ )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( true )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{}


	/// @brief Copy from non-Constant Constructor
	inline
	FArraySection( FArraySection & s ) :
		array_( s.array_ ),
		size_( s.size_ )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( false )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{}


	/// @brief Const Pointer + Size Constructor
	inline
	FArraySection( T const * array_a, size_type const size_a ) :
		array_( const_cast< T * >( array_a ) ),
		size_( size_a )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( true )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{}


	/// @brief Pointer + Size Constructor
	inline
	FArraySection( T * array_a, size_type const size_a ) :
		array_( array_a ),
		size_( size_a )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( false )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{}


	/// @brief Size + Const Pointer Constructor
	inline
	FArraySection( size_type const size_a, T const * array_a ) :
		array_( const_cast< T * >( array_a ) ),
		size_( size_a )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( true )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{}


	/// @brief Size + Pointer Constructor
	inline
	FArraySection( size_type const size_a, T * array_a ) :
		array_( array_a ),
		size_( size_a )
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
		, const_proxy_( false )
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
	{}


	/// @brief Destructor
	inline
	~FArraySection()
	{}


public: // Assignment


	/// @brief Copy Assignment
	inline
	FArraySection &
	operator =( FArraySection const & s )
	{
		if ( this != &s ) {
			array_ = s.array_;
			size_ = s.size_;
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
			const_proxy_ = s.const_proxy_;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS
		}
		return *this;
	}


public: // Inspector


	/// @brief Size
	inline
	size_type
	size() const
	{
		return size_;
	}


private: // Data


	/// @brief Pointer (non-owning) to data array
	T * array_;

	/// @brief Size of data array
	size_type size_;

	/// @brief Proxy for const data array?
#ifdef OBJEXXFCL_PROXY_CONST_CHECKS
	bool const_proxy_;
#endif // OBJEXXFCL_PROXY_CONST_CHECKS


}; // FArraySection


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_FArraySection_HH
