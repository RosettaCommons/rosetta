#ifndef INCLUDED_ObjexxFCL_ubyte_hh
#define INCLUDED_ObjexxFCL_ubyte_hh


// ubyte: Unsigned One-Byte Integer
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
#include <ObjexxFCL/ubyte.fwd.hh>

// C++ Headers
#include <cassert>
#include <cstddef>
#include <istream>


namespace ObjexxFCL {


/// @brief ubyte: One-Byte Integer
class ubyte
{


public: // Creation


	/// @brief Default Constructor
	inline
	ubyte() :
		b_( 0 )
	{}


	/// @brief Constructor
	inline
	explicit
	ubyte( unsigned short int const i ) :
		b_( i )
	{}


	/// @brief Destructor
	inline
	~ubyte()
	{}


public: // Conversion


	/// @brief short Conversion
	inline
	operator unsigned short int() const
	{
		return static_cast< unsigned short int >( b_ );
	}


public: // Assignment


	/// @brief = short
	inline
	ubyte &
	operator =( unsigned short int const i )
	{
		b_ = i;
		return *this;
	}


	/// @brief += short
	inline
	ubyte &
	operator +=( unsigned short int const i )
	{
		b_ += i;
		return *this;
	}


	/// @brief -= short
	inline
	ubyte &
	operator -=( unsigned short int const i )
	{
		b_ -= i;
		return *this;
	}


	/// @brief *= short
	inline
	ubyte &
	operator *=( unsigned short int const i )
	{
		b_ *= i;
		return *this;
	}


	/// @brief /= short
	inline
	ubyte &
	operator /=( unsigned short int const i )
	{
		assert( i != 0 );
		b_ /= i;
		return *this;
	}


public: // Incrememt/Decrement


	/// @brief ++ubyte
	inline
	ubyte &
	operator ++()
	{
		++b_;
		return *this;
	}


	/// @brief ubyte++
	inline
	ubyte const
	operator ++( int )
	{
		ubyte const old( *this );
		++b_;
		return old;
	}


	/// @brief --ubyte
	inline
	ubyte &
	operator --()
	{
		--b_;
		return *this;
	}


	/// @brief ubyte--
	inline
	ubyte const
	operator --( int )
	{
		ubyte const old( *this );
		--b_;
		return old;
	}


public: // Math


	/// @brief +ubyte
	inline
	ubyte
	operator +() const
	{
		return *this;
	}


	/// @brief -ubyte
	inline
	ubyte
	operator -() const
	{
		return ubyte( -static_cast< unsigned short int >( b_ ) );
	}


	/// @brief ubyte + ubyte
	friend
	inline
	ubyte
	operator +( ubyte const & i, ubyte const & j )
	{
		return ubyte( i.b_ + j.b_ );
	}


	/// @brief ubyte - ubyte
	friend
	inline
	ubyte
	operator -( ubyte const & i, ubyte const & j )
	{
		return ubyte( i.b_ - j.b_ );
	}


	/// @brief ubyte * ubyte
	friend
	inline
	ubyte
	operator *( ubyte const & i, ubyte const & j )
	{
		return ubyte( i.b_ * j.b_ );
	}


	/// @brief ubyte / ubyte
	friend
	inline
	ubyte
	operator /( ubyte const & i, ubyte const & j )
	{
		assert( j.b_ != 0 );
		return ubyte( i.b_ / j.b_ );
	}


public: // Bitwise Logical


	/// @brief ~ubyte
	inline
	ubyte
	operator ~() const
	{
		return ubyte( ~b_ );
	}


	/// @brief ubyte >> std::size_t
	inline
	ubyte
	operator >>( std::size_t const n ) const
	{
		return ubyte( b_ >> n );
	}


	/// @brief ubyte >> ubyte
	inline
	ubyte
	operator >>( ubyte const & n ) const
	{
		return ubyte( b_ >> static_cast< unsigned short int >( n ) );
	}


	/// @brief ubyte << std::size_t
	inline
	ubyte
	operator <<( std::size_t const n ) const
	{
		return ubyte( b_ << n );
	}


	/// @brief ubyte << ubyte
	inline
	ubyte
	operator <<( ubyte const & n ) const
	{
		return ubyte( b_ << static_cast< unsigned short int >( n ) );
	}


	/// @brief &= ubyte
	inline
	ubyte &
	operator &=( ubyte const & i )
	{
		b_ &= i.b_;
		return *this;
	}


	/// @brief |= ubyte
	inline
	ubyte &
	operator |=( ubyte const & i )
	{
		b_ |= i.b_;
		return *this;
	}


	/// @brief ^= ubyte
	inline
	ubyte &
	operator ^=( ubyte const & i )
	{
		b_ ^= i.b_;
		return *this;
	}


	/// @brief ubyte & ubyte
	friend
	inline
	ubyte
	operator &( ubyte const & i, ubyte const & j )
	{
		return ubyte( i.b_ & j.b_ );
	}


	/// @brief ubyte | ubyte
	friend
	inline
	ubyte
	operator |( ubyte const & i, ubyte const & j )
	{
		return ubyte( i.b_ | j.b_ );
	}


	/// @brief ubyte ^ ubyte
	friend
	inline
	ubyte
	operator ^( ubyte const & i, ubyte const & j )
	{
		return ubyte( i.b_ ^ j.b_ );
	}


public: // Comparison


	/// @brief ubyte == ubyte
	friend
	inline
	bool
	operator ==( ubyte const & i, ubyte const & j )
	{
		return ( i.b_ == j.b_ );
	}


	/// @brief ubyte != ubyte
	friend
	inline
	bool
	operator !=( ubyte const & i, ubyte const & j )
	{
		return ( i.b_ != j.b_ );
	}


	/// @brief ubyte < ubyte
	friend
	inline
	bool
	operator <( ubyte const & i, ubyte const & j )
	{
		return ( i.b_ < j.b_ );
	}


	/// @brief ubyte <= ubyte
	friend
	inline
	bool
	operator <=( ubyte const & i, ubyte const & j )
	{
		return ( i.b_ <= j.b_ );
	}


	/// @brief ubyte > ubyte
	friend
	inline
	bool
	operator >( ubyte const & i, ubyte const & j )
	{
		return ( i.b_ > j.b_ );
	}


	/// @brief ubyte >= ubyte
	friend
	inline
	bool
	operator >=( ubyte const & i, ubyte const & j )
	{
		return ( i.b_ >= j.b_ );
	}


public: // I/O


	/// @brief Stream Input
	friend
	inline
	std::istream &
	operator >>( std::istream & stream, ubyte & b )
	{
		if ( stream ) {
			unsigned short int s;
			stream >> s;
			b.b_ = s;
		}
		return stream;
	}


private: // Data


	/// @brief Value
	unsigned char b_;


}; // ubyte


/// @brief ubyte + ubyte
#ifndef __clang__
ubyte
operator +( ubyte const & i, ubyte const & j );
#endif


/// @brief ubyte - ubyte
ubyte
operator -( ubyte const & i, ubyte const & j );


/// @brief ubyte * ubyte
ubyte
operator *( ubyte const & i, ubyte const & j );


/// @brief ubyte / ubyte
ubyte
operator /( ubyte const & i, ubyte const & j );


/// @brief ubyte & ubyte
#ifndef __clang__
ubyte
operator &( ubyte const & i, ubyte const & j );
#endif


/// @brief ubyte | ubyte
#ifndef __clang__
ubyte
operator |( ubyte const & i, ubyte const & j );
#endif


/// @brief ubyte ^ ubyte
ubyte
operator ^( ubyte const & i, ubyte const & j );


/// @brief ubyte == ubyte
bool
operator ==( ubyte const & i, ubyte const & j );


/// @brief ubyte != ubyte
bool
operator !=( ubyte const & i, ubyte const & j );


/// @brief ubyte < ubyte
bool
operator <( ubyte const & i, ubyte const & j );


/// @brief ubyte <= ubyte
bool
operator <=( ubyte const & i, ubyte const & j );


/// @brief ubyte > ubyte
bool
operator >( ubyte const & i, ubyte const & j );


/// @brief ubyte >= ubyte
bool
operator >=( ubyte const & i, ubyte const & j );


/// @brief Stream Input
#ifndef __clang__
std::istream &
operator >>( std::istream & stream, ubyte & b );
#endif


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_ubyte_HH
