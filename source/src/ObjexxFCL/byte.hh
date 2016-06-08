#ifndef INCLUDED_ObjexxFCL_byte_hh
#define INCLUDED_ObjexxFCL_byte_hh


// byte: One-Byte Signed Integer
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
#include <ObjexxFCL/byte.fwd.hh>

// C++ Headers
#include <cassert>
#include <cstddef>
#include <istream>


namespace ObjexxFCL {


/// @brief byte: One-Byte Integer
class byte
{


public: // Creation


	/// @brief Default Constructor
	inline
	byte() :
		b_( 0 )
	{}


	/// @brief Constructor
	inline
	explicit
	byte( short int const i ) :
		b_( i )
	{}


	/// @brief Destructor
	inline
	~byte()
	{}


public: // Conversion


	/// @brief short Conversion
	inline
	operator short int() const
	{
		return static_cast< short int >( b_ );
	}


public: // Assignment


	/// @brief = short
	inline
	byte &
	operator =( short int const i )
	{
		b_ = i;
		return *this;
	}


	/// @brief += short
	inline
	byte &
	operator +=( short int const i )
	{
		b_ += i;
		return *this;
	}


	/// @brief -= short
	inline
	byte &
	operator -=( short int const i )
	{
		b_ -= i;
		return *this;
	}


	/// @brief *= short
	inline
	byte &
	operator *=( short int const i )
	{
		b_ *= i;
		return *this;
	}


	/// @brief /= short
	inline
	byte &
	operator /=( short int const i )
	{
		assert( i != 0 );
		b_ /= i;
		return *this;
	}


public: // Incrememt/Decrement


	/// @brief ++byte
	inline
	byte &
	operator ++()
	{
		++b_;
		return *this;
	}


	/// @brief byte++
	inline
	byte const
	operator ++( int )
	{
		byte const old( *this );
		++b_;
		return old;
	}


	/// @brief --byte
	inline
	byte &
	operator --()
	{
		--b_;
		return *this;
	}


	/// @brief byte--
	inline
	byte const
	operator --( int )
	{
		byte const old( *this );
		--b_;
		return old;
	}


public: // Math


	/// @brief +byte
	inline
	byte
	operator +() const
	{
		return *this;
	}


	/// @brief -byte
	inline
	byte
	operator -() const
	{
		return byte( -static_cast< short int >( b_ ) );
	}


	/// @brief byte + byte
	friend
	inline
	byte
	operator +( byte const & i, byte const & j )
	{
		return byte( i.b_ + j.b_ );
	}


	/// @brief byte - byte
	friend
	inline
	byte
	operator -( byte const & i, byte const & j )
	{
		return byte( i.b_ - j.b_ );
	}


	/// @brief byte * byte
	friend
	inline
	byte
	operator *( byte const & i, byte const & j )
	{
		return byte( i.b_ * j.b_ );
	}


	/// @brief byte / byte
	friend
	inline
	byte
	operator /( byte const & i, byte const & j )
	{
		assert( j.b_ != 0 );
		return byte( i.b_ / j.b_ );
	}


public: // Bitwise Logical


	/// @brief ~byte
	inline
	byte
	operator ~() const
	{
		return byte( ~b_ );
	}


	/// @brief byte >> std::size_t
	inline
	byte
	operator >>( std::size_t const n ) const
	{
		return byte( b_ >> n );
	}


	/// @brief byte >> byte
	inline
	byte
	operator >>( byte const & n ) const
	{
		return byte( b_ >> static_cast< short int >( n ) );
	}


	/// @brief byte << std::size_t
	inline
	byte
	operator <<( std::size_t const n ) const
	{
		return byte( b_ << n );
	}


	/// @brief byte << byte
	inline
	byte
	operator <<( byte const & n ) const
	{
		return byte( b_ << static_cast< short int >( n ) );
	}


	/// @brief &= byte
	inline
	byte &
	operator &=( byte const & i )
	{
		b_ &= i.b_;
		return *this;
	}


	/// @brief |= byte
	inline
	byte &
	operator |=( byte const & i )
	{
		b_ |= i.b_;
		return *this;
	}


	/// @brief ^= byte
	inline
	byte &
	operator ^=( byte const & i )
	{
		b_ ^= i.b_;
		return *this;
	}


	/// @brief byte & byte
	friend
	inline
	byte
	operator &( byte const & i, byte const & j )
	{
		return byte( i.b_ & j.b_ );
	}


	/// @brief byte | byte
	friend
	inline
	byte
	operator |( byte const & i, byte const & j )
	{
		return byte( i.b_ | j.b_ );
	}


	/// @brief byte ^ byte
	friend
	inline
	byte
	operator ^( byte const & i, byte const & j )
	{
		return byte( i.b_ ^ j.b_ );
	}


public: // Comparison


	/// @brief byte == byte
	friend
	inline
	bool
	operator ==( byte const & i, byte const & j )
	{
		return ( i.b_ == j.b_ );
	}


	/// @brief byte != byte
	friend
	inline
	bool
	operator !=( byte const & i, byte const & j )
	{
		return ( i.b_ != j.b_ );
	}


	/// @brief byte < byte
	friend
	inline
	bool
	operator <( byte const & i, byte const & j )
	{
		return ( i.b_ < j.b_ );
	}


	/// @brief byte <= byte
	friend
	inline
	bool
	operator <=( byte const & i, byte const & j )
	{
		return ( i.b_ <= j.b_ );
	}


	/// @brief byte > byte
	friend
	inline
	bool
	operator >( byte const & i, byte const & j )
	{
		return ( i.b_ > j.b_ );
	}


	/// @brief byte >= byte
	friend
	inline
	bool
	operator >=( byte const & i, byte const & j )
	{
		return ( i.b_ >= j.b_ );
	}


public: // I/O


	/// @brief Stream Input
	friend
	inline
	std::istream &
	operator >>( std::istream & stream, byte & b )
	{
		if ( stream ) {
			short int s;
			stream >> s;
			b.b_ = s;
		}
		return stream;
	}


private: // Data


	/// @brief Value
	signed char b_;


}; // byte


/// @brief byte + byte
byte
operator +( byte const & i, byte const & j );


/// @brief byte - byte
byte
operator -( byte const & i, byte const & j );


/// @brief byte * byte
byte
operator *( byte const & i, byte const & j );


/// @brief byte / byte
byte
operator /( byte const & i, byte const & j );


/// @brief byte & byte
byte
operator &( byte const & i, byte const & j );


/// @brief byte | byte
byte
operator |( byte const & i, byte const & j );


/// @brief byte ^ byte
byte
operator ^( byte const & i, byte const & j );


/// @brief byte == byte
bool
operator ==( byte const & i, byte const & j );


/// @brief byte != byte
bool
operator !=( byte const & i, byte const & j );


/// @brief byte < byte
bool
operator <( byte const & i, byte const & j );


/// @brief byte <= byte
bool
operator <=( byte const & i, byte const & j );


/// @brief byte > byte
bool
operator >( byte const & i, byte const & j );


/// @brief byte >= byte
bool
operator >=( byte const & i, byte const & j );


/// @brief Stream Input
#ifndef __clang__
std::istream &
operator >>( std::istream & stream, byte & b );
#endif


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_byte_HH
