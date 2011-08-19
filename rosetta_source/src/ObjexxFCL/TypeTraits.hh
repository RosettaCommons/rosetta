#ifndef INCLUDED_ObjexxFCL_TypeTraits_hh
#define INCLUDED_ObjexxFCL_TypeTraits_hh


// TypeTraits: Type Traits Template
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


namespace ObjexxFCL {


/// @brief TypeTraits: Type Traits Template
template< typename T >
struct TypeTraits
{
	typedef  T  traits_type;


	/// @brief Initial Value
	inline
	static
	traits_type
	initial_value()
	{
		return traits_type(); // Use default constructor
	}


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0; // No precision for generic types
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 0; // No minimum width for generic types
	}


	/// @brief Standard Width
	inline
	static
	int
	standard_width()
	{
		return 0; // No minimum width for generic types
	}


}; // TypeTraits


/// @brief TypeTraits: Type Traits short int Specialization
template<>
struct TypeTraits< short int >
{
	typedef  short int  traits_type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 7;
	}


	/// @brief Standard Width
	inline
	static
	int
	standard_width()
	{
		return 7;
	}


}; // TypeTraits


/// @brief TypeTraits: Type Traits int Specialization
template<>
struct TypeTraits< int >
{
	typedef  int  traits_type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 12;
	}


	/// @brief Standard Width
	inline
	static
	int
	standard_width()
	{
		return 12;
	}


}; // TypeTraits


/// @brief TypeTraits: Type Traits long int Specialization
template<>
struct TypeTraits< long int >
{
	typedef  long int  traits_type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 22; // Big enough for 64-bit LP64 representation
	}


	/// @brief Standard Width
	inline
	static
	int
	standard_width()
	{
		return 22;
	}


}; // TypeTraits


/// @brief TypeTraits: Type Traits unsigned short int Specialization
template<>
struct TypeTraits< unsigned short int >
{
	typedef  unsigned short int  traits_type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 7;
	}


	/// @brief Standard Width
	inline
	static
	int
	standard_width()
	{
		return 7;
	}


}; // TypeTraits


/// @brief TypeTraits: Type Traits unsigned int Specialization
template<>
struct TypeTraits< unsigned int >
{
	typedef  unsigned int  traits_type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 12;
	}


	/// @brief Standard Width
	inline
	static
	int
	standard_width()
	{
		return 12;
	}


}; // TypeTraits


/// @brief TypeTraits: Type Traits unsigned long int Specialization
template<>
struct TypeTraits< unsigned long int >
{
	typedef  unsigned long int  traits_type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 0;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 22; // Big enough for 64-bit LP64 representation
	}


	/// @brief Standard Width
	inline
	static
	int
	standard_width()
	{
		return 22;
	}


}; // TypeTraits


/// @brief TypeTraits: Type Traits float Specialization
template<>
struct TypeTraits< float >
{
	typedef  float  traits_type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 8;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 15;
	}


	/// @brief Standard Width
	inline
	static
	int
	standard_width()
	{
		return 16;
	}


}; // TypeTraits


/// @brief TypeTraits: Type Traits double Specialization
template<>
struct TypeTraits< double >
{
	typedef  double  traits_type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 16;
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 23;
	}


	/// @brief Standard Width
	inline
	static
	int
	standard_width()
	{
		return 26;
	}


}; // TypeTraits


/// @brief TypeTraits: Type Traits long double Specialization
template<>
struct TypeTraits< long double >
{
	typedef  long double  traits_type;


	/// @brief Precision
	inline
	static
	int
	precision()
	{
		return 33; // Big enough for 128-bit representation
	}


	/// @brief Width
	inline
	static
	int
	width()
	{
		return 40; // Big enough for 128-bit representation
	}


	/// @brief Standard Width
	inline
	static
	int
	standard_width()
	{
		return 40;
	}


}; // TypeTraits


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_TypeTraits_HH
