// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/random/mt19937.hh
/// @brief  Mersenne Twister 19937 random number generator
/// @author Ian W. Davis
///
/// @details Implementation of the Mersenne Twister random number generator
/// with "19937" parameters.  Copied and pasted into wrapper class from the
/// inventor's C code available from:
///
/// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
///
/// It's BSD licensed so we can incorporate and redistribute it freely.  I
/// choose the "dSFMT" code, which is oriented toward generating doubles
/// instead of ints and is highly optimized for modern processors, although I
/// omitted the specializations for SSE2 and Altivec that are present in the
/// original for ease of compilation.  (On my machine it's still 7% faster than
/// ran3, although random number generation is NOT the slow step in Rosetta!)
/// This is dSFMT version 1.2.1, the most current as of Nov 2007.
///
/// The Mersenne Twister is a fairly high quality, fairly fast, and widely used
/// pseudo-random number generator.  It has many parameter sets, which differ
/// in the amount of storage they require for their state, but 19937 is the
/// most used.
///
/// The major advantage over ran3 is that mt19937 can be seeded from any 32-bit
/// value, so up to ~4 billion unique simulation trajectories are possible.
/// (In fact, it can be seeded from an array of values instead, leading to even
/// more possibilities.) Although I can't find documentation, I understand from
/// David Baker that ran3 can only accept seeds in the range from approx. -1
/// million to -4 million, meaning that only ~3 million unique simulations are
/// possible.  This sometimes causes problems for jobs running on BOINC, where
/// jobs are distributed over many, many processors.  Also, ran3 has a period
/// length of (2**55)-1, whereas mt19937 has a period of (2**19937)-1, although
/// I kind of doubt we're ever really running up against that limitation.
///
/// Boost.Random has an implementation also, and says this about it:
///
/// cycle length: (2**19937) - 1
/// good uniform distribution in up to 623 dimensions
/// It is recommended as the default random number generator.
///
/// The GNU Scientific Library has an implementation also, and says this about it:
///
/// The MT19937 generator of Makoto Matsumoto and Takuji Nishimura is a variant
/// of the twisted generalized feedback shift-register algorithm, and is known
/// as the Mersenne Twister generator. It has a Mersenne prime period of
/// 2^19937 - 1 (about 10^6000) and is equi-distributed in 623 dimensions. It
/// has passed the DIEHARD statistical tests. It uses 624 words of state per
/// generator and is comparable in speed to the other generators.
///
/// See Makoto Matsumoto and Takuji Nishimura, Mersenne Twister: A
/// 623-dimensionally equidistributed uniform pseudorandom number
/// generator. ACM Transactions on Modeling and Computer Simulation, Vol.
/// 8, No. 1 (Jan. 1998), Pages 3-30


#ifndef INCLUDED_numeric_random_mt19937_hh
#define INCLUDED_numeric_random_mt19937_hh

#include <numeric/random/uniform.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <cstring>
#include <iostream>

/**
* @file dSFMT.h
*
* @brief double precision SIMD oriented Fast Mersenne Twister(dSFMT)
* pseudorandom number generator based on IEEE 754 format.
*
* @author Mutsuo Saito (Hiroshima University)
* @author Makoto Matsumoto (Hiroshima University)
*
* Copyright (C) 2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima
* University. All rights reserved.
*
* The new BSD License is applied to this software.
* see LICENSE.txt
*
* @note We assume that your system has inttypes.h.  If your system
* doesn't have inttypes.h, you have to typedef uint32_t and uint64_t,
* and you have to define PRIu64 and PRIx64 in this file as follows:
* @verbatim
typedef unsigned int uint32_t
typedef unsigned long long uint64_t
#define PRIu64 "llu"
#define PRIx64 "llx"
@endverbatim
* uint32_t must be exactly 32-bit unsigned integer type (no more, no
* less), and uint64_t must be exactly 64-bit unsigned integer type.
* PRIu64 and PRIx64 are used for printf function to print 64-bit
* unsigned int and 64-bit unsigned int in hexadecimal format.
*/


#include <stdio.h>

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
#include <inttypes.h>
#elif defined(_MSC_VER) || defined(__BORLANDC__)
typedef unsigned int uint32_t;
typedef unsigned long long uint64_t;
#define inline __inline
#else
	#include <inttypes.h>
#if defined(__GNUC__)
		#define inline __inline__
#else
		#define inline
#endif
#endif

#ifndef PRIu64
#if defined(_MSC_VER) || defined(__BORLANDC__)
#define PRIu64 "I64u"
#define PRIx64 "I64x"
#else
		#define PRIu64 "llu"
		#define PRIx64 "llx"
#endif
#endif

#ifndef UINT64_C
#define UINT64_C(v) (v ## ULL)
#endif


#define DSFMT_MEXP 19937
/*-----------------
BASIC DEFINITIONS
-----------------*/
/** Mersenne Exponent. The period of the sequence
*  is a multiple of 2^DSFMT_MEXP-1.
* #define DSFMT_MEXP 19937 */
/** DSFMT generator has an internal state array of 128-bit integers,
* and N is its size. */
#define DSFMT_N (DSFMT_MEXP / 104)
/** N32 is the size of internal state array when regarded as an array
* of 32-bit integers.*/
#define DSFMT_N32 (DSFMT_N * 4)
/** N64 is the size of internal state array when regarded as an array
* of 64-bit integers.*/
#define DSFMT_N64 (DSFMT_N * 2)

/*----------------------
the parameters of SFMT
following definitions are in dSFMT-paramsXXXX.h file.
----------------------*/
/** the pick up position of the array.
#define DSFMT_POS1 122
*/

/** the parameter of shift left as four 32-bit registers.
#define DSFMT_SL1 18
*/

/** the parameter of shift left as one 128-bit register.
* The 128-bit integer is shifted by (SL2 * 8) bits.
#define DSFMT_SL2 1
*/

/** the parameter of shift right as four 32-bit registers.
#define DSFMT_SR1 11
*/

/** the parameter of shift right as one 128-bit register.
* The 128-bit integer is shifted by (SL2 * 8) bits.
#define DSFMT_SR2 1
*/

/** A bitmask, used in the recursion.  These parameters are introduced
* to break symmetry of SIMD.
#define DSFMT_MSK1 (uint64_t)0xdfffffefULL
#define DSFMT_MSK2 (uint64_t)0xddfecb7fULL
*/

/** These definitions are part of a 128-bit period certification vector.
#define DSFMT_PCV1 UINT64_C(0x00000001)
#define DSFMT_PCV2 UINT64_C(0x00000000)
*/

#define DSFMT_LOW_MASK  UINT64_C(0x000FFFFFFFFFFFFF)
#define DSFMT_LOW_MASK32_1 0x000fffffU
#define DSFMT_LOW_MASK32_2 0xffffffffU
#define DSFMT_HIGH_CONST UINT64_C(0x3FF0000000000000)
#define DSFMT_HIGH_CONST32 0x3ff00000U


#define DSFMT_POS1 36
#define DSFMT_SL1 29
#define DSFMT_SL2 1
#define DSFMT_SR1 7
#define DSFMT_SR2 16
#define DSFMT_MSK1 UINT64_C(0x57fbfffdffff575f)
#define DSFMT_MSK2 UINT64_C(0xffff6febffffffee)
#define DSFMT_MSK32_1 0x57fbfffdU
#define DSFMT_MSK32_2 0xffff575fU
#define DSFMT_MSK32_3 0xffff6febU
#define DSFMT_MSK32_4 0xffffffeeU
#define DSFMT_PCV1 UINT64_C(0x0000000000000001)
#define DSFMT_PCV2 UINT64_C(0x000ec8f3d0b00000)
#define DSFMT_IDSTR \
	"dDSFMT-19937:36-29-1-7-16:57fbfffdffff575f-ffff6febffffffee"


/*------------------------------------------
128-bit SIMD like data type for standard C
------------------------------------------*/
/** 128-bit data structure */
union W128_T {
	uint64_t u[2];
	uint32_t u32[4];
	double d[2];
	};

/** 128-bit data type */
typedef union W128_T w128_t;

namespace numeric {
namespace random {


/*----------------
STATIC FUNCTIONS
----------------*/
inline static void lshift128(w128_t *out, const w128_t *in, int shift);
inline static uint32_t ini_func1(uint32_t x);
inline static uint32_t ini_func2(uint32_t x);
inline static int sformat_idxof(int i);


/**
* This function simulate a 32-bit array index overlapped to 64-bit
* array of LITTLE ENDIAN in BIG ENDIAN machine.
*/
//
// defined(__amd64) seems to be wrong here - commented out for now.
//
//#if (defined(__BIG_ENDIAN__) || defined(BIG_ENDIAN)) && !defined(__amd64)

#if (defined(__BIG_ENDIAN__) || defined(BIG_ENDIAN))
inline static int sformat_idxof(int i) {
	return i ^ 1;
}
#else
inline static int sformat_idxof(int i) {
		return i;
}
#endif

/**
* This function simulates SIMD 128-bit left shift by the standard C.
* The 128-bit integer given in \b in is shifted by (shift * 8) bits.
* This function simulates the LITTLE ENDIAN SIMD.
* @param out the output of this function
* @param in the 128-bit data to be shifted
* @param shift the shift value
*/
inline static void lshift128(w128_t *out, const w128_t *in, int shift) {
	out->u[0] = in->u[0] << (shift * 8);
	out->u[1] = in->u[1] << (shift * 8);
	out->u[1] |= in->u[0] >> (64 - shift * 8);
}

/**
* This function represents the recursion formula.
* @param r output
* @param a a 128-bit part of the internal state array
* @param b a 128-bit part of the internal state array
* @param c a 128-bit part of the internal state array
* @param lung a 128-bit part of the internal state array
*/
inline static void do_recursion(w128_t *r, w128_t *a, w128_t *b, w128_t *c,
	w128_t *lung) {
	w128_t x;

	lshift128(&x, a, DSFMT_SL2);
	r->u[0] = a->u[0] ^ x.u[0] ^ ((b->u[0] >> DSFMT_SR1) & DSFMT_MSK1)
		^ (c->u[0] >> DSFMT_SR2) ^ (c->u[0] << DSFMT_SL1) ^ lung->u[1];
	r->u[1] = a->u[1] ^ x.u[1] ^ ((b->u[1] >> DSFMT_SR1) & DSFMT_MSK2)
		^ (c->u[1] >> DSFMT_SR2) ^ (c->u[1] << DSFMT_SL1) ^ lung->u[0];
	r->u[0] &= DSFMT_LOW_MASK;
	r->u[1] &= DSFMT_LOW_MASK;
	lung->u[0] ^= r->u[0];
	lung->u[1] ^= r->u[1];
	r->u[0] |= DSFMT_HIGH_CONST;
	r->u[1] |= DSFMT_HIGH_CONST;
}

/**
* This function represents a function used in the initialization
* by init_by_array
* @param x 32-bit integer
* @return 32-bit integer
*/
static uint32_t ini_func1(uint32_t x) {
	return (x ^ (x >> 27)) * (uint32_t)1664525UL;
}

/**
* This function represents a function used in the initialization
* by init_by_array
* @param x 32-bit integer
* @return 32-bit integer
*/
static uint32_t ini_func2(uint32_t x) {
	return (x ^ (x >> 27)) * (uint32_t)1566083941UL;
}


/// @brief
///
/// @details
///
class mt19937_RG : public uniform_RG
{
public:

	mt19937_RG() :
		uniform_RG(),
		iseed_(0), // Clear, for consistent starting state
		sformat_idx(0), // Clear for consistent starting state
		is_sformat_initialized(0)
	{
		psformat64 = &sformat[0].d[0];
	}

	virtual ~mt19937_RG() {}

	/// @brief Set seed and state
	inline void setSeed(int const iseed)
	{
		iseed_ = iseed;
		uint32_t seed = (uint32_t) iseed;
		int i;
		uint32_t *psformat;

		psformat = &sformat[0].u32[0];
		psformat[sformat_idxof(0)] = seed;
		for ( i = 1; i < (DSFMT_N + 1) * 4; i++ ) {
			psformat[sformat_idxof(i)] = 1812433253UL
				* (psformat[sformat_idxof(i - 1)]
				^ (psformat[sformat_idxof(i - 1)] >> 30)) + i;
		}
		initial_mask();
		period_certification();
		sformat_idx = DSFMT_N64;
		is_sformat_initialized = 1;
	}

	/// @brief Set seed and state
	void setSeed( std::string const & ) { assert( false ); } //< Not implemented yet!

	inline int getSeed() {
		return iseed_;
	}

	inline double getRandom()
	{
		//This function generates and returns double precision pseudorandom
		//number which distributes uniformly in the range [1, 2).  This is
		//the primitive and faster than generating numbers in other ranges.
		//init_gen_rand() or init_by_array() must be called before this
		//function.
		//@return double precision floating point pseudorandom number

		double r;

		assert(is_sformat_initialized);

		if ( sformat_idx >= DSFMT_N * 2 ) {
			gen_rand_all();
			sformat_idx = 0;
		}
		r = psformat64[sformat_idx++];
		return r - 1.0; // to be on [0,1)
	}

	/// @brief Serializes generator state to stream losslessly.
	virtual void saveState(std::ostream & out)
	{
		out << " " << is_sformat_initialized;
		out << " " << sformat_idx;
		// psformat64 is a pointer to a memory location and is never modified
		// it should NOT be serialized or restored
		for ( int i = 0; i < DSFMT_N+1; ++i ) {
			out << " " << sformat[i].u[0] << " " << sformat[i].u[1];
		}
	}

	/// @brief Deserializes generator state from stream losslessly.
	virtual void restoreState(std::istream & in)
	{
		in >> is_sformat_initialized;
		in >> sformat_idx;
		// psformat64 is a pointer to a memory location and is never modified
		// it should NOT be serialized or restored
		for ( int i = 0; i < DSFMT_N+1; ++i ) {
			in >> sformat[i].u[0] >> sformat[i].u[1];
		}
	}

protected:
	/**
	* This function initializes the internal state array to fit the IEEE
	* 754 format.
	*/
	void initial_mask(void) {
		int i;
		uint64_t *psformat;

		psformat = &sformat[0].u[0];
		for ( i = 0; i < (DSFMT_N + 1) * 2; i++ ) {
			psformat[i] = (psformat[i] & DSFMT_LOW_MASK) | DSFMT_HIGH_CONST;
		}
	}

	/**
	* This function certificate the period of 2^{DSFMT_MEXP}-1.
	*/
	void period_certification() {
		int i, j;
		uint64_t pcv[2] = {DSFMT_PCV1, DSFMT_PCV2};
		uint64_t inner;
		uint64_t new_lung[2];
		//uint64_t work;
		uint64_t fix[2];

		fix[0] = (((DSFMT_HIGH_CONST >> DSFMT_SR1) & DSFMT_MSK2)
			^ (DSFMT_HIGH_CONST >> DSFMT_SR2)) | DSFMT_HIGH_CONST;
		fix[1] = (((DSFMT_HIGH_CONST >> DSFMT_SR1) & DSFMT_MSK1)
			^ (DSFMT_HIGH_CONST >> DSFMT_SR2)) | DSFMT_HIGH_CONST;
		fix[0] = fix[0] ^ (DSFMT_HIGH_CONST >> (64 - 8 * DSFMT_SL2));
		new_lung[0] = sformat[DSFMT_N].u[0] ^ fix[0];
		new_lung[1] = sformat[DSFMT_N].u[1] ^ fix[1];
		inner = new_lung[0] & pcv[0];
		inner ^= new_lung[1] & pcv[1];
		for ( i = 32; i > 0; i >>= 1 ) {
			inner ^= inner >> i;
		}
		inner &= 1;
		/* check OK */
		if ( inner == 1 ) {
			return;
		}
		/* check NG, and modification */
		for ( i = 0; i < 2; i++ ) {
			uint64_t work = 1;
			for ( j = 0; j < 52; j++ ) {
				if ( (work & pcv[i]) != 0 ) {
					sformat[DSFMT_N].u[i] ^= work;
					return;
				}
				work = work << 1;
			}
		}
	}

	/**
	* This function initializes the internal state array,
	* with an array of 32-bit integers used as the seeds
	* @param init_key the array of 32-bit integers, used as a seed.
	* @param key_length the length of init_key.
	*/
	void init_by_array(uint32_t init_key[], int key_length) {
		int i, j, count;
		uint32_t r;
		uint32_t *psformat32;
		int lag;
		int mid;
		int size = (DSFMT_N + 1) * 4; /* pulmonary */


		if ( size >= 623 ) {
			lag = 11;
		} else if ( size >= 68 ) {
			lag = 7;
		} else if ( size >= 39 ) {
			lag = 5;
		} else {
			lag = 3;
		}
		mid = (size - lag) / 2;

		psformat32 = &sformat[0].u32[0];
		memset(sformat, 0x8b, sizeof(sformat));
		if ( key_length + 1 > size ) {
			count = key_length + 1;
		} else {
			count = size;
		}
		r = ini_func1(psformat32[sformat_idxof(0)] ^ psformat32[sformat_idxof(mid % size)]
			^ psformat32[sformat_idxof((size - 1) % size)]);
		psformat32[sformat_idxof(mid % size)] += r;
		r += key_length;
		psformat32[sformat_idxof((mid + lag) % size)] += r;
		psformat32[sformat_idxof(0)] = r;
		i = 1;
		count--;
		for ( i = 1, j = 0; (j < count) && (j < key_length); j++ ) {
			r = ini_func1(psformat32[sformat_idxof(i)]
				^ psformat32[sformat_idxof((i + mid) % size)]
				^ psformat32[sformat_idxof((i + size - 1) % size)]);
			psformat32[sformat_idxof((i + mid) % size)] += r;
			r += init_key[j] + i;
			psformat32[sformat_idxof((i + mid + lag) % size)] += r;
			psformat32[sformat_idxof(i)] = r;
			i = (i + 1) % size;
		}
		for ( ; j < count; j++ ) {
			r = ini_func1(psformat32[sformat_idxof(i)]
				^ psformat32[sformat_idxof((i + mid) % size)]
				^ psformat32[sformat_idxof((i + size - 1) % size)]);
			psformat32[sformat_idxof((i + mid) % size)] += r;
			r += i;
			psformat32[sformat_idxof((i + mid + lag) % size)] += r;
			psformat32[sformat_idxof(i)] = r;
			i = (i + 1) % size;
		}
		for ( j = 0; j < size; j++ ) {
			r = ini_func2(psformat32[sformat_idxof(i)]
				+ psformat32[sformat_idxof((i + mid) % size)]
				+ psformat32[sformat_idxof((i + size - 1) % size)]);
			psformat32[sformat_idxof((i + mid) % size)] ^= r;
			r -= i;
			psformat32[sformat_idxof((i + mid + lag) % size)] ^= r;
			psformat32[sformat_idxof(i)] = r;
			i = (i + 1) % size;
		}
		initial_mask();
		period_certification();
		sformat_idx = DSFMT_N64;
		is_sformat_initialized = 1;
	}

	/**
	* This function fills the internal state array with double precision
	* floating point pseudorandom numbers of the IEEE 754 format.
	*/
	inline void gen_rand_all(void) {
		int i;
		w128_t lung;

		lung = sformat[DSFMT_N];
		do_recursion(&sformat[0], &sformat[0], &sformat[DSFMT_POS1], &sformat[DSFMT_N -1], &lung);
		for ( i = 1; i < DSFMT_N - DSFMT_POS1; i++ ) {
			do_recursion(&sformat[i], &sformat[i], &sformat[i + DSFMT_POS1], &sformat[i - 1],
				&lung);
		}
		for ( ; i < DSFMT_N; i++ ) {
			do_recursion(&sformat[i], &sformat[i], &sformat[i + DSFMT_POS1 - DSFMT_N],
				&sformat[i - 1], &lung);
		}
		sformat[DSFMT_N] = lung;
	}

private:
	/** The integer used to seed the RNG. */
	int iseed_;
	/** the 128-bit internal state array */
	w128_t sformat[DSFMT_N + 1];
	/** the double pointer to the 128-bit internal state array */
	double *psformat64;// = &sformat[0].d[0];
	/** index counter to the internal state array as double */
	int sformat_idx;
	/** a flag: it is 0 if and only if the internal state is not yet
	* initialized. */
	int is_sformat_initialized;// = 0;

}; // mt19937_RG


} // namespace random
} // namespace numeric

#endif // INCLUDED_numeric_random_mt19937_HH
