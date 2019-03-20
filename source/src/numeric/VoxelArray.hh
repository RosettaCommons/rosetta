// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/VoxelArray.hh
/// @brief  N-dimensional array class to store values at float offsets
/// @author William Sheffler -- original implementation
/// @author Brian Coventry (bcov@uw.edu) -- port to Rosetta


#ifndef INCLUDED_numeric_VoxelArray_HH
#define INCLUDED_numeric_VoxelArray_HH

#include <numeric/VoxelArray.fwd.hh>

#include <utility/exit.hh>

#include <boost/type_traits.hpp>
#include <boost/assert.hpp>

#include <vector>
#include <array>

#ifdef CEREAL
#include <cereal/access.hpp>
#endif
//#include <boost/serialization/access.hpp>
//#include <boost/serialization/split_member.hpp>

namespace numeric {

// DIM used to be a template parameter, but we can't link against
//  boost::multi_array which is where the magic happens.
//  If the version of boost that rosetta links against is ever
//  updated, this go can back to being N-dimensional.
template< class _Float/*=float*/, class _Value/*=float*/ >
class VoxelArray {
public:
	typedef VoxelArray<_Float,_Value> THIS;
	static size_t const DIM = 3;
	typedef _Float Float;
	typedef _Value Value;
	typedef std::array< Size, DIM > Indices;
	typedef std::array< Float, DIM > Bounds;

	typedef typename std::vector<Value>::reference data_reference_type;
	typedef typename std::vector<Value>::const_reference const_data_reference_type;
	// typedef util::SimpleArray<DIM,typename BASE::size_type> Indices;
	// typedef util::SimpleArray<DIM,Float> Bounds;

private:
	// All the things that boost::multi_array used to do...

	void
	resize( Indices const & indices ) {

		data_.clear();
		Size the_size = 1;
		for ( Size i = 0; i < DIM; i++ ) {
			the_size *= indices[i];
		}
		data_.resize( the_size );

		shape_ = indices;
	}



	data_reference_type
	operator() ( Indices const & indices ) {
		Size index = 0;
		for ( Size i = 0; i < DIM; ++i ) {
			index *= shape_[i];
			index += indices[i];
		}
		return data_[ index ];
	}

	const_data_reference_type
	operator() ( Indices const & indices ) const {
		Size index = 0;
		for ( Size i = 0; i < DIM; ++i ) {
			index *= shape_[i];
			index += indices[i];
		}
		return data_[ index ];
	}


public:
	std::vector<Value> &
	data( ) {
		return data_;
	}

	std::vector<Value> const &
	data( ) const {
		return data_;
	}

	Size
	num_elements( ) const {
		return data_.size();
	}

	Indices const &
	shape( ) const {
		return shape_;
	}

public:


	VoxelArray() {}

	template<class F1,class F2, class F3>
	VoxelArray(F1 const & lb, F2 const & ub, F3 const & cs ) {

		for ( Size i = 0; i < DIM; i++ ) {
			lb_[ i ] = lb[ i ];
			ub_[ i ] = ub[ i ];
			cs_[ i ] = cs[ i ];
		}

		Indices extents = floats_to_index(ub_);
		// std::cout << extents << std::endl;
		for ( Size i = 0; i < DIM; i++ ) extents[i] += 1; // pad by one
		this->resize(extents);
	}

	template<class Floats> Indices floats_to_index(Floats const & f) const {
		Indices ind;
		for ( Size i = 0; i < DIM; ++i ) {
			Float tmp = ((f[i]-lb_[i])/cs_[i]);
			// assert(tmp >= 0.0);
			ind[i] = tmp;
		}
		// std::cout << "floats_to_index " << f << " " << ind << std::endl;
		return ind;
	}

	Bounds indices_to_center( Indices const & idx ) const {
		Bounds c;
		for ( Size i = 0; i < DIM; ++i ) {
			c[i] = (idx[i]+0.5)*cs_[i] + lb_[i];
		}
		return c;
	}

	template<class Floats>
	const_data_reference_type
	operator[](Floats const & floats) const { return this->operator()(floats_to_index(floats)); }

	template<class Floats>
	data_reference_type
	operator[](Floats const & floats){ return this->operator()(floats_to_index(floats)); }

	Value at( Float f, Float g, Float h ) const {
		Indices idx = floats_to_index( Bounds{ { f, g, h } } );
		if ( idx[0] < this->shape()[0] && idx[1] < this->shape()[1] && idx[2] < this->shape()[2] ) {
			return this->operator()(idx);
		} else return Value(0);
	}

	template<class V>
	Value at( V const & v ) const {
		Indices idx = floats_to_index( Bounds{ { v[0], v[1], v[2] } } );
		if ( idx[0] < this->shape()[0] && idx[1] < this->shape()[1] && idx[2] < this->shape()[2] ) {
			return this->operator()(idx);
		} else return Value(0);
	}

	// void write(std::ostream & out) const {
	//  out.write( (char const*)&lb_, sizeof(Bounds) );
	//  out.write( (char const*)&ub_, sizeof(Bounds) );
	//  out.write( (char const*)&cs_, sizeof(Bounds) );
	//  for(size_t i = 0; i < DIM; ++i) out.write( (char const*)&(this->shape()[i]), sizeof() );
	//  out.write( (char const*)this->data(), this->num_elements()*sizeof(Float) );
	// }
	// void read(std::istream & in){
	//  in.read( (char*)&lb_, sizeof(Bounds) );
	//  in.read( (char*)&ub_, sizeof(Bounds) );
	//  in.read( (char*)&cs_, sizeof(Bounds) );
	//  in.read( (char*)this->data(), this->num_elements()*sizeof(Float) );
	// }
	bool operator==(THIS const & o) const {
		return lb_==o.lb_ && ub_==o.ub_ && cs_==o.cs_ && shape_ == o.shape_ && data_ == o.data_;
	}

#ifdef CEREAL
	// friend class boost::serialization::access;
    friend class cereal::access; // befriend the cereal version of access
#endif

	template<class Archive> void save(Archive & ar, const unsigned int ) const {
		BOOST_VERIFY( boost::is_pod<Float>::type::value );
		ar & lb_;
		ar & ub_;
		ar & cs_;
		for ( size_t i = 0; i < DIM; ++i ) ar & this->shape()[i];
		for ( size_t i = 0; i < this->num_elements(); ++i ) ar & this->data()[i];
	}
	template<class Archive> void load(Archive & ar, const unsigned int ){
		BOOST_VERIFY( boost::is_pod<Float>::type::value );
		ar & lb_;
		ar & ub_;
		ar & cs_;
		Indices extents;
		for ( size_t i = 0; i < DIM; ++i ) ar & extents[i];
		this->resize(extents);
		for ( size_t i = 0; i < this->num_elements(); ++i ) ar & this->data()[i];
	}
	void save( std::ostream & out ) const {
		BOOST_VERIFY( boost::is_pod<Float>::type::value );
		out.write( (char*)&lb_, sizeof(Bounds) );
		out.write( (char*)&ub_, sizeof(Bounds) );
		out.write( (char*)&cs_, sizeof(Bounds) );
		for ( size_t i = 0; i < DIM; ++i ) {
			out.write( (char*)(&(this->shape()[i])), sizeof(Size) );
		}
		for ( size_t i = 0; i < this->num_elements(); ++i ) out.write( (char*)(&(this->data()[i])), sizeof(Float) );
	}
	void load( std::istream & in ){
		BOOST_VERIFY( boost::is_pod<Float>::type::value );
		runtime_assert( in.good() );
		in.read( (char*)&lb_, sizeof(Bounds) );
		runtime_assert( in.good() );
		in.read( (char*)&ub_, sizeof(Bounds) );
		runtime_assert( in.good() );
		in.read( (char*)&cs_, sizeof(Bounds) );
		runtime_assert( in.good() );
		// todo: should I check these against the c'tor values? if not, should add default ctor?
		Indices extents;
		for ( size_t i = 0; i < DIM; ++i ) {
			in.read( (char*)(&(extents[i])), sizeof(Size) );
		}
		runtime_assert( in.good() );
		this->resize(extents);
		for ( size_t i = 0; i < this->num_elements(); ++i ) in.read( (char*)(&(this->data()[i])), sizeof(Float) );
		runtime_assert( in.good() );
	}
	// BOOST_SERIALIZATION_SPLIT_MEMBER()

private:
	Bounds lb_,ub_,cs_;

	// Multi_array stuff
	std::vector<Value> data_;
	Indices shape_;

};

template< class F, class V >
std::ostream & operator << ( std::ostream & out, VoxelArray<F,V> const & v ){
	out << "VoxelArray( lb: " << v.lb_ << " ub: " << v.ub_ << " cs: " << v.cs_ << " nelem: " << v.num_elements() << " sizeof_val: " << sizeof(V) << " )";
	return out;
}

}

#endif
