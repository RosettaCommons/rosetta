// -*- c++ -*-
/*
 * Copyright (c) 2010-2012, Jim Bosch
 * All rights reserved.
 *
 * ndarray is distributed under a simple BSD-like license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/ndarray/ndarray
 */
#ifndef NDARRAY_DETAIL_ViewBuilder_h_INCLUDED
#define NDARRAY_DETAIL_ViewBuilder_h_INCLUDED

/** 
 *  \file ndarray/detail/ViewBuilder.h @brief Implementation of arbitrary views into arrays.
 */

#include "ndarray/views.h"
#include <boost/fusion/include/push_back.hpp>
#include <boost/fusion/include/pop_back.hpp>
#include <boost/fusion/include/front.hpp>
#include <boost/fusion/include/back.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/reverse_view.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/tuple.hpp>
#include <boost/mpl/count.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/fold.hpp>

namespace ndarray {
namespace detail {

/** 
 *  @internal
 *  @brief A temporary object used in constructing a Core object in a view operation.
 */
template <typename T, int M, int N>
struct CoreTransformer {
    T * _data;
    typename Core<M>::ConstPtr _input;
    typename Core<N>::Ptr _output;
    
    CoreTransformer(
        T * data,
        typename Core<M>::ConstPtr const & input,
        typename Core<N>::Ptr const & output
    ) : _data(data), _input(input), _output(output) {}
    
    template <int M1, int N1>
    CoreTransformer(CoreTransformer<T,M1,N1> const & other) : 
        _data(other._data), _input(other._input), _output(other._output) {}
};

template <int N, int C, int I>
struct Dimensions {
    typedef boost::mpl::int_<N> ND;  // Number of dimensions in output array
    typedef boost::mpl::int_<C> RMC; // Number of contiguous dimensions in output array, from end.
    typedef boost::mpl::int_<I> IDX; // Current dimension of input array being processed.
    typedef boost::mpl::int_<N-I> N_I; 
};

template <typename Index> struct IndexTraits;

template <>
struct IndexTraits<ndarray::index::Slice> {

    template <typename D> 
    struct Append {
        typedef Dimensions< 
            D::ND::value, 
            ((D::RMC::value < D::N_I::value) ? D::RMC::value : (D::N_I::value - 1)),
            (D::IDX::value + 1)
        > Type;
    };

    /// @brief Metafunction for the result type of transform().
    template <typename T, int M, int N> struct TransformCoreResult {
        typedef CoreTransformer<T,M-1,N-1> Type;
    };

    /// @brief Apply a slice index.
    template <typename T, int M, int N>
    static CoreTransformer<T,M-1,N-1> transformCore(
        ndarray::index::Slice const & index, CoreTransformer<T,M,N> & t
    ) {
        NDARRAY_ASSERT(index.step > 0);
        NDARRAY_ASSERT(index.start <= index.stop);
        NDARRAY_ASSERT(index.start >= 0);
        NDARRAY_ASSERT(index.stop <= t._input->getSize());
        t._data += index.start * t._input->getStride();
        t._output->setSize(index.computeSize());
        t._output->setStride(t._input->getStride() * index.step);
        return t;
    }

};

template <>
struct IndexTraits<ndarray::index::Range> {

    template <typename D> 
    struct Append {
        typedef Dimensions< 
            D::ND::value, 
            ((D::RMC::value < D::N_I::value) ? D::RMC::value : D::N_I::value),
            (D::IDX::value + 1)
        > Type;
    };

    /// @brief Metafunction for the result type of transform().
    template <typename T, int M, int N> struct TransformCoreResult {
        typedef CoreTransformer<T,M-1,N-1> Type;
    };

    /// @brief Apply a range index.
    template <typename T, int M, int N>
    static CoreTransformer<T,M-1,N-1> transformCore(
        ndarray::index::Range const & index, CoreTransformer<T,M,N> & t
    ) {
        NDARRAY_ASSERT(index.start <= index.stop);
        NDARRAY_ASSERT(index.start >= 0);
        NDARRAY_ASSERT(index.stop <= t._input->getSize());
        t._data += index.start * t._input->getStride();
        t._output->setSize(index.stop - index.start);
        t._output->setStride(t._input->getStride());
        return t;
    }
};

template <>
struct IndexTraits<ndarray::index::Full> {

    template <typename D> 
    struct Append {
        typedef Dimensions< 
            D::ND::value, 
            D::RMC::value,
            (D::IDX::value + 1)
        > Type;
    };

    /// @brief Metafunction for the result type of transform().
    template <typename T, int M, int N> struct TransformCoreResult {
        typedef CoreTransformer<T,M-1,N-1> Type;
    };

    /// @brief Apply a full dimension index.
    template <typename T, int M, int N>
    static CoreTransformer<T,M-1,N-1> transformCore(
        ndarray::index::Full const &, CoreTransformer<T,M,N> & t
    ) {
        t._output->setSize(t._input->getSize());
        t._output->setStride(t._input->getStride());
        return t;
    }
};

template <>
struct IndexTraits<ndarray::index::Scalar> {

    template <typename D> 
    struct Append {
        typedef Dimensions< 
            (D::ND::value - 1),
            ((D::RMC::value < (D::N_I::value - 1)) ? D::RMC::value : (D::N_I::value - 1)),
            D::IDX::value
        > Type;
    };

    /// @brief Metafunction for the result type of transform().
    template <typename T, int M, int N> struct TransformCoreResult {
        typedef CoreTransformer<T,M-1,N> Type;
    };

    /// @brief Apply a scalar dimension index.
    template <typename T, int M, int N>
    static CoreTransformer<T,M-1,N> transformCore(
        ndarray::index::Scalar const & index, CoreTransformer<T,M,N> & t
    ) {
        NDARRAY_ASSERT(index.n >= 0);
        NDARRAY_ASSERT(index.n < t._input->getSize());
        t._data += index.n * t._input->getStride();
        return t;
    }
};

template <typename T, int M, int N, typename Index>
typename IndexTraits<Index>::template TransformCoreResult<T,M,N>::Type
transformCore(Index const & index, CoreTransformer<T,M,N> & t) {
    return IndexTraits<Index>::transformCore(index, t);
}

struct AppendIndex {

    template <typename State, typename Index>
    struct apply {
        typedef typename IndexTraits<Index>::template Append<State>::Type type;
    };

};

template <int N, int C, typename Seq_, bool isColumnMajor = (C < 0)>
struct ViewTraits;

template <int N, int C, typename Seq_>
struct ViewTraits<N,C,Seq_,false> {

    typedef typename boost::mpl::fold< Seq_, Dimensions<N,C,0>, AppendIndex >::type Dims;
    
    typedef typename Dims::ND ND;
    typedef typename Dims::RMC RMC;

};

template <int N, int C, typename Seq_>
struct ViewTraits<N,C,Seq_,true> {

    typedef typename boost::mpl::fold< 
        boost::fusion::reverse_view< typename boost::fusion::result_of::as_vector<Seq_>::type >, 
        Dimensions<N,-C,0>, AppendIndex
    >::type Dims;
    
    typedef typename Dims::ND ND;
    typedef typename boost::mpl::negate<typename Dims::RMC>::type RMC;

};

/**
 *  @internal
 *  @brief Metafunction that pads a View with extra FullDim indices to make size<Seq_>::type::value == N.
 */
template <int N, typename Seq_, bool IsNormalized=(boost::mpl::template size<Seq_>::type::value==N)>
struct ViewNormalizer {

    typedef typename boost::fusion::result_of::push_back<Seq_ const,index::Full>::type Next;

    typedef typename ViewNormalizer<N,Next>::Output Output;

    static Output apply(Seq_ const & input) {
        return ViewNormalizer<N,Next>::apply(
            boost::fusion::push_back(input, index::Full())
        );
    }
};

template <int N, typename Seq_>
struct ViewNormalizer<N,Seq_,true> {
    typedef typename boost::fusion::result_of::as_vector<Seq_>::type Output;
    static Output apply(Seq_ const & input) { return boost::fusion::as_vector(input); }
};

/**
 *  @internal @ingroup InternalGroup
 *  @brief Static function object that constructs a view Array.
 */
template <typename Array_, typename InSeq>
struct ViewBuilder {
    typedef ExpressionTraits<Array_> Traits;
    typedef typename Traits::Element Element;
    typedef typename Traits::ND InputND;
    typedef typename Traits::RMC InputRMC;
    typedef typename Traits::Core InputCore;
    typedef boost::mpl::bool_<(InputRMC::value < 0)> IsColumnMajor;

    typedef ViewNormalizer<InputND::value,InSeq> Normalizer;
    typedef typename ViewNormalizer<InputND::value,InSeq>::Output NormSeq;
    typedef ViewTraits<InputND::value,InputRMC::value,NormSeq> OutputTraits;

    typedef typename OutputTraits::ND OutputND;
    typedef typename OutputTraits::RMC OutputRMC;

    typedef ArrayRef<Element,OutputND::value,OutputRMC::value> OutputArray;
    typedef Core<OutputND::value> OutputCore;

    static OutputArray apply(Array_ const & array, InSeq const & seq) {
        CoreTransformer<Element,InputND::value,OutputND::value> initial(
            array.getData(), 
            ArrayAccess< Array_>::getCore(array),
            OutputCore::create(array.getManager())
        );
        NormSeq normSeq = Normalizer::apply(seq);
        std::pair<Element*,typename OutputCore::Ptr> final = process(normSeq, initial);
        return ArrayAccess< OutputArray >::construct(final.first, final.second);
    }

    template <int M, int N>
    static std::pair<Element*,typename OutputCore::Ptr>
    process(NormSeq const & seq, CoreTransformer<Element,M,N> t) {
        return process(seq, transformCore(boost::fusion::at_c<(InputND::value-M)>(seq), t));
    }

    static std::pair<Element*,typename OutputCore::Ptr>
    process(NormSeq const & seq, CoreTransformer<Element,0,0> t) {
        return std::make_pair(t._data, boost::static_pointer_cast<OutputCore>(t._output));
    }
    
};

/**
 *  @internal @ingroup InternalGroup
 *  @brief Wrapper function for ViewBuilder that removes the need to specify its template parameters.
 */
template <typename Array_, typename Seq_>
typename ViewBuilder<Array_, Seq_>::OutputArray
buildView(Array_ const & array, Seq_ const & seq) {
    return ViewBuilder<Array_,Seq_>::apply(array, seq);
};

} // namespace detail

} // namespace ndarray

#endif // !NDARRAY_DETAIL_ViewBuilder_h_INCLUDED
