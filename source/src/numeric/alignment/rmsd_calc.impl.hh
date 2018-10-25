// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#pragma once

#include <cstdlib>

#include <Eigen/Core>

#include "ndarray.h"
#include "ndarray/eigen.h"

#include <numeric/alignment/QCPKernel.hh>

namespace numeric
{

namespace alignment
{

/* Calculate coordinate array rmsd between two coordinate input arrays.
* @param  first_coordinates  Array of shape [n|1, c, 3].
* @param  second_coordinates Array of shape [n|1, c, 3].
* @param  out                Output array of shape [n].
*/
template <class Real>
void coordinate_array_rmsd(
	ndarray::Array<Real, 3, 1> const & first_coordinates,
	ndarray::Array<Real, 3, 1> const & second_coordinates,
	ndarray::Array<Real, 1> out
)
{
	if ( first_coordinates.template getSize<2>() != 3 ) {
		throw std::invalid_argument("first_coordinates shape[2] != 3.");
	}

	if ( second_coordinates.template getSize<2>() != 3 ) {
		throw std::invalid_argument("second_coordinates shape[2] != 3.");
	}

	if ( first_coordinates.template getSize<1>() != second_coordinates.template getSize<1>() ) {
		throw std::invalid_argument("first_coordinates and seconds_coordinates shape[1] not equal.");
	}

	int first_size = first_coordinates.template getSize<0>();
	int second_size = second_coordinates.template getSize<0>();
	int out_size = out.template getSize<0>();
	int target_size = first_size != 1 ? first_size : second_size;

	if ( first_size != 1 && first_size != target_size ) {
		throw std::invalid_argument("Input coordinate array shape mismatch.");
	}

	if ( second_size != 1 && second_size != target_size ) {
		throw std::invalid_argument("Input coordinate array shape mismatch.");
	}

	if ( out_size != target_size ) {
		throw std::invalid_argument("out buffer of incorrect shape.");
	}

	typedef Eigen::Stride<Eigen::Dynamic, 1> CoordinateStride;
	typedef Eigen::Map< Eigen::Matrix<Real, 3, Eigen::Dynamic>, Eigen::Unaligned, CoordinateStride > CoordinateMap;

	size_t system_size = first_coordinates.template getSize<1>();
	size_t first_system_stride = first_coordinates.template getStride<0>();
	size_t second_system_stride = second_coordinates.template getStride<0>();
	CoordinateStride first_coordinate_stride( first_coordinates.template getStride<1>(), 1);
	CoordinateStride second_coordinate_stride( second_coordinates.template getStride<1>(), 1);

	Eigen::Matrix<Real, 3, Eigen::Dynamic> first_coordinate_centers(3, first_coordinates.template getSize<0>());
	for ( unsigned i = 0; i < first_coordinates.template getSize<0>(); i++ ) {
		first_coordinate_centers.col(i) = CoordinateMap(
			first_coordinates.getData() + i * first_system_stride,
			3, system_size, first_coordinate_stride).rowwise().sum() / system_size;
	}

	Eigen::Matrix<Real, 3, Eigen::Dynamic> second_coordinate_centers(3, second_coordinates.template getSize<0>());
	for ( unsigned i = 0; i < second_coordinates.template getSize<0>(); i++ ) {
		second_coordinate_centers.col(i) = CoordinateMap(
			second_coordinates.getData() + i * second_system_stride,
			3, system_size, second_coordinate_stride).rowwise().sum() / system_size;
	}

	for ( int n = 0; n < target_size; n++ ) {
		int i = first_size != 1 ? n : 0;
		int j = second_size != 1 ? n : 0;

		CoordinateMap first_c(
			first_coordinates.getData() + i * first_system_stride,
			3, system_size, first_coordinate_stride);

		CoordinateMap second_c(
			second_coordinates.getData() + j * second_system_stride,
			3, system_size, second_coordinate_stride);

		out(n) = QCPKernel<Real>::calc_coordinate_rmsd(
			first_c,
			first_coordinate_centers.col(i),
			second_c,
			second_coordinate_centers.col(j));
	}
}

/* Align superimpose to onto and update superimpose_coordinates with resulting superposition transform.
* @param  src_coordinates         Array of shape [n, c, 3].
* @param  onto_coordinates        Array of shape [n|1, c, 3].
* @param  superimpose_coordinates Array of shape [n, c2, 3].
* @param  out                     Output array of shape [n].
*/
template <class Real>
void coordinate_array_superimpose(
	ndarray::Array<Real, 3, 1> const & src_coordinates,
	ndarray::Array<Real, 3, 1> const & onto_coordinates,
	ndarray::Array<Real, 3, 1> const & superimpose_coordinates,
	ndarray::Array<Real, 1> out
)
{
	if ( src_coordinates.template getSize<2>() != 3 ) {
		throw std::invalid_argument("src_coordinates shape[2] != 3.");
	}

	if ( onto_coordinates.template getSize<2>() != 3 ) {
		throw std::invalid_argument("onto_coordinates shape[2] != 3.");
	}

	if ( superimpose_coordinates.template getSize<2>() != 3 ) {
		throw std::invalid_argument("superimpose_coordinates shape[2] != 3.");
	}

	if ( src_coordinates.template getSize<1>() != onto_coordinates.template getSize<1>() ) {
		throw std::invalid_argument("src_coordinates and ontos_coordinates shape[1] not equal.");
	}

	int src_size = src_coordinates.template getSize<0>();
	int onto_size = onto_coordinates.template getSize<0>();
	int superimpose_size = superimpose_coordinates.template getSize<0>();
	int out_size = out.template getSize<0>();
	int target_size = src_size;

	if ( onto_size != 1 && onto_size != target_size ) {
		throw std::invalid_argument("src/onto coordinate array shape mismatch.");
	}

	if ( superimpose_size != target_size ) {
		throw std::invalid_argument("src/superimpose coordinate array shape mismatch.");
	}

	if ( out_size != target_size ) {
		throw std::invalid_argument("src/out buffer buffer shape mismatch.");
	}

	typedef Eigen::Stride<Eigen::Dynamic, 1> CoordinateStride;
	typedef Eigen::Map< Eigen::Matrix<Real, 3, Eigen::Dynamic>, Eigen::Unaligned, CoordinateStride > CoordinateMap;

	size_t system_size = src_coordinates.template getSize<1>();
	size_t superimpose_system_size = superimpose_coordinates.template getSize<1>();

	size_t src_system_stride = src_coordinates.template getStride<0>();
	CoordinateStride src_coordinate_stride( src_coordinates.template getStride<1>(), 1);

	size_t onto_system_stride = onto_coordinates.template getStride<0>();
	CoordinateStride onto_coordinate_stride( onto_coordinates.template getStride<1>(), 1);

	size_t superimpose_system_stride = superimpose_coordinates.template getStride<0>();
	CoordinateStride superimpose_coordinate_stride( superimpose_coordinates.template getStride<1>(), 1);

	Eigen::Matrix<Real, 3, Eigen::Dynamic> src_coordinate_centers(3, src_coordinates.template getSize<0>());
	for ( unsigned i = 0; i < src_coordinates.template getSize<0>(); i++ ) {
		src_coordinate_centers.col(i) = CoordinateMap(
			src_coordinates.getData() + i * src_system_stride,
			3, system_size, src_coordinate_stride).rowwise().sum() / system_size;
	}

	Eigen::Matrix<Real, 3, Eigen::Dynamic> onto_coordinate_centers(3, onto_coordinates.template getSize<0>());
	for ( unsigned i = 0; i < onto_coordinates.template getSize<0>(); i++ ) {
		onto_coordinate_centers.col(i) = CoordinateMap(
			onto_coordinates.getData() + i * onto_system_stride,
			3, system_size, onto_coordinate_stride).rowwise().sum() / system_size;
	}

	for ( int src_i = 0; src_i < target_size; src_i++ ) {
		int onto_i = onto_size != 1 ? src_i : 0;

		CoordinateMap src_c(
			src_coordinates.getData() + src_i * src_system_stride,
			3, system_size, src_coordinate_stride);

		CoordinateMap onto_c(
			onto_coordinates.getData() + onto_i * onto_system_stride,
			3, system_size, onto_coordinate_stride);

		CoordinateMap superimpose_c(
			superimpose_coordinates.getData() + src_i * superimpose_system_stride,
			3, superimpose_system_size, superimpose_coordinate_stride);

		Eigen::Transform<Real, 3, Eigen::Affine> superposition_transform;

		out(src_i) = QCPKernel<Real>::calc_coordinate_superposition(
			src_c,
			src_coordinate_centers.col(src_i),
			onto_c,
			onto_coordinate_centers.col(onto_i),
			superposition_transform
		);

		superimpose_c = superposition_transform * superimpose_c;
	}
}

/*
* Calculate broadcast coordinate rmsd between two coordinate input arrays.
* @param  first_coordinates  Array of shape [a, c, 3].
* @param  second_coordinates Array of shape [b, c, 3].
* @param  out                Output array of shape [a, b].
*/
template <class Real>
void coordinate_array_broadcast_rmsd(
	ndarray::Array<Real, 3, 1> const & first_coordinates,
	ndarray::Array<Real, 3, 1> const & second_coordinates,
	ndarray::Array<Real, 2> out
)
{
	if ( first_coordinates.template getSize<2>() != 3 ) {
		throw std::invalid_argument("first_coordinates shape[2] != 3.");
	}

	if ( second_coordinates.template getSize<2>() != 3 ) {
		throw std::invalid_argument("second_coordinates shape[2] != 3.");
	}

	if ( first_coordinates.template getSize<1>() != second_coordinates.template getSize<1>() ) {
		throw std::invalid_argument("first_coordinates and seconds_coordinates shape[1] not equal.");
	}

	if (
			(first_coordinates.template getSize<0>()  != out.template getSize<0>()) |
			(second_coordinates.template getSize<0>() != out.template getSize<1>())
			) {
		throw std::invalid_argument("out buffer of incorrect shape.");
	}

	typedef Eigen::Stride<Eigen::Dynamic, 1> CoordinateStride;
	typedef Eigen::Map< Eigen::Matrix<Real, 3, Eigen::Dynamic>, Eigen::Unaligned, CoordinateStride > CoordinateMap;

	size_t system_size = first_coordinates.template getSize<1>();
	size_t first_system_stride = first_coordinates.template getStride<0>();
	size_t second_system_stride = second_coordinates.template getStride<0>();
	CoordinateStride first_coordinate_stride( first_coordinates.template getStride<1>(), 1);
	CoordinateStride second_coordinate_stride( second_coordinates.template getStride<1>(), 1);

	Eigen::Matrix<Real, 3, Eigen::Dynamic> first_coordinate_centers(3, first_coordinates.template getSize<0>());
	for ( unsigned i = 0; i < first_coordinates.template getSize<0>(); i++ ) {
		first_coordinate_centers.col(i) = CoordinateMap(
			first_coordinates.getData() + i * first_system_stride,
			3, system_size, first_coordinate_stride).rowwise().sum() / system_size;
	}

	Eigen::Matrix<Real, 3, Eigen::Dynamic> second_coordinate_centers(3, second_coordinates.template getSize<0>());
	for ( unsigned i = 0; i < second_coordinates.template getSize<0>(); i++ ) {
		second_coordinate_centers.col(i) = CoordinateMap(
			second_coordinates.getData() + i * second_system_stride,
			3, system_size, second_coordinate_stride).rowwise().sum() / system_size;
	}

	for ( unsigned i = 0; i < first_coordinates.template getSize<0>(); i++ ) {
		for ( unsigned j = 0; j < second_coordinates.template getSize<0>(); j++ ) {
			CoordinateMap first_c(
				first_coordinates.getData() + i * first_system_stride,
				3, system_size, first_coordinate_stride);

			CoordinateMap second_c(
				second_coordinates.getData() + j * second_system_stride,
				3, system_size, second_coordinate_stride);

			out(i, j) = QCPKernel<Real>::calc_coordinate_rmsd(
				first_c,
				first_coordinate_centers.col(i),
				second_c,
				second_coordinate_centers.col(j));
		}
	}
}

/*
* Calculate broadcast coordinate rmsd between two coordinate input arrays, with entries given by length and start index.
* @param   coordinates_per_entry       Number of coordinates per entry to align.
* @param  first_coordinates           Array of shape [_, 3].
* @param  first_coordinate_indicies   Array of shape [a].
* @param  second_coordinates           Array of shape [_, 3].
* @param  second_coordinate_indicies   Array of shape [b].
* @param  out                Output array of shape [a, b].
*/
template <typename Index, typename Real>
void indexed_coordinate_array_broadcast_rmsd(
	Index coordinates_per_entry,
	ndarray::Array<Real, 2, 1> const & first_coordinates,
	ndarray::Array<Index, 1> const & first_coordinate_indicies,
	ndarray::Array<Real, 2, 1> const & second_coordinates,
	ndarray::Array<Index, 1> const & second_coordinate_indicies,
	ndarray::Array<Real, 2> out
)
{
	if ( first_coordinates.template getSize<1>() != 3 ) {
		throw std::invalid_argument("first_coordinates shape[1] != 3.");
	}

	if ( second_coordinates.template getSize<1>() != 3 ) {
		throw std::invalid_argument("second_coordinates shape[1] != 3.");
	}

	if (
			(first_coordinate_indicies.template getSize<0>()  != out.template getSize<0>()) |
			(second_coordinate_indicies.template getSize<0>() != out.template getSize<1>())
			) {
		throw std::invalid_argument("out buffer of incorrect shape.");
	}

	typedef Eigen::Stride<Eigen::Dynamic, 1> CoordinateStride;
	typedef Eigen::Map< Eigen::Matrix<Real, 3, Eigen::Dynamic>, Eigen::Unaligned, CoordinateStride > CoordinateMap;

	size_t system_size = coordinates_per_entry;
	size_t first_system_stride = first_coordinates.template getStride<0>();
	size_t second_system_stride = second_coordinates.template getStride<0>();
	CoordinateStride first_coordinate_stride( first_coordinates.template getStride<1>(), 1);
	CoordinateStride second_coordinate_stride( second_coordinates.template getStride<1>(), 1);

	Eigen::Matrix<Real, 3, Eigen::Dynamic> first_coordinate_centers(3, first_coordinates.template getSize<0>());
	for ( size_t i = 0; i < first_coordinates.template getSize<0>(); i++ ) {
		first_coordinate_centers.col(i) = CoordinateMap(
			first_coordinates.getData() + i * first_system_stride,
			3, system_size, first_coordinate_stride).rowwise().sum() / system_size;
	}

	Eigen::Matrix<Real, 3, Eigen::Dynamic> second_coordinate_centers(3, second_coordinates.template getSize<0>());
	for ( size_t i = 0; i < second_coordinates.template getSize<0>(); i++ ) {
		second_coordinate_centers.col(i) = CoordinateMap(
			second_coordinates.getData() + i * second_system_stride,
			3, system_size, second_coordinate_stride).rowwise().sum() / system_size;
	}

	for ( size_t i = 0; i < first_coordinates.template getSize<0>(); i++ ) {
		for ( size_t j = 0; j < first_coordinates.template getSize<0>(); j++ ) {
			CoordinateMap first_c(
				first_coordinates.getData() + i * first_system_stride,
				3, system_size, first_coordinate_stride);

			CoordinateMap second_c(
				second_coordinates.getData() + j * second_system_stride,
				3, system_size, second_coordinate_stride);

			out(i, j) = QCPKernel<Real>::calc_coordinate_rmsd(
				first_c,
				first_coordinate_centers.col(i),
				second_c ,
				second_coordinate_centers.col(j));
		}
	}
}

}
}
