// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_protocols_canonical_sampling_mc_convergence_checks_pool_util_hh
#define INCLUDED_protocols_canonical_sampling_mc_convergence_checks_pool_util_hh

#include <ObjexxFCL/FArray2D.hh>
#include <utility/vector1.hh>
#include <core/types.hh>

namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

	typedef ObjexxFCL::FArray2D<double> FArray2D_double;
	typedef utility::vector1< core::Size > Address;

	struct DataBuffer {
	public:
		DataBuffer();

		~DataBuffer();

		void
		setup( int num_slave_nodes, int nresidues, int nlevels );

		void
		address_to_buf( utility::vector1< core::Size > & address, int* buf, core::Size start_index );

		void
		buf_to_address( utility::vector1< core::Size > & address, int* buf, core::Size start_index );

		void
		farray_to_array( core::Size index,
										 FArray2D_double const& coords,
										 double* coord_buf );

		void
		farray_to_array( core::Size index,
										 core::Size num_to_add,
										 FArray2D_double const& coords,
										 double* coord_buf );

		void
		array_to_farray( core::Size index,
										 FArray2D_double & coords,
										 double* coord_buf );

		void
		array_to_farray( core::Size index,
										 core::Size num_to_add,
										 FArray2D_double & coords,
										 double* coord_buf);

		int* neighbor_addresses_;
		double* coords_transfer_buffer_;
		double* coords_receiving_buffer_;
		FArray2D_double coords_;
		FArray2D_double temp_coords_;
		core::Size num_new_neighbors_;
		int* memory_offset_;
		int* int_buf1_;
		int* winning_ranks_;
		//structures for dealing with neighbors
		FArray2D_double candidate_coords_;
		core::Size candidate_nbr_index_;
		core::Real candidate_best_rmsd_;
		double* candidate_best_rmsds_;
		Address candidate_address_;
		std::string winning_tag_;
		Address winning_address_;
		core::Size new_level_begins_;
		utility::vector1< core::Real > best_candidate_rmsds_;
		utility::vector1< bool > is_a_neighbor_;
		int* finished_;

	};


}
}
}

#endif
