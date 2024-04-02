// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_mpnn/util.cc
/// @brief Utility functions for ProteinMPNN
/// @author Brian Koepnick (koepnick@uw.edu)

#ifdef USE_TORCH

#include <protocols/protein_mpnn/util.hh>

// Core headers:
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

// Basic headers:
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <numeric/xyzVector.hh>

// Utility headers:
#include <utility/exit.hh>
#include <utility/pointer/memory.hh>

#include <torch/csrc/autograd/generated/variable_factories.h>
#include <torch/script.h>
#include <algorithm>

static basic::Tracer TR( "protocols.protein_mpnn.util" );

typedef int64_t TorchLong; // because a Win64 `long` is only 32 bits

namespace protocols {
namespace protein_mpnn {

using namespace core;
using namespace core::pose;

torch::Tensor
bb_coords_to_tensor( const core::conformation::Conformation & conf, core::Size batch_size ){

	core::Size nres = conf.size();
	std::vector< float > bb_coords_vec(batch_size * nres * 4 * 3, 0.0f);
	float * bb_coords = &(bb_coords_vec[0]);
#define BB_IDX(i,j,k,l) (((i) * nres * 4 * 3) + ((j) * 4 * 3) + ((k) * 3) + (l))

	// fill array with atom coordinates of backbone atoms
	for( core::Size batch = 1; batch <= batch_size; ++batch ){
		for( core::Size resi = 1; resi <= nres; ++resi ){

			// if not protein residue, just fill with origin coords
			if( !conf.residue_type( resi ).is_protein() ){
				for( core::Size idx = 1; idx <= 4; ++idx ){
					for( core::Size dim = 0; dim < 3; ++dim ){
						bb_coords[ BB_IDX( batch - 1, resi - 1, idx - 1, dim ) ] = 0.0;
					}
				}
			} else {
				// fill coords for first four atom indices
				for( core::Size idx = 1; idx <= 4; ++idx ){
					core::id::AtomID atm( idx, resi );
					PointPosition xyz = conf.xyz( atm );
					for( core::Size dim = 0; dim < 3; ++dim ){
						bb_coords[ BB_IDX( batch - 1, resi - 1, idx - 1, dim ) ] = static_cast<float>(xyz[dim]);
					}
				}
			}
		}
	}
#undef BB_IDX

	// get a tensor from the coordinate array
	torch::TensorOptions options = torch::TensorOptions().dtype( torch::kFloat32 );
	torch::Tensor coords_t = torch::from_blob( bb_coords, { static_cast< int >( batch_size ), static_cast< int >( nres ), 4, 3 }, options ).clone();
	return coords_t;
}

torch::Tensor
randn_to_tensor( core::Size num_res, core::Size batch_size ){
	int nres = static_cast< int >( num_res );
	int bsize = static_cast< int >( batch_size );
	torch::TensorOptions options = torch::TensorOptions().dtype( torch::kFloat32);
	return torch::randn( { bsize, nres }, options );
}

torch::Tensor
seq_to_tensor( const std::string & seq, core::Size batch_size ){
	int bsize = static_cast< int >( batch_size );

	std::vector< TorchLong > seq_aa_vec(seq.size(), 0);
	TorchLong * seq_aa = &(seq_aa_vec[0]);
	for( core::Size ii = 0; ii < seq.size(); ++ii ){
		core::Size s = AA_ALPHABET.find( seq[ ii ] );
		if( s == std::string::npos ){
			TR.Error << "Unrecognized AA: " << seq[ ii ] << std::endl;
			// Raise exception??
		} else {
			seq_aa[ ii ] = s;
		}
	}

	torch::TensorOptions options = torch::TensorOptions().dtype( torch::kLong );
	torch::Tensor seq_t = torch::from_blob( seq_aa, { bsize, static_cast< int >( seq.size() ) }, options ).clone();
	return seq_t;
}


torch::Tensor
chain_encoding_to_tensor( const core::pose::Pose & pose, core::Size batch_size ){
	int nres = static_cast< int >( pose.total_residue() );
	int bsize = static_cast< int >( batch_size );

	// assign chainID to each residue position
	std::vector< TorchLong > chain_encoding_vec(nres, 0);
	TorchLong * chain_encoding_array = &(chain_encoding_vec[0]);
	for( int ii = 0; ii < nres; ++ii ){
		chain_encoding_array[ ii ] = static_cast< int >( pose.chain( ii+1 ) );
	}

	torch::TensorOptions options = torch::TensorOptions().dtype( torch::kLong );
	return torch::from_blob( chain_encoding_array, { bsize, nres }, options ).clone();
}

torch::Tensor
residue_idx_to_tensor( const core::pose::Pose & pose, core::Size batch_size ){
	int nres = static_cast< int >( pose.total_residue() );
	int bsize = static_cast< int >( batch_size );

	// iterate over residues and add 100 to the residue index for every chainbreak
	std::vector< TorchLong > residue_idx_vec(nres, 0);
	TorchLong * residue_idx_array = &(residue_idx_vec[0]);
	for( int ii = 0; ii < nres; ++ii ){
		int chain_offset = 100 * ( pose.chain( ii+1 ) - 1 );
		residue_idx_array[ ii ] = ii + chain_offset;
	}

	torch::TensorOptions options = torch::TensorOptions().dtype( torch::kLong );
	return torch::from_blob( residue_idx_array, { bsize, nres }, options ).clone();

}

torch::Tensor
coord_mask_to_tensor( const utility::vector1_bool & mask, core::Size batch_size ){
	int nres = static_cast< int >( mask.size() );
	int bsize = static_cast< int >( batch_size );

	// convoluted conversion of vector to array to tensor
	std::vector< float > mask_vec(nres, 0.0f);
	float * mask_array = &(mask_vec[0]);
	for( int ii = 0; ii < nres; ++ii ){
		mask_array[ ii ] = mask[ ii+1 ] ? 1 : 0;
	}

	torch::TensorOptions options = torch::TensorOptions().dtype( torch::kFloat32 );
	return torch::from_blob( mask_array, { bsize, nres }, options ).clone();
}

torch::Tensor
chain_mask_to_tensor( const core::pose::Pose & pose, const utility::vector1_bool & chain_mask, core::Size batch_size ){
	int nres = static_cast< int >( pose.total_residue() );
	int bsize = static_cast< int >( batch_size );

	// convoluted conversion of vector to array to tensor
	std::vector< float > chain_mask_vec(nres, 0.0f);
	float * chain_mask_array = &(chain_mask_vec[0]);
	for( int ii = 0; ii < nres; ++ii ){
		core::Size chain_id = pose.chain( ii+1 );
		chain_mask_array[ ii ] = chain_mask[ chain_id ] ? 1 : 0;
	}

	torch::TensorOptions options = torch::TensorOptions().dtype( torch::kFloat32 );
	return torch::from_blob( chain_mask_array, { bsize, nres }, options ).clone();
}

torch::Tensor
pos_mask_to_tensor( const utility::vector1_bool & pos_mask, core::Size batch_size ){
	int nres = static_cast< int >( pos_mask.size() );
	int bsize = static_cast< int >( batch_size );

	// convoluted conversion of vector to array to tensor
	std::vector< float > pos_mask_vec(nres, 0.0f);
	float * pos_mask_array = &(pos_mask_vec[0]);
	for( int ii = 0; ii < nres; ++ii ){
		pos_mask_array[ ii ] = pos_mask[ ii+1 ] ? 1 : 0;
	}

	torch::TensorOptions options = torch::TensorOptions().dtype( torch::kFloat32 );
	return torch::from_blob( pos_mask_array, { bsize, nres }, options ).clone();
}

torch::Tensor
omit_AAs_to_tensor( utility::vector1< char > omit_AAs ){
	int alphabet_size = static_cast< int >( AA_ALPHABET.size() );
	std::vector< float > omit_AAs_vec(alphabet_size, 0.0f);
	float * omit_AAs_array = &(omit_AAs_vec[0]);

	// iterate over alphabet and check if each AA is in the omit list
	for( int ii = 0; ii < alphabet_size; ++ii ){
		char aa = AA_ALPHABET[ ii ];
		if( omit_AAs.empty() ||
		   ( std::find( omit_AAs.begin(), omit_AAs.end(), aa ) == omit_AAs.end() ) ){
			omit_AAs_array[ ii ] = 0; // aa not found in omit_AAs
		} else {
			omit_AAs_array[ ii ] = 1; // aa found in omit_AAs
		}
	}

	torch::TensorOptions options = torch::TensorOptions().dtype( torch::kFloat32 );
	return torch::from_blob( omit_AAs_array, { alphabet_size }, options ).clone();
}

torch::Tensor
omit_AAs_pos_to_tensor( const utility::vector1< utility::vector1< char > > & omit_AAs_pos, const core::Size batch_size ){
	int alphabet_size = static_cast< int >( AA_ALPHABET.size() );
	int nres = static_cast< int >( omit_AAs_pos.size() );
	int bsize = static_cast< int >( batch_size );

	std::vector< float > omit_AAs_pos_vec( bsize * nres * alphabet_size, 0.0f);
	float * omit_AAs_pos_array = &(omit_AAs_pos_vec[0]);

	int write_head = 0;
	for( int res_idx = 1; res_idx <= nres; ++res_idx ){
		utility::vector1< char > res_omit_AAs = omit_AAs_pos[ res_idx ];
		for( char aa_type : AA_ALPHABET ){
			if( res_omit_AAs.empty() ||
				( std::find( res_omit_AAs.begin(), res_omit_AAs.end(), aa_type ) == res_omit_AAs.end() ) ){
				omit_AAs_pos_array[ write_head ] = 0.0f; // AA not found in res_omit_AAs
			} else {
				omit_AAs_pos_array[ write_head ] = 1.0f; // AA found in res_omit_AAs
			}
			++write_head;
		}
	}

	torch::TensorOptions options = torch::TensorOptions().dtype( torch::kFloat32 );
	return torch::from_blob( omit_AAs_pos_array, { bsize, nres, alphabet_size }, options ).clone();
}

torch::Tensor
bias_AAs_to_tensor( utility::vector1< core::Real > bias_AAs ){
	core::Size alphabet_size = AA_ALPHABET.size();
	std::vector< float > bias_AAs_vec(alphabet_size, 0.0f);
	float * bias_AAs_array = &(bias_AAs_vec[0]);

	runtime_assert( bias_AAs.size() == alphabet_size );

	// transpose input vector to float array
	for( core::Size ii = 0; ii < alphabet_size; ++ii ){
		bias_AAs_array[ ii ] = static_cast< float >( bias_AAs[ ii+1 ] );
	}

	torch::TensorOptions options = torch::TensorOptions().dtype( torch::kFloat32 );
	return torch::from_blob( bias_AAs_array, { static_cast< int >( alphabet_size ) }, options ).clone();
}

TiedPositions
convert_tied_positions( utility::vector1< utility::vector1< core::Size > > tied_positions ){
	TiedPositions list_list;

	// transpose input vector into torch::List
	for( core::Size ii = 1; ii <= tied_positions.size(); ++ii ){
		utility::vector1< core::Size > tied_set = tied_positions[ ii ];
		torch::List< int64_t > list;
		for( core::Size jj = 1; jj <= tied_set.size(); ++jj ){
			// residue numbering of output starts at zero
			list.push_back( static_cast< int64_t >( tied_set[ jj ] ) - 1 );
		}
		list_list.push_back( list );
	}
	return list_list;
}

utility::vector1< std::string >
seqs_from_tensor( const torch::Tensor & seq_t, core::Size num_res, core::Size batch_size  ){
	int nres = static_cast< int >( num_res );
	int bsize = static_cast< int >( batch_size );
	utility::vector1< std::string > seqs;
	for( int nbatch = 0; nbatch < bsize; ++nbatch ){
		std::string seq = "";
		for( int ii = 0; ii < nres; ++ii ){
			int idx = seq_t[0][ ii ].item< int >();
			seq += AA_ALPHABET[ idx ];
		}
		seqs.push_back( seq );
	}
	return seqs;
}

} //protein_mpnn
} //protocols

#endif //USE_TORCH
