// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/AminoAcidNeighborSparseMatrix.hh
/// @brief  Amino acid neighbor sparse matrix template class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_AminoAcidNeighborSparseMatrix_hh
#define INCLUDED_core_pack_interaction_graph_AminoAcidNeighborSparseMatrix_hh

// Utility headers
#include <utility/vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.hh>

// Rosetta headers
#include <core/pack/interaction_graph/SparseMatrixIndex.hh>

namespace core {
namespace pack {
namespace interaction_graph {

// A sparse matrix for rotamer pair energies.
//
// Chris Saunders introduced the idea of reducing the rotamer pair
// energy table's memory usage by keeping a list for each residue pair
// of all their amino-acid pairs that have non-zero interaction energy.
// This sparse matrix did not allocate space for rotamer pairs that
// belonged to a pair of non-interacting-amino-acid pairs.  The sparse
// matrix representation Chris introduced when Rosetta was in Fortran
// has been extracted into its own C++ class so it can be easily reused
// in the various InteractionGraph classes.

template < typename T >
class AminoAcidNeighborSparseMatrix
{
public:

	typedef T value_type;


	/// @brief Constructor.  The amino acid neighbor sparse matrix stores the rotamer
	/// pair energies for a pair of residues - or more generally - for a pair of
	/// nodes in an interaction graph. Template type T must support:
	/// operator <
	/// operator =
	/// operator +=
	///
	/// The constructor requires references to the arrays holding state counts
	/// for each amino acid type.  Multiple sparse-matrices may refer to the same
	/// state-count-per-amino acid arrays.  By using proxy arrays instead of copying
	/// the state-counts, the sparse-matrix saves memory.
	///
	/// @param first_node_num_states_per_aa - [in] - reference to array with the
	///   number of states for the first dimension for each amino acid.
	/// @param second_node_num_states_per_aa - [in] - reference to array with
	///   the number of states for the second dimension for each amino acid.
	inline
	AminoAcidNeighborSparseMatrix
	(
		utility::vector1< int > const & first_node_num_states_per_aa,
		utility::vector1< int > const & second_node_num_states_per_aa
	) :
		num_aa_( first_node_num_states_per_aa.size() ),
		first_node_num_states_per_aatype_( first_node_num_states_per_aa ),
		second_node_num_states_per_aatype_( second_node_num_states_per_aa ),
		table_size_( 0 )
	{}

	/// @brief Method for telling the sparse matrix which amino acid pairs to allocate
	/// space for.  The input is a 2D FArray where the first dimension
	/// indexes over the amino-acids of the second node, and the second dimension
	/// indexes over the amino-acids of the first node.  This is somewhat iverted,
	/// but it is the natural consequence of a column-major 2D array.
	///
	/// @param sparse_conn_info - [in] - table of boolean values representing which
	///   amino acid pairs that neighbor.
	void
	set_sparse_aa_info
	(
		ObjexxFCL::FArray2_bool const & sparse_conn_info
	)
	{
		aa_offsets_.dimension( num_aa_, num_aa_ );
		aa_offsets_ = 0;
		int next_offset = 0;

		for ( int ii = 1; ii <= num_aa_; ++ii ) {
			int node1_num_states_for_aatype =
				first_node_num_states_per_aatype_[ ii ];
			for ( int jj = 1; jj <= num_aa_; ++jj ) {
				int node2_num_states_for_aatype =
					second_node_num_states_per_aatype_[ jj ];
				if ( ! sparse_conn_info( jj, ii) ||
						node1_num_states_for_aatype == 0 ||
						node2_num_states_for_aatype == 0 ) {
					aa_offsets_(jj, ii) = -1;
				} else {
					aa_offsets_(jj, ii) = next_offset;
					next_offset +=
						node1_num_states_for_aatype * node2_num_states_for_aatype;
				}
			}
		}

		table_size_ = next_offset;
		sparse_matrix_.dimension(table_size_);
		if ( next_offset != 0 ) {
			sparse_matrix_ = 0.0f;
		}
	}


	/// @brief returns true if node1aa and node2aa neighbor
	///
	/// @param node1aa - [in] - amino acid on node 1
	/// @param node2aa - [in] - amino acid on node 2
	inline
	bool get_sparse_aa_info( int node1aa, int node2aa) const
	{
		return (aa_offsets_(node2aa, node1aa) != -1);
	}

	/// @brief returns the number of entries in the sparse matrix
	inline
	int
	get_table_size() const
	{
		return table_size_;
	}


	/// @brief returns the number of bytes spent on the offset table
	unsigned int
	get_offset_table_size_in_bytes() const
	{
		return aa_offsets_.size() * sizeof( int );
	}


	/// @brief retrieves the value held for a pair of states.  State pairs without
	/// entries in the sparse matrix interact with zero energy.
	///
	/// @param ind1 - [in] - the SparseMatrixIndex for node1's state.
	/// @param ind2 - [in] - the SparseMatrixIndex for node2's state.
	inline
	value_type
	get( SparseMatrixIndex const & ind1, SparseMatrixIndex const & ind2 )
	const
	{
		int offset = get_offset( ind1, ind2);
		if ( offset == -1 ) return (value_type) 0;

		return sparse_matrix_( offset + get_submatrix_index(ind1, ind2) );
	}

	/// @brief cache efficient version of the energy lookup function. All private data
	/// for a particular sparse matrix has been held at a node, and now the node
	/// provides that data back to this static method.  Saves a considerable amount
	/// of time in simulated annealing - enough that this OO version of Chris's
	/// energy2b table is as fast as his was.
	///
	/// called when node2 is entertaining an alternate state during sim annealing
	///
	/// @param ind1 - [in] - the SparseMatrixIndex for node1's state.
	/// @param ind2 - [in] - the SparseMatrixIndex for node2's state.
	/// @param ind2num_states_per_aatype - [in] - for node2's state, how many other
	///   states does node2 have of the same amino acid type.
	/// @param aa_offset - [in] - offset into sparse matrix.  Node2 does not need
	///   to know what this offset represents; it only needs to know that it must
	///   provide this data to this method.
	/// @param sparse_matrix - [in] - the sparse matrix that this method reads from
	static
	inline
	value_type
	get
	(
		SparseMatrixIndex const & ind1,
		SparseMatrixIndex const & ind2,
		int ind2num_states_per_aatype,
		int aa_offset,
		ObjexxFCL::FArray1< value_type > const & sparse_matrix
	)
	{
		if ( aa_offset == -1 ) {
			return (value_type) 0;
		}

		int index = aa_offset +
			(ind2num_states_per_aatype *
			(ind1.get_state_ind_for_this_aa_type() - 1) ) +
			ind2.get_state_ind_for_this_aa_type();

		return sparse_matrix( index );
	}

	/// @brief cache efficient version of the energy lookup function. All private data
	/// for a particular sparse matrix has been held at a node, and now the node
	/// provides that data back to this static method.  Saves a considerable amount
	/// of time in simulated annealing - enough that this OO version of Chris's
	/// energy2b table is as fast as his was.
	///
	/// called when node1 is entertaining an alternate state during sim annealing
	///
	/// @param ind2 - [in] - the SparseMatrixIndex for node2's state.
	/// @param ind1_node_state_offset_minus_1 - [in] -
	/// @param ind2_num_states_per_aatype [in] -
	/// @param aa_offset - [in] - offset into sparse matrix.  Node1 does not need
	///   to know what this offset represents; it only needs to know that it must
	///   provide this data to this method.
	/// @param sparse_matrix - [in] - the sparse matrix that this method reads from
	///
	static
	inline
	value_type
	get
	(
		SparseMatrixIndex ind2,
		int ind1_node_state_offset_minus_1,
		int ind2_num_states_per_aatype,
		int aa_offset,
		ObjexxFCL::FArray1< value_type > const & sparse_matrix
	)
	{

		if ( aa_offset == -1 ) {
			return (value_type) 0;
		}

		int index = aa_offset +
			(ind2_num_states_per_aatype * ind1_node_state_offset_minus_1 ) +
			ind2.get_state_ind_for_this_aa_type();
		return sparse_matrix( index );
	}

	/// @brief how many values are stored in this sparse matrix?  Valid indices are 1 to size() for the
	/// operator [] method
	int
	size() const {
		return sparse_matrix_.size();
	}

	/// @brief accessor that does not use sparse matrix data -- ignorant of aa type for both dimensions.
	/// useful for when traversing all of the values stored in the sparse matrix.
	value_type
	operator [] ( int index ) const {
		return sparse_matrix_( index );
	}

	/// @brief non-const accessor for sparse matrix data -- ignorant of hte aa type for both dimensions.
	/// useful for when traversing all of the values stored in the sparse matrix.
	value_type &
	operator [] ( int index ) {
		return sparse_matrix_( index );
	}


	/// @brief stores a value for a pair of states.  Does not store anything for state
	/// pairs without entries in the sparse matrix.  Does not give a warning
	/// if you try to store a value where you cannot.  Maybe it should.
	///
	/// @param ind1 - [in] - the SparseMatrixIndex for node1's state.
	/// @param ind2 - [in] - the SparseMatrixIndex for node2's state.
	/// @param val - [in] - the value to store for this state pair
	inline
	void
	set
	(
		SparseMatrixIndex const & ind1,
		SparseMatrixIndex const & ind2,
		value_type const val
	)
	{
		int offset = get_offset( ind1, ind2);
		if ( offset == -1 ) return;

		sparse_matrix_( offset + get_submatrix_index( ind1, ind2) ) = val;
		return;
	}

	/// @brief
	/// cache efficient version of the energy lookup function. All private data
	/// for a particular sparse matrix has been held at a node, and now the node
	/// provides that data back to this static method.  Saves a considerable amount
	/// of time in simulated annealing - enough that this OO version of Chris's
	/// energy2b table is as fast as his was.
	///
	/// called when node2 is entertaining an alternate state during sim annealing
	///
	/// @param ind1 - [in] - the SparseMatrixIndex for node1's state.
	/// @param ind2 - [in] - the SparseMatrixIndex for node2's state.
	/// @param ind2num_states_per_aatype - [in] - for node2's state, how many other
	///   states does node2 have of the same amino acid type.
	/// @param aa_offset - [in] - offset into sparse matrix.  Node2 does not need
	///   to know what this offset represents; it only needs to know that it must
	///   provide this data to this method.
	/// @param sparse_matrix - [in] - the sparse matrix that this method reads from
	/// @param val - [in] - the value to store for this state pair
	static
	inline
	void
	set
	(
		SparseMatrixIndex const & ind1,
		SparseMatrixIndex const & ind2,
		int ind2num_states_per_aatype,
		int aa_offset,
		ObjexxFCL::FArray1< value_type > & sparse_matrix,
		value_type val
	)
	{
		if ( aa_offset == -1 ) {
			return;
		}

		int index = aa_offset +
			(ind2num_states_per_aatype *
			(ind1.get_state_ind_for_this_aa_type() - 1) ) +
			ind2.get_state_ind_for_this_aa_type();

		sparse_matrix( index ) = val;
	}

	/// @brief cache efficient version of the energy lookup function. All private data
	/// for a particular sparse matrix has been held at a node, and now the node
	/// provides that data back to this static method.  Saves a considerable amount
	/// of time in simulated annealing - enough that this OO version of Chris's
	/// energy2b table is as fast as his was.
	///
	/// called when node1 is entertaining an alternate state during sim annealing
	///
	/// @param ind2 - [in] - the SparseMatrixIndex for node2's state.
	/// @param ind1_node_state_offset_minus_1 - [in] -
	/// @param ind2_num_states_per_aatype [in] -
	/// @param aa_offset - [in] - offset into sparse matrix.  Node1 does not need
	///   to know what this offset represents; it only needs to know that it must
	///   provide this data to this method.
	/// @param sparse_matrix - [in] - the sparse matrix that this method reads from
	/// @param val - [in] - the value to store for this state pair
	static
	inline
	void
	set
	(
		SparseMatrixIndex ind2,
		int ind1_node_state_offset_minus_1,
		int ind2_num_states_per_aatype,
		int aa_offset,
		ObjexxFCL::FArray1< value_type > & sparse_matrix,
		value_type val
	)
	{

		if ( aa_offset == -1 ) {
			return;
		}

		int index = aa_offset +
			(ind2_num_states_per_aatype * ind1_node_state_offset_minus_1 ) +
			ind2.get_state_ind_for_this_aa_type();
		sparse_matrix( index ) = val;
	}

	// @brief Assigns all the entries in the sparse matrix to a single value.  Zeroed out
	// entries in the table remain zeroed out.
	inline
	void
	blanket_set
	(
		value_type val
	)
	{
		sparse_matrix_ = val;
	}

	/// @brief adds to a value for a pair of states.  Does not add anything for state
	/// pairs without entries in the sparse matrix.  Does not give a warning
	/// if you try to add to a value where you cannot.  Maybe it should.
	///
	/// @param ind1 - [in] - the SparseMatrixIndex for node1's state.
	/// @param ind2 - [in] - the SparseMatrixIndex for node2's state.
	/// @param val - [in] - the value to add for this state pair
	inline
	void
	add
	(
		SparseMatrixIndex const & ind1,
		SparseMatrixIndex const & ind2,
		value_type const val
	)
	{
		int offset = get_offset( ind1, ind2);
		if ( offset == -1 ) return;

		sparse_matrix_( offset + get_submatrix_index( ind1, ind2) ) += val;
		return;
	}

	/// @brief scale a particular energy by a constant. Danger: There's no undoing a scale-by-zero
	/// operation.
	inline
	void
	scale(
		SparseMatrixIndex const & ind1,
		SparseMatrixIndex const & ind2,
		value_type const scaler
	)
	{
		int offset = get_offset( ind1, ind2 );
		if ( offset == -1 ) return;

		sparse_matrix_( offset + get_submatrix_index( ind1, ind2 ) ) *= scaler;
	}

	/// @brief scale all energies by a constant.  Danger: There's no undoing a scale-by-zero
	/// operation.
	inline
	void
	scale(
		value_type const scaler
	)
	{
		sparse_matrix_ *= scaler;
	}

	/// @brief returns a reference to the first position in the sparse matrix
	/// to be used in a proxy array constructor.
	///
	/// used to offload the private data in this class onto a node for cache
	/// efficiency during simulated annealing.
	inline
	value_type &
	getMatrixPointer()
	{
		return sparse_matrix_(1);
	}

	/// @brief returns a reference to the first position in the amion acid pair offset
	/// matrix to be used in a proxy array constructor.
	///
	/// used to offload the private data in this class onto a node for cache
	/// efficiency during simulated annealing.
	///
	inline
	ObjexxFCL::FArray2D_int const &
	getAANeighborOffsets()
	{
		return aa_offsets_;
	}

	ObjexxFCL::FArray2D< value_type >
	get_aa_submatrix_energies(
		int node1aa,
		int node2aa
	) const
	{
		ObjexxFCL::FArray2D< value_type > submatrix(
			second_node_num_states_per_aatype_[ node2aa ],
			first_node_num_states_per_aatype_[ node1aa ]  );

		int offset = aa_offsets_( node2aa, node1aa );
		if ( offset == -1 ) {
			submatrix = value_type( 0 );
		} else {
			int const nvals = submatrix.size();
			for ( int li_src = offset, li_dest = 0; li_dest < nvals ; ++li_dest, ++li_src ) {
				submatrix[ li_dest ] = sparse_matrix_[ li_src ];
			}
		}

		return submatrix;
	}


	////////////////////////////////////////////////////////////////////////////////
	///
	/// @brief
	/// set this sparse-matrix to a completely-non-sparse state.  All entries
	/// are potentially non-zero.  Preserves the already-stored energies.
	///
	/// In some places, its easier to set the matrix to completely-non-sparse first
	/// and then to deallocate the zero-entry submatrices after computing all
	/// pair energies.
	///
	/// @global_read
	///
	/// @global_write
	///
	/// @remarks
	///
	/// @references
	///
	/// @author apl
	///
	////////////////////////////////////////////////////////////////////////////////
	inline
	void force_all_aa_neighbors()
	{
		if ( table_size_ != 0 ) {
			ObjexxFCL::FArray2D_int old_aa_offsets( aa_offsets_ ); //deep copy
			ObjexxFCL::FArray1D< value_type > old_sparse_matrix;
			old_sparse_matrix.swap( sparse_matrix_ );

			ObjexxFCL::FArray2D_bool new_aa_neighbors( num_aa_, num_aa_, true);

			set_sparse_aa_info(new_aa_neighbors);
			copy_old_data_into_new_table( old_aa_offsets, old_sparse_matrix );
		} else {
			ObjexxFCL::FArray2D_bool new_aa_neighbors( num_aa_, num_aa_, true);
			set_sparse_aa_info(new_aa_neighbors);
		}
	}

	/// @brief Sets a pair of amino acids as neighboring.  Preserves energies already
	/// stored in the table.
	inline
	void force_aa_neighbors( int node1aa, int node2aa )
	{
		if ( aa_offsets_( node2aa, node1aa ) != -1 ) return;

		ObjexxFCL::FArray2D_int old_aa_offsets( aa_offsets_ ); //deep copy
		ObjexxFCL::FArray1D< value_type > old_sparse_matrix;
		old_sparse_matrix.swap( sparse_matrix_ );

		ObjexxFCL::FArray2D_bool new_aa_neighbors( num_aa_, num_aa_, false );
		for ( int ii = 1; ii <= num_aa_; ++ii ) {
			for ( int jj = 1; jj <= num_aa_; ++jj ) {
				if ( old_aa_offsets(jj, ii) != -1 ) {
					new_aa_neighbors(jj, ii) = true;
				}
			}
		}
		new_aa_neighbors( node2aa, node1aa ) = true;

		set_sparse_aa_info( new_aa_neighbors );
		copy_old_data_into_new_table( old_aa_offsets, old_sparse_matrix );
	}

	/// @brief Deallocates amino-acid pair submatrices from the table that contain only
	/// zero-valued entries.
	void drop_zero_submatrices_where_possible()
	{
		drop_small_submatrices_where_possible( (value_type) 0 );
		return;
	}

	/// @brief Deallocates amino-acid pair submatrices from the table where the largest
	/// magnitude entry in the submatrix fails to exceed the threshold, epsilon.
	///
	/// @param epsilon - [in] - threshold for submatrix preservation
	void drop_small_submatrices_where_possible( value_type const epsilon )
	{
		ObjexxFCL::FArray2D_bool submatrix_worth_keeping( num_aa_, num_aa_, false);
		bool found_submatrix_not_worth_keeping = false;

		for ( int ii = 1; ii <= num_aa_; ++ii ) {
			int first_node_num_states_for_aa =
				first_node_num_states_per_aatype_[ ii ];
			if ( first_node_num_states_for_aa == 0 ) continue;

			for ( int jj = 1; jj <= num_aa_; ++jj ) {
				int second_node_num_states_for_aa =
					second_node_num_states_per_aatype_[ jj ];
				if ( second_node_num_states_for_aa == 0 ) continue;

				int aa_neighbor_offset = aa_offsets_( jj, ii );
				if ( aa_neighbor_offset != -1 ) {
					int submatrix_index = 1;
					bool this_submatrix_worth_keeping = false;
					for ( int kk = 1; kk <= first_node_num_states_for_aa; ++kk ) {
						for ( int ll = 1; ll <= second_node_num_states_for_aa; ++ll ) {
							int index = aa_neighbor_offset + submatrix_index;

							value_type energy_mag =
								std::abs( sparse_matrix_( index ) );
							if ( energy_mag > epsilon ) {
								this_submatrix_worth_keeping = true;
								break;
							}
							++submatrix_index;
						}
						if ( this_submatrix_worth_keeping ) break;
					}
					submatrix_worth_keeping( jj, ii ) = this_submatrix_worth_keeping;
					if ( ! this_submatrix_worth_keeping ) {
						found_submatrix_not_worth_keeping = true;
					}
				}
			}
		}

		if ( ! found_submatrix_not_worth_keeping ) return;

		ObjexxFCL::FArray1D_float old_sparse_matrix;
		old_sparse_matrix.swap( sparse_matrix_ );

		ObjexxFCL::FArray2D_int old_aa_offsets( aa_offsets_ ); //deep copy
		set_sparse_aa_info( submatrix_worth_keeping );
		copy_old_data_into_new_table( old_aa_offsets, old_sparse_matrix );

		return;
	}

	/// @brief delcare every submatrix to be not worth keeping
	/// @author Jack Maguire, jackmaguire1444@gmail.com
	void drop_all_submatrices()
	{
		ObjexxFCL::FArray2D_bool submatrix_worth_keeping( num_aa_, num_aa_, false);

		ObjexxFCL::FArray1D_float old_sparse_matrix;
		old_sparse_matrix.swap( sparse_matrix_ );

		ObjexxFCL::FArray2D_int old_aa_offsets( aa_offsets_ ); //deep copy
		set_sparse_aa_info( submatrix_worth_keeping );
		copy_old_data_into_new_table( old_aa_offsets, old_sparse_matrix );

		return;
	}

protected:

	/// @brief Lookup for offset into sparse table.  Prevents index-out-of-bounds error
	/// if either amino-acid types of the SparseMatrixIndices are 0.
	///
	/// @param ind1 - [in] - the sparse matrix index of node1
	/// @param ind2 - [in] - the sparse matrix index of node2
	inline
	int
	get_offset( SparseMatrixIndex const & ind1, SparseMatrixIndex const & ind2)
	const
	{
		if ( ind1.get_aa_type() == 0 || ind2.get_aa_type() == 0 ) return -1;

		return aa_offsets_( ind2.get_aa_type(), ind1.get_aa_type() );
	}

	/// @brief submatrix indexing function
	///
	/// @param ind1 - [in] - the sparse matrix index of node1
	/// @param ind2 - [in] - the sparse matrix index of node2
	inline
	int
	get_submatrix_index
	(
		SparseMatrixIndex const & ind1,
		SparseMatrixIndex const & ind2
	)
	const
	{
		return (ind1.get_state_ind_for_this_aa_type() - 1) *
			second_node_num_states_per_aatype_[ ind2.get_aa_type() ] +
			ind2.get_state_ind_for_this_aa_type();
	}

	/// @brief Assigns to the (freshly allocated) sparse_matrix_ member variable the
	/// values held in the input sparse matrix .  Used in both shrinking
	/// and growing the sparse matrix.
	///
	/// @param old_aa_offsets - [in] - aa offset data for the original table
	/// @param old_sparse_matrix - [in] - original table
	void copy_old_data_into_new_table
	(
		ObjexxFCL::FArray2D_int const & old_aa_offsets,
		ObjexxFCL::FArray1D< value_type > const & old_sparse_matrix
	)
	{
		for ( int ii = 1; ii <= num_aa_; ++ii ) {
			int first_node_num_states_for_aa =
				first_node_num_states_per_aatype_[ ii ];
			if ( first_node_num_states_for_aa == 0 ) continue;

			for ( int jj = 1; jj <= num_aa_; ++jj ) {
				int second_node_num_states_for_aa =
					second_node_num_states_per_aatype_[ jj ];
				if ( second_node_num_states_for_aa == 0 ) continue;

				int source_aa_offset = old_aa_offsets( jj, ii );
				int destination_aa_offset = aa_offsets_( jj, ii );

				if ( source_aa_offset != -1 && destination_aa_offset != -1 ) {
					int submatrix_index = 1;
					for ( int kk = 1; kk <= first_node_num_states_for_aa; ++kk ) {
						for ( int ll = 1; ll <= second_node_num_states_for_aa; ++ll ) {
							int index_destination = destination_aa_offset + submatrix_index;
							int index_source = source_aa_offset + submatrix_index;

							sparse_matrix_( index_destination ) =
								old_sparse_matrix( index_source );
							submatrix_index++;
						}
					}
				}
			}
		}
	}

	// Protected Member Data

	int const num_aa_;
	ObjexxFCL::FArray2D_int aa_offsets_;
	utility::vector1< int > const & first_node_num_states_per_aatype_;
	utility::vector1< int > const & second_node_num_states_per_aatype_;
	ObjexxFCL::FArray1D< value_type > sparse_matrix_;
	int table_size_;
};

} //end namespace interaction_graph
} //end namespace pack
} //end namespace core


#endif
