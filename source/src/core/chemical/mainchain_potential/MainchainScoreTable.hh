// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/mainchain_potential/MainchainScoreTable.hh
/// @brief  Headers for a general class for storing a torsional potential for mainchain resiudes.
/// @details Can be used by terms like rama, rama_prepro, p_aa_pp.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_chemical_mainchain_potential_MainchainTorsionPotential_hh
#define INCLUDED_core_chemical_mainchain_potential_MainchainTorsionPotential_hh

// Unit Headers
#include <core/chemical/mainchain_potential/MainchainScoreTable.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/pointer/deep_copy.hh>

// Numeric Headers
#include <numeric/MathNTensorBase.fwd.hh>
#include <numeric/interpolation/spline/PolycubicSplineBase.fwd.hh>
#include <numeric/interpolation/spline/BicubicSpline.fwd.hh>
#include <numeric/interpolation/spline/CubicSpline.fwd.hh>

// C++ Headers
#include <sstream>
#include <string> //getline overload

namespace core {
namespace chemical {
namespace mainchain_potential {

class MainchainScoreTable : public utility::pointer::ReferenceCount
{

public:

	/// @brief Default constructor.
	///
	MainchainScoreTable();

	/// @brief Clone function: make a copy of this object and return
	/// an owning pointer to the copy.
	MainchainScoreTableOP clone() const;

public: //Read functions:

	/// @brief Parse a Shapovalov-style rama database file and set up this MainchainScoreTable.
	/// @details Sets initialized_ to true.
	/// @param[in] filename The name of the file that was read.  (Just used for output messages -- this function does not file read).
	/// @param[in] file_contents The slurped contents of the file to parse.
	/// @param[in] res_type_name The name of the ResidueType for which we're reading data.  Data lines for other residue types will be
	/// ignored.
	/// @param[in] use_polycubic_interpolation If true, uses polycubic interpolation; if false, uses polylinear interpolation.
	void parse_rama_map_file_shapovalov(
		std::string const &filename,
		std::string const &file_contents,
		std::string const &res_type_name,
		bool const use_polycubic_interpolation
	);

public: //Functions for de novo generation

	/// @brief Reset this object completely.
	void clear();

	/// @brief Initialize this object for on-the-fly generation of a mainchain scoretable.
	/// @details Calls clear() first.
	void initialize_for_de_novo_generation( core::chemical::ResidueType const & restype, utility::vector1< core::Size > const & dimensions, utility::vector1< core::Size > const & mainchain_torsions_covered, bool const symmetrize );

	/// @brief Set an entry in the energies tensor.
	/// @note Bounds checking only in debug build!
	void set_energy( utility::vector1< core::Size > const & coords, core::Real const & energy_in );

	/// @brief Given coordinates in the energy tensor, increment the coordinates.
	/// @returns Returns "true" if increment was successful, "false" if the end of the tensor has been reached.
	bool increment_energy_coords( utility::vector1< core::Size > & coords ) const;

	/// @brief After the energies tensor has been populated, compute all derived data.
	/// @details If "normalize" is true, the probabilities are normalized to 1, and the energies adjusted accordingly.  (This
	/// has the effect of raising or lowering all of the energies by some constant.)  If "symmetrize_gly_" is true, the tensors are
	/// made symmetric.
	/// @note Sets initalized_ to true.
	void finalize_de_novo_scoretable_from_energies( core::Real const &kbt, bool const normalize );

	/// @brief After the mainchain scoretable has been finalized, write it to a scorefile.
	void write_mainchain_scoretable_to_stream( std::stringstream & outstream ) const;

public: //Accessor functions:

	/// @brief Has this MainchainScoreTable been initialized?
	/// @details False if no score table has been read yet; true otherwise.
	inline bool initialized() const { return initialized_; }

	/// @brief Given a vector of coordinates, get the corresponding vector of mainchain torsions.
	/// @details Bounds-checking only in debug mode!
	/// @note Offset is 0 if symmetrize_gly_ is false, 1/2 well if symmetrize_gly_ is true.  Offset is only added if offset_if_symmetric
	/// is set to true.
	void get_mainchain_torsions_from_coords( utility::vector1< core::Size > const & coords_in, utility::vector1< core::Real > & torsions_out, bool const offset_if_symmetric = true ) const;

	/// @brief Access an entry in the energy tensor.
	/// @details Note that bounds checking only occurs in debug mode!
	/// @note Intended for unit testing.
	core::Real const & energy_tensor( utility::vector1< core::Size > const & coords ) const;

	/// @brief Access values in this MainchainScoreTable.
	/// @details Note that the vector is deliberately not passed by reference.  The function copies the vector and ensures
	/// that all coordinates are in the range [0, 360).
	/// @note The full vector of mainchain torsions should be passed in.  If the potential is a function of fewer degrees of freedom,
	/// this function will disregard the appropraite entries in the coords vector.
	core::Real energy( utility::vector1 < core::Real > coords ) const;

	/// @brief Get the gradient with respect to x1,x2,x3,...xn for this MainchainScoreTable.
	/// @param[in] coords_in The coordinates at which to evaluate the gradient.
	/// @param[out] gradient_out The resulting gradient.
	/// @note The full vector of mainchain torsions should be passed in.  If the potential is a function of fewer degrees of freedom,
	/// this function will disregard the appropraite entries in the coords vector.
	void gradient( utility::vector1 < core::Real > coords_in, utility::vector1 < core::Real > &gradient_out ) const;

	/// @brief Set whether we should symmetrize tables for glycine.
	///
	void set_symmetrize_gly( bool const setting_in );

	/// @brief Return whether we should symmetrize tables for glycine.
	///
	inline bool symmetrize_gly() const { return symmetrize_gly_; }

	/// @brief Given the cumulative distribution function (pre-calculated), draw a random set of mainchain torsion values
	/// biased by the probability distribution.
	/// @details output is in the range (-180, 180].
	/// @note The dimensionality of the torsions vector will match the total degrees of freedom of the residue.  If the mainchain potential is
	/// a function of fewer mainchain degrees of freedom, then those torsions on which the potential does NOT depend will have values of 0 in the
	/// torsions vector.  Use the get_mainchain_torsions_covered() function to get the vector of mainchiain torsion indices that are covered.
	void draw_random_mainchain_torsion_values( utility::vector1 < core::Real > &torsions ) const;

	/// @brief Get a const reference to the vector of mainchain torsions indices that this mainchain potential covers.
	/// @details For example, for an oligourea, this would return {1, 2, 3}, since the Rama maps for oligoureas cover
	/// phi, theta, and psi (mainchain torsions 1, 2, and 3, respectively), but not mu or omega (mainchain torsions 4 and 5).
	inline utility::vector1< core::Size > const & get_mainchain_torsions_covered() const { return mainchain_torsions_covered_; }

private: //Private functions:

	/// @brief Sets the state of this MainchainScoreTable object to "initialized".
	/// @details Double-checks that it's not already initialized, and throws an error if it is.
	void set_initialized();

	/// @brief Given a residue type, set dimension_ and n_mainchain_torsions_total_, and resize all tensors appropriately.
	/// @details Assumes that all data have already been cleared.  Call clear() before invoking this method.
	/// @note The mainchain_torsions_covered vector is deliberately passed by copy.
	void set_dimension_and_ntorsions_for_restype( core::chemical::ResidueType const & restype, utility::vector1< core::Size > const & dimensions, utility::vector1< core::Size > mainchain_torsions_covered );

	/// @brief Set up the cumulative distribution function.
	/// @details The CDF is used for drawing random mainchain torsion values biased by the relative probabilities of mainchain torsion values.
	/// @note Each bin stores the probability of being in the current bin or an earlier bin (so the final bin should store a probability of 1).
	/// This differs from the convention used in Ramachandran.cc, but allows for a simpler drawing algorithm: I pick a uniformly-distributed random
	/// number from 0 to 1, loop through my bins, and stop when I get to a bin with a value greater than the value that I have drawn.
	/// @param[in] probs Tensor of probabilities.  Need not be normalized (sum to 1).
	/// @param[out] cdf Tensor for the cumulative distribution function.
	void set_up_cumulative_distribution_function( numeric::MathNTensorBaseCOP< core::Real > probs, numeric::MathNTensorBaseOP< core::Real > cdf  ) const;

	/// @brief Given a probabilities tensor, calculate the energies.
	/// @details Tensors must be the same size.  Contents of the probabilities tensor are overwritten.
	/// @param[out] energies Tensor of energies.
	/// @param[in] probs Tensor of probabilities.
	/// @param[in] kbt Boltzmann temperature (k_B*T), in Rosetta energy units.
	void
	energies_from_probs(
		numeric::MathNTensorBaseOP< core::Real > energies,
		numeric::MathNTensorBaseCOP< core::Real > probs,
		core::Real const &kbt
	) const;

	/// @brief Check that the stringstream doesn't have bad or eof status, and throw an error message if it does.
	///
	void check_linestream( std::istringstream const &linestream, std::string const &filename, bool const fail_on_eof=true ) const;

	/// @brief Initialize the energies_ and probabilities_ tensors to 0-containing N-tensors, of the
	/// dimensions given by the dimensions vector.
	void initialize_tensors( utility::vector1 < core::Size > dimensions );

	/// @brief Convert the energies from probabilities to Rosetta energy units, and add the
	/// entropic correction factor.
	void iteratively_correct_energy_tensor( core::Real const &entropy );

	/// @brief Given coordinates in the energy tensor, go to the next bin.
	/// @details As row ends are reached, the row resets to 1 and the next column is selected (and so forth down the dimensions).
	/// @returns Returns "true" if increment was successful, "false" if the end of the tensor has been reached.
	bool increment_coords( utility::vector1 < core::Size > &coords, numeric::MathNTensorBaseCOP< core::Real > tensor ) const;

	/// @brief Given a set of coordinates in a MathNTensor, get the opposite coordinates.
	/// @details For example, in a 2D 5x5 tensor, (1, 3) would yield an opposite of (3, 1).
	void get_opposite_coord( utility::vector1 < core::Size > const &coord, numeric::MathNTensorBaseCOP< core::Real > tensor, utility::vector1 <core::Size> &opposite_coord ) const;

	/// @brief Once the internal MathNTensor has been set up, set up polycubic interpolation.
	/// @details This function includes special-case logic for setting up cubic interpolation in the 1D case and
	/// bicubic interpolation in the 2D case, since these are not handled by the PolycubicSpline class.
	/// @param[in] offsets Vector of offset values, from 0 to 1 -- where centres are, as fraction of bin width.
	/// @param[in] dimensions Vector of number of entries in each dimension.  Bin widths are inferred from this: 36 entries would correspond to 10-degree bins.
	void set_up_polycubic_interpolation(
		utility::vector1 < core::Real > const &offsets,
		utility::vector1 < core::Size > const &dimensions
	);

	/// @brief Given a tensor, symmetrize it.
	/// @details Assumes that tensor stores probabilities; normalizes tensor in the process.
	void symmetrize_tensor(
		numeric::MathNTensorBaseOP< core::Real > tensor
	) const;

private: //Private variables:

	/// @brief Has this object been initialized?
	///
	bool initialized_;

	/// @brief Dimensionality of this MainchainScoreTable.
	/// @details Minimum 1, maximum 9.  A value of 0 indicates that it is uninitialized.
	core::Size dimension_;

	/// @brief N-dimensional tensor for storing energies data.
	///
	utility::pointer::DeepCopyOP< numeric::MathNTensorBase< core::Real > > energies_;

	/// @brief N-dimensional tensor for storing probabilities data.
	///
	utility::pointer::DeepCopyOP< numeric::MathNTensorBase< core::Real > > probabilities_;

	/// @brief N-dimensional tensor for storing the cumulative distribution function.
	/// @details This is used for drawing random mainchain torsion values biased by the relative
	/// probabilities of having a set of mainchain torsion values.
	utility::pointer::DeepCopyOP< numeric::MathNTensorBase< core::Real > > cdf_;

	/// @brief Is this MainchainScoreTable set up with polycubic interpolation?
	/// @details Default true.  If false, interpolation is linear.
	bool use_polycubic_interpolation_;

	/// @brief Interpolation spline for the 1D case.
	/// @details Only used if use_polycubic_interpolation_ is true and the energies_ MathNTensor is 1D.
	/// Null otherwise.
	utility::pointer::DeepCopyOP< numeric::interpolation::spline::CubicSpline > energies_spline_1D_;

	/// @brief Interpolation spline for the 2D case.
	/// @details Only used if use_polycubic_interpolation_ is true and the energies_ MathNTensor is 2D.
	/// Null otherwise.
	utility::pointer::DeepCopyOP< numeric::interpolation::spline::BicubicSpline > energies_spline_2D_;

	/// @brief Interpolation spline for the N-dimensional case, where N > 2.
	/// @details Only used if use_polycubic_interpolation_ is true and the energies_ MathNTensor is N-dimensional, where N > 2.
	/// Null otherwise.
	utility::pointer::DeepCopyOP< numeric::interpolation::spline::PolycubicSplineBase > energies_spline_ND_;

	/// @brief What is the total number of mainchain torsions for this residue type?
	/// @details Could be different than the dimension of the tensors, if certain torsions are fixed.  Note that this excludes omega (the inter-residue torsion),
	/// so for an alpha-amino acid, it should be "2".  In the case of an oligourea, it is "3" instead of "4", since mu is fixed at 180 and enforced (along with omega)
	/// by the OmegaTetherEnergy.
	core::Size n_mainchain_torsions_total_;

	/// @brief If only a subset of mainchain torsions are provided, which ones are the relevant ones?
	/// @details For oligoureas, for example, this would be {1, 2, 3}, since only phi, theta, and psi are covered by the mainchain potential, and mu and theta
	/// (torsions 4 and 5) are covered by the omega score term.
	utility::vector1< core::Size > mainchain_torsions_covered_;

	/// @brief Symmetrize glycine tables?
	/// @details Read from option system by default.
	bool symmetrize_gly_;

};

} //mainchain_potential
} //chemical
} //core

#endif //INCLUDED_core_chemical_mainchain_potential_MainchainTorsionPotential_hh
