// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/linear_algebra/GeneralizedEigenSolver.hh
/// @brief  Header file for GeneralizedEigenSolver.
/// @author Kale Kundert

#ifndef INCLUDED_numeric_linear_algebra_GeneralizedEigenSolver_HH
#define INCLUDED_numeric_linear_algebra_GeneralizedEigenSolver_HH

// Unit Headers
#include <numeric/linear_algebra/rgg.hh>

// External Headers
#include <Eigen/Dense>

// C++ headers
#include <string>
#include <vector>

namespace numeric {
namespace linear_algebra {

// Class Description {{{1

/// @class GeneralizedEigenSolver
/// @brief Solves generalized eigenvalue problems.
///
/// @tparam _MatrixType The types of the two matrices involved in the 
/// generalized eigenvalue problem.  The template parameter is expected to be 
/// an instantiation of the Eigen::Matrix class template.  Only real matrices 
/// are supported.
///
/// The generalized matrix problem is \f$Ax = \lambda Bx\f$.  Here \f$A\f$ and 
/// \f$B\f$ are real matrices, \f$x\f$ is an eigenvector, \f$\lambda\f$ is a 
/// scalar eigenvalue.  The problem is to find pairs of \f$x\f$ and 
/// \f$\lambda\f$ that satisfy the above equation.  Note that the standard 
/// eigenvalue problem \f$Ax = \lambda x\f$ is obtained by setting \f$B = I\f$.  
/// However, for this problem you should use the dedicated `EigenSolver` 
/// provided by the `Eigen` library.
///
/// Call the compute() method to solve a generalized eigenvalue problem.  
/// Alternatively, the \f$A\f$ and \f$B\f$ matrices can be simply passed to the 
/// constructor to solve the problem immediately after construction.  Once a 
/// solution has been found, call the eigenvalues() and eigenvectors() methods 
/// to access the results.  Note that even though \f$A\f$ and \f$B\f$ must be 
/// real, the eigenvalues and eigenvectors that solve the problem may be 
/// complex.  The results returned by this solver are of type `std::complex`, 
/// although the imaginary component may sometimes be null.
///
/// This solver is a thin wrapper around the RGG routine provided by EISPACK.  
/// The RGG routine itself was converted from Fortran to C++ using fable.

// }}}1

template<typename _MatrixType>
class GeneralizedEigenSolver {

// Type Definitions {{{1
public:

	/// @brief Alias for the template parameter `_MatrixType`.
	typedef _MatrixType MatrixType;

	enum {
		RowsAtCompileTime = MatrixType::RowsAtCompileTime,
		ColsAtCompileTime = MatrixType::ColsAtCompileTime,
		Options = MatrixType::Options,
		MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
		MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
	};

	/// @brief Alias for the type used to index matrices of `MatrixType`.
	typedef typename MatrixType::Index Index;

	/// @brief Alias for the scalar type used in `MatrixType` matrices.  Note 
	/// that this must be a real type for this algorithm to work.
	typedef typename MatrixType::Scalar Scalar;

	/// @brief Alias for the real component of `ScalarType`.
	typedef typename Eigen::NumTraits<Scalar>::Real RealScalar;

	/// @brief Complex type based on `RealScalar`.
	typedef std::complex<RealScalar> ComplexScalar;

	/// @brief Complex vector type used to represent the calculated eigenvalues.
	typedef Eigen::Matrix<ComplexScalar, ColsAtCompileTime, 1, Options,
					MaxColsAtCompileTime, 1> EigenvalueType;

	/// @brief Complex matrix type used to represent the calculated eigenvectors.
	typedef Eigen::Matrix<ComplexScalar, RowsAtCompileTime, ColsAtCompileTime, 
					Options, MaxRowsAtCompileTime, MaxColsAtCompileTime> EigenvectorType;

	/// @brief Vector type used to represent the real subset of calculated 
	/// eigenvalues.
	typedef Eigen::Matrix<RealScalar, Eigen::Dynamic, 1, Options, 
					MaxColsAtCompileTime, 1> RealEigenvalueType;

	/// @brief Matrix type used to represent the real subset of calculated 
	/// eigenvectors.
	typedef Eigen::Matrix<RealScalar, RowsAtCompileTime, Eigen::Dynamic, Options, 
					MaxRowsAtCompileTime, MaxColsAtCompileTime> RealEigenvectorType;

private:

	/// @brief Vector type meant to be compatible with RGG.  Sometimes used to 
	/// implicitly hold imaginary numbers, even though it is a real number type. 
	typedef Eigen::Matrix<RealScalar, ColsAtCompileTime, 1, Options & 
					~Eigen::RowMajor, MaxColsAtCompileTime, 1> ScratchVector;

	/// @brief Matrix type meant to be compatible with RGG.  Sometimes used to 
	/// implicitly hold imaginary numbers, even though it is a real number type. 
	typedef Eigen::Matrix<RealScalar, RowsAtCompileTime, ColsAtCompileTime, 
					Options & ~Eigen::RowMajor, MaxRowsAtCompileTime, 
					MaxColsAtCompileTime> ScratchMatrix;

// Constructors, Destructors, etc. {{{1
public:

	/// @brief Default constructor.
	GeneralizedEigenSolver()
			: eigenvalues_(),
				eigenvectors_(),
				real_indices_(),
				eigenvalues_were_computed_(false),
				eigenvectors_were_computed_(false),
				info_(Eigen::Success) {}

	/// @brief Default constructor with memory preallocation.
	/// @param[in] size The amount of space to preallocate for this problem.
	/// @see GeneralizedEigenSolver()
	GeneralizedEigenSolver(Index size)
			: eigenvalues_(size),
				eigenvectors_(size, size),
				real_indices_(),
				eigenvalues_were_computed_(false),
				eigenvectors_were_computed_(false),
				info_(Eigen::Success) {

		real_indices_.reserve(size);
	}

	/// @brief Construct and solve the given generalized eigenvalue problem.
	///
	/// @param[in] A Matrix on the left-hand side of the eigenvalue equation.
	/// @param[in] B Matrix on the right-hand side of the eigenvalue equation.
	/// @param[in] compute_eigenvectors Indicates whether or not the eigenvectors 
	/// should be calculated.  The eigenvalues are always calculated.
	///
	/// @details This constructor calls compute() to solve the given problem.
	/// @see GeneralizedEigenSolver()
	GeneralizedEigenSolver(
			MatrixType const & A,
			MatrixType const & B,
			bool compute_eigenvectors = true)

			: eigenvalues_(A.cols()),
				eigenvectors_(A.rows(), A.cols()),
				real_indices_(),
				eigenvalues_were_computed_(false),
				eigenvectors_were_computed_(false),
				info_(Eigen::Success) {

		real_indices_.reserve(A.cols());
		compute(A, B, compute_eigenvectors);
	}

// Public Methods {{{1
public:

	/// @brief Solve the given generalized eigenvalue problem.
	GeneralizedEigenSolver & compute(
			MatrixType const & A,
			MatrixType const & B,
			bool compute_eigenvectors = true);

	/// @brief Solve the given generalized eigenvalue problem in place.
	GeneralizedEigenSolver & compute_in_place(
			MatrixType & A,
			MatrixType & B,
			bool compute_eigenvectors = true);

	/// @brief Return the calculated (possibly complex) eigenvalues.
	/// @see compute()
	EigenvalueType eigenvalues() const {
		eigen_assert(eigenvalues_were_computed_ &&
				"No eigenvalues were computed.");
		return eigenvalues_;
	}

	/// @brief Return the calculated (possibly complex) eigenvectors.
	///
	/// @details The returned eigenvectors will be normalized.  This method can 
	/// only be called if the eigenvectors were actually calculated (i.e.  
	/// `compute_eigenvectors = true`) in the call to compute().  An assertion 
	/// will trigger if this is not the case.
	///
	/// @see compute()
	EigenvectorType eigenvectors() const {
		eigen_assert(eigenvectors_were_computed_ &&
				"No eigenvectors were computed.");
		return eigenvectors_;
	}

	/// @brief Return any real eigenvalues that were calculated.
	RealEigenvalueType real_eigenvalues() const;

	/// @brief Return any real eigenvectors that were calculated.
	RealEigenvectorType real_eigenvectors() const;

	/// @brief Return the number of real solutions that were calculated.
	/// @see compute()
	int num_real_solutions() const {
		eigen_assert(eigenvalues_were_computed_ &&
				"The compute() method was never called.");
		return real_indices_.size();
	}

	/// @brief Indicate whether or not the eigenvalue problem was successfully 
	/// solved.
	///
	/// @returns `Eigen::Success` if computation was successful, 
	/// `Eigen::NoConvergence` otherwise.
	///
	/// @see compute()
	Eigen::ComputationInfo info() const {
		eigen_assert(eigenvalues_were_computed_ &&
				"The compute() method was never called.");
		return info_;
	}


// Private Methods {{{1
private:

	/// @brief Unpacks the eigenvalues returned by RGG and repacks them into a 
	/// format that better supports complex numbers.
	void extract_eigenvalues(
			ScratchVector eigenvalue_numerators_real,
			ScratchVector eigenvalue_numerators_imaginary,
			ScratchVector eigenvalue_denominators);

	/// @brief Unpacks the eigenvectors returned by RGG and repacks them into a 
	/// format that better supports complex numbers.
	void extract_eigenvectors(
			ScratchMatrix eigenvectors);

	/// @brief Determines which eigenvalues are real (i.e. those for which the 
	/// imaginary component is zero).
	void extract_real_indices();

// Data Members {{{1
private:

	/// @brief Vector containing the calculated eigenvalues.
	EigenvalueType eigenvalues_;
	
	/// @brief Matrix containing the calculated and normalized eigenvectors.
	EigenvectorType eigenvectors_;

	/// @brief List containing the indices of any real eigenvalues.
	std::vector<Index> real_indices_;

	/// @brief Flag indicating if compute() has been called yet.
	bool eigenvalues_were_computed_;

	/// @brief Flag indicating if eigenvectors were calculated along with the 
	/// eigenvalues when compute() was called.
	bool eigenvectors_were_computed_;

	/// @brief Indicates whether or not the calculation was successful.
	Eigen::ComputationInfo info_;

// }}}1

};

// GeneralizedEigenSolver::compute {{{1

/// @param[in] A Matrix on the left-hand side of the eigenvalue equation.
/// @param[in] B Matrix on the right-hand side of the eigenvalue equation.
/// @param[in] compute_eigenvectors Indicates whether or not the eigenvectors 
/// should be calculated.  The eigenvalues are always calculated.
///
/// @details This method finds solutions to the problem \f$Ax = \lambda Bx\f$.  
/// The \f$A\f$ and \f$B\f$ matrices are not affected by this calculation.  
/// Because this behavior is intuitive, it is the default behavior.  However, 
/// it is also somewhat less efficient.  The underlying EISPACK routine uses 
/// the \f$A\f$ and \f$B\f$ matrices for scratch space, so this method has to 
/// copy the two matrices to prevent them from being changed.  Of course, this 
/// copy costs both time and memory.  If that time and memory is important to 
/// you, use the compute_in_place() method instead.
///
/// \see compute_in_place()
/// \see rgg()

template<typename MatrixType>
GeneralizedEigenSolver<MatrixType>&
GeneralizedEigenSolver<MatrixType>::compute(
		MatrixType const & A,
		MatrixType const & B,
		bool compute_eigenvectors) {

	MatrixType mutable_A (A), mutable_B(B);
	return compute_in_place(mutable_A, mutable_B, compute_eigenvectors);
}

// GeneralizedEigenSolver::compute_in_place {{{1

/// @param[in] A Matrix on the left-hand side of the eigenvalue equation.
/// @param[in] B Matrix on the right-hand side of the eigenvalue equation.
/// @param[in] compute_eigenvectors Indicates whether or not the eigenvectors 
/// should be calculated.  The eigenvalues are always calculated.
///
/// @details This method finds solutions to the problem \f$Ax = \lambda Bx\f$.  
/// The computation is performed in place, which save some resources relative 
/// to the compute() method.  However, the \f$A\f$ and \f$B\f$ matrices are 
/// trashed in the process.  Usually the compute() method is a more appropriate 
/// choice because the modest improvement in performance is not worth the 
/// counterintuitive behavior.
///
/// \see compute()
/// \see rgg()

template<typename MatrixType>
GeneralizedEigenSolver<MatrixType>&
GeneralizedEigenSolver<MatrixType>::compute_in_place(
		MatrixType & A,
		MatrixType & B,
		bool compute_eigenvectors) {

	using namespace fem::major_types;

	// Make sure the arguments are sane.
	eigen_assert(A.rows() == A.cols());
	eigen_assert(B.rows() == B.cols());
	eigen_assert(A.rows() == B.rows());
	eigen_assert(A.cols() == B.cols());

	// Note that this method has been called.
	eigenvalues_were_computed_ = true;
	eigenvectors_were_computed_ = compute_eigenvectors;

	// Define some scratch-space matrix types.
	ComplexScalar eigenvalue_numerator;
	ScratchVector eigenvalue_numerators_real;
	ScratchVector eigenvalue_numerators_imaginary;
	ScratchVector eigenvalue_denominators;
	ScratchMatrix eigenvectors;

	// Create fortran-compatible wrappers for the matrices defined above.
	fem::dims<1> vector_dimension = dimension(A.rows());
	fem::dims<2> matrix_dimension = dimension(A.rows(), A.cols());

	arr_ref<double, 2> fortran_A(
			*A.data(), matrix_dimension);

	arr_ref<double, 2> fortran_B(
			*B.data(), matrix_dimension);

	arr_ref<double, 1> fortran_eigenvalue_numerators_real(
			*eigenvalue_numerators_real.data(), vector_dimension);

	arr_ref<double, 1> fortran_eigenvalue_numerators_imaginary(
			*eigenvalue_numerators_imaginary.data(), vector_dimension);

	arr_ref<double, 1> fortran_eigenvalue_denominators(
			*eigenvalue_denominators.data(), vector_dimension);

	arr_ref<double, 2> fortran_eigenvectors(
			*eigenvectors.data(), matrix_dimension);

  int error_code = 0;

	// Invoke the EISPACK RGG routine.
	rgg(A.rows(),
			A.cols(),
			fortran_A,
			fortran_B,
			fortran_eigenvalue_numerators_real,
			fortran_eigenvalue_numerators_imaginary,
			fortran_eigenvalue_denominators,
			compute_eigenvectors,
			fortran_eigenvectors,
			error_code);

	// Check for errors.
	if (error_code) {
		info_ = Eigen::NoConvergence;
	}

	// Extract the eigenvalues.
	extract_eigenvalues(
			eigenvalue_numerators_real,
			eigenvalue_numerators_imaginary,
			eigenvalue_denominators);

	// Extract the eigenvectors.
	if (compute_eigenvectors) {
		extract_eigenvectors(eigenvectors);
	}

	// Extract the real solutions.
	extract_real_indices();

  return *this;
}

// GeneralizedEigenSolver::real_eigenvalues {{{1

/// @details The returned eigenvalues will not be in any sorted order.  
///
/// @see compute()
/// @see eigenvalues()

template<typename MatrixType>
typename GeneralizedEigenSolver<MatrixType>::RealEigenvalueType
GeneralizedEigenSolver<MatrixType>::real_eigenvalues() const {

	eigen_assert(eigenvalues_were_computed_ && "No eigenvalues were computed.");

	int num_reals = real_indices_.size();
	RealEigenvalueType real_eigenvalues(num_reals);

	for (int i = 0; i < num_reals; i++) {
		real_eigenvalues[i] = eigenvalues_(real_indices_[i]).real();
	}

	return real_eigenvalues;
}

// GeneralizedEigenSolver::real_eigenvectors {{{1

/// @details Each real eigenvector corresponds to one real eigenvalue, and the 
/// order of eigenvalues returned by real_eigenvalues() will match the order of 
/// eigenvectors returned by this method.  There is no relation between this 
/// method and the order of the complex eigenvalues returned by eigenvalues().
///
/// @see compute()
/// @see eigenvectors()

template<typename MatrixType>
typename GeneralizedEigenSolver<MatrixType>::RealEigenvectorType
GeneralizedEigenSolver<MatrixType>::real_eigenvectors() const {

	eigen_assert(eigenvectors_were_computed_ && 
			"No eigenvectors were computed.");

	int num_reals = real_indices_.size();
	int num_rows = eigenvectors_.rows();

	RealEigenvectorType real_eigenvectors(num_rows, num_reals);

	for (int i = 0; i < num_reals; i++) {
		for (int j = 0; j < num_rows; j++) {
			real_eigenvectors(j, i) = eigenvectors_(j, real_indices_[i]).real();
		}
	}

	return real_eigenvectors;
}

// GeneralizedEigenSolver::extract_eigenvalues {{{1

/// @details RGG returns information about the eigenvalues in three vectors.  
/// The first two vectors are the real and imaginary parts of the numerator, 
/// and the third vector is the (real) denominator.  The complex eigenvalues 
/// can then be calculated as: \f[\lambda = \frac{\mathrm{real\_numerator} + i 
/// \times \mathrm{imaginary\_numerator}}{\mathrm{denominator}}\f]
///
/// @see compute()

template<typename MatrixType>
void GeneralizedEigenSolver<MatrixType>::extract_eigenvalues(
		ScratchVector eigenvalue_numerators_real,
		ScratchVector eigenvalue_numerators_imaginary,
		ScratchVector eigenvalue_denominators) {

	using Eigen::internal::isMuchSmallerThan;

	Index rows = eigenvalue_denominators.rows();
	ComplexScalar eigenvalue_numerator;

	//                  real_numerator + i * imaginary_numerator
	//     eigenvalue = ----------------------------------------
	//                                denominator

	for (Index i = 0; i < rows; i++) {
		eigenvalue_numerator = ComplexScalar(
				eigenvalue_numerators_real(i),
				eigenvalue_numerators_imaginary(i));
		eigenvalues_(i) = eigenvalue_numerator / eigenvalue_denominators(i);
	}
}

// GeneralizedEigenSolver::extract_eigenvectors {{{1

/// The format of the `raw_eigenvectors` matrix produced by RGG is very dense, 
/// because it needs to represent both real and complex eigenvectors.  The role 
/// of this method is to unpack the information stored in this matrix.  The 
/// eigenvectors are represented as column vectors, but the exact meaning of 
/// each column vector depends on the corresponding eigenvalue.  If the j-th 
/// eigenvalue is real, then the j-th eigenvector is simply the j-th column of 
/// the `raw_eigenvectors` matrix.
///
/// On the other hand, if the j-th eigenvalue is complex then the j-th 
/// eigenvector must also be complex.  Complex eigenvectors require two 
/// columns to fully represent, because they have both real and imaginary 
/// components.  Fortunately, complex eigenvalues and eigenvectors always come 
/// in conjugate pairs: \f$\alpha + \beta i\f$ and \f$\alpha - \beta i\f$.  The 
/// columns in `raw_eigenvectors` explicitly give the first member of that 
/// pair, and the second member is implicitly added by the unpacking code.
///
/// @see compute()

template<typename MatrixType>
void GeneralizedEigenSolver<MatrixType>::extract_eigenvectors(
		ScratchMatrix raw_eigenvectors) {

	using Eigen::internal::isMuchSmallerThan;
	Index columns = raw_eigenvectors.cols();

	for (Index j = 0; j < columns; j++) {
		RealScalar real_eigenvalue = eigenvalues_(j).real();
		RealScalar imag_eigenvalue = eigenvalues_(j).imag();

		// This column is a real eigenvector.
		if (isMuchSmallerThan(imag_eigenvalue, real_eigenvalue)) {
			eigenvectors_.col(j) =
				raw_eigenvectors.col(j).template cast<ComplexScalar>();
			eigenvectors_.col(j).normalize();
		}

		// This column is a complex eigenvector.
		else {
			for (Index i = 0; i < columns; i++) {
				RealScalar real_coeff = raw_eigenvectors(i, j);
				RealScalar imag_coeff = raw_eigenvectors(i, j+1);

				eigenvectors_(i, j) = ComplexScalar(real_coeff, imag_coeff);
				eigenvectors_(i, j+1) = ComplexScalar(real_coeff, -imag_coeff);
			}
			eigenvectors_.col(j).normalize();
			eigenvectors_.col(j+1).normalize();
			j++;
		}
	}
}

// GeneralizedEigenSolver::extract_real_indices {{{1

/// @details This method fills `real_indices_` with the indices of all the real 
/// eigenvalues.  These indices can also be used to pick out the real 
/// eigenvectors, because there is a one-to-one correspondence between the two.  
/// The `Eigen::internal::isMuchSmallerThan` function is used to determine 
/// whether or not the imaginary component is negligible relative to the real 
/// component.
///
/// @see compute()
/// @see real_eigenvalues()
/// @see real_eigenvectors()

template<typename MatrixType>
void GeneralizedEigenSolver<MatrixType>::extract_real_indices() {

	using Eigen::internal::isMuchSmallerThan;
	Index indices = eigenvalues_.size();

	for (Index i = 0; i < indices; i++) {
		RealScalar real_eigenvalue = eigenvalues_[i].real();
		RealScalar imag_eigenvalue = eigenvalues_[i].imag();

		if (isMuchSmallerThan(imag_eigenvalue, real_eigenvalue)) {
			real_indices_.push_back(i);
		}
	}
}
// }}}1

} // end namespace linear_algebra
} // end namespace numeric

#endif

