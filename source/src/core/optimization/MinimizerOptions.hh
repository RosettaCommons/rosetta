// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/MinimizerOptions.hh
/// @brief  Minimizer options class
/// @author Phil Bradley


#ifndef INCLUDED_core_optimization_MinimizerOptions_hh
#define INCLUDED_core_optimization_MinimizerOptions_hh

#include <core/optimization/MinimizerOptions.fwd.hh>

#include <core/types.hh>

// C++ headers
#include <string>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace optimization {

class MinimizerOptions : public utility::pointer::ReferenceCount {

public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~MinimizerOptions();
	/////////////////////////////////////////////////////////////////////////////
	// c-tor's -- might want to add a number of c-tors or a string-switched
	// c-tor to replicate the dependence of, eg, xx_init on the current
	// func_switch setting

	MinimizerOptions(
		std::string const & min_type_in,
		Real const minimize_tolerance_in,
		bool const use_nblist_in, // see core/scoring/NeighborList.hh
		bool const deriv_check_in = false,
		bool const deriv_check_verbose_in = false
	);

	///
	MinimizerOptionsOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// high-level params

	// the min-type, eg "dfpmin", "linmin"
	std::string const &
	min_type() const;

	void
	min_type( std::string min_type_in );

	std::string &
	min_type();

	void
	deriv_check( bool deriv_check_in );

	void
	deriv_check_to_stdout( bool setting );

	bool
	deriv_check() const;

	bool
	deriv_check_verbose() const;

	bool
	deriv_check_to_stdout() const;

	// the tolerance for dfpmin, dfpmin_atol etc
	Real
	minimize_tolerance() const;

	Real &
	minimize_tolerance();

	void
	minimize_tolerance( Real minimize_tolerance_in );

	/// @brief Indicate whether or not a handful of optimizations regarding the 
	/// neighbor list have been enabled.
	/// @copydetails use_nblist()
	bool
	use_nblist() const;

	/// @brief Indicate whether or not a handful of optimizations regarding the 
	/// neighbor list should be enabled.
	///
	/// @details The neighbor list is based on the move map.  It includes any 
	/// atoms that can be moved by the minimizer plus their neighbors.  This list 
	/// is not updated during minimization.  All scores for atoms and atom pairs 
	/// outside the neighbor list are held fixed.  All scores for atoms and atom 
	/// pairs within the list are not cached (usually they would be), since it 
	/// assumed that they will be changing rapidly.  These optimizations are 
	/// effective when a large number of small moves are being made.  You may 
	/// prefer to enable this option when minimizing in fullatom mode, and to 
	/// disable it in centroid mode.
	///
	/// @note I wrote this docstring after I reading through a bunch of code to 
	/// figure out what this option meant, but I don't have any expertise with 
	/// the minimizer.  So the above represents my understanding of this option, 
	/// but there could very well be errors or misstatements.
	///
	/// @see core::scoring::AtomNeighbor
	void
	use_nblist( bool use_nblist_in );


	bool
	nblist_auto_update() const;

	void
	nblist_auto_update( bool setting );

	bool
	silent() const;

	void
	silent( bool silent_in );

	Real
	gmax_cutoff_for_convergence() const;

	void
	gmax_cutoff_for_convergence( Real gmax_in );

	/////////////////////////////////////////////////////////////////////////////
	// low-level params

	// bracketing params for linmin
	Real
	ax_init() const;

	Real
	xx_init() const;

	Real
	bx_init() const;


	// abs tolerance for brent
	Real
	brent_abs_tolerance() const;

	int max_iter() const;
	void max_iter(int n);

	Real ga_mutation_probability() const;
	void ga_mutation_probability(Real p);


	///////
	// data

private:
	int max_iter_;

	std::string min_type_;

	Real minimize_tolerance_;

	bool use_nblist_;
	bool nblist_auto_update_;

	bool deriv_check_;
	bool deriv_check_verbose_;
	bool deriv_check_to_stdout_;

	bool silent_;  // no minimizer output

	Real gmax_cutoff_for_convergence_;  // armijo line min: only converge if max grad < this value 
	                                    //    defaults to 1.0

	Real ax_init_;
	Real xx_init_;
	Real bx_init_;

	Real brent_abs_tolerance_;

	Real ga_mutation_probability_;

}; // MinimizerOptions


} // namespace optimization
} // namespace core


#endif // INCLUDED_core_optimization_MinimizerOptions_HH
