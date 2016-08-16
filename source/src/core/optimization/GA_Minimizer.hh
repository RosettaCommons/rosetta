// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/GA_Minimizer.hh
/// @brief  Minimizer based on Genetic Algorithm
/// @author Sergey Lyskov

#ifndef INCLUDED_core_optimization_GA_Minimizer_hh
#define INCLUDED_core_optimization_GA_Minimizer_hh

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/MinimizerOptions.hh>


#include <core/optimization/Multifunc.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace optimization {


/// Inner class for Genetic Algorithm, hold one population with some additional info
class EItem
{
public:
	EItem() {};
	EItem(const Multivec &vn):
		v( vn ),
		tag( 'X' ),
		r( 0.0 )
	{}

	Multivec v;  ///< item value
	char tag;  ///< tag for debug. (m-mutation, c-crossover, ...etc)
	double r;  ///< result of FF evaluation.

	// for sorting: compare function
	static inline bool sort_R_function(const EItem &e1, const EItem &e2) {return e1.r < e2.r;}
};


class GA_Minimizer
{
public:
	GA_Minimizer(Multifunc & func_in, MinimizerOptions const & options):
		func_( func_in ),
		allowed_time_( 0 ), // Does this make sense?
		add_original_(true),
		mutation_probability_( options.ga_mutation_probability() ),
		minimize_tolerance_( options.minimize_tolerance() )
	{}

	Real run( Multivec & phipsi_inout, int max_time); // starting position, and solution is returned here

private:
	EItem randomize(const EItem& sit, int &time);
	EItem loop(std::vector<EItem> & pop, int &time);
	void step(std::vector<EItem> &pop, int &c_time, int &mres, EItem &shift);


	/// genetic operators
	void mutate(EItem &);  // mutate given vector using normal random.
	void cross_over(EItem &V, EItem &A, EItem &B);

	Multifunc & func_;

	EItem best_;
	int allowed_time_;
	bool add_original_;  // add original best if shifting?

	Real min_error_;
	Real mutation_probability_;
	Real minimize_tolerance_;
};


} // namespace optimization
} // namespace core


#endif // INCLUDED_core_optimization_GA_Minimizer_hh
