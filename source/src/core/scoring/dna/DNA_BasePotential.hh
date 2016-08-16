// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_dna_DNA_BasePotential_hh
#define INCLUDED_core_scoring_dna_DNA_BasePotential_hh

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>

#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

namespace core {
namespace scoring {
namespace dna {


class DNA_BasePotential : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~DNA_BasePotential();
	typedef numeric::xyzMatrix< Real > Matrix;
	typedef utility::vector1< Real > Params;
	typedef conformation::Residue Residue;
	typedef ObjexxFCL::FArray3D< Real > FArray3D_Real;
	typedef ObjexxFCL::FArray4D< Real > FArray4D_Real;

public:
	/// ctor
	DNA_BasePotential();


	Real
	base_step_score(
		Residue const & rsd1,
		Residue const & rsd2
	) const;


	Real
	base_pair_score(
		Residue const & rsd1,
		Residue const & rsd2
	) const;


	void
	eval_base_step_derivative(
		Residue const & rsd1,
		Residue const & rsd2,
		Vector & F1,
		Vector & F2,
		Real const external_sign_factor // should probably be +1 or -1
	) const;


	void
	eval_base_pair_derivative(
		Residue const & rsd1,
		Residue const & rsd2,
		Vector & F1,
		Vector & F2,
		Real const sign_factor // need to think about this, see logic in pose_dna
	) const;

	void
	eval_base_pair_Z_scores(
		Residue const & rsd1,
		Residue const & rsd2,
		utility::vector1< Real > & z_scores
	) const;

	void
	eval_base_step_Z_scores(
		Residue const & rsd1,
		Residue const & rsd2,
		utility::vector1< Real > & z_scores
	) const;


private:

	enum InteractionType {
		BP_type = 1,
		BS_type
	};


	inline
	Real
	mean( InteractionType const & t, std::string const & bases, int const p ) const
	{
		int i1(0),i2(0);
		get_array_indices( t, bases, i1, i2 );
		return mean_( p, i1, i2 );
	}

	inline
	Real
	stddev( InteractionType const & t, std::string const & bases, int const p ) const
	{
		int i1(0),i2(0);
		get_array_indices( t, bases, i1, i2 );
		return stddev_( p, i1, i2 );
	}

	inline
	Real
	stiffness( InteractionType const & t, std::string const & bases, int const p1, int const p2 ) const
	{
		int i1(0),i2(0);
		get_array_indices( t, bases, i1, i2 );
		return stiffness_( p1, p2, i1, i2 );
	}

	/// "A","C","T","G"
	inline
	std::string
	base_string( Residue const & rsd ) const;

private:

	/// i1 = 1,2
	/// i2 = 1,16
	void
	get_array_indices( InteractionType const & t, std::string const & bases, int & i1, int & i2 ) const;


	void
	load_score_tables();


	void
	set_mean_and_stddev(
		InteractionType const & type,
		std::string const & bases,
		int const index,
		Real mean,
		Real stddev
	);


	void
	set_stiffness(
		InteractionType const & type,
		std::string const & bases,
		int const index1,
		int const index2,
		Real const val
	);


	Real
	base_score(
		InteractionType const & type,
		std::string const & bases,
		utility::vector1< Real > const & params
	) const;

private:

	FArray3D_Real      mean_; //(    6, 2, 16 );
	FArray3D_Real    stddev_; //(    6, 2, 16 );
	FArray4D_Real stiffness_; //( 6, 6, 2, 16 );

};


} // namespace dna
} // scoring
} // core

#endif
