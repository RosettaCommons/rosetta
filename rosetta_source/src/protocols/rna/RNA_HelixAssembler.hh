// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file loopRNA_minimizer.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_HelixAssembler_hh
#define INCLUDED_protocols_rna_RNA_HelixAssembler_hh

#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/rna/RNA_IdealCoord.fwd.hh>
#include <core/scoring/rna/RNA_IdealCoord.hh>

#include <core/types.hh>

//// C++ headers
// AUTO-REMOVED #include <cstdlib>
#include <string>

#include <utility/vector1.hh>



namespace protocols {
namespace rna {

/// @brief The RNA de novo structure modeling protocol
class RNA_HelixAssembler: public protocols::moves::Mover {
public:
	/// @brief Construct the protocol object
	RNA_HelixAssembler();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose, std::string const & sequence );

	void random_perturbation( bool const & setting ) { random_perturbation_ = setting; }

	void
	set_minimize_all( bool const & setting ) { minimize_all_ = setting; }

	void
	set_scorefxn( core::scoring::ScoreFunctionOP setting );

	void
	use_phenix_geo( bool const setting );

private:

	void
	set_Aform_torsions( core::pose::Pose & pose, Size const & n );


	void
	build_on_base_pair( core::pose::Pose & pose, Size const & n, char const & seq1, char const & seq2 );

	void
	minimize_base_step( core::pose::Pose & pose, Size const n );

	void
	put_constraints_on_base_step( core::pose::Pose & pose, Size const & n );

private:

	bool verbose_;
	bool random_perturbation_;
	bool minimize_all_;
	bool use_phenix_geo_;
	std::string const ideal_jump;

	core::chemical::ResidueTypeSetCAP rsd_set;

	core::Real const ALPHA_A_FORM;
	core::Real const BETA_A_FORM;
	core::Real const GAMMA_A_FORM;
	core::Real const DELTA_A_FORM;
	core::Real const EPSILON_A_FORM;
	core::Real const ZETA_A_FORM;
	core::Real const CHI_A_FORM;
	core::Real const NU2_A_FORM;
	core::Real const NU1_A_FORM;

	core::Real perturb_amplitude_;

	core::scoring::ScoreFunctionOP scorefxn;

	core::scoring::rna::RNA_IdealCoord ideal_coord_;

}; // class RNA_HelixAssembler



} //rna
} // protocols

#endif
