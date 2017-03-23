// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file loopRNA_minimizer.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_HelixAssembler_hh
#define INCLUDED_protocols_rna_RNA_HelixAssembler_hh

#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/types.hh>

//// C++ headers
#include <string>


namespace protocols {
namespace rna {
namespace movers {

// helper function -- which characters correspond to gaps?
bool
is_blank_seq( char const & c );

/// @brief The RNA de novo structure modeling protocol
class RNA_HelixAssembler: public protocols::moves::Mover {
public:

	/// @brief Construct the protocol object
	RNA_HelixAssembler();

	/// @brief Destroy the protocol object
	~RNA_HelixAssembler();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Apply the loop-rebuild protocol to the input pose
	using protocols::moves::Mover::apply;
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose, std::string const & full_sequence_in );

	void random_perturbation( bool const & setting ) { random_perturbation_ = setting; }

	void
	set_dump( bool const & setting ) { dump_ = setting; }

	void
	set_minimize_all( bool const & setting ) { minimize_all_ = setting; }

	void
	set_scorefxn( core::scoring::ScoreFunctionOP setting );

	void
	set_finish_scorefxn( core::scoring::ScoreFunctionOP setting );

	void
	use_phenix_geo( bool const setting );

	void set_model_and_remove_capping_residues( bool setting ){ model_and_remove_capping_residues_ = setting; }

	void
	build_helix( core::pose::Pose & pose );

	core::pose::PoseOP
	build_init_pose( std::string const & seq1, std::string const & seq2 );

private:

	void
	initialize_minimizer();

	void
	set_Aform_torsions( core::pose::Pose & pose, Size const & n ) const;

	void
	build_on_base_pair( core::pose::Pose & pose, Size const & n, std::string const & seq1, std::string const & seq2 ) const;

	void
	minimize_base_step( core::pose::Pose & pose, Size const n, core::scoring::ScoreFunctionOP scorefxn ) const;

	void
	put_constraints_on_base_step( core::pose::Pose & pose, Size const & n ) const;

	void
	add_capping_base_pairs_to_full_sequence();

	void
	get_rid_of_capping_base_pairs( core::pose::Pose & pose );

	void
	figure_out_and_remove_dangling_ends();

	void
	build_dangling_ends( core::pose::Pose & pose ) const;

	void
	build_dangle_seq1_5prime( core::pose::Pose & pose, std::string const & dangle_seq ) const;

	void
	build_dangle_seq2_5prime( core::pose::Pose & pose, std::string const & dangle_seq ) const;

	void
	build_dangle_seq1_3prime( core::pose::Pose & pose, std::string const & dangle_seq ) const;

	void
	build_dangle_seq2_3prime( core::pose::Pose & pose, std::string const & dangle_seq ) const;

	core::Size
	get_cutpoint( core::pose::Pose const & pose ) const;

	void
	append_Aform_residue( core::pose::Pose & pose, Size const & n, std::string const & nt ) const;

	void
	prepend_Aform_residue( core::pose::Pose & pose, Size const & n, std::string const & nt ) const;

	void
	perturb_torsion(
		core::pose::Pose & pose,
		utility::vector1<core::id::TorsionID> const & id_list
	) const;

	void
	minimize_append_res( core::pose::Pose & pose, Size const n ) const;

	void
	minimize_prepend_res( core::pose::Pose & pose, Size const n ) const;

	void
	fill_chain_info( core::pose::Pose & pose );

	std::string
	get_sequence( core::Size const n );

	void
	remove_first_base_pair( std::string & full_sequence,
		std::map< Size, std::string > & non_standard_residues,
		std::string & sequence_helix1,
		std::string & sequence_helix2 ) const;

	void
	remove_last_base_pair( std::string & full_sequence,
		std::map< Size, std::string > & non_standard_residues,
		std::string & sequence_helix1,
		std::string & sequence_helix2 ) const;

	core::conformation::ResidueOP
	get_residue( std::string const & nt ) const;

private:

	bool dump_;
	bool random_perturbation_;
	bool minimize_all_;
	bool minimize_jump_;
	bool use_phenix_geo_;

	core::chemical::ResidueTypeSetCOP rsd_set_;

	core::chemical::rna::RNA_FittedTorsionInfo const torsion_info_;

	core::Real perturb_amplitude_;

	core::scoring::ScoreFunctionOP scorefxn_;

	core::scoring::ScoreFunctionOP finish_scorefxn_;

	bool model_and_remove_capping_residues_;
	std::string const capping_residues_;

	core::optimization::AtomTreeMinimizerOP minimizer_;
	core::optimization::MinimizerOptionsOP minimizer_options_;

	std::string dangle_seq1_5prime_;
	std::string dangle_seq1_3prime_;
	std::string dangle_seq2_5prime_;
	std::string dangle_seq2_3prime_;

	std::string full_sequence_;
	std::map< core::Size, std::string > non_standard_residues_;

}; // class RNA_HelixAssembler


} //movers
} //rna
} //protocols

#endif
