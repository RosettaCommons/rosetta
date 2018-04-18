// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/options/StepWiseProteinModelerOptions.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_protein_StepWiseProteinModelerOptions_HH
#define INCLUDED_protocols_stepwise_modeler_protein_StepWiseProteinModelerOptions_HH

#include <protocols/stepwise/modeler/options/StepWiseProteinModelerOptions.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

#if defined(WIN32) || defined(PYROSETTA)
#include <utility/tag/Tag.hh>
#endif


namespace protocols {
namespace stepwise {
namespace modeler {
namespace options {

// multiple inheritance -- bad form -- but will replace with composition later, perhaps.
class StepWiseProteinModelerOptions: public virtual utility::pointer::ReferenceCount {

public:

	//constructor
	StepWiseProteinModelerOptions();

	StepWiseProteinModelerOptions( StepWiseProteinModelerOptions const & src );

	//destructor
	~StepWiseProteinModelerOptions();

public:

	StepWiseProteinModelerOptionsOP clone() const;

	/// @brief Describe this instance to a given output stream
	virtual
	void
	show( std::ostream & ) const{}

	/// @brief Initialize from the recursive "tag" structure.
	virtual
	void
	parse_my_tag( utility::tag::TagCOP ){}

	/// @brief The class name (its type) for a particular ResourceOptions instance.
	/// This function allows for better error message delivery.
	virtual
	std::string
	type() const{ return "StepWiseProteinModelerOptions";}

public:

	void
	initialize_from_command_line();

	void set_global_optimize( bool const & setting ){ global_optimize_ = setting; }
	bool global_optimize() const{ return global_optimize_; }

	void set_sample_beta( bool const & setting ){ sample_beta_ = setting; }
	bool sample_beta() const{ return sample_beta_; }

	void set_move_jumps_between_chains( bool const & setting ){ move_jumps_between_chains_ = setting; }
	bool move_jumps_between_chains() const{ return move_jumps_between_chains_; }

	void set_disable_sampling_of_loop_takeoff( bool const & setting ){ disable_sampling_of_loop_takeoff_ = setting; }
	bool disable_sampling_of_loop_takeoff() const{ return disable_sampling_of_loop_takeoff_; }

	void set_cart_min( bool const & setting ){ cart_min_ = setting; }
	bool cart_min() const{ return cart_min_; }

	void set_n_sample( core::Size const & setting ){ n_sample_ = setting; }
	core::Size n_sample() const{ return n_sample_; }

	void set_filter_native_big_bins( bool const & setting ){ filter_native_big_bins_ = setting; }
	bool filter_native_big_bins() const{ return filter_native_big_bins_; }

	void set_allow_virtual_side_chains( bool const & setting ){ allow_virtual_side_chains_ = setting; }
	bool allow_virtual_side_chains() const{ return allow_virtual_side_chains_; }

	void set_prepack( bool const & setting ){ prepack_ = setting; }
	bool prepack() const{ return prepack_; }

	void set_pack_protein_side_chains( bool const & setting ){ pack_protein_side_chains_ = setting; }
	bool pack_protein_side_chains() const{ return pack_protein_side_chains_; }

	void set_centroid_output( bool const & setting ){ centroid_output_ = setting; }
	bool centroid_output() const{ return centroid_output_; }

	void set_centroid_screen( bool const & setting ){ centroid_screen_ = setting; }
	bool centroid_screen() const{ return centroid_screen_; }

	void set_centroid_score_diff_cut( core::Real const & setting ){ centroid_score_diff_cut_ = setting; }
	core::Real centroid_score_diff_cut() const{ return centroid_score_diff_cut_; }

	void set_centroid_weights( std::string const & setting ){ centroid_weights_ = setting; }
	std::string centroid_weights() const{ return centroid_weights_; }

	void set_nstruct_centroid( core::Size const & setting ){ nstruct_centroid_ = setting; }
	core::Size nstruct_centroid() const{ return nstruct_centroid_; }

	void set_ghost_loops( bool const & setting ){ ghost_loops_ = setting; }
	bool ghost_loops() const{ return ghost_loops_; }

	void set_ccd_close( bool const & setting ){ ccd_close_ = setting; }
	bool ccd_close() const{ return ccd_close_; }

	void set_cluster_by_all_atom_rmsd( bool const & setting ){ cluster_by_all_atom_rmsd_ = setting; }
	bool cluster_by_all_atom_rmsd() const{ return cluster_by_all_atom_rmsd_; }

	void set_pack_weights( std::string const & setting ){ pack_weights_ = setting; }
	std::string pack_weights() const{ return pack_weights_; }

	bool const & expand_loop_takeoff() const { return expand_loop_takeoff_; }
	void set_expand_loop_takeoff( bool const & setting ){ expand_loop_takeoff_ = setting; }

	bool const & skip_coord_constraints() const { return skip_coord_constraints_; }
	void set_skip_coord_constraints( bool const & setting ){ skip_coord_constraints_ = setting; }

	void set_frag_files( utility::vector1< std::string > const & setting ){ frag_files_ = setting; }
	utility::vector1< std::string > frag_files() const { return frag_files_; }

	void set_bridge_res( utility::vector1< core::Size > const & setting ){ bridge_res_ = setting; }
	utility::vector1< core::Size > bridge_res() const { return bridge_res_; }

protected:

	void
	initialize_variables();

private:

	bool global_optimize_;
	bool sample_beta_;
	bool move_jumps_between_chains_;
	bool disable_sampling_of_loop_takeoff_;
	bool cart_min_;
	core::Size n_sample_;

	bool filter_native_big_bins_;
	bool allow_virtual_side_chains_;
	bool prepack_;
	bool pack_protein_side_chains_;

	bool centroid_output_;
	bool centroid_screen_;
	core::Real centroid_score_diff_cut_;
	std::string centroid_weights_;
	core::Size nstruct_centroid_;
	bool ghost_loops_;

	bool ccd_close_; // should this be in here, or in job parameters?
	bool cluster_by_all_atom_rmsd_;

	std::string pack_weights_;

	bool expand_loop_takeoff_;
	bool skip_coord_constraints_;
	utility::vector1< std::string > frag_files_;
	utility::vector1< core::Size > bridge_res_;

};

} //options
} //modeler
} //stepwise
} //protocols

#endif
