// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/BundleGridSampler.hh
/// @brief  Headers for BundleGridSampler.cc.  Samples conformations of a helical bundle by performing a grid
/// search in Crick parameter space.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_BundleGridSampler_hh
#define INCLUDED_protocols_helical_bundle_BundleGridSampler_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/helical_bundle/BundleGridSampler.fwd.hh>
#include <protocols/helical_bundle/BundleGridSamplerHelper.fwd.hh>
#include <protocols/helical_bundle/BundleGridSamplerHelper.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParametersSet.hh>
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <protocols/helical_bundle/PerturbBundleOptions.fwd.hh>
#include <protocols/helical_bundle/MakeBundle.fwd.hh>
#include <protocols/helical_bundle/MakeBundle.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

//JD2:
#include <protocols/jd2/JobDistributor.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>

// Project Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
	namespace helical_bundle {

		class BundleGridSampler : public protocols::moves::Mover
		{
		public: //Typedefs

				typedef core::conformation::parametric::Parameters Parameters;
				typedef core::conformation::parametric::ParametersOP ParametersOP;
				typedef core::conformation::parametric::ParametersSet ParametersSet;
				typedef core::conformation::parametric::ParametersSetOP ParametersSetOP;

				typedef protocols::helical_bundle::parameters::BundleParameters BundleParameters;
				typedef protocols::helical_bundle::parameters::BundleParametersOP BundleParametersOP;
				typedef protocols::helical_bundle::parameters::BundleParametersCOP BundleParametersCOP;
				typedef protocols::helical_bundle::parameters::BundleParametersSet BundleParametersSet;
				typedef protocols::helical_bundle::parameters::BundleParametersSetOP BundleParametersSetOP;
				typedef protocols::helical_bundle::parameters::BundleParametersSetCOP BundleParametersSetCOP;

		public:
			BundleGridSampler();
			BundleGridSampler( BundleGridSampler const &src );
			virtual ~BundleGridSampler();

			virtual protocols::moves::MoverOP clone() const;
			virtual protocols::moves::MoverOP fresh_instance() const;


			/// @brief Actually apply the mover to the pose.
			virtual void apply(core::pose::Pose & pose);

			virtual std::string get_name() const;

			virtual void parse_my_tag(
				utility::tag::TagCOP tag,
				basic::datacache::DataMap & data,
				protocols::filters::Filters_map const & filters,
				protocols::moves::Movers_map const & movers,
				core::pose::Pose const &
			);

		public:
		////////////////////////////////////////////////////////////////////////////////
		//          PUBLIC FUNCTIONS                                                  //
		////////////////////////////////////////////////////////////////////////////////
		
			/// @brief Set the reset mode.
			/// @brief If true (default), the pose is reset before generating bundles.  If false, it is not.
			void set_reset_mode( bool const val) { reset_mode_=val; return; }
			
			/// @brief Get the reset mode.
			/// @brief If true (default), the pose is reset before generating bundles.  If false, it is not.
			bool reset_mode() const { return reset_mode_; }

			/// @brief Access the r0_ BundleOptions object, by index.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsOP r0( core::Size const index ) {
				runtime_assert_string_msg(index>0 && index<=r0_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::r0() out of range.");
				return r0_[index];
			}

			/// @brief Access the r0_ BundleOptions object, by index.  This provides const access.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsCOP r0( core::Size const index ) const {
				runtime_assert_string_msg(index>0 && index<=r0_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::r0() out of range.");
				return r0_[index];
			}

			/// @brief Access the omega0_ BundleOptions object, by index.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsOP omega0( core::Size const index ) {
				runtime_assert_string_msg(index>0 && index<=omega0_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::omega0() out of range.");
				return omega0_[index];
			}

			/// @brief Access the omega0_ BundleOptions object, by index.  This provides const access.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsCOP omega0( core::Size const index ) const {
				runtime_assert_string_msg(index>0 && index<=omega0_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::omega0() out of range.");
				return omega0_[index];
			}

			/// @brief Access the delta_omega0_ BundleOptions object, by index.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsOP delta_omega0( core::Size const index ) {
				runtime_assert_string_msg(index>0 && index<=delta_omega0_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::delta_omega0() out of range.");
				return delta_omega0_[index];
			}

			/// @brief Access the delta_omega0_ BundleOptions object, by index.  This provides const access.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsCOP delta_omega0( core::Size const index ) const {
				runtime_assert_string_msg(index>0 && index<=delta_omega0_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::delta_omega0() out of range.");
				return delta_omega0_[index];
			}

			/// @brief Access the delta_omega1_ BundleOptions object, by index.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsOP delta_omega1( core::Size const index ) {
				runtime_assert_string_msg(index>0 && index<=delta_omega1_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::delta_omega1() out of range.");
				return delta_omega1_[index];
			}

			/// @brief Access the delta_omega1_ BundleOptions object, by index.  This provides const access.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsCOP delta_omega1( core::Size const index ) const {
				runtime_assert_string_msg(index>0 && index<=delta_omega1_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::delta_omega1() out of range.");
				return delta_omega1_[index];
			}

			/// @brief Access the delta_t_ BundleOptions object, by index.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsOP delta_t( core::Size const index ) {
				runtime_assert_string_msg(index>0 && index<=delta_t_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::delta_t() out of range.");
				return delta_t_[index];
			}

			/// @brief Access the delta_t_ BundleOptions object, by index.  This provides const access.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsCOP delta_t( core::Size const index ) const {
				runtime_assert_string_msg(index>0 && index<=delta_t_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::delta_t() out of range.");
				return delta_t_[index];
			}
			
			/// @brief Access the z1_offset_ BundleOptions object, by index.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsOP z1_offset( core::Size const index ) {
				runtime_assert_string_msg(index>0 && index<=z1_offset_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::z1_offset() out of range.");
				return z1_offset_[index];
			}

			/// @brief Access the z1_offset_ BundleOptions object, by index.  This provides const access.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsCOP z1_offset( core::Size const index ) const {
				runtime_assert_string_msg(index>0 && index<=z1_offset_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::z1_offset() out of range.");
				return z1_offset_[index];
			}
			
			/// @brief Access the z0_offset_ BundleOptions object, by index.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsOP z0_offset( core::Size const index ) {
				runtime_assert_string_msg(index>0 && index<=z0_offset_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::z0_offset() out of range.");
				return z0_offset_[index];
			}

			/// @brief Access the z0_offset_ BundleOptions object, by index.  This provides const access.
			/// @details This is the index in order of helices added, NOT necessarily the index of the helix.
			PerturbBundleOptionsCOP z0_offset( core::Size const index ) const {
				runtime_assert_string_msg(index>0 && index<=z0_offset_.size(), "Index passed to protocols::helical_bundle::PerturbBundle::z0_offset() out of range.");
				return z0_offset_[index];
			}

			/// @brief Access the default_r0_ BundleOptions object.
			///
			PerturbBundleOptionsOP default_r0() { return default_r0_; }

			/// @brief Access the default_r0_ BundleOptions object (const-access).
			///
			PerturbBundleOptionsCOP default_r0() const { return default_r0_; }

			/// @brief Access the default_omega0_ BundleOptions object.
			///
			PerturbBundleOptionsOP default_omega0() { return default_omega0_; }

			/// @brief Access the default_omega0_ BundleOptions object (const-access).
			///
			PerturbBundleOptionsCOP default_omega0() const { return default_omega0_; }

			/// @brief Access the default_delta_omega0_ BundleOptions object.
			///
			PerturbBundleOptionsOP default_delta_omega0() { return default_delta_omega0_; }

			/// @brief Access the default_delta_omega0_ BundleOptions object (const-access).
			///
			PerturbBundleOptionsCOP default_delta_omega0() const { return default_delta_omega0_; }

			/// @brief Access the default_delta_omega1_ BundleOptions object.
			///
			PerturbBundleOptionsOP default_delta_omega1() { return default_delta_omega1_; }

			/// @brief Access the default_delta_omega1_ BundleOptions object (const-access).
			///
			PerturbBundleOptionsCOP default_delta_omega1() const { return default_delta_omega1_; }

			/// @brief Access the default_delta_t_ BundleOptions object.
			///
			PerturbBundleOptionsOP default_delta_t() { return default_delta_t_; }

			/// @brief Access the default_delta_t_ BundleOptions object (const-access).
			///
			PerturbBundleOptionsCOP default_delta_t() const { return default_delta_t_; }

			/// @brief Access the default_z1_offset_ BundleOptions object.
			///
			PerturbBundleOptionsOP default_z1_offset() { return default_z1_offset_; }

			/// @brief Access the default_z1_offset_ BundleOptions object (const-access).
			///
			PerturbBundleOptionsCOP default_z1_offset() const { return default_z1_offset_; }

			/// @brief Access the default_z0_offset_ BundleOptions object.
			///
			PerturbBundleOptionsOP default_z0_offset() { return default_z0_offset_; }

			/// @brief Access the default_z0_offset_ BundleOptions object (const-access).
			///
			PerturbBundleOptionsCOP default_z0_offset() const { return default_z0_offset_; }

			/// @brief Add options for a new helix
			/// @details Return value is the current total number of helices after the addition.
			core::Size add_helix();

			/// @brief Set the maximum number of samples for the mover.
			/// @details If the number of gridpoints based on user options exceeds this number, an error is thrown
			/// and the mover aborts.  This is to prevent unreasonably large calculations from being attempted.
			void set_max_samples( core::Size const val ) { max_samples_=val; return; }

			/// @brief Get the maximum number of samples for the mover.
			/// @details If the number of gridpoints based on user options exceeds this number, an error is thrown
			/// and the mover aborts.  This is to prevent unreasonably large calculations from being attempted.
			core::Size max_samples() const { return max_samples_; }

			/// @brief Increments the number of helices that have been defined.
			///
			void increment_helix_count() { ++n_helices_; return; }

			/// @brief Returns the number of helices that have been defined.
			///
			core::Size n_helices() const { return n_helices_; }

			/// @brief Sets whether the selection should be for the lowest score value (true) or highest (false).
			///
			void set_selection_low( bool const val ) { select_low_=val; return; }

			/// @brief Returns whether the selection should be for the lowest score value (true) or highest (false).
			///
			bool selection_low() { return select_low_; }

			/// @brief Sets the mover that will be applied to all helical bundles generated prior to energy evaluation.
			/// @details Note: if this is used, there is no guarantee that the resulting geometry will still lie within the
			/// parameter space.  (That is, this mover could move the backbone.)
			void set_preselection_mover ( protocols::moves::MoverOP mover );

			/// @brief Sets the filter that will be applied to all helical bundles generated prior to energy evaluation.
			/// @details See the pre_selection_filter_ private member variable for details.
			void set_preselection_filter ( protocols::filters::FilterOP filter );

			/// @brief Returns "true" if and only if a preselection mover has been assigned.
			///
			bool preselection_mover_exists() const { return pre_selection_mover_exists_; }

			/// @brief Returns "true" if and only if a preselection filter has been assigned.
			///
			bool preselection_filter_exists() const { return pre_selection_filter_exists_; }

			/// @brief Set whether the mover dumps pdbs or not.
			///
			void set_pdb_output( bool const val ) { dump_pdbs_=val; return; }

			/// @brief Returns whether the mover dumps pdbs or not.
			///
			bool pdb_output() const { return dump_pdbs_; }

			/// @brief Sets the filename prefix for PDB output.
			/// @details PDB files are of the form <prefix>_#####.pdb.
			void set_pdb_prefix( std::string const &prefix ) { pdb_prefix_=prefix; return; }

			/// @brief Access the filename prefix for PDB output.
			/// @details PDB files are of the form <prefix>_#####.pdb.
			std::string pdb_prefix() { return pdb_prefix_; }

			/// @brief Sets the scorefunction for this mover.
			/// @details This must be done before calling the apply() function.
			void set_sfxn( core::scoring::ScoreFunctionOP sfxn_in ) {
				runtime_assert_string_msg( sfxn_in, "In BundleGridSampler::set_sfxn(): A null scorefunction pointer was provided." );
				sfxn_=sfxn_in;
				sfxn_set_=true;
				return;
			}

			/// @brief Returns whether the scorefunction has been set.
			///
			bool sfxn_set() const { return sfxn_set_; }
			
			/// @brief Set the nstruct mode.
			/// @details If true, each job samples one set of Crick parameters.  If false, every job samples
			/// every set of Crick parameters.  False by default.
			void set_nstruct_mode( bool const &val ) { nstruct_mode_=val; return; }
			
			/// @brief Get the nstruct mode.
			/// @details If true, each job samples one set of Crick parameters.  If false, every job samples
			/// every set of Crick parameters.  False by default.
			bool nstruct_mode( ) const { return nstruct_mode_; }
			
			/// @brief Set the nstruct repeats.
			/// @details This is set to 1 by default, which means that each nstruct number correspnds to a different
			/// set of Crick parameters.  If set greater than 1, then multiple consecutive nstruct numbers will
			/// correspond to the same Crick parameters.  This allows combinatorially combining this mover's sampling
			/// with another, similar mover's sampling.
			void set_nstruct_repeats( core::Size const val ) { nstruct_mode_repeats_=val; return; }
			
			/// @brief Get the nstruct repeats.
			/// @details This is set to 1 by default, which means that each nstruct number correspnds to a different
			/// set of Crick parameters.  If set greater than 1, then multiple consecutive nstruct numbers will
			/// correspond to the same Crick parameters.  This allows combinatorially combining this mover's sampling
			/// with another, similar mover's sampling.
			core::Size nstruct_repeats( ) const { 
				if(nstruct_mode_repeats_ < 1) return 1;
				return nstruct_mode_repeats_;
			}

		private:
		////////////////////////////////////////////////////////////////////////////////
		//          PRIVATE DATA                                                      //
		////////////////////////////////////////////////////////////////////////////////

			/// @brief Should the pose be reset before applying the GridSampler?  Default true.
			///
			bool reset_mode_;
			
			/// @brief Should the parallel sampling be done based on the job (nstruct number)?
			/// @details If true, each job samples one set of Crick parameters.  If false, every job samples
			/// every set of Crick parameters.  False by default.
			bool nstruct_mode_;
			
			/// @brief If nstruct_mode_ is true, how many times should each set of Crick parameters be repeated?
			/// @details This is set to 1 by default, which means that each nstruct number correspnds to a different
			/// set of Crick parameters.  If set greater than 1, then multiple consecutive nstruct numbers will
			/// correspond to the same Crick parameters.  This allows combinatorially combining this mover's sampling
			/// with another, similar mover's sampling.
			core::Size nstruct_mode_repeats_;

			/// @brief The selection type.
			/// @default If false, the pose with the highest score value is selected.  If true,
			/// the pose with the lowest score value is selected.  True by default.
			bool select_low_;

			/// @brief The number of helices that have been defined.
			///
			core::Size n_helices_;

			/// @brief The maximum number of gridpoints allowed.
			/// @details If the number of gridpoints based on user options exceeds this number, an error is thrown
			/// and the mover aborts.  This is to prevent unreasonably large calculations from being attempted.
			/// Default value is ten thousand (10,000).
			core::Size max_samples_;

			/// @brief Default options for sampling r0.
			/// @details May be overridden on a helix-by-helix basis.
			PerturbBundleOptionsOP default_r0_;

			/// @brief Helix-by-helix options for sampling r0.
			///
			PerturbBundleOptionsOPs r0_;

			/// @brief Default options for sampling omega0.
			/// @details May be overridden on a helix-by-helix basis.
			PerturbBundleOptionsOP default_omega0_;

			/// @brief Helix-by-helix options for sampling omega0.
			///
			PerturbBundleOptionsOPs omega0_;

			/// @brief Default options for sampling delta_omega0.
			/// @details May be overridden on a helix-by-helix basis.
			PerturbBundleOptionsOP default_delta_omega0_;

			/// @brief Helix-by-helix options for sampling delta_omega0.
			///
			PerturbBundleOptionsOPs delta_omega0_;

			/// @brief Default options for sampling delta_omega1.
			/// @details May be overridden on a helix-by-helix basis.
			PerturbBundleOptionsOP default_delta_omega1_;

			/// @brief Helix-by-helix options for sampling delta_omega1.
			///
			PerturbBundleOptionsOPs delta_omega1_;

			/// @brief Default options for sampling delta_t.
			/// @details May be overridden on a helix-by-helix basis.
			PerturbBundleOptionsOP default_delta_t_;

			/// @brief Helix-by-helix options for sampling delta_t.
			///
			PerturbBundleOptionsOPs delta_t_;
			
			/// @brief Default options for sampling z1_offset.
			/// @details May be overridden on a helix-by-helix basis.
			PerturbBundleOptionsOP default_z1_offset_;

			/// @brief Helix-by-helix options for sampling z1_offset.
			///
			PerturbBundleOptionsOPs z1_offset_;
			
			/// @brief Default options for sampling z0_offset.
			/// @details May be overridden on a helix-by-helix basis.
			PerturbBundleOptionsOP default_z0_offset_;

			/// @brief Helix-by-helix options for sampling z0_offset.
			///
			PerturbBundleOptionsOPs z0_offset_;

			/// @brief A MakeBundle mover that this mover will call.
			///
			MakeBundleOP make_bundle_;

			/// @brief Owning pointer for an (optional) pre-selection mover applied to all helical bundles before energy evaluation.
			///
			protocols::moves::MoverOP pre_selection_mover_;

			/// @brief Bool determining whether there exists a pre-selection mover that wlil be applied.
			///
			bool pre_selection_mover_exists_;

			/// @brief Owning pointer for an (optional) pre-selection filter applied to all helical bundles after the pre-selection mover but before
			/// picking the lowest-energy solution.  If PDBs are dumped, only those passing filters are dumped.
			protocols::filters::FilterOP pre_selection_filter_;

			/// @brief Bool determining whether a pre-selection filter has been set.
			///
			bool pre_selection_filter_exists_;

			/// @brief Dump a PDB file for each bundle generated?  False by default.
			///
			bool dump_pdbs_;

			/// @brief PDB filename prefix.  Filename will be of the form <prefix>_#####.pdb.
			/// @details  Defaults to "bgs_out".
			std::string pdb_prefix_;

			/// @brief Has the scorefunction been set?
			/// @details False by default.
			bool sfxn_set_;

			/// @brief The scorefunction that this mover will use to pick the lowest-energy bundle.
			/// @details Must be set prior to calling apply() function.
			core::scoring::ScoreFunctionOP sfxn_;
	
		private:
		////////////////////////////////////////////////////////////////////////////////
		//          PRIVATE FUNCTIONS                                                 //
		////////////////////////////////////////////////////////////////////////////////

			/// @brief Is a value in a list?
			///
			bool is_in_list( core::Size const val, utility::vector1 < core::Size> const &list ) const;

			/// @brief Calculate the number of grid points that will be sampled, based on the options set by the user.
			///
			core::Size calculate_total_samples() const;

		}; //BundleGridSampler class

	} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_BundleGridSampler_hh
