// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/inverse_folding/MIFST.hh
/// @brief A class for using the MIF-ST model developed by Yang et al.
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

#ifdef USE_TORCH

#ifndef INCLUDED_protocols_inverse_folding_MIFST_hh
#define INCLUDED_protocols_inverse_folding_MIFST_hh

#include <protocols/inverse_folding/MIFST.fwd.hh>

// basic headers
#include <basic/citation_manager/CitationCollectionBase.fwd.hh>

// Core headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

// Torch header
#ifdef Assert
	#undef Assert
#endif
#include <torch/script.h>

// Eigen header
#include <Eigen/Dense>

// std headers
#include <cmath>

namespace protocols {
namespace inverse_folding {

/// @brief A class for using the MIF-ST model developed Yang et al.
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)
class MIFST :  public utility::SingletonBase< MIFST >{
    friend class utility::SingletonBase< MIFST >;

public:
    /// @brief Private constructor for singleton class
    MIFST();
    MIFST( MIFST const & ) = delete;
    MIFST operator=( MIFST const & ) = delete;

    /// @brief Get the citation for MIF-ST
    static
    basic::citation_manager::CitationCollectionBaseCOP
    get_MIFST_neural_net_citation();

    ///@brief Predict amino acid probabilities for a given set of residues and reformat them into a probabilities metric
    std::map<core::Size, std::map<core::chemical::AA, core::Real>>
    sample( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const &residue_selection,
             core::select::residue_selector::ResidueSubset const &feature_selection, bool multirun, bool use_gpu);

    ///@brief Predict amino acid probabilities for a given set of residues
    torch::Tensor
    predict( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const &residue_selection,
            core::select::residue_selector::ResidueSubset const &feature_selection, bool multirun, bool use_gpu);

    /// @brief set default options
    void set_defaults();
    /// @brief get auto_download value
    bool auto_download() const { return auto_download_; }
    /// @brief set the auto_download value
    void auto_download( bool setting );


private: //methods
    /// @brief Initialize the PyTorch model used by this class.
    /// @details Called ONCE by class constructor.  Triggers read from disk!
    void init_mifst_model();

    /// @brief Downloads model from GitLab if the specified path does not exist or is empty
    /// @param[in] path_to_model Directory path where the model should be located
    /// @param[in] auto_download Whether to automatically download missing models
    static void download_model_if_not_existing( std::string const & path_to_model, bool const & auto_download );

    static std::string const GITLAB_URL_;

private: //data
    torch::jit::script::Module mifst_module_;
    // defines whether models get automatically downloaded, default value of option is false
    bool auto_download_ = false;
};

} //inverse_folding
} //protocols

#endif //INCLUDED_protocols_inverse_folding_MIFST_hh

#endif //USE_TORCH
