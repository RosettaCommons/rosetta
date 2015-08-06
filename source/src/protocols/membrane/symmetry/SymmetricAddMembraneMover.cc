// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/symmetry/SymmetricAddMembraneMover.cc
///
/// @brief      Add Membrane Representation to a Symmetric starting pose
/// @details    Given a symmetrized pose, add the membrane components,
///             ensuring all data descriptions are for the asymmetric subunit
///             (in terms of scoring) and the membrane is the root of the
///             whole system. After applying SymmetricAddMembraneMover
///             pose.conformaiton().is_membrane() AND is_symmetric( pose )
///             should both return true
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (2/9/15)

#ifndef INCLUDED_protocols_membrane_symmetry_SymmetricAddMembraneMover_cc
#define INCLUDED_protocols_membrane_symmetry_SymmetricAddMembraneMover_cc

// Unit Headers
#include <protocols/membrane/symmetry/SymmetricAddMembraneMover.hh>
#include <protocols/membrane/symmetry/SymmetricAddMembraneMoverCreator.hh>

#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/util.hh>

#include <core/pose/PDBInfo.hh>
#include <core/conformation/membrane/MembraneParams.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/symmetry/util.hh>

// Package Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static thread_local basic::Tracer TR( "protocols.membrane.symmetry.SymmetricAddMembraneMover" );

namespace protocols {
namespace membrane {
namespace symmetry {
    
using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;

/////////////////////
/// Constructors  ///
/////////////////////
    
    //increment

/// @brief Default constructor of SymmetricAddMembraneMover
/// @details Create a membrane pose, setting the membrane center
/// at center=(0, 0, 0), normal=(0, 0, 1) and loads in spans
/// and lips from the command line interface. Calls the
/// default constructor of AddMembraneMover
SymmetricAddMembraneMover::SymmetricAddMembraneMover() :
    protocols::membrane::AddMembraneMover()
{
    register_options();
    init_from_cmd();
}
    
    
/// @brief Custom Constructor for SymmetricAddMembraneMover
/// @details Creates a membrane pose, setting the membrane center
/// at center=(0, 0, 0), normal=(0, 0, 1). Uses the user provided
/// spanfile in this constructor to load the spanning topology.
/// Calls the analagous parent constructor in AddMembraneMover
SymmetricAddMembraneMover::SymmetricAddMembraneMover(
   std::string spanfile
   ) :
    protocols::membrane::AddMembraneMover( spanfile )
{
    register_options();
    init_from_cmd();
}
    
/// @brief Custom Constructor for SymmetricAddMembraneMover
/// @details Creates a membrane pose, setting the membrane center
/// at center=(0, 0, 0), normal=(0, 0, 1). Uses the user provided
/// SpanningTopology to be set in MembraneInfo. Calls the analagous
/// parent constructor in AddMembraneMover
SymmetricAddMembraneMover::SymmetricAddMembraneMover(
   SpanningTopologyOP topology
   ) :
   protocols::membrane::AddMembraneMover( topology )
{
    register_options();
    init_from_cmd();
}
    
/// @brief Custom Constructor for SymmetricAddMembraneMover
/// @details Creates a membrane pose, setting the membrane center
/// at center=(0, 0, 0), normal=(0, 0, 1). Uses the user provided
/// spanfile and lipsfile paths to load in SpanningTopology
/// objects and LipidAccInfo objects to MembraneInfo. Calls the analagous
/// parent constructor in AddMembraneMover
SymmetricAddMembraneMover::SymmetricAddMembraneMover(
    std::string spanfile,
    std::string lips_acc
    ) :
    protocols::membrane::AddMembraneMover( spanfile, lips_acc, 0 )
{
    register_options();
    init_from_cmd();
}
    
/// @brief Copy Constructor for SymmetricAddMembraneMover
/// @details Create a deep copy of SymmetricAddMembraneMover
SymmetricAddMembraneMover::SymmetricAddMembraneMover( SymmetricAddMembraneMover const & src ) :
    protocols::membrane::AddMembraneMover( src )
{}

/// @brief Destructor
SymmetricAddMembraneMover::~SymmetricAddMembraneMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
SymmetricAddMembraneMover::clone() const {
    return ( protocols::moves::MoverOP( new SymmetricAddMembraneMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
SymmetricAddMembraneMover::fresh_instance() const {
    return protocols::moves::MoverOP( new SymmetricAddMembraneMover );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
SymmetricAddMembraneMover::parse_my_tag(
    utility::tag::TagCOP,
    basic::datacache::DataMap &,
    protocols::filters::Filters_map const &,
    protocols::moves::Movers_map const &,
    core::pose::Pose const &
    ) {
}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
SymmetricAddMembraneMoverCreator::create_mover() const {
    return protocols::moves::MoverOP( new SymmetricAddMembraneMover );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
SymmetricAddMembraneMoverCreator::keyname() const {
    return SymmetricAddMembraneMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
SymmetricAddMembraneMoverCreator::mover_name() {
    return "SymmetricAddMembraneMover";
}
    

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (SymmetricAddMembraneMover)
std::string
SymmetricAddMembraneMover::get_name() const {
    return "SymmetricAddMembraneMover";
}
    
    
/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
/// @details Register mover-relevant options with JD2 - includes
/// mp seutp options: spanfiles, lipsfile, and
/// spans_from_structure.
void
SymmetricAddMembraneMover::register_options() {
    
    using namespace basic::options;
    
    option.add_relevant( OptionKeys::mp::setup::spanfiles );
    option.add_relevant( OptionKeys::mp::setup::spans_from_structure );
    option.add_relevant( OptionKeys::mp::setup::lipsfile );
    
}
    
/// @brief Initialize Mover options from the comandline
/// @details Initialize mover settings from the commandline
/// mainly in the mp, setup group: spans_from_structure,
/// spanfile and lipsfiles paths
void
SymmetricAddMembraneMover::init_from_cmd() {
    
    using namespace basic::options;
    
    // Read in User-Provided spanfile
    if ( get_spanfile().size() == 0 && option[ OptionKeys::mp::setup::spanfiles ].user() ) {
        spanfile( option[ OptionKeys::mp::setup::spanfiles ]()[1] );
    }
    
    if ( get_spanfile().size() == 0 && option[ OptionKeys::mp::setup::spans_from_structure ].user() ) {
        TR << "WARNING: Spanfile not given, topology will be created from PDB!" << std::endl;
        TR << "WARNING: Make sure your PDB is transformed into membrane coordinates!!!" << std::endl;
        spanfile( "from_structure" );
    }
    
    // Read in User-provided lipsfiles
    if ( option[ OptionKeys::mp::setup::lipsfile ].user() ) {
        
        // Set include lips to true and read in filename
        include_lips( true );
        lipsfile( option[ OptionKeys::mp::setup::lipsfile ]() );
    }
    
    // Read in user-provided membrane residue position
    if ( option[ OptionKeys::mp::setup::membrane_rsd ].user() ) {
        utility_exit_with_message( "Cannot specify a custom membrane position with symmetry code. Must setup a conformation for membranes and symmetry using SymmetricAddMembraneMover (feature coming soon)" );
    }
}

    
    
/// @brief Add a membrane virtual residue to the pose by inserting by jump
/// @details Adds a virtual residue to the pose by selecting the VRT_0_Base
/// as the anchoring residue and nres_complex+1 as tne new sequence position
/// Not equivalent to an append_by_jump!
core::Size
SymmetricAddMembraneMover::setup_membrane_virtual( Pose & pose ) {
    
    using namespace core::conformation;
    using namespace core::conformation::symmetry;
    using namespace core::chemical;
    using namespace core::kinematics;
    
    // Send a giant warning if the pose is Asymmetric!
    if ( !core::pose::symmetry::is_symmetric( pose ) ) {
        utility_exit_with_message( "Cannot create a symmetric membrane pose from an asymmetric pose. Please run SetupForSymetryMover (in protocols/simple_moves) first!" );
    }
    
    TR << "Adding a membrane residue representing the position of the membrane after residue " << pose.total_residue() << std::endl;
    
    // Get the current residue typeset of the pose and use it to determine the
    // typeset of the new membrane residue
    ResidueTypeSetCOP const & residue_set(
          ChemicalManager::get_instance()->residue_type_set( pose.is_fullatom() ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
          );
    
    // Create a new Residue from rsd typeset of type MEM
    ResidueTypeCOPs const & rsd_type_list( residue_set->name3_map("MEM") );
    ResidueType const & membrane( *rsd_type_list[1] );
    ResidueOP new_rsd( ResidueFactory::create_residue( membrane ) );
    
    // Obtain information about the symmetric setup of virtual residues
    // in the pose from the symmetry info object
    SymmetricConformation & symm_conf ( dynamic_cast< SymmetricConformation & > ( pose.conformation()) );
    
    // Grab the number of subunits and number of residues in the monomeric unit
    core::Size nsubunits( symm_conf.Symmetry_Info()->subunits() );
    core::Size nres_monomer( symm_conf.Symmetry_Info()->get_nres_subunit() );
    core::Size total_res(  pose.total_residue() ); // wants total res of the whole compelex
    core::Size num_virtuals( symm_conf.Symmetry_Info()->num_virtuals() );
    
    // Compute the target sequence position and anchoring position for the residue
    core::Size seqpos( ++total_res );
    core::Size anchor( (nres_monomer * nsubunits) + 1 ); // anchor to VRT_0 base (current root)
    
    // Insert the membrane residue by jump to the anchoring position
    TR << "Adding a membrane residue at " << seqpos << " anchored by " << anchor << std::endl;
    pose.conformation().insert_residue_by_jump( *new_rsd, seqpos, anchor, "", "", false ); // never start a new chain, use default anchoring atoms
    
    // Update number of virtuals in the symmetric pose to include the
    // membrane residue (don't change this line, ever. If MEM isnot included in
    // the virtual count, because it adds asymmetry but still accounted for in bb_independent
    // Rosetta expects a slave MEM which doesn't exist and the entire pose gets disconnected,
    // quite literally)
    symm_conf.Symmetry_Info()->num_virtuals( ++num_virtuals );
    
    // Update the foldtree such that the membrane residue is the root of the
    // symmetrized foldtree
    FoldTree ft( pose.fold_tree() );
    ft.reorder( seqpos ); // seqpos == mprsd posiiton
    pose.fold_tree( ft );
    
    // Check that the pose is still symmetric and the foldtree is still valid
    // These are very rigorous checks, leaving them here until this is better tested
    assert( core::pose::symmetry::is_symmetric( pose ) );

    return pose.total_residue();
}
    
/// @brief Helper Method - Check for Membrane residue already in the PDB
/// @details If there is an MEM residue in the PDB at the end of the pose
/// with property MEMBRANE, return a vector of all of those residues.
utility::vector1< core::SSize >
SymmetricAddMembraneMover::check_pdb_for_mem( Pose & ) {

    // initialize vector for membrane residues found in PDB
    utility::vector1< core::Size > mem_rsd;
    mem_rsd.push_back( -1 ); 
    return mem_rsd;
}

} // symmetry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_SymmetricAddMembraneMover_cc

