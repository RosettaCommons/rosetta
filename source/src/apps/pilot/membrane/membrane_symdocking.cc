// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		apps/pilot/membrane/membrane_symdocking.cc
///
/// @brief		Membrane Framework Application: Symmetric Protein-Protein Docking in Membranes
/// @details	Below is a first-pass at plugging in the membrane framework with symmetry.
///				Uses membrane representation, energ functions, and current tools for Rosetta
///				design for computing ddGs.
///				Last Modified: 12/16/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <devel/init.hh>
#include <protocols/moves/Mover.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/membrane/SetMembranePositionMover.hh> 
#include <protocols/membrane/symmetry/util.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>

#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/symmetric_docking/SymDockProtocol.hh> 
#include <core/pose/symmetry/util.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmData.hh> 
#include <core/conformation/Conformation.hh> 

// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>

#include <core/scoring/ScoreFunction.hh> 
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

static basic::Tracer TR( "apps.pilot.ralford.mp_symdocking" );

/// @brief Mover for Symmetric Docking in Membranes
class MPSymDockMover : public protocols::moves::Mover {

public:
	
	/// @brief Default Constructor
	MPSymDockMover() : Mover() {}
	
	/// @brief Required get_name method
	std::string get_name() const { return "MPSymDockMover"; }
	
	/// @brief Symmetric Docking in the Membrane
	/// @details Setup pose for symmetry, add the membrane components to the total pose,
	/// and then perform symmetric docking
	void apply( Pose & pose ) {

        using namespace core::conformation::symmetry; 
        using namespace core::pose::symmetry;
        using namespace core::scoring;
        using namespace protocols::simple_moves::symmetry; 
        using namespace protocols::membrane; 
        using namespace protocols::membrane::symmetry;
        using namespace protocols::membrane::geometry;
        using namespace protocols::symmetric_docking;
	
        // Check that the pose is neither symmetric nor a membrane protein
        // Protocol is not that advanced yet - build a clean new conf every time
        if ( is_symmetric( pose ) || pose.conformation().is_membrane() ) {
            utility_exit_with_message( "Cannot setup a new symmetric membrane protein if the conformation is already symmetric and/or a membrane pose" );
        }

        // Setup the pose for symmetry based on inputs
        SetupForSymmetryMoverOP setup_for_symm = SetupForSymmetryMoverOP( new SetupForSymmetryMover() );
        setup_for_symm->apply( pose ); 

        // Add a new membrane virtual residue to the symmetry foldtree
        setup_membrane_residue( pose ); 

        // Setup the membrane framework using the specified membrane residue
        // described above
        AddMembraneMoverOP add_memb = AddMembraneMoverOP( new AddMembraneMover( pose.total_residue() ) );
        add_memb->apply( pose ); 

        // Symmetrize the spanning topology
        SpanningTopologyOP symmetrized_topology = symmetrize_spans( pose, *pose.conformation().membrane_info()->spanning_topology() ); 
        symmetrized_topology->show();
        EmbeddingOP embeddings( new Embedding( symmetrized_topology, pose ) );
        EmbeddingDefOP final_embed = embeddings->total_embed(); 

        // Set Initial Membrane position to be the symmetrized position
        SetMembranePositionMoverOP set_initial_position( new SetMembranePositionMover( final_embed->center(), final_embed->normal() ) );
        set_initial_position->apply( pose ); 

        // Configure Symmetric MP Docking. Rosetta will symmetrize the 
        // score functions automatically
        ScoreFunctionOP sfxn_high = ScoreFunctionFactory::create_score_function( "mpframework_symdock_fa_2014" );
        ScoreFunctionOP sfxn_low = ScoreFunctionFactory::create_score_function( "mpframework_symdock_cen_2014" );

        // Setup repulsives based slide criteria (better for membranes than
        // initial contact scoring)
        SymmetricConformation & symm_conf ( dynamic_cast< SymmetricConformation & > ( pose.conformation()) );
        symm_conf.Symmetry_Info()->get_slide_info().set_SlideCriteriaType( FA_REP_SCORE );
    
        // Setup a docking protocol, don't apply filters during docking runs
        // for now
        SymDockProtocolOP symdock( new SymDockProtocol( true, false, false, sfxn_low, sfxn_high ) );
        symdock->hurry( true ); 
        //symdock->apply( pose );  

        // Done!

	}

    /// @brief Add a membrane virtual residue to the pose by inserting by jump
    /// @details Adds a virtual residue to the pose by selecting the VRT_0_Base
    /// as the anchoring residue and nres_complex+1 as tne new sequence position
    /// Not equivalent to an append_by_jump!
    void
    setup_membrane_residue( Pose & pose ) {

        using namespace core::conformation;
        using namespace core::conformation::symmetry;
        using namespace core::chemical;
        using namespace core::kinematics;
    
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
    }
};

// typedefs
typedef utility::pointer::shared_ptr< MPSymDockMover > MPSymDockMoverOP;

/// @brief Main method
int
main( int argc, char * argv [] )
{
	using namespace protocols::jd2;
	
	try {
		
		// Devel init factories
		devel::init(argc, argv);
		
		// Register JD2 options
		protocols::jd2::register_options();
		
		// Setup Membrane Symdocking & go!
		MPSymDockMoverOP symdock( new MPSymDockMover() ); 
		protocols::jd2::JobDistributor::get_instance()->go( symdock );

		return 0;
		
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
