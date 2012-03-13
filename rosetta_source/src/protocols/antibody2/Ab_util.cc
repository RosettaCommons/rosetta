// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody2/Ab_util.cc
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)

#include <protocols/antibody2/Ab_util.hh>

// Rosetta Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/id/AtomID_Map.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/import_pose/import_pose.hh>

#include <iostream>
#include <fstream>
#include <numeric/xyz.functions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/OptH.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>








static basic::Tracer TR("antibody2.Ab_util");


using namespace core;
namespace protocols{
namespace antibody2{

    
    //JQX:
    // Description (Jason contributed)
    // assuming a loop 10,11,12,13,14,15, cut_point is 13/14
    //                  or loop = Loop (10, 15, 13)
    // 1. (Rosetta Default) 
    //           set_single_loop_fold_tree()
    //           jump from 8->17 (-2 amd +2)
    // 2. (Aroop's) 
    //            simple_one_loop_fold_tree() 
    //            jump from 9->16 (-1 and +1)
    // 3. (Aroop's) 
    //            simple_fold_tree() 
    //            jump from 10->15, but you need manuly input
    // 4. (Aroop's)
    //            setup_simple_fold_tree()
    //            jump from 10->15, but you need manuly input
    
    void simple_one_loop_fold_tree(
                                   pose::Pose & pose_in,
                                   loops::Loop const & loop	) {
        using namespace kinematics;
        
        TR <<  "H3M Setting up simple one loop fold tree" << std::endl;
        
        //setup fold tree for this loop
        FoldTree f;
        f.clear();
        Size nres = pose_in.total_residue();
        Size jumppoint1 = loop.start() - 1;
        Size jumppoint2 = loop.stop() + 1;
        
        if( jumppoint1 < 1 )   jumppoint1 = 1;
        if( jumppoint2 > nres) jumppoint2 = nres;
        
        f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
        f.add_edge( jumppoint1, loop.cut(), Edge::PEPTIDE );
        f.add_edge( loop.cut() + 1, jumppoint2, Edge::PEPTIDE );
        f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
        f.add_edge( jumppoint1, jumppoint2, 1 );
        f.reorder( 1 );
        
        pose_in.fold_tree( f );
        
        TR <<  "H3M Finished setting up simple one loop fold tree" << std::endl;
        
        return;
    } // simple_one_loop_fold_tree
    
    
    
    
    
    void simple_fold_tree(
                          pose::Pose & pose_in,
                          Size jumppoint1,
                          Size cutpoint,
                          Size jumppoint2 ) {
        using namespace kinematics;
        
        TR <<  "H3M Setting up simple fold tree" << std::endl;
        
        //setup fold tree for this loop
        FoldTree f;
        f.clear();
        Size nres = pose_in.total_residue();
        
        if( jumppoint1 < 1 )   jumppoint1 = 1;
        if( jumppoint2 > nres) jumppoint2 = nres;
        
        f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
        f.add_edge( jumppoint1, cutpoint, Edge::PEPTIDE );
        f.add_edge( cutpoint + 1, jumppoint2, Edge::PEPTIDE );
        f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
        f.add_edge( jumppoint1, jumppoint2, 1 );
        f.reorder( 1 );
        
        pose_in.fold_tree( f );
        
        TR <<  "H3M Finished setting up simple fold tree" << std::endl;
        
        return;
    } // simple_fold_tree
    
    
    
    
    
    
	void setup_simple_fold_tree(
                                          Size jumppoint1,
                                          Size cutpoint,
                                          Size jumppoint2,
                                          Size nres,
                                          pose::Pose & pose_in ) {
        
		using namespace kinematics;
        
		TR << "ABM Setting up simple fold tree" << std::endl;
        
		FoldTree f;
		f.clear();
        
		f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
		f.add_edge( jumppoint1, cutpoint, Edge::PEPTIDE );
		f.add_edge( cutpoint + 1, jumppoint2, Edge::PEPTIDE );
		f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
		f.add_edge( jumppoint1, jumppoint2, 1 );
		f.reorder( 1 );
        
		pose_in.fold_tree( f );
        
		TR << "ABM Done: Setting up simple fold tree" << std::endl;
        
	} // setup_simple_fold_tree
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ///////////////////////////////////////////////////////////////////////////
    /// @begin CDR_H3_filter
    ///
    /// @brief tests if a loop has H3 like base charachteristics
    ///
    /// @detailed Uses the Shirai rules to find out if the dihedral angle
    ///           formed by CA atoms of residues n-2,n-1,n and n+1 conform to a
    ///           kinked/extended structure in accordance with the sequence. If
    ///           there is a match, a true value is returned
    ///
    /// @param[in] pose: full actual protein
    ///            loop_begin: seq numbered loop begin corresponding to pose
    ///            size: size of loop to compute loop_end
    ///
    /// @global_read reads -command line flag -base stored in dle_ns
    ///              to determine to do the complete H3 filter check or just do
    ///              a prediction of the H3 base type based on the
    ///              aforementioned dihedral angle
    ///
    /// @global_write
    ///
    /// @remarks
    ///
    /// @references Structural classification of CDR-H3 in antibodies
    ///             Hiroki Shirai, Akinori Kidera, Haruki Nakamura
    ///             FEBS Letters 399 (1996) 1-8
    ///
    /// @authors Aroop 02/04/2010
    ///
    /// @last_modified 02/04/2010
    ///////////////////////////////////////////////////////////////////////////
    bool CDR_H3_filter(
                                      const pose::Pose & pose_in,
                                      Size const loop_begin,
                                      Size const size,
                                      bool H3_filter,
                                      bool is_camelid)
    {
        
        char const light_chain = 'L';
        
        TR <<  "H3M Checking Kink/Extended CDR H3 Base Angle" << std::endl;
        
        if( !H3_filter || is_camelid )
            return( true );
        
        // Values read from plot in reference paper. Fig 1 on Page 3
        // Values adjusted to match data from antibody training set
        Real const kink_lower_bound = -10.00; // Shirai: 0
        Real const kink_upper_bound = 70.00; // Shirai: 70
        Real const extended_lower_bound = 125.00; // Shirai: ~180
        Real const extended_upper_bound = 185.00; // Shirai: ~180
        
        // Hydrogen Bond maximum value is 3.9 Angstroms - not used
        //	Real const h_bond(3.9);
        // Salt Bridge maximum value is 2.0 Angstroms - not used
        //	Real const s_bridge(4.0);
        
        // chop out the loop
        pose::Pose h3_loop( pose_in, loop_begin - 2, loop_begin + size + 1 );
        
        bool is_kinked( false );
        bool is_extended( false );
        bool is_H3( false );
        
        // extract 3 letter residue codes for the chopped loop
        std::vector <std::string> aa_name; // loop residue 3 letter codes
        for(Size ii = 1; ii <= size + 3; ii++)
            aa_name.push_back( h3_loop.residue(ii).name3() );
        
        Size const CA(2);   // CA atom position in full_coord array
        // base dihedral angle to determine kinked/extended conformation
        Real base_dihedral( numeric::dihedral_degrees(
                                                      h3_loop.residue( aa_name.size() ).xyz( CA ),
                                                      h3_loop.residue( aa_name.size() - 1).xyz( CA ),
                                                      h3_loop.residue( aa_name.size() - 2).xyz( CA ),
                                                      h3_loop.residue( aa_name.size() - 3).xyz( CA ) ) );
        
        // std::cout << "Base Dihedral: " << base_dihedral << std::endl;
        
        // setting up pseudo-periodic range used in extended base computation
        if( base_dihedral < kink_lower_bound )
            base_dihedral = base_dihedral + 360.00;
        
        
        // Rule 1a for standard kink
        if ((aa_name[aa_name.size()-3] != "ASP") &&
            (aa_name[aa_name.size()-1] == "TRP"))	{
            if( (base_dihedral > kink_lower_bound) &&
               (base_dihedral < kink_upper_bound))
            {
                // std::cout << "KINK Found" << std::endl; // aroop_temp remove
                is_kinked = true;
                is_H3 = true;
            }
        }
        
        // Rule 1b for standard extended form
        if ( ( aa_name[ aa_name.size() - 3 ] == "ASP" ) &&
            ( ( aa_name[1] != "LYS" ) && ( aa_name[1] != "ARG" ) ) &&
            ( is_H3 != true ) )     {
            
            if( ( base_dihedral > extended_lower_bound ) &&
               ( base_dihedral < extended_upper_bound) ) {
                // std::cout << "EXTENDED Found" << std::endl; // aroop_temp remove
                is_extended = true;
                is_H3 = true;
            }
            
            if(!is_H3) {
                // Rule 1b extension for special kinked form
                bool is_basic(false); // Special basic residue exception flag
                for(Size ii = 2; ii <= Size(aa_name.size() - 5); ii++) {
                    if( aa_name[ii] == "ARG" || aa_name[ii] == "LYS" ) {
                        is_basic = true;
                        break;
                    }
                }
                
                if(!is_basic) {
                    Size rosetta_number_of_L49 = pose_in.pdb_info()->pdb2pose(light_chain, 49 );
                    std::string let3_code_L49 = pose_in.residue( rosetta_number_of_L49 ).name3();
                    if( let3_code_L49 == "ARG" || let3_code_L49 == "LYS")
                        is_basic = true;
                }
                if( is_basic && ( base_dihedral > kink_lower_bound ) &&
                   ( base_dihedral < kink_upper_bound ) ) {
                    // aroop_temp remove
                    // std::cout << "KINK (special 1b) Found" << std::endl;
                    is_kinked = true;
                    is_H3 = true;
                }
            }
        }
        
        // Rule 1c for kinked form with salt bridge
        if ( ( aa_name[ aa_name.size() - 3 ] == "ASP") &&
            ( ( aa_name[1] == "LYS") || ( aa_name[1] == "ARG" ) ) &&
            ( (aa_name[0] != "LYS" ) && ( aa_name[0] != "ARG" ) ) &&
            ( is_H3 != true) ) {
            if( (base_dihedral > kink_lower_bound ) &&
               (base_dihedral < kink_upper_bound ) ) {
                // aroop_temp remove
                // std::cout << "KINK (w sb) Found" << std::endl;
                is_kinked = true;
                is_H3 = true;
            }
            if(!is_H3) {
                bool is_basic(false); // Special basic residue exception flag
                Size rosetta_number_of_L46 = pose_in.pdb_info()->pdb2pose(
                                                                          light_chain, 46 );
                std::string let3_code_L46 = pose_in.residue( rosetta_number_of_L46 ).name3();
                if( let3_code_L46 == "ARG" || let3_code_L46 == "LYS") is_basic = true;
                if( is_basic && (base_dihedral > extended_lower_bound ) &&
                   ( base_dihedral < extended_upper_bound ) ) {
                    // aroop_temp remove
                    // std::cout << "EXTENDED (special 1c) Found" << std::endl;
                    is_extended = true;
                    is_H3 = true;
                }
            }
        }
        
        // Rule 1d for extened form with salt bridge
        if ( ( aa_name[ aa_name.size() - 3 ] == "ASP") &&
            ( ( aa_name[1] == "LYS") || ( aa_name[1] == "ARG" ) ) &&
            ( ( aa_name[0] == "LYS") || ( aa_name[0] == "ARG") ) &&
            ( is_H3 != true ) ) {
            if( ( base_dihedral > extended_lower_bound ) &&
               ( base_dihedral < extended_upper_bound ) ) {
                // aroop_temp remove
                // std::cout << "EXTENDED (w sb) Found" << std::endl;
                is_extended = true;
                is_H3 = true;
            }
        }
        
        TR <<  "H3M Finished Checking Kink/Extended CDR H3 Base Angle: " << is_H3 << std::endl;
        
        return is_H3;
    } // CDR_H3_filter
    
    
    
    
    
    //TODO:
    //JQX:
    // What you need is to input the "tf" object, 
    // do something to change the value of this "tf" object
    // right now, it is OK, but JQX must come back to make sure the value of tf 
    // can be changed in this function. If not, maybe this function should return a 
    // pointer
    
    void setup_packer_task(pose::Pose & pose_in, core::pack::task::TaskFactoryOP tf ) 
    {
		using namespace pack::task;
		using namespace pack::task::operation;
        
/*		if( init_task_factory_ ) {
			tf = new TaskFactory( *init_task_factory_ );
			TR << "AbModeler Reinitializing Packer Task" << std::endl;
			return;
		}
		else{
			tf = new TaskFactory;
        }
*/
        
        tf = new TaskFactory;
        
        
		TR << "AbModeler Setting Up Packer Task" << std::endl;
        
		tf->push_back( new OperateOnCertainResidues( new PreventRepackingRLT, new ResidueLacksProperty("PROTEIN") ) );
		tf->push_back( new InitializeFromCommandline );
		tf->push_back( new IncludeCurrent );
		tf->push_back( new RestrictToRepacking );
		tf->push_back( new NoRepackDisulfides );
        
		// incorporating Ian's UnboundRotamer operation.
		// note that nothing happens if unboundrot option is inactive!
		pack::rotamer_set::UnboundRotamersOperationOP unboundrot = new pack::rotamer_set::UnboundRotamersOperation();
		unboundrot->initialize_from_command_line();
		operation::AppendRotamerSetOP unboundrot_operation = new operation::AppendRotamerSet( unboundrot );
		tf->push_back( unboundrot_operation );
		// adds scoring bonuses for the "unbound" rotamers, if any
		core::pack::dunbrack::load_unboundrot( pose_in );
        
//		init_task_factory_ = tf;
        
		TR << "AbModeler Done: Setting Up Packer Task" << std::endl;
        
	} // setup_packer_task

    

    
    
    
    
    
    
    


} // namespace antibody2
} // namespace protocols






