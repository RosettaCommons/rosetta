// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody2/AntibodyUtil.cc
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)

#include <protocols/antibody2/AntibodyUtil.hh>

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
#include <protocols/rigid/RB_geometry.hh>
#include <core/scoring/rms_util.tmpl.hh>








static basic::Tracer TR("antibody2.AntibodyUtil");


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
        
        TR <<  "Utility: Setting up simple one loop fold tree" << std::endl;
        
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
        
        TR <<  "Utility: Finished setting up simple one loop fold tree" << std::endl;
        
        return;
    } // simple_one_loop_fold_tree
    
    
    
    
    
    void simple_fold_tree(
                          pose::Pose & pose_in,
                          Size jumppoint1,
                          Size cutpoint,
                          Size jumppoint2 ) {
        using namespace kinematics;
        
        TR <<  "Utility: Setting up simple fold tree" << std::endl;
        
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
        
        TR <<  "Utility: Finished setting up simple fold tree" << std::endl;
        
        return;
    } // simple_fold_tree
    
    
    
    
    
    
	void setup_simple_fold_tree(
                                          Size jumppoint1,
                                          Size cutpoint,
                                          Size jumppoint2,
                                          Size nres,
                                          pose::Pose & pose_in ) {
        
		using namespace kinematics;
        
		TR << "Utility: Setting up simple fold tree" << std::endl;
        
		FoldTree f;
		f.clear();
        
		f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
		f.add_edge( jumppoint1, cutpoint, Edge::PEPTIDE );
		f.add_edge( cutpoint + 1, jumppoint2, Edge::PEPTIDE );
		f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
		f.add_edge( jumppoint1, jumppoint2, 1 );
		f.reorder( 1 );
        
		pose_in.fold_tree( f );
        
		TR << "Utility: Done: Setting up simple fold tree" << std::endl;
        
	} // setup_simple_fold_tree
    
    
    
    
    
    ///////////////////////////////////////////////////////////////////////////
	/// @begin all_cdr_VL_VH_fold_tree
	///
	/// @brief change to all CDR and VL-VH dock fold tree
	///
	/// @authors Aroop 07/13/2010
	///
	/// @last_modified 07/13/2010
	///////////////////////////////////////////////////////////////////////////
	void all_cdr_VL_VH_fold_tree( pose::Pose & pose_in, const loops::Loops & loops_in ) 
    {
        
		using namespace kinematics;
        
		Size nres = pose_in.total_residue();
		core::pose::PDBInfoCOP pdb_info = pose_in.pdb_info();
		char second_chain = 'H';
		Size rb_cutpoint(0);
        
		for ( Size i = 1; i <= nres; ++i ) {
			if( pdb_info->chain( i ) == second_chain) {
				rb_cutpoint = i-1;
				break;
			}
		}
        
		Size jump_pos1 ( geometry::residue_center_of_mass( pose_in, 1, rb_cutpoint ) );
		Size jump_pos2 ( geometry::residue_center_of_mass( pose_in,rb_cutpoint+1, nres ) );
        
		// make sure rb jumps do not reside in the loop region
		for( loops::Loops::const_iterator it= loops_in.begin(), it_end = loops_in.end(); it != it_end; ++it ) {
			if ( jump_pos1 >= ( it->start() - 1 ) &&
                jump_pos1 <= ( it->stop() + 1) )
				jump_pos1 = it->stop() + 2;
			if ( jump_pos2 >= ( it->start() - 1 ) &&
                jump_pos2 <= ( it->stop() + 1) )
				jump_pos2 = it->start() - 2;
		}
        
		// make a simple rigid-body jump first
		setup_simple_fold_tree(jump_pos1,rb_cutpoint,jump_pos2,nres, pose_in );
        
		// add the loop jump into the current tree,
		// delete some old edge accordingly
		FoldTree f( pose_in.fold_tree() );
        
		for( loops::Loops::const_iterator it=loops_in.begin(),
            it_end=loops_in.end(); it != it_end; ++it ) {
			Size const loop_start ( it->start() );
			Size const loop_stop ( it->stop() );
			Size const loop_cutpoint ( it->cut() );
			Size edge_start(0), edge_stop(0);
			bool edge_found = false;
			const FoldTree & f_const = f;
			Size const num_jump = f_const.num_jump();
			for( FoldTree::const_iterator it2=f_const.begin(),
                it2_end=f_const.end(); it2 !=it2_end; ++it2 ) {
				edge_start = std::min( it2->start(), it2->stop() );
				edge_stop = std::max( it2->start(), it2->stop() );
				if ( ! it2->is_jump() && loop_start > edge_start
                    && loop_stop < edge_stop ) {
					edge_found = true;
					break;
				}
			}
            
			f.delete_unordered_edge( edge_start, edge_stop, Edge::PEPTIDE);
			f.add_edge( loop_start-1, loop_stop+1, num_jump+1 );
			f.add_edge( edge_start, loop_start-1, Edge::PEPTIDE );
			f.add_edge( loop_start-1, loop_cutpoint, Edge::PEPTIDE );
			f.add_edge( loop_cutpoint+1, loop_stop+1, Edge::PEPTIDE );
			f.add_edge( loop_stop+1, edge_stop, Edge::PEPTIDE );
		}
        
		f.reorder(1);
		pose_in.fold_tree(f);
	} // all_cdr_VL_VH_fold_tree
    
    

    
    
    
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
    /// @references Structural classification of CDR-H3 in antibodies
    ///             Hiroki Shirai, Akinori Kidera, Haruki Nakamura
    ///             FEBS Letters 399 (1996) 1-8
    ///
    /// @authors Aroop 02/04/2010
    ///
    /// @last_modified 02/04/2010
    ///////////////////////////////////////////////////////////////////////////
    
    //TODO:
    //JQX:
    //work with Daisuke to put the L89 creteria into the code
    
    bool CDR_H3_filter(const pose::Pose & pose_in, loops::Loop & input_loop, bool is_camelid)
    {

        
        TR <<  "Utility: Checking Kink/Extended CDR H3 Base Angle" << std::endl;
        
        
        char const light_chain = 'L';
        
        if(is_camelid )
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
        
        // chop out the loop: 
        //JQX: 2 residues before h3, one residue after h3. Matched Rosetta2!
        Size start(input_loop.start()-2);
        Size stop(input_loop.stop()+1);
        
        
        bool is_kinked( false );
        bool is_extended( false );
        bool is_H3( false );
        
        // extract 3 letter residue codes for the chopped loop
        std::vector <std::string> aa_name; // loop residue 3 letter codes
        //JQX: pay attention here!! It is vector, not vector1! too painful to compare to R2 code
        //     just make vector, so it can match R2 code easily
        for(Size ii=start; ii<=stop;ii++){
            aa_name.push_back(pose_in.residue(ii).name3() );
//            TR<<pose_in.residue(ii).name1()<<std::endl;
        }
        
        Size const CA(2);   // CA atom position in full_coord array
        // base dihedral angle to determine kinked/extended conformation

        Real base_dihedral( numeric::dihedral_degrees(
                                                      pose_in.residue( stop ).xyz( CA ),
                                                      pose_in.residue( stop - 1).xyz( CA ),
                                                      pose_in.residue( stop - 2).xyz( CA ),
                                                      pose_in.residue( stop - 3).xyz( CA ) ) ); 

        
         TR << "Base Dihedral: " << base_dihedral << std::endl;
        
        
        // JQX: the code in the below if statement was in Rosetta 2, but Aroop did not port it into R3
        //      Maybe he has already tested that, the performance was better if using extra 
        //      sequence creteria. But for now, I still port the code in here. Maybe for some reason
        //      one still decides not to use the sequence rules.
        
        bool H3_base_only=false;
        
        if( H3_base_only) {
            std::string base;
            if((base_dihedral > kink_lower_bound) && (base_dihedral < kink_upper_bound)){
                base = "KINK";
                TR << "              " << base << std::endl;
            }
            else if((base_dihedral > extended_lower_bound) && (base_dihedral < extended_upper_bound)){
                base = "EXTENDED";
                TR << "              " << base << std::endl;
            }
            else{
                base = "NEUTRAL";
                TR << "              " << base << std::endl;
            }
            
            return( true );
        }

        
        // setting up pseudo-periodic range used in extended base computation
        if( base_dihedral < kink_lower_bound ){
            base_dihedral = base_dihedral + 360.00;
        }
        
        // Rule 1a for standard kink
        if ((aa_name[aa_name.size()-3] != "ASP") && (aa_name[aa_name.size()-1] == "TRP"))	 //aa_name.size()-3 = n-1
        {                                                                                    //aa_name.size()-1 = n+1
            if( (base_dihedral > kink_lower_bound) && (base_dihedral < kink_upper_bound))
            {
                // std::cout << "KINK Found" << std::endl; // aroop_temp remove
                is_kinked = true;
                is_H3 = true;
            }
        }
        
        // Rule 1b for standard extended form
        if (  ( aa_name[ aa_name.size() - 3 ] == "ASP" ) && 
              ( ( aa_name[1] != "LYS" ) && ( aa_name[1] != "ARG" ) ) &&   
              ( is_H3 != true )    )     //aa_name[1] = 0 position
        {
            
            if( ( base_dihedral>extended_lower_bound) && (base_dihedral<extended_upper_bound) ) {
                // std::cout << "EXTENDED Found" << std::endl; // aroop_temp remove
                is_extended = true;
                is_H3 = true;
            }
            
            if(!is_H3) {
                // Rule 1b extension for special kinked form
                bool is_basic(false); // Special basic residue exception flag
                for(Size ii = 2; ii <= Size(aa_name.size() - 5); ii++) {      //aa_name.size() - 5 = n-3
                    if( aa_name[ii] == "ARG" || aa_name[ii] == "LYS" ) {      //aa_name[2] =  0 position
                        is_basic = true;
                        break;
                    }
                }
                
                if(!is_basic) {
                    Size rosetta_number_of_L49 = pose_in.pdb_info()->pdb2pose(light_chain, 49 );
                    std::string let3_code_L49 = pose_in.residue( rosetta_number_of_L49 ).name3();
                    if( let3_code_L49 == "ARG" || let3_code_L49 == "LYS"){
                        is_basic = true;
                    }
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
                Size rosetta_number_of_L46 = pose_in.pdb_info()->pdb2pose( light_chain, 46 );
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
        if (  (  aa_name[ aa_name.size() - 3 ] == "ASP") &&
              ( ( aa_name[1] == "LYS") || ( aa_name[1] == "ARG" ) ) &&
              ( ( aa_name[0] == "LYS") || ( aa_name[0] == "ARG") ) &&
              ( is_H3 != true )       ) 
        {
            if( ( base_dihedral > extended_lower_bound ) &&
               ( base_dihedral < extended_upper_bound ) ) {
                // aroop_temp remove
                // std::cout << "EXTENDED (w sb) Found" << std::endl;
                is_extended = true;
                is_H3 = true;
            }
        }
        
        TR <<  "Utility: Finished Checking Kink/Extended CDR H3 Base Angle: " << is_H3 << std::endl;
        
        return is_H3;
    } // CDR_H3_filter
    
    
    
    
    //TODO:
    //JQX:
    // What you need is to input the "tf" object, 
    // do something to change the value of this "tf" object
    // right now, it is OK, but JQX must come back to make sure the value of tf 
    // can be changed in this function. If not, maybe this function should return a 
    // pointer
    
    void setup_packer_task(pose::Pose & pose_in, core::pack::task::TaskFactoryOP & tf ) 
    {
		using namespace pack::task;
		using namespace pack::task::operation;
        

        tf->clear();
        tf = new TaskFactory;
        
        
		TR << "Utility: Setting Up Packer Task" << std::endl;
        
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
        

        
		TR << "Utility: Done: Setting Up Packer Task" << std::endl;
        
	} // setup_packer_task



    
    
    
    
    //TODO:
    //JQX:
    // should input a variable here to let the user to adjust 1.9
    
    bool cutpoints_separation( core::pose::Pose & pose, AntibodyInfoOP & antibody_info ) 
    {
        
        bool closed_cutpoints = true;
        
        for( loops::Loops::const_iterator it=antibody_info->all_cdr_loops_.begin(),
            it_end=antibody_info->all_cdr_loops_.end(),
            it_next; it != it_end; ++it ) {
            Size cutpoint   = it->cut();
            Real separation = 10.00; // an unlikely high number
            separation = cutpoint_separation( pose, cutpoint );
            
            if( separation > 1.9 ) {
                closed_cutpoints = false;
                break;
            }
        }
        return( closed_cutpoints );
    } // cutpoints_separation
    
    Real cutpoint_separation(pose::Pose & pose_in, Size cutpoint) {
        
        Size const N ( 1 ); // N atom
        Size const C ( 3 ); // C atom
        
        // Coordinates of the C atom of cutpoint res and N atom of res cutpoint+1
        numeric::xyzVector_float peptide_C(pose_in.residue( cutpoint ).xyz( C )),
        peptide_N( pose_in.residue( cutpoint + 1 ).xyz( N ) );
        //			Real cutpoint_separation=distance(peptide_C, peptide_N);
        Real cutpoint_separation=peptide_C.distance(peptide_N);
        
        return( cutpoint_separation );
    } // cutpoint_separation

    
    
    
    
    
    
    Real global_loop_rmsd (const pose::Pose & pose_in, const pose::Pose & native_pose,loops::LoopOP current_loop ) 
    {
        using namespace scoring;
        
        Size loop_start = current_loop->start();
        Size loop_end = current_loop->stop();
        
        using ObjexxFCL::FArray1D_bool;
        FArray1D_bool superpos_partner ( pose_in.total_residue(), false );
        
        for ( Size i = loop_start; i <= loop_end; ++i ) superpos_partner(i) = true;
        
        using namespace core::scoring;
        Real rmsG = rmsd_no_super_subset( native_pose, pose_in, superpos_partner, is_protein_CA );
        return ( rmsG );
    } 
    
    
    
    
    
    
    std::string get_seq_from_a_loop(core::pose::Pose & pose_in, loops::LoopOP loop )
    {
        std::string seq="";
        for (Size it=loop->start(); it <= loop->stop(); ++it ) {
            seq+=pose_in.residue(it).name1();
        }
        return seq;
    }

    

    
    
    


} // namespace antibody2
} // namespace protocols






