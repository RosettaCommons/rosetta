#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http:#www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   create_a3b_hbs.py
## @brief  Scan the sidechain angles of a beta-AA
## @author Andrew Watkins

# Here is a sketch of the basic flow of the program...
#
# Pertubation Phase
#   +-Monte Carlo Mover---------------------------------------+
#   | +-Random Mover-( 1 / 2 / 1 / 1 )--------------------+ | |
#   | | +-Docking Mover-----------------------------------+ | |
#   | | | small rigid body movements between the peptide  | | |
#   | | | and protein for conformational diversity        | | |
#   | | +-------------------------------------------------+ | |
#   | | +-Peptide Modeling--------------------------------+ | |
#   | | | move peptide with small/shear moves to generate | | |
#   | | | conformational diversity                        | | |
#   | | +-------------------------------------------------+ | |
#   | +-----------------------------------------------------+ |
#   | +-Rotamer Trials Mover--------------------------------+ |
#   | | quick sidechain packing to find optimal rotamers    | |
#   | | before the next cycle                               | |
#   | +-----------------------------------------------------+ |
#   +---------------------------------------------------------+
#
# Design Minimization Phase
#   | +-Pack Rotamers Mover---------------------------------+ |
#   | | repack and design rotamers to explore sequence     | |
#   | | space                                              | |
#   | +-----------------------------------------------------+ |
#   | +-Minimization Mover----------------------------------+ |
#   | | energy minimize the current conformation            | |
#   | +-----------------------------------------------------+ |

# application specific options
# pert options
#RealOptionKey const mc_temp( "hddm::mc_temp" );
#RealOptionKey const pert_mc_temp( "hddm::pert_mc_temp" );
#RealOptionKey const pert_dock_rot_mag( "hddm::pert_dock_rot_mag" );
#RealOptionKey const pert_dock_trans_mag( "hddm::pert_dock_trans_mag" );
#RealOptionKey const pert_pep_small_temp( "hddm::pert_pep_small_temp" );
#RealOptionKey const pert_pep_small_H( "hddm::pert_pep_small_H" );
#RealOptionKey const pert_pep_small_L( "hddm::pert_pep_small_L" );
#RealOptionKey const pert_pep_small_E( "hddm::pert_pep_small_E" );
#RealOptionKey const pert_pep_shear_temp( "hddm::pert_pep_shear_temp" );
#RealOptionKey const pert_pep_shear_H( "hddm::pert_pep_shear_H" );
#RealOptionKey const pert_pep_shear_L( "hddm::pert_pep_shear_L" );
#RealOptionKey const pert_pep_shear_E( "hddm::pert_pep_shear_E" );

#IntegerOptionKey const pert_pep_num_rep( "hddm::pert_pep_num_rep" );
#IntegerOptionKey const pert_num( "hddm::pert_num" );
#IntegerOptionKey const dock_design_loop_num( "hddm::dock_design_loop_num" );

#BooleanOptionKey const final_design_min( "hddm::final_design_min" );
#BooleanOptionKey const use_soft_rep( "hddm::use_soft_rep" );
#BooleanOptionKey const mc_initial_pose( "hddm::mc_initial_pose" );
#BooleanOptionKey const hbs_design_first( "hddm::hbs_design_first" );

#BooleanOptionKey const pymol( "hddm::pymol" );
#BooleanOptionKey const keep_history( "hddm::keep_history" );

# design options
#RealOptionKey const desn_mc_temp( "hddm::desn_mc_temp" );
def hbs_dock_design(pose):
    nddp = NcbbDockDesignProtocol()
    protocols.ncbb.setup_filter_stats()
    nddp.apply(pose)
	#NcbbDockDesignProtocol(
	#	core::scoring::ScoreFunctionOP score_function,
	#	core::Real const mc_temp,
	#	core::Real const pert_mc_temp,
	#	core::Real const pert_dock_rot_mag,
	#	core::Real const pert_dock_trans_mag,
	#	core::Real const pert_pep_small_temp,
	#	core::Real const pert_pep_small_H,
	#	core::Real const pert_pep_small_L,
	#	core::Real const pert_pep_small_E,
	#	core::Real const pert_pep_shear_temp,
	#	core::Real const pert_pep_shear_H,
	#	core::Real const pert_pep_shear_L,
	#	core::Real const pert_pep_shear_E,
#
#		core::Size const pert_pep_num_rep,
#		core::Size const pert_num,
#		core::Size const dock_design_loop_num,
#
#		bool const no_design,
#		bool const final_design_min,
#		bool const use_soft_rep,
#		bool const mc_initial_pose,
#		bool const ncbb_design_first,
#
#		bool const pymol,
#		bool const keep_history
#
#	);

if __name__ == '__main__':
    # Argparse

    #option.add( hddm::mc_temp, "The temperature to use for the outer loop of the HDDM protocol. Defaults to 1.0." ).def( 1.0 );
    #option.add( hddm::pert_mc_temp, "The temperature to use for the pertubation phase of the HDDM protocol. Defaults to 0.8." ).def( 0.8 );
    #option.add( hddm::pert_dock_rot_mag, "The rotation magnitude for the ridged body pertubation in the pertubation phase of the HDDM protocol. Defaults to 1.0." ).def( 1 );
    #option.add( hddm::pert_dock_trans_mag, "The translation magnitude for the ridged body pertubation in the pertubation phase of the HDDM protocol. Defaults to 0.5." ).def( 0.5 );
    #option.add( hddm::pert_pep_small_temp, "" ).def( 0.8 );
    #option.add( hddm::pert_pep_shear_temp, "" ).def( 0.8 );

    #option.add( hddm::pert_pep_small_H, "" ).def( 2.0 );
    #option.add( hddm::pert_pep_small_L, "" ).def( 2.0 );
    #option.add( hddm::pert_pep_small_E, "" ).def( 2.0 );
    #option.add( hddm::pert_pep_shear_H, "" ).def( 2.0 );
    #option.add( hddm::pert_pep_shear_L, "" ).def( 2.0 );
    #option.add( hddm::pert_pep_shear_E, "" ).def( 2.0 );

    #option.add( hddm::pert_pep_num_rep, "Number of small and shear iterations for the peptide" ).def( 100 );
    #option.add( hddm::pert_num, "Number of iterations of perturbation loop per design" ).def(10);
    #option.add( hddm::dock_design_loop_num, "Number of iterations of pertubation and design" ).def(10);

    #option.add( hddm::final_design_min, "Do a final repack/design and minimization. Default true" ).def(true);
    #option.add( hddm::use_soft_rep, "Use soft repulsion for pertubation and initial design. Default false" ).def(false);
    #option.add( hddm::mc_initial_pose, "Allow initial pose to be considered as lowest energy pose. Default false" ).def(false);
    #option.add( hddm::hbs_design_first, "Design before pertubation (want when initial struct is aligned to hotspot)  Default false" ).def(false);

    #option.add( hddm::pymol, "Set up pymol mover. Default false" ).def(false);
    #option.add( hddm::keep_history, "Keep history in pymol. Requires hddm::pymol set to true. Default false" ).def(false);

    #option.add( hddm::desn_mc_temp, "The temperature to use for the design/minimization phase of the HDDM protocol. Defaults to 0.8." ).def( 0.8 );
    
    from pyrosetta import *
    from pyrosetta.rosetta import *

    init()
    
    hbs_dock_design() # args

