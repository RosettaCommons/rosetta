// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/dna/MotifLoopBuild.cc
///
/// @brief
/// @author Josh Friedman, jfried23@uw.edu, feb 2012


#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>

#include <devel/dna/MotifLoopBuild.hh>
#include <devel/dna/util_motif_loop.hh>
#include <devel/enzdes/EnzdesRemodelProtocol.hh>
#include <utility/file/FileName.hh>

#include <protocols/match/BumpGrid.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>

#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/MotifDnaPacker.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/MotifSearch.hh>
#include <protocols/motifs/motif_utils.hh>

#include <protocols/dna/DnaDesignDef.hh>
#include <protocols/dna/util.hh>
#include <protocols/dna/DnaChains.hh>
#include <protocols/dna/SeparateDnaFromNonDna.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/looprelax_protocols.hh>
//#include <protocols/loops/loop_mover/refine/LoopMover_Backrub.fwd.hh>

#include <protocols/enzdes/EnzdesFixBBProtocol.hh>
#include <protocols/enzdes/EnzdesBaseProtocol.hh>
#include <protocols/motifs/IRCollection.hh>

#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>

#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>

#include <basic/options/option.cc.gen.hh>

#include <devel/dna/MotifLoopBuildCreator.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.cc.gen.hh>
#include <core/id/types.hh>

#include <core/scoring/rms_util.hh>
#include <numeric/random/random.hh>
#include <time.h>
#include <algorithm>

static basic::Tracer tr("devel.dna.MotifLoopBuild");

//using namespace basic::options;
using namespace protocols;
using namespace devel::enzdes;
using namespace devel::dna;
using namespace basic::options;


static numeric::random::RandomGenerator RG(46117);

MotifLoopBuild::MotifLoopBuild()
{

    //Load in the score function
    scorefxnOP = ( core::scoring::getScoreFunction() );

    //rot_map = build_all_aa_rotamers();
    //Load in the full motif library
	utility::vector1< utility::file::FileName >
    listnames( basic::options::option[ basic::options::OptionKeys::motifs::list_motifs ]().vector() );
	utility::vector1< utility::file::FileName >
    motifnames( motifs::get_filenames( listnames ) );
    motif_lib = motifnames;

    tr << "******" << motif_lib.nmotifs() << " motif library members made!******\n";

}


void MotifLoopBuild::apply(core::pose::Pose & pose)
{
    using namespace core::pack::task;
	pose_=pose;
	mutate_dna( pose );
    mut_pose_=pose;
	(*scorefxnOP)(pose);

    std::string secstruct;
    core::Size num_attemps(0);
    bool pass(false);
    this->set_native_pose( new core::pose::Pose( pose ) );
    this->set_last_move_status( protocols::moves::FAIL_RETRY );

    devel::enzdes::EnzdesRemodelMover::initialize(pose);

    core::id::SequenceMappingOP seqmap( new core::id::SequenceMapping() );
    for( core::Size i(1); i <= pose.total_residue(); ++i) seqmap->push_back(i);
    set_seq_mapping(seqmap);

	setup_cached_observers(pose);

    do
    {
        ++num_attemps;

        if (num_attemps > 10) utility_exit_with_message("Too many attempts at loop closure failed. Quitting.\n");

        build_inv_rots(pose);
        secstruct = generate_secstruct_string(pose);
        pass = devel::enzdes::EnzdesRemodelMover::remodel_pose(pose, secstruct);

        if (! pass)
        {
            tr << "\n****Closure attempt "<< num_attemps << " failed! Trying again...****" << std::endl;

        }
        else
        {
            tr <<"Loop closure sucess with SS " << secstruct <<" on attempt " << num_attemps <<std::endl;
        }

    }
    while ( ! pass );


    std::cout << "\n***** Refining loop model******\n";
    //refine_pose(pose);
    (*scorefxnOP)(pose);



    utility::vector1< core::Size > design = get_flex_region();
    std::list<core::Size> design_list;
    for (utility::vector1<core::Size>::const_iterator i (design.begin()); i != design.end(); ++i)
    {
        design_list.push_back(*i);
    }

    tr << "Loop is now defined as residues " << design[1] << " - " << design[design.size()] << std::endl;

    //Changes to pose length made, must update numbering of the target residues. EnzdesRemodel takes care of the flexable region numbering for us
    if ( (mut_pose_.total_residue() != pose.total_residue() ) &&
        (dna_design_pos_[1] > design[1] ) )
    {
        core::Size offset = pose.total_residue() - mut_pose_.total_residue();

        for (core::Size i(1); i <= dna_design_pos_.size(); ++i )
        {
            dna_design_pos_[i]+=offset;
        }
    }


    tr << "Placing motifs.." << std::endl;
    place_motifs(pose, design, motif_lib);
    pose.remove_constraints();

     tr << "Begin design..." << std::endl;
    std::ostringstream designable;
    core::pack::task::PackerTaskOP seq_design_taskOP = core::pack::task::TaskFactory::create_packer_task(pose);
    seq_design_taskOP->initialize_from_command_line();

    for ( core::Size i(1); i <= pose.total_residue(); ++i)
    {
        if ( (std::find(design.begin(), design.end(), i) == design.end() ) || ( pose.residue(i).name3() != "ALA") )
        {
            seq_design_taskOP->nonconst_residue_task(i).restrict_to_repacking();
        }
        else ( designable << i << ", ");
    }
    tr << "Desiging positions: " << designable.str()  << std::endl;

    protocols::moves::MoverOP design_seqOP( new protocols::simple_moves::PackRotamersMover(scorefxnOP, seq_design_taskOP) );
    design_seqOP->apply(pose);
    pose.dump_pdb("remodel_only.pdb");

    core::pose::Pose before_refine(pose);


    protocols::loops::Loops this_loop = protocols::loops::Loops();
    this_loop.add_loop(design[1], design[ design.size()] );

    for ( core::Size ref(1); ref <= 25; ++ref )
    {
        std::ostringstream oss;
        oss << "final_model_" << ref << ".pdb";

        core::pose::Pose refine_pose( pose);

        protocols::LoopRefine lr = protocols::LoopRefine(this_loop);
        lr.apply(refine_pose);

        core::Real rmsd = core::scoring::all_atom_rmsd( refine_pose, before_refine, design_list);

        refine_pose.dump_scored_pdb( oss.str() , *scorefxnOP);

        tr << "\n RMSD of FINAL REFINEMENT: " << rmsd << "!" <<std::endl;

        core::pose::Pose seperated(refine_pose);
        protocols::dna::SeparateDnaFromNonDna seperator;
        seperator.apply(seperated);

        std::cout << "\n Energy change is " << (*scorefxnOP)(refine_pose)-(*scorefxnOP)(seperated) << " REU\n" << std::endl;
    }

}

MotifLoopBuild::~MotifLoopBuild()
{}

void
MotifLoopBuild::irc_build(core::pose::Pose & pose)
{
    using namespace protocols::motifs;

    utility::vector1<core::Size> design = get_flex_region();
    protocols::loops::LoopsOP this_loop = new protocols::loops::Loops();
    this_loop->add_loop(design[1], design[ design.size()] );
    IRCollection irc = IRCollection(pose, motif_lib, dna_design_pos_);
    irc.incorporate_motifs(pose, this_loop);

}

std::string
MotifLoopBuildCreator::keyname() const
{
	return MotifLoopBuildCreator::mover_name();
}

void MotifLoopBuild::place_motifs( core::pose::Pose & pose,
                                  utility::vector1< core::Size > & flex_pos,
                                  protocols::motifs::MotifLibrary & motif_lib)
{
    protocols::motifs::MotifSearch ms = protocols::motifs::MotifSearch();
    ms.set_motif_library(motif_lib);
    ms.run(pose, flex_pos);

}

protocols::moves::MoverOP
MotifLoopBuildCreator::create_mover() const {
	return new MotifLoopBuild();
}

std::string
MotifLoopBuildCreator::mover_name()
{
	return "MotifLoopBuild";
}


//places inverse rotamers for each design position and writes to "target_inverse_rotamers_"
//in devel::enzdes::EnzdesRemodelProtocol
void MotifLoopBuild::build_inv_rots( core::pose::Pose pose )
{

	utility::vector1< std::list < core::conformation::ResidueCOP > >  inv_rot;

    for ( core::Size restrict_pos(1); restrict_pos <= dna_design_pos_.size(); ++restrict_pos)
    {
            std::list  <core::conformation::ResidueCOP> local_list;

            do {

            motifs::MotifCOPs::const_iterator motifOP_iter;
            core::Size num;

            do {
                motifOP_iter = motif_lib.begin();
                num = RG.random_range(1,  motif_lib.nmotifs()) ;
                std::advance(motifOP_iter, num);

            }
            while (! (*motifOP_iter)->apply_check( pose, dna_design_pos_[restrict_pos] ) );


            std::string aaName( (*motifOP_iter)->restype_name1());


            core::pack::rotamer_set::RotamerSetOP aaRotsOP( get_aa_rotset(aaName));

            const core::Size num_rots( aaRotsOP->num_rotamers() );

            for ( core::Size rot(1); rot <= num_rots; ++rot)
            {
                core::conformation::Residue motif_resi (*(aaRotsOP->nonconst_rotamer(rot)) );

                (*motifOP_iter)->place_residue(pose.residue(dna_design_pos_[restrict_pos]), motif_resi);
                //Now do bump check....

                bool no_collisions = hard_sphere_check(pose, motif_resi);
                if ( no_collisions )
                    local_list.push_back(aaRotsOP->nonconst_rotamer(rot));

            }
            std::cout << "Rotamers befor bump check: " << num_rots <<"\n";
            std::cout << "Rotamers after bump check: " << local_list.size() <<"\n";

        } while ( local_list.size() < 35);

        inv_rot.push_back( local_list );
    }
    set_target_inverse_rotamers(inv_rot);
	return;
}

void MotifLoopBuild::mutate_dna(core::pose::Pose & pose)
{
	tr << "Begin call to mutate_dna()" << std::endl;
	core::scoring::dna::set_base_partner( pose );
	core::pose::Pose tmp_pose( pose );  ////a temporary copy of the pose object
  	protocols::dna::DnaDesignDefOPs design_def;
	utility::vector1< std::string > str_def( basic::options::option[ basic::options::OptionKeys::dna::design::dna_defs ]().vector() );
	protocols::dna::load_dna_design_defs_from_strings(design_def, str_def);
	motifs::make_dna_mutations(pose, design_def);
	core::scoring::dna::set_base_partner(pose);

    dna_design_pos_.clear();

    utility::vector1< core::Size > designed_pos = protocols::motifs::get_target_positions_make_dna_mutations(pose);
    protocols::dna::DnaChains dna_chains = protocols::dna::DnaChains();
    core::scoring::dna::BasePartner partner( core::scoring::dna::retrieve_base_partner_from_pose( pose ) );

    for ( core::Size i(1); i <= designed_pos.size(); ++i)
    {
        core::Size pos1, pos2;
        pos1 = (core::Size) designed_pos[i];
        pos2 = (core::Size) partner[designed_pos[i]];

        dna_design_pos_.push_back(pos1);
        //dna_design_pos_.push_back(pos2);
        tr << "Designed pos: " << pos1 << "  ::  " << pos2 << std::endl;
    }

    return;
}


