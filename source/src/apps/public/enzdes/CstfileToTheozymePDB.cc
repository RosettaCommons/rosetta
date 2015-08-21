// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/public/enzdes/CstFileToTheozymePDB.cc
/// @brief
/// @author Florian Richter, floric@u.washington.edu, june 2010


#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/match.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/util.hh>
#include <protocols/toolbox/match_enzdes_util/AlignPoseToInvrotTreeMover.hh>
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTree.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

void
create_theozyme_pdb();

int main( int argc, char * argv [] )
{
	try {
		devel::init( argc, argv );

		create_theozyme_pdb();

	} catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

void
create_theozyme_pdb()
{

	basic::Tracer tr( "apps.public.enzdes.CstfileToTheozymePDB.cc" );

	std::string cstfile_name( basic::options::option[basic::options::OptionKeys::match::geometric_constraint_file]() );

	protocols::toolbox::match_enzdes_util::EnzConstraintIOOP enz_io( new protocols::toolbox::match_enzdes_util::EnzConstraintIO( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ) );
	enz_io->read_enzyme_cstfile( cstfile_name );

	protocols::toolbox::match_enzdes_util::InvrotTreeOP invrot_tree( new protocols::toolbox::match_enzdes_util::TheozymeInvrotTree( enz_io ) );
	invrot_tree->generate_targets_and_inverse_rotamers();

	//test whether the input cstfile is in the same directory or somewhere else
	std::string::size_type const slash_loc = cstfile_name.find_last_of( '/' );
	std::string outname_base;
	if ( slash_loc == std::string::npos ) { // same directory
		outname_base = "PDB_Model_"+cstfile_name;
	} else {
		outname_base = "PDB_Model_"+cstfile_name.substr(slash_loc+1, cstfile_name.size() );
	}
	invrot_tree->dump_invrots_tree_as_multimodel_pdbs( outname_base );

	//stealth functionality: can also use this app to read in a pose
	//and align it to the theozyme
	if ( basic::options::option[basic::options::OptionKeys::in::file::s].user() ) {
		utility::vector1< std::string > input_files = basic::options::start_files();
		if ( input_files.size() == 1 ) {
			core::pose::PoseOP pose( new core::pose::Pose() );
			core::import_pose::pose_from_pdb( *pose, input_files[ 1 ] );
			protocols::toolbox::match_enzdes_util::AllowedSeqposForGeomCstOP allowed_seqpos( new protocols::toolbox::match_enzdes_util::AllowedSeqposForGeomCst() );
			allowed_seqpos->initialize_from_command_line( pose );
			protocols::toolbox::match_enzdes_util::AlignPoseToInvrotTreeMover align_pose( invrot_tree, allowed_seqpos);
			align_pose.set_add_target_to_pose( true );
			align_pose.apply( *pose );
			pose->dump_pdb("theozyme_tree_align.pdb");
		}
	}
}

