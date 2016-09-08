// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author James Thompson

// libRosetta headers

#include <core/types.hh>

#include <devel/init.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>

#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <numeric/xyzVector.hh>

#include <core/sequence/util.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>

using core::Size;
using core::Real;
using utility::vector1;

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/pose/util.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

int
main( int argc, char * argv [] ) {
	try {

		using namespace core::chemical;
		using namespace basic::options::OptionKeys;
		using namespace basic::options;
		using namespace core::import_pose::pose_stream;

		using core::id::AtomID;
		using core::id::AtomID_Map;
		using utility::vector1;
		using core::id::SequenceMapping;
		using core::sequence::SequenceAlignment;

		devel::init( argc, argv );

		// setup residue types
		ResidueTypeSetCOP rsd_set =
			ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

		PoseInputStreamOP pdb_input( new PDBPoseInputStream( option[ in::file::s ]() ) );

		core::pose::Pose pose1, pose2;
		pdb_input->fill_pose( pose1, *rsd_set );
		pdb_input->fill_pose( pose2, *rsd_set );

		vector1< std::string > align_fns = option[ in::file::alignment ]();
		vector1< SequenceAlignment > alns = core::sequence::read_aln(
			option[ cm::aln_format ](), align_fns.front()
		);

		SequenceAlignment aln = alns.front();
		SequenceMapping mapping = aln.sequence_mapping(1,2);

		vector1< Size > residues;

		if ( option[ in::target_residues ].user() ) {
			// user-specified residues
			residues = option[ in::target_residues ]();
		} else {
			// use all aligned residues
			for ( Size ii = 1; ii <= pose1.size(); ++ii ) {
				residues.push_back(ii);
			}
		}

		AtomID_Map< AtomID > atom_map;
		core::pose::initialize_atomid_map( atom_map, pose1, core::id::BOGUS_ATOM_ID );
		typedef vector1< Size >::const_iterator iter;
		for ( iter it = residues.begin(), end = residues.end(); it != end; ++it ) {
			Size const templ_ii( mapping[*it] );
			if ( templ_ii == 0 ) {
				continue;
			}
			if ( ! pose1.residue(*it).has("CA") ) continue;
			if ( ! pose2.residue(templ_ii).has("CA") ) continue;
			//std::cout << *it << " => " << templ_ii << std::endl;
			AtomID const id1( pose1.residue(*it).atom_index("CA"), *it );
			AtomID const id2( pose2.residue(templ_ii).atom_index("CA"), templ_ii );
			atom_map.set( id1, id2 );
		}

		using core::scoring::superimpose_pose;
		// rmsd numbers look wrong ...
		//core::Real const rmsd( superimpose_pose( pose1, pose2, atom_map ) );
		//std::cout << "superimposed with rmsd of " << rmsd << std::endl;
		superimpose_pose( pose1, pose2, atom_map );
		vector1< std::string > fns = option[ in::file::s ]();
		std::string output_name( fns[1] + ".super.pdb" );
		if ( option[ out::file::o ].user() ) {
			output_name = option[ out::file::o ]();
		}
		pose1.dump_pdb( output_name );
		std::cout << "wrote pdb with name " << output_name << std::endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main
