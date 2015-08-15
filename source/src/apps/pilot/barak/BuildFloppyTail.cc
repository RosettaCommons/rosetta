// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// (C) 199x-2008 Hebrew University, Jerusalem
//
/// @file   BuildPeptide.cc
//
/// @brief Application that reads in a peptides sequence file and outputs a linear peptide.
/// @author Nir London
/// @date June 01, 2009

//#define GL_GRAPHICS

// Package Headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/disulfide_util.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/scoring/Energies.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/options/option.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/run.OptionKeys.gen.hh>
#include <core/util/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>

// Numeric Headers

// C++ headers
#include <utility>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/Constraint.hh>


// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/util/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/options/util.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

// C++ headers

using core::util::T;
using core::util::Error;
using core::util::Warning;


static core::util::Tracer TR( "BuildFloppyTail" );

using namespace core;
using namespace core::options;
using namespace OptionKeys;


////////////////////////////////////////////////////////////////////////////////
/// @details Given a Pose, a protein sequence where each character represents an
/// amino acid, and a ResidueTypeSet, give the Pose a conformation of covalently
/// linked residues that match the sequence. NOTE: support making pose from a
/// fully annotated sequence now, that is, for each residue variant or ligand
/// which cannot be deduced from one letter code directly, a [] is added
/// directly following the one letter code containig the residue's fullname, e.g.
/// K[lys:NtermProteinFull]ADFGCH[HIS_D]QNVE[glu:CtermProteinFull]Z[ZN].
/// This allows a pose to be constructed with full features from a silent output
/// file, such as with distiguished HIS tautomers, various chain termini and
/// cutpoint variants etc. Currently not working with disulfide variant CYD, but
/// this is on to-do list.
void append_sequence_to_pose(
	pose::Pose & pose,
	std::string const & sequence_in,
	chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini /* true */
)
{
	using namespace core::chemical;
 	typedef core::Size Size;

	// grab residue types
	TR << "append_seq_to_pose - start" << std::endl;
	ResidueTypeCOPs requested_types = residue_types_from_sequence( sequence_in, residue_set, auto_termini );
	assert( annotated_to_oneletter_sequence( sequence_in ).length() == requested_types.size() );
	TR << "read types" << std::endl;

	// make the pose
	for ( Size i = 1, ie = requested_types.size(); i <= ie; ++i ) {
		// grab the new residue
		ResidueType const & rsd_type = *requested_types[ i ];
		core::conformation::ResidueOP new_rsd( NULL );
		new_rsd = conformation::ResidueFactory::create_residue( rsd_type );
		TR << "residue created" << std::endl;

		// yab 20090219: The following error check was in the original
		// code prior to the split into residue_types_from_sequence()
		// and this function, but it doesn't appear to be triggerable
		// because ResidueFactory always returns a residue.  I leave it
		// in for now, but consider taking it out.
		if ( !new_rsd ) {
			std::cerr << "cannot create a residue that matches the residue type "
				<< rsd_type.name1() << " " << rsd_type.name() << " at position " << i << '\n';
			utility_exit_with_message( "make_pose_from_sequence fails\n" );
		}

		TR << "append_sequence_to_pose():  seqpos: " << i << " " << new_rsd->aa() << std::endl;

		// do the actual append
		if ( rsd_type.has_variant_type( LOWER_TERMINUS ) )
					utility_exit_with_message( "did not expect a lower terminus residue\n" );
		if( new_rsd->aa() == aa_unk || new_rsd->aa() == aa_vrt ) {
				TR.Error << "found unknown aminoacid or X in sequence at position " << i <<  std::endl;
				if ( i< ie ) {
					utility_exit_with_message( "found unknown aminoacid or X in sequence - not supported\n" );
				}
		}
		pose.append_residue_by_bond( *new_rsd, true );
	}


	TR << "sequence in pose: " << pose.sequence() << std::endl;
	TR << "annotated seq: " << pose.annotated_sequence() << std::endl;

} // make_pose_from_sequence


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

	using namespace pose;
	using namespace scoring;
	using namespace conformation;
        using namespace core::chemical;

	//setup random numbers and options
	devel::init(argc, argv);

        //create a pose
        pose::Pose pose;

        //original protein pose
      	io::pdb::pose_from_pdb( pose, options::start_file() ); // gets filename from -s option
				Size origPoseLen = pose.total_residue();

        //read peptides fasta file
        std::string pepSeq = core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1]->sequence();
      	Size seqLen = pepSeq.length();

        remove_upper_terminus_type_from_pose_residue(pose,origPoseLen);
        append_sequence_to_pose(pose, pepSeq, *ChemicalManager::get_instance()->residue_type_set( "fa_standard" ), false);
        add_upper_terminus_type_to_pose_residue(pose,pose.total_residue());

				TR << "old, new, seqLen length: " << origPoseLen << ", " << pose.total_residue() << ", " << seqLen << std::endl;
				//        add_lower_terminus_type_to_pose_residue(pose,1);
				runtime_assert(seqLen + origPoseLen == pose.total_residue());

        //make peptide linear
        for (Size i=origPoseLen + 1; i <= pose.total_residue(); i++) {
            pose.set_phi(i,-135.0);
            pose.set_psi(i,135.0);
            pose.set_omega(i,180.0);
        }

        //dump pdb to output
        pose.dump_pdb("./peptide.pdb");

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
