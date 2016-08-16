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

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>


#include <core/conformation/Conformation.hh>

#include <utility/vector1.hh>

#include <sstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <devel/domain_assembly/DomainAssemblyReader.hh>
#include <devel/domain_assembly/domain_assembly_setup.hh>


// option key includes

#include <basic/options/keys/DomainAssembly.OptionKeys.gen.hh>

#include <core/pose/util.hh>

namespace devel {
namespace domain_assembly {


using namespace core;
using namespace kinematics;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace chemical;

// // protocol specific options
// namespace DomainAssembly{
//  FileOptionKey const da_setup_option_file( "DomainAssembly::da_setup_option_file" );
//  FileOptionKey const da_setup_output_pdb( "DomainAssembly::da_setup_output_pdb" );
//   //FileOptionKey const da_start_pdb( "DomainAssembly::da_start_pdb" );
//  //FileOptionKey const da_linker_file( "DomainAssembly::da_linker_file" );
//  //IntegerOptionKey const da_start_pdb_num( "DomainAssembly::da_start_pdb_num" );
//  //IntegerOptionKey const da_nruns( "DomainAssembly::da_nruns" );
// }

static THREAD_LOCAL basic::Tracer TR_da( "DomainAssemblySetup" );

/// @brief adds linkers and/or truncates a domain
///  instructions are contained in the DomainInfo member variables
void
DomainInfo::process_domain( )
{
	pose::Pose temp_pose = input_pose_;

	// trim the N-terminus
	if ( trim_nterm_ > 0 ) {
		temp_pose.conformation().delete_residue_range_slow( 1 , trim_nterm_ );
	}

	//for ( Size i = 1; i <= trim_nterm_; ++i ) {
	// temporary FoldTree manipulation - until delete_polymer fixed
	//FoldTree f( temp_pose.fold_tree() );
	//f.delete_seqpos( 1 ) ;
	//f.reorder( 2 );
	//temp_pose.fold_tree( f );
	//temp_pose.delete_polymer_residue( 1 ) ;
	//}

	// trim the c-terminus
	for ( Size i = 1; i <= trim_cterm_; ++i ) {
		temp_pose.delete_polymer_residue( temp_pose.total_residue() ) ;
	}


	ResidueTypeSetCOP residue_set ( ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );

	// add n-terminal linker if requested
	std::string linker_sequence_n = Nterm_linker_;
	if ( linker_sequence_n.length() > 0 ) {
		core::pose::remove_lower_terminus_type_from_pose_residue( temp_pose, 1 );
	}
	for ( Size pos = linker_sequence_n.length(); pos > 0; --pos ) {
		char aa = linker_sequence_n[pos-1]; // string indexing is zero-based!
		AA my_aa = aa_from_oneletter_code( aa );
		ResidueType const & rsd_type( *( residue_set->get_representative_type_aa( my_aa ) ) );
		conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( rsd_type ) );
		temp_pose.prepend_polymer_residue_before_seqpos( *new_rsd, 1, true );
		//if the residue to be added is RNA, give it the average torsion angles for A-RNA
		if ( rsd_type.is_NA() ) {
			temp_pose.set_alpha(1, -68);
			temp_pose.set_beta(1, 178);
			temp_pose.set_gamma(1, 54);
			temp_pose.set_delta(1, 82);
			temp_pose.set_epsilon(1, -153);
			temp_pose.set_zeta(1, -71);
		} else {

			// give the new residue extended chain parameters
			temp_pose.set_psi( 1, 150 );
			temp_pose.set_omega( 1, 180 );
			temp_pose.set_phi( 1, -150 );
		}
	}

	// add c-terminal linker if requested
	std::string linker_sequence = Cterm_linker_;
	//SML HACK: don't forget to remove the cterm patch from the last residue
	//using core::pose::remove_variant_type_from_pose_residue;
	//core::pose::remove_variant_type_from_pose_residue(temp_pose, core::chemical::UPPER_TERMINUS, temp_pose.total_residue());
	if ( linker_sequence.length() > 0 ) {
		core::pose::remove_upper_terminus_type_from_pose_residue( temp_pose, temp_pose.total_residue());
	}
	for ( Size pos = 0; pos < linker_sequence.length(); ++pos ) {
		char aa = linker_sequence[pos]; // string indexing is zero-based!
		AA my_aa = aa_from_oneletter_code( aa );
		ResidueType const & rsd_type( *( residue_set->get_representative_type_aa( my_aa )) );
		conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( rsd_type ) );
		temp_pose.append_residue_by_bond( *new_rsd, true );
		Size seqpos = temp_pose.total_residue( );
		//if the residue to be added is RNA, give it the average torsion angles for A-RNA
		if ( rsd_type.is_NA() ) {
			temp_pose.set_alpha(seqpos, -68);
			temp_pose.set_beta(seqpos, 178);
			temp_pose.set_gamma(seqpos, 54);
			temp_pose.set_delta(seqpos, 82);
			temp_pose.set_epsilon(seqpos, -153);
			temp_pose.set_zeta(seqpos, -71);
		} else {

			// give the new residue extended chain parameters
			//Size seqpos = temp_pose.total_residue( );
			temp_pose.set_psi( seqpos-1, 150 );
			temp_pose.set_omega( seqpos-1, 180 );
			temp_pose.set_phi( seqpos, -150 );
		}
	}
	processed_pose_ = temp_pose;
}

/// @brief calls process domain to add linkers/truncate domains
void
process_domains(
	utility::vector1< DomainInfo > & domains
)
{
	for ( Size i = 1; i <= domains.size(); ++i ) {
		domains[i].process_domain( );
	}
}

/// @brief connect the domains
///If using with pdbs files that contain RNA.  The RNA must either be at the beginning of the first pdb input
///or at the end of the last pdb input, otherwise the connect domains function will fail and crash
void
connect_domains(
	utility::vector1< DomainInfo > domains,
	pose::Pose & full_pose
)
{
	Size domain_begin(0), domain_end(0), domains_so_far_end(0);
	for ( Size i = 1; i <= domains.size(); ++i ) {
		pose::Pose domain_pose = domains[i].get_processed_pose();
		if ( full_pose.total_residue() == 0 ) {
			full_pose = domain_pose;
			domain_begin = 1;
			domain_end = full_pose.total_residue();
			//Keep track of the size of the first domain
			domains_so_far_end = full_pose.total_residue();
		} else {
			core::pose::remove_upper_terminus_type_from_pose_residue( full_pose, full_pose.total_residue( ) );
			core::pose::remove_lower_terminus_type_from_pose_residue( domain_pose, 1 );
			domain_begin = full_pose.total_residue( ) + 1;;
			//Will add sequential amino acid residues to the first domain until no more amino acids are detected
			for ( Size i = 1; i <= domain_pose.total_residue( ); ++i ) {
				if ( !((domain_pose.residue_type(i)).is_NA()) ) {
					full_pose.append_residue_by_bond( domain_pose.residue(i) );
					domain_end = i;
				}
			}
			//Adds the first nucleotide to the connected domains pose using a jump and anchoring them to the C-terminus of the connected domains
			//Only attempted if there are nucleic acid resiudes at the end of the pdb.
			if ( full_pose.total_residue() - domains_so_far_end != domain_pose.total_residue() ) {
				full_pose.append_residue_by_jump(domain_pose.residue( (domain_end + 1)), full_pose.total_residue());
			}
			//Adds the nucleotide residues at the end of the final domain pdb that were ignored previously, will crash if only one nucleotide at the end of the final domain pdb. Crude by why would you have just one RNA nucleotide at the end.
			for ( Size i = domain_end + 2; i <= domain_pose.total_residue(); ++i ) {
				full_pose.append_residue_by_bond( domain_pose.residue(i) );
			}
			full_pose.conformation().insert_ideal_geometry_at_polymer_bond( domain_begin - 1 );

			full_pose.set_phi  ( domain_begin - 1, -150.0 );
			full_pose.set_psi  ( domain_begin - 1,  150.0 );
			full_pose.set_omega( domain_begin - 1,  180.0 );
			full_pose.set_phi  ( domain_begin, -150.0 );
			full_pose.set_psi  ( domain_begin,  150.0 );
			full_pose.set_omega( domain_begin,  180.0 );
			domain_end = full_pose.total_residue();
			domains_so_far_end = full_pose.total_residue();
		}
		domains[i].set_domain_begin( domain_begin );
		domains[i].set_domain_end( domain_end );

	} // for loop over domains
}

/// @brief  A stand alone setup protocol for domain assembly.
///  Reads in a set of pdb files and connects them with specified
///  linkers.  The multidomain pose is output as a pdb file that
///  can be used as the input file for assemble_domains_optimize()
void
assemble_domains_setup()
{
	utility::vector1< DomainInfo > domains;
	pose::Pose full_pose;

	std::string option_filename = option[ DomainAssembly::da_setup_option_file ]();
	parse_da_option_file( domains, option_filename );

	TR_da << "This many domains will be connected " << domains.size() << std::endl;

	// add linkers and/or truncate domains
	process_domains( domains );


	// connect the domains
	connect_domains( domains, full_pose );

	std::string output_pdb = option[ DomainAssembly::da_setup_output_pdb ]();
	full_pose.dump_pdb( output_pdb );
}

}
}
