// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RestrictDesignToProteinDNAInterface.cc
/// @brief
/// @author ashworth

#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterfaceCreator.hh>
#include <protocols/dna/DnaChains.hh>
#include <protocols/dna/DnaDesignDef.hh>
#include <protocols/dna/util.hh> // find_basepairs
#include <protocols/dna/DnaInterfaceFinder.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
using utility::vector1;
#include <utility/string_util.hh>
using utility::string_split;

#include <ObjexxFCL/format.hh>

// option key includes

#include <basic/options/keys/dna.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>


using namespace ObjexxFCL::format;

namespace protocols {
namespace dna {

using namespace core;
using namespace chemical;
using namespace conformation;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace pose;
using namespace scoring;
using namespace constraints;
using namespace utility::tag;

using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static THREAD_LOCAL basic::Tracer TR( "protocols.dna.RestrictDesignToProteinDNAInterface", t_info );

TaskOperationOP RestrictDesignToProteinDNAInterfaceCreator::create_task_operation() const
{
	return TaskOperationOP( new RestrictDesignToProteinDNAInterface );
}

void RestrictDesignToProteinDNAInterfaceCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictDesignToProteinDNAInterface::provide_xml_schema( xsd );
}

std::string RestrictDesignToProteinDNAInterfaceCreator::keyname() const
{
	return RestrictDesignToProteinDNAInterface::keyname();
}

RestrictDesignToProteinDNAInterface::RestrictDesignToProteinDNAInterface()
: parent(),
	dna_chains_(/* 0 */),
	interface_(/* 0 */),
	reference_pose_(/* 0 */),
	base_only_( true ),
	forget_chains_and_interface_( true ),
	z_cutoff_( 0.0 ),
	close_threshold_( 10. * 10. ),
	contact_threshold_( 3.7 * 3.7 )
{}

RestrictDesignToProteinDNAInterface::~RestrictDesignToProteinDNAInterface() = default;

TaskOperationOP RestrictDesignToProteinDNAInterface::clone() const
{
	return TaskOperationOP( new RestrictDesignToProteinDNAInterface( *this ) );
}

void
RestrictDesignToProteinDNAInterface::copy_dna_chains( DnaChainsCOP dna_chains )
{
	dna_chains_ = DnaChainsOP( new DnaChains( *dna_chains ) );
}

DnaChainsCOP
RestrictDesignToProteinDNAInterface::dna_chains() const { return dna_chains_; }

void
RestrictDesignToProteinDNAInterface::copy_targeted_dna( DnaDesignDefOPs const & targeted_dna )
{
	targeted_dna_ = targeted_dna;
}

DnaDesignDefOPs const &
RestrictDesignToProteinDNAInterface::targeted_dna() const { return targeted_dna_; }

void
RestrictDesignToProteinDNAInterface::copy_interface( DnaInterfaceFinderCOP interface )
{
	interface_ = DnaInterfaceFinderOP( new DnaInterfaceFinder( *interface ) );
}

DnaInterfaceFinderCOP
RestrictDesignToProteinDNAInterface::interface() const { return interface_; }

void
RestrictDesignToProteinDNAInterface::set_reference_pose( PoseCOP pose )
{
	reference_pose_ = pose;
}

PoseCOP
RestrictDesignToProteinDNAInterface::reference_pose() const { return reference_pose_; }

void
RestrictDesignToProteinDNAInterface::parse_tag( TagCOP tag , DataMap & )
{
	typedef utility::vector1< std::string > Strings;
	if ( tag->hasOption("dna_defs") ) {
		//targeted_dna_.clear() // pros/cons?
		std::string const defs( tag->getOption< std::string >("dna_defs") );
		if ( defs == "COMMANDLINE" ) {
			TR(t_debug) << " parsing dna definitions from command line" << std::endl;
			load_dna_design_defs_from_options( targeted_dna_ );
		} else {
			TR(t_debug) << " parsing dna definitions: " << defs << std::endl;
			Strings const str_defs( string_split( defs, ',' ) );
			load_dna_design_defs_from_strings( targeted_dna_, str_defs );
		}
		TR(t_debug) << "targeted DNA is: " << targeted_dna_ << std::endl;
	}
	if ( tag->hasOption("base_only") ) base_only_ = tag->getOption< bool >("base_only");
	if ( tag->hasOption("z_cutoff") ) z_cutoff_ = tag->getOption< Real >("z_cutoff");
	if ( tag->hasOption("close_threshold") ) {
		close_threshold_ = tag->getOption< Real >("close_threshold");
	}
	if ( tag->hasOption("contact_threshold") ) {
		close_threshold_ = tag->getOption< Real >("contact_threshold");
	}
	if ( tag->hasOption("forget_chains_and_interface") ) {
		forget_chains_and_interface_ = tag->getOption< bool >("forget_chains_and_interface");
	}
}

void RestrictDesignToProteinDNAInterface::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute( "dna_defs", xs_string )
		+ XMLSchemaAttribute( "base_only", xs_boolean )
		+ XMLSchemaAttribute( "z_cutoff", xs_decimal )
		+ XMLSchemaAttribute( "close_threshold", xs_decimal )
		+ XMLSchemaAttribute( "contact_threshold", xs_decimal )
		+ XMLSchemaAttribute( "forget_chains_and_interface", xs_boolean );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}

/// @brief determines the DNA interface residues and informs a PackerTask of their appropriate packing behavior
/// @details
/// Step 2: get info about DNA chains and set up DNA packing behavior
/// Step 3: Determine protein-DNA interface
/// Step 4: apply any new restrictions to resfile pack/design settings, and any existing constraints
/// Step 5: report
void
RestrictDesignToProteinDNAInterface::apply(
	Pose const & pose,
	PackerTask & ptask
) const
{
	using namespace task;

	if ( z_cutoff_ == 0.0 ) z_cutoff_ = option[ OptionKeys::dna::design::z_cutoff ]();

	// This class does not try to automatically configure residues whose corresponding ResidueLevelTask do not have the behavior "AUTO". This behavior should be set for protein residues, either manually or by using another kind of TaskOperation.

	/*
	// if no resfile was specified,
	// assume all protein positions are "AUTO" positions (behavior automated)
	if ( ! option[ OptionKeys::packing::resfile ].user() ) {
	for ( Size i(1), end( ptask.total_residue() ); i <= end; ++i ) {
	ResidueLevelTask & rtask( ptask.nonconst_residue_task(i) );
	if ( pose.residue_type(i).is_protein() ) rtask.add_behavior("AUTO");
	}
	}
	*/

	Size const nres( pose.size() );
	runtime_assert( nres == ptask.total_residue() );


	/// Step 2: get info about DNA chains and set up DNA packing behavior

	// get basepairing info unless supplied by user
	if ( !dna_chains_ ) {
		dna_chains_ = DnaChainsOP( new DnaChains );
		find_basepairs( pose, *dna_chains_ );
	}

	// configure the packer task for any specified targeted top-stranded DNA basepair positions
	if ( ! targeted_dna_.empty() ) {
		// (by default, targeted_dna_ is an empty vector, so none of the following applies)
		 for ( auto const & def : targeted_dna_ ) {
			Size index( def->pdbpos );
			if ( pose.pdb_info() ) {
				// if pose has PDB numbering and chain info, assume DNA defs refer to them
				PDBPoseMap const & pdb_pose_map( pose.pdb_info()->pdb2pose() );
				index = pdb_pose_map.find( def->chain, def->pdbpos );
			}
			if ( ! pose.residue_type( index ).is_DNA() ) {
				std::cerr << "ERROR: DNA design def " << *def << " indicates a non-DNA position"
					<< std::endl; utility_exit();
			} else if ( ! dna_chains_->is_top( index ) ) {
				std::cerr << "ERROR: DNA design def " << *def << " is DNA but is not in the 'top' strand"
					<< std::endl; utility_exit();
			}
			DnaPosition const & pos( (*dna_chains_)[ index ] );
			runtime_assert( index == pos.top() );
			ResidueLevelTask & toptask( ptask.nonconst_residue_task( pos.top() ) );
			toptask.add_behavior("TARGET");
			if ( pos.paired() ) ptask.nonconst_residue_task( pos.bottom() ).add_behavior("TARGET");

			if ( ! def->name3.empty() ) {
				// specifying the appropriate ResidueType here is tricky, because there are multiple possible 'name3's for the nucleotides, AND because we must make sure to indicate a ResidueType that is already represented in the ResidueLevelTask
				ResidueTypeSetCOP rts( pose.residue(1).residue_type_set() );
				// a list of all existing residue types that match the input name3
				// AMW: can this be improved?
				ResidueTypeCOPs const & name3map( ResidueTypeFinder( *rts ).name3( def->name3 ).get_all_possible_residue_types() );
				// use the first ResidueType represented in the ResidueLevelTask that corresponds to one in the name3 map
				for ( auto
						allowed_type( toptask.allowed_residue_types_begin() ),
						end( toptask.allowed_residue_types_end() ); allowed_type != end; ++allowed_type ) {
					if ( std::find( name3map.begin(), name3map.end(), *allowed_type ) != name3map.end() ) {
						toptask.target_type( *allowed_type );
						break;
						//TR(t_info) << "Setting target type " << (*def)->name3
						//   << " at position " << pos.top() << std::endl;
					}
				}
				if ( ! toptask.target_type() ) {
					TR(t_info) << "Error: target type " << def->name3
						<< " does not correspond to an allowed type at position " << pos.top()
						<< std::endl; utility_exit();
				}
				if ( pos.paired() ) {
					std::string const comp_name3( dna_comp_name_str( def->name3 ) );
					ResidueTypeCOPs const & name3map_comp( ResidueTypeFinder( *rts ).name3( comp_name3 ).get_all_possible_residue_types() );
					ResidueLevelTask & bottask( ptask.nonconst_residue_task( pos.bottom() ) );
					for ( auto
							allowed_type( bottask.allowed_residue_types_begin() ),
							end( bottask.allowed_residue_types_end() ); allowed_type != end; ++allowed_type ) {
						if ( std::find( name3map_comp.begin(), name3map_comp.end(), *allowed_type ) !=
								name3map_comp.end() ) {
							bottask.target_type( *allowed_type );
							break;
							//TR(t_info) << "Setting target type " << comp_name3 << " at position " << pos.bottom()
							//   << std::endl;
						}
					}
					if ( ! bottask.target_type() ) {
						TR(t_info) << "Error: target type " << comp_name3
							<< " does not correspond to an allowed type at position " << pos.bottom()
							<< std::endl; utility_exit();
					}
				}
			}
		}

		// when targeting particular basepair positions, packing for all non-specified DNAs is disabled
		for ( Size i(1); i <= ptask.total_residue(); ++i ) {
			ResidueLevelTask & rtask( ptask.nonconst_residue_task(i) );
			if ( pose.residue_type(i).is_DNA() && !rtask.has_behavior("TARGET") ) {
				rtask.prevent_repacking();
			}
		}
	}

	/// Step 3: Determine protein-DNA interface

	// find protein-dna interface (unless user supplied one)
	if ( !interface_ ) {

		vector1< Size > dna_design_positions;
		bool limit_by_DNA(false);
		for ( Size i(1); i <= nres; ++i ) {
			if ( !pose.residue_type(i).is_DNA() ) continue;
			ResidueLevelTask const & rtask( ptask.residue_task(i) );
			// scan positions count as targeted
			if ( rtask.has_behavior("TARGET") || rtask.has_behavior("SCAN") ) {
				limit_by_DNA = true;
				break;
			}
		}
		if ( !limit_by_DNA ) {
			// there are no specifically targeted dna positions: use entire DNA interface
			for ( Size dpos(1); dpos <= nres; ++dpos ) {
				if ( !pose.residue_type(dpos).is_DNA() ) continue;
				dna_design_positions.push_back( dpos );
				// turn off all DNA repacking (...?)
				//ptask.temporarily_set_pack_residue( dpos, false );
			}
		} else {
			// targeted dna design positions exist: limit interface to within the vicinity of these
			for ( auto & it : *dna_chains_ ) {
				Size resid( it.first ); DnaPosition & dnapos( it.second );
				ResidueLevelTask & toptask( ptask.nonconst_residue_task( resid ) );
				if ( !toptask.has_behavior("TARGET") && !toptask.has_behavior("SCAN") ) continue;
				dna_design_positions.push_back( resid );
				// verbose/debug
				TR(t_info) << "\nTargeting DNA at position ";
				if ( pose.pdb_info() ) {
					TR << pose.pdb_info()->chain( dnapos.top() ) << "."
						<< pose.pdb_info()->number( dnapos.top() ) << '\n';
				} else {
					TR << pose.chain( dnapos.top() ) << "." << dnapos.top() << '\n';
				}
				TR(t_debug) << "Allowed types:\n"; toptask.print_allowed_types( TR(t_debug) );

				if ( !dnapos.paired() ) continue;
				dna_design_positions.push_back( dnapos.bottom() );
				ResidueLevelTask & bottask( ptask.nonconst_residue_task( dnapos.bottom() ) );
				// verbose/debug
				TR(t_info) << "\nTargeting DNA at position ";
				if ( pose.pdb_info() ) {
					TR << pose.pdb_info()->chain( dnapos.bottom() ) << "."
						<< pose.pdb_info()->number( dnapos.bottom() ) << '\n';
				} else {
					TR << pose.chain( dnapos.bottom() ) << "." << dnapos.bottom() << '\n';
				}
				TR(t_debug) << "Allowed types:\n"; bottask.print_allowed_types( TR(t_debug) );
			}
		}

		vector1< Size > protein_positions;
		for ( Size p_index(1); p_index <= nres; ++p_index ) {
			if ( ! ptask.pack_residue( p_index ) ) continue; // already disabled
			if ( pose.residue_type( p_index ).is_DNA() ) continue; // ignore DNA
			if ( ! pose.residue_type( p_index ).is_protein() ) {
				ptask.nonconst_residue_task( p_index ).prevent_repacking(); // not protein, disable
			} else if ( ptask.residue_task( p_index ).has_behavior("AUTO") ) {
				protein_positions.push_back(p_index); // decide how to pack/design this position below
			}
		}

		interface_ = DnaInterfaceFinderOP( new DnaInterfaceFinder( close_threshold_, contact_threshold_, z_cutoff_, base_only_ ) );
		// perform arginine rotamer sweep to decide wether or not to pack/design protein residues
		interface_->determine_protein_interface( pose, protein_positions, dna_design_positions );
	}

	/// Step 4: apply any new restrictions to resfile pack/design settings, and any existing constraints
	bool const repack_only( option[ OptionKeys::dna::design::repack_only ]() );
	ConstraintSetCOP constraint_set( pose.constraint_set() );
	 for ( auto const & itr : interface_->protein_neighbors() ) {
		Size const index( itr.first );
		DnaNeighbor const & neighbor( itr.second );

		ResidueLevelTask & restask( ptask.nonconst_residue_task( index ) );
		// residues with constraints will not be designed
		bool const cst( constraint_set && constraint_set->residue_pair_constraints_exists( index ) );
		if ( neighbor.contact() && !cst && !repack_only ) continue; // designable, no restriction
		// restriction to current and/or original amino acid type
		vector1< bool > aas( num_canonical_aas, false );
		AA current_aa( pose.residue_type(index).aa() );
		if ( neighbor.close() ) aas[ current_aa ] = true;
		if ( reference_pose_ ) {
			// intended to reflect the true native type, if desired/necessary
			// such as for 'reversion to [original] wildtype' during iterative design
			AA orig_aa( reference_pose_->residue_type(index).aa() );
			if ( current_aa != orig_aa ) aas[ orig_aa ] = true;
		}
		// aas.empty() returns true if all of its values are false
		if ( !aas.empty() ) restask.restrict_absent_canonical_aas( aas );
		else restask.prevent_repacking();
	}

	/// Step 5: report
	for ( Size i(1), end( ptask.total_residue() ); i <= end; ++i ) {
		ResidueLevelTask const & rlt( ptask.residue_task(i) );
		if ( pose.residue_type(i).is_DNA() ) {
			if ( ! rlt.being_packed() ) continue;
			TR(t_info) << "DNA types allowed at ";
			if ( pose.pdb_info() ) {
				TR << pose.pdb_info()->number(i) << " "
					<< pose.pdb_info()->chain(i);
			} else TR << i << " " << pose.chain(i);
			TR << ":";
			// set used here to avoid redundant name3's from extra adduct variant allowed_types
			std::set< std::string > name3set;
			for ( auto
					allowed_type( rlt.allowed_residue_types_begin() ),
					typesend( rlt.allowed_residue_types_end() ); allowed_type != typesend; ++allowed_type ) {
				name3set.insert( (*allowed_type)->name3() );
			}
			 for ( auto const & n3 : name3set ) {
				TR << n3 << ",";
			}
			TR << '\n';
			continue;
		}
		if ( ! rlt.being_packed() ) continue;
		if ( pose.pdb_info() ) {
			TR(t_info) << pose.pdb_info()->chain(i) << "." << pose.pdb_info()->number(i) << ".";
		} else {
			TR(t_info) << pose.chain(i) << "." << i << ".";
		}
		TR << pose.residue(i).name3();
		if ( rlt.being_designed() ) TR(t_info) << " is DESIGNABLE";
		else TR(t_info) << " is packable";
		TR(t_info) << '\n';
	}
	TR(t_info) << std::endl;
	if ( forget_chains_and_interface_ ) {
		// abandon exisiting dna chains and interface info to prevent accumulation of state
		TR.Debug << "forgetting chains and interface" << std::endl;
		dna_chains_.reset();
		interface_.reset();
	}
}

} // namespace dna
} // namespace protocols

