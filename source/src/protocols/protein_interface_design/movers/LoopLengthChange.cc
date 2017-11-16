// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/LoopLengthChange.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/LoopLengthChange.hh>
#include <protocols/protein_interface_design/movers/LoopLengthChangeCreator.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <core/chemical/ResidueType.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/kinematics/FoldTree.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
//#include <time.h> // THIS IS USED FOR DEBUGGING!
//#include <algorithm> // THIS IS USED FOR DEBUGGING!
//#include <string> // THIS IS USED FOR DEBUGGING!

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace std;
using namespace core::scoring;

static basic::Tracer TR( "protocols.protein_interface_design.movers.LoopLengthChange" );

// XRW TEMP std::string
// XRW TEMP LoopLengthChangeCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return LoopLengthChange::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP LoopLengthChangeCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new LoopLengthChange );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LoopLengthChange::mover_name()
// XRW TEMP {
// XRW TEMP  return "LoopLengthChange";
// XRW TEMP }

LoopLengthChange::LoopLengthChange() :
	Mover( LoopLengthChange::mover_name() ),
	loop_start_( 0 ), loop_end_( 0 ), delta_( 0 ), tail_segment_(false)/*, restype_char_('A') (in .hh)*/
{
}


LoopLengthChange::~LoopLengthChange() {}

void
LoopLengthChange::apply( core::pose::Pose & pose )
{
	//    TR<<"Changing loop "<<loop_start()<<"-"<<loop_end()<<" by "<<delta()<<std::endl;
	//    if ( loop_end() + delta() < loop_start() ){
	//        pose.dump_pdb("dumping_as_loop_end_greater_than_loop_start.pdb");
	//    }

	//time_t rawtime;
	//time (&rawtime);
	//std::string timerstring = ctime(&rawtime);
	//TR<<"Now dumping at time "<<timerstring<<std::endl;
	//timerstring.erase(std::remove(timerstring.begin(), timerstring.end(), '\n'), timerstring.end());
	//timerstring.erase(std::remove(timerstring.begin(), timerstring.end(), ' '), timerstring.end());
	//pose.dump_pdb(timerstring+"_dump.pdb");

	//TR<<"loop end is "<<loop_end()<<" loop_start is "<<loop_start()<<" delta is "<<delta()<<std::endl;

	runtime_assert( loop_end() >= loop_start() );
	//runtime_assert( loop_end() + delta() >= loop_start() );
	core::kinematics::FoldTree const & ft = pose.fold_tree();
	core::Size jump_count( 0 );
	if ( delta() < 0 ) {
		TR << "DEBUG: I will try to delete residues" << std::endl;
		TR << "Loop end is " << loop_end() << " delta is " << delta() << std::endl;


		for ( int del(0); del>delta(); --del ) {
			//TR<<"loop_end -delta()+ del =" <<loop_end()-delta()+ del<<std::endl;
			if ( tail_segment_ ) {
				TR << "DEBUG: I will now delete residue " << loop_end() + del << std::endl;
				//string res = static_cast<ostringstream*>( &(ostringstream() << loop_end() + del) )->str();
				//pose.dump_pdb("llc_res_before_del_res_"+res+".pdb");
				pose.delete_polymer_residue( loop_end() + del ); // there will never be jumps within the tail segment
			} else {
				if ( ft.is_jump_point( loop_end() + del - delta() - jump_count ) ) {
					TR<<"LoopLengthChange called across a jump. I'm skipping the jump residue"<<std::endl;
					jump_count++;
				}

				TR << "DEBUG: I will now delete residue " << loop_end() + del - delta() - jump_count << std::endl;
				//string res = static_cast<ostringstream*>( &(ostringstream() << loop_end() + del - delta() - jump_count) )->str();
				//pose.dump_pdb("llc_res_before_del_res_"+res+".pdb");
				pose.delete_polymer_residue( loop_end() + del - delta() - jump_count ); // was before pose.delete_polymer_residue( loop_end() + del - delta() + jump_count );
				//pose.dump_pdb("llc_res_after_del_res_"+res+".pdb");
			}
			//   pose.conformation().insert_ideal_geometry_at_polymer_bond( loop_start() );
			//   pose.conformation().insert_ideal_geometry_at_polymer_bond( loop_start() + 1 );
		}
	} else if ( delta() > 0 ) {
		using namespace core::chemical;
		using namespace core::conformation;

		ResidueTypeSetCOP residue_set( pose.residue_type_set_for_pose() );
		runtime_assert(oneletter_code_specifies_aa(restype_char())); //ensure that we aren't about to try this with aa_unk
		ResidueCOP new_res = ResidueFactory::create_residue( residue_set->name_map( name_from_aa( aa_from_oneletter_code( restype_char() ) ) ) );
		if ( tail_segment_ ) {
			/* This will be ugly: For some reason you cannot use .append_polymer_residue_after_seqpos to add residues
			*        to the very last position. Instead we are here adding n+1 residues to the second last position.  */
			for ( core::Size leng(1); leng<=(core::Size) delta() + 1; ++leng ) {
				pose.conformation().append_polymer_residue_after_seqpos( *new_res, loop_end()-2+leng, true/*build_ideal*/); // start growing the new loop from res loopEnd-1
			}
			// Now we need to delete the last residue
			core::Size unwantedLastResidueNo =  loop_end()  + delta() + 1;
			//TR << "Delta is " << delta() << std::endl;
			//pose.dump_pdb("before_delete.pdb");
			pose.conformation().delete_polymer_residue( unwantedLastResidueNo );
			//pose.dump_pdb("after_delete.pdb");
		} else {
			for ( core::Size leng(1); leng<=(core::Size) delta(); ++leng ) {
				pose.conformation().prepend_polymer_residue_before_seqpos( *new_res, loop_end()+1, true/*build_ideal*/);
				// string res = static_cast<ostringstream*>( &(ostringstream() << leng ))->str();
				// pose.dump_pdb("llc_res_"+res+".pdb");
				// pose.set_omega(loop_end()+leng-1,180.0);
			}
		}


	}
	pose.update_residue_neighbors();
	if ( pose.pdb_info() ) pose.pdb_info()->obsolete( true );
}

// XRW TEMP std::string
// XRW TEMP LoopLengthChange::get_name() const {
// XRW TEMP  return LoopLengthChange::mover_name();
// XRW TEMP }

void
LoopLengthChange::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose )
{
	loop_start( core::pose::parse_resnum( tag->getOption< std::string >( "loop_start" ), pose ) );
	loop_end( core::pose::parse_resnum( tag->getOption< std::string >( "loop_end" ), pose ) );
	delta( tag->getOption< int >( "delta" ) );

	runtime_assert( loop_end() > loop_start() );
	// runtime_assert( loop_end() + delta() >= loop_start() );

	if ( tag->hasOption("restype") ) {
		restype_char( tag->getOption< char >( "restype" ) );
		if ( !(core::chemical::oneletter_code_specifies_aa(restype_char())) ) {
			utility_exit_with_message(std::string("LoopLengthChange restype argument apparently does not specify a canonical residue type (ACDEFGHIKLMNPQRSTVWY), argument was ") + std::string( 1, restype_char()));
			//manual casting because string literal plus char doesn't work as expected; http://stackoverflow.com/questions/25812411/something-in-c
		}
		if ( delta() <= 0 ) {
			utility_exit_with_message("LoopLengthChange restype argument invalid when delta is less than zero (you are removing residues; do not specify that residue type to add!)");

		}
	}

	TR<<"LoopLengthChange with loop "<<loop_start()<<"-"<<loop_end()<<" and delta "<<delta()<<std::endl;
}

protocols::moves::MoverOP
LoopLengthChange::clone() const {
	return( protocols::moves::MoverOP( new LoopLengthChange( *this ) ));
}

void
LoopLengthChange::loop_start( core::Size const loop_start ){
	loop_start_ = loop_start;
}

core::Size
LoopLengthChange::loop_start() const{
	return( loop_start_ );
}

void
LoopLengthChange::loop_end( core::Size const l ){
	loop_end_ = l;
}

core::Size
LoopLengthChange::loop_end() const{
	return( loop_end_ );
}

void
LoopLengthChange::loop_cut( core::Size const d ){
	loop_cut_ = d;
}

core::Size
LoopLengthChange::loop_cut() const{
	return( loop_cut_ );
}

void
LoopLengthChange::restype_char( char const restype_char ){
	restype_char_ = restype_char;
}

char
LoopLengthChange::restype_char() const{
	return( restype_char_ );
}

void
LoopLengthChange::delta( int const d ){
	delta_ = d;
}

int
LoopLengthChange::delta() const{
	return( delta_ );
}
void
LoopLengthChange::tail( bool b ){
	tail_segment_ = b;
}

std::string LoopLengthChange::get_name() const {
	return mover_name();
}

std::string LoopLengthChange::mover_name() {
	return "LoopLengthChange";
}

void LoopLengthChange::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::required_attribute( "loop_start", xsct_refpose_enabled_residue_number, "Starting residue number for loop, formatted in seqpos or PDB or refpose numbering" )
		+ XMLSchemaAttribute::required_attribute( "loop_end", xsct_refpose_enabled_residue_number, "Ending residue number for loop, formatted in seqpos or PDB or refpose numbering" )
		+ XMLSchemaAttribute::required_attribute( "delta", xs_integer, "Number of residues to extend or contract the loop" )
		+ XMLSchemaAttribute("restype", xsct_canonical_res_char, "Only valid if delta is positive. Char of residue type to insert, on [ACDEFGHIKLMNPQRSTVWY] (the canonicals)." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "A Mover for changing the length of a loop; does not close the loop afterwards.", attlist );
}

std::string LoopLengthChangeCreator::keyname() const {
	return LoopLengthChange::mover_name();
}

protocols::moves::MoverOP
LoopLengthChangeCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopLengthChange );
}

void LoopLengthChangeCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopLengthChange::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
