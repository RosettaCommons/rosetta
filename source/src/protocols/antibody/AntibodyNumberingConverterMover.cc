// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/AntibodyNumberingConverterMover.cc
/// @brief Converts numbering schemes of an antibody.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/antibody/AntibodyNumberingConverterMover.hh>
#include <protocols/antibody/AntibodyNumberingConverterMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>

// Protocol headers
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/AntibodyNumberingParser.hh>
#include <protocols/antibody/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/string_constants.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.AntibodyNumberingConverterMover" );

namespace protocols {
namespace antibody {

AntibodyNumberingConverterMover::AntibodyNumberingConverterMover():
	protocols::moves::Mover( AntibodyNumberingConverterMover::mover_name() )
{
	set_defaults();
}

AntibodyNumberingConverterMover::AntibodyNumberingConverterMover(
	AntibodyNumberingSchemeEnum from,
	AntibodyNumberingSchemeEnum to
):
	protocols::moves::Mover( AntibodyNumberingConverterMover::mover_name() )
{
	set_defaults();
	from_scheme_ = from;
	to_scheme_ = to;
}

void
AntibodyNumberingConverterMover::set_defaults(){
	using namespace basic::options;

	AntibodyEnumManager manager = AntibodyEnumManager();

	//NOTE - this is to same to make sure we actually set a conversion and we fail if this is still the default.
	from_scheme_ = manager.numbering_scheme_string_to_enum( option[ OptionKeys::antibody::input_ab_scheme]() );

	if ( option[ OptionKeys::antibody::output_ab_scheme].user() ) {
		to_scheme_ = manager.numbering_scheme_string_to_enum( option[ OptionKeys::antibody::output_ab_scheme]() );
	} else {
		to_scheme_ = Chothia_Scheme;
	}
}

///@brief Set the scheme to convert the pose into.
void
AntibodyNumberingConverterMover::set_scheme_conversion(

	AntibodyNumberingSchemeEnum const from,
	AntibodyNumberingSchemeEnum const to

){
	from_scheme_ = from;
	to_scheme_ = to;

}

void
AntibodyNumberingConverterMover::set_from_scheme( AntibodyNumberingSchemeEnum const from ){
	from_scheme_ = from;
}


AntibodyNumberingConverterMover::~AntibodyNumberingConverterMover(){}

void
AntibodyNumberingConverterMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	if ( ! (tag->hasOption("from") && tag->hasOption("to") ) ) {
		utility_exit_with_message("AntibodyNumberingConverterMover: 'from' and 'to' must be specified to do the conversion!");
	}

	AntibodyEnumManager manager = AntibodyEnumManager();

	std::string from = tag->getOption<std::string>("from");
	std::string to = tag->getOption<std::string>("to");


	if ( ! manager.numbering_scheme_is_present(to) ) {
		utility_exit_with_message("Numbering scheme unknown: "+ from);
	}
	if ( ! manager.numbering_scheme_is_present(from) ) {
		utility_exit_with_message("Numbering scheme unknown: "+ to);
	}
	from_scheme_ = manager.numbering_scheme_string_to_enum(from);
	to_scheme_ = manager.numbering_scheme_string_to_enum(to);

}

void
AntibodyNumberingConverterMover::apply( core::pose::Pose & pose)
{
	using namespace core::pose;

	if ( from_scheme_ == to_scheme_ ) {
		TR << "AntibodyNumberingConverterMover: Cannot convert as schemes are identical. " << std::endl;
		return;
	}

	if ( ! ( has_chain( 'L', pose) || has_chain('H', pose) ) ) {
		utility_exit_with_message("AntibodyNumberingConverterMover only works with L and H chains. ");
	}


	AntibodyEnumManagerCOP manager = AntibodyEnumManagerCOP( new AntibodyEnumManager() );
	AntibodyNumberingParser numbering_parser = AntibodyNumberingParser(manager);
	AntibodyNumbering const numbering = numbering_parser.get_antibody_numbering(from_scheme_, North); //North is ambiguous - we really don't care here.

	utility::vector1<PDBLandmarkOP> current_landmarks = numbering.numbering_scheme_transform.at(from_scheme_);
	utility::vector1<PDBLandmarkOP> new_landmarks = numbering.numbering_scheme_transform.at(to_scheme_);

	//Renumbering logic is always so much fun.  Appologies for the confusing nature of all this.

	//core::Size last_L_resnum = core::pose::chain_end_res(pose, get_chain_id_from_chain('L', pose));
	//core::Size last_H_resnum = core::pose::chain_end_res(pose, get_chain_id_from_chain('H', pose));



	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {

		core::Size new_resnum = 0;
		char new_icode = ' ';

		core::Size chain_num = pose.chain(i);
		char chain = get_chain_from_chain_id( chain_num, pose );

		if ( ! ( chain == 'L' || chain == 'H' ) ) {
			continue;
		}

		PDBLandmark query_landmark = PDBLandmark(chain, pose.pdb_info()->number( i ), pose.pdb_info()->icode( i ) );

		bool landmark_missing = true;

		PDBLandmarkCOP new_landmark = get_matching_landmark(numbering, query_landmark, from_scheme_, to_scheme_);

		//Resnum is the same.  Nothing to be done.
		if ( *new_landmark == query_landmark ) {
			continue;
		}

		if ( new_landmark->resnum() != 0 ) {

			//TR << "Renumbering resnum "<< i << ", <" <<query_landmark.get_string() << "> to <" << new_landmark->get_string() << ">" << std::endl;

			landmark_missing = false;
			new_resnum = new_landmark->resnum();
			new_icode = new_landmark->insertion_code();

		}
		//Here, we don't have a conversion.

		//Posibilities:
		// 1) Res is part of constant region.  Here, we don't have conversion info, and we really don't care.
		//
		//     1b)  This residue does not have an insertion code,
		//         + Dealt with properly. Continue the numbering from before.
		//
		//     1a)  This residue has an insertion code, but we wouldn't want to continue the insertion operator.
		//
		//         + Dealth with properly
		//
		// 2) Res is an insertion in a CDR or framework.
		//
		//         + Dealt with properly.
		//
		// 3) Res is not found in conversion, but is the first residue of the pose.  We fail. Should never happen.
		//
		//         ! Should never happen.  Unknown how to deal with properly.  Utility_Exit, refer to PyIgClassify.
		//
		// 4) Res is not found in conversion, but is the first residue of the L or H chain.
		//
		//         ! Should never happen. Unknown how to deal with properly. Utility_exit, refer to PyIgClassify.
		//

		if ( landmark_missing ) {
			//Extend it if we are at the end.
			if ( i != 1 ) {
				core::Size last_resnum = pose.pdb_info()->number( i -1 );
				char const & last_icode = pose.pdb_info()->icode( i - 1 );

				///////////////////////////////////
				//Need to know if next replacement is the same new_resnum.  If it is, then we have to START to use insertion codes
				//So, we need to match figure out when the next match is and if that match is the same chain and same resnum as last_resnum.

				bool needs_new_icode = false;

				//We also want to know if we are really at the end of possible matches.  If so, we don't want to increment icode, even if the last resnum has an icode.

				bool downstream_matches = false;

				for ( core::Size xx = i; xx <= pose.total_residue(); ++xx ) {
					core::Size xx_chain = pose.chain(xx); //Here because for some reason we get an int for chain and I don't have time to refactor this all over Rosetta.

					if ( xx_chain != chain_num ) break; //End of chain

					PDBLandmark xx_landmark = PDBLandmark(chain, pose.pdb_info()->number( xx ), pose.pdb_info()->icode( xx ));
					PDBLandmarkCOP xx_match_landmark = get_matching_landmark(numbering, xx_landmark, from_scheme_, to_scheme_);

					if ( xx_match_landmark->resnum() == 0 ) {
						continue;
					} else {
						downstream_matches = true;
						if ( last_resnum == (xx_match_landmark->resnum() - 1) ) {
							needs_new_icode = true;
						}
						break;
					}

				}
				///////////////////////////////////

				//New Icode
				if ( needs_new_icode ) {
					new_resnum = last_resnum;
					new_icode = utility::ALPHANUMERICS.at(0);
				} else if ( last_icode == ' ' || ! downstream_matches ) {
					//Increment Resnum
					new_resnum = last_resnum +1;
					new_icode = ' ';
				} else {
					//Increment Icode
					std::size_t pos = utility::ALPHANUMERICS.find( last_icode );
					if ( pos!=std::string::npos && pos != utility::ALPHANUMERICS.size() -1 ) {
						new_icode = utility::ALPHANUMERICS.at(pos+1);
					} else {
						utility_exit_with_message( "Not enough icodes for CDR length! Last icode"+utility::to_string(last_icode));
					}
				}
			} else {
				utility_exit_with_message("First resnum not found.  Renumber into a numbering scheme first.  Try PyIgClassify. "+query_landmark.get_string());
			}
		}
		if ( new_resnum == 0 ) {
			utility_exit_with_message("ERROR in conversion. Should never reach here!");
		}

		pose.pdb_info()->number(i, new_resnum  );
		pose.pdb_info()->icode( i, new_icode );


		//Double check that we set everything properly
		core::Size resnum = pose.pdb_info()->pdb2pose(pose.pdb_info()->chain(i), new_resnum, new_icode);
		if ( resnum == 0 ) {
			std::string msg = "ERROR in conversion -  Residue not found in pose: " +
				utility::to_string(new_landmark->resnum()) + " " + new_landmark->chain() + " " + new_landmark->insertion_code() + "\n";
			utility_exit_with_message(msg);
		}

		//This does not need to be fast code, so we can check every time.
		debug_assert(resnum == i);

	}
	TR << "Renumbering complete." << std::endl;
}

protocols::moves::MoverOP
AntibodyNumberingConverterMover::clone() const
{
	return protocols::moves::MoverOP( new AntibodyNumberingConverterMover( *this ) );
}

protocols::moves::MoverOP
AntibodyNumberingConverterMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new AntibodyNumberingConverterMover );
}

// XRW TEMP std::string
// XRW TEMP AntibodyNumberingConverterMover::get_name() const
// XRW TEMP {
// XRW TEMP  return AntibodyNumberingConverterMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AntibodyNumberingConverterMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "AntibodyNumberingConverterMover";
// XRW TEMP }

void
AntibodyNumberingConverterMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, AntibodyNumberingConverterMover const & mover )
{
	mover.show(os);
	return os;
}



/////////////// Creator ///////////////

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AntibodyNumberingConverterMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new AntibodyNumberingConverterMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AntibodyNumberingConverterMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return AntibodyNumberingConverterMover::mover_name();
// XRW TEMP }

std::string AntibodyNumberingConverterMover::get_name() const {
	return mover_name();
}

std::string AntibodyNumberingConverterMover::mover_name() {
	return "AntibodyNumberingConverterMover";
}

void AntibodyNumberingConverterMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	XMLSchemaRestriction ABnumbering_enum;
	ABnumbering_enum.name("ab_numbering_schemes");
	ABnumbering_enum.base_type(xs_string);

	AntibodyEnumManager ABnumbering_manager;
	for ( auto& ABnumbering : ABnumbering_manager.get_recognized_numbering_schemes() ) {
		ABnumbering_enum.add_restriction(xsr_enumeration, ABnumbering);
	}
	xsd.add_top_level_element(ABnumbering_enum);

	attlist + XMLSchemaAttribute::required_attribute(
		"from", "ab_numbering_schemes",
		"Current numbering scheme");

	attlist + XMLSchemaAttribute::required_attribute(
		"to", "ab_numbering_schemes",
		"Numbering scheme to be converter to");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Converts numbering schemes of an antibody, independent of AntibodyInfo. "
		"By default, works on a SINGLE antibody FAB with chains L and H, as the rest of Rosetta",
		attlist );
}

std::string AntibodyNumberingConverterMoverCreator::keyname() const {
	return AntibodyNumberingConverterMover::mover_name();
}

protocols::moves::MoverOP
AntibodyNumberingConverterMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AntibodyNumberingConverterMover );
}

void AntibodyNumberingConverterMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AntibodyNumberingConverterMover::provide_xml_schema( xsd );
}


} //protocols
} //antibody

