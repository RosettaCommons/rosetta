// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/parser/FragmentReader.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/parser/FragmentReader.hh>
#include <protocols/parser/FragSetLoader.hh>

// Package Headers
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/util.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/picking_old/vall/util.hh>

#include <protocols/parser/BluePrint.hh>

#include <utility/exit.hh> // runtime_assert, utility_exit_with_message
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>

#include <core/chemical/ResidueType.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.parser.FragmentReader" );


namespace protocols {
namespace parser {

/// @brief default constructor
FragmentReader::FragmentReader():
	Parent()
{}

/// @brief value constructor
FragmentReader::FragmentReader( TagCOP const & tag ):
	Parent()
{
	parse_tag( tag );
}

/// @brief destructor
FragmentReader::~FragmentReader(){}

/// @brief parse tag
void
FragmentReader::parse_tag( TagCOP const & tag )
{
	/// the way of reading fragments: from pdbs, silent, fragfiles, or vall
	read_type_ = tag->getOption< String >( "type", "" );
	if ( read_type_.empty() ) {
		TR.Error << "No type option. " << std::endl;
		utility_exit_with_message("No type option given");
	}
	if ( read_type_ != "pdb" && read_type_ != "silent" && read_type_ != "vall" && read_type_ != "fragfile" ) {
		TR.Error << "Read type " << read_type_ << " is not registered " << std::endl;
		utility_exit_with_message("Read type not registered.");
	}

	filename_  = tag->getOption< String >( "filename", "" );
	if ( read_type_ == "pdb" || read_type_ == "silent" || read_type_ == "fragfile" ) {
		runtime_assert( !filename_.empty() );
	}

	// length of fragments
	frag_size_ = tag->getOption<Size>( "size", 0 );

	// the begin of sequence positions where fragments are stealed and inserted
	begin_ = tag->getOption<Size>( "begin", 0 );

	// the end of sequence positions where fragments are stealed and inserted
	end_   = tag->getOption<Size>( "end", 0 );

	// number of fragments per position
	nfrags_ = tag->getOption<Size>( "nfrags", 200 );

	if ( read_type_ == "vall" ) {

		runtime_assert( frag_size_ != 0 );

		/// secondary structure assignments to pick fragments from vall
		// read from blueprint
		String const blueprint( tag->getOption<String>( "blueprint", "" ) );
		if ( blueprint != "" ) {
			blueprint_ = protocols::parser::BluePrintOP( new protocols::parser::BluePrint( blueprint ) );
			ss_ = blueprint_->secstruct();
			// pick fragment using sequence information (default false)
			bool use_sequence_bias( tag->getOption<bool>( "use_sequence_bias", 0 ) );
			if ( use_sequence_bias ) {
				aa_ = blueprint_->sequence();
			}
		}

		// using abego definition which is given by blueprint file
		use_abego_ = tag->getOption<bool>( "use_abego", 0 );


		if ( ss_.empty() ) {
			ss_ = tag->getOption<String>( "ss", "" );
		}

		// make sure ss_ is not empty
		runtime_assert( !ss_.empty() );

		// amino acids to pick fragments from vall
		aa_ = tag->getOption<String>( "aa", "" );

		if ( ! aa_.empty() ) {
			runtime_assert( ss_.length() == aa_.length() );
		}
		if ( begin_ == 0 ) {
			TR << "Since option begin is emptry, Fragment is defined from the first residue. " << std::endl;
			begin_ = 1;
		}
		if ( end_ != 0 ) {
			runtime_assert( end_ == begin_ + ss_.length() - 1 );
		}
		end_ = begin_ + ss_.length() - 1;

		TR << "Picking fragments from vall for poistions " << begin_ << "-" << end_
			<< " based on ss=" << ss_ << ", aa=" << aa_ << std::endl;

	} else if ( read_type_ == "pdb" || read_type_ == "silent" ) {

		runtime_assert( frag_size_ != 0 );

		/// number of stealing times
		steal_times_ = tag->getOption<Size>( "steal_times", 1 );

		TR << "Picking Fragments from " << filename_ << " for poistions " << begin_ << "-" << end_ << std::endl;
	}
	runtime_assert( begin_ <= end_ );

}

std::string
FragmentReader::xml_element_name() { return "FragReader"; }

void FragmentReader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	typedef utility::tag::XMLSchemaAttribute Attr;

	XMLSchemaRestriction frag_reader_type;
	frag_reader_type.name( "frag_reader_type" );
	frag_reader_type.base_type( xs_string );
	frag_reader_type.add_restriction( xsr_enumeration, "pdb" );
	frag_reader_type.add_restriction( xsr_enumeration, "silent" );
	frag_reader_type.add_restriction( xsr_enumeration, "vall" );
	frag_reader_type.add_restriction( xsr_enumeration, "fragfile" );
	xsd.add_top_level_element( frag_reader_type );

	AttributeList attributes;
	attributes
		+ required_name_attribute( "The name given to the FragmentReader that will be used when defining FragmentSets" )
		+ Attr( "type", "frag_reader_type", "the source of the fragments: from pdb, silent, fragfile, or vall" )
		+ Attr( "filename", xs_string, "the file from which to read the fragments; required if using the pdb, silent or fragfile input type" )
		+ Attr( "size", xsct_non_negative_integer, "the length of the fragments" )
		+ Attr( "begin", xsct_non_negative_integer, "the beginning sequence positions where fragments are stolen and inserted" )
		+ Attr( "end", xsct_non_negative_integer, "the end of sequence positions where fragments are stolen and inserted" )
		+ Attr::attribute_w_default( "nfrags", xsct_non_negative_integer, "the number of fragments per position", "200" )
		+ Attr( "blueprint", xs_string, "optional secondary structure assignments used when picking fragments from the vall, read from a blueprint file" )
		+ Attr::attribute_w_default( "use_sequence_bias", xsct_rosetta_bool, "pick fragments using sequence information", "0" )
		+ Attr::attribute_w_default( "use_abego", xsct_rosetta_bool, "use abego definitino which is given by blueprint file when picking fragments", "0" )
		+ Attr( "ss", xs_string, "The secondary structure assignment; required if not provided by a blueprint file" )
		+ Attr( "aa", xs_string, "amino acids to pick fragments from the vall" )
		+ Attr::attribute_w_default( "steal_times", xsct_non_negative_integer, "number of times to steal a fragment when generating fragments using either the pdb or silent input types", "1" );

	XMLSchemaComplexTypeGenerator frag_reader_ct;
	frag_reader_ct.element_name( xml_element_name() )
		.complex_type_naming_func( & FragSetLoader::frag_set_loader_ct_namer )
		.add_attributes( attributes )
		.description( "Fragment readers load fragments from one of several sources; you can compose multiple fragment readers to produce a single fragment set" )
		.write_complex_type_to_schema( xsd );
}

/// @brief
void
FragmentReader::set_fragments( Pose const & pose_in, FragSetOP const & fragset )
{
	using namespace core::fragment;

	if ( begin_ == 0 ) {
		core::fragment::steal_frag_set_from_pose( pose_in, *fragset, core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new IndependentBBTorsionSRFD ), frag_size_ ) ) ) );
	} else {
		core::fragment::steal_frag_set_from_pose( pose_in, begin_, end_ , *fragset, core::fragment::FragDataCOP( core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new IndependentBBTorsionSRFD ), frag_size_ ) ) ) );
	}
}

/// @brief
void
FragmentReader::apply( FragSetOP & fragset )
{
	using core::fragment::FrameList;

	if ( read_type_ == "silent" ) {

		using core::chemical::ResidueTypeSetCOP;
		using core::chemical::ChemicalManager;
		using core::chemical::CENTROID;
		using core::import_pose::pose_stream::SilentFilePoseInputStreamOP;
		using core::import_pose::pose_stream::SilentFilePoseInputStream;

		ResidueTypeSetCOP residue_set = ChemicalManager::get_instance()->residue_type_set( CENTROID );
		SilentFilePoseInputStreamOP silent_input( new SilentFilePoseInputStream( filename_ ) );

		Size num( 0 );
		Pose pose_in;
		while ( num++ <= nfrags_ && silent_input->has_another_pose() ) {
			silent_input->fill_pose( pose_in, *residue_set );
			runtime_assert( end_ <= pose_in.size() );
			set_fragments( pose_in, fragset );
		}

	} else if ( read_type_ == "pdb" ) {

		Pose pose_in;
		utility::vector1< String > fs ( utility::string_split( filename_, ',' ) );
		for ( utility::vector1< String>::const_iterator it( fs.begin() ), end( fs.end() ); it!=end; ++it ) {
			String filename( *it );
			core::import_pose::centroid_pose_from_pdb( pose_in, filename , core::import_pose::PDB_file);
			runtime_assert( end_ <= pose_in.size() );
			for ( Size c=0; c<steal_times_; c++ ) {
				set_fragments( pose_in, fragset );
			}
		}

	} else if ( read_type_ == "fragfile" ) {

		using core::fragment::Frame;
		using core::fragment::FrameOP;
		using core::fragment::FragmentIO;
		using core::fragment::ConstFrameIterator;
		using core::fragment::ConstantLengthFragSet;
		using core::fragment::ConstantLengthFragSetOP;
		using core::fragment::OrderedFragSet;
		using core::fragment::OrderedFragSetOP;

		FragSetOP fset = FragmentIO().read_data( filename_ );
		ConstantLengthFragSetOP cf = utility::pointer::dynamic_pointer_cast< core::fragment::ConstantLengthFragSet > ( fset );
		OrderedFragSetOP of = utility::pointer::dynamic_pointer_cast< core::fragment::OrderedFragSet > ( fset );

		FrameList frames;
		if ( cf.get() != NULL && of.get() == NULL ) {
			for ( ConstFrameIterator it=cf->begin(), end( cf->end() ); it!=end; ++it ) {
				frames.push_back( core::fragment::FrameOP( new Frame( **it ) ) );
			}
		} else if ( of.get() != NULL && cf.get() == NULL ) {
			for ( ConstFrameIterator it=of->begin(), end( of->end() ); it!=end; ++it ) {
				frames.push_back( core::fragment::FrameOP( new Frame( **it ) ) );
			}
		} else {
			TR.Error << "FragmentIO returned not proper fragset. See the code." << std::endl;
			utility_exit_with_message("Error getting fragset.");
		}
		fragset->add( frames );

	} else if ( read_type_ == "vall" ) {

		using core::fragment::Frame;
		using core::fragment::FrameOP;
		using core::fragment::IndependentBBTorsionSRFD;
		using core::fragment::picking_old::vall::pick_fragments;
		using core::fragment::picking_old::vall::pick_fragments_by_ss;
		using core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa;

		FrameList frames;
		Size length = end_ - begin_ + 1;
		for ( Size j = 0, je = length; j < je; ++j ) {

			TR << "picking " << nfrags_ << " " << frag_size_ << "-mers for position " << ( begin_ + j ) << std::endl;
			String ss_sub = ss_.substr( j, frag_size_ );
			if ( ss_sub.length() < frag_size_ ) {
				ss_sub.append( frag_size_ - ss_sub.length(), 'D' );
			}

			// make fragments with sequence bias
			String aa_sub;
			if ( !aa_.empty() ) {
				aa_sub = aa_.substr( j, frag_size_ );
				if ( aa_sub.length() < frag_size_ ) {
					aa_sub.append( frag_size_ - aa_sub.length(), '.' );
				}
			} else {
				aa_sub = "";
			}

			// make fragments with abego bias
			utility::vector1< String > abego_sub;
			if ( use_abego_  ) {
				runtime_assert( ss_.length() == blueprint_->abego().size() );
				Size pos( 1 );
				abego_sub.resize( frag_size_ );
				for ( Size ii = j + 1 ; ii <= j + frag_size_; ++ii, ++pos ) {
					if ( ii > blueprint_->abego().size() ) {
						abego_sub[ pos ] = "X";
					} else {
						abego_sub[ pos ] = blueprint_->abego( ii );
					}
				}
			} else {
				abego_sub.clear();
			}

			FrameOP frame( new Frame( begin_ + j, frag_size_ ) );

			frame->add_fragment( pick_fragments( ss_sub, aa_sub, abego_sub, nfrags_, true, IndependentBBTorsionSRFD() ) );

			//  if ( !aa_.empty() ) { // make fragments with sequence bias
			//   String aa_sub = aa_.substr( j, frag_size_ );
			//   if ( aa_sub.length() < frag_size_ ) {
			//    aa_sub.append( frag_size_ - aa_sub.length(), '.' );
			//   }
			//   frame->add_fragment( pick_fragments_by_ss_plus_aa( ss_sub, aa_sub, nfrags_, true, IndependentBBTorsionSRFD() ) );
			//  } else {
			//   frame->add_fragment( pick_fragments_by_ss( ss_sub, nfrags_, true, IndependentBBTorsionSRFD() ) );
			//  }

			frames.push_back( frame );
		}
		fragset->add( frames );

	}

}

} //namespace parser
} //namespace protocols

