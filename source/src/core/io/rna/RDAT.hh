// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/rna/RDAT.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_rna_data_RDAT_HH
#define INCLUDED_core_scoring_rna_data_RDAT_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/io/rna/RDAT.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.fwd.hh>
#include <string>
////////////////////////////////////////////////////////////////////////////////////////
//
// RDAT = RNA Data Analysis Toolkit file format
//
// Described here: http://bioinformatics.oxfordjournals.org/content/28/22/3006  (Cordero et al., 2012)
// with updated specs here: http://rmdb.stanford.edu/repository/specs/
//
// Currently the standard i/o file format for RNA structure mapping data in the Das lab
//  and in EteRNA.
//
// Does this belong somewhere in core/io ?
//
////////////////////////////////////////////////////////////////////////////////////////
namespace core {
namespace io {
namespace rna {

	typedef std::pair< std::string, std::string > Annotation;

	class RDAT: public utility::pointer::ReferenceCount {

	public:

		//constructor
		RDAT();

		//constructor
		RDAT( std::string const filename );

		//constructor
		RDAT( std::string const name, std::string const sequence,
					core::Size const offset,
					utility::vector1< Size > const & seqpos,
					std::string const structure,
					utility::vector1< Annotation > const & annotations,
					utility::vector1< utility::vector1< Annotation > > const & data_annotations,
					utility::vector1< utility::vector1< core::Real > > const & reactivity,
					utility::vector1< utility::vector1< core::Real > > const & reactivity_error,
					utility::vector1< std::string > const & comments	);

		//destructor
		~RDAT();

	public:

		void
		output_rdat_to_file( std::string const filename ) const;

		void
		read_rdat_file( std::string const filename );

		void
		fill_header_information( pose::Pose & pose );

		void
		fill_data_from_features( utility::vector1< char > const & seqchars,
														 utility::vector1< core::Size > const & resnums,
														 utility::vector1< std::string > const & feature_names,
														 utility::vector1< utility::vector1< core::Real > > const & all_feature_vals );

		void
		output_data( utility::io::ozstream & out ) const;

		bool
		check_rdat() const;

		std::string const & version() const { return  version_; }
		std::string const & name() const { return  name_; }
		std::string const & sequence() const { return  sequence_; }
		std::string const & structure() const { return  structure_; }
		utility::vector1< std::string > const & sequences() const { return  sequences_; }
		utility::vector1< std::string > const & structures() const { return  structures_; }
		int const & offset() const { return  offset_; }
		utility::vector1< int > const & seqpos() const { return  seqpos_; }
		utility::vector1< Annotation > const & annotations() const { return  annotations_; }
		utility::vector1< utility::vector1< Annotation > > const & data_annotations() const { return  data_annotations_; }
		utility::vector1< utility::vector1< core::Real > > const & reactivity() const { return  reactivity_; }
		utility::vector1< utility::vector1< core::Real > > const & reactivity_error() const { return  reactivity_error_; }
		utility::vector1< std::string > const & comments() const { return  comments_; }
		utility::vector1< core::Real > const & xsel() const { return  xsel_; }
		utility::vector1< utility::vector1< core::Real > > const & xsel_refine() const { return  xsel_refine_; }
		utility::vector1< utility::vector1< core::Real > > const & trace() const { return  trace_; }

		void set_version( std::string const & setting ) {  version_ = setting; }
		void set_name( std::string const & setting ) {  name_ = setting; }
		void set_sequence( std::string const & setting ) {  sequence_ = setting; }
		void set_structure( std::string const & setting ) {  structure_ = setting; }
		void set_sequences( utility::vector1< std::string > const & setting ) {  sequences_ = setting; }
		void set_structures( utility::vector1< std::string > const & setting ) {  structures_ = setting; }
		void set_offset( int const & setting ) {  offset_ = setting; }
		void set_seqpos( utility::vector1< int > const & setting ) { seqpos_ = setting; }
		void set_annotations( utility::vector1< Annotation > const & setting ) {  annotations_ = setting; }
		void set_data_annotations( utility::vector1< utility::vector1< Annotation > > const & setting ) {  data_annotations_ = setting; }
		void set_reactivity( utility::vector1< utility::vector1< core::Real > > const & setting ) {  reactivity_ = setting; }
		void set_reactivity_error( utility::vector1< utility::vector1< core::Real > > const & setting ) {  reactivity_error_ = setting; }
		void set_comments( utility::vector1< std::string > const & setting ) {  comments_ = setting; }
		void set_xsel( utility::vector1< core::Real > const & setting ) {  xsel_ = setting; }
		void set_xsel_refine( utility::vector1< utility::vector1< core::Real > > const & setting ) {  xsel_refine_ = setting; }
		void set_trace( utility::vector1< utility::vector1< core::Real > > const & setting ) {  trace_ = setting; }

	private:

		void
		output_rdat_header( utility::io::ozstream & out ) const;

		void
		fill_sequences_and_structures();

		void
		fill_sequences_if_empty();

		void
		fill_structures_if_empty();

		void
		fill_if_empty( std::string & data_string,
									 utility::vector1< std::string > & data_strings,
									 std::string const tag ) const;

		void
		save_data( std::string const & line, utility::vector1< utility::vector1< Real > > & var  ) const;

		void
		save_data( std::string const & line, utility::vector1< Real > & var  ) const;

		void
		save_data_with_idx( utility::vector1< std::string > & var,
													Size const idx,  std::string const & value ) const;

		void
		save_data_with_idx( utility::vector1< utility::vector1< Annotation > > & var,
												Size const idx,  Annotation & value ) const;

		Annotation
		get_annotation( std::string const tag ) const;

		utility::vector1< std::string >
		str2cell( std::string const s ) const;

		std::string
		remove_tag( std::string & line, std::string const tag ) const;

		void
		fill_data_annotations_if_empty();

		void
		fill_seqpos( std::string const & seqpos_info,
								 std::string & sequence_seqpos );

		bool
		check_sequence_seqpos( std::string const & sequence_seqpos ) const;

		bool
		check_annotations( utility::vector1< Annotation > const & annotations ) const;

		bool
		check_annotation( Annotation const & annotation ) const;

	private:

		std::string version_;
		std::string name_;
		std::string sequence_;
		std::string structure_;

		// per entry
		utility::vector1< std::string > sequences_;
		utility::vector1< std::string > structures_;

		// int to add to go from 1,2,... N to conventional numbering
		int offset_;
		utility::vector1< int > seqpos_;
		utility::vector1< Annotation > annotations_;
		utility::vector1< utility::vector1< Annotation > > data_annotations_;

		// data matrices.
		utility::vector1< utility::vector1< core::Real > > reactivity_;
		utility::vector1< utility::vector1< core::Real > > reactivity_error_;

		utility::vector1< std::string > comments_;

		// to be deprecated in the future, probably.
		utility::vector1< core::Real > xsel_;
		utility::vector1< utility::vector1< core::Real > > xsel_refine_;
		utility::vector1< utility::vector1< core::Real > > trace_;

	};


	std::string
	get_tag( utility::vector1< Annotation > const & annotations, std::string const tag );

	utility::vector1< std::string >
	get_tags( utility::vector1< utility::vector1< Annotation > > const & data_annotations, std::string const tag );


} //rna
} //io
} //core

#endif
