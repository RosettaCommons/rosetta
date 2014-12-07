// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/SilentStruct.hh
///
/// @brief silent input file reader for mini
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_SilentStruct_hh
#define INCLUDED_core_io_silent_SilentStruct_hh

#include <core/io/silent/SilentStruct.hh>

// unit headers
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentEnergy.hh>

#include <core/io/silent/SilentFileData.hh>

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/sequence/AnnotatedSequence.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL headers

// C++ Headers
#include <string>
#include <map>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>


namespace core {
namespace io {
namespace silent {


	//////////////////////////////////////////////////////////////////////
	// holds all the data for a silent-structure
	class SilentStruct : public utility::pointer::ReferenceCount
#ifdef PTR_MODERN
		// New version
		, public utility::pointer::enable_shared_from_this< SilentStruct >
{
#else
{
		// Old intrusive ref-counter version
		inline SilentStructCOP shared_from_this() const { return SilentStructCOP( this ); }
		inline SilentStructOP shared_from_this() { return SilentStructOP( this ); }
#endif

		typedef std::string string;

	public:
		//constructor
		SilentStruct();

		// destructor
		virtual ~SilentStruct();

		SilentStruct( SilentStruct const & );
		SilentStruct& operator= ( SilentStruct const & );

		virtual SilentStructOP clone() const = 0;

		/// self pointers
		inline SilentStructCOP get_self_ptr() const { return shared_from_this(); }
		inline SilentStructOP get_self_ptr() { return shared_from_this(); }
		//inline SilentStructCAP get_self_weak_ptr() const { return SilentStructCAP( shared_from_this() ); }
		//inline SilentStructAP get_self_weak_ptr() { return SilentStructAP( shared_from_this() ); }

		// @brief Fill a Pose with the conformation information in this
		// SilentStruct and the FA_STANDARD ResidueTypeSet.
		virtual void fill_pose(
			core::pose::Pose & pose
		) const;

		/// @brief Fill a Pose with the conformation information in this
		/// SilentStruct and the ResidueTypeSet provided by the caller. This is
		/// a virtual method which should be implemented by classes derived from
		/// SilentStruct.
		virtual void fill_pose(
			core::pose::Pose & pose,
			core::chemical::ResidueTypeSet const & residue_set
		) const;

		/// @brief Sets the tag from the Pose DataCache.
		void set_tag_from_pose(
			const core::pose::Pose & pose
		);

		void
		precision( core::Size precision );
		core::Size precision() const;

		void
		scoreline_prefix( std::string const & prefix );
		std::string scoreline_prefix() const;

		/// @brief opposite of fill_pose -- superclass provides
    /// functionality used by most SilentStruct types, and is
    /// optionally called at the beginning of the subclass's
    /// fill_struct method. As of this writing, used only by
    /// protein SilentStructs
		virtual void fill_struct(
			core::pose::Pose const & pose,
			std::string tag = "empty_tag" );

		/// @brief calls optH if command line requests optH.
		/// must be called by derived classes.
		void finish_pose(
			core::pose::Pose & pose
		) const;

		/// @brief Do some sort of comparison between the actual RMSD of this
		/// silent-struct and the cached coordinates. Used for SilentStruct
		/// objects that are rebuild from torsions or other reduced
		/// representations of data.
		virtual Real get_debug_rmsd() = 0;

		/// @brief print out a header line to the given ozstream. In a rosetta++
		/// silent-file, this contained the lines:
		/// SEQUENCE: <protein sequence>\nSCORE: <list of score-types>.
		virtual void print_header      ( std::ostream & out ) const;
		/// @brief only print SCORE: header line
		virtual void print_score_header( std::ostream & out ) const;
		/// @brief print out a SCORE line to the given ozstream.
		virtual void print_scores      ( std::ostream & out ) const;
		/// @brief print the conformation information in SilentStruct to out.
		virtual void print_conformation( std::ostream & out ) const = 0;
		/// @brief print the comments in this SilentStruct.
		virtual void print_comments    ( std::ostream & out ) const;
		/// @brief print the resnum in this SilentStruct, if filled.
		virtual void print_residue_numbers    ( std::ostream & out ) const;

		/// @brief returns the number of residues contained by this
		/// SilentStruct.
		virtual Size nres() const {
			return nres_;
		}

		/// @brief set the tag associate with this SilentStruct
		void set_decoy_tag( std::string const & tag ) {
			decoy_tag_ = tag;
			for ( Size n = 1; n <= other_struct_list_.size(); n++ ) other_struct_list_[n]->set_decoy_tag( tag );
		}

		/// @brief returns the tag associated with this SilentStruct.
		std::string decoy_tag() const {
			return decoy_tag_;
		}

		/// @brief returns the sequence associated with this SilentStruct.
		core::sequence::AnnotatedSequence const& sequence() const {
			return sequence_;
		}

		/// @brief returns the number of residues in this SilentStruct.
		void nres( Size nres ) {
			nres_ = nres;
		}

		/// @brief sets the tag associated with this SilentStruct.
		void decoy_tag( std::string const & tag ) {
			decoy_tag_ = tag;
		}

		/// @brief sets the sequence for this SilentStruct.
		void sequence( core::sequence::AnnotatedSequence const& sequence ) {
			sequence_ = sequence;
		}

		/// @brief sets the silent_energies for this SilentStruct.
		void silent_energies( utility::vector1< SilentEnergy > const & new_se ) {
			silent_energies_ = new_se;
		}

		/// @brief sort all the silent energies by their name.
		void sort_silent_scores();

		/// @brief returns true if this SilentStruct has an energy for the given
		/// scorename, returns false otherwise.
		bool has_energy( std::string const scorename ) const;

		/// @brief Returns the energy associated with the given scorename if this
		/// SilentStruct has an energy for that scorename. Otherwise returns 0.
		core::Real get_energy( std::string const & scorename ) const;

		/// @brief Returns the energy associated with the given scorename if this
		/// SilentStruct has an energy for that scorename. Otherwise returns 0.
		std::string const & get_string_value( std::string const & scorename ) const;

		/// @brief Returns the SilentEnergy associated with this scorename.
		SilentEnergy const & get_silent_energy( std::string const & scorename ) const;

		utility::vector1< SilentEnergy > get_silent_energies(){ return silent_energies_;}

		void set_valid_energies( utility::vector1< std::string > valid );

		/// @brief Clear all of the energies in the SilentStruct. Doesn't just
		/// zero the energies, it entirely removes all knowledge of all energies
		/// from this SilentStruct.
		virtual void clear_energies() {
			silent_energies_.clear();
		}

		/// @brief Create a new SilentStruct object from the provided set of
		/// lines. This abstract method should be overwritten by derived
		/// classes. Returns false if the init_from_lines routine encounters a
		/// problem with the lines provided.
		virtual bool init_from_lines(
			utility::vector1< std::string > const & lines,
			SilentFileData & container
		) = 0;

		/// @brief add a score of a given name and value to this SilentStruct.
		/// Takes an optional weight that defaults to 1.0.
		void add_energy( std::string scorename, Real value, Real weight = 1.0 );

		/// @brief add a non-floating point score of a given name and value to this
		/// SilentStruct.
		void add_string_value( std::string scorename, std::string const & value );

		/// @brief Copy the score information in the given SilentStruct into
		/// this SilentStruct.
		void copy_scores( const SilentStruct & src_ss );

		/// @brief add a named comment to this SilentStruct object.
		/// Similar to methods for playing with energies, but
		/// mapping is string => string rather than string => Real.
		void add_comment( std::string name, std::string value );
		bool has_comment( std::string const & name ) const;
		std::string get_comment( std::string const & name ) const;
		void comment_from_line( std::string const & line );
		void erase_comment( std::string const & name );
		void clear_comments();

		std::map< std::string, std::string > get_all_comments() const;
		//void set_all_comments( std::map< std::string, std::string > );

		void parse_energies(
			std::istream & input,
			utility::vector1< std::string > const & energy_names
		);

		/// @brief Initialize this SilentStruct's energies from the given Pose.
		/// This sets energies, energy weights, and the output widths for the
		/// energies.
		void energies_from_pose( core::pose::Pose const & pose );
		/// @brief Put the energy information from this SilentStruct into the
		/// pose. Energies that correspond to a ScoreType are put into the
		/// pose.energies().total_energies() EnergyMap, all other energies are
		/// put into the ARBITRARY_FLOAT_DATA map in the pose DataCache. Also
		/// sets the scorefxn_weights in the Energies object using the
		/// information from this SilentStruct.
		void energies_into_pose( core::pose::Pose & pose ) const;

		/// @brief returns the positions of the CA atoms in this
		/// ProteinSilentStruct.  Useful for RMS calculations.
		virtual ObjexxFCL::FArray2D< Real > get_CA_xyz() const = 0;

		/// @brief Returns the vector of SilentEnergy objects associated with
		/// this SilentStruct object.
		utility::vector1< SilentEnergy > energies() const {
			return silent_energies_;
		}

		/// @brief Returns the EnergyNames that this SilentStruct contains.
		EnergyNames energy_names() const;

		// @brief Renames energies from Rosetta++ to their appropriate mini
		// names. Useful for reading older silent-files.
		void rename_energies();

		bool read_sequence(  std::string const& line );
		void read_score_headers( std::string const& line, utility::vector1< std::string > & enames, SilentFileData& container );

		/// @brief strip [...] comment from seqeunce_ and return pure
		/// one-letter sequence
		std::string one_letter_sequence() const;

		virtual core::Size mem_footprint() const { return 0; }

		//By Parin Sripakdeevong (sripakpa@stanford.edu).
		void print_parent_remarks( std::ostream & out ) const;

		//By Parin Sripakdeevong (sripakpa@stanford.edu).
		std::string get_parent_remark( std::string const & name ) const;

		//By Parin Sripakdeevong (sripakpa@stanford.edu).
		bool has_parent_remark( std::string const & name ) const;

		//By Parin Sripakdeevong (sripakpa@stanford.edu).
		void add_parent_remark( std::string const name, std::string const value );

		//By Parin Sripakdeevong (sripakpa@stanford.edu).
		void get_parent_remark_from_line( std::string const line );

		void set_residue_numbers( utility::vector1< Size > const & residue_numbers ){ residue_numbers_ = residue_numbers;}
		void set_chains( utility::vector1< char > const & chains ){ chains_ = chains;}
		void set_full_model_parameters( core::pose::full_model_info::FullModelParametersCOP setting ){ full_model_parameters_ = setting; }
		core::pose::full_model_info::FullModelParametersCOP full_model_parameters() const{ return full_model_parameters_; }

		void fill_struct_with_residue_numbers( pose::Pose const & pose );

		void
		fill_other_struct_list( pose::Pose const & pose );

		void residue_numbers_into_pose( pose::Pose & pose ) const;

		void full_model_info_into_pose( pose::Pose & pose ) const;

		void
		figure_out_residue_numbers_from_line( std::istream & line_stream );

		utility::vector1< SilentStructOP > const & other_struct_list() const { return other_struct_list_; }

		utility::vector1< SilentStructOP > & nonconst_other_struct_list() { return other_struct_list_; }

		void add_other_struct( SilentStructOP silent_struct );

		/// @brief Sets whether conversion from big-endian to little-endian (or the converse) should be forced
		/// when a binary silent structure is initialized from lines.
		virtual void set_force_bitflip( bool const setting) { force_bitflip_ = setting; return; }

		/// @brief Gets whether conversion from big-endian to little-endian (or the converse) should be forced
		/// when a binary silent structure is initialized from lines.
		virtual bool force_bitflip() const { return force_bitflip_; }

	protected:

		///@ brief helper to detect fullatom input
		// performs logical operation on fullatom and not_defined:
		// if residue is patched we cannot tell -- no change to fullatom and not_defined
		// if residue is not patched and is protein we do the following:
		// if it has 7 or 8 atoms -- no change to fullatom and not_defined
		// otherwise: well_defined = true;
		// if this residue has > 8 atoms: fullatom = true
		// if this residue has <= 6 atoms: fullatom = false
		void detect_fullatom( core::Size pos, core::Size natoms, bool &fullatom, bool& well_defined );

    /// @brief add string serialization of all WriteableCacheableData in as comments.
    /// @param pose a pose containing a datacache with WriteableCacheableData.
    void extract_writeable_cacheable_data( core::pose::Pose const& pose );

	private:

		/// @brief If true, this forces conversion from big-endian to little-endian (or the converse) when a binary silent structure
		/// is initialized from lines (init_from_lines() function).
		bool force_bitflip_;


		bool strict_column_mode_;
		Size nres_;
		std::string decoy_tag_;
		core::sequence::AnnotatedSequence sequence_;

		std::map< std::string, std::string > parent_remarks_map_; //Similar to silent_comments_ but this doesn't get outputted back to the silent_file.


		// output-related data
		std::map< std::string, std::string > silent_comments_;
		utility::vector1< SilentEnergy > silent_energies_;
    utility::vector1< std::string > cache_remarks_;

		utility::vector1< Size > residue_numbers_; // can be derived from PDB info.
		utility::vector1< char > chains_; // can be derived from PDB info.
		utility::vector1< SilentStructOP > other_struct_list_;
		core::pose::full_model_info::FullModelParametersCOP full_model_parameters_;

	private:
		/// @brief Updates the "score" entry in the silent_energies.
		void update_score();

		// number of digits after decimal place in scorefile
		core::Size precision_;
		// prefix for the SCORE: lines. Usually "SCORE:"
		std::string scoreline_prefix_;
}; // class SilentStruct

} // namespace silent
} // namespace io
} // namespace core

#endif
