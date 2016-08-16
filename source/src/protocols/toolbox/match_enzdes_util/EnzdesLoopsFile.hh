// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/match_enzdes_util/EnzdesLoopsFile.hh
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2009


#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_EnzdesLoopsFile_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_EnzdesLoopsFile_hh


#include <protocols/toolbox/match_enzdes_util/EnzdesLoopsFile.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.fwd.hh>


//#include <core/conformation
#include <core/types.hh>

#include <utility/io/izstream.fwd.hh>

#include <string>

#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

//tiny helper class for EnzdesLoopInfo
class ResInteractions{

public:

	ResInteractions();

	virtual
	~ResInteractions();

	virtual
	bool
	read_data( utility::io::izstream & data );

	virtual
	void
	write_data() const;

	core::Size targ_res() const {
		return targ_res_; }

	utility::vector1< std::string > const &
	targ_atom_names() const{
		return targ_atom_names_; }

	utility::vector1< std::string > const &
	targ_base_atom_names() const{
		return targ_base_atom_names_; }

	utility::vector1< std::string > const &
	targ_base2_atom_names() const{
		return targ_base2_atom_names_; }

	utility::vector1< std::string > const &
	loopres_atom_names() const{
		return loopres_atom_names_; }

	utility::vector1< std::string > const &
	loopres_base_atom_names() const{
		return loopres_base_atom_names_; }

	utility::vector1< std::string > const &
	loopres_base2_atom_names() const{
		return loopres_base2_atom_names_; }


	core::Size
	num_interactions() const {
		return num_interactions_; }

	toolbox::match_enzdes_util::GeomSampleInfoCOP
	dis() const;

	toolbox::match_enzdes_util::GeomSampleInfoCOP
	loop_ang() const;

	toolbox::match_enzdes_util::GeomSampleInfoCOP
	targ_ang() const;

	toolbox::match_enzdes_util::GeomSampleInfoCOP
	loop_dih() const;

	toolbox::match_enzdes_util::GeomSampleInfoCOP
	targ_dih() const;

	toolbox::match_enzdes_util::GeomSampleInfoCOP
	lt_dih() const;


protected:

	bool
	process_input_line_tokens( utility::vector1< std::string > const & tokens );

	void
	set_targ_res( core::Size targ_res) {
		targ_res_ = targ_res; }

	void
	set_targ_atom_names( utility::vector1< std::string > const & t_atom_names ) {
		targ_atom_names_ = t_atom_names; }

	void
	set_loopres_atom_names( utility::vector1< std::string > const & l_atom_names) {
		loopres_atom_names_ = l_atom_names; }

private:
	core::Size targ_res_;
	utility::vector1< std::string > targ_atom_names_;
	utility::vector1< std::string > targ_base_atom_names_;
	utility::vector1< std::string > targ_base2_atom_names_;
	core::Size num_interactions_;

	toolbox::match_enzdes_util::GeomSampleInfoOP dis_, loop_ang_, targ_ang_, loop_dih_, targ_dih_, lt_dih_;

	utility::vector1< std::string > loopres_atom_names_;
	utility::vector1< std::string > loopres_base_atom_names_;
	utility::vector1< std::string > loopres_base2_atom_names_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


//tiny helper class for EnzdesLoopInfo
class CstResInteractions : public ResInteractions{

public:

	CstResInteractions();

	virtual
	~CstResInteractions(){}

	bool
	read_data( utility::io::izstream & data );

	void
	write_data() const;

	bool
	resA() const{
		return resA_; }

	core::Size
	cst_block() const {
		return cst_block_; }

private:
	//if this is true it's resA, otherwise resB
	bool resA_;

	core::Size cst_block_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


class EnzdesLoopInfo : public utility::pointer::ReferenceCount
{


public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~EnzdesLoopInfo();

	EnzdesLoopInfo();

	bool
	read_loops_file_block( utility::io::izstream & data );

	bool
	check_data_consistency( bool report = false ) const;

	bool
	ss_strings_specified() const {
		return ss_strings_.size() != 0; }

	bool
	pdb_numb() const {
		return pdb_numb_;}

	bool
	pose_numb() const {
		return pose_numb_;}

	core::Size
	start_pdb() const {
		return loop_start_pdb_; }

	core::Size
	stop_pdb() const {
		return loop_end_pdb_; }

	char
	start_pdb_chain() const {
		return loop_start_pdb_chain_; }

	char
	stop_pdb_chain() const {
		return loop_end_pdb_chain_; }

	core::Size
	start() const {
		return loop_start_; }

	core::Size
	stop() const {
		return loop_end_; }

	core::Size
	min_length() const {
		return min_length_; }

	core::Size
	max_length() const {
		return max_length_; }


	utility::vector1< std::string > const &
	ss_strings() const {
		return ss_strings_; }

	utility::vector1< ResInteractions > const &
	res_interactions() const {
		return res_interactions_; }

	utility::vector1< CstResInteractions > const &
	cst_interactions() const {
		return cstres_interactions_; }

	//still unimplemented
	bool
	preserve_buried_contacts() const{
		return preserve_buried_contacts_; }

	//still unimplemented
	bool
	contact_buried_problematic_res() const {
		return contact_buried_problematic_res_; }


protected:

	/// @brief generate all secondary structure strings that correspond to
	/// the given blueprint. A blueprint is a string containing a succession
	/// of the following substrings: 1 ss-char followed by '(', a number and ')'
	/// or 1 ss-char followed by '(', a number, '-', a second number, and ')',
	/// i.e. H(3-5)L(2-3)H(8)L(5-10)E(4-5)L(2-3)
	/// the function interprets this string to mean a helix of length between
	/// 3 and 5, followed by a loop of length 2 or 3, followed by helix of
	/// length 8, followed by a loop of length between 5 and 10, followed
	/// by a sheet of length 4 or 5, followed by a loop of length 2 or 3.
	/// the function then generates one string for every possible combination
	/// of secondary structure elements that correspond to this blueprint.
	/// in the above example, a total of 3*2*1*6*2*2 = 144 different secondary
	/// structure strings will be generated. for every one of these strings,
	/// the function checks whether its length is within min_length_ and max_length_,
	/// and if this is the case, the string gets saved in the ss_strings_ vector.
	void
	generate_ss_strings_from_blueprint( std::string const & ss_blueprint );

	//data
private:

	//some basic information about the loop
	core::Size loop_start_, loop_end_, loop_start_pdb_, loop_end_pdb_;
	char loop_start_pdb_chain_, loop_end_pdb_chain_;
	bool pose_numb_, pdb_numb_;

	core::Size min_length_, max_length_;


	//user specified secondary structure strings
	utility::vector1< std::string > ss_strings_;

	//unimplemented at the moment
	bool preserve_buried_contacts_;

	//unimplemented at the moment
	bool contact_buried_problematic_res_;

	utility::vector1< ResInteractions > res_interactions_;

	utility::vector1< CstResInteractions > cstres_interactions_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief class to process an enzdes loops file
class EnzdesLoopsFile : public utility::pointer::ReferenceCount
{

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~EnzdesLoopsFile();

	EnzdesLoopsFile();

	bool
	read_loops_file( std::string filename );

	EnzdesLoopInfoCOP
	loop_info( core::Size l ) const;

	bool
	file_read() const {
		return file_read_; }

	core::Size num_loops() const {
		return enzloops_.size(); }

	void
	clear();

private:

	bool file_read_;

	utility::vector1< EnzdesLoopInfoOP > enzloops_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class EnzdesLoopsFile


} //match_enzdes_util
} //toolbox
} //protocols


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_match_enzdes_util_EnzdesLoopsFile )
#endif // SERIALIZATION


#endif //
