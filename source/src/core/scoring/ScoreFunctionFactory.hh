// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/conformation/ResidueFactory.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_ScoreFunctionFactory_hh
#define INCLUDED_core_scoring_ScoreFunctionFactory_hh

#include <core/scoring/ScoreFunctionFactory.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>
#include <utility/SingletonBase.hh>
#include <string>
#include <map>
#include <tuple>

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

#ifdef WIN32 //VC++ needs full class declaration
#include <core/scoring/ScoreFunction.hh> // WIN32 INCLUDE
#endif
// #include <core/chemical/ResidueType.fwd.hh>
// #include <core/conformation/Conformation.fwd.hh>
// #include <core/conformation/Residue.fwd.hh>
// #include <core/chemical/AtomTypeSet.fwd.hh>
// #include <core/chemical/MMAtomTypeSet.fwd.hh>

namespace core {
namespace scoring {

/// @brief A key for looking up previously-loaded scorefunctions (which the ScoreFunctionFactory stores
/// in a map of owning pointers indexed by keys of this type).
/// @details The key is the weights file name as a string, a comma-separated list of patches as another
/// string, and a raw pointer to the options collection.  Note that a raw pointer is appropriate in this
/// case since the goal is NOT to access the options collection, but simply to have a unique identifier
/// for whether we've seen it before.  An owning pointer or even a weak pointer (access pointer) would
/// incur unncessary overhead that we don't need for that goal.  If an OptionCollection could be hashed,
/// that might be an even better unique identifier, but for now, we'll use the memory location.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
struct ScoreFunctionKey {

public:

	/// @brief Default constructor -- explicitly deleted.
	ScoreFunctionKey() = delete;

	/// @brief Options constructor -- used to set contents of key.
	ScoreFunctionKey(
		std::string const & weights_file,
		std::string const & comma_separated_patches,
		utility::options::OptionCollection const * const options_collection_ptr
	) :
		keytuple_( std::make_tuple( weights_file, comma_separated_patches, options_collection_ptr ) )
	{}

	~ScoreFunctionKey() = default;

	/// @brief Comparison operator is needed to use this as the key for a map.
	/// We'll just use a thin wrapper for std::tuple's less than operator.
	bool
	operator< (
		ScoreFunctionKey const & other
	) const {
		return keytuple_ < other.keytuple_;
	}

private:

	/// @brief Internally, we still use a tuple since tuple::operator< is already defined.
	std::tuple < std::string, std::string, utility::options::OptionCollection const * > keytuple_;

};

/// @brief A static singleton for making a single score_function.
class ScoreFunctionFactory : public utility::SingletonBase< ScoreFunctionFactory >
{
public:
	friend class utility::SingletonBase< ScoreFunctionFactory >;

private:

	//private constructor
	ScoreFunctionFactory() = default;
	~ScoreFunctionFactory() = default;

public:

	/// @brief Returns a ScoreFunction from the database weights file  <weights_tag>
	///
	/// example(s):
	///     scorefxn = create_score_function('standard')
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.show
	///     ScoreFunction.weights
	///     ScoreType
	static
	ScoreFunctionOP
	create_score_function( std::string const & weights_tag );

	static
	ScoreFunctionOP
	create_score_function( utility::options::OptionCollection const & options, std::string const & weights_tag );

	/// @brief Returns a ScoreFunction from the database weights file  <weights_tag>
	/// with the patch <patch_tag>
	///
	/// example(s):
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.show
	///     ScoreFunction.weights
	///     ScoreType
	static
	ScoreFunctionOP
	create_score_function( std::string const & weights_tag, std::string const & patch_tag );

	static
	ScoreFunctionOP
	create_score_function( utility::options::OptionCollection const & options, std::string const & weights_tag, std::string const & patch_tag );

	/// @brief Returns a ScoreFunction from the database weights file  <weights_tag>  with patches in <patch_tags>
	static
	ScoreFunctionOP
	create_score_function( std::string const & weights_tag, utility::vector1< std::string > const & patch_tags );


	static
	ScoreFunctionOP
	create_score_function( utility::options::OptionCollection const & options, std::string const & weights_tag, utility::vector1< std::string > const & patch_tags );

	/// @brief A documentation function which reports the set of options read by the create_score_function variants
	static
	void
	list_read_options( utility::options::OptionKeyList & opts );

	/// @brief checks if the weights file is probably a talaris weights file
	/// and if it is consistent with the options system.
	/// static because C++ says it has to be; public because unit test
	/// spiritually const but you can't do that with static
	static
	bool
	validate_talaris(
		std::string const & weights_tag,
		utility::options::OptionCollection const & options );

	/// @brief checks if the weights file is probably a beta_15 weights file
	/// and if it is consistent with the options system.
	/// static because C++ says it has to be; public because unit test
	/// spiritually const but you can't do that with static
	static
	bool
	validate_beta(
		std::string const & weights_tag,
		utility::options::OptionCollection const & options ) /*const*/;

private:

	/// @brief Nonstatic version for storing the loaded scorefunction in this object.
	/// @details Loads scorefunction from disk if it has not yet been loaded (and caches it), or retrieves cached copy
	/// if it has been loaded already.  Returns clone of cached copy.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	ScoreFunctionOP
	create_score_function_nonstatic( utility::options::OptionCollection const & options, std::string const & weights_tag, utility::vector1< std::string > const & patch_tags );

	/// @brief Load a scorefunction from disk.
	/// @details TRIGGERS READ FROM DISK!  This function is needed for threadsafe creation.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	static
	core::scoring::ScoreFunctionOP
	load_score_function_from_disk(
		utility::options::OptionCollection const & options,
		std::string const & weights_tag_in,
		utility::vector1< std::string > const & patch_tags_in
	);

	/// @brief Applies user defined re-weighting from the options system. Reweights are applied as a
	/// factor of the original, so -rg_reweight 0.5 would result in half of the previously defined
	/// rg weight.
	static void apply_user_defined_reweighting_( utility::options::OptionCollection const & options, core::scoring::ScoreFunctionOP scorefxn );

	static void load_weights_file( std::string weights_tag, ScoreFunctionOP scorefxn );

private:

#ifdef MULTI_THREADED
	/// @brief A mutex for the loaded_scorefxns_ map.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	mutable utility::thread::ReadWriteMutex loaded_scorefxns_mutex_;
#endif

	/// @brief A cache of all the scorefunctions we've loaded so far, by filename, so that we don't have to
	/// load them from disk again.
	/// @details This is indexed by weights tag, comma-separated patches, pointer to options collection.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	mutable std::map < ScoreFunctionKey, core::scoring::ScoreFunctionOP > loaded_scorefxns_;

};

extern std::string const REF_2015;
extern std::string const TALARIS_2014;
extern std::string const TALARIS_2013;
extern std::string const TALARIS_2013_CART;
extern std::string const PRE_TALARIS_2013_STANDARD_WTS;
extern std::string const SCORE13;
extern std::string const CENTROID_WTS;
extern std::string const SOFT_REP_WTS;
extern std::string const SOFT_REP_DESIGN_WTS;
extern std::string const DNA_INT_WTS;
extern std::string const DNA_INT_WTS_GB;
extern std::string const MM_STD_WTS;
extern std::string const RNA_LORES_WTS;
extern std::string const RNA_HIRES_WTS;
extern std::string const RNA_LORES_PLUS_HIRES_WTS;
extern std::string const MEMB_HIGHRES_WTS; //pba

extern std::string const SCORE12_PATCH;
extern std::string const DOCK_PATCH;
extern std::string const DOCK_LOW_PATCH;

extern std::string const SCORE4_SMOOTH_CART;

extern std::string const BETA_GENPOT;
extern std::string const BETA_NOV16;
extern std::string const BETA_NOV15;
extern std::string const BETA_JULY15;

/// @brief A helper function which returns a scoring function held in an owning pointer according to the
/// user's command line parameters -score:weights and -score:patch
/// By default it returns weights=talaris2013 for fullatom,
/// and weights=cen_std and patch="" for centroid
core::scoring::ScoreFunctionOP get_score_function( bool const is_fullatom = true );

///@brief Get a ScoreFunction from cmd-line settings and dependant on the symmetrical state of the pose.
core::scoring::ScoreFunctionOP get_score_function( core::pose::Pose const & pose, bool const is_fullatom = true );

///@brief Get a ScoreFunction from cmd-line settings and dependant on the symmetrical state of the pose.
/// Local Options collection.
core::scoring::ScoreFunctionOP get_score_function( core::pose::Pose const & pose, utility::options::OptionCollection const & options, bool const is_fullatom = true );


/// @brief A helper function which creates a scoring function held in an owning pointer reading
/// from the input OptionCollection
core::scoring::ScoreFunctionOP get_score_function( utility::options::OptionCollection const & options, bool const is_fullatom = true );

/// @brief A documentation function which reports the set of options read by get_score_function.
void
list_read_options_in_get_score_function( utility::options::OptionKeyList & opts );

/// @brief A helper function that either returns a ScoreFunctionOP created by get_score_function() or
/// the one specified by the protocol which is activated by the -restore_pre_talaris_2013_behavior
/// flag.  The purpose of this function is to preserve legacy behavior for the sake of reproducibility
/// and so that a record of the old behavior is still preserved in the code to ease the process of
/// reverting the change to get_score_function if that were the wrong behavior.
core::scoring::ScoreFunctionOP get_score_function_legacy(
	std::string const & pre_talaris_2013_weight_set,
	std::string const & pre_talaris_2013_patch_file = ""
);


core::scoring::ScoreFunctionOP get_score_function_legacy(
	utility::options::OptionCollection const & options,
	std::string const & pre_talaris_2013_weight_set,
	std::string const & pre_talaris_2013_patch_file = ""
);


/// @brief A documentation function which reports the set of options read by get_score_function_legacy.
void
list_read_options_in_get_score_function_legacy( utility::options::OptionKeyList & opts );

/// @brief use the logic of get_score_function to get the name.
/// The  name format is <weights_tag>[_<patch_tag> ... ]
std::string
get_score_functionName(
	bool const is_fullatom = true );

inline
std::string
get_current_default_score_function_name(){
	return REF_2015;
}

/// @brief returns family name for a specific score function.
/// For example, ref2015_cart returns ref2015 and beta_nov16_cst returns beta_nov16
/// Returns an empty string if no match is found
std::string
basename_for_score_function( std::string const & sfxn_name );

void
apply_set_weights( ScoreFunctionOP scorefxn, utility::vector1< std::string > const & settings );

} // namespace scoring
} // namespace core

#endif
