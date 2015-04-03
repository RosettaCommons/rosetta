// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/canonical_sampling/HamiltonianExchangeMover.cc
/// @brief HamiltonianExchange methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/HamiltonianExchange.hh>
#include <protocols/canonical_sampling/HamiltonianExchangeCreator.hh>


// protocols headers
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <protocols/canonical_sampling/BiasedMonteCarlo.hh>
#include <protocols/canonical_sampling/BiasEnergy.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>

#include <protocols/rosetta_scripts/util.hh>

//#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>

// core headers
#include <core/scoring/ScoreFunction.hh>

#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/types.hh>

// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/file/file_sys_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ Headers
#include <utility/excn/Exceptions.hh>
#include <cmath>


using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer tr( "protocols.canonical_sampling.HamiltonianExchange" );


bool protocols::canonical_sampling::HamiltonianExchange::options_registered_( false );

//Mike: when you want to remove these Macros... leave them at least here as comment - since they provide documentation
void protocols::canonical_sampling::HamiltonianExchange::register_options() {
	if ( !options_registered_ ) {
		options_registered_ = true;
		Parent::register_options();
	}
}

namespace protocols {
namespace canonical_sampling {
using namespace core;

std::string
HamiltonianExchangeCreator::keyname() const {
	return HamiltonianExchangeCreator::mover_name();
}

protocols::moves::MoverOP
HamiltonianExchangeCreator::create_mover() const {
	return protocols::moves::MoverOP( new HamiltonianExchange );
}

std::string
HamiltonianExchangeCreator::mover_name() {
	return "HamiltonianExchange";
}

HamiltonianExchange::HamiltonianExchange() :
	exchange_grid_dimension_( 1 ),
	bias_energy_( /* NULL */ )
{
	set_defaults();
}

HamiltonianExchange::HamiltonianExchange(	HamiltonianExchange const & other ) :
	Parent( other ),
	hamiltonians_( other.hamiltonians_ ),
	exchange_schedules_( other.exchange_schedules_ ),
	current_exchange_schedule_( other.current_exchange_schedule_ ),
	exchange_grid_( other.exchange_grid_ ),
	exchange_grid_dimension_( other.exchange_grid_dimension_ ),
	successfully_initialized_( other.successfully_initialized_ ),
	max_coord_( other.max_coord_ ),
	bias_energy_( other.bias_energy_ )
{
	Size const nlevels( n_temp_levels() );
	runtime_assert( nlevels == hamiltonians_.size() );
}

HamiltonianExchange& HamiltonianExchange::operator=( HamiltonianExchange const& other ) {
	if ( &other == this ) return *this;
	Parent::operator=( other );
	hamiltonians_ = other.hamiltonians_;
	exchange_schedules_ = other.exchange_schedules_;
	current_exchange_schedule_ = other.current_exchange_schedule_;
	exchange_grid_ = other.exchange_grid_;
	exchange_grid_dimension_ = other.exchange_grid_dimension_;
	successfully_initialized_ = other.successfully_initialized_;
	max_coord_ = other.max_coord_;
	bias_energy_ = other.bias_energy_;
	Size const nlevels( n_temp_levels() );
	runtime_assert( nlevels == hamiltonians_.size() );
	return *this;
}

HamiltonianExchange::~HamiltonianExchange() {
}


std::string
HamiltonianExchange::get_name() const
{
	return "HamiltonianExchange";
}

protocols::moves::MoverOP
HamiltonianExchange::clone() const
{
	return protocols::moves::MoverOP( new protocols::canonical_sampling::HamiltonianExchange(*this) );
}

protocols::moves::MoverOP
HamiltonianExchange::fresh_instance() const
{
	return protocols::moves::MoverOP( new HamiltonianExchange );
}

void
HamiltonianExchange::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	pose::Pose const & pose
) {
	exchange_schedule_file_ = tag->getOption< std::string >( "exchange_schedule", "");

	Parent::parse_my_tag( tag, data, filters, movers, pose );
	if ( !successfully_initialized_ ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "Initialization of HamiltonianExchange Module failed! " );
	}
}


/// handling of options including command-line
void HamiltonianExchange::set_defaults() {
}

/// @brief Assigns user specified values to primitive members using command line options
void HamiltonianExchange::init_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	Parent::init_from_options();
}

void
HamiltonianExchange::set_monte_carlo(
  protocols::moves::MonteCarloOP monte_carlo
) {
	Parent::set_monte_carlo( monte_carlo );
	BiasedMonteCarloOP biased_mc = utility::pointer::dynamic_pointer_cast< protocols::canonical_sampling::BiasedMonteCarlo > ( monte_carlo );
	if ( biased_mc ) {
		bias_energy_ = biased_mc->bias_energy();
	}
}

void
HamiltonianExchange::initialize_simulation(
	 core::pose::Pose& pose,
		MetropolisHastingsMover const& metropolis_hastings_mover,
	 core::Size cycle //default=0; non-zero if trajectory is restarted
) {
	show( tr.Info );
	Parent::initialize_simulation( pose, metropolis_hastings_mover, cycle );
	//	Size const nlevels( n_temp_levels() );
	set_current_temp( rank()+1 );
	current_exchange_schedule_ = 0;
	next_exchange_schedule();
	monte_carlo()->reset_scorefxn( pose, *hamiltonians_[ rank()+1 ] );
}

void
HamiltonianExchange::initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & mhm,
		core::Size level,
		core::Real temp_in,
		core::Size cycle //default=0; non-zero if trajectory is restarted
) {
	show( tr.Info );
	Parent::initialize_simulation( pose, mhm, level, temp_in, cycle );
	current_exchange_schedule_ = 0;
	next_exchange_schedule();
	monte_carlo()->reset_scorefxn( pose, *hamiltonians_[ level ] );
}

Size
HamiltonianExchange::coord2key( GridCoord const& coord, GridCoord const& max_coord, Size exclude_dim ) {
	Size const n_dim( coord.size() );
	Size key( 0 );
	Size stride( 1 );
	for ( Size d=1; d<=n_dim; ++d ) {
		if ( d == exclude_dim ) continue;
		key += (coord[ d ] - 1) * stride;
		stride *= max_coord[d];
	}
	return key;
}

Size
HamiltonianExchange::nlevels_per_dim( core::Size dim ) const {
	if ( dim > exchange_grid_dimension_ ) {
		return 1;
	} else {
		return max_coord_[dim];
	}
}

void HamiltonianExchange::setup_exchange_schedule() {
	exchange_schedules_.clear();
	ExchangeSchedule list;

	//to compute unique keys we need the max coord per dimension
	int const n_dim( exchange_grid_dimension_ );
	GridCoord max_coord( n_dim, 0 );
	for ( Grid::const_iterator it=exchange_grid_.begin(); it!=exchange_grid_.end(); ++it ) {
		GridCoord const& coord( *it );
		for ( Size dim=1; (int)dim<=n_dim; ++dim ) {
			if ( coord[ dim ] > max_coord[ dim ] ) max_coord[ dim ] = coord[ dim ];
		}
	}
	max_coord_ = max_coord;

	if ( ! initialize_exchange_schedule_from_file( exchange_schedule_file_ ) ) { // if not initialize from file, do the next

	for ( Size dim=1; (int)dim<=n_dim; ++dim ) {
		//make list of uniqe groups: if dim==2 (1,1,1), (1,2,1), (1,3,1), ... (1, N, 1) should be one group
		typedef std::list< std::pair< core::Size, core::Size > > Level2GridpointList;
		typedef std::map< core::Size, Level2GridpointList > Groups;
		Groups groups;
		for ( Size i=1; i<=exchange_grid_.size(); ++i ) {
			GridCoord const& coord( exchange_grid_[ i ] );
			Size key( coord2key( coord, max_coord, dim ) );
			groups[ key ].push_back( std::make_pair( coord[ dim ], i ) );
		}

		//sort witin groups
		for ( Groups::iterator it = groups.begin(); it != groups.end(); ++it ) {
			it->second.sort();
		}

		for ( Size phase = 1; phase <= 2; ++phase ) {
			//PHASE 2 	2<->3, 4<->5...
			list.clear();
			for ( Groups::const_iterator git = groups.begin(); git != groups.end(); ++git ) {
				Level2GridpointList::const_iterator lit=git->second.begin();
				if ( phase==2 && git->second.size() > 2 ) ++lit;   //if more than 2 elements we need to switch phases
				for ( ; lit != git->second.end(); ++lit ) {
					Level2GridpointList::const_iterator first_lit = lit;
					Level2GridpointList::const_iterator second_lit = ++lit;
					if ( second_lit != git->second.end() ) {
						std::pair<int, int> elem( first_lit->second, second_lit->second );
						list.push_back(elem);
					} else break;
				}
			}
			exchange_schedules_.push_back(list);
		} // phase
	} // per dimension
	} // !initialize_exchange_schedule_from_file
}

Size
HamiltonianExchange::next_exchange_level() const {
	ExchangeSchedule const& ex( exchange_schedules_[ current_exchange_schedule_ ] );
	int self( current_temp() );
	int other( 0 );
	for ( Size i=1; i<=ex.size() && !other; i++ ) {
		if ( ex[ i ].first == self ) other = ex[ i ].second;
		if ( ex[ i ].second == self ) other = ex[ i ].first;
	}
	return other;
}


core::Real
HamiltonianExchange::temperature_move( core::Real ) {
	utility_exit_with_message( "HamiltonianExchange::temperature_move() called without pose... HamEx requires pose \
  to evaluate alternative energy function prior to switching..." );
	return -1;
}

core::Real
HamiltonianExchange::temperature_move( pose::Pose& MPI_ONLY( pose ) ) {
	using namespace ObjexxFCL::format;
	check_temp_consistency();
	//	tr.Trace << "counter: "<< temp_trial_count_+1 << " " << temperature_stride_ << " " << temp_trial_count_ % temperature_stride_ << std::endl;
	if ( !time_for_temp_move() ) return temperature();
	tr.Trace << "time for temp-move! " << std::endl;

#ifdef USEMPI
	//communication tags
	int const mpi_LEVEL_INFORM = 1;
	int const mpi_SCORE_INFORM = 2;
	int const mpi_LEVEL_DECISION = 3;

	//Size const nlevels( n_temp_levels() );
	int exchange_partner;
	bool is_master;
	find_exchange_partner( exchange_partner, is_master );
	int my_level( temperature_level() );
	int new_level;

	MPI_Request reqs[2];
	MPI_Status stats[2];
	//SEND UP ( alternative levels )
	if ( exchange_partner < 0 ) return temperature();

	Real new_bias_score( 0 );
	Real bias_score( 0 );
	if ( bias_energy_ ) {
		bias_score = bias_energy_->evaluate( pose );
		new_bias_score = bias_energy_->update_and_evaluate_replica( exchange_partner, mpi_comm(), pose );
		tr.Trace << "bias scores: " << bias_score << " " << new_bias_score << std::endl;
	}

	tr.Trace << "exchange partner: " << exchange_partner << std::endl;
	MPI_Isend( &my_level, 1, MPI_INT, exchange_partner, mpi_LEVEL_INFORM, mpi_comm(), &reqs[0]);
	MPI_Irecv( &new_level, 1, MPI_INT, exchange_partner, mpi_LEVEL_INFORM, mpi_comm(), &reqs[1]);
	MPI_Waitall(2, reqs, stats);
	tr.Trace<< "my_level: " << my_level << " other level: " << new_level << std::endl;

	std::ostringstream trial_str;
	trial_str << "HX_" << my_level << "_" << new_level;
	trial_counter().count_trial( trial_str.str() );

	Real const score( hamiltonians_[ my_level ]->score( pose ) );
	Real const new_score( hamiltonians_[ new_level ]->score( pose ) );

	float scores[ 2 ];
	int swap( 0 );
	if ( !is_master ) {
		scores[ 0 ] = new_score + new_bias_score; //the lower is other for the upper partner
		scores[ 1 ] = score + bias_score; //the higher is self for the upper partner
		tr.Trace << " send scores " << F( 8,5, scores[ 1 ]) << " " << F( 8,5, scores[ 0 ]) << std::endl;
		MPI_Send( &scores, 2, MPI_FLOAT, exchange_partner, mpi_SCORE_INFORM, mpi_comm() );
		MPI_Recv( &swap, 1, MPI_INT, exchange_partner, mpi_LEVEL_DECISION, mpi_comm(), &stats[0] );
	} else {
		MPI_Recv( &scores, 2, MPI_FLOAT, exchange_partner, mpi_SCORE_INFORM, mpi_comm(), &stats[0] );
		//x1 : our process
		//x2 : other process
		// V0 : our scorefxn, V1 : other scorefxn
		// V0(x2)-V0(x1) --> scores[1] - score
		// V1(x1)-V1(x2) --> new_score - scores[0]
		//tr.Trace << "level: " << my_level << " other level: " << new_level << " scores: "
		//<< F( 8,5, score ) << " " << F( 8,5, new_score ) << " "
		//<< F( 8,5, scores[ 0 ]) << " " << F( 8,5, scores[ 1 ]) << std::endl;
		Real const invT1( 1.0 / temperature() );
		Real const invT2( 1.0 / temperature( new_level ) );
		Real const deltaE1( scores[ 0 ] - score - bias_score );
		Real const deltaE2( scores[ 1 ] - new_score - new_bias_score );
		Real const delta( invT1*deltaE1 - invT2*deltaE2 );
		Real const r( numeric::random::rg().uniform() );
		swap = r < std::min( 1.0, std::exp( std::max(-40.0, -delta) ) ) ? 1 : 0;
		MPI_Send( &swap, 1, MPI_INT, exchange_partner, mpi_LEVEL_DECISION, mpi_comm() );
		//tr.Trace << "decision: "<< F( 4,2, invT1) << " " << F( 4,2,invT2) << " "<< F( 4,2, deltaE1 )
		//			 << " " << F( 4,2, deltaE2 ) << " " << F( 4,2, delta ) << " "
		//			 << F( 4,2, std::exp( std::max(-40.0, -delta ) ) ) <<  " " << F( 4,2, r ) << std::endl;

	}
	if ( is_master ) tr.Debug << "Exchange between level " << my_level << " and " << new_level
														<< " at ranks " << rank() << ":" << exchange_partner
														<< ( swap ? " swapped! " : " not swapped" ) << std::endl;
	if ( swap ) {
		tr.Trace<< "swap! " << std::endl;
		trial_counter().count_accepted( trial_str.str() );
		if ( bias_energy_ ) bias_energy_->swap_replicas();
		set_current_temp( new_level );
		monte_carlo()->score_function( *hamiltonians_[ new_level ] );
	}
#endif
	next_exchange_schedule();
	return temperature();
}

void HamiltonianExchange::next_exchange_schedule() {
	do {
		current_exchange_schedule_ = ( current_exchange_schedule_ + 1 ) % exchange_schedules_.size();
	}	while( !next_exchange_level() );
}

void HamiltonianExchange::clear() {
	exchange_schedules_.clear();
	exchange_grid_.clear();
	hamiltonians_.clear();
	max_coord_.clear();
	Parent::clear();
}

using namespace core::scoring;
class PatchOperation {
public:

	PatchOperation( ScoreType st, std::string const& op, Real wt ) :
		score_type_( st ),
		op_ ( op ),
		wt_ ( wt ),
		is_file_ (false)
	{}

	PatchOperation( std::string const& file ) : file_( file ) {
		is_file_ = true; //cheap Polymorphism
	}

	void apply( ScoreFunction& score ) const {
		if ( is_file_ ) {
			score.apply_patch_from_file( file_ );
		} else { //direct operation
			if ( op_ == "*=" ) {
				tr.Debug << "patching weight " << score_type_ << " with " << op_ << wt_ << std::endl;
				score.set_weight( score_type_, score.get_weight( score_type_ )*wt_ );
			} else if ( op_ == "=" ) {
				score.set_weight( score_type_, wt_ );
			} else {
				utility_exit_with_message("unrecognized scorefunction patch operation "+op_	);
			}
		} //direct operation
	} //apply

private:
	ScoreType score_type_;
	std::string op_;
	Real wt_;
	std::string file_;
	bool is_file_;
};

///brief: when an exchange schedule is given in a file, initialize from the file, instead of generate it by setup_exchange_schedule()
///each exchange pairs are specified with their temp_level, all pairs in one phase make one line in the input file
bool HamiltonianExchange::initialize_exchange_schedule_from_file( std::string const& filename ) {
	exchange_schedules_.clear();
	ExchangeSchedule list;
	utility::io::izstream in( filename );
	if ( !filename.size() ) {
		tr.Error << "no input file" << std::endl;
		return false;
	}
	if ( !in.good() ) {
		tr.Error << "cannot open file " << filename << ", will generate exchange schedule with setup_exchange_schedule()!"<<std::endl;
		return false;
	}
	std::string line;
	tr.Info << "reading exchange schedule file " << "..." << std::endl;
	while ( getline( in, line ) ) {
		list.clear();
		std::istringstream line_stream( line );
		Size grid1;
		Size grid2;
		line_stream >> grid1 >> grid2;
		std::pair<int, int> elem( grid1, grid2);
		list.push_back( elem );
		while ( line_stream.good() ) {
			line_stream >> grid1 >> grid2;
			if ( !line_stream.fail() ) {
				std::pair< int, int> elem( grid1, grid2);
				list.push_back( elem );
			}
		}
		exchange_schedules_.push_back( list );
	}
	if ( exchange_schedules_.size() ) {
		return true;
	} else { return false; }
}

bool HamiltonianExchange::initialize_from_file( std::string const& filename ) {
	typedef utility::vector1< PatchOperation > PatchOperationList;
	PatchOperationList global_patch_operations;

	clear();
	utility::vector1< core::Real > temperatures;
	utility::io::izstream in( filename );
	if ( !in.good() ) {
		tr.Error << "cannot open file " << filename << std::endl;
		utility_exit_with_message( "cannot open file "+filename );
	}

	std::string line;

	Size n_dim( 1 );
	exchange_grid_dimension_ = 0;
	// table format
	tr.Info << " reading Hamiltonian configuration file " << filename << "..." << std::endl;
	while ( getline( in, line ) ) {
		std::istringstream tag_stream( line );
		std::string tag;
		tag_stream >> tag;
		tr << "tag: " << tag << std::endl;   //zhe
		if ( tag_stream.good() && tag == "GRID_DIM" ) {
			if ( exchange_grid_dimension_ ) {
				utility_exit_with_message( "GRID_DIM statement has to come before any score-function definitions" );
			}
			tag_stream >> n_dim;
			if ( tag_stream.fail() ) {
				utility_exit_with_message( "require integer after tag GRID_DIM");
			}
			tr.Info << "using " << n_dim << " grid-dimensions" << std::endl;
			continue;
		}
		if ( tag_stream.good() && tag == "GLOBAL_PATCH" ) {
			std::string tag1, tag2;
			tag_stream >> tag1;
			if ( tag_stream.fail() ) {
				utility_exit_with_message( "GLOBAL_PATCH statement has to be followed by file-name or patch-operation" );
			}
			tag_stream >> tag2;
			if ( tag_stream.fail() ) { //this was a single file-name
				global_patch_operations.push_back( PatchOperation( tag1 ) );
			} else {
				std::istringstream line_stream( line );
				line_stream >> tag; //read GLOBAL_PATCH again
				core::scoring::ScoreType score_type;
				std::string operation;
				Real wt;
				line_stream >> score_type >> operation >> wt;
				if ( line_stream.fail() ) {
					utility_exit_with_message( "Expected score_type, operation and weight after GLOBAL_PATCH or a file-name" );
				}
				global_patch_operations.push_back( PatchOperation( score_type, operation, wt ) );
			}
			continue;
		}
		if ( tag_stream.good() && tag[0]=='#' ) {
			tr.Debug << "read comment line: " << line << std::endl;
			continue;
		}

		if ( !tag_stream.good() ) {
			tr.Debug << "read empy line: " << line << std::endl;
			continue;
		}

		// first proper line ( not GRID_DIM or comment )
		exchange_grid_dimension_ = n_dim;
		std::istringstream line_stream( line+" " ); //start over reading the same line
		Real temp;
		std::string score_name;
		std::string patch_name;
		GridCoord coord( n_dim );
		for ( Size d = 1; d<=n_dim; ++d ) {
			line_stream >> coord[ d ];
		}
		if ( !line_stream.good() ) {
			tr.Error << "format error in hamiltonian exchange file : " << filename << " at line " << line << std::endl;
			tr.Error << "expected " << n_dim << " integer values for the grid-coordinate" << std::endl;
			utility_exit_with_message( "format error in "+filename );
			return false;
		}
		line_stream >> temp >> score_name;
		if ( !line_stream.good() ) {
			tr.Error << "format error in hamiltonian exchange file : " << filename << " at line " << line << std::endl;
			tr.Error << "expected real value for temperature and string for score-name" << std::endl;
			utility_exit_with_message( "format error in "+filename );
			return false;
		}
		line_stream >> patch_name;
		using namespace core::scoring;
		ScoreFunctionOP score;
		if ( !line_stream.good() ) { // no patch
			score = ScoreFunctionFactory::create_score_function( score_name );
		} else { // patch or patch = "NO_PATCH"
			score = ScoreFunctionFactory::create_score_function( score_name, patch_name );
		}
		for ( PatchOperationList::const_iterator it = global_patch_operations.begin(); it != global_patch_operations.end(); ++it ) {
			it->apply( *score );
		}
		tr.Debug << "line_stream still good after PATCH reading" << std::endl;
		while ( line_stream.good() ) {
			std::string tag;
			core::scoring::ScoreType score_type;
			std::string operation;
			Real wt;
			line_stream >> tag;
			if ( line_stream.fail() ) {  //nothing else to read
				break;
			}
			if ( tag == "ETABLE" ) {
				line_stream >> tag;
				score->set_etable( tag );
				continue;
			} else {
				std::istringstream tag_stream( tag );
				tag_stream >> score_type;
				line_stream >> operation >> wt;
			}
			if ( line_stream.fail() ) {
				tr.Debug << "tried to read a X op wt triple and failed.. " << std::endl;
				break;
			}
			PatchOperation patch( score_type, operation, wt );
			patch.apply( *score );
		}

		tr.Debug << "adding gridd cell with temperature: " << temp << " and score " << score << " at grid-coord: ";
		for ( GridCoord::const_iterator it=coord.begin(); it!=coord.end(); ++it ) {
			tr.Debug << *it << " ";
		}
		tr.Debug << std::endl;

		temperatures.push_back( temp );
		hamiltonians_.push_back( score );
		exchange_grid_.push_back( coord );
	}
	set_temperatures( temperatures );
	runtime_assert( temperatures.size() == n_temp_levels() );
	runtime_assert( n_temp_levels() == hamiltonians_.size() );
	runtime_assert( n_temp_levels() == exchange_grid_.size() );

	setup_exchange_schedule();
	successfully_initialized_ = true;
	return true; //succesfully initialized
}


void HamiltonianExchange::show( std::ostream& os ) const {
	using namespace ObjexxFCL::format;
	// All osput will be 80 characters - 80 is a nice number, don't you think?
	std::string line_marker = "///";
	os << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	os << line_marker << A( 47, "HamiltonianExchange Module" ) << space( 27 ) << line_marker << std::endl;
	os << line_marker << space( 74 ) << line_marker << std::endl;
	// Display the movable jumps that will be used in docking
	os << line_marker << A( 20, "Hamiltonian Cells: " ) << I( 5, n_temp_levels() )
		 << A( 40, "Hamiltonian Grid Dimension: " ) << I( 5, exchange_grid_dimension_ ) << space( 4 ) << line_marker << std::endl;
	os << line_marker << repeat( 74, '-' ) << line_marker << std::endl;
	for ( Size level=1; level<= n_temp_levels(); ++level ) {
		os << line_marker << A( 20, "Grid Cell: " );
		for ( Size d=1; d<=exchange_grid_dimension_; d++ ) {
			os << I( 2, exchange_grid_[ level ][ d ] );
		}
		os << A( 15, " Temperature: " ) << F( 5, 3, temperature( level ) )
			 << space( 74-20-20-3*exchange_grid_dimension_ ) << line_marker << std::endl;
		if ( tr.Debug.visible() ) hamiltonians_[ level ]->show( os );
		os << std::endl;
		os << line_marker << repeat( 74, '-' ) << line_marker << std::endl;
	}

	os << line_marker << repeat( 74, '=' ) << line_marker << std::endl;
	os << line_marker << A( 40, "Exchange Schedules " ) << std::endl;
	for ( utility::vector0< ExchangeSchedule >::const_iterator ex_it = exchange_schedules_.begin(); ex_it != exchange_schedules_.end(); ++ex_it ) {
		ExchangeSchedule const& ex( *ex_it );
		os << line_marker << repeat( 74, '-' ) << line_marker << std::endl;
		for ( ExchangeSchedule::const_iterator it = ex.begin(); it != ex.end(); ++it ) {
			GridCoord const& coord1( exchange_grid_[ it->first ] );
			GridCoord const& coord2( exchange_grid_[ it->second ] );
			Size const n_dim( exchange_grid_dimension_ );
			os << line_marker << A( 20, "( ");
			for ( Size d=1; d<=n_dim; ++d ) {
				os << I( 2, coord1[ d ] ) << (d != n_dim ? ", " : " )  " );
			}
			os << A( 10, " <-->   ( ");
			for ( Size d=1; d<=n_dim; ++d ) {
				os << I( 2, coord2[ d ] ) << (d != n_dim ? ", " : " )  " );
			}
			os << space( 74-( 20+2+4*n_dim+2+10+8+4*n_dim+2 ) ) << line_marker << std::endl;
		}
	}
	// Close the box I have drawn
	os << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
}

std::ostream& operator << ( std::ostream & os, HamiltonianExchange const& obj ) {
	obj.show( os );
	return os;
}

} //moves
} //protocols

