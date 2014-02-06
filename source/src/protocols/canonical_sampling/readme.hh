// This header is for documentation purposes only; it doesn't contain any code.

/// @namespace protocols::canonical_sampling
///
/// @brief Algorithms for rigorously sampling native state ensembles.  See 
/// disclaimer in detailed description.
///
/// @details When I started using the algorithms in this namespace, I had a 
/// difficult time orienting myself due to the complete lack of documentation.  
/// I decided that it would be helpful to write documentation as part of the 
/// process of figuring out how the code works.  I mention all this to get to 
/// my point, which is a disclaimer:
/// 
/// I didn't necessarily understand how everything worked when I was writing  
/// this documentation.  As I went along, I certainly became more and more 
/// comfortable with the code.  I even successfully refactored parts of it.  
/// But some of the documentation is speculative, so take everything with a 
/// grain of salt.  And please fix any errors you find!  With that said, and 
/// without any further ado, here is the description of the 
/// protocols::canonical_sampling namespace:
///
/// The central class in this namespace is MetropolisHastingsMover.  Its 
/// responsibility is to run a Monte Carlo loop.  This class is a Mover, and 
/// like all movers it has an apply() method that does all the work.  It also 
/// provides a number of getters and setters so you can control exactly what 
/// the algorithm does.  In particular, there are three classes that most 
/// influence how the MetropolisHastingsMover behaves.  These are 
/// ThermoynamicMover, TemperatureController, and MonteCarlo.
///
/// Subclasses of ThermodynamicMover are responsible for making moves which 
/// obey detailed balance.  Some classes which fulfill this interface include 
/// BackrubMover, LoopClosureMover, ShearMover, SmallMover, SidechainMover, and 
/// SidechainMCMover.  Note that both sidechain and backbone moves are
/// available.  Related to ThermodynamicMover is ThermodynamicObserver.  
/// Observer subclasses are meant to observe the progress of the simulation.  
/// The most important of these derive from TrajectoryRecorder and are used to 
/// keep track of the poses being sampled.  You can add a mover to your 
/// simulation using MetropolisHastingsMover.add_mover() or any of the 
/// similarly named convenience methods, and you can add an observer using 
/// MetropolisHastingsMover.add_observer().
///
/// Subclasses of TemperatureController are responsible for setting the 
/// temperature during the simulation.  The most basic controller is 
/// FixedTemperatureControler, which simply maintains a constant temperature 
/// and is used by default.  More sophisticated temperature protocols allow for 
/// faster equilibration.  These include SimulatedTempering, ParallelTempering, 
/// and HamiltonianExchange.  The former two algorithms can both be used to 
/// sample a defined thermodynamic ensemble, which is a nice feature.
///
/// The actual simulation itself is delegated to the standard MonteCarlo class.  
/// The temperature and the score function used in the simulation are taken 
/// from a MonteCarlo object provided via 
/// MetropolisHastingsMover.set_monte_carlo().  To smartly keep track of 
/// acceptance rates during simulations where the temperature changes, setup 
/// your MonteCarlo object to use the MultiTempTrialCounter by calling 
/// MonteCarlo.set_counter().  This custom counter separately keeps track of 
/// moves for each temperature level, which is useful information.
