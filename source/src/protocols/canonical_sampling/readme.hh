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
/// All of the documentation in the module was written by someone who didn't 
/// really understand how the code works.  I was more just taking notes on how 
/// I felt particular pieces of code should work.  The idea was to draw 
/// attention to pieces of code that didn't work like I expected.  Therefore, 
/// even though the documentation is phrased in an authoritative way, it is 
/// based primarily on my feelings and secondarily on the code.
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
/// available.  To add a particular move to your simulation, use 
/// MetropolisHastingsMover.add_mover() or any of the similarly-named 
/// convenience methods.
/// 
/// Related to the ThermodynamicMover class is the ThermodynamicObserver class.  
/// Observer subclasses are meant to observe and change the parameters of this 
/// simulation itself.  There are three major subclasses: TrialCounterObserver, 
/// TrajectoryRecorder, and TemperatureController.  The trial counters keep 
/// track of move statistics and the trajectory recorders keep track of the 
/// poses being sampled.  TemperatureController is an especially important 
/// class because it is setup to manage the temperature of the Metropolis 
/// Hastings simulation.  Protocols like simulated annealing and parallel 
/// tempering are implemented by children of this class.  You can add an 
/// observer to your simulation using MetropolisHastingsMover.add_observer().
///
/// The actual simulation itself is delegated to the standard MonteCarlo class.  
/// You can use a custom MonteCarlo object in a canonical simulation by calling 
/// MetropolisHastingsMover.set_monte_carlo().
