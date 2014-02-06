// #error Don't include; for documentation only.

/// @namespace protocols::kinematic_closure
///
/// @brief Implementation of the kinematic closure algorithm for sampling 
/// protein backbone conformations.
///
/// @details The most important classes in this namespace are KicMover, 
/// ClosureProblem, and ClosureSolution.  KicMover is meant for use in other 
/// protocols, while ClosureProblem and ClosureSolution drive that kinematic 
/// closure algorithm under the hood.
///
/// KicMover objects are capable of sampling different backbone conformations 
/// in a pose.  Exactly how these objects work can be heavily customized using 
/// the perturber and picker classes.  Some of these classes are very general, 
/// while others are specific to certain applications.  Feel free to create new 
/// sampler classes for your applications.  For example, if you decide to write 
/// custom pivot picking and perturbing algorithms that are meant to be used 
/// together, make a new sampler class that wraps around KicMover and uses 
/// your custom algorithms by default.
///
/// The ClosureProblem and ClosureSolution classes are used by the KicMover 
/// as follows.  First, a ClosureProblem object is created.  This object is 
/// used to pick the pivot residues, to sample the nonpivot torsions, and to 
/// find backbone conformations consistent with all that.  All possible 
/// solutions are found, and each solution is returned as ClosureSolution 
/// object.  Then KicMover then picks a solution, using a handful of filters, 
/// and applies it to the pose.
