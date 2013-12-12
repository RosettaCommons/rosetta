// #error Don't include; for documentation only.

/// @namespace protocols::kinematic_closure
///
/// @brief Implementation of the kinematic closure algorithm for sampling 
/// protein backbone conformations.
///
/// @details The most important classes in this namespace are 
/// samplers::KicSampler, ClosureProblem, and ClosureSolution.  KicSampler is 
/// meant for use in other protocols, while ClosureProblem and ClosureSolution 
/// drive that kinematic closure algorithm under the hood.
///
/// KicSampler objects are capable of sampling different backbone conformations 
/// in a pose.  Exactly how these objects work can be heavily customized using 
/// the perturber and picker classes.  Some of these classes are very general, 
/// while others are specific to certain applications.
///
/// The ClosureProblem and ClosureSolution classes are used by the KicSampler 
/// as follows.  First, a ClosureProblem object is created.  This object is 
/// used to pick the pivot residues, to sample the nonpivot torsions, and to 
/// find backbone conformations consistent with all that.  All possible 
/// solutions are found, and each solution is returned as ClosureSolution 
/// object.  Then KicSampler then picks a solution, using a handful of filters, 
/// and applies it to the pose.
