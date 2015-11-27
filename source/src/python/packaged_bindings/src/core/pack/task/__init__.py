from __task_all_at_once_ import *

def standard_task_factory():
    """Create default task factory."""
    from rosetta.core.pack.task.operation import InitializeFromCommandline, NoRepackDisulfides

    tf = TaskFactory()
    tf.push_back(InitializeFromCommandline())
    tf.push_back(NoRepackDisulfides())
    return tf

def standard_packer_task(pose):
    """Create default packer task over given pose."""
    tf = standard_task_factory()
    task = tf.create_task_and_apply_taskoperations(pose)
    return task
