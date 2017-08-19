import os, sys, json
import logging

import pyrosetta.rosetta as rosetta

_logger = logging.getLogger("rosetta")

def initialize_logging():
    """Initialize rosetta logging if not already initialized."""

    # Check if rosetta or root logger have had handlers initialized.
    # Preceeding calls to logging.basicConfig will initialize root handlers.

    # If not, add an output handler to logger 'rosetta' targeting stdout and
    # set log level to INFO.
    if (not logging.root.handlers) and (not _logger.handlers):
        _logger.addHandler(logging.StreamHandler(sys.stdout))
        _logger.setLevel(logging.INFO)


def _notebook_logging_callback(logstring):
    """Logging callback which pring Rosetta output on Jupyter notenook."""
    sys.stdout.write(logstring)


# Added because the new bindings have an issues setting the _logging_tracer.output_callback to a python function
# 

class CustomJupyterTracer( rosetta.basic.PyTracer ):
    def __init__(self):
        super( CustomJupyterTracer, self).__init__()
        self.set_ios_hook = rosetta.basic.Tracer.set_ios_hook
        
    def output_callback(self, logstring):
        return _notebook_logging_callback(logstring)

    def __del__(self):
        self.set_ios_hook(None,'')
            
def set_logging_handler(notebook):
    """Redirect all rosetta trace output through the logging.Logger 'rosetta'."""
    global _logging_tracer
    _logging_tracer = CustomJupyterTracer()

    # Set the tracer hook, do not enable 'raw' in order to
    # allow rosetta's tracing settings to have priority.
    rosetta.basic.Tracer.set_ios_hook(_logging_tracer, rosetta.basic.Tracer.get_all_channels_string(), False)

    # Mute rosetta's standard logging output
    rosetta.basic.Tracer.super_mute(True)
