# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   TestBindings.py
## @brief  Run bindings test and demo scrips
## @author Sergey Lyskov


import os, sys, commands, datetime


def execute(message, commandline, return_=False, untilSuccesses=False):
    print message, commandline
    while True:
        (res, output) = commands.getstatusoutput(commandline)
        print output

        if res and untilSuccesses: pass  # Thats right - redability COUNT!
        else: break

        print "Error while executing %s: %s\n" % (message, output)
        print "Sleeping 60s... then I will retry..."
        time.sleep(60)

    if res:
        print "\nEncounter error while executing: " + commandline
        if return_==True: return True
        else: sys.exit(1)

    if return_ == 'output': return output
    else: return False



def main(args):


    #tests = execute('getting list of tests...', 'ls test/T*.py', return_= 'output').split()
    #tests = filter(lambda x: x.startswith('T') and x.endswith('.py'), sorted( os.listdir('test/') ) )

    #tests = filter(lambda x: x.endswith('.py'), sorted( os.listdir('test/') + os.listdir('demos/') ) )

    def get_py_files(dir_):
	if not os.path.exists(dir_): return []
	return [dir_ + '/' + f for f in os.listdir(dir_) if f.endswith('.py')]

    tests = sorted( get_py_files('test') + get_py_files('demo') )

    print 'Preparingn to run:\n%s\n' % '\n'.join(tests)

    #for t in map(lambda x: x.replace('/', '.'), args[1:]) or tests:
    for t in args[1:] or tests:
        #print '\nRunning %s...' % t
        #__import__( t[:-3] )

        started = datetime.datetime.today()
        execute('\n\nExecuting %s...' % t, 'export PYTHONPATH=`pwd`:$PYTHONPATH && %s %s' % (sys.executable, t) )
        print '\nFinished {} in {}'.format(t, datetime.datetime.today() - started)

        #__import__( 'test.' + t[:-3] )
        #execute('Executing %s...' % t, 'export PYTHONPATH=`pwd`:$PYTHONPATH && python %s' % t)
        #execute('Executing %s...' % t, 'source SetPyRosettaEnvironment.sh && python %s' % t)

    print '\nAll PyRosetta Tests passed!\n'



if __name__ == "__main__":
    print '%s::__main__, started at %s...' % (sys.argv[0], datetime.datetime.now())
    main(sys.argv)
