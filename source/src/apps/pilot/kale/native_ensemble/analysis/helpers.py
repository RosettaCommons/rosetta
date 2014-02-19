from __future__ import division

import sys, os, re
import numpy, pylab
import sqlalchemy, schema

plural = lambda x: (x, '' if x == 1 else 's')

# Rosetta/command helpers

def find_rosetta_root():
    import subprocess
    try: 
        with open(os.devnull, 'w') as devnull:
            git = 'git', 'rev-parse', '--show-toplevel'
            stdout = subprocess.check_output(git, stderr=devnull)
            return stdout.strip()

    except subprocess.CalledProcessError:
        raise RosettaNotFound()

def rosetta_command(command):
    return os.path.join(
            find_rosetta_root(), 'source', 'bin', command)

def analysis_script(command):
    return os.path.join(
            find_rosetta_root(), 'source', 'src', 'apps', 'pilot',
            'kale', 'native_ensemble', 'analysis', command)

def db_to_rosetta(database):
    if database.startswith('sqlite'):
        return database[10:]

def database_arguments(url):
    if url.drivername == 'mysql':
        return [
                '-inout:dbms:mode', url.drivername,
                '-inout:dbms:database_name', url.database,
                '-inout:dbms:user', url.username,
                '-inout:dbms:password', url.password,
                '-inout:dbms:host', url.host,
                '-inout:dbms:port', str(url.port),
        ]
    elif url.drivername == 'sqlite':
        return [
                '-inout:dbms:database_name', url.database,
        ]
    else:
        message = "Unknown database driver '{0.database}'."
        raise ValueError(message.format(url))


# Database helpers

from cStringIO import StringIO
URL, SESSION = None, None

def connect_to_database(url):
    from contextlib import contextmanager

    @contextmanager  # (no fold)
    def session_manager():
        global URL, SESSION
        try:
            engine = sqlalchemy.create_engine(url)
            URL = sqlalchemy.engine.url.make_url(url)
            SESSION = sqlalchemy.orm.sessionmaker(bind=engine)()
            schema.Base.metadata.create_all(engine)
            yield
            SESSION.commit()
        except:
            if SESSION: SESSION.rollback()
            raise
        finally:
            if SESSION: SESSION.close()
            URL, SESSION = None, None

    return session_manager()

def cache_trajectory(job, force=False):
    import sys, subprocess, shlex
    from collections import defaultdict

    assert URL is not None

    command = 'query_trajectory -mute all -query:job {} '
    command = rosetta_command(command.format(job))
    command = shlex.split(command) + database_arguments(URL)
    process = subprocess.Popen(command, stdout=subprocess.PIPE)

    num_frames = 0
    frame = -1

    while process.poll() is None:
        for line in iter(process.stdout.readline, ''):
            try: type, fields = line.split(' ', 1)
            except: continue

            if type == 'header':
                num_frames, loop_begin, loop_end = map(int, fields.split())
                residues = loop_end - loop_begin + 1;

                iterations = numpy.zeros(num_frames)
                scores = numpy.zeros(num_frames)
                rmsds = numpy.zeros(num_frames)
                torsions = defaultdict(lambda: numpy.zeros(num_frames))

            if type == 'iter':
                frame += 1
                iterations[frame] = int(fields)

                status = '\rCaching trajectory for job {} [{}/{}]'
                sys.stdout.write(status.format(job, frame+1, num_frames))
                sys.stdout.flush()

            if type == 'pose':
                score, rmsd = fields.split()
                scores[frame] = float(score)
                rmsds[frame] = float(rmsd)

            if type == 'residue':
                index, phi, psi, omega, chis = fields.split(' ', 4)

                torsions['phi/%s' % index][frame] = float(phi)
                torsions['psi/%s' % index][frame] = float(psi)
                torsions['omega/%s' % index][frame] = float(omega)

                for x, chi in enumerate(chis.split()):
                    torsions['chi%d/%s' % (x+1, index)][frame] = float(chi)

    cache_array(job, 'iterations', iterations, force=force)
    cache_array(job, 'scores', scores, force=force)
    cache_array(job, 'rmsds', rmsds, force=force)
    cache_array(job, 'torsions', torsions, force=force)

def cache_array(job, type, array, force=False):
    query = SESSION.query(schema.NumpyCache).filter_by(job_id=job, type=type)
    exists = bool(query.count())

    if exists and force:
        print "Replacing the '{}' cache.".format(type)
        query.delete()
    if exists and not force:
        print "Skipping the '{}' cache.".format(type)
        return

    output = StringIO()

    if isinstance(array, numpy.ndarray):
        numpy.save(output, array)
        bytes = array.nbytes
    elif isinstance(array, dict):
        numpy.savez(output, **array)
        bytes = sum(x.nbytes for x in array.values())
    else:
        raise TypeError("array must be either numpy.ndarray or dict")

    # This issue here is that the MySQL binary types have maximum sizes.  The 
    # type used to cache numpy arrays is `MEDIUMBLOB', which has a maximum size 
    # of 2**24 bytes (16MB).  If the array we're trying to cache is bigger than 
    # this, this script will just fail gracefully.  The schema may need to be 
    # adjusted to use a `LONGBLOB' to get around the problem.

    blob_column = schema.NumpyCache.data.property.columns[0]
    max_size = blob_column.type.length

    if bytes > max_size:
        error = "The '{}' cache cannot be saved because it is too big ({}MB)."
        print error.format(type, bytes // 2**20)
        return

    cache = schema.NumpyCache(job_id=job, type=type, data=output.getvalue())
    SESSION.add(cache)

def load_array(job, type):
    try:
        query = SESSION.query(schema.NumpyCache)
        record = query.filter_by(job_id=job, type=type).one()
        input = StringIO(record.data)
        return numpy.load(input)

    except sqlalchemy.orm.exc.NoResultFound:
        cache_trajectory(job)
        return load_array(job, type)

def load_arrays(job, type):
    with load_array(job, type) as npz_manager:
        return dict(npz_manager)

def delete_job(job):
    wipe_cache(job)

    query = SESSION.query(schema.Trajectory).filter_by(job_id=job)
    print "Removing %d record%s from table 'trajectories'..." % plural(query.count())
    query.delete()

    query = SESSION.query(schema.Move).filter_by(job_id=job)
    print "Removing %d record%s from table 'moves'..." % plural(query.count())
    query.delete()

    query = SESSION.query(schema.Job).filter_by(id=job)
    print "Removing %d record%s from table 'jobs'..." % plural(query.count())
    query.delete()

def wipe_cache(job):
    query = SESSION.query(schema.NumpyCache).filter_by(job_id=job)
    print "Removing %d record%s from table 'numpy_cache'..." % plural(query.count())
    query.delete()


# Script helpers

def define_args(*specializations):
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('job', type=int)
    parser.add_argument('--database', '-d', default='sqlite:///sandbox.db')

    if 'torsion' in specializations:
        parser.add_argument('--include', '-i', nargs='*', default=[])
        parser.add_argument('--exclude', '-x', nargs='*', default=[])

    if 'plot' in specializations:
        parser.add_argument('--output', '-o')
        parser.add_argument('--quiet', '-q', action='store_true')
        parser.add_argument('--foreground', '-f', action='store_true')

    return parser

def job_title(arguments):
    return '{0.database}/job-{0.job}'.format(arguments)

def pick_torsions(arguments, torsions):
    results = set()
    include = set(arguments.include)
    exclude = set(arguments.exclude)

    abbreviations = {
            'bb': ('phi', 'psi')
    }

    for find, replace in abbreviations.items():
        if find in include:
            include.remove(find)
            include.update(replace)
        if find in exclude:
            exclude.remove(find)
            exclude.update(replace)

    for key in torsions.keys():
        name, index = key.split('/')
        type = name.strip('0123456789')

        # Pass on keys that aren't in the 'include' list.
        if not include: pass
        elif (key not in include) and \
             (name not in include) and \
             (type not in include) and \
             (index not in include): continue

        # Pass on keys that are in the 'exclude' list.
        if key in exclude: continue
        if name in exclude: continue
        if type in exclude: continue
        if index in exclude: continue

        results.add(key)

    return results

def make_plot(arguments):
    from mpldatacursor import datacursor

    finish_script()
    datacursor(formatter='{label}'.format)

    if not arguments.quiet:
        if arguments.foreground:
            pylab.show()
        elif not os.fork():
            pylab.show()
    if arguments.output:
        pylab.savefig(arguments.output)

def finish_script():
    if not sys.stdout.isatty():
        print "Done."
        sys.stdout.flush()

