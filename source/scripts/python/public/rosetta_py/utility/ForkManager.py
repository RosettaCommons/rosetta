import datetime
import fcntl
import os
import pty
import re
import sys
import time

class ForkManager:
	"""If the master python process does not output to a terminal (i.e. it is
	run from a script or is redirected to a file), then you must run it in
	unbuffered mode (i.e. python -u) so that the forked slave processes can
	have their output captured correctly."""

	def __init__(self, max_slaves = 1, incremental_output = False,
	             only_output_bad = False, time_limit = 0, error_output = []):

		self.max_slaves = max_slaves
		self.slaves = {}
		self.bad = []
		self.ismaster = True

		self.incremental_output = incremental_output
		if only_output_bad:
			self.incremental_output = False
			self.only_output_bad = True

		self.time_limit = time_limit

		default_error_output = [
			'python: double free or corruption',
			'python: corrupted double-linked list',
			'python: free\(\): invalid pointer',
			'python: free\(\): invalid next size \(normal\)',
			'python: free\(\): invalid next size \(fast\)',
		]

		self.errorre = re.compile("("+"|".join(default_error_output+error_output)+")")

	def fork(self, name):
		"""Fork the current process. Returns true if we are now in a slave
		process or if max_slaves is set to 0 (in which case no forking
		happens). Otherwise, returns false to the master process after there is
		at least one less than the maximum number of slaves. Note that for code
		brevity, the output of this command is the inverse of the usual fork()
		that returns true to the master and false to the slave."""

		if self.max_slaves == 0:
			return True

		pid, fd = pty.fork()

		if pid:

			self.add(Slave(name, pid, fd))
			self.wait(self.max_slaves - 1)
			return False

		else:

			self.ismaster = False
			return True

	def add(self, slave):
		"""Add a slave, this is usually called by fork() and not by user code"""

		self.slaves[slave.pid] = slave

	def wait(self, num_slaves = 0):
		"""Pause execution until the number of slaves running is less than or
		equal to num_slaves"""

		# wait for one to finish
		while len(self.slaves) > num_slaves:

			time.sleep(.1)

			for pid in self.slaves.keys():

				self.slaves[pid].update()
				if self.incremental_output:
					self.slaves[pid].print_new_output()
					sys.stdout.flush()

				if self.slaves[pid].status != None:
					if not self.incremental_output and (not self.only_output_bad or self.slaves[pid].status):
						self.slaves[pid].print_all_output()
						sys.stdout.flush()
					print self.slaves[pid].name + " finished with status " + str(self.slaves[pid].status) + " duration " + str(datetime.timedelta(seconds = self.slaves[pid].time_end - self.slaves[pid].time_start))
					sys.stdout.flush()
					if self.slaves[pid].status != 0:
						self.bad.append(self.slaves[pid])
					del self.slaves[pid]
					if len(self.slaves) > num_slaves:
						break
					continue

				if self.slaves[pid].killed:
					continue

				if self.time_limit:
					if time.time() - self.slaves[pid].time_start > self.time_limit:
						print self.slaves[pid].name + " exceeded time limit"
						sys.stdout.flush()
						self.slaves[pid].kill()
						continue

				if self.slaves[pid].new_output:
					if self.errorre.search(self.slaves[pid].new_output):
						print self.slaves[pid].name + " output an error"
						sys.stdout.flush()
						self.slaves[pid].kill()
						continue

	def exit_slave(self, status = 0):
		"""Exit from a slave process. Does nothing when called from the master
		process"""

		if not self.ismaster:
			sys.exit(status)

class Slave:

	def __init__(self, name, pid, fd):

		self.name = name
		self.pid = pid
		self.fd = fd
		self.status = None
		self.killed = False
		self.output = []
		self.next_output_line = 0
		self.new_output = False
		self.prepend_name = True
		self.time_start = time.time()
		self.time_end = None

		# turn on non-blocking I/O
		fcntl.fcntl(fd, fcntl.F_SETFL, os.O_NONBLOCK)

	def update(self, status = None):
		"""Update the status and output of a slave"""

		#print "Updating", self.name

		if status != None:
			self.status = status
			# turn off non-blocking I/O
			fcntl.fcntl(fd, fcntl.F_SETFL, 0)

		if self.status == None:
			pid, status = os.waitpid(self.pid, os.WNOHANG)
			if pid != 0:
				self.status = status
				self.time_end = time.time()

		new_output = ""
		try:
			readoutput = os.read(self.fd, 1000)
			while readoutput:
				new_output += readoutput
				readoutput = os.read(self.fd, 1000)
		except:
			None

		self.new_output = ""
		if new_output:
			new_output = new_output.splitlines(True)
			if len(self.output) == 0 or self.output[-1].endswith("\n"):
				self.output += new_output
				self.new_output = "".join(new_output)
			else:
				self.output[-1] += new_output[0]
				self.output += new_output[1:]
				self.new_output = "".join([self.output[-1]+new_output[0]] + new_output[1:])

		if self.status != None and self.fd != None:
			os.close(self.fd)
			self.fd = None

		#if self.output:
		#	print self.name + ": " + self.output.strip()

	def kill(self, signal = 15):

		os.kill(self.pid, signal)
		self.killed = True

	def print_new_output(self):

		next_output_line = self.next_output_line
		if len(self.output) > self.next_output_line:
			next_output_line = len(self.output)
			if not self.output[-1].endswith("\n"):
				next_output_line -= 1
		if next_output_line > self.next_output_line:
			output = self.output[self.next_output_line:next_output_line]
			if self.prepend_name:
				output = [self.name + ": " + line for line in output]
			print "".join(output),
		self.next_output_line = next_output_line

	def print_all_output(self):

		if len(self.output):
			output = self.output
			if self.prepend_name:
				output = [self.name + ": " + line for line in output]
			print "".join(output),
