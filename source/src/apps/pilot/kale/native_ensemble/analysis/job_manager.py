#!/usr/bin/env python

from __future__ import division

import pygtk; pygtk.require('2.0'); import gtk, gobject
import sqlalchemy, schema, helpers

# Wanted Features
# ===============
# 1. Copy command-line to clipboard.
# 2. Right click on jobs to launch scripts.

class AppWindow (gtk.Window):

    def __init__(self, *urls):
        gtk.Window.__init__(self, gtk.WINDOW_TOPLEVEL)
        self.set_border_width(2)
        self.set_default_size(644, 365)
        self.add_events(gtk.gdk.KEY_PRESS_MASK)
        self.connect('delete_event', self.quit)
        self.connect('key_press_event', self.on_hotkey)

        self.vbox = gtk.VBox()
        self.vbox.set_spacing(2)
        self.add(self.vbox)

        self.jobs = DatabaseWindow(*urls)
        self.vbox.pack_start(self.jobs)

        self.commands = CommandWindow(self.jobs)
        self.vbox.pack_start(self.commands, expand=False)

        self.show_all()

    def main(self, foreground=False):
        import os
        if foreground: gtk.main()
        elif not os.fork(): gtk.main()

    def quit(self, widget, event=None):
        gtk.main_quit()
        return False

    def push_status(self, widget):
        self.vbox.pack_start(widget, expand=False)
        self.vbox.reorder_child(widget, 1)

    def pop_status(self, widget):
        if widget in self.vbox.children():
            self.vbox.remove(widget)

    def on_hotkey(self, widget, event):
        key = gtk.gdk.keyval_name(event.keyval)
        control = (event.state & gtk.gdk.CONTROL_MASK)

        if control and key == 's': 
            self.jobs.show_move_stats()
            return True


class DatabaseWindow (gtk.Notebook):

    def __init__(self, *urls):
        gtk.Notebook.__init__(self)
        self.pages = []

        hbox = gtk.HBox()
        self.set_action_widget(hbox, gtk.PACK_END)

        # Make a page for each database.

        for url in urls:
            self.add_page(url)

        # Make an add button.
        
        image = gtk.Image()
        image.set_from_stock(gtk.STOCK_ADD, gtk.ICON_SIZE_BUTTON)
        button = gtk.Button()
        button.add(image)
        button.set_relief(gtk.RELIEF_NONE)
        button.connect('clicked', self.open_page)
        hbox.add(button)

        # Make a refresh button.

        image = gtk.Image()
        image.set_from_stock(gtk.STOCK_REFRESH, gtk.ICON_SIZE_BUTTON)
        button = gtk.Button()
        button.add(image)
        button.set_relief(gtk.RELIEF_NONE)
        button.connect('clicked', self.refresh)
        hbox.add(button)

        hbox.show_all()

    def refresh(self, *args):
        for page in self.pages:
            page.refresh()

    def add_page(self, url):
        page = JobWindow(url)
        page.show_all()
        label = page.get_title()
        self.append_page(page, label)
        self.pages.append(page)

    def open_page(self, *args):
        dialog = gtk.Dialog(
                buttons=(
                    gtk.STOCK_CANCEL, False,
                    gtk.STOCK_OK, True))
        width = 75

        # Ask for a file name.

        label = gtk.Label("Protocol")
        label.set_size_request(width, -1)

        combo = gtk.combo_box_new_text()
        combo.append_text('sqlite')
        combo.append_text('mysql')
        combo.set_active(0)

        hbox = gtk.HBox()
        hbox.pack_start(label, expand=False, padding=5)
        hbox.pack_start(combo, padding=5)
        dialog.vbox.pack_start(hbox, padding=5)

        # Ask how many frames to record.

        label = gtk.Label("Path")
        label.set_size_request(width, -1)

        entry = gtk.Entry()

        hbox = gtk.HBox()
        hbox.pack_start(label, expand=False, padding=5)
        hbox.pack_start(entry, padding=5)
        dialog.vbox.add(hbox)

        dialog.show_all()
        proceed = dialog.run()
        if proceed:
            protocol = combo.get_model()[combo.get_active()][0]
            path = entry.get_text()

        dialog.destroy()
        if not proceed:
            return 

        if protocol == 'sqlite':
            path = 'sqlite:///{}'.format(path)
        if protocol == 'mysql':
            path = 'mysql:///{}'.format(path)

        self.add_page(path)

    def get_page(self):
        index = self.get_current_page()
        return self.get_nth_page(index)

    def show_move_stats(self, *args):
        page = self.get_page()
        popup = MoveStatsWindow(page)


class JobWindow (gtk.ScrolledWindow):

    def __init__(self, url):
        gtk.ScrolledWindow.__init__(self)

        self.url = url
        engine = sqlalchemy.create_engine(url)
        self.session_factory = sqlalchemy.orm.sessionmaker(bind=engine)

        self.model = gtk.ListStore(
                object, int, str, str, str, str, int, str, str)
        self.view = gtk.TreeView(self.model)
        self.view.set_rubber_banding(True)
        self.view.set_rules_hint(True)
        self.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.add_with_viewport(self.view)

        selection = self.view.get_selection()
        selection.set_mode(gtk.SELECTION_MULTIPLE)

        right_align = gtk.CellRendererText()
        left_align = gtk.CellRendererText()
        left_align.set_property('xalign', 1.0)
        editable = gtk.CellRendererText()
        editable.set_property('editable', True)
        editable.connect('edited', self.edit_description)

        columns = [
                gtk.TreeViewColumn("Job", right_align, text=1),
                gtk.TreeViewColumn("Algorithm", right_align, text=2),
                gtk.TreeViewColumn("Protein", right_align, text=3),
                gtk.TreeViewColumn("Residues", right_align, text=4),
                gtk.TreeViewColumn("Duration", left_align, text=5),
                gtk.TreeViewColumn("Iterations", left_align, text=6),
                gtk.TreeViewColumn("Efficiency", left_align, text=7),
                gtk.TreeViewColumn("Description", editable, text=8),
        ]

        for column in columns:
            self.view.append_column(column)

        self.refresh()

    def refresh(self):
        session = self.session_factory()
        selection = self.view.get_selection()
        ids_selected = [job.id for job in self.get_selected_jobs()]

        self.model.clear()

        for record in session.query(schema.Job).all():
            query = session.query(schema.Move).filter_by(job_id=record.id)
            moves = query.all()

            try:
                trials = sum(x.num_trials for x in moves)
                accepted = sum(x.num_accepted for x in moves)
                efficiency = '%.2f%%' % (100 * accepted / trials)
            except ZeroDivisionError:
                continue
            
            row = self.model.append()
            self.model.set(row, 0, record)
            self.model.set(row, 1, record.id)
            self.model.set(row, 2, record.algorithm)
            self.model.set(row, 3, record.protein)
            self.model.set(row, 4, record.loop)
            self.model.set(row, 5, record.duration)
            self.model.set(row, 6, record.iterations)
            self.model.set(row, 7, efficiency)
            self.model.set(row, 8, record.description)

            if record.id in ids_selected:
                selection.select_iter(row)
            if not ids_selected and not selection.count_selected_rows():
                selection.select_iter(row)

    def get_title(self):
        return gtk.Label(self.url)

    def get_selected_jobs(self):
        jobs = []
        def append_result(model, path, iter):
            job = self.model.get_value(iter, 0)
            jobs.append(job)

        selection = self.view.get_selection()
        selection.selected_foreach(append_result)
        return jobs

    def edit_description(self, *args):
        cell, column, description = args
        iter = self.model.get_iter(column)
        id = self.model.get_value(iter, 1)

        try:
            session = self.session_factory()
            record = session.query(schema.Job).filter_by(id=id).one()
            record.description = description
            session.add(record)
            session.commit()
        except:
            session.rollback()
            raise
        finally:
            session.close()

        self.refresh()


class CommandWindow (gtk.HBox):
    
    def __init__(self, database_window):
        gtk.HBox.__init__(self)

        self.database_window = database_window
        self.status_window = StatusWindow()
        vbox = gtk.VBox()

        for script in CommandManager.scripts:
            button = gtk.Button(script.title)
            button.set_size_request(125, -1)
            button.connect('clicked', self.launch, script)

            hbox = gtk.HBox()
            hbox.pack_start(button, expand=False)
            script.options(hbox, self.launch)
            vbox.add(hbox)

        self.pack_start(vbox)
        self.pack_start(self.status_window, padding=5)
        
    def launch(self, widget, script):
        page = self.database_window.get_page()
        jobs = page.get_selected_jobs()

        for runner in script.run(page, jobs):
            runner.run(self.status_window)


class CommandRunner (gtk.EventBox):

    def __init__(self, title, command):
        gtk.EventBox.__init__(self)

        self.add_events(gtk.gdk.BUTTON_PRESS_MASK)
        self.connect('button-press-event', self.menu)

        self.title = title
        self.command = command
        self.echo_mode = False

    def __str__(self):
        return '%s: %s' % (self.title, self.command)

    def run(self, status_window):
        import subprocess, os, fcntl
        from subprocess import PIPE, STDOUT

        self.status_window = status_window

        # Run the command such that it doesn't block.

        self.buffer = ''
        self.process = subprocess.Popen(
                self.command, stdout=PIPE, stderr=STDOUT, shell=True)

        fileno = self.process.stdout.fileno()
        flags = fcntl.fcntl(fileno, fcntl.F_GETFL) | os.O_NONBLOCK
        fcntl.fcntl(fileno, fcntl.F_SETFL, flags)
        gobject.io_add_watch(fileno, gobject.IO_IN | gobject.IO_HUP, self.poll)

        # Create a dialog to show the status of the command.

        self.progress = gtk.ProgressBar()
        self.progress.set_text("Running %s..." % self.title.lower())

        self.add(self.progress)
        self.status_window.push(self)
        self.show_all()

    def poll(self, fileno, condition):
        import os, re, sys

        job_finished = (condition == gobject.IO_HUP)

        # Read as much as possible out of the pipe.

        try:
            while True:
                output = os.read(fileno, 4096); self.buffer += output
                progress_update = re.search(r'(.*)\[(\d+)/(\d+)\]', output)
                script_finished = re.match(r'^Done\.$', output)

                if self.status_window.always_show_output:
                    sys.stdout.write(output)
                    sys.stdout.flush()

                if not output:
                    break

                if progress_update:
                    message, current, total = progress_update.groups()
                    fraction = min(1, int(current) / int(total))
                    self.progress.set_fraction(fraction)

                    message = message.strip()
                    if fraction == 1: message = "Running %s" % self.title
                    if message: self.progress.set_text(message.strip())

                if script_finished:
                    self.progress.set_fraction(1)
                    self.progress.set_text("Running %s" % self.title)

        except OSError:
            pass

        # Update the gui if the job has finished.

        if job_finished:
            self.process.wait()
            self.buffer += self.process.stdout.read()

            app.jobs.refresh()

            if self.process.poll() == 0:    # Success
                self.dismiss()
            else:                           # Failure
                message = "Failed with error code {}."
                message = message.format(self.process.poll())
                self.progress.set_text(message)

        return not job_finished

    def menu(self, widget, event):
        if event.button != 3: return False

        menu = gtk.Menu()

        if self.process.poll() is None:
            abort_item = gtk.MenuItem('Abort')
            abort_item.connect('activate', self.abort)
            menu.append(abort_item)

        dismiss_item = gtk.MenuItem('Dismiss')
        dismiss_item.connect('activate', self.dismiss)
        menu.append(dismiss_item)

        dump_item = gtk.MenuItem('Show output')
        dump_item.connect('activate', self.dump)
        menu.append(dump_item)

        menu.show_all()
        menu.popup(None, None, None, event.button, event.time)

        return True

    def dismiss(self, *args):
        self.status_window.pop(self)

    def abort(self, *args):
        self.dismiss()
        if self.process.poll() is None:
            self.process.kill()

    def dump(self, *args):
        return_code = self.process.poll()
        header = self.command

        if return_code is not None:
            self.dismiss()

        print header
        print '=' * min(len(header), 80)
        print self.buffer
        if return_code > 0:
            print "Command failed with error code %d." % self.process.poll()
        print

    def echo(self, value):
        self.echo_mode = value


class StatusWindow (gtk.VBox):

    def __init__(self):
        gtk.VBox.__init__(self)
        self.always_show_output = False

        # Setup the status window.

        self.view = gtk.Viewport()
        self.empty_label = gtk.Label("Waiting for commands...")
        self.empty_label.set_sensitive(False)
        self.status_list = gtk.VBox()

        window = gtk.ScrolledWindow()
        window.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
        window.add(self.view)

        # Setup the control widgets.

        def dismiss_all(button):
            for runner in self.status_list.get_children():
                runner.dismiss()

        def abort_all(button):
            for runner in self.status_list.get_children():
                runner.abort()

        def toggle_output(toggle):
            self.always_show_output = toggle.get_active()

        self.dismiss_button = gtk.Button('Dismiss All')
        self.dismiss_button.connect('clicked', dismiss_all)
        self.abort_button = gtk.Button('Abort All')
        self.abort_button.connect('clicked', abort_all)
        self.output_toggle = gtk.CheckButton('Always show output')
        self.output_toggle.set_alignment(1, 0.5)
        self.output_toggle.connect('toggled', toggle_output)

        row_1 = gtk.HBox()
        row_2 = gtk.HBox()

        row_1.pack_end(self.dismiss_button, expand=False)
        row_1.pack_end(self.abort_button, expand=False)
        row_2.pack_end(self.output_toggle, expand=False)

        self.pack_start(window, padding=5)
        self.pack_end(row_2, expand=False)
        self.pack_end(row_1, expand=False)

        self.refresh()

    def refresh(self):
        if not self.status_list.get_children():
            active = False
            widget = self.empty_label
        else:
            active = True
            widget = self.status_list

        self.dismiss_button.set_sensitive(active)
        self.abort_button.set_sensitive(active)

        if self.view.get_child():
            self.view.remove(self.view.get_child())

        self.view.add(widget)
        self.view.show_all()

    def push(self, runner):
        self.status_list.pack_start(runner, expand=False)
        self.status_list.reorder_child(runner, 0)
        self.refresh()

    def pop(self, runner):
        if runner in self.status_list.get_children():
            self.status_list.remove(runner)
            self.refresh()


class MoveStatsWindow (gtk.Window):

    def __init__(self, page):
        gtk.Window.__init__(self)
        self.set_type_hint(gtk.gdk.WINDOW_TYPE_HINT_DIALOG)

        session = page.session_factory()
        jobs = page.get_selected_jobs()

        model = gtk.ListStore(object, int, str, int, int, str)
        view = gtk.TreeView(model)
        view.get_selection().set_mode(gtk.SELECTION_NONE)

        left_align = gtk.CellRendererText()
        right_align = gtk.CellRendererText()
        right_align.set_property('xalign', 1)

        columns = [
                gtk.TreeViewColumn('Job', left_align, text=1),
                gtk.TreeViewColumn('Move', left_align, text=2),
                gtk.TreeViewColumn('Accepted', right_align, text=3),
                gtk.TreeViewColumn('Proposed', right_align, text=4),
                gtk.TreeViewColumn('Efficiency', right_align, text=5),
        ]
        for column in columns:
            view.append_column(column)

        columns[0].set_sort_column_id(1)
        columns[1].set_sort_column_id(2)

        self.add(view)
        self.show_all()

        for job in jobs:
            moves = session.query(schema.Move).filter_by(job_id=job.id).all()
            for move in moves:
                row = model.append()
                efficiency = '%.2f%%' % (100 * move.efficiency)
                model.set(row, 0, move)
                model.set(row, 1, job.id)
                model.set(row, 2, move.type)
                model.set(row, 3, move.num_accepted)
                model.set(row, 4, move.num_trials)
                model.set(row, 5, efficiency)



class CommandManager (type):
    scripts = []

    def __init__(cls, name, bases, dict):
        super(CommandManager, cls).__init__(name, bases, dict)
        if name != 'CommandLauncher':
            CommandManager.scripts.append(cls())


class CommandLauncher (object):
    __metaclass__ = CommandManager

    title = ''
    default = ''
    command = 'echo {job} --database {url} {options}'

    def options(self, hbox, callback):
        self.entry = gtk.Entry()
        self.entry.set_text(self.default)
        self.entry.connect('activate', callback, self)
        hbox.pack_start(self.entry, expand=True, fill=True)

    def run(self, page, jobs):
        url = page.url
        options = self.entry.get_text()
        commands = []

        for job in jobs:
            command = self.command.format(url=url, job=job.id, options=options)
            commands.append(CommandRunner(self.title, command))

        return commands



class DeleteJob (CommandLauncher):
    title = "Delete job"
    command = helpers.analysis_script(
            'delete_job.py {job} --database {url} {options}')

    def run(self, page, jobs):
        ids = [str(x.id) for x in jobs]

        if not jobs:
            raise AssertionError
        elif len(jobs) == 1:
            message = 'job {}'.format(ids[0])
        elif len(jobs) == 2:
            message = 'jobs {} and {}'.format(*ids)
        else:
            parts = str.join(', ', ids[:-1]), ids[-1]
            message = 'jobs {} and {}'.format(*parts)

        # Confirm that the user really means to delete something.

        dialog = gtk.MessageDialog(
            flags=gtk.DIALOG_DESTROY_WITH_PARENT,
            type=gtk.MESSAGE_QUESTION, 
            buttons=gtk.BUTTONS_OK_CANCEL,
            message_format="Delete %s?" % message)

        dialog.format_secondary_text(
                "Please confirm that you mean to take this action.  "
                "Deleted jobs cannot be recovered.")

        result = dialog.run()
        dialog.destroy()

        # Run the script, if the user gave confirmation.

        if result == gtk.RESPONSE_CANCEL:
            return []

        return CommandLauncher.run(self, page, jobs)
    

class WipeCache (CommandLauncher):
    title = "Wipe cache"
    command = helpers.analysis_script(
            'wipe_cache.py {job} --database {url} {options}')
    
class ScoreVsRmsd (CommandLauncher):
    title = "Score vs rmsd"
    command = helpers.analysis_script(
            'score_vs_rmsd.py {job} --database {url} --foreground {options}')

class ScoreVsTime (CommandLauncher):
    title = "Score vs time"
    command = helpers.analysis_script(
            'score_vs_time.py {job} --database {url} --foreground {options}')

class TorsionVsTime (CommandLauncher):
    title = "Torsion vs time"
    default = '-i chi1'
    command = helpers.analysis_script(
            'torsion_vs_time.py {job} --database {url} --foreground {options}')

class CorrelationTime (CommandLauncher):
    title = "Autocorrelation"
    default = '-i chi1'
    command = helpers.analysis_script(
            'autocorrelation.py {job} --database {url} --foreground {options}')

class SaveMovie (CommandLauncher):
    title = "Save movie"

    def run(self, page, jobs):
        commands = []

        for job in jobs:
            path = self.entry.get_text().format(job)
            frames = self.combo.get_model()[self.combo.get_active()][0]

            title = "trajectory_movie on job %s" % job.id
            database = helpers.db_to_rosetta(page.url)
            command = 'trajectory_movie ' \
                    '-movie:job {} -movie:out {} -movie:frames {} -database_name {}'
            command = helpers.rosetta_command(command)
            command = command.format(job.id, path, frames, database)
            commands.append(CommandRunner(self.title, command))

        return commands

    def options(self, hbox, callback):
        self.entry = gtk.Entry()
        self.entry.set_text('{0.id}.movie.pdb')
        self.entry.connect('activate', callback, self)

        self.combo = gtk.combo_box_new_text()
        for i in range(50, 500, 50):
            self.combo.append_text(str(i))
        self.combo.set_active(4)

        hbox.pack_start(self.entry, expand=True, fill=True)
        hbox.pack_start(self.combo, expand=False)



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('databases', nargs='*', default=['all'])
    parser.add_argument('--foreground', '-f', action='store_true')
    arguments = parser.parse_args()

    # Pick which databases to show.

    urls = set()
    for db in arguments.databases:
        if db == 'all':
            urls.add('sqlite:///sandbox.db')
        elif db == 'sqlite':
            urls.add('sqlite:///sandbox.db')
        elif db == 'mysql':
            pass
        else:
            urls.add(db)

    # Run the GUI.

    app = AppWindow(*urls)
    app.main(arguments.foreground)

