"""
Interface for MDEngine objects. New MDEngines can be implemented by constructing a new class that inherits from MDEngine
and implements its abstract methods.
"""

import abc
import os
import re
import math
import glob
import numpy
import shutil
import django
import random
import pytraj
import mdtraj
import mdtraj.formats
from atesa import utilities

class MDEngine(abc.ABC):
    """
    Abstract base class for molecular dynamics engines.

    Implements methods for all of the engine-specific tasks that ATESA might need.

    """

    @abc.abstractmethod
    def get_frame(self, trajectory, frame, settings):
        """
        Return a new file containing just the frame'th frame of a trajectory in Amber .rst7 format

        Parameters
        ----------
        trajectory : str
            Name of trajectory file to obtain last frame from
        frame : int
            Index of frame to return; 1-indexed, -1 gives last frame, 0 is invalid
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        last_frame : str
            Name of .rst7 format coordinate file corresponding to desired frame of trajectory, if it exists; an empty
            string otherwise

        """
        pass

    @abc.abstractmethod
    def get_input_file_aimless_shooting(self, settings, thread, job_index, **kwargs):
        """
        Build from template with kwargs and return the input file for an aimless_shooting simulation

        Parameters
        ----------
        settings : argparse.Namespace
            Settings namespace object
        thread : Thread
            The relevant Thread object on which to act
        job_index : int
            Index corresponding to the relevant step within thread.current_type
        kwargs : dict
            Keyword arguments used to to fill the templated input file, as needed

        Returns
        -------
        input_file : str
            Name of the simulation input file

        """
        pass

    @abc.abstractmethod
    def get_input_file_committor_analysis(self, settings, thread, job_index, **kwargs):
        """
        Build from template with kwargs and return the input file for a committor analysis simulation

        Parameters
        ----------
        settings : argparse.Namespace
            Settings namespace object
        thread : Thread
            The relevant Thread object on which to act
        job_index : int
            Index corresponding to the relevant step within thread.current_type
        kwargs : dict
            Keyword arguments used to to fill the templated input file, as needed

        Returns
        -------
        input_file : str
            Name of the simulation input file

        """
        pass

    @abc.abstractmethod
    def get_input_file_find_ts(self, settings, thread, job_index, **kwargs):
        """
        Build from template with kwargs and return the input file for a find_ts simulation

        Parameters
        ----------
        settings : argparse.Namespace
            Settings namespace object
        thread : Thread
            The relevant Thread object on which to act
        job_index : int
            Index corresponding to the relevant step within thread.current_type
        kwargs : dict
            Keyword arguments used to to fill the templated input file, as needed

        Returns
        -------
        input_file : str
            Name of the simulation input file

        """
        pass

    @abc.abstractmethod
    def get_input_file_umbrella_sampling(self, settings, thread, job_index, **kwargs):
        """
        Build from template with kwargs and return the input file for an umbrella sampling simulation

        Parameters
        ----------
        settings : argparse.Namespace
            Settings namespace object
        thread : Thread
            The relevant Thread object on which to act
        job_index : int
            Index corresponding to the relevant step within thread.current_type
        kwargs : dict
            Keyword arguments used to to fill the templated input file, as needed

        Returns
        -------
        input_file : str
            Name of the simulation input file

        """
        pass

    @abc.abstractmethod
    def get_input_file_equilibrium_path_sampling(self, settings, thread, job_index, **kwargs):
        """
        Build from template with kwargs and return the input file for an equilibrium path sampling simulation

        Parameters
        ----------
        settings : argparse.Namespace
            Settings namespace object
        thread : Thread
            The relevant Thread object on which to act
        job_index : int
            Index corresponding to the relevant step within thread.current_type
        kwargs : dict
            Keyword arguments used to to fill the templated input file, as needed

        Returns
        -------
        input_file : str
            Name of the simulation input file

        """
        pass

    @abc.abstractmethod
    def control_rst_format(self, restart_file):
        """
        Handle building an Amber-style .rst7 file from whatever the restart_file actually was.

        This method is used whenever there is a need specifically for an .rst7 restart file with velocity information.
        When velocity information is optional or not required, it's easier to use mdengine.get_frame.

        Parameters
        ----------
        restart_file : str
            Path to restart file directly produced by md engine

        Returns
        -------
        rst7 : str
            Path to newly created .rst7 file

        """
        pass

    @abc.abstractmethod
    def simulation_cleanup(self, thread_history, settings):
        """
        Handle any cleanup needed after an simulation job completes but before it is interpreted.

        Parameters
        ----------
        thread_history : argparse.Namespace
            Thread history object
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        None

        """
        pass


class AdaptAmber(MDEngine):
    """
    Adapter class for Amber MDEngine.

    """

    def get_frame(self, trajectory, frame, settings):
        new_restart_name = trajectory + '_frame_' + str(frame) + '.rst7'
        if not os.path.exists(trajectory):
            return ''   # since it's possible to call this before the trajectory file has been initialized
        if frame >= 1:
            shift_frame = frame - 1     # because write_traj is 0-indexed but get_frame is 1-indexed
        elif frame == -1:
            shift_frame = -1
        else:
            raise IndexError('invalid frame index for get_frame: ' + str(frame) + ' (must be >= 1, or exactly -1)')

        traj = mdtraj.load(trajectory, top=settings.topology)
        if traj.n_frames == 0:
            return ''

        try:
            traj[shift_frame].save_amberrst7(new_restart_name, force_overwrite=True)
            # pytraj.write_traj(new_restart_name, traj, format='rst7', frame_indices=[shift_frame], options='multi', overwrite=True, velocity=True)
        except IndexError:
            raise IndexError('frame index ' + str(frame) + ' is out of range for trajectory: ' + trajectory)
        if not os.path.exists(new_restart_name):
            raise OSError('expected mdtraj to write ' + new_restart_name + ' but it was not found.\nThe most likely '
                          'explanation is that the ATESA process is unable to write to the working directory (' +
                          settings.working_directory + '). You may have run out of storage space.')

        return new_restart_name

    def get_input_file_aimless_shooting(self, settings, thread, job_index, **kwargs):
        return settings.path_to_input_files + '/' + settings.job_type + '_' + thread.current_type[job_index] + '_' + settings.md_engine.lower() + '.in'

    def get_input_file_committor_analysis(self, settings, thread, job_index, **kwargs):
        return settings.path_to_input_files + '/' + settings.job_type + '_' + thread.current_type[job_index] + '_' + settings.md_engine.lower() + '.in'

    def get_input_file_umbrella_sampling(self, settings, thread, job_index, **kwargs):
        input_file_name = 'umbrella_sampling_' + str(thread.history.window) + '_' + str(thread.history.index) + '.in'

        # If the appropriate file already exists, simply return it
        if os.path.exists(settings.working_directory + '/' + input_file_name):
            return input_file_name

        # Build restraint file that implements us_pathway_restraints_file, if applicable
        def add_pathway_restraints(thread, settings):
            """
            Build appropriate pathway restraints DISANG file for the window in the specified thread

            Parameters
            ----------
            settings : argparse.Namespace
                Settings namespace object

            Returns
            -------
            None

            """
            if not True in ['nmropt=1' in line for line in
                            open(settings.path_to_input_files + '/umbrella_sampling_prod_' + settings.md_engine.lower() + '.in',
                                 'r').readlines()]:
                raise RuntimeError('did not find \'nmropt=1\' in input file: ' + settings.path_to_input_files +
                                   '/umbrella_sampling_prod_' + settings.md_engine.lower() + '.in, which is required when a'
                                   ' us_pathway_restraints_file is specified in the config file. Make sure that '
                                   'exactly that string is present.')
            if not os.path.exists(settings.working_directory + '/us_pathway_restraints_' + str(thread.history.window) +
                                  '.DISANG'):
                if not os.path.exists(settings.us_pathway_restraints_file):
                    raise FileNotFoundError('Attempted to implement us_pathway_restraints for umbrella sampling, but '
                                            'could not find the indicated file: ' + settings.us_pathway_restraints_file)

                # Define a helper function
                def closest(lst, K):
                    # Return index of closest value to K in list lst
                    return min(range(len(lst)), key=lambda i: abs(lst[i] - float(K)))

                # Define function to clean up some floating point precision issues;
                # e.g.,  numpy.arange(-6,12,0.1)[-1] = 11.899999999999935,
                # whereas safe_arange(-6,12,0.1)[-1] = 11.9
                def safe_arange(start, stop, step):
                    return step * numpy.arange(start / step, stop / step)

                window_centers = safe_arange(settings.us_rc_min, settings.us_rc_max, settings.us_rc_step)

                # Build min_max.pkl file if it doesn't already exist
                if not os.path.exists('min_max.pkl'):
                    # Partition each line in as_full_cvs.out into the nearest US window
                    as_full_cvs = open(settings.us_pathway_restraints_file, 'r')  # load file as an iterator
                    as_full_cvs_line = next(as_full_cvs)  # get first line

                    # Initialize results list; first index is window, second index is CV, third index is 0 for min and 1 for max
                    # E.g., min_max[12][5][1] is the max value of CV6 in the 13th window
                    min_max = [[[None, None] for null in range(len(as_full_cvs_line.split()))] for null in
                               range(len(window_centers))]

                    rc_minmax = [[], []]  # this minmax is for passing to utilities.get_cvs.reduce_cv
                    if settings.rc_reduced_cvs:
                        # Prepare cv_minmax list
                        asout_lines = [[float(item) for item in
                                        line.replace('A <- ', '').replace('B <- ', '').replace(' \n', '').replace('\n',
                                        '').split(' ')] for line in open(settings.as_out_file, 'r').readlines()]
                        open(settings.as_out_file, 'r').close()
                        mapped = list(map(list, zip(*asout_lines)))
                        rc_minmax = [[numpy.min(item) for item in mapped], [numpy.max(item) for item in mapped]]

                        def reduce_cv(unreduced_value, local_index, rc_minmax):
                            # Returns a reduced value for a CV given an unreduced value and the index within as.out corresponding to that CV
                            this_min = rc_minmax[0][local_index]
                            this_max = rc_minmax[1][local_index]
                            return (float(unreduced_value) - this_min) / (this_max - this_min)

                    while True:  # loop can only end via 'break'
                        try:
                            # Split line into cvs
                            this_line = as_full_cvs_line.split()

                            # Evaluate the reaction coordinate value for this line
                            if settings.rc_reduced_cvs:
                                this_line_temp = []
                                for cv_index in range(len(this_line)):
                                    if 'cv' + str(
                                            int(cv_index + 1)) in settings.rc_definition.lower():  # this catches some errant CVs (e.g., CV22 gets 2 and 22) but it's okay
                                        this_line_temp.append(reduce_cv(this_line[cv_index], cv_index, rc_minmax))
                                    else:
                                        this_line_temp.append(None)
                                rc = utilities.evaluate_rc(settings.rc_definition, this_line_temp)
                            else:
                                rc = utilities.evaluate_rc(settings.rc_definition, this_line)
                            window_index = closest(window_centers, rc)

                            # Add to min and max as appropriate
                            for cv_index in range(len(this_line)):
                                if (min_max[window_index][cv_index][0] is None or float(
                                        min_max[window_index][cv_index][0]) > float(
                                        this_line[cv_index])) and not math.isnan(float(this_line[cv_index])):
                                    min_max[window_index][cv_index][0] = this_line[cv_index]  # set new min
                                if (min_max[window_index][cv_index][1] is None or float(
                                        min_max[window_index][cv_index][1]) < float(
                                        this_line[cv_index])) and not math.isnan(float(this_line[cv_index])):
                                    min_max[window_index][cv_index][1] = this_line[cv_index]  # set new max

                            as_full_cvs_line = next(as_full_cvs)  # get next line for next iteration
                        except StopIteration:  # raised by next(as_full_cvs) when trying to read past the end
                            break

                    # Dump as a .pkl file so we don't have to do this againZ
                    pickle.dump(min_max, open('min_max.pkl', 'wb'), protocol=2)

                # Now build the actual DISANG file.
                min_max = pickle.load(open('min_max.pkl', 'rb'))  # load the min_max pickle file
                window_index = closest(window_centers, thread.history.window)  # identify the appropriate window_index
                open(settings.working_directory + '/us_pathway_restraints_' + str(thread.history.window) + '.DISANG',
                     'w').close()
                with open(settings.working_directory + '/us_pathway_restraints_' + str(thread.history.window) +
                          '.DISANG', 'a') as f:
                    f.write(
                        'ATESA automatically generated restraint file implementing us_pathway_restraints_file option\n')
                    for cv_index in range(len(min_max[window_index])):
                        atoms, optype, nat = utilities.interpret_cv(cv_index + 1, settings)  # get atom indices and type for this CV
                        try:
                            this_min = float(min_max[window_index][cv_index][0])
                            this_max = float(min_max[window_index][cv_index][1])
                        except TypeError:  # no min and/or max for this window_index, so skip it
                            continue
                        f.write(' &rst iat=' + ', '.join([str(atom) for atom in atoms]) + ', r1=' +
                                str(this_min - (this_max - this_min)) + ', r2=' + str(this_min) + ', r3=' +
                                str(this_max) + ', r4=' + str(
                            this_max + (this_max - this_min)) + ', rk2=100, rk3=100, /\n')

        if settings.us_pathway_restraints_file:
            add_pathway_restraints(thread, settings)

        # Break down settings.rc_definition into component terms
        # rc_definition must be a linear combination of terms for US anyway, so our strategy here is to condense it
        # into just numbers, CV terms, and +/-/* characters, and then split into a list to iterate over.
        condensed_rc = settings.rc_definition.replace(' ', '').replace('--', '+').replace('+-', '-').replace('-', '+-')
        condensed_rc = [item for item in condensed_rc.split('+') if not item in ['', ' ']]
        alp0 = sum([float(item) for item in condensed_rc if not 'CV' in item])  # extract constant terms
        condensed_rc = [item for item in condensed_rc if 'CV' in item]  # remove constant terms
        if len(condensed_rc) == 0:
            raise RuntimeError(
                'it appears that the provided reaction coordinate contains no CVs or it otherwise '
                'improperly formatted. Only linear combinations of constant terms and CVs (with '
                'coefficients) are permitted. The offending RC definition is: ' + settings.rc_definition)

        # Prepare cv_minmax list for evaluating min and max terms later, if appropriate
        if settings.rc_reduced_cvs:
            asout_lines = [[float(item) for item in
                            line.replace('A <- ', '').replace('B <- ', '').replace(' \n', '').replace('\n', '').split(
                                ' ')] for line in open(settings.as_out_file, 'r').readlines()]
            open(settings.as_out_file, 'r').close()
            mapped = list(map(list, zip(*asout_lines)))
            rc_minmax = [[numpy.min(item) for item in mapped], [numpy.max(item) for item in mapped]]

        if settings.us_implementation.lower() == 'amber_rxncor':
            # Make sure template input file includes 'irxncor=1'
            if not True in ['irxncor=1' in line for line in open(settings.path_to_input_files + '/umbrella_sampling_prod_' + settings.md_engine.lower() + '.in', 'r').readlines()]:
                raise RuntimeError('did not find \'irxncor=1\' in input file: ' + settings.path_to_input_files +
                                   '/umbrella_sampling_prod_' + settings.md_engine.lower() + '.in, but it is required for '
                                   'umbrella sampling with us_implementation = \'amber_rxncor\'. Make sure that exactly'
                                   ' that string is present.')
            shutil.copy(settings.path_to_input_files + '/umbrella_sampling_prod_' + settings.md_engine.lower() + '.in',
                        settings.working_directory + '/' + input_file_name)

            # Here, we'll write the &rxncor and &rxncor_order_parameters namelists manually
            with open(settings.working_directory + '/' + input_file_name, 'a') as file:
                # if not open(settings.working_directory + '/' + input_file_name, 'r').readlines()[-1] == '':    # if not last line of template empty
                #     file.write('\n')
                file.write(' &rxncor\n')
                file.write('  rxn_dimension=' + str(settings.rc_definition.count('CV')) + ',\n')
                file.write('  rxn_kconst=' + str(settings.us_restraint) + ',\n')
                file.write('  rxn_c0=' + str(thread.history.window) + ',\n')
                file.write('  rxn_out_fname=\'rcwin_' + str(thread.history.window) + '_' + str(thread.history.index) + '_us.dat\',\n')
                file.write('  rxn_out_frq=1,\n')
                file.write(' &end\n')
                file.write(' &rxncor_order_parameters\n')
                file.write('  alp0=' + str(alp0) + ',\n')
                ordinal = 1
                for term in condensed_rc:
                    # First, obtain the type of CV (optype) and the number of atoms involved (nat)
                    cv_index = int(re.findall('CV[0-9]+', term)[0].replace('CV', ''))
                    atoms, optype, nat = utilities.interpret_cv(cv_index, settings)  # get atom indices and type for this CV

                    # Next, obtain the value of "alp" for this CV (coefficient / (max - min))
                    coeff = term.replace('CV' + str(cv_index), '').replace('*','')
                    try:
                        null = float(coeff)
                    except ValueError:
                        raise RuntimeError('unable to cast coefficient of CV' + str(cv_index) + ' to float. It must be '
                                           'specified in a non-standard way. Offending term is: ' + term)
                    if settings.rc_reduced_cvs:
                        this_min = rc_minmax[0][cv_index - 1]
                        this_max = rc_minmax[1][cv_index - 1]

                        if optype in ['angle', 'dihedral']:     # convert from angles to radians for irxncor
                            this_min = this_min * numpy.pi / 180
                            this_max = this_max * numpy.pi / 180
                    else:   # to effectively turn off reduction of variables, we set...
                        this_min = 0
                        this_max = 1

                    alp = float(coeff)/(this_max - this_min)

                    # Finally, write it out and increment ordinal
                    file.write('\n  optype(' + str(ordinal) + ')=\'' + optype + '\',\n')
                    file.write('  alp(' + str(ordinal) + ')=' + str(alp) + ',\n')
                    file.write('  factnorm(' + str(ordinal) + ')=1.0,\n')
                    file.write('  offnorm(' + str(ordinal) + ')=' + str(this_min) + ',\n')
                    file.write('  nat(' + str(ordinal) + ')=' + str(nat) + ',\n')
                    if not optype == 'diffdistance':
                        file.write('  nat1(' + str(ordinal) + ')=' + str(nat) + ',\n')
                        for nat_index in range(nat):
                            at = str(atoms[nat_index])
                            file.write('  at(' + str(nat_index + 1) + ',' + str(ordinal) + ')=' + str(at) + ',\n')
                    else:
                        file.write('  nat1(' + str(ordinal) + ')=2,\n')
                        for nat_index in [0, 1]:
                            at = str(atoms[nat_index])
                            file.write('  at(' + str(nat_index + 1) + ',' + str(ordinal) + ')=' + str(at) + ',\n')
                        file.write('  nat2(' + str(ordinal) + ')=2,\n')
                        for nat_index in [2, 3]:
                            at = str(atoms[nat_index])
                            file.write('  at(' + str(nat_index + 1) + ',' + str(ordinal) + ')=' + str(at) + ',\n')

                    ordinal += 1

                file.write(' &end\n')

        elif settings.us_implementation.lower() == 'plumed':
            # We want to add specify a plumedfile in the Amber input file, and then make that file
            # First, deal with the Amber input file.
            # Check to ensure that the string 'plumed=1' appears in the input file
            if not True in ['plumed=1' in line.lower().replace(' ','') for line in open(settings.path_to_input_files + '/umbrella_sampling_prod_' + settings.md_engine.lower() + '.in', 'r').readlines()]:
                raise RuntimeError('Required option plumed=1 not found in the Amber input file: ' +
                                   settings.path_to_input_files + '/umbrella_sampling_prod_' + settings.md_engine.lower() + '.in')

            # Then replace the template slot we asked the user to prepare with the appropriate plumedfile name
            plumedfile = 'plumed_' + str(thread.history.window) + '_' + str(thread.history.index) + '.in'
            lines = open(settings.path_to_input_files + '/umbrella_sampling_prod_' + settings.md_engine.lower() + '.in', 'r').readlines()
            with open(settings.working_directory + '/' + input_file_name, 'w') as f:
                for line in lines:
                    if '{{ plumedfile }}' in line:
                        line = line.replace('{{ plumedfile }}', '\'' + plumedfile + '\'')
                    f.write(line)

            # Next, build the plumedfile
            if not os.path.exists(plumedfile):
                with open(plumedfile, 'w') as f:
                    f.write('# PLUMED file to define restraints for umbrella sampling centered at: ' + str(thread.history.window) + '\n')
                    f.write('UNITS LENGTH=A ENERGY=kcal/mol\n')
                    # Iterate through each term in the RC to add definitions for them to the plumed file
                    labels = []
                    RCfunc = str(alp0)
                    for term in condensed_rc:
                        # First, obtain the type of CV (optype) and the number of atoms involved (nat)
                        cv_index = int(re.findall('CV[0-9]+', term)[0].replace('CV', ''))
                        atoms, optype, nat = utilities.interpret_cv(cv_index, settings)  # get atom indices and type for this CV

                        coeff = term.replace('CV' + str(cv_index), '').replace('*', '')
                        try:
                            null = float(coeff)
                        except ValueError:
                            raise RuntimeError('unable to cast coefficient of CV' + str(cv_index) + ' to float. It must be '
                                               'specified in a non-standard way. Offending term is: ' + term)

                        if settings.rc_reduced_cvs:
                            this_min = rc_minmax[0][cv_index - 1]
                            this_max = rc_minmax[1][cv_index - 1]

                        else:  # to effectively turn off reduction of variables, we set...
                            this_min = 0
                            this_max = 1

                        alp = float(coeff) / (this_max - this_min)

                        # 'distance', 'angle', 'dihedral', or 'diffdistance'
                        if optype == 'distance':
                            f.write('DISTANCE LABEL=CV' + str(cv_index) + ' ATOMS=' + str(atoms[0]) + ',' + str(atoms[1]) + '\n')
                            labels.append('CV' + str(cv_index))
                            cv_str = '(cv' + str(cv_index) + '-' + str(this_min) + ')'
                        elif optype == 'angle':
                            f.write('ANGLE LABEL=CV' + str(cv_index) + ' ATOMS=' + str(atoms[0]) + ',' + str(atoms[1]) + ',' + str(atoms[2]) + '\n')
                            labels.append('CV' + str(cv_index))
                            cv_str = '((cv' + str(cv_index) + '*180/(pi))' + '-' + str(this_min) + ')'
                        elif optype == 'dihedral':
                            f.write('TORSION LABEL=CV' + str(cv_index) + ' ATOMS=' + str(atoms[0]) + ',' + str(
                                atoms[1]) + ',' + str(atoms[2]) + ',' + str(atoms[3]) + '\n')
                            labels.append('CV' + str(cv_index))
                            cv_str = '((cv' + str(cv_index) + '*180/(pi))' + '-' + str(this_min) + ')'
                        elif optype == 'diffdistance':
                            f.write('DISTANCE LABEL=CV' + str(cv_index) + 'A ATOMS=' + str(atoms[0]) + ',' + str(atoms[1]) + '\n')
                            f.write('DISTANCE LABEL=CV' + str(cv_index) + 'B ATOMS=' + str(atoms[2]) + ',' + str(atoms[3]) + '\n')
                            labels.append('CV' + str(cv_index) + 'A')
                            labels.append('CV' + str(cv_index) + 'B')
                            cv_str = '((cv' + str(cv_index) + 'a-cv' + str(cv_index) + 'b)' + '-' + str(this_min) + ')'
                        else:
                            raise RuntimeError('unrecognized CV type: ' + optype)

                        # Construct RC function
                        RCfunc += '+' + str(alp) + '*' + cv_str

                    f.write('CUSTOM ...\n')
                    f.write('  LABEL=RC\n')
                    f.write('  ARG=' + ','.join(labels) + '\n')
                    f.write('  VAR=' + ','.join([item.lower() for item in labels]) + '\n')
                    f.write('  FUNC=' + RCfunc + '\n')
                    f.write('  PERIODIC=NO\n')
                    f.write('... CUSTOM\n')
                    f.write('restraint-rc: RESTRAINT ARG=RC KAPPA=' + str(2 * float(settings.us_restraint)) + ' AT=' + str(thread.history.window) + '\n')
                    f.write('PRINT ARG=RC FILE=rcwin_' + str(thread.history.window) + '_' + str(thread.history.index) + '_us.dat')

        else:
            raise RuntimeError('unrecognized us_implementation option: ' + settings.us_implementation)

        with open(settings.working_directory + '/' + input_file_name, 'a') as file:
            if settings.us_pathway_restraints_file:
                if True in ['type="end"' in line.lower() for line in open(settings.path_to_input_files + '/umbrella_sampling_prod_' + settings.md_engine.lower() + '.in', 'r').readlines()]:
                    raise RuntimeError('The umbrella sampling input file ' + settings.path_to_input_files +
                                       '/umbrella_sampling_prod_' + settings.md_engine.lower() + '.in appears to contain an'
                                       ' &wt namelist with \'type="END"\', which must be added by ATESA when using '
                                       'the us_pathway_restraints_file option. Please remove it and try again.')
                file.write(' &wt\n  type="END",\n &end\n')
                file.write('DISANG=us_pathway_restraints_' + str(thread.history.window) + '.DISANG')
            else:
                if not True in ['type="end"' in line.lower() for line in open(settings.path_to_input_files + '/umbrella_sampling_prod_' + settings.md_engine.lower() + '.in', 'r').readlines()]:
                    file.write(' &wt\n  type="END",\n &end\n')
                    if not settings.suppress_us_warning:
                        print('Did not find an &wt namelist with \'type="END"\' in the umbrella sampling input file' +
                              settings.path_to_input_files + '/umbrella_sampling_prod_' + settings.md_engine.lower() + '.in, so'
                              ' ATESA added it automatically.')
                        settings.suppress_us_warning = True     # so that this is only printed once

        return input_file_name

    def get_input_file_equilibrium_path_sampling(self, settings, thread, job_index, **kwargs):
        if not settings.eps_n_steps % settings.eps_out_freq == 0:
            raise RuntimeError('eps_n_steps must be evenly divisible by eps_out_freq')

        if thread.current_type == ['init']:
            return settings.path_to_input_files + '/' + settings.job_type + '_' + thread.current_type[
                job_index] + '_' + settings.md_engine.lower() + '.in'
        else:
            if job_index == 0:  # have to roll to determine the number of fwd and bwd steps
                roll = random.randint(1, int(settings.eps_n_steps / settings.eps_out_freq)) * settings.eps_out_freq
                thread.history.prod_lens.append([roll, settings.eps_n_steps - roll])

            input_file_name = 'eps_' + str(thread.history.prod_lens[-1][job_index]) + '.in'

            if not os.path.exists(input_file_name):
                template = settings.env.get_template(settings.path_to_input_files + '/' + settings.job_type + '_' +
                                                     thread.current_type[job_index] + '_' + settings.md_engine.lower()
                                                     + '.in')
                filled = template.render(django.template.Context(
                    {'nstlim': str(thread.history.prod_lens[-1][job_index]), 'ntwx': str(settings.eps_out_freq)}))
                with open(input_file_name, 'w') as newfile:
                    newfile.write(filled)
                    newfile.close()

            return input_file_name

    def get_input_file_find_ts(self, settings, thread, job_index, **kwargs):
        # In find TS, not only do we want to get the input filename, we also want to store the initial basin and write
        # the restraint file to push it into the other basin. First, get the basin...
        commit = utilities.check_commit(thread.history.prod_inpcrd[0], settings)
        if commit == '':
            raise RuntimeError('the coordinates provided during a find_ts run must represent a structure in either the '
                               'fwd or bwd basin, but the coordinate file ' + thread.history.prod_inpcrd[0] + ' is '
                               'in neither.\nIf it is a transition state guess, you should set jobtype = '
                               '\'aimless_shooting\' instead to begin aimless shooting.')
        elif commit in ['fwd', 'bwd']:
            thread.history.init_basin = commit
        else:
            raise RuntimeError('internal error in utilities.check_commit(); did not return valid output with coordinate'
                               ' file: ' + thread.history.prod_inpcrd[0])

        def write_find_ts_restraint(basin, inp_file):
            # Writing an amber-style restraint file from scratch
            open('find_ts_restraints.disang', 'w').write('')  # initialize file
            with open('find_ts_restraints.disang', 'a') as f:
                f.write('DISANG restraint file produced by ATESA with find_ts > 0 to produce transition state guesses\n')
                for def_index in range(len(basin[0])):
                    extra = 0  # additional distance to add to basin definition to push *into* basin rather than to its edge
                    if basin[3][def_index] == 'lt':
                        extra = -0.1 * basin[2][def_index]  # minus 10% # todo: this doesn't work for negative angles/dihedrals (moves them towards zero instead of "more negative"), which is fine for now since angles and dihedrals are not yet supported
                    elif basin[3][def_index] == 'gt':
                        extra = 0.1 * basin[2][def_index]   # plus 10%
                    else:
                        raise RuntimeError('entries in the last list in commitment definitions (commit_fwd and commit_bwd) '
                                           'must be either \'lt\' (less than) or \'gt\' (greater than)')

                    f.write(' &rst\n')
                    f.write('  iat=' + str(basin[0][def_index]) + ',' + str(basin[1][def_index]) + ',\n')
                    f.write('  r1=0, r2=' + str(basin[2][def_index] + extra) + \
                            ', r3=' + str(basin[2][def_index] + extra) + \
                            ', r4=' + str(basin[2][def_index] + extra + 2) + ',\n')
                    f.write('  rk2=500, rk3=500,\n')
                    f.write(' &end\n')

            return inp_file     # unmodified, as find_ts_amber already includes the reference to find_ts_restraints.disang

        # Then, write the restraint file
        if thread.history.init_basin == 'fwd':
            other_basin = settings.commit_bwd
        else:  # == 'bwd', the only other valid option
            other_basin = settings.commit_fwd
        input_file = write_find_ts_restraint(other_basin, settings.path_to_input_files + '/find_ts_prod_' + settings.md_engine.lower() + '.in')

        return input_file

    def control_rst_format(self, restart_file):
        return restart_file     # Amber restart files already match the desired format by default

    def simulation_cleanup(self, thread_history, settings):
        pass    # don't need to do anything for Amber


class AdaptCP2K(MDEngine):
    """
    Adapter class for CP2K MDEngine.

    """

    def get_frame(self, trajectory, frame, settings):
        new_restart_name = trajectory + '_frame_' + str(frame) + '.rst7'
        if not os.path.exists(trajectory):
            return ''  # since it's possible to call this before the trajectory file has been initialized
        if frame >= 1:
            shift_frame = frame - 1  # because write_traj is 0-indexed but get_frame is 1-indexed
        elif frame == -1:
            shift_frame = -1
        else:
            raise IndexError('invalid frame index for get_frame: ' + str(frame) + ' (must be >= 1, or exactly -1)')

        traj = mdtraj.load(trajectory, top=settings.topology)
        if traj.n_frames == 0:
            return ''

        try:
            traj[shift_frame].save_amberrst7(new_restart_name, force_overwrite=True)
        except IndexError:
            raise IndexError('frame index ' + str(frame) + ' is out of range for trajectory: ' + trajectory)
        if not os.path.exists(new_restart_name):
            raise OSError(
                'expected mdtraj to write ' + new_restart_name + ' but it was not found.\nThe most likely '
                'explanation is that the ATESA process is unable to write to the working directory (' +
                settings.working_directory + '). You may have run out of storage space.')

        # ATESA expects velocity information in the .rst7 file but CP2K doesn't have it

        return new_restart_name

    def get_input_file_aimless_shooting(self, settings, thread, job_index, **kwargs):
        templ_file = settings.path_to_input_files + '/' + settings.job_type + '_' + thread.current_type[job_index] + '_' + settings.md_engine.lower() + '.in'
        template = settings.env.get_template(templ_file)

        # As appropriate, need to extract box information and velocities from inpcrd and add them to the kwargs
        mtraj = mdtraj.load(kwargs['inpcrd'], top=kwargs['prmtop'])
        box_xyz = ' '.join([str(item * 10) for item in mtraj.unitcell_lengths[0]])    # * 10 to convert nm to Å
        box_abc = ' '.join([str(item) for item in mtraj.unitcell_angles[0]])

        # Assume that the input coordinate file contains velocities
        # Velocities begin on line ceil(0.5 * mtraj.n_atoms) + first_coord_line and end on the penultimate non-empty line
        velocities = ''
        lines = open(kwargs['inpcrd'], 'r').readlines()
        # First, identify first line containing coordinates
        first_coord_line = 0
        for line in lines:
            if len(line.split()) == 6:  # first coord line is always six values
                try:    # ensure they're all numbers
                    null = [float(item) for item in line.split()]
                except TypeError:
                    continue
                break   # we're done, first_coord_line is the first coord line
            first_coord_line += 1
        # Now identify the last line containing velocities
        last_vel_line = -1
        while len(lines[last_vel_line].split()) == 0:   # ignore all blank lines at the end
            last_vel_line -= 1
        last_vel_line -= 1  # since it's the penultimate line
        for ii in range((math.ceil(0.5 * mtraj.n_atoms) + first_coord_line), len(lines) + last_vel_line + 1):
            velocities += ' '.join(lines[ii].split()[0:3]) + '\n' + ' '.join(lines[ii].split()[3:6]) + '\n'
        velocities = velocities[:-1]    # remove trailing newline

        # cast kwargs['inpcrd'] to .pdb file? As an alternative to writing elements into input files directly
        mtraj.save_pdb(kwargs['inpcrd'] + '.pdb', force_overwrite=True)
        kwargs['inpcrd'] = kwargs['inpcrd'] + '.pdb'

        kwargs.update({ 'box_xyz': box_xyz, 'box_abc': box_abc, 'velocities': velocities})

        filled = template.render(django.template.Context(kwargs))
        newfilename = kwargs['name'] + '.inp'
        with open(newfilename, 'w') as newfile:
            newfile.write(filled)
            newfile.close()
        return newfilename

    def get_input_file_committor_analysis(self, settings, thread, job_index, **kwargs):
        pass

    def get_input_file_umbrella_sampling(self, settings, thread, job_index, **kwargs):
        pass

    def get_input_file_equilibrium_path_sampling(self, settings, thread, job_index, **kwargs):
        pass

    def get_input_file_find_ts(self, settings, thread, job_index, **kwargs):
        templ_file = settings.path_to_input_files + '/' + settings.job_type + '_' + thread.current_type[job_index] + '_' + settings.md_engine.lower() + '.in'
        template = settings.env.get_template(templ_file)

        # As appropriate, need to extract box information and velocities from inpcrd and add them to the kwargs
        mtraj = mdtraj.load(kwargs['inpcrd'], top=kwargs['prmtop'])
        box_xyz = ' '.join([str(item * 10) for item in mtraj.unitcell_lengths[0]])  # * 10 to convert nm to Å
        box_abc = ' '.join([str(item) for item in mtraj.unitcell_angles[0]])

        # Need to define template fillers for colvars (definitions of dimensions to constrain) and collective
        # (definitions of what the constraints actually are)
        if thread.history.init_basin == 'fwd':
            other_basin = settings.commit_bwd
        else:  # == 'bwd', the only other valid option
            other_basin = settings.commit_fwd
        colvars = ''
        collective = ''
        for def_index in range(len(other_basin[0])):
            extra = 0  # additional distance to add to basin definition to push *into* basin rather than to its edge
            if other_basin[3][def_index] == 'lt':
                extra = -0.1 * other_basin[2][def_index]  # minus 10% # todo: this doesn't work for negative angles/dihedrals (moves them towards zero instead of "more negative"), which is fine for now since angles and dihedrals are not yet supported
            elif other_basin[3][def_index] == 'gt':
                extra = 0.1 * other_basin[2][def_index]  # plus 10%
            else:
                raise RuntimeError('entries in the last list in commitment definitions (commit_fwd and commit_bwd) '
                                   'must be either \'lt\' (less than) or \'gt\' (greater than)')

            init_value = mdtraj.compute_distances(mtraj, numpy.array([[other_basin[0][def_index] - 1, other_basin[1][def_index] - 1]]))[0][0] * 10

            # todo: find a solution to Django escaping ampersands that isn't just adding "|safe" in the template files

            colvars += '&COLVAR\n'
            colvars += '  &DISTANCE\n'
            colvars += '    ATOMS ' + str(other_basin[0][def_index]) + ' ' + str(other_basin[1][def_index]) + '\n'
            colvars += '  &END DISTANCE\n'
            colvars += '&END COLVAR\n'

            collective += '&COLLECTIVE\n'
            collective += '  INTERMOLECULAR TRUE\n'
            collective += '  COLVAR ' + str(def_index + 1) + '\n'
            collective += '  TARGET [angstrom] ' + str(init_value) + '\n'
            collective += '  TARGET_LIMIT [angstrom] ' + str(other_basin[2][def_index] + extra) + '\n'
            collective += '  TARGET_GROWTH [angstrom/fs] ' + str(((other_basin[2][def_index] + extra) - init_value)/100) + '\n'
            collective += '&END COLLECTIVE\n'

        # Remove trailing newlines
        colvars = colvars[:-1]
        collective = collective[:-1]

        # cast kwargs['inpcrd'] to .pdb file? As an alternative to writing elements into input files directly
        mtraj.save_pdb(kwargs['inpcrd'] + '.pdb', force_overwrite=True)
        kwargs['inpcrd'] = kwargs['inpcrd'] + '.pdb'

        kwargs.update({'box_xyz': box_xyz, 'box_abc': box_abc, 'colvars': colvars, 'collective': collective})

        filled = template.render(django.template.Context(kwargs))
        newfilename = kwargs['name'] + '.inp'
        with open(newfilename, 'w') as newfile:
            newfile.write(filled)
            newfile.close()
        return newfilename

    def control_rst_format(self, restart_file):
        # CP2K restart files are actually input files with atomic coordinates and velocities in them. This method parses
        # those sections and then writes them directly out to a new .rst7 file from scratch.

        if os.path.exists(restart_file + '.bak'):
            os.rename(restart_file + '.bak', restart_file)

        os.rename(restart_file, restart_file + '.bak')  # make a backup of the original restart file
        open(restart_file, 'w').close()     # make a new placeholder for the .rst7 file

        lines = open(restart_file + '.bak', 'r').readlines()
        coords = []
        velocs = []
        cell = [0, 0, 0]
        coords_yet = False
        velocs_yet = False
        cell_yet = False
        for line in lines:
            if '&COORD' in line:
                coords_yet = True
                continue
            elif '&VELOCITY' in line:
                velocs_yet = True
                continue
            elif '&END COORD' in line:
                coords_yet = False
                continue
            elif '&END VELOCITY' in line:
                velocs_yet = False
                continue
            elif '&CELL' in line:
                cell_yet = True
                continue
            elif '&END CELL' in line:
                cell_yet = False
                continue
            if coords_yet:
                coords.append(['%.7f' % float(item) for item in line.split()[1:4]])
            elif velocs_yet:
                velocs.append(['%.7f' % float(item) for item in line.split()[0:3]])
            elif cell_yet:
                if line.split()[0] == 'A':
                    cell[0] = '%.7f' % float(line.split()[1])
                if line.split()[0] == 'B':
                    cell[1] = '%.7f' % float(line.split()[2])
                if line.split()[0] == 'C':
                    cell[2] = '%.7f' % float(line.split()[3])


        try:
            assert len(coords) == len(velocs)
        except AssertionError:
            raise RuntimeError('Error when parsing CP2K restart file: ' + restart_file + '.bak\n'
                               'Mismatch in number of atoms and number of atomic velocities')

        with open(restart_file, 'a') as f:
            f.write('.rst7 format restart file built by ATESA based on contents of ' + restart_file + '.bak\n')
            f.write(str(len(coords)) + '\n')
            for ii in range(len(coords)):
                f.write(''.join([string.rjust(12) for string in coords[ii]]))
                if ii % 2 == 1 or ii + 1 == len(coords):
                    f.write('\n')
            for ii in range(len(velocs)):
                f.write(''.join([string.rjust(12) for string in velocs[ii]]))
                if ii % 2 == 1 or ii + 1 == len(velocs):
                    f.write('\n')
            f.write(''.join([string.rjust(12) for string in cell]))
            f.write('  90.0000000  90.0000000  90.0000000\n')   # todo: any way to support other box geometries?

    def simulation_cleanup(self, thread_history, settings):
        # todo: this whole concept is pretty kludgey, can I formalize this better?

        if type(thread_history.prod_trajs[-1]) == str:
            iterable = [thread_history.prod_trajs[-1]]
        elif type(thread_history.prod_trajs[-1]) == list:
            iterable = thread_history.prod_trajs[-1]
        else:
            raise RuntimeError('Unexpected data type for thread_history.prod_trajs[-1]: ' + str(type(thread_history.prod_trajs[-1])))

        for filename in iterable:
            if not os.path.exists(filename):
                files = glob.glob('*-' + filename + '-*')
                if len(files) == 1:
                    real_name = files[0]
                else:
                    raise RuntimeError(
                        'Could not identify file: ' + filename + '\nIt does not exist and there is no single'
                        ' unique file matching the pattern \'*-' + filename + '-*\'.')
            else:
                continue
            os.rename(real_name, filename)