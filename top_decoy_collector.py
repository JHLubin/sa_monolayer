#!/usr/bin/python
"""
Selects top decoys from folders made by submit_script_gen.py on MARCC for 
LK peptide SAM surface docks. If solution state decoys are being analyzed, a 
new score list will be generated.

Created by Joseph Lubin
"""

import argparse
import os
from shutil import copyfile, move
# Further import of PyRosetta is only necessary to analyze solution models


def parse_args():
	""" Collecting inputs """
	parser = argparse.ArgumentParser()
	parser.add_argument("folder", help="What folder below the current \
		 						working directory do you want to use?")
	parser.add_argument("-sol", "--solution_state", action="store_true", 
						help="Do ou want the solution state decoys, rather \
								than the surface docked ones?")
	parser.add_argument("-d", "--decoy_count", type=int, default=100, 
							help="Collect the top [how many] decoys?")
	parser.add_argument("-l", "--location", type=str, choices=['m', 'l'], 
						default='m', help="Are PDB's on MARCC (m, default) \
											or downloaded to local (l)?")
	parser.add_argument("-s", "--silence", action="store_true", 
						help="Silence output from calculations")
	args = parser.parse_args()

	return args


def get_folders(direc, f_name, args):
	""" 
	This function returns an input folder with PDB files and an output folder
	to which top decoys, a score file, and the flags file will be copied
	"""
	# Getting path for folder containing decoys
	if args.location == 'm': 	# Running from MARCC-formatted folders
		input_folder = os.path.join(direc, f_name + '_output')
	else: 	# Running from folders already created by this script
		input_folder = direc

	# Making top decoys folder
	td_folder_name = f_name + '.top' + str(args.decoy_count)
	top_decoys_folder = os.path.join(input_folder, td_folder_name)
	if not os.path.isdir(top_decoys_folder):
		os.makedirs(top_decoys_folder)

	return input_folder, [td_folder_name, top_decoys_folder]


def display_time(start, elapsed, run_count, total_count):
	""" Displays the time for scoring PDBs """
	# Giving an initial estimate of computation time
	if run_count == 1: 
		total_time = elapsed / 2 * total_count / 3600
		out = "Rough estimated time to score this folder: {} hours"
		print out.format(round(total_time,2))
		return

	# Giving a total at completion
	if (run_count + 1) == total_count:
		out = "Scoring completed for this folder. Completion time: {} hours"
		print out.format(round(elapsed,2))
		return

	# Updating estimate every 100 decoys
	if (run_count + 1) % 100 == 0:
		est = elapsed * (float(total_count / (run_count + 1)) - 1) / 3600
		out = "\t{} PDBs scored. Estimated time remaining: {} hours"
		print out.format(run_count + 1, round(est, 2))
		return


def rescore_sol(folder, args):
	"""
	This function takes an input folder and scores all the solution state
	decoys in the folder. It returns a sorted list of names and scores.
	This requires the importation of PyRosetta. There is a limitation that
	this function will only use the default Rosetta score function. In this
	case, it is using talaris2014.
	"""
	import rosetta 
	import rosetta.core.scoring.solid_surface
	opts = '-include_surfaces -mute basic -mute core -mute protocols'
	rosetta.init(extra_options = opts)
	if not args.silence:
		from time import time

	# Getting list of all files in the folder
	f_name = os.path.basename(folder).replace('_output', '')
	if not args.silence:
		print '\n\nFolder:\t{}'.format(f_name)

	folder_list = os.listdir(folder) 	# Full folder

	# Narrowing list to only solution models
	sol_pdbs = []
	for i in folder_list:
		if 'Sol' in i:
			sol_pdbs.append(i)
	sol_pdbs.sort()
	count = len(sol_pdbs)

	if not args.silence:
		print "Scoring {} PDBs".format(count)

	# Writing unsorted scores file
	scoresc = os.path.join(folder, 'sol_score.sc')
	header = ('\t' * 6).join(['Description', 'Score'])
	with open(scoresc, 'w') as s:
		s.write(header)

	# Scoring solution PDBs and listing scores
	score_erors = {}
	sf = rosetta.get_fa_scorefxn()
	sol_scores = []
	start = time()
	for i in range(len(sol_pdbs)):
		try:
			pdb = sol_pdbs[i]
			p = rosetta.pose_from_pdb(os.path.join(folder, pdb))
			score = sf(p)
			sol_scores.append(score)

			# Adding score to unsorted list
			with open(scoresc, 'a') as s:
				s.write('\n{}\t{}'.format(pdb, score))

			if not args.silence:
				elapsed = time() - start
				display_time(start, elapsed, i, count)

		except RuntimeError:
			print "Unable to read PDB: {}".format(pdb)
			if score_erors.has_key(f_name):
				score_erors[f_name].append(pdb)
			else:
				score_erors.update({f_name: [pdb]})

	# Combining files names and scores, sorting by scores
	s_name_scores = sorted(zip(sol_pdbs, sol_scores), key=lambda x:x[1])

	return s_name_scores, header, score_erors


def read_scorefile(scorefile):
	"""
	This function uses the output score.sc file from the MARCC run to find the
	top docked decoys based on Rosetta score. It requires a score file input.
	Outputs a sorted scores list
	"""
	# Reading score file
	score_text = open(scorefile).readlines()
	score_header = score_text[1]
	raw_scores = score_text[2:]

	# Getting total scores and filenames
	total_scores = []
	for line in raw_scores:
		total_scores.append(line.split())

	# Converting string numbers to floats and sorting
	for line in range(len(total_scores)):
		# Removing "SCORE:" from each line
		total_scores[line] = total_scores[line][1:]
		for i in range(len(total_scores[line])-1):
			# Last item is file name
			float_val = float(total_scores[line][i])
			total_scores[line][i] = float_val
		
	total_scores.sort()

	# Reorganizing data
	row_size = len(total_scores[0])
		# Cleaning header
	headlist = score_header.split()[-row_size:]    # Strips all before "score"
	header = str(headlist[-1]) 	# Convert to str with file name before scores
	for column in range(len(headlist)-1):
		header += '\t' + str(headlist[column]).replace(':','')

		# Re-ordering scores to match header
	for line in range(len(total_scores)):
		total_scores[line] = [total_scores[line][-1]] + \
								total_scores[line][1:-1]

	return total_scores, header


def pdb_copy(pdb_name, source, destination, args):
	""" Copy a PDB file from one folder to another """
	if args.solution_state:
		dec = pdb_name
	else:
		dec = pdb_name + '.pdb.gz'

	decoy_source = os.path.join(source, dec)
	decoy_destination = os.path.join(destination, dec)
	copyfile(decoy_source, decoy_destination)
	
	if not args.silence:
		print '\t' + dec


def score_and_flags(name, dir_in, dir_out, score_list, header, args):
	""" Output a sorted scores file and copy over the simulation flags """
	if not args.silence:
		print 'Writing sorted scores file and copying flags\n\n'

	# Making output scores file
	out_scores = os.path.join(dir_out, 'score.sc')
	with open(out_scores, 'w') as sc:
		sc.write(header + '\n')
		for line in score_list:
			line_text = '\t'.join([str(i) for i in line])
			sc.write(line_text + '\n')

	#Copying flags file
	flag_name = name + '.flags'
	flag_source = os.path.join(dir_in, flag_name)
	flag_destination = os.path.join(dir_out, flag_name)
	copyfile(flag_source, flag_destination)


def take_top_decoys(folder, name, args):
	"""
	Given a folder containing PDB files, this function will get the 
	appropriate directory or subdirectory containing PDB files. Then it will
	either rescore solution state PDBs or read and sort the score.sc file to
	get the lowest scoring decoys. The specified number of best scorers will 
	be copied to a top decoys folder, along with the flags file that generated 
	them, and a sorted scores file.
	"""
	folder_full = os.path.join(folder, name)
	input_folder, out_folder = get_folders(folder_full, name, args)

	# Getting scores list
	if args.solution_state: 	# For solution models, need to re-score
		scores, header, errors = rescore_sol(input_folder, args)

	else: 	# Otherwise, can read scores directly from score.sc
		if not args.silence:
			print 'Folder:\t' + name

		scorefile = os.path.join(input_folder, 'score.sc')
		scores, header = read_scorefile(scorefile)
		errors = {}

	# Getting top models
	top_models = []
	for i in range(args.decoy_count):
		top_models.append(scores[i][0])

	# Copying top decoys
	print '\nCopying top {} decoys'.format(args.decoy_count)
	for decoy in top_models:
		pdb_copy(decoy, input_folder, out_folder[1], args)
	
	# Making output scores file and copying over flags for reference
	score_and_flags(name, folder_full, out_folder[1], 
						scores, header, args)

	return out_folder, errors


def main():
	args = parse_args()

	# Setting path
	main_directory = os.getcwd() + "/" + args.folder

	# Getting list of subfolders
	subdirectories = next(os.walk(main_directory))[1]
	subdirectories.sort()

	# Picking out top decoys in each subfolder
	out_folders = {}
	pdb_exclusions = {}
	folder_exclusions = []

	for folder in subdirectories:
		#verify that folder has PDBs (exclude results folders, etc.)
		run = False
		for file in os.listdir(os.path.join(main_directory, folder)):
			if any(x in file for x in ['.pdb.gz', '_output']):
				run = True
				break

		if run:
			out_folder, errors = take_top_decoys(main_directory, folder, args)
			out_folders.update({out_folder[0]: out_folder[1]})
			pdb_exclusions.update(errors)
			
		else:
			print 'Excluding:\t\t' + folder
			folder_exclusions.append(folder)

	# Making results directory
	results_folder = os.path.join(main_directory, 'results')
	if not os.path.isdir(results_folder):
		os.makedirs(results_folder)
	for folder in out_folders:
		move(out_folders[folder], os.path.join(results_folder, folder))

	# Displaying results directory to copy-paste for scp/rsync
	print '\n\n' + results_folder + '/*'

	# Displaying excluded folders and PDBs
	if len(folder_exclusions) > 0:
		print '\nFOLDERS EXCLUDED:'
		for exclusion in folder_exclusions:
			print exclusion

	if len(pdb_exclusions) > 0:
		print '\nPDBS EXCLUDED:'
		for folder in pdb_exclusions:
			print folder
			for pdb in pdb_exclusions[folder]:
				print '\t' + pdb

if __name__ == '__main__':
	main()