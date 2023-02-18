import os
import requests
import time
import bs4
from bs4 import BeautifulSoup
import concurrent.futures
import pandas as pd
import numpy as np
import jpredapi
from time import sleep
import time
import warnings
import uuid
from requests_toolbelt import MultipartEncoder
from decimal import Decimal
import argparse

# Algos used:
# (gor1, dpm, gor3, phd, predator, hnn, mlrc, sopm, dsc) via PRABI, (JPred) via JPred, (prof, sspro, yaspin, JNet, PSIPRED, sympred) via Sympred

# Use conda environment called 'PyVkabat'
# $ conda activate PyVkabat

################################ Configuration Area ################################ 
# all
#seq_name = "test_sequence_1"
#sequence = "NYLDLFSHKNMKLKERVLIPVKQYPKFNFVGKILGPQGNTIKRLQEETGAKISVLGKGSMRDKAKEEELRKGGDPKYAHLNMDLHVFIEVFGPPCEAYALMAHAMCEVKKFLVPDMMD"
alignment_width = "100"

# gor1
constants = "0"
dch = "40"
dce = "35"
dct = "0"
dcc = "0"

# predator
use_dssp_or_stride = "dssp" # "dssp" or "stride"

# sopm
states = "3" # "3" is Helix, Sheet, Coil and "4" is Helix, Sheet, Turn, Coil
threshold = "8" # similarity threshold (default is 8)
width = "17" # window width (default is 17)

# JPred
email = ''
skipPDB = True
jpred_mode = "single" # options are single, batch, and msa
jpred_user_format = "raw" # options are raw, fasta, msf, blc
jpred_seq_or_file = "sequence" # options are sequence or file
jpred_input_file = ''
check_jpred_status = int(60) # Number of minutes to keep checking the status before giving up.

#############################################################################

def parse_arguments():
	
	parser = argparse.ArgumentParser(description='Calculates secondary structure variability (vkabat = k * N / n1) for each resiue in a protein sequence.')
	
	# positional arguments
	parser.add_argument('sequence', metavar='<sequence>', type=str, help='Enter sequence (using 1 letter amino acid abbreviations)')
	
	# optional arguments
	parser.add_argument('--name', metavar='<protein name>', type=str, help='Enter the name of the protein sequence or job name. This is a way to keep track of the submission and will affect the csv file names deposited in the project folder.')
	parser.add_argument('--dir', metavar='<output directory path>', type=str, help='Enter the path to the output directory. This is where the output files will be sent.')
	parser.add_argument('--email', metavar='<email>', type=str, help='Enter your email. Some servers will email you the results.')
	parser.add_argument('--jpred_timeout', metavar='<time>', type=int, help='Enter the maximum allowable time for JPred data retreival in seconds.')
	parser.add_argument('--yaspin_timeout', metavar='<time>', type=int, help='Enter the maximum allowable time for Yaspin data retreival in seconds.')
	parser.add_argument('--sympred_timeout', metavar='<time>', type=int, help='Enter the maximum allowable time for Sympred data retreival in seconds.')

	args = parser.parse_args()

	# seq_name
	global seq_name
	if args.name == None:
		seq_name = "test"
	else:
		seq_name = args.name
	print(f'Sequence Name: {seq_name}')

	# sequence
	global sequence
	sequence = args.sequence
	print(f'Input Sequence: {sequence}')
	print(f'Sequence length: {len(sequence)}')

	# jpred_timeout
	global jpred_timeout
	if args.jpred_timeout == None:
		jpred_timeout = 1300 # 30 minutes (1300 seconds)
	else:
		jpred_timeout = args.jpred_timeout
	print(f'Maximum time allowed for JPred: {jpred_timeout} seconds.')

	# yaspin_timeout
	global yaspin_timeout
	if args.yaspin_timeout == None:
		yaspin_timeout = 1300 # 30 minutes (1300 seconds)
	else:
		yaspin_timeout = args.yaspin_timeout
	print(f'Maximum time allowed for Yaspin: {yaspin_timeout} seconds.')

	# sympred_timeout
	global sympred_timeout
	if args.sympred_timeout == None:
		sympred_timeout = 1300 # 30 minutes (1300 seconds)
	else:
		sympred_timeout = args.sympred_timeout
	print(f'Maximum time allowed for Sympred: {sympred_timeout} seconds.')

	# output_directory
	global output_directory
	if args.dir == None:
		output_directory = os.getcwd()
	else:
		output_directory = args.dir
	print(f'Output Directory: {output_directory}')

def runPrabi():
	prabi_start_time = time.time()
	print('Running PRABI')
	
	prabi_algos = {
					'gor1':{
							'url':'https://npsa-prabi.ibcp.fr/cgi-bin/secpred_gor.pl', 
							'data':{
									'title':str(seq_name), 
									'notice':str(sequence), 
									'ali_width':str(alignment_width), 
									'constants':str(constants), 
									'dch':str(dch), 
									'dce':str(dce), 
									'dct':str(dct), 
									'dcc':str(dcc)}},
					
					'gor3':{
							'url':'https://npsa-prabi.ibcp.fr/cgi-bin/secpred_gib.pl',
							'data':{
									'title':str(seq_name), 
									'notice':str(sequence), 
									'ali_width':str(alignment_width)}},
					
					'dpm':{
							'url':'https://npsa-prabi.ibcp.fr/cgi-bin/secpred_dpm.pl',
							'data':{
									'title':str(seq_name), 
									'notice':str(sequence), 
									'ali_width':str(alignment_width)}},
					
					'predator':{
								'url':'https://npsa-prabi.ibcp.fr/cgi-bin/secpred_preda.pl',
								'data':{
										'title':str(seq_name), 
										'notice':str(sequence), 
										'ali_width':str(alignment_width),
										'predatorssmat':str(use_dssp_or_stride)}},
					
					'hnn':{
							'url':'https://npsa-prabi.ibcp.fr/cgi-bin/secpred_hnn.pl',
							'data':{
									'title':str(seq_name), 
									'notice':str(sequence), 
									'ali_width':str(alignment_width)}},
					
					'sopm':{
							'url':'https://npsa-prabi.ibcp.fr/cgi-bin/secpred_sopm.pl',
							'data':{
									'title':str(seq_name), 
									'notice':str(sequence), 
									'ali_width':str(alignment_width),
									'states':str(states),
									'threshold':str(threshold),
									'width':str(width)}},
									
					'mlrc':{
							'url':'https://npsa-prabi.ibcp.fr/cgi-bin/secpred_mlr.pl',
							'data':{
									'title':str(seq_name), 
									'notice':str(sequence), 
									'ali_width':str(alignment_width)}},
								
					'dsc':{
							'url':'https://npsa-prabi.ibcp.fr/cgi-bin/secpred_dsc.pl',
							'data':{
									'title':str(seq_name), 
									'notice':str(sequence), 
									'ali_width':str(alignment_width)}}
					
			}
			
			
	# Define a function to run each algorithm in a thread
	def run_algorithm(key):
		print(f'Submitting request to PRABI for {key}')
		response = requests.post(prabi_algos[key]['url'], data=prabi_algos[key]['data'])
		soup = BeautifulSoup(response.text, features='html.parser')

		output_list = []
		for x in soup.find('code').findChildren('font'):
			output_list.append(x.text)

		modified_output_list = []
		for letter in output_list:
			letter = letter.upper()
			if letter == 'T':
				modified_output_list.append(letter.replace('T', 'C'))
			else:
				modified_output_list.append(letter)

		print(f'{key}: {modified_output_list}')
		print(f'{key} completed in {time.time()-prabi_start_time} seconds')
		return {key: modified_output_list}

	# Run each algorithm in a separate thread
	output = {}
	with concurrent.futures.ThreadPoolExecutor() as executor:
		futures = [executor.submit(run_algorithm, key) for key in prabi_algos.keys()]
		for future in concurrent.futures.as_completed(futures):
			out = future.result()
			output.update(out)

	prabi_end_time = time.time()
	prabi_execution_time = prabi_end_time - prabi_start_time
	print(f'PRABI algos took {prabi_execution_time} secs to complete.')

	return output
	
def run_alt_JPred():
	alt_JPred_start_time = time.time()
	
	print('Running alt JPred (no API)')

	# Define the form data to be submitted
	data = {
		'seq': sequence,
		'fileup': ("", ""),
		"input": "seq",
		"pdb": "on",
		'email': email,
		"queryName": seq_name
	}
	
	multipart_data = MultipartEncoder(fields=data, boundary="--------------------------124327508021125478525198")
	
	# Set the headers for the request
	headers = {
		"Content-Type": multipart_data.content_type
	}
	
	# Define the URL of the JPred web form
	url = 'https://www.compbio.dundee.ac.uk/jpred/cgi-bin/jpred_form'

	# Submit the form data using a POST request
	response = requests.post(url, data=multipart_data, headers=headers)
	
	soup = BeautifulSoup(response.text, features='html.parser')
	
	output_list = []
	for x in soup.find('div', {'id':'content'}).find('a'):
		output_list.append(x.text)
	
	job_id = output_list[0].split('chklog?')[-1]

	# Check if the POST request was successful
	if response.status_code == 200:
		print(f'JPred Job ID: {job_id}')
		
		# Define the URL of the job results page
		# Get result (not using the jpredapi)
		jpred_result_base = 'http://www.compbio.dundee.ac.uk/jpred4/results'
		jpred_simple_result_url = (f'{jpred_result_base}/{job_id}/{job_id}.simple.html')

		# Wait for the job to complete and the results to become available
		print('JPred: Waiting for results...')
		
		while True:
			
			# check elapsed time
			if (time.time() - alt_JPred_start_time) > jpred_timeout:
				print(f'JPred timed out. Exceeded {jpred_timeout} seconds.')
				break
			
			# Send a GET request to the results page
			results_page = requests.get(jpred_simple_result_url)
				
			if results_page.status_code == 404:
				#print(f'Request failed with status code {results_page.status_code}')
				sleep(5)
			elif results_page.status_code == 200:
				# Get result (not using the jpredapi)

				jpred_results = requests.get(jpred_simple_result_url)
				jpred_soup = BeautifulSoup(jpred_results.text, 'html.parser')
				#print(jpred_soup)
				
				jpred_file_string = ''
				for x in jpred_soup.find('code'):
					jpred_file_string = jpred_file_string + x.text
					
				#print(jpred_file_string)
				
				string_list = jpred_file_string.split('\n')
				
				out = ''
				for i, string in enumerate(string_list):
					if (i != 0) and (i % 2 != 0):
						out = out + string
				
				out_list = []	
				for letter in out:
					if letter == '-':
						out_list.append(letter.replace('-', 'C'))
					else:
						out_list.append(letter)
						
				out = {'JPred': out_list}
				print(out)
				
				alt_JPred_end_time = time.time()
				alt_JPred_execution_time = alt_JPred_end_time - alt_JPred_start_time
				print(f'alt_JPred took {alt_JPred_execution_time} sec')
				
				return out
				break
			else:
				print(f'JPred: Request failed with status code {results_page.status_code}')
				print(f'JPred: Aborting operation.')
				return None
				break
		
	else:
		print('JPred: Unable to get job ID')
		return None

def runJPred():
	
	jpred_start_time = time.time()
	
	print('Running JPred')
			
	# Run JPred API using the Python version of the perl API
	# See https://github.com/MoseleyBioinformaticsLab/jpredapi for details
	# Some Documentation: https://jpredapi.readthedocs.io/en/latest/tutorial.html#using-jpredapi-as-a-library
	
	# submit the request to jpredapi
	print(f'Submitting request to JPred via jpredapi.')
	if jpred_seq_or_file == 'sequence':
		job = BeautifulSoup(jpredapi.submit(mode=str(jpred_mode), user_format=str(jpred_user_format), seq=str(sequence), skipPDB=skipPDB, email=str(email), silent=True).text, features='html.parser')
	elif jpred_seq_or_file == 'file':
		job = BeautifulSoup(jpredapi.submit(mode=str(jpred_mode), user_format=str(jpred_user_format), file=str(jpred_input_file), skipPDB=skipPDB, email=str(email), silent=True).text, features='html.parser')
	else:
		print('Something went wrong with the JPred submission. Check the jpred_seq_or_file variable.')
		print(f'Currently, jpred_seq_or_file = {jpred_seq_or_file}. Valid options are "sequence" or "file"')
		print('Aborting JPred.')
		return
		
	# get the job id of the submission
	link = []
	for x in job.find('a'):
		link.append(x.text)
		
	job_id = str(link[0].split('/chklog?')[-1])
	
	# Check the status of the job
	i = int(1)
	while i < (check_jpred_status + 1):
		print(f'JPred status check: {str(i)}')
		
		check = BeautifulSoup(jpredapi.status(jobid=job_id, silent=True).text, 'html.parser')
		#print(check)
		
		info = []
		for idx, line in enumerate(check):
			
			if 'Your job is 100% complete...' in line:
				print('JPred job is 100% complete. Fetching results now.')
				i = check_jpred_status + 1
				
			elif 'Results available at the following URL:' in line:
				info.append(line)
				
				status = info[0].split('.')[0].split(' ')[-1]
				
				if status == 'finished':
					print(status)
					i = check_jpred_status + 1
				
			elif 'complete...' in line:
				i += 1
				sleep(30)
			
			else:
				print('Breaking out of JPred while loop...')
				break
			
	# Get result (not using the jpredapi)
	jpred_result_base = 'http://www.compbio.dundee.ac.uk/jpred4/results'
	jpred_simple_result_url = (f'{jpred_result_base}/{job_id}/{job_id}.simple.html')

	jpred_results = requests.get(jpred_simple_result_url)
	jpred_soup = BeautifulSoup(jpred_results.text, 'html.parser')
	print(jpred_soup)
	
	jpred_file_string = ''
	for x in jpred_soup.find('code'):
		jpred_file_string = jpred_file_string + x.text
		
	#print(jpred_file_string)
	
	string_list = jpred_file_string.split('\n')
	
	out = ''
	for i, string in enumerate(string_list):
		if (i != 0) and (i % 2 != 0):
			out = out + string
	
	out_list = []	
	for letter in out:
		if letter == '-':
			out_list.append(letter.replace('-', 'C'))
		else:
			out_list.append(letter)
			
	out = {'JPred': out_list}
	print(out)
	
	jpred_end_time = time.time()
	jpred_execution_time = jpred_end_time - jpred_start_time
	print(f'JPred took {jpred_execution_time} sec')
	
	return out
	
def runSympred5():
	
	sympred_start_time = time.time()
	
	print('Running Sympred')
	
	boundary = str(uuid.uuid4())
	
	sequence_with_defline = '>abcd\n' + sequence

	# Define the form data to be submitted
	data = {
		'seq': sequence_with_defline,
		'seq_file': ('', ''),
		'email': email,
		'sympred': 'Do prediction',
		'conmethod': 'D',
		'conweight': 'N',
		'window': '21',
		'database': 'nr',
		'pred1': '-phdpsi',
		'pred2': '-prof',
		'pred3': '-sspro',
		'pred6': '-jnet',
		'pred7': '-psipred',
		'MB': '',
		'mbjob[description]': seq_name
	}
	
	# Set the headers for the request
	headers = {
		'Content-Type': 'multipart/form-data; boundary=' + boundary
	}
	
	payload = ''
	for key, value in data.items():
		payload += f'--{boundary}\r\n'
		if isinstance(value, tuple):  # File upload
			payload += f'Content-Disposition: form-data; name="{key}"; filename="{value[0]}"\r\n'
			payload += 'Content-Type: application/octet-stream\r\n\r\n'
			payload += value[1] + '\r\n'
		else:
			payload += f'Content-Disposition: form-data; name="{key}"\r\n\r\n'
			payload += value + '\r\n'

	payload += f'--{boundary}--\r\n'

	# Define the URL of the Sympred web form
	url = "https://www.ibi.vu.nl/programs/sympredwww/"

	# Submit the form data using a POST request
	response = requests.post(url, headers=headers, data=payload, allow_redirects=True)

	# Check if the POST request was successful
	if response.status_code == 202:
		# Extract the job ID from the redirect URL
		job_id = response.url.split('/')[-2]
		print(f'Sympred Job ID: {job_id}')
		
		# Define the URL of the job results page
		results_url = f'http://zeus.few.vu.nl/jobs/{job_id}/result.hpred'

		# Wait for the job to complete and the results to become available
		print('Sympred: Waiting for results...')
		
		while True:
			
			# check elapsed time
			if (time.time() - sympred_start_time) > sympred_timeout:
				print(f'Sympred timed out. Exceeded {sympred_timeout} seconds.')
				break
			
			# Send a GET request to the results page
			results_page = requests.get(results_url)
				
			if results_page.status_code == 202:
				sleep(5)
				
			elif results_page.status_code == 404:
				sleep(5)
				
			elif results_page.status_code == 200:

				text_lines = results_page.text.split('\n')
				text_lines_no_comments = []
				for line in text_lines:
					if "#" in line:
						continue
					elif line == '':
						continue
					else:
						text_lines_no_comments.append(line)
						
				AA_lines = []
				PHD_lines = []
				PROF_lines = []
				SSPRO_lines = []
				JNET_lines = []
				PSIPRED_lines = []
				#SYMPRED_lines = []
				
				for line in text_lines_no_comments:
					if 'AA    ' in line:
						AA_lines.append(line)
					elif 'PHD    ' in line:
						PHD_lines.append(line)
					elif 'PROF    ' in line:
						PROF_lines.append(line)
					elif 'SSPRO    ' in line:
						SSPRO_lines.append(line)
					elif 'JNET    ' in line:
						JNET_lines.append(line)
					elif 'PSIPRED    ' in line:
						PSIPRED_lines.append(line)
					#elif 'SYMPRED    ' in line:
						#SYMPRED_lines.append(line)
					else:
						pass
				
				space_count = 0
				for letter in AA_lines[0]:
					if letter == ' ':
						space_count += 1
				space_count += 2
				
				aa = ''
				for line in AA_lines:
					aa += line.replace(line, line[space_count:])
					
				phd = ''
				for line in PHD_lines:
					phd += line.replace(line, line[space_count:])
					
				prof = ''
				for line in PROF_lines:
					prof += line.replace(line, line[space_count:])
					
				sspro = ''
				for line in SSPRO_lines:
					sspro += line.replace(line, line[space_count:])
					
				jnet = ''
				for line in JNET_lines:
					jnet += line.replace(line, line[space_count:])
					
				psipred = ''
				for line in PSIPRED_lines:
					psipred += line.replace(line, line[space_count:])
					
				#sympred = ''
				#for line in SYMPRED_lines:
					#sympred += line.replace(line, line[space_count:])
				
				def replace_space_with_c(letter):
					if letter == ' ':
						out = letter.replace(letter, "C")
					else:
						out = letter
						
					return(out)
				
				aa_list = []
				for letter in aa:
					aa_list.append(letter)
				aa_dict = {'AA':aa_list}
				
				phd_list = []
				for letter in phd:
					phd_list.append(replace_space_with_c(letter))
				phd_dict = {'PHD':phd_list}
				
				prof_list = []
				for letter in prof:
					prof_list.append(replace_space_with_c(letter))
				prof_dict = {'PROF':prof_list}
				
				sspro_list = []
				for letter in sspro:
					sspro_list.append(replace_space_with_c(letter))
				sspro_dict = {'SSPRO':sspro_list}
				
				jnet_list = []
				for letter in jnet:
					jnet_list.append(replace_space_with_c(letter))
				jnet_dict = {'JNET':jnet_list}
				
				psipred_list = []
				for letter in psipred:
					psipred_list.append(replace_space_with_c(letter))
				psipred_dict = {'PSIPRED':psipred_list}
				
				# ~ sympred_list = []
				# ~ for letter in sympred:
					# ~ sympred_list.append(replace_space_with_c(letter))
				# ~ sympred_dict = {'SYMPRED':sympred_list}
				
					
				output_data = dict()
				for dictionary in [phd_dict, prof_dict, sspro_dict, jnet_dict, psipred_dict]:
					output_data.update(dictionary)
					
				print(output_data)
				
				sympred_end_time = time.time()
				sympred_execution_time = sympred_end_time - sympred_start_time
				print(f'Sympred took {sympred_execution_time} sec')
				
				
				return output_data
				break
			else:
				print(f'Sympred: Request failed with status code {results_page.status_code}')
				print(f'Sympred: Aborting operation.')
				break
		
	else:
		print('Sympred: Unable to get job ID')
	
def runYaspin():
	yaspin_start_time = time.time()
	
	print('Running Yaspin')
	
	boundary = str(uuid.uuid4())

	# Define the form data to be submitted
	data = {
		'seq': sequence,
		'seq_file': ('', ''),
		'pssm_file': ('', ''),
		'smethod': 'nr',
		'nnmethod': 'dssp',
		'email': email,
		'mbjob[description]': seq_name,
		'yaspin_align': 'YASPIN prediction',
	}
	
	# Set the headers for the request
	headers = {
		'Content-Type': 'multipart/form-data; boundary=' + boundary
	}
	
	# Build the request payload
	payload = ''
	for key, value in data.items():
		payload += '--' + boundary + '\r\n'
		if isinstance(value, tuple):
			payload += 'Content-Disposition: form-data; name="' + key + '"; filename="' + value[0] + '"\r\n'
			payload += 'Content-Type: ' + value[1] + '\r\n\r\n'
		else:
			payload += 'Content-Disposition: form-data; name="' + key + '"\r\n\r\n'
		payload += str(value) + '\r\n'

	payload += '--' + boundary + '--\r\n'

	# Define the URL of the Yaspin web form
	url = 'https://www.ibi.vu.nl/programs/yaspinwww/'

	# Submit the form data using a POST request
	response = requests.post(url, headers=headers, data=payload, allow_redirects=True)

	# Check if the POST request was successful
	if response.status_code == 202:
		# Extract the job ID from the redirect URL
		job_id = response.url.split('/')[-2]
		print(f'Yaspin Job ID: {job_id}')
		
		# Define the URL of the job results page
		results_url = f'http://zeus.few.vu.nl/jobs/{job_id}/results.out'

		# Wait for the job to complete and the results to become available
		print('Yaspin: Waiting for results...')
		
		while True:
			
			# check elapsed time
			if (time.time() - yaspin_start_time) > yaspin_timeout:
				print(f'Yaspin timed out. Exceeded {yaspin_timeout} seconds.')
				break
			
			# Send a GET request to the results page
			results_page = requests.get(results_url)
				
			if results_page.status_code == 404:
				sleep(5)
			elif results_page.status_code == 200:
				text_lines = results_page.text.split('\n')
				yaspin_string = ''
				for line in text_lines:
					if '*' in line:
						continue
					elif 'Pred: ' in line:
						yaspin_string += line.split('Pred: ')[1].split('\n')[0]
					else:
						continue
				
				yaspin_out = []
				for letter in yaspin_string:
					if letter == '-':
						yaspin_out.append(letter.replace('-', 'C'))
					else:
						yaspin_out.append(letter)
					
				yaspin = {'YASPIN': yaspin_out}
				print(yaspin)
				
				yaspin_end_time = time.time()
				yaspin_execution_time = yaspin_end_time - yaspin_start_time
				print(f'Yaspin completed in {yaspin_execution_time} sec')
				return yaspin
				break
				
			else:
				print(f'Yaspin: Request failed with status code {results_page.status_code}')
				print(f'Yaspin: Aborting operation.')
				return None
				break
		
	else:
		print('Yaspin: Unable to get job ID')
		return None
		
def process_data(data_list):
	pd_start = time.time()
	print(f'Processing data...')
	
	# Make a new dictionary that collects all of the other dictionary information
	all_algos_dict = dict()
	for each in data_list:
		all_algos_dict.update(each)
			
	for key, value in zip(all_algos_dict.keys(), all_algos_dict.values()):
		print(f'{key}: {value}\n')
		
	def calc_vkabat(ss_assignment_list):
		# k = number of secondary structure classes predicted for a residue (1, 2, or 3)
		# N = total number of predictions (15)
		# n1 = how often most predicted ss is observed for a residue
		
		N = len(ss_assignment_list)
		
		unique_ss_assignment_list = np.unique(np.array(ss_assignment_list))
		
		k = len(unique_ss_assignment_list)
		
		counts_list = []
		E_COUNT = 0
		H_COUNT = 0
		C_COUNT = 0
		T_COUNT = 0
		for ss in ss_assignment_list:
				if ss == 'E':
					E_COUNT += 1
				elif ss == 'H':
					H_COUNT += 1
				elif ss == 'C':
					C_COUNT += 1
				elif ss == 'T':
					T_COUNT += 1
				else:
					print('Encountered secondary structure assignment other than E, H, C, or T!')
					
		counts_list.append(E_COUNT)
		counts_list.append(H_COUNT)
		counts_list.append(C_COUNT)
		counts_list.append(T_COUNT)
		
		total_counts = 0
		for count in counts_list:
			total_counts += count

		E_perc = round(Decimal(str(E_COUNT / total_counts * 100)), 2)
		H_perc = round(Decimal(str(H_COUNT / total_counts * 100)), 2)
		C_perc = round(Decimal(str(C_COUNT / total_counts * 100)), 2)
		T_perc = round(Decimal(str(T_COUNT / total_counts * 100)), 2)
		
		n1 = max(counts_list)
		
		vkabat = k * N / n1
		
		calc_dict = {
						'E_COUNT': E_COUNT, 
						'H_COUNT': H_COUNT, 
						'C_COUNT': C_COUNT, 
						'T_COUNT': T_COUNT, 
						'total_counts': total_counts,
						'E_perc': E_perc,
						'H_perc': H_perc,
						'C_perc': C_perc,
						'T_perc': T_perc,
						'k': k,
						'N': N,
						'n1': n1,
						'vkabat':vkabat
						} 
		
		return calc_dict
		
	i = 0
	list_2d = []
	while i < len(sequence):
		list_1d = []
		for x in all_algos_dict.values():
			list_1d.append(x[i])
		list_2d.append(list_1d)
		i += 1
		
	E_COUNT_list = []
	H_COUNT_list = []
	C_COUNT_list = []
	T_COUNT_list = []
	total_counts_list = []
	E_perc_list = []
	H_perc_list = []
	C_perc_list = []
	T_perc_list = []
	k_list = []
	N_list = []
	n1_list = []
	vkabat_list = []
	for each in list_2d:
		E_COUNT_list.append(calc_vkabat(each)['E_COUNT'])
		H_COUNT_list.append(calc_vkabat(each)['H_COUNT'])
		C_COUNT_list.append(calc_vkabat(each)['C_COUNT'])
		T_COUNT_list.append(calc_vkabat(each)['T_COUNT'])
		total_counts_list.append(calc_vkabat(each)['total_counts'])
		E_perc_list.append(calc_vkabat(each)['E_perc'])
		H_perc_list.append(calc_vkabat(each)['H_perc'])
		C_perc_list.append(calc_vkabat(each)['C_perc'])
		T_perc_list.append(calc_vkabat(each)['T_perc'])
		k_list.append(calc_vkabat(each)['k'])
		N_list.append(calc_vkabat(each)['N'])
		n1_list.append(calc_vkabat(each)['n1'])
		vkabat_list.append(calc_vkabat(each)['vkabat'])
	
	vkabat_out_dict = {
						'E_COUNT': E_COUNT_list, 
						'H_COUNT': H_COUNT_list, 
						'C_COUNT': C_COUNT_list, 
						'T_COUNT': T_COUNT_list, 
						'total_counts': total_counts_list,
						'E_perc': E_perc_list,
						'H_perc': H_perc_list,
						'C_perc': C_perc_list,
						'T_perc': T_perc_list,
						'k': k_list,
						'N': N_list,
						'n1': n1_list,
						'vkabat': vkabat_list
						}
	
	print(f'vkabat: {vkabat_out_dict["vkabat"]}')
	
	all_algos_dict.update(vkabat_out_dict)
	df = pd.DataFrame(data=all_algos_dict)
	
	print('')
	print(df)
	print('')
	
	# Write DataFrame csv file
	vkabat_dataframe_file_name = str(seq_name) + '_vkabat_dataframe.csv'
	vkabat_dataframe_file_name_path = os.path.join(output_directory, vkabat_dataframe_file_name)
	print(f'Writing vkabat (data frame) csv file {vkabat_dataframe_file_name_path}')
	df.to_csv(vkabat_dataframe_file_name_path)
	
	# Write vkabat csv file
	vkabat_only_file_name = str(seq_name) + '_vkabat.csv'
	vkabat_only_file_name_path = os.path.join(output_directory, vkabat_only_file_name)
	print(f'Writing vkabat (only) csv file {vkabat_only_file_name_path}')
	df2 = pd.DataFrame(data={'vkabat': vkabat_out_dict["vkabat"]})
	df2.to_csv(vkabat_only_file_name_path, index=False)
	
	
	pd_end = time.time()
	print(f'Processing Data completed in {pd_end - pd_start} seconds')
	
	return vkabat_out_dict
	
def print_banner():
	print('\n                                 ')
	print('   R   U   N   N   I   N   G   :	  ')
	print('   _   _   _   _   _   _   _   _   ')
	print('  / \ / \ / \ / \ / \ / \ / \ / \  ')
	print(' ( P | y | V | k | a | b | a | t ) ')
	print('  \_/ \_/ \_/ \_/ \_/ \_/ \_/ \_/  ')
	print('   Version: 0.0.1                  ')
	print('                                   ')
	print('   Code available at: "https://github.com/1000000000000000000000000000000/PyVkabat"\n\n')

def main():
	
	# start clock
	prog_start_time = time.time()
	
	# print banner
	print_banner()
	
	# argparse code
	parse_arguments()

	# supress irrelevant warnings in bs4
	warnings.filterwarnings("ignore", category=UserWarning, module='bs4')
	
	with concurrent.futures.ThreadPoolExecutor() as executor:
		prabi_future = executor.submit(runPrabi)
		jpred_future = executor.submit(run_alt_JPred)
		yaspin_future = executor.submit(runYaspin)
		sympred_future = executor.submit(runSympred5)
		
		prabi_output = prabi_future.result()
		jpred_output = jpred_future.result()
		yaspin_output = yaspin_future.result()
		sympred_output = sympred_future.result()
    
	print(f'Total elapsed time to retrieve data: {time.time() - prog_start_time} seconds\n')
	
	algo_list = [prabi_output, yaspin_output, jpred_output, sympred_output]
	successful_algo_list = []
	failed_algo_list = []
	for algo in algo_list:
		if algo == None:
			continue
		else:
			successful_algo_list.append(algo)
	
	process_data(successful_algo_list)
	
	print(f'Total running time: {time.time() - prog_start_time} seconds.')
	print('Done.')
	
if __name__ == '__main__':
	main()
