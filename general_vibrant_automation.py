import os

'''
Automation script to run VIBRANT on a directory of FASTA files
'''

directory = '/path/to/directory/with/FASTA/files'
destination = '/path/to/output/directory'

for filename in os.listdir(directory):
    if filename.startswith('GCA'): # makes sure its a fasta file it is accessing just in case
        accession = filename.split('.fna') # based on how my files were named, isolated the accession number
        folder_name = accession[0]
        upload_command =  f'VIBRANT_run.py -i {directory}/{filename} -folder {destination}/{folder_name}_VIBRANT_OUTPUT'
        os.system(upload_command)