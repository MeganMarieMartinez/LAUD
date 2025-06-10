import glob
import os
import csv 


'''
Parses through a directory of VIBRANT outputs to extract number of phages found, number of high and medium quality phage drafts 
and writes them to files
'''

# change source to be the directory with VIBRANT outputs
source = '/path/to/directory/with/vibrant/outputs'


# function to grab the accession number
def accession_finder(file):
    x = file.split(f'{source}/')
    y = x[1]
    x = y.split('_ASM')
    accession = x[0]
    return accession

# function to extract the filename
def filename_finder(file):
    x = file.split('/')
    filename = x[-1]
    return filename

# function to write high quality drafts to a fasta file
def high_quality_writer(high_records): # this takes a list
    with open('high_quality.txt','a') as f:
        f.write(str(high_records)+ '\n')
        f.close()

# function to write the medium quality drafts to a fasta file
def med_quality_writer(med_records): # this takes a list
    with open('med_quality.txt','a') as f:
        #csvwriter = csv.writer(csvfile)
        f.write(str(med_records) + '\n')
        f.close()

no_phages = [] # list to just hold accessions that did not find a phage
rows = [] # this will hold ALL THE ROWS
med_records = []
high_records = []
re_run = []
genus_species = ''
UMB = ''

# This for loop just determines if there is a putative phage or not 
phage_flag = False
for directory in glob.glob(f'{source}/*'):
    info_holder = []
    directory_size = os.path.getsize(directory)
    accession = accession_finder(directory)
    if directory_size == 16384:
        for log_file in glob.glob(f'{source}/{directory}/VIBRANT_{accession}_genomic/VIBRANT_log_run*'):
            with open(log_file,"r") as f:
                data = f.readlines()
                # finds the line with the # of phages, checks how many there are 
                for line in data:
                    if "putative phages were identified." in line:
                        if line.startswith('0'):
                            no_phages.append(accession)
                        else:
                            phage_flag = True
                            phages_found = int(line[0])
                    if 'No phages' in line:
                        no_phages.append(accession)
    elif directory_size == 4096:
        re_run.append(directory) 
    else:
        for log_file in glob.glob(f'{directory}/VIBRANT_{accession}_*_genomic/VIBRANT_log_run*'): # accesses just the log files
            with open(log_file,"r") as f:
                data = f.readlines()
                # finds the line with the # of phages, checks how many there are 
                for line in data:
                    if "putative phages were identified." in line:
                        if line.startswith('0'):
                            no_phages.append(accession)
                        else:
                            phage_flag = True
                            if line[1].isdigit():
                                phages_found = int(line[0]+line[1])
                            else:
                                phages_found = int(line[0])
                    if 'No phages' in line:
                        no_phages.append(accession)
    if phage_flag is True:
        medium = 0 # medium quality draft
        high = 0 # high quality draft
        medium_scaffolds = []
        high_scaffolds = []
        quality_file = str(glob.glob(f'{directory}/VIBRANT_{accession}_*_genomic/VIBRANT_results_*_genomic/VIBRANT_genome_quality*'))
        quality_file = quality_file[2:-2]
        print(quality_file)
        if quality_file != '':
            with open(quality_file) as fd:
                rd = csv.reader(fd, delimiter="\t", quotechar='"')
                next(rd)
                for row in rd:
                    scaffold = row[0]
                    s = scaffold.split(',')
                    s = s[0]
                    s = s.split(' ')
                    genus_species = s[1] + ' ' + s[2]
                    for word in s:
                        if word.startswith('UMB'):
                            UMB = word
                    if row[-1] == 'medium quality draft':
                        medium += 1
                        medium_scaffolds.append(scaffold)
                    if row[-1] == 'high quality draft':
                        high += 1
                        high_scaffolds.append(scaffold)
                    if UMB == '':
                        UMB = s[-2]
                high_quality_writer(high_scaffolds)
                med_quality_writer(medium_scaffolds)
    info_holder = info_holder + [accession,genus_species,UMB,phages_found,medium,high]


    rows.append(info_holder)


                            
with open('sample_master.txt', 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
    fields = ['Accession Number','Genus and species','Intraspecies ID','# of Phages Found','Medium Quality Drafts','High Quality Drafts']
    # writing the fields
    csvwriter.writerow(fields)
    # writing the data rows
    csvwriter.writerows(rows)

