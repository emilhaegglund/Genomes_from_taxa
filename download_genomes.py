#!/usr/bin/env python

import requests
import hashlib
import argparse
import multiprocessing as mp
from NCBI_taxonomy_tree import ncbiTaxonomyTree
import os
import tarfile

def parse_arguments():
    """Parse command line arguments.
    Return:
        args: Object with command line arguments.
    """
    parser = argparse.ArgumentParser(description='Download genomes from NCBI')

    parser.add_argument('--taxid', '-t', required=True,
                        help='Taxonomical ID')

    parser.add_argument('--update', '-u',
                        help='Update')

    parser.add_argument('--output', '-o', required=True,
                        help='Output directory to store downloaded genomes')

    parser.add_argument('--filter', '-f',
                        help='''Filter on assembly status. Allowed values are
                        complete, scaffold, contig. Example use -f complete,contig
                        will download genomes with Complete genome or contig as assembly status''')
    parser.add_argument('--workers', '-w',
                        help='Number of parallel downloads')

    args = parser.parse_args()


    return args

def download_assembly_summary(basedir, store=True):
    """
    Download NCBI assembly summary

    Arguments:
    store - store the assembly summary on file for faster access

    Returns
    assembly_summary - list containing each entry in the assembly summary file
    """
    print('Download assembly summary')
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt'
    assembly_summary = []
    assembly_summary_file = os.path.join(basedir, 'assembly_summary.txt')
    # open assembly summary file
    if store:
        f = open(assembly_summary_file, 'w')

    ret = requests.get(url, stream=True)

    # handle the request
    for line in ret.iter_lines():
        if store:
            f.write(line.decode('utf-8', 'ignore') + '\n')
        if line.decode('utf-8')[0] != '#':  # skip comment lines
            line = line.decode('utf-8', 'ignore').strip('\n')
            line = line.split('\t')
            assembly_summary.append(line)

    # close assembly summary file
    if store:
        f.close()

    return assembly_summary

def read_assembly_summary(basedir):
    """
    Read local copy of NCBI assembly summary file.

    Returns
    assembly_summary - list containing each entry in the assembly summary file
    """
    assembly_summary = []
    assembly_summary_file = os.path.join(basedir, 'assembly_summary.txt')
    with open(assembly_summary_file, 'r') as f:
        for line in f:
            if line[0] != '#':  # skip comment lines
                line = line.strip('\n')
                line = line.split('\t')
                assembly_summary.append(line)

    print(len(assembly_summary))
    return assembly_summary

def download_taxonomy(basedir):
    """
    Download NCBI taxonomy data
    """
    print('Download taxonomy')

    taxonomy_file = os.path.join(basedir, 'taxdump.tar.gz')

    url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
    ret = requests.get(url, stream=True)
    with open(taxonomy_file, 'wb') as f:
        for chunk in ret.iter_content(4096):
            f.write(chunk)

    t = tarfile.open(taxonomy_file, 'r:gz')

    t.extract('nodes.dmp', path=basedir)
    t.extract('names.dmp', path=basedir)

    tree = create_ncbi_tree(basedir)

    return tree

def check_taxonomy_files(basedir):
    """
    Check if nodes.dmp and names.dmp exists

    Returns:
    True if both file exists
    False otherwise
    """

    names = False
    nodes = False

    names_file = os.path.join(basedir, 'names.dmp')
    nodes_file = os.path.join(basedir, 'nodes.dmp')

    if os.path.isfile(names_file):
        names = True

    if os.path.isfile(nodes_file):
        nodes = True

    if nodes == True and names == True:
        return True
    else:
        return False


def create_ncbi_tree(basedir, store=True):
    """
    Create a tree object from taxonomy files
    TODO: Try to store the tree to file for faster accessing.

    Returns:
    tree - ncbi taxonomy tree
    """
    print('Build tree')
    names_file = os.path.join(basedir, 'names.dmp')
    nodes_file = os.path.join(basedir, 'nodes.dmp')

    tree = ncbiTaxonomyTree.NcbiTaxonomyTree(nodes_filename=nodes_file, names_filename=names_file)

    return tree

def download(url):
    """
    Download files using requests.get.

    Returns:
    ret -
    """
    ret = requests.get(url, stream=True)

    return ret

def download_genome(assembly_line):
    """
    Download NCBI entry for taxid

    Arguments:
    taxid - Taxnomic identifier
    Returns:
    """
    print('Download: ' + assembly_line[7] + ' ' + assembly_line[0])
    # create output dir
    base_dir = assembly_line[-1]
    basename = assembly_line[7] + '_' + assembly_line[0]
    basename = basename.replace(' ', '_')
    basename = basename.replace('/', '_')
    #else:
    #   basename = assembly_line[0]

    output_dir = os.path.join(base_dir, basename)

    os.mkdir(output_dir)

    # download hash file
    ftp_url = assembly_line[19]
    https_url = ftp_to_https(ftp_url)
    hash_url = https_url + '/md5checksums.txt'
    hash_dict = download_hash_file(hash_url, output_dir)

    # download all other files in directory
    for file in hash_dict.keys():

        # check for and create subfolders, most of it is renaming of folders.
        if len(file.split('/')) > 2:
            outfile = file[2:].split('/')
            subfolder = outfile[0]
            subfolder = subfolder.split('_')[3:]
            subfolder = assembly_line[7].replace(' ', '_') + '_' + assembly_line[0] + '_' + '_'.join(subfolder)
            end = '/'.join(outfile[1:-1])
            path = os.path.join(output_dir, subfolder, end)
            if not os.path.isdir(path):
                os.makedirs(path)

        # create url
        file_url = https_url + file[1:]

        # download file
        #print('Download: ' + file_url)
        ret = download(file_url)

        # handle download
        if file[2:] == 'annotation_hashes.txt':
            out_file = file[2:]
        else:
            out_file = file[2:].split('_')
            out_file = basename + '_' + '_'.join(out_file[3:])
        output_path = os.path.join(output_dir, out_file)
        with open(output_path, 'wb') as f:
            for chunk in ret.iter_content():
                f.write(chunk)

        # check the integrity of the file
        local_checksum = hash_file(output_path)
        remote_checksum = hash_dict[file]
        if not local_checksum == remote_checksum:
            print('File checksum disagree, removing file')

def download_hash_file(url, output_dir):
    """
    Download the md5checksums.txt to check the integrity of downloaded files.
    Also use this file to find the structure of the directory.

    Arguments:
    url - url to md5checksums.txt
    output_dir - path to outpur directory for specie

    Returns:
    hash_dict - key: path to file in directory, value: md5sum for file
    """
    hash_dict = {}
    hash_path = os.path.join(output_dir, 'md5checksums.txt')

    # download md5checksum.txt
    ret = download(url)

    # handle the return
    with open(hash_path, 'w') as f:
        for line in ret.iter_lines():
            f.write(line.decode('utf-8') + '\n')
            line = line.decode('utf-8').strip('\n')
            line = line.split(' ')
            hash_dict[line[2]] = line[0]

    return hash_dict


def ftp_to_https(ftp_url):
    """
    Convert ftp urls to https urls

    Argument:
    ftp_url - ftp://link/to/something

    Returns:
    url - https://link/to/something
    """
    url = ftp_url.replace('ftp://', 'https://')

    return url

def species_to_download(taxid, tree, assembly_summary):
    """
    Find all descendants to the taxid that are in the assembly_summary.txt and return
    a list with the corresponding entries

    Arguments:
    taxid - Taxonomical id to
    tree - NCBI taxonomy tree
    assembly_summary - list of entries from the assembly_summary.txt file

    Returns:
    species_list - entries from assembly_summary.txt that are descendats from taxid
    """
    species_list = []
    descendants = tree.getLeaves(taxid)

    # find descendants that are also in the assembly summary
    for specie in descendants:
        for assembly_line in assembly_summary:
            if int(specie) == int(assembly_line[5]):
                species_list.append(assembly_line)

    return species_list

def hash_file(path):
    """
    Caclulate the checksum of a file

    Arguments:
    path -- path to the file

    Returns:
    local_checksum -- the md5 checksum of the local file
    """
    local_checksum = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            local_checksum.update(chunk)

    return local_checksum.hexdigest()

def main():

    args = parse_arguments()
    update = False
    taxid = int(args.taxid)
    base_dir = args.output
    filter_values = args.filter.split(',')
    assembly_filter = []
    filter_value_map = {'complete':'Complete Genome', 'scaffold':'Scaffold', 'contig':'Contig'}
    for value in filter_values:
        if value in filter_value_map.keys():
            assembly_filter.append(filter_value_map[value])
        else:
            print(value + ' is not a valid filter and will not be used')

    # Set number of workers to use
    if not args.workers:
        j = 1
    else:
        j = int(args.workers)

    if not os.path.exists(base_dir):
        os.mkdir(base_dir)
    assembly_summary_file = os.path.join(base_dir, 'assembly_summary.txt')
    if os.path.isfile(assembly_summary_file) and not update:
        assembly_summary = read_assembly_summary(base_dir)
    else:
        assembly_summary = download_assembly_summary(base_dir)

    if check_taxonomy_files(base_dir) and not update:
        tree = create_ncbi_tree(base_dir)
    else:
        tree = download_taxonomy(base_dir)

    species_list = species_to_download(taxid, tree, assembly_summary)
    for assembly_line in species_list:
        assembly_line.append(base_dir)

    print('Total number of genomes in ' + str(taxid) + ':' + str(len(species_list)))

    filtered_species_list = []
    for assembly_line in species_list:
        if assembly_line[11] in assembly_filter:
            filtered_species_list.append(assembly_line)
    print('Will download ' + str(len(filtered_species_list)) + ' genomes')

    pool = mp.Pool(processes=j)
    pool.map(download_genome, filtered_species_list)

    #for assembly_line in species_list:
    #    print(assembly_line[7])
    #    download_genome(3, assembly_line, base_dir)


    directory_list_path = base_dir + '/list_of_dirs.txt'
    directory_list = open(directory_list_path, 'w')
    for assembly_line in species_list:
        directory = os.path.join(base_dir, assembly_line[0])
        directory_list.write(directory + '\t' + assembly_line[7])

main()
