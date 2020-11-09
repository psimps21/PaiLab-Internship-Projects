import os
import pwd
from argparse import ArgumentParser
from pathlib import Path
import datetime
import csv
import multiprocessing
from multiprocessing import Pool


# Traverses a directory to find the total size, size of first level directories, size of second level
# directories, and a list of unacceptable files
# All data collected is written into a csv
def multiprocessing_monitor(args):
    start = datetime.datetime.now()
    root = Path(args.root)
    tree_dict = {
        'layer1': [str(root.joinpath(i)) for i in os.listdir(str(root))]
    }
    user_unaccpt_files = {}
    name_dir_size = sum([os.path.getsize(i) for i in tree_dict['layer1']])

    pool = Pool(processes=multiprocessing.cpu_count())
    futures = pool.map(worker_process, [(i, args.filetypes) for i in tree_dict['layer1']])
    produce_dictionary_csv(futures, args.outdir, user_unaccpt_files, start, name_dir_size)


# Produces a cvs file from data in a dictionary
def produce_dictionary_csv(futures_list, outdir, user_unaccpt_files, start, name_dir_size):
    logtime = datetime.datetime.now()
    total_size = 0
    path = os.path.join(outdir, 'Multiprocessing_User_Data_Size.csv')
    with open(path, 'w') as total_data:
        data_doc = csv.writer(total_data)
        data_doc.writerow(['owner', 'path', 'totalsizedata(Gb)'])
        for dict in futures_list:
            if len(dict['subdir_dict_list']) != 0:
                for i in dict['subdir_dict_list']:
                    data_doc.writerow([i['owner'], i['path'], round(bytes_to_gigabytes(i['size']), 3)])

            total_size += dict['size']

        total_size += name_dir_size
        data_doc.writerow(['RUNTIME', logtime-start])
        data_doc.writerow(['LOG TIME', logtime])
        data_doc.writerow(['TOTAL SIZE', round(bytes_to_gigabytes(total_size), 3)])

    path2 = os.path.join(outdir, 'Unaccepted_Files.csv')
    with open(path2, 'w') as unaccpt_files:
        file_doc = csv.writer(unaccpt_files)
        file_doc.writerow(['owner', 'filesize(Gb)', 'filetype', 'filepath'])
        for dict in futures_list:
            for i in dict['unaccpted_files']:
                file_doc.writerow([i[3], round(bytes_to_gigabytes(i[2]), 3), i[1], i[0]])



# Executes os.walk on a directory in order to get the size of all internal files and subdirectories.
# Also keeps track of all unacceptable file types that a use has
# Returns a dictionary containing size and a list of unaccepted file paths
def worker_process(root):
    unaccpt_file_list = [
        '.sam',
        '.fastq',
        '.fq']
    if root[1]:
        unaccpt_file_list.extend(root[1])
    data_dict = {
        'size': 0,
        'name': os.path.basename(root[0]),
        'unaccpted_files': [],
        'subdir_dict_list': [],
        'owner': pwd.getpwuid(os.stat(root[0]).st_uid)[0]
    }

    if os.path.isdir(root[0]):
        for dir_path, subdirs, files in os.walk(root[0]):
            for subdir in subdirs:
                if dir_path == root[0]:
                    data_dict['subdir_dict_list'].append({
                        'path': os.path.join(root[0], subdir),
                        'name': subdir,
                        'size': 0,
                        'owner': pwd.getpwuid(os.stat(os.path.join(root[0], subdir)).st_uid)[0]
                    })

                for i in data_dict['subdir_dict_list']:
                    if i['name'] in dir_path.split('/'):
                        i['size'] += os.path.getsize(os.path.join(dir_path, subdir))

                data_dict['size'] += os.path.getsize(os.path.join(dir_path, subdir))

            for file in files:
                file_path = Path(dir_path).joinpath(file)
                suffixes = file_path.suffixes
                data_dict['size'] += os.path.getsize(str(file_path))

                for i in data_dict['subdir_dict_list']:
                    if i['name'] in dir_path.split('/') and file[0] != '.':
                        i['size'] += os.path.getsize(str(file_path))

                if len(suffixes) != 0:
                    if suffixes[-1] in unaccpt_file_list:
                        file_tuple = (str(Path(dir_path).joinpath(file)), suffixes[-1], os.path.getsize(str(file_path)),
                                      pwd.getpwuid(os.stat(str(file_path)).st_uid)[0])
                        data_dict['unaccpted_files'].append(file_tuple)

    elif os.path.isfile(root[0]):
        suffixes = Path(root[0]).suffixes
        if len(suffixes) != 0:
            if suffixes[-1] in unaccpt_file_list:
                file_tuple = (str(root[0]), suffixes[-1], os.path.getsize(str(root[0])),
                              pwd.getpwuid(os.stat(str(root[0])).st_uid)[0])
                data_dict['unaccpted_files'].append(file_tuple)

    return data_dict


def get_args():
    parser = ArgumentParser(
        description='Determine the size of all subdirectories and identify the presence of specified file types within'
                    ' a given directory.')
    parser.add_argument('root',
                        help='Path to the desired directory.')
    parser.add_argument('-o', '--outdir',
                        help='The directory in which the output files will be created.',
                        default='.')
    parser.add_argument('-f', '--filetypes',
                        help='List of all file types to search for.',
                        default=None,
                        nargs='*')

    return parser.parse_args()


def bytes_to_kilobytes(bytes):
    return bytes/1024

def bytes_to_megabytes(bytes):
    return bytes_to_kilobytes(bytes)/1024

def bytes_to_gigabytes(bytes):
    return bytes_to_megabytes(bytes)/1024


if __name__ == '__main__':
    args = get_args()
    multiprocessing_monitor(args)
