## This script contains multiple python functions for general use

import os
import logging
logger = logging.getLogger(__name__)

# Checks if the directory exists, if not creates a new directory
def mkdirs(output_dir):
    try:
        os.listdir(output_dir)
    except FileNotFoundError:
        logger.info('Making directory for ' + output_dir)
        os.mkdir(output_dir)

# Count the number of lines in the provided file
def count_lines(input_file):
    logger.info('Counting the number of lines in {}'.format(input_file))
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
        no_of_lines = len(lines)
    return no_of_lines

# Concatenates files within a specified folder using a prefix of the filename (all files should be compressed)
def concat_prefix_gz(output_zipped, input_dir, prefix):
    import os
    import gzip
    logger.info('Concatenating all compressed files with prefix {} into {}'.format(prefix, output_zipped))
    with gzip.open(output_zipped, 'wb') as outfile:
        for file in sorted(os.listdir(input_dir)):
            if file.startswith(prefix):
                with gzip.open(input_dir+file, 'rt') as infile:
                    for line in infile:
                        binline = str.encode(line)
                        outfile.write(binline)

# Concatenate two compressed files
def concat_gz(output_zipped, input1_zipped, input2_zipped):
    import os
    import gzip
    logger.info('Concatenating compressed {} and {}'.format(input1_zipped, input2_zipped))
    with gzip.open(output_zipped, 'wb') as outfile:
        with gzip.open(input1_zipped, 'rt') as infile_1:
            for line in infile_1:
                binline = str.encode(line)
                outfile.write(binline)
        with gzip.open(input2_zipped, 'rt') as infile_2:
            for line in infile_2:
                binline = str.encode(line)
                outfile.write(binline)

# Concatenates files within a specified folder using a prefix of the filename (Files should not be compressed)
def concat_prefix_file(output_file, input_dir, prefix):
    import os
    logger.info('Concatenating all files with prefix {} into {}'.format(prefix, output_file))
    with open(output_file, 'w') as outfile:
        for file in sorted(os.listdir(input_dir)):
            if file.startswith(prefix):
                with open(input_dir+file, 'r') as infile:
                    for line in infile:
                        outfile.write(line)

# Concatenate two files (Files should not be compressed)
def concat_files(output_file, input1_file, input2_file):
    logger.info('Concatenating {} and {}'.format(input1_file, input2_file))
    with open(output_file, 'w') as outfile:
        with open(input1_file, 'r') as infile_1:
            for line in infile_1:
                outfile.write(line)
        with open(input2_file, 'r') as infile_2:
            for line in infile_2:
                outfile.write(line)

# Move files from a specified directory using a prefix of the filename
def move_with_prefix(dir_with_files, prefix_of_files, new_path):
    logger.info('Moving all files with prefix {} to {}'.format(prefix_of_files, new_path))
    for file in os.listdir(dir_with_files):
        if file.startswith(prefix_of_files):
            new_file_path = os.path.join(new_path, file)
            old_file_path = os.path.join(dir_with_files, file)
            os.rename(old_file_path, new_file_path)

# Find a pattern within a file and print those lines + optional lines before and after the pattern to a new compressed file (files should be compressed)
def find_pattern_gz(output_file, input_file, pattern, prev_lines=0, after_lines=0):
    import gzip
    logger.info('Finding pattern {} in compressed file {} and printing {} lines before the pattern and {} lines after the pattern'.format(pattern, input_file, prev_lines, after_lines))
    with gzip.open(output_file, 'wb') as outfile:
        with gzip.open(input_file, 'rt') as infile:
            lines = infile.readlines()
            for index, line in enumerate(lines):
                if pattern in line:
                    pattern_list = lines[min(index, index-prev_lines): max(index, index+after_lines)]
                    for line in pattern_list:
                        binline = str.encode(line)
                        outfile.write(binline)

# Find a pattern within a file and print those lines + optional lines before and after the pattern to a new file (files should not be compressed)
def find_pattern(output_file, input_file, pattern, prev_lines=0, after_lines=0):
    logger.info('Finding pattern {} in file {} and printing {} lines before the pattern and {} lines after the pattern'.format(pattern, input_file, prev_lines, after_lines))
    with open(output_file, 'w') as outfile:
        with open(input_file, 'r') as infile:
            lines = infile.readlines()
            for index, line in enumerate(lines):
                if pattern in line:
                    pattern_list = lines[min(index, index-prev_lines): max(index, index+after_lines)]
                    for line in pattern_list:
                        outfile.write(line)
