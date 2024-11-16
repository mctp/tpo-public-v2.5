#!/usr/bin/env python3

# split intervals generated from ScatterIntervalsByNs to segments in separated files
# usage: prepare_intervals.py [-h] -i <intervals> -o <output_directory>

import argparse
ap = argparse.ArgumentParser()
ap.add_argument('-i', required=True, help='intervals generated from ScatterIntervalsByNs')
ap.add_argument('-o', required=True, help='output directory storing the splits')
args = ap.parse_args()

with open(args.i) as reader:
    total = 0
    file_count = 0

    header = ''
    regions = ''

    for line in reader:
        if line.startswith('@'):
            header += line
            continue
        line = line.rstrip()
        chrname = line.split('\t')[0]
        start = int(line.split('\t')[1])
        end = int(line.split('\t')[2])

        length = end - start +1
        total = total + length
        regions = regions + f'{chrname}\t{start}\t{end}\t+\t.\n'

        if total >= 1000000:
            prefix = str(file_count).zfill(3)
            with open(args.o + f'/{prefix}.interval_list', 'w') as writer:
                writer.write(header + regions)

            file_count += 1
            total = 0
            regions = ''

    # the last output
    prefix = str(file_count).zfill(3)
    with open(args.o + f'/{prefix}.interval_list', 'w') as writer:
        writer.write(header + regions)

