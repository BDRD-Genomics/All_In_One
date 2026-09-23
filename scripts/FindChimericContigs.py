#!/usr/bin/env python
# Logan Voegtly
# Edited by Anthony Cicalo
# Naval Medical Research Center, Biological Defense Research Directorate

import os
import os.path as path
import sys
import re
import json
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import argparse
import ete3
import csv 

#MAX_TARGET_SEQS = 10
taxid_index = 5
start_index = 2
end_index = 3
ETE3_TAXA_DB = os.path.join(os.sep, '/export', 'nextflow', 'bdrd', 'all_in_one_pipeline', 'bin', 'etetoolkit', 'taxa.sqlite')
# ETE3_TAXA_DB = os.path.join(os.sep, 'Users', 'lvoegtly', '.etetoolkit', 'taxa.sqlite')
ncbi = ete3.NCBITaxa(dbfile=ETE3_TAXA_DB)
desired_ranks = ['superkingdom','genus','order','family','phylum', 'class','species']
#desired_rank = 'family'
#ncbi = ete3.NCBITaxa()
NO_ORDER_ID = -1
MIN_GAP_LEN = 200
FIND_CONTIG_AND_OFFSET = re.compile(r'(.*)_start_pos:(\d+)')
UNCULTURED_BACTERIUM_TAXID = '77133'
SYNTHETIC_CONSTRUCT_TAXID = '32630'
PHAGE_ORDER_TAXID = '28883'
LEN_BETWEEN_TAXA = 300
### New Function ### Testing started 05132025 (acicalo)
#def has_data(filepath):
#    with open(filepath, newline='') as tsvfile:
#        reader = csv.reader(tsvfile, delimiter='\t')
#        try:
#            headers = next(reader)
#        except StopIteration:
#            return False, 0  # Empty file
#        rows = list(reader)
#        return (len(rows) > 0), len(rows)

#def write_summary(result, count, prefix, output_path="chimeric_summary.txt"):
#    with open(output_path, "w") as out:
#        out.write(f"ChimericEvidence: {'Yes' if result else 'No'}\n\n")
#        if result:
#            out.write(f"Alert!! Sample {prefix} is found to have some evidence of genetic engineering.\n\n")
#            out.write("1. Access the final report.\n")
#            out.write(f"2. View the evidence of genetic engineering (Number of Chimeric Regions: {count})\n")
#            out.write("3. Download all the outputs from this analysis.\n")
#        else:
#            out.write(f"Sample {prefix} does NOT have any evidence of genetic engineering.\n\n")
#            out.write("1. Assess the final report.\n")
#            out.write("2. Download all the outputs from the analysis.\n")

##### Testing updated function (05/23/2025
#def has_data(filepath):
#    with open(filepath, newline='') as tsvfile:
#        reader = csv.reader(tsvfile, delimiter='\t')
#        try:
#            headers = next(reader)
#        except StopIteration:
#            return False, 0  # Empty file
#        rows = list(reader)
#        return (len(rows) > 0), len(rows)

#def count_lines(filepath):
#    """Return number of non-header lines (excluding blank lines)."""
#    if not os.path.exists(filepath) or os.stat(filepath).st_size == 0:
#        return 0
#    with open(filepath, 'r') as f:
#        lines = [line for line in f if line.strip()]
#        return len(lines) - 1 if lines else 0  # Exclude header

#def write_summary(chimera_result, chimera_count, prefix, vf_file_contigs, ar_file_contigs, vf_file_reads, ar_file_reads, output_path="chimeric_summary.txt"):
#    vf_count_contigs = count_lines(vf_file_contigs)
#    ar_count_contigs = count_lines(ar_file_contigs)
#    vf_count_reads = count_lines(vf_file_reads)
#    ar_count_reads = count_lines(ar_file_reads)
#    with open(output_path, "w") as out:
#        out.write(f"ChimericEvidence: {'Yes' if chimera_result else 'No'}\n\n")

#        if chimera_result:
#            out.write(f"Alert!! Sample {prefix} is found to have some evidence of genetic engineering.\n\n")
#            out.write("1. Access the final report.\n")
#            out.write(f"2. View the evidence of genetic engineering (Number of Chimeric Regions: {chimera_count})\n")
#            out.write("3. Download all the outputs from this analysis.\n")
#        else:
#            out.write(f"Sample {prefix} does NOT have any evidence of genetic engineering.\n\n")
#            out.write("1. Assess the final report.\n")
#            out.write("2. Download all the outputs from the analysis.\n")

#        out.write("\nAdditional Evidence Summary:\n")
#        out.write(f"- Virulence Factors (VF) found in Contigs: {vf_count_contigs} found\n")
#        out.write(f"- Antimicrobial Resistance (AR) Genes found in Contigs: {ar_count_contigs} found\n")
#        out.write(f"- Virulence Factors (VF) found in Reads: {vf_count_reads} found\n")
#        out.write(f"- Antimicrobial Resistance (AR) Genes found in Reads: {ar_count_reads} found\n")

def insert_hit(hit_range, coverage_list:list):
    if len(coverage_list) == 0:
        coverage_list.append(hit_range)
        return coverage_list
    # Loop through coverage ranges for current order_id
    hit_start = hit_range[0]
    hit_end = hit_range[1]
    list_index = 0
    list_length = len(coverage_list)
    last_loop = False
    # for current_range in coverage_list:
    while not last_loop:
        current_range = coverage_list[list_index]
        # Last in coverage list range
        if list_length == list_index + 1:
            last_loop = True
        current_start = current_range[0]
        current_end = current_range[1]
        #---|hit|---|cur_range|---
        if hit_start < hit_end < current_start - 1:
            coverage_list.insert(list_index, hit_range)
            break
        # ---|cur_range|---|hit|---
        if current_end + 1 < hit_start < hit_end:
            if last_loop:
                coverage_list.append(hit_range)
                break
            list_index += 1
            continue
        # ---- | hit  === cur_range | ----
        if hit_start < current_start - 1 <= hit_end <= current_end:
            new_range = [hit_start, current_end]
            coverage_list.pop(list_index)
            coverage_list.insert(list_index, new_range)
            break
        # ---- | cur_range  === hit | ----
        if current_start <= hit_start <= current_end + 1 < hit_end:
            new_range = [current_start, hit_end]
            coverage_list.pop(list_index)
            if last_loop:
                coverage_list.append(new_range)
            # Make new coordinates be the hit coordinates
            hit_range = new_range
            hit_start = hit_range[0]
            hit_end = hit_range[1]
            continue
        # ---- | hit_start === cur_range === hit_end | ----
        if hit_start < current_start and hit_end > current_end:
            coverage_list.pop(list_index)
            if last_loop:
                coverage_list.append(hit_range)
            continue
        if current_start <= hit_start < hit_end <= current_end:
            break
        print("Not matched range %s" % hit_range)

    return coverage_list

class Contig:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.length = len(self.sequence)
        self.order_hits = {}
        self.order_ids = []
        self.family_hits = {}
        self.family_ids = []
        self.genus_hits = {}
        self.genus_ids = []
        self.superkingdom_hits = {}
        self.superkingdom_ids = []
        self.phylum_hits = {}
        self.phylum_ids = []
        self.class_hits = {}
        self.class_ids = []
        self.species_hits = {}
        self.species_ids = []

    def add_blast_hit(self, blast_hit, start_offset=0):
        hit_start = int(blast_hit[start_index]) + start_offset
        hit_end = int(blast_hit[end_index]) + start_offset
        hit_range = [hit_start, hit_end]
        hit_taxid = int(blast_hit[taxid_index])
        desired_taxid_list = get_desired_ranks(hit_taxid, desired_ranks)

        try:
            hit_order_id = desired_taxid_list['order_id']
            if hit_order_id == '<not present>':
                hit_order_id = hit_taxid
        except KeyError as e:# no order id for hit
            hit_order_id = hit_taxid
        try:
            hit_family_id = desired_taxid_list['family_id']
            if hit_family_id == '<not present>':
                hit_family_id = hit_taxid
        except KeyError as e:# no order id for hit
            hit_family_id = hit_taxid
        try:
            hit_genus_id = desired_taxid_list['genus']
            if hit_genus_id == '<not present>':
                hit_genus_id = hit_taxid
        except KeyError as e:# no order id for hit
            hit_genus_id = hit_taxid
        try:
            hit_superkingdom_id = desired_taxid_list['superkingdom']
            if hit_superkingdom_id == '<not present>':
                hit_superkingdom_id = hit_taxid
        except KeyError as e:# no order id for hit
            hit_superkingdom_id = hit_taxid
        try:
            hit_phylum_id = desired_taxid_list['phylum']
            if hit_phylum_id == '<not present>':
                hit_phylum_id = hit_taxid
        except KeyError as e:# no order id for hit
            hit_phylum_id = hit_taxid
        try:
            hit_class_id = desired_taxid_list['class']
            if hit_class_id == '<not present>':
                hit_class_id = hit_taxid
        except KeyError as e:# no order id for hit
            hit_class_id = hit_taxid
        try:
            hit_species_id = desired_taxid_list['species']
            if hit_species_id == '<not present>':
                hit_species_id = hit_taxid
        except KeyError as e:# no order id for hit
            hit_species_id = hit_taxid
        # Order Hits 
        if hit_order_id not in self.order_ids:
            self.order_ids.append(hit_order_id)
            print(ncbi.translate_to_names([hit_order_id]))
            self.order_hits[hit_order_id] = {'first_hit': hit_start, 'last_hit': hit_end, 'coverage': [[hit_start, hit_end]]}
            #return
        order_hits = self.order_hits[hit_order_id]
        first_hit = order_hits['first_hit']
        last_hit = order_hits['last_hit']
        coverage_list = order_hits['coverage']
        if hit_start < first_hit:
            order_hits['first_hit'] = hit_start
        if last_hit < hit_end:
            order_hits['last_hit'] = hit_end
        order_hits['coverage'] = insert_hit(hit_range, coverage_list)
        order_hits = self.order_hits[hit_order_id]
        first_hit = order_hits['first_hit']
        last_hit = order_hits['last_hit']
        coverage_list = order_hits['coverage']
        if hit_start < first_hit:
            order_hits['first_hit'] = hit_start
        if last_hit < hit_end:
            order_hits['last_hit'] = hit_end
        order_hits['coverage'] = insert_hit(hit_range, coverage_list)
        if hit_family_id not in self.family_ids:
            self.family_ids.append(hit_family_id)
            print(ncbi.translate_to_names([hit_family_id]))
            self.family_hits[hit_family_id] = {'first_hit': hit_start, 'last_hit': hit_end, 'coverage': [[hit_start, hit_end]]}
            #return
        family_hits = self.family_hits[hit_family_id]
        first_hit = family_hits['first_hit']
        last_hit = family_hits['last_hit']
        coverage_list = family_hits['coverage']
        if hit_start < first_hit:
            family_hits['first_hit'] = hit_start
        if last_hit < hit_end:
            family_hits['last_hit'] = hit_end
        family_hits['coverage'] = insert_hit(hit_range, coverage_list)
        if hit_genus_id not in self.genus_ids:
            self.genus_ids.append(hit_genus_id)
            print(ncbi.translate_to_names([hit_genus_id]))
            self.genus_hits[hit_genus_id] = {'first_hit': hit_start, 'last_hit': hit_end, 'coverage': [[hit_start, hit_end]]}
        if hit_genus_id not in self.genus_ids:
            self.genus_ids.append(hit_genus_id)
            print(ncbi.translate_to_names([hit_genus_id]))
            self.genus_hits[hit_genus_id] = {'first_hit': hit_start, 'last_hit': hit_end, 'coverage': [[hit_start, hit_end]]}
            #return
        genus_hits = self.genus_hits[hit_genus_id]
        first_hit = genus_hits['first_hit']
        last_hit = genus_hits['last_hit']
        coverage_list = genus_hits['coverage']
        if hit_start < first_hit:
            genus_hits['first_hit'] = hit_start
        if last_hit < hit_end:
            genus_hits['last_hit'] = hit_end
        genus_hits['coverage'] = insert_hit(hit_range, coverage_list)
        if hit_superkingdom_id not in self.superkingdom_ids:
            self.superkingdom_ids.append(hit_superkingdom_id)
            print(ncbi.translate_to_names([hit_superkingdom_id]))
            self.superkingdom_hits[hit_superkingdom_id] = {'first_hit': hit_start, 'last_hit': hit_end, 'coverage': [[hit_start, hit_end]]}
            #return
        superkingdom_hits = self.superkingdom_hits[hit_superkingdom_id]
        first_hit = superkingdom_hits['first_hit']
        last_hit = superkingdom_hits['last_hit']
        coverage_list = superkingdom_hits['coverage']
        if hit_start < first_hit:
            superkingdom_hits['first_hit'] = hit_start
        if last_hit < hit_end:
            superkingdom_hits['last_hit'] = hit_end
        superkingdom_hits['coverage'] = insert_hit(hit_range, coverage_list)
        if hit_phylum_id not in self.phylum_ids:
            self.phylum_ids.append(hit_phylum_id)
            print(ncbi.translate_to_names([hit_phylum_id]))
            self.phylum_hits[hit_phylum_id] = {'first_hit': hit_start, 'last_hit': hit_end, 'coverage': [[hit_start, hit_end]]}
            #return
        phylum_hits = self.phylum_hits[hit_phylum_id]
        first_hit = phylum_hits['first_hit']
        last_hit = phylum_hits['last_hit']
        coverage_list = phylum_hits['coverage']
        if hit_start < first_hit:
            phylum_hits['first_hit'] = hit_start
        if last_hit < hit_end:
            phylum_hits['last_hit'] = hit_end
        phylum_hits['coverage'] = insert_hit(hit_range, coverage_list)
        if hit_class_id not in self.class_ids:
            self.class_ids.append(hit_class_id)
            print(ncbi.translate_to_names([hit_class_id]))
            self.class_hits[hit_class_id] = {'first_hit': hit_start, 'last_hit': hit_end, 'coverage': [[hit_start, hit_end]]}
            #return
        class_hits = self.class_hits[hit_class_id]
        first_hit = class_hits['first_hit']
        last_hit = class_hits['last_hit']
        coverage_list = class_hits['coverage']
        if hit_start < first_hit:
            class_hits['first_hit'] = hit_start
        if last_hit < hit_end:
            class_hits['last_hit'] = hit_end
        class_hits['coverage'] = insert_hit(hit_range, coverage_list)
        if hit_species_id not in self.species_ids:
            self.species_ids.append(hit_species_id)
            print(ncbi.translate_to_names([hit_species_id]))
            self.species_hits[hit_species_id] = {'first_hit': hit_start, 'last_hit': hit_end, 'coverage': [[hit_start, hit_end]]}
            #return
        species_hits = self.species_hits[hit_species_id]
        first_hit = species_hits['first_hit']
        last_hit = species_hits['last_hit']
        coverage_list = species_hits['coverage']
        if hit_start < first_hit:
            species_hits['first_hit'] = hit_start
        if last_hit < hit_end:
            species_hits['last_hit'] = hit_end
        species_hits['coverage'] = insert_hit(hit_range, coverage_list)

    def getChimericRegions(self):
        if desired_rank == "order":
            # List of taxid1, region end, taxid2, region start
            chimeric_region_list = []
            current_taxid_index_dict = {}
            # Initialize index dict
            for order_id in self.order_ids:
                current_taxid_index_dict[order_id] = 0
            # Contigs with no hits
            if not self.order_hits:
                return chimeric_region_list
            first_hit = 0
            previous_hit = 0
            last_hit = 0
            # Determine first hit
            for order_id, order_hits in self.order_hits.items():
                current_first_hit = [order_id, order_hits['coverage'][current_taxid_index_dict[order_id]]]
                current_last_hit = [order_id, order_hits['coverage'][-1]]
                if not first_hit:
                    # current hit is order_id, [start stop] is found by order_hits{'coverage': [hit list index]}
                    first_hit = current_first_hit
                elif current_first_hit[1][0] < first_hit[1][0] :
                    first_hit = current_first_hit
                if not last_hit:
                    last_hit = current_last_hit
                elif current_last_hit[1][1] > last_hit[1][1]:
                    last_hit = current_last_hit
            is_last_hit = False

            not_past_first_hit = True
            while not is_last_hit:
                next_hit = 0
                for order_id, order_hits in self.order_hits.items():
                    try:
                        current_hit = [order_id, order_hits['coverage'][current_taxid_index_dict[order_id]]]
                    except IndexError:
                        # Keep going to the next hit if the it fails
                        continue
                    # Do nothing till the first hit is found
                    if not_past_first_hit:
                        if current_hit[1][0] == first_hit[1][0]:
                            # increment the current hit taxa index
                            # current_taxid_index_dict[current_hit[0]] += 1
                            previous_hit = current_hit
                            not_past_first_hit = False
                            next_hit = current_hit
                            break
                        else:
                            continue
                    # Searching for the next hit in the sequence
                    if not next_hit:
                        next_hit = current_hit
                    elif current_hit[1][0] < next_hit[1][0]:
                        next_hit = current_hit
                # Check if the next and previous hits are of the same order
                try:
                    if next_hit[0] != previous_hit[0]:
                        # length between hits of different taxa
                        length_between_hits = int(next_hit[1][0]) - int(previous_hit[1][1])
                        chimeric_region_list.append([previous_hit[0], previous_hit[1][1], next_hit[0], next_hit[1][0]])
                    current_taxid_index_dict[next_hit[0]] += 1
                    # if the last hit to get out of the while loop
                    if next_hit[1][1] == last_hit[1][1]:
                        is_last_hit = True
                    previous_hit = next_hit
                except TypeError as e:
                    print("ERROR: %s" % e)
                    break
            return chimeric_region_list

        elif desired_rank == "family":
            # List of taxid1, region end, taxid2, region start
            chimeric_region_list = []
            current_taxid_index_dict = {}
            # Initialize index dict
            for family_id in self.family_ids:
                current_taxid_index_dict[family_id] = 0
            # Contigs with no hits
            if not self.family_hits:
                return chimeric_region_list
            first_hit = 0
            previous_hit = 0
            last_hit = 0
            # Determine first hit
            for family_id, family_hits in self.family_hits.items():
                current_first_hit = [family_id, family_hits['coverage'][current_taxid_index_dict[family_id]]]
                current_last_hit = [family_id, family_hits['coverage'][-1]]
                if not first_hit:
                    # current hit is family_id, [start stop] is found by family_hits{'coverage': [hit list index]}
                    first_hit = current_first_hit
                elif current_first_hit[1][0] < first_hit[1][0] :
                    first_hit = current_first_hit
                if not last_hit:
                    last_hit = current_last_hit
                elif current_last_hit[1][1] > last_hit[1][1]:
                    last_hit = current_last_hit
            is_last_hit = False

            not_past_first_hit = True
            while not is_last_hit:
                next_hit = 0
                for family_id, family_hits in self.family_hits.items():
                    try:
                        current_hit = [family_id, family_hits['coverage'][current_taxid_index_dict[family_id]]]
                    except IndexError:
                        # Keep going to the next hit if the it fails
                        continue
                    # Do nothing till the first hit is found
                    if not_past_first_hit:
                        if current_hit[1][0] == first_hit[1][0]:
                            # increment the current hit taxa index
                            # current_taxid_index_dict[current_hit[0]] += 1
                            previous_hit = current_hit
                            not_past_first_hit = False
                            next_hit = current_hit
                            break
                        else:
                            continue
                    # Searching for the next hit in the sequence
                    if not next_hit:
                        next_hit = current_hit
                    elif current_hit[1][0] < next_hit[1][0]:
                        next_hit = current_hit
                # Check if the next and previous hits are of the same family
                try:
                    if next_hit[0] != previous_hit[0]:
                        # length between hits of different taxa
                        length_between_hits = int(next_hit[1][0]) - int(previous_hit[1][1])
                        chimeric_region_list.append([previous_hit[0], previous_hit[1][1], next_hit[0], next_hit[1][0]])
                    current_taxid_index_dict[next_hit[0]] += 1
                    # if the last hit to get out of the while loop
                    if next_hit[1][1] == last_hit[1][1]:
                        is_last_hit = True
                    previous_hit = next_hit
                except TypeError as e:
                    print("ERROR: %s" % e)
                    break
            return chimeric_region_list
        elif desired_rank == "genus":
            # List of taxid1, region end, taxid2, region start
            chimeric_region_list = []
            current_taxid_index_dict = {}
            # Initialize index dict
            for genus_id in self.genus_ids:
                current_taxid_index_dict[genus_id] = 0
            # Contigs with no hits
            if not self.genus_hits:
                return chimeric_region_list
            first_hit = 0
            previous_hit = 0
            last_hit = 0
            # Determine first hit
            for genus_id, genus_hits in self.genus_hits.items():
                current_first_hit = [genus_id, genus_hits['coverage'][current_taxid_index_dict[genus_id]]]
                current_last_hit = [genus_id, genus_hits['coverage'][-1]]
                if not first_hit:
                    # current hit is genus_id, [start stop] is found by genus_hits{'coverage': [hit list index]}
                    first_hit = current_first_hit
                elif current_first_hit[1][0] < first_hit[1][0] :
                    first_hit = current_first_hit
                if not last_hit:
                    last_hit = current_last_hit
                elif current_last_hit[1][1] > last_hit[1][1]:
                    last_hit = current_last_hit
            is_last_hit = False

            not_past_first_hit = True
            while not is_last_hit:
                next_hit = 0
                for genus_id, genus_hits in self.genus_hits.items():
                    try:
                        current_hit = [genus_id, genus_hits['coverage'][current_taxid_index_dict[genus_id]]]
                    except IndexError:
                        # Keep going to the next hit if the it fails
                        continue
                    # Do nothing till the first hit is found
                    if not_past_first_hit:
                        if current_hit[1][0] == first_hit[1][0]:
                            # increment the current hit taxa index
                            # current_taxid_index_dict[current_hit[0]] += 1
                            previous_hit = current_hit
                            not_past_first_hit = False
                            next_hit = current_hit
                            break
                        else:
                            continue
                    # Searching for the next hit in the sequence
                    if not next_hit:
                        next_hit = current_hit
                    elif current_hit[1][0] < next_hit[1][0]:
                        next_hit = current_hit
                # Check if the next and previous hits are of the same genus
                try:
                    if next_hit[0] != previous_hit[0]:
                        # length between hits of different taxa
                        length_between_hits = int(next_hit[1][0]) - int(previous_hit[1][1])
                        chimeric_region_list.append([previous_hit[0], previous_hit[1][1], next_hit[0], next_hit[1][0]])
                    current_taxid_index_dict[next_hit[0]] += 1
                    # if the last hit to get out of the while loop
                    if next_hit[1][1] == last_hit[1][1]:
                        is_last_hit = True
                    previous_hit = next_hit
                except TypeError as e:
                    print("ERROR: %s" % e)
                    break
            return chimeric_region_list
        elif desired_rank == "superkingdom":
            # List of taxid1, region end, taxid2, region start
            chimeric_region_list = []
            current_taxid_index_dict = {}
            # Initialize index dict
            for superkingdom_id in self.superkingdom_ids:
                current_taxid_index_dict[superkingdom_id] = 0
            # Contigs with no hits
            if not self.superkingdom_hits:
                return chimeric_region_list
            first_hit = 0
            previous_hit = 0
            last_hit = 0
            # Determine first hit
            for superkingdom_id, superkingdom_hits in self.superkingdom_hits.items():
                current_first_hit = [superkingdom_id, superkingdom_hits['coverage'][current_taxid_index_dict[superkingdom_id]]]
                current_last_hit = [superkingdom_id, superkingdom_hits['coverage'][-1]]
                if not first_hit:
                    # current hit is superkingdom_id, [start stop] is found by superkingdom_hits{'coverage': [hit list index]}
                    first_hit = current_first_hit
                elif current_first_hit[1][0] < first_hit[1][0] :
                    first_hit = current_first_hit
                if not last_hit:
                    last_hit = current_last_hit
                elif current_last_hit[1][1] > last_hit[1][1]:
                    last_hit = current_last_hit
            is_last_hit = False

            not_past_first_hit = True
            while not is_last_hit:
                next_hit = 0
                for superkingdom_id, superkingdom_hits in self.superkingdom_hits.items():
                    try:
                        current_hit = [superkingdom_id, superkingdom_hits['coverage'][current_taxid_index_dict[superkingdom_id]]]
                    except IndexError:
                        # Keep going to the next hit if the it fails
                        continue
                    # Do nothing till the first hit is found
                    if not_past_first_hit:
                        if current_hit[1][0] == first_hit[1][0]:
                            # increment the current hit taxa index
                            # current_taxid_index_dict[current_hit[0]] += 1
                            previous_hit = current_hit
                            not_past_first_hit = False
                            next_hit = current_hit
                            break
                        else:
                            continue
                    # Searching for the next hit in the sequence
                    if not next_hit:
                        next_hit = current_hit
                    elif current_hit[1][0] < next_hit[1][0]:
                        next_hit = current_hit
                # Check if the next and previous hits are of the same superkingdom
                try:
                    if next_hit[0] != previous_hit[0]:
                        # length between hits of different taxa
                        length_between_hits = int(next_hit[1][0]) - int(previous_hit[1][1])
                        chimeric_region_list.append([previous_hit[0], previous_hit[1][1], next_hit[0], next_hit[1][0]])
                    current_taxid_index_dict[next_hit[0]] += 1
                    # if the last hit to get out of the while loop
                    if next_hit[1][1] == last_hit[1][1]:
                        is_last_hit = True
                    previous_hit = next_hit
                except TypeError as e:
                    print("ERROR: %s" % e)
                    break
        elif desired_rank == "phylum":
            # List of taxid1, region end, taxid2, region start
            chimeric_region_list = []
            current_taxid_index_dict = {}
            # Initialize index dict
            for phylum_id in self.phylum_ids:
                current_taxid_index_dict[phylum_id] = 0
            # Contigs with no hits
            if not self.phylum_hits:
                return chimeric_region_list
            first_hit = 0
            previous_hit = 0
            last_hit = 0
            # Determine first hit
            for phylum_id, phylum_hits in self.phylum_hits.items():
                current_first_hit = [phylum_id, phylum_hits['coverage'][current_taxid_index_dict[phylum_id]]]
                current_last_hit = [phylum_id, phylum_hits['coverage'][-1]]
                if not first_hit:
                    # current hit is phylum_id, [start stop] is found by phylum_hits{'coverage': [hit list index]}
                    first_hit = current_first_hit
                elif current_first_hit[1][0] < first_hit[1][0] :
                    first_hit = current_first_hit
                if not last_hit:
                    last_hit = current_last_hit
                elif current_last_hit[1][1] > last_hit[1][1]:
                    last_hit = current_last_hit
            is_last_hit = False

            not_past_first_hit = True
            while not is_last_hit:
                next_hit = 0
                for phylum_id, phylum_hits in self.phylum_hits.items():
                    try:
                        current_hit = [phylum_id, phylum_hits['coverage'][current_taxid_index_dict[phylum_id]]]
                    except IndexError:
                        # Keep going to the next hit if the it fails
                        continue
                    # Do nothing till the first hit is found
                    if not_past_first_hit:
                        if current_hit[1][0] == first_hit[1][0]:
                            # increment the current hit taxa index
                            # current_taxid_index_dict[current_hit[0]] += 1
                            previous_hit = current_hit
                            not_past_first_hit = False
                            next_hit = current_hit
                            break
                        else:
                            continue
                    # Searching for the next hit in the sequence
                    if not next_hit:
                        next_hit = current_hit
                    elif current_hit[1][0] < next_hit[1][0]:
                        next_hit = current_hit
                # Check if the next and previous hits are of the same phylum
                try:
                    if next_hit[0] != previous_hit[0]:
                        # length between hits of different taxa
                        length_between_hits = int(next_hit[1][0]) - int(previous_hit[1][1])
                        chimeric_region_list.append([previous_hit[0], previous_hit[1][1], next_hit[0], next_hit[1][0]])
                    current_taxid_index_dict[next_hit[0]] += 1
                    # if the last hit to get out of the while loop
                    if next_hit[1][1] == last_hit[1][1]:
                        is_last_hit = True
                    previous_hit = next_hit
                except TypeError as e:
                    print("ERROR: %s" % e)
                    break
        elif desired_rank == "class":
            # List of taxid1, region end, taxid2, region start
            chimeric_region_list = []
            current_taxid_index_dict = {}
            # Initialize index dict
            for class_id in self.class_ids:
                current_taxid_index_dict[class_id] = 0
            # Contigs with no hits
            if not self.class_hits:
                return chimeric_region_list
            first_hit = 0
            previous_hit = 0
            last_hit = 0
            # Determine first hit
            for class_id, class_hits in self.class_hits.items():
                current_first_hit = [class_id, class_hits['coverage'][current_taxid_index_dict[class_id]]]
                current_last_hit = [class_id, class_hits['coverage'][-1]]
                if not first_hit:
                    # current hit is class_id, [start stop] is found by class_hits{'coverage': [hit list index]}
                    first_hit = current_first_hit
                elif current_first_hit[1][0] < first_hit[1][0] :
                    first_hit = current_first_hit
                if not last_hit:
                    last_hit = current_last_hit
                elif current_last_hit[1][1] > last_hit[1][1]:
                    last_hit = current_last_hit
            is_last_hit = False

            not_past_first_hit = True
            while not is_last_hit:
                next_hit = 0
                for class_id, class_hits in self.class_hits.items():
                    try:
                        current_hit = [class_id, class_hits['coverage'][current_taxid_index_dict[class_id]]]
                    except IndexError:
                        # Keep going to the next hit if the it fails
                        continue
                    # Do nothing till the first hit is found
                    if not_past_first_hit:
                        if current_hit[1][0] == first_hit[1][0]:
                            # increment the current hit taxa index
                            # current_taxid_index_dict[current_hit[0]] += 1
                            previous_hit = current_hit
                            not_past_first_hit = False
                            next_hit = current_hit
                            break
                        else:
                            continue
                    # Searching for the next hit in the sequence
                    if not next_hit:
                        next_hit = current_hit
                    elif current_hit[1][0] < next_hit[1][0]:
                        next_hit = current_hit
                # Check if the next and previous hits are of the same class
                try:
                    if next_hit[0] != previous_hit[0]:
                        # length between hits of different taxa
                        length_between_hits = int(next_hit[1][0]) - int(previous_hit[1][1])
                        chimeric_region_list.append([previous_hit[0], previous_hit[1][1], next_hit[0], next_hit[1][0]])
                    current_taxid_index_dict[next_hit[0]] += 1
                    # if the last hit to get out of the while loop
                    if next_hit[1][1] == last_hit[1][1]:
                        is_last_hit = True
                    previous_hit = next_hit
                except TypeError as e:
                    print("ERROR: %s" % e)
                    break
        elif desired_rank == "species":
            # List of taxid1, region end, taxid2, region start
            chimeric_region_list = []
            current_taxid_index_dict = {}
            # Initialize index dict
            for species_id in self.species_ids:
                current_taxid_index_dict[species_id] = 0
            # Contigs with no hits
            if not self.species_hits:
                return chimeric_region_list
            first_hit = 0
            previous_hit = 0
            last_hit = 0
            # Determine first hit
            for species_id, species_hits in self.species_hits.items():
                current_first_hit = [species_id, species_hits['coverage'][current_taxid_index_dict[species_id]]]
                current_last_hit = [species_id, species_hits['coverage'][-1]]
                if not first_hit:
                    # current hit is species_id, [start stop] is found by species_hits{'coverage': [hit list index]}
                    first_hit = current_first_hit
                elif current_first_hit[1][0] < first_hit[1][0] :
                    first_hit = current_first_hit
                if not last_hit:
                    last_hit = current_last_hit
                elif current_last_hit[1][1] > last_hit[1][1]:
                    last_hit = current_last_hit
            is_last_hit = False

            not_past_first_hit = True
            while not is_last_hit:
                next_hit = 0
                for species_id, species_hits in self.species_hits.items():
                    try:
                        current_hit = [species_id, species_hits['coverage'][current_taxid_index_dict[species_id]]]
                    except IndexError:
                        # Keep going to the next hit if the it fails
                        continue
                    # Do nothing till the first hit is found
                    if not_past_first_hit:
                        if current_hit[1][0] == first_hit[1][0]:
                            # increment the current hit taxa index
                            # current_taxid_index_dict[current_hit[0]] += 1
                            previous_hit = current_hit
                            not_past_first_hit = False
                            next_hit = current_hit
                            break
                        else:
                            continue
                    # Searching for the next hit in the sequence
                    if not next_hit:
                        next_hit = current_hit
                    elif current_hit[1][0] < next_hit[1][0]:
                        next_hit = current_hit
                # Check if the next and previous hits are of the same species
                try:
                    if next_hit[0] != previous_hit[0]:
                        # length between hits of different taxa
                        length_between_hits = int(next_hit[1][0]) - int(previous_hit[1][1])
                        chimeric_region_list.append([previous_hit[0], previous_hit[1][1], next_hit[0], next_hit[1][0]])
                    current_taxid_index_dict[next_hit[0]] += 1
                    # if the last hit to get out of the while loop
                    if next_hit[1][1] == last_hit[1][1]:
                        is_last_hit = True
                    previous_hit = next_hit
                except TypeError as e:
                    print("ERROR: %s" % e)
                    break
            return chimeric_region_list
        else:
            print("Desired Rank Doesn't Exist")
# from https://bioinformatics.stackexchange.com/questions/5540/python-scripting-with-ete3-to-query-ncbis-taxonomy-sqlite3-warning-can-only
def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    return{'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

def writeToGFF(gff_file, contig_name, feature_type, start, stop, feature_id, match_count=0):
    feature_name = '%s_%s' % (feature_id, match_count)
    metadata = 'ID=%s' % feature_name

    # S3_Cedar_Virus	.	Primer	4641	4661	.	-	.	Name=R1
    if stop < start:
        # Always have the smaller number first
        line = (contig_name, feature_type, stop, start, '-', metadata)
    else:
        line = (contig_name, feature_type, start, stop, '+', metadata)
    gff_line = '%s\tEDGE\t%s\t%s\t%s\t.\t%s\t.\t%s\n' % line
    gff_file.write(gff_line)

def runBlast(query_fasta, blast_db, output_file, num_cpu=4):
    blast_command = NcbiblastnCommandline(
        query=query_fasta,
        db=blast_db,
        outfmt='6 qseqid qlen qstart qend sacc staxid scomnames salltitles slen sstart send pident evalue bitscore length',
        num_threads=num_cpu,
        out=output_file,
        max_target_seqs=MAX_TARGET_SEQS,
        culling_limit=1
    )
    print('running blast: %s' % blast_command)
    if not path.isfile(output_file):
        blast_command()
    # Put contigs into a seq list
    return

def processBlastResults(blast_results_file, contigs_dict, is_gap=False):
    try:
        blast_output_file = open(blast_results_file, 'r')
    except IOError as e:
        exit('ERROR could not open blast output file %s. %s' % (blast_results_file, e))

    try:
        for blast_line in blast_output_file:
            blast_line = blast_line.rstrip()
            # 0:qseqid 1:qlen 2:qstart 3:qend 4:sacc 5:staxid 6:scomnames 7:salltitles 8:slen 9:sstart 10:send 11:pident 12:evalue 13:bitscore 14:length
            split_blast_line = blast_line.split('\t')
            blast_query_contig = split_blast_line[0]
            start_offset = 0
            contigs_dict[blast_query_contig].add_blast_hit(split_blast_line, start_offset=start_offset)
    except KeyError as e:
        exit('result not in contig')

    return contigs_dict

def main(args):
    contigs_file_name = args.contigs
    prefix = args.prefix
    out_dir = args.outDir
    blast_db1 = args.blast_db1
    num_cpu = args.num_cpu
    global desired_rank
    desired_rank = args.desired_rank
    global proj_dir
    proj_dir = args.project_dir
    #desired_ranks = ' '.join(desired_ranks)
    #desired_ranks = desired_ranks.split(',')
    global MAX_TARGET_SEQS
    MAX_TARGET_SEQS = args.max_target_seqs
    # Check if output dir exists if not create it
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    contigs_dict = {}
    try:
        for sequence in SeqIO.parse(contigs_file_name, 'fasta'):
            contigs_dict[sequence.id] = Contig(sequence.id, sequence.seq)
    except IOError as e:
        exit('ERROR: could not open contig fasta file %s. %s' % (contigs_file_name, e)) 
    # Run the initial BLAST
    print("Running Initial BLAST against RefSeq")
    blast_output_file_name = path.join(out_dir, ('%s_blast_results.tsv' % prefix))
    runBlast(contigs_file_name, blast_db1, blast_output_file_name, num_cpu)
    print("Done running initial BLAST against RefSeq")
    print("Processing initial BLAST results")
    # Add hit coverage to Contig
    contigs_dict = processBlastResults(blast_output_file_name, contigs_dict)
    print("Done processing initial BLAST results")
    if desired_rank == "order":
        # Output Files
        hits_gff_file_name = path.join(out_dir, ('%s_hits.gff' % prefix))
        chimeric_fasta_file_name = path.join(out_dir, ('%s_chimeric_regions.fasta' % prefix))
        chimeric_gff_file_name = path.join(out_dir, ('%s_chimeric_regions.gff' % prefix))
        chimeric_regions_table_file_name = path.join(out_dir, ('%s_chimeric_regions.tsv' % prefix))
        try:
            hits_gff_file = open(hits_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open gff file for writing %s: %s' % (hits_gff_file_name, e))
        hits_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_fasta_file = open(chimeric_fasta_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric fasta file for writing %s: %s' % (chimeric_fasta_file_name, e))
        try:
            chimeric_gff_file = open(chimeric_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric gff file for writing %s: %s' % (chimeric_gff_file_name, e))
        chimeric_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_regions_table_file = open(chimeric_regions_table_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric regions file for writing %s: %s' % (chimeric_regions_table_file_name, e))
        chimeric_regions_table_file.write('chimeric_count\tcontig_name\tfirst_order_id\tsecond_order_id\tregion_start\t'
                                          'region_end\n')
    elif desired_rank == "family":
        # Output Files
        hits_gff_file_name = path.join(out_dir, ('%s_hits.gff' % prefix))
        chimeric_fasta_file_name = path.join(out_dir, ('%s_chimeric_regions.fasta' % prefix))
        chimeric_gff_file_name = path.join(out_dir, ('%s_chimeric_regions.gff' % prefix))
        chimeric_regions_table_file_name = path.join(out_dir, ('%s_chimeric_regions.tsv' % prefix))
        try:
            hits_gff_file = open(hits_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open gff file for writing %s: %s' % (hits_gff_file_name, e))
        hits_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_fasta_file = open(chimeric_fasta_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric fasta file for writing %s: %s' % (chimeric_fasta_file_name, e))
        try:
            chimeric_gff_file = open(chimeric_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric gff file for writing %s: %s' % (chimeric_gff_file_name, e))
        chimeric_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_regions_table_file = open(chimeric_regions_table_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric regions file for writing %s: %s' % (chimeric_regions_table_file_name, e))
        chimeric_regions_table_file.write('chimeric_count\tcontig_name\tfirst_family_id\tsecond_family_id\tregion_start\t'
                                          'region_end\n')
    elif desired_rank == "genus":
        # Output Files
        hits_gff_file_name = path.join(out_dir, ('%s_hits.gff' % prefix))
        chimeric_fasta_file_name = path.join(out_dir, ('%s_chimeric_regions.fasta' % prefix))
        chimeric_gff_file_name = path.join(out_dir, ('%s_chimeric_regions.gff' % prefix))
        chimeric_regions_table_file_name = path.join(out_dir, ('%s_chimeric_regions.tsv' % prefix))
        try:
            hits_gff_file = open(hits_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open gff file for writing %s: %s' % (hits_gff_file_name, e))
        hits_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_fasta_file = open(chimeric_fasta_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric fasta file for writing %s: %s' % (chimeric_fasta_file_name, e))
        try:
            chimeric_gff_file = open(chimeric_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric gff file for writing %s: %s' % (chimeric_gff_file_name, e))
        chimeric_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_regions_table_file = open(chimeric_regions_table_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric regions file for writing %s: %s' % (chimeric_regions_table_file_name, e))
        chimeric_regions_table_file.write('chimeric_count\tcontig_name\tfirst_genus_id\tsecond_genus_id\tregion_start\t'
                                          'region_end\n')
    elif desired_ranks == "superkingdom":
        # Output Files
        hits_gff_file_name = path.join(out_dir, ('%s_hits.gff' % prefix))
        chimeric_fasta_file_name = path.join(out_dir, ('%s_chimeric_regions.fasta' % prefix))
        chimeric_gff_file_name = path.join(out_dir, ('%s_chimeric_regions.gff' % prefix))
        chimeric_regions_table_file_name = path.join(out_dir, ('%s_chimeric_regions.tsv' % prefix))
        try:
            hits_gff_file = open(hits_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open gff file for writing %s: %s' % (hits_gff_file_name, e))
        hits_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_fasta_file = open(chimeric_fasta_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric fasta file for writing %s: %s' % (chimeric_fasta_file_name, e))
        try:
            chimeric_gff_file = open(chimeric_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric gff file for writing %s: %s' % (chimeric_gff_file_name, e))
        chimeric_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_regions_table_file = open(chimeric_regions_table_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric regions file for writing %s: %s' % (chimeric_regions_table_file_name, e))
        chimeric_regions_table_file.write('chimeric_count\tcontig_name\tfirst_superkingdom_id\tsecond_superkingdom_id\tregion_start\t'
                                          'region_end\n')
    elif desired_rank == "phylum":
        # Output Files
        hits_gff_file_name = path.join(out_dir, ('%s_hits.gff' % prefix))
        chimeric_fasta_file_name = path.join(out_dir, ('%s_chimeric_regions.fasta' % prefix))
        chimeric_gff_file_name = path.join(out_dir, ('%s_chimeric_regions.gff' % prefix))
        chimeric_regions_table_file_name = path.join(out_dir, ('%s_chimeric_regions.tsv' % prefix))
        try:
            hits_gff_file = open(hits_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open gff file for writing %s: %s' % (hits_gff_file_name, e))
        hits_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_fasta_file = open(chimeric_fasta_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric fasta file for writing %s: %s' % (chimeric_fasta_file_name, e))
        try:
            chimeric_gff_file = open(chimeric_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric gff file for writing %s: %s' % (chimeric_gff_file_name, e))
        chimeric_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_regions_table_file = open(chimeric_regions_table_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric regions file for writing %s: %s' % (chimeric_regions_table_file_name, e))
        chimeric_regions_table_file.write('chimeric_count\tcontig_name\tfirst_phylum_id\tsecond_phylum_id\tregion_start\t'
                                          'region_end\n')
    elif desired_rank == "class":
        # Output Files
        hits_gff_file_name = path.join(out_dir, ('%s_hits.gff' % prefix))
        chimeric_fasta_file_name = path.join(out_dir, ('%s_chimeric_regions.fasta' % prefix))
        chimeric_gff_file_name = path.join(out_dir, ('%s_chimeric_regions.gff' % prefix))
        chimeric_regions_table_file_name = path.join(out_dir, ('%s_chimeric_regions.tsv' % prefix))
        try:
            hits_gff_file = open(hits_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open gff file for writing %s: %s' % (hits_gff_file_name, e))
        hits_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_fasta_file = open(chimeric_fasta_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric fasta file for writing %s: %s' % (chimeric_fasta_file_name, e))
        try:
            chimeric_gff_file = open(chimeric_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric gff file for writing %s: %s' % (chimeric_gff_file_name, e))
        chimeric_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_regions_table_file = open(chimeric_regions_table_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric regions file for writing %s: %s' % (chimeric_regions_table_file_name, e))
        chimeric_regions_table_file.write('chimeric_count\tcontig_name\tfirst_class_id\tsecond_class_id\tregion_start\t'
                                          'region_end\n')
    elif desired_rank == "species":
        # Output Files
        hits_gff_file_name = path.join(out_dir, ('%s_hits.gff' % prefix))
        chimeric_fasta_file_name = path.join(out_dir, ('%s_chimeric_regions.fasta' % prefix))
        chimeric_gff_file_name = path.join(out_dir, ('%s_chimeric_regions.gff' % prefix))
        chimeric_regions_table_file_name = path.join(out_dir, ('%s_chimeric_regions.tsv' % prefix))
        try:
            hits_gff_file = open(hits_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open gff file for writing %s: %s' % (hits_gff_file_name, e))
        hits_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_fasta_file = open(chimeric_fasta_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric fasta file for writing %s: %s' % (chimeric_fasta_file_name, e))
        try:
            chimeric_gff_file = open(chimeric_gff_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric gff file for writing %s: %s' % (chimeric_gff_file_name, e))
        chimeric_gff_file.write('##gff-version 3.2.1\n')
        try:
            chimeric_regions_table_file = open(chimeric_regions_table_file_name, 'w')
        except IOError as e:
            exit('ERROR: Could not open chimeric regions file for writing %s: %s' % (chimeric_regions_table_file_name, e))
        chimeric_regions_table_file.write('chimeric_count\tcontig_name\tfirst_species_id\tsecond_species_id\tregion_start\t'
                                          'region_end\n')
    else:
        print("Taxonomic Rank does not exist.")
    if desired_rank == "order":
        print("Identifying chimeric regions based on taxonomic order")
        for contig_name, contig in contigs_dict.items():
            contig_length = len(contig.sequence)
            hits_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_region_list = contig.getChimericRegions()

              # Print hits gff file
            for taxid, hits in contig.order_hits.items():
                hit_count = 0
                tax_order = ncbi.translate_to_names([taxid])[0]
                for hit in hits['coverage']:
                    hit_start = hit[0]
                    hit_end = hit[1]
                    writeToGFF(hits_gff_file, contig_name, tax_order, hit_start, hit_end, tax_order, hit_count)
                    hit_count += 1
            # Write chimeric regions
            # List of taxid1, region end, taxid2, region start
            chimeric_count = 0
            for chimeric_region in chimeric_region_list:
                first_order_id = ncbi.translate_to_names([chimeric_region[0]])[0]
                region_start = int(chimeric_region[1])
                second_order_id = ncbi.translate_to_names([chimeric_region[2]])[0]
                region_end = int(chimeric_region[3])
                region_length = region_end - region_start
                # Goal is to have 300bp regions
                if region_length < LEN_BETWEEN_TAXA:
                    adjustment = 150 - int(region_length/2)
                    region_start -= adjustment
                    region_end += adjustment
                    header = '>%s_order_%s_%s_position_%s_%s' % (
                    contig_name, first_order_id, second_order_id, region_start, region_end)
                    gff_feature_id = '%s_order_%s_%s_position_%s_%s' % (
                    contig_name, first_order_id, second_order_id, region_start, region_end)
                    region_sequence = contig.sequence[region_start:region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', region_start, region_end, gff_feature_id,
                               chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_order_id, second_order_id, region_start, region_end))
                    chimeric_count += 1
                  # Divide into two regions, one at each junction
                else:
                    adjustment = int(LEN_BETWEEN_TAXA / 2)
                    first_region_start = region_start - adjustment
                    first_region_end = region_start + adjustment
                    second_region_start = region_end - adjustment
                    second_region_end = region_end + adjustment

                    header = '>%s_order_%s_%s_position_%s_%s' % (
                        contig_name, first_order_id, second_order_id, first_region_start, first_region_end)
                    gff_feature_id = '%s_order_%s_%s_position_%s_%s' % (
                        contig_name, first_order_id, second_order_id, first_region_start, first_region_end)
                    region_sequence = contig.sequence[first_region_start:first_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', first_region_start, first_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_order_id, second_order_id, region_start, region_end))
                    chimeric_count += 1
                    header = '>%s_order_%s_%s_position_%s_%s' % (
                        contig_name, first_order_id, second_order_id, second_region_start, second_region_end)
                    gff_feature_id = '%s_order_%s_%s_position_%s_%s' % (
                        contig_name, first_order_id, second_order_id, second_region_start, second_region_end)
                    region_sequence = contig.sequence[second_region_start:second_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', second_region_start, second_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_order_id, second_order_id, region_start, region_end))
                    chimeric_count += 1

    elif desired_rank == "family":
        print("Identifying chimeric regions based on taxonomic family")
        for contig_name, contig in contigs_dict.items():
            contig_length = len(contig.sequence)
            hits_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_region_list = contig.getChimericRegions()

              # Print hits gff file
            for taxid, hits in contig.family_hits.items():
                hit_count = 0
                tax_family = ncbi.translate_to_names([taxid])[0]
                for hit in hits['coverage']:
                    hit_start = hit[0]
                    hit_end = hit[1]
                    writeToGFF(hits_gff_file, contig_name, tax_family, hit_start, hit_end, tax_family, hit_count)
                    hit_count += 1
            # Write chimeric regions
            # List of taxid1, region end, taxid2, region start
            chimeric_count = 0
            for chimeric_region in chimeric_region_list:
                first_family_id = ncbi.translate_to_names([chimeric_region[0]])[0]
                region_start = int(chimeric_region[1])
                second_family_id = ncbi.translate_to_names([chimeric_region[2]])[0]
                region_end = int(chimeric_region[3])
                region_length = region_end - region_start
                # Goal is to have 300bp regions
                if region_length < LEN_BETWEEN_TAXA:
                    adjustment = 150 - int(region_length/2)
                    region_start -= adjustment
                    region_end += adjustment
                    header = '>%s_family_%s_%s_position_%s_%s' % (
                    contig_name, first_family_id, second_family_id, region_start, region_end)
                    gff_feature_id = '%s_family_%s_%s_position_%s_%s' % (
                    contig_name, first_family_id, second_family_id, region_start, region_end)
                    region_sequence = contig.sequence[region_start:region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', region_start, region_end, gff_feature_id,
                               chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_family_id, second_family_id, region_start, region_end))
                    chimeric_count += 1
                  # Divide into two regions, one at each junction
                else:
                    adjustment = int(LEN_BETWEEN_TAXA / 2)
                    first_region_start = region_start - adjustment
                    first_region_end = region_start + adjustment
                    second_region_start = region_end - adjustment
                    second_region_end = region_end + adjustment

                    header = '>%s_family_%s_%s_position_%s_%s' % (
                        contig_name, first_family_id, second_family_id, first_region_start, first_region_end)
                    gff_feature_id = '%s_family_%s_%s_position_%s_%s' % (
                        contig_name, first_family_id, second_family_id, first_region_start, first_region_end)
                    region_sequence = contig.sequence[first_region_start:first_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', first_region_start, first_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_family_id, second_family_id, region_start, region_end))
                    chimeric_count += 1
                    header = '>%s_family_%s_%s_position_%s_%s' % (
                        contig_name, first_family_id, second_family_id, second_region_start, second_region_end)
                    gff_feature_id = '%s_family_%s_%s_position_%s_%s' % (
                        contig_name, first_family_id, second_family_id, second_region_start, second_region_end)
                    region_sequence = contig.sequence[second_region_start:second_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', second_region_start, second_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_family_id, second_family_id, region_start, region_end))
                    chimeric_count += 1

    elif desired_rank == "genus":
        print("Identifying chimeric regions based on taxonomic genus")
        for contig_name, contig in contigs_dict.items():
            contig_length = len(contig.sequence)
            hits_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_region_list = contig.getChimericRegions()

              # Print hits gff file
            for taxid, hits in contig.genus_hits.items():
                hit_count = 0
                tax_genus = ncbi.translate_to_names([taxid])[0]
                for hit in hits['coverage']:
                    hit_start = hit[0]
                    hit_end = hit[1]
                    writeToGFF(hits_gff_file, contig_name, tax_genus, hit_start, hit_end, tax_genus, hit_count)
                    hit_count += 1
            # Write chimeric regions
            # List of taxid1, region end, taxid2, region start
            chimeric_count = 0
            for chimeric_region in chimeric_region_list:
                first_genus_id = ncbi.translate_to_names([chimeric_region[0]])[0]
                region_start = int(chimeric_region[1])
                second_genus_id = ncbi.translate_to_names([chimeric_region[2]])[0]
                region_end = int(chimeric_region[3])
                region_length = region_end - region_start
                # Goal is to have 300bp regions
                if region_length < LEN_BETWEEN_TAXA:
                    adjustment = 150 - int(region_length/2)
                    region_start -= adjustment
                    region_end += adjustment
                    header = '>%s_genus_%s_%s_position_%s_%s' % (
                    contig_name, first_genus_id, second_genus_id, region_start, region_end)
                    gff_feature_id = '%s_genus_%s_%s_position_%s_%s' % (
                    contig_name, first_genus_id, second_genus_id, region_start, region_end)
                    region_sequence = contig.sequence[region_start:region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', region_start, region_end, gff_feature_id,
                               chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_genus_id, second_genus_id, region_start, region_end))
                    chimeric_count += 1
                  # Divide into two regions, one at each junction
                else:
                    adjustment = int(LEN_BETWEEN_TAXA / 2)
                    first_region_start = region_start - adjustment
                    first_region_end = region_start + adjustment
                    second_region_start = region_end - adjustment
                    second_region_end = region_end + adjustment

                    header = '>%s_genus_%s_%s_position_%s_%s' % (
                        contig_name, first_genus_id, second_genus_id, first_region_start, first_region_end)
                    gff_feature_id = '%s_genus_%s_%s_position_%s_%s' % (
                        contig_name, first_genus_id, second_genus_id, first_region_start, first_region_end)
                    region_sequence = contig.sequence[first_region_start:first_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', first_region_start, first_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_genus_id, second_genus_id, region_start, region_end))
                    chimeric_count += 1
                    header = '>%s_genus_%s_%s_position_%s_%s' % (
                        contig_name, first_genus_id, second_genus_id, second_region_start, second_region_end)
                    gff_feature_id = '%s_genus_%s_%s_position_%s_%s' % (
                        contig_name, first_genus_id, second_genus_id, second_region_start, second_region_end)
                    region_sequence = contig.sequence[second_region_start:second_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', second_region_start, second_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_genus_id, second_genus_id, region_start, region_end))
                    chimeric_count += 1

    elif desired_rank == "superkingdom":
        print("Identifying chimeric regions based on taxonomic superkingdom")
        for contig_name, contig in contigs_dict.items():
            contig_length = len(contig.sequence)
            hits_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_region_list = contig.getChimericRegions()

              # Print hits gff file
            for taxid, hits in contig.superkingdom_hits.items():
                hit_count = 0
                tax_superkingdom = ncbi.translate_to_names([taxid])[0]
                for hit in hits['coverage']:
                    hit_start = hit[0]
                    hit_end = hit[1]
                    writeToGFF(hits_gff_file, contig_name, tax_superkingdom, hit_start, hit_end, tax_superkingdom, hit_count)
                    hit_count += 1
            # Write chimeric regions
            # List of taxid1, region end, taxid2, region start
            chimeric_count = 0
            for chimeric_region in chimeric_region_list:
                first_superkingdom_id = ncbi.translate_to_names([chimeric_region[0]])[0]
                region_start = int(chimeric_region[1])
                second_superkingdom_id = ncbi.translate_to_names([chimeric_region[2]])[0]
                region_end = int(chimeric_region[3])
                region_length = region_end - region_start
                # Goal is to have 300bp regions
                if region_length < LEN_BETWEEN_TAXA:
                    adjustment = 150 - int(region_length/2)
                    region_start -= adjustment
                    region_end += adjustment
                    header = '>%s_superkingdom_%s_%s_position_%s_%s' % (
                    contig_name, first_superkingdom_id, second_superkingdom_id, region_start, region_end)
                    gff_feature_id = '%s_superkingdom_%s_%s_position_%s_%s' % (
                    contig_name, first_superkingdom_id, second_superkingdom_id, region_start, region_end)
                    region_sequence = contig.sequence[region_start:region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', region_start, region_end, gff_feature_id,
                               chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_superkingdom_id, second_superkingdom_id, region_start, region_end))
                    chimeric_count += 1
                  # Divide into two regions, one at each junction
                else:
                    adjustment = int(LEN_BETWEEN_TAXA / 2)
                    first_region_start = region_start - adjustment
                    first_region_end = region_start + adjustment
                    second_region_start = region_end - adjustment
                    second_region_end = region_end + adjustment

                    header = '>%s_superkingdom_%s_%s_position_%s_%s' % (
                        contig_name, first_superkingdom_id, second_superkingdom_id, first_region_start, first_region_end)
                    gff_feature_id = '%s_superkingdom_%s_%s_position_%s_%s' % (
                        contig_name, first_superkingdom_id, second_superkingdom_id, first_region_start, first_region_end)
                    region_sequence = contig.sequence[first_region_start:first_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', first_region_start, first_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_superkingdom_id, second_superkingdom_id, region_start, region_end))
                    chimeric_count += 1
                    header = '>%s_superkingdom_%s_%s_position_%s_%s' % (
                        contig_name, first_superkingdom_id, second_superkingdom_id, second_region_start, second_region_end)
                    gff_feature_id = '%s_superkingdom_%s_%s_position_%s_%s' % (
                        contig_name, first_superkingdom_id, second_superkingdom_id, second_region_start, second_region_end)
                    region_sequence = contig.sequence[second_region_start:second_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', second_region_start, second_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_superkingdom_id, second_superkingdom_id, region_start, region_end))
                    chimeric_count += 1

    elif desired_rank == "phylum":
        print("Identifying chimeric regions based on taxonomic phylum")
        for contig_name, contig in contigs_dict.items():
            contig_length = len(contig.sequence)
            hits_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_region_list = contig.getChimericRegions()

              # Print hits gff file
            for taxid, hits in contig.phylum_hits.items():
                hit_count = 0
                tax_phylum = ncbi.translate_to_names([taxid])[0]
                for hit in hits['coverage']:
                    hit_start = hit[0]
                    hit_end = hit[1]
                    writeToGFF(hits_gff_file, contig_name, tax_phylum, hit_start, hit_end, tax_phylum, hit_count)
                    hit_count += 1
            # Write chimeric regions
            # List of taxid1, region end, taxid2, region start
            chimeric_count = 0
            for chimeric_region in chimeric_region_list:
                first_phylum_id = ncbi.translate_to_names([chimeric_region[0]])[0]
                region_start = int(chimeric_region[1])
                second_phylum_id = ncbi.translate_to_names([chimeric_region[2]])[0]
                region_end = int(chimeric_region[3])
                region_length = region_end - region_start
                # Goal is to have 300bp regions
                if region_length < LEN_BETWEEN_TAXA:
                    adjustment = 150 - int(region_length/2)
                    region_start -= adjustment
                    region_end += adjustment
                    header = '>%s_phylum_%s_%s_position_%s_%s' % (
                    contig_name, first_phylum_id, second_phylum_id, region_start, region_end)
                    gff_feature_id = '%s_phylum_%s_%s_position_%s_%s' % (
                    contig_name, first_phylum_id, second_phylum_id, region_start, region_end)
                    region_sequence = contig.sequence[region_start:region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', region_start, region_end, gff_feature_id,
                               chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_phylum_id, second_phylum_id, region_start, region_end))
                    chimeric_count += 1
                  # Divide into two regions, one at each junction
                else:
                    adjustment = int(LEN_BETWEEN_TAXA / 2)
                    first_region_start = region_start - adjustment
                    first_region_end = region_start + adjustment
                    second_region_start = region_end - adjustment
                    second_region_end = region_end + adjustment

                    header = '>%s_phylum_%s_%s_position_%s_%s' % (
                        contig_name, first_phylum_id, second_phylum_id, first_region_start, first_region_end)
                    gff_feature_id = '%s_phylum_%s_%s_position_%s_%s' % (
                        contig_name, first_phylum_id, second_phylum_id, first_region_start, first_region_end)
                    region_sequence = contig.sequence[first_region_start:first_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', first_region_start, first_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_phylum_id, second_phylum_id, region_start, region_end))
                    chimeric_count += 1
                    header = '>%s_phylum_%s_%s_position_%s_%s' % (
                        contig_name, first_phylum_id, second_phylum_id, second_region_start, second_region_end)
                    gff_feature_id = '%s_phylum_%s_%s_position_%s_%s' % (
                        contig_name, first_phylum_id, second_phylum_id, second_region_start, second_region_end)
                    region_sequence = contig.sequence[second_region_start:second_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', second_region_start, second_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_phylum_id, second_phylum_id, region_start, region_end))
                    chimeric_count += 1

    elif desired_rank == "class":
        print("Identifying chimeric regions based on taxonomic class")
        for contig_name, contig in contigs_dict.items():
            contig_length = len(contig.sequence)
            hits_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_region_list = contig.getChimericRegions()

              # Print hits gff file
            for taxid, hits in contig.class_hits.items():
                hit_count = 0
                tax_class = ncbi.translate_to_names([taxid])[0]
                for hit in hits['coverage']:
                    hit_start = hit[0]
                    hit_end = hit[1]
                    writeToGFF(hits_gff_file, contig_name, tax_class, hit_start, hit_end, tax_class, hit_count)
                    hit_count += 1
            # Write chimeric regions
            # List of taxid1, region end, taxid2, region start
            chimeric_count = 0
            for chimeric_region in chimeric_region_list:
                first_class_id = ncbi.translate_to_names([chimeric_region[0]])[0]
                region_start = int(chimeric_region[1])
                second_class_id = ncbi.translate_to_names([chimeric_region[2]])[0]
                region_end = int(chimeric_region[3])
                region_length = region_end - region_start
                # Goal is to have 300bp regions
                if region_length < LEN_BETWEEN_TAXA:
                    adjustment = 150 - int(region_length/2)
                    region_start -= adjustment
                    region_end += adjustment
                    header = '>%s_class_%s_%s_position_%s_%s' % (
                    contig_name, first_class_id, second_class_id, region_start, region_end)
                    gff_feature_id = '%s_class_%s_%s_position_%s_%s' % (
                    contig_name, first_class_id, second_class_id, region_start, region_end)
                    region_sequence = contig.sequence[region_start:region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', region_start, region_end, gff_feature_id,
                               chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_class_id, second_class_id, region_start, region_end))
                    chimeric_count += 1

                  # Divide into two regions, one at each junction
                else:
                    adjustment = int(LEN_BETWEEN_TAXA / 2)
                    first_region_start = region_start - adjustment
                    first_region_end = region_start + adjustment
                    second_region_start = region_end - adjustment
                    second_region_end = region_end + adjustment

                    header = '>%s_class_%s_%s_position_%s_%s' % (
                        contig_name, first_class_id, second_class_id, first_region_start, first_region_end)
                    gff_feature_id = '%s_class_%s_%s_position_%s_%s' % (
                        contig_name, first_class_id, second_class_id, first_region_start, first_region_end)
                    region_sequence = contig.sequence[first_region_start:first_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', first_region_start, first_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_class_id, second_class_id, region_start, region_end))
                    chimeric_count += 1
                    header = '>%s_class_%s_%s_position_%s_%s' % (
                        contig_name, first_class_id, second_class_id, second_region_start, second_region_end)
                    gff_feature_id = '%s_class_%s_%s_position_%s_%s' % (
                        contig_name, first_class_id, second_class_id, second_region_start, second_region_end)
                    region_sequence = contig.sequence[second_region_start:second_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', second_region_start, second_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_class_id, second_class_id, region_start, region_end))
                    chimeric_count += 1

    elif desired_rank == "species":
        print("Identifying chimeric regions based on taxonomic species")
        for contig_name, contig in contigs_dict.items():
            contig_length = len(contig.sequence)
            hits_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_gff_file.write('##sequence-region %s 1 %s\n' % (contig_name, contig_length))
            chimeric_region_list = contig.getChimericRegions()

              # Print hits gff file
            for taxid, hits in contig.species_hits.items():
                hit_count = 0
                tax_species = ncbi.translate_to_names([taxid])[0]
                for hit in hits['coverage']:
                    hit_start = hit[0]
                    hit_end = hit[1]
                    writeToGFF(hits_gff_file, contig_name, tax_species, hit_start, hit_end, tax_species, hit_count)
                    hit_count += 1
            # Write chimeric regions
            # List of taxid1, region end, taxid2, region start
            chimeric_count = 0
            for chimeric_region in chimeric_region_list:
                first_species_id = ncbi.translate_to_names([chimeric_region[0]])[0]
                region_start = int(chimeric_region[1])
                second_species_id = ncbi.translate_to_names([chimeric_region[2]])[0]
                region_end = int(chimeric_region[3])
                region_length = region_end - region_start
                # Goal is to have 300bp regions
                if region_length < LEN_BETWEEN_TAXA:
                    adjustment = 150 - int(region_length/2)
                    region_start -= adjustment
                    region_end += adjustment
                    header = '>%s_species_%s_%s_position_%s_%s' % (
                    contig_name, first_species_id, second_species_id, region_start, region_end)
                    gff_feature_id = '%s_species_%s_%s_position_%s_%s' % (
                    contig_name, first_species_id, second_species_id, region_start, region_end)
                    region_sequence = contig.sequence[region_start:region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', region_start, region_end, gff_feature_id,
                               chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_species_id, second_species_id, region_start, region_end))
                    chimeric_count += 1
                  # Divide into two regions, one at each junction
                else:
                    adjustment = int(LEN_BETWEEN_TAXA / 2)
                    first_region_start = region_start - adjustment
                    first_region_end = region_start + adjustment
                    second_region_start = region_end - adjustment
                    second_region_end = region_end + adjustment

                    header = '>%s_species_%s_%s_position_%s_%s' % (
                        contig_name, first_species_id, second_species_id, first_region_start, first_region_end)
                    gff_feature_id = '%s_species_%s_%s_position_%s_%s' % (
                        contig_name, first_species_id, second_species_id, first_region_start, first_region_end)
                    region_sequence = contig.sequence[first_region_start:first_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', first_region_start, first_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_species_id, second_species_id, region_start, region_end))
                    chimeric_count += 1
                    header = '>%s_species_%s_%s_position_%s_%s' % (
                        contig_name, first_species_id, second_species_id, second_region_start, second_region_end)
                    gff_feature_id = '%s_species_%s_%s_position_%s_%s' % (
                        contig_name, first_species_id, second_species_id, second_region_start, second_region_end)
                    region_sequence = contig.sequence[second_region_start:second_region_end]
                    chimeric_fasta_file.write('%s\n%s\n' % (header, region_sequence))
                    writeToGFF(chimeric_gff_file, contig_name, 'chimera', second_region_start, second_region_end,
                               gff_feature_id, chimeric_count)
                    chimeric_regions_table_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                        chimeric_count, contig_name, first_species_id, second_species_id, region_start, region_end))
                    chimeric_count += 1

    else:
        print("Desired Rank Doesn't Exist")
    
    print("Done identifying chimeric regions.")
    # === Remove Chimeric Regions that do not span the full 300 bp's  ===
    #with open(chimeric_regions_table_file, newline='') as infile, tempfile.NamedTemporaryFile('w', delete=False, newline='') as tmpfile:
	#reader = csv.reader(infile, delimiter='\t')
	#writer = csv.writer(tmpfile, delimiter='\t')
	#header = next(reader)
	#writer.writerow(header)
	#for row in reader:
	#	has_negative = any(cell.lstrip('-').replace('.', '', 1).isdigit() and float(cell) < 0 for cell in row)
	#	if not has_negative:
	#		writer.writerow(row)
    #shutil.move(tmpfile.name, input_file)
    #chimeric_regions_table = path.join(out_dir, f"{prefix}_chimeric_regions.tsv")
    #chim_result, chim_count  = has_data(chimeric_regions_table)
    # Files from GFA
    #vf_file_contigs = path.join(proj_dir,"AssemblyBasedAnalysis","SpecialtyGenes",f"{prefix}_VF_genes_ShortBRED_krona_list.txt")
    #ar_file_contigs = path.join(proj_dir,"AssemblyBasedAnalysis","SpecialtyGenes",f"{prefix}_AR_genes_rgi.txt")
    #vf_file_reads   = path.join(proj_dir,"ReadsBasedAnalysis","SpecialtyGenes",f"{prefix}_VF_genes_ShortBRED_krona_list.txt")
    #ar_file_reads   = path.join(proj_dir,"ReadsBasedAnalysis","SpecialtyGenes",f"{prefix}_AR_genes_ShortBRED_table.txt")
    # Output File
    #summary_path = path.join(out_dir, "chimeric_summary.txt")
    #write_summary(chim_result, chim_count, prefix, vf_file_contigs, ar_file_contigs, vf_file_reads, ar_file_reads, output_path=summary_path)
    #write_summary(result, count, prefix,summary_path)
    #print(f"[INFO] Summary written to {summary_path}")
def run():
    parser = argparse.ArgumentParser(description="Process contigs to find chimeric regions.")
    parser.add_argument('-c', '--contigs', required=True, dest='contigs', help='Input Contigs to be processed for '
                                                                               'chimeric regions')
    parser.add_argument('-p', '--prefix', required=True, dest='prefix', help='Prefix for output files and headers.')
    parser.add_argument('-o', '--outDir', default=os.path.curdir, dest='outDir', help='Directory to write output to, '
                                                                                      'default is current dirctory.')
    parser.add_argument('--blast_db1', required=True, dest='blast_db1', help='Full path to BLAST database to be used '
                                                                             'for first BLAST')
    parser.add_argument('-t', '--threads', default=4, dest='num_cpu', help='Number of cpu threads to use for BLAST.')    
    parser.add_argument('--desired_rank',default='order', help='Desired Ranks for taxonomic classification.')
    parser.add_argument('-m', '--max_target_seqs', default=5, dest='max_target_seqs', help='Number of max target sequences from BLAST.')
    parser.add_argument('--project_dir',default="/home/", help='Project Directory')
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
