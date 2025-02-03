#!/usr/bin/env python

''' This script takes in the full list of human genetic variants downloaded from the Uniprot FTP site (/project/data/VariantP/variant_data/homo_sapiens_variation.txt)
and filters it according to a list of canonical proteins (/project/data/VariantP/variant_data/canonical_filtered.txt) for all missense variants.
'''

import Bio
import re
import numpy as np
import pandas as pd
from Bio.SeqUtils import seq1
import os
from datetime import datetime

# Function to get variant file of current protein
def get_variant_file(uniprot_id, output_path):
    """
    Uses the current uniprot_id to retrieve the file containing all human variants
    mapped to the protein. 
    """
    # File path to all variant files
    variants_dir = '/project/data/VariantP/split_data'
    output_path = '/project/data/VariantP/variant_data/variant_table_files'


    # Current protein
    protein = uniprot_id.rstrip()

    # Get first letter of accession code
    letter_dir = protein[0]

    # Get variant file using first letter
    try:
        variant_file = open(f'{variants_dir}/{letter_dir}/{protein}.txt', 'r')
        print('Current Protein:', protein)
        return variant_file

    except FileNotFoundError:
        print('Protein ', protein, 'has no reported variants. Moving to next protein.')

        # Catch any additional canonical proteins that don't contain any variants in a .txt file
        with open(f'{output_path}/no_variants.txt', 'a') as none_file:
            none_file.write(protein + '\n')

        return FileNotFoundError

# Function to extract phenotype information
def get_phenotype_columns(variant):
    """
    Filters extracted data of current variant for relevant phenotypic information.

    :param str, variant: string of variant information separated by commas from original homo_sapiens_variation.txt file.
    """
    # Initiate list to store data
    phenotype_data = []

    # Get snp id
    snp = variant[3]
    phenotype = variant[5:8]
    phenotype_data.append(snp)
    for item in phenotype:
        phenotype_data.append(item)

    # Returns phenotype data columns
    return phenotype_data

# Function to extract general variant information
def get_variant_columns(variant):
    """
    Filters extracted data of current variant for variant table information.

    :param str, variant: string of variant information separated by commas from original homo_sapiens_variation.txt file.
    """

    # Get snp id
    snp = variant[3]

    # Regular expression to get wt, position, and mutant variant information
    pattern = re.compile(r'([a-zA-Z\W]*)(\d+)([a-zA-Z\W]*)')

    # Extract groups
    mutant_data = variant[2].replace('p.', '')
    match = re.findall(pattern, mutant_data)

    # Get variant information
    wt = match[0][0]
    position = match [0][1]
    mutant = match [0][2]

    # Translate three 3 code to 1 letter code
    if len(wt) == 3:
        wt = seq1(wt)
    if len(mutant) == 3:
        mutant = seq1(mutant)

    # Format columns for table
    variant_data = [variant[1], wt, position, mutant, snp, variant[12]]
    
    # Returns filtered variant table columns
    return variant_data

# Function to extract database information
def get_database_columns(variant):
    """
    Filters extracted data of current variant for variant source (database) information.
    Sorts any unique databases into a new column.

    :param str, variant: string of variant information separated by commas from original homo_sapiens_variation.txt file.
    """
    # Initialize list for storing databases
    database_data = []

    # Get databases
    databases = variant[len(variant)-1].split(',')

    # Append the db
    for db in databases:
        database_data.append(db)
        
    return database_data

# Function to filter phenotype data
def filter_phenotype_data(phenotype_cols):
    """
    Cleans up phenotype columns and sorts source information into individual columns for OMIM id,
    ECO id, and PubMed id(s).

    :param list, phenotype_cols: list of phenotypic information associated to a single variant extracted from homo_sapiens_variantion.txt file.
    """
    # Variable to store phenotype sources
    source_cols = []


    # Check if variant has PubMed ids with no tag associated
    match2 = re.findall(r'\b\d{8}\b', phenotype_cols[3])
    if match2:
        for pubmed_id in match2:
            source_cols.append(pubmed_id)

    # Check if variant has a MIM id either in phenotype column ...
    match4 = re.search(r'MIM:(\d+)', phenotype_cols[2])
    match5 = re.search(r'MIM:(\d+)', phenotype_cols[3])

    if match4:
        # Extract only MIM id and update column in phenotype table
        phenotype_cols[3] = match4.group(1)
        # Update phenotype col to only contain phenotype info
        phenotype_cols[2] = phenotype_cols[2][:-len(match4.group())-3]
    # ... Or in source column
    elif match5:
        phenotype_cols[3] = match5.group(1)
    else:
        phenotype_cols[3] = '-'

    # Check if source is only SNP id
    match6 = re.fullmatch(r'(RCV)\d{9}', phenotype_cols[3])
    if match6:
        # Update as space holder
        phenotype_cols[3] = '-'
    else:
        pass

    return phenotype_cols, source_cols

# Function to get all filtered data
def variant_filter(proteins, output_path):
    """
    Uses list of all canonical proteins (~20,000) to filter out information for all missense variants on all
    proteins.

    :param str, proteins: List of all proteins in current batch to be processed. If no increment was specified, list contains all proteins.
    :param str, output_path: Directory path for output files containing extracted and filtered variant data
    """
    # Initialize counters as check
    var_counter = 0

    # Initialize dict of lists to organize each table's information for all missense variants
    variant_table = []
    pheno_table = []
    db_table = []
    source_table = []
    all_data = {'Variants': variant_table, 'Phenotype': pheno_table, 'Databases': db_table, 'Sources': source_table}


    # To keep track of non-repeating variant ids
    ids = []

    for uniprot_id in proteins:
        # Get variant file
        all_variants = get_variant_file(uniprot_id, output_path) ### This is when current protein check is printed

        if all_variants == FileNotFoundError:
            continue
        else:
            # Extract variants by line
            for line in all_variants:
                # Variable to store current variant information
                line = line.rstrip()
                variant = line.split('\t')

                # Get snp id
                var_id = variant[3]

                # Check if variant is missense variant
                if variant[4] == 'missense variant' or variant[4] == 'Missense':
                    
                    # # Extract columns of data for each table (variant, phenotype, database)
                    variant_cols = get_variant_columns(variant)
                    phenotype_cols = get_phenotype_columns(variant)
                    database_cols = get_database_columns(variant)

                    # Update variant table list
                    variant_table.append(variant_cols)

                    # Update variant count
                    var_counter += 1

                    # If variant has phenotype data
                    if any(data != '-' for data in phenotype_cols[1:4]):
                        # Cleaning up the phenotype data
                        filtered_pcols, source_cols = filter_phenotype_data(phenotype_cols)


                        # Update phenotype table, skipping repeats
                        if any(data != '-' for data in filtered_pcols[1:3]):
                            pheno_table.append(filtered_pcols)

                        # Update source table
                        for source in source_cols:
                            new_row = [var_id, source]
                            source_table.append(new_row)


                    # If repeat variant
                    if var_id in ids:
                            # Check for database redundancy
                            for db in database_cols:
                                # Checking the last ten variants for repetition
                                # Append this db to the variants pre-existing row of db's
                                db_table.append([var_id, db])

                    else:
                        # If novel variant, append db information as new rows to table list
                        for db in database_cols:
                            new_row = [var_id, db]
                            db_table.append(new_row)


                    # Update ids list with new variant id

                    ids.append(var_id)

                else:
                    pass
    
        
    # Update dictionary with all three table lists
    all_data['Variants'] = variant_table
    all_data['Phenotype'] = pheno_table
    all_data['Databases'] = db_table
    all_data['Sources'] = source_table


    return all_data

# Function to print variant data to outfiles
def output_variant_files(output_path, curr_var, curr_pheno, curr_db, curr_sources, append=True):

    """
    Generates outfiles for variant information respective to each of the three tables in the database.

    :param str, output_path: Desired file path for output text files.
    :param list, curr_var: Nested list of variant table data, with a list per table row, for current protein batch.
    :param list, curr_pheno: Nested list of phenotype table data for current protein batch.
    :param list, curr_db: Nested list of database table data for current protein batch.
    :param bool, append: Boolean parameter to prevent overwriting of files; Default is True, but if first batch append=False and empty files will be generated with 'w'.
    """
    mode = 'a' if append else 'w'
    # Write missense variant information for current protein to outfile
    with open(f'{output_path}/variant_table.txt', mode) as var_outfile:
        
        header = 'uniprot_id\twt\tposition\tmutant\tsnp_id\tensembl_id\n'
        if mode == 'w':
            var_outfile.write(header)

        # Append variant information
        for var in curr_var:
            row = '\t'.join(var)
            var_outfile.write(row + '\n')

    # Write variant phenotype information to outfile
    with open(f'{output_path}/phenotype_table.txt', mode) as pheno_outfile:

        # Generate column names
        header = ['snp_id', 'clinical_sig', 'phenotype', 'omim_id']

        header = '\t'.join(header) + '\n'
        if mode == 'w':
            pheno_outfile.write(header)

        # Append phenotype information
        for var in curr_pheno:
            row = '\t'.join(var)
            pheno_outfile.write(row + '\n')

    # Write variant phenotype information to outfile
    with open(f'{output_path}/sources_table.txt', mode) as pheno_source_outfile:

        # Generate column names
        header = ['snp_id', 'pubmed_id']

        header = '\t'.join(header) + '\n'
        if mode == 'w':
            pheno_source_outfile.write(header)

        # Append phenotype information
        for pubmed in curr_sources:
            row = '\t'.join(pubmed)
            pheno_source_outfile.write(row + '\n')

    # Write variant source information to outfile
    with open(f'{output_path}/db_table.txt', mode) as db_outfile:
        
        header = 'snp_id\tdatabase\n'
        if mode == 'w':
            db_outfile.write(header)

        # Append database information
        for var in curr_db:
            row = '\t'.join(var)
            db_outfile.write(row + '\n')

# Function to get count of missense variants extracted
def get_variant_count(variant_table):
    """
    Counts number of unique missense variants extracted for database population.

    :param list, variant_table: Nested list containing all rows of variants extracted for current batch of proteins.
    """
    variant_df = pd.DataFrame(variant_table)
    # Get unique variants
    unique = variant_df.iloc[:, 0:4].drop_duplicates()
    var_count = unique.shape[0]
    # Return number of rows in df
    return var_count


# Function to get all filtered variant data
def get_variant_data(canonical_path, output_path, increment=None, print_files=False):
    """
    Returns a dictionary with 3 key/value pairs. Values are nested lists of filtered missense variant data of each table with associated
     keys to identify the table ('Variants', 'Phenotype', 'Databases').

    :param str, canonical_path: File path to list of canonical proteins
    :param str, output_path: Desired file path for output text files
    :param int, increment: Defines the number of proteins to be scanned for variants at a time. If not specified,
     default is None and all proteins will be scanned for variants in a single batch.
    :param bool, print_files: Specifies if variant information should be printed to .txt outfiles. If not specified, default is False and outfiles will not be produced.

    """
    # Initialize run time, variant counter, and dictionary to store extracted variant data
    startTime = datetime.now()
    var_count = 0
    total_data = {'Variants': [], 'Phenotype': [], 'Databases': [], 'Sources': []}

    # Clear file for proteins with no variants if it exists already
    try:
        os.remove(f'{output_path}/no_variants.txt')
    except FileNotFoundError:
        pass

    # For each canonical protein, get variant file
    with open(f'{canonical_path}', 'r') as in_file:
        # Get all canonical proteins
        all_proteins = in_file.readlines()

        # If entire protein list being scanned
        if increment == None:
            total_data = variant_filter(all_proteins, output_path)
            var_count = get_variant_count(all_proteins['Variants'])
            
            # Check number of proteins and total missense variants found
            print('Proteins processed:', len(all_proteins), 'Variant table rows:', len(total_data['Variants']), 'Unique variants recorded:', var_count)


        # If proteins being scanned in increments
        elif increment is not None:
            # Split proteins into batches
            num_batches = len(all_proteins) // increment 

            # Print the number of total proteins
            print('Number of proteins:', len(all_proteins))

            # Get batch start and end
            for batch in range(num_batches+1):
                batch_start = batch * increment
                batch_end = (batch+1) * increment
                processed = batch_end

                # If last batch
                if batch == num_batches:
                    # Batch processes up to last protein
                    batch_end = (len(all_proteins) + 1)
                    processed = len(all_proteins)
                
                # Get current batch of proteins
                curr_batch = all_proteins[batch_start:batch_end]

                # Filter variants for all proteins in batch (CALLING VARIANT FILTER)
                curr_batch_variants = variant_filter(curr_batch, output_path)

                # Check counts
                var_count += get_variant_count(curr_batch_variants['Variants'])
                print('Proteins processed:', processed, 'Unique variants recorded:', var_count)

                # If filtered variant data is being exported to outfiles
                if print_files==True:
                    # Get variant data for current batch
                    curr_var = curr_batch_variants['Variants']
                    curr_pheno = curr_batch_variants['Phenotype']
                    curr_db = curr_batch_variants['Databases']
                    curr_sources = curr_batch_variants['Sources']
                    
                    # If this is the first batch of proteins, no need to append
                    if batch_start == 0:
                        output_variant_files(output_path, curr_var, curr_pheno, curr_db, curr_sources, append=False)
                    # If any other batch, append to current extracted data
                    else:
                        output_variant_files(output_path, curr_var, curr_pheno, curr_db, curr_sources, append=True)
                
                # Append current batch dictionary to total dictionary
                for table_type in curr_batch_variants.keys():
                    total_data[table_type].extend(curr_batch_variants[table_type])
    
    # Print time script took to run
    print('Time of processing: ', datetime.now() - startTime)


def main():

    # Define file paths for canonical proteins and variant files
    # Note: Make sure output_path is set to desired path for variant text files
    output_path = '/project/data/VariantP/variant_data/variant_table_files'
    canonical_test = '/project/data/VariantP/variant_data/canonical_test20.txt'
    canonical_filtered = '/project/data/VariantP/variant_data/canonical_filtered.txt'

    # Run filter for desired proteins at desired batch size
    all_data = get_variant_data(canonical_filtered, output_path, increment=500, print_files=True)

    return all_data

if __name__ == '__main__':
    main()


