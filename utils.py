def query_ptm(params):
    with connection.cursor() as cursor:
        # Query to get all ptms
        cursor.execute("""
        SELECT
            CASE
                WHEN gene_names.uniprot_id = %(user_protein)s THEN %(user_protein)s
                WHEN gene_names.gene_name = %(user_protein)s THEN gene_names.uniprot_id
            END AS uniprot_id,
            gene_names.gene_name,
            alternative_gene_names.alternative_gene_name,
            protein_names.protein_name,
            ptm.ptm_type,
            ptm.start AS ptm_start,
            ptm.end AS ptm_end,
            ptm.residue AS wt_residue,
            CONCAT(ptm_source.source, ": ", ptm_source.id) AS ptm_source,
            CASE
                WHEN ptm.residue = %(user_wt)s THEN 'True'
                ELSE 'False'
            END AS wt_check
        FROM
            gene_names
        LEFT JOIN 
            alternative_gene_names ON gene_names.uniprot_id = alternative_gene_names.uniprot_id
        LEFT JOIN 
            protein_names ON gene_names.uniprot_id = protein_names.uniprot_id
        LEFT JOIN
            ptm ON protein_names.uniprot_id = ptm.uniprot_id
        LEFT JOIN 
            ptm_source ON ptm.uniprot_id = ptm_source.uniprot_id
            AND ptm.ptm_type = ptm_source.ptm_type
            AND ptm.start = ptm_source.start
            AND ptm.end = ptm_source.end
            AND ptm.residue = ptm_source.residue
        WHERE
            gene_names.uniprot_id = %(user_protein)s OR gene_names.gene_name = %(user_protein)s;
        """, params)

        # get results from query and return as a DataFrame
        rows = cursor.fetchall()
        # checks if valid Uniprot ID or gene name exits in database
        if len(rows) == 0:
            return(None)

        ptm_df = pd.DataFrame(rows, columns = [col[0] for col in cursor.description])
        # ptm_df.sort_values(by = "ptm_start", inplace=True)

        # Condense the DataFrame such that:
        #   alternative gene names are in a list 
        #   ptm_info is a list of dictionaries with {ptm_type, ptm_start, ptm_end, wt_residue, ptm_source, wt_check}
        grouped = ptm_df.groupby(['uniprot_id', 'gene_name', 'protein_name'])

        # Initialize an empty dictionary to store the condensed data
        condensed_dict = {}

        # Iterate through each group
        for group_name, group_data in grouped:
            # extract unique identifiers
            uniprot_id, gene_name, protein_name = group_name
            
            ptm_info_set = set() # create a set to store unique ptm_info dictionaries
            alternative_gene_names = set()  # to store unique alternative gene names
            
            for _, row in group_data.iterrows():
                # Create a hashable version of ptm_info dictionary
                ptm_info_hashable = tuple(row[['ptm_type', 'ptm_start', 'ptm_end', 'wt_residue', 'ptm_source', 'wt_check']])
                
                # Add ptm_info to set if it's not already present
                if ptm_info_hashable not in ptm_info_set:
                    ptm_info_set.add(ptm_info_hashable)
                
                # Collect unique alternative gene names
                alternative_gene_names.add(row['alternative_gene_name'])
            
            # Convert set of ptm_info back to list of dictionaries
            ptm_info_list = [dict(zip(['ptm_type', 'ptm_start', 'ptm_end', 'wt_residue', 'ptm_source', 'wt_check'], item)) for item in ptm_info_set]
            
            # Create the final dictionary entry
            entry = {
                'uniprot_id': uniprot_id,
                'gene_name': gene_name,
                'alternative_gene_names': list(alternative_gene_names),
                'protein_name': protein_name,
                'ptm_info': ptm_info_list
            }
            
            # Use uniprot_id as the dictionary key
            condensed_dict[uniprot_id] = entry

    return entry
