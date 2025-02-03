from django.shortcuts import render
from .forms import QueryForm
from django.http import HttpResponse
from django.db import connection

# Create your views here.
# def index(request):
#     return HttpResponse("Hello, world. You're at the dbQuery index.")

def process_form(request):
    if request.method == "POST":
        form = QueryForm(request.POST)
        if form.is_valid():
            # print("Form is valid")
            uniprot_id = form.cleaned_data["uniprot_id"]
            gene_name = form.cleaned_data["gene_name"]
            position = form.cleaned_data['position']
            residue = form.cleaned_data["residue"]
            if uniprot_id:
                input_data = {'user_protein':uniprot_id, 'user_position': position, 'user_wt': residue}
            elif gene_name:
                input_data = {'user_protein':gene_name, 'user_position': position, 'user_wt': residue}


            # Empty dict to store data
            returned_data = {}

            with connection.cursor() as cursor:
                    # Query to get topology and pdb file path
                    cursor.execute("""
                    SELECT
                        CASE
                            WHEN gene_names.uniprot_id = %(user_protein)s THEN %(user_protein)s
                            WHEN gene_names.gene_name = %(user_protein)s THEN gene_names.uniprot_id
                        END AS uniprot_id,
                        gene_names.gene_name,
                        topology.region_type,
                        topology.region_number,
                        topology.start,
                        topology.end,
                        structure.structure AS structure_path
                    FROM
                        gene_names
                    LEFT JOIN
                        topology ON gene_names.uniprot_id = topology.uniprot_id
                    LEFT JOIN
                        structure ON gene_names.uniprot_id = structure.uniprot_id
                    WHERE
                        gene_names.uniprot_id = %(user_protein)s OR gene_names.gene_name = %(user_protein)s;

                    """, input_data)

                    topology = cursor.fetchall()
                    #topology_columns = [column[0] for column in cursor.description]


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
                        ptm_source.source AS ptm_source,
                        ptm_source.id AS ptm_id,
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
                    """, input_data)

                    ptms = cursor.fetchall()
                    ptms_df = pd.DataFrame(ptms, col_names = [col[0] for col in cursor.description])

                    # Condense df to put alt gene names are in list
                    ptms_df_grouped = ptm_df.groupby(['uniprot_id', 'gene_name', 'protein_name'])

                    # Gets all additional ptms
                    cursor.execute("""
                    SELECT
                        CASE
                            WHEN gene_names.uniprot_id = %(user_protein)s THEN %(user_protein)s
                            WHEN gene_names.gene_name = %(user_protein)s THEN gene_names.uniprot_id
                        END AS uniprot_id,
                        gene_names.gene_name,
                        protein_names.protein_name,
                        additional_ptm.description AS addtl_ptm,
                        CONCAT (additional_ptm_sources.source, '-', additional_ptm_sources.id) AS addtl_ptm_source
                    FROM
                        gene_names
                    LEFT JOIN 
                        protein_names ON gene_names.uniprot_id = protein_names.uniprot_id
                    LEFT JOIN
                        additional_ptm on protein_names.uniprot_id = additional_ptm.uniprot_id
                    LEFT JOIN
                        additional_ptm_sources on additional_ptm.uniprot_id = additional_ptm_sources.uniprot_id
                    WHERE
                        gene_names.uniprot_id = %(user_protein)s OR gene_names.gene_name = %(user_protein)s;
                    """, input_data)

                    # Get results
                    results = cursor.fetchall()
                    # Get column names
                    col_names = [column[0] for column in cursor.description]

                    # Post processing data
                    resultsdf = pd.DataFrame(results, columns=col_names)
                    addtl_ptms = resultsdf.groupby(['addtl_ptm'])['addtl_ptm_source'].agg(list).reset_index()



                    # Query to get all variants at user input position
                    cursor.execute("""
                    SELECT
                        CASE
                            WHEN gene_names.uniprot_id = %(user_protein)s THEN %(user_protein)s
                            WHEN gene_names.gene_name = %(user_protein)s THEN gene_names.uniprot_id
                        END AS uniprot_id,
                        variants.snp_id, variants.position, variants.wt, variants.mutant, variants.ensembl_id, variant_databases.database_name AS db,
                        phenotype.clinical_significance, phenotype.phenotype AS disease, phenotype.omim_id, phenotype_sources.pubmed_id AS pheno_source,
                        topology.region_type, topology.region_number, topology.start, topology.end,
                        CASE
                            WHEN variants.wt = %(user_wt)s THEN 'True'
                            ELSE 'False'
                        END AS wt_check
                    FROM
                        gene_names
                    LEFT JOIN
                        variants ON gene_names.uniprot_id = variants.uniprot_id
                        AND (%(user_position)s = variants.position)
                    LEFT JOIN
                        variant_databases ON variants.snp_id = variant_databases.snp_id
                    LEFT JOIN
                        phenotype ON variants.snp_id = phenotype.snp_id
                    LEFT JOIN
                        phenotype_sources ON variants.snp_id = phenotype_sources.snp_id
                    LEFT JOIN
                        topology ON variants.uniprot_id = topology.uniprot_id
                        AND (variants.position BETWEEN topology.start AND topology.end)
                    WHERE
                        gene_names.uniprot_id = %(user_protein)s OR gene_names.gene_name = %(user_protein)s;
                    """, input_data)

                    results = cursor.fetchall()
                    
                    # Convert to dataframe
                    col_names = [column[0] for column in cursor.description]
                    resultsdf = pd.DataFrame(results, columns=col_names)

                    # Get unique variant rows
                    variants = resultsdf.drop_duplicates(['snp_id', 'mutant', 'clinical_significance', 'disease', 'region_type', 'region_number', 'start', 'end', 'wt_check']).reset_index(drop=True)
                    
                    # Get lists for repeating columns
                    repeats = resultsdf.groupby(['uniprot_id', 'snp_id', 'mutant']).agg({
                        'ensembl_id': lambda x: list(set(x)),
                        'db': lambda x: list(set(x)),
                        'pheno_source': lambda x: list(set(x)),
                    }).reset_index()
                    
                    # Update to final formatting
                    variants.update(repeats)



                    # Collect all data to be returned
                    returned_data = {'topology':topology, 'ptms':ptms, 'addtl_ptms': addtl_ptms_df, 'variants':variants}

            return render(returned_data, "dbQuery/results.html", {"form_data": form.cleaned_data})
    else:
        form = QueryForm()

    return render(request, 'dbQuery/home.html', {'form': form})



def home(request):
    return render(request, 'dbQuery/home.html')

def documentation(request):
    return render(request, 'dbQuery/documentation.html')

def dataset(request):
    return render(request, 'dbQuery/dataset.html')

def contact(request):
    return render(request, 'dbQuery/contact.html')

def formFill(request):
    return HttpResponse("Function executed successfully")




