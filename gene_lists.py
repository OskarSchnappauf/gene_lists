import streamlit as st
import pandas as pd
import re
import requests
from bioservices import KEGG

# Function for matching genes from text areas
def find_matching_genes(list1, list2):
    set1 = set(gene.lower() for gene in list1.split())
    set2 = set(gene.lower() for gene in list2.split())
    matches = set1.intersection(set2)
    return matches

# Function for KEGG pathway analysis
def kegg_gene_list():
    st.header("KEGG Gene List Extractor")

    # User input for KEGG pathway ID
    user_input_id = st.text_input("Enter KEGG Pathway ID (e.g., map04064 or hsa04064)", "map04064")
    organism_prefix = st.text_input("Enter organism code (e.g., hsa for human, mmu for mouse)", "hsa")

    if st.button("Fetch KEGG Genes"):
        try:
            # --- Normalize Pathway ID ---
            if user_input_id.startswith("map"):
                pathway_id = user_input_id.replace("map", organism_prefix)
            elif user_input_id.startswith(organism_prefix):
                pathway_id = user_input_id
            else:
                st.error("Invalid KEGG pathway ID. Use 'hsa04810' or 'map04810'.")
                return

            kegg = KEGG()
            raw_kegg_data = kegg.get(pathway_id)
            parsed_kegg = kegg.parse(raw_kegg_data)

            # --- Get Pathway Name ---
            pathway_name_list = parsed_kegg.get("NAME", ["unknown_pathway"])
            pathway_name_raw = pathway_name_list[0] if isinstance(pathway_name_list, list) else pathway_name_list
            st.subheader(f"Pathway Name: {pathway_name_raw}")

            # --- Extract Genes ---
            gene_entries = parsed_kegg.get("GENE", {})
            genes = []

            if isinstance(gene_entries, dict):
                for gene_id, desc in gene_entries.items():
                    full_description = desc[0] if isinstance(desc, list) else desc
                    gene_symbol = full_description.split(";")[0].strip()
                    genes.append({
                        "GeneID": gene_id,
                        "GeneSymbol": gene_symbol,
                        "FullDescription": full_description
                    })

            gene_df = pd.DataFrame(genes)

            if not gene_df.empty:
                st.success(f"Extracted {len(gene_df)} genes.")
                st.dataframe(gene_df)

                # Show as comma-separated list
                gene_list_str = ",".join(gene_df["GeneSymbol"].dropna().astype(str))
                st.text_area("Comma-separated gene list", value=gene_list_str, height=100)

                # Download button for gene symbols
                csv = gene_df[["GeneSymbol"]].to_csv(index=False).encode("utf-8")
                st.download_button("Download Gene Symbols as CSV", csv, f"{pathway_id}_genes.csv", "text/csv")
            else:
                st.warning("No genes found in the pathway.")

            # --- Display KEGG Pathway Image ---
            png_url = f"https://www.kegg.jp/kegg/pathway/{organism_prefix}/{pathway_id}.png"
            st.image(png_url, caption=f"KEGG Pathway Map: {pathway_id}", use_container_width=True)

        except Exception as e:
            st.error(f"Error fetching or parsing KEGG data: {e}")

# Main function for the Streamlit app
def main():
    st.title('Gene Matcher and Gene List Analysis App')

    # Sidebar menu for selecting functionality
    menu = ["Text Area Gene Matcher", "File-based Gene Matcher", "KEGG Gene List"]
    choice = st.sidebar.selectbox("Choose the functionality", menu)

    if choice == "Text Area Gene Matcher":
        # Text Area Gene Matcher
        st.header('Text Area Gene Matcher')

        genes_list_1 = st.text_area("Paste your first list of genes here:", height=200)
        genes_list_2 = st.text_area("Paste your second list of genes here:", height=200)

        if st.button('Find Matches'):
            if genes_list_1 and genes_list_2:
                matching_genes = find_matching_genes(genes_list_1, genes_list_2)
                if matching_genes:
                    matching_genes_upper = [gene.upper() for gene in matching_genes]
                    st.success(f"Matching Genes: {', '.join(matching_genes_upper)}")
                else:
                    st.warning("No matching genes found.")
            else:
                st.error("Please paste both lists of genes to find matches.")

    elif choice == "File-based Gene Matcher":
        # File-based Gene Matcher
        st.header('File-based Gene Matcher')

        # Hardcoded paths for the gene list CSVs
        said_genes_path = 'SAIDs_book.csv'
        said_candidates_path = 'candidate_partner_genes.csv'
        innate_genes_path = 'innate_genes.csv'

        # Loading all gene lists from hardcoded CSVs
        try:
            said_genes_df = pd.read_csv(said_genes_path)
            said_genes_set = set(said_genes_df['gene'].dropna().unique())

            said_candidates_df = pd.read_csv(said_candidates_path)
            said_candidates_set = set(said_candidates_df['gene'].dropna().unique())

            innate_genes_df = pd.read_csv(innate_genes_path)
            innate_genes_set = set(innate_genes_df['gene'].dropna().unique())
        except Exception as e:
            st.error(f"Failed to load gene lists: {e}")
            return

        # File uploader for the large CSV, TSV, or TXT file
        large_file = st.file_uploader("Upload your large file (CSV, TSV, or TXT)", type=['csv', 'tsv', 'txt'])
        if large_file is not None:
            file_extension = large_file.name.split('.')[-1]
            delimiter = ',' if file_extension.lower() == 'csv' else '\t'

            try:
                large_df = pd.read_csv(large_file, delimiter=delimiter)

                # Check if 'gene' or 'ann.gene' column exists
                gene_column = None
                if 'gene' in large_df.columns:
                    gene_column = 'gene'
                elif 'ann.gene' in large_df.columns:
                    gene_column = 'ann.gene'
                elif 'SYMBOL' in large_df.columns:
                    gene_column = 'SYMBOL'
                else:
                    st.error("The uploaded file does not contain a 'gene' or 'ann.gene' column.")
                    return

                # Apply the matching logic for all gene lists
                large_df['SAID genes'] = large_df[gene_column].apply(lambda x: x if x in said_genes_set else '')
                large_df['SAID candidates'] = large_df[gene_column].apply(lambda x: x if x in said_candidates_set else '')
                large_df['Innate Genes'] = large_df[gene_column].apply(lambda x: x if x in innate_genes_set else '')

                # Display the DataFrame
                st.write(large_df)

                # Convert DataFrame to CSV and let user download it
                csv = large_df.to_csv(index=False).encode('utf-8')
                st.download_button("Download modified file", csv, "modified_large_file.csv", "text/csv", key='download-csv')
            except Exception as e:
                st.error(f"Failed to read the large file: {e}")

    elif choice == "KEGG Gene List":
        kegg_gene_list()

if __name__ == "__main__":
    main()