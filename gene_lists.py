import streamlit as st
import pandas as pd
import re
import requests
from bioservices import KEGG

# ------------------ HELPER FUNCTIONS ------------------ #

# Function for matching genes from text areas
def find_matching_genes(list1, list2):
    set1 = set(gene.lower() for gene in list1.split())
    set2 = set(gene.lower() for gene in list2.split())
    matches = set1.intersection(set2)
    return matches

def guess_gene_column(df, reference_sets):
    """
    Try to guess which column contains gene symbols by:
    - Checking for gene-like strings
    - Maximizing overlap with known reference gene sets
    """
    # unified reference set of known genes (uppercased)
    ref_genes = set()
    for s in reference_sets:
        ref_genes |= {str(g).upper() for g in s}

    best_col = None
    best_score = -1

    # simple pattern: letters/numbers/-/_/. and at least one letter
    pattern = re.compile(r"^[A-Za-z0-9_.-]+$")

    for col in df.columns:
        series = df[col].astype(str).str.strip()
        if series.empty:
            continue

        # Filter to entries that look like plausible gene symbols
        mask = series.str.match(pattern) & series.str.contains("[A-Za-z]")
        if mask.sum() == 0:
            continue

        # Uppercase candidate values
        vals = {v.upper() for v in series[mask] if v and v.upper() != "NAN"}
        overlap = len(vals & ref_genes)

        if overlap > best_score:
            best_score = overlap
            best_col = col

    return best_col

def normalize_key(x):
    """Normalize a value to an uppercase gene key."""
    return str(x).strip().upper()

# ------------------ KEGG FUNCTION ------------------ #

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
            st.image(png_url, caption=f"KEGG Pathway Map: {pathway_id}", use_column_width=True)

        except Exception as e:
            st.error(f"Error fetching or parsing KEGG data: {e}")

# ------------------ MAIN APP ------------------ #

def main():
    st.title('Gene Matcher and Gene List Analysis App')

    # Sidebar menu for selecting functionality
    # File-based Gene Matcher is now the DEFAULT (index=0)
    menu = ["File-based Gene Matcher", "Text Area Gene Matcher", "KEGG Gene List"]
    choice = st.sidebar.selectbox("Choose the functionality", menu, index=0)

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

        # Paths to master lists
        saids_path = 'SAIDs_2025.csv'                # columns: Disease, Gene, Inheritance
        immune_genes_path = 'Immune_genes_2025.csv'  # columns: Disease, Gene, Inheritance

        try:
            # ----- SAIDs_2025 -----
            saids_df = pd.read_csv(saids_path)

            saids_df['Gene_upper'] = (
                saids_df['Gene']
                .astype(str)
                .str.strip()
                .str.upper()
            )

            # gene -> SAID disease
            said_gene_to_disease = (
                saids_df
                .groupby('Gene_upper')['Disease']
                .apply(lambda x: '; '.join(sorted(set(map(str, x)))))
                .to_dict()
            )

            # gene -> SAID inheritance
            said_gene_to_inheritance = (
                saids_df
                .groupby('Gene_upper')['Inheritance']
                .apply(lambda x: '; '.join(sorted(set(map(str, x)))))
                .to_dict()
            )

            said_gene_set = set(said_gene_to_disease.keys())

            # ----- Immune_genes_2025 -----
            immune_df = pd.read_csv(immune_genes_path)

            immune_df['Gene_upper'] = (
                immune_df['Gene']
                .astype(str)
                .str.strip()
                .str.upper()
            )

            # gene -> Immune disease
            immune_gene_to_disease = (
                immune_df
                .groupby('Gene_upper')['Disease']
                .apply(lambda x: '; '.join(sorted(set(map(str, x)))))
                .to_dict()
            )

            # gene -> Immune inheritance
            immune_gene_to_inheritance = (
                immune_df
                .groupby('Gene_upper')['Inheritance']
                .apply(lambda x: '; '.join(sorted(set(map(str, x)))))
                .to_dict()
            )

            immune_gene_set = set(immune_gene_to_disease.keys())

        except Exception as e:
            st.error(f"Failed to load master gene lists: {e}")
            return

        # For auto-detection of gene column: union of both sets
        reference_sets = [said_gene_set, immune_gene_set]

        # File uploader for the large CSV, TSV, or TXT file
        large_file = st.file_uploader(
            "Upload your large file (CSV, TSV, or TXT)",
            type=['csv', 'tsv', 'txt']
        )

        if large_file is not None:
            file_extension = large_file.name.split('.')[-1]
            delimiter = ',' if file_extension.lower() == 'csv' else '\t'

            try:
                large_df = pd.read_csv(large_file, delimiter=delimiter)

                # --- Auto-detect gene column ---
                gene_column = guess_gene_column(large_df, reference_sets)

                # Fallback to common names if needed
                if gene_column is None:
                    if 'gene' in large_df.columns:
                        gene_column = 'gene'
                    elif 'ann.gene' in large_df.columns:
                        gene_column = 'ann.gene'
                    elif 'SYMBOL' in large_df.columns:
                        gene_column = 'SYMBOL'

                if gene_column is None:
                    st.error(
                        "Could not automatically detect a gene column. "
                        "Please ensure there is a column with gene symbols."
                    )
                    st.write("Available columns:", list(large_df.columns))
                    return
                else:
                    st.info(f"Using column **'{gene_column}'** as gene column.")

                # Normalize gene column
                large_df[gene_column] = large_df[gene_column].astype(str).str.strip()

                # ---- Annotate uploaded file using SAIDs_2025 and Immune_genes_2025 ----

                def said_gene_symbol(x):
                    key = normalize_key(x)
                    return key if key in said_gene_to_disease else ''

                def immune_gene_symbol(x):
                    key = normalize_key(x)
                    return key if key in immune_gene_to_disease else ''

                large_df['SAID gene'] = large_df[gene_column].apply(said_gene_symbol)
                large_df['SAID disease'] = large_df[gene_column].apply(
                    lambda x: said_gene_to_disease.get(normalize_key(x), '')
                )
                large_df['SAID disease inheritance'] = large_df[gene_column].apply(
                    lambda x: said_gene_to_inheritance.get(normalize_key(x), '')
                )

                large_df['Immune gene'] = large_df[gene_column].apply(immune_gene_symbol)
                large_df['Immune disease'] = large_df[gene_column].apply(
                    lambda x: immune_gene_to_disease.get(normalize_key(x), '')
                )
                large_df['Immune disease inheritance'] = large_df[gene_column].apply(
                    lambda x: immune_gene_to_inheritance.get(normalize_key(x), '')
                )

                # Display the DataFrame
                st.write(large_df)

                # Convert DataFrame to CSV and let user download it
                csv = large_df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    "Download modified file",
                    csv,
                    "modified_large_file.csv",
                    "text/csv",
                    key='download-csv'
                )

            except Exception as e:
                st.error(f"Failed to read or process the uploaded file: {e}")

    elif choice == "KEGG Gene List":
        kegg_gene_list()

if __name__ == "__main__":
    main()