import streamlit as st
import pandas as pd

# Function for matching genes from text areas
def find_matching_genes(list1, list2):
    set1 = set(gene.lower() for gene in list1.split())
    set2 = set(gene.lower() for gene in list2.split())
    matches = set1.intersection(set2)
    return matches

# Main function for the Streamlit app
def main():
    st.title('Gene Matcher and Gene List Analysis App')

    # Sidebar menu for selecting functionality
    menu = ["Text Area Gene Matcher", "File-based Gene Matcher"]
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
        said_genes_path = '/Users/oskarschnappauf/Desktop/python/streamlit/gene_lists/SAIDs_book.csv'
        said_candidates_path = '/Users/oskarschnappauf/Desktop/python/streamlit/gene_lists/candidate_partner_genes.csv'
        innate_genes_path = '/Users/oskarschnappauf/Desktop/python/streamlit/gene_lists/innate_genes.csv'

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

if __name__ == "__main__":
    main()
