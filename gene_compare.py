import streamlit as st

def find_matching_genes(list1, list2):
    # Split the input strings into lists, convert them to lowercase, and convert them to sets for comparison
    set1 = set(gene.lower() for gene in list1.split())
    set2 = set(gene.lower() for gene in list2.split())
    # Find the intersection of both sets
    matches = set1.intersection(set2)
    return matches

def main():
    st.title('Gene Matcher App')

    # Create text areas for input
    genes_list_1 = st.text_area("Paste your first list of genes here:", height=200)
    genes_list_2 = st.text_area("Paste your second list of genes here:", height=200)

    # Button to find matches
    if st.button('Find Matches'):
        if genes_list_1 and genes_list_2:
            matching_genes = find_matching_genes(genes_list_1, genes_list_2)
            if matching_genes:
                # Convert matching genes to uppercase
                matching_genes_upper = [gene.upper() for gene in matching_genes]
                st.success(f"Matching Genes: {', '.join(matching_genes_upper)}")
            else:
                st.warning("No matching genes found.")
        else:
            st.error("Please paste both lists of genes to find matches.")

if __name__ == "__main__":
    main()
