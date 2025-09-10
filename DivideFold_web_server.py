import streamlit as st
import streamlit.components.v1 as components
from st_forna_component import forna_component
import os
import numpy as np
import colorsys

# Session state initialization
for key, default in [("rna_sequence", ""), ("use_example_clicked", False)]:
    if key not in st.session_state:
        st.session_state[key] = default
st.set_page_config(layout="wide")

# Title
st.title("DivideFold web server")

EXAMPLE_SEQUENCE = "UUGAUGGAGAGUUUGAUCCUGGCUCAGGACGAACGCUGGCGGCGUGCUUAACACAUGCAAGUCGAGCGGGGAGCUUCGGCUCUCAGCGGCGAACGGGUGAGUAACACGUGGGCAACCUGCCCCGAGCUCUGGAAUAACCUCGGGAAACCGGGGCUAAUGCCGGAUAUGACAUUGCCGGGCAUCUGGUGGUGUGGAAAGAUUUAUCGGCUCGGGAUGGGCCCGCGGCCUAUCAGCUUGUUGGUGGGGUGAUGGCCUACCAAGGCGACGACGGGUAGCCGGCCUGAGAGGGCGACCGGCCACACUGGGACUGAGACACGGCCCAGACUCCUACGGGAGGCAGCAGUGGGGAAUAUUGCGCAAUGGGCGGAAGCCUGACGCAGCGACGCCGCGUGGGGGAUGACGGCCUUCGGGUUGUAAACCUCUUUCAGCAGGGACGAAGCGAGAGUGACGGUACCUGCAGAAGAAGCACCGGCCAACUACGUGCCAGCAGCCGCGGUAAUACGUAGGGUGCAAGCGUUGUCCGGAAUUAUUGGGCGUAAAGAGCUCGUAGGCGGCCUGUUGCGUCGGCUGUGAAAACCCGGGGCUCAACUCCGGGCCUGCAGUCGAUACGGGCAGGCUAGAGUCCGGCAGGGGAGACUGGAAUUCCUGGUGUAGCGGUGAAAUGCGCAGAUAUCAGGAGGAACACCGGUGGCGAAGGCGGGUCUCUGGGCCGGAACUGACGCUAAGGAGCGAAAGCGUGGGGAGCGAACAGGAUUAGAUACCCUGGUAGUCCACGCCGUAAACGUUGGGCGCUAGGUGUGGGGGACCUUCCACGGCCUCCGUGCCGCAGCUAACGCAUUAAGCGCCCCGCCUGGGGAGUACGGCCGCAAGGCUAAAACUCAAAGGAAUUGACGGGGGCCCGCACAAGCGGCGGAGCAUGUGGCUUAAUUCGAUGCAACGCGAAGAACCUUACCAGGGCUUGACAUGCAGGGAAAUCUCGUAGAGAUACGGGGUCCGUAAGGGUCCUGCACAGGUGGUGCAUGGCUGUCGUCAGCUCGUGUCGUGAGAUGUUGGGUUAAGUCCCGCAACGAGCGCAACCCUCGUCCUAUGUUGCCAGCGAGUCGUGUCGGGGACUCAUAGGAGACUGCCGGGGUCAACUCGGAGGAAGGUGGGGAUGACGUCAAGUCAUCAUGCCCCUUACGUCCUGGGCUGCACACAUGCUACAAUGGCCGGUACAAAGGGCUGCGAUGCCGUGAGGUGGAGCGAAUCCCAAAAAGCCGGUCUCAGUUCGGAUCGGGGUCUGCAACUCGACCCCGUGAAGUCGGAGUCGCUAGUAAUCGCAGAUCAGCAAUGCUGCGGUGAAUACGUUCCCGGGCCUUGUACACACCGCCCGUCACGUCACGAAAGUCGGUAACACCCGAAGCCGGUGGCCUAACCCUUGUGGGGGGAGCCGUCGAAGGUGGGACCGGCGAUUGGGACGAAGUCGUAACAAGGUAGCCGUACCGGAAGGUGCGGCUGGAUCACCUCCUUUCU"

# Get params
token = st.query_params.get("token", None)
if token is None:
    st.warning("No token provided. Please access this app through a proper job link.")
    st.stop()

with st.spinner("Loading libraries..."):
    from dividefold.predict import dividefold_predict, knotfold_predict, ipknot_predict, probknot_predict, pkiss_predict, rnafold_predict, linearfold_predict, mxfold2_predict


def main_page():
    # Set example sequence
    if st.session_state.get("use_example_clicked"):
        st.session_state.rna_sequence = EXAMPLE_SEQUENCE
        st.session_state.use_example_clicked = False

    # Sequence input
    sequence_input = st.text_area("RNA sequence:", value="", height=200, key="rna_sequence")

    # Button to set example sequence
    if st.button("Use example sequence"):
        st.session_state.use_example_clicked = True
        st.rerun()

    left, right = st.columns(2)

    # Maximum fragment length
    max_fragment_length = left.text_input("Maximum fragment length (nt):", value=1000, placeholder="e.g., 1000", key="max_fragment_length")

    # Second parameter
    predict_fnc = right.selectbox("Structure prediction model:", options=["KnotFold", "IPknot", "ProbKnot", "pKiss", "RNAfold", "LinearFold", "MXfold2"], index=0, key="predict_fnc")

    # Submit button
    if st.button("Predict"):
        # Prepare parameters
        predict_fnc = {"KnotFold": knotfold_predict,
                    "IPknot": ipknot_predict,
                    "ProbKnot": probknot_predict,
                    "pKiss": pkiss_predict,
                    "RNAfold": rnafold_predict,
                    "LinearFold": linearfold_predict,
                    "MXfold2": mxfold2_predict}[predict_fnc]
        max_fragment_length = int(max_fragment_length)

        # Loading screen spinner
        with st.spinner("Running DivideFold..."):
            seq = st.session_state.get("rna_sequence")
            pred, frags = dividefold_predict(seq,
                                             predict_fnc=predict_fnc,
                                             max_fragment_length=max_fragment_length,
                                             return_fragments=True)

            # Write results
            if not os.path.exists("results"):
                os.mkdir("results")
            results_path = f"results/{token}.txt"
            with open(results_path, "w") as f:
                f.write(f"{seq}\n{pred}\n{[f.tolist() for f in frags]}\n")

        # Show results
        st.success("Prediction complete!")
        results_page(token)


def results_page(token):
    # Read results
    results_path = f"results/{token}.txt"
    with open(results_path, "r") as f:
        results_data = f.read()
    seq, pred, frags_str = results_data.strip().split("\n")

    # Parse frags back to list
    frags = []
    for this_frag_str in frags_str[3:-3].split("]], [["):
        this_frag = []
        for subfrag_str in this_frag_str.split("], ["):
            start, end = subfrag_str.split(", ")
            this_frag.append([int(start), int(end)])
        frags.append(this_frag)

    # Write prediction
    st.markdown("### 🧬 Predicted secondary structure")
    st.code(f"{seq}\n{pred}", language="text")

    # Write fragments
    st.markdown("### 🧬 Predicted fragments")
    st.code(frags, language="text")

    # Plot structure
    plot_structure(seq, pred, frags)


def plot_structure(seq, pred, frags):
    st.markdown("### Visualization")

    # Sample colors
    def sample_color():
        i = 0
        while True:
            # Use golden ratio to spread hues evenly across [0,1]
            hue = (i * 0.618033988749895) % 1
            r, g, b = colorsys.hsv_to_rgb(hue, 0.65, 0.95)
            r, g, b = tuple(round(x * 255) for x in (r, g, b))
            i += 1
            yield (r, g, b)

    color_sampler = sample_color()
    colors = [next(color_sampler) for _ in range(len(frags))]
    hexadecimal_converter = "#%02x%02x%02x"
    colors = [hexadecimal_converter % c for c in colors]

    # Reorder and index fragments
    frag_idx = {}
    frag_lens = []
    for k, f in enumerate(frags):
        frag_lens.append(0)
        for start, end in f:
            frag_lens[-1] += end - start + 1
            for i in range(start, end + 1):
                frag_idx[i] = k
    frag_order = np.argsort(frag_lens)[::-1]
    frag_rank = np.argsort(frag_order)

    # Build color string for Forna
    colors_str = " ".join([f"{i}:{colors[frag_rank[frag_idx[i]]]}" \
                                    for i in range(1, len(seq) + 1)])

    # Forna
    forna_component(structure = pred, # RNA structure
                    sequence =  seq, # RNA sequence
                    height = 630, # height of the component
                    animation = True, # boolean to enable animation/interactivity
                    zoomable = True, # boolean to enable zooming
                    # label_interval = label_interval, # interval for numbering the nucleotides
                    node_label = True, # boolean to enable node labeling
                    editable = False, # boolean to enable editing
                    color_scheme = "custom", # color scheme for the nucleotides: options are 'structure', 'sequence', 'positions', 'custom'
                    colors = colors_str, # custom colors for the nucleotides: string with space separated index:color pairs
                    )


main_page()
