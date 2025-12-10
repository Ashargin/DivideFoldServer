import streamlit as st
import streamlit.components.v1 as components
from st_forna_component import forna_component
import os
import numpy as np
import colorsys
import uuid

# Session state initialization
for key, default in [("rna_sequence", ""), ("use_example_clicked", False)]:
    if key not in st.session_state:
        st.session_state[key] = default
st.set_page_config(layout="wide")

# Title
st.title("DivideFold+ web server")
st.markdown("##### DivideFold+ partitions a long RNA sequence into smaller fragments while preserving as much of the structure as possible and then uses a secondary structure prediction model on the fragments.")

EXAMPLE_SEQUENCE = "UUGAUGGAGAGUUUGAUCCUGGCUCAGGACGAACGCUGGCGGCGUGCUUAACACAUGCAAGUCGAGCGGGGAGCUUCGGCUCUCAGCGGCGAACGGGUGAGUAACACGUGGGCAACCUGCCCCGAGCUCUGGAAUAACCUCGGGAAACCGGGGCUAAUGCCGGAUAUGACAUUGCCGGGCAUCUGGUGGUGUGGAAAGAUUUAUCGGCUCGGGAUGGGCCCGCGGCCUAUCAGCUUGUUGGUGGGGUGAUGGCCUACCAAGGCGACGACGGGUAGCCGGCCUGAGAGGGCGACCGGCCACACUGGGACUGAGACACGGCCCAGACUCCUACGGGAGGCAGCAGUGGGGAAUAUUGCGCAAUGGGCGGAAGCCUGACGCAGCGACGCCGCGUGGGGGAUGACGGCCUUCGGGUUGUAAACCUCUUUCAGCAGGGACGAAGCGAGAGUGACGGUACCUGCAGAAGAAGCACCGGCCAACUACGUGCCAGCAGCCGCGGUAAUACGUAGGGUGCAAGCGUUGUCCGGAAUUAUUGGGCGUAAAGAGCUCGUAGGCGGCCUGUUGCGUCGGCUGUGAAAACCCGGGGCUCAACUCCGGGCCUGCAGUCGAUACGGGCAGGCUAGAGUCCGGCAGGGGAGACUGGAAUUCCUGGUGUAGCGGUGAAAUGCGCAGAUAUCAGGAGGAACACCGGUGGCGAAGGCGGGUCUCUGGGCCGGAACUGACGCUAAGGAGCGAAAGCGUGGGGAGCGAACAGGAUUAGAUACCCUGGUAGUCCACGCCGUAAACGUUGGGCGCUAGGUGUGGGGGACCUUCCACGGCCUCCGUGCCGCAGCUAACGCAUUAAGCGCCCCGCCUGGGGAGUACGGCCGCAAGGCUAAAACUCAAAGGAAUUGACGGGGGCCCGCACAAGCGGCGGAGCAUGUGGCUUAAUUCGAUGCAACGCGAAGAACCUUACCAGGGCUUGACAUGCAGGGAAAUCUCGUAGAGAUACGGGGUCCGUAAGGGUCCUGCACAGGUGGUGCAUGGCUGUCGUCAGCUCGUGUCGUGAGAUGUUGGGUUAAGUCCCGCAACGAGCGCAACCCUCGUCCUAUGUUGCCAGCGAGUCGUGUCGGGGACUCAUAGGAGACUGCCGGGGUCAACUCGGAGGAAGGUGGGGAUGACGUCAAGUCAUCAUGCCCCUUACGUCCUGGGCUGCACACAUGCUACAAUGGCCGGUACAAAGGGCUGCGAUGCCGUGAGGUGGAGCGAAUCCCAAAAAGCCGGUCUCAGUUCGGAUCGGGGUCUGCAACUCGACCCCGUGAAGUCGGAGUCGCUAGUAAUCGCAGAUCAGCAAUGCUGCGGUGAAUACGUUCCCGGGCCUUGUACACACCGCCCGUCACGUCACGAAAGUCGGUAACACCCGAAGCCGGUGGCCUAACCCUUGUGGGGGGAGCCGUCGAAGGUGGGACCGGCGAUUGGGACGAAGUCGUAACAAGGUAGCCGUACCGGAAGGUGCGGCUGGAUCACCUCCUUUCU"

# Get params
if "token" not in st.session_state:
    st.session_state["token"] = str(uuid.uuid4())
token = st.session_state["token"]

with st.spinner("Loading libraries..."):
    from dividefold.predict import dividefold_predict, knotfold_predict, ipknot_predict, probknot_predict, rnafold_predict, linearfold_predict, mxfold2_predict


def main_page():
    # Settings
    parameter_font_size = 18

    # Set example sequence
    if st.session_state.get("use_example_clicked"):
        st.session_state.rna_sequence = EXAMPLE_SEQUENCE
        st.session_state.use_example_clicked = False

    # Sequence input
    st.markdown(f'<h3 style="margin-bottom:-100px; font-size:{parameter_font_size}px;">RNA sequence:</h3>', unsafe_allow_html=True)
    sequence_input = st.text_area("", value="", height=200, key="rna_sequence")

    # Button to set example sequence
    if st.button("Use example sequence"):
        st.session_state.use_example_clicked = True
        st.rerun()

    left, right = st.columns(2)

    # Maximum fragment length
    left.markdown(f'<h3 style="margin-bottom:-100px; font-size:{parameter_font_size}px;">Maximum fragment length (nt):</h3>', unsafe_allow_html=True)
    max_fragment_length = left.number_input("", min_value=50, max_value=1000, step=1, value=1000, placeholder="e.g., 1000", key="max_fragment_length", help="The maximum length allowed for the fragments, in nucleotides.\nShould be between 50 and 1,000.")

    # Second parameter
    right.markdown(f'<h3 style="margin-bottom:-100px; font-size:{parameter_font_size}px;">Structure prediction model:</h3>', unsafe_allow_html=True)
    predict_fnc_selectbox = right.selectbox("", options=["KnotFold (best accuracy, slower)", "IPknot (good accuracy, faster)", "LinearFold", "RNAfold", "MXfold2", "ProbKnot"], index=0, key="predict_fnc", help="The secondary structure prediction model used on the fragments. We recommend KnotFold for higher accuracy or IPknot for faster prediction.")

    # Submit button
    if st.button("Predict"):
        # Prepare parameters
        seq = st.session_state.get("rna_sequence")
        predict_fnc = {"KnotFold (best accuracy, slower)": knotfold_predict,
                    "IPknot (good accuracy, faster)": ipknot_predict,
                    "ProbKnot": probknot_predict,
                    "RNAfold": rnafold_predict,
                    "LinearFold": linearfold_predict,
                    "MXfold2": mxfold2_predict}[predict_fnc_selectbox]

        # Validate inputs
        if not seq:
            st.error("Please input an RNA sequence.")
        elif len(seq) <= max_fragment_length:
            st.error("The RNA sequence must be longer than the maximum fragment length.")
        else:

            # Loading screen spinner
            with st.spinner("Running DivideFold+..."):
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
    frag_lens = []
    for this_frag_str in frags_str[3:-3].split("]], [["):
        this_frag = []
        frag_lens.append(0)
        for subfrag_str in this_frag_str.split("], ["):
            start, end = subfrag_str.split(", ")
            this_frag.append((int(start), int(end)))
            frag_lens[-1] += int(end) - int(start) + 1
        frags.append(this_frag)

    # Reorder fragments
    frag_order = np.argsort(frag_lens)[::-1]
    frags = [frags[i] for i in frag_order]

    # Sample colors
    def sample_color():
        i = 0
        while True:
            # Use golden ratio to spread hues evenly across [0,1]
            hue = (i * 0.618033988749895) % 1
            r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 0.75)
            r, g, b = tuple(round(x * 255) for x in (r, g, b))
            i += 1
            yield (r, g, b)

    color_sampler = sample_color()
    colors = [next(color_sampler) for _ in range(len(frags))]
    hexadecimal_converter = "#%02x%02x%02x"
    colors = [hexadecimal_converter % c for c in colors]

    # Write prediction
    write_prediction(seq, pred, frags, colors)

    # Write fragments
    write_fragments(frags, colors)

    # Plot structure
    plot_structure(seq, pred, frags, colors)


def write_colored_text(texts, colors, with_comma=False):
    if with_comma:
        formatted_texts = [[] for _ in range(len(texts))]
        formatted_colors = []
        for i in range(len(colors)):
            formatted_colors.append(colors[i])
            if i < len(colors) - 1:
                formatted_colors.append("#000000")
            for j in range(len(texts)):
                formatted_texts[j].append(texts[j][i])
                if i < len(colors) - 1:
                    formatted_texts[j].append(", ")
        texts = formatted_texts
        colors = formatted_colors

    html_blocks = [[f'<span style="color:{c};">{t}</span>' for t, c in zip(row_texts, colors)] for row_texts in texts]
    row_sep = "\n" + " " * 4
    html_color_text = row_sep.join(["".join(row_hb) for row_hb in html_blocks])
    st.markdown(f"""
    <div style="
        background-color: #f0f0f0;
        padding: 15px;
        border-radius: 5px;
        font-family: monospace;
        white-space: pre;       /* Prevent wrapping */
        overflow-x: auto;       /* Horizontal scroll */
        overflow-y: hidden;     /* Prevent vertical scroll */
    ">{html_color_text}</div>
    """, unsafe_allow_html=True)


def write_prediction(seq, pred, frags, colors):
    st.markdown("### 🧬 Predicted secondary structure")
    subfrag_colors = [(subf, c) for f, c in zip(frags, colors) for subf in f]
    subfrag_colors = sorted(subfrag_colors, key=lambda x: x[0][0])
    subfrags, subcolors = zip(*subfrag_colors)
    subtxts = [[txt[subf[0] - 1:subf[1] - 1] for subf in subfrags] for txt in [seq, pred]]
    write_colored_text(subtxts, subcolors)

def write_fragments(frags, colors):
    st.markdown("### ✂️ Predicted fragments")
    write_colored_text([frags], colors, with_comma=True)

def plot_structure(seq, pred, frags, colors):
    st.markdown("### 🔎 Visualization")

    # Index fragments
    frag_idx = {}
    for k, f in enumerate(frags):
        for start, end in f:
            for i in range(start, end + 1):
                frag_idx[i] = k

    # Build color string for Forna
    colors_str = " ".join([f"{i}:{colors[frag_idx[i]]}" for i in range(1, len(seq) + 1)])

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
