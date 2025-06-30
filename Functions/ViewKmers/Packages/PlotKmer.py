import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import itertools
import statistics

fixed_colors = [
    "rgb(141, 211, 199)",
    "rgb(255, 255, 179)",
    "rgb(190, 186, 218)",
    "rgb(251, 128, 144)",
    "rgb(128, 177, 211)",
    "rgb(253, 180, 98)",
    "rgb(179, 222, 105)",
    "rgb(0, 255, 255)",
    "rgb(252, 205, 229)",
    "rgb(217, 217, 217)"
] 

def generate_colors(kmers):
    kmer_count = len(kmers)
    if kmer_count > 10:
        raise ValueError("Inputed more than 10 Kmers. Please lower the number of Kmers (<=10).")
    colors = fixed_colors[:kmer_count]
    return {kmers[i]: colors[i] for i in range(kmer_count)}

def plot_kmers(SeqencesDict, KmerList, selected_genes, up_stream, K2M_Dict):
    kmer_colors = generate_colors(KmerList)
    fig = go.Figure()
    y_offset = 0
    gene_data_list = []

    max_seq_len = max(len(seq[0]) for seq in SeqencesDict.values())
    x_range = [-up_stream, max_seq_len - up_stream]

    # 1. Plot gene bar
    for selected_gene in selected_genes:
        NT_sequence = SeqencesDict[selected_gene][0]
        T_sequence = SeqencesDict[selected_gene][1]

        # Non-template strand
        fig.add_trace(go.Scatter(
            x=[i - up_stream for i in range(len(NT_sequence))],
            y=[y_offset] * len(NT_sequence),
            mode='lines',
            line=dict(color="darkgray", width=20),
            name="Non-Template",
            showlegend=False
        ))

        # Template strand
        fig.add_trace(go.Scatter(
            x=[i - up_stream for i in range(len(T_sequence))],
            y=[y_offset + 0.66] * len(T_sequence),
            mode='lines',
            line=dict(color="lightgray", width=20),
            name="Non-Template",
            showlegend=False
        ))

        # Gene label
        fig.add_annotation(
            x=-up_stream - 5,
            y=y_offset + 0.33,
            text=selected_gene,
            showarrow=False,
            xanchor="right",
            font=dict(size=16, color="white"),
            bgcolor="#2C3E50",
            bordercolor="white",
            borderwidth=1,
            borderpad=4,
        )

        gene_data_list.append((selected_gene, NT_sequence, T_sequence, y_offset))
        y_offset += 2

    # 2. Plot Kmer hit
    legend_drawn = set()

    for kmer in KmerList:
        color = kmer_colors[kmer]

        for selected_gene, NT_seq, T_seq, offset in gene_data_list:
            # Non-template strand
            pos = 0
            while pos < len(NT_seq):
                pos = NT_seq.find(kmer, pos)
                if pos == -1:
                    break
                x = list(range(pos - up_stream, pos - up_stream + len(kmer)))
                y = [offset] * len(x)
                fig.add_trace(go.Scatter(
                    x=x,
                    y=y,
                    mode='lines',
                    line=dict(color=color, width=40),
                    name=f"Kmer: {kmer}" if kmer not in legend_drawn else "",
                    showlegend=kmer not in legend_drawn,
                    hovertemplate=(
                        f"Kmer: {kmer}<br>"
                        f"Position: {pos-up_stream} ~ {pos-up_stream + len(kmer) - 1}<br>"
                        f"Kmer2Motif: {', '.join(K2M_Dict[kmer])}<extra></extra>"
                        if K2M_Dict and kmer in K2M_Dict
                        else f"Kmer: {kmer}<br>"
                            f"Position: {pos-up_stream} ~ {pos-up_stream + len(kmer) - 1}<extra></extra>"
                    )
                ))
                legend_drawn.add(kmer)
                pos += 1

            # Template strand
            pos = 0
            while pos < len(T_seq):
                pos = T_seq.find(kmer, pos)
                if pos == -1:
                    break
                x = list(range(pos - up_stream, pos - up_stream + len(kmer)))
                y = [offset + 0.66] * len(x)
                fig.add_trace(go.Scatter(
                    x=x,
                    y=y,
                    mode='lines',
                    line=dict(color=color, width=40),
                    name=f"Kmer: {kmer}" if kmer not in legend_drawn else "",
                    showlegend=kmer not in legend_drawn,
                    hovertemplate=(
                        f"Kmer: {kmer}<br>"
                        f"Position: {pos-up_stream} ~ {pos-up_stream + len(kmer) - 1}<br>"
                        f"Kmer2Motif: {', '.join(K2M_Dict[kmer])}<extra></extra>"
                        if K2M_Dict and kmer in K2M_Dict
                            else f"Kmer: {kmer}<br>"
                                f"Position: {pos-up_stream} ~ {pos-up_stream + len(kmer) - 1}<extra></extra>"
                        )
                    ))
                legend_drawn.add(kmer)
                pos += 1

    # TSS vertical line
    fig.add_shape(
        type="line",
        x0=0, x1=0,
        y0=-1, y1=y_offset - 1 + 1,
        line=dict(color="#2C3E50", width=2, dash="dot"),
        layer="below"
    )

    # Layout
    fig.update_layout(
        title=dict(
            text=(
                "<span style='color:lightgray;'>■: Template strand</span> | "
                "<span style='color:darkgray;'>■: Non-template strand</span>"
            ),
            x=0.5,
            xanchor='center'),
        xaxis_title="Sequence Position",
        yaxis=dict(showticklabels=False),
        showlegend=True,
        height=100 + y_offset * 50,
        margin=dict(r=100)
    )

    return fig


def Show_df(selected_kmers, selected_genes, SeqencesDict, Orientation):
    kmer_colors = generate_colors(selected_kmers)

    if len(selected_kmers) > 1:
        df = pd.DataFrame(index=selected_kmers, columns=selected_kmers)

        for kmer1, kmer2 in itertools.combinations(selected_kmers, 2):
            distances = []
            for gene in selected_genes:
                if Orientation == "NT":
                    sequence = SeqencesDict[gene][0]
                elif Orientation == "T":
                    sequence = SeqencesDict[gene][1]
                pos1 = [sequence.find(kmer1, i) for i in range(len(sequence)) if sequence.find(kmer1, i) != -1]
                pos2 = [sequence.find(kmer2, i) for i in range(len(sequence)) if sequence.find(kmer2, i) != -1]
                if pos1 and pos2:
                    for p1 in pos1:
                        for p2 in pos2:
                            distances.append(abs(p1 - p2))
            if distances:
                avg_distance = sum(distances) / len(distances)
                std_distance = statistics.stdev(distances) if len(distances) > 1 else 0
                value = f"{avg_distance:.2f} ± {std_distance:.2f}"
            else:
                value = "N/A"

            df.loc[kmer1, kmer2] = value
            df.loc[kmer2, kmer1] = value

        for kmer in selected_kmers:
            df.loc[kmer, kmer] = "--"

        def highlight(cell_value, row_kmer, col_kmer):
            if cell_value == "--":
                return "background-color: #eeeeee; color: black;"
            if cell_value == "N/A":
                return "background-color: #aaaaaa; color: white;"
            color1 = kmer_colors[row_kmer]
            color2 = kmer_colors[col_kmer]
            return f"background-color: {color1}; color: white;"

        # apply style
        def style_func(val, row_idx, col_idx):
            row_kmer = df.index[row_idx]
            col_kmer = df.columns[col_idx]
            return highlight(val, row_kmer, col_kmer)

        styled_df = df.style.apply(
            lambda df_: pd.DataFrame(
                [[style_func(df_.iloc[i, j], i, j) for j in range(df_.shape[1])] for i in range(df_.shape[0])],
                index=df_.index,
                columns=df_.columns
            ),
            axis=None
        )

        st.dataframe(styled_df)

def main(new_folder, SeqencesDict, KmerList, up_stream, K2M_Dict):
    if 'KmerList' not in st.session_state:
        st.session_state.KmerList = KmerList

    # 1. select Genes
    selected_genes = st.multiselect(
        "Select gene (Mutiple allowed)", 
        list(SeqencesDict.keys()),
        default=st.session_state.get('selected_genes', SeqencesDict.keys())
    )

    # 2. select Kmers
    selected_kmers = st.multiselect(
        "Select Kmer (Mutiple allowed)", 
        st.session_state.KmerList,
        default=st.session_state.get('selected_kmers', st.session_state.KmerList),
        key="kmer_multiselect"
    )
    st.session_state.selected_kmers = selected_kmers

    if selected_genes:
        st.header("Kmer Distribution in Selected Genes")
        st.text("Click 'Fullscreen' to have better experience")
        fig = plot_kmers(SeqencesDict, selected_kmers, selected_genes, up_stream, K2M_Dict)
        st.plotly_chart(fig)

        # Add Average distance information
        st.header("Average Kmer Spacing in Genes")
        st.subheader("In Non-Template(Coding) Strands")
        Show_df(selected_kmers, selected_genes, SeqencesDict, "NT")
        st.subheader("In Template Strands")
        Show_df(selected_kmers, selected_genes, SeqencesDict, "T")
    else:
        st.warning("Select at least one gene to show the chart.")