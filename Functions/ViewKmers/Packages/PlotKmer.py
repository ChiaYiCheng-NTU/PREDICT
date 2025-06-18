import streamlit as st
import plotly.graph_objects as go

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

def plot_kmers(SeqencesDict, KmerList, selected_genes, up_stream):
    kmer_colors = generate_colors(KmerList)
    
    fig = go.Figure()
    y_offset = 0
    legend_added = set()
    
    max_seq_len = max(len(seq) for seq in SeqencesDict.values())
    x_range = [-up_stream, max_seq_len - up_stream]  # count from upstream (negative)

    for selected_gene in selected_genes:
        sequence = SeqencesDict[selected_gene]

        # gene body gray bar
        fig.add_trace(go.Scatter(
            x=[i - up_stream for i in range(len(sequence))],
            y=[y_offset] * len(sequence),
            mode='lines',
            line=dict(color="lightgray", width=20),
            name=f"Gene: {selected_gene}",
            showlegend=False
        ))

        # gene label
        fig.add_annotation(
            x=-up_stream - 5,
            y=y_offset,
            text=selected_gene,
            showarrow=False,
            xanchor="right",
            font=dict(size=16, color="white"),
            bgcolor="#2C3E50",
            bordercolor="white",
            borderwidth=1,
            borderpad=4,
        )

        # plot kmer distribution
        for kmer in KmerList:
            color = kmer_colors[kmer]
            start = 0
            while start < len(sequence):
                start = sequence.find(kmer, start)
                if start == -1:
                    break
                show_legend = kmer not in legend_added
                fig.add_trace(go.Scatter(
                    x=[start - up_stream, start - up_stream + len(kmer) - 1],
                    y=[y_offset, y_offset],
                    mode='lines',
                    line=dict(color=color, width=40),
                    name=f"Kmer: {kmer}" if show_legend else "",
                    showlegend=show_legend
                ))
                if show_legend:
                    legend_added.add(kmer)
                start += 1

        y_offset += 2

    # TSS Vertical line
    fig.add_shape(
        type="line",
        x0=0, x1=0,
        y0=-1, y1=y_offset - 1 + 1,
        line=dict(color="#2C3E50", width=2, dash="dot"),
        layer="below"
    )

    # Layout
    fig.update_layout(
        title=f"Kmer Distribution in Selected Genes",
        xaxis_title="Sequence Position",
        yaxis=dict(showticklabels=False),
        showlegend=True,
        height=100 + y_offset * 50,
        margin=dict(r=100)
    )

    return fig


def main(new_folder, SeqencesDict, KmerList, up_stream):
    if 'KmerList' not in st.session_state:
        st.session_state.KmerList = KmerList

    # 1. select Genes
    selected_genes = st.multiselect(
        "Select gene (Mutiple allowed)", 
        list(SeqencesDict.keys())
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
        fig = plot_kmers(SeqencesDict, selected_kmers, selected_genes, up_stream)
        st.plotly_chart(fig)

        # Add Average distance information
        kmer_colors = generate_colors(selected_kmers)
        st.subheader("Average Kmer Spacing in Genes")
        if len(selected_kmers) > 1:
            import itertools
            import statistics
            for kmer_pair in itertools.combinations(selected_kmers, 2):
                kmer1, kmer2 = kmer_pair
                color1, color2 = kmer_colors[kmer1], kmer_colors[kmer2]
                distances = []
                for gene in selected_genes:
                    sequence = SeqencesDict[gene]
                    pos1 = [sequence.find(kmer1, i) for i in range(len(sequence)) if sequence.find(kmer1, i) != -1]
                    pos2 = [sequence.find(kmer2, i) for i in range(len(sequence)) if sequence.find(kmer2, i) != -1]
                    if pos1 and pos2:
                        for p1 in pos1:
                            for p2 in pos2:
                                distances.append(abs(p1 - p2))
                if distances:
                    avg_distance = sum(distances) / len(distances)
                    std_distance = statistics.stdev(distances) if len(distances) > 1 else 0
                    st.markdown(
                        f"<span style='color:{color1};'>■</span>  Avg distance of {kmer1} and "
                        f"<span style='color:{color2};'>■</span> {kmer2}:{avg_distance:.2f}, "
                        f"Distance std:{std_distance:.2f}",
                        unsafe_allow_html=True
                    )
                else:
                    st.markdown(
                        f"<span style='color:{color1};'>■</span> Avg distance of {kmer1} and "
                        f"<span style='color:{color2};'>■</span> {kmer2}: N/A, "
                        f"Distance std: N/A",
                        unsafe_allow_html=True
                    )
    else:
        st.warning("Select at least one gene to show the chart.")