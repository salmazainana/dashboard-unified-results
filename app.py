import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Unified Results Dashboard",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 Unified Results Dashboard")

# ── Data loading ──────────────────────────────────────────────────────────────
DATA_DIR = Path(__file__).parent / "data"

@st.cache_data
def load_file(path: str, _mtime: float) -> pd.DataFrame:
    """_mtime is passed so the cache is busted whenever the file changes."""
    df = pd.read_csv(path, sep="\t", compression="gzip" if path.endswith(".gz") else None)
    df.replace("", np.nan, inplace=True)
    for col in df.columns:
        if col != "gene":
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df

def _make_label(filename: str) -> str:
    base = filename.replace(".tsv.gz", "")
    parts = base.split(".unified_results")
    prefix = parts[0].replace("_", " ").upper()
    if len(parts) > 1 and parts[1]:
        filter_map = {"n_cohorts_ge4": "n≥4 cohorts", "MAF_le_001": "MAF≤1%"}
        filter_parts = [filter_map.get(f, f) for f in parts[1].lstrip(".").split(".")]
        return f"{prefix} ({', '.join(filter_parts)})"
    return prefix

result_files = sorted(DATA_DIR.glob("*.tsv.gz"))
file_labels = {_make_label(f.name): str(f) for f in result_files}

if not file_labels:
    st.error("No `.tsv.gz` result files found in the `data/` folder.")
    st.stop()

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("⚙️ Settings")

    selected_label = st.selectbox("Dataset", list(file_labels.keys()))
    selected_path = file_labels[selected_label]
    mtime = Path(selected_path).stat().st_mtime
    df = load_file(selected_path, mtime)

    st.markdown("---")
    p_options = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05]
    p_threshold = st.select_slider(
        "p-value significance threshold",
        options=p_options,
        value=0.05,
        format_func=lambda x: f"{x:.0e}" if x < 0.01 else str(x),
    )

    n_sig = int((df["p_model"] < p_threshold).sum())
    n_total = int(df["p_model"].notna().sum())
    st.metric("Significant genes", f"{n_sig} / {n_total}", help=f"Genes with p_model < {p_threshold}")

    st.markdown("---")
    st.markdown("**Predictors to show**")
    show_constraint    = st.checkbox("Log constraint",     value=True)
    show_pathogenicity = st.checkbox("Log pathogenicity",  value=True)
    show_plof          = st.checkbox("pLoF indicator",     value=True)
    show_missense      = st.checkbox("Missense indicator", value=True)

    st.markdown("---")
    hide_missing_pmodel = st.checkbox("Hide genes without p_model", value=False)

    st.markdown("---")
    if st.button("🔄 Clear cache & reload"):
        st.cache_data.clear()
        st.rerun()

# ── Derive useful columns ─────────────────────────────────────────────────────
df_display = df.copy()
if hide_missing_pmodel:
    df_display = df_display.dropna(subset=["p_model"])

df_display["-log10_p_model"] = -np.log10(df_display["p_model"].replace(0, np.nan))

# Build predictor map based on sidebar toggles
predictors = {}
if show_constraint:
    predictors["Log constraint"] = ("coef_log_constraint", "p_log_constraint")
if show_pathogenicity:
    predictors["Log pathogenicity"] = ("coef_log_pathogenicity", "p_log_pathogenicity")
if show_plof:
    predictors["pLoF indicator"] = ("coef_pLoF_indicator", "p_pLoF_indicator")
if show_missense:
    predictors["Missense indicator"] = ("coef_missense_indicator", "p_missense_indicator")

# ── Tabs ──────────────────────────────────────────────────────────────────────
tab1, tab2, tab3, tab4 = st.tabs(["📋 Table", "📊 Rankings", "🟥 Heatmap", "🔬 Gene Detail"])

# ─────────────────────────────────────────────────────────────────────────────
# TAB 1 — Interactive table
# ─────────────────────────────────────────────────────────────────────────────
with tab1:
    st.subheader("Gene Results Table")
    st.caption(f"Significance threshold: p < {p_threshold}  |  {len(df_display)} genes")

    # Build display df with rounded values
    table_cols = ["gene", "p_model"] + [col for pred in predictors.values() for col in pred]
    table_df = df_display[["gene", "p_model", "-log10_p_model"] + [c for p in predictors.values() for c in p]].copy()
    table_df = table_df.sort_values("p_model", na_position="last")

    def highlight_sig(val, col):
        if "p_" in col and pd.notna(val):
            if val < p_threshold:
                return "background-color: #d4edda; color: #155724"
            elif val < 0.25:
                return "background-color: #fff3cd; color: #856404"
        return ""

    # Style p-value columns
    p_cols = [c for c in table_df.columns if c.startswith("p_")]
    def color_pval(v):
        try:
            if pd.notna(v) and float(v) < p_threshold:
                return "background-color: #d4edda; color: #155724"
            elif pd.notna(v) and float(v) < 0.25:
                return "background-color: #fff3cd; color: #856404"
        except (TypeError, ValueError):
            pass
        return ""

    styled = table_df.style.format(
        {c: "{:.3g}" for c in table_df.columns if c != "gene"}
    ).map(color_pval, subset=p_cols)

    st.dataframe(styled, use_container_width=True, height=500)

    # Download button
    csv = table_df.to_csv(index=False)
    st.download_button("⬇️ Download table as CSV", csv, file_name=f"{selected_label}_results.csv", mime="text/csv")


# ─────────────────────────────────────────────────────────────────────────────
# TAB 2 — Gene rankings
# ─────────────────────────────────────────────────────────────────────────────
with tab2:
    st.subheader("Gene Rankings by Model p-value")

    rank_df = df_display.dropna(subset=["p_model"]).sort_values("-log10_p_model", ascending=True)
    if rank_df.empty:
        st.warning("No genes with a p_model value in this dataset.")
    else:
        rank_df["significant"] = rank_df["p_model"] < p_threshold
        rank_df["color"] = rank_df["significant"].map({True: "#2ecc71", False: "#95a5a6"})

        fig = go.Figure()

        fig.add_trace(go.Bar(
            x=rank_df["-log10_p_model"],
            y=rank_df["gene"],
            orientation="h",
            marker_color=rank_df["color"],
            hovertemplate="<b>%{y}</b><br>-log₁₀(p) = %{x:.3f}<br>p = %{customdata:.3g}<extra></extra>",
            customdata=rank_df["p_model"],
        ))

        # Significance line
        sig_line = -np.log10(p_threshold)
        fig.add_vline(x=sig_line, line_dash="dash", line_color="red",
                      annotation_text=f"p={p_threshold}", annotation_position="top right")

        fig.update_layout(
            height=max(400, len(rank_df) * 22),
            xaxis_title="-log₁₀(p_model)",
            yaxis_title="Gene",
            margin=dict(l=10, r=10, t=30, b=40),
            showlegend=False,
        )
        st.plotly_chart(fig, use_container_width=True)


# ─────────────────────────────────────────────────────────────────────────────
# TAB 3 — Heatmap
# ─────────────────────────────────────────────────────────────────────────────
with tab3:
    st.subheader("Predictor Heatmap — Genes × Predictors")

    if not predictors:
        st.warning("Select at least one predictor in the sidebar.")
    else:
        heatmap_mode = st.radio(
            "Color by:",
            ["-log₁₀(p-value)", "Coefficient value"],
            horizontal=True,
        )

        genes_sorted = df_display.sort_values("p_model", na_position="last")["gene"].tolist()
        pred_labels = list(predictors.keys())

        if heatmap_mode == "-log₁₀(p-value)":
            matrix = np.full((len(genes_sorted), len(pred_labels)), np.nan)
            for j, (label, (coef_col, p_col)) in enumerate(predictors.items()):
                vals = df_display.set_index("gene").reindex(genes_sorted)[p_col].values
                with np.errstate(divide="ignore"):
                    matrix[:, j] = -np.log10(vals.astype(float))

            colorscale = "Reds"
            zmin = 0
            finite = matrix[np.isfinite(matrix)]
            zmax = float(np.percentile(finite, 99)) if len(finite) > 0 else None
            if zmax is not None:
                matrix = np.where(np.isfinite(matrix), matrix, zmax)
            colorbar_title = "-log₁₀(p)"
        else:
            matrix = np.full((len(genes_sorted), len(pred_labels)), np.nan)
            for j, (label, (coef_col, p_col)) in enumerate(predictors.items()):
                vals = df_display.set_index("gene").reindex(genes_sorted)[coef_col].values
                matrix[:, j] = vals.astype(float)

            # Clip extreme outliers for readability
            finite = matrix[np.isfinite(matrix)]
            if len(finite) > 0:
                vmax = np.percentile(np.abs(finite), 95)
                matrix = np.clip(matrix, -vmax, vmax)
            colorscale = "RdBu_r"
            zmin, zmax = None, None
            colorbar_title = "Coefficient"

        fig = go.Figure(data=go.Heatmap(
            z=matrix,
            x=pred_labels,
            y=genes_sorted,
            colorscale=colorscale,
            zmin=zmin,
            zmax=zmax,
            colorbar=dict(title=colorbar_title),
            hoverongaps=False,
            hovertemplate="<b>%{y}</b> — %{x}<br>Value: %{z:.3g}<extra></extra>",
        ))

        fig.update_layout(
            height=max(400, len(genes_sorted) * 22),
            margin=dict(l=10, r=10, t=30, b=60),
            xaxis=dict(side="top"),

        )
        st.plotly_chart(fig, use_container_width=True)


# ─────────────────────────────────────────────────────────────────────────────
# TAB 4 — Gene detail
# ─────────────────────────────────────────────────────────────────────────────
with tab4:
    st.subheader("Gene Detail View")

    all_genes = df_display["gene"].tolist()
    selected_gene = st.selectbox("Select a gene", all_genes)

    if selected_gene:
        row = df_display[df_display["gene"] == selected_gene].iloc[0]

        # Summary metrics
        col1, col2 = st.columns(2)
        with col1:
            p_val = row.get("p_model")
            if pd.notna(p_val):
                st.metric("Model p-value", f"{p_val:.4g}")
            else:
                st.metric("Model p-value", "N/A")
        with col2:
            log_p = row.get("-log10_p_model")
            if pd.notna(log_p):
                st.metric("-log₁₀(p_model)", f"{log_p:.3f}")
            else:
                st.metric("-log₁₀(p_model)", "N/A")

        st.markdown("---")

        if not predictors:
            st.warning("Select at least one predictor in the sidebar.")
        else:
            # Build per-predictor table
            detail_rows = []
            for label, (coef_col, p_col) in predictors.items():
                coef = row.get(coef_col)
                p    = row.get(p_col)
                sig  = "✅" if pd.notna(p) and p < p_threshold else ("⚠️" if pd.notna(p) and p < 0.25 else "")
                detail_rows.append({
                    "Predictor": label,
                    "Coefficient": coef,
                    "p-value": p,
                    "Sig": sig,
                })

            detail_df = pd.DataFrame(detail_rows)

            # Grouped bar chart: coefficients
            fig = go.Figure()
            colors = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12"]
            for i, r in detail_df.iterrows():
                coef = r["Coefficient"]
                p    = r["p-value"]
                bar_color = colors[i % len(colors)]
                opacity = 1.0 if (pd.notna(p) and p < p_threshold) else 0.4
                coef_str = f"{coef:.4g}" if pd.notna(coef) else "N/A"
                p_str = f"{p:.4g}" if pd.notna(p) else "N/A"
                fig.add_trace(go.Bar(
                    name=r["Predictor"],
                    x=[r["Predictor"]],
                    y=[coef] if pd.notna(coef) else [0],
                    marker_color=bar_color,
                    opacity=opacity,
                    hovertemplate=(
                        f"<b>{r['Predictor']}</b><br>"
                        f"Coefficient: {coef_str}<br>"
                        f"p-value: {p_str}<extra></extra>"
                    ),
                ))

            fig.update_layout(
                title=f"Predictor Coefficients — {selected_gene}",
                yaxis_title="Coefficient",
                xaxis_title="Predictor",
                showlegend=False,
                height=400,
            )
            fig.add_hline(y=0, line_color="black", line_width=0.8)
            st.plotly_chart(fig, use_container_width=True)
            st.caption(f"🟢 Opaque bars = p < {p_threshold} (significant). Faded bars = not significant or missing p-value.")

            # Detail table
            st.dataframe(
                detail_df.style.format({"Coefficient": "{:.4g}", "p-value": "{:.4g}"}, na_rep="—"),
                use_container_width=True,
                hide_index=True,
            )
