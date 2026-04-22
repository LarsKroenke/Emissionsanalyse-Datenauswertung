import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

# table part


def ex2df(file):
    return pd.read_csv(file, sep='\t', engine='python')


def concat_df_7(probe):
    file_1 = f'2k-klebstoffe/component-tables/7-days/{probe}_1.txt'
    file_2 = f'2k-klebstoffe/component-tables/7-days/{probe}_2.txt'
    file_3 = f'2k-klebstoffe/component-tables/7-days/{probe}_3.txt'

    df_1 = ex2df(file_1)
    df_2 = ex2df(file_2)
    try:
        df_3 = ex2df(file_3)
    except FileNotFoundError:
        df_3 = pd.DataFrame()

    return pd.concat([df_1, df_2, df_3], ignore_index=True)


def concat_df_28(probe):
    file_1 = f'2k-klebstoffe/component-tables/28-days/{probe}_1.txt'
    file_2 = f'2k-klebstoffe/component-tables/28-days/{probe}_2.txt'
    file_3 = f'2k-klebstoffe/component-tables/28-days/{probe}_3.txt'

    df_1 = ex2df(file_1)
    df_2 = ex2df(file_2)
    df_3 = ex2df(file_3)

    return pd.concat([df_1, df_2, df_3], ignore_index=True)


def group_dfs(df_all):

    df_all["Retentionszeit"] = (df_all["Retention Time"].astype(str).str.replace(
        ',', '.', regex=False).pipe(pd.to_numeric, errors='coerce') / 0.11).round() * 0.11

    df_all["Amount"] = (df_all["Amount "].astype(str).str.replace(
        ',', '.', regex=False).pipe(pd.to_numeric, errors='coerce'))

    df_all["Prob"] = (df_all["Prob "].astype(str).str.replace(
        ',', '.', regex=False).pipe(pd.to_numeric, errors='coerce'))

    df_all["Est. Amt."] = (df_all["Est. Amt."].astype(str).str.replace(
        ',', '.', regex=False).pipe(pd.to_numeric, errors='coerce'))

    df_all["RSI "] = (df_all["RSI "].astype(str).str.replace(
        ',', '.', regex=False).pipe(pd.to_numeric, errors='coerce'))

    df_all["Peak Name"] = (
        df_all["Peak Name "]
    )

    summary = (
        df_all
        .groupby(["Peak Name", "Retentionszeit"], dropna=False)
        .agg(
            n=("Amount", "count"),
            mean=("Amount", "mean"),
            td8=("Est. Amt.", "mean"),
            sx=("Amount", "std"),
            P=("Prob", "mean"),
            RSI=("RSI ", "mean")
        )
        .reset_index()
    )
    summary[["mean", "td8", "sx", "P", "RSI"]] = summary[[
        "mean", "td8", "sx", "P", "RSI"]].round(1)
    summary.rename(columns={'Peak Name': 'Komponente', 'Retentionszeit': 'Retentionszeit', 'mean': 'Toluol Äq.',
                   'sx': 'Standardabweichung', 'P': 'Probability', 'td8': 'int. Std. Äq.'}, inplace=True)
    return summary.sort_values('Retentionszeit')


# plot the data
def chrom2df(file_path):
    # Read file
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # Find start of chromatogram numeric data
    start_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Time (min)"):
            start_idx = i + 1
            break

    # Read numeric section
    df = pd.read_csv(
        file_path,
        sep="\t",
        # skiprows=start_idx,
        names=["Time", "Step", "Value"],
        engine="python"
    )

    # Remove commas and convert to float
    df["Time"] = pd.to_numeric(df["Time"], errors="coerce")
    df["Value"] = (
        df["Value"]
        .astype(str)
        .str.replace(",", "", regex=False)
    )
    df["Value"] = pd.to_numeric(df["Value"], errors="coerce")

    df = df.dropna()
    return df


def read_plot_data_7(kleber):
    file_7 = f'2k-klebstoffe/spectra/7-days/{kleber}.txt'
    return chrom2df(file_7)


def read_plot_data_28(kleber):
    file_28 = f'2k-klebstoffe/spectra/28-days/{kleber}.txt'
    return chrom2df(file_28)


def plot_spectra(df_7, df_28, slider, on_7, on_28, x_offset, y_offset):
    x_7 = df_7["Time"]
    y_7 = df_7["Value"]
    x_28 = df_28["Time"] + x_offset
    y_28 = df_28["Value"] - y_offset

    fig, ax = plt.subplots(figsize=(10, 4))
    if on_7:
        ax.plot(x_7, y_7, linewidth=0.25, color='black', label='7 Tage')
    if on_28:
        ax.plot(x_28, y_28, linewidth=0.25, color='red', label='28 Tage')
    ax.set_xlabel('Retentionszeit in min')
    ax.set_ylabel('Intensität')
    ax.legend()
    # ax.spines.right.set_color(None)
    # ax.spines.top.set_color(None)
    # ax.spines.bottom.set_bounds(slider)
    # ax.spines.left.set_bounds(0, 5.5*10**9)
    ax.set_xticks(np.arange(0, 60, 1), minor=True)
    ax.set_xticks(np.arange(0, 61, 5))
    ax.set_xlim(slider)
    return fig


# streamlit code
option = st.sidebar.selectbox(
    'Analyten auswählen',
    ('pn7444', 'lord', '9075l', 'ep6055', 'dp6330')
)

# plot integration
slider = st.sidebar.slider(
    'Retentionszeitbereich (min)',
    0, 60, (0, 60)
)

st.sidebar.write('Spektrum ein- ausblenden')
on_7 = st.sidebar.toggle('7 Tage', value=True)
on_28 = st.sidebar.toggle('28 Tage', value=True)

st.sidebar.write('Chromatogramm offset justieren')
x_offset = st.sidebar.slider('X offset', 0.0, 2.0, step=0.05, value=0.5)
y_offset = st.sidebar.slider(
    'Y offset', 0.0*10**9, 0.5*10**9, step=0.01*10**9, value=0.1*10**9)

df_7 = read_plot_data_7(option)
df_28 = read_plot_data_28(option)

plot = plot_spectra(df_7, df_28, slider, on_7, on_28, x_offset, y_offset)

st.write(plot)

# 7 days calculations
df_all_7 = concat_df_7(option)
component_table_7 = group_dfs(df_all_7)

# 28 days calculations
df_all_28 = concat_df_28(option)
component_table_28 = group_dfs(df_all_28)

st.write('Komponententabelle (7 Tage)')
st.write(component_table_7)

st.write('Komponententabelle (28 Tage)')
st.write(component_table_28)
