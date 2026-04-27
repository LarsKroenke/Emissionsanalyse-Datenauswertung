import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
import numpy as np

st.set_page_config(
    page_title="Full Width App",
    page_icon="🖥️",
    layout="wide",  # Enables full-width layout
    initial_sidebar_state="expanded"
)


def main():
    option = st.sidebar.selectbox(
        'Analyten auswählen',
        ('pn7444', 'lord', '9075l', 'ep6055', 'dp6330')
    )

    option_hit = st.sidebar.selectbox(
        'Hit der Spektrenbibliothek auswählen', (1, 2, 3))

    hit_argument = ['Peak Name', 'Peak Name Hit 2', 'Peak Name Hit 3']
    df_summary = read_file(option, '7')
    df_stat = prepare_data(df_summary, hit_argument[option_hit-1])

    df_summary_28 = read_file(option, '28')
    df_stat_28 = prepare_data(df_summary_28, hit_argument[option_hit-1])

    # plot integration
    x, y, x_28, y_28 = read_raw_spectra(option)

    st.sidebar.write('Spektrum ein- ausblenden')
    on_7 = st.sidebar.toggle('7 Tage', value=True)
    on_28 = st.sidebar.toggle('28 Tage', value=False)

    st.sidebar.write('Toggle Annotations')
    annotation = st.sidebar.toggle('Annotation', value=True)

    slider = st.sidebar.slider(
        'Retentionszeitbereich (min)',
        0, 60, (0, 60)
    )
    if on_28:
        slider_y = st.sidebar.slider(
            'Intensitätsbereich',
            0.0, max(max(y), max(y_28)), (0.0, max(max(y), max(y_28)))
        )
    else:
        slider_y = st.sidebar.slider(
            'Intensitätsbereich',
            0.0, max(y), (0.0, max(y))
        )

    st.sidebar.write('Chromatogramm offset justieren')
    x_offset = st.sidebar.slider('X offset', 0.0, 2.0, step=0.05, value=0.0)
    y_offset = st.sidebar.slider(
        'Y offset', 0.0*10**9, 0.5*10**9, step=0.01*10**9, value=0.0)

    plot = plot_spectra(x, y, x_28, y_28, on_7, on_28,
                        slider, slider_y, df_stat, x_offset, y_offset, annotation)

    # Inhalte schreiben
    # st.header('Auswertung der Dreifachbestimmung')
    st.write(plot)
    st.write(
        '7 Tage Messung (gruppiert nach Retentionszeit und übereinstimmendem Hit)')
    st.write(df_stat)
    st.write(
        '28 Tage Messung (gruppiert nach Retentionszeit und übereinstimmendem Hit)')
    st.write(df_stat_28)


def read_file(probe, days):
    base_path = f'2k-klebstoffe/component-tables/{days}-days'
    dfs = []

    for i in range(1, 4):
        try:
            df = pd.read_csv(f"{base_path}/{probe}_{i}.txt",
                             sep='\t', engine='python')
            dfs.append(df)
        except FileNotFoundError:
            continue

    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def prepare_data(df_all, hit):

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
        .groupby([hit, "Retentionszeit"], dropna=False)
        .agg(
            n=("Amount", "count"),
            mean=("Amount", "mean"),
            td8=("Est. Amt.", "mean"),
            sx=("Amount", "std"),
            P=("Prob", "mean"),
            RSI=("RSI ", "mean")
        )
    )
    summary[["mean", "td8", "sx", "P", "RSI"]] = summary[[
        "mean", "td8", "sx", "P", "RSI"]].round(1)
    summary.rename(columns={'Peak Name': 'Komponente', 'Retentionszeit': 'Retentionszeit', 'mean': 'Toluol Äq.',
                   'sx': 'Standardabweichung', 'P': 'Probability', 'td8': 'int. Std. Äq.'}, inplace=True)
    return summary.sort_values('Retentionszeit').reset_index()


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


def read_raw_spectra(kleber):
    file_7 = f'2k-klebstoffe/spectra/7-days/{kleber}.txt'
    df_7 = chrom2df(file_7)

    file_28 = f'2k-klebstoffe/spectra/28-days/{kleber}.txt'
    df_28 = chrom2df(file_28)

    x_7 = df_7["Time"]
    y_7 = df_7["Value"]
    x_28 = df_28["Time"]  # + x_offset
    y_28 = df_28["Value"]  # - y_offset

    return x_7, y_7, x_28, y_28


def plot_spectra(x_7, y_7, x_28, y_28, on_7, on_28, slider, slider_y, df_stat, x_offset, y_offset, annotate):
    from scipy.interpolate import interp1d

    fig, ax = plt.subplots(figsize=(10, 4))
    if on_7:
        ax.plot(x_7, y_7, linewidth=0.25, color='black', label='7 Tage')
    if on_28:
        ax.plot(x_28 + x_offset, y_28 + y_offset,
                linewidth=0.25, color='red', label='28 Tage')
    ax.set_xlabel('Retentionszeit in min')
    ax.set_ylabel('Intensität')
    ax.legend()
    ax.spines.right.set_color(None)
    ax.spines.top.set_color(None)
    ax.spines.bottom.set_bounds(slider)
    ax.spines.left.set_bounds(min(slider_y), max(slider_y))
    ax.set_xticks(np.arange(min(slider), max(slider), 1), minor=True)
    ax.set_xticks(np.arange(min(slider), max(slider)+1, 5))
    ax.set_xlim(slider)
    ax.set_ylim(min(slider_y)-(max(slider_y)/25), max(slider_y))

    # df_measurement_1
    # df_stat["Retention Time"] = (df_stat["Retention Time"].astype(str).str.replace(
    # ',', '.', regex=False).pipe(pd.to_numeric, errors='coerce'))

    # df_stat
    df_stat["Retentionszeit"] = (df_stat["Retentionszeit"].astype(str).str.replace(
        ',', '.', regex=False).pipe(pd.to_numeric, errors='coerce'))

    # for i in range(len(df_stat)):
    #     target_time = df_stat.iloc[i, 2]
    #     # intensity = np.interp(target_time, x_7, y_7)

    #     f = interp1d(x_7, y_7, kind='cubic', fill_value="extrapolate")
    #     intensity = f(target_time)
    #     ax.annotate(f'{i}', xy=(target_time, intensity), fontsize=6)
    #     num += 1

    if annotate & on_7:
        window = 0.05  # adjust

        x = np.asarray(x_7)
        y = np.asarray(y_7)

        for i in range(len(df_stat)):
            target_time = df_stat.iloc[i, 1]

            # mask safely
            mask = (x >= target_time - window) & (x <= target_time + window)

            if not np.any(mask):
                continue  # skip if nothing found

            local_x = x[mask]
            local_y = y[mask]

            if len(local_y) == 0:
                continue

            idx = np.argmax(local_y)

            peak_time = local_x[idx]
            peak_intensity = local_y[idx]

            ax.annotate(
                f'{i}',
                xy=(peak_time, peak_intensity),
                xytext=(0, 5),
                textcoords='offset points',
                ha='center',
                fontsize=5
            )
    return fig


if __name__ == "__main__":
    main()
