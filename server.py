import os
import pandas as pd
import numpy as np
import base64
import io
from dash import Dash, html, dcc, dash_table, Input, Output, State
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import scipy.stats as stats
from flask import Flask, send_file

# ── Statistical helpers ───────────────────────────────────────────────────────

def benjamini_hochberg(p_values):
    p_values = np.array(p_values)
    n = len(p_values)
    ranked = np.argsort(p_values)
    sorted_p = p_values[ranked]
    adj = np.empty_like(p_values)
    for i, idx in enumerate(ranked):
        adj[idx] = min(1, sorted_p[i] * n / (i + 1))
    for i in range(n - 2, -1, -1):
        adj[ranked[i]] = min(adj[ranked[i]], adj[ranked[i + 1]])
    return adj

# ── App setup ─────────────────────────────────────────────────────────────────

port = int(os.environ.get("PORT", 10000))

app = Dash(
    __name__,
    external_stylesheets=[
        dbc.themes.MATERIA,
        "https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&display=swap",
    ],
)
server = app.server
app.config.suppress_callback_exceptions = True

# ── Custom CSS (minimal, non-conflicting) ─────────────────────────────────────

app.index_string = """
<!DOCTYPE html>
<html>
  <head>
    {%metas%}
    <title>Module Gene Overlap Analysis</title>
    {%favicon%}
    {%css%}
    <style>
      body, .dash-dropdown, input, label, p, h1, h2, h3, h4, h5, h6 {
        font-family: 'Roboto', sans-serif !important;
      }
      body { background: #ECEFF1; }

      /* App bar */
      #appbar {
        background: linear-gradient(135deg, #004D40 0%, #00897B 100%);
        box-shadow: 0 2px 6px rgba(0,0,0,.30);
        padding: 18px 24px;
        border-radius: 3px;
      }
      #appbar h1 {
        color: #fff;
        font-size: 1.25rem;
        font-weight: 500;
        letter-spacing: .01em;
        margin: 0;
      }
      #appbar p {
        color: rgba(255,255,255,.65);
        font-size: .78rem;
        font-weight: 400;
        margin: 2px 0 0;
      }

      /* Control card */
      #control-card {
        background: #fff;
        border: none !important;
        border-radius: 2px !important;
        box-shadow: 0 1px 3px rgba(0,0,0,.12), 0 1px 2px rgba(0,0,0,.24);
      }
      #control-card .card-body { padding: 0; }

      /* Sections inside control card */
      .ctrl-section {
        padding: 16px 20px;
        border-bottom: 1px solid #ECEFF1;
      }
      .ctrl-section:last-child { border-bottom: none; }
      .ctrl-label {
        font-size: .7rem;
        font-weight: 500;
        letter-spacing: .1em;
        text-transform: uppercase;
        color: #78909C;
        margin-bottom: 8px;
      }

      /* Upload zone */
      .upload-zone {
        width: 100%;
        border: 2px dashed #80CBC4 !important;
        border-radius: 4px !important;
        background: #F0FAF9;
        transition: all .2s ease;
        cursor: pointer;
      }
      .upload-zone:hover {
        border-color: #00897B !important;
        background: #E0F2F1 !important;
      }

      /* Run button */
      #run-analysis {
        background: #00897B !important;
        border: none !important;
        border-radius: 2px !important;
        font-weight: 500 !important;
        letter-spacing: .08em !important;
        text-transform: uppercase !important;
        font-size: .85rem !important;
        padding: 12px !important;
        width: 100% !important;
        box-shadow: 0 2px 4px rgba(0,137,123,.35) !important;
        transition: background .15s, box-shadow .15s !important;
      }
      #run-analysis:hover {
        background: #00695C !important;
        box-shadow: 0 4px 8px rgba(0,137,123,.4) !important;
      }

      /* Output cards */
      .out-card {
        background: #fff;
        border: none !important;
        border-radius: 2px !important;
        box-shadow: 0 1px 3px rgba(0,0,0,.12), 0 1px 2px rgba(0,0,0,.24);
      }
      .out-card-header {
        padding: 14px 16px 12px;
        background: #37474F;
        color: #fff;
        font-size: .72rem;
        font-weight: 500;
        letter-spacing: .1em;
        text-transform: uppercase;
        border-radius: 2px 2px 0 0;
      }
      .out-card-body { padding: 0; background: #fff; border-radius: 0 0 2px 2px; }

      /* Stat chips */
      .stat-chips-row { display: flex; gap: 12px; margin-bottom: 16px; }
      .stat-chip {
        flex: 1;
        display: flex;
        flex-direction: column;
        align-items: center;
        padding: 14px 16px;
        background: #fff;
        border-radius: 2px;
        box-shadow: 0 1px 3px rgba(0,0,0,.12), 0 1px 2px rgba(0,0,0,.24);
      }
      .stat-chip .num  { font-size: 1.6rem; font-weight: 300; line-height: 1.1; }
      .stat-chip .lbl  { font-size: .68rem; font-weight: 500; letter-spacing: .06em;
                         text-transform: uppercase; color: #90A4AE; margin-top: 2px; }

      /* Gene chips */
      .gene-tag {
        display: inline-block;
        padding: 1px 9px;
        margin: 2px 2px;
        background: #E0F2F1;
        border: 1px solid #80CBC4;
        border-radius: 12px;
        font-size: .75rem;
        color: #00695C;
        font-family: 'Courier New', monospace;
      }

      /* Table */
      .dash-table-container .dash-spreadsheet-container .dash-spreadsheet-inner td {
        font-size: 13px !important;
        font-family: 'Roboto', sans-serif !important;
      }
      .download-btn {
        background: #00897B !important;
        border: none !important;
        border-radius: 2px !important;
        font-size: .75rem !important;
        font-weight: 500 !important;
        letter-spacing: .07em !important;
        text-transform: uppercase !important;
      }
    </style>
  </head>
  <body>
    {%app_entry%}
    <footer>{%config%}{%scripts%}{%renderer%}</footer>
  </body>
</html>
"""

# ── Data loading ──────────────────────────────────────────────────────────────

def load_builtin_datasets():
    try:
        xlsx = pd.ExcelFile("data/data.xlsx")
        out = {}
        for sheet in xlsx.sheet_names:
            df = pd.read_excel(xlsx, sheet_name=sheet)
            df = df[~df['Module'].str.contains('grey', case=False)]
            out[sheet] = df
        return out
    except Exception as e:
        print(f"Error loading built-in datasets: {e}")
        return {}

BUILTIN_DATASETS = load_builtin_datasets()
BUILTIN_DATASET_OPTIONS = [{'label': name, 'value': name} for name in BUILTIN_DATASETS]

def create_example_csv():
    if BUILTIN_DATASETS and not os.path.exists('data/example.csv'):
        os.makedirs('data', exist_ok=True)
        first = next(iter(BUILTIN_DATASETS))
        BUILTIN_DATASETS[first].to_csv('data/example.csv', index=False)

create_example_csv()

# ── Layout ────────────────────────────────────────────────────────────────────

app.layout = dbc.Container([

    # App bar
    dbc.Row(dbc.Col(
        html.Div([
            html.H1("Module Gene Overlap Analysis"),
            html.P("Statistical comparison of co-expression modules · Fisher's Exact & Hypergeometric · BH correction"),
        ], id="appbar"),
    width=12), className="mb-4"),

    # Main row
    dbc.Row([

        # ── Left: control panel ───────────────────────────────────────────────
        dbc.Col([
            dbc.Card([
                dbc.CardBody([

                    # 1. Upload
                    html.Div([
                        html.P("Upload Your Data", className="ctrl-label"),
                        dcc.Upload(
                            id='upload-data',
                            children=html.Div([
                                html.Span("☁  ", style={"fontSize": "1.4rem", "verticalAlign": "middle",
                                                         "color": "#00897B"}),
                                html.Span("Drag & drop  or  ", style={"color": "#78909C", "fontSize": ".85rem"}),
                                html.A("browse", style={"color": "#00897B", "fontWeight": "500",
                                                         "fontSize": ".85rem"}),
                                html.Br(),
                                html.Span("Excel (.xlsx) or CSV (.csv)",
                                          style={"fontSize": ".72rem", "color": "#90A4AE"}),
                            ], className="text-center py-3"),
                            className="upload-zone",
                            multiple=False,
                        ),
                        html.Div(id='upload-output', className="mt-2"),
                        html.A(
                            dbc.Button("Download Example File", size="sm", color="light",
                                       style={"color": "#00695C", "fontWeight": "500",
                                              "fontSize": ".78rem", "width": "100%",
                                              "marginTop": "8px", "border": "1px solid #B2DFDB"}),
                            href="/download-example", target="_blank",
                            id="download-example-button",
                            style={"textDecoration": "none"},
                        ),
                    ], className="ctrl-section"),

                    # 2. Reference dataset
                    html.Div([
                        html.P("Reference Dataset", className="ctrl-label"),
                        dcc.Dropdown(
                            id='builtin-dataset-dropdown',
                            options=BUILTIN_DATASET_OPTIONS,
                            value=BUILTIN_DATASET_OPTIONS[0]['value'] if BUILTIN_DATASET_OPTIONS else None,
                            clearable=False,
                        ),
                    ], className="ctrl-section"),

                    # 3. Parameters
                    html.Div([
                        html.P("Analysis Parameters", className="ctrl-label"),
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Universe Size", style={"fontSize": ".78rem", "color": "#546E7A",
                                                                   "fontWeight": "500"}),
                                dbc.Input(id='universe-size', type='number',
                                          value=20000, min=1000, step=1000, size="sm"),
                            ], width=6),
                            dbc.Col([
                                dbc.Label("FDR (α)", style={"fontSize": ".78rem", "color": "#546E7A",
                                                             "fontWeight": "500"}),
                                dbc.Input(id='alpha-threshold', type='number',
                                          value=0.05, min=0.001, max=0.1, step=0.001, size="sm"),
                            ], width=6),
                        ], className="g-2 mb-3"),
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Min. Overlap", style={"fontSize": ".78rem", "color": "#546E7A",
                                                                   "fontWeight": "500"}),
                                dbc.Input(id='min-overlap', type='number',
                                          value=5, min=1, step=1, size="sm"),
                            ], width=6),
                            dbc.Col([
                                dbc.Label("Test Method", style={"fontSize": ".78rem", "color": "#546E7A",
                                                                  "fontWeight": "500"}),
                                dcc.RadioItems(
                                    id='test-type',
                                    options=[
                                        {'label': ' Fisher Exact', 'value': 'fisher'},
                                        {'label': ' Hypergeometric', 'value': 'hypergeometric'},
                                    ],
                                    value='fisher',
                                    inputStyle={"marginRight": "5px"},
                                    labelStyle={"display": "block", "fontSize": ".82rem",
                                                "marginBottom": "4px", "color": "#546E7A"},
                                ),
                            ], width=6),
                        ], className="g-2"),
                        html.P("BH-adjusted p-values will be computed.",
                               style={"fontSize": ".7rem", "color": "#90A4AE", "marginTop": "10px",
                                      "marginBottom": "0"}),
                    ], className="ctrl-section"),

                    # 4. Run button
                    html.Div([
                        dbc.Button("Run Analysis", id="run-analysis"),
                    ], className="ctrl-section"),

                ]),
            ], id="control-card"),
        ], width=4),

        # ── Right: analysis output ────────────────────────────────────────────
        dbc.Col([
            dcc.Loading(
                type="circle", color="#00897B",
                children=html.Div(id='analysis-output'),
            ),
        ], width=8),

    ], className="g-3"),

    # Results table row
    dbc.Row(dbc.Col(
        dcc.Loading(type="circle", color="#00897B",
                    children=html.Div(id='result-tables')),
        width=12,
    ), className="mt-3"),

    # Footer
    html.Hr(style={"borderColor": "#CFD8DC", "marginTop": "3rem"}),
    html.P(
        "Module Gene Overlap Analysis  ·  Fisher's Exact & Hypergeometric  ·  BH correction",
        className="text-center",
        style={"fontSize": ".72rem", "color": "#B0BEC5", "marginBottom": "2rem"},
    ),

    # Stores
    dcc.Store(id='user-dataset'),
    dcc.Store(id='overlap-results'),
    dcc.Store(id='significant-overlaps'),

], fluid=True, style={"paddingLeft": "24px", "paddingRight": "24px"})

# ── Helpers ───────────────────────────────────────────────────────────────────

def parse_upload(contents, filename):
    _, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        if 'xlsx' in filename:
            df = pd.read_excel(io.BytesIO(decoded))
        elif 'csv' in filename:
            df = pd.read_csv(io.BytesIO(decoded))
        else:
            return None, "Please upload an Excel (.xlsx) or CSV (.csv) file."
        if not all(c in df.columns for c in ['Gene', 'Module']):
            return None, "File is missing required columns: Gene, Module."
        df = df[~df['Module'].str.contains('grey', case=False)]
        return df, ""
    except Exception as e:
        return None, f"Parse error: {e}"


def out_card(header, body):
    """Reusable output card with dark header."""
    return html.Div([
        html.Div(header, className="out-card-header"),
        html.Div(body,   className="out-card-body"),
    ], className="out-card")

# ── Callbacks ─────────────────────────────────────────────────────────────────

@app.callback(
    [Output('upload-output', 'children'),
     Output('user-dataset', 'data')],
    Input('upload-data', 'contents'),
    State('upload-data', 'filename'),
)
def cb_upload(contents, filename):
    if contents is None:
        return "", None
    df, err = parse_upload(contents, filename)
    if df is None:
        return dbc.Alert(err, color="danger",
                         style={"fontSize": ".82rem", "padding": "8px 12px",
                                "marginBottom": "0", "borderRadius": "2px"}), None
    n_mod = len(df['Module'].unique())
    return dbc.Alert(
        [html.Strong("✓ "), f"{filename}",
         html.Span(f"  ·  {n_mod} modules  ·  {len(df)} genes",
                   style={"opacity": ".8"})],
        color="success",
        style={"fontSize": ".82rem", "padding": "8px 12px",
               "marginBottom": "0", "borderRadius": "2px"},
    ), df.to_dict('records')


def run_overlap_analysis(user_data, builtin_name, universe_size, min_overlap, test_type):
    if user_data is None:
        return None, "Please upload your data file first."
    if builtin_name not in BUILTIN_DATASETS:
        return None, "Selected reference dataset not found."

    user_df = pd.DataFrame(user_data)
    ref_df  = BUILTIN_DATASETS[builtin_name].copy()
    user_df['Gene'] = user_df['Gene'].str.upper()
    ref_df['Gene']  = ref_df['Gene'].str.upper()

    results = []
    for um in user_df['Module'].unique():
        u_genes = set(user_df[user_df['Module'] == um]['Gene'])
        for rm in ref_df['Module'].unique():
            r_genes = set(ref_df[ref_df['Module'] == rm]['Gene'])
            overlap = u_genes & r_genes
            n_ov = len(overlap)
            if n_ov < min_overlap:
                continue
            nu, nr = len(u_genes), len(r_genes)
            ct = np.array([
                [n_ov,            nu - n_ov],
                [nr - n_ov, universe_size - nu - nr + n_ov],
            ])
            if test_type == 'fisher':
                _, pval = stats.fisher_exact(ct)
            else:
                pval = stats.hypergeom.sf(n_ov - 1, universe_size, nu, nr)
            results.append({
                'UserModule': um, 'RefDataset': builtin_name, 'RefModule': rm,
                'Overlap': n_ov, 'GenesInUserModule': nu, 'GenesInRefModule': nr,
                'P-value': pval, 'OverlappingGenes': list(overlap),
            })

    if not results:
        return None, ("No overlaps found above the minimum threshold. "
                      "Try lowering the minimum overlap count or selecting a different reference.")

    rdf = pd.DataFrame(results)
    rdf['BH-adjusted P-value'] = benjamini_hochberg(rdf['P-value'])
    return rdf, ""


@app.callback(
    [Output('analysis-output', 'children'),
     Output('result-tables', 'children'),
     Output('overlap-results', 'data'),
     Output('significant-overlaps', 'data')],
    Input('run-analysis', 'n_clicks'),
    [State('user-dataset', 'data'),
     State('builtin-dataset-dropdown', 'value'),
     State('universe-size', 'value'),
     State('alpha-threshold', 'value'),
     State('min-overlap', 'value'),
     State('test-type', 'value')],
    prevent_initial_call=True,
)
def perform_analysis(n_clicks, user_data, builtin_name,
                     universe_size, alpha, min_overlap, test_type):
    if not n_clicks:
        return None, None, None, None

    results_df, err = run_overlap_analysis(
        user_data, builtin_name, universe_size, min_overlap, test_type
    )
    if results_df is None:
        return dbc.Alert(err, color="warning",
                         style={"borderRadius": "2px"}), None, None, None

    sig_df = results_df[results_df['BH-adjusted P-value'] < alpha].sort_values('BH-adjusted P-value')

    # ── Heatmap ──────────────────────────────────────────────────────────────
    user_modules = sorted(pd.DataFrame(user_data)['Module'].unique())
    ref_modules  = sorted(BUILTIN_DATASETS[builtin_name]['Module'].unique())
    hm = pd.DataFrame(1.0, index=user_modules, columns=ref_modules)
    for _, row in results_df.iterrows():
        hm.loc[row['UserModule'], row['RefModule']] = row['BH-adjusted P-value']
    hm  = -np.log10(hm)
    hmr = np.round(hm.values, 2)

    keep_r = ~np.all(hmr <= 0.01, axis=1)
    keep_c = ~np.all(hmr <= 0.01, axis=0)
    fdata  = hmr[keep_r][:, keep_c]
    fy     = np.array(user_modules)[keep_r].tolist()
    fx     = np.array(ref_modules)[keep_c].tolist()

    fig = go.Figure(go.Heatmap(
        z=fdata, x=fx, y=fy,
        colorscale=[[0, '#f5f5f5'], [0.3, '#80CBC4'],
                    [0.6, '#00897B'], [0.85, '#E91E63'], [1, '#880E4F']],
        hoverongaps=False,
        hovertemplate='Your Module: %{y}<br>Reference: %{x}<br>-log10(adj.p): %{z:.2f}<extra></extra>',
    ))
    fig.update_layout(
        xaxis_title=f'{builtin_name} Modules',
        yaxis_title='Your Data Modules',
        height=800, autosize=True,
        margin=dict(t=50, b=70, l=160, r=40),
        xaxis={'side': 'top'},
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        font=dict(family='Roboto, sans-serif', color='#37474F', size=12),
    )
    fig.update_xaxes(tickangle=45)
    fig.add_annotation(
        text="-log10(BH-adjusted p-value)  ·  brighter = more significant overlap",
        xref="paper", yref="paper", x=0.5, y=-0.09,
        showarrow=False, font=dict(size=11, color='#90A4AE'),
    )

    # ── Stat chips ───────────────────────────────────────────────────────────
    n_sig = len(sig_df)

    def chip(val, lbl, color):
        return html.Div([
            html.Span(str(val), className="num", style={"color": color}),
            html.Span(lbl,      className="lbl"),
        ], className="stat-chip")

    stat_row = html.Div([
        chip(len(results_df), "Pairs Tested",             "#00897B"),
        chip(n_sig,           f"Significant (α={alpha})", "#E91E63" if n_sig else "#90A4AE"),
        chip(len(fy),         "Modules in Heatmap",       "#3949AB"),
    ], className="stat-chips-row")

    heatmap_div = html.Div([
        stat_row,
        out_card(
            f"Overlap Heatmap  ·  Your Data  ×  {builtin_name}",
            dcc.Graph(figure=fig, config={'displayModeBar': True, 'displaylogo': False},
                      style={"padding": "8px"}),
        ),
    ])

    # ── Significant results table ─────────────────────────────────────────────
    if n_sig > 0:
        tdf = sig_df.copy()
        tdf['Top Overlapping Genes'] = tdf['OverlappingGenes'].apply(
            lambda x: ', '.join(x) if len(x) <= 5
                      else ', '.join(x[:5]) + f'  … +{len(x) - 5}'
        )
        tdata = tdf.drop(columns=['OverlappingGenes']).to_dict('records')

        sig_table = out_card(
            [html.Span("Significant Module Overlaps"),
             html.Span(f" ({n_sig})",
                       style={"background": "#E91E63", "borderRadius": "10px",
                              "padding": "1px 8px", "fontSize": ".7rem",
                              "marginLeft": "8px", "verticalAlign": "middle"})],
            html.Div([
                html.P("Select a row to see overlapping genes.",
                       style={"fontSize": ".78rem", "color": "#90A4AE",
                              "padding": "12px 16px 4px", "margin": 0}),
                dash_table.DataTable(
                    id='significant-table',
                    columns=[
                        {"name": "Your Module",          "id": "UserModule"},
                        {"name": "Reference Module",     "id": "RefModule"},
                        {"name": "Overlap",              "id": "Overlap"},
                        {"name": "Genes (Yours)",        "id": "GenesInUserModule"},
                        {"name": "Genes (Ref)",          "id": "GenesInRefModule"},
                        {"name": "P-value",              "id": "P-value",
                         "type": "numeric", "format": {"specifier": ".2e"}},
                        {"name": "BH adj. P-value",      "id": "BH-adjusted P-value",
                         "type": "numeric", "format": {"specifier": ".2e"}},
                        {"name": "Top Overlapping Genes","id": "Top Overlapping Genes"},
                    ],
                    data=tdata,
                    sort_action="native",
                    sort_by=[{"column_id": "BH-adjusted P-value", "direction": "asc"}],
                    page_size=15,
                    style_table={"overflowX": "auto"},
                    style_cell={
                        "textAlign": "left",
                        "padding": "10px 14px",
                        "fontFamily": "Roboto, sans-serif",
                        "fontSize": "13px",
                        "border": "none",
                        "borderBottom": "1px solid #ECEFF1",
                        "color": "#37474F",
                    },
                    style_header={
                        "backgroundColor": "#37474F",
                        "color": "#fff",
                        "fontWeight": "500",
                        "fontSize": "11px",
                        "letterSpacing": ".06em",
                        "textTransform": "uppercase",
                        "border": "none",
                        "padding": "12px 14px",
                    },
                    style_data_conditional=[
                        {"if": {"row_index": "odd"}, "backgroundColor": "#FAFAFA"},
                        {"if": {"state": "selected"},
                         "backgroundColor": "#E0F2F1",
                         "border": "1px solid #00897B !important"},
                        {"if": {"column_id": "BH-adjusted P-value",
                                "filter_query": "{BH-adjusted P-value} < 0.001"},
                         "color": "#E91E63", "fontWeight": "500"},
                    ],
                    row_selectable="single",
                ),
                html.Div(id='selected-overlap-output',
                         style={"padding": "0 16px"}),
                html.Div([
                    dbc.Button("Download Gene List", id="download-genes",
                               className="download-btn mt-3",
                               style={"display": "none"}),
                    dcc.Download(id="download-gene-list"),
                ], style={"padding": "0 16px 16px"}),
            ]),
        )
    else:
        sig_table = out_card(
            "Significant Module Overlaps",
            html.Div([
                html.Div("🔍", style={"fontSize": "2rem"}),
                html.P("No significant overlaps found.",
                       style={"fontWeight": "500", "color": "#546E7A", "margin": "8px 0 4px"}),
                html.P("Try adjusting the FDR threshold or minimum overlap count.",
                       style={"fontSize": ".83rem", "color": "#90A4AE", "margin": 0}),
            ], className="text-center py-5"),
        )

    results_json = results_df.to_json(orient='records')
    sig_json     = sig_df.to_json(orient='records') if n_sig > 0 else None
    return heatmap_div, sig_table, results_json, sig_json


@app.callback(
    [Output('selected-overlap-output', 'children'),
     Output('download-genes', 'style')],
    Input('significant-table', 'selected_rows'),
    State('significant-overlaps', 'data'),
)
def show_selected_overlap(selected_rows, significant_data):
    hidden  = {"display": "none"}
    visible = {"display": "inline-block", "background": "#00897B", "border": "none",
               "borderRadius": "2px", "fontSize": ".75rem", "fontWeight": "500",
               "letterSpacing": ".07em", "textTransform": "uppercase"}

    if not selected_rows or significant_data is None:
        return None, hidden

    try:
        sig_df = pd.read_json(significant_data, orient='records')
    except Exception:
        sig_df = pd.DataFrame(significant_data)

    row   = sig_df.iloc[selected_rows[0]]
    genes = sorted(row['OverlappingGenes'])
    shown = genes[:80]

    chips = [html.Span(g, className="gene-tag") for g in shown]
    if len(genes) > 80:
        chips.append(html.Span(f"  +{len(genes) - 80} more",
                               style={"fontSize": ".75rem", "color": "#90A4AE"}))

    panel = html.Div([
        html.Div([
            html.Span(f"{row['UserModule']}  ∩  {row['RefModule']}",
                      style={"fontWeight": "500", "color": "#37474F", "fontSize": ".88rem"}),
            html.Span(f"  ·  {len(genes)} genes",
                      style={"fontSize": ".78rem", "color": "#78909C"}),
        ], style={"padding": "10px 14px", "background": "#E0F2F1",
                  "borderBottom": "1px solid #B2DFDB",
                  "borderRadius": "2px 2px 0 0"}),
        html.Div(chips, style={"padding": "12px 14px", "maxHeight": "180px",
                                "overflowY": "auto", "lineHeight": "2.2"}),
    ], style={"border": "1px solid #B2DFDB", "borderRadius": "2px",
              "marginTop": "12px", "background": "#fff"})

    return panel, visible


@app.callback(
    Output("download-gene-list", "data"),
    Input("download-genes", "n_clicks"),
    [State('significant-table', 'selected_rows'),
     State('significant-overlaps', 'data')],
    prevent_initial_call=True,
)
def download_gene_list(n_clicks, selected_rows, significant_data):
    if not n_clicks or not selected_rows or significant_data is None:
        return None
    try:
        sig_df = pd.read_json(significant_data, orient='records')
    except Exception:
        sig_df = pd.DataFrame(significant_data)
    row  = sig_df.iloc[selected_rows[0]]
    text = '\n'.join(sorted(row['OverlappingGenes']))
    return dict(
        content=text,
        filename=(f"overlapping_genes_YourData_{row['UserModule']}_"
                  f"{row['RefDataset']}_{row['RefModule']}.txt"),
    )

# ── Flask route ───────────────────────────────────────────────────────────────

@server.route('/download-example')
def download_example():
    try:
        return send_file('data/example.csv', mimetype='text/csv',
                         as_attachment=True, download_name='example.csv')
    except Exception as e:
        return f"Error: {e}", 404


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=port, debug=False)
