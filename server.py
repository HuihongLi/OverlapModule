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

# Benjamini-Hochberg correction implementation
def benjamini_hochberg(p_values):
    p_values = np.array(p_values)
    n = len(p_values)
    ranked_indices = np.argsort(p_values)
    sorted_p = p_values[ranked_indices]
    adjusted_p = np.empty_like(p_values)
    
    for i, idx in enumerate(ranked_indices):
        adjusted_p[idx] = min(1, sorted_p[i] * n / (i + 1))
    
    # Ensure monotonicity
    for i in range(n-2, -1, -1):
        adjusted_p[ranked_indices[i]] = min(adjusted_p[ranked_indices[i]], 
                                          adjusted_p[ranked_indices[i+1]])
    
    return adjusted_p

# Initialize Dash app with Bootstrap
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server
app.config.suppress_callback_exceptions = True

# Load built-in datasets
def load_builtin_datasets():
    try:
        excel_file = "data/data.xlsx"
        builtin_datasets = {}
        
        xlsx = pd.ExcelFile(excel_file)
        for sheet_name in xlsx.sheet_names:
            df = pd.read_excel(xlsx, sheet_name=sheet_name)
            df = df[~df['Module'].str.contains('grey', case=False)]
            builtin_datasets[sheet_name] = df
            
        return builtin_datasets
    except Exception as e:
        print(f"Error loading built-in datasets: {str(e)}")
        return {}

# Load datasets and set options
BUILTIN_DATASETS = load_builtin_datasets()
BUILTIN_DATASET_OPTIONS = [{'label': name, 'value': name} for name in BUILTIN_DATASETS.keys()]

# Create example CSV if needed
def create_example_csv():
    if BUILTIN_DATASETS and not os.path.exists('data/example.csv'):
        first_dataset_name = next(iter(BUILTIN_DATASETS))
        BUILTIN_DATASETS[first_dataset_name].to_csv('data/example.csv', index=False)
        print("Created data/example.csv from built-in dataset")

create_example_csv()  # Create example during startup

# App layout
app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H1("Module Gene Overlap Analysis", className="text-center mb-4"),
            html.Hr(),
        ], width=12)
    ]),
    
    dbc.Row([
        dbc.Col([
            html.H4("Upload Your Data", className="mb-2"),
            dbc.Row([
                dbc.Col([
                    dcc.Upload(
                        id='upload-data',
                        children=html.Div(['Drag and Drop or ', html.A('Select Excel/CSV File')]),
                        style={
                            'width': '100%', 'height': '60px', 'lineHeight': '60px',
                            'borderWidth': '1px', 'borderStyle': 'dashed',
                            'borderRadius': '5px', 'textAlign': 'center',
                        },
                        multiple=False
                    ),
                ], width=9, className="pe-0"),
                dbc.Col([
                    html.Div([
                        html.A(
                            dbc.Button("Download Example", color="secondary", className="w-100", style={'height': '60px'}),
                            href="/download-example",
                            target="_blank",
                            id="download-example-button"
                        ),
                    ], style={'height': '100%'}, className="d-flex align-items-center"),
                ], width=3, className="ps-0"),
            ], className="align-items-center mb-3"),
            html.Div(id='upload-output'),
            
            html.Hr(),
            
            html.H4("Reference Dataset Selection", className="mb-2"),
            html.P("Select a built-in dataset to compare with your data:"),
            dcc.Dropdown(
                id='builtin-dataset-dropdown',
                options=BUILTIN_DATASET_OPTIONS,
                value=BUILTIN_DATASET_OPTIONS[0]['value'] if BUILTIN_DATASET_OPTIONS else None,
                clearable=False
            ),
            
            html.Hr(),
            
            html.H4("Analysis Parameters", className="mb-2"),
            dbc.Row([
                dbc.Col([
                    html.P("Total gene universe size:"),
                    dbc.Input(id='universe-size', type='number', value=20000, min=1000, step=1000),
                ], width=6),
                dbc.Col([
                    html.P("Significance threshold (FDR):"),
                    dbc.Input(id='alpha-threshold', type='number', value=0.05, min=0.001, max=0.1, step=0.001),
                    html.Small("Benjamini-Hochberg adjusted p-values will be used", className="text-muted"),
                ], width=6),
            ]),
            
            dbc.Row([
                dbc.Col([
                    html.P("Minimum overlap count:"),
                    dbc.Input(id='min-overlap', type='number', value=5, min=1, step=1),
                ], width=6),
                dbc.Col([
                    html.P("Statistical test:"),
                    dcc.RadioItems(
                        id='test-type',
                        options=[
                            {'label': 'Fisher Exact Test', 'value': 'fisher'},
                            {'label': 'Hypergeometric Test', 'value': 'hypergeometric'}
                        ],
                        value='fisher',
                        inline=True
                    ),
                ], width=6),
            ]),
            
            html.Br(),
            dbc.Button("Run Analysis", id="run-analysis", color="primary", className="mt-2"),
            
        ], width=4),
        
        dbc.Col([
            html.Div(id='analysis-output')
        ], width=8),
    ]),
    
    dbc.Row([
        dbc.Col([
            html.Div(id='result-tables'),
        ], width=12),
    ]),
    
    # Store components
    dcc.Store(id='user-dataset'),
    dcc.Store(id='overlap-results'),
    dcc.Store(id='significant-overlaps')
], fluid=True)

# File upload handler
def parse_upload(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    
    try:
        if 'xlsx' in filename:
            df = pd.read_excel(io.BytesIO(decoded))
        elif 'csv' in filename:
            df = pd.read_csv(io.BytesIO(decoded))
        else:
            return None, "Please upload an Excel (.xlsx) or CSV (.csv) file"
        
        if not all(col in df.columns for col in ['Gene', 'Module']):
            return None, "Uploaded file is missing required columns (Gene, Module)"
        
        df = df[~df['Module'].str.contains('grey', case=False)]
        
        return df, ""
    except Exception as e:
        return None, f"Error parsing file: {str(e)}"

@app.callback(
    [Output('upload-output', 'children'),
     Output('user-dataset', 'data')],
    [Input('upload-data', 'contents')],
    [State('upload-data', 'filename')]
)
def update_output(contents, filename):
    if contents is None:
        return "", None
    
    user_df, error_message = parse_upload(contents, filename)
    
    if user_df is None:
        return html.Div([error_message], style={'color': 'red'}), None
    
    return html.Div([
        html.P(f"Successfully uploaded: {filename}"),
        html.P(f"Found {len(user_df['Module'].unique())} modules and {len(user_df)} gene entries")
    ], style={'color': 'green'}), user_df.to_dict('records')

# Run analysis between datasets
def run_overlap_analysis(user_data, builtin_dataset_name, universe_size, min_overlap, test_type):
    if user_data is None:
        return None, "Please upload your data file first."
    
    if builtin_dataset_name not in BUILTIN_DATASETS:
        return None, "Selected reference dataset not found."
    
    # Convert user data to DataFrame
    user_df = pd.DataFrame(user_data)
    ref_df = BUILTIN_DATASETS[builtin_dataset_name]
    
    # Get unique modules and convert genes to uppercase
    user_modules = user_df['Module'].unique().tolist()
    ref_modules = ref_df['Module'].unique().tolist()
    user_df['Gene'] = user_df['Gene'].str.upper()
    ref_df['Gene'] = ref_df['Gene'].str.upper()
    
    results = []
    
    # Compare modules
    for user_module in user_modules:
        user_genes = set(user_df[user_df['Module'] == user_module]['Gene'].tolist())
        
        for ref_module in ref_modules:
            ref_genes = set(ref_df[ref_df['Module'] == ref_module]['Gene'].tolist())
            
            # Calculate overlap
            overlap = user_genes.intersection(ref_genes)
            overlap_count = len(overlap)
            
            if overlap_count < min_overlap:
                continue
            
            # Statistical test setup
            n_user_genes = len(user_genes)
            n_ref_genes = len(ref_genes)
            
            contingency_table = np.array([
                [overlap_count, n_user_genes - overlap_count],
                [n_ref_genes - overlap_count, universe_size - n_user_genes - n_ref_genes + overlap_count]
            ])
            
            # Run appropriate test
            if test_type == 'fisher':
                odds_ratio, p_value = stats.fisher_exact(contingency_table)
            else:  # hypergeometric
                p_value = stats.hypergeom.sf(
                    overlap_count - 1,
                    universe_size,
                    n_user_genes,
                    n_ref_genes
                )
            
            # Store result
            results.append({
                'UserModule': user_module,
                'RefDataset': builtin_dataset_name,
                'RefModule': ref_module,
                'Overlap': overlap_count,
                'GenesInUserModule': n_user_genes,
                'GenesInRefModule': n_ref_genes,
                'P-value': p_value,
                'OverlappingGenes': list(overlap)
            })
    
    if not results:
        return None, "No significant overlaps found. Try reducing the minimum overlap threshold or selecting a different reference dataset."
    
    # Convert to DataFrame and add adjustments
    results_df = pd.DataFrame(results)
    results_df['BH-adjusted P-value'] = benjamini_hochberg(results_df['P-value'])
    
    return results_df, ""

@app.callback(
    [Output('analysis-output', 'children'),
     Output('result-tables', 'children'),
     Output('overlap-results', 'data'),
     Output('significant-overlaps', 'data')],
    [Input('run-analysis', 'n_clicks')],
    [State('user-dataset', 'data'),
     State('builtin-dataset-dropdown', 'value'),
     State('universe-size', 'value'),
     State('alpha-threshold', 'value'),
     State('min-overlap', 'value'),
     State('test-type', 'value')],
    prevent_initial_call=True
)
def perform_analysis(n_clicks, user_data, builtin_dataset_name, universe_size, alpha, min_overlap, test_type):
    if n_clicks is None or not n_clicks or user_data is None or not builtin_dataset_name:
        return None, None, None, None
    
    # Run the analysis
    results_df, error_message = run_overlap_analysis(user_data, builtin_dataset_name, universe_size, min_overlap, test_type)
    
    if results_df is None:
        return html.Div([error_message], style={'color': 'red'}), None, None, None
    
    # Filter significant results
    significant_df = results_df[results_df['BH-adjusted P-value'] < alpha].sort_values('BH-adjusted P-value')
    
    # Create heatmap data
    user_df = pd.DataFrame(user_data)
    user_modules = sorted(user_df['Module'].unique().tolist())
    ref_df = BUILTIN_DATASETS[builtin_dataset_name]
    ref_modules = sorted(ref_df['Module'].unique().tolist())
    
    # Initialize heatmap with high p-values (1.0 = not significant)
    heatmap_data = pd.DataFrame(1.0, index=user_modules, columns=ref_modules)
    
    # Fill in adjusted p-values from results
    for _, row in results_df.iterrows():
        heatmap_data.loc[row['UserModule'], row['RefModule']] = row['BH-adjusted P-value']
    
    # Apply -log10 transformation and round values
    heatmap_data = -np.log10(heatmap_data)
    heatmap_data_rounded = np.round(heatmap_data.values, 2)
    
    # Remove empty rows and columns
    non_empty_rows = ~np.all(heatmap_data_rounded <= 0.01, axis=1)
    non_empty_cols = ~np.all(heatmap_data_rounded <= 0.01, axis=0)
    
    # Filter data for non-empty rows/columns
    filtered_data = heatmap_data_rounded[non_empty_rows][:, non_empty_cols]
    filtered_y = np.array(user_modules)[non_empty_rows].tolist()
    filtered_x = np.array(ref_modules)[non_empty_cols].tolist()
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=filtered_data,
        x=filtered_x,
        y=filtered_y,
        colorscale='YlOrRd',
        hoverongaps=False,
        hovertemplate='Module: %{y}<br>Reference: %{x}<br>-log10(adj.p): %{z:.2f}<extra></extra>'
    ))
    
    fig.update_layout(
        xaxis_title=f'{builtin_dataset_name} Modules',
        yaxis_title='Your Data Modules',
        autosize=True,
        height=800,
        width=1000,
        margin=dict(t=100, b=50, l=150, r=50),
        xaxis={'side': 'top'},
    )
    
    fig.update_xaxes(tickangle=45)
    fig.add_annotation(
        text="-log10(adjusted p-value) scale: darker = more significant overlap",
        xref="paper", yref="paper",
        x=0.5, y=-0.07,
        showarrow=False
    )
    
    heatmap_div = html.Div([
        html.H4(f"Overlap Analysis: Your Data vs {builtin_dataset_name}", className="text-center mt-4"),
        dcc.Graph(figure=fig)
    ])
    
    # Create DataTable for significant results
    if len(significant_df) > 0:
        table_df = significant_df.copy()
        table_df['OverlappingGenes_str'] = table_df['OverlappingGenes'].apply(
            lambda x: ', '.join(x) if len(x) <= 5 else ', '.join(x[:5]) + '...'
        )
        
        table_data = table_df.drop(columns=['OverlappingGenes']).rename(
            columns={'OverlappingGenes_str': 'Top Overlapping Genes'}
        )
        
        significant_table = html.Div([
            html.H4("Significant Module Overlaps", className="text-center mt-4"),
            dash_table.DataTable(
                id='significant-table',
                columns=[
                    {"name": "Your Module", "id": "UserModule"},
                    {"name": "Reference Module", "id": "RefModule"},
                    {"name": "Overlap Count", "id": "Overlap"},
                    {"name": "Genes in Your Module", "id": "GenesInUserModule"},
                    {"name": "Genes in Ref Module", "id": "GenesInRefModule"},
                    {"name": "P-value", "id": "P-value", "type": "numeric", "format": {"specifier": ".2e"}},
                    {"name": "BH-adjusted P-value", "id": "BH-adjusted P-value", "type": "numeric", "format": {"specifier": ".2e"}},
                    {"name": "Top Overlapping Genes", "id": "Top Overlapping Genes"},
                ],
                data=table_data.to_dict('records'),
                sort_action="native",
                sort_by=[{"column_id": "BH-adjusted P-value", "direction": "asc"}],
                page_size=15,
                style_table={'overflowX': 'auto'},
                style_cell={'textAlign': 'left'},
                style_header={
                    'backgroundColor': 'rgb(230, 230, 230)',
                    'fontWeight': 'bold'
                },
                style_data_conditional=[
                    {
                        'if': {'row_index': 'odd'},
                        'backgroundColor': 'rgb(248, 248, 248)'
                    }
                ],
                row_selectable="single"
            ),
            html.Div(id='selected-overlap-output'),
            dbc.Button(
                "Download Gene List", 
                id="download-genes", 
                color="info", 
                className="mt-3",
            ),
            dcc.Download(id="download-gene-list")
        ])
    else:
        significant_table = html.Div([
            html.H4("No Significant Module Overlaps Found", className="text-center mt-4"),
            html.P("Try adjusting the significance threshold or minimum overlap count.")
        ])
    
    # Store as JSON for easier handling
    results_json = results_df.to_json(orient='records')
    significant_json = significant_df.to_json(orient='records') if len(significant_df) > 0 else None
    
    return heatmap_div, significant_table, results_json, significant_json

@app.callback(
    [Output('selected-overlap-output', 'children'),
     Output('download-genes', 'style')],
    [Input('significant-table', 'selected_rows')],
    [State('significant-overlaps', 'data')]
)
def show_selected_overlap(selected_rows, significant_data):
    if not selected_rows or significant_data is None:
        return None, {'display': 'none'}

    try:
        significant_df = pd.read_json(significant_data, orient='records')
    except:
        significant_df = pd.DataFrame(significant_data)

    row_idx = selected_rows[0]
    selected_overlap = significant_df.iloc[row_idx]
    overlapping_genes = selected_overlap['OverlappingGenes']

    gene_list = [html.Li(gene) for i, gene in enumerate(sorted(overlapping_genes)) if i < 30]

    return (
        html.Div([
            html.H5(f"Overlapping Genes: {selected_overlap['UserModule']} (Your Data) ∩ "
                    f"{selected_overlap['RefModule']}", className="mt-3"),
            html.Div([
                html.P(f"Number of overlapping genes: {len(overlapping_genes)}"),
                html.Div([
                    html.Ul(gene_list, style={'columnCount': '3', 'columnGap': '20px'})
                ]),
                html.P("..." if len(overlapping_genes) > 30 else "")
            ], style={'maxHeight': '300px', 'overflowY': 'auto'})
        ]),
        {'display': 'inline-block'}  # show button
    )

@app.callback(
    Output("download-gene-list", "data"),
    Input("download-genes", "n_clicks"),
    [State('significant-table', 'selected_rows'),
     State('significant-overlaps', 'data')],
    prevent_initial_call=True
)
def download_gene_list(n_clicks, selected_rows, significant_data):
    if not n_clicks or not selected_rows or significant_data is None:
        return None
    
    try:
        significant_df = pd.read_json(significant_data, orient='records')
    except:
        significant_df = pd.DataFrame(significant_data)
    
    row_idx = selected_rows[0]
    selected_overlap = significant_df.iloc[row_idx]
    overlapping_genes = selected_overlap['OverlappingGenes']
    gene_text = '\n'.join(sorted(overlapping_genes))
    
    return dict(
        content=gene_text,
        filename=f"overlapping_genes_YourData_{selected_overlap['UserModule']}_"
                 f"{selected_overlap['RefDataset']}_{selected_overlap['RefModule']}.txt"
    )
# Get port from environment variable (Render sets this)
port = int(os.environ.get("PORT", 10000))

app = Dash(__name__)
server = app.server

# Server route for example file download
@server.route('/download-example')
def download_example():
    try:
        return send_file('data/example.csv', 
                        mimetype='text/csv',
                        as_attachment=True,
                        download_name='data/example.csv')
    except Exception as e:
        print(f"Error serving data/example.csv: {str(e)}")
        return "Error: Example file not available", 404

if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=port, debug=False)