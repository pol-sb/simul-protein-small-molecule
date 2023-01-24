from dash import Dash, dcc, html, Input, Output
import dash_bootstrap_components as dbc

import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np

import pandas as pd

app = Dash(__name__)


SUBP_TITL = ("IDP Condensate Average", "DRG Condensate Average")

df = pd.read_pickle("results/simulations_df.pkl")


colors = {
    "background": "#fbf1c7",
    "bg0": "#fbf1c7",
    "bg1": "#ebdbb2",
    "bg2": "#d5c4a1",
    "bg3": "#bdae93",
    "bg4": "#a89984",
    "fg0": "#282828",
    "fg1": "#3c3836",
    "fg2": "#504945",
    "fg3": "#665c54",
    "fg4": "#7c6f64",
    "gray1": "#928374",
    "blue1": "#458588",
    "blue2": "#076678",
    "yellow1": "#d79921",
    "yellow2": "#b57614",
    "green1": "#98971a",
    "green2": "#79740e",
    "aqua1": "#689d6a",
    "aqua2": "#427b58",
    "purple1": "#b16286",
    "purple2": "#8f3f71",
    "orange1": "#d65d0e",
    "orange2": "#af3a03",
    "red1": "#cc241d",
    "red2": "#9d0006",
    "text": "#111111",
}

fig = make_subplots(rows=1, cols=2, subplot_titles=SUBP_TITL)
fig.add_trace(
    row=1,
    col=1,
    trace=go.Scatter(
        x=df["idp_plat_avg"],
        y=df["temp"],
        mode="markers",
        name="Average IDP Density",
    ),
)

fig.add_trace(
    row=1,
    col=2,
    trace=go.Scatter(
        x=df["drg_plat_avg"], y=df["temp"], mode="markers", name="Average DRG Density"
    ),
)

fig.update_traces(yaxis="y1")

fig.update_layout(
    plot_bgcolor=colors["background"],
    paper_bgcolor=colors["background"],
    font_color=colors["text"],
    hovermode="y unified",
    hoverlabel=dict(bgcolor="white", font_size=16, font_family="Sans Serif"),
    showlegend=False,
)

app.layout = html.Div(
    style={"backgroundColor": colors["bg1"]},
    children=[
        html.Div(
            [
                html.Div(
                    [
                        html.H1(
                            children="IDP Condensates",
                            style={
                                "textAlign": "center",
                                "backgroundcolor": colors["background"],
                                "color": colors["text"],
                            },
                        ),
                    ]
                ),
                html.Div(
                    [
                        html.Div(
                            children="Testing dash web apps using simulation data for MD simulations of IDP condensates.",
                            style={
                                "color": colors["text"],
                                "textAlign": "center",
                            },
                        ),
                    ]
                ),
                html.Div(
                    [
                        html.Div(
                            dcc.Graph(id="idp-graph", figure=fig),
                            style={
                                "display": "inline-block",
                                "width": "95%",
                                "marginLeft": "20px",
                            },
                        ),
                        html.Div(
                            children=[
                                html.Div(
                                    [
                                        "Lambda:",
                                        dcc.Dropdown(
                                            df["lambd"].unique(),
                                            0.7,
                                            placeholder="Lambda Value",
                                            id="xaxis-column",
                                        ),
                                    ],
                                    style={
                                        "width": "30%",
                                        "verticalAlign": "baseline",
                                    },
                                ),
                                html.Div(
                                    [
                                        "Concentration (mM):",
                                        dcc.Dropdown(
                                            df["conc"].unique(),
                                            20,
                                            placeholder="Concentration Value",
                                            id="concval",
                                        ),
                                    ],
                                    style={
                                        "width": "30%",
                                        "verticalAlign": "baseline",
                                    },
                                ),
                                html.Div(
                                    [
                                        "Sigma:",
                                        dcc.Dropdown(
                                            df["sigma"].unique(),
                                            0.45,
                                            placeholder="Sigma Value",
                                            id="sigmaval",
                                        ),
                                    ],
                                    style={
                                        "width": "30%",
                                        "verticalAlign": "baseline",
                                    },
                                ),
                            ],
                            style={
                                "display": "inline-block",
                                "width": "40%",
                                "margin-left": "50px",
                                "marginTop": "50px",
                            },
                        ),
                    ],
                    style={
                        "display": "flex",
                        "width": "100%",
                    },
                ),
            ]
        ),
    ],
)


@app.callback(
    Output("idp-graph", "figure"),
    Input("xaxis-column", "value"),
    Input("concval", "value"),
    Input("sigmaval", "value"),
)
def update_graph(xaxis_column_name, concval_name, sigmaval_name):
    dff = df[df["lambd"] == xaxis_column_name]
    dff = dff[dff["conc"] == concval_name]
    dff = dff[dff["sigma"] == sigmaval_name]

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=SUBP_TITL,
    )
    fig.add_trace(
        row=1,
        col=1,
        trace=go.Scatter(
            x=dff["idp_plat_avg"],
            y=dff["temp"],
            mode="markers",
            name="Average IDP Density",
            marker=dict(
                autocolorscale=True,
                colorscale="viridis",
                symbol="square-dot",
                color=dff["idp_plat_avg"],
                cauto=True,
                colorbar=dict(),
            ),
            customdata=np.stack((dff["conc"], dff["sigma"]), axis=-1),
            hovertemplate="<br>T(K) = %{y}<br>IDP avg. conc = %{x}<br>conc = %{customdata[0]} <br>sig = %{customdata[1]}",
        ),
    )

    fig.add_trace(
        row=1,
        col=2,
        trace=go.Scatter(
            x=dff["drg_plat_avg"],
            y=dff["temp"],
            mode="markers",
            marker=dict(
                autocolorscale=True,
                colorscale="viridis",
                symbol="square-dot",
                color=dff["drg_plat_avg"],
                cauto=True,
                colorbar=dict(),
            ),
            name="Average DRG Density",
        ),
    )

    fig.update_traces(yaxis="y1")
    fig.update_layout(
        margin={"l": 40, "b": 40, "t": 40, "r": 0},
        hovermode="y unified",
        hoverlabel=dict(
            bgcolor="white", font_size=16, font_family="Computer Modern Roman"
        ),
        showlegend=False,
    )

    return fig


if __name__ == "__main__":
    app.run_server(debug=True, host="0.0.0.0")
