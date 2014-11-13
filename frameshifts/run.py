import plotly.plotly as py
import sys
from plotly.graph_objs import *

def main():

    tc = {
    "mafft" : [0.958,   0.958,   0.938,   0.942,   0.945,   0.961,   0.934,   0.956,   0.933,   0.943],
    "clustalo" : [0.754,   0.794,   0.721,   0.847,   0.776,   0.737,   0.815,   0.815,   0.758,   0.733],
    "macse" : [0.94,    0.977,   0.974,   0.978,   0.979,   0.93,    0.975,   0.919,   0.977,   0.922],
    "noah" : [0.922,   0.908,   0.917,   0.918,   0.906,   0.92,    0.918,   0.921,   0.91,    0.913]
    }


    tc["mafft"].extend([0.64,    0.607,   0.551,   0.563,   0.617,   0.584,   0.61,    0.593,   0.628,   0.596])
    tc["clustalo"].extend([0.0484,  0.0702,  0.0906,  0.0915,  0.113,   0.189,   0.107,   0.0421,  0.0386,  0.0842])
    tc["macse"].extend([0.318,   0.14,    0.178,   0.162,   0.16,    0.196,   0.117,   0.133,   0.189,   0.249])
    tc["noah"].extend([0.398,   0.411,   0.39,    0.356,   0.319,   0.402,   0.366,   0.361,   0.358,   0.407])

    labels = ["SP"] * 10
    labels.extend(["TC"] * 10)

    traces = []
    for name, vals in tc.iteritems():

        traces.append(Box(
                    y=vals,
                    x=labels,
                    name=name,
                    opacity=0.8 
                )
            )
        layout = Layout(
            showlegend=True,
            boxmode="group",
            title="Sum of pairs and total column scores",
            titlefont = Font(
                size = 20),
            yaxis=YAxis(
                title='Score',
                showgrid=False,
                autorange = True,
                titlefont = Font(
                    size = 20)
            ),
            xaxis=XAxis(
                showticklabels=False),
            legend = Legend(
                font = Font(
                    size = 18
                ),
                yanchor = "bottom",
                xanchor = "left")
        )
    data = Data(traces)
    fig = Figure(data=data, layout=layout)
    py.sign_in("imogen", "mtf1tawct2")
    py.plot(fig)

if __name__ == "__main__":
    main()