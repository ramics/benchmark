import os
import pickle
import sys
from scipy.stats.mstats import rankdata
from scipy import stats
import plotly.plotly as py
from plotly.graph_objs import *
import numpy
import pprint as pp

def main(folder):

    read_scores = {}
    read_ranks = {}
    read_lengths = {}
    read_info_content = {}

    final_scores = {}
    final_ranks = {}

    num_tests = {}

    filenames = []

    num_ranks = 0

    # get general info
    tool_0 = os.listdir(folder)[0]
    pkl_file = open(os.path.join(folder, tool_0), 'rb')
    (reads, results) = pickle.load(pkl_file)
    for filename, values in reads.iteritems():
        read_scores[filename] = {}
        read_ranks[filename] = {}
        filenames.append(filename)

    tools = os.listdir(folder)

    # get specific info
    for tool in tools:
        pkl_file = open(os.path.join(folder, tool), 'rb')
        (reads, results) = pickle.load(pkl_file)
        final_scores[tool] = [0] * 2
        final_ranks[tool] = [0] * 2

        num_tests[tool] = len(results)
        for filename, values in results.iteritems():

            read_scores[filename][tool] = values["scores"]
            for i in range(len(values["scores"])):
                final_scores[tool][i] += values["scores"][i]


    #for filename, tool_scores in read_scores.iteritems():
    #    print filename, tool_scores 

    for filename, tool_scores in read_scores.iteritems():

        # We only rank on the ones that didn't error for anyone
        if len(tool_scores) < len(tools):
            continue
        num_ranks += 1

        raw_scores = [(tool, scores) for tool, scores in tool_scores.iteritems()]
        ranks = [
            rankdata([1 - r[1][score] for r in raw_scores]
            ) for score in range(2)]
        for p, (tool, _) in enumerate(raw_scores):
            for score_type in range(2):    
                final_ranks[tool][score_type] += ranks[score_type][p]

    results = []
    
    print "Tool\ttests\tsp\ttc\tsp\ttc"
    for tool in tools:
        results = [tool, num_tests[tool]]
        results.extend([s/float(num_tests[tool]) for s in final_scores[tool]])
        results.extend([s/float(num_ranks) for s in final_ranks[tool]])
        print "%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f" % (
                tuple(results) 
        )

    traces = []
    tool_names = {"clustalo.sh":"Clustal Omega",
        "clustalw.sh":"Clustal W",
        "muscle.sh":"MUSCLE",
        "mafft.sh":"MAFFT L-INS-I",
        "prank.sh":"PRANK",
        "prank_f.sh":"PRANK+F",
        "noah_full.sh":"NOAH (full algorithm)",
        "noah_basic.sh":"NOAH (chained guide tree, no FB traversal)",
        "noah_no_fb.sh":"NOAH (guide arc, no FB traversal)",
        "noah_no_arc.sh":"NOAH (chained guide tree, FB traversal)"
    }

    for tool in tools:
        chosen_files = [f for f in filenames if tool in read_scores[f].keys()]
        qscore_data = [read_scores[f][tool][0] for f in chosen_files]
        #tc_data = [read_scores[f][tool][1] for f in chosen_files]
        labels = ["Sum of pairs" for _ in range(len(qscore_data))]
        #labels.extend(["Total column" for _ in range(len(tc_data))])
        all_data = qscore_data
        
        traces.append(Box(
                y=all_data,
                name=tool_names[tool],
                opacity=0.8 
            )
        )
        layout = Layout(
            showlegend=True,
            title="Sum of pairs (SP) scores",
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
    #py.plot(fig)



if __name__ == "__main__":
    main(sys.argv[1])