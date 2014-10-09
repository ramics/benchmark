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
        read_lengths[filename] = values["length"]
        filenames.append(filename)

    tools = os.listdir(folder)

    # get specific info
    for tool in tools:
        pkl_file = open(os.path.join(folder, tool), 'rb')
        (reads, results) = pickle.load(pkl_file)
        final_scores[tool] = [0] * 4
        final_ranks[tool] = [0] * 4

        num_tests[tool] = len(results)
        for filename, values in results.iteritems():

            read_scores[filename][tool] = values["scores"]
            for i in range(len(values["scores"])):
                final_scores[tool][i] += values["scores"][i]

    for filename, tool_scores in read_scores.iteritems():

        # We only rank on the ones that didn't error for anyone
        if len(tool_scores) < len(tools):
            continue
        num_ranks += 1

        raw_scores = [(tool, scores) for tool, scores in tool_scores.iteritems()]
        ranks = [
            rankdata([1 - r[1][score] for r in raw_scores]
            ) for score in range(4)]
        for p, (tool, _) in enumerate(raw_scores):
            for score_type in range(4):    
                final_ranks[tool][score_type] += ranks[score_type][p]

    results = []
    
    print "Tool\ttests\tq\ttc\tcline\tmodeler\tq\ttc\tcline\tmodeler"
    for tool in tools:
        results = [tool, num_tests[tool]]
        results.extend([s/float(num_tests[tool]) for s in final_scores[tool]])
        results.extend([s/float(num_ranks) for s in final_ranks[tool]])
        print "%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (
                tuple(results) 
        )

    length_bins = [round(x) for x in stats.mstats.mquantiles([read_lengths[f] for f in filenames], [0, 0.2, 0.4, 0.6, 0.8, 1])]
    
    information_traces = []
    length_traces = []
    #pp.pprint(read_scores)
    for tool in tools:
        chosen_files = [f for f in filenames if tool in read_scores[f].keys()]
        length_trend_range = [read_lengths[f] for f in chosen_files]
        trend_data = [read_scores[f][tool][2] for f in chosen_files]

        #print read_scores

        trend_data_length = [x for (y, x) in sorted(zip(length_trend_range, trend_data))]
       
        length_trend_range = sorted(length_trend_range)
        


        #length_bins = [round(x) for x in stats.mstats.mquantiles(length_trend_range, [0, 0.2, 0.4, 0.6, 0.8, 1])]
        #info_bins = [round(x, 2) for x in stats.mstats.mquantiles(info_trend_range, [0, 0.2, 0.4, 0.6, 0.8, 1])]

        digitized_length = numpy.digitize(length_trend_range, length_bins)
       

        digitized_length_strings = ["<" + str(i) for i in length_bins]
        digitized_length_strings.append(">" + str(length_bins[len(length_bins) - 1]))



        binned_lengths = [digitized_length_strings[digitized_length[i]] for i in range(len(digitized_length))]

        length_traces.append(Box(
                y=trend_data_length,
                x=binned_lengths,
                name=tool,
                opacity=0.75 
            )
        )

    data2 = Data(length_traces)
    layout2 = Layout(
            boxmode='group',
            showlegend=True,
            title="Relationship between QScore and number of sequences in MSA",
            titlefont = Font(
                size = 28),
            yaxis=YAxis(
                title='Score',
                showgrid=False,
                autorange = True,
                titlefont = Font(
                    size = 20)
            ),
            xaxis = XAxis(
                title = "Number of sequences",
                titlefont = Font(
                    size = 20)
            ),
            legend = Legend(
                font = Font(
                    size = 18
                ),
                yanchor = "bottom",
                xanchor = "left")
    )
    fig2 = Figure(data=data2, layout=layout2)

    py.sign_in("imogen", "mtf1tawct2")
    #py.plot(fig2)


if __name__ == "__main__":
    main(sys.argv[1])