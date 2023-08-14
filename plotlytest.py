
import dash
import dash_html_components as html
import dash_core_components as dcc
import numpy as np
import plotly.graph_objs as go

from dash.dependencies import Input, Output
import fastaparser
import model

Model = model.def_model()
Model.load_weights("metalogs/logs/fit/20230813-104958-300-300/my_model_weights.h5")
FILE_PATH =  "./data/Saccharomyces_cerevisiae/GCF_000146045.2_R64_genomic.fna"
#FILE_PATH = "./Val_Data/AllValPromotorSeq1000.fsa"
RANGE = [17845,19589]
SEQUENCE_NAME = "NC_001137.3"
#SEQUENCE_NAME = "Saccharomyces_cerevisiae-10"

def get_sequence(FILE_PATH, RANGE, SEQUENCE_NAME):
    
    Range_start = RANGE[0]-1
    Range_end = RANGE[1] + 1
     

    with open(FILE_PATH, 'r') as FastaFile:
        parser = fastaparser.Reader(FastaFile, parse_method='quick')
    
        for seq in parser:
            if SEQUENCE_NAME not in seq.header: continue
            print("Sequence:", seq.header)
            sequence = seq.sequence
            sequence_len = len(sequence)
            print("Sequence:", sequence_len)
            if sequence_len < 300 :
                print("SEKVENCE JE MOC KRATKA, VELIKOST JEN: ", sequence_len)

            if Range_start < 149 or sequence_len < Range_end + 151:
                print("ZMENTE RANGE, SEKVENCE JE MOC NA KRAJICH")

            return sequence[Range_start-149: Range_end + 151]

        
seq = get_sequence(FILE_PATH, RANGE, SEQUENCE_NAME)
print(len(seq))

def IterateOverSeq(seq, i):
    window = 300
    stride = 1
    len(seq)
    while window <= len(seq):
        i += stride
        window += stride 
        return seq[i+150:i+150+300]
    

def toNumbers(letter):
    switcher = {
        "A": 0,
        "T": 1,
        "C": 2,
        "G": 3,
        "a": 0,
        "t": 1,
        "c": 2,
        "g": 3,
        }   
    return switcher.get(letter,0)
        

def GetPredictions(seq):
    x = [toNumbers(x) for x in seq]
    X = []
    X.append(x)
    X= np.asarray(X, dtype=np.int32)

    predictions = Model.predict(X, verbose=0)
    return predictions




X =  list(range(RANGE[0]-2, RANGE[1]+1))
OFFGENE = [0]
CDS = [0] 
PROMOTOR = [0]
# Example app.
figure = dict(data=[{'x': [], 'y': [], 'name' : 'CDS'}, {'x': [], 'y': []} ], layout=dict(xaxis=dict(range=[RANGE[0], RANGE[1]+1]), yaxis=dict(range=[0, 1])))
app = dash.Dash(__name__, update_title=None)  # remove "Updating..." from title
app.layout = html.Div([dcc.Graph(id='graph', animate=True), dcc.Interval(id="interval", interval=50)])


@app.callback(Output('graph', 'figure'), [Input('interval', 'n_intervals')])
def update_data(n_intervals):
    
    if type(n_intervals) != int :
        index = 0    
    else :
        index = n_intervals

    seqtopred = IterateOverSeq(seq, index)
    pred = GetPredictions(seqtopred)
    OFFGENE.append(list(pred[0])[0])
    CDS.append(list(pred[0])[1])
    PROMOTOR.append(list(pred[0])[2])
    # tuple is (dict of new data, target trace index, number of points to keep)


    data = go.Scatter(
        x = list(X),
        y=list(OFFGENE),
        name='OFFgene',
    )

    data2 = go.Scatter(
        x = list(X),
        y=list(CDS),
        name='CDS',
    )

    data3 = go.Scatter(
        x = list(X),
        y=list(PROMOTOR),
        name='Promotor',
    )
    return {'data': [data, data2, data3],'layout' : go.Layout(xaxis=dict(range=[RANGE[0], RANGE[1]+1]),
                                                yaxis=dict(range=[0,1]),)}
    return dict(x=[[x[index]]], y=[[y[index]]]) , [0]



if __name__ == '__main__':
    app.run_server()