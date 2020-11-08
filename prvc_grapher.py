import json
import plotly
import plotly.graph_objs as go
import plotly.io
from plotly.offline import iplot, init_notebook_mode
# Using plotly + cufflinks in offline mode
# import cufflinks
# cufflinks.go_offline(connected=True)
# init_notebook_mode(connected=True)
import numpy as np
import pandas as pd


def read_from_json(fname):
    with open(fname, 'r') as jfile:
        data = json.load(jfile)

    return data


def read_from_multiple_csv(file_list):
    datalist = []
    for file in file_list:
        datalist.append(pd.read_csv(file))

    return datalist


def read_from_multiple_josns(file_list):
    datadict = {}
    for fname in file_list:
        with open(fname, 'r') as jfile:
            datadict[fname.split('_')[2]] = json.load(jfile)

    return datadict


def bar_graph(data):
    plotly.offline.plot({
        'data': [go.Bar(x=list(data['_total_snps'].keys()), y=list(data['_total_snps'].values()))],
        'layout': go.Layout(title='Distribution of SNPs by Type (4sU-100)', xaxis={'title':'Type of SNP (Ref/Read)'},
                            yaxis={'title':'Number of Occurrences'})

    }, auto_open=True)


def percent_snps_tc_graph(datadict):
    plotly.offline.plot({
        'data': [go.Bar(x=['Zero', 'Ten', 'Fifty', 'One Hundred'],
                        y=[datadict[i]['_percent_snps_in_file']['T/C'][1] for i in datadict.keys()])],
        'layout': go.Layout(title='Percent of SNPs that are T/C per Dataset', xaxis={'title': 'Dataset'},
                            yaxis={'title': 'Percentage of SNPs that are T/C'})

    }, auto_open=True)


def percent_snps_not_tc_graph(datadict):
    plotly.offline.plot({
        'data': [go.Bar(x=['Zero', 'Ten', 'Fifty', 'One Hundred'],
                        y=[sum([v[1] for k, v in datadict[i]['_percent_snps_in_file'].items() if k != 'T/C'])
                           for i in datadict.keys()])],
        'layout': go.Layout(title='Percent of SNPs that are not T/C per Dataset', xaxis={'title': 'Dataset'},
                            yaxis={'title': 'Percentage of SNPs that are not T/C'})

    }, auto_open=True)


def percent_bases_tc_graph(datadict):
    plotly.offline.plot({
        'data': [go.Bar(x=['Zero', 'Ten', 'Fifty', 'One Hundred'],
                        y=[(datadict[i]['_percent_snps_in_file']['T/C'][0] / datadict[i]['_total_read_len']) * 100
                           for i in datadict.keys()])],
        'layout': go.Layout(title='Percent of Bases that are T/C SNPs per Dataset', xaxis={'title': 'Dataset'},
                            yaxis={'title': 'Percent of Bases that are T/C SNPs'})

    }, auto_open=True)


def percent_bases_not_tc_graph(datadict):
    plotly.offline.plot({
        'data': [go.Bar(x=['Zero', 'Ten', 'Fifty', 'One Hundred'],
                        y=[(sum(v[0] for k, v in datadict[i]['_percent_snps_in_file'].items() if k != 'T/C') /
                           datadict[i]['_total_read_len']) * 100 for i in datadict.keys()])],
        'layout': go.Layout(title='Percent of Bases that are not T/C SNPs per Dataset', xaxis={'title': 'Dataset'},
                            yaxis={'title': 'Percent of Bases that are not T/C SNPs'})

    }, auto_open=True)


def percent_snp_bases_graph(datadict):
    plotly.offline.plot({
        'data': [go.Bar(x=['Zero', 'Ten', 'Fifty', 'One Hundred'],
                        y=[(sum(v[0] for k, v in datadict[i]['_percent_snps_in_file'].items()) /
                           datadict[i]['_total_read_len']) * 100 for i in datadict.keys()])],
        'layout': go.Layout(title='Percent of Bases that are SNPs per Dataset', xaxis={'title': 'Dataset'},
                            yaxis={'title': 'Percent of Bases that are SNPs'})

    }, auto_open=True)


if __name__ == '__main__':
    file_list2 = ['ERCC_PRVC_Data_and_File_Monitoring/ercc_00043_0_perc_4sU_Concise_PRVC_Data.csv',
                  'ERCC_PRVC_Data_and_File_Monitoring/ercc_00043_10_perc_4sU_Concise_PRVC_Data.csv',
                  'ERCC_PRVC_Data_and_File_Monitoring/ercc_00043_50_perc_4sU_Concise_PRVC_Data.csv',
                  'ERCC_PRVC_Data_and_File_Monitoring/ercc_00043_100_perc_4sU_Concise_PRVC_Data.csv']

    data_list = read_from_multiple_csv(file_list2)
