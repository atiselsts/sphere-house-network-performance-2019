#!/usr/bin/python3

import os, sys, time, re, copy, json, datetime, math
import numpy as np
import pylab as pl                                                
import matplotlib.legend_handler as lh

LABELS = [
    "Single gateway",
    "Two gateways, layer-3 load balancing",
    "Two gateways, layer-2 load balancing",
]

OUT_DIR = "../plots"

NUM_SOURCES = 8
POINTS = [6, 10, 14, 18]
TOTAL_PPS = [x * NUM_SOURCES for x in POINTS]

COLORS = ["red", "orange", "blue"]
MARKERS = ["s", "s", "s", "s"]

data = {"l2": [[100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0], [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0], [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0], [100.0, 99.91787571858747, 99.77888336097291, 100.0, 100.0, 100.0, 100.0, 100.0]], "baseline": [[100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0], [99.97989141363362, 91.94752774974774, 96.87373941105284, 99.89983974358974, 97.55511022044088, 99.09837707874173, 99.33854479855682, 99.81967541574835], [73.07749077490774, 61.69800059683676, 72.32537577365163, 71.01449275362319, 71.15232207956888, 71.34911806769426, 71.80190353282788], [57.594862813777, 51.28564749883123, 55.552782630396806, 55.30236634531113, 56.158615042695054, 55.45361875637105, 54.9519586104952]], "l3": [[100.0, 99.8623537508603, 100.0, 99.79296066252589, 99.40284795590262, 100.0, 99.72426470588235, 99.9769159741459], [100.0, 100.0, 100.0, 100.0, 99.92774566473989, 100.0, 100.0, 100.0], [99.9844744604875, 99.96894892097501, 99.92243251628918, 100.0, 100.0, 100.0, 99.98446963814257, 99.39412769923878], [98.68269008550959, 95.56844547563806, 91.37931034482759, 99.21387283236994, 91.92532088681448, 97.05677867902665, 89.81503160852259, 97.12696941612604]]}

def symlog(x):
    if x <= 1:
        return 0
    return math.log(x, 10)

def plot(allData, filenameOut):
    pl.figure(figsize=(4.5, 3))

    # turn off vertical grid                                                                                                                                       
    pl.grid(axis="x")

    width = 0.5

    i = 0
    for data, label in zip(allData, LABELS):
        x = TOTAL_PPS
        print("mean=", [np.mean(u) for u in data])
        packets_lost_per_100thousand = []
        for node in data:
            packets_lost_per_100thousand.append([(100 - u) * 1000 for u in node])
        y = [symlog(np.mean(u)) for u in packets_lost_per_100thousand]
        yerr1 = [symlog(min(u)) for u in packets_lost_per_100thousand]
        yerr2 = [symlog(max(u)) for u in packets_lost_per_100thousand]

        # change the values from absolute to relative error values 
        yerr1 = [u - v for u,v in zip(y, yerr1)]
        yerr2 = [v - u for u,v in zip(y, yerr2)]

        pl.errorbar(x, y, [yerr1, yerr2], label=label, color=COLORS[i], marker=MARKERS[i], ) #markersize=10, markeredgewidth=2 )
        i += 1

    pl.ylabel("PDR, %")
    pl.xlabel("Packets per second, total")

    pl.xticks(TOTAL_PPS, [str(x) for x in TOTAL_PPS])

    pl.xlim(4 * NUM_SOURCES, 20 * NUM_SOURCES)

    pl.ylim(5, -0.5)
    pl.yticks(range(6), ["100%", "99.99%", "99.9%", "99%", "90%", "0%"])

    legend = pl.legend(
        bbox_to_anchor=(0.5, 1.33), loc='upper center', ncol=1,
        handler_map={lh.Line2D: lh.HandlerLine2D(numpoints=1)})

    pl.savefig(OUT_DIR + "/" + filenameOut, format='pdf',
               bbox_extra_artists=(legend,),
               bbox_inches='tight')
    pl.close()


def main():
    plot([                                                                                                                                                         
        data["baseline"],
        data["l3"],
        data["l2"]], "load_balancing.pdf")

main()
