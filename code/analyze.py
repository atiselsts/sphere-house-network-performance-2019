#!/usr/bin/python3

import os, re, sys, json, math
import pylab as pl
import matplotlib.legend_handler as lh
from matplotlib import ticker
import numpy as np

################################################

DATA_DIR = "../"
OUT_DIR = "../plots"

PACKET_RATES = [2, 4, 6, 8, 10, 12]

SPHERE_NODES = ["3", "4", "5", "6", "7"]
NUM_NODES = len(SPHERE_NODES)

NODE_FILENAMES = [
    "gw1-logF.log--2019-08-22",
    "gw4-logF.log--2019-08-22",
    "gw7-logF.log--2019-08-22",
    "gw9-logF.log--2019-08-22",    
    "gw10-logF.log--2019-08-22",
]

experiments = [
    ["orchestra-6", "Orchestra (6 slots)"],
    ["orchestra-17", "Orchestra (17 slots)"],
    ["orchestra-99", "Orchestra (99 slots)"],
    ["oldmsf", "MSF v2"],
    ["newmsf", "MSF v4"],
    ["sphere", "SPHERE"]
]

################################################

# regexp helper
class Matcher:
    def __init__(self, pattern, flags=0):
        self._pattern = re.compile(pattern, flags)
        self._hit = None
    def match(self, line):
        self._hit = re.match(self._pattern, line)
        return self
    def search(self, line):
        self._hit = re.search(self._pattern, line)
        return self._hit
    def matched(self):
        return self._hit != None
    def group(self, idx):
        return self._hit.group(idx)
    def as_int(self, idx):
        return int(self._hit.group(idx))

################################################

# host logs:
# 1. "Publishing "
# 2. "1533912951.0 > "

MATCHER_PUBLISHING = Matcher(r"Publishing (.*)$")
MATCHER_TIMESTAMP = Matcher(r"([0-9]+)\.[0-9] > .*$")


# node logs:
# 1. "1533909566.4 > LINK STATS to 1: "
# 2. "1533909566.5 > drop: "
# 3. "1566502790.2 > *** starting packetgen: 6 Hz"

MATCHER_PACKETGEN = Matcher(r"([0-9]+)\.[0-9] > \*\*\* starting packetgen: ([0-9]+) Hz$")
MATCHER_LINK_STATS = Matcher(r"([0-9]+)\.[0-9] > LINK STATS to 1: ([0-9]+) ([0-9]+) ([0-9]+)$")
MATCHER_DROP_CAUSE = Matcher(r"([0-9]+)\.[0-9] > drop: ([0-9]+) ([0-9]+) ([0-9]+)$")

################################################

def safediv(a, b):
    if b == 0:
        return 0
    return a / b

################################################

class ExperimentRun:
    def __init__(self, rate):
        self.rate = rate
        self.seqnums = {}
        self.link_sent = {}
        self.link_acked = {}
        for n in SPHERE_NODES:
            self.seqnums[n] = set()
            self.link_sent[n] = 0
            self.link_acked[n] = 0
        self.sent = 0
        self.acked = 0
        self.dropped_rtrx = 0
        self.dropped_queue = 0
        self.dropped_other = 0

    def add_seqnum(self, node, seqnum):
        self.seqnums[node].add(seqnum)

    def get_pdr(self):
        expected = 0
        received = 0
        expected_per_node = []
        received_per_node = []

        for n in SPHERE_NODES:
            if len(self.seqnums[n]) == 0:
                print("no packets for ", n)
            else:
                e = max(self.seqnums[n]) - min(self.seqnums[n]) + 1
                r = len(self.seqnums[n])
                expected += e
                received += r
                expected_per_node.append(e)
                received_per_node.append(r)
        if expected == 0:
            print("no packets at all")
        return safediv(received, expected), \
            [safediv(r, u) for r, u in zip(received_per_node, expected_per_node)]

    def add_link_stats(self, node, sent, acked):
        self.link_sent[node] += sent
        self.link_acked[node] += acked

    def add_drop_cause(self, rtrx, queue, other):
        self.dropped_rtrx += rtrx
        self.dropped_queue += queue
        self.dropped_other += other

    def get_prr(self):
        sent = 0
        acked = 0

        per_node = []
        for n in SPHERE_NODES:
            sent += self.link_sent[n]
            acked += self.link_acked[n]
            per_node.append(safediv(self.link_acked[n], self.link_sent[n]))

        return safediv(acked, sent), per_node

################################################

class Experiment:
    def __init__(self, name, desc):
        self.name = name
        self.desc = desc
        self.rates = {}
        for r in PACKET_RATES:
            self.rates[r] = ExperimentRun(r)

all_experiments = {}

################################################

def is_readable(filename):
    return os.path.isfile(filename) and os.access(filename, os.R_OK)

################################################

def import_results(experiment):

    dirname = os.path.join(DATA_DIR, experiment.name)
    print("Import from " + dirname)

    rate_start_ts = {}

    for i, n in enumerate(SPHERE_NODES):
        filename = NODE_FILENAMES[i]
        full_filename = os.path.join(dirname, filename)
        if not is_readable(full_filename):
            print("Not accessible: ", full_filename)
            continue

        previous_rate = -1
        current_rate = -1
        current_rate_start_ts = -1

        with open(full_filename, "r") as f:
            for line in f.readlines():
                line = line.strip()

                m = MATCHER_PACKETGEN.match(line)
                if m.matched():
                    ts = m.as_int(1)
                    previous_rate = current_rate
                    current_rate = m.as_int(2)
                    current_rate_start_ts = ts + 30
                    rate_start_ts[current_rate] = ts
                    continue

                m = MATCHER_LINK_STATS.match(line)
                if m.matched():
                    ts = m.as_int(1)

                    # decide which rate to use
                    if ts >= current_rate_start_ts:
                        rate = current_rate
                    else:
                        rate = previous_rate
                    if rate != -1:
                        # valid rate, add to the results
                        sent = m.as_int(2)
                        acked = m.as_int(3)
                        experiment.rates[rate].add_link_stats(n, sent, acked)
                    continue

                m = MATCHER_DROP_CAUSE.match(line)
                if m.matched():
                    ts = m.as_int(1)
                    # decide which rate to use
                    if ts >= current_rate_start_ts:
                        rate = current_rate
                    else:
                        rate = previous_rate
                    if rate != -1:
                        # valid rate, add to the results
                        rtrx = m.as_int(2)
                        queue = m.as_int(3)
                        other = m.as_int(4)
                        experiment.rates[current_rate].add_drop_cause(rtrx, queue, other)
                    continue

    filename = os.path.join(dirname, "host.log")
    if not is_readable(filename):
        print("Not accessible: ", filename)
        return

    current_rate = 0
    next_rate = 2

    with open(filename, "r") as f:
        for line in f.readlines():
            line = line.strip()

            m = MATCHER_TIMESTAMP.match(line)
            if m.matched():
                ts = m.as_int(1)
                if next_rate <= PACKET_RATES[-1] and ts >= rate_start_ts.get(next_rate, 0):
                    # start a new rate
                    current_rate = next_rate
                    next_rate += 2
                continue

            if current_rate:
                m = MATCHER_PUBLISHING.match(line)
                if m.matched():
                    #print("load publishing")
                    try:
                        obj = json.loads(m.group(1))
                    except:
                        print("failed to parse json:", line)
                        continue

                    if "gw" in obj:
                        stat = obj["gw"][0]
                        seqnum = stat["mc"]
                        node = stat["uid"][-1]
                        experiment.rates[current_rate].add_seqnum(node, seqnum)

################################################

def read_stats():
    for name, description in experiments:
        all_experiments[name] = Experiment(name, description)
        import_results(all_experiments[name])
        print("")


################################################

TOTAL_PPS = [x * NUM_NODES for x in PACKET_RATES]

MARKERS = ["o", "o", "o", "v", "v", "."]
COLORS = ["orange", "red", "brown", "lightgreen", "darkgreen", "blue"]

def symlog(x):
    if x <= 1:
        return 0
    return math.log(x, 10)

def plot(all_data, all_data_minmax, ylabel, filename_out):
    pl.figure(figsize=(4, 2.5))

    # turn off vertical grid
    pl.grid(axis="x")

    width = 0.5
    i = 0

    for i, data in enumerate(all_data):
        # drop the zeros: refer to experiments
        d = data #[u for u in data if u != 0.0]
        x = TOTAL_PPS[:len(d)]

        e = all_data_minmax[i]
        label = experiments[i][1]

        if "PDR" in ylabel:
            # use symlog scale
            packets_lost_per_100thousand = [(1.0 - u) * 100000 for u in d]
            y = [symlog(u) for u in packets_lost_per_100thousand]
            print("lost=", packets_lost_per_100thousand)
            print("     ", y)

            packets_lost_per_100thousand = [(1.0 - u) * 100000 for u,v in e]
            yerr1 = [symlog(u) for u in packets_lost_per_100thousand]
            print("max lost=", packets_lost_per_100thousand)            

            packets_lost_per_100thousand = [(1.0 - v) * 100000 for u,v in e]
            yerr2 = [symlog(u) for u in packets_lost_per_100thousand]
            print("min lost=", packets_lost_per_100thousand)            

        else:
            # use linear scale
            y = [u * 100 for u in d]
            yerr1 = [u * 100 for u,v in e]
            yerr2 = [v * 100 for u,v in e]

        print(y)
        print(yerr1)
        print(yerr2)
        print()

        # change the values from absolute to relative error values 
        yerr1 = [u - v for u,v in zip(y, yerr1)]
        yerr2 = [v - u for u,v in zip(y, yerr2)]

        pl.errorbar(x, y, [yerr1, yerr2], label=label, # ls="None",
                    alpha=(0.5 if "SPHERE" in label else 1.0),
                    marker=MARKERS[i],
                    markersize=4, markeredgewidth=2,  color=COLORS[i])
        i += 1

    pl.ylabel(ylabel)
    pl.xlabel("Packets per second, total")

    pl.xticks(TOTAL_PPS, [str(x) for x in TOTAL_PPS])

    if "PDR" in ylabel:
        pl.ylim(5, -0.5)
        pl.yticks(range(6), ["100%", "99.99%", "99.9%", "99%", "90%", "0%"])
    else:
        pl.ylim([0, 100])

    legend = pl.legend(
        bbox_to_anchor=(0.5, 1.42), loc='upper center', ncol=2,
        handler_map={lh.Line2D: lh.HandlerLine2D(numpoints=1)})

    pl.savefig(OUT_DIR + "/" + filename_out, format='pdf',
               bbox_extra_artists=(legend,),
               bbox_inches='tight')
    pl.close()

################################################

def show_stats():
    all_pdr = []
    all_prr = []

    all_pdr_minmax = []
    all_prr_minmax = []

    for name, desc in experiments:
        print(name)
        exp = all_experiments[name]

        pdr = []
        pdr_minmax = []

        prr = []
        prr_minmax = []

        for r in PACKET_RATES:
            all, per_node = exp.rates[r].get_pdr()
            pdr.append(all)
            pdr_minmax.append((min(per_node, default=0), max(per_node, default=0)))

            all, per_node = exp.rates[r].get_prr()
            prr.append(all)
            prr_minmax.append((min(per_node, default=0), max(per_node, default=0)))

            print("{} : {:.6f}% PDR {:.2f}% PRR {} {} {} dropped".format(
                r, 100.0 * exp.rates[r].get_pdr()[0],
                100.0 * exp.rates[r].get_prr()[0],
                exp.rates[r].dropped_rtrx,
                exp.rates[r].dropped_queue,
                exp.rates[r].dropped_other
            ))
            print("  PRR:", exp.rates[r].get_prr()[1])
        print("")

        all_pdr.append(pdr)
        all_prr.append(prr)

        all_pdr_minmax.append(pdr_minmax)
        all_prr_minmax.append(prr_minmax)

    plot(all_pdr, all_pdr_minmax, "End-to-end PDR, %", "sphere_vs_orchestra_6top_pdr.pdf")
    plot(all_prr, all_prr_minmax, "Link-layer PAR, %", "sphere_vs_orchestra_6top_prr.pdf")

################################################

def plot_energy_estimate():
    num_rx_slots_per_second = {
        "orchestra-6" : 100 / 6.0 + 100 / 99.0,
        "orchestra-17" : 100 / 17.0 + 100 / 99.0,
        "orchestra-99" : 100 / 99.0 + 100 / 99.0,
        "oldmsf" : 100 / 99.0 + 100 / 99.0,
        "newmsf" : 100 / 99.0 + 100 / 99.0,
        "sphere" : 100 / 99.0 + 100 / 99.0,
    }

    voltage = 3.6
    current_consumption_for_rx_ma = 7.22
    guard_time_seconds = 0.001 # 1ms

    charge_for_rx_slot_mc = current_consumption_for_rx_ma * guard_time_seconds
    energy_for_rx_slot_mj = voltage * charge_for_rx_slot_mc
    print("mc=", charge_for_rx_slot_mc)
    print("mj=", energy_for_rx_slot_mj)

    fig, ax1 = pl.subplots()
    fig.set_size_inches(4, 3)
    ax2 = ax1.twinx()
    ax1.grid(b=False, axis="x")
    ax2.grid(b=False, axis="both")
    width = 0.5

    labels = [x[1] for x in experiments]
    x = np.arange(0, len(labels), 1)

    y = []
    for name, desc in experiments:
        y.append(num_rx_slots_per_second[name])

    ax1.bar(x, y)

    ax1.set_ylim([0, 20])
    ax2.set_ylim([0, 20 * charge_for_rx_slot_mc])

    ax1.set_xlabel('Scheduler')
    ax1.set_ylabel('Slots per sec')
    ax2.set_ylabel('Charge per sec, mC')
    
    ax1.set_xticklabels([""] + labels, rotation=90)

#   legend = pl.legend()
    pl.savefig(OUT_DIR + "/energy_estimate.pdf", format='pdf',
#               bbox_extra_artists=(legend,),
               bbox_inches='tight')
    pl.close()

 
################################################

def main():
    try:
        os.mkdir(OUT_DIR)
    except:
        pass
    read_stats()
    show_stats()
    plot_energy_estimate()

###########################################

if __name__ == '__main__':
    main()
    print("all done!")
