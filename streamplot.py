#  Name: 
#  Author: rotem.tal
#  Description:
#
import pandas as pd
from matplotlib.pyplot import *
from matplotlib.lines import Line2D
import numpy as np
import numbers
import csv
import ast
from collections import Counter
import os
import warnings
import fire
from colour import Color
import argparse
import json
from io import StringIO
import pickle
from argparse import ArgumentParser


class StreamPlot:
    UNOCCURRING_STR = 'NOTHINGNEW'
    DEF_SUB_N = "DEFAULT_SUB_NM"

    def __init__(self, subject_name_format='sub_$R_$N_$T', subject_name_sep='_', cmap=None, consts=None,
                 print_other=False, other_n='Other', other_color='k', bacteria_taxo_foramt='$S$B__',
                 bacteria_taxonomic_sep='|', depth='p', depth_dict=None, rec_col='record_id'):
        """
        Initilize a stream plot object
        :param subject_name_format: Format of samples name, $R denotes the ID, $N the sample number, $T sample classification (time stamp, sick visits etc..)
        :param subject_name_sep: Separator used in sample names
        :param cmap: Color map (as dictionary)
        :param consts: Mapping from sample number to time as dictionary (e.g. sample 4 is 6 months, than the entry '4': 6 should be inserted)
        :param print_other: Print bacteria not present in the color map
        :param other_n: Name for bacteria not present in the color map
        :param other_color: Color for bacteria not present in the color map
        :param bacteria_taxo_foramt: Bacteria name format, %S denots separator, $B denotes taxa name
        :param bacteria_taxonomic_sep: Separator between bacteria taxa ranking
        :param depth: Depth at which gradient would be printed
        :param depth_dict: Dictionary for possible taxonomic rankings (default is from kingdom to species)
        :param rec_col: Meta-Data record ID column name
        """
        _cmap = cmap
        if isinstance(cmap, str):
            with open(cmap) as f:
                if cmap.endswith('json'):
                    _cmap = json.load(f)
                elif cmap.endswith("csv"):
                    reader = csv.reader(f)
                    _cmap = {r[0]: r[1] for r in reader}
                else:
                    try:
                        _cmap = ast.literal_eval(cmap)
                    except:
                        warnings.warn(f"Color map parameter should be a valid: json file, csv file, string representing python dict{os.linesep}Using default coloring")
                        _cmap = None
        self.cmap = _cmap if _cmap is not None else {'p__Actinobacteria': '#ffff00', 'p__Proteobacteria': '#009800', 'p__Bacteroidetes': '#0000ff', 'p__Firmicutes': '#cc00cc', 'p__Verrucomicrobia': "#ff1a75"}
        self.next_depth = depth_dict if depth_dict is not None else {'k': 'p', 'p': 'c', 'c': 'o', 'o': 'f', 'f': 'g', 'g': 's', 's': self.UNOCCURRING_STR}
        self.bacteria_abundance_sum = None
        self.other = print_other
        self.other_n = other_n
        if self.other:
            self.cmap[self.other_n] = other_color
        self.consts = consts if consts else {'0': 0, '1': 0.5, '2': 1, '3': 2, '4': 4, '5': 6, '6': 9, '7': 12}
        self.sub_nm_format = subject_name_format
        if depth[0].lower() in self.next_depth.keys():
            self.depth = depth[0].lower()
        else:
            self.depth = 'p'
            warnings.warn('Invalid bacteria depth parameter, defaulting to phylum level', RuntimeWarning)
        self._bacteria_limiter = bacteria_taxo_foramt.replace("$B", self.next_depth[self.depth]).replace('$S', bacteria_taxonomic_sep)
        self.depth = bacteria_taxo_foramt.replace('$B', self.depth).replace('$S', "")
        self._bacteria_sep = bacteria_taxonomic_sep
        self.sep = subject_name_sep
        self.extract_from_name_format()
        self._subjects = {}
        self._events = {}
        self._ind_len = 0
        self._meta = None
        self._rec = rec_col
        self._noise = None
        self.complete_color_array = {}

    def extract_from_name_format(self):
        """
        finds and stores the indices (denoted by the subject name separator) of ID ($R), Number ($N) and tag ($T) in the name format
        """
        try:
            self._r = np.argwhere([0 if i != '$R' else 1 for i in self.sub_nm_format.split(self.sep)]).ravel()[0]
        except:
            self._r = None
            warnings.warn("Name format doesn't contain $R as an indicator for subject ID, cannot insert "
                             "more than 1 subject at a time")
        try:
            self._t = np.argwhere([0 if i != '$T' else 1 for i in self.sub_nm_format.split(self.sep)]).ravel()[0]
        except:
            self._t = None
        try:
            self._n = np.argwhere([0 if i != '$N' else 1 for i in self.sub_nm_format.split(self.sep)]).ravel()[0]
        except:
            self._n = None
            warnings.warn("Format does not contain sample numbering  mentioned as $N, software behaviour may be unexpected")

    def insert_subject(self, sub_csv, header=0, index_col=0, delim=','):
        """
        Insert a single subject from csv file
        :param sub_csv: Path to csv
        :param header: Header row
        :param index_col: Index column
        :param delim: Delimiter (for csv)
        :return:
        """
        if not self._r and len(self._subjects):
            warnings.warn("Name format did not contain indicator or ID, cannot insert more than 1 subject at a time",
                          RuntimeWarning)
            return False
        df = pd.read_csv(sub_csv, header=header, index_col=index_col, delimiter=delim)
        if self._r:
            nm = df.columns[0].split(self.sep)[self._r]
        else:
            nm = StreamPlot.DEF_SUB_N
        df = df.reindex(sorted(df.columns, key=lambda x: int(x.split(self.sep)[self._n])), axis=1) if self._n else df
        self.bacteria_abundance_sum = (len(self._subjects) * self.bacteria_abundance_sum + df.sum(1)) / (len(self._subjects) + 1) if self.bacteria_abundance_sum is not None else df.sum(1)
        self._subjects[nm] = df
        self._ind_len = len(df.index)
        self.color_by_abundance()
        return self

    def load_subjects(self, sub_csv, header=0, index_col=0, delim=','):
        """
        Insert a csv of several subjects
        :param sub_csv: Path to csv
        :param header: Header row
        :param index_col: Index column
        :param delim: csv delimiter
        :return:
        """
        if not self._r:
            warnings.warn("Name format did not contain indicator or ID, cannot insert more than 1 subject at a time",
                          RuntimeWarning)
            return False
        df = pd.read_csv(sub_csv,header=header, index_col=index_col, delimiter=delim)
        unique = np.unique([i.split(self.sep)[self._r] for i in df.columns])
        self.bacteria_abundance_sum = (len(self._subjects) * self.bacteria_abundance_sum + df.sum(1)) / (len(self._subjects) + len(unique)) if self.bacteria_abundance_sum is not None else df.sum(1)
        ind = [(i, [j for j in df.columns if j.split(self.sep)[self._r] == i]) for i in unique]
        self._ind_len = len(df.index)
        if self._n:
            self._subjects.update({i: df[cols].reindex(sorted(df[cols].columns,
                                                              key=lambda x: int(x.split(self.sep)[self._n])),
                                                       axis=1) for i, cols in ind})
        else:
            self._subjects.update({i: df[cols] for i, cols in ind})
        self.color_by_abundance()
        return self

    def color_by_abundance(self):
        """
        Defines the color mapping of the precision depth based on the overall abundance of the taxa in the data
        :return:
        """
        for fam in self.cmap.keys():
            if fam == self.other_n:
                continue
            cur_fam = self.bacteria_abundance_sum.drop(list(filter(lambda x: fam not in x, self.bacteria_abundance_sum.index))).sort_values(ascending=False)
            fam_col = Color(self.cmap[fam])
            colors_gradient = list(fam_col.range_to(Color('black'), len(cur_fam.index) + 42))
            i = 0
            for taxa in cur_fam.index:
                if cur_fam[taxa] == 0:
                    self.complete_color_array[taxa] = '#000000'
                else:
                    self.complete_color_array[taxa] = str(colors_gradient[i])
                    i += 1

    def plot_stream(self, sub_names=None, folder=None, title_format="sub %s streamplot", plot_type='wiggle',
                    suffix='png', x_lab='Time (m)', y_lab='Relative abundance', event_color='#380814'):
        """
        Plot streamsplots
        :param sub_names: List of subjects to plot
        :param folder: Folder to save the plots
        :param title_format: Plot title format (%s where the subject ID should be inserted)
        :param plot_type: Plot baseline style
        :param suffix: suffix to save the plots by (file type)
        :param x_lab: Label for X axis
        :param y_lab: Label for y axis
        :param event_color: Color to print events at
        :return:
        """
        if "%s" not in title_format:
            title_format += "%s"
        if sub_names:
            plot_arr = {sub: df for sub, df in self._subjects.items() if sub in sub_names}
        else:
            plot_arr = self._subjects
        if len(plot_arr) == 0:
            warnings.warn("No subjects to plots")
            return
        for rec, sub in plot_arr.items():
            filt_sub = sub.drop([i for i in sub.index if i not in self.complete_color_array.keys()])
            colors = [self.complete_color_array[i] if len(self.complete_color_array[i]) == 7 or not self.complete_color_array[i].startswith('#')
                      else "".join([x*2 if x != '#' else x for x in self.complete_color_array[i]]) for i in filt_sub.index]
            if self.other:
                d = pd.Series(sub.drop(self.complete_color_array.keys()).sum(0), name=self.other_n, index=sub.columns)
                filt_sub = filt_sub.append(d)
                colors.append(self.cmap[self.other_n])
            x_axis = [self.consts[i.split(self.sep)[self._n]] for i in filt_sub.columns] if self._n is not None and \
                                                                                            self.consts is not None else \
                    np.arange(1, len(sub.columns)+1)
            stackplot(x_axis, filt_sub, baseline=plot_type, colors=colors)
            if rec in self._events.keys():
                ymin, ymax = gca().get_ylim()
                axvline(self._events[rec][0], color=event_color)
                gca().text(self._events[rec][0]+0.07, ymax - 0.2, self._events[rec][1])
            leg_ord = np.unique(filt_sub.rename(lambda x: x[x.find(self.depth):x.find(self._bacteria_limiter)]).index)
            k,b = zip(*[(Line2D([0], [0], color=self.cmap[i], lw=5), i) for i in np.flip(leg_ord, 0)])
            legend(k,b, bbox_to_anchor=(1.04, 0.6),loc='upper left')
            title(title_format%rec)
            xlabel(x_lab)
            ylabel(y_lab)
            if folder is not None:
                if not os.path.isdir(folder):
                    os.makedirs(folder)
                if not folder.endswith(os.sep):
                    folder += os.sep
            else: folder = ""
            savefig(folder+plot_type+rec+f'.{suffix}', bbox_inches="tight")
            show()
        return self

    def add_event(self, date, txt, sub_name=None):
        """
        Add an event to be printed on a subjects plot
        :param date: Location on the stream graph
        :param txt: Text for the event (short)
        :param sub_name: Subject to add event to
        :return:
        """
        nm = sub_name if sub_name else list(self._subjects.keys())[0]
        self._events[nm] = (date, txt)
        return self

    def events_from_meta(self,event_name='Onset', rec_col='record_id', event_col='time_sx', time_col=None,
                         event_ide=None,event_thr=None,meta_csv='data/GMAPMetaClean.xlsx', delim=',',sheet_name='Metadata', header=0, index=0):
        """
        Automatically add events from the meta data
        :param event_name: General name for the event (to be printed)
        :param rec_col: Meta data ID column name
        :param event_col: Event column name
        :param time_col: Time column name (to extract event from)
        :param event_ide: Match to indicate event
        :param event_thr: Threshold (for numerical events) to be considered an event
        :param meta_csv: Path to meta data csv
        :param delim: meta data delimiter
        :param sheet_name: sheet name for xls files
        :param header: meta data header row
        :param index: meta data index column
        :return:
        """
        if self._meta is None:
            try:
                self.add_meta(meta_csv, delim,sheet_name, header, index)
            except:
                warnings.warn("No meta data found, events would not be added")
                return False
        if self._r:
            rec_id = self._meta[rec_col].unique()
        else:
            return False
        for i in rec_id:
            c_m = self._meta[self._meta[rec_col] == i]
            e_time = c_m[event_col].iloc[0] if (time_col is None and event_ide is None) else c_m[time_col][c_m[event_col] == event_ide]
            if not isinstance(e_time, numbers.Number):
                try:
                    if len(e_time) > 1:
                        e_time = e_time.iloc[0]
                except:
                    warnings.warn("Ignoring invalid event time input for record id %s, event not added"%i,
                                  RuntimeWarning)
                    continue
            if not (event_thr is None):
                if (c_m[event_col] <= event_thr).any():
                    self._events[str(i)] = (e_time, event_name)
            else:
                self._events[str(i)] = (e_time, event_name)
        return self


    def __empty(self):
        """
        Clears saved data
        :return:
        """
        self._subjects = {}
        self._events = {}
        self._meta = None

    def remove_subject(self, sub_name):
        """
        Delete a certain subject from streamer memory
        :param sub_name: subject name
        :return:
        """
        if sub_name in self._subjects.keys():
            self._subjects.pop(sub_name)
            return True
        return False

    def set_other(self, print_other, n, c):
        """
        Change the Others (bacteria not in the cmap) settings
        :param print_other: Boolean, indicate to print others or not
        :param n: Others printed name
        :param c: Color
        :return:
        """
        self.other = print_other
        self.other_n = n
        self.cmap[self.other_n] = c


    def binary_ordered_cluster(self, k, thr, t_thr=None, clean=False, folder='Plotter'):
        """
        Clusters subjects, ignoring the abundance level
        :param k: Top bacteria
        :param thr: Number of corresponding bacteria per time point
        :param t_thr: Number of time points for a subject to be clustered
        :param clean: Ignore bacteria not in the cmap
        :param folder: Folder path to save to
        :return:
        """
        if self._n is None:
            warnings.warn("Name format must contain sample chronology (marked $N) in order to preform clustering")
            return
        if clean:
            sub_m = {sub: df.drop(list(filter(lambda x: all([fam not in x for fam in self.cmap.keys()]),
                                              df.index))) for sub, df in self._subjects.items()}
        else:
            sub_m = self._subjects
        t_thr = thr if t_thr is None else t_thr
        top_gen = {sub: frozenset(map(lambda x: (int(x.split(self.sep)[self._n]), frozenset(df[x].nlargest(k).index)),
                                      df.columns)) for sub, df in sub_m.items()}
        cluster_sub = list(sub_m.keys())[0]
        clusters = {top_gen[cluster_sub]: [cluster_sub]}
        for sub in top_gen.keys():
            sub_d = dict(top_gen[sub])
            for key in clusters.keys():
                if sub in clusters[key]:
                    continue
                key_d = dict(key)
                inters = {(t, sub_d[t] | key_d[t]) for t in key_d.keys() if t in sub_d.keys() and
                          len(sub_d[t] & key_d[t]) > thr}
                if len(inters) > t_thr:
                    clusters[key].append(sub)
                    inters.update([(t, key_d[t]) for t in key_d.keys() if t not in sub_d.keys()])
                    inters.update([(t, sub_d[t]) for t in sub_d.keys() if t not in key_d.keys()])
                    clusters[frozenset(inters)] = clusters.pop(key)
                    break
            clusters[top_gen[sub]] = [sub]
        printr=[v for v in clusters.values() if len(v)>1]
        if not os.path.isdir(folder):
            os.makedirs(folder)
        for i, c in enumerate(printr):
            c_pa = folder+os.sep+ 'Streamer_Clusteres_k=%s_thr=%s_t_thr=%s'%(k, thr, t_thr)+os.sep+'cluster_%s'%i
            self.plot_stream(c, c_pa)
            if self._meta is not None:
                c_m = self._meta.drop([i for i in self._meta.index if str(self._meta[self._rec][i]) not in c])
                for u in c:
                    c_m = c_m.drop([i for i in c_m.index if str(c_m[self._rec][i])==u][:-1])
                c_n = len(c_m.index)
                with open(c_pa+os.sep+"cluster_stats.txt", 'w') as f:
                    f.write(f"{os.linesep}".join([f"{x}: {os.linesep}\t" + f"{os.linesep}\t".join([f'{t[0]} - {t[1]/c_n}'
                                                                   for t in Counter(c_m[x]).items()])
                                       for x in c_m.columns]))
        return self



    def cluster(self, k, thr, t_thr, delta, folder):
        """
        Clusters subjects according to graphical clustering algorithm
        :param k: Top bacteria
        :param thr: Number of corresponding bacteria per time point
        :param t_thr: Number of time points for a subject to be clustered
        :param clean: Ignore bacteria not in the cmap
        :param delta: allowed abundance level differencefor bacteria
        :param folder: Folder path to save to
        """
        if self._n is None:
            warnings.warn("Name format must contain sample chronology (marked $N) in order to preform clustering")
            return
        top_gen = {sub: dict(map(lambda x: (int(x.split(self.sep)[self._n]), df[x].nlargest(k)),
                                 df.columns)) for sub, df in self._subjects.items()}
        cluster_sub = list(top_gen.keys())[0]
        clusters = {tuple([cluster_sub]): top_gen[cluster_sub]}
        for sub, df_dict in top_gen.items():
            free = True
            for cluster_tup, cluster_df_dict in clusters.items():
                if sub in cluster_tup:
                    continue
                t_inters = {t: df.index.intersection(cluster_df_dict[t].index) for t, df in df_dict.items() if t in cluster_df_dict.keys()}
                inters = {t: (df_dict[t][i] + len(cluster_tup)*cluster_df_dict[t][i])/(len(cluster_tup)+1)
                          for t, i in t_inters.items() if len(i) >= thr and
                          ((df_dict[t][i] - cluster_df_dict[t][i]).abs() < delta).sum() >= thr}
                if len(inters) >= t_thr:
                    inters.update({t: df for t, df in cluster_df_dict.items() if t not in df_dict.keys()})
                    inters.update({t: df for t, df in df_dict.items() if t not in cluster_df_dict.keys()})
                    clusters[cluster_tup + tuple([sub])] = inters
                    clusters.pop(cluster_tup)
                    free = False
                    break
            if free: clusters[tuple([sub])] = df_dict
        print(len(clusters))
        printr, cluster_size = zip(*[(v,len(v)) for v in clusters.keys() if len(v) > 1])
        folder = folder if folder.endswith(os.sep) else folder + os.sep
        for i, sub_tup in enumerate(printr):
            cur_fold = folder + 'Streamer_Clusteres_k=%s_thr=%s_t_thr=%s_delta=%s' % (k, thr, t_thr, delta) + os.sep
            cur_path = cur_fold + 'cluster_%s' % i
            self.plot_stream(sub_tup, cur_path)
            if self._meta is not None:
                c_m = self._meta.drop([i for i in self._meta.index if str(self._meta[self._rec][i]) not in sub_tup])
                # TODO: add option to print only certain columns from the
                # meta data
                for u in sub_tup:
                    c_m = c_m.drop([i for i in c_m.index if str(c_m[self._rec][i]) == u][:-1])
                c_n = len(c_m.index)
                with open(cur_path + os.sep + "cluster_stats.txt", 'w') as f:
                    f.write(f"{os.linesep}".join([f"{x}: {os.linesep}\t"+ f"{os.linesep}\t".join(['%s - %f' % (t[0], t[1] / c_n)
                                                                     for t in Counter(c_m[x]).items()])
                                       for x in c_m.columns]))
        cluster_med_size, num_clustered = np.median(cluster_size), np.sum(cluster_size)
        num_single = len(top_gen) - num_clustered
        with open(cur_fold + 'summary.txt', 'w') as s:
            s.write(f"Number of clusters: %s{os.linesep}Total subjects clustered: %s{os.linesep}Cluster size median: %s{os.linesep}Number of singletons: %s"%(len(cluster_size), num_clustered, cluster_med_size, num_single))
        return self

    def add_meta(self, meta_csv='data/GMAPMetaClean.xlsx', delim=',',sheet_name='Metadata', header=0, index=0):
        """
        Add meta data
        :param meta_csv: Path to csv, xls or xlsx
        :param delim: data delimiter
        :param sheet_name: for xlsx files
        :param header: Header row
        :param index: Index column
        :return:
        """
        if isinstance(meta_csv, str):
            if meta_csv.endswith("xlsx"):
                self._meta = pd.read_excel(meta_csv, sheet_name=sheet_name)
                self._meta.set_index('sample', inplace=True)
            elif meta_csv.endswith('xls'):
                self._meta = pd.read_excel(meta_csv)
            elif meta_csv.endswith("csv"):
                self._meta = pd.read_csv(meta_csv,header=header, index_col=index, delimiter=delim)
            else:
                warnings.warn(f"Meta data must be one of the following types: \".csv\" \".xls\" \".xlas\"{os.linesep}")
        else:
            self._meta = meta_csv
        return self

    def save(self, path):
        """
        Save current stackplot object
        :param path: Path to save to
        :return:
        """
        if os.sep in path and not os.path.isdir(path[:path.rfind(os.sep)]):
            os.mkdir(path[:path.rfind(os.sep)])
        path = path + ".SCG" if not path.endswith(".SCG") else path
        with open(path, '+wb') as of:
            pickle.dump(self, of)



def read_dict(some_dict):
    """
    Reads a colormap dictionary from json, csv or a string representing a
    dictionary
    :param some_dict: Path to .json or .csv file, or string
    :return:
    """
    if os.path.isfile(some_dict):
        with open(some_dict) as f:
            if some_dict.endswith('json'):
                n_dict = json.load(f)
            elif some_dict.endswith("csv"):
                reader = csv.reader(f)
                n_dict = {r[0]: r[1] for r in reader}
            else:
                print("Error: Invalid string, expected be a valid: json file, csv file or string representing python dict")
                return None
    else:
        try:
            n_dict = ast.literal_eval(some_dict)
            assert isinstance(n_dict, dict)
        except:
            print(
                "Error: Invalid string, expected be a valid: json file, csv file or string representing python dict")
            return None
    return n_dict

def get_name(streamer):
    """
    Get subject name format from the user
    :param streamer: Current streamer object
    :return:
    """
    while True:
        streamer.sub_nm_format = input(f"Enter name format:{os.linesep}")
        streamer.sep = input(f"Enter name separator:{os.linesep}")
        accept = input(
            f"new subject name format is {streamer.sub_nm_format} where {streamer.sep} is the separator{os.linesep}Ok? (y/n/x){os.linesep}")
        if accept == 'y':
            streamer.extract_from_name_format()
            return True
        if accept == 'x':
            return False


def construct_streamer(streamer):
    """
    Construct the streamer object
    :param streamer: object to be constructed or modified
    :return:
    """
    bin_opt = ['y', 'n']
    print("Please consult the README file for clarity. We will now modify the configurations step by step")
    print("You can save the configuration using \'s\' or start again using x (Will keep modifications)")
    print(f"Current subject name format is {streamer.sub_nm_format} where {streamer.sep} is the separator")
    while True:
        action = input(f"Change subject name format? (y/n/x/s){os.linesep}")
        if action == 'n':
            break
        if action == 'y':
            if get_name(streamer):
                break
        if action == 'x':
            run(streamer)
            return -1
        if action == 's':
            streamer.save(input("Enter save path:"))
        print("Invalid key")
    print(
        f"You can add a mapping from sample number ($N in the subject name format) to time{os.linesep}For example the current mapping is {json.dumps(streamer.consts, indent=4)} and the time is in motnh.")
    while True:
        action = input(f"Change the mapping? (y/n/s/x){os.linesep}")
        if action == 'n':
            break
        if action == 'y':
            d = read_dict(input(f"Enter mapping dictionary (path or string):{os.linesep}"))
            if d == None:
                print("Could not load mapping")
            else:
                streamer.consts = d
                break
        if action == 'x':
            run(streamer)
            return -1
        if action == 's':
            streamer.save(input("Enter save path:"))
    print(f"Default color map is {json.dumps(streamer.cmap, indent=4)}")
    while True:
        print("Change color map? (y/n/s/x)")
        action = input()
        if action == 'n':
            break
        if action == 'y':
            cmap_d = read_dict(input(f"Enter color map (path or string):{os.linesep}"))
            if cmap_d is None:
                print("Could not load color map")
            else:
                streamer.cmap = cmap_d
                break
        else:
            print("Invalid key")
        if action == 'x':
            run(streamer)
            return -1
        if action == 's':
            streamer.save(input("Enter save path:"))
    while True:
        print("Change color map? (y/n)")
        action = input()
        if action == 'n':
            break
        if action == 'y':
            pass



def run(streamer):
    """
    Run the configuration editor
    :param streamer: Streamer object to modify
    :return:
    """
    print("Streamer - Stream plot printing and clustering tool")
    cur_opt = ['l','c']
    bin_opt = ['y', 'n']
    print("This program will help you create a configuration file for later use")
    while True:
        action = input("to load a streamer configuration (\".SCG\") file using press \'l\', Alternatively press \'c\' to construct a new one")
        if action.lower() == 'l':
            while True:
                path = input("Enter file path or 'x' to return:")
                if path.lower() == 'x':
                    break
                if not os.path.isfile(path.strip()):
                    print("Could not find path")
                    continue
                if not path.endswith("SCG"):
                    print("Invalid file type")
                    continue
                with open(path, 'rb') as f:
                    streamer = pickle.load(f)
                    if construct_streamer(streamer) < 0:
                        run(streamer)
                        return
                    else:
                        return
        if action.lower == 'c':
            construct_streamer(streamer)
            break
        if action.lower() == 'x':
            break



def main():
    """
    Run this program without arguments to configure a new streamer object
    :return:
    """
    if len(sys.argv) <= 1:
        print("No arguments passed, running the StreamPlot configuration edtior")
        streamer = StreamPlot()
        run(streamer)
    else:
        parser = ArgumentParser()
        parser.add_argument(['-cfg', '-configuration'], help='Path to StreamPlot Configuration file (.SCG)', required=True, dest='cfg')
        parser.add_argument(['-s','-sub','-subject_matrix'], help='Path to subject microbiome reads', dest='sub')
        parser.add_argument(['-d','-delim'], help='Optional: Subject matrix delimiter (default is \",\")', dest='delim1')
        parser.add_argument(['-m', '-meta'], help='Path to Meta-data file (as csv, xls, xlsx)', dest='meta')
        parser.add_argument(['-m_d, -meta_delimiter'], dest='delim2', help='Metadata delimiter (default is \",\")')
        parser.add_argument(['-m_s', '-meta_sheet'], dest='sheet', help='Metadata sheet name (for Excel files)')
        parser.add_argument(['-c','-cluster'], nargs=4, help=f'Cluster subjects and plot clusters, parameters:{os.linesep}Number of dominant bacteria to consider{os.linesep}Matching threshold{os.linesep}Time threshold{os.linesep}Delta', dest='c')
        parser.add_argument(['-p','-plot'], help='Plot all inserted subjects', dest='p')
        parser.add_argument(['-f','-folder'], dest='f', help='Folder to save to, default is \"Plots\"')
        parser.add_argument(['-e', '-events'], dest='events', nargs='*', help=f'Add events from the meta data. Optional arguments (as ordered):{os.linesep}Event name (default \"Onset\"){os.linesep}Meta data subject ID column (default \"record_id\"){os.linesep}Event columns, Time column and events identifier (e.g. \"Symptomatic\", \"Case\" etc.)')
        args = parser.parse_args()
        if args.cfg.endswith(".SCG"):
            try:
                streamer = pickle.load(args.cfg)
            except:
                raise RuntimeError("File not found")
        else:
            raise RuntimeError("Cfg file must be of type .SCG")
        if args.sub:
            delim  = args.delim1 if args.delim1 else ","
            streamer.load_subjects(args.sub, delim)
        if args.meta:
            delim = args.delim2 if args.delim2 else ","
            sheet = args.sheet if args.sheet else "Metadata"
            streamer.add_meta(args.meta, delim, sheet)
        folder = args.f if args.f else 'Plots'
        if args.c:
            streamer.cluster(int(args.c[0]),int(args.c[1]),int(args.c[2]),
                            float(args.c[3]), folder)
        if args.p:
            streamer.plot_stream(folder=folder)
        if args.events:
            # TODO use streamer.events_from_meta() with the args
            # from args.events to implement this
            pass


if __name__ == '__main__':
    main()