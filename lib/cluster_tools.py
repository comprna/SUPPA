# -*- coding: utf-8 -*-
"""
Created on Wed May 25 04:20:00 CEST 2016

@authors: Juan C Entizne, Juan L. Trincado
@email: juancarlos.entizne01[at]estudiant.upf.edu,
        juanluis.trincado[at]upf.edu
"""

import os
import logging
from itertools import islice
from sklearn.cluster import DBSCAN
from collections import defaultdict
from sklearn.metrics import silhouette_score
from lib.optics import *



def sliding_windows(seq, n=2):
    '''
    Returns a sliding sliding_windows (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
    '''

    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def is_index_valid(indexes, lst):

    # Setting logging preferences
    logger = logging.getLogger(__name__)

    last_el = int(indexes[-1][-1])
    lst_len = int(len(lst))
    if last_el < len(lst):
        logger.error("Invalid index. Index %d is smaller than the number of columns in the file (%d)."
                     % (last_el, lst_len))
        sys.exit(1)
    elif last_el > lst_len:
        logger.error("Invalid index. Index %d is larger than the number of columns in the file (%d)."
                     % (last_el, lst_len))
        sys.exit(1)
    else:
        pass

    # If the groups of indexes are continuous the diff between the value of the last element of the first
    # indexes_group minus the first element of the next index_group must be equal to 1
    for indx_pair in sliding_windows(indexes):
        diff = indx_pair[1][0]-indx_pair[0][-1]

        if diff <= 0:
            msg = str(indx_pair[0][0])+"-"+str(indx_pair[0][-1])+", "+str(indx_pair[1][0])+"-"+str(indx_pair[1][-1])
            logger.error("Indexes %s are not valid. Cluster groups must not overlap." % msg)
            sys.exit(1)
        elif diff > 1:
            msg = str(indx_pair[0][0])+"-"+str(indx_pair[0][-1])+", "+str(indx_pair[1][0])+"-"+str(indx_pair[1][-1])
            logger.error("Indexes %s are not valid. Cluster groups must be continuous." % msg)
            sys.exit(1)
        else:
            pass

    return True


def average_psi_per_condition(psi_vals, indexes):

    # Setting logging preferences
    logger = logging.getLogger(__name__)

    new_vector = []

    # Check if there is empty elements in index list
    if '' in indexes.split(","):
        logger.error("One of the indexes group is not valid.")
        sys.exit(1)

    grouped_ranges = [range_pair.split("-") for range_pair in indexes.split(",")]
    group_indexes = [list(map(int, ranges)) for ranges in grouped_ranges]

    if is_index_valid(group_indexes, psi_vals):
        for indx in group_indexes:
            frst = indx[0]-1
            last = indx[-1]
            group = np.mean(psi_vals[frst:last])

            new_vector.append(group)

    return new_vector


def process_cluster_input(dpsi, psivec, sig_threshold, dpsi_threshold, indexes):

    sig_events = []
    with open(dpsi) as fh:
        next(fh)
        for ln in fh:
            sl = ln.strip("\n").split("\t")
            ev_id = sl[0]
            ev_dpsi_vals = list(map(abs, map(float, sl[1::2])))
            ev_pvals = list(map(float, sl[2::2]))
            for dpsi_val, p_val in zip(ev_dpsi_vals, ev_pvals):
                if dpsi_val > dpsi_threshold and p_val < sig_threshold:
                    sig_events.append(ev_id)

    psi_matrix, ev_id_lst = ([] for _ in range(2))
    with open(psivec) as fh:
        next(fh)
        for ln in fh:
            sl = ln.strip("\n").split("\t")
            ev_id = sl[0]
            psi_vals_str = sl[1:]
            if ev_id in sig_events:
                # Strict form to filter events: event excluded if -1.0 in any sample
                if "-1.0" not in psi_vals_str:
                    psi_vals = [float(val) for val in psi_vals_str]
                    if indexes:
                        vals = average_psi_per_condition(psi_vals, indexes)
                    else:
                        vals = psi_vals

                    psi_matrix.append(vals)
                    ev_id_lst.append(ev_id)

    #psi_matrix is the average psi per condition (only for events that at least one of the comparisons yielded a significant comparison)
    #ev_id_lst is just the list of the events included in the previous analysis
    return psi_matrix, ev_id_lst


def DBSCAN_cluster(psi_matrix, eventid_lst, dist, minpts, metric):

    # Setting logging preferences
    logger = logging.getLogger(__name__)

    # The metric is "cosine" works only with the algorithm "brute"
    if metric == "cosine":
        alg = 'brute'
    else:
        alg = 'auto'

    try:
        db = DBSCAN(eps=dist, min_samples=minpts, metric=metric, algorithm=alg).fit(psi_matrix)
        labels = db.labels_
    except:
        logger.error("Unknown error: {}".format(sys.exc_info()))
        sys.exit(1)

    eventid_labels_dict = {k: v for k, v in zip(eventid_lst, labels)}

    return eventid_labels_dict, labels


def write_averaged_cluster_output(psi_matrix, eventid_lst, eventid_labels_dict, output):

    res = []
    for ev_id, avg in zip(eventid_lst, psi_matrix):
        avg_vals = "\t".join(list(map(str, avg)))
        cluster_id = eventid_labels_dict[ev_id]
        el = (ev_id, cluster_id, avg_vals)
        res.append(el)

    res.sort(key=lambda x: x[1], reverse=True)

    with open("%s.clustvec" % output, "w+") as fh:
        for el in res:
            ln = "%s\t%d\t%s\n" % (el[0], el[1], el[2])
            fh.write(ln)


def calculate_rmsstd(cluster_vectors):

    n_cond = len(cluster_vectors[0])
    # Degrees of Freedom
    dof = n_cond-1

    return np.sqrt(sum(np.std(cluster_vectors, axis=0))/dof)


def calculate_cluster_scores(x, cluster_labels, output):

    with open("%s_scores.log" % output, "w+") as fh:
        # Filter out singleton "cluster" (labeled as -1)
        filtered_x, filtered_cluster_labels, singletons = ([] for _ in range(3))
        cluster_groups = defaultdict(list)
        for vec, lab in zip(x, cluster_labels):
            if not lab == -1:
                filtered_x.append(vec)
                filtered_cluster_labels.append(lab)

                cluster_groups[lab].append(vec)
            else:
                singletons.append(vec)

        ln = "Number of clustered events: %d/%d (%f%%)\n" % (len(filtered_x), len(filtered_x)+len(singletons),
                                                           (len(filtered_x)/(len(filtered_x)+len(singletons)))*100)
        print(ln.strip("\n"))
        fh.write(ln)

        for group in cluster_groups:
                n_events = len(cluster_groups[group])
                ln = "Cluster %d contains %d events\n" % (group, n_events)
                print(ln.strip("\n"))
                fh.write(ln)

        rmsstd_scores = []
        for group in cluster_groups:
            rmsstd = calculate_rmsstd(np.array(cluster_groups[group]))
            ln = "The RMSSTD score for cluster %d is %f\n" % (group, rmsstd)
            print(ln.strip("\n"))
            fh.write(ln)

            rmsstd_scores.append(rmsstd)

        try:
            silhouette_avg = silhouette_score(np.array(filtered_x), np.array(filtered_cluster_labels))
            ln = "The average silhouette score is : %f\n" % silhouette_avg
            print(ln.strip("\n"))
            fh.write(ln)
        except:
            silhouette_avg = float("nan")
            ln = "Impossible to calculate silhouette score. Only 1 cluster group identified.\n"
            print(ln.strip("\n"))
            fh.write(ln)

    return silhouette_avg, rmsstd_scores


def create_points_list(psi_matrix, eventid_lst):
    '''Generate a list of points from the psi_matrix'''
    points_list = []
    for i,x in enumerate(psi_matrix):
        points_list.append(Point(eventid_lst[i],x))
    return points_list

def generate_labels(clusters, eventid_lst):
    '''From the clusters generated from OPTICS, label the event_ids'''
    eventid_labels_dict = {}
    # Initialize to -1 both structures
    for x in eventid_lst : eventid_labels_dict[x] = -1
    labels = [-1 for x in range(len(eventid_lst))]

    #Read clusters, annotating the cluster number for each event
    if(len(clusters)!=0):
        for number, cluster in enumerate(clusters):
            for point in cluster.points:
                eventid_labels_dict[point.id] = number
                labels[eventid_lst.index(point.id)] = number

    return eventid_labels_dict, labels


def cluster_analysis(dpsi, psivec, sig_threshold, dpsi_threshold, eps, minpts, metric, indexes, clustering,
                     separation, output):

    path = os.path.dirname(os.path.realpath(dpsi))
    os.chdir(path)

    psi_matrix, eventid_lst = process_cluster_input(dpsi, psivec, sig_threshold, dpsi_threshold, indexes)

    if(clustering=="DBSCAN"):
        eventid_labels_dict, labels = DBSCAN_cluster(psi_matrix, eventid_lst, eps, minpts, metric)
        #eventid_labels_dict are the labels of the clustering for eacg event

        write_averaged_cluster_output(psi_matrix, eventid_lst, eventid_labels_dict, output)
        calculate_cluster_scores(psi_matrix, labels, output)

    else:
        #OPTICS
        points_list = create_points_list(psi_matrix, eventid_lst) #Transform the points on psi_matrix to Points from optics.py
        optics = Optics(points_list, eps, minpts)  # Maximum radius to be considered, cluster size >= 2 points
        optics.run()  # run the algorithm
        clusters = optics.cluster(separation)  # minimum threshold for clustering (upper limit to separate the clusters)
        eventid_labels_dict, labels = generate_labels(clusters, eventid_lst)
        write_averaged_cluster_output(psi_matrix, eventid_lst, eventid_labels_dict, output)
        calculate_cluster_scores(psi_matrix, labels, output)


