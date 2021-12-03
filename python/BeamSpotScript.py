import ROOT
import os
import re
import argparse
import numpy as np
import csv

parser = argparse.ArgumentParser("")
parser.add_argument('-rB', '--reduceBrilCalc', action='store_true', default=False, help="")
parser.add_argument('-rS', '--reduceSummary',  action='store_true', default=False, help="")
parser.add_argument('-rM', '--reduceMatch',    action='store_true', default=False, help="")
parser.add_argument('-rF', '--reduceFinal',    action='store_true', default=False, help="")
parser.add_argument('-t', '--tree', action='store_true', default=False, help="")
args = parser.parse_args()
reduceB = args.reduceBrilCalc
reduceS = args.reduceSummary
reduceM = args.reduceMatch
reduceF = args.reduceFinal
makeTree = args.tree

def reduce_brilcalc(tag=''):
    print 'reduce_brilcalc()...'
    fout = open('2016_brilcalc_byls_reduced'+tag+'.csv', 'w')
    ls_file = open('2016_brilcalc_byls'+tag+'.txt', 'r')
    ls_lines = [ x.strip() for x in ls_file.readlines()]
    for l in ls_lines:
        #print l
        if '|' not in l: continue
        if 'STABLE' not in l: continue
        l2 = l.split('|')
        run  = l2[1].rstrip().lstrip().split(':')[0]
        ls   = l2[2].rstrip().lstrip().split(':')[0]
        lumi = l2[7].rstrip().lstrip()
        avg  = l2[8].rstrip().lstrip()
        #print('('+run+','+ls+') '+lumi)
        fout.write(run+','+ls+','+lumi+','+avg+'\n')
    fout.close()
    return

def reduce_summary():
    print 'reduce_summary()...'
    fout = open('summary_reduced.csv', 'w')
    bs_file = open('summary.txt', 'r')
    bs_lines = [ x.strip() for x in bs_file.readlines() ]
    for ib,b in enumerate(bs_lines):
        if 'hash' in b:
            b1   = bs_lines[ib].split(',')
            run  = b1[0]
            ls = b1[1].split(' ')[0]
            target = bs_lines[ib+5:ib+13]
            X0,X0err     = target[0].split(' ')[6], target[0].split(' ')[8]
            Y0,Y0err     = target[1].split(' ')[6], target[1].split(' ')[8]
            Z0,Z0err     = target[2].split(' ')[6], target[2].split(' ')[8]
            sZ0,sZ0err   = target[3].split(' ')[7], target[3].split(' ')[9]
            dxdz,dxdzerr = target[4].split(' ')[10], target[4].split(' ')[12]
            dydz,dydzerr = target[5].split(' ')[10], target[5].split(' ')[12]
            fout.write(run+','+ls+','+X0+','+X0err+','+Y0+','+Y0err+','+Z0+','+Z0err+','+sZ0+','+sZ0err+','+dxdz+','+dxdzerr+','+dydz+','+dydzerr+'\n')
    fout.close()
    return

def reduce_match(tag=''):
    print 'reduce_match()...'
    with open('summary_reduced.csv', 'rb') as master:
        master_indices = dict( (r[0]+'_'+r[1], [ r[2],r[3],   #X0
                                                 r[4],r[5],   #Y0
                                                 r[6],r[7],   #Z0
                                                 r[8],r[9],   #sZ0
                                                 r[10],r[11], #dxdz
                                                 r[12],r[13]  #dydz                                             
                                             ] ) for i,r in enumerate(csv.reader(master)))
    with open('2016_brilcalc_byls_reduced'+tag+'.csv', 'rb') as hosts:
        with open('2016_brilcalc_byls_matched'+tag+'.csv', 'wb') as results:
            reader = csv.reader(hosts)
            writer = csv.writer(results)
            for row in reader:
                key = row[0]+'_'+row[1]
                index = master_indices.get(key)
                if index is not None:
                    writer.writerow(row + 
                                    [master_indices[key][0]] + [master_indices[key][1]] + #X0 
                                    [master_indices[key][2]] + [master_indices[key][3]] + #Y0
                                    [master_indices[key][4]] + [master_indices[key][5]] + #Z0
                                    [master_indices[key][6]] + [master_indices[key][7]] + #sZ0
                                    [master_indices[key][8]] + [master_indices[key][9]] + #dxdz
                                    [master_indices[key][10]] + [master_indices[key][11]] #dydz
                    )
                else:
                    writer.writerow(row)
    return

def reduce_final(tag=''):
    print 'reduce_final()...'
    fout = open('summary_final'+tag+'.csv', 'w')
    ls_file = open('2016_brilcalc_byls_matched'+tag+'.csv', 'r')
    ls_lines = [ x.strip() for x in ls_file.readlines() ]
    n_ls_lines = len(ls_lines)
    Ltot = 0.
    for il in range(0,len(ls_lines)):
        l = ls_lines[il]
        if l.count(',')>3:
            dL = float( l.split(',')[2] )
            avg = float( l.split(',')[3] )*dL
            il2 = il+1
            found_next = False
            while (not found_next) and il2<n_ls_lines:
                lpost = ls_lines[il2]
                if lpost.count(',') == 3:
                    lumi = lpost.split(',')[2]
                    avgi = float(lpost.split(',')[3])
                    dL += float(lumi)
                    avgi *= float(lumi)
                    avg += avgi
                    il2 += 1
                else: found_next = True
            avg /= (dL if dL>0. else 1.0) 
            Ltot += dL
            l_split = l.split(',')
            l_split[2] = str(dL)
            l_split[3] = str(avg)
            newline = ','.join(l_split)
            #print(newline)
            fout.write(newline+'\n')
    print Ltot,'/fb in summary_final.csv'
    fout.close()
    return


def check(tag=''):
    ls_file = open('2016_brilcalc_byls_matched'+tag+'.csv', 'rb')
    bs_file = open('summary_reduced.csv', 'rb')
    ls_lines = [ x.strip() for x in ls_file.readlines() ]
    bs_lines = [ x.strip() for x in bs_file.readlines() ]
    unmatched = 0
    dL = 0.
    for l in ls_lines:
        dL += float(l.split(',')[2])
    print "Total lumi in 2016_brilcalc_byls_matched.csv:", dL
    for b in bs_lines:
        run,ls = b.split(',')[0], b.split(',')[1]
        match = False
        for l in ls_lines:
            run2,ls2 = l.split(',')[0], l.split(',')[1]
            if run2==run and ls2==ls:
                match = True
                break
        if not match:
            print run,ls, "NOT MATCHED"
            unmatched += 1
    print unmatched, 'unmatched'
    return

def make_tree(tag=''):
    print 'make_tree()...'
    f = ROOT.TFile('out'+tag+'.root','RECREATE')
    tree = ROOT.TTree('tree', 'tree')
    dL = np.empty((1), dtype="float32")
    avg = np.empty((1), dtype="float32")
    X0 = np.empty((1), dtype="float32")
    X0err = np.empty((1), dtype="float32")
    Y0 = np.empty((1), dtype="float32")
    Y0err = np.empty((1), dtype="float32")
    Z0 = np.empty((1), dtype="float32")
    Z0err = np.empty((1), dtype="float32")
    sZ0 = np.empty((1), dtype="float32")
    sZ0err = np.empty((1), dtype="float32")
    dxdz = np.empty((1), dtype="float32")
    dxdzerr = np.empty((1), dtype="float32")
    dydz = np.empty((1), dtype="float32")
    dydzerr = np.empty((1), dtype="float32")

    tree.Branch("dL", dL, "dL/F")
    tree.Branch("avg", avg, "avg/F")
    tree.Branch("X0", X0, "X0/F")
    tree.Branch("X0err", X0err, "X0err/F")
    tree.Branch("Y0", Y0, "Y0/F")
    tree.Branch("Y0err", Y0err, "Y0err/F")
    tree.Branch("Z0", Z0, "Z0/F")
    tree.Branch("Z0err", Z0err, "Z0err/F")
    tree.Branch("sZ0", sZ0, "sZ0/F")
    tree.Branch("sZ0err", sZ0err, "sZ0err/F")
    tree.Branch("dxdz", dxdz, "dxdz/F")
    tree.Branch("dxdzerr", dxdzerr, "dxdzerr/F")
    tree.Branch("dydz", dydz, "dydz/F")
    tree.Branch("dydzerr", dydzerr, "dydzerr/F")

    ls_file = open('summary_final'+tag+'.csv', 'r')
    ls_lines = [ x.strip() for x in ls_file.readlines() ]
    for l in ls_lines:
        dL[0] = l.split(',')[2]
        avg[0] = l.split(',')[3]
        X0[0] = l.split(',')[4]
        X0err[0] = l.split(',')[5]
        Y0[0] = l.split(',')[6]
        Y0err[0] = l.split(',')[7]
        Z0[0] = l.split(',')[8]
        Z0err[0] = l.split(',')[9]
        sZ0[0] = l.split(',')[10]
        sZ0err[0] = l.split(',')[11]
        dxdz[0] = l.split(',')[12]
        dxdzerr[0] = l.split(',')[13]
        dydz[0] = l.split(',')[14]
        dydzerr[0] = l.split(',')[15]
        tree.Fill()
    tree.Write()
    f.Close()
    return

def make_tree_avg(tag=''):
    print 'make_tree_avg()...'
    f = ROOT.TFile('out_avg'+tag+'.root','RECREATE')
    tree = ROOT.TTree('tree', 'tree')
    dL = np.empty((1), dtype="float32")
    avg = np.empty((1), dtype="float32")
    tree.Branch("dL", dL, "dL/F")
    tree.Branch("avg", avg, "avg/F")
    ls_file = open('2016_brilcalc_byls_matched'+tag+'.csv', 'r')
    ls_lines = [ x.strip() for x in ls_file.readlines() ]
    for l in ls_lines:
        dL[0] = l.split(',')[2]
        avg[0] = l.split(',')[3]
        tree.Fill()
    tree.Write()
    f.Close()
    return

tag = '_69200ub'

if reduceB:
    reduce_brilcalc(tag)
if reduceS:
    reduce_summary()
if reduceM:
    reduce_match(tag)
if reduceF:
    reduce_final(tag)
if makeTree:
    make_tree(tag)

#check()
make_tree_avg(tag)

#make_tree(firstLS, lastLS)
