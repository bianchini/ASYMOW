import ROOT
import os
import re
import argparse
import numpy as np
import csv
import math


'''
# Dump BeamSpot infos from DB:
> wget https://raw.githubusercontent.com/cms-sw/cmssw/master/CondTools/BeamSpot/test/BeamSpotRcdPrinter_cfg.py .
> cmsRun BeamSpotRcdPrinter_cfg.py  # for 2016: --startIOV 1173204676640769 --endIOV 1219959690625054
The first and last IOV are determined by looking at the json:
/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt
using the unopack functions. E.g.:
root[0] #include "CondFormats/Common/interface/TimeConversions.h"
root[1] cond::time::pack({284044, 30})

# Produce a brilcalc output file containing a summary of lumi and avg pileup LS-by-LS:
> setenv PATH $HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
> pip install --user --upgrade brilws
> ~/.local/bin/brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt --byls --minBiasXsec 69200 --xingTr 0.1 > 2016_brilcalc_byls_69200ub.txt
# FIXME: brilcalc supports the -o option to output a csv file directly
'''

parser = argparse.ArgumentParser("")
parser.add_argument('-e', '--era', type=str, default='2016', help="")
parser.add_argument('-t', '--tag', type=str, default='69200ub', help="")
parser.add_argument('-rS', '--reduceSummary',  action='store_true', default=False, help="")
parser.add_argument('-rB', '--reduceBrilCalc', action='store_true', default=False, help="")
parser.add_argument('-rM', '--reduceMatch',    action='store_true', default=False, help="")
parser.add_argument('-rF', '--reduceFinal',    action='store_true', default=False, help="")
parser.add_argument('-mT', '--makeTree', action='store_true', default=False, help="")
parser.add_argument('-mP', '--makePdf', action='store_true', default=False, help="")
parser.add_argument('-mPP','--makePlot', action='store_true', default=False, help="")
parser.add_argument('-mPP2','--makePlot2', action='store_true', default=False, help="")
args = parser.parse_args()
era = args.era
tag = '_'+args.tag
reduceS = args.reduceSummary
reduceB = args.reduceBrilCalc
reduceM = args.reduceMatch
reduceF = args.reduceFinal
makeTree= args.makeTree
makePdf= args.makePdf
makePlot= args.makePlot
makePlot2= args.makePlot2

#tag = '_69200ub'
#tag = '_72400ub'
#tag = '_66000ub'

Run2016F_nonHIPM = ['278769', '278801', '278802', '278803', '278804', '278805', '278808']
Run2016G_first   = ['278820']


def reduce_summary():
    print 'reduce_summary()...'
    fout = open('summary_2016_reduced.csv', 'w')
    bs_file = open('summary_2016.txt', 'r')
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

def check_era(run,era):
    is_postVFP = (int(run)>=int(Run2016G_first[0]) or (run in Run2016F_nonHIPM))
    if era=='2016': 
        return True
    elif era=='2016preVFP' and is_postVFP: return False
    elif era=='2016postVFP' and not is_postVFP: return False
    return True

def reduce_brilcalc(tag='', era=''):
    print 'reduce_brilcalc()...'
    fout = open(era+'_brilcalc_byls_reduced'+tag+'.csv', 'w')
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
        if(check_era( run, era)):
            fout.write(run+','+ls+','+lumi+','+avg+'\n')
    fout.close()
    return


def reduce_match(tag='',era=''):
    print 'reduce_match()...'
    with open('summary_2016_reduced.csv', 'rb') as master:
        master_indices = dict( (r[0]+'_'+r[1], [ r[2],r[3],   #X0
                                                 r[4],r[5],   #Y0
                                                 r[6],r[7],   #Z0
                                                 r[8],r[9],   #sZ0
                                                 r[10],r[11], #dxdz
                                                 r[12],r[13]  #dydz                                             
                                             ] ) for i,r in enumerate(csv.reader(master)))
    with open(era+'_brilcalc_byls_reduced'+tag+'.csv', 'rb') as hosts:
        with open(era+'_brilcalc_byls_matched'+tag+'.csv', 'wb') as results:
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

def reduce_final(tag='',era=''):
    print 'reduce_final()...'
    fout = open('summary_'+era+'_final'+tag+'.csv', 'w')
    ls_file = open(era+'_brilcalc_byls_matched'+tag+'.csv', 'r')
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
    print Ltot,'/fb in summary_'+era+'_final.csv'
    fout.close()
    return


def check(tag='',era=''):
    ls_file = open(era+'_brilcalc_byls_matched'+tag+'.csv', 'rb')
    bs_file = open('summary_2016_reduced.csv', 'rb')
    ls_lines = [ x.strip() for x in ls_file.readlines() ]
    bs_lines = [ x.strip() for x in bs_file.readlines() ]
    unmatched = 0
    dL = 0.
    for l in ls_lines:
        dL += float(l.split(',')[2])
    print "Total lumi in "+era+"_brilcalc_byls_matched.csv:", dL
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

def make_tree(tag='',era=''):
    print 'make_tree()...'
    f = ROOT.TFile('out_'+era+tag+'.root','RECREATE')
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
 
    ls_file = open('summary_'+era+'_final'+tag+'.csv', 'r')
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

def make_tree_avg(tag='',era=''):
    print 'make_tree_avg()...'
    f = ROOT.TFile('out_'+era+'_avg'+tag+'.root','RECREATE')
    tree = ROOT.TTree('tree', 'tree')
    dL = np.empty((1), dtype="float32")
    avg = np.empty((1), dtype="float32")
    tree.Branch("dL", dL, "dL/F")
    tree.Branch("avg", avg, "avg/F")
    ls_file = open(era+'_brilcalc_byls_matched'+tag+'.csv', 'r')
    ls_lines = [ x.strip() for x in ls_file.readlines() ]
    for l in ls_lines:
        dL[0] = l.split(',')[2]
        avg[0] = l.split(',')[3]
        tree.Fill()
    tree.Write()
    f.Close()
    return

def make_pdfs(tag='',era=''):
    print 'make_pdfs()...'

    zbins, zmin, zmax = 55, -15., 15.

    mc = ROOT.TFile('./NanoAOD.root','READ')
    tmc = mc.Get('demo/Events')
    hmc = ROOT.TH1F('mc_vtxZ','', zbins, zmin, zmax)
    tmc.Draw('GenVertex_z>>mc_vtxZ') 
    fitmc = ROOT.TF1('fitmc', '[0]/TMath::Sqrt(2*TMath::Pi())/[1]*TMath::Exp( -0.5*(x-[2])*(x-[2])/[1]/[1] )', zmin, zmax)    
    fitmc.SetParameter(0, hmc.Integral())
    fitmc.SetParameter(1, hmc.GetRMS())
    fitmc.SetParameter(2, hmc.GetMean())
    hmc.Fit(fitmc, 'R')
    fitmc.SetParameter(0, 1.0)
    hmc.Scale(1./hmc.Integral())

    f = ROOT.TFile('out_'+era+tag+'.root','READ')
    tree = f.Get('tree')

    hdict = {}
    keys = ['0_50', '5_10', '10_15', '15_20', '20_25',  '25_30', '30_35', '35_40', '40_45']
    for k in keys:
        hdict[k] = ROOT.TH1F('vtxZ_'+k,k+';vtx_z;pdf',zbins, zmin,zmax)
    h2dict = {}
    h2dict['sZ0_avg'] = ROOT.TH2F('sZ0_avg',';avg;sZ0', 7, 5, 40, 20, 3.0, 4.0)
    h2dict['Z0_avg']  = ROOT.TH2F('Z0_avg', ';avg;Z0',  7, 5, 40, 20, -1.0, 1.0)

    dL = np.empty((1), dtype="float32")
    avg = np.empty((1), dtype="float32")
    Z0 = np.empty((1), dtype="float32")
    sZ0 = np.empty((1), dtype="float32")
    tree.SetBranchAddress('dL', dL)
    tree.SetBranchAddress('avg', avg)
    tree.SetBranchAddress('Z0', Z0)
    tree.SetBranchAddress('sZ0', sZ0)
    for i in range(0, tree.GetEntries()):
        tree.GetEntry(i)
        for k,v in hdict.items():
            for b in range(1, v.GetNbinsX()+1):
                z = v.GetBinCenter(b)
                val = ROOT.TMath.Gaus(z, Z0[0], sZ0[0], True )
                avg_min, avg_max = float(k.split('_')[0]),float(k.split('_')[1]) 
                if avg[0]>=avg_min and avg[0]<avg_max:
                    v.Fill( z, val*dL[0] )
        h2dict['sZ0_avg'].Fill(avg[0],sZ0[0],dL[0])
        h2dict['Z0_avg'].Fill(avg[0], Z0[0],dL[0])

    for k,v in hdict.items():
        norm = v.Integral('width')
        if norm>0.:
            v.Scale(1./norm)
      
    fout = ROOT.TFile('ratio_'+era+tag+'.root', 'RECREATE')
    fout.cd()
    hratiodict = {}
    for k,v in hdict.items():
        v.Write()
        hratiodict[k] = v.Clone('ratio_'+k)        
        for b in range(1, v.GetNbinsX()+1):
            #hratiodict[k].SetBinContent(b, v.GetBinContent(b)/hmc.GetBinContent(b) if hmc.GetBinContent(b)>0. else 1.0)    
            vmc = fitmc.Eval(v.GetBinCenter(b))
            hratiodict[k].SetBinContent(b, v.GetBinContent(b)/vmc if vmc>0. else 1.0)    
    for k,v in hratiodict.items():
        v.Write()
    for k,v in h2dict.items():
        v.Write()
    hmc.Write()
    fout.Close()


def plot(tag,era=''):    
    from ROOT import kRed, kDashed
    print 'plot()...'
    fin = ROOT.TFile('ratio_'+era+tag+'.root', 'READ')
    keys = ['0_50', '5_10', '10_15', '15_20', '20_25',  '25_30', '30_35', '35_40', '40_45']
    c = ROOT.TCanvas()
    leg = ROOT.TLegend(0.7,0.5,0.9,0.9)
    leg.SetHeader('Average pileup')
    hmc = fin.Get('mc_vtxZ')
    h0 = None
    for ik,k in enumerate(keys):
        h = fin.Get('ratio_'+k)
        # compute mean and RMS of weight
        mu,rms2 = 0., 0.        
        for ib in range(1, hmc.GetNbinsX()+1):            
            mu += h.GetBinContent(ib)*hmc.GetBinContent(ib)
        for ib in range(1, hmc.GetNbinsX()+1):            
            rms2 += ((h.GetBinContent(ib)-mu)**2)*hmc.GetBinContent(ib)          
        print mu, math.sqrt(rms2)
        
        h.SetLineColor(ik+1)
        h.SetLineWidth(2)
        h.SetStats(0)
        leg.AddEntry(h, '['+k.split('_')[0]+', '+k.split('_')[1]+') #Delta='+'{:.1f}'.format(rms2*0.5*100)+'%', "L")
        if ik==0: 
            h0 = h
            h.SetLineWidth(4)
            h.SetLineStyle(kDashed)
            h.SetMaximum(3.0)
            h.SetMinimum(0.0)
            h.SetTitle("")
            h.GetXaxis().SetTitle("Z_{vtx} [cm]")
            h.GetXaxis().SetTitleSize(0.05)            
            h.GetYaxis().SetTitleSize(0.05)            
            h.GetXaxis().SetTitleOffset(0.9)
            h.GetYaxis().SetTitle("Data/MC")
            h.Draw("HIST")
        else: 
            h.Draw("HISTSAME")
    h0.SetTitle(era)
    h0.Draw("SAME")
    c.Update()
    leg.Draw()
    c.Draw()
    c.SaveAs('ratioVsPU_'+era+tag+'.png')
    #raw_input()


def plot2D(tag,era=''):    
    print 'plot2D()...'
    fin = ROOT.TFile('ratio_'+era+tag+'.root', 'READ')
    keys = ['sZ0_avg','Z0_avg']
    c = ROOT.TCanvas()
    for k in keys:
        h2 = fin.Get(k)
        h2.SetTitle(era)
        c.cd()
        h2.Draw('COLZ')
        c.SaveAs(k+'_'+era+'.png')
    fin.Close()
    return
    

if reduceS:
    reduce_summary()
if reduceB:
    reduce_brilcalc(tag,era)
if reduceM:
    reduce_match(tag,era)
if reduceF:
    reduce_final(tag,era)
if makeTree:
    make_tree(tag,era)
#check()
#make_tree_avg(tag,era)
if makePdf:
    make_pdfs(tag,era)
if makePlot:
    plot(tag,era)
if makePlot2:
    plot2D(tag,era)
