"""
use case:
jkml_v2.synthetic_data(50,100,8,lambda x: x.find('TGCTGCTA'),lambda: 'TGCTGCTA',freq=0.5,signal=2)
"""

import os, time, math
import random as rd
import numpy as np
import scipy
import re
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
from pylab import *


####################################################################################################################
# Kmer model
####################################################################################################################

bases = ['A','T','C','G']

def add_base(li):
        """Used in make_mer_list to add one more base to list"""
        new_li = []
        for s in li:
            for b in bases:
                new_li.append(s+b)
        return new_li

def make_mer_list(mer_len):
    """Makes a list of all n-mers"""
    li = bases
    for i in range(mer_len-1):
        li = add_base(li)
    return li

def kmer_matrix(seqs,mer_len):

    mer_dict = dict(zip(make_mer_list(mer_len),range(4**mer_len)))
    rows,cols = [],[]
    r = 0
    for i in xrange(len(seqs)):
        cur_seq = seqs[i]
        for b in range(len(cur_seq)-mer_len+1):
            rows.append(r)
            cols.append(mer_dict[cur_seq[b:b+mer_len]])
        # if(r%10000)==0:
        #     print(r),
        r+=1
    vals = np.ones_like(cols)
    rows.append(r-1)
    cols.append(4**mer_len-1)
    vals = np.append(vals,0)
    X = scipy.sparse.csr_matrix((vals,(rows,cols)),dtype=np.float64)
    return X

def kmer_matrix_pos(seqs,mer_len,pos_win):

    seq_series = pd.Series(seqs)

    result = []

    for i in range(mer_len-pos_win+1):
        result.append(kmer_matrix(seq_series.str.slice(i,i+pos_win).values,pos_win))

    return scipy.sparse.hstack(result)

####################################################################################################################
# Synthetic data generation
####################################################################################################################

def gen_random_seq(length, rna=False):
    if rna:
        return ''.join(rd.choice('CGUA') for _ in xrange(length))
    else:
        return ''.join(rd.choice('CGTA') for _ in xrange(length))


def random_dna_list(length, count=1, rna=False):

    if type(length) == int:
        return [gen_random_seq(length, rna) for i in range(count)]
    else:
        return [gen_random_seq(rd.randint(length[0],length[1]), rna) for i in range(count)]


def decision(probability):
    return rd.random() < probability


def spike_in_motif(seqL, motif_len, motif_detector, motif_generator, freq, rna=False):
    resultL = []

    for s in seqL:

        while 1:
            idx = motif_detector(s)
            if idx == -1:
                break
            s = s[:idx] + gen_random_seq(motif_len, rna) + s[(idx + motif_len):]

        if decision(freq):
            s_new = bytearray(s)
            posSta = rd.randint(0, len(s) - motif_len)
            s_new[posSta:posSta + motif_len] = motif_generator()
            s = str(s_new)

        resultL.append(s)

    return resultL


def synthetic_data(seq_len, seq_count, motif_len, motif_detector, motif_generator, freq=0.5):
    """
    only supports:
    - 0 or 1 motifs per sequence
    - fixed length motif
    """

    print('synthesizing random sequences')
    seqL_raw = random_dna_list(seq_len, seq_count)

    print('spiking in designated motif')
    seqL = spike_in_motif(seqL_raw, motif_len, motif_detector, motif_generator, freq)

    return seqL


def get_sec_structure(seqL):

    tempfilepath = '/tmp/%s' % time.time()

    with open(tempfilepath,'w') as tempfile:
        tempfile.write('\n'.join(seqL))

    secL = []

    with os.popen('/home/jinkuk/.local/bin/RNAfold --noPS -4 --maxBPspan=24 --noconv < %s' % tempfilepath, 'r') as secF:

        while 1:

            l = secF.readline()

            if not l:
                break

            l = secF.readline()
            secL.append(l.split()[0])

    return secL


def gen_simul_data(seqL,secL,motif_detector,signal,noise,tt_bound):

    if seqL == []: seqL = [''] * len(secL)
    if secL == []: secL = [''] * len(seqL)

    print('spiking in signal and noise')

    countA = np.array([motif_detector(p, s) != -1 for p,s in zip(seqL,secL)]).astype(np.float32)[:, np.newaxis]
    targetA = np.random.randn(*countA.shape) * noise + countA * signal

    print('freq: %.5f (%d/%d, %d/%d)' % (sum(countA)/len(countA),sum(countA[:tt_bound]),len(countA[:tt_bound]),sum(countA[tt_bound:]),len(countA[tt_bound:])))

    print('evaluating correlation')

    model = linear_model.Ridge(alpha=0)
    _ = model.fit(countA[:tt_bound], targetA[:tt_bound])

    print('R2:', model.score(countA[:tt_bound], targetA[:tt_bound]), model.score(countA[tt_bound:], targetA[tt_bound:]))

    return countA, targetA


####################################################################################################################
# Input encoding
####################################################################################################################

# def format_rnafold(fold):
#
#    matrixL = []
#
#    for f in fold:
#
#        if f == '(':
#        elif f == ')':
#            arr = [0., 1., 0.]
#        elif f == '.':
#            arr = [0., 0., 1.]
#        elif f == '*':
#            arr = [0.33, 0.33, 0.33]
#        else:
#            raise Exception
#
#        matrixL.append(np.array(arr, dtype=np.float32))
#
#    return np.vstack(matrixL).reshape(1,len(fold),1,3)


def encode(seq='', fold='', dim=2, ambiguous='zero'):

    if seq and not fold:
        channels = 4
    elif not seq and fold:
        channels = 3
    elif seq and fold:
        channels = 7
    else:
        raise Exception

    matrixL = []

    seq_len = max(len(seq), len(fold))

    for i in range(seq_len):

        arr = []

        if seq:

            c = seq[i]

            if c == 'A':
                arr += [1., 0., 0., 0.]
            elif c == 'C':
                arr += [0., 1., 0., 0.]
            elif c == 'G':
                arr += [0., 0., 1., 0.]
            elif c == 'T' or c == 'U':
                arr += [0., 0., 0., 1.]
            elif c == 'N':
                if ambiguous == 'zero':
                    arr += [0, 0, 0, 0]
                elif ambiguous == 'avg':
                    arr += [.25, .25, .25, .25]
                else:
                    raise Exception
            else:
                raise Exception

        if fold:

            f = fold[i]

            if f == '(':
                arr += [1., 0., 0.]
            elif f == ')':
                arr += [0., 1., 0.]
            elif f == '.':
                arr += [0., 0., 1.]
                # arr += [-1., -1., 1.]
                # arr += [-1000., -1000., 1.]
            elif f == '?':
                if ambiguous == 'zero':
                    arr += [0, 0, 0]
                elif ambiguous == 'avg':
                    arr += [.333, .333, .333]
                else:
                    raise Exception
            else:
                raise Exception

        matrixL.append(np.array(arr, dtype=np.float32))

    if dim==2:
        return np.vstack(matrixL).reshape(1, seq_len, 1, channels)
    elif dim==1:
        return np.vstack(matrixL).reshape(1, seq_len, channels)
    else:
        raise Exception


def unify_len(seqL, seq_type, max_seq_len, align='left'):

    resultL = []

    for seq in seqL:

        if len(seq) < max_seq_len:

            deficit = max_seq_len - len(seq)

            if align=='left':
                left_n = 0
            elif align=='center':
                left_n = deficit / 2  # random.randint(0,deficit)
            else:
                raise Exception

            if seq_type == 'nt':
                seq = 'N'*left_n + seq + 'N'*(deficit-left_n)
            elif seq_type == 'fold':
                seq = '?'*left_n + seq + '?'*(deficit-left_n)
            else:
                raise Exception

        elif len(seq) > max_seq_len:

            #seq = seq[:max_seq_len]
            print(seq, max_seq_len)
            raise Exception

        resultL.append(seq)

    return resultL


def pad_seq(seqL, n):
    return ['N'*n + s + 'N'*n for s in seqL]

def pad_fold(foldL, n):
    return ['?'*n + f + '?'*n for f in foldL]

def encode_batch(seqL, foldL, max_seq_len, padding=0, dim=2, align='left', ambiguous='zero'):

    if len(seqL)>0:
        seqL = pad_seq(unify_len(seqL,'nt',max_seq_len, align),padding)
        num_seq = len(seqL)

    if len(foldL)>0:
        foldL = pad_fold(unify_len(foldL,'fold',max_seq_len, align),padding)
        num_seq = len(foldL)

    matrixL = []

    for i in range(num_seq):

        if len(seqL)>0 and len(foldL)==0:
            matrixL.append(encode(seq=seqL[i],dim=dim, ambiguous=ambiguous))
        elif len(seqL)==0 and len(foldL)>0:
            matrixL.append(encode(fold=foldL[i], dim=dim, ambiguous=ambiguous))
        elif len(seqL)>0 and len(foldL)>0:
            matrixL.append(encode(seq=seqL[i],fold=foldL[i],dim=dim, ambiguous=ambiguous))
        else:
            raise Exception

    return np.vstack(matrixL)


def encode_target_track(tPosL, signalL, seqLenL, maxSeqLen, delay):

    matrixL = []

    for i in range(len(tPosL)):

        tPos = eval(tPosL[i])
        signal = eval(signalL[i])

        arr = [0.] * seqLenL[i] + [-1.] * (maxSeqLen-seqLenL[i])

        for j in range(len(tPos)):
            assert(j<seqLenL[i])
            if tPos[j] + delay < maxSeqLen:
                arr[tPos[j] + delay] = signal[j]

        matrixL.append(np.array(arr, dtype=np.float32))

    return np.vstack(matrixL).reshape(len(tPosL), maxSeqLen, 1)


def encode_mask(seq_lenL, max_seq_len, nFilters):

    mat = np.ones((len(seq_lenL), max_seq_len, nFilters), dtype=np.bool_)

    for i, seq_len in enumerate(seq_lenL):
        mat[i, seq_len:max_seq_len, :] = 0

    return mat


def encode_mask_N(seqL, max_seq_len, nFilters): # mask N's in the sequence as well

    mat = np.ones((len(seqL), max_seq_len, nFilters), dtype=np.bool_)

    for i, seq in enumerate(seqL):
        mat[i, len(seq):max_seq_len, :] = 0
        mat[i, [m.start() for m in re.finditer("N",seq)]] = 0

    return mat


####################################################################################################################
# Visualize filter
####################################################################################################################

ntL = ['A', 'C', 'G', 'T']
pairingL = ['(', ')', '.']


def filter_string(filter_weights,char_map,cutoff=0.05,min_max='max'):

    abs_max = np.amax(np.abs(filter_weights))

    if min_max == 'max':
        idx_selectL = np.argmax(filter_weights,axis=1)
    else:
        idx_selectL = np.argmin(filter_weights,axis=1)

    resultL = []

    for pos in range(filter_weights.shape[0]):

        pos_abs_max = np.amax(abs(filter_weights[pos,:]))

        if pos_abs_max <= abs_max * cutoff :
            resultL.append('*')
        else:
            resultL.append(char_map[idx_selectL[pos]])

    return ''.join(resultL)


def filter_summary(W_conv, W_reg=None, cutoff=0.05, min_max='max'):

    filters = W_conv.swapaxes(0, 1).swapaxes(2, 3).swapaxes(1, 2)[0]

    if W_reg == None:
        W_reg = np.ones((filters.shape[0],1))

    num_channels = filters.shape[2]

    if num_channels in (4,7):
        motifL = []
        for i in range(len(filters)):
            motifL.append(filter_string(filters[i][:,:4],ntL,cutoff,min_max))

    if num_channels in (3,7):
        foldL = []
        for i in range(len(filters)):
            foldL.append(filter_string(filters[i][:, -3:], pairingL,cutoff,min_max))

    filterRangeA = np.array([max(abs(np.min(filters[i])), abs(np.max(filters[i]))) for i in range(len(filters))])
    filterPriorityL = list(filterRangeA * W_reg[:,0])
    filterRegrWeightL = list(W_reg[:,0])

    if num_channels == 3:
        summaryL = zip(foldL, filterPriorityL, filterRegrWeightL, range(len(filters)))
        summaryL.sort(lambda x, y: cmp(y[1], x[1]))
    elif num_channels == 4:
        summaryL = zip(motifL, filterPriorityL, filterRegrWeightL, range(len(filters)))
        summaryL.sort(lambda x, y: cmp(y[1], x[1]))
    elif num_channels == 7:
        summaryL = zip(motifL, foldL, filterPriorityL, filterRegrWeightL, range(len(filters)))
        summaryL.sort(lambda x, y: cmp(y[2], x[2]))
    else:
        raise Exception

    return summaryL


def filter_plot(W_conv,filtIdx,style='bwr',cutoff=0.05,figSize=20):

    filters = W_conv.swapaxes(0,1).swapaxes(2,3).swapaxes(1,2)[0]

    fig = plt.figure()

    filter_weights = filters[filtIdx]
    abs_max = np.amax(abs(filter_weights))

    ax = fig.add_subplot(1,1,1)
    cax = ax.matshow(filter_weights.transpose(),cmap=plt.get_cmap(style),vmin=-abs_max,vmax=abs_max)

    if filters.shape[2] == 3:
        ticklabelL = pairingL
        print(filter_string(filter_weights, pairingL, cutoff))
    elif filters.shape[2] == 4:
        ticklabelL = ntL
        print(filter_string(filter_weights, ntL, cutoff))
    elif filters.shape[2] == 7:
        ticklabelL = ntL + pairingL
        print(filter_string(filter_weights[:,:4], ntL, cutoff))
        print(filter_string(filter_weights[:,4:], pairingL, cutoff))
    else:
        raise Exception

    print('range: %.1e, %.1e' % (np.amin(filter_weights),np.amax(filter_weights)))

    ax.set_yticklabels(['']+ticklabelL)
    cbar = fig.colorbar(cax, orientation='horizontal')
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation='vertical')

    fig.set_figwidth(figSize)


####################################################################################################################
# Batch calculation
####################################################################################################################

def batch_loss(tf_sess, tf_loss, tf_feed_dict, nSamples, blockSize):

    loss_tmp = 0.
    nDivide = int(math.ceil(nSamples / float(blockSize)))

    for jj in range(nDivide):

        myIdx1 = jj * blockSize

        if jj == nDivide - 1:
            myIdx2 = nSamples
        else:
            myIdx2 = (jj + 1) * blockSize

        feed_dict_tmp = {}

        for k, v in tf_feed_dict.items():
            if type(v) == np.ndarray:
                feed_dict_tmp[k] = v[myIdx1:myIdx2]
            else:
                feed_dict_tmp[k] = v

        loss_tmp += tf_sess.run(tf_loss, feed_dict=feed_dict_tmp)

    return loss_tmp / nDivide


def batch_pred(tf_sess, tf_Y_pred, tf_feed_dict, nSamples, blockSize):

    Y_pred_tmp = []
    nDivide = int(math.ceil(nSamples / float(blockSize)))

    for jj in range(nDivide):

        myIdx1 = jj * blockSize

        if jj == nDivide - 1:
            myIdx2 = nSamples
        else:
            myIdx2 = (jj + 1) * blockSize

        feed_dict_tmp = {}

        for k, v in tf_feed_dict.items():
            if type(v) == np.ndarray:
                feed_dict_tmp[k] = v[myIdx1:myIdx2]
            else:
                feed_dict_tmp[k] = v

        Y_pred_tmp.append(tf_sess.run(tf_Y_pred, feed_dict=feed_dict_tmp))

    return np.vstack(Y_pred_tmp)

####################################################################################################################
# Performance visualization
####################################################################################################################

from sklearn.metrics import r2_score

def r2_scatter(y_obs,y_pred,test_fold_idx):

    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    fsize = 15

    if test_fold_idx != 0:
        labelL = ['train', 'test']
    else:
        labelL = ['all']

    fig = figure(figsize=(6*len(labelL),6))

    for i,label in enumerate(labelL):

        ax = fig.add_subplot(1,len(labelL),i+1)
        ax.scatter(y_obs[i],y_pred[i],s=10,alpha=0.5,edgecolor='None')
        # ax.set_xlim((-0.02,1.02))
        # ax.set_ylim((-0.02,1.02))
        ax.set_xlabel('obs',fontsize=fsize-3)
        ax.set_ylabel('pre',fontsize=fsize-3)
        ax.set_title('%s set (%s), R2=%.3f, N=%d' % (label,test_fold_idx,r2_score(y_obs[i],y_pred[i]),len(y_obs[i])),fontsize=fsize)
        ax.tick_params(labelsize=fsize-5)

        # ax.text(0.19,0.47,'$R^2=%0.2f$'% model.score(csr_count[idx_train,:],df_info[bool_idx_train].decay),fontsize=fsize)

####################################################################################################################
# KERAS callbacks
####################################################################################################################


from keras.callbacks import Callback

class R2Callback(Callback):

    def __init__(self, target_cutoff, *test_data):
        self.target_cutoff = target_cutoff
        self.test_data = test_data

    def on_epoch_end(self, epoch, logs={}):

        x, y = self.test_data
        y_pred = self.model.predict(x)
        y_flat = y.reshape(y.shape[0]*y.shape[1],1)
        y_pred_flat = y_pred.reshape(y_pred.shape[0]*y_pred.shape[1], 1)

        y_idx = y_flat > self.target_cutoff

        print('Epoch %d R2: %.3f' % (epoch,r2_score(y_flat[y_idx], y_pred_flat[y_idx])))


class BlankCallback(Callback):

    from keras.callbacks import Callback

    def on_epoch_end(self, epoch, logs={}):
        print


####################################################################################################################
# Dynamic input model
####################################################################################################################

def dynamic_input_predict(_model, _x, _l, _y):

    y_L = []
    y_pred_L = []

    for j in range(_y.shape[0]):
        seq_len = _l[j]
        y_test_pred = _model.predict_on_batch(_x[j:j + 1, :seq_len, :])

        y_L.append(_y[j, :seq_len, 0])
        y_pred_L.append(y_test_pred[0, :, 0])

    return np.hstack(y_L), np.hstack(y_pred_L)



