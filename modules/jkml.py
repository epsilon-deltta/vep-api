import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#import tensorflow as tf

nts = ['A','C','G','T']
pairing = ['(',')','.']

chanel7 = nts+pairing

## this should go to a module: jkml

def format_seq_2d(seq):
    
    matrixL = []
    
    for c in seq:
        
        if c == 'A':
            arr = np.array([1., 0., 0., 0.], dtype=np.float32)
        elif c == 'C':
            arr = np.array([0., 1., 0., 0.], dtype=np.float32)
        elif c == 'G':
            arr = np.array([0., 0., 1., 0.], dtype=np.float32)
        elif c == 'T' or c == 'U':
            arr = np.array([0., 0., 0., 1.], dtype=np.float32)
        elif c == 'N':
            arr = np.array([.25, .25, .25, .25], dtype=np.float32)
        else:
            raise Exception

        matrixL.append(arr)
    
    return np.vstack(matrixL).reshape(len(seq),4)

def format_seq(seq):
    
    matrixL = []
    
    for c in seq:
        
        if c == 'A':
            arr = np.array([1., 0., 0., 0.], dtype=np.float32)
        elif c == 'C':
            arr = np.array([0., 1., 0., 0.], dtype=np.float32)
        elif c == 'G':
            arr = np.array([0., 0., 1., 0.], dtype=np.float32)
        elif c == 'T' or c == 'U':
            arr = np.array([0., 0., 0., 1.], dtype=np.float32)
        elif c == 'N':
            arr = np.array([.25, .25, .25, .25], dtype=np.float32)
        else:
            raise Exception

        matrixL.append(arr)
    
    return np.vstack(matrixL).reshape(1,len(seq),1,4)

def format_seq_batch(seqL,max_seq_len):

    matrixL = []
    
    for seq in seqL:

        if len(seq) < max_seq_len:
            seq += 'N'*(max_seq_len-len(seq))
        
        if len(seq) > max_seq_len:
            seq = seq[:max_seq_len]

        matrixL.append(format_seq(seq))
        
    return np.vstack(matrixL)

def format_seq_rnafold(seq,fold):
    
    matrixL = []
    
    for c,f in zip(seq,fold):
        
        if c == 'A':
            arr = [1., 0., 0., 0.]
        elif c == 'C':
            arr = [0., 1., 0., 0.]
        elif c == 'G':
            arr = [0., 0., 1., 0.]
        elif c == 'T' or c == 'U':
            arr = [0., 0., 0., 1.]
        elif c == 'N':
            arr = [.25, .25, .25, .25]
        else:
            raise Exception

        if f == '(':
            arr += [1., -1., -1.]
        elif f == ')':
            arr += [-1., 1., -1.]
        elif f == '.':
            #arr += [-1., -1., 1.]
            arr += [-1000., -1000., 1.]
        elif f == '*':
            arr += [0., 0., 0.]
        else:
            raise Exception

        matrixL.append(np.array(arr, dtype=np.float32))

    return np.vstack(matrixL).reshape(1,len(seq),1,7)

def format_seq_rnafold_batch(seqL,foldL,max_seq_len):

    matrixL = []
    
    for seq,fold in zip(seqL,foldL):

        if len(seq) < max_seq_len:
            seq += 'N'*(max_seq_len-len(seq))
            fold += '*'*(max_seq_len-len(fold))

        elif len(seq) > max_seq_len:
            seq = seq[:max_seq_len]
            fold = fold[:max_seq_len]

        matrixL.append(format_seq_rnafold(seq,fold))
        
    return np.vstack(matrixL)

#def format_rnafold(fold):
#    
#    matrixL = []
#    
#    for f in fold:
#        
#        if f == '(':
#            arr = [1., 0., 0.]
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

def format_rnafold(fold):
    
    matrixL = []
    
    for f in fold:
        
        if f == '(':
            arr = [1., -1., -1.]
        elif f == ')':
            arr = [-1., 1., -1.]
        elif f == '.':
            arr = [-100., -100., 1.]
        elif f == '*':
            arr = [0, 0, 0]
        else:
            raise Exception

        matrixL.append(np.array(arr, dtype=np.float32))

    return np.vstack(matrixL).reshape(1,len(fold),1,3)

def format_rnafold_batch(foldL,max_seq_len):

    matrixL = []
    
    for fold in foldL:

        if len(fold) < max_seq_len:
            fold += '*'*(max_seq_len-len(fold))

        elif len(fold) > max_seq_len:
            fold = fold[:max_seq_len]

        matrixL.append(format_rnafold(fold))
        
    return np.vstack(matrixL)

def show_filter(W_conv,W_reg,showImage=False,figWidth=20):

    filters = W_conv.swapaxes(0,1).swapaxes(2,3).swapaxes(1,2)[0]

    motifL = []

    for i in range(len(filters)):
        motifL.append(''.join(map(lambda x: nts[x], np.argmax(filters[i],axis=1))))

    motif_regL = zip(motifL,W_reg[:,0],range(len(motifL)))
    motif_regL.sort(lambda x,y: cmp(y[1],x[1]))
    
    if showImage:

        fig = plt.figure()

        numPlots = min(len(filters),5)
        
        for i in range(numPlots):

            filter_values = filters[motif_regL[i][2]].transpose()

            abs_max = max(abs(np.min(filter_values)),abs(np.max(filter_values))) * 0.5

            ax = fig.add_subplot(1,numPlots+1,i+1)
            cax = ax.matshow(filter_values,cmap=plt.get_cmap('gray_r'),vmin=-abs_max,vmax=abs_max)

            ax.set_yticklabels(['']+nts)
            cbar = fig.colorbar(cax, orientation='horizontal')
            cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation='vertical')

        fig.set_figwidth(figWidth)

    return motif_regL

def show_filter_rnafold(W_conv,W_reg,showImage=False,figWidth=20):

    filters = W_conv.swapaxes(0,1).swapaxes(2,3).swapaxes(1,2)[0]

    motifL = []

    for i in range(len(filters)):
        motifL.append(''.join(map(lambda x: nts[x], np.argmax(filters[i][:,:4],axis=1))))

    foldL = []

    for i in range(len(filters)):
        foldL.append(''.join(map(lambda x: pairing[x], np.argmax(filters[i][:,4:],axis=1))))

    motif_regL = zip(motifL,foldL,W_reg[:,0],range(len(motifL)))
    motif_regL.sort(lambda x,y: cmp(y[2],x[2]))
    
    if showImage:

        fig = plt.figure()

        numPlots = min(len(filters),5)
        
        for i in range(numPlots):

            filter_values = filters[motif_regL[i][2]].transpose()

            abs_max = max(abs(np.min(filter_values)),abs(np.max(filter_values))) * 0.5

            ax = fig.add_subplot(1,numPlots+1,i+1)
            cax = ax.matshow(filter_values,cmap=plt.get_cmap('gray_r'),vmin=-abs_max,vmax=abs_max)

            ax.set_yticklabels(['']+nts)
            cbar = fig.colorbar(cax, orientation='horizontal')
            cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation='vertical')

        fig.set_figwidth(figWidth)

    return motif_regL

def show_filter_rnafold_only(W_conv,W_reg,figWidth=20):

    filters = W_conv.swapaxes(0,1).swapaxes(2,3).swapaxes(1,2)[0]

    foldL = []

    for i in range(len(filters)):
        foldL.append(''.join(map(lambda x: pairing[x], np.argmax(filters[i],axis=1))))

    motif_regL = zip(foldL,W_reg[:,0],range(len(foldL)))
    motif_regL.sort(lambda x,y: cmp(y[1],x[1]))

    return motif_regL

def show_filter_rnafold2(W_conv,W_reg,showImage=False,figWidth=20):

    filters = W_conv.swapaxes(0,1).swapaxes(2,3).swapaxes(1,2)[0]

    motifL = []

    for i in range(len(filters)):
        motifL.append(''.join(map(lambda x: chanel7[x], np.argmax(filters[i],axis=1))))

    motif_regL = zip(motifL,W_reg[:,0],range(len(motifL)))
    motif_regL.sort(lambda x,y: cmp(y[1],x[1]))
    
    if showImage:

        fig = plt.figure()

        numPlots = min(len(filters),5)
        
        for i in range(numPlots):

            filter_values = filters[motif_regL[i][2]].transpose()

            abs_max = max(abs(np.min(filter_values)),abs(np.max(filter_values))) * 0.5

            ax = fig.add_subplot(1,numPlots+1,i+1)
            cax = ax.matshow(filter_values,cmap=plt.get_cmap('gray_r'),vmin=-abs_max,vmax=abs_max)

            ax.set_yticklabels(['']+nts)
            cbar = fig.colorbar(cax, orientation='horizontal')
            cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation='vertical')

        fig.set_figwidth(figWidth)

    return motif_regL

def tabulize_filter(W_conv,W_reg,outFilePath=None):

    filters = W_conv.swapaxes(0,1).swapaxes(2,3).swapaxes(1,2)[0]

    motifL = []

    for i in range(len(filters)):
        motifL.append(''.join(map(lambda x: nts[x], np.argmax(filters[i],axis=1))))

    df = pd.DataFrame({'motif':motifL,'weight':W_reg[:,0]})
    df.sort_values(by='weight',ascending=False,inplace=True)
    df.reset_index(inplace=True)

    if outFilePath: df.to_csv(outFilePath)

    return df

def show_one_filter(W_conv,filtIdx):

    filters = W_conv.swapaxes(0,1).swapaxes(2,3).swapaxes(1,2)[0]

    fig = plt.figure()

    filter_values = filters[filtIdx].transpose()

    abs_max = max(abs(np.min(filter_values)),abs(np.max(filter_values))) * 0.5

    ax = fig.add_subplot(1,1,1)
    cax = ax.matshow(filter_values,cmap=plt.get_cmap('Greys'),vmin=-abs_max,vmax=abs_max)

    ax.set_yticklabels(['']+nts)
    cbar = fig.colorbar(cax, orientation='horizontal')
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation='vertical')

    fig.set_figwidth(8)

def show_one_filter_rnafold(W_conv,filtIdx):

    filters = W_conv.swapaxes(0,1).swapaxes(2,3).swapaxes(1,2)[0]

    fig = plt.figure()

    filter_values = filters[filtIdx].transpose()

    abs_max = max(abs(np.min(filter_values)),abs(np.max(filter_values))) * 0.5

    ax = fig.add_subplot(1,1,1)
    cax = ax.matshow(filter_values,cmap=plt.get_cmap('bwr'),vmin=-abs_max,vmax=abs_max)

    ax.set_yticklabels(['']+nts+pairing)
    cbar = fig.colorbar(cax, orientation='horizontal')
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation='vertical')

    fig.set_figwidth(20)

def show_one_filter_rnafold_only(W_conv,filtIdx):

    filters = W_conv.swapaxes(0,1).swapaxes(2,3).swapaxes(1,2)[0]

    fig = plt.figure()

    filter_values = filters[filtIdx].transpose()

    abs_max = max(abs(np.min(filter_values)),abs(np.max(filter_values))) * 0.5

    ax = fig.add_subplot(1,1,1)
    cax = ax.matshow(filter_values,cmap=plt.get_cmap('gray_r'),vmin=-abs_max,vmax=abs_max)

    ax.set_yticklabels(['']+pairing)
    cbar = fig.colorbar(cax, orientation='horizontal')
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation='vertical')

    fig.set_figwidth(10)

def softmax(x):
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()

def get_filter_text(W_conv,filtIdx,outFilePath):

    filters = W_conv.swapaxes(0,1).swapaxes(2,3).swapaxes(1,2)[0]

    outFile = open(outFilePath,'w')
    outFile.write('\t%s\n' % '\t'.join(nts))

    for i,posValues in enumerate(filters[filtIdx]):
        outFile.write("%s\t%s\n" % (i+1,'\t'.join(map(str,posValues))))

def get_filter_idx_by_rank(W_conv,W_reg,filt_rank):

    filters = W_conv.swapaxes(0,1).swapaxes(2,3).swapaxes(1,2)[0]

    motifL = []

    for i in range(len(filters)):
        motifL.append(''.join(map(lambda x: nts[x], np.argmax(filters[i],axis=1))))

    motif_regL = zip(motifL,W_reg[:,0],range(len(motifL)))
    motif_regL.sort(lambda x,y: cmp(y[1],x[1]))

    return motif_regL[filtIdx][2]

def get_filter_mat(W_conv,W_reg,filtIdx=0):

    filters = W_conv.swapaxes(0,1).swapaxes(2,3).swapaxes(1,2)[0]

    motifL = []

    for i in range(len(filters)):
        motifL.append(''.join(map(lambda x: nts[x], np.argmax(filters[i],axis=1))))

    motif_regL = zip(motifL,W_reg[:,0],range(len(motifL)))
    motif_regL.sort(lambda x,y: cmp(y[1],x[1]))

    return filters[motif_regL[filtIdx][2]]


def batch_loss(tf_sess,tf_loss,tf_feed_dict,nSamples,nDivide):

    loss_tmp = 0.
    nBatch = nSamples/nDivide

    for jj in range(nDivide):
        
        myIdx1 = jj*nBatch
        myIdx2 = (jj+1)*nBatch

        feed_dict_tmp = {}
        for k,v in tf_feed_dict.items():
            if type(v) == np.ndarray:
                feed_dict_tmp[k] = v[myIdx1:myIdx2]
            else:
                feed_dict_tmp[k] = v
            
        loss_tmp += tf_sess.run(tf_loss, feed_dict=feed_dict_tmp)

    return loss_tmp/float(nDivide)
    
def batch_pred(tf_sess,tf_Y_pred,tf_feed_dict,nSamples,nDivide):

    Y_pred_tmp = []
    nBatch = nSamples/nDivide

    for jj in range(nDivide):
        
        myIdx1 = jj*nBatch
        myIdx2 = (jj+1)*nBatch

        feed_dict_tmp = {}
        for k,v in tf_feed_dict.items():
            if type(v) == np.ndarray:
                feed_dict_tmp[k] = v[myIdx1:myIdx2]
            else:
                feed_dict_tmp[k] = v
            
        Y_pred_tmp.append(tf_sess.run(tf_Y_pred, feed_dict=feed_dict_tmp))

    return np.vstack(Y_pred_tmp)

def match_seq_fold(seq_regex_compiled,fold_pat,seq,fold):

    if len(seq) != len(fold):
        print len(seq), len(fold)
        raise Exception

    target_len = len(seq)
    pat_len = len(fold_pat)

    for i in range(target_len):

        if fold[i:i+pat_len]==fold_pat and seq_regex_compiled.match(seq[i:i+pat_len]):
            return True

    return False
