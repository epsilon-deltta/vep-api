#!/usr/bin/python

import sys, os, datetime, re
from glob import glob

# MISC

parseParam = lambda l: {k:v for (k,v) in l}

## SYSTEM CONFIGURATION

Z_base = '/home/jinkuk/'

#storageBase = '/dg-shared/for_jinkuk/pipeline/'
#apacheBase = '/dg-shared/for_jinkuk/pipeline/html/'

storageBase = '/dg-local/for_jinkuk/pipeline/'
apacheBase = '/dg-local/for_jinkuk/pipeline/html/'
sraBase = '/dg-local/for_jinkuk/ncbi/sra/'
resultBase = '/mnt/gtex/results/'
#resultBase = '/dg-shared/gtex_align_2/'
originalBase = '/mnt/gtex/dbGaP-11537/sra/'

init_file_ext = 'gsnap'

def prepare_baseDir(projectN, transfer=True):

    global storageBase, apacheBase

    if glob(storageBase+projectN):
        print ('File directory: already exists')
    else:
        os.system('mkdir %s/%s; chmod a+w %s/%s' % (storageBase,projectN, storageBase,projectN))
        print ('File directory: created')

    if glob(apacheBase+projectN):
        print ('Log directory: already exists')
    else:
        os.system('mkdir %s/%s; chmod a+w %s/%s' % (apacheBase,projectN, apacheBase,projectN))
        print ('Log directory: created')

#    if transfer:
#        os.system('mkdir -p %s/%s' % (resultBase,projectN))

    return(storageBase+projectN)

def logOK(baseDir,logFileN,logExistsFn):

    if glob('%s/%s' % (baseDir,logFileN)):
        return logExistsFn(open('%s/%s' % (baseDir,logFileN)).readlines())
    else:
        return False

def resultOK(baseDir,outFilePostFix):

    resultFileOK = True

    for postFix in outFilePostFix:

        outFileNL = glob('%s/*%s' % (baseDir,postFix))

        if len(outFileNL)<1:
            resultFileOK = False
            break

        for outFileN in outFileNL:
            if os.path.getsize(outFileN)==0:
                resultFileOK = False
                break

    return resultFileOK

def fn_mkdir(htmlF,baseDir):

    htmlF.write('<p><b>Directory: </b>')

    htmlF.write('<br>' + baseDir.replace('//','/'))

    if glob(baseDir):
        htmlF.write(' already exists')
    else:
        ret = os.system('mkdir '+baseDir )

        if ret != 0:
            htmlF.write(' failed to be created')
            sys.exit(1)

        htmlF.write(' created')

    os.system('chmod a+rw '+baseDir)

def fn_ln(htmlF,baseDir,inputFilePathL,sampN):

    htmlF.write('<p><b>Input link files: </b>')

    for inputFileP in inputFilePathL:

        inputFileN = inputFileP.split('/')[-1]
        htmlF.write('<br>' + inputFileN)

        if glob('%s/%s' % (baseDir,inputFileN)):
            htmlF.write(' already exists')
        else:
            ret = os.system('ln -f -s %s %s/%s' % (inputFileP,baseDir,inputFileN))

            if ret!=0:
                htmlF.write(' failed to be created')
                sys.exit(1)

            htmlF.write(' created')

    htmlF.write('</p>')

def fn_execute(htmlF, fn, paramL, paramH, baseDir, logFileN, logExistsFn, execute_onward=False, rerun=False):

    htmlF.write('<p><b>Execution: </b>')

    if execute_onward:
        htmlF.write('Running (previous steps were updated)</p>')
    elif rerun:
        htmlF.write('Running (user enforced execution)</p>')
    else:
        
        logFileOK = logOK(baseDir, logFileN, logExistsFn)

        if not logFileOK:
            htmlF.write('Running (not-integral log file)</p>')
        else:
            htmlF.write('Skipping (previously completed)</p>')

    if execute_onward or rerun or not logFileOK:
        apply(fn,paramL,paramH)
        execute_onward = True

    return execute_onward

def fn_content(htmlF,baseDir,logFileN):

    htmlF.write('<p><b>Log:</b>')
    htmlF.write('<div class="log_box" style="height:400px;width:65%;border:1px solid #ccc;overflow:auto;">')
    htmlF.write('<pre>' + ''.join(open('%s/%s' % (baseDir,logFileN)).readlines()) + '</pre></p>')
    htmlF.write('</div>')

def fn_results(htmlF, baseDir, outFilePostFix):
    
    htmlF.write('<p><b>Result files:</b><br>')

    for postFix in outFilePostFix:
        outFileNL = glob('%s/*%s' % (baseDir, postFix))
#        if len(outFileNL) == -1 or os.path.getsize(outFileNL[0]) != 0:
        if len(outFileNL) > 0:
            for outFileN in outFileNL:
                sizeF = (float(os.path.getsize(outFileN)))/(1024*1024)
                creationD = datetime.datetime.fromtimestamp(os.path.getmtime(outFileN)).replace(microsecond=0)
                htmlF.write('-- %s , %s (%.3f MB) <br>' % (creationD, outFileN.split('/')[-1], sizeF))

def fn_links(htmlF, projectN, baseDir, outLinkPostFix):

    for postFix in outLinkPostFix:
        outLinkNL = glob('%s/*%s' % (baseDir, postFix))
        for outLinkN in outLinkNL:
            if os.path.getsize(outLinkN) != 0:
                creationD = datetime.datetime.fromtimestamp(os.path.getmtime(outLinkN)).replace(microsecond=0)
                htmlF.write('-- %s, <a href="./%s">%s</a> <br>' % (creationD, outLinkN[len(storageBase + projectN)+1:], outLinkN[len(baseDir)+1:]))
#                htmlF.write('-- %s, <a href="./%s">%s</a> <br>' % (creationD, outLinkN[len(apacheBase+projectN)+1:], outLinkN[len(apacheBase+projectN)+1:] ))

def fn_files(htmlF,baseDir,prevFileS):

    htmlF.write('<p><b>New files:</b><br>')

    allFileS = set(glob(baseDir+'/*'))

    newFileL = list(allFileS.difference(prevFileS))
    newFileL.sort(lambda x,y: cmp(x,y))

    for i in range(len(newFileL)):
        sizeF = (float(os.path.getsize(newFileL[i])))/(1024*1024)
        htmlF.write('-- %s. %s (%.3f MB) <br>' % (i+1, newFileL[i].split('/')[-1], sizeF))

    htmlF.write('</p>')

    return allFileS

def fn_clean(baseDir, logFileN, logExistsFn, outFilePostFix, deleteFilePostFix, htmlF):
    
    if not logOK(baseDir,logFileN,logExistsFn): #failed execution
        print 'stopping pipeline: log file not ok'
        htmlF.write('<br><b>Stopping pipeline due to not-integral log file</b></br>')
        sys.exit(1)
        return

    if not resultOK(baseDir,outFilePostFix): # failed execution or previously successfully run
        print 'skipping cleaning: result not ok'
        return

    for fN in glob(baseDir+'/*'):
        
        if sum(map(fN.endswith, deleteFilePostFix)) > 0:
            print 'deleting %s' % (fN,)
            os.system('rm %s' % (fN,))

def main(inputFilePathL, genSpecFn, sampN, projectN='test_yn', server='smc1', genome='hg19', sra_dump=False, clean=False, transfer=False, delete_original=False):

    ## adjust storageBase/apacheBase depending on project(CNA, CNA_corr, Purity, Clonality)
    prepare_baseDir(projectN, transfer)

    # HTML log file initiation
    htmlFileN = '%s/%s/%s.html' % (apacheBase,projectN,sampN)
    htmlF = open(htmlFileN, 'w', 0)
    os.system('chmod a+rw %s' % htmlFileN)

    htmlF.write('<DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd"><html><head></head><body>')

    # creating sample data directory and linking input file

    startTime = datetime.datetime.now().replace(microsecond=0)

    baseDir = '%s/%s/%s' % (storageBase,projectN,sampN)

    htmlF.write('<h2>%s, %s</h2>' % (projectN,sampN))
    htmlF.write('<hr><b>Step 0: Set-up</b><hr>')

    fn_mkdir(htmlF,baseDir)

    inputFilePre = inputFilePathL[0].split('/')[-1].split('.')[0]

    if sra_dump:

        if len(inputFilePathL) != 1:
            raise Exception

        htmlF.write('<p/><b>SRA to FASTQ dump: </b><br>')


        htmlF.write('Creating link in SRA baseDir: ln -f -s %s %s<br>' % (inputFilePathL[0],sraBase,))
        os.system('ln -f -s %s %s' % (inputFilePathL[0],sraBase,))

        htmlF.write('Extracting FASTQ from SRA: cd %s;/home/jinkuk/D/tools/sratoolkit.2.7.0-ubuntu64/bin/fastq-dump %s.sra --split-3<br>' % (sraBase,inputFilePre))
        os.system('cd %s;/home/jinkuk/D/tools/sratoolkit.2.7.0-ubuntu64/bin/fastq-dump %s.sra --split-3' % (sraBase,inputFilePre))


        htmlF.write('Renaming the extracted fastq files')
        os.system('cd %s;mv -f %s_1.fastq %s/%s.1.%s.fq' % (sraBase,inputFilePre,baseDir,inputFilePre,init_file_ext))
        os.system('cd %s;mv -f %s_2.fastq %s/%s.2.%s.fq' % (sraBase,inputFilePre,baseDir,inputFilePre,init_file_ext))

        if glob('%s/%s.fastq' % (sraBase,inputFilePre)):
            os.system('cd %s;mv -f %s.fastq %s/%s.U.%s.fq' % (sraBase,inputFilePre,baseDir,inputFilePre,init_file_ext))
        else:
            os.system('touch %s/%s.U.%s.fq' % (baseDir,inputFilePre,init_file_ext))

        inputFilePathL = ('%s/%s.1.%s.fq' % (sraBase,inputFilePre,init_file_ext),'%s/%s.2.%s.fq' % (sraBase,inputFilePre,init_file_ext),'%s/%s.U.%s.fq' % (sraBase,inputFilePre,init_file_ext))

    else:

        fn_ln(htmlF,baseDir,inputFilePathL,sampN)

    prevFileS = fn_files(htmlF,baseDir,set([]))

    endTime = datetime.datetime.now().replace(microsecond=0)
    elapsedT = (endTime - startTime)
    htmlF.write('<b> Step 0 elapsed time : %s </b><br><br>' % (elapsedT,))

    # Step 1-N

    execute_onward = False
    specL = genSpecFn(baseDir, server, genome)

    for i in range(len(specL)):

        startTime = datetime.datetime.now().replace(microsecond=0)

        logFileN = '%s%s' % (sampN,specL[i]['logPostFix'])

        htmlF.write('<hr><b>Step %s: %s: %s</b><hr>' % (i+1,specL[i]['name'],specL[i]['desc']))
        
        execute_onward = fn_execute(htmlF, specL[i]['fun'], specL[i]['paramL'], specL[i]['paramH'], baseDir, logFileN, specL[i]['logExistsFn'], execute_onward, specL[i]['rerun']) # execute if conditions are met

        fn_content(htmlF,baseDir,logFileN)

        fn_results(htmlF, baseDir, specL[i]['outFilePostFix'])

        if 'outLinkPostFix' in specL[i]:
            fn_links(htmlF, projectN, baseDir, specL[i]['outLinkPostFix'])

        prevFileS = fn_files(htmlF,baseDir,prevFileS)
        
        endTime = datetime.datetime.now().replace(microsecond=0)
        elapsedT = (endTime - startTime)
        htmlF.write('<b> Step %s elapsed time : %s </b><br><br>' % (i+1, elapsedT))

        if specL[i]['clean']: # clean if user-specified
             fn_clean(baseDir, logFileN, specL[i]['logExistsFn'], specL[i]['outFilePostFix'], specL[i]['deleteFilePostFix'], htmlF) # clean only if both log and results are integral

#    htmlF.write('<hr><b>Step %s: Final data file transfer to S3</b><hr>' % (i+1,))
#    htmlF.write('<p>')
#
#    for t in transfer:
#
#        os.system('aws s3 mb s3://jinkuk-%s' % (projectN))
#
#        for p in glob('%s/*%s' % (baseDir,t)):
#
#            os.system('aws s3 cp %s s3://jinkuk-%s' % (p,projectN))
#            htmlF.write('aws s3 cp %s s3://jinkuk-%s' % (p,projectN))
#
#    htmlF.write('</p>')

    if clean:

        # deleting sra link/file
        os.system('rm -f %s/%s.sra' % (sraBase,inputFilePre))
        os.system('rm -f %s/%s.sra.vdbcache.cache' % (sraBase,inputFilePre))

        # deleting fq file
        os.system('rm -f %s/%s.1.%s.fq' % (sraBase,inputFilePre,init_file_ext))
        os.system('rm -f %s/%s.2.%s.fq' % (sraBase,inputFilePre,init_file_ext))

        htmlF.write('<p>sra cleanup completed</p>')
        print 'sra cleanup completed'

    if transfer:

        os.system('mkdir -p %s/%s' % (resultBase,projectN))
        os.system('mv -f %s %s/%s' % (baseDir,resultBase,projectN))
        htmlF.write('<p>result transfer completed</p>')
        print 'result transfer completed'

    if delete_original:

        os.system('rm -f %s/%s.sra' % (originalBase,inputFilePre))
        os.system('rm -f %s/%s.sra.vdbcache.cache' % (originalBase,inputFilePre))

    htmlF.close()
