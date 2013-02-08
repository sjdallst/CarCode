#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      2013sdallstr
#
# Created:     05/10/2012
# Copyright:   (c) 2013sdallstr 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
import wave
import struct
import numpy
import scipy
import matplotlib
from pylab import *
from math import *
from itertools import *
import os
from pprint import pprint

#slice spectrogrma to get psd vector for a spectrum of frequencies at each period of time

#take two signals slice them by equal time intervals cross correlate the fft of each chunk
def main():
    resultlist= [0,0,0,0]
    refrenceffts = generaterefFFT()
    for item in os.listdir("mldata/back"):
        result = matchwords("mldata/back/"+item,refrenceffts)[0]
        resultlist[result] +=1



    print("--BACK RESULTS--")
    print("")
    print("Number Correct: " + str(resultlist[0]) + " Accuracy: " + str(resultlist[0]/59.0*100) )
    print("Number Wrong: " + str(sum(resultlist[1:])) + "  Percent Wrong: " + str(sum(resultlist[1:])/59.0*100))
    print("Predicted as LEFT: " + str(resultlist[2])  + " Percent: " + str(100*float(resultlist[2])/(sum(resultlist[1:]))))
    print("Predicted as RIGHT: " + str(resultlist[3])  + " Percent: " + str(100*float(resultlist[3])/(sum(resultlist[1:]))))
    print("Predicted as FORWARD: " + str(resultlist[1])  + " Percent: " + str(100*float(resultlist[1])/(sum(resultlist[1:]))))
    print("")
    resultlist= [0,0,0,0]

    for item in os.listdir("mldata/forward"):
        result = matchwords("mldata/forward/"+item,refrenceffts)[0]
        resultlist[result] +=1



    print("--FORWARD RESULTS--")
    print("")
    print("Number Correct: " + str(resultlist[1]) + " Accuracy: " + str(resultlist[1]/59.0*100) )
    print("Number Wrong: " + str(59-resultlist[1]) + "  Percent Wrong: " + str((59-resultlist[1])/59.0*100))
    print("Predicted as LEFT: " + str(resultlist[2])  + " Percent: " + str(100*float(resultlist[2])/(59-resultlist[1])))
    print("Predicted as RIGHT: " + str(resultlist[3])  + " Percent: " + str(100*float(resultlist[3])/(59-resultlist[1])))
    print("Predicted as BACK: " + str(resultlist[0])  + " Percent: " + str(100*float(resultlist[0])/(59-resultlist[1])))
    print("")
    resultlist= [0,0,0,0]

    for item in os.listdir("mldata/left"):

        result = matchwords("mldata/left/"+item,refrenceffts)

        resultlist[result[0]] +=1



    print("--LEFT RESULTS--")
    print("")
    print("Number Correct: " + str(resultlist[2]) + " Accuracy: " + str(resultlist[2]/59.0*100) )
    print("Number Wrong: " + str(59-resultlist[2]) + "  Percent Wrong: " + str((59-resultlist[2])/59.0*100))
    print("Predicted as FORWARD: " + str(resultlist[1])  + " Percent: " + str(100*float(resultlist[1])/(59-resultlist[2])))
    print("Predicted as RIGHT: " + str(resultlist[3])  + " Percent: " + str(100*float(resultlist[3])/(59-resultlist[2])))
    print("Predicted as BACK: " + str(resultlist[0])  + " Percent: " + str(100*float(resultlist[0])/(59-resultlist[2])))
    print("")

    resultlist= [0,0,0,0]

    for item in os.listdir("mldata/right"):
        result = matchwords("mldata/right/"+item,refrenceffts)[0]
        resultlist[result] +=1



    print("--RIGHT RESULTS--")
    print("")
    print("Number Correct: " + str(resultlist[3]) + " Accuracy: " + str(resultlist[3]/59.0*100) )
    print("Number Wrong: " + str(59-resultlist[3]) + "  Percent Wrong: " + str((59-resultlist[3])/59.0*100))
    print("Predicted as FORWARD: " + str(resultlist[1])  + " Percent: " + str(100*float(resultlist[1])/(59-resultlist[3])))
    print("Predicted as LEFT: " + str(resultlist[2])  + " Percent: " + str(100*float(resultlist[2])/(59-resultlist[3])))
    print("Predicted as BACK: " + str(resultlist[0])  + " Percent: " + str(100*float(resultlist[0])/(59-resultlist[3])))
    print("")

def matchwords(testwavname,refdata):
    testwave = getwave(testwavname)

    backfref_fft,forwardref_fft,leftref_fft,rightref_fft,refbackwave,refforwardwave,refleftwave,refrightwave = refdata

    test_backfft = normalizedfft(truncatewaves(refbackwave,testwave))
    test_backfft = test_backfft/norm(test_backfft)
    test_forwardfft = normalizedfft(truncatewaves(refforwardwave,testwave))
    test_forwardfft =test_forwardfft/norm(test_forwardfft)
    test_leftfft = normalizedfft(truncatewaves(refleftwave,testwave))
    test_leftfft = test_leftfft/norm(test_leftfft)
    test_rightfft = normalizedfft(truncatewaves(refrightwave,testwave))
    test_rightfft = test_rightfft/norm(test_rightfft)

    forwardsimilarity = euclidean(forwardref_fft,test_forwardfft)

    backsimilarity = euclidean(backfref_fft,test_backfft)

    rightsimilarity = euclidean(leftref_fft,test_leftfft)

    leftsimilarity = euclidean(rightref_fft,test_rightfft)

    """
    cosinesimilarityforwardvec = getsimilarityvector(specgram(refforwardwave),specgram(truncatewaves(refforwardwave,testwave)))
    cosinesimilaritybackwardvec = getsimilarityvector(specgram(refbackwardwave),specgram(truncatewaves(refbackwardwave,testwave)))
    cosinesimilarityrightvec = getsimilarityvector(specgram(refrightwave),specgram(truncatewaves(refrightwave,testwave)))
    cosinesimilarityleftvec = getsimilarityvector(specgram(refleftwave),specgram(truncatewaves(refleftwave,testwave)))
    """

    """
    forwardsimilarity = numpy.mean(cosinesimilarityforwardvec)
    backsimilarity = numpy.mean(cosinesimilaritybackwardvec)
    rightsimilarity = numpy.mean(cosinesimilarityrightvec)
    leftsimilarity = numpy.mean(cosinesimilarityleftvec)
    """
    """
    refbackwave = getwave("backsam30.wav")
    refforwardwave = getwave("forwardsam30.wav")
    refleftwave = getwave("leftsam30.wav")
    refrightwave = getwave("rightsam30.wav")





    backfref_fft = normalizedfft(refbackwave)
    backfref_fft = backfref_fft/norm(backfref_fft)

    forwardref_fft = normalizedfft(refforwardwave)
    forwardref_fft = forwardref_fft/norm(forwardref_fft)
    leftref_fft = normalizedfft(refleftwave)
    leftref_fft = leftref_fft/norm(leftref_fft)
    rightref_fft = normalizedfft(refrightwave)

    rightref_fft = rightref_fft/norm(rightref_fft )
    """

    similaritydict = {backsimilarity: 0,forwardsimilarity:1,leftsimilarity:2,rightsimilarity:3 }

    closest = min([backsimilarity,forwardsimilarity,rightsimilarity,leftsimilarity ])
    #print("Closest Match: ", similaritydict[closest])
    return (similaritydict[closest],similaritydict)

def generaterefFFT():
    refbackwave = getwave("backsam30.wav")
    refforwardwave = getwave("forwardsam30.wav")
    refleftwave = getwave("leftsam60.wav")
    refrightwave = getwave("rightsam30.wav")





    backfref_fft = normalizedfft(refbackwave)
    backfref_fft = backfref_fft/norm(backfref_fft)

    forwardref_fft = normalizedfft(refforwardwave)
    forwardref_fft = forwardref_fft/norm(forwardref_fft)
    leftref_fft = normalizedfft(refleftwave)
    leftref_fft = leftref_fft/norm(leftref_fft)
    rightref_fft = normalizedfft(refrightwave)

    rightref_fft = rightref_fft/norm(rightref_fft )

    return(backfref_fft,forwardref_fft,leftref_fft,rightref_fft,refbackwave,refforwardwave,refleftwave,refrightwave)



def getwave(filename):
    wave1 = wave.open(filename,"r")
    wave1length =  wave1.getnframes()
    soundlist = [(ord(x)-128)for x in wave1.readframes(wave1length)]
    soundlist= filtersound(soundlist)
    return soundlist

def euclidean(vec1,vec2):
    sumval = 0
    for index in range(len(vec1)):
        sumval += (vec1[index]-vec2[index])**2

    return sqrt(sumval)

def cosinesimilarity(testvector,refvector):
    dot = sum([refvector[i] * testvector[i] for i in range(len(refvector))])
    sumrefvector = sum([refvector[i]**2 for i in range(len(refvector))])
    magrefvector = sqrt(sumrefvector)
    sumtestvector = sum([testvector[i]**2 for i in range(len(refvector))])
    magtestvector = sqrt(sumtestvector)

    return (dot/(magrefvector*magtestvector))


def truncatewaves(ref,test):

    if(len(ref)<len(test)):
        return test[:len(ref)]
    else:
        test.extend([0]*(len(ref)-len(test)))
        return test

def getsimilarityvector(refgram,testgram):
    powerref = refgram[0]
    powertest = testgram[0]
    cosinesimilarityvec = []
    for x in range(powerref.shape[0]):

        powertestslice = numpy.take(powertest,[x],axis = 0)
        powerrefslice = numpy.take(powerref,[x],axis = 0)
        simscore = cosinesimilarity(powerrefslice,powertestslice)
        cosinesimilarityvec.append(simscore)

    return cosinesimilarityvec





"""
def rms(wavearray):
    return sqrt(1.0/len(wavearray)*sum([x*x for x in wavearray]))
def rmsnormalize(refrms,wavearray):
    poslist = [x for x in wavearray if x >0]
    print(sum(poslist),sum(wavearray))
    a = len(wavearray)
    b = 2*((-1*sum(poslist)*len(poslist))+((a- len(poslist))*(sum(wavearray)-sum(poslist))))
    c = sum([x*x for x in wavearray]) - refrms**2*len(wavearray)
    offset = (-b+sqrt(b**2-4*a*c))/2*a
    for index,item in enumerate(wavearray[:]):
        if item == 0:
            continue
        if item < 0:
            wavearray[index] += -1*offset
        else:
            wavearray[index] += offset
"""
def filtersound(wavearray):
    startindex = 0
    for index,val in enumerate(wavearray):
        if val != 0 and val != -1:
            startindex = index
            break
    stopindex = 0
    for index,val in reversed(list(enumerate(wavearray))):
        if val != 128 and val != 127:
            stopindex = index
            break

    return wavearray[startindex:stopindex+1]

def normalizedfft(framelist):
    nUniquePts = ceil((len(framelist)+1)/2.0)


    fftframelist = fft(framelist)
    fftframelist = fftframelist[0:nUniquePts]
    fftframelist = abs(fftframelist)
    fftframelist=fftframelist/float(len(framelist))
    fftframelist = fftframelist**2
    p = fftframelist
    if len(framelist) % 2 > 0: # we've got odd number of points fft
      p[1:len(p)] = p[1:len(p)] * 2
    else:
      p[1:len(p) -1] = p[1:len(p) - 1] * 2 # we've got even number of points fft
    p = [10*log(x) for x in p]
    return p
    """
    nUniquePts = ceil((wave1length+1)/2.0)
    wave2samplefreq = wave2.getframerate()n
    wave2length =  wave2.getnframes()
    fftframelist = fft(framelist)
    fftframelist = fftframelist[0:nUniquePts]
    fftframelist = abs(fftframelist)
    fftframelist=fftframelist/float(wave1length)
    fftframelist = fftframelist**2
    p = fftframelist
    if wave1length % 2 > 0: # we've got odd number of points fft
      p[1:len(p)] = p[1:len(p)] * 2
    else:
      p[1:len(p) -1] = p[1:len(p) - 1] * 2 # we've got even number of points fft
    freqArray = arange(0, nUniquePts, 1.0) * (wave1samplefreq / float(wave1length))
    p = [10*log(x) for x in p]
    plot(freqArray, p)
    xlabel('Frequency (kHz)')
    ylabel('Power (dB)')
    show()
    """
"""
    backfref_fft = normalizedfft(refbackwave)

    forwardref_fft = normalizedfft(refforwardwave)
    leftref_fft = normalizedfft(refleftwave)
    rightref_fft = normalizedfft(refrightwave)

    test_backfft = normalizedfft(truncatewaves(refbackwave,testwave))
    test_forwardfft = normalizedfft(truncatewaves(refforwardwave,testwave))
    test_leftfft = normalizedfft(truncatewaves(refleftwave,testwave))
    test_rightfft = normalizedfft(truncatewaves(refrightwave,testwave))

    backsimilarity = cosinesimilarity(test_backfft,   backfref_fft)
    forwardsimilarity = cosinesimilarity(test_forwardfft,   forwardref_fft)
    leftsimilarity = cosinesimilarity(test_leftfft,   leftref_fft)
    rightsimilarity = cosinesimilarity(test_rightfft,   rightref_fft)
    def normalizedfft(framelist):
    nUniquePts = ceil((len(framelist)+1)/2.0)


    fftframelist = fft(framelist)
    fftframelist = fftframelist[0:nUniquePts]
    fftframelist = abs(fftframelist)
   # fftframelist=fftframelist/float(len(framelist))
    fftframelist = fftframelist**2
    p = fftframelist
    if len(framelist) % 2 > 0: # we've got odd number of points fft
      p[1:len(p)] = p[1:len(p)] * 2
    else:
      p[1:len(p) -1] = p[1:len(p) - 1] * 2 # we've got even number of points fft
    p = [10*log(x) for x in p]
    return p
    """




main()

