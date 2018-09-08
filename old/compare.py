'''
Utility functions to compare different bounds 
old code (~2008), not cleaned-up.
'''

from curve import *

def compareList(D):
    '''
    For a given list L of n bounds(in table form), output a n by n table with
    T[i][j]=#codes where L[i]>L[j]
    '''
    L=[]
    K=[]
    for i in range(len(D)):
        New=[]
        Newimp=[]
        for j in range(len(D)):
            cnt=0
            maximp=0
            for deg in range(len(D[0])):
                for k in range(len(D[0][0])):
                    maximp = max (maximp,D[j][deg][k]-D[i][deg][k])
                    cnt=cnt+(D[j][deg][k]>D[i][deg][k])
            New.append(cnt)
            Newimp.append(maximp)
        L.append(New)
        K.append(Newimp)
    return L,K

def degOptimal(L,showDivisors=False,fullList=False):
    '''
    Maximal distance per degree according to all bounds
    bounds in the passed list that achieve it
    NOTE:Since the name of the bound is only stored in the variable name
    the output just gives the place in the original list
    Ex. degOptimal([DLM,DDP,DMT]) will return at 13 
    13,15,[2]
    This means at deg=13, the maximal code has distance 15, and only the DMT bound reaches 15.
    '''
    LDEG=[]
    for deg in range(len(L[0])):
        maxD=[max(L[i][deg]) for i in range(len(L))]      
        maxmaxD=max(maxD)
        if(fullList):
            List=[]
            for i in range(len(maxD)):
                if(maxD[i]==maxmaxD):
                    List.append(i)
        else:
            for i in range(len(maxD)):
                if(maxD[i]==maxmaxD):
                    List=i
                    break
        if(showDivisors):
            LDEG.append([deg,maxmaxD,[i for i,val in enumerate(L[-1][deg]) if maxmaxD==val],List])        
        else:
            LDEG.append([deg,maxmaxD,List])
    return LDEG

def Compare2(D1,D2,curve,sdiff=-1):
    '''
    The fuction requires that D1 <= D2, otherwise returns an error
    The output is 1D tables L and M with
    L[i]=#codes with D2-D1=i
    M[i]=#maximum deg code with D2-D1=i
    '''
    m=curve.m
    lD=min([len(D1),len(D2)])
    Diff=[[D2[deg][i]-D1[deg][i] for i in range(m)] for deg in range(lD)]
    mx=max([max(i) for i in Diff])
    print 'max diff between D1 and D2 ',mx
    L=[0 for i in range(mx+1)]
    M=[0 for i in range(mx+1)]
    for deg in range(lD):
        for i in range(m):
            k=Diff[deg][i]
            if k>=0:
                L[k]=L[k]+1
                M[k]=deg
            else:
                print deg,i,"problem D1>D2"
    if(sdiff!=-1):
        select=[]
        for deg in range(lD):
            for i in range(m):
                if(Diff[deg][i]==sdiff):
                    select.append([deg,i])
        print select   
    return L,M

def Compare3(D1,D2,D3,sdiff1=-1,sdiff2=-1):  
    '''
    The fuction requires that D1 <= D2 <= D3, otherwise returns an error
    The output is 2D tables L and M with
    L[i][j]=#codes with D2-D1=i and D3-D2=j
    M[i][j]=#maximum deg code with D2-D1=i and D3-D2=j
    use sdiff1=i and sdiff2=j to see the codes which have exactly D1+i=D2 and D2+j=D3
    '''
    lD=min([len(D1),len(D2),len(D3)])
    Diff1=[[D2[deg][i]-D1[deg][i] for i in range(m)] for deg in range(lD)]
    Diff2=[[D3[deg][i]-D2[deg][i] for i in range(m)] for deg in range(lD)]
    mx1=max([max(i) for i in Diff1])
    mx2=max([max(i) for i in Diff2])
    print 'max diff between D1 and D2 ',mx1
    print 'max diff between D2 and D3 ',mx2
    L=[[0 for j in range(mx2+1)] for i in range(mx1+1)]
    M=[[0 for j in range(mx2+1)] for i in range(mx1+1)]
    for deg in range(lD):
        for i in range(m):
            k=Diff1[deg][i]
            l=Diff2[deg][i]
            if(k<0 or l<0):
                print deg,i," problem D1>D2 or D2>D3"
                returnrkirov
            L[k][l]=L[k][l]+1
            M[k][l]=deg
    if(sdiff1!=-1 or sdiff2!=-1):
        select=[]
        for deg in range(lD):
            for i in range(m):
                if((sdiff1==-1 or Diff1[deg][i]==sdiff1) and (Diff2[deg][i]==sdiff2 or sdiff2==-1)):
                    select.append([deg,i])
        print select    
    return L,M
