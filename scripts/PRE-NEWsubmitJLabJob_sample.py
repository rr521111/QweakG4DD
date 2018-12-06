#////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//                                                                                                        //
#//      Script to generate polarization in files for a set of jobs. Should be followed by NEW-submit...   //
#//      These scripts are separated because the polarization step takes a long time with interpolation.   //
#//                                                                                                        //
#//      Robert Radloff, Ohio University, 2018                                                             //
#//                                                                                                        //
#////////////////////////////////////////////////////////////////////////////////////////////////////////////

#!/apps/python/PRO/bin/python
from subprocess import call
import sys,os,time

def main():
    
    #center, x,y,z=0,335,560
    _xP=0.
    _yP=335.0
    _zP=560.
    _Px=0.#deg
    _Py=0.
    _beamE=1160#MeV
    _email="rradloff@jlab.org"
    _source="/u/home/rradloff/dd/QweakG4DD"
    _directory="/volatile/hallc/qweak/rradloff"
    _tracking=0 #0=primary only | 1=prim + opt photon | 2=no optical ph and 10x faster than 3=full
    _stpSize=0.02
    _nEv=10000
    _nrStop=9999
    _nrStart=0
    _pol="V"
    modTrj=0 ## 0:standard G4 propagation(wght sims) 1:debug print == big NONO! 2: modifyTraj
    submit=1
    nDist=208
    sample=1
    useSWIF=1 #0: uses jsub 1: uses SWIF+jsub
    LeadExists=1 #0: removes preradiator
    
    idRoot= _pol+'_sampled_'+str(nDist)+'_%03dk_'% (_nEv/1000)
    _workflowID="rradloff_"+idRoot+"_polmake"

    for nr in range(_nrStart,_nrStop): # repeat for nr jobs
        _idN= idRoot+'_%04d'% (nr) 
        print _idN

	if sample == 1:
            call(["mkdir",_directory+"/"+_idN])

    createXMLfile(_source,_directory,idRoot,_nrStart,_nrStop,_email,sample,_nEv,nDist,_pol)


    if submit==1:
        if useSWIF==1:
            print "submitting position sampled with id",_idN," between ",_nrStart,_nrStop," using designated SWIF workflow "+str(_workflowID)
            call(["swif","add-jsub","-workflow",str(_workflowID),"-create","-script",_directory+"/"+idRoot+".xml"])
        elif useSWIF==0:
            print "submitting position sampled with id",_idN," between ",_nrStart,_nrStop
            call(["jsub","-xml",_directory+"/"+idRoot+".xml"])
    else:
        print "NOT submitting position sampled with id",_idN," between ",_nrStart,_nrStop
        
    print "I am all done"

def createXMLfile(source,writeDir,idRoot,nStart,nStop,email,sample,nEv,_nDist,pol):
    
    if not os.path.exists(source+"/scripts/jobs"):
        os.makedirs(source+"/scripts/jobs")

    f=open(writeDir+"/"+idRoot+".xml","w")
    f.write("<Request>\n")
    f.write("  <Email email=\""+email+"\" request=\"false\" job=\"true\"/>\n")
    f.write("  <Project name=\"qweak\"/>\n")
#    f.write("  <Track name=\"debug\"/>\n")
    f.write("  <Track name=\"simulation\"/>\n")
    f.write("  <Name name=\""+idRoot+"\"/>\n")
    f.write("  <OS name=\"centos7\"/>\n")
    f.write("  <Memory space=\"2000\" unit=\"MB\"/>\n")

    for nr in range(nStart,nStop): # repeat for nr jobs
        f.write("  <Job>\n")
        idName= writeDir+"/"+idRoot+'_%04d'%(nr)
        f.write("    <Command><![CDATA[\n")
        f.write(""+source+"/build/QweakSimRoot -l -b .L "+source+"/rootScripts/samplePrimaryDist3.C\\("+str(nEv)+",1,"+str(_nDist)+",\\\""+str(pol)+"\\\",\\\""+idName+"\\\"\\)")
        f.write("    ]]></Command>\n")
        f.write("    <Stdout dest=\""+idName+"/log/log.out\"/>\n")
        f.write("    <Stderr dest=\""+idName+"/log/log.err\"/>\n")
        f.write("  </Job>\n\n")

    f.write("</Request>\n")
    f.close()
    return 0

                    
if __name__ == '__main__':
    main()
                            
