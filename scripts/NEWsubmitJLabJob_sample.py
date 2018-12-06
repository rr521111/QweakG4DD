#////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//                                                                                                        //
#//      Script run jobs using interpolated polarization maps.                                             //
#//      These scripts are separated because the polarization step takes a long time with interpolation.   //
#//                                                                                                        //
#//      Robert Radloff, Ohio University, 2018                                                             //
#//                                                                                                        //
#////////////////////////////////////////////////////////////////////////////////////////////////////////////

#!/apps/python/PRO/bin/python
from subprocess import call,Popen
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
    _workflowID="rradloff_"+idRoot

    for nr in range(_nrStart,_nrStop): # repeat for nr jobs
        _idN= idRoot+'_%04d'% (nr) 
        print _idN
        createMacFile(_directory,_idN,_xP,_yP,_zP,_Px,_Py,_tracking,_beamE,_nEv,nr,modTrj,sample,_pol,_stpSize,nDist,LeadExists)

        call(["cp",_source+"/build/QweakSimG4",_directory+"/"+_idN+"/QweakSimG4"])
        call(["cp",_source+"/myQweakCerenkovOnly.mac",_directory+"/"+_idN+"/myQweakCerenkovOnly.mac"])
	#os.chdir(_directory+"/"+_idN)
	#call(["cd",_directory+"/"+_idN])
	#call(["QweakSimG4","myRun.mac"])

    createXMLfile(_source,_directory,idRoot,_nrStart,_nrStop,_email,sample)


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


def createMacFile(directory,idname,
                  xPos,yPos,zPos,
                  Px,Py,tracking,
                  beamE,nEv,nr,modTrj,sample,pol,stpSize,_nDist,_LeadExists):
    if not os.path.exists(directory+"/"+idname+"/log"):
        os.makedirs(directory+"/"+idname+"/log")
   
    f=open(directory+"/"+idname+"/myRun.mac",'w')
    f.write("/control/execute myQweakCerenkovOnly.mac\n")
    f.write("/PrimaryEvent/SetBeamPositionX "+str(xPos)+" cm\n")
    f.write("/PrimaryEvent/SetBeamPositionY "+str(yPos)+" cm\n")
    f.write("/PrimaryEvent/SetBeamPositionZ "+str(zPos)+" cm\n")
    f.write("/PrimaryEvent/SetBeamDirectionX "+str(Px)+" deg\n")
    f.write("/PrimaryEvent/SetBeamDirectionY "+str(Py)+" deg\n")
    if sample==1:
        f.write("/PrimaryEvent/SetFixedPosMom false\n")
        f.write("/PrimaryEvent/SetPolarization f\n")
    else:
        f.write("/PrimaryEvent/SetFixedPosMom true\n")
        f.write("/PrimaryEvent/SetPolarization "+str(pol)+"\n")
    f.write("/PhysicsProcesses/settingFlag "+str(modTrj)+"\n")    
    f.write("/EventGen/SetBeamEnergy    "+str(beamE)+" MeV\n")
    f.write("/TrackingAction/TrackingFlag "+str(tracking)+"\n")
    f.write("/EventGen/SelectOctant "+str(_nDist%100)+"\n")
    if _LeadExists == 0:
        f.write("/Cerenkov/SetPreradiatorMaterial Air\n")
    seedA=int(time.time()/2000.)+   100*nr+nr
    seedB=int(time.time()/300. ) +10000*nr+nr
    f.write("/Cerenkov/SetPbStepSize "+str(stpSize)+" mm\n");
    f.write("/HallC/GeometryUpdate\n");
    f.write("/random/setSeeds "+str(seedA)+" "+str(seedB)+"\n")
    f.write("/run/beamOn "+str(nEv)+"\n")
    f.close()
    return 0

def createXMLfile(source,writeDir,idRoot,nStart,nStop,email,sample):
    
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
    f.write("  <Command><![CDATA[\n")
    f.write("QweakSimG4 myRun.mac\n")
    f.write("  ]]></Command>\n")
    f.write("  <Memory space=\"2000\" unit=\"MB\"/>\n")

    for nr in range(nStart,nStop): # repeat for nr jobs
        f.write("  <Job>\n")
        idName= writeDir+"/"+idRoot+'_%04d'%(nr)
        f.write("    <Input src=\""+idName+"/QweakSimG4\" dest=\"QweakSimG4\"/>\n")
        f.write("    <Input src=\""+idName+"/myQweakCerenkovOnly.mac\" dest=\"myQweakCerenkovOnly.mac\"/>\n")
        f.write("    <Input src=\""+idName+"/myRun.mac\" dest=\"myRun.mac\"/>\n")
        if sample==1:
            f.write("    <Input src=\""+idName+"/positionMomentum.in\" dest=\"positionMomentum.in\"/>\n")
            f.write("    <Input src=\""+idName+"/polarization.in\" dest=\"polarization.in\"/>\n")
        f.write("    <Output src=\"QwSim_0.root\" dest=\""+idName+"/QwSim_0.root\"/>\n")
        f.write("    <Output src=\"o_tuple.root\" dest=\""+idName+"/o_tuple.root\"/>\n")
        f.write("    <Stdout dest=\""+idName+"/log/log.out\"/>\n")
        f.write("    <Stderr dest=\""+idName+"/log/log.err\"/>\n")
        f.write("  </Job>\n\n")

    f.write("</Request>\n")
    f.close()
    return 0

                    
if __name__ == '__main__':
    main()
                            
