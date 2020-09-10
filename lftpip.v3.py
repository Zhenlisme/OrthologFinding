import os,re,time,subprocess,sys
from Bio import SeqIO
import argparse

class Doliftover():
    def __init__(self,DBdict,querygenome,targetgenome,bedfile,mainpath="",runmem=20,nprocess=30,splitq=30,pslfile=""):
        self.begaintime=time.time()
        self.psl=pslfile
        self.runmem=runmem
        self.splitq=splitq
        self.process=nprocess
        self.anotation_bed=bedfile
        self.qgnom = querygenome
        self.tgnom = targetgenome
        self.qmarker=self.qgnom.split("/")[-1].split(".")[0]
        self.tmarker=self.tgnom.split("/")[-1].split(".")[0]
        self.mainpath="myliftoverdir" if not mainpath else mainpath
        self.newbed="".join([self.qmarker,"_new.bed"])
        self.DBdict=DBdict

        DBfiles=["".join([self.tmarker,'.',i]) for i in ['bck','des','prj','sds','ssp','suf','tis']]
        if not os.path.exists(os.path.abspath(self.DBdict)):
            os.mkdir(os.path.abspath(self.DBdict))
            os.chdir(os.path.abspath(self.DBdict))
            subprocess.run("lastdb %s %s"%(self.tmarker,self.tgnom),shell=True)
        else:
            os.chdir(os.path.abspath(self.DBdict))
            if sum([os.path.isfile(i) for i in DBfiles])!=7:
                os.system("rm -rf *")
                subprocess.run("lastdb %s %s"%(self.tmarker,self.tgnom),shell=True)
        self.DBmarker=self.DBdict.rstrip("/")+"/"+self.tmarker
        sys.stdout.write("The data base name is %s\n" % (self.DBmarker))
        os.chdir("../")            

        if os.path.exists(os.path.abspath(self.mainpath)):
            os.system("rm -rf %s"%(os.path.abspath(self.mainpath)))
            os.mkdir(os.path.abspath(self.mainpath))
            os.chdir(os.path.abspath(self.mainpath))
            sys.stdout.write("The directory '%s' has been re-created.\n"%(os.getcwd()))
        else:
            os.mkdir(os.path.abspath(self.mainpath))
            os.chdir(os.path.abspath(self.mainpath))
        sys.stdout.write("The work directory has been changed into %s\n" %(os.getcwd()))
        
        os.system("faToTwoBit %s %s" % (self.tgnom,self.tmarker+".2bit"))
        os.system("faToTwoBit %s %s" % (self.qgnom,self.qmarker+".2bit"))
        subprocess.run("twoBitInfo %s %s" % (self.tmarker + ".2bit", self.tmarker + ".chromsize"),shell=True)
        subprocess.run("twoBitInfo %s %s" % (self.qmarker + ".2bit", self.qmarker + ".chromsize"),shell=True)
        sys.stdout.write("Two bit files and chromesize files have been created.\n")

    def runsubprocess(self,pbsfile):
        if self.runmem==0:
            pidnumber=subprocess.check_output(["qsub","-l","nodes=1:ppn=1",pbsfile])
        else:
            pidnumber=subprocess.check_output(["qsub","-l","nodes=1:ppn=1,mem=%dgb"%(self.runmem),pbsfile])
       # pidnumber=subprocess.check_output(["qsub","-q","big","-l","nodes=smp1:ppn=1,mem=%dgb"%(self.runmem),pbsfile])
        sys.stdout.write("Task %s was submited.\n"%str(pidnumber).strip("b'\\\\n"))
        return str(pidnumber).strip("b'\\\\n")

    def BeforeBlastzEnd(self):
        os.mkdir("splitbycount")
        os.system("faSplit sequence %s %d splitbycount/"%(self.qgnom,self.splitq))
        pbscontent=(self.doblastz(query="splitbycount/%s"%(qf)) for qf in os.listdir("splitbycount/") if qf.endswith(".fa"))
        return pbscontent

    def Maftopsl(self):
        def trans(maf):
            marker=re.split('[/.]',maf)[-2]
            pbsfile="".join(["lastz_script/",marker,"_topsl.pbs"])
            with open(pbsfile,'w') as F:
                F.write("cd %s\n"%(os.getcwd()))
                F.write("maf-convert psl %s > %s"%(maf,"".join(["lastzf/",marker,".psl"])))
            return pbsfile
        mafiles=["lastzf/"+i for i in os.listdir("lastzf") if i.endswith(".maf")]
        return (trans(f) for f in mafiles)
        

    def doblastz(self,query):
        if not os.path.exists(os.path.abspath("lastz_script")):
            os.mkdir(os.path.abspath("lastz_script"))
        if not os.path.exists(os.path.abspath("lastzf")):
            os.mkdir(os.path.abspath("lastzf"))
        pbsmarker=re.match('splitbycount/(.+?).fa',query).group(1)
        pbsfile="".join(["lastz_script/",self.qmarker,"_",pbsmarker,'_lastz.pbs'])
        sys.stdout.write("The pbsfile is %s\n"%(pbsfile))
        opfile="".join(["lastzf/",self.qmarker,"_",pbsmarker,'.maf'])
        commandline="lastal %s %s > %s"%(self.DBmarker,query,opfile)
        with open(pbsfile,'w') as F:
            F.write("cd %s\n"%(os.getcwd()))
            F.write(commandline)
        return pbsfile

    def pbsmanage(self,pbscontent):
        Totaltime=0
        pidcontent=[]
        while True:
            try:
                pidcontent.append(self.runsubprocess(next(pbscontent)))
                while len(pidcontent)==self.process:
                    [pidcontent.remove(i) for i in pidcontent if not re.findall("\n%s\s.+?\s[QR]\s.+" % (i), subprocess.getoutput("qstat"))]
                    time.sleep(120)
                    Totaltime+=120
            except:
                break
        while pidcontent:
            time.sleep(120)
            Totaltime+=120
            [pidcontent.remove(i) for i in pidcontent if not re.findall("\n%s\s.+?\s[QR]\s.+"%(i),subprocess.getoutput("qstat"))]
        return Totaltime

    def AfterBlastzEnd(self):
        if not self.psl:
            pslfile="".join([self.tmarker,"_",self.qmarker,".psl"])
            subprocess.run("cat lastzf/*.psl > %s"%(pslfile),shell=True)
        else:
            pslfile=self.psl
        ##convert to net
        chainfile="".join([self.qmarker,"_",self.tmarker,".chain"])
        subprocess.run("axtChain -linearGap=medium -psl %s %s %s %s"%\
                      (pslfile,self.tmarker+".2bit",self.qmarker+".2bit",chainfile),shell=True)
        ##merge chain and splited them
        if not os.path.exists("splited_chain/"):
            os.mkdir("splited_chain/")
        subprocess.run("chainMergeSort %s|chainSplit %s stdin"%(chainfile,os.path.abspath("splited_chain/")),shell=True)
        ##make net chain
        overchainf="To".join([self.tmarker,self.qmarker])+".over.chain"
        [subprocess.run("chainNet %s %s %s stdout /dev/null|netChainSubset stdin %s stdout >> %s" %
        ("splited_chain/"+chainf,self.tmarker + ".chromsize",self.qmarker + ".chromsize","splited_chain/"+chainf,overchainf),shell=True)
         for chainf in os.listdir("splited_chain")]
        ##swap chain to be query-referenced
        tbestchainf = "_".join([self.tmarker, self.qmarker]) + ".tBest.chain"
        subprocess.run("chainStitchId %s stdout|chainSwap stdin stdout|chainSort stdin %s" % (overchainf,tbestchainf),shell=True)
        ##Net those on query(non-xxx) to get query-refd reciprocal best net
        rbestnetf="_".join([self.tmarker, self.qmarker])+".rbest.net"
        subprocess.run("chainPreNet %s %s %s stdout|chainNet -minSpace=1 stdin %s %s stdout /dev/null |netSyntenic stdin %s"%(tbestchainf,self.qmarker + ".chromsize",self.tmarker + ".chromsize",self.qmarker + ".chromsize",self.tmarker + ".chromsize",rbestnetf),shell=True)
        ##Extract query-refd reciprocal best chain
        rbestchainf="_".join([self.tmarker, self.qmarker]) + ".rBest.chain"
        subprocess.run("netChainSubset %s %s stdout|chainStitchId stdin %s"%(rbestnetf,tbestchainf,rbestchainf),shell=True)
        ##Swap to get xxx-refd reciprocal best chain
        sortrbestchainf="_".join([self.tmarker, self.qmarker])+".sorted.rbest.chain"
        subprocess.run("chainSwap %s stdout|chainSort stdin %s"%(rbestchainf,sortrbestchainf),shell=True)
        ##Net those on xxx to get xxx-refd recipical best net
        rbestnetf2 = "_".join([self.tmarker, self.qmarker]) + ".rbest.net2"
        subprocess.run("chainPreNet %s %s %s stdout|chainNet -minSpace=1 -minScore=0 stdin %s %s stdout /dev/null |netSyntenic stdin %s"%\
        (sortrbestchainf,self.tmarker + ".chromsize",self.qmarker + ".chromsize",self.tmarker + ".chromsize",self.qmarker + ".chromsize",rbestnetf2),shell=True)
        os.system("mkdir recip_chain/")
        subprocess.run("chainMergeSort %s|chainSplit %s stdin"%(sortrbestchainf,"recip_chain/"),shell=True)
        opchainfile="To".join([self.tmarker,self.qmarker])+"_recipl.over.chain"
        [subprocess.run("chainNet %s %s %s stdout /dev/null|netChainSubset stdin %s stdout >> %s" %
                        ("recip_chain/" + chainf, self.tmarker + ".chromsize", self.qmarker + ".chromsize",
                         "recip_chain/" + chainf, opchainfile),shell=True)
        for chainf in os.listdir("recip_chain/")]
        return opchainfile
    def liftover(self,chainf):
        subprocess.run("liftOver %s %s %s /dev/null"%(self.anotation_bed,chainf,self.newbed),shell=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="To Run Liftover Program Automatically. Writen by Zhen Li at Jan, 4, 2020")
    parser.add_argument("-q", "--querygenome", required=True, type=str,
                        help="The query genome sequence which should be fasta format.")
    parser.add_argument("-db","--dbdict",required=True,type=str,
                        help="The db directory name.")
    parser.add_argument("-psl","--psl",type=str,default="",help="Whether use the psl format file directly.")
    parser.add_argument("-t", "--targetgenome", required=True, type=str,
                        help="The target genome sequence which should be fasta format.")
    parser.add_argument("-b", "--anotationfile", required=True, type=str,
                        help="The anotation file which should be bed format.")
    parser.add_argument("-sq","--splitq",type=int,default=30,help="The number that the query genome should be splited into")
    parser.add_argument("-d","--optdirect",type=str,
                        help="The destinated directory in where all of your output will be stored.",default="")
    parser.add_argument("-p","--process",default=0,type=int,help="The maximum process number that could run at the same time.")
    parser.add_argument("-m","--runmem",type=int,help="The memory size(gb) that each task will be used.",default=20)
    Args=parser.parse_args()

    if not Args.psl:
        Mainprogram=Doliftover(querygenome=os.path.abspath(Args.querygenome),DBdict=os.path.abspath(Args.dbdict),targetgenome=os.path.abspath(Args.targetgenome),bedfile=os.path.abspath(Args.anotationfile),mainpath=Args.optdirect,runmem=Args.runmem,nprocess=Args.process,splitq=Args.splitq)
        pbscontent=Mainprogram.BeforeBlastzEnd()
        checkpoint=Mainprogram.pbsmanage(pbscontent)
        sys.stdout.write("%.3f hours used by lastz.\n"%round(checkpoint/3600,3))
        mafiles=Mainprogram.Maftopsl()
        checkpoint=Mainprogram.pbsmanage(mafiles)
        sys.stdout.write("%.3f hours used by maf-convert.\n"%round(checkpoint/3600,3))
    else:
        Mainprogram=Doliftover(querygenome=os.path.abspath(Args.querygenome),DBdict=os.path.abspath(Args.dbdict),targetgenome=os.path.abspath(Args.targetgenome),bedfile=os.path.abspath(Args.anotationfile),mainpath=Args.optdirect,runmem=Args.runmem,nprocess=Args.process,splitq=Args.splitq,pslfile=os.path.abspath(Args.psl))

    chainfile=Mainprogram.AfterBlastzEnd()
    Mainprogram.liftover(chainf=chainfile)
    sys.stdout.write("%.3f hours used totally.\n"%round((time.time()-Mainprogram.begaintime)/3600,3))

