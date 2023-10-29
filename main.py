import copy

########################################################################################################################################
########################################################################################################################################
# Amino acid monomers in polypeptides
testnode = ['Cys', 'Cys', 'Glu', 'Tyr', 'Cys', 'Cys', 'Asn', 'Pro', 'Ala', 'Cys', 'Thr', 'Gly', 'Cys', 'Tyr']
# Bond in a polypeptide where the first two digits of each list are the sequence number of the amino acid at either end of the bond, and the third digit is the bond type.
#1:peptide bond, 2:disulfide bond
testbond = [[1, 2, 1], [2, 3, 1], [3, 4, 1], [4, 5, 1], [5, 6, 1], [6, 7, 1], [7, 8, 1], [8, 9, 1], [9, 10, 1],
                    [10, 11, 1], [11, 12, 1], [12, 13, 1], [13, 14, 1], [1, 6, 2], [2, 10, 2], [5, 13, 2]]
########################################################################################################################################
########################################################################################################################################

Piptideweight={'Ala':89.0477, 'Arg':174.1117, 'Asn':132.0535, 'Asp':133.0375, 'Cys':121.0198, 'Gln':146.0691, 'Glu':147.0532, 'Gly':75.032, 'His':155.0695, 'Ile':131.0946,
               'Leu':131.0946, 'Lys':146.1055, 'Met':149.0511, 'Phe':165.079, 'Pro':115.0633, 'Ser':105.0426, 'Thr':119.0582, 'Trp':204.0899, 'Tyr':181.0739, 'Val':117.079}
BondTypeWeight={1:-18.01056, 2:-2.01565}

class PiptedStructure:
    def __init__(self,testnode,testbond):
        self.PiptideNode=testnode
        self.PiptideBond=testbond

class SuperNode:
    def __init__(self, NodeGraph):
        self.aa = NodeCount(NodeGraph)
        self.a=[[],self.aa[1]]
        for ai in self.aa[0]:
            self.a[0].append(ai[0])
        self.NodeNumber = len(self.a[1])
        self.FragmentList = []
        self.FragmentMonomerNumberList = []
        self.FragmentEndNodeList = []
        for i in range(0,len(self.a[0])):
            self.FragmentList.append([i+1])
            self.FragmentMonomerNumberList.append(len(self.a[0][i])-1)
            self.FragmentEndNodeList.append([0,0])
            for j in range(0,len(self.a[1])):
                if self.a[0][i][0]==self.a[1][j]:
                    self.FragmentEndNodeList[-1][0]=j+1
                if self.a[0][i][-1]==self.a[1][j]:
                    self.FragmentEndNodeList[-1][1] = j+1
        self.FragmentLinkList = []
        for ii in range(0,len(self.FragmentList)):
            for jj in range(0,ii):
                if self.FragmentEndNodeList[ii][0]== self.FragmentEndNodeList[jj][0] or self.FragmentEndNodeList[ii][0]== self.FragmentEndNodeList[jj][1] or self.FragmentEndNodeList[ii][1]== self.FragmentEndNodeList[jj][0] or self.FragmentEndNodeList[ii][1]== self.FragmentEndNodeList[jj][1]:
                    self.FragmentLinkList.append([ii+1,jj+1])
        self.IsFragmentCircle=IsFragmentCircle(self.FragmentEndNodeList)

def NodeCount(NodeGraph):
    PreNodeList=[]
    for i in range(0,len(NodeGraph.PiptideNode)):
        PreNodeList.append([i+1])
    for j in NodeGraph.PiptideBond:
        if len(PreNodeList[j[0]-1])==1:
            PreNodeList[j[0]-1].append(1)
        else:
            PreNodeList[j[0] - 1][1]+=1
        if len(PreNodeList[j[1]-1])==1:
            PreNodeList[j[1]-1].append(1)
        else:
            PreNodeList[j[1] - 1][1]+=1
    NodeNumber = 0
    NodeReflection=[]
    for i in PreNodeList:
        if i[1]!=2:
            NodeNumber+=1
            NodeReflection.append(i[0])
    PreFragmentlist=[]
    for j in NodeReflection:
        NodePath=[[j],[]]
        for k in range(0,len(P.PiptideBond)):
            if P.PiptideBond[k][0]==NodePath[0][0]:
                NodePath[0].append(P.PiptideBond[k][1])
                NodePath[1].append(k)
                if NodePath[0][-1] not in NodeReflection:
                    PreFragmentlist.append(NodePathFind(NodePath, P.PiptideBond, NodeReflection))
                else:
                    PreFragmentlist.append(NodePath)
            if P.PiptideBond[k][1]==NodePath[0][0]:
                NodePath[0].append(P.PiptideBond[k][0])
                NodePath[1].append(k)
                if NodePath[0][-1] not in NodeReflection:
                    PreFragmentlist.append(NodePathFind(NodePath, P.PiptideBond, NodeReflection))
                else:
                    PreFragmentlist.append(NodePath)
            NodePath=[[j],[]]
    c=len(PreFragmentlist)
    for ii in range(0,c):
        if PreFragmentlist[c-ii-1][0][0]>PreFragmentlist[c-ii-1][0][-1]:
            PreFragmentlist.remove(PreFragmentlist[c-ii-1])
        elif PreFragmentlist[c-ii-1][0][0]==PreFragmentlist[c-ii-1][0][-1]:
            if PreFragmentlist[c-ii-1][0][1]>PreFragmentlist[c-ii-1][0][-2]:
                PreFragmentlist.remove(PreFragmentlist[c - ii - 1])

    NewPreFragmentlist=[]
    for jj in PreFragmentlist:
        jjj=copy.deepcopy(jj)
        jjjj=jjj[1][0]
        jjj[1][0]=jjj[1][-1]
        jjj[1][-1]=jjjj
        if jj not in NewPreFragmentlist and jjj not in NewPreFragmentlist:
            NewPreFragmentlist.append(jj)
    PreFragmentlist=NewPreFragmentlist
    return [PreFragmentlist,NodeReflection]

def NodePathFind(NodePath,link,NodeReflection):
    EndJudge = True
    while EndJudge == True:
        for k in range(0,len(link)):
            if k not in NodePath[1]:
                if link[k][1] not in NodePath[0][1:] and link[k][0]==NodePath[0][-1]:
                    NodePath[0].append(link[k][1])
                    NodePath[1].append(k)
                    break
                if link[k][0] not in NodePath[0][1:] and link[k][1]==NodePath[0][-1]:
                    NodePath[0].append(link[k][0])
                    NodePath[1].append(k)
                    break
        if NodePath[0][-1] in NodeReflection:
            EndJudge=False
    return NodePath

def IsFragmentCircle(FragmentEndNodeList):
    IsFragmentCircle=[]
    for c in FragmentEndNodeList:
        if c[0] == c[1]:
            IsFragmentCircle.append(1)
        else:
            IsFragmentCircle.append(0)
    return IsFragmentCircle

def CreatAllFragment(FragmentList):
    FragmentSubgraph = []
    for i in range(0,len(FragmentList)):
        FragmentSubgraph.append([])
        if i == 0:
            FragmentSubgraph[i]=copy.deepcopy(FragmentList)
        else:
            for j in FragmentSubgraph[i-1]:
                TemporaryAddNodeList=[]
                for k in j:
                    for addnode in PeptideHydrolyze.FragmentLinkList:
                        if k==addnode[0]:
                            if addnode[1] not in j and addnode[1] not in TemporaryAddNodeList:
                                TemporaryAddNodeList.append(addnode[1])
                        if k==addnode[1]:
                            if addnode[0] not in j and addnode[0] not in TemporaryAddNodeList:
                                TemporaryAddNodeList.append(addnode[0])
                j.append(TemporaryAddNodeList)
                for TemporaryAddNode in TemporaryAddNodeList:
                    jcopy=j.copy()
                    jcopy.pop()
                    jcopy.append(TemporaryAddNode)
                    jadd = jcopy
                    jadd.sort()
                    if jadd not in FragmentSubgraph[i]:
                        FragmentSubgraph[i].append(jadd)
    return FragmentSubgraph

def GetAllPreStucture(FragmentSubgraph):
    PreStuctureList=[]
    for i in range(0,len(FragmentSubgraph)-1):
        for j in FragmentSubgraph[i]:
            x=[]
            y=[]
            z=j[-1]
            for k in range (0,len(j)-1):
                x.append(j[k])
                y.append(PeptideHydrolyze.FragmentEndNodeList[j[k]-1][0])
                y.append(PeptideHydrolyze.FragmentEndNodeList[j[k] - 1][1])
            y=list(set(y))
            PreStuctureList.append([x,y,z])
    PreStuctureList.append([FragmentSubgraph[-1][0],[],[]])
    for f in PeptideHydrolyze.FragmentList:
        PreStuctureList.append([[],[],f])
    for i in range(0, PeptideHydrolyze.NodeNumber):
        x=[]
        y=[i+1]
        z=[]
        for zi in range(0,len(PeptideHydrolyze.FragmentEndNodeList)):
            if PeptideHydrolyze.FragmentEndNodeList[zi][0]==y[0] or PeptideHydrolyze.FragmentEndNodeList[zi][1]==y[0]:
                z.append(zi+1)
        PreStuctureList.append([x,y,z])
    return PreStuctureList

def RoundFragmentLinkType(TypeSum,FragmentSubgraph):
    for i in range(0,len(FragmentSubgraph)-1):
        for j in FragmentSubgraph[i]:
            ExistNodeList=[]
            typelist=[]
            for k in range(0,len(j)-1):
                if PeptideHydrolyze.FragmentEndNodeList[j[k]-1][0] not in ExistNodeList:
                    ExistNodeList.append(PeptideHydrolyze.FragmentEndNodeList[j[k]-1][0])
                if PeptideHydrolyze.FragmentEndNodeList[j[k]-1][1] not in ExistNodeList:
                    ExistNodeList.append(PeptideHydrolyze.FragmentEndNodeList[j[k]-1][1])
            num=1
            for f in j[-1]:
                if PeptideHydrolyze.FragmentMonomerNumberList[f-1]==0:
                    typelist.append(0)
                    numx=1
                else:
                    if PeptideHydrolyze.FragmentEndNodeList[f-1][0] in ExistNodeList and PeptideHydrolyze.FragmentEndNodeList[f-1][1] in ExistNodeList:
                        typelist.append(2)
                        numx = MultipleNumber(2,PeptideHydrolyze.FragmentMonomerNumberList[f-1])
                    else:
                        typelist.append(1)
                        numx = MultipleNumber(1, PeptideHydrolyze.FragmentMonomerNumberList[f-1])
                num*=numx
            TypeSum += num
    TypeSum+=1
    return TypeSum
def MultipleNumber(a,b):
    if a==2:
        return int(b*(b+1)/2)
    else:
        return int(b)

def SigleFragmentCal(TypeSum):
    for i in range(0,len(PeptideHydrolyze.FragmentList)):
        TypeSum+=int((PeptideHydrolyze.FragmentMonomerNumberList[i]*(PeptideHydrolyze.FragmentMonomerNumberList[i]-1))/2)
    return TypeSum

def SigleNodeCal(TypeSum):
    for i in range(0,PeptideHydrolyze.NodeNumber):
        PossibleFragment = []
        for j in range(0,len(PeptideHydrolyze.FragmentEndNodeList)):
            if PeptideHydrolyze.FragmentEndNodeList[j][0]-1==i and j not in PossibleFragment:
                PossibleFragment.append(j)
            if PeptideHydrolyze.FragmentEndNodeList[j][1]-1==i and j not in PossibleFragment:
                PossibleFragment.append(j)
        num=1
        for k in PossibleFragment:
            if PeptideHydrolyze.IsFragmentCircle[k]==0:
                num*=(PeptideHydrolyze.FragmentMonomerNumberList[k])
            else:
                num*=int((PeptideHydrolyze.FragmentMonomerNumberList[k])*int(PeptideHydrolyze.FragmentMonomerNumberList[k]+1)/2)
        TypeSum+=num
    return TypeSum

def TotalFragmentCount(FragmentSubgraph):
    TypeSum = 0
    TypeSum = RoundFragmentLinkType(TypeSum,FragmentSubgraph)
    TypeSum = SigleFragmentCal(TypeSum)
    TypeSum = SigleNodeCal(TypeSum)
    return TypeSum

def GetPeptideCombination(PreStuctureList):
    FinalI=[]
    MiddleStuctureList=[]
    for i in PreStuctureList:
        if i[0]==[]:
            if i[1]==[]:
                a=SignalFrengmentGenerate(PeptideHydrolyze.a[0][i[2][0]-1],PeptideHydrolyze.aa[0][i[2][0]-1][1])
                for aa in a:
                    if aa!=[]:
                        MiddleStuctureList.append(aa)
            else:
                a=SignalNodeFrengmentGenerate(i[1],i[2])
                for aa in a:
                    if aa!=[]:
                        MiddleStuctureList.append(aa)
        else:
            a=MultipleFrengmentGenerate(i[0],i[1],i[2])
            for aa in a:
                if aa != []:
                    MiddleStuctureList.append(aa)
    for end in MiddleStuctureList:
        end.append(WeightCalculation(end,P.PiptideNode,P.PiptideBond))
    for ei in range(0,len(MiddleStuctureList)):
        for ej in range(ei,len(MiddleStuctureList)):
            if MiddleStuctureList[ei][-1]<MiddleStuctureList[ej][-1]:
                ek=MiddleStuctureList[ej]
                MiddleStuctureList[ej]=MiddleStuctureList[ei]
                MiddleStuctureList[ei]=ek
    for endtest in MiddleStuctureList:
        endtest[0].sort()
        endtest[1].sort()
    asb=[]
    for asbi in MiddleStuctureList:
        # print(asbi)
        I=GetIUPAC(asbi)
        if I not in FinalI:
            FinalI.append(I)
    for ff in FinalI:
        print(ff[0],ff[1])
        #the result

def GetIUPAC(I):
    Nodee=I[0]
    Bonde=I[1]
    w=I[2]
    NewBond=[]
    for bi in Bonde:
        NewBond.append(P.PiptideBond[bi])
    SSbond=[]
    for si in NewBond:
        if si[2]==2:
            SSbond.append(si)
    WorkList=[]
    for i in range(0,len(Nodee)):
        p=P.PiptideNode[Nodee[i]-1]
        for j in range(0,len(SSbond)):
            if SSbond[j][0]==Nodee[i] or SSbond[j][1]==Nodee[i]:
                p=p+'('+str(j+1)+')'
        if i == 0:
            WorkList.append([p])
        else:
            Link=False
            for bi in NewBond:
                if bi[2]==1:
                    if bi[0]==Nodee[i-1] and bi[1]==Nodee[i]:
                        Link=True
            if Link==True:
                WorkList[-1].append(p)
            else:
                WorkList.append([])
                WorkList[-1].append(p)
    ei=''
    for wi in WorkList:
        ei=ei+'H-'
        for wii in wi:
            ei=ei+wii+'-'
        ei+='OH.'
    ei=ei[:-1]
    return [ei,w]

def SignalFrengmentGenerate(a,abond):
    aresult=[]
    if len(a)<=2:
        aresult.append([])
    else:
        for i in range(1,len(a)):
            for j in range(i+1,len(a)):
                w=a[i:j]
                wb = abond[i:j-1]
                aresult.append([w, wb])
    return aresult

def SignalNodeFrengmentGenerate(Node,Frengment):
    aresult=[]
    if len(Frengment)==0:
        pass
    else:
        MultipleCombination=[]
        for mi in Frengment:
            MultipleNodeJudgment=1
            MultipleCombinationi=[]
            a = copy.deepcopy(PeptideHydrolyze.a[0][mi - 1])
            n = PeptideHydrolyze.a[1][Node[0] - 1]
            abond = copy.deepcopy(PeptideHydrolyze.aa[0][mi - 1][1])
            if n == a[0] and n==a[-1]:
                MultipleNodeJudgment=0
            if n!=a[0]:
                MultipleNodeJudgment=-1
            if MultipleNodeJudgment==1:
                for msi in range(1, PeptideHydrolyze.FragmentMonomerNumberList[mi - 1] + 1):
                    w = a[0:msi]
                    wb = abond[0:msi - 1]
                    MultipleCombinationi.append([w, wb])
                if MultipleCombinationi!=[]:
                    MultipleCombination.append(MultipleCombinationi)
            elif MultipleNodeJudgment==-1:
                a.reverse()
                abond.reverse()
                for msi in range(1, PeptideHydrolyze.FragmentMonomerNumberList[mi - 1] + 1):
                    w = a[0:msi]
                    wb = abond[0:msi - 1]
                    MultipleCombinationi.append([w, wb])
                if MultipleCombinationi!=[]:
                    MultipleCombination.append(MultipleCombinationi)
            else:

                for msi in range(0, PeptideHydrolyze.FragmentMonomerNumberList[mi - 1]):
                    for msj in range(0,msi+1):
                        w = []
                        wb = []
                        waa=a[msi+1:len(a)]
                        wab=a[0:msj+1]
                        wba=abond[msi + 1:len(abond)]
                        wbb=abond[0:msj]
                        for waai in waa:
                            if waai not in w:
                                w.append(waai)
                        for wabi in wab:
                            if wabi not in w:
                                w.append(wabi)
                        for wbai in wba:
                            if wbai not in wb:
                                wb.append(wbai)
                        for wbbi in wbb:
                            if wbbi not in wb:
                                wb.append(wbbi)
                        MultipleCombinationi.append([w, wb])
                if MultipleCombinationi!=[]:
                    MultipleCombination.append(MultipleCombinationi)
        MultipleControl=[]
        MultipleControlLim=[]
        for mm in MultipleCombination:
            MultipleControl.append(0)
            MultipleControlLim.append(len(mm))
        MultipleControlLim.append(1)
        MultipleControlLimMultiplication=[]
        la = 1
        for li in range(0,len(MultipleControlLim)-1):
            la*=MultipleControlLim[len(MultipleControlLim)-li-1]
            MultipleControlLimMultiplication.append(la)
        SumType = 1
        for si in MultipleControlLim:
            SumType *= si
        for wi in range(0,SumType):
            for wii in range(0,len(MultipleControl)):
                MultipleControl[wii]=int(((wi - (wi % MultipleControlLimMultiplication[len(MultipleControl) - wii - 1]))/MultipleControlLimMultiplication[len(MultipleControl) - wii - 1])%MultipleControlLim[wii])
            apreresult=[]
            for ma in range(0,len(MultipleControl)):
                apreresult.append(MultipleCombination[ma][MultipleControl[ma]])
            newapr=[[],[]]
            for aa in apreresult:
                if aa[0]!=[]:
                    for aaa in aa[0]:
                        if aaa not in newapr[0]:
                            newapr[0].append(aaa)
                if aa[1]!=[]:
                    for aaa in aa[1]:
                        if aaa not in newapr[1]:
                            newapr[1].append(aaa)
            aresult.append(newapr)
    return aresult

def MultipleFrengmentGenerate(FullFrengment,Node,SemiFrengment):
    aresult=[]
    fsatomic=[]
    fbound=[]
    nlist=[]
    for ni in Node:
        nlist.append(PeptideHydrolyze.a[1][ni - 1])
    for ff in FullFrengment:
        fa = copy.deepcopy(PeptideHydrolyze.a[0][ff - 1])
        fabond = copy.deepcopy(PeptideHydrolyze.aa[0][ff - 1][1])
        for fai in fa:
            if fai not in fsatomic:
                fsatomic.append(fai)
        for fabondi in fabond:
            if fabondi not in fbound:
                fbound.append(fabondi)
    MultipleCombination = []
    for mi in SemiFrengment:
        MultipleNodeJudgment = 1
        MultipleCombinationi = []
        a = copy.deepcopy(PeptideHydrolyze.a[0][mi - 1])
        abond = copy.deepcopy(PeptideHydrolyze.aa[0][mi - 1][1])
        if a[0] in nlist and a[-1] in nlist:
            MultipleNodeJudgment = 0
        if a[0] not in nlist:
            MultipleNodeJudgment = -1
        if MultipleNodeJudgment == 1:
            for msi in range(1, PeptideHydrolyze.FragmentMonomerNumberList[mi - 1] + 1):
                w = a[0:msi]
                wb = abond[0:msi - 1]
                MultipleCombinationi.append([w, wb])
            if MultipleCombinationi != []:
                MultipleCombination.append(MultipleCombinationi)
        elif MultipleNodeJudgment == -1:
            a.reverse()
            abond.reverse()
            for msi in range(1, PeptideHydrolyze.FragmentMonomerNumberList[mi - 1] + 1):
                w = a[0:msi]
                wb = abond[0:msi - 1]
                MultipleCombinationi.append([w, wb])
            if MultipleCombinationi != []:
                MultipleCombination.append(MultipleCombinationi)
        else:
            for msi in range(0, PeptideHydrolyze.FragmentMonomerNumberList[mi - 1]):
                for msj in range(0, msi + 1):
                    w = []
                    wb = []
                    waa = a[msi + 1:len(a)]
                    wab = a[0:msj + 1]
                    wba = abond[msi + 1:len(abond)]
                    wbb = abond[0:msj]
                    for waai in waa:
                        if waai not in w:
                            w.append(waai)
                    for wabi in wab:
                        if wabi not in w:
                            w.append(wabi)
                    for wbai in wba:
                        if wbai not in wb:
                            wb.append(wbai)
                    for wbbi in wbb:
                        if wbbi not in wb:
                            wb.append(wbbi)
                    MultipleCombinationi.append([w, wb])
            if MultipleCombinationi != []:
                MultipleCombination.append(MultipleCombinationi)
    MultipleControl = []
    MultipleControlLim = []
    for mm in MultipleCombination:
        MultipleControl.append(0)
        MultipleControlLim.append(len(mm))
    MultipleControlLim.append(1)
    MultipleControlLimMultiplication = []
    la = 1
    for li in range(0, len(MultipleControlLim) - 1):
        la *= MultipleControlLim[len(MultipleControlLim) - li - 1]
        MultipleControlLimMultiplication.append(la)
    SumType = 1
    for si in MultipleControlLim:
        SumType *= si
    for wi in range(0, SumType):
        for wii in range(0, len(MultipleControl)):
            MultipleControl[wii] = int(((wi - (wi % MultipleControlLimMultiplication[len(MultipleControl) - wii - 1]))/MultipleControlLimMultiplication[len(MultipleControl) - wii - 1])%MultipleControlLim[wii])
        apreresult = []
        for ma in range(0, len(MultipleControl)):
            apreresult.append(MultipleCombination[ma][MultipleControl[ma]])
        newapri=[fsatomic,fbound]
        newapr=copy.deepcopy(newapri)
        for aa in apreresult:
            if aa[0] != []:
                for aaa in aa[0]:
                    if aaa not in newapr[0]:
                        newapr[0].append(aaa)
            if aa[1] != []:
                for aaa in aa[1]:
                    if aaa not in newapr[1]:
                        newapr[1].append(aaa)
        aresult.append(newapr)
    return aresult

def WeightCalculation(structure,PiptideNode,PiptideBond):
    satomin=structure[0]
    bond=structure[1]
    m=0.0
    for si in satomin:
        m+=Piptideweight[PiptideNode[si-1]]
    for bi in bond:
        m+=BondTypeWeight[PiptideBond[bi][2]]
    return round(m,4)
#################################
#################################
P=PiptedStructure(testnode,testbond)
PeptideHydrolyze=SuperNode(P)
GetPeptideCombination(GetAllPreStucture(CreatAllFragment(PeptideHydrolyze.FragmentList)))
at=(PeptideHydrolyze.FragmentEndNodeList)
bt=(PeptideHydrolyze.a)
ct=(PeptideHydrolyze.FragmentMonomerNumberList)
print('There are '+str(TotalFragmentCount(CreatAllFragment(PeptideHydrolyze.FragmentList)))+' different structure.')
#################################
#################################


