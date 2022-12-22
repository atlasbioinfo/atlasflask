from ast import arg
from copyreg import constructor
import math
import json,os,time
from multiprocessing import connection
from flask_restful import Resource,reqparse
from sqlalchemy.orm import load_only
from flask import jsonify
from g4atlasflask.models import ResearchDetails,Species,SpeStatistics,ATCGStats,DownloadTable
from g4atlasflask.models import HsQuickSearchGene,MmQuickSearchGene,SyQuickSearchGene,PpQuickSearchGene,EcQuickSearchGene,ScQuickSearchGene,PaQuickSearchGene,PfQuickSearchGene,AtQuickSearchGene,OsQuickSearchGene
from g4atlasflask.models import HsAvaliableData,MmAvaliableData,SyAvaliableData,PpAvaliableData,EcAvaliableData,ScAvaliableData,PaAvaliableData,PfAvaliableData,AtAvaliableData,OsAvaliableData
from g4atlasflask.models import HsGqsDetails,MmGqsDetails,SyGqsDetails,PpGqsDetails,EcGqsDetails,ScGqsDetails,PaGqsDetails,PfGqsDetails,AtGqsDetails,OsGqsDetails
from g4atlasflask.models import HsMrnaSeq,MmMrnaSeq,SyMrnaSeq,PpMrnaSeq,EcMrnaSeq,ScMrnaSeq,PaMrnaSeq,PfMrnaSeq,AtMrnaSeq,OsMrnaSeq
from g4atlasflask.models import HsRTCountCov,MmRTCountCov,SyRTCountCov,PpRTCountCov,EcRTCountCov,ScRTCountCov,PaRTCountCov,PfRTCountCov,AtRTCountCov,OsRTCountCov
from g4atlasflask.models import HsChemicalProbGini,MmChemicalProbGini,SyChemicalProbGini,PpChemicalProbGini,EcChemicalProbGini,ScChemicalProbGini,PaChemicalProbGini,PfChemicalProbGini,AtChemicalProbGini,OsChemicalProbGini
from g4atlasflask.models import IdentifiedRG4,IdentifiedRG4Dataset,TypeSpecificStatic

from g4atlasflask.models import HsRG4RSR,MmRG4RSR,SyRG4RSR,PpRG4RSR,EcRG4RSR,ScRG4RSR,PaRG4RSR,PfRG4RSR,AtRG4RSR,OsRG4RSR

gqsNameSeach={
    "GAv1AT":AtGqsDetails,
    "GAv1EC":EcGqsDetails,
    "GAv1HS":HsGqsDetails,
    "GAv1MM":MmGqsDetails,
    "GAv1PA":PaGqsDetails,
    "GAv1PF":PfGqsDetails,
    "GAv1PP":PpGqsDetails,
    "GAv1OS":OsGqsDetails,
    "GAv1SY":SyGqsDetails,
    "GAv1SC":ScGqsDetails
}

gqsNameSeachRes=[
    {"value":["GAv1AT","Arabidopsis thaliana","TAIR10"],"short":"Search for rG4 ID"},
    {"value":["GAv1EC","Escherichia coli","ASM584v2"],"short":"Search for rG4 ID"},
    {"value":["GAv1HS","Homo sapiens","GRCh38"],"short":"Search for rG4 ID"},
    {"value":["GAv1MM","Mus musculus","GRCm39"],"short":"Search for rG4 ID"},
    {"value":["GAv1PA","Pseudomonas aeruginosa","ASM676v1"],"short":"Search for rG4 ID"},
    {"value":["GAv1PF","Plasmodium falciparum","ASM276v2"],"short":"Search for rG4 ID"},
    {"value":["GAv1PP","Pseudomonas putida","ASM273612v1"],"short":"Search for rG4 ID"},
    {"value":["GAv1OS","Oryza sativa","IRGSP-1.0"],"short":"Search for rG4 ID"},
    {"value":["GAv1SY","Synechococcussp.WH8102","ASM19597v1"],"short":"Search for rG4 ID"},
    {"value":["GAv1SC","Saccharomy cescerevisiae","R64-1-1"],"short":"Search for rG4 ID"}
]

quickSearchGeneArr=[HsQuickSearchGene,MmQuickSearchGene,SyQuickSearchGene,PpQuickSearchGene,EcQuickSearchGene,ScQuickSearchGene,PaQuickSearchGene,PfQuickSearchGene,AtQuickSearchGene,OsQuickSearchGene]

avaliableDataArr=[HsAvaliableData,MmAvaliableData,SyAvaliableData,PpAvaliableData,EcAvaliableData,ScAvaliableData,PaAvaliableData,PfAvaliableData,AtAvaliableData,OsAvaliableData]

gqsDetailsArr=[HsGqsDetails,MmGqsDetails,SyGqsDetails,PpGqsDetails,EcGqsDetails,ScGqsDetails,PaGqsDetails,PfGqsDetails,AtGqsDetails,OsGqsDetails]

mrnaArr=[HsMrnaSeq,MmMrnaSeq,SyMrnaSeq,PpMrnaSeq,EcMrnaSeq,ScMrnaSeq,PaMrnaSeq,PfMrnaSeq,AtMrnaSeq,OsMrnaSeq]

rtArr=[HsRTCountCov,MmRTCountCov,SyRTCountCov,PpRTCountCov,EcRTCountCov,ScRTCountCov,PaRTCountCov,PfRTCountCov,AtRTCountCov,OsRTCountCov]

giniArr=[HsChemicalProbGini,MmChemicalProbGini,SyChemicalProbGini,PpChemicalProbGini,EcChemicalProbGini,ScChemicalProbGini,PaChemicalProbGini,PfChemicalProbGini,AtChemicalProbGini,OsChemicalProbGini]

rsrArr=[HsRG4RSR,MmRG4RSR,SyRG4RSR,PpRG4RSR,EcRG4RSR,ScRG4RSR,PaRG4RSR,PfRG4RSR,AtRG4RSR,OsRG4RSR]

class QuickSearchResource(Resource):
    #med=71
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        self.reqparse.add_argument('input', type=str, default="",
                                   location='json')
        super(QuickSearchResource,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 71):
            return "Identification error"
        args.input=args.input.strip()
        res=set()
        jsRes=[]
        if (args.input[0:4].upper()=="GAV1"):
            if (len(args.input)<=6):
                return jsonify({"res":gqsNameSeachRes})
            tname="GAv1"+args.input[4:6]
            if (tname not in gqsNameSeach):
                return jsonify({"res":gqsNameSeachRes})
            res=gqsNameSeach[tname].query.filter(gqsNameSeach[tname].gqsid.ilike('%'+args.input[6:]+'%')).limit(30).all()
            for i in res:
                value=[i.gqsid,i.mrnaid,i.beg,i.end,i.speInfo.spename]
                tdis=i.pqslabel
                jsRes.append({"value":value,"short":tdis})
        else:
            args.input=args.input.upper()
            for subt in quickSearchGeneArr:    
                res = set(subt.query.filter(subt.gene_name.ilike('%'+args.input+'%')).limit(10).all()).union(res)
                res = set(subt.query.filter(subt.geneid.ilike('%'+args.input+'%')).limit(10).all()).union(res)
            for i in res:
                value=[i.geneid]
                if (i.gene_name!=""):
                    value.append(i.gene_name)
                value.append(i.speInfo.spename)
                tdis=i.short_disc
                if (len(i.short_disc)>60):
                    tdis=i.short_disc[0:60]+"..."
                jsRes.append({"value":value,"short":tdis})
        return jsonify({"res":jsRes})

class SubSearchResource(Resource):
    #med=55
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        self.reqparse.add_argument('search', type=str, default="",
                                   location='json')
        super(SubSearchResource,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 55):
            return "Identification error"
        args.search=args.search.strip()
        # if (args.search[4:]=="GAv1"):
        #     print(args.search)         
        res=set()
        
        for subt in quickSearchGeneArr:    
            res = set(subt.query.filter(subt.gene_name.ilike('%'+args.search+'%')).limit(20).all()).union(res)
            res = set(subt.query.filter(subt.geneid.ilike('%'+args.search+'%')).limit(20).all()).union(res)
        jsRes=[]
        for i in res:
            tavaRes=avaliableDataArr[i.speindex-1].query.filter(avaliableDataArr[i.speindex-1].mrnaid == i.mrnaid).all()
            if (len(tavaRes)!=0):
                ava=[]
                hasAva=False
                for tava in tavaRes:
                    if (tava.seqfpkm>1):
                        hasAva=True
                        ava.append({
                        "seqDataName":tava.resdetails.exptype,
                        "treatment":tava.resdetails.treatment,
                        "fpkm":round(tava.seqfpkm,3)
                        })    
                if (hasAva):            
                    jsRes.append({"gene":i.geneid,"mrna":i.mrnaid,"name":i.gene_name,"gtype":i.gene_type,"short":i.short_disc,"spe":i.speInfo.spename,"PQS":i.gqscount,"Ava":ava,"hasAva":hasAva,"speindex":i.speindex})
                    continue
            jsRes.append({"gene":i.geneid,"mrna":i.mrnaid,"name":i.gene_name,"gtype":i.gene_type,"short":i.short_disc,"spe":i.speInfo.spename,"PQS":i.gqscount,"Ava":[],"hasAva":False,"speindex":i.speindex})
        return jsonify({"res":jsRes})

class GetSubName(Resource):
    #med=72
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        self.reqparse.add_argument('input', type=str, default="",
                                   location='json')
        super(GetSubName,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 72):
            return "Identification error"
        searchJSON=json.loads(args.input.replace("\'","\""))
        res = quickSearchGeneArr[int(searchJSON["speindex"])-1].query.filter(quickSearchGeneArr[int(searchJSON["speindex"])-1].geneid == searchJSON["gname"]).first()
        
        pqs=gqsDetailsArr[int(searchJSON["speindex"])-1].query.filter(gqsDetailsArr[int(searchJSON["speindex"])-1].mrnaid == res.mrnaid).all()

        ave=avaliableDataArr[int(searchJSON["speindex"])-1].query.filter(avaliableDataArr[int(searchJSON["speindex"])-1].mrnaid == res.mrnaid).all()

        mrna=mrnaArr[int(searchJSON["speindex"])-1].query.filter(mrnaArr[int(searchJSON["speindex"])-1].mrnaid == res.mrnaid).first()
        # gores=res.goinfo
        # # get go
        # goinfo=[]
        # for go in gores:
        #     goinfo.append({"goacc":go.goacc,"goname":go.goname,"godefine":go.godefine})
        # get pqs info
        pqsinfo=[]
        for tpqs in pqs:
            pqsinfo.append({
                "pqsid":tpqs.gqsid,
                "pqsbeg":tpqs.beg,
                "pqsend":tpqs.end,
                "pqslabel":tpqs.pqslabel,
                "g4hscore":tpqs.g4hscore
                })
        ava=[]
        for tava in ave:
            ava.append({
                "index":tava.researchdetailsid,
                "pubid":tava.resdetails.pubid,
                "spe":tava.resdetails.spe,
                "expType":tava.resdetails.exptype,
                "treatment":tava.resdetails.treatment,
                "citation":tava.resdetails.research.citation,
                "title":tava.resdetails.research.title,
                "resindex":tava.resdetails.index,
                "fpkm":round(tava.seqfpkm,3)
            })                
        return jsonify({
            "mrna":res.mrnaid,
            "name":res.gene_name,
            "gtype":res.gene_type,
            "short":res.short_disc,
            "spe":res.speInfo.spename,
            "geneID":res.geneid,
            # "goinfo":goinfo,
            "seq":mrna.seq,
            "pqs":pqsinfo,
            "ava":ava
            })

class GetPQSinfo(Resource):
    #med=959
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        self.reqparse.add_argument('input', type=str, default="",
                                   location='json')
        super(GetPQSinfo,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 959):
            return "Identification error"
        args.input=args.input.strip()
        res = QuickSearchGqs.query.filter(QuickSearchGqs.mrnaid == args.input).all()
        pqsinfo=[]
        for tpqs in res:
            pqsinfo.append({
                "pqsid":tpqs.gqsid,
                "pqsbeg":tpqs.gqsdetail.beg,
                "pqsend":tpqs.gqsdetail.end,
                "pqslabel":tpqs.gqsdetail.pqslabel,
                "g4hscore":tpqs.gqsdetail.g4hscore
                })
        return jsonify({
            "pqsinfo":pqsinfo
            })

class GetSpecificPQSinfo(Resource):
    #med=377
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        self.reqparse.add_argument('speindex', type=int, default="",
                                   location='json')
        self.reqparse.add_argument('gqsid', type=str, default="",
                                   location='json')
        super(GetSpecificPQSinfo,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 377):
            return "Identification error"
        gqsRes = gqsDetailsArr[args["speindex"]-1].query.filter(gqsDetailsArr[args["speindex"]-1].gqsid == args["gqsid"]).first()
        
        return jsonify({
            "pqsid":gqsRes.gqsid,
            "mrnaid":gqsRes.mrnaid,
            "pqsbeg":gqsRes.beg,
            "pqsend":gqsRes.end,
            "gqsflankseq":gqsRes.pqsseq,
            "pqslabel":gqsRes.pqslabel,
            "g4hscore":gqsRes.g4hscore,
            "beglen":gqsRes.beglen,
            "gqslen":gqsRes.gqslen,
            "endlen":gqsRes.endlen,
            "pqsstrudb":gqsRes.pqsstrudb,
            "pqsmfe":gqsRes.pqsmfe,
            "pqsbpp":gqsRes.pqsbpp,
            "pqsconsstrudb":gqsRes.pqsconsstrudb,
            "pqsconsmfe":gqsRes.pqsconsmfe,
            "pqsconsbpp":gqsRes.pqsconsbpp,
            "spename":gqsRes.speInfo.spename,
            "speindex":gqsRes.speindex
            })

class GetAraSeq(Resource):
    #med=295
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        self.reqparse.add_argument('input', type=str, default="",
                                   location='json')
        super(GetAraSeq,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 295):
            return "Identification error"
        args.input=args.input.strip()
        res = MrnaSeq.query.filter(MrnaSeq.mrnaid == args.input).first()
        return jsonify({
            "mrna":res.mrnaid,
            "seq":res.seq
            })

class GetAllSpecies(Resource):
    #med=366
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        super(GetAllSpecies,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 366):
            return "Identification error"
        spes = Species.query.all()
        
        speData=[{} for i in range(10)]
        for tspe in spes:
            tava_exp=set()
            for tres in tspe.researches:
                tava_exp.add(tres.exptype)
            speData[tspe.index-1]={
                "spename":tspe.spename,
                "speimage":tspe.speimage,
                "exps":list(tava_exp),
                "speversion":tspe.speversion,
                "spelink":tspe.spelink,
                "speindex":tspe.index
            }
        [speData[2],speData[8]]=[speData[8],speData[2]]
        [speData[3],speData[9]]=[speData[9],speData[3]]
        return jsonify({
            "speinfo":speData
            })

class GetAllResearchDetails(Resource):
    #med=708
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        super(GetAllResearchDetails,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 708):
            return "Identification error"
        ress = ResearchDetails.query.all()

        allres=[]
        for tsrr in ress:
            allres.append({
                "pubid":tsrr.pubid,
                "spe":tsrr.spe,
                "exptype":tsrr.exptype,
                "treatment":tsrr.treatment,
                "srrqc":"static/home/QC_reports/research"+str(tsrr.index)+"_report/multiqc_report.html",
                "citation":tsrr.research.citation,
                "srrs":tsrr.srrs,
                "resindex":tsrr.index
            })

        return jsonify({
            "allres":allres
            })

class GetExpData(Resource):
    #med=532
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        self.reqparse.add_argument('mrnaid', type=str, default="")
        self.reqparse.add_argument('pos', type=str, default="")
        self.reqparse.add_argument('input', type=str, default="",
                                   location='json')
        self.reqparse.add_argument('speindex', type=int, default="",
                                   location='json')
        self.reqparse.add_argument('gqsid', type=str, default="",
                                   location='json')
                                   
        super(GetExpData,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 532):
            return "Identification error"
        targetExp={}
        tpos=args.pos.split("-")
        flankBeg=int(tpos[0])-int(tpos[2])-1
        flankEnd=int(tpos[1])+int(tpos[3])
        if (flankEnd-flankBeg+1 > 161):
            return "Please do not do this. If you need more data please contact yiliang.ding@jic.ac.uk"
        for targs in args.input.split(","):
            ttmp=targs.split("_")
            if (ttmp[0] not in targetExp):
                targetExp[ttmp[0]]=set()
            targetExp[ttmp[0]].add(ttmp[1])
            
        expinfo=[]
        for texp in targetExp:
            res=rtArr[args["speindex"]-1].query.filter(rtArr[args["speindex"]-1].mrnaid == args.mrnaid, rtArr[args["speindex"]-1].researchdetailsid == int(texp)).first()
            if ("rt" in targetExp[texp] or "cov" in targetExp[texp]):
                if ("rt" in targetExp[texp]):
                    trt=",".join(res.rt.split(",")[flankBeg+1:flankEnd+1])
                    expinfo.append({
                        "resid":res.researchdetailsid,
                        "cite":" - ".join([res.resdetails.exptype, res.resdetails.treatment]),
                        "subcite": res.resdetails.research.citation,
                        "datatype":"rt",
                        "data":trt
                    })
                if ("cov" in targetExp[texp]):
                    tcov=",".join(res.cov.split(",")[flankBeg:flankEnd])
                    expinfo.append({
                        "resid":res.researchdetailsid,
                        "cite":" - ".join([res.resdetails.exptype, res.resdetails.treatment]),
                        "subcite": res.resdetails.research.citation,
                        "datatype":"cov",
                        "data":tcov
                    })
            if ("deltaRSR" in targetExp[texp] or "BnormTest" in targetExp[texp]):
                    rsrAve=rsrArr[int(args["speindex"])-1].query.filter(rsrArr[int(args["speindex"])-1].gqsid == args["gqsid"]).all()
                    for trsr in rsrAve:
                        # print(texp+"\t"+trsr.researchKId)
                        if (trsr.researchKId != int(texp)):
                            continue
                        if ("deltaRSR" in targetExp[texp]):
                            trsrK=trsr.rsrK.split(",")
                            trsrLi=trsr.rsrLi.split(",")
                            tdelta=[]
                            for i in range(len(trsrK)):
                                if (trsrK[i]!="N" and trsrLi[i]!="N"):
                                    tdelta.append(str(float(trsrK[i])-float(trsrLi[i])))
                                else:
                                    tdelta.append("N")
                            expinfo.append({
                                "resid":res.researchdetailsid,
                                "cite":" - ".join([res.resdetails.exptype, res.resdetails.treatment]),
                                "subcite": res.resdetails.research.citation,
                                "datatype":"deltaRSR",
                                "data":",".join(tdelta)
                            })
                        if ("BnormTest" in targetExp[texp]):
                            tpvalue=trsr.bnormP.split(",")
                            for i in range(len(tpvalue)):
                                if (tpvalue[i]!="N"):
                                    if (float(tpvalue[i])==0):
                                        tpvalue[i]=0.00001
                                    tpvalue[i]=str(-math.log10(float(tpvalue[i])))
                            
                            expinfo.append({
                                "resid":res.researchdetailsid,
                                "cite":" - ".join([res.resdetails.exptype, res.resdetails.treatment]),
                                "subcite": res.resdetails.research.citation,
                                "datatype":"BnormTest",
                                "data":",".join(tpvalue)
                            })

        return jsonify({
            "expinfo":expinfo
            })

class GetSpeInfo(Resource):
    #med=621
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        self.reqparse.add_argument('input', type=str, default="",
                                   location='json')
        super(GetSpeInfo,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 621):
            return "Identification error"
        speQuery=Species.query.filter(Species.index == int(args.input)).first()
        speStatQuery=SpeStatistics.query.filter(SpeStatistics.speindex == int(args.input)).all()
        speATCGStat=ATCGStats.query.filter(ATCGStats.speindex == int(args.input)).first()
        validatedDatasets=IdentifiedRG4Dataset.query.filter(IdentifiedRG4Dataset.speindex == int(args.input)).all()
        typeSpeATCGStat=TypeSpecificStatic.query.filter(TypeSpecificStatic.speindex == int(args.input)).all()
        # validatedDatasets=set(IdentifiedRG4.query.filter(IdentifiedRG4.speindex == int(args.input)).with_entities(IdentifiedRG4.pubid, IdentifiedRG4.exptype, IdentifiedRG4.treatment).all())
        valideRG4Data=[]
        for val in validatedDatasets:
            valideRG4Data.append({
                "index":val.index,
                "pubid":val.pubid,
                "speindex":val.speindex,
                "exptype":val.exptype,
                "treatment":val.treatment,
                "cite":val.research.citation
            })

        typeSpeArr=[]
        for tstat in typeSpeATCGStat:
            typeSpeArr.append({
                "biotype":tstat.biotype,
                "freqA":tstat.freqA,
                "freqT":tstat.freqT,
                "freqC":tstat.freqC,
                "freqG":tstat.freqG,
                "freqGQS":tstat.freqGQS
            })

        srrQueryInfo=[]
        for tsrr in speQuery.researches:
            if (tsrr.treatment=="vivo" or tsrr.pubid=="30413703"):
                srrQueryInfo.append({
                    "pubid":tsrr.pubid,
                    "spe":tsrr.spe,
                    "exptype":tsrr.exptype,
                    "treatment":tsrr.treatment,
                    "srrqc":"static/home/QC_reports/research"+str(tsrr.index)+"_report/multiqc_report.html",
                    "citation":tsrr.research.citation,
                    "vivo":1,
                    "srrs":tsrr.srrs,
                    "resindex":tsrr.index
                })
                continue
            srrQueryInfo.append({
                    "pubid":tsrr.pubid,
                    "spe":tsrr.spe,
                    "exptype":tsrr.exptype,
                    "treatment":tsrr.treatment,
                    "srrqc":"static/home/QC_reports/research"+str(tsrr.index)+"_report/multiqc_report.html",
                    "citation":tsrr.research.citation,
                    "vivo":0,
                    "srrs":tsrr.srrs,
                    "resindex":tsrr.index
                })

        speStatInfo=[]
        for tstat in speStatQuery:
            speStatInfo.append({
                "name":tstat.typename,
                "value":tstat.count
            })

        def checkZero(number):
            if (number==0):
                return "NA"
            return number

        speATCG=[
            ["Nuc","full","utr5","cds","utr3"],
            ["A",checkZero(speATCGStat.full_A),checkZero(speATCGStat.utr5_A),checkZero(speATCGStat.cds_A),checkZero(speATCGStat.utr3_A)],
            ["T",checkZero(speATCGStat.full_T),checkZero(speATCGStat.utr5_T),checkZero(speATCGStat.cds_T),checkZero(speATCGStat.utr3_T)],
            ["C",checkZero(speATCGStat.full_C),checkZero(speATCGStat.utr5_C),checkZero(speATCGStat.cds_C),checkZero(speATCGStat.utr3_C)],
            ["G",checkZero(speATCGStat.full_G),checkZero(speATCGStat.utr5_G),checkZero(speATCGStat.cds_G),checkZero(speATCGStat.utr3_G)],
        ]
        speGQSStat={
            "full_gqs":speATCGStat.full_gqs,
            "utr5_gqs":speATCGStat.utr5_gqs,
            "cds_gqs":speATCGStat.cds_gqs,
            "utr3_gqs":speATCGStat.utr3_gqs
        }
        return jsonify({
            "expdata":srrQueryInfo,
            "spename":speQuery.spename,
            "speversion":speQuery.speversion,
            "spelink":speQuery.spelink,
            "cdnastat":speStatInfo,
            "atcgstat":speATCG,
            "gqsstat":speGQSStat,
            "validatedRG4Data":valideRG4Data,
            "typeSpeArr":typeSpeArr
            })

class GetExperimentList(Resource):
    #med=1309
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        self.reqparse.add_argument('speindex', type=int, default="",
                                   location='json')
        self.reqparse.add_argument('mrnaid', type=str, default="",
                                   location='json')
        self.reqparse.add_argument('gqsid', type=str, default="",
                                   location='json')
        self.reqparse.add_argument('beg', type=int, default="",
                                   location='json')
        self.reqparse.add_argument('end', type=int, default="",
                                   location='json')

        super(GetExperimentList,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 1309):
            return "Identification error"
        ave=avaliableDataArr[int(args["speindex"])-1].query.filter(avaliableDataArr[int(args["speindex"])-1].mrnaid == args["mrnaid"]).all()
        ava=[]
        avaExpDataCount=0

        for tava in ave:
            if (tava.resdetails.pubid == "30413703"):
                continue
            if (round(tava.seqfpkm,3)>1):
                avaExpDataCount+=1
            ava.append({
                "index":tava.researchdetailsid,
                "pubid":tava.resdetails.pubid,
                "spe":tava.resdetails.spe,
                "expType":tava.resdetails.exptype,
                "treatment":tava.resdetails.treatment,
                "citation":tava.resdetails.research.citation,
                "title":tava.resdetails.research.title,
                "fpkm":round(tava.seqfpkm,3)
            })
        avaExpDataCount=avaExpDataCount*2

        tdata=IdentifiedRG4Dataset.query.filter(IdentifiedRG4Dataset.speindex== int(args["speindex"])).all()
        texptype=set()
        for tt in tdata:
            texptype.add(tt.exptype)

        rsrRes=[]
        if ("rG4-seq" in texptype):
            rsrAve=rsrArr[int(args["speindex"])-1].query.filter(rsrArr[int(args["speindex"])-1].gqsid == args["gqsid"]).all()
            for trsr in rsrAve:
                avaExpDataCount+=1
                tresK=ResearchDetails.query.filter(ResearchDetails.index == int(trsr.researchKId)).first()
                tresLi=ResearchDetails.query.filter(ResearchDetails.index == int(trsr.researchLiId)).first()

                rsrRes.append({
                    "researchKId":trsr.researchKId,
                    "researchLiId":trsr.researchLiId,
                    "rsrK_cite":tresK.research.citation,
                    "rsrK_treat":tresK.treatment,
                    "rsrLi_treat":tresLi.treatment,
                })
        
        return jsonify({
            "ava_exp_info":ava,
            "avaExpDataCount":avaExpDataCount,
            "rsrinfo":rsrRes
            })

class GetDownloadTable(Resource):
    #med=2731
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        super(GetDownloadTable,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 2731):
            return "Identification error"
        res=DownloadTable.query.all()
        predictRG4=[]
        rg4fold=[]
        giniinfo=[]
        annotation=[]
        for r in res:
            if (r.datatype=="PredictedRG4s"):
                predictRG4.append({
                    "spename":r.speinfo.spename,
                    "spelink":r.speinfo.spelink,
                    "datatype":r.datatype,
                    "datasource":r.datasource,
                    "datasourcelink":r.datasourcelink,
                    "note":r.note,
                    "datalink":r.datalink
                })
                continue
            if (r.datatype=="RG4sRNAStructure"):
                rg4fold.append({
                    "spename":r.speinfo.spename,
                    "spelink":r.speinfo.spelink,
                    "datatype":r.datatype,
                    "datasource":r.datasource,
                    "datasourcelink":r.datasourcelink,
                    "note":r.note,
                    "datalink":r.datalink
                })
                continue
            if (r.datatype=="Gini index of RT-stop"):
                giniinfo.append({
                    "spename":r.speinfo.spename,
                    "spelink":r.speinfo.spelink,
                    "datatype":r.datatype,
                    "datasource":r.datasource,
                    "datasourcelink":r.datasourcelink,
                    "note":r.note,
                    "datalink":r.datalink
                })
                continue
            annotation.append({
                "spename":r.speinfo.spename,
                "spelink":r.speinfo.spelink,
                "datatype":r.datatype,
                "datasource":r.datasource,
                "datasourcelink":r.datasourcelink,
                "note":r.note,
                "datalink":r.datalink
            })

        return jsonify({
            "predictRG4":predictRG4,
            "rg4fold":rg4fold,
            "giniinfo":giniinfo,
            "annotation":annotation
        })

class GetSpeciesExpList(Resource):
    #med=187
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        self.reqparse.add_argument('speindex', type=int, default="",
                                   location='json')
        super(GetSpeciesExpList,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 187):
            return "Identification error"
        
        ave=ResearchDetails.query.filter(ResearchDetails.speindex == args.speindex).all()
        ava=[] 
        for tava in ave:
            if (tava.exptype=="rG4-seq"):
                continue
            ava.append({
                "index":tava.index,
                "pubid":tava.pubid,
                "spe":tava.spe,
                "expType":tava.exptype,
                "treatment":tava.treatment,
                "citation":tava.research.citation,
                "title":tava.research.title
            }) 
        
        return jsonify({
            "ava_exp_info":ava
            })

class GetGiniIndex(Resource):
    #med=152
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        self.reqparse.add_argument('speindex', type=int, default="",
                                   location='json')
        self.reqparse.add_argument('dataindex', type=str, default="",
                                   location='json')
        self.reqparse.add_argument('overlap', type=float, default="",
                                   location='json')
        
        super(GetGiniIndex,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 152):
            return "Identification error"
        researchArr=args["dataindex"].split(",")
        res=set()
        identiRG4Info=[]
        for id in researchArr:
            tdata=IdentifiedRG4Dataset.query.filter(IdentifiedRG4Dataset.index== int(id)).first()
            res=set(IdentifiedRG4.query.filter(IdentifiedRG4.speindex == int(args["speindex"]),IdentifiedRG4.exptype == tdata.exptype,IdentifiedRG4.treatment == tdata.treatment).limit(100).all()).union(res)
        
        for r in res:
            if (int(args["overlap"])):
                if (r.eorn=="N"):
                    continue
            tlabel="N"
            if (r.eorn == "E"):
                tlabel="E"

            identiRG4Info.append({
                "exptype":r.exptype,
                "treatment":r.treatment,
                "gqsid":r.gqsid,
                "length":int(r.end)-int(r.beg)+1,
                "gene":r.gene,
                "eorn":tlabel,
            })
        
        return jsonify({
            "identifiedRG4List":identiRG4Info
            })

class GetContactMessage(Resource):
    #med=78
    def __init__(self):
        self.reqparse = reqparse.RequestParser()
        self.reqparse.add_argument('med',type=int,required=True,location='json')
        self.reqparse.add_argument('input', type=str, default="",
                                   location='json')
        super(GetContactMessage,self).__init__()

    def post(self):
        args = self.reqparse.parse_args(strict=True)
        if (args.med != 78):
            return "Identification error"
        with open("msg.tsv","a") as out:
            out.write(args.input+"\n")

        return "Cheers"
