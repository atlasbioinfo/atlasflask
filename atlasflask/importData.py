import imp
import os
from g4atlasflask.models import Species,Researches,ResearchDetails,SpeStatistics,UtrInfo,ATCGStats,DownloadTable

from g4atlasflask.models import HsGqsDetails,MmGqsDetails,SyGqsDetails,PpGqsDetails,EcGqsDetails,ScGqsDetails,PaGqsDetails,PfGqsDetails,AtGqsDetails,OsGqsDetails,TypeSpecificStatic

from g4atlasflask.models import HsQuickSearchGene,MmQuickSearchGene,SyQuickSearchGene,PpQuickSearchGene,EcQuickSearchGene,ScQuickSearchGene,PaQuickSearchGene,PfQuickSearchGene,AtQuickSearchGene,OsQuickSearchGene

from g4atlasflask.models import HsAvaliableData,MmAvaliableData,SyAvaliableData,PpAvaliableData,EcAvaliableData,ScAvaliableData,PaAvaliableData,PfAvaliableData,AtAvaliableData,OsAvaliableData

from g4atlasflask.models import HsMrnaSeq,MmMrnaSeq,SyMrnaSeq,PpMrnaSeq,EcMrnaSeq,ScMrnaSeq,PaMrnaSeq,PfMrnaSeq,AtMrnaSeq,OsMrnaSeq

from g4atlasflask.models import HsRTCountCov,MmRTCountCov,SyRTCountCov,PpRTCountCov,EcRTCountCov,ScRTCountCov,PaRTCountCov,PfRTCountCov,AtRTCountCov,OsRTCountCov

from g4atlasflask.models import HsChemicalProbGini,MmChemicalProbGini,SyChemicalProbGini,PpChemicalProbGini,EcChemicalProbGini,ScChemicalProbGini,PaChemicalProbGini,PfChemicalProbGini,AtChemicalProbGini,OsChemicalProbGini

from g4atlasflask.models import HsRG4RSR,PaRG4RSR,EcRG4RSR,PfRG4RSR,AtRG4RSR

from g4atlasflask.models import IdentifiedRG4,IdentifiedRG4Dataset,IdentifiedG4RPSeq

from g4atlasflask import rdb

quickSearchGeneArr=["NA",HsQuickSearchGene,MmQuickSearchGene,SyQuickSearchGene,PpQuickSearchGene,EcQuickSearchGene,ScQuickSearchGene,PaQuickSearchGene,PfQuickSearchGene,AtQuickSearchGene,OsQuickSearchGene]

gqsDetailsArr=[HsGqsDetails,MmGqsDetails,SyGqsDetails,PpGqsDetails,EcGqsDetails,ScGqsDetails,PaGqsDetails,PfGqsDetails,AtGqsDetails,OsGqsDetails]

avaliableDataArr=[HsAvaliableData,MmAvaliableData,SyAvaliableData,PpAvaliableData,EcAvaliableData,ScAvaliableData,PaAvaliableData,PfAvaliableData,AtAvaliableData,OsAvaliableData]

mrnaArr=[HsMrnaSeq,MmMrnaSeq,SyMrnaSeq,PpMrnaSeq,EcMrnaSeq,ScMrnaSeq,PaMrnaSeq,PfMrnaSeq,AtMrnaSeq,OsMrnaSeq]

rtArr=[HsRTCountCov,MmRTCountCov,SyRTCountCov,PpRTCountCov,EcRTCountCov,ScRTCountCov,PaRTCountCov,PfRTCountCov,AtRTCountCov,OsRTCountCov]

giniArr=[HsChemicalProbGini,MmChemicalProbGini,SyChemicalProbGini,PpChemicalProbGini,EcChemicalProbGini,ScChemicalProbGini,PaChemicalProbGini,PfChemicalProbGini,AtChemicalProbGini,OsChemicalProbGini]

rsrDict={"research_2_1.rsr":HsRG4RSR,"research_39_38.rsr":HsRG4RSR,"research_3_1.rsr":HsRG4RSR,"research_41_40.rsr":HsRG4RSR,"research_62_61.rsr":PaRG4RSR,"research_64_63.rsr":EcRG4RSR,"research_66_65.rsr":PfRG4RSR,"research_67_65.rsr":PfRG4RSR,"research_72_71.rsr":AtRG4RSR,"research_73_71.rsr":AtRG4RSR}

def importQuickSearchGene(tspe, file, drop = False,create = False):
    print("Updating QuichSearch for "+str(tspe))
    if (drop):
        quickSearchGeneArr[tspe].__table__.drop(rdb.engine)
    if (create):
        quickSearchGeneArr[tspe].__table__.create(rdb.engine)

    with open(file, 'r') as f:
        f.readline()
        for line in f:
            tmp = line.strip().split("\t")
            tsp = quickSearchGeneArr[tspe](
                geneid=tmp[0],
                mrnaid=tmp[1],
                gene_name=tmp[2],
                short_disc=tmp[3],
                gene_type=tmp[4],
                gqscount=int(tmp[6]),
                speindex=tspe
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success")

# def importGOInfo(file, drop = False,create = False):
#     print("Updating GO info")
#     if (drop):
#         GOInfomation.__table__.drop(rdb.engine)
#     if (create):
#         GOInfomation.__table__.create(rdb.engine)
#     count=0
#     with open(file, 'r') as f:
#         f.readline()
#         for line in f:
#             tmp = line.strip().split("ATLAS")
#             count+=1
#             tsp = GOInfomation(
#                 geneid=tmp[0],
#                 goacc=tmp[1],
#                 goname=tmp[2],
#                 godefine=tmp[3]
#             )
#             rdb.session.add(tsp)
#         rdb.session.commit()
#     print("Success! Updated "+str(count)+" records")

# def importQuickSearchGQS(file, drop = False,create = False):
#     print("Updating QuickSearchGQS "+file)
#     if (drop):
#         QuickSearchGqs.__table__.drop(rdb.engine)
#     if (create):
#         QuickSearchGqs.__table__.create(rdb.engine)

#     with open(file,"r") as f:
#         for line in f:
#             tmp=line.strip().split("\t")
#             tsp = QuickSearchGqs(
#                 gqsid=tmp[0],
#                 mrnaid=tmp[1]
#             )
#             rdb.session.add(tsp)
#         rdb.session.commit()
#     print("Success!")

def importPQSDetails(file, index, drop = False,create = False):
    print("Updating PQSDetails "+file)
    if (drop):
        gqsDetailsArr[index].__table__.drop(rdb.engine)
    if (create):
        gqsDetailsArr[index].__table__.create(rdb.engine)

    with open(file,"r") as f:
        for line in f:
            tmp=line.strip().split("\t")
            tsp = gqsDetailsArr[index](
                gqsid=tmp[0],
                mrnaid=tmp[1],
                beg=int(tmp[2]),
                end=int(tmp[3]),
                pqsseq=tmp[11],
                pqslabel=tmp[5],
                g4hscore=float(tmp[6]),
                beglen=int(tmp[8]),
                gqslen=int(tmp[9]),
                endlen=int(tmp[10]),
                pqsstrudb=tmp[12],
                pqsmfe=tmp[13],
                pqsbpp=tmp[14],
                pqsconsstrudb=tmp[15],
                pqsconsmfe=tmp[16],
                pqsconsbpp=tmp[17],
                speindex=index+1
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importResearch(file, drop = False,create = False):
    print("Updating Researches")
    if (drop):
        Researches.__table__.drop(rdb.engine)
    if (create):
        Researches.__table__.create(rdb.engine)

    with open(file,"r") as f:
        for line in f:
            tmp=line.strip().split("\t")
            tsp = Researches(
                pubid=tmp[0],
                title=tmp[1],
                citation=tmp[2]
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importResearchDetails(file, drop = False,create = False):
    print("Updating Researches Details")
    if (drop):
        ResearchDetails.__table__.drop(rdb.engine)
    if (create):
        ResearchDetails.__table__.create(rdb.engine)

    speDict={
        "H. sapiens":1,
        "M. musculus (mESC)":2,
        "Synechococcus (WH8102)":3,
        "P. putida":4,
        "E. coli":5,
        "S. cerevisiae ":6,
        "H. sapiens (Hela)":1,
        "H. sapiens (HEK293T)":1,
        "M. musculus":2,
        "P. aeruginosa":7,
        "P. falciparum":8,
        "A. thaliana":9,
        "O. sativa ssp. japonica":10
    }
    with open(file,"r") as f:
        for line in f:
            tmp=line.strip().split("\t")
            tsp = ResearchDetails(
                index=int(tmp[0]),
                pubid=tmp[1],
                speindex=int(speDict[tmp[2]]),
                spe=tmp[2],
                exptype=tmp[3],
                treatment=tmp[4],
                srrs=tmp[5]
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importAvaliabledata(index,file,researchidArr=[], drop = False,create = False):
    print("Updating Avaliable data "+file)
    if (drop):
        avaliableDataArr[index].__table__.drop(rdb.engine)
    if (create):
        avaliableDataArr[index].__table__.create(rdb.engine)
    with open(file,"r") as f:
        f.readline()
        for line in f:
            line=line.replace("\"","")
            tmp=line.strip().split("\t")
            for i in range(len(researchidArr)):
                tsp= avaliableDataArr[index](
                    mrnaid=tmp[0],
                    researchdetailsid=researchidArr[i],
                    seqfpkm=float(tmp[i+1])
                )
                rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importRTcov(index,file, researchid,drop = False,create = False):
    print("Updating RTcov "+file)
    if (drop):
        rtArr[index].__table__.drop(rdb.engine)
    if (create):
        rtArr[index].__table__.create(rdb.engine)
    with open(file,"r") as f:
        for line in f:
            if (line[0]=='*'):
                continue
            tmp=line.strip().split("\t")
            tsp = rtArr[index](
                mrnaid=tmp[0],
                researchdetailsid=researchid,
                rt=tmp[1],
                cov=tmp[2]
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importMRNASeq(index,file, drop = False,create = False):
    print("Updating mrnaseq "+file)
    if (drop):
        mrnaArr[index].__table__.drop(rdb.engine)
    if (create):
        mrnaArr[index].__table__.create(rdb.engine)
    with open(file,"r") as f:
        fcont=f.read().split(">")
        for cont in fcont[1:]:
            tmp=cont.split("\n")
            name=tmp[0].split(" ")[0]
            tsp = mrnaArr[index](
                mrnaid=name,
                seq="".join(tmp[1:])
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importSpecies(file, drop = False,create = False):
    print("Updating Species list "+file)
    if (drop):
        Species.__table__.drop(rdb.engine)
    if (create):
        Species.__table__.create(rdb.engine)
    with open(file,"r") as f:
        title=f.readline()
        for line in f:
            tmp=line.strip().split("\t")
            tsp = Species(
                index=int(tmp[0]),
                spename=tmp[1],
                speversion=tmp[2],
                spelink=tmp[3],
                speimage=tmp[4]
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importSpeStatistics(file, drop = False,create = False):
    print("Updating Species list "+file)
    if (drop):
        SpeStatistics.__table__.drop(rdb.engine)
    if (create):
        SpeStatistics.__table__.create(rdb.engine)

    with open(file,"r") as f:
        title=f.readline()
        for line in f:
            tmp=line.strip().split("\t")
            tsp = SpeStatistics(
                speindex=int(tmp[0]),
                typename=tmp[1],
                count=int(tmp[2])
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importUTRinfo(file, drop = False,create = False):
    print("utrinfo "+file)
    if (drop):
        UtrInfo.__table__.drop(rdb.engine)
    if (create):
        UtrInfo.__table__.create(rdb.engine)

    with open(file,"r") as f:
        for line in f:
            tmp=line.strip().split("\t")
            tsp = UtrInfo(
                mrnaid=tmp[0],
                utrfive=int(tmp[1]),
                utrcds=int(tmp[2]),
                utrthree=int(tmp[3])
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!") 

def importATCGStats(file, drop = False,create = False):
    print("atcgStats "+file)
    if (drop):
        ATCGStats.__table__.drop(rdb.engine)
    if (create):
        ATCGStats.__table__.create(rdb.engine)
    file2spe={
        "Paerugi_cdna.atcg_gqs": 6 ,
        "Syne_cdna.atcg_gqs":  2,
        "Pfalciparum_cdna.atcg_gqs": 7 ,
        "Yeast_cdna.atcg_gqs": 5 ,
        "Homo_sapiens_cdna.atcg_gqs": 0 ,
        "Oryza_sativa_cdna.atcg_gqs": 9 ,
        "Mus_musculus_cdna.atcg_gqs": 1 ,
        "Ecoli_cdna.atcg_gqs":  4,
        "Pputida_cdna.atcg_gqs":  3,
        "Arabidopsis_thaliana_cdna.atcg_gqs": 8 
    }
    with open(file,"r") as f:
        f.readline()
        for line in f:
            tmp=line.strip().split("\t")
            tsp = ATCGStats(
                speindex=file2spe[tmp[0]]+1,
                full_A=float(tmp[1]),
                full_T=float(tmp[2]),
                full_C=float(tmp[3]),
                full_G=float(tmp[4]),
                utr5_A=float(tmp[5]),
                utr5_T=float(tmp[6]),
                utr5_C=float(tmp[7]),
                utr5_G=float(tmp[8]),
                cds_A=float(tmp[9]),
                cds_T=float(tmp[10]),
                cds_C=float(tmp[11]),
                cds_G=float(tmp[12]),
                utr3_A=float(tmp[13]),
                utr3_T=float(tmp[14]),
                utr3_C=float(tmp[15]),
                utr3_G=float(tmp[16]),
                full_gqs=float(tmp[17])*1000,
                utr5_gqs=float(tmp[18])*1000,
                cds_gqs=float(tmp[19])*1000,
                utr3_gqs=float(tmp[20])*1000
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!") 

def importDownloadResource(file, drop = False,create = False):
    print("Updating download resource "+file)
    if (drop):
        DownloadTable.__table__.drop(rdb.engine)
    if (create):
        DownloadTable.__table__.create(rdb.engine)
    with open(file,"r") as f:
        f.readline()
        for line in f:
            tmp=line.strip().split("\t")
            tsp = DownloadTable(
                speindex=int(tmp[0]),
                datatype=tmp[1],
                datasource=tmp[2],
                datasourcelink=tmp[3],
                note=tmp[4],
                datalink=tmp[5]
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importGini(speindex,file,resid, drop = False,create = False):
    print("Updating download resource "+file)
    if (drop):
        giniArr[speindex].__table__.drop(rdb.engine)
    if (create):
        giniArr[speindex].__table__.create(rdb.engine)
    with open(file,"r") as f:
        f.readline()
        for line in f:
            tmp=line.strip().split("\t")
            if (float(tmp[2])==0 or tmp[1]=="NA"):
                continue
            tsp = giniArr[speindex](
                speindex=speindex+1,
                gqsid=tmp[0],
                coverageg=float(tmp[2]),
                giniindex=float(tmp[1]),
                researchdetailsid=resid
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importRSR(file, drop = False,create = False):
    print("Updating mrnaseq "+file)
    ids=file[:-4].split("_")
    if (drop):
        rsrDict[file].__table__.drop(rdb.engine)
    if (create):
        rsrDict[file].__table__.create(rdb.engine)
    with open("../../G4AtlasData/RG4ExpInfo/"+file,"r") as f:
        for line in f:
            tmp=line.strip().split("\t")
            tsp = rsrDict[file](
                gqsid=tmp[0],
                researchKId=int(ids[1]),
                researchLiId=int(ids[2]),
                rsrK=tmp[1],
                rsrLi=tmp[2],
                bnormP=tmp[3]
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importIdentifiedRG4(index,file, drop = False,create = False):
    print("Updating identified rG4s "+file)
    if (drop):
        IdentifiedRG4.__table__.drop(rdb.engine)
    if (create):
        IdentifiedRG4.__table__.create(rdb.engine)
    with open(file,"r") as f:
        f.readline()
        for line in f:
            tmp=line.strip().split("\t")
            tsp = IdentifiedRG4(
                pubid=tmp[0],
                speindex=tmp[1],
                exptype=tmp[2],
                treatment=tmp[3],
                gqsid=tmp[4],
                gene=tmp[5],
                beg=int(tmp[6]),
                end=tmp[7],
                eorn=tmp[8]
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importIdentifiedRG4Data(file, drop = False,create = False):
    print("Updating identified Data "+file)
    if (drop):
        IdentifiedRG4Dataset.__table__.drop(rdb.engine)
    if (create):
        IdentifiedRG4Dataset.__table__.create(rdb.engine)
    with open(file,"r") as f:
        f.readline()
        for line in f:
            tmp=line.strip().split("\t")
            tsp = IdentifiedRG4Dataset(
                pubid=tmp[0],
                speindex=int(tmp[1]),
                exptype=tmp[2],
                treatment=tmp[3]
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")

def importTypeSpecificSTAT(file,drop = False,create = False):
    print("Updating ATCG Stat "+file)
    if (drop):
        TypeSpecificStatic.__table__.drop(rdb.engine)
    if (create):
        TypeSpecificStatic.__table__.create(rdb.engine)
    with open(file,"r") as f:
        
        for line in f:
            tmp=line.strip().split("\t")
            tsp = TypeSpecificStatic(
                speindex=int(tmp[0]),
                biotype=tmp[1],
                freqA=float(tmp[2]),
                freqT=float(tmp[3]),
                freqC=float(tmp[4]),
                freqG=float(tmp[5]),
                freqGQS=float(tmp[6])
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")


def importIdentifiedLigandData(file,drop = False,create = False):
    print("Updating identified ligand "+file)
    if (drop):
        IdentifiedG4RPSeq.__table__.drop(rdb.engine)
    if (create):
        IdentifiedG4RPSeq.__table__.create(rdb.engine)
    with open(file,"r") as f:
        f.readline()
        for line in f:
            tmp=line.strip().split("\t")
            tsp = IdentifiedG4RPSeq(
                allInfo=";".join(tmp)
            )
            rdb.session.add(tsp)
        rdb.session.commit()
    print("Success!")