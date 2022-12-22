import click,os
from atlasflask import rdb
# from atlasflask.importData import importQuickSearchGene

def register_commands(app):
    @app.cli.command()
    @click.option('--drop', is_flag=True, help="Creat after drop")
    def initdb(drop):
        if drop:
            rdb.drop_all()
            click.echo("Dropdb db")  
        else:  
            rdb.create_all()
            click.echo("Initialized db")

    @app.cli.command()
    def add():
        
        print("added")
# importGOInfo('../../G4AtlasData/GOInfomation/Arago_Biomart_fix_final.tsv',False,True)