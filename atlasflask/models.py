from atlasflask import rdb

class Species(rdb.Model):
    index=rdb.Column(rdb.Integer, primary_key=True, index=True)
    spename=rdb.Column(rdb.String(20))
    speversion=rdb.Column(rdb.Text)
    spelink=rdb.Column(rdb.Text)
    speimage=rdb.Column(rdb.Text)
    researches=rdb.relationship("ResearchDetails",backref="species")
