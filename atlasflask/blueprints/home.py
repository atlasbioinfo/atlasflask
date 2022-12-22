from flask import render_template,abort,Blueprint,request,jsonify
from jinja2 import TemplateNotFound


home=Blueprint("home",__name__)
@home.route("/",methods=['GET','POST'])
def show():
    try:
        return render_template('home/index.html')
    except TemplateNotFound:
        abort(404)


