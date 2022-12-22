from flask import Flask,render_template,Blueprint,request
from atlasflask.blueprints.home import home
from atlasflask.settings import config
from atlasflask.extensions import rdb
from flask_restful import Api
# from atlasflask.g4atlas_api import QuickSearchResource
from flask_cors import CORS
from atlasflask.command import register_commands

import os

def create_app(config_name=None):
    if config_name is None:
        config_name = os.getenv('FLASK_CONFIG', 'development')
    
    app=Flask('atlasflask')
    app.config.from_object(config[config_name])
    CORS(app,resources={r'/aq/*':{"origins":"*"}})
    register_blueprints(app)
    register_extensions(app)
    register_commands(app)
    register_api(app)

    return app

def register_api(app):
    api = Api(app)


def register_errors(app):
    # @app.errorhandler(400)
    # def bad_request(e):
    #     return render_template('errors/400.html'), 400
    @app.errorhandler(404)
    def pagenotfound(e):
        return render_template('errors/404.html'), 404
    # @app.errorhandler(404)
    # def page_not_found(e):
    #     return render_template('errors/404.html'), 404

    # @app.errorhandler(500)
    # def internal_server_error(e):
    #     return render_template('errors/500.html'), 500


def register_blueprints(app):
    app.register_blueprint(home)
    
def register_extensions(app):
    rdb.init_app(app)


