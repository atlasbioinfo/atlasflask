import os
import sys

class BaseConfig(object):
    SECRET_KEY = os.getenv('SECRET_KEY')
    DEBUG_TB_INTERCEPT_REDIRECTS = False
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    SQLALCHEMY_RECORD_QUERIES = True

class DevelopmentConfig(BaseConfig):
    SQLALCHEMY_DATABASE_URI = 'postgresql://hyu::bioinfo@localhost:5433/RNAstructure'

class TestingConfig(BaseConfig):
    TESTING = True
    WTF_CSRF_ENABLED = False
    SQLALCHEMY_DATABASE_URI = 'postgresql://hyu::bioinfo@localhost:5433/RNAstructure'


class ProductionConfig(BaseConfig):
    SQLALCHEMY_DATABASE_URI = 'postgresql://hyu::bioinfo@localhost:5433/RNAstructure'


config = {
    'development': DevelopmentConfig,
    'testing': TestingConfig,
    'production': ProductionConfig
}
