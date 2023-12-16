"""__init__.py"""

from flask import Flask
from flask_sqlalchemy import SQLAlchemy


app = Flask(__name__)

app.config.from_object("webapp.config")

db = SQLAlchemy(app)


from webapp.blue.enrichapp.models import *
from webapp.blue.enrichapp.views import enrich_app

app.register_blueprint(enrich_app)

with app.app_context():
    db.create_all()
