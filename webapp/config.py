"""config.py"""

from datetime import timedelta

PERMANENT_SESSION_LIFETIME = timedelta(days=1)

SECRET_KEY = "mysecretkey"
SQLALCHEMY_DATABASE_URI = "sqlite:///GeneSet.db"
