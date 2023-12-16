"""model.py"""
from flask_wtf import FlaskForm
from wtforms import SelectField, TextAreaField, SubmitField
from wtforms.validators import DataRequired
from sqlalchemy import Column, String, Integer

# from sqlalchemy.orm import declarative_base
from webapp import db


class InputGenes(FlaskForm):
    """Input genes"""

    input_species = SelectField(
        choices=[("human"), ("mouse"), ("rat")],
        # default="human",
        validators=[DataRequired()],
    )
    input_genes = TextAreaField("Input genes")
    submit = SubmitField("Submit")


# Base = declarative_base()


class GOBP(db.Model):
    """GO_Biological_Process_2023"""

    __tablename__ = "GOBP"
    id = Column(Integer, primary_key=True)
    Term = Column(String)
    Gene = Column(String)


class GOCC(db.Model):
    """GO_Cellular_Component_2023"""

    __tablename__ = "GOCC"
    id = Column(Integer, primary_key=True)
    Term = Column(String)
    Gene = Column(String)


class GOMF(db.Model):
    """GO_Molecular_Function_2023"""

    __tablename__ = "GOMF"
    id = Column(Integer, primary_key=True)
    Term = Column(String)
    Gene = Column(String)
