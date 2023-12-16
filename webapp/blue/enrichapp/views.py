"""views.py"""
import re
from math import log10, log2

from flask import Blueprint, render_template, request, session, redirect, url_for

import pandas as pd
from flatsplode import flatsplode

import plotly.graph_objects as go
from plotly.offline import plot

from biothings_client import get_client
import gseapy as gp

from webapp import db
from webapp.blue.enrichapp.models import InputGenes, GOBP, GOCC, GOMF

enrich_app = Blueprint("enrich_app", __name__, template_folder="templates")


@enrich_app.route("/", methods=["GET", "POST"])
def settings():
    """Settings page"""

    form = InputGenes()

    mg = get_client("gene")

    if "input_species" in session:
        form.input_species.data = session["input_species"]
        species = form.input_species.data
    else:
        species = "human"

    if "input_genes" in session:
        form.input_genes.data = session["input_genes"]

    # if "gene_table_fix" in session:
    #     gene_table_fix = pd.read_json(session["gene_table_fix"])
    # else:
    #     gene_table_fix = pd.DataFrame()

    # if "human_gene_table_fix" in session:
    #     human_gene_table_fix = pd.read_json(session["human_gene_table_fix"])
    # else:
    #     human_gene_table_fix = pd.DataFrame()

    gene_table_fix = pd.DataFrame()
    human_gene_table_fix = pd.DataFrame()

    error_message = None

    if (
        request.method == "POST"
        and form.validate_on_submit()
        and form.input_genes.data != ""
    ):
        gene_list = [
            i for i in re.split(r"\n|\r|\t|,|;| ", form.input_genes.data) if i != ""
        ]

        species = form.input_species.data

        try:
            mygene_table = mg.querymany(
                gene_list,
                scopes="symbol",
                fields="symbol,name,entrezgene,ensembl.gene,homologene.id",
                species=species,
                returnall=True,
            )

        except Exception as e:
            error_message = f"An error occurred: {e}"

        else:
            gene_table_list = []

            for gene_info in mygene_table["out"]:
                gene_table_list.append(
                    (
                        pd.DataFrame(list(flatsplode(gene_info))).filter(
                            [
                                "query",
                                "symbol",
                                "name",
                                "entrezgene",
                                "ensembl.gene",
                                "homologene.id",
                            ]
                        )
                    )
                )

            gene_table = pd.DataFrame()

            gene_table = pd.concat(
                [gene_table, pd.concat(gene_table_list, ignore_index=True)]
            )

            gene_table = (
                gene_table.rename(
                    columns={
                        "ensembl.gene": "ensemblgene",
                        "homologene.id": "homologene",
                    }
                )
                .assign(homologene=lambda df: df.homologene.astype("Int64"))
                .assign(homologene=lambda df: df.homologene.astype("object"))
            )

            gene_table_fix = (
                pd.DataFrame({"query": gene_list})
                .merge(
                    (
                        gene_table.filter(["query", "symbol", "name"])
                        .query("symbol.notna()")
                        .drop_duplicates()
                    ),
                    how="left",
                    on=["query"],
                )
                .merge(
                    (
                        gene_table.filter(["query", "symbol", "name", "entrezgene"])
                        .query("entrezgene.notna()")
                        .drop_duplicates()
                    ),
                    how="left",
                    on=["query", "symbol", "name"],
                )
                .merge(
                    (
                        gene_table.filter(["query", "symbol", "name", "ensemblgene"])
                        .query("ensemblgene.notna()")
                        .drop_duplicates()
                    ),
                    how="left",
                    on=["query", "symbol", "name"],
                )
                .merge(
                    (
                        gene_table.filter(["query", "symbol", "name", "homologene"])
                        .query("homologene.notna()")
                        .drop_duplicates()
                    ),
                    how="left",
                    on=["query", "symbol", "name"],
                )
                .drop_duplicates()
            )

            human_gene_table = (
                gene_table_fix.query("homologene.notna()")
                .filter(["query", "homologene"])
                .drop_duplicates()
            )

            gene_table_fix = gene_table_fix.drop(
                columns=["homologene"]
            ).drop_duplicates()

            if species != "human":
                try:
                    mygene_table = mg.querymany(
                        human_gene_table.homologene,
                        scopes="homologene.id",
                        fields="symbol,name,entrezgene,ensembl.gene",
                        species="human",
                        returnall=True,
                    )

                except Exception as e:
                    error_message = f"An error occurred: {e}"

                else:
                    gene_table_list = []

                    for gene_info in mygene_table["out"]:
                        gene_table_list.append(
                            (
                                pd.DataFrame(list(flatsplode(gene_info))).filter(
                                    [
                                        "query",
                                        "symbol",
                                        "name",
                                        "entrezgene",
                                        "ensembl.gene",
                                    ]
                                )
                            )
                        )

                    gene_table = pd.DataFrame()

                    gene_table = pd.concat(
                        [gene_table, pd.concat(gene_table_list, ignore_index=True)]
                    )

                    gene_table = (
                        gene_table.rename(
                            columns={
                                "query": "homologene",
                                "ensembl.gene": "ensemblgene",
                            }
                        )
                        .assign(homologene=lambda df: df.homologene.astype("Int64"))
                        .assign(homologene=lambda df: df.homologene.astype("object"))
                    )

                    human_gene_table_fix = (
                        pd.DataFrame({"query": gene_list})
                        .merge(
                            (
                                human_gene_table.query(
                                    "homologene.notna()"
                                ).drop_duplicates()
                            ),
                            how="left",
                            on="query",
                        )
                        .merge(
                            (
                                gene_table.filter(["homologene", "symbol", "name"])
                                .query("symbol.notna()")
                                .drop_duplicates()
                            ),
                            how="left",
                            on="homologene",
                        )
                        .merge(
                            (
                                gene_table.filter(
                                    ["homologene", "symbol", "name", "entrezgene"]
                                )
                                .query("entrezgene.notna()")
                                .drop_duplicates()
                            ),
                            how="left",
                            on=["homologene", "symbol", "name"],
                        )
                        .merge(
                            (
                                gene_table.filter(
                                    ["homologene", "symbol", "name", "ensemblgene"]
                                )
                                .query("ensemblgene.notna()")
                                .drop_duplicates()
                            ),
                            how="left",
                            on=["homologene", "symbol", "name"],
                        )
                        .drop_duplicates()
                    )

        gene_table_fix.index = range(1, len(gene_table_fix) + 1)
        gene_table_fix = gene_table_fix.fillna("")

        human_gene_table_fix.index = range(1, len(human_gene_table_fix) + 1)
        human_gene_table_fix = human_gene_table_fix.fillna("")

        if form.input_species.data == "human":
            gp_gene_list = (
                gene_table_fix.query("symbol.notna()")
                .drop_duplicates("symbol")
                .symbol.to_list()
            )
        else:
            gp_gene_list = (
                human_gene_table_fix.query("symbol.notna()")
                .drop_duplicates("symbol")
                .symbol.to_list()
            )

        session["input_species"] = species
        session["input_genes"] = form.input_genes.data
        # session["gene_table_fix"] = gene_table_fix.to_json()
        # session["human_gene_table_fix"] = human_gene_table_fix.to_json()
        session["gp_gene_list"] = gp_gene_list

    return render_template(
        "settings.html",
        gene_table=gene_table_fix.to_html(classes="display", header="true", index=True),
        human_gene_table=human_gene_table_fix.to_html(
            classes="display", header="true", index=True
        ),
        species=species,
        form=form,
        error_message=error_message,
    )


@enrich_app.route("/clear_session", methods=["POST"])
def clear_session():
    """Clear session"""
    session.clear()
    return redirect(url_for("settings"))


def replace_points(text):
    """Replace points in text"""
    spaces = [m.start() for m in re.finditer(" ", text)]

    if len(spaces) >= 8:
        return text[: spaces[7]] + "..."
    else:
        return text


@enrich_app.route("/GeneOntology", methods=["GET", "POST"])
def gene_ontology():
    """Gene Ontology page"""

    if "gp_gene_list" not in session:
        gp_df_1 = pd.DataFrame()
        fig_1 = go.Figure()
        gp_df_2 = pd.DataFrame()
        fig_2 = go.Figure()
        gp_df_3 = pd.DataFrame()
        fig_3 = go.Figure()
    else:
        gp_gene_list = session["gp_gene_list"]

        # GOBP
        gp_records = db.session.query(GOBP).all()
        db.session.close()

        geneset_dict = {}

        for record in gp_records:
            if record.Term in geneset_dict:
                geneset_dict[record.Term].append(record.Gene)
            else:
                geneset_dict[record.Term] = [record.Gene]

        enr = gp.enrichr(
            gene_list=gp_gene_list,
            gene_sets=geneset_dict,
            # gene_sets=["GO_Biological_Process_2023"],
            # organism="human",
            # outdir=None,
        )

        gp_results = (
            enr.results.filter(
                [
                    "Gene_set",
                    "Term",
                    "Overlap",
                    "P-value",
                    "Adjusted P-value",
                    "Odds Ratio",
                    "Combined Score",
                    "Genes",
                ]
            )
            .rename(
                columns={
                    "P-value": "p_val",
                    "Adjusted P-value": "adj_p",
                    "Odds Ratio": "OddsRatio",
                    "Combined Score": "Score",
                }
            )
            .sort_values("p_val")
            # .query(f'adj_p < {cutoff}')
            .query("p_val < 0.05")
            .assign(
                Term_ed=lambda df: df.Term.str.replace(" \\(GO:.*\\)$", "", regex=True)
            )
            .assign(mLog10_p_val=lambda df: df.p_val.apply(log10))
            .assign(mLog10_p_val=lambda df: df.mLog10_p_val * -1)
            .assign(mLog10_adj_p=lambda df: df.adj_p.apply(log10))
            .assign(mLog10_adj_p=lambda df: df.mLog10_adj_p * -1)
            .assign(num_Score=lambda df: df.Score)
            .assign(Log2_Score=lambda df: df.Score.apply(log2))
            .assign(Genes=lambda df: df.Genes.str.replace(";", "; "))
            .assign(p_val=lambda df: df.p_val.apply(lambda x: f"{x:.3e}"))
            .assign(adj_p=lambda df: df.adj_p.apply(lambda x: f"{x:.3e}"))
            .assign(OddsRatio=lambda df: df.OddsRatio.apply(lambda x: f"{x:.3f}"))
            .assign(Score=lambda df: df.Score.apply(lambda x: f"{x:.3f}"))
        )

        gp_df1 = gp_results.filter(
            ["Term", "Overlap", "p_val", "adj_p", "OddsRatio", "Score", "Genes"]
        )

        gp_df1.index = range(1, len(gp_df1) + 1)

        gp_df2 = gp_results.assign(
            Term_ed=lambda df: df.Term_ed.apply(replace_points)
        ).head(30)

        fig = go.Figure(
            go.Scatter(
                x=gp_df2.num_Score,
                y=gp_df2.Term_ed,
                mode="markers",
                orientation="h",
                marker=dict(
                    color=gp_df2.mLog10_p_val,
                    colorbar=dict(title="-log10<br>(p_val)", len=0.4),
                    size=gp_df2.Log2_Score * 2,
                    line=dict(color="DarkSlateGrey", width=1),
                ),
                customdata=gp_df2[["Term", "p_val", "adj_p", "Score"]],
                hovertemplate=(
                    "<b>Term:</b> %{customdata[0]}<br>"
                    + "<b>p_val:</b> %{customdata[1]}<br>"
                    + "<b>adj_p:</b> %{customdata[2]}<br>"
                    + "<b>Score:</b> %{customdata[3]}<br>"
                    + "<extra></extra>"
                ),
            )
        )

        fig.update_layout(
            autosize=True,
            height=600,
            # width=800,
            title="GO Biological Process 2023 (TOP 30)",
            xaxis_title="Score",
            yaxis_title="Term",
            margin=dict(l=10, r=10, t=40, b=10),
        )

        fig.update_yaxes(autorange="reversed")

        gp_df_1 = gp_df1
        fig_1 = fig

        # GOCC
        gp_records = db.session.query(GOCC).all()
        db.session.close()

        geneset_dict = {}

        for record in gp_records:
            if record.Term in geneset_dict:
                geneset_dict[record.Term].append(record.Gene)
            else:
                geneset_dict[record.Term] = [record.Gene]

        enr = gp.enrichr(
            gene_list=gp_gene_list,
            gene_sets=geneset_dict,
            # gene_sets=["GO_Biological_Process_2023"],
            # organism="human",
            # outdir=None,
        )

        gp_results = (
            enr.results.filter(
                [
                    "Gene_set",
                    "Term",
                    "Overlap",
                    "P-value",
                    "Adjusted P-value",
                    "Odds Ratio",
                    "Combined Score",
                    "Genes",
                ]
            )
            .rename(
                columns={
                    "P-value": "p_val",
                    "Adjusted P-value": "adj_p",
                    "Odds Ratio": "OddsRatio",
                    "Combined Score": "Score",
                }
            )
            .sort_values("p_val")
            # .query(f'adj_p < {cutoff}')
            .query("p_val < 0.05")
            .assign(
                Term_ed=lambda df: df.Term.str.replace(" \\(GO:.*\\)$", "", regex=True)
            )
            .assign(mLog10_p_val=lambda df: df.p_val.apply(log10))
            .assign(mLog10_p_val=lambda df: df.mLog10_p_val * -1)
            .assign(mLog10_adj_p=lambda df: df.adj_p.apply(log10))
            .assign(mLog10_adj_p=lambda df: df.mLog10_adj_p * -1)
            .assign(num_Score=lambda df: df.Score)
            .assign(Log2_Score=lambda df: df.Score.apply(log2))
            .assign(Genes=lambda df: df.Genes.str.replace(";", "; "))
            .assign(p_val=lambda df: df.p_val.apply(lambda x: f"{x:.3e}"))
            .assign(adj_p=lambda df: df.adj_p.apply(lambda x: f"{x:.3e}"))
            .assign(OddsRatio=lambda df: df.OddsRatio.apply(lambda x: f"{x:.3f}"))
            .assign(Score=lambda df: df.Score.apply(lambda x: f"{x:.3f}"))
        )

        gp_df1 = gp_results.filter(
            ["Term", "Overlap", "p_val", "adj_p", "OddsRatio", "Score", "Genes"]
        )

        gp_df1.index = range(1, len(gp_df1) + 1)

        gp_df2 = gp_results.assign(
            Term_ed=lambda df: df.Term_ed.apply(replace_points)
        ).head(30)

        fig = go.Figure(
            go.Scatter(
                x=gp_df2.num_Score,
                y=gp_df2.Term_ed,
                mode="markers",
                orientation="h",
                marker=dict(
                    color=gp_df2.mLog10_p_val,
                    colorbar=dict(title="-log10<br>(p_val)", len=0.4),
                    size=gp_df2.Log2_Score * 2,
                    line=dict(color="DarkSlateGrey", width=1),
                ),
                customdata=gp_df2[["Term", "p_val", "adj_p", "Score"]],
                hovertemplate=(
                    "<b>Term:</b> %{customdata[0]}<br>"
                    + "<b>p_val:</b> %{customdata[1]}<br>"
                    + "<b>adj_p:</b> %{customdata[2]}<br>"
                    + "<b>Score:</b> %{customdata[3]}<br>"
                    + "<extra></extra>"
                ),
            )
        )

        fig.update_layout(
            autosize=True,
            height=600,
            # width=800,
            title="GO Cellular Component 2023 (TOP 30)",
            xaxis_title="Score",
            yaxis_title="Term",
            margin=dict(l=10, r=10, t=40, b=10),
        )

        fig.update_yaxes(autorange="reversed")

        gp_df_2 = gp_df1
        fig_2 = fig

        # GOMF
        gp_records = db.session.query(GOMF).all()
        db.session.close()

        geneset_dict = {}

        for record in gp_records:
            if record.Term in geneset_dict:
                geneset_dict[record.Term].append(record.Gene)
            else:
                geneset_dict[record.Term] = [record.Gene]

        enr = gp.enrichr(
            gene_list=gp_gene_list,
            gene_sets=geneset_dict,
            # gene_sets=["GO_Biological_Process_2023"],
            # organism="human",
            # outdir=None,
        )

        gp_results = (
            enr.results.filter(
                [
                    "Gene_set",
                    "Term",
                    "Overlap",
                    "P-value",
                    "Adjusted P-value",
                    "Odds Ratio",
                    "Combined Score",
                    "Genes",
                ]
            )
            .rename(
                columns={
                    "P-value": "p_val",
                    "Adjusted P-value": "adj_p",
                    "Odds Ratio": "OddsRatio",
                    "Combined Score": "Score",
                }
            )
            .sort_values("p_val")
            # .query(f'adj_p < {cutoff}')
            .query("p_val < 0.05")
            .assign(
                Term_ed=lambda df: df.Term.str.replace(" \\(GO:.*\\)$", "", regex=True)
            )
            .assign(mLog10_p_val=lambda df: df.p_val.apply(log10))
            .assign(mLog10_p_val=lambda df: df.mLog10_p_val * -1)
            .assign(mLog10_adj_p=lambda df: df.adj_p.apply(log10))
            .assign(mLog10_adj_p=lambda df: df.mLog10_adj_p * -1)
            .assign(num_Score=lambda df: df.Score)
            .assign(Log2_Score=lambda df: df.Score.apply(log2))
            .assign(Genes=lambda df: df.Genes.str.replace(";", "; "))
            .assign(p_val=lambda df: df.p_val.apply(lambda x: f"{x:.3e}"))
            .assign(adj_p=lambda df: df.adj_p.apply(lambda x: f"{x:.3e}"))
            .assign(OddsRatio=lambda df: df.OddsRatio.apply(lambda x: f"{x:.3f}"))
            .assign(Score=lambda df: df.Score.apply(lambda x: f"{x:.3f}"))
        )

        gp_df1 = gp_results.filter(
            ["Term", "Overlap", "p_val", "adj_p", "OddsRatio", "Score", "Genes"]
        )

        gp_df1.index = range(1, len(gp_df1) + 1)

        gp_df2 = gp_results.assign(
            Term_ed=lambda df: df.Term_ed.apply(replace_points)
        ).head(30)

        fig = go.Figure(
            go.Scatter(
                x=gp_df2.num_Score,
                y=gp_df2.Term_ed,
                mode="markers",
                orientation="h",
                marker=dict(
                    color=gp_df2.mLog10_p_val,
                    colorbar=dict(title="-log10<br>(p_val)", len=0.4),
                    size=gp_df2.Log2_Score * 2,
                    line=dict(color="DarkSlateGrey", width=1),
                ),
                customdata=gp_df2[["Term", "p_val", "adj_p", "Score"]],
                hovertemplate=(
                    "<b>Term:</b> %{customdata[0]}<br>"
                    + "<b>p_val:</b> %{customdata[1]}<br>"
                    + "<b>adj_p:</b> %{customdata[2]}<br>"
                    + "<b>Score:</b> %{customdata[3]}<br>"
                    + "<extra></extra>"
                ),
            )
        )

        fig.update_layout(
            autosize=True,
            height=600,
            # width=800,
            title="GO Molecular Function 2023 (TOP 30)",
            xaxis_title="Score",
            yaxis_title="Term",
            margin=dict(l=10, r=10, t=40, b=10),
        )

        fig.update_yaxes(autorange="reversed")

        gp_df_3 = gp_df1
        fig_3 = fig

    return render_template(
        "gene_ontology.html",
        gp_df_1=gp_df_1.to_html(classes="display", header="true", index=True),
        gp_plot_1=plot(
            fig_1,
            output_type="div",
            include_plotlyjs=False,
            config={"responsive": True},
        ),
        gp_df_2=gp_df_2.to_html(classes="display", header="true", index=True),
        gp_plot_2=plot(
            fig_2,
            output_type="div",
            include_plotlyjs=False,
            config={"responsive": True},
        ),
        gp_df_3=gp_df_3.to_html(classes="display", header="true", index=True),
        gp_plot_3=plot(
            fig_3,
            output_type="div",
            include_plotlyjs=False,
            config={"responsive": True},
        ),
    )
