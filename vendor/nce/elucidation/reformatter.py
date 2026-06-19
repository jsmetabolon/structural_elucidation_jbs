"""Styled HTML and Excel formatting for elucidation summary DataFrames.

Applies column-specific renderers (formulas, SMILES links, spectra plots) and
exposes reformat as the main entry used by UserFeatureReporter and elucidator.
"""

import base64
import io
import logging
from collections import defaultdict
from urllib.parse import quote as urlquote

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import regex as re
from elucidation.htmlUtils import wrapit
from elucidation.postgresSmilesDb import PostgresSmilesDb
from elucidation.rdkitUtils import draw_smiles_image_for_html, canonicalize_smiles as canonicalize
from elucidation.conveniences import qw, lmap, is_empty_df


def fancy_fmla(s):
    """
    This takes a formula string and turns it into HTML with subscript numbers.
    :param s: plain formula string: CH4
    :return: HTML formula string: CH<sub>4</sub>
    """
    s = str(s)
    if s == "nan":
        return ""
    if s == "" or not s[0] in "QWERTYUIOPASDFGHJKLZXCVBNM":
        return s
    s = re.sub(r"([A-Z][a-z]?[']?)(\d+)", r"\g<1><sub>\g<2></sub>", s)
    s = re.sub(r"[']", "&#8217;", s)
    if s.endswith("(-1)") or s.endswith("(1-)"):
        s = s[:-4] + "<sup>-</sup>"
    if s.endswith("(+1)") or s.endswith("(1+)"):
        s = s[:-4] + "<sup>+</sup>"
    if s.endswith("(-2)") or s.endswith("(2-)"):  # who knows what the future may bring.
        s = s[:-4] + "<sup>2-</sup>"
    if s.endswith("(+2)") or s.endswith("(2+)"):
        s = s[:-4] + "<sup>2+</sup>"
    return s


def fancy_column(s):
    s = fancy_fmla(s)
    if len(s) > 0 and s[0] in "QWERTYUIOPASDFGHJKLZXCVBNM":
        return f"<p>&nbsp;{s}&nbsp;</p>"  # Formulas are horizontal
    else:
        return f'<p style="writing-mode:vertical-rl;">{s}</p>'


def hover(hover_color="#ffff99"):
    """
    Causes a row in a table to be highlghted when the cursor is over it.
    :param hover_color: highlighting color, default is a yellow.
    :return: A dictionary included in styles for dataframe formatter.
    """
    return dict(selector="tr:hover", props=[("background-color", "%s" % hover_color)])


def highlight_values(data, vals=None, color="lightgreen"):
    """
    Will highlight values in a dataframe column.
    :param data: dataframe to be formatted.
    :param vals: values that will be highlighted
    :param color: highlight color.
    :return:
    """
    if vals is None:
        vals = {}
    attr = "background-color: {}".format(color)
    is_val = data.isin(vals)
    if data.ndim == 1:  # Series from .apply(axis=0) or axis=1
        return [attr if v else "" for v in is_val]
    else:  # from .apply(axis=None)
        return pd.DataFrame(np.where(is_val, attr, ""), index=data.index, columns=data.columns)


def dotf(k):
    """
    returns a function that can be used to format reals to the speicified number of decimal places using f format.
    :param k: number of decimal places.
    :return: a function of one argument (float to be formatted)
    """
    s = "{:,." + str(k) + "f}"
    return lambda x: s.format(x) if not np.isnan(x) else ""


def dote(k):
    """
    returns a function that can be used to format reals to the speicified number of decimal places using e format.
    :param k: number of decimal places.
    :return: a function of one argument (float to be formatted)
    """
    s = "{:,." + str(k) + "e}"
    return lambda x: s.format(x) if not np.isnan(x) else ""


# f'<a href="https://cactus.nci.nih.gov/chemical/structure/{ss}/image" target="_blank">{s}</a>'
# %http://zinc15.docking.org/substances/?ecfp4_fp=%20O=C(O)CNC(=O)c1ccc[nH]1%20
# urlquote('O=CCCC(=O)O')


def smiles_link(s):
    """
    Given a SMILES string, it creates an html box containing SMILES, structure image, link to ZINC database (for, e.g.,
    supplier information), and a name (though name is currently unimplemented).
    @param s: SMILES string
    @type s: str
    @return: html for the boxed table.
    @rtype: str
    """
    s = str(s)  # might be NaN
    # don't build structures for these little guys because they are obvious and take up too much space.
    if s in qw("*O *C *N *CC *CCC *CCCC *CCCCC *CCCCCC *I *Br *Cl"):
        return s
    if s == "nan":
        return ""
    try:
        pic = draw_smiles_image_for_html(s)
        img = '<img align="left" src="data:image/png;base64,%s">' % pic
    except Exception as e:
        img = ""
    zinc_smiles = s.replace("*", "[H]")
    try:
        zinc_smiles = canonicalize(zinc_smiles)
    except Exception as e:
        pass
    url = f"http://zinc15.docking.org/substances/?ecfp4_fp={urlquote(zinc_smiles)}"
    zinc = f'<a href="{url}" target="_blank"><b>ZINC</b></a>'
    name = get_a_smiles_name(zinc_smiles)
    lis = [t for t in re.split(r"(.{30})", s) if t != ""]
    t = "<br>".join(lis)
    return f"<table><tr><td>{t}</td></tr><tr><td>{img}</td></tr><tr><td>{zinc}</td></tr><tr><td>{name}</td></tr></table>"


def plot_matched_spectrum(spectrum_match, show=False):
    """
    Plots match of an observed spectrum to a SMILES string.
    There is a gray dot, red triangle, or green triangle for each ion in the MS2 spectrum.
    :param spectrum_match: a data structure describing the three groups.
    :return: html string for the image.
    """

    if not spectrum_match:
        return ""
    try:
        fig, ax = plt.subplots(nrows=1, ncols=1)
        fig.set_size_inches(4, 3)
        ax.set_facecolor("white")
        ax.grid(color="gray", linewidth=1)
        successes, failures, unmodeled = spectrum_match
        if unmodeled:  # don't unzip something empty!
            masses, ratios = zip(*unmodeled)
            # for m, i in unmodeled:
            #     plt.plot([m, m], [0, i], lw=0.5, color="gray")
            plt.plot(masses, ratios, ".", color="black", alpha=0.5)

        if failures:
            masses, ratios = zip(*failures)
            # for m, i in failures:
            #     plt.plot([m, m], [0, i], lw=0.5, color="red")
            plt.plot(masses, ratios, "v", color="deeppink")

        if successes:
            masses, ratios = zip(*successes)
            # for m, i in successes:
            #     plt.plot([m, m], [0, i], lw=0.5, color="green")
            plt.plot(masses, ratios, "^", color="lime")

        s = io.BytesIO()
        plt.savefig(s, format="png", bbox_inches="tight")
        plt.close("all")

        s = base64.b64encode(s.getvalue()).decode("utf-8").replace("\n", "")
        s = '<img align="left" src="data:image/png;base64,%s">' % s
        s = f"<table><tr><td>{s}</td></tr></table>"
        return s
    except Exception as e:
        print(f"Exception in plot_matched_spectrum: {e}")


def get_a_smiles_name(smiles):  # TODO make this try MLIMS first.
    return ""
    print("get_a_smiles_name")
    sm = PostgresSmilesDb.get_instance()  #  This doesn't depend on library.
    print(sm)
    name = sm.get_a_smiles_name(smiles)
    print(name)
    return name


def smiles_link_list(lis):  # TODO I think this might be dead.
    if type(lis) is list:
        links = lmap(smiles_link, lis)
        links = "\n".join([f"<tr><td>{link}</td></tr>" for link in links])
        return f"<table>\n{links}\n</table>"
    return ""


def fmt_percent(x):
    """
    function to format a percentage.
    :param x: float to be turned into percentage.
    :return: string with percent sign.
    """
    return "%.1f%%" % (100 * x)


def forgiving_int(x):
    """
    Formatter to format integers, but it's forgiving if argument is not an integer.
    :param x: an object
    :return: a string for the object.
    """
    if np.isnan(x):
        return ""
    try:
        return str(int(x))
    except Exception as e:
        return str(x)


def fancy_list(lis):
    """
    Given a list of plain formulas, return HTML for a stacked list of subscripted formulas.
    :param lis: a list of plain formulas
    :return: HTML element for stacked formulas.
    """

    return "<br>".join([fancy_fmla(f) for f in lis])


def build_format_dictionary():
    """
    Builds a dictionary holding formaters for different columns names.  The current scheme is that columns should have
    unique names even across all tables.  It's getting a little awkward.
    :return: a dictionary mapping column names to formattting functions.
    """

    format_dict = defaultdict(lambda: lambda x: x)

    # These columns should be in float format with 1 decimal place.
    format_dict.update(
        {c: dotf(1) for c in qw("C_predict S_predict Br_predict Cl_predict O_predict")}
    )

    format_dict.update({c: dotf(3) for c in qw("recall precision prob probability mass_frac")})

    format_dict.update({c: dotf(0) for c in qw("ri pred_ri")})

    # These columns have 4 decimal places.
    format_dict.update(
        {
            c: dotf(4)
            for c in qw(
                "mass target_mass ion_mass dev_m prec_m adduct_m diff_m parent_mass mz_feature mz_buddy inferred_m"
                + " neutral_mass neutral_m delta theor_delta mz m_z"
            )
        }
    )

    # These are in scientific form with 2 decimal places.
    format_dict.update(
        {
            c: dote(2)
            for c in qw(
                "area height med_int mean_int wtd_mean wtd_std wtd_stderr conf_lo conf_hi"
                + " ratio_mean ratio_std ratio_25% ratio_50% ratio_75% ratio_max ratio_min"
                + " log_prec_int log_adduct_int theor_ratio max_ratio med_ratio min_ratio"
                + " mean_ratio std_ratio wgt_mean_ratio wgt_med_ratio"
                + " prec_int adduct_int log_ratio slope"
            )
        }
    )

    # Floating point, 2 decimal places.
    format_dict.update(
        {
            c: dotf(2)
            for c in qw("start_rt apex_rt end_rt ppmsd_ratio std_errors  r2 pearsonr ppm c_o_v")
        }
    )

    # Floating point, 3 decimal places.
    format_dict.update(
        {
            c: dotf(3)
            for c in qw(
                "cv ratio cv_ratio wgt_cv norm_int norm_pct_dev rel_intens adj_prob recall score "
            )
            + qw("precision match_score old_score rd_score my_score")
        }
    )

    # These columns contains formulas that should have subscripting added.
    format_dict.update({c: fancy_fmla for c in qw("cand_fmla neutral _ fmla mol_formula formula")})

    # These have lists of things that should be joined by commas.
    format_dict.update(
        {
            c: lambda lis: ", ".join(lmap(str, lis))
            for c in qw("has_misc is_misc has_isotopes is_isotopes")
        }
    )

    # These are floats that should be displayed as percentages.
    format_dict.update({c: fmt_percent for c in qw("wgt_dev_med wgt_dev_theor %C_lo %C_hi ")})

    # These should be displayed as integers, but sometimes they are missing.
    format_dict.update(
        {
            c: forgiving_int
            for c in qw("dbe ms2 ri pred_ri rank_sum chemical_id start_scan apex_scan end_scan")
        }
    )

    # These contain SMILES string that should link to NIH images.
    format_dict.update(
        {c: smiles_link for c in qw("canonical_smiles subsmiles smiles zapped my_smiles rd_smiles")}
    )

    format_dict.update(
        {
            c: plot_matched_spectrum
            for c in qw("spectrum_match spectrum_match_1 my_spectrum_match rd_spectrum_match")
        }
    )

    format_dict.update({c: smiles_link_list for c in qw("smiles_db")})

    # These are lists of formulas that should be stacked and have subscripts added.
    format_dict.update({c: fancy_list for c in qw("fmlas adduct_fmlas")})
    return format_dict


default_styles = [
    hover(),
    dict(
        selector="th",
        props=[
            ("font-size", "14pt"),
            ("font-weight", "bold"),
            ("background-color", "lightgray"),
            ("font-family", "Calibri"),
            ("min-width", "10px"),
        ],
    ),
    dict(
        selector="td",
        props=[
            ("font-size", "12pt"),
            ("text-align", "left"),
            ("min-width", "10px"),
            ("border", "2px solid black"),
            ("padding", "5px"),
            ("font-family", "Calibri"),
        ],
    ),
    dict(selector="caption", props=[("caption-side", "bottom")]),
]

format_dict = build_format_dictionary()


def reformat(df_in, index=None, alt=None, gripe=False, styles=default_styles):

    """
    Turns dataframe into HTML form for inclusion in reports.
    :param df_in: The dataframe.
    :param index: The column to use as an index.  If None, then supress index.
    :param alt: Text to display if dataframe is empty.
    :param gripe: If True, complain about formatting problems.  Formatting problems generally don't "break" anything.
    :param styles: a list of dictionaries specifying styles of various elements.
    :return:
    """
    # Handle empty dataframe.
    if alt is not None and is_empty_df(df_in):
        return wrapit(alt, tag="h3")
    if is_empty_df(df_in):
        return None

    stack_cols = qw("db_id ri_obs ri_pred ri_X name source fmla")
    # This is for columns with lists that need to be stacked.
    df = df_in.copy()
    for c in df.columns:
        try:
            if c in format_dict:
                df[c] = df[c].apply(format_dict[c])
            elif set("CHONPSIB") & set(
                str(c)
            ):  # It's the ion consistency report. The column heading is a precursor formula.
                # but c might be an integer; can't make a set from it then.
                df[c] = df[c].apply(fancy_list)
        except Exception as e:
            if gripe:  # ignore this one for now.  Sometimes it's a string.
                logging.warning(str(e))
                logging.warning(f"Trouble reformatting column {c} and not reformatting.\nContents:")
                logging.warning(list(df[c]))
        if c in stack_cols and isinstance(df[c].iloc[0], list):
            df[c] = df[c].apply(lambda lis: "<br>".join(map(str, lis)))
    if len(df.columns) == 1:
        df.columns = [""]
    df.columns = lmap(fancy_column, df.columns)

    if False and type(index) is str:  # maybe later.
        df.index = df[index]
        df = df.drop(index, axis=1)
    styled = df.style
    # highlight YES, NO, and MAYBE
    styled = styled.apply(highlight_values, vals={"Y", "yes", "Yes", "YES"}, color="palegreen")
    styled = styled.apply(highlight_values, vals={"N", "no", "No", "NO"}, color="hotpink")
    styled = styled.apply(highlight_values, vals={"maybe", "Maybe", "MAYBE", "?"}, color="#ffff99")

    if index is not True:
        styled = styled.hide_index()
    styled = styled.set_table_styles(styles)
    return styled.render()


# THESE ARE USEFUL PATTERNS FOUND ON THE WEB.
#
# def highlight_too_big(threshold=1):
#     """
#     Returns a formatting function that highlights numbers greater than the threshold in pink.
#     Not in use, but keeping it as an example of this technique.
#     :param threshold: The threshold for highlighting
#     :return: A function that highlights.
#     """
#     def highlighter(data):
#         if data.ndim == 1:  # Series from .apply(axis=0) or axis=1
#             colors = ["hotpink", "white"]
#             print(data)
#             attrs = [f'"background-color: {colors[int(n>threshold)]};"' for n in data]
#             print(attrs)
#             return attrs
#         else:
#             print(data)
#             raise Exception("highlight_too_big will only highlight one column")
#
#     return highlighter
#
#
# def highlight_max(data, color="yellow"):  # example from web
#     """
#         highlight the maximum in a Series or DataFrame
#         """
#     attr = "background-color: {}".format(color)
#     if data.ndim == 1:  # Series from .apply(axis=0) or axis=1
#         is_max = data == data.max()
#         return [attr if v else "" for v in is_max]
#     else:  # from .apply(axis=None)
#         is_max = data == data.max().max()
#         return pd.DataFrame(np.where(is_max, attr, ""), index=data.index, columns=data.columns)
#
#
# df = pd.DataFrame([[1, 2], [1, 3], [4, 6]], columns=['A', 'B'])
# t = df.style.apply(highlight_max, color='darkorange', axis=0)
#
# myhtml = df.style.set_properties(**{'font-size': '11pt', 'font-family': 'Calibri', 'border-collapse': 'collapse',
#                                     'border': '1px solid black'}).render()
