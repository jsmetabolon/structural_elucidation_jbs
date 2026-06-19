# coding: utf-8
"""HTML report primitives for per-feature elucidation output.

Provides heading/box helpers, embedded PNG encoding, and FileReporter for
writing standalone HTML reports consumed by UserFeatureReporter and MS1 clustering.
"""

import base64
import datetime
import logging
import os

import matplotlib.pyplot as plt


def heading(s):
    """
    Generate an HTML heading: a red bar over blue text.
    :param s: text for header
    :return: HTML string.
    """
    return f'<hr height="10px" color="red"/><br/><h1 style="color:blue">{s}</h1><br>'


def build_tds(markups):
    """
    Creates a row of a table (a bunch of <td> items) from a list of cell contents.
    :param markups: List of stuff to go in the cells.
    :return: Long string of <td>s.
    """
    return "".join([f"<td>{markup}</td>\n" for markup in markups])


def today_string():
    """
    return a string representation of the first day of the rest of your life.
    :return: a string
    """
    today = datetime.date.today()
    today = "%4d-%02d-%02d" % (today.year, today.month, today.day)
    return today


def encode_png(png_file):
    """
    Encodes a png file in binary for dump into HTML file.
    :param png_file: Name of png file to encode.
    """
    return base64.b64encode(open(png_file, "rb").read()).decode("utf-8").replace("\n", "")


class NullReporter:
    def __init__(self, *argv, **kwargs):
        return

    def close_report(self, *argv, **kwargs):
        return

    def img(self, *argv, **kwargs):
        return

    def comment(self, *argv, **kwargs):
        return

    def start_accordion(self, *argv, **kwargs):
        return

    def start_small_accordion(self, *argv, **kwargs):
        return

    def end_accordion(self, *argv, **kwargs):
        return

    def htmls(self, *argv, **kwargs):
        return

    def flush(self, *argv, **kwargs):
        return


class FileReporter(object):
    """Accumulate HTML fragments and flush them to a report file on disk."""

    _text = ""
    _fp = None
    _title = None
    _file = None
    _reminder = ""

    def __init__(self, filename):
        """
        creates FileReporter object. Writes preamble to the report.
        :param filename: the name of the file to which the report is to be directed.
        """
        self._text = ""
        self._file = filename
        self._fp = open(filename, "w", encoding="utf-8")

        self._text = """<!DOCTYPE html>
        <head>
                <meta charset="UTF-8">
                <meta name="viewport" content="width=device-width, initial-scale=1">          
                <!-- BOOTSTRAP Latest compiled and minified CSS -->
                <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">
                <!-- Optional theme -->
                <!--link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css" integrity="sha384-rHyoN1iRsVXV4nD0JutlnGaslCJuC7uwjduW9SVrLvRYooPp2bWYgmgJQIXwl/Sp" crossorigin="anonymous"-->
                <!-- Latest compiled and minified JavaScript -->
                <!--script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script-->
                <style type="text/css"> div.wrapper {display: inline-block; page-break-inside: avoid;
                                                    font-size:11px;} </style>         
     <style>
    .accordion {
      background-color: #eee;
      color: olive;
      cursor: pointer;
      padding: 12px;
      width: 60%;
      border: 4px black
      text-align: left;
      outline: none;
      font-size: 16px;
      transition: 0.4s;
      position: relative;
    overflow: hidden; 
    }
    
    .active, .accordion:hover {
      background-color: #ccc; 
    }
    
    .panel {
      padding: 0 18px;
      display: none;
      background-color: white;
      overflow: hidden;
    }
    .accordion:after {
  content: '+';
  font-size: 13px;
  color: red;
  float: right;
  margin-left: 5px;
  position: absolute;
  transform: translateY(50%);
}

.active:after {
  content: '-';
  color: red;
}
.div {
    overflow-x:scroll;
}
    </style> 

              </head>
                        """
        return

    def close_report(self):
        """
        flushes the text for the report and closes the file.
        :return: file name of report.
        """
        self._text += """
<script>
var acc = document.getElementsByClassName("accordion");
var i;

for (i = 0; i < acc.length; i++) {
  acc[i].addEventListener("click", function() {
    this.classList.toggle("active");
    var panel = this.nextElementSibling;
    if (panel.style.display === "block") {
      panel.style.display = "none";
    } else {
      panel.style.display = "block";
    }
  });
}
</script> 
            </body>
            </html>
            """
        # print(self._text)
        self._fp.write(self._text)
        self._fp.close()
        self._text = ""
        return self._file

    @staticmethod
    def img(_=None):
        """
        Adds a plot in line in the report. Saves the plot to a file, encodes it, and puts the encoding into the report.
        :param _: A bit of a kludge, it allows us to write FileReporter.img(plt.plot(...)), for example,
        instead of plt.plot(); FileReporter.img()
        :return: an html element to add to the report.
        """
        import io
        import base64

        s = io.BytesIO()
        logging.debug("about to save figure")
        plt.savefig(s, format="png", bbox_inches="tight")
        logging.debug("saved figure")
        plt.close("all")
        s = base64.b64encode(s.getvalue()).decode("utf-8").replace("\n", "")
        return '<img align="left" src="data:image/png;base64,%s">' % s

    def comment(self, *texts):
        """
        Adds a simple text (or texts) to the report as a series of lines.
        :param texts:
        :return:  None; updates self._text.
        """
        text = " ".join([str(t) for t in texts])
        self._text += "<p>%s</p>\n" % text.replace("\n", "<br />\n")
        return

    def start_accordion(self, title):
        self._text += f"""<button class="accordion">{title}</button><div class="panel">"""

    def start_small_accordion(self, title):
        self._text += f"""<button style="width:51%;padding:12px;border4px blue" class="accordion">{title}</button><div class="panel">"""

    def end_accordion(self):
        self._text += "</div>"
        self.flush()  # what a great time to flush the file!

    def htmls(self, *markups):
        """
        Adds HTML elements to the report as a row in a table.
        :param markups: the items to be added.
        :return:  None; updates self._text.
        """
        html = """
            <div style="overflow-x:scroll;">
            <br/>
            <p align="left">
            <table>
            <tr>
            %s
            </tr>
            </table>
                        <br/>
            </p>
            <br/>
            </div>
                        <br/>
            """ % build_tds(
            markups
        )
        self._text += html
        return

    def flush(self):
        self._fp.write(self._text)
        self._text = ""
        self._fp.flush()
        os.fsync(self._fp.fileno())


def wrapit(text, tag="p", color="lightgray"):
    return f"""<table style="border: 3px; border-style:solid; border-color:{color}; padding: 15px;">
        <tr style="padding:15px"><td><{tag}>{text}</{tag}></td></tr></table>"""


def boxit(text, tag="p", color="lightgray"):
    """
    Puts a text in a little box with a gray (default) border.
    :param text: text for box
    :param tag: "p" is the default; not sure what else might work.
    :param color: border color
    :return:  an HTML string for the box.
    """
    return f'<{tag} style="border: 3px; border-style:solid; border-color:{color}; padding: 1em;">{text}</{tag}>'


def vboxit(text, tag="p", color="lightgray"):
    """
    :param text: the text for the box.
    :param tag: "p" is the default; not sure what else might work.
    :param color: border color.
    :return: HTML string for the box.
    """
    return f'<{tag} style="writing-mode:vertical-rl;border: 3px; border-style:solid; border-color:{color}; padding: 0.5em;">{text}</{tag}>'
