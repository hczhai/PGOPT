# -*- coding: utf-8 -*-

import os, re, numpy as np, time, copy
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Table, TableStyle
from reportlab.platypus.tableofcontents import TableOfContents
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch, cm
from reportlab.lib.styles import ParagraphStyle, StyleSheet1
from reportlab.platypus import PageBreak, KeepTogether
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.lib.colors import Color
from reportlab.lib.enums import TA_CENTER
from reportlab.pdfgen import canvas

from meta.projection import DrawCluster
from meta.energetics import DrawEnergetics, httoev, httok, multi_name, multi_color, mini_name
from meta.analysis import OptimizationData
from meta.ana_utils.base import time_span_str, time_span_short
from meta.plot import PGPlot
from meta.ana_utils.thermal import gas_phase_partition, gas_phase_partition_d, surface_partition
from meta.ana_utils.thermal import gas_phase_partition_dd, gas_phase_zpe, surface_zpe
from meta.ana_utils.thermal import ktoht
from meta.ana_utils.plot_rpath import DrawRGraph
from cluster.base import Cluster
from surface.base import ClusterAtSurface

colors_base = {"blue": [40, 130, 220], "black": [20, 20, 20]}
colors = {}
for k, v in colors_base.items():
    colors[k] = Color(*list(np.array(v) / 255.0))

def freq_str(freql, width):
    rstr = ""
    rrstr = ""
    freqt = []
    for i, f in enumerate(freql):
        ft = "%.0f" % f
        if ft == "-0": ft = "0"
        if len(freqt) == 0 or ft != "0":
            freqt.append(ft)
        elif freqt[-1] == "0" and ft == "0":
            freqt[-1] = 2
        elif isinstance(freqt[-1], int) and ft == "0":
            freqt[-1] += 1
        elif freqt[-1] != "0" and ft == "0":
            freqt.append("0")
    for i, f in enumerate(freqt):
        if isinstance(f, int):
            freqt[i] = "0^%d" % f
    for i, f in enumerate(freqt):
        tstrb = tstr = rstr + ("" if i == 0 else ", ") + f
        if i != len(freql) - 1: tstr += ", ..."
        w = stringWidth(tstr, fontName="Helvetica", fontSize=10)
        if w > width: break
        rstr = tstrb
        rrstr = tstr
    return rrstr

class PGTable(Table, object):
    def __init__(self, data, colWidths=None, dx=0.0, dy=0.0, **args):
        super(PGTable, self).__init__(data, colWidths, **args)
        self.dx = dx
        self.dy = dy

    def split(self, availWidth, availHeight):
        sp = super(PGTable, self).split(availWidth, availHeight)
        for s in sp:
            s.dx = self.dx
            s.dy = self.dy
        return sp

    def draw(self):
        self.canv.translate(self.dx, self.dy)
        super(PGTable, self).draw()

def creation_method_name(mn):
    mn = mn.upper()
    x = { "BLDA": "Bond Length Distribution Algorithm",
        "CK": "Coalescence Kick method",
        "RV": "Random Vector method",
        "BC": "Bounds Checking method" }
    return "%s (%s)" % (mn, x[mn] if mn in x else "Unknown")

def basis_short(ba):
    return ba.replace("(", "").replace(")", "").replace("def2-", "")

class NumberedCanvas(canvas.Canvas):
    def __init__(self, *args, **kwargs):
        canvas.Canvas.__init__(self, *args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._doc.pageCounter += 1
        self._startPage()

    def save(self):
        """add page info to each page (page x of y)"""
        self._doc.pageCounter = 1
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.draw_page_number(num_pages)
            canvas.Canvas.showPage(self)
        canvas.Canvas.save(self)

    def draw_page_number(self, page_count):
        self.setFont("Helvetica", 11)
        self.drawRightString(letter[0] - 0.5 * inch, 0.5 * inch,
            "%d / %d" % (self._pageNumber, page_count))

class TOCDocTemplate(SimpleDocTemplate):
    def afterFlowable(self, flowable):
        if hasattr(flowable, 'toc_key'):
            text = flowable.getPlainText()
            style = flowable.style.name
            if style == "heading1":
                # self.canv.bookmarkPage(flowable.toc_key)
                self.notify('TOCEntry', (0, text, self.page, flowable.toc_key))
                self.canv.addOutlineEntry(flowable.toc_cont, flowable.toc_key, 0, True)
            elif style == "heading2":
                # self.canv.bookmarkPage(flowable.toc_key)
                self.notify('TOCEntry', (1, text, self.page, flowable.toc_key))
                self.canv.addOutlineEntry(flowable.toc_cont, flowable.toc_key, 1, True)

class TwoColumnTOCTemplate(TOCDocTemplate):
    def build(self, flowables, canvasmaker=canvas.Canvas, column_gap=None):
        self._calc()
        if column_gap is None: column_gap = 1 * cm
        self.addPageTemplates(
            [ PageTemplate(frames=
                [ Frame( self.leftMargin, self.bottomMargin, self.width / 2,
                        self.height, id='left', rightPadding=column_gap / 2,
                        showBoundary=0),
                    Frame( self.leftMargin + self.width / 2, self.bottomMargin,
                        self.width / 2, self.height, id='right', 
                        leftPadding=column_gap / 2, showBoundary=0),
                ])])
        BaseDocTemplate.build(self, flowables, canvasmaker=canvasmaker)
    
    def handle_pageBegin(self):
        self._handle_pageBegin()

def init_fonts():
    pass

class SimpleDraw(object):
    def __init__(self, filename, finals):
        self.filename = filename
        self.finals = finals
        init_fonts()
        self.doc = TwoColumnTOCTemplate(self.filename, pagesize=letter,
            leftMargin=0.5*inch, rightMargin=0.5*inch, topMargin=0.6*inch, bottomMargin=0.5*inch)

    def build(self, surface_depth=3.0, ratio=1.0, plot_charge=False, rotmat=None, zdepth=1.0,
        title="", plot_force=False, force_factor=1.0, force_image=3, perspective=True):
        if len(self.finals) == 0 or self.finals[0].surfnum == 0:
            rt1, _ = Report.test_width(self.finals)
        else:
            rt1 = Report.test_width_surf(self.finals, surface_depth=surface_depth,
                                         plot_charge=plot_charge, rotmat=rotmat)
        rt1 *= ratio
        story = []
        pstyle = ParagraphStyle(name="body", fontName="Helvetica",
            fontSize=11, leading=15, spaceAfter=10, textColor=colors["black"])
        tstyle = ParagraphStyle(name="tit", fontName="Helvetica",
            fontSize=13, leading=15, spaceAfter=10, textColor=colors["black"])
        emin = self.finals[0].energy if len(self.finals) != 0 else 0.0
        if title != "":
            story.append(Paragraph(title, tstyle))
        for ic, c in enumerate(self.finals):
            print (ic)
            xg = DrawCluster(c, simple=False, ratio=rt1, surface_depth=surface_depth,
                plot_charge=plot_charge, rotmat=rotmat, zdepth=zdepth, plot_force=plot_force,
                force_factor=force_factor, force_image=force_image, perspective=perspective)
            if c.energy is not None:
                e, ex = c.energy, httoev * (c.energy - emin)
                story.append(Paragraph("[#%d] " % (ic + 1) + c.label + 
                    "<br /> E = %.4f a.u. (%.3f eV)" % (e, ex), pstyle))
            else:
                story.append(Paragraph("[#%d] " % (ic + 1) + c.label, pstyle))
            if plot_charge:
                ichg = 0.0
                for iat in range(c.surfnum, c.n):
                    ichg += c.props["charges"][iat]
                story.append(Paragraph("Cluster charges = %.3f" % ichg, pstyle))
            story.append(xg)
        self.doc.build(story, canvasmaker=NumberedCanvas)

class Report(object):

    def __init__(self, filename, ipff=None, input_file=None, energy_cutoff=[ 0.40, 1.00 ], 
        energy_cutoff_hard=None, no_freq=False, no_path=False, nebmins_ref=None):
        opt_data = OptimizationData(input_file, ipff, nebmins_ref=nebmins_ref) \
            if input_file is not None else None
        self.ipff = ipff
        self.no_freq = no_freq # do not include TS/MIN in summary
        self.filename = filename
        self.d = opt_data
        self.energy_cutoff = energy_cutoff
        self.energy_cutoff_hard = energy_cutoff_hard
        self.no_path = no_path
        # init fonts
        init_fonts()
        self.doc = TOCDocTemplate(self.filename, pagesize=letter,
            leftMargin=0.5*inch, rightMargin=0.5*inch, topMargin=0.6*inch, bottomMargin=0.5*inch)
        self.stylesheet = StyleSheet1()
        self.stylesheet.add(ParagraphStyle(name="title", fontName="Helvetica",
            fontSize=18, leading=22, spaceAfter=22, textColor=colors["blue"]))
        self.stylesheet.add(ParagraphStyle(name="heading1", parent=self.stylesheet["title"],
            fontSize=15, leading=18, spaceAfter=10, spaceBefore=18, textColor=colors["black"]))
        self.stylesheet.add(ParagraphStyle(name="heading2", parent=self.stylesheet["title"],
            fontSize=12, leading=15, spaceAfter=10, spaceBefore=10, textColor="#343434", fontName="Helvetica"))
        self.stylesheet.add(ParagraphStyle(name="body", fontName="Helvetica",
            fontSize=11, leading=15, spaceAfter=10, textColor=colors["black"]))
        self.stylesheet.add(ParagraphStyle(name="cellbody", parent=self.stylesheet["body"],
            fontSize=10))
        self.stylesheet.add(ParagraphStyle(name="cellcenter", parent=self.stylesheet["body"],
            fontSize=10, alignment=TA_CENTER))
        self.stylesheet.add(ParagraphStyle(name="celltitle", parent=self.stylesheet["cellbody"],
            fontName="Helvetica", alignment=TA_CENTER))
        self.stylesheet.add(ParagraphStyle(name="celltitleleft", parent=self.stylesheet["cellbody"],
            fontName="Helvetica"))
        self.stylesheet.add(ParagraphStyle(name="indbody", parent=self.stylesheet["body"],
            leftIndent=15))
        self.tablestyles = {}
        self.tablestyles["main"] = TableStyle([
            ('FONT',    (1, 0), (-1, -1), "Helvetica"),
            ('FONT',    (0, 0), (0, -1), "Helvetica"),
            ('ALIGN',    (0, 0), (0, -1), "RIGHT"),
            ('FONTSIZE',    (0, 0), (-1, -1), 10),
            ("BACKGROUND", (0, 0), (0, -1), "#F4F4F4"),
            ('INNERGRID',    (0, 0), (-1, -1), 0.8, "#E1E1E1"),
            ('LINEABOVE',    (0, 0), (-1, 0),1, "#A5A5A5"),
            ('LINEBELOW',    (0, -1), (-1, -1),1, "#A5A5A5"),
            ('VALIGN',    (0, 0), (-1, -1), "MIDDLE")
        ])
        self.tablestyles["mainx"] = copy.deepcopy(self.tablestyles["main"])
        self.tablestyles["mainx"]._cmds += [
            ('ALIGN', (0, 0), (-1, -1), "CENTER"),
            ('LINEBELOW',    (0, 0), (-1, 0), 1, "#C5C5C5"),
            ('FONT',    (0, 0), (-1, 0), "Helvetica"),
        ]
        self.tablestyles["mainxx"] = copy.deepcopy(self.tablestyles["mainx"])
        self.tablestyles["mainxx"]._cmds += [
            ('FONT',    (1, 0), (1, -1), "Helvetica"),
        ]
    
    @staticmethod
    def test_width(finals, clip=False, clip_circle=False):
        if len(finals) == 0: return 5.0, 2.7
        tests = finals[:10]
        rwid = (letter[0] / 2 - 0.6 * inch) * 0.52
        mwid = (letter[0] / 2 - 0.6 * inch) * 0.48
        rots = [ [90.0, 1.0, 0.0, 0.0], [90.0, 0.0, 1.0, 0.0], [45.0, 1.0, 1.0, 0.0] ]
        widm = np.array([DrawCluster(t, ratio=5.0).width for t in tests]).mean() + 20
        widr = np.array([sum([DrawCluster(t, simple=True, rotmat=r, clip=clip,
            clip_circle=clip_circle, ratio=2.7).width for r in rots])
            for t in tests]).mean() + 20
        rt1 = 5.0 * mwid / widm
        rt2 = 2.7 * rwid / widr
        return rt1, rt2
    
    @staticmethod
    def test_width_surf(finals, surface_depth=3.0, plot_charge=True, rotmat=None):
        if len(finals) == 0:
            return 8.0
        tests = finals[:10]
        awid = letter[0] - 1.2 * inch
        # rots = [ [0.0, 1.0, 0.0, 0.0], "best", [90.0, -1.0, 0.0, 0.0] ]
        if rotmat is None:
            rots = [ [0.0, 1.0, 0.0, 0.0], [90.0, -1.0, 0.0, 0.0] ]
        else:
            rots = [ rotmat ]
        wida = np.array([sum([DrawCluster(t, rotmat=r, 
            ratio=8.0, surface_depth=surface_depth, plot_charge=plot_charge).width
            for r in rots]) for t in tests]).mean() + 50
        if awid / wida > 0.4:
            awid = wida * 0.4
        return 8.0 * awid / wida
    
    def draw_cluster_at_surface(self, finals):
        if len(finals) == 0: return None
        cdr = []
        mener = finals[0].energy if len(finals) != 0 else 0.0
        efactor = 27.21139
        awid = letter[0] - 1.2 * inch
        qwid = stringWidth("Freqs: ", fontName="Helvetica", fontSize=10)
        fwid = awid * 2.0 / 3.0 - 12 - 5 - qwid
        # rots = [ [0.0, 1.0, 0.0, 0.0], "best", [90.0, -1.0, 0.0, 0.0] ]
        rots = [ [0.0, 1.0, 0.0, 0.0], [90.0, -1.0, 0.0, 0.0] ]
        rt = Report.test_width_surf(finals, self.surface_depth)
        for c in finals:
            xa, xc = [ DrawCluster(c, simple=False, ratio=rt, rotmat=r, 
                surface_depth=self.surface_depth) for r in rots ]
            cdr.append(Table([[xa, xc]], (xa.width + 12, xc.width), 
                style=TableStyle([ ('LEFTPADDING', (0, 0), (-1, -1), 0)]), hAlign='LEFT'))
        data = []
        style = copy.deepcopy(self.tablestyles["main"])
        for i in range(len(cdr)):
            c = finals[i]
            pyname = "<a name=\"X%s\" />[<font face=\"Helvetica\">%d</font>] "    % (c.aname, i + 1)
            if c.similar is not None:
                pyname += " ~ [<a href=\"#X%s\" >%d</a>] (d = %.2f) " % (finals[c.similar[0]].aname, 
                    c.similar[0] + 1, c.similar[1])
            c.details = [
                [ Paragraph(pyname + "<font face=\"Helvetica\">#%s</font> %s (<font color=\"%s\">%s</font>)" % 
                    (c.tname, mini_name(c.minimum_type), multi_color(c.multiplicity), 
                    multi_name(c.multiplicity)), self.stylesheet['cellbody']), 
                    Paragraph("""<font face=\"Helvetica\">E</font> = %s (%.3f eV)""" % 
                    (("%.5f" if c.energy <= -1000 else "%.6f") % c.energy, 
                    httoev * (c.energy - mener)), self.stylesheet['cellbody']), 
                    Paragraph("""<font face=\"Helvetica\">[Relax] Time: </font> %s 
                    <font face=\"Helvetica\">Step: </font> %d""" % (time_span_short(c.props["time"]), 
                    int(c.props["step"])), self.stylesheet['cellbody']) ], 
                ]
            if "freqs" in c.props:
                c.detailsx = [
                    [ Paragraph("""<font face=\"Helvetica\">Freqs: </font> %s""" % 
                        freq_str(c.props["freqs"], fwid), self.stylesheet['cellbody']), "", 
                        Paragraph("""<font face=\"Helvetica\">[Freq] Time: </font> %s""" % 
                        time_span_short(c.props["time.freqs"]), self.stylesheet['cellbody']) ], 
                    [ cdr[i], "", "" ] ]
            else:
                c.detailsx = [ [ cdr[i], "", "" ] ]
            data += c.details + c.detailsx
            len_ed = len(data)
            style._cmds += [
                ('SPAN',    (0, len_ed - 1), (-1, len_ed - 1)), 
            	('LINEBELOW',    (0, len_ed - 1), (-1, len_ed - 1), 1, "#A5A5A5")
            ]
            if "freqs" in c.props:
                style._cmds += [ ('SPAN',    (0, len_ed - 2), (1, len_ed - 2)) ]
        style._cmds += [
            ('FONT',    (0, 0), (-1, -1), "Helvetica"), 
            ('FONTSIZE',    (0, 0), (-1, -1), 10), 
            ("BACKGROUND", (0, 0), (0, -1), "#FFFFFF"), 
            ('ALIGN',    (0, 0), (0, -1), "CENTER"), 
        ]
        return Table(data, (awid * 4 / 10, awid * 2.8 / 10, awid * 3.2 / 10), style=style, hAlign="LEFT")

    def draw_cluster(self, finals):
        if len(finals) == 0: return None
        cdx = []
        cdr = []
        mwid = 0
        mener = finals[0].energy if len(finals) != 0 else 0.0
        efactor = 27.21139
        rwid = (letter[0] / 2 - 0.6 * inch) * 0.52
        mwid = (letter[0] / 2 - 0.6 * inch) * 0.48
        qwid = stringWidth("Freqs: ", fontName="Helvetica", fontSize=10)
        fwid = rwid - 12 - 5 - qwid
        rots = [ [90.0, 1.0, 0.0, 0.0], [90.0, 0.0, 1.0, 0.0], [45.0, 1.0, 1.0, 0.0] ]
        rotns = [ "TOP", "LEFT", "TILTED" ]
        rt1, rt2 = Report.test_width(finals)
        for c in finals:
            xg = DrawCluster(c, simple=False, ratio=rt1)
            cdx.append(xg)
            xa, xb, xc = [ DrawCluster(c, simple=True, ratio=rt2, rotmat=r, title=n)
                for r, n in zip(rots, rotns) ]
            cdr.append(Table([[xa, xb, xc]], (xa.width, xb.width, xc.width), 
                style=TableStyle([ ('LEFTPADDING', (0, 0), (-1, -1), 0)]), hAlign='LEFT'))
        data = []
        for i in range(len(cdx)):
            c = finals[i]
            pyname = "<a name=\"X%s\" />[<font face=\"Helvetica\">%d</font>]"    % (c.aname, i + 1)
            if c.similar is not None:
                pyname += " ~ [<a href=\"#X%s\" >%d</a>] (d = %.2f)" % (finals[c.similar[0]].aname, 
                    c.similar[0] + 1, c.similar[1])
            cpt = c.pointgroup()
            ptp = Paragraph("<font face=\"HelveticaItalic\">%s<sub rise=2 size=8>%s</sub></font>" 
                % (cpt[0], cpt[1:]), self.stylesheet['cellbody'])
            c.details = [
                [ Paragraph(pyname, self.stylesheet['cellbody']), 
                    Paragraph("<font face=\"Helvetica\">#%s</font> %s (<font color=\"%s\">%s</font>)" % 
                    (c.tname, mini_name(c.minimum_type), multi_color(c.multiplicity), 
                    multi_name(c.multiplicity)), self.stylesheet['cellbody']), ptp], 
                [cdx[i], Paragraph("""<font face=\"Helvetica\">E</font> = %s (%.3f eV)""" % 
                    (("%.5f" if c.energy <= -1000 else "%.6f") % c.energy, 
                    httoev * (c.energy - mener)), self.stylesheet['cellbody']), ""], 
                ["", Paragraph("""<font face=\"Helvetica\">[Relax] Time: </font> %s 
                    <font face=\"Helvetica\">Step: </font> %d""" % (time_span_short(c.props["time"]), 
                    int(c.props["step"])), self.stylesheet['cellbody']), ""] ]
            if "vip.energy" in c.props:
                c.details += [
                    [cdx[i], Paragraph("""<font face=\"Helvetica\">E<sub rise=2 size=8>vip</sub></font> = %s (%.3f eV)""" % 
                        (("%.5f" if c.props["vip.energy"] <= -1000 else "%.6f") % c.props["vip.energy"], 
                        httoev * c.props["vip.energy"]), self.stylesheet['cellbody']), ""], 
                ]
            if "freqs" in c.props:
                c.detailsx = [
                ["", Paragraph("""<font face=\"Helvetica\">Freqs: </font> %s""" % 
                    freq_str(c.props["freqs"], fwid), self.stylesheet['cellbody']), "" ], 
                ["", Paragraph("""<font face=\"Helvetica\">[Freq] Time: </font> %s""" % 
                    time_span_short(c.props["time.freqs"]), self.stylesheet['cellbody']), "" ], ["", cdr[i], ""] ]
            else:
                c.detailsx = [ ["", "", ""], ["", "", ""], ["", cdr[i], ""] ]
            c.detailsy = [ ["", cdr[i], ""] ]
        style = copy.deepcopy(self.tablestyles["main"])
        data = []
        for i in range(0, len(finals), 2):
            if i + 1 >= len(finals):
                long_mode = "freqs" in finals[i].props
            elif "freqs" in finals[i].props or "freqs" in finals[i + 1].props:
                long_mode = True
            else:
                long_mode = False
            xa = finals[i].details
            xb = finals[i].detailsx if long_mode else finals[i].detailsy
            xa += xb
            if i + 1 >= len(finals):
                ya = [ ["", "", ""] ] * len(xa)
            else:
                ya = finals[i + 1].details
                yb = finals[i + 1].detailsx if long_mode else finals[i + 1].detailsy
                ya += yb
            len_st = len(data)
            for j in range(len(xa)):
                data += [ xa[j] + ya[j] ]
            len_ed = len(data) - 1
            style._cmds += [
                ('SPAN',    (0, len_st + 1), (0, len_ed)), 
                ('SPAN',    (3, len_st + 1), (3, len_ed)), 
            	('LINEBELOW',    (0, len_ed), (-1, len_ed), 1, "#A5A5A5")
            ]
            for ix in range(len_st + 1, len_ed + 1):
                style._cmds += [
                    ('SPAN',    (1, ix), (2, ix)), ('SPAN',    (4, ix), (5, ix)), 
                ]
        style._cmds += [
            ('FONT',    (0, 0), (-1, -1), "Helvetica"), 
            ('FONTSIZE',    (0, 0), (-1, -1), 10), 
            ("BACKGROUND", (0, 0), (0, -1), "#FFFFFF"), 
            ('ALIGN',    (0, 0), (0, -1), "CENTER"), 
            ('ALIGN',    (3, 0), (3, -1), "CENTER"), 
        ]
        # if not self.no_freq:
            # style._cmds += [
            #     ("BACKGROUND", (0, 0), (0, -1), "#F4F4F4"), 
            #     ("BACKGROUND", (3, 0), (3, -1), "#F4F4F4"), 
            # ]
        cwid = 27
        return Table(data, (mwid, rwid - cwid, cwid, mwid, rwid - cwid, cwid), 
            style=style, hAlign="LEFT")

    def draw_path_at_surface(self, paths, surface_depth=3.0):
        if len(paths) == 0:
            return None
        style = copy.deepcopy(self.tablestyles["main"])
        awid = letter[0] - 1.2 * inch
        data = []
        _, rto = Report.test_width(paths.values()[0].structs)
        for i, (xi, p) in enumerate(sorted(paths.items(), key=lambda x: int(x[0]))):
            print ('Path %d / %d' % (i + 1, len(paths)))
            tid = int(xi)
            imx = ""
            if p.focus:
                imx += "*"
            if p.ismax:
                imx += " (max)"
            if not p.isdirect:
                imx += " (int)"
            imx += " H = %.3f / %.3f eV" % tuple(httoev * p.barriers)
            data += [[Paragraph("""[<font face=\"Helvetica\">%d</font>] <font face=\"Helvetica\">
                                #%d</font>%s""" % (i, tid, imx), self.stylesheet['cellbody']),
                      "From %d to %d [d = %.2f]" % (p.from_to[0], p.from_to[1], p.diffs[0]),
                      Paragraph("""<font face=\"Helvetica\">[Path] Time: </font> %s
                                <font face=\"Helvetica\">Step: </font> %d""" %
                                (time_span_short(p.props["time"]), int(p.props["step"])),
                                self.stylesheet['cellbody'])]]
            data += [[PGPlot.path_plot(p.energies, p.structs, p.diffs, surface_depth,
                                       ratio=rto), "", ""]]
            len_ed = len(data)
            style._cmds += [
                ('SPAN',      (0, len_ed - 1), (-1, len_ed - 1)), 
            	('LINEBELOW', (0, len_ed - 1), (-1, len_ed - 1), 1, "#A5A5A5")
            ]
        style._cmds += [
            ("FONT",       (0, 0), (-1, -1), "Helvetica"),
            ("FONTSIZE",   (0, 0), (-1, -1), 10),
            ("BACKGROUND", (0, 0), ( 0, -1), "#FFFFFF"),
            ("ALIGN",      (0, 0), ( 0, -1), "CENTER"),
        ]
        return Table(data, (awid * 0.4, awid * 0.28, awid * 0.32), style=style, hAlign="LEFT")
    
    def draw_path_min_at_surface(self, mins, rto_scale = 1.5, surface_depth=3.0, lcol=5):
        if len(mins) == 0:
            return None
        style = copy.deepcopy(self.tablestyles["main"])
        awid = self.doc.width - 1.2 * inch
        dataog = []
        _, rto = Report.test_width(mins.values()[0:5], clip=True)
        rto *= rto_scale
        emin = min([p.energy for p in mins.values()])
        from meta.ana_utils.plot_rpath import to_letters
        for i, (xi, p) in enumerate(sorted(mins.items(), key=lambda x: int(x[0]))):
            tid = int(xi)
            dataog += [[Paragraph("""[<font face=\"Helvetica\">%s</font>] E = %.3f eV"""
                      % (to_letters(tid), (p.energy - emin) * httoev), self.stylesheet['cellbody'])]]
            dataog += [[DrawCluster(p, simple=True, surface_depth=surface_depth,
                        clip=True, ratio=rto)]]
        data = []
        empty = Paragraph("", self.stylesheet['cellbody'])
        while (len(dataog) / 2) % lcol != 0:
            dataog.extend([[empty], [empty]])
        for i in range(0, len(dataog), lcol * 2):
            data += [dataog[i:i+lcol*2:2], dataog[i+1:i+lcol*2:2]]
            len_ed = len(data)
            style._cmds += [
            	('LINEBELOW', (0, len_ed - 1), (-1, len_ed - 1), 1, "#A5A5A5")
            ]
        style._cmds += [
            ("FONT",       (0, 0), (-1, -1), "Helvetica"),
            ("FONTSIZE",   (0, 0), (-1, -1), 10),
            ("BACKGROUND", (0, 0), ( 0, -1), "#FFFFFF"),
            ("ALIGN",      (0, 0), ( 0, -1), "CENTER"),
        ]
        return Table(data, (awid * (1.0 / lcol), ) * lcol, style=style, hAlign="LEFT")

    def prop_graph(self, desp, idesp, name, pdesp, rot=False, vib=True, max_count=None, tmax=200.0):
        sub_title = self.heading(desp, key=idesp, level=1)
        idx = 0
        story = []
        trange = np.arange(0.5, tmax, 2.0)
        for dx in self.d.d:
            for g in dx.groups:
                if len(g.finals) != 0 and isinstance(g.finals[0], Cluster):
                    if name not in g.finals[0].props: continue
                    if vib and "freqs" not in g.finals[0].props: continue
                    sub_story = []
                    if idx == 0:
                        sub_story = [ sub_title ]
                        idx = 1
                    gname = "%s.%s" % (dx.name, g.name)
                    sub_story.append(self.heading(g.description, cont="Group %s" % gname, 
                        key="%s-gp-%s" % (idesp, gname), level=2))
                    style = copy.deepcopy(self.tablestyles["mainxx"])
                    mener = g.finals[0].energy
                    tdata = []
                    for t in trange:
                        sumz = 0.0
                        sump = 0.0
                        for c in (g.finals if max_count is None else g.finals[:max_count]):
                            if vib and "freqs" not in c.props: continue
                            if name not in c.props: continue
                            zc = gas_phase_partition(c, t, rel_energy=mener, rot=rot, vib=vib)
                            sumz += zc
                            sump += zc * c.props[name]
                        pv = sump / sumz * httoev
                        if np.isfinite(pv): tdata.append([t, pv])
                    sub_story.append(PGPlot.prop_plot(tdata, ymin=True, ylabel=pdesp))
                    story.append(KeepTogether(sub_story))
        return story

    def heat_capacity(self, desp, idesp, rot=False, vib=True, max_count=None, tmax=200.0):
        sub_title = self.heading(desp, key=idesp, level=1)
        idx = 0
        story = []
        trange = np.arange(0.5, tmax, 2.0)
        for dx in self.d.d:
            for g in dx.groups:
                if len(g.finals) != 0 and isinstance(g.finals[0], Cluster):
                    if vib and "freqs" not in g.finals[0].props: continue
                    sub_story = []
                    if idx == 0:
                        sub_story = [ sub_title ]
                        idx = 1
                    gname = "%s.%s" % (dx.name, g.name)
                    sub_story.append(self.heading(g.description, cont="Group %s" % gname, 
                        key="%s-gp-%s" % (idesp, gname), level=2))
                    style = copy.deepcopy(self.tablestyles["mainxx"])
                    mener = g.finals[0].energy
                    tdata = []
                    for t in trange:
                        beta = 1 / (t * ktoht)
                        suma = 0.0
                        sumb = 0.0
                        sumz = 0.0
                        for c in (g.finals if max_count is None else g.finals[:max_count]):
                            if vib and "freqs" not in c.props: continue
                            zc = gas_phase_partition(c, t, rel_energy=mener, rot=rot, vib=vib)
                            wc = gas_phase_partition_d(c, t, rel_energy=mener, rot=rot, vib=vib)
                            vc = gas_phase_partition_dd(c, t, rot=rot, vib=vib)
                            suma += wc * zc
                            sumz += zc
                            sumb += wc**2 * zc + vc * zc
                        cv = beta**2 * (-(suma / sumz)**2 + sumb / sumz)
                        if np.isfinite(cv): tdata.append([t, cv])
                    sub_story.append(PGPlot.prop_plot(tdata))
                    story.append(KeepTogether(sub_story))
        return story

    def probabilities(self, desp, idesp, rot=False, vib=True, max_prop_count=10):
        sub_title = self.heading(desp, key=idesp, level=1)
        idx = 0
        story = []
        for dx in self.d.d:
            for g in dx.groups:
                if len(g.finals) == 0:
                    continue
                is_surf = isinstance(g.finals[0], ClusterAtSurface)
                if vib and "freqs" not in g.finals[0].props:
                    continue
                sub_story = []
                if idx == 0:
                    sub_story = [ sub_title ]
                    idx = 1
                gname = "%s.%s" % (dx.name, g.name)
                sub_story.append(self.heading(g.description, cont="Group %s" % gname, 
                    key="%s-gp-%s" % (idesp, gname), level=2))
                style = copy.deepcopy(self.tablestyles["mainxx"])
                mener = g.finals[0].energy
                ets = [ [ "name", "multiplicity", 
                    "ZPE(eV)" if vib else "E (eV)", "P (100K)", 
                    "P (300K)", "P (400K)", "P (500K)", "P (700K)", "P (1000K)" ] ]
                ets[0] = [Paragraph(i, self.stylesheet['celltitle']) for i in ets[0]]
                ets[0][0].style = self.stylesheet['celltitleleft']
                props, propsum = [], []
                for itt, t in enumerate([ 100, 200, 300, 500, 700, 1000 ]):
                    if is_surf:
                        propx = [ surface_partition(c, t, rel_energy=mener, vib=vib, 
                            iprint=(itt == 0)) for c in g.finals ]
                    else:
                        propx = [ gas_phase_partition(c, t, rel_energy=mener, rot=rot, vib=vib, 
                            iprint=(itt == 0)) for c in g.finals ]
                    props.append(propx)
                    propsum.append(sum(propx))
                for ic, c in enumerate(g.finals[:max_prop_count]):
                    if ic < len(g.finalsx):
                        pyname = "<a href=\"#X%s\" >[<font face=\"Helvetica\">%d</font>]</a>"    % (c.aname, ic + 1)
                    else:
                        pyname = "[<font face=\"Helvetica\">%d</font>]"    % (ic + 1)
                    pyname += " <font face=\"Helvetica\">#%s</font> %s"    % (c.tname, mini_name(c.minimum_type)[0]
                        if mini_name(c.minimum_type) != "UNK" else "")
                    if c.props["nrep"] != 1:
                        pyname += " <font color=\"#AAAAAA\" face=\"HelveticaItalic\">* %d</font>" % c.props["nrep"]
                    muname = "<font color=\"%s\">%s</font>" % (multi_color(c.multiplicity), multi_name(c.multiplicity))
                    if not is_surf:
                        cpt = c.pointgroup()
                        if cpt != "C1":
                            muname += " &middot; <font face=\"HelveticaItalic\">%s<sub rise=2 size=8>%s</sub></font>" \
                                % (cpt[0], cpt[1:])
                    if not vib:
                        ee = httoev * (c.energy - mener)
                    elif not is_surf:
                        ee = httoev * gas_phase_zpe(c)
                    else:
                        ee = httoev * surface_zpe(c)
                    etsf = [ pyname, muname, "%.5f" % ee ]
                    for itt, t in enumerate([ 100, 200, 300, 500, 700, 1000 ]):
                        etsf.append("%.5f" % (np.float64(props[itt][ic] / propsum[itt])))
                    etsf[0] = Paragraph(etsf[0], self.stylesheet['cellbody'])
                    etsf[1] = Paragraph(etsf[1], self.stylesheet['cellcenter'])
                    ets.append(etsf)
                tnrep = sum([c.props["nrep"] for c in g.finals[:max_prop_count]])
                ttnrep = sum([c.props["nrep"] for c in g.finals])
                etsf = [ "<font face=\"Helvetica\">total</font> (%d of %d)" % 
                    (len(g.finals[:max_prop_count]), len(g.finals)), 
                    "<font size=\"8\" color=\"#AAAAAA\" face=\"HelveticaItalic\">(* %d of %d)</font>"
                    % (tnrep, ttnrep), "" ]
                for itt, t in enumerate([ 100, 200, 300, 500, 700, 1000 ]):
                    etsf.append("%.5f" % (sum(props[itt][:max_prop_count]) / propsum[itt]))
                etsf[0] = Paragraph(etsf[0], self.stylesheet['cellbody'])
                etsf[1] = Paragraph(etsf[1], self.stylesheet['cellbody'])
                ets.append(etsf)
                style._cmds += [ ('ALIGN', (1, 1), (-1, -1), "CENTER"), 
                    ('ALIGN', (0, 0), (0, 0), "LEFT"), 
                    ('LINEBELOW',    (0, -2), (-1, -2), 1, "#C5C5C5"), 
                    ('LINEAFTER',    (2, 0), (2, -1), 1, "#C5C5C5"), ]
                sub_story.append(self.draw_simple_table(ets, style=style, 
                    widths=(90, 65, 50, 52, 52, 52, 52, 52, 57), span=range(3, 8)))
                story.append(KeepTogether(sub_story))
        return story

    def simple_wrap(self, data):
        return Table([[data]], style=[("LEFTPADDING", (0, 0), (0, 0), 15)], hAlign="LEFT")

    def draw_simple_table(self, data, style=None, widths=None, span=[], noacross=None):
        if style is None: style = self.tablestyles["main"]
        if len(span) != 0:
            style = copy.deepcopy(style)
            otsg, ostc = None, None
            l = 0
            for iix, ix in enumerate(data):
                if otsg is None:
                    otsg, otsc = [], []
                    for c in span:
                        otsg.append([[l, l]])
                        otsc.append(ix[c])
                else:
                    for i, c in enumerate(span):
                        if ix[c] != otsc[i] or (noacross is not None and i > noacross and iix != 0 and 
                            ix[span[noacross]] != data[iix - 1][span[noacross]]):
                            otsc[i] = ix[c]
                            otsg[i].append([l, l])
                        else:
                            otsg[i][-1][1] = l
                l += 1
            for i, c in enumerate(span):
                for x, y in otsg[i]:
                    if x != y:
                        style._cmds += [ ('SPAN', (c, x), (c, y)) ]
        if widths is None:
            table = PGTable(data, dx=15, style=style, hAlign="LEFT")
        else:
            table = PGTable(data, widths, dx=15, style=style, hAlign="LEFT")
        return table

    def heading(self, name, key, cont=None, level=1):
        if cont is None: cont = name
        p = Paragraph("<a name=\"%s\" />%s" % (key, name), self.stylesheet['heading%d' % level])
        p.toc_key = key
        p.toc_cont = cont
        return p
    
    def build_single(self, minx, surface_depth=3.0, max_prop_count=10, ratio=2.0, lcol=1, wid=4, hei=4):
        print ('building single pdf ...')
        self.surface_depth = surface_depth
        self.doc = TOCDocTemplate(self.filename, pagesize=(wid*inch, hei*inch),
            leftMargin=0.2*inch, rightMargin=0.2*inch, topMargin=0.2*inch, bottomMargin=0.2*inch)
        if isinstance(minx, list):
            mcc = { str(i): j for i, j in enumerate(minx) }
        else:
            mcc = { '0': minx }
        print (len(mcc))
        dcx = self.draw_path_min_at_surface(mcc, ratio, surface_depth, lcol=lcol)
        self.doc.build([dcx])

    def build_graph(self, surface_depth=3.0, max_prop_count=10):
        print ('building pdf-graph ...')
        self.surface_depth = surface_depth
        self.doc = TOCDocTemplate("graph.pdf", pagesize=(11*inch, 11*inch),
            leftMargin=0.5*inch, rightMargin=0.5*inch, topMargin=0.6*inch, bottomMargin=0.5*inch)
        story = []
        for dx in self.d.d:
            for g in dx.groups:
                if g.nebgraph is not None:
                    if len(g.nebgraph.vs) == 0:
                        dc = None
                    else:
                        dc = DrawRGraph(g.nebgraph, width=(self.doc.width - 2.2 * inch))
                    if len(g.nebmins) == 0:
                        dcx = None
                    else:
                        dcx = self.draw_path_min_at_surface(g.nebmins, surface_depth=surface_depth,
                            lcol=7)
                    if dc is not None:
                        story.append(dc)
                    if dcx is not None:
                        story.append(dcx)
        if len(story) != 0:
            self.doc.build(story)

    def build(self, surface_depth=3.0, max_prop_count=10):
        print ('building pdf ...')
        self.surface_depth = surface_depth
        story = []
        title = Paragraph('<b>' + self.d.title + '</b>', self.stylesheet['title'])
        story.append(title)
        story.append(Paragraph("""This is a computer automatically generated document. <br/> 
            Copyright &copy; %s Parallel Global Optimization Toolkit<br/>
            <font face="Helvetica">Program Author:</font> Huanchen Zhai (Alexandrova Research Group, UCLA)<br/>
            <font face="Helvetica">Date:</font> %s""" % (time.strftime("%Y"), 
            time.strftime("%b %d, %Y %H:%M:%S")), self.stylesheet['body']))
        
        story.append(self.heading('Contents', key="contents", level=1))
        toc = TableOfContents()
        toc.levelStyles = [ ParagraphStyle(fontName='Helvetica', fontSize=12, name='TOCHeading1', 
            leftIndent=40, firstLineIndent=-20, spaceBefore=0, leading=12),
            ParagraphStyle(fontName='Helvetica', fontSize=10, name='TOCHeading2',
                leftIndent=60, firstLineIndent=-20, spaceBefore=0, leading=10),
        ]
        story.append(toc)
        story.append(PageBreak())

        if self.energy_cutoff_hard is not None and self.energy_cutoff_hard > 0.0:
            for dx in self.d.d:
                for g in dx.groups:
                    min_e = g.finals[0].energy if len(g.finals) != 0 else 0.0
                    cutoff = self.energy_cutoff_hard
                    g.finalsx = [c for c in g.finals if c.energy < (cutoff / httoev) + min_e]
        else:
            for dx in self.d.d:
                for g in dx.groups:
                    g.finalsx = g.finals

        sub_title = self.heading('Energetics Summary', key="summary", level=1)
        for idx, dx in enumerate(self.d.d):
            sub_story = []
            if idx == 0:
                sub_story = [ sub_title ]
            for g in dx.groups:
                gfxs = [ ('', g.finalsx) ]
                _, rto = Report.test_width(g.finalsx)
                for ifx, gfx in gfxs:
                    if len(ifx) != 0: ifxx = ifx.upper() + ' '
                    else: ifxx = ''
                    if len(gfx) == 0: continue
                    gname = "%s.%s" % (dx.name, g.name)
                    sub_story.append(self.heading(ifxx + g.description, cont="%sGroup %s" % (ifxx, gname), 
                        key="summary-gp%s-%s" % (ifx, gname), level=2))
                    es_text = "<font face=\"Helvetica\">Multiplicity color: </font>"
                    mults = sorted(list(set([d.multiplicity for d in gfx])))
                    ex_texx = [ "<font color=\"%s\">%s</font>" % (multi_color(m), multi_name(m)) for m in mults ]
                    es_text += ', '.join(ex_texx)
                    if len(gfx) != 0:
                        from cluster.coval import AtomicWeight, AtomicColor
                        es_text += ".&nbsp;&nbsp;<font face=\"Helvetica\">Element color: </font>"
                        xelemcs = list(gfx[0].elems)
                        if hasattr(gfx[0], "surf"): xelemcs += list(gfx[0].surf.elems)
                        elemcs = sorted(list(set([d for d in xelemcs])), key=lambda x:AtomicWeight.x[x])
                        ex_texx = [ "<font color=\"%s\">%s</font>" % (AtomicColor.x[x], x) for x in elemcs ]
                        es_text += ', '.join(ex_texx)
                    if gfx[0].aname is not None:
                        es_text += (".<br/>All energies are relative to <a href=\"X%s\" fontName=\"Helvetica\"> " + 
                            "the putative global minimum: %.6f [1] (#%s), <font color=\"%s\">%s</font></a>.") % \
                            (gfx[0].aname, gfx[0].energy, gfx[0].tname, 
                            multi_color(gfx[0].multiplicity),    multi_name(gfx[0].multiplicity))
                    else:
                        es_text += (".<br/>All energies are relative to <font face=\"Helvetica\"> " + 
                            "the putative global minimum: %.6f [1], <font color=\"%s\">%s</font></font>.") % \
                            (gfx[0].energy, multi_color(gfx[0].multiplicity),    multi_name(gfx[0].multiplicity))
                    sub_story.append(Paragraph(es_text, self.stylesheet['indbody']))
                    for ii, i in enumerate(self.energy_cutoff):
                        ener = DrawEnergetics(clus=gfx, cutoff=i, ratio=rto, no_freq=self.no_freq, 
                            surface_depth=surface_depth)
                        sub_story.append(Paragraph(
                            "%s<font face=\"Helvetica\">Range:</font> 0.00 ~ %.2f eV (%d structures)" % 
                            (("<br/>" if ii != 0 else ""), i, len(ener.clus)), self.stylesheet['indbody']))
                        sub_story.append(self.simple_wrap(ener))
            story.append(KeepTogether(sub_story))

        # cmd history
        sub_title = self.heading('Command History', key="hist", level=1)
        idx = 0
        for dx in self.d.d:
            sub_story = []
            if idx == 0:
                sub_story = [ sub_title ]
                idx = 1
            dxname = "[%s] %s" % (dx.name, dx.filename)
            sub_story.append(self.heading(dxname, cont=dxname, key="hist-i-%s" % dx.name, 
                level=2))
            sub_story.append(Paragraph("<br/>".join(dx.cmd_hist), self.stylesheet['body']))
            story.append(KeepTogether(sub_story))

        # global information
        story.append(self.heading('Global Information', key="global", level=1))
        global_info_data = [ ]
        style = copy.deepcopy(self.tablestyles["main"])
        for d in self.d.d:
            x = "File %s" % d.name
            statx = len(global_info_data)
            global_info_data += [
                [x, "operating system", "", d.xsystem["os"]], 
                [x, "CPU", "", "%s @ %s MHz" % (d.xsystem["CPU"], d.xsystem["MHz"])], 
                [x, "host name", "", d.xsystem["hostname"]], 
                [x, "cores per node", "", d.xsystem["cores"]], 
                [x, "file name", "", d.filename], 
                [x, "total CPU hours", "", d.xglobal["total_cpu_time"] / 3600], 
                [x, "times", "start", time.strftime("%b %d, %Y %H:%M:%S", d.xglobal["start_time"])], 
                [x, "times", "finish", time.strftime("%b %d, %Y %H:%M:%S", d.xglobal["finish_time"])], 
                [x, "times", "total", time_span_str(d.xglobal["total_time"])], 
                [x, "times", "running", time_span_str(d.xglobal["running_time"])], 
                [x, "times", "idle", time_span_str(d.xglobal["idle_time"])], 
            ]
            edx = len(global_info_data)
            for i in range(6):
                style._cmds += [ ('SPAN',    (1, i + statx), (2, i + statx) ) ]
            style._cmds += [ ('LINEBELOW',    (0, edx - 1), (-1, edx - 1), 1, "#A5A5A5") ]
        style._cmds += [
            ('ALIGN', (1, 0), (2, -1), "RIGHT"), 
            ('FONT',    (1, 0), (2, -1), "Helvetica"), 
        ]
        story.append(self.draw_simple_table(global_info_data, span=[0, 1], style=style))
        story = story[:-2] + [ KeepTogether(story[-2:]) ]

        # energy computation
        story.append(self.heading('Energy Computation Parameters', key="comp-param", level=1))
        dft_params = [
            [ "group id", "program", "basis", "functional", "charge", 
                "mini step", "max step", "disp step" ]
        ]
        for d in self.d.d:
            for g in d.groups:
                dft_params += [
                    [ "%s.%s" % (d.name, g.name), g.param["program"], g.param["basis"], 
                        g.param["functional"], g.param["charge"], 
                        g.param["step"], g.param["max_step"], g.param["freq_step"] ]
                ]
        story.append(self.draw_simple_table(dft_params, style=self.tablestyles["mainx"], 
            span=range(1, 8)))
        story = story[:-2] + [ KeepTogether(story[-2:]) ]

        # Local Optimization Options
        story.append(self.heading('Local Optimization Options', key="local-opts", level=1))
        opts_texts = []
        for d in self.d.d:
            for r in d.runs:
                opts_texts += [
                    """<font face="Helvetica">%s-%s (%s): </font> %s""" % (
                        "%s.%s" % (d.name, r.parent.parent.name), r.pname, 
                        r.task_type + ("*" if r.focus == r else ""), r.param["others"])
                ]
        story.append(Paragraph("<br />".join(opts_texts), self.stylesheet['body']))
        story = story[:-2] + [ KeepTogether(story[-2:]) ]

        # creation
        story.append(self.heading('Creation Parameters', key="creation-param", level=1))
        cre_names = []
        cre_data = []
        for d in self.d.d:
            for s in d.stages:
                tname = "%s.%s.%s" % (d.name, s.parent.name, s.stage)
                found = None
                for i, c in enumerate(cre_data):
                    if s.creation_param == c:
                        found = i
                        break
                if found is None:
                    cre_names.append([tname])
                    cre_data.append(s.creation_param)
                else:
                    cre_names[i].append(tname)
        creation_params = [ ]
        style = copy.deepcopy(self.tablestyles["main"])
        style._cmds += [
            ('ALIGN', (1, 0), (1, -1), "RIGHT"), 
            ('FONT',    (1, 0), (1, -1), "Helvetica"), 
        ]
        for n, d in zip(cre_names, cre_data):
            x = Paragraph("<br/>".join(n), self.stylesheet['celltitle'])
            st_stat = len(creation_params)
            creation_params += [
                [x, "name", Paragraph(d["fname"], self.stylesheet['cellbody']), 1], 
                [x, "method", creation_method_name(d["method"]), 2], 
                [x, "number", str(d["number"]), 3], 
            ]
            for k in range(3):
                style._cmds += [ ('SPAN', (2, st_stat + k), (3, st_stat + k)) ]
            if d["method"] == "blda":
                st_stat = len(creation_params)
                creation_params += [[x,    "order", d["order"], "" ]]
                style._cmds += [ ('SPAN', (2, st_stat), (3, st_stat)) ]
                for k in sorted(d["dist"]["mu"].keys()):
                    mu, sig = d["dist"]["mu"][k], d["dist"]["sigma"][k]
                    sta = "statistics"
                    if len(mu) == 2:
                        creation_params += [[x, sta, k, Paragraph("""<font face="Helvetica" color="#5577AA">[1]</font> 
                            %.3f Â± %.3f&nbsp;&nbsp;<font face="Helvetica" color="#5577AA">[2]</font>
                            %.3f Â± %.3f""" % (mu[0], sig[0], mu[1], sig[1]), self.stylesheet['cellbody'])]]
                    elif len(mu) == 1:
                        creation_params += [[x, sta, k, Paragraph("""<font face="Helvetica" color="#5577AA">[1]</font> 
                            %.3f Â± %.3f""" % (mu[0], sig[0]), self.stylesheet['cellbody'])]]
            elif d["method"] == "ck":
                st_stat = len(creation_params)
                creation_params += [[x,    "drel", d["drel"], "" ]]
                style._cmds += [ ('SPAN', (2, st_stat), (3, st_stat)) ]
            edx = len(creation_params)
            style._cmds += [ ('LINEBELOW',    (0, edx - 1), (-1, edx - 1), 1, "#A5A5A5") ]
        story.append(self.draw_simple_table(creation_params, style=style, widths=(60, 60, 120, 180), 
            span=[0, 1]))
        story = story[:-2] + [ KeepTogether(story[-2:]) ]

        # Filtering
        story.append(self.heading('Filtering Parameters', key="filter-param", level=1))
        filtering_params = [
            [ "group id", "creation", "", "runtime", "", "final", "" ], 
            [ "", "max diff", "report", "max diff", "report", "max diff", "report" ]
        ]
        for d in self.d.d:
            for g in d.groups:
                filtering_params += [
                    [ "%s.%s" % (d.name, g.name), 
                        "%.2f" % g.param["filtering"]["create"]["max_diff"], 
                        "%.2f" % g.param["filtering"]["create"]["max_diff_report"], 
                        "%.2f" % g.param["filtering"]["run"]["max_diff"], 
                        "%.2f" % g.param["filtering"]["run"]["max_diff_report"], 
                        "%.2f" % g.param["filtering"]["final"]["max_diff"], 
                        "%.2f" % g.param["filtering"]["final"]["max_diff_report"]
                    ]
                ]
        style = copy.deepcopy(self.tablestyles["main"])
        style._cmds += [
            ('ALIGN', (0, 0), (-1, -1), "CENTER"), 
            ('LINEBELOW',    (0, 1), (-1, 1), 1, "#C5C5C5"), 
            ('FONT',    (0, 0), (-1, 1), "Helvetica"), 
            ('SPAN',    (1, 0), (2, 0) ), ('SPAN',    (3, 0), (4, 0) ), ('SPAN',    (5, 0), (6, 0) ), 
            ('SPAN',    (0, 0), (0, 1) ), 
        ]
        story.append(self.draw_simple_table(filtering_params, style=style, span=range(1, 7)))
        story = story[:-2] + [ KeepTogether(story[-2:]) ]

        # Job Parallel Information
        story.append(self.heading('Job Parallelism Information', key="job-para", level=1))
        jp_info = [
            [ "group <br/> id", "run id", "task", "average available <br/> processors", 
                "average running<br/> processors", 
                "maximum available <br/> processors", "maximum running <br/> processors", 
                "cores per <br/>processors" ]
        ]
        jp_info[0] = [Paragraph(i, self.stylesheet['celltitle']) for i in jp_info[0]]
        style = copy.deepcopy(self.tablestyles["mainxx"])
        for d in self.d.d:
            for r in d.runs:
                jp_info += [
                    [ "%s.%s" % (d.name, r.parent.parent.name), r.pname, r.task_type, 
                    int(r.avail_procs.mean()), int(r.run_procs.mean()), 
                    int(r.avail_procs.max()), int(r.run_procs.max()), int(r.factor) ]
                ]
            edx = len(jp_info)
            style._cmds += [ ('LINEBELOW',    (0, edx - 1), (-1, edx - 1), 1, "#A5A5A5") ]
        story.append(self.draw_simple_table(jp_info, style=style, 
            widths=(40, 40, 40, 70, 70, 70, 70, 70), span=range(0, 3), noacross=0))
        story = story[:-2] + [ KeepTogether(story[-2:]) ]

        # Job Time Information
        story.append(self.heading('Job Time Information', key="job-time", level=1))
        jt_info = [ [ "group id", "run id", "task", "total wall time", "total node hours", 
                "total CPU hours", "total cost CPU hours" ] ]
        jt_info[0] = [Paragraph(i, self.stylesheet['celltitle']) for i in jt_info[0]]
        style = copy.deepcopy(self.tablestyles["mainxx"])
        for d in self.d.d:
            for r in d.runs:
                jt_info += [
                    [ "%s.%s" % (d.name, r.parent.parent.name), r.pname, r.task_type, 
                    time_span_short(r.xtimes["walltime"]), int(r.xtimes["nodetime"] / 3600),
                    int(r.xtimes["cputime"] / 3600), int(r.xtimes["costtime"] / 3600) ] ]
            edx = len(jt_info)
            style._cmds += [ ('LINEBELOW',    (0, edx - 1), (-1, edx - 1), 1, "#A5A5A5") ]
        story.append(self.draw_simple_table(jt_info, style=style, 
            widths=(40, 40, 40, 60, 60, 60, 60), span=range(0, 3), noacross=0))
        story = story[:-2] + [ KeepTogether(story[-2:]) ]

        # Primary Configurations
        story.append(self.heading('Primary Configurations', key="prim-config", level=1))
        pcf = [ [ "group id", "run id", "sources", 
                "total", "failed", "max", "converged", "minima", "transition states", "duplicates"]]
        pcf[0] = [Paragraph(i, self.stylesheet['celltitle']) for i in pcf[0]]
        style = copy.deepcopy(self.tablestyles["mainxx"])
        for d in self.d.d:
            for r in d.runs:
                if r.sour_name == "unknown": continue
                pcf += [ [ "%s.%s" % (d.name, r.parent.parent.name), r.pname, 
                        r.sour_name, len(r.xdata.runtime["states"].values()), 
                        len([x for x in r.xdata.runtime["states"].values() if x in [0, 2, 3]]), 
                        len([x for x in r.xdata.runtime["states"].values() if x == 4]), 
                        len([x for x in r.xdata.runtime["states"].values() if x == 1]), 
                        len([x for x in r.optima.values() if x.ttype == "local" and x.minimum_type == 1]), 
                        len([x for x in r.optima.values() if x.ttype == "local" and x.minimum_type == 2]), 
                        len([x for x in r.xdata.runtime["states"].values() if x in [5, 6]])
                    ] ]
            edx = len(pcf)
            style._cmds += [ ('LINEBELOW',    (0, edx - 1), (-1, edx - 1), 1, "#A5A5A5") ]
        story.append(self.draw_simple_table(pcf, style=style, 
            widths=(40, 40, 80, 35, 40, 35, 45, 50, 40, 38), span=range(0, 3), noacross=1))
        story = story[:-2] + [ KeepTogether(story[-2:]) ]

        # Primary Times
        story.append(self.heading('Primary Times', key="prim-time", level=1))
        pcf = [ [ "group id", "run id", "sources", 
                "total", "failed", "converged", "minima", "transition states", "duplicates"]]
        pcf[0] = [Paragraph(i, self.stylesheet['celltitle']) for i in pcf[0]]
        style = copy.deepcopy(self.tablestyles["mainxx"])
        for d in self.d.d:
            for r in d.runs:
                if r.sour_name == "unknown": continue
                pp = np.array([ 
                    sum([x.props["time"] for x in r.optima.values()]), 
                    sum([x.props["time"] for x in r.optima.values() if r.xdata.runtime["states"][str(x.tid)] == 0]), 
                    sum([x.props["time"] for x in r.optima.values() if r.xdata.runtime["states"][str(x.tid)] == 1]), 
                    sum([x.props["time"] for x in r.optima.values() if x.ttype == "local" and x.minimum_type == 1]), 
                    sum([x.props["time"] for x in r.optima.values() if x.ttype == "local" and x.minimum_type == 2]), 
                    sum([x.props["time"] for x in r.optima.values() if r.xdata.runtime["states"][str(x.tid)] in [5, 6]])
                ])
                pp = [int(x) for x in pp / 3600 * r.factor]
                pcf += [ [ "%s.%s" % (d.name, r.parent.parent.name), r.pname, r.sour_name ] + pp ]
            edx = len(pcf)
            style._cmds += [ ('LINEBELOW',    (0, edx - 1), (-1, edx - 1), 1, "#A5A5A5") ]
        story.append(self.draw_simple_table(pcf, style=style, 
            widths=(40, 40, 80, 35, 40, 70, 50, 40, 70), span=range(0, 3), noacross=1))
        story = story[:-2] + [ KeepTogether(story[-2:]) ]

        # Filtered Optima
        story.append(self.heading('Filtered Configurations', key="filtered-config", level=1))
        pfc = [
            [ "group id", "multiplicity",    "minima", "transition states", "unknowns", 
                "filtered minima", "filtered transition states", "filtered unknowns"]]
        pfc[0] = [Paragraph(i, self.stylesheet['celltitle']) for i in pfc[0]]
        style = copy.deepcopy(self.tablestyles["mainxx"])
        for d in self.d.d:
            for g in d.groups:
                for im, m in sorted(g.multi_stat.items(), key=lambda x: int(x[0])):
                    pfc += [ [ "%s.%s" % (d.name, g.name), multi_name(int(im)) ] + m ]
                line = len(pfc) - 1
                style._cmds += [ ('LINEBELOW',    (0, line), (-1, line), 1, "#C5C5C5"), ]
                if len(g.multi_stat) != 1:
                    pfc += [ [ "%s.%s" % (d.name, g.name), "total" ] + g.stat ]
                    line = len(pfc) - 1
                    style._cmds += [ ('LINEBELOW',    (0, line), (-1, line), 1, "#C5C5C5"), ]
            edx = len(pfc)
            style._cmds += [ ('LINEBELOW',    (0, edx - 1), (-1, edx - 1), 1, "#A5A5A5") ]
        story.append(self.draw_simple_table(pfc, style=style, 
            widths=(40, 55, 50, 60, 60, 50, 60, 60), span=range(0, 1)))
        story = story[:-2] + [ KeepTogether(story[-2:]) ]

        # Optimization Time and Steps
        story.append(self.heading('Optimization Time and Steps', key="time-step", level=1))
        ots = [
            [ "group id", "run id", "task", "sources", "number", 
                "steps", "steps", "steps", "time", "time", "time"], 
            [ "", "", "", "", "", "mean", "max", "min", "mean", "max", "min"]]
        ots[0] = [Paragraph(i, self.stylesheet['celltitle']) for i in ots[0]]
        ots[1] = [Paragraph(i, self.stylesheet['celltitle']) for i in ots[1]]
        style = copy.deepcopy(self.tablestyles["mainxx"])
        for d in self.d.d:
            for r in d.runs:
                if r.sour_name == "unknown": continue
                xsteps = np.array([ x.props["step"] for x in r.optima.values() ])
                xtimes = np.array([ x.props["time"] for x in r.optima.values() ])
                if len(xtimes) == 0:
                    xsteps = np.array([0.0])
                if len(xtimes) == 0:
                    xtimes = np.array([0.0])
                ots += [ [ "%s.%s" % (d.name, r.parent.parent.name), r.pname, 
                    r.task_type, r.sour_name, len(r.xdata.runtime["states"].values()) ] + 
                    map(int, [ xsteps.mean(), xsteps.max(), xsteps.min() ]) + 
                    map(time_span_short, [ xtimes.mean(), xtimes.max(), xtimes.min() ]) ]
            edx = len(ots)
            style._cmds += [ ('LINEBELOW',    (0, edx - 1), (-1, edx - 1), 1, "#A5A5A5") ]
        style._cmds += [ 
            ('SPAN', (5, 0), (7, 0)), ('SPAN', (8, 0), (10, 0)), 
            ('SPAN', (0, 0), (0, 1)), ('SPAN', (1, 0), (1, 1)), ('SPAN', (2, 0), (2, 1)), 
            ('SPAN', (3, 0), (3, 1)), ('SPAN', (4, 0), (4, 1)), 
            ('LINEBELOW',    (0, 1), (-1, 1), 1, "#C5C5C5"), 
            ('LINEBELOW',    (0, 0), (-1, 0), 0.8, "#E1E1E1"), 
        ]
        story.append(self.draw_simple_table(ots, style=style, 
            widths=(40, 40, 50, 85, 40, 37, 37, 37, 45, 45, 45), span=range(0, 3), noacross=0))
        story = story[:-2] + [ KeepTogether(story[-2:]) ]

        # Energetics
        story.append(self.heading('Energetics', key="energetics", level=1))
        ets = [ [ "stage", "run id", "sources", "minimum", 
                "delta (eV)", "maximum", "delta (eV)" ] ]
        ets[0] = [Paragraph(i, self.stylesheet['celltitle']) for i in ets[0]]
        style = copy.deepcopy(self.tablestyles["mainxx"])
        for d in self.d.d:
            for g in d.groups:
                xeners, xfinals = {}, {}
                gm = None
                for s in g.stages:
                    for r in s.runs:
                        if r.sour_name == "unknown": continue
                        xepre = np.array([ x.props["initial_energy"] for x in r.optima.values() 
                            if "initial_energy" in x.props ])
                        if len(xepre) == 0: continue
                        xeners[r.pname] = xepre
                    for r in s.runs:
                        if r.sour_name == "unknown": continue
                        xepre = np.array([ x.props["energy"] for x in r.optima.values() 
                            if "energy" in x.props and x.props["energy"] is not None ])
                        if len(xepre) == 0: continue
                        xfinals[r.pname] = xepre
                        if gm is None or xfinals[r.pname].min() < gm:
                            gm = xfinals[r.pname].min()
                for s in g.stages:
                    for r in s.runs:
                        if r.pname not in xeners: continue
                        ets += [ [ "%s.%s initial" % (d.name, r.parent.parent.name), r.pname, 
                            r.sour_name, "%12.6f" % xeners[r.pname].min(), 
                            "%8.3f" % ((xeners[r.pname].min() - gm) * httoev), 
                            "%12.6f" % xeners[r.pname].max(), 
                            "%8.3f" % ((xeners[r.pname].max() - gm) * httoev) ] ]
                line = len(ets) - 1
                style._cmds += [ ('LINEBELOW',    (0, line), (-1, line), 1, "#C5C5C5"), ]
                for s in g.stages:
                    for r in s.runs:
                        if r.pname not in xfinals: continue
                        ets += [ [ "%s.%s final" % (d.name, r.parent.parent.name), r.pname,    
                            r.sour_name, "%12.6f" % xfinals[r.pname].min(), 
                            "%8.3f" % ((xfinals[r.pname].min() - gm) * httoev), 
                            "%12.6f" % xfinals[r.pname].max(), 
                            "%8.3f" % ((xfinals[r.pname].max() - gm) * httoev) ] ]
                line = len(ets) - 1
                style._cmds += [ ('LINEBELOW',    (0, line), (-1, line), 1, "#C5C5C5"), ]
            line = len(ets) - 1
            style._cmds += [ ('LINEBELOW',    (0, line), (-1, line), 1, "#A5A5A5"), ]
        story.append(self.draw_simple_table(ets, style=style,
            widths=(65, 40, 80, 75, 50, 75, 50), span=range(0, 3), noacross=0))
        story = story[:-2] + [ KeepTogether(story[-2:]) ]

        # Configuration Convergence
        sub_title = self.heading('Configuration Convergence', key="convergence", level=1)
        sub_note = Paragraph("""<font face="Helvetica">Thick lines:</font> 
             converged structures with energy up to %.2f eV 
            with respect to global minimum of each multiplicity. <br/>
            <font face="Helvetica">Thin lines:</font> all converged structrues. """ 
            % self.d.conv_cut, self.stylesheet['body'])
        sub_title = KeepTogether([sub_title, sub_note])
        idx = 0
        for dx in self.d.d:
            for g in dx.groups:
                for s in g.stages:
                    if s.conv_maxs == [0, 0]: continue
                    sub_story = []
                    if idx == 0:
                        sub_story = [ sub_title ]
                        idx = 1
                    sname = "%s.%s.%s" % (dx.name, g.name, s.stage)
                    sub_story.append(self.heading(s.description, cont="Stage %s" % sname, 
                        key="conv-sta-%s" % sname, level=2))
                    sub_story.append(PGPlot.conv_plot(s.conv_data, s.conv_maxs))
                    story.append(KeepTogether(sub_story))

        # Energy and Step Distribution
        sub_title = self.heading('Energy and Step Distribution', key="distri", level=1)
        idx = 0
        for dx in self.d.d:
            for g in dx.groups:
                if g.xsour_structi.keys() == [ "final" ] and len(g.xsour_structi["final"]) == 0 and \
                    len(g.locals) == 0: continue
                sub_story = []
                if idx == 0:
                    sub_story = [ sub_title ]
                    idx = 1
                gname = "%s.%s" % (dx.name, g.name)
                sub_story.append(self.heading(g.description, cont="Group %s" % gname, 
                    key="distr-gp-%s" % gname, level=2))
                sub_story.append(PGPlot.ener_plot(g.xsour_structi, g.locals))
                story.append(KeepTogether(sub_story))
            if len(dx.groups) > 1:
                sub_story = []
                stt = {}
                for g in dx.groups: stt.update(g.xsour_structi)
                sub_story.append(self.heading("[Total]", cont="File %s" % dx.name, 
                    key="distr-total-%s" % dx.name, level=2))
                sub_story.append(PGPlot.ener_plot(stt, []))
                story.append(KeepTogether(sub_story))
        
        # Monte-Carlo Energetics
        sub_title = self.heading('Monte-Carlo Energetics', key="mcener", level=1)
        idx = 0
        for dx in self.d.d:
            for r in dx.runs:
                if r.mc is None: continue
                sub_story = []
                if idx == 0:
                    sub_story = [ sub_title ]
                    idx = 1
                rname = "%s.%s.%s.%d" % (dx.name, r.parent.parent.name, r.stage, r.m)
                sub_story.append(self.heading(r.description, cont="Run %s" % rname, 
                    key="mcener-run-%s" % rname, level=2))
                story.append(KeepTogether(sub_story))
                yvalues = np.array([xx[1] for x in r.mce.values() for xx in x])
                if len(yvalues) == 0:
                    continue
                tmaxy, tminy = yvalues.max(), yvalues.min()
                for imck, imc in sorted(r.mce.items(), key=lambda x: x[0]):
                    story.append(PGPlot.mcener_plot(imc, int(imck), len(r.mce), miny=tminy, maxy=tmaxy))

        # Probabilities
        story += self.probabilities('Boltzmann Distribution Probabilities', "boltz",
            vib=False, rot=False, max_prop_count=max_prop_count)
        story += self.probabilities('Vibrational Entropy Included Probabilities', "vibentro", 
            vib=True, rot=False, max_prop_count=max_prop_count)
        story += self.probabilities('Vib-Rot Entropy Included Probabilities', "vrentro", 
            vib=True, rot=True, max_prop_count=max_prop_count)

        # heat capacities at 200K
        story += self.heat_capacity('Electronic Heat Capacity (0-200 K)', "hc", vib=False, rot=False, tmax=200.0)
        story += self.heat_capacity('Elec-Vib Heat Capacity (0-200 K)', "vibhc", vib=True, rot=False, tmax=200.0)
        story += self.heat_capacity('Elec-Vib-Rot Heat Capacity (0-200 K)', "vrhc", vib=True, rot=True, tmax=200.0)

        # heat capacities at 1000K
        story += self.heat_capacity('Electronic Heat Capacity (0-1000 K)', "hct", vib=False, rot=False, tmax=1000.0)
        story += self.heat_capacity('Elec-Vib Heat Capacity (0-1000 K)', "vibhct", vib=True, rot=False, tmax=1000.0)
        story += self.heat_capacity('Vib-Rot Heat Capacity (0-1000 K)', "vrhct", vib=True, rot=True, tmax=1000.0)

        # VIP at 1000K
        story += self.prop_graph('Electronic Vertical Ionization Potential', "vip", "vip.energy", 
            "VIP Energy (eV)", vib=False, rot=False, tmax=1000.0)
        story += self.prop_graph('Elec-Vib Vertical Ionization Potential', "vibvip", "vip.energy", 
            "VIP Energy (eV)", vib=True, rot=False, tmax=1000.0)
        story += self.prop_graph('Elec-Vib-Rot Vertical Ionization Potential', "vrvip", "vip.energy", 
            "VIP Energy (eV)", vib=True, rot=True, tmax=1000.0)

        # Details
        sub_title = self.heading('Details', key="details", level=1)
        if self.energy_cutoff_hard is not None and self.energy_cutoff_hard > 0.0:
            sub_note = Paragraph("""<font face="Helvetica">Note: </font> 
                Only structures with energy up to %.2f eV 
                with respect to global minimum are listed.""" % self.energy_cutoff_hard, 
                self.stylesheet['body'])
            sub_title = KeepTogether([sub_title, sub_note])
        idx = 0
        for dx in self.d.d:
            for g in dx.groups:
                if len(g.finalsx) == 0: continue
                sub_story = []
                if idx == 0:
                    sub_story = [ sub_title ]
                    idx = 1
                gname = "%s.%s" % (dx.name, g.name)
                sub_story.append(self.heading(g.description, cont="Group %s" % gname, 
                    key="detail-gp-%s" % gname, level=2))
                if hasattr(g.finalsx[0], "surf"):
                    dc = self.draw_cluster_at_surface(g.finalsx)
                else:
                    dc = self.draw_cluster(g.finalsx)
                if dc is not None: sub_story.append(dc)
                story.append(KeepTogether(sub_story))
        
        if not self.no_path:
            print ('building paths ...')
            # Reaction Paths
            sub_title = self.heading('Reaction Paths', key="rpaths", level=1)
            idx = 0
            for dx in self.d.d:
                for g in dx.groups:
                    if len(g.nebpaths) == 0:
                        continue
                    sub_story = []
                    if idx == 0:
                        sub_story = [sub_title]
                        idx = 1
                    gname = "%s.%s" % (dx.name, g.name)
                    sub_story.append(self.heading(g.description, cont="Group %s" % gname, 
                        key="rpaths-gp-%s" % gname, level=2))
                    if g.nebpaths.values()[0].surf is not None:
                        dc = self.draw_path_at_surface(g.nebpaths, surface_depth)
                    else:
                        dc = None
                    if dc is not None:
                        sub_story.append(dc)
                    story.append(KeepTogether(sub_story))

        sub_title = self.heading('Reaction Path Minima', key="rpmins", level=1)
        idx = 0
        for dx in self.d.d:
            for g in dx.groups:
                if len(g.nebmins) == 0:
                    continue
                sub_story = []
                if idx == 0:
                    sub_story = [sub_title]
                    idx = 1
                gname = "%s.%s" % (dx.name, g.name)
                sub_story.append(self.heading(g.description, cont="Group %s" % gname, 
                    key="rpmins-gp-%s" % gname, level=2))
                if g.nebmins.values()[0].surf is not None:
                    dc = self.draw_path_min_at_surface(g.nebmins, surface_depth)
                else:
                    dc = None
                if dc is not None:
                    sub_story.append(dc)
                if g.nebgraph is not None:
                    dc = DrawRGraph(g.nebgraph)
                    sub_story.append(dc)
                story.append(KeepTogether(sub_story))

        self.doc.multiBuild(story, canvasmaker=NumberedCanvas)
