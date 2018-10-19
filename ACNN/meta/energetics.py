
from reportlab.platypus.flowables import Flowable
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle
from meta.projection import DrawCluster
import numpy as np
from reportlab.lib.enums import TA_CENTER
from cluster.base import Cluster

httoev = 27.21138505
httok = 315774.646

def multi_name(mi):
    listm = ["mixed", "singlet", "doublet", "triplet", "quartet", "quintet", 
        "sextet", "septet", "octet", "nonet", "dectet", 
        "undectet", "duodectet", "tredectet", "quatuordectet", "quindectet", 
        "sedectet", "septendectet", "decennoctet", "decennovtet", "vigetettet"]
    if mi >= 0 and mi <= 20:
        return listm[mi]
    else:
        return str(mi) + '-tet'

def multi_color(mi):
    listm = ["#236A45", "#236A45", "#AC5A23", "#A954CE", "#AA528C", "#5477E5", 
        "#77AE42", "#EC8659", "#44B883", "#764161", "#AC62D2", 
        "#E24AA9", "#54B3B1", "#7A7C97", "#42DDEC", "#C466C7", 
        "#9A6BA4", "#6C3ABA", "#379AC2", "#D4825A", "#82DD72"]
    if mi >=0 and mi <= 20:
        return listm[mi]
    else:
        return "#111111"

def mini_name(mi):
    x = ["UNK", "MIN", "TS"]
    return x[mi]

def draw_flowable(flow, canvas, x, y, valign='TOP'):
    hh = flow.canv if hasattr(flow, "canv") else None
    canvas.saveState()
    flow.canv = canvas
    if valign == 'TOP':
        flow.canv.translate(x, y - flow.height)
    else:
        flow.canv.translate(x, y)
    flow.draw()
    canvas.restoreState()
    flow.canv = hh

class BigCanvas(Flowable, object):
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.items = []
        self.coords = []

    def add(self, draw, x, y):
        self.items.append(draw)
        self.coords.append((x, y))

    def wrap(self, *opts):
        return (self.width, self.height)

    def draw(self):
        for i, (x, y) in zip(self.items, self.coords):
            if x is not None:
                draw_flowable(i, self.canv, x, y, valign=None)
            else:
                i(self.canv)

class GraphItem(object):
    def __init__(self, width, height, line_pos, idx):
        self.width = width
        self.height = height
        self.selected = -1
        self.line_pos = line_pos
        self.idx = idx
        self.put_pos = None

    def get_diff(self):
        return np.abs(self.put_pos - self.line_pos)

    @staticmethod
    def select(items, width):
        xpos = [0.0, 0.0]
        while True:
            fmin = [[], []]
            for k in [0, 1]:
                for i in items:
                    if i.selected == -1:
                        if xpos[k] + i.width > width: continue
                        i.put_pos = xpos[k] + i.width / 2
                        dd = i.get_diff()
                        if len(fmin[k]) == 0 or dd < fmin[k][1]:
                            fmin[k] = [i.idx, dd]
                if len(fmin[k]) != 0:
                    items[fmin[k][0]].selected = k
                    xpos[k] += items[fmin[k][0]].width
            if len(fmin[0]) == 0 and len(fmin[1]) == 0:
                break
        ksel = [[], []]
        khei = [0, 0]
        for k in [0, 1]:
            for i in items:
                if i.selected == k:
                    ksel[k].append(i)
            xpos = 0.0
            for i in ksel[k]:
                i.put_pos = xpos + i.width / 2
                xpos += i.width
            if len(ksel[k]) == 0: continue
            if ksel[k][-1].put_pos < width - ksel[k][-1].width / 2:
                st = ksel[k][0].width / 2
                ratio = width - ksel[k][-1].width / 2 - st
                if ksel[k][-1].put_pos != st:
                    ratio = ratio / (ksel[k][-1].put_pos - st)
                for i in ksel[k]:
                    i.put_pos = (i.put_pos - st) * ratio + st
            ma = None
            for ii, i in enumerate(ksel[k]):
                i.left = i.right = False
                if ii == 0: i.left = True
                if ii == len(ksel[k]) - 1: i.right = True
                i.idxk = ii
                if ma is None or i.get_diff() < ma:
                    mai = ii
                    ma = i.get_diff()
            for i in ksel[k][:mai][::-1] + ksel[k][mai:]:
                if i.left: ltouch = 0
                else:
                    nl = ksel[k][i.idxk - 1]
                    ltouch = nl.put_pos + nl.width / 2
                if i.right: rtouch = width
                else:
                    nr = ksel[k][i.idxk + 1]
                    rtouch = nr.put_pos - nr.width / 2
                if i.put_pos < i.line_pos and rtouch > i.put_pos + i.width / 2 + 1:
                    i.put_pos += min(i.line_pos - i.put_pos, rtouch - (i.put_pos + i.width / 2))
                elif i.put_pos > i.line_pos and ltouch < i.put_pos - i.width / 2 - 1:
                    i.put_pos -= min(i.put_pos - i.line_pos, i.put_pos - i.width / 2 - ltouch)
            khei[k] = max([i.height for i in ksel[k]])
            ksel[k] = [i.idx for i in ksel[k]]
        return ksel[0], ksel[1], khei[0], khei[1]

class DrawEnergetics(Flowable, object):
    def __init__(self, clus, width=500, cutoff=0.2, ratio=2.7, no_freq=False, surface_depth=3.0):
        self.canv = None
        super(DrawEnergetics, self).__init__()
        self.width = width
        self.clus = clus
        self.no_freq = no_freq
        self.min_e = min_e = self.clus[0].energy
        self.clus = [c for c in self.clus if c.energy < (cutoff / httoev) + min_e]
        self.es = [(c.energy - min_e) * httoev for c in self.clus]
        self.maxe = cutoff
        self.scale = 10
        self.fontsize = 7
        self.ratio = ratio
        self.graphs = [ DrawCluster(c, simple=True, ratio=self.ratio) if isinstance(c, Cluster)
            else DrawCluster(c, simple=True, ratio=self.ratio * 1.25, rotmat=[0.0, 1.0, 0.0, 0.0], 
            surface_depth=surface_depth) for c in self.clus ]
        self.margin_ratio = 0.05
        self.space = 5
        self.vspace = 5
        self.pstyle = ParagraphStyle(name="body", fontName="Arial",
            fontSize=self.fontsize, leading=self.fontsize , textColor="#323232", alignment=TA_CENTER)
        for i, g in enumerate(self.graphs):
            ener_xpos = self.es[i] / self.maxe * (self.width * (1 - self.margin_ratio * 2)) + \
                self.width * self.margin_ratio
            g.dpos = ener_xpos
        self.graph_items = [ GraphItem(g.width + self.space * 2, g.height, g.dpos, ig) 
            for ig, g in enumerate(self.graphs) ]
        self.measure()
    
    def wrap(self, *opts):
        return (self.width, self.height)
    
    def draw_connection(self, i, heis):
        self.canv.setStrokeColor(multi_color(self.clus[i.idx].multiplicity))
        self.canv.setFillColor(multi_color(self.clus[i.idx].multiplicity))
        new_xpos = i.put_pos
        ener_xpos = i.line_pos
        self.canv.line(ener_xpos, heis[0], ener_xpos, heis[1])
        self.canv.line(ener_xpos, heis[1], new_xpos, heis[2])
        self.canv.circle(new_xpos, heis[2], 0.8, stroke=0, fill=1)
        self.canv.circle(ener_xpos, heis[0], 0.8, stroke=0, fill=1)
    
    def measure(self):
        self.upp, self.dwp, maxha, maxhb = GraphItem.select(self.graph_items, self.width)
        self.height = 15 + self.fontsize + 15 + self.vspace * 4 + 45 + maxha + maxhb \
             + self.fontsize * 4
        self.hfheight = 15 + self.fontsize + 15 + self.vspace * 2 + maxhb \
             + self.fontsize * 2

    def draw(self):
        if [i.minimum_type for i in self.clus].count(0) == len(self.clus):
            self.no_freq = True
        self.canv.setStrokeColor("#C8C8C8")
        self.canv.setLineWidth(0.25)
        # if not self.no_freq:
        self.canv.rect(0, 0, self.width, self.height, fill=0, stroke=1)
        self.canv.setFillColor("#E2E2E2")
        self.canv.rect(self.width * self.margin_ratio, self.hfheight - 5, 
            self.width * (1 - self.margin_ratio * 2), 5, fill=1, stroke=0)
        self.canv.setLineWidth(0.5)
        self.canv.setFont("Arial", self.fontsize)
        self.canv.setStrokeColor("#575757")
        self.canv.setFillColor("#575757")
        for i in range(self.scale + 1):
            text = "%.2f eV" % (i * self.maxe / self.scale)
            xpos = float(i) / self.scale * (self.width * (1 - self.margin_ratio * 2)) + \
                self.width * self.margin_ratio
            self.canv.line(xpos, self.hfheight, xpos, self.hfheight - 5 - self.fontsize)
            if i != self.scale:
                self.canv.drawString(xpos + 2, self.hfheight - 5 - self.fontsize, text)
        for ii, i in enumerate(self.clus):
            self.canv.setStrokeColor(multi_color(i.multiplicity))
            xpos = self.es[ii] / self.maxe * (self.width * (1 - self.margin_ratio * 2)) + \
                self.width * self.margin_ratio
            self.canv.line(xpos, self.hfheight - 5, xpos, self.hfheight + 18)
        for i in self.dwp:
            heis = [
                self.hfheight - 5 - self.fontsize - 2 - self.vspace, 
                self.hfheight - 5 - self.fontsize - 10 - self.vspace, 
                self.hfheight - 5 - self.fontsize - 15 - self.vspace
            ]
            gpos = self.graph_items[i].put_pos - self.graphs[i].width / 2
            draw_flowable(self.graphs[i], self.canv, gpos, heis[2])
            self.draw_connection(self.graph_items[i], heis)
            if self.no_freq:
                text1 = "[%d] (%s)" % (i + 1, multi_name(self.clus[i].multiplicity))
            else:
                text1 = "[%d] (%s)<sup rise=2>%d</sup>" % (i + 1, 
                    mini_name(self.clus[i].minimum_type), self.clus[i].multiplicity)
            text2 = "%.2f eV" % self.es[i]
            cpt = self.clus[i].pointgroup()
            if cpt != "C1":
                text2 += " (%s<sub rise=2>%s</sub>)" % (cpt[0], cpt[1:])
            self.pstyle.textColor = multi_color(self.clus[i].multiplicity)
            if self.clus[i].aname is not None:
                para = Paragraph("<a href=\"#X%s\">%s</a>" % (self.clus[i].aname, 
                    text1), self.pstyle)
            else:
                para = Paragraph("%s" % text1, self.pstyle)
            para.wrap(self.graphs[i].width * 2, self.fontsize * 2)
            para.drawOn(self.canv, gpos - self.graphs[i].width / 2, 
                heis[2] - self.graphs[i].height - self.fontsize)
            para2 = Paragraph(text2, self.pstyle)
            para2.wrap(self.graphs[i].width * 2, self.fontsize * 2)
            para2.drawOn(self.canv, gpos - self.graphs[i].width / 2, 
                heis[2] - self.graphs[i].height - self.fontsize * 2)
        for i in self.upp:
            heis = [ self.hfheight + 20 + self.vspace, 
                self.hfheight + 28 + self.vspace, self.hfheight + 35 + self.vspace ]
            gpos = self.graph_items[i].put_pos - self.graphs[i].width / 2
            draw_flowable(self.graphs[i], self.canv, gpos, heis[2], valign='BOTTOM')
            self.draw_connection(self.graph_items[i], heis)
            if self.no_freq:
                text1 = "[%d] (%s)" % (i + 1, multi_name(self.clus[i].multiplicity))
            else:
                text1 = "[%d] (%s)<sup rise=2>%d</sup>" % (i + 1, 
                    mini_name(self.clus[i].minimum_type), self.clus[i].multiplicity)
            text2 = "%.2f eV" % self.es[i]
            cpt = self.clus[i].pointgroup()
            if cpt != "C1":
                text2 += " (%s<sub rise=2>%s</sub>)" % (cpt[0], cpt[1:])
            self.pstyle.textColor = multi_color(self.clus[i].multiplicity)
            if self.clus[i].aname is not None:
                para = Paragraph("<a href=\"#X%s\">%s</a>" % (self.clus[i].aname, 
                    text1), self.pstyle)
            else:
                para = Paragraph("%s" % text1, self.pstyle)
            para.wrap(self.graphs[i].width * 2, self.fontsize * 2)
            para.drawOn(self.canv, gpos - self.graphs[i].width / 2, 
                heis[2] + self.graphs[i].height + self.fontsize * 2 - 3)
            para2 = Paragraph(text2, self.pstyle)
            para2.wrap(self.graphs[i].width * 2, self.fontsize * 2)
            para2.drawOn(self.canv, gpos - self.graphs[i].width / 2, 
                heis[2] + self.graphs[i].height + self.fontsize - 3)
