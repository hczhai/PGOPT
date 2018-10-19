# -*- coding: utf-8 -*-

from reportlab.graphics.shapes import Drawing, Line, String, Group, Rect
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.widgets.markers import makeMarker
from meta.energetics import httoev, multi_name, multi_color, BigCanvas
from reportlab.lib.colors import HexColor
from reportlab.lib.attrmap import AttrMap, AttrMapValue
from reportlab.pdfbase.pdfmetrics import stringWidth
import numpy as np

def source_color(sour):
    listm = ["#236A45", "#57AA22", "#A954CE", "#AA528C", "#5477E5", 
        "#77AE42", "#EC8659", "#44B883", "#764161", "#CC4242", 
        "#A24AB9", "#54B3B1", "#7A7C97", "#42DDEC", "#C466C7", 
        "#9A6BA4", "#6C3ABA", "#379AC2", "#D4825A", "#82DD72"]
    listx = [ "local", "structs", "trans", "create/rv", 
        "create/blda", "create/ck", "final", "create/bc", 
        "step", "nn", "partial" ]
    if sour in listx:
        return listm[listx.index(sour)]
    else:
        return "#111111"

class AddedPlot(LinePlot, object):

    _attrMap = AttrMap(BASE=LinePlot, g = AttrMapValue(None, desc='Groups.'))

    def __init__(self):
        super(AddedPlot, self).__init__()
        self.g = Group()
    
    def draw(self):
        g = super(AddedPlot, self).draw()
        for i in self.g.contents:
            if isinstance(i, Rect):
                x = self.xValueAxis.scale(i.x)
                width = self.xValueAxis.scale(i.width) - self.xValueAxis.scale(0)
                y = self.yValueAxis.scale(i.y)
                height = self.yValueAxis.scale(i.height) - self.yValueAxis.scale(0)
                g.add(Rect(x, y, width, height, strokeColor=i.strokeColor, 
                    strokeWidth=i.strokeWidth, fillColor=i.fillColor))
            else:
                g.add(i)
        return g

class PGPlot(object):

    @staticmethod
    def blank_plot(maxx, maxy, minx=0.0, miny=0.0):
        if maxx == 0.0: maxx = 1.0
        if maxy < 1E-8 and miny == 0.0: maxy = 1.0
        lpt = AddedPlot()
        lpt.strokeColor = HexColor("#454545")
        lpt.joinedLines = 1
        lpt.xValueAxis.labels.fontName = 'Arial'
        lpt.xValueAxis.tickDown = 0
        lpt.xValueAxis.tickUp = 5
        lpt.yValueAxis.tickLeft = 0
        lpt.yValueAxis.tickRight = 5
        lpt.yValueAxis.labels.fontName = 'Arial'
        lpt.xValueAxis.valueMin = minx
        lpt.xValueAxis.valueMax = maxx
        lpt.yValueAxis.valueMin = miny
        lpt.yValueAxis.valueMax = maxy
        lpt.data = [ [(0, 0), (0, 0)], [(0, 0), (0, 0)] ]
        rto = 7
        if maxx - minx > 10 or isinstance(minx, int):
            lpt.xValueAxis.labelTextFormat = '%.0f'
        elif maxx - minx > 1: lpt.xValueAxis.labelTextFormat = '%.1f'
        elif maxx - minx > 0.1: lpt.xValueAxis.labelTextFormat = '%.2f'
        else: lpt.xValueAxis.labelTextFormat = '%.3f'
        if maxy - miny > 10: lpt.yValueAxis.labelTextFormat = '%.0f'
        elif maxy - miny > 1: lpt.yValueAxis.labelTextFormat = '%.1f'
        elif maxy - miny > 0.1: lpt.yValueAxis.labelTextFormat = '%.2f'
        else: lpt.yValueAxis.labelTextFormat = '%.3f'
        xsp = float(lpt.xValueAxis.labelTextFormat % ((maxx - minx) / rto))
        ysp = float(lpt.yValueAxis.labelTextFormat % ((maxy - miny) / rto))
        if xsp > 10: xsp = int(xsp) / 10 * 10
        if ysp > 10: ysp = int(ysp) / 10 * 10
        if xsp == 0.0:
            xsp = 1.0
        if ysp == 0.0:
            ysp = 1.0
        if np.abs(miny - maxy) < 1E-10:
            maxy = miny + 1.0
            ysp = 1.0
        lpt.xValueAxis.valueSteps = list(np.arange(minx, maxx + xsp / 2, xsp))
        lpt.yValueAxis.valueSteps = list(np.arange(miny, maxy + ysp / 2, ysp))
        if lpt.xValueAxis.valueSteps[-1] > maxx:
            lpt.xValueAxis.valueMax = lpt.xValueAxis.valueSteps[-1]
        if lpt.yValueAxis.valueSteps[-1] > maxy:
            lpt.yValueAxis.valueMax = lpt.yValueAxis.valueSteps[-1]
        return lpt

    @staticmethod
    def line_plot(data, colors, widths, maxx, maxy, minx=0.0, miny=0.0,
                  marker=False, intpol=False):
        ym, yn = maxy, miny
        if not intpol:
            lpt = PGPlot.blank_plot(maxx, maxy, minx=minx, miny=miny)
            if len(data) != 0 and len(data[0]) != 0:
                lpt.data = data
            mm = makeMarker('FilledCircle', size=2)
            for ii, i in enumerate(zip(colors, widths)):
                lpt.lines[ii].strokeColor = i[0]
                lpt.lines[ii].strokeWidth = i[1]
                if marker:
                    lpt.lines[ii].symbol = mm
        else:
            from scipy.interpolate import CubicSpline
            ldata = []
            for d in data:
                dd = np.array(d)
                ipx = CubicSpline(dd[:, 0], dd[:, 1], bc_type='clamped')
                datax = np.arange(minx, max(dd[:, 0]), (maxx - minx) / 500.0 - 1E-6)
                datay = ipx(datax)
                ldata += [zip(datax, datay), d]
                ym = max(max(datay), ym)
                yn = min(min(datay), yn)
            lpt = PGPlot.blank_plot(maxx, ym, minx=minx, miny=yn)
            lpt.data = ldata
            mm = makeMarker('FilledCircle', size=3, fillColor=HexColor("#FFFFFF"),
                            strokeColor=colors[0], strokeWidth=widths[0])
            for ii, i in enumerate(zip(colors, widths)):
                lpt.lines[ii * 2].strokeColor = i[0]
                lpt.lines[ii * 2].strokeWidth = i[1]
                lpt.lines[ii * 2 + 1].strokeColor = None
                lpt.lines[ii * 2 + 1].strokeWidth = 0.0
                if marker:
                    lpt.lines[ii * 2 + 1].symbol = mm
        return lpt, ym, yn

    @staticmethod
    def mcener_plot(imc, imck, imct, miny, maxy):
        xlabel = "Step Index"
        ylabel = "Accepted Relaxed Energy (eV)"
        data = [(x[0], x[1] * httoev) for x in imc]
        lpt, _, _ = PGPlot.line_plot([data], [HexColor(source_color("trans"))], [1.0],
                               max([d[0] for d in data]), maxy * httoev,
                               0.0, miny * httoev, marker=True)
        drw = Drawing(width=650, height=230)
        lpt.x, lpt.y = 20, 15
        lpt.width, lpt.height = 450, 180
        drw.add(lpt)
        drw.add(String(480, 2, xlabel, fontName="Arial", fontSize=12,
                       fillColor=HexColor("#151515")))
        drw.add(String(2, drw.height - 25, ylabel,
                       fontName="Arial", fontSize=12, fillColor=HexColor("#151515")))
        drw.add(String(380, drw.height - 48, "Walker %d of %d" % (imck + 1, imct),
                       fontName="Arial", fontSize=12,
                       fillColor=HexColor(source_color("trans"))))
        return drw

    @staticmethod
    def path_plot(ipe, cs, diffs, surface_depth, ratio):
        xlabel = "Image Index"
        ylabel = "Energy (eV)"
        data = np.array([(ix, (x - ipe[0]) * httoev) for ix, x in enumerate(ipe)])
        mimages = int(max([d[0] for d in data]))
        lpt, ym, yn = PGPlot.line_plot([data], [HexColor(source_color("local"))], [1.0],
                               mimages, data[:, 1].max(),
                               0, data[:, 1].min(), marker=True, intpol=True)
        ym, yn = lpt.yValueAxis.valueMax, lpt.yValueAxis.valueMin
        from cluster.base import Cluster
        from meta.projection import DrawCluster
        graphs = []
        for c in cs:
            if isinstance(c, Cluster):
                graphs.append(DrawCluster(c, simple=True, ratio=ratio))
            else:
                graphs.append(DrawCluster(c, simple=True, ratio=ratio * 1.25,
                                          rotmat=[0.0, 1.0, 0.0, 0.0],
                                          surface_depth=surface_depth))
        grh = irh = 190
        irw = 500
        ilw = 80
        grleft = []
        igr = 0
        igrh = 0.0
        while igr < len(graphs):
            grleft.append(graphs[igr])
            igrh += graphs[igr].height
            igr += 1
            if igrh > grh:
                break
        gwid, ghei = max([g.width for g in grleft]), sum([g.height for g in grleft])
        grtop = []
        igrw = 0.0
        irdw = 10.0
        while igr < len(graphs):
            grtop.append(graphs[igr])
            igrw += graphs[igr].width + irdw
            igr += 1
        bhei = 0.0 # height below plot
        blw = 1
        if len(grtop) > 6:
            bhei = max([gb.height for gb in grtop[1:-1:2]])
            blw = 2
        gthei = max([g.height for g in grtop[0:-1:blw]]) + grh
        ghei = max(ghei, irh, gthei)
        drw = Drawing(width=irw, height=irh)
        lpt.x, lpt.y = 35, 20
        lpt.width, lpt.height = 230, 150
        drw.add(lpt)
        drw.add(String(270, 12, xlabel, fontName="Arial", fontSize=12,
                       fillColor=HexColor("#151515")))
        drw.add(String(0, irh - 10, ylabel, fontName="Arial", fontSize=12,
                       fillColor=HexColor("#151515")))
        for i in range(0, mimages):
            iiy = lpt.y + lpt.height * ((data[i, 1] + data[i + 1, 1]) / 2 - yn) / (ym - yn)
            iix = lpt.x + lpt.width * (i + 0.5) / mimages
            drw.add(String(iix, iiy, "%.2f" % diffs[i + 1], fontName="Arial",
                           fontSize=8, fillColor=HexColor("#4A4A4A"), textAnchor="middle"))
        bc = BigCanvas(width=irw + gwid + ilw, height=ghei + bhei)
        bc.add(drw, gwid + ilw, bhei)
        blines = []
        igrh = 0.0
        for ig, g in enumerate(grleft):
            bc.add(g, gwid - g.width + ilw, igrh + bhei)
            blines.append([gwid + ilw, igrh + g.height / 2 + bhei,
                           gwid + ilw + lpt.width * ig / mimages + lpt.x,
                           lpt.y + lpt.height * (data[ig, 1] - yn) / (ym - yn) + bhei])
            igrh += g.height
        ign = len(grleft)
        igrw = ilw + gwid
        for ig, g in enumerate(grtop[0:-1:blw]):
            bc.add(g, igrw, grh + bhei)
            blines.append([igrw + g.width / 2, grh + bhei,
                           gwid + ilw + lpt.width * (ig * blw + ign) / mimages + lpt.x,
                           lpt.y + lpt.height * (data[ig * blw + ign, 1] - yn) / (ym - yn) + bhei])
            igrw += g.width + irdw
        if len(grtop) != 0:
            hh = drw.height - grtop[-1].height
            if hh < 0:
                hh = 0
            bc.add(grtop[-1], ilw + gwid + lpt.width + irdw + lpt.x, hh + bhei)
            blines.append([ilw + gwid + lpt.width + irdw + lpt.x, hh + grtop[-1].height / 2 + bhei,
                           gwid + ilw + lpt.width + lpt.x,
                           lpt.y + lpt.height * (data[-1, 1] - yn) / (ym - yn) + bhei])
        if blw == 2:
            igrw = ilw + gwid
            for ig, g in enumerate(grtop[1:-1:blw]):
                bc.add(g, igrw, 0.0)
                blines.append([igrw + g.width / 2, g.height,
                               gwid + ilw + lpt.width * (ig * blw + 1 + ign) / mimages + lpt.x,
                               lpt.y + lpt.height * \
                               (data[ig * blw + 1 + ign, 1] - yn) / (ym - yn) + bhei])
                igrw += g.width + irdw

        def dd(canv, blx):
            canv.setStrokeColor("#B1B1B1")
            canv.setLineWidth(0.5)
            canv.line(*blx)

        for bl in blines:
            bc.add(lambda x, bl=bl: dd(x, bl), None, None)
        return bc

    @staticmethod
    def ener_plot(xcour_str, xlocals, ecut=15.0):
        lmin = None
        data, colors, names = [], [], []
        for k, v in xcour_str.items():
            nv = np.array([c.energy for c in v if c.energy is not None]) * httoev
            nmin = nv.min() if nv.shape[0] != 0 else 0.0
            if lmin is None or lmin > nmin:
                lmin = nmin
            data.append(nv)
            if len(nv) != 0:
                colors.append(HexColor(source_color(k)))
        for iv, nv in enumerate(data):
            nv = nv - lmin
            data[iv] = [x for x in nv if x < ecut]
        for ix, (k, v) in enumerate(xcour_str.items()):
            nv = np.array(data[ix])
            if len(nv) != 0:
                names.append([u"%s (%.1f ± %.1f)" % (k, nv.mean(), nv.std()), 
                    HexColor(source_color(k))])
        data = [x for x in data if len(x) != 0]
        nv = np.array([c.props["step"] for c in xlocals], dtype=float)
        datas = [nv]
        lpt = PGPlot.histogram(data, colors)
        drw = Drawing(width=500, height=200)
        lpt.x, lpt.y = 25, 30
        lpt.width, lpt.height = 230, 150
        drw.add(lpt)
        if len(xlocals) != 0:
            lpt2 = PGPlot.histogram(datas, [ HexColor(source_color('local')) ])
            lpt2.x, lpt2.y = 300, 30
            lpt2.width, lpt2.height = 230, 150
            drw.add(lpt2)
        sh = drw.height - 30
        dsh = 10
        sxx = 0.0
        sxr = [stringWidth(na, fontName="Arial", fontSize=10) for na, cl in names]
        for ix, (na, cl) in enumerate(names):
            drw.add(Line(lpt.x + 202 - sxr[ix], sh, lpt.x + 222 - sxr[ix], sh, 
                strokeColor=cl, strokeWidth=2.0))
            drw.add(String(lpt.x + 225 - sxr[ix], sh - 3, na, 
                fontName="Arial", fontSize=10, fillColor=cl))
            sh -= dsh
        if len(xlocals) != 0:
            smean = datas[0].mean()
            sstd = datas[0].std()
            snames = [ u"%.0f ± %.0f" % (smean, sstd) ]
            sc = [stringWidth(x, fontName="Arial", fontSize=10) for x in snames]
            sh = drw.height - 30
            for x, xsc in zip(snames, sc):
                drw.add(String(lpt2.x + 225 - xsc, sh - 3, x, 
                    fontName="Arial", fontSize=10, 
                    fillColor=HexColor(source_color('local'))))
                sh -= dsh
        drw.add(String(200, 0, "Energy (eV)", fontName="Arial", fontSize=12, 
            fillColor=HexColor("#151515")))
        drw.add(String(0, drw.height - 10, "Probability", fontName="Arial", fontSize=12, 
            fillColor=HexColor("#151515")))
        if len(xlocals) != 0:
            drw.add(String(482, 0, "# of Steps", fontName="Arial", fontSize=12, 
                fillColor=HexColor("#151515")))
            drw.add(String(275, drw.height - 10, "Probability", fontName="Arial", fontSize=12, 
                fillColor=HexColor("#151515")))
        return drw

    @staticmethod
    def histogram(data, colors, alpha=0.6, bins=40):
        hist = []
        maxx, maxy = None, None
        for d in data:
            dx = np.array(d)
            dmax = dx.max() if dx.shape[0] != 0 else 0.0
            if maxx is None or dmax > maxx: maxx = dmax
        for d in data:
            dx = np.array(d)
            cnts = np.zeros(bins)
            for x in dx:
                if maxx == 0.0:
                    idx = 0
                else:
                    idx = int(np.floor(bins * x / maxx))
                if idx == bins: idx = bins - 1
                cnts[idx] += 1
            cnts = cnts / float(cnts.sum())
            hist.append(cnts)
            if maxy is None or cnts.max() > maxy: maxy = cnts.max()
        if maxx is None:
            maxx, maxy = 1.0, 1.0
        lpt = PGPlot.blank_plot(maxx, maxy)
        wid = maxx / bins
        stkc = HexColor("#111111")
        stkc.alpha = alpha
        for h, c in zip(hist, colors):
            c.alpha = alpha
            for ii, i in enumerate(h):
                if i > 0:
                    rc = Rect(ii * wid, 0, wid, i, 
                        strokeColor=stkc, strokeWidth=1.0, fillColor=c)
                    lpt.g.add(rc)
        return lpt
    
    @staticmethod
    def prop_plot(data, xmin=False, ymin=False, 
        xlabel="Temperatrure / K", ylabel="Heat Capacity / a.u."):
        lpt, _, _ = PGPlot.line_plot([ [tuple(d) for d in data ]], 
            [ HexColor(source_color("trans")) ], [ 1.0 ], 
            max([d[0] for d in data]), max([d[1] for d in data]), 
            min([d[0] for d in data]) if xmin else 0.0, 
            min([d[1] for d in data]) if ymin else 0.0)
        drw = Drawing(width=500, height=220)
        lpt.x, lpt.y = 100, 15
        lpt.width, lpt.height = 300, 180
        drw.add(lpt)
        drw.add(String(410, 2, xlabel, fontName="Arial", fontSize=12, 
            fillColor=HexColor("#151515")))
        drw.add(String(82, drw.height - 13, ylabel, 
            fontName="Arial", fontSize=12, fillColor=HexColor("#151515")))
        return drw

    @staticmethod
    def conv_plot(conv_data, conv_maxs):
        data = []
        colors = []
        widths = []
        names = []
        for k, u in sorted(conv_data.items(), key=lambda x: int(x[0])):
            data.append([(float(d[0]), float(d[1][0])) for d in u])
            data.append([(float(d[0]), float(d[1][1])) for d in u])
            colors += [HexColor(multi_color(int(k)))] * 2
            widths += [2.0, 0.2]
            names += [[ multi_name(int(k)), HexColor(multi_color(int(k)))]]
        lpt, _, _ = PGPlot.line_plot(data, colors, widths, conv_maxs[0], conv_maxs[1])
        drw = Drawing(width=500, height=220)
        lpt.x, lpt.y = 100, 15
        lpt.width, lpt.height = 300, 180
        drw.add(lpt)
        sh = drw.height - 35
        dsh = 10
        sx = [ 10, 30 ]
        for na, cl in names:
            drw.add(Line(lpt.x + sx[0], sh, lpt.x + sx[1], sh, 
                strokeColor=cl, strokeWidth=2.0))
            drw.add(String(lpt.x + sx[1] + 3, sh - 3, na, 
                fontName="Arial", fontSize=10, fillColor=cl))
            sh -= dsh
        drw.add(String(410, 2, "Total # of Config.", fontName="Arial", fontSize=12, 
            fillColor=HexColor("#151515")))
        drw.add(String(82, drw.height - 13, "# of Converged Config.", 
            fontName="Arial", fontSize=12, fillColor=HexColor("#151515")))
        return drw
