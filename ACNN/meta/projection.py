
import numpy as np
from reportlab.lib.colors import Color, HexColor
from reportlab.platypus.flowables import Flowable
from cluster.coval import Covalent, AtomicColor
from meta.ana_utils.plot_charge import ChargePlotter, print_charge

# AtomicColor.x['B'] = '#20A4F3'
# AtomicColor.x['Pt'] = '#515052'

def rotate(cla, ang, x, y, z):
    ang = ang / 180.0 * np.pi
    norm = np.linalg.norm(np.array([x, y, z]))
    x, y, z = np.array([x, y, z]) / norm
    x = np.array([[np.cos(ang) + x**2 * (1 - np.cos(ang)),
                   x * y * (1 - np.cos(ang)) - z * np.sin(ang),
                   x * z * (1 - np.cos(ang)) + y * np.sin(ang)],
                  [y * x * (1 - np.cos(ang)) + z * np.sin(ang),
                   np.cos(ang) + y**2 * (1 - np.cos(ang)),
                   y * z * (1 - np.cos(ang)) - x * np.sin(ang)],
                  [z * x * (1 - np.cos(ang)) - y * np.sin(ang),
                   z * y * (1 - np.cos(ang)) + x * np.sin(ang),
                   np.cos(ang) + z**2 * (1 - np.cos(ang))]])
    return np.dot(cla, x.T)

def hex_color(c):
    cc = HexColor(c)
    return np.array([cc.red, cc.green, cc.blue]) * 255.0

class ProjectedShape(object):
    write = np.array([255.0, 255.0, 255.0])
    line_color = np.array([50.0, 50.0, 50.0])
    dot_color = np.array([100.0, 115.0, 155.0])
    atom_radius = 1.3
    line_radius = 0.3
    atom_glow = np.array([atom_radius * 0.35, -atom_radius * 0.35])
    contrast = 0.5
    def __init__(self, atom=True, elem=None, coords=None, ends=None, colors=None,
                 islight=False, isdark=False):
        self.atom = atom
        self.islight = islight
        self.isdark = isdark
        if self.atom:
            self.z = coords[0, 2]
            self.coord = coords[0]
            self.elem = elem
            self.color = hex_color(AtomicColor.x[elem])
            self.dcolor = hex_color(AtomicColor.x[elem])
        else:
            self.coords = coords
            self.ends = ends
            self.z = coords[:, 2].mean()
            self.colors = colors
            self.color = colors.mean(axis=0)
            self.dcolor = self.dot_color

    def calculate_colorx(self, zmax, zmin):
        x = 1.0 - (self.z - zmin) / (zmax - zmin)
        self.colorx = self.color * (1 - x) + ProjectedShape.write * x
        self.linex = ProjectedShape.line_color * (1 - x) + ProjectedShape.write * x
        self.dotx = self.dcolor * (1 - x) + ProjectedShape.write * x

    def divide(self, ind=20):
        rr = []
        if not self.atom:
            if self.coords[0, 2] == self.coords[1, 2]:
                self.ends = [True, True]
                return [self]
            cor = self.coords[1] - self.coords[0]
            col = self.colors[1] - self.colors[0]
            for i in range(ind):
                cor1 = cor * i / ind + self.coords[0]
                cor2 = cor * (i + 1) / ind + self.coords[0]
                col1 = col * i / ind + self.colors[0]
                col2 = col * (i + 1) / ind + self.colors[0]
                ends = [i == 0, i == ind - 1]
                rr.append(ProjectedShape(atom=False, coords=np.array([cor1, cor2]),
                                         ends=ends, colors=np.array([col1, col2]),
                                         isdark=self.isdark, islight=self.islight))
        return rr

    def draw_atom_simple(self, c):
        lstr, fillc, strc = 0, self.dotx, self.dotx * self.contrast
        if self.isdark:
            lstr = 1
        elif self.islight:
            fillc = self.dotx * self.contrast + ProjectedShape.write * (1 - self.contrast)
        c.setFillColor(Color(*list(fillc / 255.0)))
        c.setStrokeColor(Color(*list(strc / 255.0)))
        c.circle(self.coord[0], self.coord[1], ProjectedShape.atom_radius * 0.5,
                 stroke=lstr, fill=1)

    def draw_atom(self, c):
        fillc, strc = self.colorx, self.linex
        if self.islight:
            fillc = self.colorx * self.contrast + ProjectedShape.write * (1 - self.contrast)
            strc = self.linex * self.contrast + ProjectedShape.write * (1 - self.contrast)
        c.setStrokeColor(Color(*list(strc / 255.0)))
        c.setFillColor(Color(*list(fillc / 255.0)))
        glow = ProjectedShape.atom_glow
        c.circle(self.coord[0], self.coord[1], ProjectedShape.atom_radius,
                 stroke=1, fill=1)
        c.setFillColor(Color(*list(fillc / 255.0 * 0.5 +
                                   ProjectedShape.write / 255.0 * 0.5)))
        c.circle(self.coord[0] - glow[0], self.coord[1] - glow[1],
                 ProjectedShape.atom_radius * 0.2, stroke=0, fill=1)

    def draw_line_simple(self, c):
        strc = self.linex
        if self.islight:
            strc = self.linex * self.contrast + ProjectedShape.write * (1 - self.contrast)
        c.setStrokeColor(Color(*list(strc / 255.0)))
        c.setLineWidth(0.3)
        c.setLineCap(1)
        c.line(*(list(self.coords[0, :2]) + list(self.coords[1, :2])))

    def draw_line(self, c):
        fillc, strc = self.colorx, self.linex
        if self.islight:
            fillc = self.colorx * self.contrast + ProjectedShape.write * (1 - self.contrast)
            strc = self.linex * self.contrast + ProjectedShape.write * (1 - self.contrast)
        pv = self.coords[1, :2] - self.coords[0, :2]
        pg = np.dot([[np.cos(np.pi/2), np.sin(np.pi/2)],
                     [-np.sin(np.pi/2), np.cos(np.pi/2)]], pv)
        pg = pg / np.linalg.norm(pg) * ProjectedShape.line_radius
        c.setStrokeColor(Color(*list(strc / 255.0)))
        c.setFillColor(Color(*list(fillc / 255.0)))
        p = c.beginPath()
        p.moveTo(*list(self.coords[0, :2] + pg))
        p.lineTo(*list(self.coords[1, :2] + pg))
        p.lineTo(*list(self.coords[1, :2] - pg))
        p.lineTo(*list(self.coords[0, :2] - pg))
        p.lineTo(*list(self.coords[0, :2] + pg))
        c.drawPath(p, stroke=0, fill=1)
        c.line(*(list(self.coords[0, :2] + pg) + list(self.coords[1, :2] + pg)))
        c.line(*(list(self.coords[0, :2] - pg) + list(self.coords[1, :2] - pg)))
        lr = ProjectedShape.line_radius
        pv = pv / np.linalg.norm(pv) * lr
        for i, e in enumerate(self.ends):
            if e:
                if (self.coords[i, 2] > self.coords[1 - i, 2]) == (i == 0):
                    pvx = pv
                else:
                    pvx = -pv
                if self.coords[i, 2] < self.coords[1 - i, 2]:
                    p = c.beginPath()
                    p.moveTo(*list(self.coords[i, :2] + pg))
                    p.curveTo(*(list(self.coords[i, :2] + pg + pvx)
                                + list(self.coords[i, :2] - pg + pvx)
                                + list(self.coords[i, :2] - pg)))
                    p.lineTo(*list(self.coords[i, :2] + pg))
                    c.drawPath(p, stroke=0, fill=1)
                if self.coords[i, 2] < self.coords[1 - i, 2]:
                    p = c.beginPath()
                    p.moveTo(*list(self.coords[i, :2] + pg))
                    p.curveTo(*(list(self.coords[i, :2] + pg + pvx)
                                + list(self.coords[i, :2] - pg + pvx)
                                + list(self.coords[i, :2] - pg)))
                    c.drawPath(p, stroke=1, fill=0)

    def draw(self, c, simple):
        if not simple:
            if self.atom:
                self.draw_atom(c)
            else:
                self.draw_line(c)
        else:
            if self.atom:
                self.draw_atom_simple(c)
            else:
                self.draw_line_simple(c)

class ChargeShape(object):
    plus = np.array([200.0, 100.0, 100.0])
    minus = np.array([100.0, 100.0, 200.0])
    fontsize = 1.8
    def __init__(self, apos, cpos, charge, wid, ele, z, islight):
        self.apos, self.cpos, self.charge, self.wid = apos, cpos, charge, wid
        self.ele, self.z, self.islight = ele, z, islight
        self.atom = True
        self.cc = np.array(self.plus) if self.charge > 0 else np.array(self.minus)

    def calculate_colorx(self, zmax, zmin):
        x = 1.0 - (self.z - zmin) / (zmax - zmin)
        # if x > 0.3: x = 0.3
        self.cc = self.cc * (1 - x) + ProjectedShape.write * x

    def shift(self, d):
        self.apos += d
        self.cpos += d

    def draw(self, c, simple):
        assert not simple
        cc = self.cc
        if self.islight:
            cc = self.cc * ProjectedShape.contrast + \
                ProjectedShape.write * (1 - ProjectedShape.contrast)
        c.setStrokeColor(Color(*list(cc / 255.0)))
        c.line(*(list(self.apos) + list(self.cpos)))
        c.line(*(list(self.cpos) + [self.cpos[0] + self.wid, self.cpos[1]]))
        c.setFillColor(Color(*list(cc / 255.0)))
        c.setFont("Arial", self.fontsize)
        xpos, ypos = self.cpos
        text = print_charge(self.charge) + "|%s" % self.ele
        if self.wid < 0:
            xpos += self.wid
            text = "%s|" % self.ele + print_charge(self.charge)
        if self.cpos[1] < self.apos[1]:
            ypos -= self.fontsize
        c.drawString(xpos, ypos + 0.2, text)
        c.setFillColor(Color(*list(cc / 255.0)))
        c.circle(self.apos[0], self.apos[1],
                 ProjectedShape.atom_radius * 0.2, stroke=0, fill=1)

class Projection(object):
    def __init__(self, f, n=0, aspect=0, fov=0, perspective=False):
        self.x = np.eye(4)
        if not perspective:
            self.x[2, 2] = 1.0 / f
        else:
            fov = fov / 180.0 * np.pi
            self.x[0, 0] = 1.0 / (np.tan(fov * 0.5) * aspect)
            self.x[1, 1] = 1.0 / np.tan(fov * 0.5)
            self.x[2, 2] = f / (f - n)
            self.x[2, 3] = 1.0
            self.x[3, 2] = f * n / (n - f)
            self.x[3, 3] = 0.0

    def translate(self, d):
        x = np.eye(4)
        x[:, 3] = np.array([d[0], d[1], d[2], 1.0])
        self.x = np.dot(self.x, x)

    def scale(self, d):
        x = np.eye(4)
        for i in range(3):
            x[i, i] = d[i]
        self.x = np.dot(self.x, x)

    def rotate(self, ang, x, y, z):
        ang = ang / 180.0 * np.pi
        norm = np.linalg.norm(np.array([x, y, z]))
        x, y, z = np.array([x, y, z]) / norm
        x = np.array([[np.cos(ang) + x**2 * (1 - np.cos(ang)),
                       x * y * (1 - np.cos(ang)) - z * np.sin(ang),
                       x * z * (1 - np.cos(ang)) + y * np.sin(ang), 0],
                      [y * x * (1 - np.cos(ang)) + z * np.sin(ang),
                       np.cos(ang) + y**2 * (1 - np.cos(ang)),
                       y * z * (1 - np.cos(ang)) - x * np.sin(ang), 0],
                      [z * x * (1 - np.cos(ang)) - y * np.sin(ang),
                       z * y * (1 - np.cos(ang)) + x * np.sin(ang),
                       np.cos(ang) + z**2 * (1 - np.cos(ang)), 0],
                      [0, 0, 0, 1]])
        self.x = np.dot(self.x, x)

    @staticmethod
    def extend(vs):
        vs = np.array([[v[0], v[1], v[2], 1] for v in vs])
        return vs

    def make(self, vs):
        if len(vs.shape) == 1:
            vs.reshape(1, vs.shape[0])
        if vs.shape[1] == 3:
            vs = Projection.extend(vs)
        r = np.tensordot(self.x, vs, axes=[1, 1]).T
        r[:, 3] = np.array([(x if x != 0.0 else 1.0) for x in r[:, 3]])
        r = np.array([rx[:3] / np.abs(rx[3]) for rx in r])
        return r[:, :3]

class DrawCluster(Flowable, object):
    def __init__(self, clus, ratio=5.0, simple=False, rotmat=None,
                 surface_depth=3.0, ind=2, title="", clip=False, clip_circle=False,
                 plot_charge=True, zdepth=1.0, plot_force=False, force_factor=1.0,
                 force_image=3, perspective=True):
        if rotmat is None:
            rotmat = [0.0, 1.0, 0.0, 0.0]
        self.canv = None
        super(DrawCluster, self).__init__()
        p = Projection(fov=1.5, aspect=1.0, f=50.0, n=0.3, perspective=perspective)
        p.translate([0, 0, -200])
        p.scale([2.8] * 3)
        self.fontsize = 3.0 * ratio
        self.charge_fontsize = 1.8
        self.fhratio = 1.0
        self.p = p
        self.clus = clus
        if hasattr(clus, "surf"):
            self.clus = clus.to_cluster()
            if hasattr(clus, "props"):
                self.clus.props = clus.props
        self.ratio = ratio
        self.simple = simple
        self.rotmat = rotmat
        self.title = title
        self.surface_depth = surface_depth
        self.ind = ind
        self.clip = clip
        self.clip_circle = clip_circle
        self.plot_charge = plot_charge
        self.zdepth = zdepth
        self.draws = []
        self.shaped = None
        self.std_mean = None
        self.load_cluster()
        if plot_force:
            self.draws = []
            dx = np.array(range(1, force_image + 1))
            zp = list(dx * force_factor / force_image)
            zm = list(-dx * force_factor / force_image)[::-1]
            zz = zm + zp
            xatoms = np.array(self.clus.atoms)
            xcharge = self.plot_charge
            self.plot_charge = False
            for z in zz:
                self.clus.atoms = xatoms + z * self.clus.forces
                self.load_cluster(def_depth=-1)
            self.clus.atoms = xatoms
            self.plot_charge = xcharge
            self.load_cluster()

    def load_cluster(self, def_depth=0):
        pt = self.clus.atoms
        pte = self.clus.elems
        ptis = [def_depth] * self.clus.n # color depth for atoms
        new_surfnum = self.clus.surfnum
        ptchar = [None] * self.clus.n # plotting charges
        if hasattr(self.clus, "props") and "charges" in self.clus.props and self.plot_charge:
            ptchar = list(self.clus.props["charges"])
        if self.surface_depth is not None and self.clus.surfnum != 0:
            ptx = []
            ptex = []
            ptis = []
            ptcharx = []
            high_z = pt[0:self.clus.surfnum].max(axis=0)[2]
            for i in range(0, self.clus.surfnum):
                if pt[i, 2] >= high_z - self.surface_depth:
                    ptx.append(pt[i])
                    ptex.append(pte[i])
                    ptis.append(-1)
                    ptcharx.append(ptchar[i])
            new_surfnum = len(ptx)
            ptx += list(pt[self.clus.surfnum:])
            ptex += list(pte[self.clus.surfnum:])
            ptis += [1] * (self.clus.n - self.clus.surfnum)
            ptcharx += list(ptchar[self.clus.surfnum:])
            pt, pte, ptis = np.array(ptx), np.array(ptex), np.array(ptis)
            ptchar = ptcharx
        if self.rotmat == "best":
            import cluster.align, cluster.base
            ali = cluster.align.Align()
            ali.clus = [cluster.base.Cluster(n=len(pt), atoms=pt, elems=pte)]
            ali.clus[0].surfnum = len(pt) - self.clus.n + self.clus.surfnum
            ali.align()
            pt = ali.clus[0].atoms
            pt = rotate(pt, *[80.0, -1.0, 0.0, 0.0])
        else:
            pt = rotate(pt, *self.rotmat)
        if self.std_mean is None:
            self.std_mean = pt.mean(axis=0).reshape(1, 3)
        pt = pt - self.std_mean
        xr = pt[:, 0].max() - pt[:, 0].min()
        yr = pt[:, 1].max() - pt[:, 1].min()
        shapex = self.p.make(np.array([[0, 0, -max(xr, yr) / 2],
                                       [0, 0, max(xr, yr) / 2]]))
        zd = np.abs(shapex[0, 2] - shapex[1, 2])
        shape = self.p.make(pt)
        ptg = []
        ptgis = [] # color depth for lines
        if not self.simple:
            rdsl = ProjectedShape.atom_radius / 1.8 / 10.0
            rdsr = ProjectedShape.atom_radius / 1.8 / 3.6
        else:
            rdsl = ProjectedShape.atom_radius / 1.8 / 1.8
            rdsr = ProjectedShape.atom_radius / 1.8 / 1.8
        ptchb = [False] * new_surfnum
        for i in range(0, len(pt)):
            for j in range(i):
                ll = np.linalg.norm(pt[j] - pt[i])
                xlen = Covalent.x[pte[i]] + Covalent.x[pte[j]]
                if ll < rdsl + rdsr or ll > xlen + 0.45:
                    continue
                svl = (pt[j] - pt[i]) / ll * rdsl
                svr = (pt[j] - pt[i]) / ll * rdsr
                if shape[i, 2] > shape[j, 2]:
                    ptg += [pt[i] + svl, pt[j] - svr]
                else:
                    ptg += [pt[i] + svr, pt[j] - svl]
                if ptis[i] == 1 and ptis[j] == 1:
                    ptgis.append(1)
                elif ptis[i] == -1 or ptis[j] == -1:
                    ptgis.append(-1)
                else:
                    ptgis.append(0)
                if i >= new_surfnum and j < new_surfnum:
                    ptchb[j] = True
                elif j >= new_surfnum and i < new_surfnum:
                    ptchb[i] = True
        for i in range(0, new_surfnum):
            if not ptchb[i]:
                ptchar[i] = None
        ptg = np.array(ptg)
        shapel = self.p.make(ptg).reshape((ptg.shape[0] / 2, 2, ptg.shape[1])) if len(ptg) != 0 else np.zeros((0, ))
        cdraws = []
        if not self.simple:
            cp = ChargePlotter(shape, ptchar, self.charge_fontsize, pte)
            for ic, cpx in enumerate(cp.cpts):
                if cpx is not None:
                    cdraws.append(ChargeShape(apos=cpx[2], cpos=cpx[1], charge=cpx[0],
                                              wid=cpx[3], ele=pte[ic], z=shape[ic][2],
                                              islight=(ptis[ic] == -1)))
        rd = ProjectedShape.atom_radius * 1.5
        mx, my = shape[:, 0].min(), shape[:, 1].min()
        nx, ny = shape[:, 0].max(), shape[:, 1].max()
        if len(cdraws) != 0:
            hei = self.charge_fontsize
            mx = min(mx, np.array([[cd.cpos[0], cd.cpos[0] + cd.wid] for cd in cdraws]).min())
            nx = max(nx, np.array([[cd.cpos[0], cd.cpos[0] + cd.wid] for cd in cdraws]).max())
            my = min(my, np.array([[cd.cpos[1], cd.cpos[1] - hei] for cd in cdraws]).min())
            ny = max(ny, np.array([[cd.cpos[1], cd.cpos[1] + hei] for cd in cdraws]).max())
        cmx, cmy = mx, my
        cnx, cny = nx, ny
        if new_surfnum != len(shape):
            cmx, cmy = shape[new_surfnum:, 0].min(), shape[new_surfnum:, 1].min()
            cnx, cny = shape[new_surfnum:, 0].max(), shape[new_surfnum:, 1].max()
        if self.clip:
            mx, my = cmx, cmy
            nx, ny = cnx, cny
        if self.clip_circle:
            ccave = shape[new_surfnum:].mean(axis=0)
            ccr = np.linalg.norm(shape[new_surfnum:] - ccave, axis=1).max()
            mx, my = ccave[0] - ccr, ccave[1] - ccr
            nx, ny = ccave[0] + ccr, ccave[1] + ccr
        dx, dy = nx - mx, ny - my
        if self.shaped is None:
            shaped = -np.array([[mx - rd, my - rd, 0]])
            self.shaped = shaped
        else:
            shaped = np.array(self.shaped)
        self.width = self.ratio * (dx + rd * 2)
        self.height = self.ratio * (dy + rd * 2)
        shape += shaped
        if len(shapel) != 0:
            shapel += shaped
        for cd in cdraws:
            cd.shift(shaped[0,:2])
        self.cwh = np.array([[cmx - rd, cmy - rd], [cnx + rd, cny + rd]]) + shaped[:, :2]
        self.cwh *= self.ratio
        if len(self.title) != 0:
            self.height += self.fontsize / self.fhratio
        self.shs = shape[new_surfnum:] * self.ratio

        gray = np.array([220, 220, 220])
        draws = []
        for ikx, s in enumerate(shape):
            x = ProjectedShape(atom=True, elem=pte[ikx], isdark=(ptis[ikx] == 1),
                               islight=(ptis[ikx] == -1), coords=np.array([s]))
            draws.append(x)
        for iks, sl in enumerate(shapel):
            x = ProjectedShape(atom=False, coords=np.array([sl[0], sl[1]]),
                               colors=np.array([gray, gray]),
                               isdark=(ptgis[iks] == 1), islight=(ptgis[iks] == -1))
            draws += x.divide(ind=self.ind)
        draws += cdraws

        draws.sort(key=lambda x: [x.z, 1 if x.atom else 0])
        zs = [x.z for x in draws]
        zmax, zmin = max(zs), min(zs)
        if not self.simple:
            if zmax - zmin < zd * 0.8:
                zmin = zmax - zd * 0.8
            zmax, zmin = zmax + (zmax - zmin) * 0.3, zmin - (zmax - zmin) * 0.5
        else:
            if zmax - zmin < zd * 0.5:
                zmin = zmax - zd * 0.5
            zmax, zmin = zmax + (zmax - zmin) * 0.2, zmin - (zmax - zmin) * 0.3
        zmin = zmax - self.zdepth * (zmax - zmin)
        for d in draws:
            d.calculate_colorx(zmax, zmin)
        self.draws += draws

    def wrap(self, *_):
        return (self.width, self.height)

    def draw(self):
        self.canv.setLineWidth(0.1)
        # self.canv.setStrokeColor("#3E3E3E")
        # self.canv.rect(0, 0, self.width, self.height, fill=0)
        if self.clip:
            self.canv.setLineWidth(0.5)
            self.canv.setStrokeColor("#333333")
            p = self.canv.beginPath()
            p.moveTo(0, 0)
            p.lineTo(0, self.height)
            p.lineTo(self.width, self.height)
            p.lineTo(self.width, 0)
            p.lineTo(0, 0)
            self.canv.clipPath(p, stroke=0, fill=0)
        elif self.clip_circle:
            self.canv.setLineWidth(1.0)
            self.canv.setFillColor("#FFFFFF")
            p = self.canv.beginPath()
            p.moveTo(self.width, self.height / 2)
            p.arcTo(0, 0, self.width, self.height, startAng=0, extent=360)
            self.canv.clipPath(p, stroke=0, fill=1)
        if len(self.title) != 0:
            self.canv.setFillColor("#8E8E8E")
            self.canv.setFont("Arial", self.fontsize)
            self.canv.drawCentredString(self.width / 2, self.height - self.fontsize, self.title)
        self.canv.scale(self.ratio, self.ratio)
        for d in self.draws:
            d.draw(self.canv, self.simple)
