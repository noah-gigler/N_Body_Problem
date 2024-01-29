# -*- coding: utf-8 -*-
"""
    Animated 3D sinc function
"""
from PyQt5 import QtWidgets 
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import pyqtgraph as pg
import numpy as np
import sys

def read_xyz(filename):
    # line format x y z vx vy vz
    with open(filename, 'r') as f:
        lines = f.readlines()
        lines = [line.split() for line in lines]
        lines = np.array(lines)
        lines = lines.astype(np.float64)
        return lines[:, 0], lines[:, 1], lines[:, 2]

class Visualizer(object):
    def __init__(self):
        self.app = QtWidgets.QApplication(sys.argv)
        self.w = gl.GLViewWidget()
        self.w.opts['distance'] = 40
        self.w.setWindowTitle('pyqtgraph example: GLLinePlotItem')
        self.w.setGeometry(0, 110, 1920, 1080)
        self.w.show()

        # create the background grids
        gx = gl.GLGridItem()
        gx.rotate(90, 0, 1, 0)
        gx.translate(-10, 0, 0)
        self.w.addItem(gx)
        gy = gl.GLGridItem()
        gy.rotate(90, 1, 0, 0)
        gy.translate(0, -10, 0)
        self.w.addItem(gy)
        gz = gl.GLGridItem()
        gz.translate(0, 0, -10)
        self.w.addItem(gz)

        self.i = 0

        self.x, self.y, self.z = read_xyz(path + 'frame{}.txt'.format(self.i))

        self.plotter = gl.GLScatterPlotItem(pos=np.array([self.x, self.y, self.z]).transpose(), color=pg.glColor('w'),
                                          size=0.1, pxMode=True)
        self.w.addItem(self.plotter)

    def start(self):
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QtWidgets.QApplication.instance().exec_()

    def set_plotdata(self, points, color):
        self.plotter.setData(pos=points, color=color)

    def update(self):
        self.i += 1
        if self.i == 60:
            self.i = 0
        self.x, self.y, self.z = read_xyz(path + 'xv_{}.txt'.format(self.i))

        self.set_plotdata(points=np.array([self.x, self.y, self.z]).transpose(), color=pg.glColor('w'))
        self.w.grabFramebuffer().save('frames/img_{}.png'.format(self.i))   

    def animation(self):
        timer = QtCore.QTimer()
        timer.timeout.connect(self.update)
        timer.start(5)
        self.start()

path = 'output/sim_basic/'

# Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    v = Visualizer()
    v.animation()
