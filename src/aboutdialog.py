'''
@package aboutdialog.py

aboutdialog.py was writen by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
and is distributed under LGPL version 3

Geneve February 2015
'''

from PyQt5 import *

from gui_aboutdialog import Ui_AboutDialog

class AboutDialog(QtWidgets.QDialog, Ui_AboutDialog):
    def __init__(self, parent=None):
        QtWidgets.QDialog.__init__(self,parent)
        self.setupUi(self)
        self.closeButton.clicked.connect(self.close_)

    def close_(self):
        self.reject()
