'''
@package main.py

main.py was writen by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
and is distributed under LGPL version 3

Geneve November 2016
'''

import os
import shutil
import sys
import threading
from PyQt5 import *

import gui_qabsorb as qabs
#from icons_rc import *
from molutils import *
from tabmodel import *
from aboutdialog import *

# pyinstaller bug under windows
#import Tkinter
#import FileDialog

mpath = None
#try:
#    mpath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
#except NameError:
#    mpath = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
#mpath = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
mpath = os.getcwd()+"/"
mpath += "../models"

#if not mpath in sys.path:
#    sys.path.insert(1, mpath)
#del mpath

class QAbsorb(QtWidgets.QWidget, qabs.Ui_QAbsorb):
    def __init__(self, parent=None):
        """ init method """
        super(QAbsorb, self).__init__(parent)
        self.setupUi(self)
        self.progressBar.hide()
        self.openSMIButton.clicked.connect(self.openSMIlst)
        self.quitButton.clicked.connect(self.close)
        self.aboutButton.clicked.connect(self.about)
        self.calculateButton.clicked.connect(self.computemodels)
        #self.pampamodel.stateChanged.connect(self.computemodels)

        self.tablemodel = TableModel(self)
        self.tableView.setModel(self.tablemodel)
        self.tableView.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.tableView.customContextMenuRequested.connect(self.openTableMenu)

        #self.lineEdit.setText("/Users/marco/tubastatinA.smi")


    def openSMIlst(self):
        """ open smiles list method """
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open smile list', '', 'SMILES files ( *.smi *.smile)')
        if fname and os.path.isfile(fname):
            self.lineEdit.setText(fname)

    def close(self):
        """ close method """
        QtWidgets.QApplication.quit()

    def about(self):
        """ about method """
        adialog = AboutDialog()
        adialog.exec_()
        return

    def openTableMenu(self, position):
        """ context menu event """
        menu = QMenu(self)
        exportAction = menu.addAction("Export table as CSV")
        action = menu.exec_(self.tableView.viewport().mapToGlobal(position))

        if action == exportAction:
            fname = QFileDialog.getSaveFileName(self, "Save File", "CSV (*.csv)");
            self.tablemodel.SaveTable(fname)
            return
        else:
            return

    def keyPressEvent(self, e):
        """ key press method """
        if (e.modifiers() & QtCore.Qt.ControlModifier):
            if e.key() == QtCore.Qt.Key_C: #copy
                if len(self.tableView.selectionModel().selectedIndexes()) > 0:
                    previous = self.tableView.selectionModel().selectedIndexes()[0]
                    columns = []
                    rows = []
                    for index in self.tableView.selectionModel().selectedIndexes():
                        if previous.column() != index.column():
                            columns.append(rows)
                            rows = []
                        rows.append(str(index.data().toPyObject()))
                        previous = index
                    columns.append(rows)

                    # add rows and columns to clipboard
                    clipboard = ""
                    nrows = len(columns[0])
                    ncols = len(columns)
                    for r in xrange(nrows):
                        for c in xrange(ncols):
                            clipboard += columns[c][r]
                            if c != (ncols-1):
                                clipboard += '\t'
                        clipboard += '\n'

                    # copy to the system clipboard
                    sys_clip = QApplication.clipboard()
                    sys_clip.setText(clipboard)


    def getdesc_(self, mtoget, mpath, tabres):
        GetDescValues(self.lineEdit.text(), mtoget, mpath, tabres)

    def computemodels(self):
        if len(self.lineEdit.text()) > 0:
            self.progressBar.show()
            self.progressBar.setMinimum(0)
            self.progressBar.setMaximum(0)
            del self.tablemodel.arraydata[:]
            del self.tablemodel.header[:]
            self.tablemodel.clean()
            self.tableView.model().layoutChanged.emit()
            header = ["Depiction", "Name"]
            mtoget = []
            if self.pampamodel.isChecked():
                header.append("HDM-PAMPA (GIT+: Permeable; GIT-: Not permeable)")
                mtoget.append("m3")
            self.tablemodel.setHeader(header)
            header.append("SMILE")

            tabres = []
            th = threading.Thread(target = self.getdesc_, args=(mtoget, mpath, tabres,))
            th.start()

            while True:
                if len(tabres) > 0 or th.isAlive() == False:
                    break
                else:
                    QtWidgets.QApplication.processEvents()


            #self.getdesc_(mtoget, mpath, tabres)

            for i in range(len(tabres)):
                row = [GenQPixmapDepiction(tabres[i][1]), tabres[i][0]]
                for j in range(2, len(tabres[i])):
                    try:
                      row.append(round(tabres[i][j], 2))
                    except:
                      row.append(tabres[i][j])
                row.append(tabres[i][1])
                self.tablemodel.addRow(row)

            self.tableView.model().layoutChanged.emit()
            self.tableView.resizeRowsToContents()
            self.tableView.resizeColumnsToContents()
            self.progressBar.hide()


def main():
    app = QtWidgets.QApplication(sys.argv)
    form = QAbsorb()
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()
