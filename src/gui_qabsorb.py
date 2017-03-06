# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'qabsorb.ui'
#
# Created by: PyQt5 UI code generator 5.7.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_QAbsorb(object):
    def setupUi(self, QAbsorb):
        QAbsorb.setObjectName("QAbsorb")
        QAbsorb.resize(924, 573)
        self.gridLayout_3 = QtWidgets.QGridLayout(QAbsorb)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.tableView = QtWidgets.QTableView(QAbsorb)
        self.tableView.setObjectName("tableView")
        self.gridLayout_3.addWidget(self.tableView, 2, 0, 1, 5)
        self.aboutButton = QtWidgets.QPushButton(QAbsorb)
        self.aboutButton.setMinimumSize(QtCore.QSize(110, 32))
        self.aboutButton.setMaximumSize(QtCore.QSize(110, 32))
        self.aboutButton.setObjectName("aboutButton")
        self.gridLayout_3.addWidget(self.aboutButton, 3, 3, 1, 1)
        self.quitButton = QtWidgets.QPushButton(QAbsorb)
        self.quitButton.setMinimumSize(QtCore.QSize(110, 32))
        self.quitButton.setMaximumSize(QtCore.QSize(110, 32))
        self.quitButton.setObjectName("quitButton")
        self.gridLayout_3.addWidget(self.quitButton, 3, 4, 1, 1)
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.openSMIButton = QtWidgets.QPushButton(QAbsorb)
        self.openSMIButton.setObjectName("openSMIButton")
        self.horizontalLayout.addWidget(self.openSMIButton)
        self.lineEdit = QtWidgets.QLineEdit(QAbsorb)
        self.lineEdit.setObjectName("lineEdit")
        self.horizontalLayout.addWidget(self.lineEdit)
        self.gridLayout_2.addLayout(self.horizontalLayout, 0, 0, 1, 1)
        self.label = QtWidgets.QLabel(QAbsorb)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 0, 1, 2, 1)
        self.groupBox = QtWidgets.QGroupBox(QAbsorb)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName("gridLayout")
        self.pampamodel = QtWidgets.QCheckBox(self.groupBox)
        self.pampamodel.setChecked(True)
        self.pampamodel.setObjectName("pampamodel")
        self.gridLayout.addWidget(self.pampamodel, 0, 0, 1, 1)
        self.gridLayout_2.addWidget(self.groupBox, 1, 0, 1, 1)
        self.gridLayout_3.addLayout(self.gridLayout_2, 0, 0, 1, 5)
        self.calculateButton = QtWidgets.QPushButton(QAbsorb)
        self.calculateButton.setMaximumSize(QtCore.QSize(97, 32))
        self.calculateButton.setObjectName("calculateButton")
        self.gridLayout_3.addWidget(self.calculateButton, 1, 0, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_3.addItem(spacerItem, 1, 2, 1, 3)
        spacerItem1 = QtWidgets.QSpacerItem(383, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_3.addItem(spacerItem1, 3, 0, 1, 3)
        self.progressBar = QtWidgets.QProgressBar(QAbsorb)
        self.progressBar.setProperty("value", 24)
        self.progressBar.setObjectName("progressBar")
        self.gridLayout_3.addWidget(self.progressBar, 1, 1, 1, 1)

        self.retranslateUi(QAbsorb)
        QtCore.QMetaObject.connectSlotsByName(QAbsorb)

    def retranslateUi(self, QAbsorb):
        _translate = QtCore.QCoreApplication.translate
        QAbsorb.setWindowTitle(_translate("QAbsorb", "QAbsorb"))
        self.aboutButton.setText(_translate("QAbsorb", "About"))
        self.quitButton.setText(_translate("QAbsorb", "Quit"))
        self.openSMIButton.setText(_translate("QAbsorb", "Open SMILES list"))
        self.label.setText(_translate("QAbsorb", "<html><head/><body><p><img src=\":/icons/qabsorb.128x128.png\"/></p></body></html>"))
        self.groupBox.setTitle(_translate("QAbsorb", "Models"))
        self.pampamodel.setText(_translate("QAbsorb", "HDM-PAMPA (GIT)"))
        self.calculateButton.setText(_translate("QAbsorb", "Calculate"))

import icons_rc
